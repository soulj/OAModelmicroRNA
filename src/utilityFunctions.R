
#utility function to get results data frame
getResultsDataFrame <- function(dds, condition, numerator, 
                                denominator) {
   data <- as.data.frame(lfcShrink(dds, contrast = c(condition,numerator, denominator),type = "normal",BPPARAM = MulticoreParam(3)))
  
  data <- data[, c("log2FoldChange", "padj","lfcSE")]
  #colnames(data) <- paste(paste(numerator, denominator,sep = "vs"), colnames(data), sep = "_")
  
  colnames(data) <- c("logFC","adj.P.Val","lfcSE")
  return(data)
}

#have to specify the host to a particular release to work around a biomart problem
#https://support.bioconductor.org/p/9143874/
convertHumanGeneList <- function(x){
  human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",version=103)
  mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl",version=103)
  genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol",values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=TRUE)
  return(genesV2)
}


annotateTables <- function(diffExpTable,en2gene,species="Human"){
  diffExpTable$ID <- rownames(diffExpTable)
  diffExpTable <- na.omit(merge(diffExpTable, en2gene, by.x = "ID", 
                                by.y = "gene_id"))
  
  if(species=="Mouse"){
  #get a table of the human and gene symbols
  human2geneTable <- convertHumanGeneList(diffExpTable$gene_name)
  
  #merge the annotation table to this table
  diffExpTable <- merge(diffExpTable,human2geneTable,by.x="gene_name",by.y="MGI.symbol")
  }
  
  return(diffExpTable)
  
  
}


plotVolcano <- function(diffExpTable,en2gene=NULL,fcThreshold=1.5,nDown=7,nUp=16){
  
  diffExpTable$ID <- rownames(diffExpTable)
  diffExpTable <- na.omit(merge(diffExpTable, en2gene, by.x = "ID", 
                                by.y = "gene_id"))
  
  diffExpTable$logpval <- -log10(diffExpTable$adj.P.Val)
  diffExpTable$test <- ifelse(diffExpTable$adj.P.Val<=0.05 & abs(diffExpTable$logFC)>=log2(fcThreshold),TRUE,FALSE)
  
  
  diffExpTable$Pi <- diffExpTable$logFC * -log10(diffExpTable$adj.P.Val)
  chosen <- diffExpTable[order(abs(diffExpTable$Pi),decreasing = TRUE),]
  chosen <- chosen[!grepl("Gm|E03",chosen$gene_name),]
  chosen <- chosen[chosen$gene_biotype=="protein_coding" & chosen$adj.P.Val<=0.05,]
  chosen.up <- chosen[ chosen$logFC>=log2(fcThreshold),"gene_name"]
  chosen.down <- chosen[ chosen$logFC <= log2(1/fcThreshold),"gene_name"]
  chosen <- c(chosen.up[1:nUp],chosen.down[1:nDown])
  
  g <- ggplot(diffExpTable, aes(x=logFC, y=logpval)) +
    geom_point(aes(colour=test), size=3, alpha=1) +
    scale_fill_carto_d( palette = "Safe") +
    geom_vline(xintercept=log2(fcThreshold), colour='blue',linetype=2) +
    geom_vline(xintercept=-log2(fcThreshold), colour='blue',linetype=2) +
    geom_hline(yintercept=-log10(0.05), colour='blue',linetype=2) +
    cowplot::theme_cowplot(font_size = 28) +
    xlab("log2 fold change") + ylab("-log10 adj p-value") + theme(legend.position = "none") +
    ggrepel::geom_text_repel(data = subset(diffExpTable,gene_name %in% chosen), aes(label=gene_name),size=7,max.overlaps = Inf,force = 2,min.segment.length = 0,seed = 42) +
    theme(axis.line = element_line(size = 1))

  
  return(g)
}

plotOATargets <- function(diffExpTable,en2gene,OATargets,fcThreshold){
  
  diffExpTable$ID <- rownames(diffExpTable)
  diffExpTable <- na.omit(merge(diffExpTable, en2gene, by.x = "ID", 
                                by.y = "gene_id"))

  diffExpTable$logpval <- -log10(diffExpTable$adj.P.Val)
  
  diffExpTable <- merge(diffExpTable,OATargets,by.x="gene_name",by.y="MGI.symbol")
  
  chosen <- diffExpTable[order(abs(diffExpTable$logFC),decreasing = TRUE),]
  chosen <- chosen[chosen$adj.P.Val<=0.05,]
  chosen.up <- chosen[ chosen$logFC >= log2(fcThreshold) ,"gene_name"]
  chosen.down <- chosen[ chosen$logFC <= log2(1/fcThreshold),"gene_name"]
  chosen <- c(chosen.up[1:16],chosen.down[1:7])
  
  g <- ggplot(diffExpTable, aes(x=logFC, y=logpval)) +
    geom_point(aes(colour=effectConsensus), size=3, alpha=1) +
    scale_colour_manual(values = c("#88CCEE", "#CC6677" ,"#DDCC77")) +
    geom_vline(xintercept=log2(1.5), colour='blue',linetype=2) +
    geom_vline(xintercept=-log2(1.5), colour='blue',linetype=2) +
    geom_hline(yintercept=-log10(0.05), colour='blue',linetype=2) +
    cowplot::theme_cowplot(font_size = 28) +
    xlab("log2 fold change") + ylab("-log10 adj p-value") + theme(legend.position = "none") +
    ggrepel::geom_text_repel(data = subset(diffExpTable,gene_name %in% chosen), aes(label=gene_name),size=6,max.overlaps = Inf,force = 2,min.segment.length = 0,seed = 42) +
    theme(axis.line = element_line(size = 1))

  
  return(list(g,diffExpTable))
}



plotOATargetsMir <- function(diffExpTable,mirs){
  
  diffExpTable$ID <- sapply(strsplit(as.character(diffExpTable$ID),"_"),"[[",1)
  diffExpTable$ID <- gsub("\\.","-",diffExpTable$ID)
  diffExpTable$ID <- gsub("mmu-","",diffExpTable$ID)
  diffExpTable$mir <- gsub("-3p|-5p","",diffExpTable$ID)
  diffExpTable$mir <- tolower(diffExpTable$mir)
  diffExpTable$mir <- gsub("a|b|c|d","",diffExpTable$mir)

  diffExpTable <- merge(diffExpTable,mirs,by="mir")
  diffExpTable$logpval <- -log10(diffExpTable[,4])
  
  colnames(diffExpTable)[3] <- "log2FC"
  
  chosen <- diffExpTable[order(diffExpTable[,4],decreasing = FALSE),]
  chosen.up <- chosen[ chosen[,3]>0,"ID"]
  chosen.down <- chosen[ chosen[,3]<0,"ID"]
  chosen <- c(chosen.up[1:10],chosen.down[1])
  
  g <- ggplot(diffExpTable, aes(x=log2FC, y=logpval)) +
    geom_point(aes(colour=effectConsensus), size=3, alpha=1) +
    scale_colour_manual(values=c('grey', 'red','blue')) +
    geom_vline(xintercept=log2(1.5), colour='blue',linetype=2) +
    geom_vline(xintercept=-log2(1.5), colour='blue',linetype=2) +
    geom_hline(yintercept=-log10(0.05), colour='blue',linetype=2) +
    cowplot::theme_cowplot(font_size = 22) +
    xlab("log2 fold change") + ylab("-log10 adj p-value") + theme(legend.position = "none") +
    ggrepel::geom_text_repel(data = subset(diffExpTable,ID %in% chosen), aes(label=ID),max.overlaps = Inf,force = 2,min.segment.length = 0,seed = 42)
  
  return(g)
}


#get the cosine similarity to all the datasets
cos.sim <- function(i,X,query) {
  
  A = query[,grep("log|GeneSymbol|HGNC.symbol|gene_name",colnames(query))]
  B = X[,i,drop=F]
  merged <- merge(A,B,by.x=1,by.y="row.names")
  
  A=merged[,2]
  B=merged[,3]
  
  return( sum(A*B,na.rm = T)/sqrt(sum(A^2,na.rm = T)*sum(B^2,na.rm = T)) )
}


getSim <- function(dataset,foldchangeTable,accessions,datasetName="query"){
  sim <- sapply(seq_along(foldchangeTable)[-ncol(foldchangeTable)],cos.sim,foldchangeTable,dataset)
  sim <- data.frame(colnames(foldchangeTable)[-ncol(foldchangeTable)],sim)
  
  sim$dataset <- datasetName
  sim <- merge(sim,accessions,by.x=1,by.y="combined")
  sim$zscore <- scale(sim$sim)
  colnames(sim)[1] <- "ID"
  sim
  
}

tidySims <- function(results,name){
  
  r <- results[,c(1,7)]
  colnames(r) <- c("ID",name)
  return(r)
}


# function to get the results table for a limma contrast
getLimmaResultsDataFrame <- function(contrast, fit) {
  
  print(contrast)
  data <- topTable(fit, coef = contrast, number = Inf, 
                   sort.by = "P",confint = TRUE)
  #data <- independentFiltering(data, filter = data$AveExpr, objectType = "limma")
  #colnames(data) <- paste(paste(numerator, denominator, sep = "vs"), colnames(data), sep = "_")
  return(data)
}

runDiffExp <- function(txiData,sampleTable,comparisonsTable,OATargets,variables,method,block=NULL,PCAColumns,foldChangeTable,species="Human",saveCounts=FALSE,nDown=10,nUp=10,fcThreshold=1.5){
  
  #read in and prepare the gene counts for DESeq2
  txi <- readRDS(txiData)
  
  #make the the colData
  sampleTable <- read.delim(sampleTable)
  comparisonsTable <- read.delim(comparisonsTable)
  sampleTable <- sampleTable[, !grepl("File",colnames(sampleTable))]
  
  #make sure the sample metadata and order of the counts match
  txi$counts <- txi$counts[, match(as.character(sampleTable[,1]), colnames(txi$counts))]
  txi$abundance <- txi$abundance[, match(as.character(sampleTable[,1]), colnames(txi$abundance))]
  txi$length <- txi$length[, match(as.character(sampleTable[,1]), colnames(txi$length))]
  
  #specify the design formula
  designFormula <- as.formula(paste("~", paste(variables, collapse = " + ")))
  
  #make the deseq2 object
  dds <- DESeqDataSetFromTximport(txi = txi, colData = sampleTable, 
                                  design = designFormula)
  #run DESeq2
  dds <- DESeq(dds,parallel=TRUE,BPPARAM = MulticoreParam(4))
  
  #optionally save the counts matrix
  if(saveCounts){
    mrnaCounts <- counts(dds,normalized=TRUE)
    saveRDS(mrnaCounts,file="data/mir199Counts.RDS")
  }
  
  
  if(method=="limma"){
  
  designFormula <- as.formula(paste("~ 0 + ", paste(variables, collapse = " + ")))
  design  <- model.matrix(designFormula, colData(dds))
  
  #use limma voom for the replicated samples experiments
  countData <- DGEList(txi$counts)
  countData <- calcNormFactors( countData )
  
  #remove very lowly expressed genes
  keep <- filterByExpr(countData,design)
  countData <- countData[keep, ]
  
  #run voom the first time
  vobj <- voom(countData, design, plot=FALSE)
  
  if(!is.null(block)){
    
    dupcor <- duplicateCorrelation(vobj,design,block=colData(dds)[,block])
    
    # run voom again considering the duplicateCorrelation results
    vobj <- voom( countData, design, plot=FALSE, block=colData(dds)[,block], correlation=dupcor$consensus)
    
    # Estimate linear mixed model with a single variance component
    # Fit the model for each gene, 
    dupcor <- duplicateCorrelation(vobj, design, block=colData(dds)[,block])
    
    # This step uses only the genome-wide average for the random effect
    fit <- lmFit(vobj, design, block=colData(dds)[,block], correlation=dupcor$consensus)
    
    contrasts <- apply(comparisonsTable, 1, function(x) paste0(make.names(x[1]), 
                                                               make.names(x[2]), "-", make.names(x[1]), make.names(x[3])))
    contrast.matrix <- makeContrasts(contrasts = contrasts, 
                                     levels = colnames((coef(fit))))
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    resultsTable <- lapply(seq_along(contrasts), getLimmaResultsDataFrame,fit2)
    names(resultsTable) <- contrasts
    
  } else{
    fit <- lmFit(vobj, design)
    fit <- eBayes(fit)
    coefs <- paste0(comparisonsTable$Factor, comparisonsTable$Numerator)
    resultsTable <- lapply(coefs, getLimmaResultsDataFrame,fit)
    
    names(resultsTable) <- coefs
  }
  
  
  } else {
    
    #make one results table from all the comparisons
    resultsTable <- lapply(1:nrow(comparisonsTable), 
                           function(x) getResultsDataFrame(dds, comparisonsTable[x,1], comparisonsTable[x, 2], comparisonsTable[x,3]))
    
    
  }
  
  
  #map the ensembl gene to the gene symbols
  if(species=="Human"){
    endf <- GenomicFeatures::genes(EnsDb.Hsapiens.v103, return.type = "DataFrame")
  } else {
    endf <- GenomicFeatures::genes(EnsDb.Mmusculus.v103, return.type = "DataFrame")
  }
  en2gene <- as.data.frame(endf[, c("gene_id", "gene_name","gene_biotype")])  
  
  #make volcano plots for the comparisons
  volcanoPlots <- lapply(resultsTable,plotVolcano,en2gene,fcThreshold,nDown,nUp)
  
  #plot the differential expression of known OA regulators from OATargets
  if(!is.null(OATargets)){
  OATargetPlots <- lapply(resultsTable,plotOATargets,en2gene,OATargets,fcThreshold)
  } else {
    OATargetPlots <- NA
  }
  
  #annotate individual results tables for downstream analysis
  cutTables <- lapply(resultsTable,annotateTables,en2gene)
  
  #save one big table of results
  resultsTable <- do.call(cbind, resultsTable)
  resultsTable <- merge(resultsTable, en2gene, by.x = "row.names", 
                        by.y = "gene_id")
  geneNameInd <- which(colnames(resultsTable)=="gene_name")
  
  resultsTable <- resultsTable[, c(1,geneNameInd:ncol(resultsTable), 2:(geneNameInd-1))]
  colnames(resultsTable)[1] <- "ID"
  write.table(resultsTable, file = foldChangeTable, col.names = TRUE, 
              row.names = F, sep = "\t", quote = FALSE)
  
  #plot the pca of the data
  vsd <- varianceStabilizingTransformation(dds)
  g <- plotCustomPCA(vsd,PCAColumns) + theme_cowplot(font_size = 28)
  
  g <- g +  theme(axis.line = element_line(size = 1.2))
  
  return(list(volcano=volcanoPlots,pca=g,oatargets=OATargetPlots,cutTables=cutTables))
}

targetPlot <- function(diffExp,targetScan,TPMs){

#combine target scan with gene exp data
mir199Diff.targetscan <- merge(diffExp,targetScan,by.x="gene_name",by.y="Target.gene")

#add TPMs
mir199Diff.targetscan <- merge(mir199Diff.targetscan,TPMs,by.x="ID",by.y="ind")
mir199Diff.targetscan <- mir199Diff.targetscan[ mir199Diff.targetscan$values>=5,]
mir199Diff.targetscan <- mir199Diff.targetscan %>% group_by(gene_name) %>% dplyr::slice(which.min(Predicted.occupancy..low.miRNA.)) %>% as.data.frame()
g <- ggplot(mir199Diff.targetscan,aes(x=Predicted.occupancy..low.miRNA.,logFC)) +
  theme_cowplot(font_size = 21) + geom_point(colour="#EE8879FF")  + xlab("Targetscan score") + ylab("log2FC mir199a-5p inhibition vs Ctrl") +
ggrepel::geom_text_repel(data = subset(mir199Diff.targetscan, adj.P.Val <= 0.05 & logFC >= log2(1)), aes(label=gene_name,fontface ="bold"),size=6,color="black",max.overlaps = Inf,force = 35,min.segment.length = 0,seed = 42) +
 theme(legend.position = "none") + 
  theme(axis.line = element_line(size = 1))


return(list(targets=mir199Diff.targetscan,plot=g))
}



#function to plot the fold change by mean expression for the microRNA data
plotMA <- function(diffExpTable,meanCounts){
  
  
  diffExpTable <- na.omit(diffExpTable)
  
  #tidy up microrna name
  diffExpTable$ID <- sapply(strsplit(as.character(diffExpTable$ID),"_"),"[[",1)
  diffExpTable$ID <- gsub("\\.","-",diffExpTable$ID)
  diffExpTable$ID <- gsub("mmu-","",diffExpTable$ID)
  
  #prepare the labels for the plot
  diffExpTable$test <- ifelse(diffExpTable$adj.P.Val<=0.05,TRUE,FALSE)

  chosen <- diffExpTable[order(diffExpTable$meanCounts,decreasing = TRUE),]
  chosen <- chosen[ chosen$adj.P.Val<=0.05,]
  chosen.up <- chosen[ chosen[,2]>=log2(1.5),"ID"]
  chosen.down <- chosen[ chosen[,2]<=log2(1/1.5),"ID"]
  chosen <- c(chosen.up[1:8],chosen.down[1:5])
  
  chosen <- chosen[chosen != "miR-199a-3p" ]
  
  #make a composite label for mir199a/b -3p as sequence is the same
  diffExpTable$label <- diffExpTable$ID
  diffExpTable[diffExpTable$label == "miR-199b-3p","label" ] <- "miR-199a/b-3p"
  
  #make the plot
  ggplot(diffExpTable,aes(x=log2(meanCounts),y=logFC)) + geom_point(aes(color=test)) +
    scale_colour_manual(values=c('grey', 'red')) +
    geom_hline(yintercept=log2(1.5), colour='blue',linetype=2) +
    geom_hline(yintercept=-log2(1.5), colour='blue',linetype=2) +theme_cowplot(font_size = 28) +
    ylab("log2 fold change") + xlab("log2 Normalised Counts") +
    ggrepel::geom_text_repel(data = subset(diffExpTable,ID %in% chosen), aes(label=label),size=5,max.overlaps = Inf,force = 9,max.time=2,min.segment.length = 0,seed = 123,nudge_y = 0.2) +
    theme(legend.position = "none") + theme(axis.line = element_line(size = 1))

  
  
  
}

plotMirVolcano <- function(diffExpTable){
  
  #tidy the mir names for plotting
  diffExpTable$mir <- sapply(strsplit(diffExpTable[,1],"_"),"[[",1)
  diffExpTable$mir <- gsub("\\.","-",diffExpTable$mir)
  diffExpTable$mir <- gsub("mmu-","",diffExpTable$mir)
  
  #prepare the y axis and colour
  diffExpTable$logpval <- -log10(diffExpTable$adj.P.Val)
  diffExpTable$test <- ifelse(diffExpTable$adj.P.Val<=0.05,TRUE,FALSE)
  
  #get the mirs to label
  chosen <- diffExpTable[order(diffExpTable$adj.P.Val,decreasing = FALSE),]
  chosen.up <- chosen[ chosen[,2]>0,"mir"]
  chosen.down <- chosen[ chosen[,2]<0,"mir"]
  chosen <- c(chosen.up[1:8],chosen.down[1:8])
  colnames(diffExpTable)[2] <- "log2FC"
  
  diffExpTable$label=diffExpTable$mir
  
  chosen <- chosen[!chosen %in% c("miR-486b-3p","miR-486b-5p") ]
  diffExpTable[diffExpTable$label == "miR-486a-3p","label" ] <-"miR-486a/b-3p"
  diffExpTable[diffExpTable$label == "miR-486a-5p","label" ] <-"miR-486a/b-5p"
  
  #make the volcano plot
  g <- ggplot(diffExpTable, aes(x=log2FC, y=logpval)) +
    geom_point(aes(colour=test), size=3, alpha=1) +
    scale_colour_manual(values=c('grey', 'red','blue')) +
    geom_vline(xintercept=log2(1.5), colour='blue',linetype=2) +
    geom_vline(xintercept=-log2(1.5), colour='blue',linetype=2) +
    geom_hline(yintercept=-log10(0.05), colour='blue',linetype=2) +
    cowplot::theme_cowplot(font_size = 28) +
    xlab("log2 fold change") + ylab("-log10 adj p-value") + theme(legend.position = "none") +
    ggrepel::geom_text_repel(data = subset(diffExpTable,mir %in% chosen), aes(label=label),max.overlaps = Inf,size=5,force = 9,min.segment.length = 0,seed = 42) +
    theme(axis.line = element_line(size = 1))

  
  return(g)
}

#function to plot the PCA with samples grouped by two variables - based on the DESeq2 plotPCA function
plotCustomPCA <- function(object, intgroup, ntop = 30000) {
  #get the variable genes
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  
  #get the principle components and the variance explained by them
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  
  #build the table for plotting
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], intgroup.df, 
                  name = colnames(object))
  
  #plot
      if (length(intgroup) > 1) {
    color <- sym(colnames(d)[3])
    shape <- sym(colnames(d)[4])
    g <- ggplot(data = d, aes(x = PC1, y = PC2, 
                              color = !!color, shape = !!shape)) + geom_point(size = 5) + xlab(paste0("PC1: ", 
                                                                                                  round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] * 100), 
                  "% variance")) + coord_fixed() + scale_fill_carto_d( palette = "Safe") +  
      scale_shape_discrete(name = colnames(intgroup.df)[2])
      } else {
        
        color <- sym(colnames(d)[3])
        g <- ggplot(data = d, aes(x = PC1, y = PC2,
                                  color = !!color)) + geom_point(size = 5) + xlab(paste0("PC1: ", 
                                                                                                      round(percentVar[1] * 100), "% variance")) + 
          ylab(paste0("PC2: ", round(percentVar[2] * 100), 
                      "% variance")) + coord_fixed() + scale_fill_carto_d( palette = "Safe")
        
  }
  
  g <- g+theme_cowplot()
  
  return(g)
}

makeMap <- function(res) {
  
  ENTRY_DEL = "\001"
  KEY_DEL = "\002"
  
  response = content(res, "text")
  
  arr = unlist(strsplit(response, ENTRY_DEL, fixed = TRUE))
  
  list_map <- list("")
  vec_map_names <- c("");
  
  for (str in arr) {
    arrKeyValue = unlist(strsplit(str, KEY_DEL, fixed = TRUE));
    
    if (length(arrKeyValue) > 1) {
      list_map[length(list_map) + 1] <- arrKeyValue[2]
      vec_map_names[length(vec_map_names) + 1] <- arrKeyValue[1]
    }
  }
  
  names(list_map) <- vec_map_names
  
  list_map
}

getMirDipTargets <- function(microRNAs,minimumScore="High") {
  
  # const values
  url <- "http://ophid.utoronto.ca/mirDIP"
  
  mapScore <- list("0", "1", "2", "3");
  names(mapScore) <- c("Very High", "High", "Medium", "Low")
  
  
  
    parameters <- list(
      genesymbol = "",
      microrna = microRNAs,
      scoreClass = mapScore[minimumScore]
    )
    
    # ... send http POST
    res <- POST(paste(url, "/Http_U", sep = ""), body = parameters, encode = "form", verbose())#
    res <- makeMap(res)$results
    tmp <- tempfile()
    writeLines(res,tmp)
    res <- read.delim(tmp)
    return(res)
  }
  
plotSkeletalVis <- function(similarityScores){
  
  
  #get the rank of the similarity
  similarityScores$rank <- rank(-similarityScores$zscore)
  similarityScores$ID <- make.names(similarityScores$accession,unique = TRUE)
  
  chosen <- similarityScores[ similarityScores$zscore>2,]
  chosen <- chosen[    grepl("surgery|injured|Sham|ACLT",chosen$comparisonsText),"ID"]
  
  similarityScores$test <- ifelse(similarityScores$ID %in% chosen,"PTOA","no")
  similarityScores <- similarityScores %>% mutate(`PTOA Model` = case_when(
    accession == "GSE112641" ~ "ACL rupture",
    accession %in% c("GSE110268_3","GSE121033","GSE125320") ~ "ACLT",
    accession %in% c("GSE147529","GSE41342","GSE53857","GSE147529") ~ "DMM",
    accession %in% c("GSE103416") ~ "MIA",
    TRUE ~ "NA"))
  
  similarityScores[!similarityScores$ID %in% chosen, "PTOA Model"] <- "NA"
  
  
  ggplot(similarityScores %>% arrange(desc(`PTOA Model`)),aes(x=zscore,y=rank)) +
    geom_point(aes(colour=`PTOA Model`), size=3, alpha=1) +
    scale_y_continuous(trans='log10')+
    scale_colour_manual(values = c("#88CCEE", "#CC6677","#DDCC77","#DDCC77"),limits = c('ACL rupture', 'ACLT',"DMM")) +
    theme_cowplot(font_size = 28) +
    ylab("Rank of similarity") +
    xlab("Z-Score of cosine similarity") +
    geom_vline(xintercept=2, colour='blue',linetype=2) +
    theme(axis.line = element_line(size = 1)) +
    theme(legend.position = c(0.1, 0.2))
  
  
  
}

