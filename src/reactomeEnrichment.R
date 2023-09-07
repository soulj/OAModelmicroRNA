
getReactomePathways <- function(species){
  
  
  #the reactome pathways are only provided with entrez IDs
  #so we need to convert them to gene symbols to match the data we have
  
  symbol2eg <- dplyr::case_when(
    species == "Human" ~ "org.Hs.egSYMBOL",
    species == "Mouse" ~ "org.Mm.egSYMBOL",
    species == "Pig" ~ "org.Ss.egSYMBOL",
    species == "Cow" ~ "org.Bt.egSYMBOL",
    species == "Rat" ~ "org.Rn.egSYMBOL",
  )
  
  symbol2eg <- as.list(get(symbol2eg))
  
  #get eg to reactome pathway
  reactome2eg <- as.list(reactomePATHID2EXTID)
  
  speciesID <- dplyr::case_when(
    species == "Human" ~ "R-HSA",
    species == "Mouse" ~ "R-MMU",
    species == "Pig" ~ "R-SSC",
    species == "Cow" ~ "R-BTA",
    species == "Rat" ~"R-RNO",
  )
  
  #filter to just species pathways
  reactome2eg <- reactome2eg[grep(speciesID,names(reactome2eg))]
  
  #function to search 
  grepREACTOME <- function(id,mapkeys){
    unique(unlist(mapkeys[id],use.names=FALSE))
  }
  
  #convert the entrez ids to gene symbols for each pathway
  reactome <- lapply(reactome2eg,grepREACTOME,symbol2eg)
  
  
  #get the pathway names rather than the ids
  reactome2name <- as.list(reactomePATHID2NAME)
  reactomeNames <- sapply(names(reactome),grepREACTOME,reactome2name)
  names(reactome) <- reactomeNames
  
  return(reactome)
}

getEnrichedPathways <- function(data, sigData, pathway,geneExonLengths, 
                                species, threshold = 0.05) {
  data <- as.data.frame(data)
  sigData <- as.data.frame(sigData)
  genes <- ifelse(as.data.frame(data)[, 1] %in% as.data.frame(sigData)[, 1], 1, 0)
  names(genes) <- data[, 1]
  
  geneExonLengths <- geneExonLengths[match(data[, 1], geneExonLengths$gene_name),]
  x <- nullp(genes, bias.data = geneExonLengths$length, plot.fit = T)
  PATHWAY = goseq(x, gene2cat = pathway)
  PATHWAY$padj = p.adjust(PATHWAY$over_represented_pvalue, method = "BH")
  PATHWAY$Coverage = PATHWAY$numDEInCat/PATHWAY$numInCat * 100
  PATHWAY$Adjpvaluelog = -log10(PATHWAY$padj)
  PATHWAY$Term <- unlist(PATHWAY$Term)
  PATHWAY.sig = PATHWAY[PATHWAY$padj <= threshold, ]
  if (nrow(PATHWAY.sig) == 0) {
    print("no significant pathways")
    quit(status=3)
  }
  pathwayResults = list()
  for (i in 1:nrow(PATHWAY.sig)) {
    pathwayTerm = PATHWAY.sig$category[i]
    termIDs = pathway[[pathwayTerm]]
    sig = sigData[sigData[, 1] %in% termIDs, ]
    pathwayResults[[pathwayTerm]] = sig[, 1]
  }
  names(pathwayResults) = PATHWAY.sig$category
  pathwayResults = lapply(pathwayResults, function(x) paste(x, 
                                                            sep = "", collapse = " "))
  
  print(head(pathwayResults))
  
  PATHWAY.sig$enrichment <- PATHWAY.sig$numDEInCat/PATHWAY.sig$numInCat
  PATHWAY.sig$Fold_Enrichment <- PATHWAY.sig$numDEInCat/nrow(sigData)/(PATHWAY.sig$numInCat/nrow(data))
  
  print(head(PATHWAY.sig))
  print(colnames(PATHWAY.sig))
  
  pathwayResults <- data.frame(Term = names(pathwayResults), numGenes=PATHWAY.sig$numDEInCat,
             Genes = unlist(pathwayResults), Adj.pvalue = PATHWAY.sig$padj, 
             PercentageCoverage = PATHWAY.sig$enrichment*100, FoldEnrichment=PATHWAY.sig$Fold_Enrichment)
  
  
  pathwayResults$Term <- gsub("Homo sapiens: ","",pathwayResults$Term )
  return(pathwayResults)
}

pathwayEnrichment <- function (differentialExpression, reactomePathways, geneExonLengths, foldChangeOnly=FALSE, 
                               foldchange = 1.5, padj = 0.05,direction="both") {    
  
  colnames(differentialExpression)[1:2] <- c("GeneSymbol","log2FC")
  if (foldChangeOnly == TRUE) {
    differentialExpression <- differentialExpression %>% 
      group_by(GeneSymbol) %>% dplyr::slice(which.max(abs(log2FC))) %>% 
      as.data.frame
    
    if(direction=="both"){
      differentialExpression.sig <- na.omit(differentialExpression[abs(differentialExpression[, 
                                                                                              2]) >= log2(foldchange), ])
    } else if(direction=="up"){
      differentialExpression.sig <- na.omit(differentialExpression[differentialExpression[, 
                                                                                          2] >= log2(foldchange), ])
    } else {
      differentialExpression.sig <- na.omit(differentialExpression[differentialExpression[, 
                                                                                          2] <= -log2(1/foldchange), ])
    }
    
    
  }    else {
    colnames(differentialExpression)[3] <- "padj"
    differentialExpression <- differentialExpression %>% 
      group_by(GeneSymbol) %>% dplyr::slice(which.min(padj)) %>% 
      as.data.frame
    
    
    if(direction=="both"){
      differentialExpression.sig <- na.omit(differentialExpression[abs(differentialExpression[, 
                                                                                              2]) >= log2(foldchange) & differentialExpression[, 
                                                                                                                                               3] <= padj, ])
    } else if(direction=="up"){
      differentialExpression.sig <- na.omit(differentialExpression[differentialExpression[,2] >= log2(foldchange) & differentialExpression[, 3] <= padj, ])
    } else {
      differentialExpression.sig <- na.omit(differentialExpression[differentialExpression[, 
                                                                                          2] <= -log2(1/foldchange) & differentialExpression[, 
                                                                                                                                             3] <= padj, ])
    }
  }
  
   pathwayResults <- getEnrichedPathways(differentialExpression, 
                                        differentialExpression.sig, reactomePathways,geneExonLengths)

}



plotEnrichment <- function(results,n=10,IDs=NULL) {
  if(is.null(IDs)){
    results <- results[ order(results$Adj.pvalue),][1:n,]
  } else {
    results <- results[ results$Term %in% IDs,]
  }
  
  results$log_p <- -log10(results$Adj.pvalue)
  results$Term <- str_wrap(results$Term,width=25)
  results$Term <- fct_reorder(results$Term,results$Adj.pvalue,.desc = TRUE)
  g <- ggplot(results,aes(x = FoldEnrichment,y=Term)) +
    geom_point(aes(color = log_p,size=numGenes)) + theme_cowplot(font_size = 38) +
    xlab("Fold Enrichment") + theme(axis.title.y = element_blank()) +
    theme(axis.text.y =element_text(size = 30)) +
    labs(size = "# genes", color = expression(-log[10](p))) +
    scale_size(range=c(8,25),breaks=c(2,5,10),guide = guide_legend(order = 1)) +
    theme(axis.line = element_line(size = 1.4)) +
    scale_color_gradient(low = "#EE8879FF", high = "#E64B35FF") +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)))
  
  
  return(g)
}
