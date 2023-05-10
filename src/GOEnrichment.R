
GOEnrichment <- function (differentialExpression, 
                          geneLengths,species,foldChangeOnly, 
                          foldchange, 
                          padj,direction="both",decreaseSize=TRUE) 
{
  suppressPackageStartupMessages(library("goseq"))
  
  
  
  
  suppressPackageStartupMessages(library("GO.db"))
  suppressPackageStartupMessages(library("GOSemSim"))
  suppressPackageStartupMessages(library("tidyr"))
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("org.Hs.eg.db"))

  getEnrichedGOTerms <- function(diffExp, sigDiffExp, species,geneLengths, 
                                 threshold = 0.05,decreaseSize) {
    
    if (species == "Human") {
      suppressPackageStartupMessages(library("org.Hs.eg.db"))
      gene2GO <- AnnotationDbi::select(org.Hs.eg.db, keys(org.Hs.egGO2EG), 
                                       c("ENTREZID", "SYMBOL"), "GOALL")
    }        else if (species == "Mouse") {
      suppressPackageStartupMessages(library("org.Mm.eg.db"))
      gene2GO <- AnnotationDbi::select(org.Mm.eg.db, keys(org.Mm.egGO2EG), 
                                       c("ENTREZID", "SYMBOL"), "GOALL")
    }
    gene2GO <- unstack(gene2GO[, c(1, 5)])
    genes <- ifelse(diffExp[, 1] %in% sigDiffExp[, 1], 1, 
                    0)
    names(genes) <- diffExp[, 1]
    geneLengths <- geneLengths[match(diffExp[, 1], geneLengths$gene_name), ]
    x <- nullp(genes, bias.data = geneLengths$length, plot.fit = F)
    GOTerms <- goseq(x, gene2cat = gene2GO)
    GOTerms <- GOTerms[GOTerms$ontology == "BP", ]
    GOTerms$padj <- p.adjust(GOTerms$over_represented_pvalue, 
                             method = "BH")
    GOTerms.sig <- GOTerms[GOTerms$padj <= threshold, ]
    if (nrow(GOTerms.sig) == 0) {
      print("no significant GO terms!")
      return(NULL)
    }
    GOTerms.sig$enrichment <- GOTerms.sig$numDEInCat/GOTerms.sig$numInCat
    
    GOTerms.sig$Fold_Enrichment <- GOTerms.sig$numDEInCat/nrow(sigDiffExp)/(GOTerms.sig$numInCat/nrow(diffExp))
    
    
    GOResults = list()
    for (i in 1:nrow(GOTerms.sig)) {
      GOTerm <- GOTerms.sig$category[i]
      index <- sapply(gene2GO, function(x) GOTerm %in% 
                        x)
      termIDs <- names(index[index == "TRUE"])
      sig <- sigDiffExp[sigDiffExp[, 1] %in% termIDs, ]
      GOResults[[GOTerm]] = sig[, 1]
    }
    names(GOResults) = GOTerms.sig$term
    GOResults <- lapply(GOResults, function(x) paste(x, sep = "", 
                                                     collapse = " "))
    GOResults <- data.frame(Term = names(GOResults), ID = GOTerms.sig$category, numGenes=GOTerms.sig$numDEInCat,
                            Genes = unlist(GOResults), Adj.pvalue = GOTerms.sig$padj, 
                            PercentageCoverage = GOTerms.sig$enrichment*100, FoldEnrichment=GOTerms.sig$Fold_Enrichment)
    print("simplify")
    if(nrow(GOResults)<3) return(list(GOResults = GOResults,GOResults.reduced=NULL))
    GOResults.reduced <- try(simplify(GOResults, gene2GO,decreaseSize))
    if(inherits(GOResults.reduced,"try-error")) return(list(GOResults = GOResults,GOResults.reduced=NULL))
    return(list(GOResults = GOResults, GOResults.reduced = GOResults.reduced))
  }
  simplify <- function(GORes, gene2GO,decreaseSize=TRUE) {
    if (species == "Human") {
      semData <- godata("org.Hs.eg.db", "SYMBOL", "BP")
    }        else if (species == "Mouse") {
      semData <- godata("org.Mm.eg.db", "SYMBOL", "BP")
    }
    sim <- mgoSim(GORes$ID, GORes$ID, semData = semData, 
                  measure = "Rel", combine = NULL)
    sim[is.na(sim)] <- 0
    
    if(decreaseSize){
      print("decreasing")
    
    go1 <- go2 <- similarity <- NULL
    sim.df <- as.data.frame(sim)
    sim.df$go1 <- row.names(sim.df)
    sim.df <- gather(sim.df, go2, similarity, -go1)
    sim.df <- sim.df[!is.na(sim.df$similarity), ]
    sim.df <- sim.df[order(sim.df$similarity, decreasing = T), 
                     ]
    sim.df <- sim.df[sim.df$go1 != sim.df$go2, ]
    sim.df <- sim.df[sim.df$similarity > 0.4, ]
    GO2Gene <- unstack(stack(gene2GO)[2:1])
    freq <- sapply(GO2Gene, length)
    freqCutOff <- length(gene2GO) * 0.05
    highFreqTerms <- names(freq[freq > freqCutOff])
    sim.df$remove <- apply(sim.df, 1, function(x) {
      if (x[1] %in% highFreqTerms) {
        return(x[1])
      }
      if (x[2] %in% highFreqTerms) {
        return(x[2])
      }            else {
        return(NA)
      }
    })
    remove <- na.omit(sim.df$remove)
    sim.df <- sim.df[is.na(sim.df$remove), ]
    sim.df$go1.pval <- GORes$Adj.pvalue[match(sim.df$go1, 
                                              GORes$ID)]
    sim.df$go2.pval <- GORes$Adj.pvalue[match(sim.df$go2, 
                                              GORes$ID)]
    childTerms <- as.list(GOBPCHILDREN)
    print(sim.df)
    for (i in 1:nrow(sim.df)) {
      if (sim.df[i, "go1"] %in% remove) {
        next
      }
      if (sim.df[i, "go2"] %in% remove) {
        next
      }
      go1.pval <- sim.df[i, "go1.pval"]
      go2.pval <- sim.df[i, "go2.pval"]
      if (go1.pval == go2.pval) {
        go1 <- sim.df[i, "go1"]
        go2 <- sim.df[i, "go2"]
        if (go2 %in% childTerms[[go1]]) {
          remove <- c(remove, go2)
          next
        }                else if (go1 %in% childTerms[[go2]]) 
          remove <- c(remove, go1)
        next
      }
      remove <- c(remove, sim.df[i, which.max(c(go1.pval, 
                                                go2.pval))])
    }
    GORes.filt <- GORes[!GORes$ID %in% remove, ]
    sim.filt <- sim[as.character(GORes.filt$ID), as.character(GORes.filt$ID)]
    } else {
      GORes.filt <- GORes
      sim.filt <- sim[as.character(GORes$ID), as.character(GORes$ID)]
    }
    
    fit <- cmdscale(1 - sim.filt, eig = TRUE, k = 2)
    x <- fit$points[, 1]
    y <- fit$points[, 2]
    GORes.filt.plot <- GORes.filt
    GORes.filt.plot$x <- x
    GORes.filt.plot$y <- y
    GORes.filt.plot$log10Adjpvalue <- -log10(GORes.filt.plot$Adj.pvalue)
    return(GORes.filt.plot)
  }
  
  differentialExpression <- na.omit(differentialExpression)
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
  GOResultsList <- getEnrichedGOTerms(differentialExpression, 
                                      differentialExpression.sig, species, geneLengths,0.05,decreaseSize)
  
  
  
  return(GOResultsList)
}


plotGOMDS <- function(plotData,k=4){
  
  plotData$group <- kmeans(plotData[,c("x","y")],k)$cluster
  plotData$Term <- str_wrap(plotData$Term,width = 10)
  terms <- plotData %>% group_by(group) %>% dplyr::slice_max(log10Adjpvalue,n=1)
  
  
  
  #plot the top GO terms grouped by semantic similarity
  goPlot <- plotData %>% ggplot(aes(x=x, y=y, size = PercentageCoverage,fill=log10Adjpvalue)) +
    scale_colour_viridis_c(option="A") +
    theme_cowplot(font_size = 36) +
    geom_point(alpha=0.5, shape=21, color="black") +
    scale_size(range = c(.5, 24), name="Coverage",guide=F) +
    ylab("Dim1") +
    xlab("Dim2") + ggrepel::geom_text_repel(data = subset(plotData,Term %in% terms$Term),inherit.aes = FALSE,aes(x=x,y=y,label=Term),size=7,min.segment.length = 0,max.overlaps = 5) +
   theme(axis.line = element_line(size = 1.3))

  
  
  return(goPlot)
}


plotEnrichment <- function(results,n=10,IDs=NULL) {
  if(is.null(IDs)){
    results <- results[ order(results$Adj.pvalue),][1:n,]
  } else {
    results <- results[ results$ID %in% IDs,]
  }
  
  
  
  results$log_p <- -log10(results$Adj.pvalue)
  results$Term <- str_wrap(results$Term,width=25)
  results$Term <- fct_reorder(results$Term,results$Adj.pvalue,.desc = TRUE)
  g <- ggplot(results,aes(x = FoldEnrichment,y=Term)) +
    geom_point(aes(color = log_p,size=numGenes)) + theme_cowplot(font_size = 38) +
    xlab("Fold Enrichment") + theme(axis.title.y = element_blank()) +
    theme(axis.text.y =element_text(size = 30)) +
    labs(size = "# genes", color = expression(-log[10](p))) +
    scale_size(range=c(8,25)) +
    theme(axis.line = element_line(size = 1.4)) +
    scale_color_gradient(low = "#EE8879FF", high = "#E64B35FF") +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)))
  
  
  return(g)
}
