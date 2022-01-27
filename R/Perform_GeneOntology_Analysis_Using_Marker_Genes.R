#' A Perform_GeneOntology_Analysis_Using_Marker_Genes Function
#' Except EnrichR, all others uses DEGs
#'
#' This function allows you to perform differential gene analysis.
#' @param Temp.object A list of Seurat objects between which to find anchors for downstream integration.
#' @param Temp.MarkerList A data.frame containg marker genes.
#' @param saveDIR Path to save generated data.
#' @param SuffixName Suffix. to be added in the directory name as 
#' @param methods ClusterProfiler (enrichGO)
#' @param GroupNameDEGs ClusterProfiler (enrichGO)
#' @param OrgDB ClusterProfiler (enrichGO)
#' @param OrgName ClusterProfiler (enrichGO)
#' @param EnrichDB EnrichR DB. For list of DBs: listEnrichrDbs()... https://maayanlab.cloud/Enrichr/#libraries
#' @param topnumber Max #pathways to be plotted per group
#' @param PDFheight PDFheight
#' @param PDFwidth PDFwidth
#' @keywords Temp.object, saveDIR, msigdbrCategoryList="H", SubmsigdbrCategorySubset = "", minSizeGeneSet=5, topnumber=10, PDFheight=10, PDFwidth=12
#' @export
#' @examples
#' Perform_GeneOntology_Analysis_Using_Marker_Genes()



Perform_GeneOntology_Analysis_Using_Marker_Genes <- function(Temp.object, Temp.MarkerList, saveDIR, SuffixName="GO",
                                 methods=c("ClusterProfilerM","ClusterProfilerAlt1", "ClusterProfilerAlt2"), EnrichDB="GO_Biological_Process_2021",
                                 GroupNameDEGs="cluster", OrgDB=org.Mm.eg.db, OrgName="mmu",
                                 topnumber=10, PDFheight=12, PDFwidth=18){
  
  #Temp.object=Temp.object
  #Temp.MarkerList=markers.Info
  #saveDIR=GOdirMain
  #methods=c("ClusterProfiler","ClusterProfilerAlt1", "ClusterProfilerAlt2")
  #GroupNameDEGs="cluster"
  #OrgDB=org.Mm.eg.db ### eval(as.name(OrgDB))
  #EnrichDB="GO_Biological_Process_2021"
  #OrgName="mmu"
  #topnumber=25
  #PDFheight=12
  #PDFwidth=18
  
  setwd(saveDIR)
  GOdirMain <- paste(getwd(),paste0("GeneOntology_Analysis"),sep="/"); print(GOdirMain)
  dir.create(file.path(getwd(),paste0("GeneOntology_Analysis")), showWarnings = FALSE)
  setwd(GOdirMain)
  
  
  print(paste0("Top Genes to plot:",topnumber))
  print(table(Temp.MarkerList[,GroupNameDEGs]))
  
  ## EnrighGO
  # https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/deg-go-enrichment.html#load-deg
  # https://github.com/cotneylab/scRNA_EnamelKnot/blob/main/2_Seurat_tooth_scRNA_for_publication.Rscript
  if("ClusterProfilerM" %in% methods){
    
    ## DEGs GO
    library(Seurat)
    library(tidyverse)
    library(magrittr)
    library(ReactomePA)
    library(clusterProfiler)
    library(enrichplot)
    library(org.Mm.eg.db)
    library(org.Hs.eg.db)
    library(DOSE)
    
    print("Processing GO using ClusterProfiler")
    
    setwd(GOdirMain)
    GOdirSub <- paste(getwd(),paste0("GO_ClusterProfiler"),sep="/"); print(GOdirSub)
    dir.create(file.path(getwd(),paste0("GO_ClusterProfiler")), showWarnings = FALSE)
    
    setwd(GOdirSub)
    pdf(paste0(SuffixName,"_",GroupNameDEGs,"_Based_using_ClusterProfiler_top_",topnumber,".pdf"),height=PDFheight, width=PDFwidth)
    
    
    ## Get gene names per cluster
    deg.ls <- split(rownames(Temp.MarkerList), f = Temp.MarkerList[,GroupNameDEGs]); head(deg.ls)
    
    #Transfer gene symbol into entrez id
    geneid.ls <- deg.ls %>% map(~{
      
      # here for macaque
      gene.df <- AnnotationDbi::select(OrgDB,
                                       keys = .x,
                                       columns = c("ENTREZID", "SYMBOL"),
                                       keytype = "SYMBOL")
      
      gene <- gene.df$ENTREZID
      gene <- gene[which(!is.na(gene))]
      gene <- unique(gene)
      
      return(gene)
    })
    
    
    ## GO for gene list
    ## go enrichment for all the variable
    gene.ls <- geneid.ls
    
    # enrichGO
    compGO <- compareCluster(geneCluster   = gene.ls,
                             fun           = "enrichGO",
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH", 
                             OrgDb = OrgDB, 
                             ont = 'BP')
    
    # enrichKEGG
    compKEGG <- compareCluster(geneCluster   = gene.ls,
                               fun           = "enrichKEGG",
                               pvalueCutoff  = 0.05,
                               pAdjustMethod = "BH", 
                               organism = OrgName)
    
    #compPathway <- compareCluster(geneCluster   = gene.ls,
    #                              fun           = "enrichPathway",
    #                              pvalueCutoff  = 0.05,
    #                              pAdjustMethod = "BH")
    
    ## dot plot
    g1 <- dotplot(compGO, showCategory = topnumber, title = "GO Enrichment Analysis using ClusterProfiler (enrichGO)")
    g2 <- dotplot(compKEGG, showCategory = topnumber, title = "KEGG Pathway Enrichment Analysis using ClusterProfiler (enrichKEGG)")
    #g3 <- dotplot(compPathway, showCategory = 10, title = "REACTOME Pathway Enrichment Analysis")
    print(g1)
    print(g2)
    
    ## go enrichment per cluster (individual)
    go.ls <- geneid.ls %>% map(~{
      
      eGO <- enrichGO(
        gene          = .x,
        OrgDb         = OrgDB,
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        readable      = TRUE
      )
      
      return(eGO)
      
    })
    
    
    kegg.ls <- gene.ls %>% map(~{
      eKEGG <- enrichKEGG(
        gene = .x,
        pvalueCutoff = 0.05, 
        organism = OrgName
      )
      return(eKEGG)
    })
    
    
    #pathway.ls <- gene.ls %>% map(~{
    #  ePathway <- enrichPathway(
    #    gene = .x,
    #    pvalueCutoff = 0.05,
    #    readable = TRUE, 
    #  )
    #  return(ePathway)
    #})
    
    
    ## enrichment visu
    barplotTerm <- function(object,
                            x = "Count",
                            color = 'p.adjust',
                            showCategory = topnumber,
                            font.size = 12,
                            title = "") {
      ## use *height* to satisy barplot generic definition
      ## actually here is an enrichResult object.
      colorBy <- color
      
      df <- fortify(object, showCategory = showCategory, by = x)
      df$p.adjust <- -log10(df$p.adjust)
      #df <- df[c(1:3,9:12,15,16),]
      if (colorBy %in% colnames(df)) {
        p <-
          ggplot(df, aes_string(x = x, y = "Description", fill = colorBy)) +
          theme_dose(font.size) +
          scale_fill_continuous(
            low = "red",
            high = "blue",
            name = color,
            guide = guide_colorbar(reverse = TRUE)
          )
      } else {
        p <- ggplot(df, aes_string(x = x, y = "Description")) +
          theme_dose(font.size) +
          theme(legend.position = "none")
      }
      
      
      p + geom_col(fill = color) + ggtitle(title) + xlab('-log10 FDR') + ylab(NULL)
    }
    
    
    
    lapply(1:length(gene.ls), function(x){
      
      name <- names(gene.ls)[[x]]
      g1 = barplotTerm(go.ls[[x]], showCategory = topnumber, title = paste0(name, " GO using ClusterProfiler (enrichGO)"), color = 'blue', x = 'p.adjust')
      g2 = barplotTerm(kegg.ls[[x]], showCategory = topnumber, title = paste0(name, " KEGG using ClusterProfiler (enrichKEGG)"), color = 'blue', x = 'p.adjust')
      print(plot_grid(g1, g2, ncol = 2))
      
    })
    
    dev.off()
    
    
  }
  
  
  ## EnrighGO
  # https://github.com/cotneylab/scRNA_EnamelKnot/blob/main/2_Seurat_tooth_scRNA_for_publication.Rscript
  # 
  if("ClusterProfilerAlt1" %in% methods){
  
    print("Processing GO using ClusterProfilerAlt1")
    
    setwd(GOdirMain)
    GOdirSub <- paste(getwd(),paste0("GO_ClusterProfiler_Alt1"),sep="/"); print(GOdirSub)
    dir.create(file.path(getwd(),paste0("GO_ClusterProfiler_Alt1")), showWarnings = FALSE)
    
    setwd(GOdirSub)
    pdf(paste0(SuffixName,"_",GroupNameDEGs,"_Based_using_ClusterProfiler_Alt1_top_",topnumber,".pdf"),height=PDFheight, width=PDFwidth)
    
    TotalGroups=sort(unique(Temp.MarkerList[,GroupNameDEGs]))
    for(GroupNameDEGsInd in TotalGroups){
      #GroupNameDEGsInd="NZO_LF"
    
    print(paste0("Processing ",GroupNameDEGsInd," (",which(TotalGroups %in% GroupNameDEGsInd)," of total ",length(TotalGroups),")"))
      
    Temp.MarkerList.Ind <- Temp.MarkerList[Temp.MarkerList[,GroupNameDEGs] == GroupNameDEGsInd,]
        
    f<-c()
    n<-tail(Temp.MarkerList.Ind[order(Temp.MarkerList.Ind$avg_log2FC),], 50)
    f<-rbind(f, n); print(dim(f))
    
    value_bp <- enrichGO(gene = f$gene,
                         universe = rownames(Temp.object@assays$RNA@data),
                         OrgDb = OrgDB, 
                         keyType = 'SYMBOL',
                         readable = F,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.10,
                         pAdjustMethod = "BH")
    
    
    
    value_bp1<-simplify(value_bp, cutoff=0.75, by="qvalue")
    value_bp1<-as.data.frame(value_bp)
    print(value_bp1)
    assign(paste(GroupNameDEGsInd, "ontology", sep="_"), value_bp1)
    write.table(value_bp1, paste0("GO_Enriched_",GroupNameDEGs,"_Based_using_ClusterProfiler_Alt1_top_",topnumber,"_for_",GroupNameDEGsInd,".txt"), quote = FALSE, sep = "\t", row.names = TRUE)
    p1 <- dotplot(value_bp, showCategory=topnumber, title = paste0("enrichGO: ",GroupNameDEGsInd))
    p2 <- barplot(value_bp, showCategory=topnumber, title = paste0("enrichGO: ",GroupNameDEGsInd))
    print(plot_grid(p1, p2, ncol = 2))
    
    }
    
    dev.off()
    
  }
  
  
  
  ## EnrighGO
  # https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/
  # 
  if("ClusterProfilerAlt2" %in% methods){
    
    print("Processing GO using ClusterProfilerAlt2")
    
    library(clusterProfiler)
    library(wordcloud)
    
    setwd(GOdirMain)
    GOdirSub <- paste(getwd(),paste0("GO_ClusterProfiler_Alt2"),sep="/"); print(GOdirSub)
    dir.create(file.path(getwd(),paste0("GO_ClusterProfiler_Alt2")), showWarnings = FALSE)
    
    setwd(GOdirSub)
    pdf(paste0(SuffixName,"_",GroupNameDEGs,"_Based_using_ClusterProfiler_Alt2_top_",topnumber,".pdf"),height=PDFheight, width=PDFwidth)
    
    TotalGroups=sort(unique(Temp.MarkerList[,GroupNameDEGs]))
    for(GroupNameDEGsInd in TotalGroups){
      #GroupNameDEGsInd="B6_LF"
      
      print(paste0("Processing ",GroupNameDEGsInd," (",which(TotalGroups %in% GroupNameDEGsInd)," of total ",length(TotalGroups),")"))
      
      Temp.MarkerList.Ind <- Temp.MarkerList[Temp.MarkerList[,GroupNameDEGs] == GroupNameDEGsInd,]
      
    #### THIS NEEDS TO BE CHANGED FOR UNIVERSE GENES
    ## we want the log2 fold change 
    #original_gene_list <- marker.LF.Beta$avg_log2FC; length(gene_list)
    ## name the vector
    #names(original_gene_list) <- marker.LF.Beta$gene
    ## omit any NA values 
    #gene_list<-na.omit(original_gene_list)
    ## sort the list in decreasing order (required for clusterProfiler)
    ##gene_list = sort(gene_list, decreasing = TRUE); length(gene_list)
    #### TILL HERE.. THIS NEEDS TO BE CHANGED FOR UNIVERSE GENES
    
    # Exctract significant results (padj < 0.05)
    sig_genes_df = subset(Temp.MarkerList.Ind, p_val_adj < 0.05); dim(sig_genes_df)
    # From significant results, we want to filter on log2fold change
    genes <- sig_genes_df$avg_log2FC
    
    # Name the vector
    names(genes) <- sig_genes_df$gene
    # omit NA values
    genes <- na.omit(genes)
    genes <- head(genes,50); print(genes)
    # filter on min log2fold change (log2FoldChange > 2)
    #genes <- names(genes)[abs(genes) > 2]
    genesNames <- names(genes)
    
    go_enrich1 <- enrichGO(gene = genesNames,
                          universe = rownames(Temp.object@assays$RNA@data),
                          OrgDb = OrgDB, 
                          keyType = 'SYMBOL',
                          readable = F,
                          ont = "BP",
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.10,
                          pAdjustMethod = "BH")
    
    go_enrich<-simplify(go_enrich1, cutoff=0.75, by="qvalue")
    
    library(enrichplot)
    p1 <- dotplot(go_enrich, showCategory=topnumber, title = paste0("enrichGO: ",GroupNameDEGsInd), font.size = 10)
    p2 <- barplot(go_enrich, showCategory = topnumber, title = paste0("enrichGO: ",GroupNameDEGsInd), drop = TRUE,  font.size = 10)
    print(plot_grid(p1, p2, ncol = 2))
    
    p3 <- upsetplot(go_enrich)
    print(p3)
    
    wcdf<-read.table(text=go_enrich$GeneRatio, sep = "/")[1]
    wcdf$term<-go_enrich[,2]
    wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)
    
    
    #emapplot(go_enrich)
    
    print(goplot(go_enrich, showCategory = topnumber, title = paste0("enrichGO: ",GroupNameDEGsInd)))
    
    # categorySize can be either 'pvalue' or 'geneNum'
    print(cnetplot(go_enrich, categorySize="pvalue", foldChange=genes, title = paste0("enrichGO: ",GroupNameDEGsInd)))
    
    #library(pathview)
    ## Produce the native KEGG plot (PNG)
    #dme <- pathview(gene.data=gene_list, pathway.id="dme04080", species = kegg_organism, gene.idtype=gene.idtype.list[3])
    ## Produce a different plot (PDF) (not displayed here)
    #dme <- pathview(gene.data=gene_list, pathway.id="dme04080", species = kegg_organism, gene.idtype=gene.idtype.list[3], kegg.native = F)
    #knitr::include_graphics("dme04080.pathview.png")
    
    }
    
  dev.off()  
  }
  
  
  
  
  ## https://ucdavis-bioinformatics-training.github.io/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/scrnaseq_analysis/scRNA_Workshop-PART6.html
  ## https://ucdavis-bioinformatics-training.github.io/2021-March-Single-Cell-RNA-Seq-Analysis/data_analysis/scRNA_Workshop-PART6_fixed
  ## TopGO
  if("TopGO" %in% methods){
    
    print("Processing GO using TopGO")
    
    library(topGO)
    
    setwd(GOdirMain)
    GOdirSub <- paste(getwd(),paste0("GO_TopGO"),sep="/"); print(GOdirSub)
    dir.create(file.path(getwd(),paste0("GO_TopGO")), showWarnings = FALSE)
    
    setwd(GOdirSub)
    pdf(paste0(SuffixName,"_",GroupNameDEGs,"_Based_using_TopGO_top_",topnumber,".pdf"),height=PDFheight, width=PDFwidth)
    
    TotalGroups=sort(unique(Temp.MarkerList[,GroupNameDEGs]))
    for(GroupNameDEGsInd in TotalGroups){
      #GroupNameDEGsInd="B6_LF"
      
      print(paste0("Processing ",GroupNameDEGsInd," (",which(TotalGroups %in% GroupNameDEGsInd)," of total ",length(TotalGroups),")"))
      
      Temp.MarkerList.Ind <- Temp.MarkerList[Temp.MarkerList[,GroupNameDEGs] == GroupNameDEGsInd,]
      
      # install org.Mm.eg.db if not already installed (for mouse only)
      # install.packages("org.Mm.eg.db") 
      print(head(Idents(Temp.object)))
      
      #SCdata.Temp <- SubsetData(Temp.object, ident.use = GroupNameUse)
      # install org.Mm.eg.db from Bioconductor if not already installed (for mouse only)
      SCdata.Temp <- subset(Temp.object, idents = GroupNameDEGsInd)
      expr <- as.matrix(GetAssayData(SCdata.Temp))
      # Select genes that are expressed > 0 in at least 75% of cells (somewhat arbitrary definition)
      n.gt.0 <- apply(expr, 1, function(x)length(which(x > 0)))
      expressed.genes <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.75)]
      all.genes <- rownames(expr)
      
      # define geneList as 1 if gene is in expressed.genes, 0 otherwise
      geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
      names(geneList) <- all.genes
      
      if(OrgName == "mmu"){topGOdbName = "org.Mm.eg.db"
      } else if(OrgName == "hsa"){topGOdbName = "org.Hs.eg.db"}
      
      print(paste0("TopGO DB Name: ",topGOdbName))
      
      # Create topGOdata object
      GOdata <- new("topGOdata",
                    ontology = "BP", # use biological process ontology
                    allGenes = geneList,
                    geneSelectionFun = function(x)(x == 1),
                    annot = annFUN.org, mapping = topGOdbName, ID = "symbol")
      
      # Test for enrichment using Fisher's Exact Test
      ## First, we perform a classical enrichment analysis by testing the over-representation of GO terms within the
      ## group of differentially expressed genes. For the method classic each GO category is tested independently.
      resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
      resultFisher
      
      #Next we will test the enrichment using the Kolmogorov-Smirnov test. We will use the both the classic and the elim method.
      resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
      resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
      
      #In the following example, we list the top 10 significant GO terms identified by the elim method. At
      #the same time we also compare the ranks and the p-values of these GO terms with the ones obtained by the classic method.
      allRes <- GenTable(GOdata, classicFisher = resultFisher,
                         classicKS = resultKS, elimKS = resultKS.elim,
                         orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 50, numChar = 100)
      head(allRes)
      
      write.table(allRes, paste0("GO_Enriched_",GroupNameDEGs,"_Based_using_TopGO_top_",topnumber,"_for_",GroupNameDEGsInd,".txt"), quote = FALSE, sep = "\t", row.names = TRUE)
      
      require(ggplot2)
      library(scales)
      
      ntop <- 30
      goEnrichment <- allRes
      goEnrichment$elimKS <- as.numeric(goEnrichment$elimKS)
      goEnrichment <- goEnrichment[goEnrichment$elimKS < 0.05,] # filter terms for KS p<0.05
      goEnrichment <- goEnrichment[,c("GO.ID","Term","elimKS")]
      goEnrichment
      ggdata <- goEnrichment[1:ntop,]
      ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
      gg1 <- ggplot(ggdata,
                    aes(x = Term, y = -log10(elimKS), size = -log10(elimKS), fill = -log10(elimKS))) +
        
        expand_limits(y = 1) +
        geom_point(shape = 21) +
        scale_size(range = c(2.5,12.5)) +
        scale_fill_continuous(low = 'royalblue', high = 'red4') +
        
        xlab('') + ylab('Enrichment score') +
        labs(
          title = paste0('GO Biological processes: ',GroupNameDEGsInd),
          subtitle = 'Top 30 terms ordered by Kolmogorov-Smirnov p-value',
          caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
        
        geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
                   linetype = c("dotted", "longdash", "solid"),
                   colour = c("black", "black", "black"),
                   size = c(0.5, 1.5, 3)) +
        
        theme_bw(base_size = 24) +
        theme(
          legend.position = 'right',
          legend.background = element_rect(),
          plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
          plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
          plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
          
          axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
          axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
          axis.title = element_text(size = 12, face = 'bold'),
          axis.title.x = element_text(size = 12, face = 'bold'),
          axis.title.y = element_text(size = 12, face = 'bold'),
          axis.line = element_line(colour = 'black'),
          
          #Legend
          legend.key = element_blank(), # removes the border
          legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
          legend.text = element_text(size = 14, face = "bold"), # Text size
          title = element_text(size = 14, face = "bold")) +
        
        coord_flip()
      
      print(gg1)
      
      par(cex = 0.5)
      print(showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = 2, useInfo = 'all'))
    }
    
    dev.off()  
  }
  
  
  
  ## https://ucdavis-bioinformatics-training.github.io/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/scrnaseq_analysis/scRNA_Workshop-PART6.html
  ## https://ucdavis-bioinformatics-training.github.io/2021-March-Single-Cell-RNA-Seq-Analysis/data_analysis/scRNA_Workshop-PART6_fixed
  ## EnrichR
  if("EnrichR" %in% methods){
    
    print("Processing GO using EnrichR")
    
    ### EnrichR Database list
    ## https://maayanlab.cloud/Enrichr/#libraries
    #DEenrichRPlot(object = Temp.object, ident.1 = "Basal_1",ident.2 = "Basal_2",enrich.database = "GO_Biological_Process_2021",max.genes = 50, num.pathway = 20)
    library(enrichR)
    dbs <- listEnrichrDbs()
    
    setwd(GOdirMain)
    GOdirSub <- paste(getwd(),paste0("GO_EnrichR"),sep="/"); print(GOdirSub)
    dir.create(file.path(getwd(),paste0("GO_EnrichR")), showWarnings = FALSE)
    
    setwd(GOdirSub)
    pdf(paste0(SuffixName,"_",GroupNameDEGs,"_Based_using_EnrichR_top_",topnumber,".pdf"),height=PDFheight, width=PDFwidth)
    
    TotalGroups=sort(unique(Temp.MarkerList[,GroupNameDEGs]))
    for(GroupNameDEGsInd in TotalGroups){
      #GroupNameDEGsInd="B6_LF"
      
      # install org.Mm.eg.db if not already installed (for mouse only)
      # install.packages("org.Mm.eg.db") 
      print(head(Idents(Temp.object)))
      print(paste0("Processing EnrichR for: ",GroupNameDEGsInd))
      
      EnrichDGEs <- DEenrichRPlot(object = Temp.object, ident.1 = GroupNameDEGsInd,enrich.database = EnrichDB,max.genes = 50, num.pathway = 50)
      head(EnrichDGEs$data)
      
      write.table(EnrichDGEs$data, paste0("GO_Enriched_",GroupNameDEGs,"_Based_using_EnrichR_top_",topnumber,"_for_",GroupNameDEGsInd,".txt"), quote = FALSE, sep = "\t", row.names = TRUE)
      print(EnrichDGEs)
      
    }
    
    dev.off()  
  }
  
  
  
  
#Function
}





