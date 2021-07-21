#' A Doublet_Detection_SCDShybrid Function
#'
#' This function allows you to detect doublets using both Doublet_Finder algorithm.
#' @param SeuratObject Seurat Object.
#' @param saveDIR Path to save Quality plots and RDS data.
#' @param Sample Sample Name.
#' @param Species Species Name. Valid options are hsa or mmu.
#' @param FeatureUseCount Number of features to select as top variable features; only used when selection.method is set to 'dispersion' or 'vst'
#' @param PCAnum Number of PCs to be used
#' @param resClus Resolution to be used for clustering
#' @param ClusPallette Color pallete to be used for clusters
#' @param SCDShybrid Run SCDS hybrid Doublet Detection Algorithm
#' @param plotCCgene Plots TOP2A gene or not
#' @keywords SeuratObject, saveDIR, Sample, Species, FeatureUseCount, PCAnum, resClus, ClusPallette, SCDS plotCCgene
#' @export
#' @examples
#' Doublet_Detection_SCDShybrid()



Doublet_Detection_SCDShybrid <- function(SeuratObject, saveDIR, Sample, Species="hsa", FeatureUseCount=2500, PCAnum=10, resClus = 0.5, ClusPallette=ClusPallette, SCDShybrid = TRUE, plotCCgene = TRUE){
  
  #SeuratObject=SCdata
  #saveDIR=saveDIR
  #Sample=Sample
  #FeatureUseCount=2000
  #Species="mmu"
  #FeatureUseCount=2500
  #PCAnum=10
  #resClus = 0.5
  #ClusPallette=ClusPallette
  #SCDShybrid = TRUE
  #plotCCgene = TRUE
  
  #library(DoubletDecon)
  
  
  library(scds)
  library(scater)
  library(rsvd)
  library(Rtsne)
  library(cowplot)
  
  
  setwd(saveDIR)
  DDdir <- paste(getwd(),paste0("Doublet_Detection_",Sample),sep="/"); print(DDdir)
  dir.create(file.path(getwd(),paste0("Doublet_Detection_",Sample)), showWarnings = FALSE)
  
  print(paste0("PrePorocessing data before running Doublet Detection Algorithms for: ",Sample))
  RUNProcessData="YES"
  if(RUNProcessData == "YES"){
    
    setwd(DDdir)
    pdf(file=paste0("PreProcess_Doublet_Detection_",Sample,"_using_PCA_",PCAnum,"_res_",resClus,".pdf"),height = 10,width = 12)
    print(paste0("Loading Seurat object for: ",Sample))
    
    ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
    #SeuratObject <- NormalizeData(SeuratObject)
    # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
    
    if(Species=="hsa"){
      print("Counting MT % for Human")
      SeuratObject[["percent.mt"]] <- PercentageFeatureSet(SeuratObject, pattern = "^MT-")
      SeuratObject[["percent.rb"]] <- PercentageFeatureSet(SeuratObject, pattern = "^RP[SL]")
      
    } else {
      print("Counting MT % for Mouse")
      SeuratObject[["percent.mt"]] <- PercentageFeatureSet(SeuratObject, pattern = "^mt-")
      SeuratObject[["percent.rb"]] <- PercentageFeatureSet(SeuratObject, pattern = "^Rp[sl]")
    }
    
    print("Finished Counting")
    print(head(SeuratObject@meta.data))
    
    # Visualize QC metrics as a violin plot
    #p1 <- VlnPlot(SeuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), cols = ClusPallette, pt.size = 0.00, ncol = 1)
    #print(p1)
    #plot1 <- FeatureScatter(SeuratObject, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = ClusPallette)
    #plot2 <- FeatureScatter(SeuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = ClusPallette)
    #print(CombinePlots(plots = list(plot1, plot2)))
    
    print("Skipped VlnPlot")
    
    SeuratObject <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = FeatureUseCount)
    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(SeuratObject), 10)
    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(SeuratObject)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    print(CombinePlots(plots = list(plot1, plot2)))
    
    SeuratObject <- ScaleData(SeuratObject)
    SeuratObject <- RunPCA(SeuratObject)
    print(ElbowPlot(SeuratObject))
    print(ElbowPlot(SeuratObject, ndims = 50))
    print(DimHeatmap(SeuratObject, dims = 1:12, cells = 500, balanced = TRUE))
    
    ## Clusters
    SeuratObject <- FindNeighbors(SeuratObject, dims = 1:PCAnum)
    SeuratObject <- FindClusters(SeuratObject, resolution = resClus)
    
    Idents(object = SeuratObject) <- "seurat_clusters"
    #print(VlnPlot(SeuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), cols = ClusPallette, pt.size = 0.01, ncol = 2))
    
    print("Skipped VlnPlot")
    
    SeuratObject <- RunUMAP(SeuratObject, dims = 1:PCAnum)
    #print(DimPlot(SeuratObject, reduction = "umap", label=TRUE, label.size = 8, pt.size = 0.5, cols = ClusPallette) + ggtitle(paste0(Sample, " (",nrow(SeuratObject@meta.data)," cells)")))
    dev.off()
    
    print("Skipped DimPlot")
  }
  
  
  if (SCDS == TRUE) {
    
    setwd(DDdir)
    SCDSDir <- paste(getwd(),paste0("SCDS_",Sample),sep="/"); print(SCDSDir)
    dir.create(file.path(getwd(),paste0("SCDS_",Sample)), showWarnings = FALSE)
    
    setwd(SCDSDir)
    pdf(file=paste0("Plots_SCDS_",Sample,"_using_PCA_",PCAnum,"_res_",resClus,".pdf"),height = 10,width = 12)
    
    print(paste0("Started SCDS process for ",Sample, " using PCA ",PCAnum))
    
    sce <- as.SingleCellExperiment(SeuratObject, assay = "RNA")
    sce = cxds(sce,retRes = TRUE)
    sce = bcds(sce,retRes = TRUE,verb=TRUE)
    sce = cxds_bcds_hybrid(sce)
    
    CD  = as.data.frame(colData(sce))
    head(CD)
    head(cbind(CD$cxds_score,CD$bcds_score, CD$hybrid_score))
    
    par(mfcol=c(1,3))
    print(boxplot(sce$cxds_score   ~ sce$Project, main="cxds"))
    print(boxplot(sce$bcds_score   ~ sce$Project, main="bcds"))
    print(boxplot(sce$hybrid_score ~ sce$Project, main="hybrid"))
    #boxplot(sce$hybrid_score ~ sce$hybrid_call, main="hybrid")
    
    doublet_cutoff <- median(sce$hybrid_score) + 2*(sd(sce$hybrid_score))
    CD$DoubletSCDShybrid <- "Singlet"
    CD[CD$hybrid_score > doublet_cutoff, "DoubletSCDShybrid"] = "Doublet"
    head(CD)
    print("Estimated Singlets and Doublets")
    table(CD$DoubletSCDShybrid)
    
    print("Finished Doublet Detection steps")
    
    ## Run SCDS with varying classification stringencies ----------------------------------------------------------------
    head(SeuratObject@meta.data)
    SeuratObject@meta.data$cxds_score <- CD$cxds_score[match(rownames(SeuratObject@meta.data), rownames(CD))]
    SeuratObject@meta.data$bcds_score <- CD$bcds_score[match(rownames(SeuratObject@meta.data), rownames(CD))]
    SeuratObject@meta.data$hybrid_score <- CD$hybrid_score[match(rownames(SeuratObject@meta.data), rownames(CD))]
    SeuratObject@meta.data$DoubletSCDShybrid <- CD$DoubletSCDShybrid[match(rownames(SeuratObject@meta.data), rownames(CD))]
    head(SeuratObject@meta.data)
    
    print(table(SeuratObject@meta.data$DoubletSCDShybrid))
    
    head(SeuratObject@meta.data)
    cutoff.df <- data.frame(Doublets = table(SeuratObject@meta.data$DoubletSCDShybrid)); print(cutoff.df)
    TableDF <- cutoff.df
    FontsDF <- c(2.5,2.5,2.5)
    titleDF <- paste0(Sample,": Doublets Detected")
    func.PlotTable.General(TableDF, FontsDF, titleDF, 20)
    
    print(Create_Table(SeuratObject, BeforeFilter=0))
    
    print("Plotted Table")
    
    #p1 <- DimPlot(object = SeuratObject, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, cols = ClusPallette)
    #p2 <- DimPlot(object = SeuratObject, reduction = "umap", group.by = "SCDS", pt.size=1.5, cols=c("red", "dodgerblue", "black"))
    #print(plot_grid(p1, p2, NULL))
    dev.off()
    
    write.table(SeuratObject@meta.data,file=paste0("Doublets_Detected_",Sample,"_using_PCA_",PCAnum,"_res_",resClus,".txt"),quote=F,sep="\t")
    
    DoubletCells <- SeuratObject@meta.data[SeuratObject@meta.data$DoubletSCDShybrid != "Singlet",]; dim(DoubletCells)
    head(DoubletCells)
    write.table(DoubletCells,file=paste0("Discarded_Cells_",OutName,Sample,"_using_PCA_",PCAnum,"_res_",resClus,".txt"),quote=F,sep="\t")
    
  }
  
  
  
  #setwd(plotWD)
  setwd(SCDSDir)
  DoubletfromSCDShybrid <- read.table(file = paste0("Discarded_Cells_", 
                                                    OutName, Sample, "_using_PCA_", PCAnum, "_res_", 
                                                    resClus, ".txt"))
  head(DoubletfromSCDShybrid)
  dim(DoubletfromSCDShybrid)
  
  common <- rownames(DoubletfromSCDShybrid)
  length(common)
  setwd(DDdir)
  SeuratObject@meta.data$DoubletfromSCDShybrid <- "Singlet"
  SeuratObject@meta.data[rownames(DoubletfromSCDShybrid), "DoubletfromSCDShybrid"] <- "Doublet"
  print(table(SeuratObject@meta.data$DoubletfromSCDShybrid, SeuratObject@meta.data$DoubletSCDShybrid))
  write.table(SeuratObject@meta.data, file = paste0("Comparison_SCDS_Doublets_Detected_", 
                                                    Sample, "_using_PCA_", PCAnum, "_res_", resClus, 
                                                    ".txt"), quote = F, sep = "\t")
  
  print("Doing Final Plotting")
  
  setwd(DDdir)
  pdf(file = paste0("Comparison_Plots_Doublets_Detected_", 
                    Sample, "_using_PCA_", PCAnum, "_res_", resClus, 
                    ".pdf"), height = 10, width = 12)
  table(SeuratObject@meta.data$DoubletfromSCDShybrid)
  p1 <- DimPlot(SeuratObject, label.size = 6, group.by = "DoubletfromSCDShybrid") + 
    ggtitle(paste0("DoubletfromSCDShybrid"))
  table(SeuratObject@meta.data$DoubletfromSCDShybrid)
  
  
  if(plotCCgene==TRUE){
    if(Species=="hsa"){
      p2 <- FeaturePlot(SeuratObject, features = "TOP2A")
    } else {
      p2 <- FeaturePlot(SeuratObject, features = "Top2a")
    }
    
    
    print(plot_grid(p1, p2, NULL, NULL, ncol = 2))
  } else if(plotCCgene==FALSE){
    
    print(plot_grid(p1, NULL, NULL, NULL, ncol = 2))
  }
  
  head(SeuratObject@meta.data)
  cutoff.df <- data.frame(Doublets = table(SeuratObject@meta.data$DoubletSCDShybrid)); print(cutoff.df)
  TableDF <- cutoff.df
  FontsDF <- c(2.5,2.5,2.5)
  titleDF <- paste0(Sample,": Doublets Detected")
  func.PlotTable.General(TableDF, FontsDF, titleDF, 20)
  
  Create_Table(SeuratObject, BeforeFilter=0)
  
  dev.off()
  
  print("Done")
  print(Sys.time())
}



