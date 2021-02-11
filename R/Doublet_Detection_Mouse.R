#' A Read10X_Norm_Mouse Function
#'
#' This function allows you to express your love of cats.
#' @param SeuratObject Seurat Object.
#' @param saveDIR Path to save Quality plots and RDS data.
#' @param Sample Sample Name.
#' @param FeatureUseCount Number of features to select as top variable features; only used when selection.method is set to 'dispersion' or 'vst'
#' @param PCAnum Number of PCs to be used
#' @param resClus Resolution to be used for clustering
#' @param ClusPallette Color pallete to be used for clusters
#' @param DoubletFinder Run DoubletFinder
#' @param DoubletDecon Run DoubletDecon
#' @keywords SeuratObject, saveDIR, Sample, FeatureUseCount, PCAnum, resClus, ClusPallette, DoubletFinder, DoubletDecon
#' @export
#' @examples
#' Read10X_Norm_Mouse()



Doublet_Detection_Mouse <- function(SeuratObject, saveDIR, Sample, FeatureUseCount=2500, PCAnum=10, resClus = 0.5, ClusPallette=ClusPallette, DoubletFinder = TRUE, DoubletDecon = TRUE){
  
  #SeuratObject=TempSCdata
  #saveDIR=saveDIR
  #Sample=Sample
  
  library(DoubletDecon)
  library(DoubletFinder)
  
  setwd(saveDIR)
  DDdir <- paste(getwd(),paste0("Doublet_Detection_",Sample),sep="/"); print(DDdir)
  dir.create(file.path(getwd(),paste0("Doublet_Detection_",Sample)), showWarnings = FALSE)
  
  print(paste0("PrePorocessing data before running Doublet Detection Algorithms for: ",Sample))
  RUNProcessData="YES"
  if(RUNProcessData == "YES"){
    
    setwd(DDdir)
    pdf(file=paste0("PreProcess_Doublet_Detection_",Sample,"_using_PCA_",PCAnum,"_res_",resClus,".pdf"),height = 10,width = 12)
    print(paste0("Loading Seurat object for: ",name))
    
    ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
    #SeuratObject <- NormalizeData(SeuratObject)
    # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
    SeuratObject[["percent.mt"]] <- PercentageFeatureSet(SeuratObject, pattern = "^mt-")
    # Visualize QC metrics as a violin plot
    p1 <- VlnPlot(SeuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), cols = ClusPallette, pt.size = 0.00, ncol = 1)
    print(p1)
    plot1 <- FeatureScatter(SeuratObject, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = ClusPallette)
    plot2 <- FeatureScatter(SeuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = ClusPallette)
    print(CombinePlots(plots = list(plot1, plot2)))
    
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
    p1 <- VlnPlot(SeuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), cols = ClusPallette, pt.size = 0.00, ncol = 1)
    print(p1)
    
    SeuratObject <- RunUMAP(SeuratObject, dims = 1:PCAnum)
    print(DimPlot(SeuratObject, reduction = "umap", label=TRUE, label.size = 8, pt.size = 0.5, cols = ClusPallette) + ggtitle(paste0(Sample, " (",nrow(SeuratObject@meta.data)," cells)")))
    dev.off()
    
  }
  
  
  if (DoubletFinder == TRUE) {
    
    setwd(DDdir)
    DoubletFinderDir <- paste(getwd(),paste0("DoubletFinder_",Sample),sep="/"); print(DoubletFinderDir)
    dir.create(file.path(getwd(),paste0("DoubletFinder_",Sample)), showWarnings = FALSE)
    
    setwd(DoubletFinderDir)
    pdf(file=paste0("Plots_DoubletFinder_",Sample,"_using_PCA_",PCAnum,"_res_",resClus,".pdf"),height = 10,width = 12)
    
    print(paste0("Started DoubletFinder process for ",Sample, " using PCA ",PCAnum))
    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.list_SeuratObject <- paramSweep_v3(SeuratObject, PCs = 1:PCAnum, sct = FALSE)
    sweep.stats_SeuratObject <- summarizeSweep(sweep.res.list_SeuratObject, GT = FALSE)
    bcmvn_SeuratObject <- find.pK(sweep.stats_SeuratObject)
    
    
    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    homotypic.prop <- modelHomotypic(SeuratObject$seurat_clusters)           ## ex: annotations <- SeuratObject@meta.data$ClusteringResults
    nExp_poi <- round(0.075*length(colnames(x = SeuratObject))); nExp_poi  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)); nExp_poi.adj
    
    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    SeuratObject <- doubletFinder_v3(SeuratObject, PCs = 1:PCAnum, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    head(SeuratObject@meta.data)
    colnames(SeuratObject@meta.data)[colnames(SeuratObject@meta.data) %like% 'pANN_'] <- "pANNcomputed"
    colnames(SeuratObject@meta.data)[colnames(SeuratObject@meta.data) %like% 'DF.classifications_'] <- "Doublet_LowConf"
    head(SeuratObject@meta.data)
    
    SeuratObject <- doubletFinder_v3(SeuratObject, PCs = 1:PCAnum, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANNcomputed", sct = FALSE)
    head(SeuratObject@meta.data)
    colnames(SeuratObject@meta.data)[colnames(SeuratObject@meta.data) %like% 'DF.classifications_'] <- "Doublet_HighConf"
    head(SeuratObject@meta.data)
    table(SeuratObject@meta.data$Doublet_LowConf)
    table(SeuratObject@meta.data$Doublet_HighConf)
    table(SeuratObject@meta.data$Doublet_LowConf, SeuratObject@meta.data$Doublet_HighConf)
    
    SeuratObject@meta.data$DoubletFinder <- SeuratObject@meta.data$Doublet_LowConf
    SeuratObject@meta.data[SeuratObject@meta.data$Doublet_LowConf == "Doublet", "DoubletFinder"] <- "Doublet_LowConf"
    table(SeuratObject@meta.data$DoubletFinder)
    SeuratObject@meta.data[SeuratObject@meta.data$Doublet_HighConf == "Doublet", "DoubletFinder"] <- "Doublet_HighConf"
    print(table(SeuratObject@meta.data$DoubletFinder))
    
    head(SeuratObject@meta.data)
    cutoff.df <- data.frame(Doublets = table(SeuratObject@meta.data$DoubletFinder)); cutoff.df
    TableDF <- cutoff.df
    FontsDF <- c(2.5,2.5,2.5)
    titleDF <- paste0(Sample,": Doublets Detected")
    func.PlotTable.General(TableDF, FontsDF, titleDF, 20)
    
    Create_Table(SeuratObject, BeforeFilter=0)
    
    
    p1 <- DimPlot(object = SeuratObject, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, cols = ClusPallette)
    p2 <- DimPlot(object = SeuratObject, reduction = "umap", group.by = "DoubletFinder", pt.size=1.5, cols=c("red", "dodgerblue", "black"))
    print(plot_grid(p1, p2, NULL))
    dev.off()
    
    write.table(SeuratObject@meta.data,file=paste0("Doublets_Detected_",Sample,"_using_PCA_",PCAnum,"_res_",resClus,".txt"),quote=F,sep="\t")
    
    DoubletCells <- SeuratObject@meta.data[SeuratObject@meta.data$DoubletFinder != "Singlet",]; dim(DoubletCells)
    head(DoubletCells)
    write.table(DoubletCells,file=paste0("Discarded_Cells_",OutName,Sample,"_using_PCA_",PCAnum,"_res_",resClus,".txt"),quote=F,sep="\t")
    
  }
  
  
  
  if (DoubletDecon == TRUE) {
    
    setwd(DDdir)
    DoubletDeconDir <- paste(getwd(),paste0("DoubletDecon_",Sample),sep="/"); print(DoubletDeconDir)
    dir.create(file.path(getwd(),paste0("DoubletDecon_",Sample)), showWarnings = FALSE)
    
    setwd(DoubletDeconDir)
    pdf(file=paste0("Plots_DoubletDecon_",Sample,"_using_PCA_",PCAnum,"_res_",resClus,".pdf"),height = 10,width = 12)
    
    #Seurat_Pre_Process(expressionFile, genesFile, clustersFile)
    #expressionFile: Normalized expression matrix or counts file as a .txt file (expression from Seurat's NormalizeData() function)
    #genesFile: Top marker gene list as a .txt file from Seurat's top_n() function
    #clustersFile: Cluster identities as a .txt file from Seurat object @ident
    
    ## find markers for every cluster compared to all remaining cells, report only the positive ones
    Idents(object = SeuratObject) <- "seurat_clusters"
    #SeuratObject <- subset(SeuratObject, downsample=100)
    #SeuratObject.markers <- FindAllMarkers(SeuratObject, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    #top10 <- SeuratObject.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    ##DoHeatmap(SeuratObject, features = top10$gene) + NoLegend()
    #top50 <- SeuratObject.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
    #write.table(top50,paste0("Top50_DGEs_",OutName,NameInpdf,"_Doublets_Detected_",Sample,".txt"),quote=F,sep="\t")
    print(paste0("Completed Marker identification for sample ",Sample, " using PCA ",PCAnum))
    
    newFiles <- Improved_Seurat_Pre_Process(SeuratObject, num_genes=50, write_files=FALSE)
    
    filename=paste0(Sample,"_PCA",PCAnum)
    setwd(DoubletDeconDir)
    write.table(newFiles$newExpressionFile, paste0(filename, "_expression.txt"), sep="\t")
    #write.table(newFiles$newFullExpressionFile, paste0(filename, "_fullExpression.txt"), sep="\t")
    write.table(newFiles$newGroupsFile, paste0(filename , "_groups.txt"), sep="\t", col.names = F)
    #closeAllConnections()
    
    results=Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile, 
                               groupsFile=newFiles$newGroupsFile, 
                               filename=filename, 
                               location=DoubletDeconDir,
                               fullDataFile=NULL, 
                               removeCC=FALSE, 
                               species="mmu", 
                               rhop=1.1, 
                               write=TRUE, 
                               PMF=TRUE, 
                               useFull=FALSE, 
                               heatmap=FALSE,
                               centroids=TRUE,
                               num_doubs=100, 
                               only50=FALSE,
                               min_uniq=4,
                               nCores=1)
    
    
    head(results$Final_doublets_groups); dim(results$Final_doublets_groups)
    
    SeuratObject@meta.data$DoubletDecon <- "Singlet"
    #SeuratObject@meta.data[rownames(results$Final_doublets_groups),"DoubletDecon"] <- "Doublet"
    SeuratObject@meta.data[gsub("\\.","-",rownames(results$Final_doublets_groups)),"DoubletDecon"] <- "Doublet"
    head(SeuratObject@meta.data)
    table(SeuratObject@meta.data$DoubletDecon)
    table(SeuratObject@meta.data$DoubletDecon, SeuratObject@meta.data$seurat_clusters)
    
    temp <- results$DRS_doublet_table
    rownames(temp) <- gsub("\\.","-",rownames(temp))
    temp$DoubletDecon <- temp$isADoublet
    temp$DoubletDecon <- gsub("FALSE","Singlet",temp$DoubletDecon)
    temp$DoubletDecon <- gsub("TRUE","Doublet",temp$DoubletDecon)
    
    write.table(temp,file=paste0("Doublets_Detected_",Sample,"_using_PCA_",PCAnum,"_res_",resClus,".txt"),quote=F,sep="\t")
    DoubletCells <- temp[temp$DoubletDecon != "Singlet",]; dim(DoubletCells)
    head(DoubletCells)
    write.table(DoubletCells,file=paste0("Discarded_Cells_",OutName,Sample,"_using_PCA_",PCAnum,"_res_",resClus,".txt"),quote=F,sep="\t")
    

    TableDF <- as.data.frame.matrix(table(SeuratObject@meta.data$seurat_clusters, SeuratObject@meta.data$DoubletDecon))
    TableDF$DoubPerc <- round((TableDF$Doublet/(TableDF$Doublet+TableDF$Singlet)*100),2)
    FontsDF <- c(2.1,2.1,2.1)
    titleDF <- paste0(Sample)
    func.PlotTable.withRowNames.General(TableDF, FontsDF, titleDF, 15)
    
    p1 <- DimPlot(SeuratObject, label = TRUE, label.size = 6, cols = ClusPallette)
    p2 <- DimPlot(SeuratObject, label = TRUE, label.size = 5, group.by = "DoubletDecon")
    p3 <- FeaturePlot(SeuratObject, features = "Ooep")
    p4 <- FeaturePlot(SeuratObject, features = "Zp3")
    print(plot_grid(p1, p2, p3, p4, ncol = 2))
    dev.off()
    
    
    
  }
  
  
  CompareAlgos="YES"
  if(CompareAlgos == "YES"){
    
    #setwd(plotWD)
    setwd(DoubletFinderDir)
    DoubletfromFinder <- read.table(file = paste0("Discarded_Cells_",OutName,Sample,"_using_PCA_",PCAnum,"_res_",resClus,".txt")); head(DoubletfromFinder); dim(DoubletfromFinder)
    setwd(DoubletDeconDir)
    DoubletfromDecon <- read.table(file = paste0("Discarded_Cells_",OutName,Sample,"_using_PCA_",PCAnum,"_res_",resClus,".txt")); head(DoubletfromDecon); dim(DoubletfromDecon)
    head(DoubletfromFinder); dim(DoubletfromDecon)
    common <- intersect(rownames(DoubletfromFinder), rownames(DoubletfromDecon)); length(common)
    
    
    setwd(DDdir)
    SeuratObject@meta.data$DoubletfromFinder <- "Singlet"
    SeuratObject@meta.data[rownames(DoubletfromFinder),"DoubletfromFinder"] <- "Doublet"
    SeuratObject@meta.data$DoubletfromDecon <- "Singlet"
    SeuratObject@meta.data[rownames(DoubletfromDecon),"DoubletfromDecon"] <- "Doublet"
    SeuratObject@meta.data$DoubletDeconFinder <- "Singlet"
    SeuratObject@meta.data[common,"DoubletDeconFinder"] <- "Doublet"
    write.table(SeuratObject@meta.data,file=paste0("Comparison_DoubletDeconFinder_Doublets_Detected_",Sample,"_using_PCA_",PCAnum,"_res_",resClus,".txt"),quote=F,sep="\t")
    
    setwd(DDdir)
    pdf(file=paste0("Comparison_Plots_Doublets_Detected_",Sample,"_using_PCA_",PCAnum,"_res_",resClus,".pdf"),height = 10,width = 12)
    
    table(SeuratObject@meta.data$DoubletfromFinder)
    p1 <- DimPlot(SeuratObject, label.size = 6, group.by = "DoubletfromFinder") + ggtitle(paste0("DoubletfromFinder"))
    table(SeuratObject@meta.data$DoubletfromDecon)
    p2 <- DimPlot(SeuratObject, label.size = 6, group.by = "DoubletfromDecon") + ggtitle(paste0("DoubletfromDecon"))
    table(SeuratObject@meta.data$DoubletDeconFinder)
    p3 <- DimPlot(SeuratObject, label.size = 6, group.by = "DoubletDeconFinder") + ggtitle(paste0("Doublet_Decon_and_Finder"))
    p4 <- FeaturePlot(SeuratObject, features = "Ooep")
    print(plot_grid(p1, p2, p3, p4, ncol = 2))
    dev.off()
    
  }
  
  print("Done")
  print(Sys.time())
}

