#' A Test_Cluster_Resolution_Harmony Function
#'
#' This function allows you to compare different cluster resolutions.
#' @param Temp.object Seurat objects to be used for the QC plots.
#' @param saveDIR Directory to save the plots.
#' @param IdentToBatchCorrect Identity to be used for Harmony batch correction
#' @param ThetaToBatchCorrect Theta values to be used while running harmony. Diversity clustering penalty parameter. Specify for each variable in group.by.vars. Default theta=2. theta=0 does not encourage any diversity. Larger values of theta result in more diverse clusters.
#' @param SuffixName Suffix. to be added in the directory name as 
#' @param ColToPlot Mention one ident only
#' @param ColPaletteToPlot Mention one palette only for the ident
#' @param mingenes minimum gene number that will be mentioned in the output file name
#' @param PCAuse PCA dim to use
#' @param ClusResList Vector of the PCAs to be used for testing
#' @param DownSamplePCA #Cells to plot in UMAP for PCA testing
#' @keywords Temp.object, saveDIR, IdentToBatchCorrect, ThetaToBatchCorrect, SuffixName, ColToPlot, ColPaletteToPlot mingenes PCAuse ClusResList
#' @export
#' @examples
#' Test_Cluster_Resolution_Harmony()



Test_Cluster_Resolution_Harmony <- function(Temp.object, saveDIR, IdentToBatchCorrect="orig.ident", ThetaToBatchCorrect=2, SuffixName="Testing_Cluster_Resolution",  
                                 ColToPlot=ColToPlot, ColPaletteToPlot=ColPaletteToPlot,mingenes=500, PCAuse=25,
                                 ClusResList=c(0.1, 0.15, 0.2, 0.3), ClusOrder = ClusOrderFrom1, DownSamplePCA=1000){
  
  
  print(paste0("Covariates being used:",IdentToBatchCorrect))
  print(paste0("Theta for Covariates being used:",ThetaToBatchCorrect))
  print(paste0("PCA being used :",PCAuse))
  
  print(paste0("Testing Cluster Resolutions:"))
  print(paste0(ClusResList))
  LastClusRes=rev(ClusResList)[1]
  print(paste0("Last ClusRes to test is:",LastClusRes))
  ClusResinpdf = paste(ClusResList, collapse = '_'); ClusResinpdf
  
  #Temp.object=SCdata
  #saveDIR=plotWD.Subset
  #IdentToBatchCorrect=c("DietLibrary")
  #ThetaToBatchCorrect=THETAtest
  #SuffixName=paste0(OutputName,"_",SubsetName,"_theta_",THETAinpdf)
  #ColToPlot=ColToPlot
  #ColPaletteToPlot=ColPaletteToPlot
  #mingenes=500
  #PCAuse=10
  #ClusResList = c(0.1, 0.15, 0.2, 0.3)
  #ClusOrder = ClusOrderFrom1
  #DownSamplePCA=1000
  
  ClusOrder <- ClusOrderFrom1
  
  
  setwd(saveDIR) 
  for(PCAdim in PCAuse){
    
    #PCAdim=10
    
    Temp.object <- RunPCA(object = Temp.object, npcs = PCAdim, verbose = FALSE)
    Temp.object <- RunHarmony(Temp.object, group.by.vars = IdentToBatchCorrect, theta=ThetaToBatchCorrect)
    Temp.object <- RunUMAP(Temp.object, reduction = "harmony", dims = 1:PCAdim)
    
    pdf(file=paste0("Testing_Cluster_Resolution_",SuffixName,"_minGenes_",mingenes,"_PCA_",PCAuse,"_ClusRes_",ClusResinpdf,".pdf"),height = 9,width = 15)
    
    ClusRes.list = list()
    for(resUse in ClusResList){
      #resUse=0.1
      print(paste0("Processing res:",resUse))
      Temp.object <- FindNeighbors(Temp.object, reduction = "harmony", dims = 1:PCAdim)
      Temp.object <- FindClusters(Temp.object, resolution=resUse)
      Temp.object$seurat_clusters <- MakeClustersFrom1(Temp.object$seurat_clusters, ClusOrder)
      
      #sort(unique(Temp.object$seurat_clusters))
      #paste0(ToUseColumn,"_Based")
      table(Temp.object@meta.data[,paste0("RNA_snn_res.",resUse)])
      Temp.object@meta.data[,paste0("RNA_snn_res.",resUse)] <- Temp.object$seurat_clusters
      table(Temp.object@meta.data[,paste0("RNA_snn_res.",resUse)])
      
      # Projecting singlet identities on TSNE visualization
      #DimPlot(Temp.object, group.by = "HTO_classification")
      Temp.object@meta.data$seurat_clusters <- MakeClustersFrom1(Temp.object@meta.data$seurat_clusters, ClusOrder)
      Temp.object.DS <- subset(Temp.object,downsample=DownSamplePCA)
      ClusRes.list[[as.character(resUse)]] <- DimPlot(Temp.object.DS, group.by = "seurat_clusters", pt.size = 0.5, reduction = "umap", cols = ClusPallette, label = T, label.size = 8) + ggtitle(paste0("Res:",resUse, ", PCA:",PCAdim))
      
    }
    
    
    ClusRescols=ceiling(length(ClusResList)/2); ClusRescols
    print(paste0("Columns being plotted for ClusRes: ",ClusRescols))
    
    print(plot_grid(plotlist = ClusRes.list, ncol = ClusRescols))
    
    
    Idents(Temp.object) <- ColToPlot
    q1 <- DimPlot(Temp.object, pt.size = 0.5, reduction = "umap", cols = ColPaletteToPlot, 
                  label = F, label.size = 7) + ggtitle(paste0("PCA:", PCAdim))
    
    TopPanel <- plot_grid(q1, ncol = ClusRescols)
    print(plot_grid(TopPanel, nrow = 2))
    dev.off()
  }
  
  return(Temp.object)
  
  
  print("Done")
  print(Sys.time())
}

