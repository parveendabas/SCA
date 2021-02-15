#' A Perform_CCA Function
#'
#' This function allows you to perform CCA. For the selected PC, different resoultion ("0.1", "0.2", "0.3", "0.5", "0.8") will be computed for clustering
#' @param TempAll.object A list of Seurat objects between which to find anchors for downstream integration.
#' @param NameInpdf Path to save Quality plots and RDS data.
#' @param saveDIR Suffix to be added to CCA files.
#' @param FeatureUseCount A numeric value. This will call SelectIntegrationFeatures to select the provided number of features to be used in anchor finding
#' @param CCAdimchosen Which dimensions to use from the CCA to specify the neighbor search space
#' @param res Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.
#' @param res k.filter	How many neighbors (k) to use when filtering anchors
#' @param plots Save CCA plots
#' @param save Save integrated CCA RDS Seurat object
#' @keywords TempAll.object, NameInpdf, saveDIR, FeatureUseCount, CCAdimchosen, res, k.filter, plots, save
#' @export
#' @examples
#' Perform_CCA()



Perform_CCA <- function(TempAll.object, NameInpdf, saveDIR, FeatureUseCount=2500, CCAdimchosen=30, res = 0.5, k.filter = 200, plots = TRUE, save = TRUE){
  
  #TempAll.object=Sampleall.object
  
  print(paste0("Performing CCA for ",length(TempAll.object)," Objects"))
  
  reference.list <- c(TempAll.object)
  print(paste0("Finding IntegrationAnchors"))
  sample.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:CCAdimchosen, anchor.features = FeatureUseCount, k.filter = k.filter)
  
  ## We then pass these anchors to the IntegrateData function, which returns a Seurat object.
  ## The returned object will contain a new Assay, which holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
  print(paste0("Integrating Data"))
  temp.integrated <- IntegrateData(anchorset = sample.anchors, dims = 1:CCAdimchosen)
  ##After running IntegrateData, the Seurat object will contain a new Assay with the integrated expression matrix. Note that the original (uncorrected values) 
  ## are still stored in the object in the “RNA” assay, so you can switch back and forth.
  
  #saveRDS(object = temp.integrated, file = paste0(NameInpdf,".rds"))
  
  ## We can then use this new integrated matrix for downstream analysis and visualization. Here we scale the integrated data, run PCA, and visualize the results with UMAP. 
  ## The integrated datasets cluster by cell type, instead of by technology.
  # switch to integrated assay. The variable features of this assay are automatically set during
  # IntegrateData
  DefaultAssay(object = temp.integrated) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  temp.integrated <- ScaleData(object = temp.integrated, verbose = FALSE)
  temp.integrated <- RunPCA(object = temp.integrated, npcs = CCAdimchosen, verbose = FALSE)
  temp.integrated <- RunUMAP(object = temp.integrated, reduction = "pca", dims = 1:CCAdimchosen)
  
  p=q=list()
  temp.integrated1 <- temp.integrated
  Tempreslist=c("0.1", "0.2", "0.3", "0.5", "0.8")
  for(Tempres in Tempreslist){
    #for(Tempres in c("0.1", "0.2", "0.3")){
    #Tempres="0.1"
    
    print(paste0("Processing for Temp resolution ",Tempres))
    temp.integrated1 <- FindNeighbors(temp.integrated1, dims = 1:CCAdimchosen)
    temp.integrated1 <- FindClusters(temp.integrated1, resolution = as.numeric(Tempres))
    
    
    if (plots == TRUE) {
      print(paste0("Generating different resolution plots for PCA ",CCAdimchosen))
      
      p[[Tempres]]  <- DimPlot(temp.integrated1, reduction = "umap", group.by = "seurat_clusters", cols = ClusPalette, label = T, label.size = 5)
      q[[Tempres]]  <- DimPlot(temp.integrated1, reduction = "umap", group.by = "orig.ident", cols = Dark.Pallette)
      #print(plot_grid(p1, p2, NULL, NULL, nrow = 2))
    }
    
  }
  rm(temp.integrated1)
  
  pdf(file=paste0("Plots_CCA_",NameInpdf,"_PC",CCAdimchosen,"_Different_Resolutions.pdf"),height = 10,width = 14)
  print(ElbowPlot(temp.integrated, ndims = CCAdimchosen) + ggtitle(paste("PCs ",CCAdimchosen)))
  print(plot_grid(q[["0.1"]]))
  #print(plot_grid(p[["0.1"]], p[["0.2"]], p[["0.3"]], p[["0.4"]], p[["0.5"]], p[["0.6"]], p[["0.7"]], p[["0.8"]], p[["0.9"]], ncol = 3),
  blank <- ggplot() + theme_bw() + ggtitle(paste0("PCs ",CCAdimchosen)) + theme(plot.title = element_text(size=50, colour = "red"))
  print(plot_grid(p[["0.1"]], p[["0.2"]], p[["0.3"]], p[["0.5"]], p[["0.8"]], blank, ncol = 3, labels=c("0.1", "0.2", "0.3", "0.5", "0.8")))
  dev.off()
  
  temp.integrated <- FindNeighbors(temp.integrated, dims = 1:CCAdimchosen)
  temp.integrated <- FindClusters(temp.integrated, resolution = res)
  
  setwd(saveDIR)
  if (plots == TRUE) {
    print("Generating quality plots")
    pdf(file=paste0("Plots_CCA_",NameInpdf,".pdf"),height = 10,width = 14)
    print(ElbowPlot(temp.integrated, ndims = CCAdimchosen) + ggtitle(paste("PCs ",CCAdimchosen)))
    p1 <- DimPlot(temp.integrated, reduction = "umap", group.by = GroupName, cols = Dark.Pallette)
    p2 <- DimPlot(temp.integrated, reduction = "umap", group.by = "seurat_clusters", cols = ClusPalette, label = T, label.size = 5)
    print(plot_grid(p1, p2, NULL, NULL, nrow = 2))
    print(DimPlot(temp.integrated, reduction = "umap", group.by = "seurat_clusters", split.by = GroupName, ncol = 3, cols = ClusPalette))
    print(DimPlot(temp.integrated, reduction = "umap", group.by = GroupName, split.by = "seurat_clusters", ncol = 5, cols = Dark.Pallette))
    #Create_Table(temp.integrated)
    print(VlnPlot(temp.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), pt.size = 0.5, ncol = 2)) 
    dev.off()
  }
  
  
  if (save == TRUE) {
  print("Saving Seurat RDS object and meta data")
  setwd(saveDIR)
  write.table(temp.integrated@meta.data,file=paste0("Meta_Data_",NameInpdf,".txt"),quote=F,sep="\t")
  saveRDS(object = temp.integrated, file = paste0(NameInpdf,".rds"))
  }
  
  return(temp.integrated)
  rm(TempAll.object)
  
  print("Done")
  print(Sys.time())
}

