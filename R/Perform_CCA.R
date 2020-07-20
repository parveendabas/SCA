#' A Perform_CCA Function
#'
#' This function allows you to express your love of cats.
#' @param Sampleall.object A list of Seurat objects between which to find anchors for downstream integration.
#' @param NameInpdf Path to save Quality plots and RDS data.
#' @param saveDIR Suffix to be added to CCA files.
#' @param FeatureUseCount A numeric value. This will call SelectIntegrationFeatures to select the provided number of features to be used in anchor finding
#' @param CCAdimchosen Which dimensions to use from the CCA to specify the neighbor search space
#' @param res Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.
#' @param plots Save CCA plots
#' @param save Save integrated CCA RDS Seurat object
#' @keywords Sampleall.object, NameInpdf, saveDIR, FeatureUseCount, CCAdimchosen, res, plots, save
#' @export
#' @examples
#' Perform_CCA()



Perform_CCA <- function(Sampleall.object, NameInpdf, saveDIR, FeatureUseCount=2500, CCAdimchosen=30, res = 0.5, plots = TRUE, save = TRUE){
  
  TempAll.object=Sampleall.object
  
  print(paste0("Performing CCA for ",length(TempAll.object)," Objects"))
  
  reference.list <- c(TempAll.object)
  print(paste0("Finding IntegrationAnchors"))
  sample.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:CCAdimchosen, anchor.features = FeatureUseCount)
  
  ## We then pass these anchors to the IntegrateData function, which returns a Seurat object.
  ## The returned object will contain a new Assay, which holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
  print(paste0("Integrating Data"))
  sample.integrated <- IntegrateData(anchorset = sample.anchors, dims = 1:CCAdimchosen)
  ##After running IntegrateData, the Seurat object will contain a new Assay with the integrated expression matrix. Note that the original (uncorrected values) 
  ## are still stored in the object in the “RNA” assay, so you can switch back and forth.
  
  #saveRDS(object = sample.integrated, file = paste0(NameInpdf,".rds"))
  
  ## We can then use this new integrated matrix for downstream analysis and visualization. Here we scale the integrated data, run PCA, and visualize the results with UMAP. 
  ## The integrated datasets cluster by cell type, instead of by technology.
  # switch to integrated assay. The variable features of this assay are automatically set during
  # IntegrateData
  DefaultAssay(object = sample.integrated) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  sample.integrated <- ScaleData(object = sample.integrated, verbose = FALSE)
  sample.integrated <- RunPCA(object = sample.integrated, npcs = CCAdimchosen, verbose = FALSE)
  sample.integrated <- RunUMAP(object = sample.integrated, reduction = "pca", dims = 1:CCAdimchosen)
  
  sample.integrated <- FindNeighbors(sample.integrated, dims = 1:CCAdimchosen)
  sample.integrated <- FindClusters(sample.integrated, resolution = res)
  
  setwd(saveDIR)
  if (plots == TRUE) {
    print("Generating quality plots")
    pdf(file=paste0("Plots_CCA_",Sample,".pdf"),height = 10,width = 12)
    print(ElbowPlot(sample.integrated, ndims = CCAdimchosen) + ggtitle(paste("PCs ",CCAdimchosen)))
    p1 <- DimPlot(sample.integrated, reduction = "umap", group.by = GroupName)
    p2 <- DimPlot(sample.integrated, reduction = "umap", group.by = "seurat_clusters")
    print(plot_grid(p1, p2, NULL, NULL, nrow = 2))
    print(DimPlot(sample.integrated, reduction = "umap", group.by = "seurat_clusters", split.by = GroupName, ncol = 3))
    print(DimPlot(sample.integrated, reduction = "umap", group.by = GroupName, split.by = "seurat_clusters", ncol = 5))
    Create_Table(sample.integrated)
    print(VlnPlot(SCdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), pt.size = 0.5, ncol = 2)) 
    dev.off()
  }
  
  
  if (save == TRUE) {
  print("Saving Seurat RDS object and meta data")
  setwd(saveDIR)
  write.table(sample.integrated@meta.data,file=paste0("Meta_Data_",NameInpdf,".txt"),quote=F,sep="\t")
  saveRDS(object = sample.integrated, file = paste0(NameInpdf,".rds"))
  }
  
  return(sample.integrated)
  rm(TempAll.object)
  
  print("Done")
  print(Sys.time())
}

