#' A Perform_CCA_Diff_PCAdim_res Function
#'
#' This function allows you to express your love of cats.
#' @param Sampleall.object A list of Seurat objects between which to find anchors for downstream integration.
#' @param NameInpdf Path to save Quality plots and RDS data.
#' @param saveDIR Suffix to be added to CCA files.
#' @param FeatureUseCount A numeric value. This will call SelectIntegrationFeatures to select the provided number of features to be used in anchor finding
#' @param plots Save CCA plots
#' @param save Save integrated CCA RDS Seurat object
#' @keywords Sampleall.object, NameInpdf, saveDIR, FeatureUseCount, TempCCAdimchosen, res, plots, save
#' @export
#' @examples
#' Perform_CCA_Diff_PCAdim_res()



Perform_CCA_Diff_PCAdim_res <- function(Sampleall.object, NameInpdf, saveDIR, TempCCAdimchosenlist=c(15, 20, 30), FeatureUseCount=2500, plots = TRUE, save = FALSE){
  
  TempAll.object=Sampleall.object
  print(paste0("Performing CCA for ",length(TempAll.object)," Objects"))
  
  pdf(file=paste0("Plots_CCA_",Sample,"_Different_PCs_and_Resolutions.pdf"),height = 16,width = 18)
  print(ElbowPlot(sample.integrated, ndims = 30) + ggtitle(paste((unlist(TempCCAdimchosenlist)), collapse=",")))
  
  for(TempCCAdimchosen in TempCCAdimchosenlist){
  #TempCCAdimchosen=20
  print(paste0("Processing for Temp PCA ",TempCCAdimchosen))
    
  reference.list <- c(TempAll.object)
  print(paste0("Finding IntegrationAnchors"))
  sample.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:TempCCAdimchosen, anchor.features = FeatureUseCount)
  
  ## We then pass these anchors to the IntegrateData function, which returns a Seurat object.
  ## The returned object will contain a new Assay, which holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
  print(paste0("Integrating Data"))
  sample.integrated <- IntegrateData(anchorset = sample.anchors, dims = 1:TempCCAdimchosen)
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
  sample.integrated <- RunPCA(object = sample.integrated, npcs = TempCCAdimchosen, verbose = FALSE)
  sample.integrated <- RunUMAP(object = sample.integrated, reduction = "pca", dims = 1:TempCCAdimchosen)
  
  setwd(saveDIR)
  p=q=list()
  
  Tempreslist=c("0.1", "0.2", "0.3", "0.5", "0.8")
  for(Tempres in Tempreslist){
  #for(Tempres in c("0.1", "0.2", "0.3")){
  #Tempres="0.1"
    
  print(paste0("Processing for Temp resolution ",Tempres))
  sample.integrated <- FindNeighbors(sample.integrated, dims = 1:TempCCAdimchosen)
  sample.integrated <- FindClusters(sample.integrated, resolution = as.numeric(Tempres))
  
  
  if (plots == TRUE) {
  print(paste0("Generating different resolution plots for PCA ",TempCCAdimchosen))
    
  p[[Tempres]]  <- DimPlot(sample.integrated, reduction = "umap", group.by = "seurat_clusters")
  q[[Tempres]]  <- DimPlot(sample.integrated, reduction = "umap", group.by = GroupName)
  #print(plot_grid(p1, p2, NULL, NULL, nrow = 2))
  }
  
  }
  
  print(plot_grid(q[["0.1"]]))
  #print(plot_grid(p[["0.1"]], p[["0.2"]], p[["0.3"]], p[["0.4"]], p[["0.5"]], p[["0.6"]], p[["0.7"]], p[["0.8"]], p[["0.9"]], ncol = 3),
  blank <- ggplot() + theme_bw() + ggtitle(paste0("PCs ",TempCCAdimchosen)) + theme(plot.title = element_text(size=50, colour = "red"))
  print(plot_grid(p[["0.1"]], p[["0.2"]], p[["0.3"]], p[["0.5"]], p[["0.8"]], blank, ncol = 3, labels=c("0.1", "0.2", "0.3", "0.5", "0.8")))
  
  
  if (save == TRUE) {
    print("Saving Seurat RDS object and meta data")
    setwd(saveDIR)
    write.table(sample.integrated@meta.data,file=paste0("Meta_Data_",NameInpdf,".txt"),quote=F,sep="\t")
    saveRDS(object = sample.integrated, file = paste0(NameInpdf,"_PCA",TempCCAdimchosen,".rds"))
  }
  
  }
  
  dev.off()
  
  print("Done")
  print(Sys.time())
}

