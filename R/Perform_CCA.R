#' A Perform_CCA Function
#'
#' This function allows you to express your love of cats.
#' @param matrix.DIR Path to 10X directory.
#' @param saveDIR Path to save Quality plots and RDS data.
#' @param Sample Sample Name.
#' @param mincells Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
#' @param mingenes Include cells where at least this many features are detected.
#' @param mtpercent Include cells reporting at most this much mitochondrial transcript percentage.
#' @param rbpercent Include cells reporting at most this much ribosomal transcript percentage.
#' @param FeatureUseCount Number of features to select as top variable features; only used when selection.method is set to 'dispersion' or 'vst'
#' @keywords matrix.DIR, saveDIR, Sample, mincells, mingenes, mtpercent, rbpercent, FeatureUseCount
#' @export
#' @examples
#' Perform_CCA()



Perform_CCA <- function(Sampleall.object, FeatureUseCount=2500, plots = TRUE, save = TRUE){
  
  print(paste0("Performing CCA",Sample))
  print("Reading 10X Dir:")
  print(matrix.DIR)
  print(saveDIR)
  data <- Read10X(data.dir = matrix.DIR)
  ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
  
  print("Creating Seurat Object")
  print(paste0("min.cells:",mincells, ", min.features:",mingenes))
  SCdata <- CreateSeuratObject(counts = data, project = Sample, min.cells = mincells, min.features = mingenes)
  SCdata@meta.data$Cells <- rownames(SCdata@meta.data)
  SCdata@meta.data$Project <- Sample
  head(SCdata@meta.data)
  
  reference.list <- c(Sampleall.object)
  sample.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:CCAdimchosen, anchor.features = FeatureUseInCCA)
  
  ## We then pass these anchors to the IntegrateData function, which returns a Seurat object.
  ## The returned object will contain a new Assay, which holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
  sample.integrated <- IntegrateData(anchorset = sample.anchors, dims = 1:CCAdimchosen)
  ##After running IntegrateData, the Seurat object will contain a new Assay with the integrated expression matrix. Note that the original (uncorrected values) 
  ## are still stored in the object in the “RNA” assay, so you can switch back and forth.
  
  #saveRDS(object = sample.integrated, file = paste0(OutName,".rds"))
  
  ## We can then use this new integrated matrix for downstream analysis and visualization. Here we scale the integrated data, run PCA, and visualize the results with UMAP. 
  ## The integrated datasets cluster by cell type, instead of by technology.
  # switch to integrated assay. The variable features of this assay are automatically set during
  # IntegrateData
  DefaultAssay(object = sample.integrated) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  sample.integrated <- ScaleData(object = sample.integrated, verbose = FALSE)
  sample.integrated <- RunPCA(object = sample.integrated, npcs = CCAdimchosen, verbose = FALSE)
  sample.integrated <- RunUMAP(object = sample.integrated, reduction = "pca", dims = 1:CCAdimchosen)
  saveRDS(object = sample.integrated, file = paste0(OutName,".rds"))
  
  if (save == TRUE) {
  print("Saving Seurat RDS object and meta data")
  setwd(saveDIR)
  write.table(SCdata@meta.data,file=paste0("Meta_Data_",Sample,".txt"),quote=F,sep="\t")
  saveRDS(SCdata, file = paste0(Sample,".rds"))
  }
  
  return(SCdata)
  
  print("Done")
  print(Sys.time())
}

