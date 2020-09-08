#' A Read10X_RemoveDoublets_Norm_Mouse Function
#'
#' This function allows you to express your love of cats.
#' @param matrix.DIR Path to 10X directory
#' @param saveDIR Path to save Quality plots and RDS data.
#' @param Sample Sample Name.
#' @param Doublet.DIR Path to Doublet.DIR directory (List of cells already called doublets)
#' @param mincells Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
#' @param mingenes Include cells where at least this many features are detected.
#' @param mtpercent Include cells reporting at most this much mitochondrial transcript percentage.
#' @param rbpercent Include cells reporting at most this much ribosomal transcript percentage.
#' @param FeatureUseCount Number of features to select as top variable features; only used when selection.method is set to 'dispersion' or 'vst'
#' @param PCAnum Number of PCs to be used
#' @param resClus Resolution to be used for clustering
#' @param plots Save QC plots
#' @param save Save RDS Seurat object
#' @keywords matrix.DIR, saveDIR, Sample, mincells, mingenes, mtpercent, rbpercent, FeatureUseCount, plots, save
#' @export
#' @examples
#' Read10X_RemoveDoublets_Norm_Mouse()



Read10X_RemoveDoublets_Norm_Mouse <- function(matrix.DIR, saveDIR, Sample, Doublet.DIR, mincells=3, mingenes=500, mtpercent=20, rbpercent=50, FeatureUseCount=2500, PCAnum=10, resClus = 0.5, plots = TRUE, save = TRUE){
  
  #matrix.DIR=matrix.DIR
  #Doublet.DIR=Doublet.DIR
  #saveDIR=saveDIR
  #Sample=Sample
  
  print(paste0("Processing Sample:",Sample))
  print("Reading 10X Dir:")
  print(matrix.DIR)
  print(saveDIR)
  data <- Read10X(data.dir = matrix.DIR)
  ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
  
  print("Creating Seurat Object")
  print(paste0("min.cells:",mincells, ", min.features:",mingenes))
  BeforeFilter=ncol(data); print(paste0(Sample, " cells BEFORE filtering:",BeforeFilter))
  SCdata <- CreateSeuratObject(counts = data, project = Sample, min.cells = mincells, min.features = mingenes)
  SCdata@meta.data$Cells <- rownames(SCdata@meta.data)
  SCdata@meta.data$Project <- Sample
  head(SCdata@meta.data)
  
  ## Remove Doublets
  print(Doublet.DIR); setwd(Doublet.DIR)
  Doubletfile <- dir(path = Doublet.DIR, pattern = paste0("Comparison_DoubletDeconFinder_Doublets_Detected_",Sample,"_using_PCA_",PCAnum,"_res_",resClus,".txt"), full.names = T, recursive = F)
  DoubletsDetected <- read.table(file = Doubletfile, header = T, sep = "\t"); head(DoubletsDetected); dim(DoubletsDetected)
  CellsToUse <- rownames(DoubletsDetected[DoubletsDetected$DoubletDeconFinder == "Singlet",]); length(CellsToUse)
  DiscardedDoublets <- rownames(DoubletsDetected[!DoubletsDetected$DoubletDeconFinder == "Singlet",]); length(DiscardedDoublets)
    
  SCdata1 <- SCdata
  SCdata <- subset(SCdata1, cells=CellsToUse); SCdata
  rm(SCdata1)
  
  print("Counting MT and Ribosomal %")
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  mt.genes <- rownames(SCdata@assays$RNA@counts)[rownames(SCdata@assays$RNA@counts) %like% "^mt-"]; mt.genes; length(mt.genes)
  SCdata[["percent.mt"]] <- PercentageFeatureSet(SCdata, pattern = "^mt-")
  rb.genes <- rownames(SCdata@assays$RNA@counts)[rownames(SCdata@assays$RNA@counts) %like% "^Rp[sl]"]; rb.genes; length(rb.genes)
  SCdata[["percent.rb"]] <- PercentageFeatureSet(SCdata, pattern = "^Rp[sl]")
  #print(VlnPlot(SCdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 2))  
  
  print("Filtering Data based on specified filters")
  print(paste0("mtpercent:",mtpercent, ", rbpercent:",rbpercent))
  SCdata <- subset(SCdata, subset = percent.mt < mtpercent & percent.rb < rbpercent)
  print(paste0(Sample, " cells AFTER filtering:",nrow(SCdata@meta.data)))
  
  print("Normalizing Data")
  print("Normalizing Method: LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.")
  SCdata <- NormalizeData(SCdata)
  SCdata <- FindVariableFeatures(object = SCdata, selection.method = "vst", nfeatures = FeatureUseCount, verbose = FALSE)
  
  if (plots == TRUE) {
  print("Generating quality plots")
    setwd(saveDIR)
    pdf(file=paste0("QC_",Sample,".pdf"),height = 8,width = 10)
    Create_Table_Doublets(SCdata, BeforeFilter=BeforeFilter, Doublets=length(DiscardedDoublets), SampleName=Sample, mincells=mincells, mingenes=mingenes, mtpercent=mtpercent, rbpercent=rbpercent)
    #Create_Table(SCdata, BeforeFilter=BeforeFilter, SampleName=Sample, mincells=mincells, mingenes=mingenes, mtpercent=mtpercent, rbpercent=rbpercent)
    print(VlnPlot(SCdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), pt.size = 0.5, ncol = 2)) 
    dev.off()
  } else {
      print("Skipping plotting table, plots == FALSE")
    }
  
  if (save == TRUE) {
  print("Saving Seurat RDS object and meta data")
  setwd(saveDIR)
  write.table(SCdata@meta.data,file=paste0("Meta_Data_",Sample,".txt"),quote=F,sep="\t")
  saveRDS(SCdata, file = paste0(Sample,".rds"))
  } else {
    print("Skipping saving RDS object, save == FALSE")
  }
  
  return(SCdata)
  
  print("Done")
  print(Sys.time())
}

