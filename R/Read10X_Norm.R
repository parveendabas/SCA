#' A Read10X_Norm Function
#'
#' This function allows you to express your love of cats.
#' @param matrix.DIR Path to 10X directory.
#' @param saveDIR Path to save Quality plots and RDS data.
#' @param Sample Sample Name.
#' @param Species Species Name. Valid options are hsa or mmu
#' @param mincells Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
#' @param mingenes Include cells where at least this many features are detected.
#' @param mtpercent Include cells reporting at most this much mitochondrial transcript percentage.
#' @param rbpercent Include cells reporting at most this much ribosomal transcript percentage.
#' @param FeatureUseCount Number of features to select as top variable features; only used when selection.method is set to 'dispersion' or 'vst'
#' @param plots Save QC plots
#' @param save Save RDS Seurat object
#' @keywords matrix.DIR, saveDIR, Sample, mincells, mingenes, mtpercent, rbpercent, FeatureUseCount, plots, save
#' @export
#' @examples
#' Read10X_Norm()



Read10X_Norm <- function(matrix.DIR, saveDIR, Sample, Species="hsa", mincells=3, mingenes=500, mtpercent=20, rbpercent=50, FeatureUseCount=2500, plots = TRUE, save = TRUE){
  
  print(paste0("Processing Sample:",Sample))
  print("Reading 10X Dir:")
  print(matrix.DIR)
  print(saveDIR)
  data <- Read10X(data.dir = matrix.DIR)
  ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
  
  print("Creating Seurat Object")
  print(paste0("min.cells:",mincells, ", min.features:",mingenes))
  BeforeFilter=ncol(data)
  print(paste0(Sample, " cells BEFORE filtering:",BeforeFilter))
  SCdata <- CreateSeuratObject(counts = data, project = Sample, min.cells = mincells, min.features = mingenes)
  SCdata@meta.data$Cells <- rownames(SCdata@meta.data)
  SCdata@meta.data$Project <- Sample
  head(SCdata@meta.data)
  
  if(Species=="hsa"){
  print("Counting MT and Ribosomal % for Human")
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  mt.genes <- rownames(SCdata@assays$RNA@counts)[rownames(SCdata@assays$RNA@counts) %like% "^MT-"]; mt.genes; length(mt.genes)
  SCdata[["percent.mt"]] <- PercentageFeatureSet(SCdata, pattern = "^MT-")
  rb.genes <- rownames(SCdata@assays$RNA@counts)[rownames(SCdata@assays$RNA@counts) %like% "^RP[SL]"]; rb.genes; length(rb.genes)
  SCdata[["percent.rb"]] <- PercentageFeatureSet(SCdata, pattern = "^RP[SL]")
  #print(VlnPlot(SCdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 2))  
  } else if (Species=="mmu"){
    print("Counting MT and Ribosomal % for Mouse")
    # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
    mt.genes <- rownames(SCdata@assays$RNA@counts)[rownames(SCdata@assays$RNA@counts) %like% "^mt-"]; mt.genes; length(mt.genes)
    SCdata[["percent.mt"]] <- PercentageFeatureSet(SCdata, pattern = "^mt-")
    rb.genes <- rownames(SCdata@assays$RNA@counts)[rownames(SCdata@assays$RNA@counts) %like% "^Rp[sl]"]; rb.genes; length(rb.genes)
    SCdata[["percent.rb"]] <- PercentageFeatureSet(SCdata, pattern = "^Rp[sl]")
    #print(VlnPlot(SCdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 2))  
  } else {
    print("Valid options for species are Human or Mouse only")
    break
  }
  
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
    Create_Table(SCdata, BeforeFilter=BeforeFilter, mincells=mincells, mingenes=mingenes, mtpercent=mtpercent, rbpercent=rbpercent)
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

