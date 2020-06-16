#' A Read10X_Norm_Mouse Function
#'
#' This function allows you to express your love of cats.
#' @param matrix.DIR Path to 10X directory.
#' @param mincells Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
#' @param mingenes Include cells where at least this many features are detected.
#' @param mtpercent Include cells reporting at most this much mitochondrial transcript percentage.
#' @param rbpercent Include cells reporting at most this much ribosomal transcript percentage.
#' @param Sample Sample Name.
#' @keywords matrix.DIR, mincells, mingenes, mtpercent, rbpercent, Sample
#' @export
#' @examples
#' Read10X_Norm_Mouse()



Read10X_Norm_Mouse <- function(matrix.DIR, Sample, mincells=3, mingenes=500, mtpercent=20, rbpercent=50){
  
  print(paste0("Processing Sample:",Sample))
  print("Reading 10X Dir:")
  print(matrix.DIR)
  data <- Read10X(data.dir = matrix.DIR)
  ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
  
  print("Creating Seurat Object")
  print(paste0("min.cells:",mincells, ", min.features:",mingenes))
  SCdata <- CreateSeuratObject(counts = data, project = Sample, min.cells = mincells, min.features = mingenes)
  SCdata@meta.data$Cells <- rownames(SCdata@meta.data)
  SCdata@meta.data$Project <- Sample
  head(SCdata@meta.data)
  
  print("Counting MT and Ribosomal %")
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  mt.genes <- rownames(SCdata@assays$RNA@counts)[rownames(SCdata@assays$RNA@counts) %like% "^mt-"]; mt.genes; length(mt.genes)
  SCdata[["percent.mt"]] <- PercentageFeatureSet(SCdata, pattern = "^mt-")
  rb.genes <- rownames(SCdata@assays$RNA@counts)[rownames(SCdata@assays$RNA@counts) %like% "^Rp[sl]"]; rb.genes; length(rb.genes)
  SCdata[["percent.rb"]] <- PercentageFeatureSet(SCdata, pattern = "^Rp[sl]")
  print(VlnPlot(SCdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 2))  
  
  print("Filtering Data based on specified filters")
  print(paste0("mtpercent:",mtpercent, ", rbpercent:",rbpercent))
  SCdata <- subset(SCdata, subset = percent.mt < mtpercent & percent.rb < rbpercent)
  
  print("Normalizing Data")
  print("Normalizing Method: LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.")
  SCdata <- NormalizeData(SCdata)
  
  return(SCdata)
  
  print("Done")
  print(Sys.time())
}

