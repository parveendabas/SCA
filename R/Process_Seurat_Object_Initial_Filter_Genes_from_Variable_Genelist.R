#' A Process_Seurat_Object_Initial_Filter_Genes_from_Variable_Genelist Function
#'
#' This function performs a) mt, ribo % calculation, b) Normalization c) Filtering from variable genes d) Scaling of the data
#' @param Temp.object Seurat Object
#' @param FeatureNum number of Features for variable genes
#' @param Species Species Name. Valid options are hsa or mmu
#' @param mt mitochondrial genes to be removed from Seurat object 
#' @param rb ribosomal genes to be removed from Seurat object 
#' @param cc CellCycle genes (Converted from Human to Mouse in needed) to be removed from Seurat object 
#' @keywords Temp.object, mt, rb, cc
#' @export
#' @examples
#' Process_Seurat_Object_Initial_Filter_Genes_from_Variable_Genelist()

Process_Seurat_Object_Initial_Filter_Genes_from_Variable_Genelist <- function(Temp.object, FeatureNum=2500, Species="hsa", mt = TRUE, rb = TRUE, cc = FALSE){
  
  #Temp.object=Temp.object
  #saveDIR=pkWD
  #Species="mmu"
  #FeatureNum=2500
  
  if(Species=="hsa"){
    print("Counting MT and Ribosomal % for Human")
    # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
    mt.genes <- rownames(Temp.object@assays$RNA@counts)[rownames(Temp.object@assays$RNA@counts) %like% "^MT-"]; mt.genes; length(mt.genes)
    Temp.object[["percent.mt"]] <- PercentageFeatureSet(Temp.object, pattern = "^MT-")
    rb.genes <- rownames(Temp.object@assays$RNA@counts)[rownames(Temp.object@assays$RNA@counts) %like% "^RP[SL]"]; rb.genes; length(rb.genes)
    Temp.object[["percent.rb"]] <- PercentageFeatureSet(Temp.object, pattern = "^RP[SL]")
    #print(VlnPlot(Temp.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 2))  
  } else if (Species=="mmu"){
    print("Counting MT and Ribosomal % for Mouse")
    # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
    mt.genes <- rownames(Temp.object@assays$RNA@counts)[rownames(Temp.object@assays$RNA@counts) %like% "^mt-"]; mt.genes; length(mt.genes)
    Temp.object[["percent.mt"]] <- PercentageFeatureSet(Temp.object, pattern = "^mt-")
    rb.genes <- rownames(Temp.object@assays$RNA@counts)[rownames(Temp.object@assays$RNA@counts) %like% "^Rp[sl]"]; rb.genes; length(rb.genes)
    Temp.object[["percent.rb"]] <- PercentageFeatureSet(Temp.object, pattern = "^Rp[sl]")
    #print(VlnPlot(Temp.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 2))  
  } else {
    print("Valid options for species are Human or Mouse only")
    break
  }
  
  Temp.object <- NormalizeData(Temp.object) %>% FindVariableFeatures(nfeatures=FeatureNum)
  genes.filter.list <- Filter_Genes_from_Variable_Genelist(Temp.object, Species = Species, mt = mt, rb = rb, cc = cc)
  
  print(paste0("Variable genes length before filtering: ",length(Temp.object@assays$RNA@var.features)))
  length(intersect(Temp.object@assays$RNA@var.features, genes.filter.list))
  Temp.object@assays$RNA@var.features <- Temp.object@assays$RNA@var.features[!Temp.object@assays$RNA@var.features %in% genes.filter.list]
  nrow(Temp.object@meta.data)
  print(paste0("Variable genes length after filtering: ",length(Temp.object@assays$RNA@var.features)))
  
  print(paste0("Scaling now"))
  Temp.object <- ScaleData(object = Temp.object) 
  
  return(Temp.object)
  
  print("Done")
  print(Sys.time())
}



