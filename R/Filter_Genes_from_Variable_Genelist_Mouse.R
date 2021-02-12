#' A Filter_Genes_from_Variable_Genelist_Mouse Function
#'
#' This function allows you to filter genes for mt, ribosomal and cell cycle list.
#' @param SCdata Seurat Object
#' @param Species Species Name. Valid options are hsa or mmu
#' @param mt mitochondrial genes (MOUSE) to be removed from Seurat object 
#' @param rb ribosomal genes (MOUSE) to be removed from Seurat object 
#' @param cc CellCycle genes (Converted from Human to Mouse) to be removed from Seurat object 
#' @keywords SCdata, mt, rb, cc
#' @export
#' @examples
#' Filter_Genes_from_Variable_Genelist_Mouse()

Filter_Genes_from_Variable_Genelist_Mouse <- function(SCdata, Species="hsa", mt = TRUE, rb = TRUE, cc = FALSE){
  
  
  if(Species=="hsa"){
    
    
    print(paste0("Fetching mitochondrial and ribosomal genes"))
    mt.genes <- rownames(SCdata@assays$RNA@counts)[rownames(SCdata@assays$RNA@counts) %like% "^MT-"]; mt.genes; length(mt.genes)
    rb.genes <- rownames(SCdata@assays$RNA@counts)[rownames(SCdata@assays$RNA@counts) %like% "^RP[SL]"]; rb.genes; length(rb.genes)
    
    if(cc == TRUE){
      print(paste0("Adding cell cycle genes to filter list"))  
      s.genes <- cc.genes$s.genes; head(s.genes); length(s.genes)
      g2m.genes <- cc.genes$g2m.genes; head(g2m.genes); length(g2m.genes)
      cellcycle.genes <- c(s.genes, g2m.genes); length(cellcycle.genes)
      genes.filter <- c(mt.genes, rb.genes, cellcycle.genes); head(genes.filter); print(length(genes.filter))
    } else {
      genes.filter <- c(mt.genes, rb.genes); head(genes.filter); print(length(genes.filter))
    }
    
  } else if (Species=="mmu"){
  
  print(paste0("Fetching mitochondrial and ribosomal genes"))
  mt.genes <- rownames(SCdata@assays$RNA@counts)[rownames(SCdata@assays$RNA@counts) %like% "^mt-"]; mt.genes; length(mt.genes)
  rb.genes <- rownames(SCdata@assays$RNA@counts)[rownames(SCdata@assays$RNA@counts) %like% "^Rp[sl]"]; rb.genes; length(rb.genes)
  
  if(cc == TRUE){
  print(paste0("Adding cell cycle genes to filter list"))  
  s.genes <- cc.genes$s.genes; head(s.genes); length(s.genes)
  s.genes <- stringr::str_to_title((s.genes)); head(s.genes); length(s.genes)
  g2m.genes <- cc.genes$g2m.genes; head(g2m.genes); length(g2m.genes)
  g2m.genes <- stringr::str_to_title((g2m.genes)); head(g2m.genes); length(g2m.genes)
  cellcycle.genes <- c(s.genes, g2m.genes); length(cellcycle.genes)
  genes.filter <- c(mt.genes, rb.genes, cellcycle.genes); head(genes.filter); print(length(genes.filter))
  } else {
    genes.filter <- c(mt.genes, rb.genes); head(genes.filter); print(length(genes.filter))
  }
  
  } else {
    print("Valid options for species are Human or Mouse only")
    break
  }
  return(genes.filter)
  
  print("Done")
  print(Sys.time())
}



