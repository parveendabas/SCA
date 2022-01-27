#' A load_Packages_SeuratAnalysis Function
#'
#' This function load all the required packages for Seurat Analysis.
#' @param load_Packages_SeuratAnalysis load_Packages_SeuratAnalysis
#' @keywords load_Packages_SeuratAnalysis
#' @export 
#' @examples
#' load_Packages_SeuratAnalysis()



load_Packages_SeuratAnalysis <- function(){

#library(easypackages)
library(lubripack)

 
	lubripack(c("Seurat", "DoubletFinder", "Matrix", "stringr", "useful", "ggplot2", "gridExtra", "grid", "gtable", "data.table", "cowplot", "plyr", "reshape2", "tidyr", "pheatmap", "RColorBrewer", "tidyverse", "colorspace", "splines", "AnnotationDbi", "dplyr", "gridBase", "ggrepel", "calibrate", "ggplotify")) 

	lubripack(c("RColorBrewer", "ggpubr", "harmony", "SingleCellExperiment", "SummarizedExperiment", "easypackages", "SCA"))

	
  
}






