#' A Create_Table_Doublets Function
#'
#' This function create table for Seurat Normalzied data.
#' @param SCdata Seurat object to count number of cells.
#' @param BeforeFilter Total number of cells before filtering
#' @param Doublets Total number of Doublets detected
#' @param SampleName Sample Name
#' @param mincells Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
#' @param mingenes Include cells where at least this many features are detected.
#' @param mtpercent Include cells reporting at most this much mitochondrial transcript percentage.
#' @param rbpercent Include cells reporting at most this much ribosomal transcript percentage.
#' @param nRows Maximum number of rows to be prined.
#' @param coreFont core Font Size.
#' @param colheadFont colhead Font Size.
#' @param rowheadFont rowhead Font Size.
#' @param titlesize titlesize Font Size.
#' @keywords matrix.DIR, mincells, mingenes, mtpercent, rbpercent, Sample
#' @export
#' @examples
#' Create_Table_Doublets()



Create_Table_Doublets <- function(SCdata, BeforeFilter, Doublets, SampleName=NULL, mincells=3, mingenes=500, mtpercent=20, rbpercent=50, nRows=40, coreFont=2.5, colheadFont=2.5, rowheadFont=2.5, titlesize=25){
  
  cutoff.df <- data.frame(Value = t(data.frame(Cells = nrow(SCdata@meta.data), min.genes=mingenes, min.cells = mincells, MTpercent=mtpercent, RBpercent = rbpercent, Doublets=Doublets, BeforeFilter = BeforeFilter))); cutoff.df
  cutoff.df$Sample <- rownames(cutoff.df)
  print(head(cutoff.df))
  #colnames(cutoff.df) <- c("Sample", "Value")
  cutoff.df <- cutoff.df[,c("Sample", "Value")]
  TableDF <- cutoff.df
  FontsDF <- c(coreFont,colheadFont,rowheadFont)
  #titleDF <- paste0(Sample)
  titleDF <- paste0(SampleName)
  #func.PlotTable.General(TableDF, FontsDF, titleDF, 25)
  
  TableGI <- head(TableDF, nRows)
  rownames(TableGI) <- NULL
  tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
  mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = FontsDF[1])),
    colhead = list(fg_params=list(cex = FontsDF[2])),
    rowhead = list(fg_params=list(cex = FontsDF[3])))
  
  t1 <- tableGrob(TableGI, theme = mytheme)
  title <- textGrob(paste0(titleDF), gp = gpar(fontsize = titlesize))
  padding <- unit(5,"mm")
  
  table <- gtable_add_rows( t1, heights = grobHeight(title) + padding, pos = 0)
  table <- gtable_add_grob(table, title, 1, 1, 1, ncol(table))
  
  grid.newpage()
  grid.draw(table)
  
  
  print("Done")
  print(Sys.time())
}

