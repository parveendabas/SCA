#' A load_Colors Function
#'
#' This function load all the required packages for Seurat Analysis.
#' @param load_Colors load_Colors
#' @keywords load_Colors
#' @export 
#' @examples
#' load_Colors()



load_Colors <- function(){
 
	### Colors
## https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html#npg
libraries(c("easypackages","ggsci", "ggplot2", "gridExtra", "scales"))
mypal_npg = pal_npg("nrc", alpha = 0.7)(10); mypal_npg ## 10 Colors
mypal_nejm = pal_nejm("default", alpha = 0.7)(8); mypal_nejm ## 10 Colors
mypal_jco = pal_jco("default", alpha = 0.7)(10); mypal_jco ## 10 Colors
mypal_ucscgb = pal_ucscgb("default", alpha = 0.7)(26); mypal_ucscgb ## 10 Colors
mypal_startrek = pal_startrek("uniform", alpha = 0.7)(7); mypal_startrek ## 10 Colors
mypal_d3 = pal_d3("category20")(20); mypal_d3 ## 10 Colors

#show_col(c(mypal_ucscgb, mypal_npg, mypal_jco, mypal_nejm))
#show_col(c(mypal_jco[1:5], pal_d3("category20")(20), pal_lancet("lanonc")(9)))
#show_col(pal_d3("category20")(20))
#show_col(pal_d3("category10")(10))
#show_col(pal_d3("category20b")(20))
#show_col(pal_d3("category20c")(20))

ClusPallette <- c(mypal_ucscgb, mypal_npg, mypal_jco, mypal_nejm); ClusPalette
PairedPalette12 <- brewer.pal(n = 12, name = "Paired")
OtherPalette <- c(mypal_jco[1:5], pal_d3("category20")(20), pal_lancet("lanonc")(9))
 
  
}






