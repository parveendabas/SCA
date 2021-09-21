#' A Final_ChosenPCA_and_Res_Harmony Function
#'
#' This function allows you to perform differential gene analysis.
#' @param Temp.object Seurat objects to be used for the QC plots.
#' @param saveDIR Directory to save the plots.
#' @param IdentToBatchCorrect Identity to be used for Harmony batch correction (Max 1)
#' @param ThetaToBatchCorrect Theta values to be used while running harmony. Diversity clustering penalty parameter. Specify for each variable in group.by.vars. Default theta=2. theta=0 does not encourage any diversity. Larger values of theta result in more diverse clusters.
#' @param SuffixName Suffix. to be added in the directory name as 
#' @param ChosenPCA ChosenPCA
#' @param ChosenRes ChosenRes
#' @param ColNamesToPlot Vector of identities to plot (Max = 3)
#' @param ColPaletteToPlot list of color palettes to plot for the identities (Max = 3)
#' @param mingenes minimum gene number that will be mentioned in the output file name
#' @keywords Temp.object, saveDIR, IdentToBatchCorrect, ThetaToBatchCorrect, SuffixName, ChosenPCA, ChosenRes, ColNamesToPlot, ColPaletteToPlot mingenes
#' @export
#' @examples
#' Final_ChosenPCA_and_Res_Harmony()



Final_ChosenPCA_and_Res_Harmony <- function(Temp.object, saveDIR, IdentToBatchCorrect="orig.ident", ThetaToBatchCorrect=2, SuffixName="Testing_PCA_Dim", 
                               ColNamesToPlot=ColNamesToPlot, ColPaletteToPlot=ColPaletteToPlot,mingenes=500,
                               ChosenPCA=25, ChosenRes=0.1,
                               ClusOrder = ClusOrderFrom1){
  
  print(paste0("Covariates being used:",IdentToBatchCorrect))
  print(paste0("Theta for Covariates being used:",ThetaToBatchCorrect))
  print(paste0("Chosen PCA is: ",ChosenPCA))
  print(paste0("Chosen Cluster Resolution is: ",ChosenRes))
  
  
  #Temp.object=SCdata
  #saveDIR=plotWD.Subset
  #IdentToBatchCorrect=c("DietLibrary")
  #ThetaToBatchCorrect=THETAtest
  #SuffixName=paste0(OutputName,"_",SubsetName,"_theta_",THETAinpdf)
  #ColNamesToPlot=ColNamesToPlot
  #ColPaletteToPlot=ColPaletteToPlot
  #mingenes=500
  #ChosenPCA=10
  #ChosenRes=0.2
  #ClusOrder = ClusOrderFrom1
  
  #Temp.object=SCdata
  #saveDIR=pkWD
  #SuffixName=paste0("Final_ChosenPCA_and_Resolution")
  #ColNamesToPlot=c("seurat_clusters", "DietLibrary", "DietStrain")
  #ColPaletteToPlot=list(ClusPallette, OtherPalette, PairedPalette12)
  #IdentToBatchCorrect="DietLibrary"
  #ChosenPCA=50
  #ChosenRes=0.25
  #mingenes=500
  
  setwd(saveDIR) 
  
  ClusOrder <- ClusOrderFrom1
  
    # We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
    Temp.object <- RunPCA(object = Temp.object, npcs = ChosenPCA, verbose = FALSE)
    Temp.object <- RunHarmony(Temp.object, group.by.vars = IdentToBatchCorrect, theta=ThetaToBatchCorrect)
    Temp.object <- RunUMAP(Temp.object, reduction = "harmony", dims = 1:ChosenPCA)
    
    Temp.object <- FindNeighbors(Temp.object, reduction = "harmony", dims = 1:ChosenPCA, graph.name = c("test", "cR1"))
    Temp.object <- FindClusters(Temp.object, resolution=ChosenRes, graph.name = "cR1")
    graphResName <- colnames(Temp.object@meta.data)[colnames(Temp.object@meta.data) %like% "cR1_res"]
    #graphResName <- graphResName[graphResName %like% ChosenRes]
    graphResName <- graphResName[str_detect(graphResName, paste0(ChosenRes,"$"))]; print(graphResName)
    
    
    Temp.object$seurat_clusters <- Temp.object@meta.data[,graphResName]
    table(Temp.object@meta.data[,graphResName])
    table(Temp.object@meta.data[,"seurat_clusters"])
    Temp.object$seurat_clusters <- MakeClustersFrom1(Temp.object$seurat_clusters)
    #sort(unique(Temp.object$seurat_clusters))
    
    TempColumn=paste0("cR1_res.",ChosenRes); TempColumn
    paste0(TempColumn,"_Based")
    table(Temp.object@meta.data[,paste0("cR1_res.",ChosenRes)])
    Temp.object@meta.data[,paste0("cR1_res.",ChosenRes)] <- Temp.object$seurat_clusters
    table(Temp.object@meta.data[,paste0("cR1_res.",ChosenRes)])
    
    QCcols.list=list()
    for(i in 1:length(ColNamesToPlot)){
      #i=2
      ident1=ColNamesToPlot[i]
      ident1Palette=ColPaletteToPlot[[i]]
      print(paste0("ident1 for UMAP plotting is: ",ident1, " (",i,")"))
      
      Idents(Temp.object) <- ident1
      Idents(Temp.object) <- factor(Idents(Temp.object), levels = sort(unique(Temp.object@meta.data[,ident1])))
      
      if(ident1 == "seurat_clusters"){
        QCcols.list[[ident1]] <- DimPlot(Temp.object, reduction = "umap", cols = ident1Palette, label = T, label.size = 6)   +  ggtitle(paste0("PCA:",ChosenPCA, ", res:",ChosenRes))
      } else if (ident1 == "CT"){
        QCcols.list[[ident1]] <- DimPlot(Temp.object, reduction = "umap", cols = ident1Palette, label = T, label.size = 6)
      } else {
        QCcols.list[[ident1]] <- DimPlot(Temp.object, reduction = "umap", cols = ident1Palette, label = F, label.size = 6)
      }
      
    }
    
    Idents(Temp.object) <- "seurat_clusters"
    p1 <- VlnPlot(Temp.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.00, cols = ClusPallette)
    p2 <- VlnPlot(Temp.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.01, cols = ClusPallette)
    
    
    TopPanel <- plot_grid(plotlist = QCcols.list, ncol = length(ColNamesToPlot))
    BottomPanel <- plot_grid(p1)
    
    pdf(file=paste0("FINAL_",SuffixName,"_minGenes_",mingenes,"_PCA_",ChosenPCA,"_Res_",ChosenRes,".pdf"),height = 9,width = 15)
    print(plot_grid(TopPanel, BottomPanel, nrow = 2))
    print(plot_grid(p2, nrow = 2))
    dev.off()
    
    
  return(Temp.object)
  
  
  print("Done")
  print(Sys.time())
}



