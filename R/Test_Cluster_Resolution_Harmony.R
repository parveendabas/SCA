#' A Test_Cluster_Resolution_Harmony Function
#'
#' This function allows you to perform differential gene analysis.
#' @param Temp.object Seurat objects to be used for the QC plots.
#' @param saveDIR Directory to save the plots.
#' @param IdentToBatchCorrect Identity to be used for Harmony batch correction (Max 1)
#' @param SuffixName Suffix. to be added in the directory name as 
#' @param ColNamesToPlot Vector of identities to plot (Max = 3)
#' @param ColPaletteToPlot list of color palettes to plot for the identities (Max = 3)
#' @param mingenes minimum gene number that will be mentioned in the output file name
#' @param ClusResList Vector of the PCAs to be used for testing
#' @keywords Temp.object, saveDIR, IdentToBatchCorrect, SuffixName, ColNamesToPlot, ColPaletteToPlot mingenes PCAdimList
#' @export
#' @examples
#' Test_Cluster_Resolution_Harmony()



Test_Cluster_Resolution_Harmony <- function(Temp.object, saveDIR, IdentToBatchCorrect="orig.ident", SuffixName="Testing_Cluster_Resolution",  
                                 ColNamesToPlot=ColNamesToPlot, ColPaletteToPlot=ColPaletteToPlot,mingenes=500,
                                 ClusResList=c(0.1, 0.15, 0.2, 0.3), ClusOrder = ClusOrderFrom1){
  
  
  #Temp.object=SCdata
  #saveDIR=pkWD
  #SuffixName=paste0("Testing_Cluster_Resolution")
  #ColNamesToPlot="DietStrain"
  #ColPaletteToPlot=Dark.Pallette
  #IdentToBatchCorrect="orig.ident"
  #ClusResList=c(0.1, 0.15, 0.2, 0.3)
  #mingenes=500
  
  ClusOrder <- ClusOrderFrom1
  
  print(paste0("Testing Cluster Resolutions:"))
  print(paste0(ClusResList))
  
  setwd(saveDIR) 
  for(PCAdim in c(25)){
    
    #PCAdim=25
    
    Temp.object <- RunPCA(object = Temp.object, npcs = PCAdim, verbose = FALSE)
    Temp.object <- RunHarmony(Temp.object, group.by.vars = IdentToBatchCorrect)
    Temp.object <- RunUMAP(Temp.object, reduction = "harmony", dims = 1:PCAdim)
    
    p1=q1=list()
    
    pdf(file=paste0("Cluster_Resolution_Testing_Analysis_",SuffixName,"_minGenes_",mingenes,".pdf"),height = 14,width = 28)
    for(resUse in ClusResList){
      #resUse=0.1
      print(paste0("Processing res:",resUse))
      Temp.object <- FindNeighbors(Temp.object, reduction = "harmony", dims = 1:PCAdim)
      Temp.object <- FindClusters(Temp.object, resolution=resUse)
      Temp.object$seurat_clusters <- MakeClustersFrom1(Temp.object$seurat_clusters)
      
      #sort(unique(Temp.object$seurat_clusters))
      #paste0(ToUseColumn,"_Based")
      table(Temp.object@meta.data[,paste0("RNA_snn_res.",resUse)])
      Temp.object@meta.data[,paste0("RNA_snn_res.",resUse)] <- Temp.object$seurat_clusters
      table(Temp.object@meta.data[,paste0("RNA_snn_res.",resUse)])
      
      # Projecting singlet identities on TSNE visualization
      #DimPlot(Temp.object, group.by = "HTO_classification")
      p1[[as.character(resUse)]] <- DimPlot(Temp.object, group.by = "seurat_clusters", pt.size = 0.5, reduction = "umap", cols = ClusPalette, label = T, label.size = 8) + ggtitle(paste0("PCA:",PCAdim, ", Res:",resUse))
      Idents(Temp.object) <- ColNamesToPlot
      #Idents(Temp.object) <- factor(Idents(Temp.object), levels = GroupOrder)
      q1[[as.character(resUse)]] <- DimPlot(Temp.object, pt.size = 0.5, reduction = "umap", cols = ColPaletteToPlot, label = F, label.size = 7) + ggtitle(paste0("PCA:",PCAdim, ", Res:",resUse))
      #DimPlot(Temp.object, group.by = "LinDay", cols = NeuronPalette, pt.size = 1, reduction = "tsne")
      
    }
    
    TopPanel <- plot_grid(plotlist = p1, ncol = length(ClusResList))
    BottomPanel <- plot_grid(plotlist = q1, ncol = length(ClusResList))
    print(plot_grid(TopPanel, BottomPanel, nrow = 2))
    #print(plot_grid(p1[["0.1"]], p1[["0.2"]], p1[["0.3"]], p1[["0.4"]], p1[["0.5"]], q1[["0.1"]], q1[["0.2"]], q1[["0.3"]], q1[["0.4"]], q1[["0.5"]], nrow = 2))
    
    dev.off()
  }
  
  return(Temp.object)
  
  
  print("Done")
  print(Sys.time())
}

