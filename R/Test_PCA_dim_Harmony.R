#' A Test_PCA_dim_Harmony Function
#'
#' This function  allows you to compare different no. of PCAs for Harmony. Cluster Resolution is set to 0.1.
#' @param Temp.object Seurat objects to be used for the QC plots.
#' @param saveDIR Directory to save the plots.
#' @param IdentToBatchCorrect Identity to be used for Harmony batch correction (Max 1)
#' @param ThetaToBatchCorrect Theta values to be used while running harmony. Diversity clustering penalty parameter. Specify for each variable in group.by.vars. Default theta=2. theta=0 does not encourage any diversity. Larger values of theta result in more diverse clusters.
#' @param SuffixName Suffix. to be added in the directory name as 
#' @param ColNamesToPlot Vector of identities to plot (Max = 3)
#' @param ColPaletteToPlot list of color palettes to plot for the identities (Max = 3)
#' @param mingenes minimum gene number that will be mentioned in the output file name
#' @param PCAdimList Vector of the PCAs to be used for testing
#' @param DownSamplePCA #Cells to plot in UMAP for PCA testing
#' @keywords Temp.object, saveDIR, IdentToBatchCorrect, ThetaToBatchCorrect, SuffixName, ColNamesToPlot, ColPaletteToPlot mingenes PCAdimList
#' @export
#' @examples
#' Test_PCA_dim_Harmony()



Test_PCA_dim_Harmony <- function(Temp.object, saveDIR, IdentToBatchCorrect="orig.ident", ThetaToBatchCorrect=2, SuffixName="Testing_PCA_Dim",  
                               ColNamesToPlot=ColNamesToPlot, ColPaletteToPlot=ColPaletteToPlot,mingenes=500,
                               PCAdimList=c(20, 25, 30, 40, 50), ClusOrder = ClusOrderFrom1, DownSamplePCA=1000){
  
  
  #Temp.object=SCdata
  #saveDIR=plotWD.Subset
  #IdentToBatchCorrect=c("DietLibrary")
  #ThetaToBatchCorrect=THETAtest
  #SuffixName=paste0(OutputName,"_",SubsetName,"_theta_",THETAinpdf)
  #ColNamesToPlot=ColNamesToPlot
  #ColPaletteToPlot=ColPaletteToPlot
  #mingenes=500
  #PCAdimList = c(10, 15, 20, 25, 30, 40, 50)
  #ClusOrder = ClusOrderFrom1
  #DownSamplePCA=1000
  
  print(paste0("Covariates being used:",IdentToBatchCorrect))
  print(paste0("Theta for Covariates being used:",ThetaToBatchCorrect))
  
  ClusOrder <- ClusOrderFrom1
  
  print(paste0("Testing PCA dim:"))
  print(paste0(PCAdimList))
  LastPCAdim=rev(PCAdimList)[1]
  print(paste0("Last PCAdim to test is:",LastPCAdim))
  PCADiminpdf = paste(PCAdimList, collapse = '_'); PCADiminpdf
  
  setwd(saveDIR) 
  #for(resPlot in c(0.1, 0.2, 0.3, 0.4)){
  for(resPlot in c(0.1)){
    #resPlot=0.1
    print(paste0("Processing resolution:",resPlot))
    
    pdf(file=paste0("Testing_PCA_Dimensions_",SuffixName,"_minGenes_",mingenes,"_Res_",resPlot,"_PCAs_",PCADiminpdf,".pdf"),height = 9,width = 15)
    
    PCAumap.list=list()
    if(length(ColNamesToPlot) == 2){ print("Total ColNamesToPlot length: 2"); PCAumap.list.ident2=list()}
    if(length(ColNamesToPlot) == 3){ print("Total ColNamesToPlot length: 3"); PCAumap.list.ident2=list(); PCAumap.list.ident3=list()}
    if(length(ColNamesToPlot) == 4){ print("Total ColNamesToPlot length: 4"); PCAumap.list.ident2=list(); PCAumap.list.ident3=list(); PCAumap.list.ident4=list()}

    for(PCAdim in PCAdimList){
      #for(PCAdim in c(10)){
      
      #PCAdim=25
      print(paste0("Processing PCA:",PCAdim))
      
      Temp.object <- RunPCA(object = Temp.object, npcs = PCAdim, verbose = FALSE)
      Temp.object <- RunHarmony(Temp.object, group.by.vars = IdentToBatchCorrect, theta=ThetaToBatchCorrect)
      Temp.object <- RunUMAP(Temp.object, reduction = "harmony", dims = 1:PCAdim)
      Temp.object <- FindNeighbors(Temp.object, reduction = "harmony", dims = 1:PCAdim)
      Temp.object <- FindClusters(Temp.object, resolution=resPlot)
      
      ident1=ColNamesToPlot[1]
      ident1Palette=ColPaletteToPlot[[1]]
      
      Temp.object@meta.data$seurat_clusters <- MakeClustersFrom1(Temp.object@meta.data$seurat_clusters)
      Temp.object.DS <- subset(Temp.object,downsample=DownSamplePCA)
      PCAumap.list[[as.character(PCAdim)]] <- DimPlot(Temp.object.DS, reduction = "umap", cols = ident1Palette, label = T, label.size = 6)   +  ggtitle(paste0("PCA:",PCAdim, ", res:",resPlot))

	if(length(ColNamesToPlot) == 2){ print("Saving plots of ident 2"); 
	PCAumap.list.ident2[[as.character(PCAdim)]] <- DimPlot(Temp.object.DS, reduction = "umap", group.by=ColNamesToPlot[2], cols = ColPaletteToPlot[[2]], label = T, label.size = 6)
	} else if(length(ColNamesToPlot) == 3){ print("Saving plots of idents 2 & 3");  
	PCAumap.list.ident2[[as.character(PCAdim)]] <- DimPlot(Temp.object.DS, reduction = "umap", group.by=ColNamesToPlot[2], cols = ColPaletteToPlot[[2]], label = T, label.size = 6)
	PCAumap.list.ident3[[as.character(PCAdim)]] <- DimPlot(Temp.object.DS, reduction = "umap", group.by=ColNamesToPlot[3], cols = ColPaletteToPlot[[3]], label = T, label.size = 6)
	} else if(length(ColNamesToPlot) == 4){ print("Saving plots of idents 2, 3 & 4");  
        PCAumap.list.ident2[[as.character(PCAdim)]] <- DimPlot(Temp.object.DS, reduction = "umap", group.by=ColNamesToPlot[2], cols = ColPaletteToPlot[[2]], label = T, label.size = 6)
        PCAumap.list.ident3[[as.character(PCAdim)]] <- DimPlot(Temp.object.DS, reduction = "umap", group.by=ColNamesToPlot[3], cols = ColPaletteToPlot[[3]], label = T, label.size = 6)
        PCAumap.list.ident4[[as.character(PCAdim)]] <- DimPlot(Temp.object.DS, reduction = "umap", group.by=ColNamesToPlot[4], cols = ColPaletteToPlot[[4]], label = T, label.size = 6)
        }
    }
    
    
    PCAcols=ceiling(length(PCAdimList)/2); PCAcols
    print(paste0("Columns being plotted for PCAdims: ",PCAcols))
    
    print(ElbowPlot(Temp.object, ndims = LastPCAdim))
    print(plot_grid(plotlist = PCAumap.list, ncol = PCAcols))
    
	if(length(ColNamesToPlot) == 2){ print("Plotting ident 2"); 
	print(plot_grid(plotlist = PCAumap.list.ident2, ncol = PCAcols))
        } else if(length(ColNamesToPlot) == 3){ print("Plotting idents 2 & 3");
	print(plot_grid(plotlist = PCAumap.list.ident2, ncol = PCAcols))
	print(plot_grid(plotlist = PCAumap.list.ident3, ncol = PCAcols))
        } else if(length(ColNamesToPlot) == 4){ print("Plotting idents 2, 3 & 4");
	print(plot_grid(plotlist = PCAumap.list.ident2, ncol = PCAcols))
	print(plot_grid(plotlist = PCAumap.list.ident3, ncol = PCAcols))
	print(plot_grid(plotlist = PCAumap.list.ident4, ncol = PCAcols))
        }


      
        QCcols.list=list()
        for(i in 1:length(ColNamesToPlot)){
          #i=2
          ident1=ColNamesToPlot[i]
          ident1Palette=ColPaletteToPlot[[i]]
          print(paste0("ident1 for UMAP plotting is: ",ident1, " (",i,")"))
          
          
          Idents(Temp.object) <- ident1
          Idents(Temp.object) <- factor(Idents(Temp.object), levels = sort(unique(Temp.object@meta.data[,ident1])))
          
          if(ident1 == "seurat_clusters"){
            QCcols.list[[ident1]] <- DimPlot(Temp.object, reduction = "umap", cols = ident1Palette, label = T, label.size = 6)   +  ggtitle(paste0("PCA:",PCAdim, ", res:",resPlot))
          } else if (ident1 == "CT"){
            QCcols.list[[ident1]] <- DimPlot(Temp.object, reduction = "umap", cols = ident1Palette, label = T, label.size = 6)
          } else {
            QCcols.list[[ident1]] <- DimPlot(Temp.object, reduction = "umap", cols = ident1Palette, label = F, label.size = 6)
          }
          
        }
        
        ident1=ColNamesToPlot[1]
        ident1Palette=ColPaletteToPlot[[1]]
        Idents(Temp.object) <- ident1
        p1 <- VlnPlot(Temp.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.00, cols = ClusPallette)
        p2 <- VlnPlot(Temp.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.01, cols = ClusPallette)
        
        
        TopPanel <- plot_grid(plotlist = QCcols.list, ncol = length(ColNamesToPlot))
        BottomPanel <- plot_grid(p1)
        print(plot_grid(TopPanel, BottomPanel, nrow = 2))
        print(plot_grid(p2, nrow = 2))
        
        
      
      
    #print(plot_grid(plotlist = QCcols.list, ncol = ceiling(length(ColNamesToPlot)/2)))
    
    dev.off()
  }
  
  return(Temp.object)
  
  
  print("Done")
  print(Sys.time())
}

