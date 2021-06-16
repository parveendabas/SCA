#' A Make_Regular_Plots Function
#'
#' This function allows you to perform differential gene analysis.
#' @param Temp.object Seurat objects to be used for the QC plots.
#' @param saveDIR Directory to save the plots.
#' @param SuffixName Suffix. to be added in the directory name as 
#' @param Species Species Name. Valid options are hsa or mmu
#' @param ColNamesToPlot Vector of identities to plot
#' @param ColPaletteToPlot list of color palettes to plot for the identities
#' @param HistogramPlot Plot Histograms or not
#' @param VlnPlot Plot Violins or not
#' @param QCumapPlots Plot QC UMAP or not
#' @param TablePlot Plot Tables or not
#' @param UMAPPlot Plot UMAP or not
#' @keywords Temp.object, saveDIR, SuffixName, Species, ColNamesToPlot, ColOrdersToPlot
#' @export
#' @examples
#' Make_Regular_Plots()



Make_Regular_Plots <- function(Temp.object, saveDIR, 
                               ColNamesToPlot=ColNamesToPlot, ColPaletteToPlot=ColPaletteToPlot,
                               SuffixName="QC_Plots", Species="hsa", 
                               HistogramPlot = TRUE, VlnPlot = TRUE, QCumapPlots = TRUE, TablePlot = TRUE, UMAPPlot = TRUE){
  
  
  #Temp.object=SCdata
  #saveDIR=pkWD
  #SuffixName=paste0("Regular_Plots_","Cube_Harmony")
  #ColNamesToPlot=c("seurat_clusters", "CT", "DietStrain", "Project", "Diet", "Strain")
  #ColPaletteToPlot=list(ClusPalette, Light.Pallette, Dark.Pallette, Dark.Pallette, Dark.Pallette, Dark.Pallette)
  
  
  setwd(saveDIR)
  plotWD <- paste(getwd(),paste0(SuffixName),sep="/"); print(plotWD)
  dir.create(file.path(getwd(),paste0(SuffixName)), showWarnings = FALSE)
  
  setwd(plotWD)
  HistogramPlot="YES"
  if(HistogramPlot=="YES"){
    
    pdf(file=paste0("Histogram_Plots_",SuffixName,".pdf"),height = 10,width = 12)
    for(i in 1:(length(ColNamesToPlot)-1)){
      #i=1
      ident1=ColNamesToPlot[i]
      ident1Palette=ColPaletteToPlot[[i]]
      print(paste0("ident1 for Histogram plotting is: ",ident1, " (",i,")"))
      
      for(j in (1:length(ColNamesToPlot))[!1:length(ColNamesToPlot) %in% i]){
        #j=2
        ident2=ColNamesToPlot[j]
        ident2Palette=ColPaletteToPlot[j]
        print(paste0("ident2 for plotting is: ",ident2, " (",j,")"))
        
    setwd(plotWD)
    print(paste0("Plotting Histogram Plot"))
    Hist <- melt(Temp.object@meta.data[,c(ident1, ident2)])
    head(Hist)
    #Hist[,ToUseCol] <- factor(Hist[,ToUseCol],levels = ToUseOrder)
    p1 <- ggplot(Hist, aes_string(ident2)) + geom_bar(aes_string(fill = ident1))  + scale_fill_manual(values=ident1Palette) +
      labs(x = ident2, fill=ident1) +
      theme(strip.background = element_rect(color="black", fill="grey", size=1.5, linetype="solid"),
            strip.text.x = element_text(size = 15, colour = "black"),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=15), legend.text=element_text(size=15), legend.title=element_text(size=15),
            legend.key.size = unit(0.4, "cm"), plot.title = element_text(size = 20)) +
      ggtitle(paste0(ident1," vs ",ident2))
    
    p2 <- ggplot(Hist, aes_string(x = ident2, fill = ident1)) +
      geom_bar(position="fill") + labs(x = ident2, y = "Percentage %", fill=ident1)  + scale_fill_manual(values=ident1Palette) +
      theme(strip.background = element_rect(color="black", fill="grey", size=1.5, linetype="solid"),
            strip.text.x = element_text(size = 15, colour = "black"),
            axis.title.x = element_text(color="black", size=15, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=15), legend.text=element_text(size=15), legend.title=element_text(size=15),
            legend.key.size = unit(0.4, "cm"), plot.title = element_text(size = 20)) +
      ggtitle(paste0(ident1," vs ",ident2))
    print(plot_grid(p1, p2, nrow = 2))
    
    #Bracket j
      }
      #Bracket i
    }
    dev.off()
  }
  
  
  
  VlnPlot="YES"
  if(VlnPlot=="YES"){
    
    pdf(file=paste0("QC_Violin_Plots_",SuffixName,".pdf"),height = 10,width = 12)
    for(i in 1:length(ColNamesToPlot)){
      #i=2
      ident1=ColNamesToPlot[i]
      ident1Palette=ColPaletteToPlot[[i]]
      print(paste0("ident1 for VlnPlot plotting is: ",ident1, " (",i,")"))
      
    Idents(Temp.object) <- ident1
    Idents(Temp.object) <- factor(Idents(Temp.object), levels = sort(unique(Temp.object@meta.data[,ident1])))
    print(VlnPlot(Temp.object, features = c("nFeature_RNA", "percent.mt", "percent.rb"), cols = ident1Palette, pt.size = 0.00, ncol = 1))
    }
    
    dev.off()
  }
  
  
  QCumapPlots="YES"
  if(QCumapPlots=="YES"){
    
    library(cmocean)
    pdf(file=paste0("QC_umap_Plots_",SuffixName,".pdf"),height = 9,width = 13)
    QCcols=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
    titelsUse = c("UMI (nCount_RNA)", "Genes (nFeature_RNA)", "mt%", "ribo%")
    names(titelsUse) <- QCcols
    QCcols.list=list()
    for(QCcolsUse in QCcols){
      #QCcolsUse="nCount_RNA"
      
      QCcols.list[[QCcolsUse]] <- FeaturePlot(Temp.object, features = QCcolsUse) +
        scale_colour_viridis_c(option = "viridis", direction = -1) + 
        ggtitle(paste0(titelsUse[QCcolsUse]))
      
      
      #scale_colour_viridis_c(option = "plasma", direction = -1)
      #scale_colour_viridis_c(option = "viridis", direction = -1)
      #scale_colour_viridis_c(option = "viridis")
      #scale_color_gradientn(colours = cmocean('amp')(50))
      #scale_color_gradientn(colours = cmocean('balance')(50))
    }
    
    if(Species=="hsa"){
      print(paste0("Plotting cell cycle gene TOP2A for ",Species))
      QCcols.list[["CC"]] <- FeaturePlot(Temp.object, features = "TOP2A") +
        scale_colour_viridis_c(option = "viridis", direction = -1) +
        ggtitle(paste0("Cell Cycle Gene"))
    } else {
      print(paste0("Plotting cell cycle gene Top2a for ",Species))
      QCcols.list[["CC"]] <- FeaturePlot(Temp.object, features = "Top2a") +
        scale_colour_viridis_c(option = "viridis", direction = -1) +
        ggtitle(paste0("Cell Cycle Gene"))
    }
    
    Idents(Temp.object) <- "seurat_clusters"
    print(paste0("Plotting Clusters"))
    QCcols.list[["seurat_clusters"]] <- DimPlot(Temp.object, reduction = "umap", cols = ClusPalette, label = T, label.size = 7, pt.size = 0.4)
    
    print(plot_grid(QCcols.list[[1]], QCcols.list[[2]], QCcols.list[["CC"]], QCcols.list[[3]], QCcols.list[[4]], QCcols.list[["seurat_clusters"]], ncol = 3))
    
    dev.off()
  }
  
  
  
  TablePlot="YES"
  if(TablePlot=="YES"){
    
    pdf(file=paste0("Tables_",SuffixName,".pdf"),height = 8,width = 16)
    for(i in 1:(length(ColNamesToPlot)-1)){
      #i=1
      ident1=ColNamesToPlot[i]
      ident1Palette=ColPaletteToPlot[[i]]
      print(paste0("ident1 for Table plotting is: ",ident1, " (",i,")"))
      
      PlotTable.withRowNames.OneValue(Temp.object, ident1, ident1, 15, 2, 2.5, 2.5)
      
      z=i+1
      for(j in z:length(ColNamesToPlot)){
        #j=2
        ident2=ColNamesToPlot[j]
        ident2Palette=ColPaletteToPlot[j]
        print(paste0("ident2 for plotting is: ",ident2, " (",j,")"))
        
        PlotTable.withRowNames.TwoValues(Temp.object, ident1, ident2, paste0(ident1," vs ",ident2), 15, 1.7, 2, 2)
    
      }
      
    }
  dev.off()  
  }
  
  
  
  
  UMAPPlot="YES"
  if(UMAPPlot=="YES"){
    
    pdf(file=paste0("UMAP_Plots_",SuffixName,".pdf"),height = 10,width = 16)
    QCcols.list=list()
    for(i in 1:length(ColNamesToPlot)){
      #i=2
      ident1=ColNamesToPlot[i]
      ident1Palette=ColPaletteToPlot[[i]]
      print(paste0("ident1 for UMAP plotting is: ",ident1, " (",i,")"))
      
      Idents(Temp.object) <- ident1
      Idents(Temp.object) <- factor(Idents(Temp.object), levels = sort(unique(Temp.object@meta.data[,ident1])))
      
      if(ident1 == "seurat_clusters" | ident1 == "CT"){
      QCcols.list[[ident1]] <- DimPlot(Temp.object, reduction = "umap", cols = ident1Palette, label = T, label.size = 6)
      } else {
        QCcols.list[[ident1]] <- DimPlot(Temp.object, reduction = "umap", cols = ident1Palette, label = F, label.size = 6)
      }
    
    }
    
    print(plot_grid(plotlist = QCcols.list, ncol = ceiling(length(ColNamesToPlot)/2)))
    dev.off()
  }
  
  
  print("Done")
  print(Sys.time())
}

