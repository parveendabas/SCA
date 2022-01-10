#' A Perform_DGE_ONEvsALL Function
#'
#' This function allows you to perform differential gene analysis.
#' @param Temp.object A list of Seurat objects between which to find anchors for downstream integration.
#' @param saveDIR Path to save generated data.
#' @param SuffixName Suffix. to be added in the directory name as 
#' @param MainCol Column to be used (Default is seurat_cluster)
#' @param MainColOrder Factor order to be used (Default is numerical sorted order)
#' @param MainColPallete color palette to be used
#' @param ColNamesToPlot Vector of identities to plot
#' @param ColPaletteToPlot list of color palettes to plot for the identities
#' @param Species Species Name. Valid options are hsa or mmu
#' @param FDR FDR cutoff
#' @param FCcutoff Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25 Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param topnumber Max #genes to be plotted per identity
#' @param min.pctagecells only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.25
#' @param test.use Denotes which test to use.
#' @param latent.vars Variables to test, used only when test.use is one of 'LR', 'negbinom', 'poisson', or 'MAST'
#' @param min.diff.pct only test genes that show a minimum difference in the fraction of detection between the two groups. Set to -Inf by default
#' @param only.pos.Genes Only return positive markers (TRUE by default)
#' @param downsampleHeatmap Max #cells to be plotted per identity
#' @param plots Save CCA plots
#' @param save Save integrated CCA RDS Seurat object
#' @param ComputeDGEs To run FindAllMarkers or not
#' @keywords Temp.object, saveDIR, SuffixName="ALLcells", MainCol, MainColOrder, MainColPallete, ColNamesToPlot, ColPaletteToPlot, Species, FDR, FCcutoff, topnumber, downsampleHeatmap, plots, save
#' @export
#' @examples
#' Perform_DGE_ONEvsALL()



Perform_DGE_ONEvsALL <- function(Temp.object, saveDIR, SuffixName="ALLcells", 
                                 MainCol="seurat_clusters", MainColOrder=ClusOrder, MainColPallete=ClusPallette,
                                 ColNamesToPlot=ColNamesToPlot, ColPaletteToPlot=ColPaletteToPlot, Species="hsa", 
                                 FDR = 0.1, FCcutoff = 1.5, topnumber = 5, min.pctagecells = 0.25, test.use = "wilcox", only.pos.Genes = TRUE, latent.vars = NULL, min.diff.pct = -Inf, downsampleHeatmap = 300, plots = TRUE, save = TRUE, ComputeDGEs="YES"){
  
  #Temp.object <- SCdata
  #saveDIR <- pkWD
  
  DgeNameInpdf <- paste0(MainCol,"_Based_",SuffixName,"_",test.use,"_min.pct_",min.pctagecells); DgeNameInpdf
  
  setwd(saveDIR)
  plotWD <- paste(getwd(),paste0("DGEs"),sep="/"); print(plotWD)
  dir.create(file.path(getwd(),paste0("DGEs")), showWarnings = FALSE)
  
  setwd(plotWD)
  plotWD1 <- paste(getwd(),paste0("DGEs_",DgeNameInpdf),sep="/"); print(plotWD1)
  dir.create(file.path(getwd(),paste0("DGEs_",DgeNameInpdf)), showWarnings = FALSE)
  
  Idents(Temp.object) <- MainCol
  DefaultAssay(Temp.object) <- "RNA"
  
  print(paste0("Differential Genes Positive only: ",only.pos.Genes))
  
  setwd(plotWD1)
  RUNheatmap="YES"
  if(RUNheatmap=="YES"){
  
    Idents(Temp.object) <- MainCol
    ##MainColOrder <- sort(unique(Temp.object@meta.data[,MainCol])); MainColOrder
    Idents(Temp.object) <- factor(Idents(Temp.object), levels= MainColOrder[MainColOrder %in% unique(Temp.object@meta.data[,MainCol])])


    if(ComputeDGEs=="YES"){
    markers <- FindAllMarkers(Temp.object, min.pct = min.pctagecells, logfc.threshold = log2(FCcutoff), assay = "RNA", test.use = test.use, latent.vars = latent.vars, min.diff.pct = min.diff.pct, only.pos = only.pos.Genes)
    markers <- markers[markers$p_val_adj < FDR,]; dim(markers)
    
    if(Species=="hsa"){
      print("Removing MT and Ribosomal genes of Human")
    markers <- Remove_Genes_Rp_mt_Rna_Human(markers)
    } else {
      print("Removing MT and Ribosomal genes of Mouse")
    markers <- Remove_Genes_Rp_mt_Rna_Mouse(markers)
    }
    } else{

    print(paste0("DGEs step already completed, reading marker info from the directory"))
    setwd(plotWD1)
    markers <- read.table(file = paste0("DEGs_Heatmap_",DgeNameInpdf,".txt"), header = T, row.names = 1)
    print(head(markers)); print(dim(markers))

    }
    
    head(markers, n = 15); print(dim(markers))
    
    if (nrow(markers) == 0) {
      
      print("Skipping Plotting")
      
      } else if(nrow(markers) == 1) {
        print("Plotting only DotPlot")

	pdf(file = paste0("DEGs_DotPlot_", DgeNameInpdf, 
                          ".pdf"), height = 6, width = 8)
        p1 <- DotPlot(Temp.object, features = (markers$gene), 
                      cols = c("gray80", "red")) + theme(axis.title.x = element_blank(), 
                                                         axis.title.y = element_blank(), axis.text = element_text(size = 13), 
                                                         legend.text = element_text(size = 13), legend.title = element_text(size = 15), 
                                                         legend.key.size = unit(0.4, "cm")) + RotatedAxis() + 
                scale_colour_viridis_c(option = "plasma") + ggtitle(paste0(MainCol))
        p2 <- VlnPlot(SeuratObj.CT.Strain, features = markers$gene, pt.size = 0)  + geom_boxplot(width=0.1, outlier.shape=NA)
        
        print(plot_grid(p1, p2, NULL, nrow = 2))
        
        dev.off()
        
      } else if (nrow(markers) >1 ) {

      top <- markers %>% group_by(cluster) %>% top_n(n = topnumber, wt = avg_log2FC); top <- top[!duplicated(top$gene),]; dim(top); 
      topFDR <- markers %>% dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>% dplyr::group_by(cluster) %>% dplyr::slice(1:topnumber); topFDR <- topFDR[!duplicated(topFDR$gene),]; dim(topFDR);
      write.table(markers, file = paste0("DEGs_Heatmap_",DgeNameInpdf,".txt"),quote=F,sep="\t")
      SCdata.temp.Heatmap <- subset(Temp.object, downsample=downsampleHeatmap)
      SCdata.temp.Heatmap <- ScaleData(object = SCdata.temp.Heatmap, verbose = FALSE, features = markers$gene, scale.max = 2)
      dtype="scale.data"
      nor.exp <- GetAssayData(object = SCdata.temp.Heatmap, slot = dtype); print(dim(nor.exp))
      UseGenes <- intersect(markers$gene, rownames(nor.exp)); length(UseGenes)
      nor.exp <- nor.exp[UseGenes,,drop=FALSE]; dim(nor.exp)
      
      meta.data.plot <- SCdata.temp.Heatmap@meta.data[,ColNamesToPlot]; head(meta.data.plot)
      meta.data.plot[,MainCol] <- factor(meta.data.plot[,MainCol], levels = MainColOrder[MainColOrder %in% unique(Temp.object@meta.data[,MainCol])])
      #meta.data.plot[,ColNamesToPlot[2]] <- factor(meta.data.plot[,ColNamesToPlot[2]], levels = sort(unique(meta.data.plot[,ColNamesToPlot[2]])))
      #meta.data.plot <- meta.data.plot[order(factor(meta.data.plot[,MainCol], levels = MainColOrder[MainColOrder %in% unique(Temp.object@meta.data[,MainCol])]), factor(meta.data.plot[,ColNamesToPlot[2]], levels = sort(unique(meta.data.plot[,ColNamesToPlot[2]])))),]; head(meta.data.plot)
      meta.data.plot <- meta.data.plot[order(factor(meta.data.plot[,MainCol], levels = MainColOrder[MainColOrder %in% unique(Temp.object@meta.data[,MainCol])])),]; head(meta.data.plot)
      nor.exp <- nor.exp[,rownames(meta.data.plot), drop=FALSE]
      Group.list <- c("Group1", "Group2", "Group3", "Group4", "Group5", "Group6", "Group7")
      colnames(meta.data.plot) <- c("Cluster", Group.list[1:(length(ColNamesToPlot)-1)])
      head(meta.data.plot)
      print(dim(nor.exp))
      
      Cluster.FULL <- MainColPallete[MainColOrder %in% unique(meta.data.plot[,"Cluster"])]; names(Cluster.FULL) <- MainColOrder[MainColOrder %in% unique(meta.data.plot[,"Cluster"])]; Cluster.FULL
      Group1.FULL <- ColPaletteToPlot[[2]][1:length(unique(meta.data.plot[,"Group1"]))]; names(Group1.FULL) <- sort(unique(meta.data.plot[,"Group1"])); Group1.FULL
      
      if(length(ColNamesToPlot)>2){ Group2.FULL <- ColPaletteToPlot[[3]][1:length(unique(meta.data.plot[,"Group2"]))]; names(Group2.FULL) <- sort(unique(meta.data.plot[,"Group2"])); Group2.FULL 
      ann_colors = list(
        Cluster = Cluster.FULL[names(Cluster.FULL) %in% as.character(unique(meta.data.plot$Cluster))],
        #CT = CT.FULL[names(CT.FULL) %in% as.character(unique(meta.data.plot$CT))],
        #Project = Project.FULL[names(Project.FULL) %in% as.character(unique(meta.data.plot$Project))],
        Group1 = Group1.FULL[names(Group1.FULL) %in% as.character(unique(meta.data.plot$Group1))],
        Group2 = Group2.FULL[names(Group2.FULL) %in% as.character(unique(meta.data.plot$Group2))]
      )}
      
      if(length(ColNamesToPlot)>3){ Group3.FULL <- ColPaletteToPlot[[4]][1:length(unique(meta.data.plot[,"Group3"]))]; names(Group3.FULL) <- sort(unique(meta.data.plot[,"Group3"])); Group3.FULL 
      ann_colors = list(
        Cluster = Cluster.FULL[names(Cluster.FULL) %in% as.character(unique(meta.data.plot$Cluster))],
        #CT = CT.FULL[names(CT.FULL) %in% as.character(unique(meta.data.plot$CT))],
        #Project = Project.FULL[names(Project.FULL) %in% as.character(unique(meta.data.plot$Project))],
        Group1 = Group1.FULL[names(Group1.FULL) %in% as.character(unique(meta.data.plot$Group1))],
        Group2 = Group2.FULL[names(Group2.FULL) %in% as.character(unique(meta.data.plot$Group2))],
        Group3 = Group3.FULL[names(Group3.FULL) %in% as.character(unique(meta.data.plot$Group3))]
        )}

      if(length(ColNamesToPlot)>4){ Group4.FULL <- ColPaletteToPlot[[5]][1:length(unique(meta.data.plot[,"Group4"]))]; names(Group4.FULL) <- sort(unique(meta.data.plot[,"Group4"])); Group4.FULL 
      ann_colors = list(
        Cluster = Cluster.FULL[names(Cluster.FULL) %in% as.character(unique(meta.data.plot$Cluster))],
        #CT = CT.FULL[names(CT.FULL) %in% as.character(unique(meta.data.plot$CT))],
        #Project = Project.FULL[names(Project.FULL) %in% as.character(unique(meta.data.plot$Project))],
        Group1 = Group1.FULL[names(Group1.FULL) %in% as.character(unique(meta.data.plot$Group1))],
        Group2 = Group2.FULL[names(Group2.FULL) %in% as.character(unique(meta.data.plot$Group2))],
        Group3 = Group3.FULL[names(Group3.FULL) %in% as.character(unique(meta.data.plot$Group3))],
        Group4 = Group4.FULL[names(Group4.FULL) %in% as.character(unique(meta.data.plot$Group4))]
        )}
      
      
      if(length(ColNamesToPlot)>5){ Group5.FULL <- ColPaletteToPlot[[6]][1:length(unique(meta.data.plot[,"Group5"]))]; names(Group5.FULL) <- sort(unique(meta.data.plot[,"Group5"])); Group5.FULL 
      ann_colors = list(
        Cluster = Cluster.FULL[names(Cluster.FULL) %in% as.character(unique(meta.data.plot$Cluster))],
        #CT = CT.FULL[names(CT.FULL) %in% as.character(unique(meta.data.plot$CT))],
        #Project = Project.FULL[names(Project.FULL) %in% as.character(unique(meta.data.plot$Project))],
        Group1 = Group1.FULL[names(Group1.FULL) %in% as.character(unique(meta.data.plot$Group1))],
        Group2 = Group2.FULL[names(Group2.FULL) %in% as.character(unique(meta.data.plot$Group2))],
        Group3 = Group3.FULL[names(Group3.FULL) %in% as.character(unique(meta.data.plot$Group3))],
        Group4 = Group4.FULL[names(Group4.FULL) %in% as.character(unique(meta.data.plot$Group4))],
        Group5 = Group5.FULL[names(Group5.FULL) %in% as.character(unique(meta.data.plot$Group5))]
      )}
      
      if(length(ColNamesToPlot)>6){ Group6.FULL <- ColPaletteToPlot[[7]][1:length(unique(meta.data.plot[,"Group6"]))]; names(Group6.FULL) <- sort(unique(meta.data.plot[,"Group6"])); Group6.FULL 
      ann_colors = list(
        Cluster = Cluster.FULL[names(Cluster.FULL) %in% as.character(unique(meta.data.plot$Cluster))],
        #CT = CT.FULL[names(CT.FULL) %in% as.character(unique(meta.data.plot$CT))],
        #Project = Project.FULL[names(Project.FULL) %in% as.character(unique(meta.data.plot$Project))],
        Group1 = Group1.FULL[names(Group1.FULL) %in% as.character(unique(meta.data.plot$Group1))],
        Group2 = Group2.FULL[names(Group2.FULL) %in% as.character(unique(meta.data.plot$Group2))],
        Group3 = Group3.FULL[names(Group3.FULL) %in% as.character(unique(meta.data.plot$Group3))],
        Group4 = Group4.FULL[names(Group4.FULL) %in% as.character(unique(meta.data.plot$Group4))],
        Group5 = Group5.FULL[names(Group5.FULL) %in% as.character(unique(meta.data.plot$Group5))],
        Group6 = Group6.FULL[names(Group6.FULL) %in% as.character(unique(meta.data.plot$Group6))]
      )}
      


      names(ann_colors) <- ColNamesToPlot      
      colnames(meta.data.plot) <- ColNamesToPlot
      print(ann_colors)
      print(head(meta.data.plot))
      
      colors <- c(seq(-2,2,by=0.01))
      my_palette <- c(colorRampPalette(colors = c("darkblue", "#a7c5f2", "#e6f0f5", "gray97", "darksalmon", "orangered3", "darkred")) (n = length(colors)))
      
      pdf(file=paste0("Counts_DEGs_Heatmap_",DgeNameInpdf,".pdf"),height = 8,width = 10)
      maxylim.temp <- max(table(markers$cluster)); maxylim <- maxylim.temp+ceiling(maxylim.temp*0.1)
      p1 <- ggplot(markers, aes(cluster)) + geom_bar(aes(fill=cluster)) + ylim(0,maxylim) + scale_fill_manual(values=ColPaletteToPlot[[1]]) + geom_text(stat='count', aes(label=..count..), vjust=-1)
      
      # https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
      library(UpSetR)
      library(ComplexHeatmap)
      VarNames <- sort(unique(markers$cluster)); VarNames
      UpsetSet <- lapply(X = VarNames, FUN = function(x){ (markers[markers$cluster == x,"gene"])}); UpsetSet
      names(UpsetSet) <- VarNames
      m <- make_comb_mat(UpsetSet)
      m.df <- as.data.frame(m)
      #UpSet(m)
      UpSetPlotOrder <- levels(Idents(Temp.object)); UpSetPlotOrder <- UpSetPlotOrder[UpSetPlotOrder %in% rownames(m.df)]
      print(paste0("Upset Order is:",UpSetPlotOrder))
      p2 <- grid.grabExpr(draw(UpSet(m, comb_col = "#0000FF", bg_col = "#F0F0FF", bg_pt_col = "#CCCCFF", top_annotation = upset_top_annotation(m, add_numbers = TRUE),
                                     right_annotation = upset_right_annotation(m, add_numbers = TRUE), set_order = UpSetPlotOrder)))
      
      print(plot_grid(p1, p2, ncol = 1, rel_heights = c(1, 1.5)))
      dev.off()

      pdf(file=paste0("DEGs_Heatmap_",DgeNameInpdf,".pdf"),height = 14,width = 20)
      print(DotPlot(Temp.object, features = (topFDR$gene), cols= c("gray80", "red"))  + 
              theme(axis.title.x=element_blank(), axis.title.y = element_blank(), axis.text=element_text(size=13), legend.text=element_text(size=13), legend.title=element_text(size=15),
                    legend.key.size = unit(0.4, "cm")) + RotatedAxis() + scale_colour_viridis_c(option = "plasma") + ggtitle(paste0(MainCol)))
      print(DotPlot(Temp.object, features = sort(topFDR$gene), cols= c("gray80", "red"))  + 
              theme(axis.title.x=element_blank(), axis.title.y = element_blank(), axis.text=element_text(size=13), legend.text=element_text(size=13), legend.title=element_text(size=15),
                    legend.key.size = unit(0.4, "cm")) + RotatedAxis() + scale_colour_viridis_c(option = "plasma") + ggtitle(paste0(MainCol)))
      print(pheatmap(nor.exp[topFDR$gene,],annotation_col=meta.data.plot,show_colnames=F,show_rownames=T,  color=my_palette, breaks=colors, 
                     annotation_colors = ann_colors, cluster_rows = FALSE, cluster_cols = FALSE, main = paste0("Heatmap All markers, Max cells plotted per main category:",downsampleHeatmap), fontsize = 13))
      print(pheatmap(nor.exp[top$gene,],annotation_col=meta.data.plot,show_colnames=F,show_rownames=T,  color=my_palette, breaks=colors, 
                     annotation_colors = ann_colors, cluster_rows = FALSE, cluster_cols = FALSE, main = paste0("Heatmap All markers, Max cells plotted per main category:",downsampleHeatmap), fontsize = 13))
      
      
      dev.off()
      
      rm(SCdata.temp.Heatmap)

        print(table(markers$cluster))
	return(markers)
    }
  }
  
  
  print("Done")
  print(Sys.time())
}

