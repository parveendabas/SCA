#' A Perform_DGE_ONEvsALL Function
#'
#' This function allows you to perform differential gene analysis.
#' @param Temp.object A list of Seurat objects between which to find anchors for downstream integration.
#' @param SuffixName Suffix. to be added in the directory name as 
#' @param saveDIR Path to save generated data.
#' @param GroupColOrder Factor order for Sample group column name
#' @param Species Species Name. Valid options are hsa or mmu
#' @param GroupColPalette Sample group color palette
#' @param GroupCol Sample group column name
#' @param ToUseCol Column to be used (Default is seurat_cluster)
#' @param ToUseOrder Factor order to be used (Default is numerical sorted order)
#' @param ToUsePallete color palette to be used
#' @param FDR FDR cutoff
#' @param FCcutoff Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25 Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param topnumber Max #genes to be plotted per identity
#' @param downsampleHeatmap Max #cells to be plotted per identity
#' @param plots Save CCA plots
#' @param save Save integrated CCA RDS Seurat object
#' @keywords Temp.object, SuffixName="ALLcells", saveDIR, GroupColOrder , GroupColPalette, GroupCol, ToUseCol=, ToUseOrder, ToUsePallete, FDR, FCcutoff, topnumber, downsampleHeatmap, plots, save
#' @export
#' @examples
#' Perform_DGE_ONEvsALL()



Perform_DGE_ONEvsALL <- function(Temp.object, saveDIR, GroupColOrder, Species="hsa",
                                 SuffixName="ALLcells", GroupColPalette=Dark.Pallette, GroupCol="orig.ident",
                                 ToUseCol="seurat_clusters", ToUseOrder=ClusOrder, ToUsePallete=ClusPallette,
                                 FDR = 0.1, FCcutoff = 1.5, topnumber = 5, downsampleHeatmap = 300, plots = TRUE, save = TRUE){
  
  #Temp.object <- Seurat.object
  DgeNameInpdf <- paste0(ToUseCol,"_Based_",SuffixName); DgeNameInpdf
  
  setwd(saveDIR)
  plotWD <- paste(getwd(),paste0("DGEs"),sep="/"); print(plotWD)
  dir.create(file.path(getwd(),paste0("DGEs")), showWarnings = FALSE)
  
  setwd(plotWD)
  plotWD1 <- paste(getwd(),paste0("DGEs_",DgeNameInpdf),sep="/"); print(plotWD1)
  dir.create(file.path(getwd(),paste0("DGEs_",DgeNameInpdf)), showWarnings = FALSE)
  
  GroupColOrder = GroupColOrder
  
  Idents(Temp.object) <- ToUseCol
  DefaultAssay(Temp.object) <- "RNA"
  
  
  HistogramPlot="YES"
  if(HistogramPlot=="YES"){
    print(paste0("Plotting Histogram Plot"))
    Hist <- melt(Temp.object@meta.data[,c(GroupCol, ToUseCol)])
    head(Hist)
    #Hist[,ToUseCol] <- factor(Hist[,ToUseCol],levels = ToUseOrder)
    pdf(file=paste0("Count_Per_Cluster_",DgeNameInpdf,".pdf"),height = 8,width = 10)
    p1 <- ggplot(Hist, aes_string(ToUseCol)) + geom_bar(aes_string(fill = GroupCol))  + scale_fill_manual(values=GroupColPalette) +
      theme(strip.background = element_rect(color="black", fill="grey", size=1.5, linetype="solid"),
            strip.text.x = element_text(size = 20, colour = "black"),
            axis.title.x = element_text(color="black", size=20, face="bold"),
            axis.title.y = element_text(color="black", size=20, face="bold"),
            axis.text=element_text(size=20),
            legend.text=element_text(size=20), legend.title=element_text(size=20))
    
    
    p2 <- ggplot(Hist, aes_string(x = ToUseCol, fill = GroupCol)) +
      geom_bar(position="fill") + labs(x = "Cluster", y = "Percentage %", fill="Group")  + scale_fill_manual(values=GroupColPalette) +
      theme(strip.background = element_rect(color="black", fill="grey", size=1.5, linetype="solid"),
            strip.text.x = element_text(size = 20, colour = "black"),
            axis.title.x = element_text(color="black", size=20, face="bold"),
            axis.title.y = element_text(color="black", size=20, face="bold"),
            axis.text=element_text(size=20),
            legend.text=element_text(size=20), legend.title=element_text(size=20))
    print(plot_grid(p1, p2, rows = 2))
    dev.off()
    
    
  }
  
  
  TablePlot="YES"
  if(TablePlot=="YES"){
    print(paste0("Plotting Histogram Plot"))
    pdf(file=paste0("Table_",DgeNameInpdf,".pdf"),height = 8,width = 10)
    TableDF <- as.data.frame.matrix(table(Temp.object@meta.data[,GroupCol], Temp.object@meta.data[,ToUseCol]))
    TableDF <- TableDF[,colSums(TableDF) > 0 ]
    TableDF <- TableDF[rowSums(TableDF) > 0,]
    #FontsDF <- c(core,colhead,rowhead)
    #FontsDF <- c(1.2,1.2,1.2)
    FontsDF <- c(1.3,1.3,1.3)
    #FontsDF <- c(2.5,2.5,2.5)
    titleDF <- paste0("#Cells")
    func.PlotTable.withRowNames.General(TableDF, FontsDF, titleDF, 13)
    
    TableDF <- as.data.frame(table(Temp.object@meta.data[,GroupCol]))
    FontsDF <- c(1.3,1.3,1.3)
    #FontsDF <- c(2.5,2.5,2.5)
    titleDF <- paste0(GroupCol)
    func.PlotTable.withRowNames.General(TableDF, FontsDF, titleDF, 13)
    
    TableDF <- as.data.frame(table(Temp.object@meta.data[,ToUseCol]))
    FontsDF <- c(1.3,1.3,1.3)
    #FontsDF <- c(2.5,2.5,2.5)
    titleDF <- paste0(ToUseCol)
    func.PlotTable.withRowNames.General(TableDF, FontsDF, titleDF, 13)
    
    dev.off()
    
    
    
  }
  
  
  RUNVlnPlot="YES"
  if(RUNVlnPlot=="YES"){
    pdf(file=paste0("QC_Violin_Plots_",DgeNameInpdf,".pdf"),height = 10,width = 14)
    
    Idents(Temp.object) <- ToUseCol
    print(VlnPlot(Temp.object, features = c("nFeature_RNA", "percent.mt", "percent.rb"), cols = ToUsePallete, pt.size = 0.00, ncol = 1))
    
    Idents(Temp.object) <- GroupCol
    print(VlnPlot(Temp.object, features = c("nFeature_RNA", "percent.mt", "percent.rb"), cols = GroupColPalette, pt.size = 0.00, ncol = 1))
    
    dev.off()
  }
  
  
  RUNumapPlot="YES"
  if(RUNumapPlot=="YES"){
    pdf(file=paste0("UMAP_Plots_",DgeNameInpdf,".pdf"),height = 10,width = 14)
    
    Idents(Temp.object) <- GroupCol
    Idents(Temp.object) <- factor(Idents(Temp.object), levels = GroupColOrder)
    p1 <- DimPlot(Temp.object, reduction = "umap", cols = GroupColPalette, label = F, label.size = 6)
    
    Idents(Temp.object) <- ToUseCol
    ##Idents(Temp.object) <- factor(Idents(Temp.object), levels = ToUseOrder)
    p2 <- DimPlot(Temp.object, reduction = "umap", cols = ToUsePallete, label = T, label.size = 7, pt.size = 0.7)
    
    
    print(plot_grid(p1, p2, NULL, NULL, nrow = 2))
    print(plot_grid(p2))
    
    dev.off()
  }
  
  
  setwd(plotWD1)
  RUNheatmap="YES"
  if(RUNheatmap=="YES"){
    Idents(Temp.object) <- ToUseCol
    ##ToUseOrder <- sort(unique(Temp.object@meta.data[,ToUseCol])); ToUseOrder
    Idents(Temp.object) <- factor(Idents(Temp.object), levels= ToUseOrder[ToUseOrder %in% unique(Temp.object@meta.data[,ToUseCol])])
    markers <- FindAllMarkers(Temp.object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log2(FCcutoff), assay = "RNA")
    markers <- markers[markers$p_val_adj < FDR,]; dim(markers)
    
    if(Species=="hsa"){
      print("Removing MT and Ribosomal genes of Human")
    markers <- Remove_Genes_Rp_mt_Rna_Human(markers)
    } else {
      print("Removing MT and Ribosomal genes of Mouse")
    markers <- Remove_Genes_Rp_mt_Rna_Mouse(markers)
    }
    
    head(markers, n = 15); print(dim(markers))
    if(nrow(markers) > 0){
      top <- markers %>% group_by(cluster) %>% top_n(n = topnumber, wt = avg_logFC); top <- top[!duplicated(top$gene),]; dim(top); 
      topFDR <- markers %>% dplyr::arrange(p_val_adj, desc(avg_logFC)) %>% dplyr::group_by(cluster) %>% dplyr::slice(1:topnumber); topFDR <- topFDR[!duplicated(topFDR$gene),]; dim(topFDR);
      write.table(markers, file = paste0("DEGs_Heatmap_",DgeNameInpdf,".txt"),quote=F,sep="\t")
      SCdata.temp.Heatmap <- subset(Temp.object, downsample=downsampleHeatmap)
      SCdata.temp.Heatmap <- ScaleData(object = SCdata.temp.Heatmap, verbose = FALSE, features = markers$gene)
      dtype="scale.data"
      nor.exp <- GetAssayData(object = SCdata.temp.Heatmap, slot = dtype); print(dim(nor.exp))
      UseGenes <- intersect(markers$gene, rownames(nor.exp)); length(UseGenes)
      nor.exp <- nor.exp[UseGenes,,drop=FALSE]; dim(nor.exp)
      table(SCdata.temp.Heatmap@meta.data[,ToUseCol], SCdata.temp.Heatmap@meta.data[,GroupName])
      #check <- intersect(GroupOrder, unique(SCdata.temp.Heatmap@meta.data[,GroupName])); check
      meta.data.plot <- SCdata.temp.Heatmap@meta.data[,c(ToUseCol, GroupCol)]; head(meta.data.plot)
      meta.data.plot[,ToUseCol] <- factor(meta.data.plot[,ToUseCol], levels = ToUseOrder[ToUseOrder %in% unique(Temp.object@meta.data[,ToUseCol])])
      meta.data.plot[,GroupCol] <- factor(meta.data.plot[,GroupCol], levels = GroupColOrder)
      meta.data.plot <- meta.data.plot[order(factor(meta.data.plot[,ToUseCol], levels = ToUseOrder[ToUseOrder %in% unique(Temp.object@meta.data[,ToUseCol])]), factor(meta.data.plot[,GroupCol], levels = GroupColOrder)),]; head(meta.data.plot)
      nor.exp <- nor.exp[,rownames(meta.data.plot), drop=FALSE]
      colnames(meta.data.plot) <- c("Cluster", "Group")
      head(meta.data.plot)
      print(dim(nor.exp))
      
      #Cluster.FULL = c(`1` = "#0000ee", `2` = "#27408B", `3` = "#56B4E9", `4` = "#00ffff", `5` = "#ff0000", `6` = "#8b0000", `7` = "#CD5C5C", `8` = "#ff80bf", `9` = "#008000", `10` = "#00ff00", `11` = "#CAFF70", `12` = "#E69F00", `13` = "#ffb90f", `14` = "#DEB887", `15` = "#F0E442", `16` = "#ffff00", `17` = "#8B008B", `18` = "#bf3eff", `19` = "#7FFFD4", `20` = "#b5d4ff", `21` = "#FF1493", `22` = "#ff8247", `23` = "#836fff", `24` = "#787878", `25` = "#999999", `26` = "#D3D3D3", `27` = "#000000")
      #Group.FULL = c(`B6` = "#00008b", `S100A4-/-` = "#cd5b45")
      
      Cluster.FULL <- ToUsePallete[ToUseOrder %in% unique(meta.data.plot[,"Cluster"])]; names(Cluster.FULL) <- ToUseOrder[ToUseOrder %in% unique(meta.data.plot[,"Cluster"])]; Cluster.FULL
      Group.FULL <- GroupColPalette[1:length(GroupColOrder)]; names(Group.FULL) <- GroupColOrder; Group.FULL
      
      ann_colors = list(
        Cluster = Cluster.FULL[names(Cluster.FULL) %in% as.character(unique(meta.data.plot$Cluster))],
        #CT = CT.FULL[names(CT.FULL) %in% as.character(unique(meta.data.plot$CT))],
        #Project = Project.FULL[names(Project.FULL) %in% as.character(unique(meta.data.plot$Project))],
        Group = Group.FULL[names(Group.FULL) %in% as.character(unique(meta.data.plot$Group))]
        
      )
      colors <- c(seq(-2,4,by=0.01),max(nor.exp))
      my_palette <- c(colorRampPalette(colors = c("#a7c5f2", "#e6f0f5", "gray97", "darksalmon", "orangered3", "darkred")) (n = length(colors)-3), "firebrick4")
      
      pdf(file=paste0("DEGs_Heatmap_",DgeNameInpdf,".pdf"),height = 14,width = 20)
      print(pheatmap(nor.exp[topFDR$gene,],annotation_col=meta.data.plot,show_colnames=F,show_rownames=T,  color=my_palette, breaks=colors, 
                     annotation_colors = ann_colors, cluster_rows = FALSE, cluster_cols = FALSE, main = paste0("Heatmap All markers, Max cells plotted per group:",downsampleHeatmap), fontsize = 13))
      print(pheatmap(nor.exp[top$gene,],annotation_col=meta.data.plot,show_colnames=F,show_rownames=T,  color=my_palette, breaks=colors, 
                     annotation_colors = ann_colors, cluster_rows = FALSE, cluster_cols = FALSE, main = paste0("Heatmap All markers, Max cells plotted per group:",downsampleHeatmap), fontsize = 13))
      #print(pheatmap(nor.exp,annotation_col=meta.data.plot,show_colnames=F,show_rownames=T,  color=my_palette, breaks=colors, 
      #               cluster_rows = FALSE, cluster_cols = FALSE, main = paste0("Heatmap All markers"), fontsize = 13))
      dev.off()
      
      
      rm(SCdata.temp.Heatmap)
    }
  }
  
  
  setwd(plotWD1)
  SCdata.temp <- Temp.object
  Idents(SCdata.temp) <- ToUseCol
  #ToUseOrder <- sort(unique(SCdata.temp@meta.data[,ToUseCol])); ToUseOrder
  Idents(SCdata.temp) <- factor(Idents(SCdata.temp), levels= ToUseOrder)
  for(ID1 in sort(unique(Temp.object@meta.data[,ToUseCol]))){
    #for(ID1 in c(7,  8,  9,  10, 11, 12)){
    #ID1="1"
    print(paste0("Processing ",ID1))
    
    setwd(plotWD1)
    DEGcompdir <- paste(getwd(),paste0("DGEs_",ToUseCol,"_",ID1,"_VS_Rest"),sep="/"); print(DEGcompdir)
    dir.create(file.path(getwd(),paste0("DGEs_",ToUseCol,"_",ID1,"_VS_Rest")), showWarnings = FALSE)
    
    print(paste0("Finding out differential genes of ---- ",ID1, " ---- VS rest of the cells"))
    cond.markers <- FindMarkers(SCdata.temp, slot = "data", ident.1 = ID1, min.pct=0.25, logfc.threshold=log2(FCcutoff), verbose = TRUE, assay = "RNA")
    cond.markers$gene <- rownames(cond.markers)
    print(head(cond.markers, n = 15)); print(dim(cond.markers))
    cond.markers <- cond.markers[cond.markers$p_val_adj < FDR,]; dim(cond.markers)
    
    if(Species=="hsa"){
      print("Removing MT and Ribosomal genes of Human")
    cond.markers <- Remove_Genes_Rp_mt_Rna_Human(cond.markers)
    } else {
      print("Removing MT and Ribosomal genes of Mouse")
      cond.markers <- Remove_Genes_Rp_mt_Rna_Mouse(cond.markers)
    }
    
    sigGenes <- rownames(cond.markers)
    print(paste0("DONE: DGEs of ",ID1, " VS rest of the cells"))
    
    setwd(DEGcompdir)
    write.table(cond.markers, file = paste0("DGEs_",ID1,"_VS_Rest_FDR_",FDR,"_and_Fold_",FCcutoff,"_Genes",nrow(cond.markers),".txt"),quote=F,sep="\t")
    
    DGEcondName <- c("UP", "DOWN"); DGEcondName
    DGEcondGenes <- list()
    temp.DGE <- cond.markers[cond.markers$avg_logFC > 0,]; head(temp.DGE); dim(temp.DGE)
    DGEcondGenes[[1]] <- rownames(temp.DGE); head(DGEcondGenes[[1]]); length(DGEcondGenes[[1]])
    write.table(temp.DGE, file = paste0("DGEs_",ID1,"_VS_Rest_FDR_",FDR,"_and_Fold_",FCcutoff,"_UP_Genes",nrow(temp.DGE),".txt"),quote=F,sep="\t")
    
    temp.DGE <- cond.markers[cond.markers$avg_logFC < 0,]; head(temp.DGE); dim(temp.DGE)
    DGEcondGenes[[2]] <- rownames(temp.DGE); head(DGEcondGenes[[2]]); length(DGEcondGenes[[2]])
    write.table(temp.DGE, file = paste0("DGEs_",ID1,"_VS_Rest_FDR_",FDR,"_and_Fold_",FCcutoff,"_DOWN_Genes",nrow(temp.DGE),".txt"),quote=F,sep="\t")
    
    if(nrow(cond.markers[cond.markers$avg_logFC > 0,]) > 0){
      
      ## Discarding negative genes  
      cond.markers  <- cond.markers[cond.markers$avg_logFC > 0,]; dim(cond.markers)
      cond.markers <- head(cond.markers,30)
      
      ### Key plots
      RUNViolinNoPoints="YES"
      if(RUNViolinNoPoints=="YES"){
        setwd(DEGcompdir)
        print(paste0("ViolinPlot of positive to 30 genes by FC -- ",ID1, " VS rest of the cells"))
        pdf(file=paste0("DEGs_Violin_plots_",ID1,"_VS_Rest_FDR_",FDR,"_and_Fold_",FCcutoff,".pdf"),height = 24,width = 35)
        genesTo.Plot <- rownames(cond.markers); genesTo.Plot
        if(length(genesTo.Plot) < 7){
          set <- split(genesTo.Plot, ceiling(seq_along(genesTo.Plot)/6)); length(set)
        } else {
          set <- split(genesTo.Plot, ceiling(seq_along(genesTo.Plot)/15)); length(set)
        }
        for(i in 1:length(set)){
          #i=2
          print(paste0("Set no. ",i," of total ",length(set)," sets"))
          #plots <- VlnPlot(SCdata.temp, features = set[[i]], split.by = "Patient", pt.size = 0, combine = FALSE)
          plots <- VlnPlot(SCdata.temp, features = set[[i]], pt.size = 0, combine = FALSE, cols = ToUsePallete[ToUseOrder %in% unique(SCdata.temp@meta.data[,ToUseCol])])
          print(CombinePlots(plots = plots, ncol = 3))
          
        }
        dev.off()
      }
      
      RUNViolinPoints="YES"
      if(RUNViolinPoints=="YES"){
        setwd(DEGcompdir)
        print(paste0("ViolinPlot with points of positive to 30 genes by FC -- ",ID1, " VS rest of the cells"))
        pdf(file=paste0("DEGs_Violin_plots_",ID1,"_VS_Rest_FDR_",FDR,"_and_Fold_",FCcutoff,"_Points.pdf"),height = 24,width = 35)
        genesTo.Plot <- rownames(cond.markers); genesTo.Plot
        if(length(genesTo.Plot) < 7){
          set <- split(genesTo.Plot, ceiling(seq_along(genesTo.Plot)/6)); length(set)
        } else {
          set <- split(genesTo.Plot, ceiling(seq_along(genesTo.Plot)/15)); length(set)
        }
        for(i in 1:length(set)){
          #i=2
          print(paste0("Set no. ",i," of total ",length(set)," sets"))
          #plots <- VlnPlot(SCdata.temp, features = set[[i]], split.by = "Patient", pt.size = 0, combine = FALSE)
          plots <- VlnPlot(SCdata.temp, features = set[[i]], pt.size = 0.1, combine = FALSE, cols =  ToUsePallete[ToUseOrder %in% unique(SCdata.temp@meta.data[,ToUseCol])])
          print(CombinePlots(plots = plots, ncol = 3))
          
        }
        dev.off()
      }
      
      
      RUNDotPlot="YES"
      if(RUNDotPlot=="YES"){
        setwd(DEGcompdir)
        print(paste0("Dotplot of positive to 30 genes by FC -- ",ID1, " VS rest of the cells"))
        pdf(file=paste0("DEGs_DotPlot_plots_",ID1,"_VS_Rest_FDR_",FDR,"_and_Fold_",FCcutoff,".pdf"),height = 14,width =20)
        genesTo.Plot <- rownames(cond.markers); genesTo.Plot
        set <- split(genesTo.Plot, ceiling(seq_along(genesTo.Plot)/30)); length(set)
        for(i in 1:length(set)){
          #i=1
          print(paste0("Set no. ",i," of total ",length(set)," sets"))
          print(DotPlot(SCdata.temp, features = rev(set[[i]]), dot.scale = 10, cols = c("dodgerblue", "red")) + RotatedAxis())
          #print(DotPlot(SCdata.temp, features = rev(set[[i]]), dot.scale = 10, cols = c("blue", "red", "green", "pink"), split.by = GroupName) + RotatedAxis())
          
        }
        dev.off()
      }  
      
      
      
      
      RUNthis="YES"
      if(RUNthis=="YES"){
        setwd(DEGcompdir)
        print(paste0("Featureplot of positive to 30 genes by FC -- ",ID1, " VS rest of the cells"))
        SCdata.temp.Feature <- subset(SCdata.temp, downsample=downsampleHeatmap)
        SCdata.temp.Feature
        pdf(file=paste0("DEGs_Feature_plots_",ID1,"_VS_Rest_FDR_",FDR,"_and_Fold_",FCcutoff,".pdf"),height = 16,width =18)
        genesTo.Plot <- rownames(cond.markers); genesTo.Plot
        if(length(genesTo.Plot) < 7){
          set <- split(genesTo.Plot, ceiling(seq_along(genesTo.Plot)/6)); length(set)
        } else {
          set <- split(genesTo.Plot, ceiling(seq_along(genesTo.Plot)/15)); length(set)
        }
        for(i in 1:length(set)){
          #i=1
          print(paste0("Set no. ",i," of total ",length(set)," sets"))
          print(FeaturePlot(SCdata.temp.Feature, features = set[[i]], ncol=3, #max.cutoff = 4, 
                            cols = c("grey", "blue")))
          
        }
        
        dev.off()
      }
      
      
      RUNthis="YES"
      if(RUNthis=="YES"){
        
        if(nrow(cond.markers) > 1){
          print(paste0("Heatmap of positive to 30 genes by FC -- ",ID1, " VS rest of the cells"))
          sigGenes <- rownames(cond.markers); sigGenes
          ### Plotting Heatmap
          ### ReCopying Meta Data into SampleInfo.temp after Corrections
          table(SCdata.temp@meta.data[,ToUseCol])
          SCdata.temp.Heatmap <- subset(Temp.object, downsample=downsampleHeatmap)
          SCdata.temp.Heatmap <- ScaleData(object = SCdata.temp.Heatmap, verbose = FALSE, features = cond.markers$gene)
          dtype="scale.data"
          nor.exp <- GetAssayData(object = SCdata.temp.Heatmap, slot = dtype); print(dim(nor.exp))
          UseGenes <- intersect(rownames(cond.markers), rownames(nor.exp)); length(UseGenes)
          nor.exp <- nor.exp[UseGenes,,drop=FALSE]; dim(nor.exp)
          
          table(SCdata.temp.Heatmap@meta.data[,ToUseCol], SCdata.temp.Heatmap@meta.data[,GroupName])
          meta.data.plot <- SCdata.temp.Heatmap@meta.data[,c(ToUseCol, GroupCol)]; head(meta.data.plot)
          meta.data.plot[,ToUseCol] <- factor(meta.data.plot[,ToUseCol], levels = ToUseOrder[ToUseOrder %in% unique(Temp.object@meta.data[,ToUseCol])])
          meta.data.plot[,GroupCol] <- factor(meta.data.plot[,GroupCol], levels = GroupColOrder)
          meta.data.plot <- meta.data.plot[order(factor(meta.data.plot[,ToUseCol], levels = ToUseOrder[ToUseOrder %in% unique(Temp.object@meta.data[,ToUseCol])]), factor(meta.data.plot[,GroupCol], levels = GroupColOrder)),]; head(meta.data.plot)
          nor.exp <- nor.exp[,rownames(meta.data.plot), drop=FALSE]
          colnames(meta.data.plot) <- c("Cluster", "Group")
          head(meta.data.plot)
          print(dim(nor.exp))
          
          setwd(DEGcompdir)
          ## Plot only CTfinal
          colors <- c(seq(-2,4,by=0.01),max(nor.exp))
          my_palette <- c(colorRampPalette(colors = c("#a7c5f2", "#e6f0f5", "gray97", "darksalmon", "orangered3", "darkred")) (n = length(colors)-3), "firebrick4")
          pdf(file=paste0("DEGs_Heatmap_plots_",ID1,"_VS_Rest_FDR_",FDR,"_and_Fold_",FCcutoff,"_Points.pdf"),height = 14,width = 16)
          print(pheatmap(nor.exp,annotation_col=meta.data.plot,show_colnames=F,show_rownames=T,  color=my_palette, breaks=colors, 
                         annotation_colors = ann_colors, cluster_rows = FALSE, cluster_cols = FALSE, main = paste0("Heatmap, ",ID1,"_VS_Rest, Max cells plotted per group:",downsampleHeatmap), fontsize = 13))
          print(pheatmap(nor.exp,annotation_col=meta.data.plot,show_colnames=F,show_rownames=T,  color=my_palette, breaks=colors, 
                         annotation_colors = ann_colors, cluster_rows = TRUE, cluster_cols = FALSE, main = paste0("Heatmap, ",ID1,"_VS_Rest, Max cells plotted per group:",downsampleHeatmap), fontsize = 13))
          print(pheatmap(nor.exp,annotation_col=meta.data.plot,show_colnames=F,show_rownames=T,  color=my_palette, breaks=colors, 
                         annotation_colors = ann_colors, cluster_rows = TRUE, cluster_cols = TRUE, main = paste0("Heatmap, ",ID1,"_VS_Rest, Max cells plotted per group:",downsampleHeatmap), fontsize = 13))
          
          dev.off()
        }
      }
    }
  }
  
  print("Done")
  print(Sys.time())
}

