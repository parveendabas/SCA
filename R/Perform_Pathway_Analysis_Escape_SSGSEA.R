#' A Perform_Pathway_Analysis_Escape_SSGSEA Function
#'
#' This function allows you to perform differential gene analysis.
#' @param Temp.object A list of Seurat objects between which to find anchors for downstream integration.
#' @param saveDIR Path to save generated data.
#' @param SuffixName Suffix. to be added in the directory name as 
#' @param msigdbrCategoryList Msigdb Category List. H, C1:C8. https://www.gsea-msigdb.org/gsea/msigdb/
#' @param SubmsigdbrCategorySubset SubCategoryUse list based on MsigDB selected
#' @param minSizeGeneSet minimum Size of GeneSet for filtering
#' @param topnumber Max #pathways to be plotted per group
#' @param PDFheight PDFheight
#' @param PDFwidth PDFwidth
#' @keywords Temp.object, saveDIR, msigdbrCategoryList="H", SubmsigdbrCategorySubset = "", minSizeGeneSet=5, topnumber=10, PDFheight=10, PDFwidth=12
#' @export
#' @examples
#' Perform_Pathway_Analysis_Escape_SSGSEA()



Perform_Pathway_Analysis_Escape_SSGSEA <- function(Temp.object, saveDIR, SuffixName="Pathway_FGSEA", msigdbrCategoryList="H", SubmsigdbrCategorySubset = "", 
                                 minSizeGeneSet=5, topnumber=10, PDFheight=10, PDFwidth=12){
  
  #Temp.object=SeuratObj
  #msigdbrCategoryList="H"
  #SubmsigdbrCategorySubset = ""
  #minSizeGeneSet=5
  #topnumber=10
  #PDFheight=10
  #PDFwidth=12
  
  
  library(escape)
  library(Seurat)
  library(dittoSeq)
  library(ggplot2)
  library(EnrichmentBrowser)
  library(GSEABase)
  library(GSVA)
  
  colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                              "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                              "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))
  
  print(paste0("Top Genes to plot:",topnumber))
  
        
        for(CategoryUse in c(msigdbrCategoryList)){
          #CategoryUse="H"
          print(paste0("Processing category: ",CategoryUse))
          
          print(paste0("Check Point:  ----------->>>>>>>>>>>>>> MAIN Cat"))
          
          gene.sets <- getGeneSets(library = "H")  ## escape Packgae;
          gene.sets
          
          
          ES <- enrichIt(obj = SeuratObj, 
                         gene.sets = gene.sets, 
                         groups = 1000, cores = -1, 
                         min.size = NULL)
          
          
          
          
          
          
          for(SubCategoryUse in SubmsigdbrCategoryList){
            #SubCategoryUse="H"
            print(paste0("Processing category: ",SubCategoryUse))
            
            if(SubCategoryUse == "H"){
              print(paste0("No SubCategoryUse for H"))
              m_df_H <- msigdbr(species = "Mus musculus", category = CategoryUse) 
              
            } else {
              print(paste0("Processing Category:",CategoryUse," and SubCategoryUse:",SubCategoryUse))
              m_df_H <- msigdbr(species = "Mus musculus", category = CategoryUse, subcategory = SubCategoryUse)
            }
            
            Annot <- Temp.object@meta.data
            
            Annot.pathway2=as.data.frame(levels(as.factor(m_df_H$gs_name)))
            names(Annot.pathway2)="pathway"
            
            fgsea_sets<- m_df_H %>% split(x = .$gene_symbol, f = .$gs_name)
            print(paste0("Length of gene sets:",length(fgsea_sets)))
            
            Groups=levels(as.factor(Annot[,GroupName]))
            print(paste0("Total Groups:",length(Groups)))
            print(paste0(Groups))
            
            #run test
            Groups.genes <- wilcoxauc(Temp.object , GroupName)
            dplyr::count(Groups.genes, group)
            #name=paste0(CategoryUse,"_",SubCategoryUse); name
            
            name=paste0(CategoryUse,SubCategoryUse); name
            if(CategoryUse=="H"){name=paste0(CategoryUse); name}
            
            for (i in 1:length(Groups)){
              #i=1
              X=  Groups[i]
              
              cat("\n"); print(Sys.time()); 
              print(paste0("Processing Group number:",i," -> ",X))
              
              Genes<-Groups.genes %>%
                dplyr::filter(group == Groups[[i]]) %>%
                arrange(desc(auc)) %>% 
                dplyr::select(feature, auc)
              
              ranks<- deframe(Genes)
              print(head(ranks))
              
              
              fgseaRes<- fgsea(fgsea_sets, stats = ranks,minSize=minSizeGeneSet)
              
              fgseaResTidy <- fgseaRes %>%
                as_tibble() %>%
                arrange(desc(NES))
              fgseaResTidy %>% 
                dplyr::select(-leadingEdge, -ES, -log2err) %>% 
                arrange(padj) %>% 
                head()
              fgseaResTidy$Enrichment = ifelse(fgseaResTidy$NES > 0, "Up-regulated", "Down-regulated")
              Filtidy<-fgseaResTidy %>% filter(padj < 0.05) 
              filtRes = rbind(Filtidy %>% filter(Enrichment == "Up-regulated") %>% head(n= 20),
                              Filtidy %>% filter(Enrichment == "Down-regulated") %>% head(n= 20))
              Y<-fgseaResTidy
              hmc=as.data.frame(Y)
              hmc=apply(Y,2,as.character) 
              #write.csv(hmc,paste(name, X,".csv"))
              names(Y)[names(Y)=="NES"]=paste(X)
              Annot.pathway<-Y[,c("pathway",paste(X))]
              
              if(i == 1){
                print(paste0("Processing Entry First:",i))
                Annot.pathway.Merge=Annot.pathway
              } else {
                print(paste0("Processing Entry :",i))
                Annot.pathway.Merge<-merge(Annot.pathway, Annot.pathway.Merge, by.x="pathway", by.y="pathway")
              }
              print(head(Annot.pathway.Merge,1))
              
            }
            
            print(paste0("Check Point: ----------->>>>>>>>>>>>>> 1"))
            
            Annot.pathway2 <- Annot.pathway.Merge
            
            print(head(Annot.pathway2))
            
            ##make a heatmap
            rownames(Annot.pathway2)=Annot.pathway2$pathway
            #Annot.pathway2=Annot.pathway2[,-1]
            Annot.pathway2$pathway <- NULL
            positions=Groups
            
            Annot.pathway2=Annot.pathway2[,positions]
            Annot.pathway2[is.na(Annot.pathway2)]=0
            
            print(paste0("Check Point: ----------->>>>>>>>>>>>>>  2"))
            
            pathnames=c(row.names(Annot.pathway2 %>% top_n(n = topnumber, wt = eval(parse(text=names(Annot.pathway2)[1])))),row.names(Annot.pathway2 %>% top_n(n = -(topnumber), wt = eval(parse(text=names(Annot.pathway2)[1])))))
            for (i in 2: length(names(Annot.pathway2))){
              pathnames1=c(row.names(Annot.pathway2 %>% top_n(n = topnumber, wt = eval(parse(text=names(Annot.pathway2)[i])))),row.names(Annot.pathway2 %>% top_n(n = -(topnumber), wt = eval(parse(text=names(Annot.pathway2)[i])))))
              pathnames=c(pathnames,pathnames1)
            }
            pathnames =pathnames[!duplicated(pathnames)]
            Annot.pathway4=Annot.pathway2[pathnames,]
            Annot.pathway4[is.na(Annot.pathway4)]=0
            redblackblu=CustomPalette(low = "blue", high = "red", mid = "white", k = 1000)
            
            print(paste0("Check Point:  ----------->>>>>>>>>>>>>> 3"))
            
            setwd(plotWD)
            pdf(paste0(GroupName,"_Based_Pathway_ALLcells_Category_",name,"_cells",downsampleCells,"_top_",topnumber,".pdf"),height=PDFheight, width=PDFwidth)
            #pheatmap(Annot.pathway4,scale = "none",cluster_rows = T,cluster_cols = F,border_color = "black",treeheight_row = 0,legend_breaks = c(-5,-2,0,2,5),
            #         fontsize =8, cellwidth = 15,cellheight =10,color= viridis(8),angle_col = 90,
            #         main = paste("                                          ",name,objname," (NES)"))
            
            print(pheatmap(as.matrix(Annot.pathway4),scale = "none",cluster_rows = T,cluster_cols = F,border_color = "black",treeheight_row = 0,legend_breaks = c(-5,-2,0,2,5),
                           fontsize =8, cellwidth = 15,cellheight =10,color= c(brewer.pal(n = 5, name = "Blues"), plasma(8)), #angle_col = 90,
                           main = paste("                                          ",name," (NES)")))
            print(pheatmap(as.matrix(Annot.pathway4),scale = "none",cluster_rows = T,cluster_cols = T,border_color = "black",treeheight_row = 0,legend_breaks = c(-5,-2,0,2,5),
                           fontsize =8, cellwidth = 15,cellheight =10,color= c(brewer.pal(n = 5, name = "Blues"), plasma(8)), #angle_col = 90,
                           main = paste("                                          ",name," (NES)")))
            
            
            dev.off()
            
            print(paste0("Check Point: ----------->>>>>>>>>>>>>>  4"))
            saveRDS(Annot.pathway2, file=paste0(GroupName,"_Based_Pathway_ALLcells_Category_",name,"_cells",downsampleCells,".rds"))
          
            print(paste0("Check Point: ----------->>>>>>>>>>>>>>  5"))
            
            #SubCategoryUse
          }
          
          print(paste0("Check Point:  ----------->>>>>>>>>>>>>> 6"))
          
          #CategoryUse
        }
        
        print(paste0("Check Point:  ----------->>>>>>>>>>>>>> 7"))
        
        #Function
}





