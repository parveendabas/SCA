#' A Perform_Pathway_Analysis_FGSEA Function
#'
#' This function allows you to perform differential gene analysis. 
#' Category: H, SubCategory: None
#' Category: C2, SubCategory: CP=29, CP:KEGG=186, CP:PID=196, CP:BIOCARTA=292, CP:WIKIPATHWAYS=615, CP:REACTOME=1604, CGP=3368
#' Category: C5, SubCategory: GO:BP=7480, GO:CC=996, GO:MF=1707
#' @param Temp.object A list of Seurat objects between which to find anchors for downstream integration.
#' @param saveDIR Path to save generated data.
#' @param SuffixName Suffix. to be added in the directory name as 
#' @param msigdbrCategoryList Msigdb Category List. H, C1:C8. https://www.gsea-msigdb.org/gsea/msigdb/
#' @param SubmsigdbrCategorySubset SubCategoryUse list based on MsigDB selected
#' @param minSizeGeneSet minimum Size of GeneSet for filtering
#' @param DownsamplePerIdent Number of cells to use per Idents
#' @param GroupName Column to use for Identity
#' @param FDRcutoff FDRcutoff
#' @param topnumber Max #pathways to be plotted per group
#' @param PDFheight PDFheight
#' @param PDFwidth PDFwidth
#' @keywords Temp.object, saveDIR, msigdbrCategoryList="H", SubmsigdbrCategorySubset = "", minSizeGeneSet=5, topnumber=10, PDFheight=10, PDFwidth=12
#' @export
#' @examples
#' Perform_Pathway_Analysis_FGSEA()



Perform_Pathway_Analysis_FGSEA <- function(Temp.object, saveDIR, SuffixName="Pathway_FGSEA", msigdbrCategoryList="H", SubmsigdbrCategorySubset = "", 
				DownsamplePerIdent="ALL", GroupName="Cluster", FDRcutoff=0.1,
                                 minSizeGeneSet=5, topnumber=10, PDFheight=10, PDFwidth=12){
  
  #Temp.object=SCdata.temp
  #saveDIR=plotWD.subset
  #SuffixName="Beta_Both_LF_and_HF"
  #msigdbrCategoryList="H"
  #SubmsigdbrCategorySubset=""
  
  #DownsamplePerIdent="ALL"
  #GroupName=IdentGroup
  #FDRcutoff=0.1
  
  #minSizeGeneSet=3
  #topnumber=10
  #PDFheight=10
  #PDFwidth=12
  
  
  library(presto)
  library(msigdbr)
  library(dplyr)
  library(fgsea)
  library(tibble)
  library(pheatmap)
  library(viridis)
  
  setwd(saveDIR)
  PathwaydirMain <- paste(getwd(),paste0("Pathway_Analysis"),sep="/"); print(PathwaydirMain)
  dir.create(file.path(getwd(),paste0("Pathway_Analysis")), showWarnings = FALSE)
  setwd(PathwaydirMain)
  
  Temp.object.1 <- Temp.object 
  if(DownsamplePerIdent == "ALL"){
    print("Using All Cells for the Analysis")
    } else {
    DownsamplePerIdent=as.numeric(DownsamplePerIdent)
    print(paste0("Downsampling the data, per ident cells = ",DownsamplePerIdent))
    Temp.object=subset(Temp.object, downsample=DownsamplePerIdent) 
    }

    print(Temp.object.1)
    print(Temp.object)

 
  print(paste0("Top Genes to plot:",topnumber))
        
        for(CategoryUse in c(msigdbrCategoryList)){
          #CategoryUse="H"
          print(paste0("Processing category: ",CategoryUse))
          
          print(paste0("Check Point:  ----------->>>>>>>>>>>>>> MAIN Cat"))
          
          if(length(SubmsigdbrCategorySubset) == 1 & unique(SubmsigdbrCategorySubset=="")){
            print(paste0("SubCategoryUse is not provided, using default options for categories"))
            
            print(paste0("Check Point:  ----------->>>>>>>>>>>>>> MAIN Sub"))
            
            if(CategoryUse == "H"){
              
              print(paste0("Assigning Empty for Category H"))
              SubmsigdbrCategoryList="H"
              #50 = 50
              
            } else if(CategoryUse == "C2"){
              
              print(paste0("Assigning Empty for Category H"))
              SubmsigdbrCategoryList=c("CP", "CP:KEGG", "CP:PID", "CP:BIOCARTA", "CP:WIKIPATHWAYS", "CP:REACTOME")
              ## GeneSet Length
              #"CP" = 29
              #"CP:KEGG" = 186
              #"CP:PID" = 196
              #"CP:BIOCARTA" = 292
              #"CP:WIKIPATHWAYS" = 615
              #"CP:REACTOME" = 1604
              #"CGP" = 3368
              
            } else if(CategoryUse == "C5"){
              
              print(paste0("Assigning Empty for Category H"))
              SubmsigdbrCategoryList=c("GO:CC", "GO:MF", "GO:BP")
              ### Gene Sets
              #GO:BP=7480
              #GO:CC=996
              #GO:MF=1707
            }
          } else {
            print(paste0("User Define sub category"))
            SubmsigdbrCategoryList=SubmsigdbrCategorySubset
            print(paste0(SubmsigdbrCategoryList))
          }
          
          print(paste0("Using following SubmsigdbrCategoryList"))
          print(paste0(SubmsigdbrCategoryList))
          
          
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
                arrange(dplyr::desc(auc)) %>% 
                dplyr::select(feature, auc)
              
              ranks<- deframe(Genes)
              print(head(ranks))
              
              
              fgseaRes<- fgsea(fgsea_sets, stats = ranks,minSize=minSizeGeneSet)
              
              fgseaResTidy <- fgseaRes %>%
                as_tibble() %>%
                arrange(dplyr::desc(NES))
		

              fgseaResTidy %>% 
                dplyr::select(-leadingEdge, -ES, -log2err) %>% 
                arrange(padj) %>% 
                head()


              fgseaResTidy$Enrichment = ifelse(fgseaResTidy$NES > 0, "Up-regulated", "Down-regulated")
              Filtidy<-fgseaResTidy %>% filter(padj < FDRcutoff) 
	      fgseaResTidy$cluster <- X
	       
	      print(paste0("Pathways passed at FDR cutoff of ",FDRcutoff," are: ",dim(Filtidy)))

              #filtRes = rbind(Filtidy %>% filter(Enrichment == "Up-regulated") %>% head(n=250),
              #                Filtidy %>% filter(Enrichment == "Down-regulated") %>% head(n=250))

              Y<-fgseaResTidy
              hmc=as.data.frame(Y)
              hmc=apply(Y,2,as.character) 
              #write.csv(hmc,paste(name, X,".csv"))
              names(Y)[names(Y)=="NES"]=paste(X)
              Annot.pathway<-Y[,c("pathway",paste(X))]
              print(head(Annot.pathway,1))
              
              if(i == 1){
                print(paste0("Processing Entry First:",i))
                Annot.pathway.Merge=Annot.pathway
                Combined.df <- fgseaResTidy; print(tail(Combined.df)); print(dim(Combined.df))
              } else {
                print(paste0("Processing Entry :",i))
                Annot.pathway.Merge<-merge(Annot.pathway, Annot.pathway.Merge, by.x="pathway", by.y="pathway")
                Combined.df <- rbind(Combined.df,fgseaResTidy); print(tail(Combined.df)); print(dim(Combined.df))
              }
              print(head(Annot.pathway.Merge,1))
              
            }
            
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
            
            setwd(PathwaydirMain)
            pdf(paste0(SuffixName,"_",GroupName,"_Based_Pathway_ALLcells_Category_",name,"_cells",DownsamplePerIdent,"_top_",topnumber,".pdf"),height=PDFheight, width=PDFwidth)
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
            
            saveRDS(Annot.pathway2, file=paste0(SuffixName,"_",GroupName,"_Based_Pathway_ALLcells_Category_",name,"_cells",DownsamplePerIdent,"_NES.rds"))
            saveRDS(Combined.df, file=paste0(SuffixName,"_",GroupName,"_Based_Pathway_ALLcells_Category_",name,"_cells",DownsamplePerIdent,".rds"))
          
            
            #SubCategoryUse
          }
          
          
          #CategoryUse
        }
        
        
        #Function
}





