library(fgsea)
library(data.table)
library(ggplot2)
library(tibble)
library(tidyr)
library(autoimage)
library(dplyr)
library(broom)

###################### GSEA for log2FC Translational Efficiency ####################
#note: ribosomal proteins were removed from the data frame
Ribo_df$log2FC_TE<-RPF27to31$RPF_log2FC-RPF27to31$Total_log2FC

Ribo_df<-Ribo_df[order(Ribo_df$log2FC_TE,decreasing=TRUE),]
Ribo_df<-subset(Ribo_df,Gene_ID!='CNOT1')
ranks <- deframe(Ribo_df[,c('Gene_ID','log2FC_TE')])
head(ranks, 20)


# Biological Process
pathwaystouse <- gmtPathways("\\\\john-doe/gw/Systems/Sarah/CNOTpaper/gene_ontology/c5.bp.v7.1.symbols.gmt")
xlab<-'Biological Process'

fgseaRes <- fgsea(pathways=pathwaystouse, stats=ranks, nperm=5000,minSize=5,maxSize=500)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

ff<-as.data.frame(fgseaResTidy[,c(1:5,7)])
write.table(x=subset(ff,padj<0.05),file = paste0(xlab,'_log2FC.csv'),col.names=TRUE,row.names=FALSE,sep=',')



# Molecular Function
pathwaystouse <- gmtPathways("\\\\john-doe/gw/Systems/Sarah/CNOTpaper/gene_ontology/c5.mf.v7.1.symbols.gmt")
xlab<-'Molecular Function'

fgseaRes <- fgsea(pathways=pathwaystouse, stats=ranks, nperm=5000,minSize=5,maxSize=500)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

ff<-as.data.frame(fgseaResTidy[,c(1:5,7)])
write.table(x=subset(ff,padj<0.05),file = paste0(xlab,'_log2FC.csv'),col.names=TRUE,row.names=FALSE,sep=',')



#Cellular Component
pathwaystouse <- gmtPathways("\\\\john-doe/gw/Systems/Sarah/CNOTpaper/gene_ontology/c5.cc.v7.1.symbols.gmt")
xlab<-'Cellular Component'

fgseaRes <- fgsea(pathways=pathwaystouse, stats=ranks, nperm=5000,minSize=5,maxSize=500)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

ff<-as.data.frame(fgseaResTidy[,c(1:5,7)])
write.table(x=subset(ff,padj<0.05),file = paste0(xlab,'_log2FC.csv'),col.names=TRUE,row.names=FALSE,sep=',')







