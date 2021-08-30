library(fgsea)
library(data.table)
library(ggplot2)
library(tibble)
library(tidyr)
library(autoimage)
library(dplyr)

######################## mRNA half life GSEA analysis ############################

halflife_df<-halflife_df[order(halflife_df$log2FC_halflife,decreasing=TRUE),]
halflife_df<-subset(halflife_df,gene_name!='CNOT1')
ranks <- deframe(halflife_df[,c('gene_name','log2FC_halflife')])
head(ranks, 20)

#Biological Process
# Load the pathways into a named list
pathwaystouse <- gmtPathways("c5.bp.v7.1.symbols.gmt")
xlab<-'Biological Process'

#run analysis and plot significant terms
fgseaRes <- fgsea(pathways=pathwaystouse, stats=ranks, nperm=5000,minSize=50,maxSize=500)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

#plot data
ggplot(subset(fgseaResTidy,padj<0.05), aes(reorder(pathway, NES), NES))+scale_fill_gradient2(low='red3',mid='white',high='green4',breaks=c(-3,0,3),labels=c(-3,0,3))+
  geom_col(aes(fill=NES))+coord_flip()+theme(panel.grid.minor=element_blank(),panel.grid.minor.x=element_blank())+theme(axis.title=element_text(size=12),axis.text=element_text(size=12))+
labs(x=xlab, y="Normalized Enrichment Score")+ggtitle(xlab)+theme_minimal()+theme(axis.title.y=element_text(size=18),plot.title = element_text(size=16))






