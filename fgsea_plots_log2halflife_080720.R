library(fgsea)
library(data.table)
library(ggplot2)
library(tibble)
library(tidyr)
library(autoimage)
library(dplyr)


#get data
halflives<-read.table('Z:/Systems/Sarah/CNOTdataAnalysis/mRNAhalflife/halflife_oneTperG_counts/nls_halflife_010620.txt',stringsAsFactors = FALSE,header=TRUE)
names(halflives)
nrow(halflives)
halflives<-subset(halflives,startsWith(halflives$gene_name,'RPL')==FALSE)
nrow(halflives)
halflives<-subset(halflives,startsWith(halflives$gene_name,'RPS')==FALSE)
nrow(halflives)
halflives<-subset(halflives,startsWith(halflives$gene_name,'MRP')==FALSE)
nrow(halflives)
halflives<-subset(halflives,startsWith(halflives$gene_name,'RP11-')==FALSE)
nrow(halflives)
#

test<-merge(halflives,humanDF[,c(5,10:12)],by.x='gene_name',by.y='Gene_ID')
nrow(test) 
test<-subset(test,nchar(CDS_sequence)>200)
nrow(test)
test<-subset(test,nchar(CDS_sequence)<10000)
nrow(test)


#####################

#could try rank by p-value (giving +ve and -ve to p-val of direction of change??)
setwd('\\\\john-doe/gw/Systems/Sarah/CNOTpaper/gene_ontology')


################### mRNA half life ############################

test_gsea<-test
test_gsea<-test_gsea[order(test_gsea$log2FC_halflife,decreasing=TRUE),]
test_gsea<-subset(test_gsea,gene_name!='CNOT1')
ranks <- deframe(test_gsea[,c('gene_name','log2FC_halflife')])
head(ranks, 20)


# Load the pathways into a named list
pathwaystouse <- gmtPathways("c5.bp.v7.1.symbols.gmt")
xlab<-'Biological Process'

pathwaystouse <- gmtPathways("c5.mf.v7.1.symbols.gmt")
xlab<-'Molecular Function'

pathwaystouse <- gmtPathways("c5.cc.v7.1.symbols.gmt")
xlab<-'Cellular Component'

pathwaystouse <- gmtPathways("h.all.v7.1.symbols.gmt")
xlab<-'Hallmark Pathways'

#run analysis and plot significant terms
fgseaRes <- fgsea(pathways=pathwaystouse, stats=ranks, nperm=5000,minSize=50,maxSize=500)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

#
ggplot(subset(fgseaResTidy,padj<0.05), aes(reorder(pathway, NES), NES))+scale_fill_gradient2(low='red3',mid='white',high='green4',breaks=c(-3,0,3),labels=c(-3,0,3)) +
  geom_col(aes(fill=NES))+coord_flip() +labs(x=xlab, y="Normalized Enrichment Score")+ggtitle(xlab)+theme_minimal()+theme(axis.title.y=element_text(size=18),plot.title = element_text(size=16))

#ff<-as.data.frame(fgseaResTidy[,c(1:5,7)])
#write.table(x=subset(ff,padj<0.05),file = paste0(xlab,'_log2FChalflife.csv'),col.names=TRUE,row.names=FALSE,sep=',')
#rm(ff)

ff<-read.table('Biological Process_lo2FChalflife.csv',sep=',',stringsAsFactors = FALSE,header=TRUE)
fgseaResTidy<-ff

###########################################################################
###########################################################################
###########################################################################

ff<-as.data.frame(fgseaResTidy[,c(1:5,7)])
nrow(ff)

#RPFcomp CC
remove_terms<-c('GO_BIOLOGICAL_ADHESION','GO_ION_TRANSMEMBRANE_TRANSPORT','GO_SYNAPSE_ORGANIZATION','GO_SECRETORY_GRANULE','GO_ANCHORING_JUNCTION','GO_INTRINSIC_COMPONENT_OF_PLASMA_MEMBRANE','GO_CELL_SUBSTRATE_JUNCTION','GO_POSTSYNAPSE','GO_ENDOPLASMIC_RETICULUM_LUMEN')
ff<-subset(ff,(pathway %in% remove_terms)==FALSE)
nrow(ff)

#
#RPFcomp BP
remove_terms<-c('GO_BIOLOGICAL_ADHESION','GO_ION_TRANSMEMBRANE_TRANSPORT','GO_CELL_MORPHOGENESIS_INVOLVED_IN_DIFFERENTIATION','GO_CATION_TRANSMEMBRANE_TRANSPORT','GO_CELLULAR_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND','GO_CELLULAR_RESPONSE_TO_NITROGEN_COMPOUND','GO_EXOCYTOSIS','GO_CELLULAR_ION_HOMEOSTASIS','GO_CELL_CELL_SIGNALING_BY_WNT','GO_CELL_MORPHOGENESIS_INVOLVED_IN_NEURON_DIFFERENTIATION','GO_EMBRYONIC_MORPHOGENESIS','GO_TRANSMEMBRANE_RECEPTOR_PROTEIN_TYROSINE_KINASE_SIGNALING_PATHWAY','GO_AXON_DEVELOPMENT','GO_NEURON_DEVELOPMENT')
ff<-subset(ff,(pathway %in% remove_terms)==FALSE)
nrow(ff)



#
ggplot(subset(ff,padj<0.05), aes(reorder(pathway, NES), NES))+scale_fill_gradient2(low='red3',mid='white',high='green4',breaks=c(-3,0,3),labels=c(-3,0,3))+
  geom_col(aes(fill=NES))+coord_flip()+theme(panel.grid.minor=element_blank(),panel.grid.minor.x=element_blank())+theme(axis.title=element_text(size=12),axis.text=element_text(size=12)) +labs(x=xlab, y="Normalized Enrichment Score")+ggtitle(xlab)+theme_minimal()+theme(axis.title.y=element_text(size=18),plot.title = element_text(size=16))



