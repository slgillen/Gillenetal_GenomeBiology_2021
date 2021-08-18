library(fgsea)
library(data.table)
library(ggplot2)
library(tibble)
library(tidyr)
library(autoimage)
library(dplyr)
library(broom)

#get data
RPF27to31<-read.table('~/cnotdata/paperplots/SILAC/RPFCDSfilt_27to31_withgroups.txt',stringsAsFactors = FALSE,header=TRUE,sep='\t')
nrow(RPF27to31)
notRPGs<-read.delim('Z:/Systems/Sarah/CNOTdataAnalysis/mRNAs_noRPGsnoHistones.txt',header=TRUE,stringsAsFactors = FALSE,sep='\t')
RPF27to31<-subset(RPF27to31, (Gene_ID %in% notRPGs$Gene_ID)==TRUE)
nrow(RPF27to31)
RPF27to31<-subset(RPF27to31,startsWith(RPF27to31$Gene_ID,'RPL')==FALSE)
nrow(RPF27to31)
RPF27to31<-subset(RPF27to31,startsWith(RPF27to31$Gene_ID,'RPS')==FALSE)
nrow(RPF27to31)
RPF27to31<-subset(RPF27to31,startsWith(RPF27to31$Gene_ID,'MRP')==FALSE)
nrow(RPF27to31)
RPF27to31<-subset(RPF27to31,startsWith(RPF27to31$Gene_ID,'RP11-')==FALSE)
nrow(RPF27to31)
RPF27to31$RPF_TOT_DIFF<-RPF27to31$RPFsCDSfilt_27to31_log2FoldChange-RPF27to31$Total_log2FoldChange
#

RPF27to31<-merge(RPF27to31,humanDF[,c(5,10:12)],by.x='Gene_ID',by.y='Gene_ID')
nrow(RPF27to31) 
RPF27to31<-subset(RPF27to31,nchar(CDS_sequence)>200)
nrow(RPF27to31)
RPF27to31<-subset(RPF27to31,nchar(CDS_sequence)<10000)
nrow(RPF27to31)

#focus on groupType2 TE up and TE down groups to streamline number of significant terms
RPF27to31<-RPF27to31[,c(1:6)]
nrow(RPF27to31)
names(RPF27to31)<-c('Gene_ID','RNA_log2FC','RNA_padj','RPF_log2FC','RPF_padj','Group')

RPF27to31$RPF_TOT_DIFF<-RPF27to31$RPF_log2FC-RPF27to31$RNA_log2FC
RPF27to31$GroupType2<-rep('all other mRNAs',nrow(RPF27to31))

#three group version
for(i in 1:nrow(RPF27to31)){
  if((RPF27to31[i,'RPF_TOT_DIFF']<0.1) & (RPF27to31[i,'RPF_TOT_DIFF']>(-0.1))){ #or 0.05
    RPF27to31[i,'GroupType2']<-'no TE change'
  }
  if((RPF27to31[i,'RPF_TOT_DIFF']>0.2)){
    RPF27to31[i,'GroupType2']<-'increased TE'
  }
  if((RPF27to31[i,'RPF_TOT_DIFF']<(-0.2))){
    RPF27to31[i,'GroupType2']<-'decreased TE'
  }
}


#####################

#could try rank by p-value (giving +ve and -ve to p-val of direction of change??)
setwd('\\\\john-doe/gw/Systems/Sarah/CNOTpaper/gene_ontology_100920')



test_gsea<-subset(RPF27to31,GroupType2=='increased TE' | GroupType2=='decreased TE')
test_gsea<-RPF27to31

nrow(test_gsea)
test_gsea<-test_gsea[order(test_gsea$RPF_TOT_DIFF,decreasing=TRUE),]
test_gsea<-subset(test_gsea,Gene_ID!='CNOT1')
ranks <- deframe(test_gsea[,c('Gene_ID','RPF_TOT_DIFF')])
head(ranks, 20)


# Load the pathways into a named list
pathwaystouse <- gmtPathways("\\\\john-doe/gw/Systems/Sarah/CNOTpaper/gene_ontology/c5.bp.v7.1.symbols.gmt")
xlab<-'Biological Process'

pathwaystouse <- gmtPathways("\\\\john-doe/gw/Systems/Sarah/CNOTpaper/gene_ontology/c5.mf.v7.1.symbols.gmt")
xlab<-'Molecular Function'

pathwaystouse <- gmtPathways("\\\\john-doe/gw/Systems/Sarah/CNOTpaper/gene_ontology/c5.cc.v7.1.symbols.gmt")
xlab<-'Cellular Component'



#run analysis and plot significant terms
fgseaRes <- fgsea(pathways=pathwaystouse, stats=ranks, nperm=5000,minSize=5,maxSize=500)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

#
ggplot(subset(fgseaResTidy,padj<0.05), aes(reorder(pathway, NES), NES))+scale_fill_gradient2(low='red3',mid='white',high='green4',breaks=c(-3,0,3),labels=c(-3,0,3)) +
  geom_col(aes(fill=NES))+coord_flip() +labs(x=xlab, y="Normalized Enrichment Score")+ggtitle(xlab)+theme_minimal()+theme(axis.title.y=element_text(size=18),plot.title = element_text(size=16))

ff<-as.data.frame(fgseaResTidy[,c(1:5,7)])
#write.table(x=subset(ff,padj<0.05),file = paste0(xlab,'_RPFTOTDIFF_X.csv'),col.names=TRUE,row.names=FALSE,sep=',')
rm(ff)


ff<-read.table('Biological Process_RPFTOTDIFF_filtered_v2.csv',sep=',',stringsAsFactors = FALSE,header=TRUE)
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



