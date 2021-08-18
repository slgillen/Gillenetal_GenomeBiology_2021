#TPMs and spike-in calculations


spikeins<-readDNAStringSet('//john-doe/gw/Systems/Sarah/LabData/CNOTpaperRevisions/ERcyto_data/ERCC_spikein.fa')
seq_name = names(spikeins)
sequence = paste(spikeins)
spikeindf<-data.frame(seq_name, sequence,stringsAsFactors = FALSE)

#sort first part of ENSG before . in gene ID merge

#try CPM

#try TPM

#try where to add/remove the spike-in in the adjustments

#check Chiara RNA concentrations between conditions



###################################
setwd("//john-doe/gw/Systems/Sarah/LabData/CNOTpaperRevisions/ERcyto_data/fcounts")

sample_names<-c('A_Ctrl_Cyto','C_Ctrl_Cyto','E_Ctrl_Cyto','A_KD_Cyto','C_KD_Cyto','E_KD_Cyto','A_Ctrl_ER','C_Ctrl_ER','E_Ctrl_ER','A_KD_ER','C_KD_ER','E_KD_ER')

samples<-list()
for(i in 1:length(sample_names)){
  xx<-read.table(paste0('fcounts_',sample_names[i],'_matrix.txt'),stringsAsFactors = FALSE,header=TRUE)
  names(xx)<-c('ENSG',sample_names[i])
  samples[[i]]<-xx
  rm(xx)
}

ERcyto_readcounts<-samples[[1]]
for(i in 2:length(samples)){
  ERcyto_readcounts<-merge(ERcyto_readcounts,samples[[i]],by='ENSG')
}
dim(ERcyto_readcounts)
ERcyto_readcounts<-subset(ERcyto_readcounts,A_Ctrl_Cyto>0 & C_Ctrl_Cyto>0 & E_Ctrl_Cyto>0)
nrow(ERcyto_readcounts)
names(ERcyto_readcounts)
apply(ERcyto_readcounts[,2:13],2,function(x) sum(x))



#separate the spike-ins v2
spikeins<-subset(ERcyto_readcounts,startsWith(ENSG,'ERCC-')==TRUE)
nrow(spikeins)
apply(spikeins[,2:13],2,function(x) sum(x))
spikein_normfactors<-apply(spikeins[,2:13],2,function(x) sum(x))
spikein_normfactors2<-apply(ERcyto_readcounts[,2:13],2,function(x) sum(x))
spikein_normfactors
spikein_normfactors2
spikein_normfactors<-(spikein_normfactors/spikein_normfactors2)*100 #get proportional counts
spikein_normfactors #these are the numbers to divide the read counts by to get quantitative between ER and cyto within a condition

#cyto/ER ratio of spike-in
ratio_normfactors<-spikein_normfactors[1:6]/spikein_normfactors[7:12]
ratio_normfactors

#remove spike-ins
ERcyto_readcounts2<-subset(ERcyto_readcounts,startsWith(ENSG,'ERCC-')==FALSE)
apply(ERcyto_readcounts2[,2:13],2,function(x) sum(x))

#change to CPM
for(i in 2:ncol(ERcyto_readcounts2)){
  ERcyto_readcounts2[,i]<-(ERcyto_readcounts2[,i]/sum(ERcyto_readcounts2[,i]))*1000000
}
apply(ERcyto_readcounts2[,2:13],2,function(x) sum(x))

#normalise CYTO samples to ER/cyto ratio spike-in
ERcyto_readcounts_SI<-ERcyto_readcounts2
for(i in 2:7){
  ERcyto_readcounts_SI[,i]<-ERcyto_readcounts_SI[,i]/ratio_normfactors[i-1]
}
apply(ERcyto_readcounts_SI[,2:13],2,function(x) sum(x))

#MDS plot of data - sanity check of the data spike-in ratio norm
plotMDS(as.matrix(ERcyto_readcounts_SI[,c(2:13)])) #all samples
plotMDS(as.matrix(ERcyto_readcounts_SI[,c(2:7)])) #cyto
plotMDS(as.matrix(ERcyto_readcounts_SI[,c(8:13)])) #ER
plotMDS(as.matrix(ERcyto_readcounts_SI[,c(2,3,4,8,9,10)])) #control
plotMDS(as.matrix(ERcyto_readcounts_SI[,c(5,6,7,11,12,13)])) #KD




###############################################################################################
################################### based on spike-in #########################################
###############################################################################################


#per replicate log2FC change
ERcyto_log2FCs_SI<-ERcyto_readcounts_SI
#names(ERcyto_log2FCs_SI)<-c('ENSG','A_Ctrl_ERcyto','C_Ctrl_ERcyto','E_Ctrl_ERcyto','A_KD_ERcyto','C_KD_ERcyto','E_KD_ERcyto','Ctrl_ERcyto_average','KD_ERcyto_average')

ERcyto_log2FCs_SI$ENSG<-ERcyto_readcounts_SI$ENSG
ERcyto_log2FCs_SI$A_Ctrl_ERcyto<-log(ERcyto_readcounts_SI$A_Ctrl_ER/ERcyto_readcounts_SI$A_Ctrl_Cyto,2)
ERcyto_log2FCs_SI$C_Ctrl_ERcyto<-log(ERcyto_readcounts_SI$C_Ctrl_ER/ERcyto_readcounts_SI$C_Ctrl_Cyto,2)
ERcyto_log2FCs_SI$E_Ctrl_ERcyto<-log(ERcyto_readcounts_SI$E_Ctrl_ER/ERcyto_readcounts_SI$E_Ctrl_Cyto,2)

ERcyto_log2FCs_SI$A_KD_ERcyto<-log(ERcyto_readcounts_SI$A_KD_ER/ERcyto_readcounts_SI$A_KD_Cyto,2)
ERcyto_log2FCs_SI$C_KD_ERcyto<-log(ERcyto_readcounts_SI$C_KD_ER/ERcyto_readcounts_SI$C_KD_Cyto,2)
ERcyto_log2FCs_SI$E_KD_ERcyto<-log(ERcyto_readcounts_SI$E_KD_ER/ERcyto_readcounts_SI$E_KD_Cyto,2)

ERcyto_log2FCs_SI$A_cyto_KDtoCtrl<-log(ERcyto_readcounts_SI$A_KD_Cyto/ERcyto_readcounts_SI$A_Ctrl_Cyto,2)
ERcyto_log2FCs_SI$C_cyto_KDtoCtrl<-log(ERcyto_readcounts_SI$C_KD_Cyto/ERcyto_readcounts_SI$C_Ctrl_Cyto,2)
ERcyto_log2FCs_SI$E_cyto_KDtoCtrl<-log(ERcyto_readcounts_SI$E_KD_Cyto/ERcyto_readcounts_SI$E_Ctrl_Cyto,2)

ERcyto_log2FCs_SI$A_ER_KDtoCtrl<-log(ERcyto_readcounts_SI$A_KD_ER/ERcyto_readcounts_SI$A_Ctrl_ER,2)
ERcyto_log2FCs_SI$C_ER_KDtoCtrl<-log(ERcyto_readcounts_SI$C_KD_ER/ERcyto_readcounts_SI$C_Ctrl_ER,2)
ERcyto_log2FCs_SI$E_ER_KDtoCtrl<-log(ERcyto_readcounts_SI$E_KD_ER/ERcyto_readcounts_SI$E_Ctrl_ER,2)

for(i in 1:nrow(ERcyto_log2FCs_SI)){
  ERcyto_log2FCs_SI[i,'Ctrl_ER_average']<-mean(c(ERcyto_log2FCs_SI[i,'A_Ctrl_ER'],ERcyto_log2FCs_SI[i,'C_Ctrl_ER'],ERcyto_log2FCs_SI[i,'E_Ctrl_ER']),na.rm=TRUE)
}

for(i in 1:nrow(ERcyto_log2FCs_SI)){
  ERcyto_log2FCs_SI[i,'Ctrl_cyto_average']<-mean(c(ERcyto_log2FCs_SI[i,'A_Ctrl_Cyto'],ERcyto_log2FCs_SI[i,'C_Ctrl_Cyto'],ERcyto_log2FCs_SI[i,'E_Ctrl_Cyto']),na.rm=TRUE)
}

for(i in 1:nrow(ERcyto_log2FCs_SI)){
  ERcyto_log2FCs_SI[i,'KD_ER_average']<-mean(c(ERcyto_log2FCs_SI[i,'A_KD_ER'],ERcyto_log2FCs_SI[i,'C_KD_ER'],ERcyto_log2FCs_SI[i,'E_KD_ER']),na.rm=TRUE)
}

for(i in 1:nrow(ERcyto_log2FCs_SI)){
  ERcyto_log2FCs_SI[i,'KD_cyto_average']<-mean(c(ERcyto_log2FCs_SI[i,'A_KD_Cyto'],ERcyto_log2FCs_SI[i,'C_KD_Cyto'],ERcyto_log2FCs_SI[i,'E_KD_Cyto']),na.rm=TRUE)
}

##
for(i in 1:nrow(ERcyto_log2FCs_SI)){
  ERcyto_log2FCs_SI[i,'Ctrl_ERcyto_average']<-mean(c(ERcyto_log2FCs_SI[i,'A_Ctrl_ERcyto'],ERcyto_log2FCs_SI[i,'C_Ctrl_ERcyto'],ERcyto_log2FCs_SI[i,'E_Ctrl_ERcyto']),na.rm=TRUE)
}

for(i in 1:nrow(ERcyto_log2FCs_SI)){
  ERcyto_log2FCs_SI[i,'KD_ERcyto_average']<-mean(c(ERcyto_log2FCs_SI[i,'A_KD_ERcyto'],ERcyto_log2FCs_SI[i,'C_KD_ERcyto'],ERcyto_log2FCs_SI[i,'E_KD_ERcyto']),na.rm=TRUE)
}

for(i in 1:nrow(ERcyto_log2FCs_SI)){
  ERcyto_log2FCs_SI[i,'cyto_KDtoCtrl_average']<-mean(c(ERcyto_log2FCs_SI[i,'A_cyto_KDtoCtrl'],ERcyto_log2FCs_SI[i,'C_cyto_KDtoCtrl'],ERcyto_log2FCs_SI[i,'E_cyto_KDtoCtrl']),na.rm=TRUE)
}

for(i in 1:nrow(ERcyto_log2FCs_SI)){
  ERcyto_log2FCs_SI[i,'ER_KDtoCtrl_average']<-mean(c(ERcyto_log2FCs_SI[i,'A_ER_KDtoCtrl'],ERcyto_log2FCs_SI[i,'C_ER_KDtoCtrl'],ERcyto_log2FCs_SI[i,'E_ER_KDtoCtrl']),na.rm=TRUE)
}

ERcyto_log2FCs_SI$KDtoCtrl<-ERcyto_log2FCs_SI$KD_ERcyto_average-ERcyto_log2FCs_SI$Ctrl_ERcyto_average
ERcyto_log2FCs_SI$ERtoCYTO<-ERcyto_log2FCs_SI$ER_KDtoCtrl_average-ERcyto_log2FCs_SI$cyto_KDtoCtrl_average
nrow(ERcyto_log2FCs_SI)
ERcyto_log2FCs_SI_filt<-subset(ERcyto_log2FCs_SI,is.finite(KDtoCtrl)==TRUE)
nrow(ERcyto_log2FCs_SI_filt)


checks<-ERcyto_log2FCs_SI[,c(1,2)]
checks$group<-rep('A_Ctrl_ERcyto',nrow(checks))
names(checks)[2]<-'ERcyto'
for(i in 3:7){
  checks2<-ERcyto_log2FCs_SI[,c(1,i)]
  checks2$group<-rep(names(ERcyto_log2FCs_SI[i]),nrow(checks2))
  names(checks2)[2]<-'ERcyto'
  checks<-rbind(checks,checks2)
  rm(checks2)
}

ggplot(checks,aes(x=group,y=ERcyto))+geom_violin(draw_quantiles=c(0.25,0.5,0.75))+theme_classic()



ggplot(ERcyto_readcounts_SI,aes(x=A_Ctrl_Cyto,y=A_Ctrl_ER))+geom_point()+theme_classic()+
  scale_x_log10()+scale_y_log10()

ggplot(ERcyto_readcounts_SI,aes(x=C_Ctrl_Cyto,y=C_Ctrl_ER))+geom_point()+theme_classic()+
  scale_x_log10()+scale_y_log10()

ggplot(ERcyto_readcounts_SI,aes(x=E_Ctrl_Cyto,y=E_Ctrl_ER))+geom_point()+theme_classic()+
  scale_x_log10()+scale_y_log10()


#k-means clustering with 2 clusters to identify ER-target mRNAs
##----------------------------------------------------------------------
## K-means CLUSTERING --------------------------------------------------
##----------------------------------------------------------------------
kmeans_test<-ERcyto_log2FCs_SI_filt[,c(14,15,16)]
rownames(kmeans_test)<-ERcyto_log2FCs_SI_filt$ENSG
nrow(kmeans_test)
kmeans_test<-subset(kmeans_test,(rownames(kmeans_test) %in% ERcyto$ENSG)==TRUE)
nrow(kmeans_test)

values<-1:10
sum_of_squares<-c()
for (k in values){
  kmeans_model<-kmeans(kmeans_test,k)
  within_cluster<-sum(kmeans_model$withinss)
  sum_of_squares[k]<-within_cluster
  
}
x=1:10
y=sum_of_squares
plot(x,y)

kmeans_model_2<-kmeans(kmeans_test,2)

kmeans_test$cluster_k2<-kmeans_model_2$cluster

kmeans_test$cluster_k2<-factor(kmeans_test$cluster_k2,levels=c('1','2'))

kmeans_test$ENSG<-rownames(kmeans_test)
kmeans_merge<-merge(kmeans_test[,c(4,7)],ERcyto_log2FCs_SI,by='ENSG')
kmeans_merge<-kmeans_merge[order(kmeans_merge$cluster_k2),]

f1<-ggplot(kmeans_merge,aes(x=Ctrl_cyto_average,y=Ctrl_ER_average,col=cluster_k2))+geom_point(size=0.3)+theme_classic()+
  scale_x_log10()+scale_y_log10()+xlab('cyto CPM')+ylab('ER CPM')+scale_colour_manual(values=c('grey68','royalblue3'))+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
  plot.title=element_text(size=16),legend.title=element_blank())
ggsave(plot=f1,filename='control_ERcyto.png',height=3.3,width=4.1)

GID<-cyto[,c(1,2)]
RPF27to31<-merge(GID,RPF27to31,by='Gene_ID')
kmeans_merge<-merge(kmeans_merge,RPF27to31,by='ENSG')

#which of TE groups are ER-targets?
nrow(subset(kmeans_merge,GroupType2=='increased TE' & cluster_k2=='2'))
nrow(subset(kmeans_merge,GroupType2=='increased TE' & cluster_k2=='1'))

nrow(subset(kmeans_merge,GroupType2=='no TE change' & cluster_k2=='2'))
nrow(subset(kmeans_merge,GroupType2=='no TE change' & cluster_k2=='1'))

nrow(subset(kmeans_merge,GroupType2=='decreased TE' & cluster_k2=='2'))
nrow(subset(kmeans_merge,GroupType2=='decreased TE' & cluster_k2=='1'))

nrow(subset(kmeans_merge,GroupType2=='all other mRNAs' & cluster_k2=='2'))
nrow(subset(kmeans_merge,GroupType2=='all other mRNAs' & cluster_k2=='1'))

#####################################################################################################################################
plotdata<-merge(M[,c(1,10)],kmeans_merge[,c(37,2,39)],by.x='gene_name',by.y='Gene_ID')

f1<-ggplot(plotdata,aes(x=cluster_k2,y=mean_FR,fill=cluster_k2))+geom_boxplot(outlier.size=0.5)+theme_classic()+
  ylab('log2FC protein production (siCNOT1/siControl)')+scale_fill_manual(values=c('grey68','royalblue3'))+coord_trans(ylim=c(-1.25,1.25))+
  theme(axis.text.x=element_blank(),axis.title=element_blank(),axis.ticks.x=element_blank())
ggsave(plot = f1,filename='ER_SILAC.png',height=3,width=2.5)

#stats
KT<-kruskal.test(mean_FR~cluster_k2,data=plotdata) 


plotdata<-subset(plotdata,GroupType2!='all other mRNAs')
plotdata$GroupType2<-factor(plotdata$GroupType2,levels=c('decreased TE','no TE change','increased TE'))
plotdata<-plotdata[order(plotdata$GroupType2),]
f1<-ggplot(plotdata,aes(x=GroupType2,y=mean_FR,fill=GroupType2))+geom_boxplot(outlier.size=0.5)+theme_classic()+
  scale_fill_manual(values=c('royalblue3','goldenrod2','red3'))+coord_trans(ylim=c(-1.25,1.25))+
  theme(axis.text.x=element_blank(),axis.title=element_blank(),axis.ticks.x=element_blank())
ggsave(plot = f1,filename='TEgroups_SILAC.png',height=3,width=2.5)


setwd("//john-doe/gw/Systems/Sarah/LabData/CNOTpaperRevisions/ERcyto_data")

cyto<-read.table('Cyto_KDvCtrl_DESeq2output.csv',stringsAsFactors = FALSE,header=TRUE,sep=',')
names(cyto)[3:7]<-paste('cyto_',names(cyto)[3:7],sep='')

ER<-read.table('ER_KDvCtrl_DESeq2output.csv',stringsAsFactors = FALSE,header=TRUE,sep=',')
names(ER)[3:7]<-paste('ER_',names(ER)[3:7],sep='')


cyto<-merge(cyto,kmeans_merge[,c(1,2,37,44)],by='ENSG')
dim(cyto)
cyto<-subset(cyto, GroupType2!='all other mRNAs')
cyto$GroupType2<-factor(cyto$GroupType2,levels=c('decreased TE','no TE change','increased TE'))
cyto<-cyto[order(cyto$GroupType2),]

f1<-ggplot(cyto,aes(x=cluster_k2,y=cyto_log2FoldChange,fill=cluster_k2))+geom_boxplot(outlier.size=0.5)+theme_classic()+
  scale_fill_manual(values=c('grey68','royalblue3'))+coord_trans(ylim=c(-1.25,1.25))+
  theme(axis.text.x=element_blank(),axis.title=element_blank(),axis.ticks.x=element_blank())
ggsave(plot = f1,filename='cyto_RNAseq.png',height=3,width=2.5)

#stats
KT<-kruskal.test(cyto_log2FoldChange~cluster_k2,data=cyto) 

f1<-ggplot(cyto,aes(x=GroupType2,y=cyto_log2FoldChange,fill=GroupType2))+geom_boxplot(outlier.size=0.5)+theme_classic()+
  scale_fill_manual(values=c('royalblue3','goldenrod2','red3'))+coord_trans(ylim=c(-1.25,1.25))+
  theme(axis.text.x=element_blank(),axis.title=element_blank(),axis.ticks.x=element_blank())
ggsave(plot = f1,filename='cyto_RNAseq_TEgroups.png',height=3,width=2.5)


ER<-merge(ER,kmeans_merge[,c(1,2,37,44)],by='ENSG')
dim(ER)
ER<-subset(ER, GroupType2!='all other mRNAs')
ER$GroupType2<-factor(ER$GroupType2,levels=c('decreased TE','no TE change','increased TE'))
ER<-ER[order(ER$GroupType2),]

f1<-ggplot(ER,aes(x=cluster_k2,y=ER_log2FoldChange,fill=cluster_k2))+geom_boxplot(outlier.size=0.5)+theme_classic()+
  scale_fill_manual(values=c('royalblue3','goldenrod2','red3'))+coord_trans(ylim=c(-1.25,1.25))+
  theme(axis.text.x=element_blank(),axis.title=element_blank(),axis.ticks.x=element_blank())
ggsave(plot = f1,filename='ER_RNAseq.png',height=3,width=2.5)

#stats
KT<-kruskal.test(ER_log2FoldChange~cluster_k2,data=ER) 

f1<-ggplot(ER,aes(x=GroupType2,y=ER_log2FoldChange,fill=GroupType2))+geom_boxplot(outlier.size=0.5)+theme_classic()+
  scale_fill_manual(values=c('royalblue3','goldenrod2','red3'))+coord_trans(ylim=c(-1.25,1.25))+
  theme(axis.text.x=element_blank(),axis.title=element_blank(),axis.ticks.x=element_blank())
ggsave(plot = f1,filename='ER_RNAseq_TEgroups.png',height=3,width=2.5)




######################
plotdata<-ERcyto_log2FCs_SI_filt
nrow(plotdata)
plotdata$ENSG_short<-sapply(plotdata$ENSG,function(x) strsplit(x,split='[.]')[[1]][1])

GeneIDs$ENSG_short<-sapply(GeneIDs$ENSG,function(x) strsplit(x,split='[.]')[[1]][1])
GeneIDs<-GeneIDs[!duplicated(GeneIDs$ENSG_short),]

plotdata<-merge(GeneIDs[,2:3],plotdata,by='ENSG_short') #check how much CPM loss in Gene ID loss - go back to fcounts to sort? g27 v g28?
nrow(plotdata)

plotdata<-merge(plotdata,kmeans_test[,c('ENSG','cluster_k2')],by='ENSG')
nrow(plotdata)
#ERcluster_genes<-subset(plotdata,cluster_k2=='2')$Gene_ID
#length(ERcluster_genes)
plotdata$cluster_k2<-factor(plotdata$cluster_k2,levels=c('2','1'))
plotdata<-subset(plotdata,startsWith(plotdata$Gene_ID,'MT-')==FALSE)
nrow(plotdata)


##
f1<-ggplot(kmeans_merge,aes(x=Ctrl_cyto_average,y=Ctrl_ER_average,col=cluster_k2))+geom_point(size=0.3)+theme_classic()+
  scale_x_log10()+scale_y_log10()+xlab('cyto CPM')+ylab('ER CPM')+scale_colour_manual(values=c('grey68','royalblue3'))+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
ggsave(plot=f1,filename='control_ERcyto.png',height=4,width=4.5)

f1<-ggplot(kmeans_merge,aes(x=Ctrl_ER_average,y=KD_ER_average,col=cluster_k2))+geom_point(size=0.75)+theme_classic()+
  scale_x_log10()+scale_y_log10()+xlab('Control ER CPM')+ylab('KD ER CPM')+scale_colour_manual(values=c('grey68','royalblue3'))+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
ggsave(plot=f1,filename='ER_KDvCtrl.png',height=4,width=4.5)

f1<-ggplot(kmeans_merge,aes(x=Ctrl_cyto_average,y=KD_cyto_average,col=cluster_k2))+geom_point(size=0.75)+theme_classic()+
  scale_x_log10()+scale_y_log10()+xlab('Control cyto CPM')+ylab('KD cyto CPM')+scale_colour_manual(values=c('grey68','royalblue3'))+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
ggsave(plot=f1,filename='cyto_KDvCtrl.png',height=4,width=4.5)

f1<-ggplot(kmeans_merge,aes(x=KD_cyto_average,y=KD_ER_average,col=cluster_k2))+geom_point(size=0.75)+theme_classic()+
  scale_x_log10()+scale_y_log10()+xlab('Control cyto CPM')+ylab('KD cyto CPM')+scale_colour_manual(values=c('grey68','royalblue3'))+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
ggsave(plot=f1,filename='KD_ERcyto.png',height=4,width=4.5)

##
colourset<-c('grey68','royalblue3')
f1<-ggplot(plotdata, aes(y=cyto_KDtoCtrl_average,x=cluster_k2,fill=cluster_k2))+geom_violin(draw_quantiles=c(0.25,0.5,0.75))+scale_fill_manual(values=colourset)+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text.x = element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank(),
                        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12),
                        plot.title=element_text(size=16),legend.title=element_blank())+coord_trans(limy=c(-4,4))
ggsave(plot=f1,filename='cyto_KDtoCtrl_average.png',height=4,width=3)

colourset<-c('grey68','royalblue3')
f1<-ggplot(plotdata, aes(y=ER_KDtoCtrl_average,x=cluster_k2,fill=cluster_k2))+geom_violin(draw_quantiles=c(0.25,0.5,0.75))+scale_fill_manual(values=colourset)+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text.x = element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank(),
                        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12),
                        plot.title=element_text(size=16),legend.title=element_blank())+coord_trans(limy=c(-4,4))
ggsave(plot=f1,filename='ER_KDtoCtrl_average.png',height=4,width=3)

colourset<-c('grey68','royalblue3')
f1<-ggplot(plotdata, aes(y=Ctrl_ERcyto_average,x=cluster_k2,fill=cluster_k2))+geom_violin(draw_quantiles=c(0.25,0.5,0.75))+scale_fill_manual(values=colourset)+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text.x = element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank(),
                        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12),
                        plot.title=element_text(size=16),legend.title=element_blank())+coord_trans(limy=c(-4,4))
ggsave(plot=f1,filename='Ctrl_ERcyto_average.png',height=4,width=3)

colourset<-c('grey68','royalblue3')
f1<-ggplot(plotdata, aes(y=KD_ERcyto_average,x=cluster_k2,fill=cluster_k2))+geom_violin(draw_quantiles=c(0.25,0.5,0.75))+scale_fill_manual(values=colourset)+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text.x = element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank(),
                        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12),
                        plot.title=element_text(size=16),legend.title=element_blank())+coord_trans(limy=c(-4,4))
ggsave(plot=f1,filename='KD_ERcyto_average.png',height=4,width=3)


######## now do boxplots for gene groups ###########


#Reid 2011 data
Reid2011<-read.table('//john-doe/gw/Systems/Sarah/CNOTpaper/SRPinvestigations/Reid_2011_ER_Cyto.csv',sep=',',stringsAsFactors = FALSE,header=TRUE)
unique(Reid2011$Type)

plotdata<-merge(plotdata,Reid2011[,c(2,16)],by.x='Gene_ID',by.y='Gene')
plotdata$Type<-factor(plotdata$Type,levels=c("","Cytosolic","ER-targeted"))
plotdata<-plotdata[order(plotdata$Type),]

colourset<-c('grey68','darkorange','royalblue3')
ggplot(plotdata, aes(y=Ctrl_ERcyto_average,x=Type,fill=Type))+geom_boxplot(alpha=0.75,outlier.size=0.75)+scale_fill_manual(values=colourset)+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text.x = element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank(),
                        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12),
                        plot.title=element_text(size=16),legend.title=element_blank())+coord_trans(limy=c(-4,4))

colourset<-c('grey68','darkorange','royalblue3')
ggplot(plotdata, aes(y=KD_ERcyto_average,x=Type,fill=Type))+geom_boxplot(alpha=0.75,outlier.size=0.75)+scale_fill_manual(values=colourset)+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text.x = element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank(),
                        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12),
                        plot.title=element_text(size=16),legend.title=element_blank())+coord_trans(limy=c(-4,4))






############### try % in each fraction for group assignment #############

plotdata$fractionofregion<-plotdata$Ctrl_ER_average/(plotdata$Ctrl_ER_average+plotdata$Ctrl_cyto_average)
plotdata$fractiongroup<-rep('other',nrow(plotdata))

for(i in 1:nrow(plotdata)){
  if(plotdata[i,'fractionofregion']>0.75){
    plotdata[i,'fractiongroup']<-'ER-targeted'
  }
  if(plotdata[i,'fractionofregion']<0.25){
    plotdata[i,'fractiongroup']<-'cyto localised'
  }
}

for(g in unique(plotdata$fractiongroup)){
  print(g)
  print(nrow(subset(plotdata,fractiongroup==g)))
}

SILAC<-merge(M,RPF27to31,by.x='gene_name',by.y='Gene_ID')
dim(SILAC)
SILAC<-merge(SILAC,plotdata[,c(3,38,40)],by.x='gene_name',by.y='Gene_ID')
nrow(SILAC)

f1<-ggplot(SILAC,aes(y=as.numeric(SILAC$mean_FR),x=cluster_k2,fill=cluster_k2))+geom_boxplot(outlier.size=0.5)+theme_classic()+scale_fill_manual(values=c('grey68','royalblue3'))
f1<-f1+ylab('average protein production')+coord_trans(limy=c(-1,1))
f1<-f1+theme(legend.text=element_text(size=12),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank(),
             axis.text.y = element_text(size = 14), axis.title = element_text(size = 16),
             plot.title=element_text(size=16),legend.title=element_blank())
f1
ggsave(plot=f1,filename='SILAC_newERcytogroups.png',width=3,height=3.5)


f1<-ggplot(SILAC,aes(y=as.numeric(SILAC$mean_FR),x=fractiongroup,fill=fractiongroup))+geom_boxplot(outlier.size=0.5)+theme_classic()+scale_fill_manual(values=c('darkorange','royalblue3','grey68'))
f1<-f1+ylab('average protein production')+coord_trans(limy=c(-1,1))
f1<-f1+theme(legend.text=element_text(size=12),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank(),
             axis.text.y = element_text(size = 14), axis.title = element_text(size = 16),
             plot.title=element_text(size=16),legend.title=element_blank())
f1
ggsave(plot=f1,filename='SILAC_newERcytogroup_v2.png',width=4,height=3.5)



#################################

#
notRPGs<-read.delim('\\\\john-doe/gw/Systems/Sarah/CNOTdata/finaldraftdata_171019/mRNAs_noRPGsnoHistones.txt',header=TRUE,stringsAsFactors = FALSE,sep='\t')
RPF27to31<-read.delim('\\\\john-doe/gw/Systems/Sarah/CNOTdata/finaldraftdata_171019/RPFCDSfilt_27to31_withgroups.txt',header=TRUE,stringsAsFactors = FALSE)
nrow(notRPGs)
nrow(RPF27to31)
RPF27to31<-subset(RPF27to31, (Gene_ID %in% notRPGs$Gene_ID)==TRUE)
nrow(RPF27to31)
RPF27to31<-subset(RPF27to31, startsWith(RPF27to31$Gene_ID,'RP11-')==FALSE)
nrow(RPF27to31)
RPF27to31$RPF_TOT_DIFF<-RPF27to31$RPFsCDSfilt_27to31_log2FoldChange-RPF27to31$Total_log2FoldChange
names(RPF27to31)
names(RPF27to31)[2]<-'RNA_log2FC'
names(RPF27to31)[4]<-'RPF_log2FC'



GeneIDs$ENSG_short<-sapply(GeneIDs$ENSG,function(x) strsplit(x,split='[.]')[[1]][1])
GeneIDs<-GeneIDs[!duplicated(GeneIDs$ENSG_short),]
ERplotdata<-merge(GeneIDs,kmeans_merge[,c(1,2,36)],by='ENSG')

ERplotdata<-merge(RPF27to31,ERplotdata,by='Gene_ID')
dim(ERplotdata)

ERplotdata$cluster_k2<-factor(ERplotdata$cluster_k2,levels=c('2','1'))
ERplotdata<-ERplotdata[order(ERplotdata$cluster_k2),]                            
colourset<-c('grey68','royalblue3')
f1<-ggplot(ERplotdata,aes(x=RNA_log2FC,y=RPF_log2FC,col=cluster_k2))+geom_point(size=0.75)+theme_classic()+xlab('log2FC RNA (siCNOT1/siControl)')+
  ylab('log2FC RPF (siCNOT1/siControl)')+geom_abline(col='black',slope=1)+
  scale_colour_manual(values=colourset)+coord_trans(ylim=c(-1.75,1.75),xlim=c(-1.75,1.75))+
  theme(legend.position='none',axis.text = element_text(size = 12),axis.title = element_text(size = 12))
f1
ggsave(plot=f1,filename='newER_inRPFTOT.png',width=3.5,height=3.5)
#stats
KT<-kruskal.test(RPF_TOT_DIFF~cluster_k2,data=ERplotdata) 



################### supplemental table - clusters, log2 cyto RNA, log2 ER RNA ######################



setwd("//john-doe/gw/Systems/Sarah/LabData/CNOTpaperRevisions/ERcyto_data")
cytoDEseq<-read.table('Cyto_KDvCtrl_DESeq2output.csv',stringsAsFactors = FALSE,header=TRUE,sep=',')
ERDEseq<-read.table('ER_KDvCtrl_DESeq2output.csv',stringsAsFactors = FALSE,header=TRUE,sep=',')

names(cytoDEseq)[3:7]<-paste0('cyto_',names(cytoDEseq)[3:7])
names(ERDEseq)[3:7]<-paste0('ER_',names(ERDEseq)[3:7])

ercytoDEseq<-merge(cytoDEseq[,c(1,2,4,7)],ERDEseq[,c(2,4,7)],by='Gene_ID')
dim(ercytoDEseq)

ercytoDEseq<-merge(ercytoDEseq,kmeans_test[,c(4,7)],by='ENSG')
dim(ercytoDEseq)

f1<-ggplot(ercytoDEseq,aes(x=cluster_k2,y=cyto_log2FoldChange,fill=cluster_k2))+geom_boxplot(outlier.size=0.5)+theme_classic()+
  scale_fill_manual(values=c('grey68','royalblue3'))+coord_trans(ylim=c(-1.25,1.25))+
  theme(axis.text.x=element_blank(),axis.title=element_blank(),axis.ticks.x=element_blank())

f1<-ggplot(ercytoDEseq,aes(x=cluster_k2,y=ER_log2FoldChange,fill=cluster_k2))+geom_boxplot(outlier.size=0.5)+theme_classic()+
  scale_fill_manual(values=c('grey68','royalblue3'))+coord_trans(ylim=c(-1.25,1.25))+
  theme(axis.text.x=element_blank(),axis.title=element_blank(),axis.ticks.x=element_blank())

forsupp<-ercytoDEseq[,c(2,7,3:6)]
forsupp$cluster<-as.character(forsupp$cluster)
forsupp$cluster2<-rep('',nrow(forsupp))
names(forsupp)<-c('gene_name','cluster','cyto_RNA_log2FC','cyto_padj','ER_RNA_log2FC','ER_padj')
for(i in 1:nrow(forsupp)){
  if(forsupp[i,'cluster']=='1'){
    forsupp[i,'cluster2']<-'ER-targeted'
  }else{
    forsupp[i,'cluster2']<-'cytosolic'
  }
}

forsupp<-forsupp[,c(1,9,3:6)]
names(forsupp)[2]<-'cluster'
write.csv(forsupp,file='Supplemental_Table_5.csv',quote=FALSE,sep=',',col.names = TRUE,row.names=FALSE)





