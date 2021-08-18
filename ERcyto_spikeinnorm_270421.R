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
#ratios_nanodrop<-c(0.38395,0.38026,0.44511,0.29480,0.31888,0.39918) #for keeping cyto constant
ratios_nanodrop<-c(2.604517342,2.629801892,2.246651786,3.392073996,3.136024986,2.505110926) #for keeping ER constant & * cyto counts by
ratio_normfactors
ratios_nanodrop

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


#normalise CYTO samples to ER/cyto ratio nanodrop
ERcyto_readcounts_nano<-ERcyto_readcounts2
for(i in 2:7){
  ERcyto_readcounts_nano[,i]<-ERcyto_readcounts_nano[,i]*ratios_nanodrop[i-1]
}
apply(ERcyto_readcounts_nano[,2:13],2,function(x) sum(x))




#MDS plot of data - sanity check of the data spike-in ratio norm
plotMDS(as.matrix(ERcyto_readcounts_SI[,c(2:13)])) #all samples
plotMDS(as.matrix(ERcyto_readcounts_SI[,c(2:7)])) #cyto
plotMDS(as.matrix(ERcyto_readcounts_SI[,c(8:13)])) #ER
plotMDS(as.matrix(ERcyto_readcounts_SI[,c(2,3,4,8,9,10)])) #control
plotMDS(as.matrix(ERcyto_readcounts_SI[,c(5,6,7,11,12,13)])) #KD


#MDS plot of data - sanity check of the data nanodrop ratio norm
plotMDS(as.matrix(ERcyto_readcounts_nano[,c(2:13)])) #all samples
plotMDS(as.matrix(ERcyto_readcounts_nano[,c(2:7)])) #cyto
plotMDS(as.matrix(ERcyto_readcounts_nano[,c(8:13)])) #ER
plotMDS(as.matrix(ERcyto_readcounts_nano[,c(2,3,4,8,9,10)])) #control
plotMDS(as.matrix(ERcyto_readcounts_nano[,c(5,6,7,11,12,13)])) #KD



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
kmeans_model_3<-kmeans(kmeans_test,3)
kmeans_model_4<-kmeans(kmeans_test,4)

kmeans_test$cluster_k2<-kmeans_model_2$cluster
kmeans_test$cluster_k3<-kmeans_model_3$cluster
kmeans_test$cluster_k4<-kmeans_model_4$cluster

kmeans_test$cluster_k2<-factor(kmeans_test$cluster_k2,levels=c('2','1'))
kmeans_test$cluster_k3<-factor(kmeans_test$cluster_k3,levels=c('1','2','3'))
kmeans_test$cluster_k4<-factor(kmeans_test$cluster_k4,levels=c('1','2','3','4'))

kmeans_test$ENSG<-rownames(kmeans_test)
kmeans_merge<-merge(kmeans_test[,4:7],ERcyto_log2FCs_SI,by='ENSG')
kmeans_merge<-kmeans_merge[order(kmeans_merge$cluster_k2),]

f1<-ggplot(kmeans_merge,aes(x=Ctrl_cyto_average,y=Ctrl_ER_average,col=cluster_k2))+geom_point(size=0.5)+theme_classic()+
  scale_x_log10()+scale_y_log10()+xlab('cyto CPM')+ylab('ER CPM')+scale_colour_manual(values=c('grey68','royalblue3'))+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
  plot.title=element_text(size=16),legend.title=element_blank())
ggsave(plot=f1,filename='control_ERcyto.png',height=4,width=4.5)



#####################################################################################################################################


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




