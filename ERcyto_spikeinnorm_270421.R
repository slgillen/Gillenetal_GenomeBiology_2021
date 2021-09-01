
spikeins<-readDNAStringSet('ERCC_spikein.fa')
seq_name = names(spikeins)
sequence = paste(spikeins)
spikeindf<-data.frame(seq_name, sequence,stringsAsFactors = FALSE)


sample_names<-c('Rep1_siControl_cyto','Rep2_siControl_cyto','Rep3_siControl_cyto','Rep1_siCNOT1_cyto','Rep2_siCNOT1_cyto','Rep3_siCNOT1_cyto','Rep1_siControl_ER','Rep2_siControl_ER','Rep3_siControl_ER','Rep1_siCNOT1_ER','Rep2_siCNOT1_ER','Rep3_siCNOT1_ER')

samples<-list()
for(i in 1:length(sample_names)){
 xx<-read.table(paste0('fcounts_',sample_names[i],'_matrix.txt'),stringsAsFactors = FALSE,header=TRUE) #inputs are formatted outputs from featureCounts
 names(xx)<-c('ENSG',sample_names[i])
 samples[[i]]<-xx
 rm(xx)
}

ERcyto_readcounts<-samples[[1]]
for(i in 2:length(samples)){
  ERcyto_readcounts<-merge(ERcyto_readcounts,samples[[i]],by='ENSG')
}
dim(ERcyto_readcounts)
ERcyto_readcounts<-subset(ERcyto_readcounts,Rep1_siControl_cyto>5 & Rep2_siControl_cyto>5 & Rep3_siControl_cyto>5) #pre-filtering to remove genes with no/very low reads
nrow(ERcyto_readcounts)
names(ERcyto_readcounts)


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



