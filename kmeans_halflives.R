setwd("//john-doe/gw/Systems/Sarah/CNOTpaper/halflife_clusters")


#k-means clustering with 2 clusters to identify differential half-life changes
##----------------------------------------------------------------------
## K-means CLUSTERING --------------------------------------------------
##----------------------------------------------------------------------
kmeans_test<-halflife_TE[,c(10,13)]
rownames(kmeans_test)<-halflife_TE$Gene_ID

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

##
kmeans_test<-halflife_TE[,c(10,14)]
rownames(kmeans_test)<-halflife_TE$Gene_ID

head(kmeans_test)
kmeans_test <- scale(kmeans_test)

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

kmeans_model_4<-kmeans(kmeans_test,4)
kmeans_model_5<-kmeans(kmeans_test,5)
kmeans_model_6<-kmeans(kmeans_test,6)
kmeans_model_7<-kmeans(kmeans_test,7)
kmeans_model_8<-kmeans(kmeans_test,8)

kmeans_test2<-data.frame(kmeans_test)
kmeans_test2$cluster_k4<-kmeans_model_4$cluster
kmeans_test2$cluster_k4<-factor(kmeans_test2$cluster_k4,levels=c('1','2','3','4'))

kmeans_test<-halflife_TE[,c(1,10,14)]
kmeans_test$cluster_k4<-kmeans_model_4$cluster
kmeans_test$cluster_k5<-kmeans_model_5$cluster
kmeans_test$cluster_k6<-kmeans_model_6$cluster
kmeans_test$cluster_k7<-kmeans_model_7$cluster
kmeans_test$cluster_k8<-kmeans_model_8$cluster

kmeans_test$cluster_k4<-factor(kmeans_test$cluster_k4,levels=c('1','2','3','4'))
kmeans_test$cluster_k5<-factor(kmeans_test$cluster_k5,levels=c('1','2','3','4','5'))
kmeans_test$cluster_k6<-factor(kmeans_test$cluster_k6,levels=c('1','2','3','4','5','6'))
kmeans_test$cluster_k7<-factor(kmeans_test$cluster_k8,levels=c('1','2','3','4','5','6','7'))
kmeans_test$cluster_k8<-factor(kmeans_test$cluster_k8,levels=c('1','2','3','4','5','6','7','8'))


#kmeans_test$Gene_ID<-rownames(kmeans_test)
kmeans_merge<-merge(kmeans_test,halflife_TE[,c(1,13)],by='Gene_ID')
kmeans_merge<-kmeans_merge[order(kmeans_merge$cluster_k4),]

f1<-ggplot(kmeans_merge,aes(x=cluster_k4,y=log2FC_halflife,fill=cluster_k4))+geom_boxplot()+theme_classic()+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
f1
ggsave(plot=f1,filename='log2hf_clusters_k4.png',height=4,width=4.5)

f1<-ggplot(kmeans_merge,aes(x=cluster_k4,y=siControl_half_life,fill=cluster_k4))+geom_boxplot()+theme_classic()+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
f1
ggsave(plot=f1,filename='controlhf_clusters_k4.png',height=4,width=4.5)

f1<-ggplot(kmeans_merge,aes(x=cluster_k4,y=siCNOT1_half_life,fill=cluster_k4))+geom_boxplot()+theme_classic()+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
f1
ggsave(plot=f1,filename='nothf_clusters_k4.png',height=4,width=4.5)


f1<-ggplot(kmeans_merge,aes(x=cluster_k5,y=log2FC_halflife,fill=cluster_k5))+geom_boxplot()+theme_classic()+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
f1
ggsave(plot=f1,filename='log2hf_clusters_k5.png',height=4,width=4.5)

f1<-ggplot(kmeans_merge,aes(x=cluster_k5,y=siControl_half_life,fill=cluster_k5))+geom_boxplot()+theme_classic()+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
f1
ggsave(plot=f1,filename='controlhf_clusters_k5.png',height=4,width=4.5)

f1<-ggplot(kmeans_merge,aes(x=cluster_k5,y=siCNOT1_half_life,fill=cluster_k5))+geom_boxplot()+theme_classic()+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
f1
ggsave(plot=f1,filename='nothf_clusters_k5.png',height=4,width=4.5)


f1<-ggplot(kmeans_merge,aes(x=cluster_k6,y=log2FC_halflife,fill=cluster_k6))+geom_boxplot()+theme_classic()+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
f1
ggsave(plot=f1,filename='log2hf_clusters_k6.png',height=4,width=4.5)

f1<-ggplot(kmeans_merge,aes(x=cluster_k6,y=siControl_half_life,fill=cluster_k6))+geom_boxplot()+theme_classic()+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
f1
ggsave(plot=f1,filename='controlhf_clusters_k6.png',height=4,width=4.5)

f1<-ggplot(kmeans_merge,aes(x=cluster_k6,y=siCNOT1_half_life,fill=cluster_k6))+geom_boxplot()+theme_classic()+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
f1
ggsave(plot=f1,filename='nothf_clusters_k6.png',height=4,width=4.5)



clusters<-c('1','2','3','4','5','6')

for(i in 1:6){
  print(i)
  print(mean(subset(kmeans_merge,cluster_k6==clusters[i])$log2FC_halflife))
  print(mean(subset(kmeans_merge,cluster_k6==clusters[i])$siControl_half_life))
  print(mean(subset(kmeans_merge,cluster_k6==clusters[i])$siCNOT1_half_life))
}


###
kmeans_test$Gene_ID<-rownames(kmeans_test)
kmeans_merge<-merge(kmeans_test,halflife_TE[,c(1,14)],by='Gene_ID')
kmeans_merge<-kmeans_merge[order(kmeans_merge$cluster_k8),]

f1<-ggplot(kmeans_merge,aes(x=cluster_k8,y=log2FC_halflife,fill=cluster_k8))+geom_boxplot()+theme_classic()+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
f1
ggsave(plot=f1,filename='log2hf_clusters_k8.png',height=4,width=4.5)

f1<-ggplot(kmeans_merge,aes(x=cluster_k8,y=siControl_half_life,fill=cluster_k8))+geom_boxplot()+theme_classic()+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
f1
ggsave(plot=f1,filename='conhf_clusters_k8.png',height=4,width=4.5)


f1<-ggplot(kmeans_merge,aes(x=cluster_k8,y=siCNOT1_half_life,fill=cluster_k8))+geom_boxplot()+theme_classic()+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
f1
ggsave(plot=f1,filename='nothf_clusters_k8.png',height=4,width=4.5)


clusters<-c('1','2','3','4','5','6','7','8')

for(i in 1:8){
  print(i)
  print(mean(subset(kmeans_merge,cluster_k8==clusters[i])$log2FC_halflife))
  print(mean(subset(kmeans_merge,cluster_k8==clusters[i])$siControl_half_life))
  print(mean(subset(kmeans_merge,cluster_k8==clusters[i])$siCNOT1_half_life))
}




###
kmeans_test$Gene_ID<-rownames(kmeans_test)
kmeans_merge<-merge(kmeans_test,halflife_TE[,c(1,14)],by='Gene_ID')
kmeans_merge<-kmeans_merge[order(kmeans_merge$cluster_k7),]

f1<-ggplot(kmeans_merge,aes(x=cluster_k7,y=log2FC_halflife,fill=cluster_k7))+geom_boxplot()+theme_classic()+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
f1
ggsave(plot=f1,filename='log2hf_clusters_k7.png',height=4,width=4.5)

f1<-ggplot(kmeans_merge,aes(x=cluster_k7,y=siControl_half_life,fill=cluster_k7))+geom_boxplot()+theme_classic()+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
f1
ggsave(plot=f1,filename='conhf_clusters_k7.png',height=4,width=4.5)


f1<-ggplot(kmeans_merge,aes(x=cluster_k7,y=siCNOT1_half_life,fill=cluster_k7))+geom_boxplot()+theme_classic()+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.title=element_text(size=16),legend.title=element_blank())
f1
ggsave(plot=f1,filename='nothf_clusters_k7.png',height=4,width=4.5)


clusters<-c('1','2','3','4','5','6','7','8')

for(i in 1:7){
  print(i)
  print(mean(subset(kmeans_merge,cluster_k7==clusters[i])$log2FC_halflife))
  print(mean(subset(kmeans_merge,cluster_k7==clusters[i])$siControl_half_life))
  print(mean(subset(kmeans_merge,cluster_k7==clusters[i])$siCNOT1_half_life))
}



##
clusters<-c('1','2','3','4','5','6','7','8')

for(i in 1:5){
  print(i)
  print(mean(subset(kmeans_merge,cluster_k5==clusters[i])$log2FC_halflife))
  print(mean(subset(kmeans_merge,cluster_k5==clusters[i])$siControl_half_life))
  print(mean(subset(kmeans_merge,cluster_k5==clusters[i])$siCNOT1_half_life))
}


##
clusters<-c('1','2','3','4','5','6','7','8')

for(i in 1:4){
  print(i)
  print(mean(subset(kmeans_merge,cluster_k4==clusters[i])$log2FC_halflife))
  print(mean(subset(kmeans_merge,cluster_k4==clusters[i])$siControl_half_life))
  print(mean(subset(kmeans_merge,cluster_k4==clusters[i])$siCNOT1_half_life))
}

for(i in 1:4){
  print(i)
  print(median(subset(kmeans_merge,cluster_k4==clusters[i])$log2FC_halflife))
  print(median(subset(kmeans_merge,cluster_k4==clusters[i])$siControl_half_life))
  print(median(subset(kmeans_merge,cluster_k4==clusters[i])$siCNOT1_half_life))
}

for(i in 1:4){
  print(i)
  print(median(subset(kmeans_merge,cluster_k4_2==clusters[i])$log2FC_halflife))
  print(median(subset(kmeans_merge,cluster_k4_2==clusters[i])$siControl_half_life))
  print(median(subset(kmeans_merge,cluster_k4_2==clusters[i])$siCNOT1_half_life))
}

##
kmeans_test$Gene_ID<-rownames(kmeans_test)
kmeans_merge<-merge(kmeans_test,halflife_TE[,c(1,14)],by='Gene_ID')
kmeans_merge<-kmeans_merge[order(kmeans_merge$cluster_k4),]

kmeans_merge$cluster_k4_2<-as.character(kmeans_merge$cluster_k4_2)
kmeans_merge$cluster_k4<-as.character(kmeans_merge$cluster_k4)

ggplot(kmeans_merge,aes(x=siControl_half_life,y=siCNOT1_half_life,col=cluster_k4_2))+geom_point()+theme_classic()
ggplot(kmeans_merge,aes(x=siControl_half_life,y=siCNOT1_half_life,col=cluster_k8))+geom_point()+theme_classic()


heatmap.2(as.matrix(kmeans_merge[,c(2,3)]),RowSideColours=kmeans_merge$cluster_k4)

library(RColorBrewer)
library(gplots)

coul <- colorRampPalette(brewer.pal(8, "YlOrRd"))(50)

png(file = "halflifeclustersX.png",width=350,height=850)
heatmap.2(as.matrix(kmeans_merge[,c(2,3)]),scale='col',col = coul,Rowv=FALSE,Colv=TRUE,dendrogram='none',RowSideColors = as.character(kmeans_merge$cluster_k4),trace="none",density.info='none', hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(20,12))
dev.off()


coul <- colorRampPalette(brewer.pal(8, "PuRd"))(50)

png(file = "halflifeclusterY.png",width=350,height=850)
heatmap.2(as.matrix(kmeans_merge[,c(2,3,9)]),scale='col',col = coul,Rowv=FALSE,Colv=TRUE,dendrogram='none',RowSideColors = as.character(kmeans_merge$cluster_k4),trace="none",density.info='none', hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(20,12))
dev.off()

png(file = "halflifeclusterY.png",width=350,height=850)
heatmap.2(as.matrix(kmeans_merge[1:5000,c(2,3,9)]),scale='col',col = coul,Rowv=FALSE,Colv=TRUE,dendrogram='none',RowSideColors = as.character(kmeans_merge[1:5000,'cluster_k4']),trace="none",density.info='none', hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(20,12),labRow=FALSE)
dev.off()


png(file = "halflifeclusterY.png",width=350,height=850)
heatmap.2(as.matrix(kmeans_merge[1:5000,c(2,3,9)]),scale='col',col = coul,Rowv=FALSE,Colv=FALSE,dendrogram='none',RowSideColors = as.character(kmeans_merge[1:5000,'cluster_k4']),trace="none",density.info='none',labRow=FALSE)
dev.off()


library(heatmap3)
plotdata<-kmeans_merge[,c(2,3,9,4)]

png(file = "halflifeclusterZ.png",width=350,height=850)
heatmap3(as.matrix(kmeans_merge[,c(2,3,9)]),scale='col',hclustfun = hclust(x,method = 'centroid'),col = coul,Rowv=FALSE,Colv=FALSE,showColDendro=F, showRowDendro = F,RowSideColors = as.character(kmeans_merge[,'cluster_k4']),labRow=FALSE)
dev.off()

coul <- colorRampPalette(brewer.pal(8, "Spectral"))(50)
png(file = "halflifeclusterZ.png",width=350,height=850)
heatmap3(as.matrix(plotdata[,c(1,2,3)]),scale='none',col = coul,Rowv=FALSE,Colv=FALSE,showColDendro=F, showRowDendro = F, distfun = function(x) dist(x,method = 'euclidean'),RowSideColors = as.character(plotdata[,4]),labRow=FALSE)
dev.off()


png(file = "halflifeclusterZ2.png",width=350,height=850)
heatmap3(as.matrix(kmeans_test2[,c(1,2)]),scale='none',col = coul,Rowv=FALSE,Colv=FALSE,showColDendro=F, showRowDendro = F, distfun = function(x) dist(x,method = 'euclidean'),RowSideColors = as.character(kmeans_test2[,3]),labRow=FALSE)
dev.off()

kmeans_test2<-kmeans_test2[order(kmeans_test2$cluster_k4),]
png(file = "halflifeclusterZ3.png",width=350,height=850)
heatmap3(as.matrix(kmeans_test2[,c(1,2)]),scale='col',col = coul,Rowv=FALSE,Colv=FALSE,showColDendro=F, showRowDendro = F, distfun = function(x) dist(x,method = 'euclidean'),RowSideColors = as.character(kmeans_test2[,3]),labRow=FALSE)
dev.off()


png(file = "halflifeclusterZ3.png",width=350,height=850)
heatmap3(as.matrix(kmeans_test2[,c(1,2)]),scale='col',col = coul,Rowv=FALSE,Colv=FALSE,showColDendro=F, showRowDendro = F,RowSideColors = as.character(kmeans_test2[,3]),labRow=FALSE)
dev.off()


