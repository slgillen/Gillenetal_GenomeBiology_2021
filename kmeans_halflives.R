library('factoextra')
library(ggplot2)
library(scales)

#k-means clustering of mRNA half-life changes with CNOT1 knockdown
halflife_data$siControl_half_life_log10<-log(halflife_data$siControl_half_life,10)
halflife_data$siCNOT1_half_life_log10<-log(halflife_data$siCNOT1_half_life,10)

halflife_df<-halflife_data[,c('siControl_halflife_log10','siCNOT1_halflife_log10')]
rownames(halflife_df)<-halflife_data$Gene_ID


#determine optimal number of clusters
fviz_nbclust(halflife_df, kmeans, method='silhouette')
png('kmeans_clusters_elbowmethod.png',width=400,height=250)

#plot data
res.k3<-hkmeans(halflife_df,3,hc.method='centroid')
png('kmeans_clusters_3.png',width=350,height=300)
fviz_cluster(res.k3, ellipse=FALSE,labelsize=6,ggtheme = theme_classic(),geom=c('point'),pointsize=1,stand=TRUE,shape='circle')+scale_colour_manual(values = c("tan2", "springgreen3","violetred3"))
dev.off()

kmeans_test$hk_clusters3<-res.hk3$cluster
kmeans_test$hk_clusters3<-factor(as.character(kmeans_test$hk_clusters3),levels=c('3','1','2'))
kmeans_test<-kmeans_test[order(kmeans_test$hk_clusters3),]
f1<-ggplot(kmeans_test,aes(x=hk_clusters3,y=siControl_half_life,fill=hk_clusters3))+geom_boxplot(outlier.size=0.5)+theme_classic()+scale_y_log10()+
  scale_fill_manual(values=c('violetred3','tan2','springgreen3'))+ylab('mRNA half-life (control)')+coord_trans(limy=c(0.5,150))+
  theme(axis.text.y=element_text(size=10),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_text(size=12),axis.title.x=element_blank())
f1
ggsave(f1, filename = 'k3_siControl_halflife_boxplot.png',width=3,height=3.5)

f1<-ggplot(kmeans_test,aes(x=hk_clusters3,y=siCNOT1_half_life,fill=hk_clusters3))+geom_boxplot(outlier.size=0.5)+theme_classic()+scale_y_log10()+
  scale_fill_manual(values=c('violetred3','tan2','springgreen3'))+ylab('mRNA half-life (siCNOT1)')+coord_trans(limy=c(0.5,150))+
  theme(axis.text.y=element_text(size=10),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_text(size=12),axis.title.x=element_blank())
f1
ggsave(f1, filename = 'k3_siCNOT1_halflife_boxplot.png',width=3,height=3.5)
