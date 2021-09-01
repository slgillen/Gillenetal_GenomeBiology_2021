#delta from average of TPM normalsied replicates, by group - outlier removal
library(seqinr)
library(Biostrings)
library(zoo)
library(ggplot2)
library(dplyr)
library(gplots)
library(colorRamps)

#humanDF sequences - most abundant transcript per gene used as the sequence

################ functions ##################
#Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#Count up all the RPK values in a sample and divide this number by 1,000,000. This is your per million scaling factor.
#Divide the RPK values by the per million scaling factor. This gives you TPM.
calc_TPM<-function(Total,HDF){
  print(nrow(Total))
  Total<-merge(Total,HDF[,c(5,3,2,7:12)],by.x='gene_name',by.y='Gene_ID')
  print(nrow(Total))
  
  #calculate TPM
  Total$length_kb<-nchar(Total$sequence)/1000
  Total$RPK<-Total[,3]/Total$length_kb
  Total$TPM<-Total$RPK/(sum(Total$RPK)/1000000)
  Total$TPM<-Total$TPM*20
  return(Total)
}

searchFeatures<-function(x,feature){
  section<-strsplit(x,split='|',fixed=TRUE)[[1]][[grep(feature,strsplit(x,split='|',fixed=TRUE)[[1]])]]
  return(section)
}

########## Total RNA sample - sort TPM ############
sample_names<-c('siControl_rep1','siControl_rep2','siControl_rep3','siCNOT1_rep1','siCNOT1_rep2','siCNOT1_rep3')
Total_gene_counts<-NULL
for(i in 1:length(sample_names)){
  Total_gene_counts[[i]]<-read.table(paste0(countsdir,'geneID_reverse_counts_TotalRNA_',sample_names[i],'.txt'),stringsAsFactors = FALSE,header=FALSE)
  names(Total_gene_counts[[i]])<-c('ENSG','gene_name','Total_read_counts')
  print(sum(Total_gene_counts[[i]][,3]))
  Total_gene_counts[[i]]<-calc_TPM(Total_gene_counts[[i]],humanDF)
  print(nrow(Total_gene_counts[[i]]))
  Total_gene_counts[[i]]<-subset(Total_gene_counts[[i]],TPM>0)
  print(nrow(Total_gene_counts[[i]]))
  Total_gene_counts[[i]]<-subset(Total_gene_counts[[i]],TPM>1)
  print(nrow(Total_gene_counts[[i]]))
}


lapply(Total_gene_counts,function(x) nrow(x))
lapply(Total_gene_counts,function(x) length(unique(x$gene_name)))
lapply(Total_gene_counts,function(x) sum(x$TPM))

########## RPF data - normalise for library size ############
sample_names<-c('siControl_rep1','siControl_rep2','siControl_rep3','siCNOT1_rep1','siCNOT1_rep2','siCNOT1_rep3')
RPF_counts<-NULL
for(i in 1:length(sample_names)){
  RPF_counts[[i]]<-read.table(paste0(countsdir,'RPFcountsTable_RPF_',sample_names[i],'.txt'),stringsAsFactors = FALSE,header=TRUE,sep='\t')
  RPF_sum<-sum(as.numeric(RPF_counts[[i]]$Total_RPF_counts))
  print(RPF_sum)
  RPF_counts[[i]]$Total_RPF_counts<-(RPF_counts[[i]]$Total_RPF_counts/RPF_sum)*20000000
  SF<-20000000/RPF_sum
  print(SF)
  for(j in 1:nrow(RPF_counts[[i]])){
    datax<-as.numeric(strsplit(as.character(RPF_counts[[i]][j,3]),split=' ')[[1]])
    RPF_counts[[i]][j,3]<-paste(c(datax*SF),sep=',',collapse=',')
  }
  print(nrow(RPF_counts[[i]]))
  RPF_counts[[i]]<-subset(RPF_counts[[i]], Total_RPF_counts>20)
  print(nrow(RPF_counts[[i]]))
}

lapply(RPF_counts, function(x) nrow(x))
lapply(RPF_counts, function(x) sum(x$Total_RPF_counts))


for(k in 1:length(RPF_counts)){
  s<-RPF_counts[[k]]
  s$ENST<-sapply(strsplit(as.character(s[,'Transcript_ID']), "|",fixed=TRUE), "[[", 1)
  s$ENSG<-sapply(strsplit(as.character(s[,'Transcript_ID']), "|",fixed=TRUE), "[[", 2)
  s$Gene_ID<-sapply(strsplit(as.character(s[,'Transcript_ID']), "|",fixed=TRUE), "[[", 6)
  s$length<-sapply(strsplit(as.character(s[,'Transcript_ID']), "|",fixed=TRUE), "[[", 7)
  
  
  s$UTR5pos<-rep(NA,nrow(s))
  s$CDSpos<-rep(NA,nrow(s))
  s$UTR3pos<-rep(NA,nrow(s))
  for(i in 1:nrow(s)){
    info<-s[i,1]
    
    if(grepl('CDS:',info)==TRUE){
      CDSp<-searchFeatures(info,'CDS:') 
      s[i,'CDSpos']<-strsplit(CDSp,split='-',fixed=TRUE)[[1]][2]
    }
    if(grepl('UTR5',info)==TRUE){
      UTR5p<-searchFeatures(info,'UTR5')
      s[i,'UTR5pos']<-strsplit(UTR5p,split='-',fixed=TRUE)[[1]][2]
    }
    if(grepl('UTR3',info)==TRUE){
      UTR3p<-searchFeatures(info,'UTR3')
      s[i,'UTR3pos']<-strsplit(UTR3p,split='-',fixed=TRUE)[[1]][2]
    }
  }
  RPF_counts[[k]]<-s
}
lapply(RPF_counts, function(x) nrow(x))

########## RPF data - normalise for TPM ############
sample_names<-c('siControl_rep1','siControl_rep2','siControl_rep3','siCNOT1_rep1','siCNOT1_rep2','siCNOT1_rep3')
RPF_Total<-NULL
for(i in 1:length(sample_names)){
  #print(nrow(Total_gene_counts[[i]]))
  #Total_gene_counts[[i]]<-Total_gene_counts[[i]][!duplicated(Total_gene_counts[[i]]$gene_name),]
  print(nrow(Total_gene_counts[[i]]))
  print(nrow(RPF_counts[[i]]))
  RPF_Total[[i]]<-merge(RPF_counts[[i]][,c(7,5,3,4)],Total_gene_counts[[i]][,c(1,6:11,14)],by.x='Gene_ID',by.y='gene_name')
  #print(nrow(RPF_Total[[i]]))
  #RPF_Total[[i]]<-RPF_Total[[i]][!duplicated(RPF_Total[[i]]$Gene_ID),]
  #print(nrow(RPF_Total[[i]]))
}

lapply(RPF_Total, function(x) nrow(x))
lapply(RPF_Total, function(x) length(unique(x$Gene_ID)))


#TPM normalise
for(i in 1:length(sample_names)){
  for(j in 1:nrow(RPF_Total[[i]])){
    datay<-as.numeric(strsplit(as.character(RPF_Total[[i]][j,'pos_RPF_counts']),split=',')[[1]])
    RPF_Total[[i]][j,3]<-paste(c(datay/RPF_Total[[i]][j,'TPM']),sep=',',collapse=',')
  }
}

lapply(RPF_Total, function(x) nrow(x))
lapply(RPF_Total, function(x) length(unique(x$Gene_ID)))

########## keep only genes detected in all samples ############
keep_genes<-RPF_Total[[1]]$Gene_ID
length(keep_genes)
for(i in 2:length(sample_names)){
  keep_genes<-keep_genes[(keep_genes %in% RPF_Total[[i]]$Gene_ID)==TRUE]
  print(length(keep_genes))
}


########## get in-frame reads + whole CDS ############
#consider looking at the 5'UTR also
for(i in 1:length(sample_names)){
  RPF_Total[[i]]$CDS_section_pos_counts<-rep('',nrow(RPF_Total[[i]]))
  RPF_Total[[i]]$CDS_section_InF_pos_counts<-rep('',nrow(RPF_Total[[i]]))
  for(j in 1:nrow(RPF_Total[[i]])){
    CDSsection_p1<-strsplit(RPF_Total[[i]][j,'pos_RPF_counts'],split=',')[[1]]
    UTR5pos<-as.numeric(RPF_Total[[i]][j,'UTR5pos'])
    CDSpos<-as.numeric(RPF_Total[[i]][j,'CDSpos'])
    CDSsection_p2<-CDSsection_p1[(UTR5pos-11):(CDSpos-11)]# account for offset
    RPF_Total[[i]][j,'CDS_section_pos_counts']<-paste0(CDSsection_p2,collapse=',')
    CDSsection_p2_InF<-CDSsection_p2[seq(1,length(CDSsection_p2),3)]
    RPF_Total[[i]][j,'CDS_section_InF_pos_counts']<-paste0(CDSsection_p2_InF,collapse=',')
  }
}
names(RPF_Total[[1]])
lapply(RPF_Total, function(x) nrow(x))
lapply(RPF_Total, function(x) length(unique(x$Gene_ID)))

lapply(RPF_Total,function(x) nrow(x))
RPF_Total<-lapply(RPF_Total,function(x) subset(x,nchar(x$CDS_sequence)>300))
lapply(RPF_Total,function(x) nrow(x))

###### conduct delta for each biological replicate ########
delta_replicates<-NULL

for(i in 1:3){
  genes<-unique(RPF_Total[[1]]$Gene_ID)
  delta_data<-data.frame(matrix(nrow=length(genes),ncol=2,0),stringsAsFactors = FALSE)
  names(delta_data)<-c('gene_name','CDS_delta')
  
  print(length(genes))
  for(j in 1:length(genes)){
    
    C<-subset(RPF_Total[[i]],Gene_ID==genes[j])$CDS_section_InF_pos_counts #can also just use InF
    N<-subset(RPF_Total[[(i+3)]],Gene_ID==genes[j])$CDS_section_InF_pos_counts # can also just use InF
    
    Cx<-as.numeric(strsplit(C,split=',')[[1]])
    Nx<-as.numeric(strsplit(N,split=',')[[1]])
    gene_delta<-Nx-Cx
    
    delta_data[j,'gene_name']<-genes[j]
    delta_data[j,'CDS_delta']<-paste(gene_delta,collapse=',')

  }
  delta_replicates[[i]]<-delta_data
  rm(delta_data)
  
}

lapply(delta_replicates, function(x) nrow(x))
lapply(delta_replicates, function(x) length(unique(x$gene_name)))


######### average from replicates #######
delta_average<-data.frame(matrix(nrow=nrow(delta_replicates[[1]]),ncol=2,0),stringsAsFactors = FALSE)
names(delta_average)<-c('gene_name','average_CDS_delta')
genes<-unique(delta_replicates[[1]]$gene_name)
for(i in 1:length(genes)){
  g1<-as.numeric(strsplit(subset(delta_replicates[[1]],gene_name==genes[i])$CDS_delta,split=',')[[1]])
  g2<-as.numeric(strsplit(subset(delta_replicates[[2]],gene_name==genes[i])$CDS_delta,split=',')[[1]])
  g3<-as.numeric(strsplit(subset(delta_replicates[[3]],gene_name==genes[i])$CDS_delta,split=',')[[1]])
  
  delta_average[i,'gene_name']<-genes[i]
  delta_average[i,'average_CDS_delta']<-paste(rowMeans(cbind(g1,g2,g3)),collapse=',')
}


###### plot delta data ########

delta_average_groups<-merge(delta_average,RPFdata[,c(1,6)],by.x='gene_name',by.y='Gene_ID')

window_sort<-function(num_vec,nwindows){
  num_vec<-as.numeric(num_vec)
  windowlen<-length(num_vec)/nwindows
  windowlen<-floor(windowlen)
  starts<-seq(1,length(num_vec),windowlen)
  starts<-starts[1:nwindows]
  ends<-c(seq(windowlen+1,length(num_vec),windowlen))
  ends<-ends[1:nwindows]
  ends[nwindows]<-length(num_vec)
  sum_windows<-apply(data.frame(starts,ends),1, function(x) (sum(num_vec[x[1]:x[2]],na.rm=TRUE)/length(num_vec[x[1]:x[2]])))
  return(sum_windows)
}

plot_all_metagene<-function(deltadata,sample_names,nwindows,plotname,colours,outliers){
  
  to_plot<-NULL
  to_plot2<-NULL
  for(k in 1:length(sample_names)){
    ss<-lapply(subset(deltadata,(Group==sample_names[k]))$average_CDS_delta, function(x) window_sort(as.numeric(strsplit(x,split=',')[[1]]),nwindows))
    CDSes<-do.call(rbind,ss)
    rownames(CDSes)<-subset(deltadata,(Group %in% sample_names[k])==TRUE)$gene_name
    print(nrow(CDSes))
    CDSes<-subset(CDSes, is.finite(CDSes[,1])==TRUE)
    print(nrow(CDSes))
    CDSes<-subset(CDSes,(rownames(CDSes) %in% outliers)==FALSE)
    print(nrow(CDSes))
    #return(CDSes)#for random peak check
    toplot<-data.frame(matrix(nrow=nwindows,ncol=3),stringsAsFactors = FALSE)
    names(toplot)<-c('window','delta','Group')
    toplot$window<-seq(1,nwindows,1)
    toplot$delta<-apply(CDSes,2,function(x) median(x))
    toplot$Group<-rep(sample_names2[k],nrow(toplot))
    to_plot<-rbind(to_plot,toplot)
    
    toplot2<-data.frame(matrix(nrow=nwindows-2,ncol=3),stringsAsFactors = FALSE)
    names(toplot2)<-c('window','delta','Group')
    toplot2$window<-seq(1,nwindows-2,1)
    toplot2$delta<-rollapply(toplot$delta,3,mean)
    toplot2$Group<-rep(sample_names[k],nrow(toplot2))
    to_plot2<-rbind(to_plot2,toplot2)
    
    rm(toplot)
    rm(toplot2)
    
  }
  
  to_plot$Group<-factor(to_plot$Group,levels=sample_names2)
  to_plot<-to_plot[order(to_plot$Group),]
  
  fx<-ggplot(to_plot, aes(x=window,y=delta,col=Group))+geom_line(size=2)+scale_colour_manual(values=c(colours))+theme_classic()+xlab('% CDS')+ylab('median normalised ribosome occupancy change (siCNOT1-siControl)')
  fx<-fx+theme_classic()+theme(legend.text=element_text(size=16),axis.text.x = element_text(size = 16),
                               axis.text.y = element_text(size = 16), axis.title = element_text(size = 16),
                               plot.title=element_text(size=16),legend.title=element_blank(),panel.grid = element_blank())+scale_x_continuous(labels=c(0,25,50,75,100),breaks=c(0,10,20,30,40))
  
  ggsave(plot=fx,filename=paste0(plotname,'_CDS.png'),width=12,height=7)
  
  
  to_plot2$Group<-factor(to_plot2$Group,levels=sample_names2)
  to_plot2<-to_plot2[order(to_plot2$Group),]
  
  fx<-ggplot(to_plot2, aes(x=window,y=delta,col=Group))+geom_line(size=2)+scale_colour_manual(values=c(colours))+theme_classic()+xlab('% CDS')+ylab('median normalised ribosome occupancy change (siCNOT1-siControl)')
  fx<-fx+theme_classic()+theme(legend.text=element_text(size=16),axis.text.x = element_text(size = 16),
                               axis.text.y = element_text(size = 16), axis.title = element_text(size = 16),
                               plot.title=element_text(size=16),legend.title=element_blank(),panel.grid = element_blank())+scale_x_continuous(labels=c(0,25,50,75,100),breaks=c(0,9,18,28,38))
  
  ggsave(plot=fx,filename=paste0(plotname,'_CDS_rollapply.png'),width=8,height=3.5)
  
  
}




group_names<-c('TE down','no TE change','TE up')
group_names2<-c('TE down','no TE change','TE up')
group_colours<-c('dodgerblue3','goldenrod2','red3')
outliers_toremove<-c()


deltaplot<-plot_all_metagene(delta_average_groups,group_names[c(1,2,3,6,7)],group_names2[c(1,2,3,6,7)],40,'delta_TEgroups',group_colours,outliers_remove) 




