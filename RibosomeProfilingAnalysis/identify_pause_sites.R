################ functions ##################
#get Transcripts per million for Total RNA
calc_TPM<-function(Total,HDF){
  print(nrow(Total))
  Total<-merge(Total,HDF[,c(5,3,2,7:12)],by.x='gene_name',by.y='Gene_ID')
  print(nrow(Total))
  
  #calculate TPM
  Total$length_kb<-nchar(Total$sequence)/1000
  Total$RPK<-Total[,3]/Total$length_kb
  Total$TPM<-Total$RPK/(sum(Total$RPK)/1000000)
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
}

lapply(Total_gene_counts,function(x) nrow(x))
lapply(Total_gene_counts,function(x) sum(x$TPM))


########## RPF data - normalise for library size ############
sample_names<-c('siControl_rep1','siControl_rep2','siControl_rep3','siCNOT1_rep1','siCNOT1_rep2','siCNOT1_rep3')
RPF_counts<-NULL
for(i in 1:length(sample_names)){
  RPF_counts[[i]]<-read.table(paste0(countsdir,'RPFcountsTable_sorted_RPFaligned_',sample_names[i],'.txt'),stringsAsFactors = FALSE,header=TRUE,sep='\t')
  
  print(sum(as.numeric(RPF_counts[[i]]$Total_RPF_counts)))
  
  RPF_sum<-sum(as.numeric(RPF_counts[[i]]$Total_RPF_counts))
  print(RPF_sum)
  RPF_counts[[i]]$Total_RPF_counts<-(RPF_counts[[i]]$Total_RPF_counts/RPF_sum)*1000000
  SF<-1000000/RPF_sum
  print(SF)
  for(j in 1:nrow(RPF_counts[[i]])){
    datax<-as.numeric(strsplit(as.character(RPF_counts[[i]][j,3]),split=' ')[[1]])
    RPF_counts[[i]][j,3]<-paste(c(datax*SF),sep=',',collapse=',')
  }
  print(nrow(RPF_counts[[i]]))
  RPF_counts[[i]]<-subset(RPF_counts[[i]], Total_RPF_counts>10)
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

########## RPF data - normalise for mRNA abundance (TPM from Total RNA-seq) ############
sample_names<-c('siControl_rep1','siControl_rep2','siControl_rep3','siCNOT1_rep1','siCNOT1_rep2','siCNOT1_rep3')
RPF_Total<-NULL
for(i in 1:length(sample_names)){
  RPF_Total[[i]]<-merge(RPF_counts[[i]][,c(7,5,3,4)],Total_gene_counts[[i]][,c(1,6:11,14)],by.x='Gene_ID',by.y='gene_name')
}

for(i in 1:length(sample_names)){
  for(j in 1:nrow(RPF_Total[[i]])){
    datay<-as.numeric(strsplit(as.character(RPF_Total[[i]][j,'pos_RPF_counts']),split=',')[[1]])
    datay<-datay/RPF_Total[[i]][j,'TPM']
    RPF_Total[[i]][j,3]<-paste(datay,sep=',',collapse=',')
    
  }
  print(nrow(RPF_Total[[i]]))
}


##### sort reads and combine replicates #####
for(i in 1:length(sample_names)){
  RPF_Total[[i]]$CDS_section_InF_pos_counts<-rep('',nrow(RPF_Total[[i]]))
  for(j in 1:nrow(RPF_Total[[i]])){
    CDSsection_p1<-as.numeric(strsplit(RPF_Total[[i]][j,'pos_RPF_counts'],split=',')[[1]]) 
    UTR5pos<-as.numeric(RPF_Total[[i]][j,'UTR5pos'])
    CDSpos<-as.numeric(RPF_Total[[i]][j,'CDSpos'])
    CDSsection_p2<-CDSsection_p1[(UTR5pos+34):(CDSpos-27)]
    CDSsection_p2_InF<-CDSsection_p2[seq(1,length(CDSsection_p2),3)]
    RPF_Total[[i]][j,'CDS_section_InF_pos_counts']<-paste0(CDSsection_p2_InF,collapse=',')
  }
}


control_data<-RPF_Total[1:3]
geneids<-control_data[[1]]$Gene_ID

control_combined<-data.frame(matrix(nrow=length(geneids),ncol=2,0),stringsAsFactors = FALSE)
names(control_combined)<-c('Gene_ID','CDS_counts')

for(i in 1:length(geneids)){
  control_combined[i,'Gene_ID']<-geneids[i]
  
  r1<-as.numeric(strsplit(subset(control_data[[1]], Gene_ID==geneids[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r2<-as.numeric(strsplit(subset(control_data[[2]], Gene_ID==geneids[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r3<-as.numeric(strsplit(subset(control_data[[3]], Gene_ID==geneids[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  
  control_combined[i,'CDS_counts']<-paste(colMeans(rbind(r1,r2,r3)),collapse=',')
  rm(r1)
  rm(r2)
  rm(r3)
  
}

not_data<-RPF_Total[4:6]
geneids<-not_data[[1]]$Gene_ID

not_combined<-data.frame(matrix(nrow=length(geneids),ncol=2,0),stringsAsFactors = FALSE)
names(not_combined)<-c('Gene_ID','CDS_counts')

for(i in 1:length(geneids)){
  not_combined[i,'Gene_ID']<-geneids[i]
  
  r1<-as.numeric(strsplit(subset(not_data[[1]], Gene_ID==geneids[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r2<-as.numeric(strsplit(subset(not_data[[2]], Gene_ID==geneids[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r3<-as.numeric(strsplit(subset(not_data[[3]], Gene_ID==geneids[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  
  not_combined[i,'CDS_counts']<-paste(colMeans(rbind(r1,r2,r3)),collapse=',')
  rm(r1)
  rm(r2)
  rm(r3)
  
}

##### get delta per biological replicate and then average #####
geneids<-unique(control_combined$Gene_ID)
delta_combined<-data.frame(matrix(nrow=length(geneids),ncol=2,0),stringsAsFactors = FALSE)
names(delta_combined)<-c('Gene_ID','delta_CDS_counts')

for(i in 1:length(geneids)){
  delta_combined[i,'Gene_ID']<-geneids[i]
  
  r1C<-as.numeric(strsplit(subset(control_data[[1]], Gene_ID==geneids[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r1N<-as.numeric(strsplit(subset(not_data[[1]], Gene_ID==geneids[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r1D<-r1N-r1C

  
  r2C<-as.numeric(strsplit(subset(control_data[[2]], Gene_ID==geneids[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r2N<-as.numeric(strsplit(subset(not_data[[2]], Gene_ID==geneids[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r2D<-r2N-r2C

  
  r3C<-as.numeric(strsplit(subset(control_data[[3]], Gene_ID==geneids[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r3N<-as.numeric(strsplit(subset(not_data[[3]], Gene_ID==geneids[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r3D<-r3N-r3C

  
  delta_combined[i,'delta_CDS_counts']<-paste(colMeans(rbind(r1D,r2D,r3D)),collapse=',')
  rm(r1C,r1N,r1D)
  rm(r2C,r2N,r2D)
  rm(r3C,r3N,r3D)
}



################## identify peaks in control condition, siCNOT1 condition and if the peak is resolved/induced/sustained after CNOT1 depletion ###########        

names(control_combined)[2]<-'con_CDS_counts'
names(not_combined)[2]<-'not_CDS_counts'

combined_data<-merge(control_combined,not_combined,by='Gene_ID')
combined_data<-merge(combined_data,delta_combined[,c(1,2)],by='Gene_ID')

combined_data$Cpauseincrease_mean<-numeric(nrow(combined_data))
combined_data$Cpauseincrease_mean_pos<-numeric(nrow(combined_data))

for(i in 1:nrow(combined_data)){
  d<-as.numeric(strsplit(combined_data[i,'con_CDS_counts'],split=',')[[1]])
  dup<-d[d>0]
  
  dup_mean<-mean(dup)
  combined_data[i,'Cpauseincrease_mean']<-length(dup[which(dup>(10*dup_mean))])
  combined_data[i,'Cpauseincrease_mean_pos']<-paste0(which(d>(10*dup_mean)),collapse=',')
  
  rm(d,dup,dup_mean)
}

combined_data$Npauseincrease_mean<-numeric(nrow(combined_data))
combined_data$Npauseincrease_mean_pos<-numeric(nrow(combined_data))

for(i in 1:nrow(combined_data)){
  d<-as.numeric(strsplit(combined_data[i,'not_CDS_counts'],split=',')[[1]])
  dup<-d[d>0]
 
  dup_mean<-mean(dup)
  combined_data[i,'Npauseincrease_mean']<-length(dup[which(dup>(10*dup_mean))])
  combined_data[i,'Npauseincrease_mean_pos']<-paste0(which(d>(10*dup_mean)),collapse=',')
  
  rm(d,dup,dup_mean)
}

combined_data$Dpauseincrease_mean<-numeric(nrow(combined_data))
combined_data$Dpausedecrease_mean<-numeric(nrow(combined_data))
combined_data$Dpauseincrease_mean_pos<-numeric(nrow(combined_data))
combined_data$Dpausedecrease_mean_pos<-numeric(nrow(combined_data))

for(i in 1:nrow(combined_data)){
  d<-as.numeric(strsplit(combined_data[i,'delta_CDS_counts'],split=',')[[1]])
  dup<-d[d>0]
  dup_mean<-mean(dup)
  combined_data[i,'Dpauseincrease_mean']<-length(dup[which(dup>(10*dup_mean))])
  combined_data[i,'Dpauseincrease_mean_pos']<-paste0(which(d>(10*dup_mean)),collapse=',')
  
  ddown<-d[d<0]
  ddown_mean<-mean(ddown)
  combined_data[i,'Dpausedecrease_mean']<-length(ddown[which(ddown<(10*ddown_mean))])
  combined_data[i,'Dpausedecrease_mean_pos']<-paste0(which(d<(10*ddown_mean)),collapse=',')
  
  rm(d,dup,ddown,dup_mean,ddown_mean)
}

       
########## categorise pause types ############
       
pauseinfo<-combined_data[,c(1,7,10,16,15)]

sort_line<-function(listx){
  linesx<-NULL
  vecs<-unique(as.numeric(c(strsplit(listx[[2]],split=',')[[1]],strsplit(listx[[3]],split=',')[[1]],strsplit(listx[[4]],split=',')[[1]],strsplit(listx[[5]],split=',')[[1]])))
  if(length(vecs)>0){
    for(i in 1:length(vecs)){  
      linex<-c(listx[[1]],vecs[i])
      if((vecs[i] %in% as.numeric(strsplit(listx[[2]],split=',')[[1]]))==TRUE){
        linex<-c(linex,'Y')
      }else{
        linex<-c(linex,'N')
      }
      
      if((vecs[i] %in% as.numeric(strsplit(listx[[3]],split=',')[[1]]))==TRUE){
        linex<-c(linex,'Y')
      }else{
        linex<-c(linex,'N')
      }
      
      if((vecs[i] %in% as.numeric(strsplit(listx[[4]],split=',')[[1]]))==TRUE){
        linex<-c(linex,'Y')
      }else{
        linex<-c(linex,'N')
      }
      
      if((vecs[i] %in% as.numeric(strsplit(listx[[5]],split=',')[[1]]))==TRUE){
        linex<-c(linex,'Y')
      }else{
        linex<-c(linex,'N')
      }
      linesx<-rbind(linesx,linex)
      rm(linex)
    }
    return(linesx)
  }
  
}


pauses<-NULL
for(i in 1:nrow(pauseinfo)){
  pauses<-rbind(pauses,sort_line(pauseinfo[i,]))
}
pauses<-data.frame(pauses,stringsAsFactors = FALSE)

names(pauses)<-c('gene_name','position','con_peak','not_peak','delta_decrease','delta_increase')
pauses$position<-as.numeric(pauses$position)

sustainedpause<-subset(pauses,con_peak=='Y' & not_peak=='Y' & delta_decrease=='N' & delta_increase=='N')
inducedpause<-subset(pauses,con_peak=='N' & not_peak=='Y' & delta_decrease=='N' & delta_increase=='Y')
resolvedpause<-subset(pauses,con_peak=='Y' & not_peak=='N' & delta_decrease=='Y' & delta_increase=='N')

resolvedonly<-subset(resolvedpause,(gene_name %in% inducedpause$gene_name)==FALSE)
resolvedonly<-subset(resolvedonly,(gene_name %in% sustainedpause$gene_name)==FALSE)

inducedonly<-subset(inducedpause,(gene_name %in% resolvedpause$gene_name)==FALSE)
inducedonly<-subset(inducedonly,(gene_name %in% sustainedpause$gene_name)==FALSE)

sustainedonly<-subset(sustainedpause,(gene_name %in% resolvedpause$gene_name)==FALSE)
sustainedonly<-subset(sustainedonly,(gene_name %in% inducedpause$gene_name)==FALSE)











