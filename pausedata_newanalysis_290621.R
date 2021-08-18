setwd("//john-doe/gw/Systems/Sarah/CNOTpaper/peakcalling")
countsdir<-'//john-doe/gw/Systems/Sarah/CNOTpaper/peakcalling/CountsFiles/'

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

#countsdir2<-'//john-doe/gw/Systems/Sarah/CNOTpaper/peakcalling/CountsFiles/tousefrom231220/'

########## RPF data - normalise for library size ############
sample_names<-c('siControl_rep1','siControl_rep2','siControl_rep3','siCNOT1_rep1','siCNOT1_rep2','siCNOT1_rep3')
#sample_names2<-c('con1','con2','con3','not1','not2','not3')
RPF_counts<-NULL
for(i in 1:length(sample_names)){
  #RPF_counts[[i]]<-read.table(paste0(countsdir2,'RPFcountsTable_sorted_RPFaligned_',sample_names2[i],'_27to31_oneTperG.txt'),stringsAsFactors = FALSE,header=TRUE,sep='\t')
  RPF_counts[[i]]<-read.table(paste0(countsdir,'RPFcountsTable_RPF_',sample_names[i],'_282930.txt'),stringsAsFactors = FALSE,header=TRUE,sep='\t')
  
  print(sum(as.numeric(RPF_counts[[i]]$Total_RPF_counts)))
  
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

########## keep only genes present in all samples ############
keep_genes<-RPF_Total[[i]]$Gene_ID
length(keep_genes)
for(i in 2:length(sample_names)){
  keep_genes<-keep_genes[(keep_genes %in% RPF_Total[[i]]$Gene_ID)==TRUE]
  print(length(keep_genes))
}
#also filter for RPF group list - i.e. no RPGs etc
length(keep_genes)
keep_genes<-keep_genes[(keep_genes %in% RPF27to31$Gene_ID)==TRUE] #changed from halflife_TE
length(keep_genes)

for(i in 1:length(sample_names)){
  print(nrow(RPF_Total[[i]]))
  RPF_Total[[i]]<-subset(RPF_Total[[i]],(Gene_ID %in% keep_genes)==TRUE)
  print(nrow(RPF_Total[[i]]))
}

for(i in 1:length(RPF_Total)){
  print(nrow(RPF_Total[[i]]))
  RPF_Total[[i]]<-subset(RPF_Total[[i]],startsWith(Gene_ID,'RP11-')==FALSE)
  print(nrow(RPF_Total[[i]]))
}

long_genes<-subset(RPF_Total[[2]], nchar(CDS_sequence)>8000)$Gene_ID

for(i in 1:length(RPF_Total)){
  print(nrow(RPF_Total[[i]]))
  RPF_Total[[i]]<-subset(RPF_Total[[i]],(Gene_ID %in% long_genes)==FALSE)
  print(nrow(RPF_Total[[i]]))
}


########## get in-frame reads of CDS ############
for(i in 1:length(sample_names)){
  RPF_Total[[i]]$CDS_section_pos_counts<-rep('',nrow(RPF_Total[[i]]))
  RPF_Total[[i]]$CDS_section_InF_pos_counts<-rep('',nrow(RPF_Total[[i]]))
  for(j in 1:nrow(RPF_Total[[i]])){
    CDSsection_p1<-as.numeric(strsplit(RPF_Total[[i]][j,'pos_RPF_counts'],split=',')[[1]]) #added as.numeric()
    UTR5pos<-as.numeric(RPF_Total[[i]][j,'UTR5pos'])
    CDSpos<-as.numeric(RPF_Total[[i]][j,'CDSpos'])
    #CDSsection_p2<-CDSsection_p1[(UTR5pos-11):(CDSpos-11)]# includes start and stop
    CDSsection_p2<-CDSsection_p1[(UTR5pos+34):(CDSpos-27)]
    RPF_Total[[i]][j,'CDS_section_pos_counts']<-paste0(CDSsection_p2,collapse=',')
    CDSsection_p2_InF<-CDSsection_p2[seq(1,length(CDSsection_p2),3)]
    RPF_Total[[i]][j,'CDS_section_InF_pos_counts']<-paste0(CDSsection_p2_InF,collapse=',')
  }
}

min100CDSlen<-subset(humanDF,nchar(CDS_sequence)>150)$Gene_ID

control_data<-RPF_Total[1:3]
unique_genes<-unique(control_data[[1]]$Gene_ID)
length(unique_genes)
unique_genes<-unique_genes[(unique_genes %in% control_data[[2]]$Gene_ID)==TRUE]
length(unique_genes)
unique_genes<-unique_genes[(unique_genes %in% control_data[[3]]$Gene_ID)==TRUE]
length(unique_genes)
unique_genes<-unique_genes[(unique_genes %in% min100CDSlen)==TRUE]
length(unique_genes)
control_combined<-data.frame(matrix(nrow=length(unique_genes),ncol=2,0),stringsAsFactors = FALSE)
names(control_combined)<-c('Gene_ID','CDS_counts')

for(i in 1:length(unique_genes)){
  control_combined[i,'Gene_ID']<-unique_genes[i]
  
  r1<-as.numeric(strsplit(subset(control_data[[1]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r2<-as.numeric(strsplit(subset(control_data[[2]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r3<-as.numeric(strsplit(subset(control_data[[3]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  
  control_combined[i,'CDS_counts']<-paste(colMeans(rbind(r1,r2,r3)),collapse=',')
  rm(r1)
  rm(r2)
  rm(r3)
  
}

not_data<-RPF_Total[4:6]
unique_genes<-unique(not_data[[1]]$Gene_ID)
length(unique_genes)
unique_genes<-unique_genes[(unique_genes %in% not_data[[2]]$Gene_ID)==TRUE]
length(unique_genes)
unique_genes<-unique_genes[(unique_genes %in% not_data[[3]]$Gene_ID)==TRUE]
length(unique_genes)
unique_genes<-unique_genes[(unique_genes %in% min100CDSlen)==TRUE]
length(unique_genes)
not_combined<-data.frame(matrix(nrow=length(unique_genes),ncol=2,0),stringsAsFactors = FALSE)
names(not_combined)<-c('Gene_ID','CDS_counts')

for(i in 1:length(unique_genes)){
  not_combined[i,'Gene_ID']<-unique_genes[i]
  
  r1<-as.numeric(strsplit(subset(not_data[[1]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r2<-as.numeric(strsplit(subset(not_data[[2]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r3<-as.numeric(strsplit(subset(not_data[[3]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  
  not_combined[i,'CDS_counts']<-paste(colMeans(rbind(r1,r2,r3)),collapse=',')
  rm(r1)
  rm(r2)
  rm(r3)
  
}

#and average of deltas that were done separately for each replicate
#try delta and ratio? #ratio would not working for TPM adjustment??

unique_genes<-unique(control_combined$Gene_ID)
length(unique_genes)
unique_genes<-unique_genes[(unique_genes %in% not_combined$Gene_ID)==TRUE]
length(unique_genes)
delta_combined<-data.frame(matrix(nrow=length(unique_genes),ncol=2,0),stringsAsFactors = FALSE)
names(delta_combined)<-c('Gene_ID','delta_CDS_counts')

for(i in 1:length(unique_genes)){
  delta_combined[i,'Gene_ID']<-unique_genes[i]
  
  r1C<-as.numeric(strsplit(subset(control_data[[1]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r1N<-as.numeric(strsplit(subset(not_data[[1]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r1D<-r1N-r1C

  
  r2C<-as.numeric(strsplit(subset(control_data[[2]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r2N<-as.numeric(strsplit(subset(not_data[[2]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r2D<-r2N-r2C

  
  r3C<-as.numeric(strsplit(subset(control_data[[3]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r3N<-as.numeric(strsplit(subset(not_data[[3]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r3D<-r3N-r3C

  
  delta_combined[i,'delta_CDS_counts']<-paste(colMeans(rbind(r1D,r2D,r3D)),collapse=',')
  rm(r1C,r1N,r1D)
  rm(r2C,r2N,r2D)
  rm(r3C,r3N,r3D)
}


###########################################################################################################
###########################################################################################################
###########################################################################################################


#get dataframe of positions of RPFs
get_plotdf_individual<-function(countsdf,readtype){
  #readcounts<-c(strsplit(countsdf[1,'UTR5_counts'],split=',')[[1]],strsplit(countsdf[1,'CDS_counts'],split=',')[[1]])
  readcounts<-strsplit(countsdf[1,readtype],split=',')[[1]]
  readcounts<-as.numeric(readcounts)
  
  output<-data.frame(matrix(ncol=2,nrow=length(readcounts),0))
  #count<-count
  output[,2]<-readcounts
  output[,1]<-seq(15,length(readcounts),1)
  
  names(output)<-c('position','value')
  output$position<-as.numeric(output$position)
  output$value<-as.numeric(output$value)
  return(output)
}


#need to add counts for 5'UTR, CDS - split plots?
setwd("//john-doe/gw/Systems/Sarah/CNOTpaper/cotranslationalassembly/individual_pauses")
setwd("//john-doe/gw/Systems/Sarah/CNOTpaper/pausesites_revisions")

genelist<-c('SLC7A5','NEK7')
colours<-c('grey58','violetred2')

genelist<-c('NAMPT')
colours<-c('grey68','darkmagenta')

genelist<-c('NEK7')
colours<-c('grey68','darkorange')

genelist<-c('USP13','PHB','RPTOR')
colours<-c('grey68','green4')

for(k in 1:length(genelist)){ #might be good to try with total mRNA TPM normalisation
  gene<-genelist[k]
  
  if((gene %in% delta_combined$Gene_ID)==TRUE){
    z<-subset(control_combined,Gene_ID==gene)
    print(nrow(z))
    RPFdata<-get_plotdf_individual(z,'con_CDS_counts') 
    RPFdata$group<-rep('mRNA',nrow(RPFdata))  
    
    
    gene<-genelist[k]
    z<-subset(not_combined,Gene_ID==gene)
    print(nrow(z))
    RPFdata2<-get_plotdf_individual(z,'not_CDS_counts') 
    RPFdata2$group<-rep('mRNA',nrow(RPFdata2))  
    
    RPFdata$condition<-rep('siControl',nrow(RPFdata))
    RPFdata2$condition<-rep('siCNOT1',nrow(RPFdata2))
    findmax<-max(RPFdata$value,RPFdata2$value)
    
    
    f1<-ggplot(data=RPFdata, aes(x=position,y=value,fill=group,col=group))+geom_bar(position = 'dodge2',stat='identity')+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+ggtitle(gene)+xlab('position in CDS (codons)')+ylab('ribosome occupancy')
    f1<-f1+scale_colour_manual(values=colours[1])+theme_bw()+theme(legend.text=element_text(size=10),axis.line = element_line(colour = "black"),axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8), axis.title = element_text(size = 8), panel.grid.major = element_blank(),panel.grid.minor=element_blank(),plot.title=element_text(size=10),legend.title=element_blank(),legend.position='none',panel.border=element_blank())+coord_trans(ylim=c(0,findmax))
    
    filenamex<-paste(gene,'_siControl_TPMnorm.png',sep='',collapse='')
    ggsave(filename=filenamex,plot=f1,width=3.4, height=1.4)
    
    
    ##################################
    
    
    f1<-ggplot(data=RPFdata2, aes(x=position,y=value,fill=group,col=group))+geom_bar(position = 'dodge2',stat='identity')+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+ggtitle(gene)+xlab('position in CDS (codons)')+ylab('ribosome occupancy')
    f1<-f1+scale_colour_manual(values=colours[2])+theme_bw()+theme(legend.text=element_text(size=10),axis.line = element_line(colour = "black"),axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8), axis.title = element_text(size = 8), panel.grid.major = element_blank(),panel.grid.minor=element_blank(),plot.title=element_text(size=10),legend.title=element_blank(),legend.position='none',panel.border=element_blank())+coord_trans(ylim=c(0,findmax))
    
    
    filenamex<-paste(gene,'_siCNOT1_TPMnorm.png',sep='',collapse='')
    ggsave(filename=filenamex,plot=f1,width=3.4, height=1.4)
    
    
    ##################################
    print('delta')
    z<-subset(delta_combined,Gene_ID==gene)
    print(nrow(z))
    RPFdata3<-get_plotdf_individual(z,'delta_CDS_counts') 
    RPFdata3$group<-rep('mRNA',nrow(RPFdata3))
    findmin<-min(RPFdata3$value)
    findmax<-max(RPFdata3$value)
    
    f1<-ggplot(data=RPFdata3, aes(x=position,y=value,fill=group,col=group))+geom_bar(position = 'dodge2',stat='identity')+scale_x_continuous(expand=c(0.0,0))+scale_y_continuous(expand=c(0,0))+ylab('delta')+ggtitle(gene)+xlab('position in CDS (codons)')+ylab('delta')
    f1<-f1+scale_colour_manual(values=colours[2])+theme_bw()+theme(legend.text=element_text(size=10),axis.line = element_line(colour = "black"),axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8), axis.title = element_text(size = 8), panel.grid.major = element_blank(),panel.grid.minor=element_blank(),plot.title=element_text(size=10),legend.title=element_blank(),legend.position='none',panel.border=element_blank())+coord_trans(ylim=c(findmin,findmax))
    
    
    filenamex<-paste(gene,'_delta_TPMnorm.png',sep='',collapse='')
    ggsave(filename=filenamex,plot=f1,width=3.4, height=1.4)

    ##################################
    #print('ratio')
    #z<-subset(delta_combined,Gene_ID==gene)
    #(nrow(z))
    #RPFdata4<-get_plotdf_individual(z,'ratio_CDS_counts') 
    #RPFdata4$group<-rep('mRNA',nrow(RPFdata4))
    #findmin<-min(RPFdata4$value)
    #findmax<-max(RPFdata4$value)
    
    #f1<-ggplot(data=RPFdata4, aes(x=position,y=value,col=group))+geom_bar(position = 'dodge2',stat='identity')+scale_x_continuous(expand=c(0.025,0))+ylab('ratio')+ggtitle(gene)+xlab('position in CDS (nt)')
    #<-f1+scale_colour_manual(values=colours[2])+theme_bw()+theme(legend.text=element_text(size=14),axis.line = element_line(colour = "black"),axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), panel.grid.major = element_blank(),panel.grid.minor=element_blank(),plot.title=element_text(size=16),legend.title=element_blank(),legend.position='none',panel.border=element_blank())+coord_trans(limy=c(findmin,findmax))
    
    #filenamex<-paste(gene,'_ratio_TPMnorm.png',sep='',collapse='')
    #ggsave(filename=filenamex,plot=f1,width=4, height=1.9)
    
    ##
    rm(RPFdata, RPFdata2, RPFdata3) #, RPFdata4
    
  }else{
    print(gene)
    print('data not available')
  }
}



######################################################################################################################
######################################################################################################################
######################################################################################################################

setwd("//john-doe/gw/Systems/Sarah/CNOTpaper/pausesites_revisions")

names(control_combined)[2]<-'con_CDS_counts'
names(not_combined)[2]<-'not_CDS_counts'

combined_data<-merge(control_combined,not_combined,by='Gene_ID')
combined_data<-merge(combined_data,delta_combined[,c(1,2)],by='Gene_ID')

combined_data$Cpauseincrease_med<-numeric(nrow(combined_data))
#combined_data$Cpausedecrease_med<-numeric(nrow(combined_data))
combined_data$Cpauseincrease_mean<-numeric(nrow(combined_data))
#combined_data$Cpausedecrease_mean<-numeric(nrow(combined_data))

for(i in 1:nrow(combined_data)){
  d<-as.numeric(strsplit(combined_data[i,'con_CDS_counts'],split=',')[[1]])
  dup<-d[d>0]
  
  dup_med<-median(dup)
  dup_mean<-mean(dup)
  combined_data[i,'Cpauseincrease_med']<-length(dup[which(dup>(10*dup_med))])
  combined_data[i,'Cpauseincrease_mean']<-length(dup[which(dup>(10*dup_mean))])
  
  #ddown<-d[d<0]
  #ddown_med<-median(ddown)
  #ddown_mean<-mean(ddown)
  #combined_data[i,'Cpausedecrease_med']<-length(ddown[which(ddown<(10*ddown_med))])
  #combined_data[i,'Cpausedecrease_mean']<-length(ddown[which(ddown<(10*ddown_mean))])
  
  rm(d,dup,dup_med,dup_mean)#,ddown,ddown_med,ddown_mean)
}



combined_data$Npauseincrease_med<-numeric(nrow(combined_data))
#combined_data$Npausedecrease_med<-numeric(nrow(combined_data))
combined_data$Npauseincrease_mean<-numeric(nrow(combined_data))
#combined_data$Npausedecrease_mean<-numeric(nrow(combined_data))

for(i in 1:nrow(combined_data)){
  d<-as.numeric(strsplit(combined_data[i,'not_CDS_counts'],split=',')[[1]])
  dup<-d[d>0]
  
  dup_med<-median(dup)
  dup_mean<-mean(dup)
  combined_data[i,'Npauseincrease_med']<-length(dup[which(dup>(10*dup_med))])
  combined_data[i,'Npauseincrease_mean']<-length(dup[which(dup>(10*dup_mean))])
  
  #ddown<-d[d<0]
  #ddown_med<-median(ddown)
  #ddown_mean<-mean(ddown)
  #combined_data[i,'Npausedecrease_med']<-length(ddown[which(ddown<(10*ddown_med))])
  #combined_data[i,'Npausedecrease_mean']<-length(ddown[which(ddown<(10*ddown_mean))])
  
  rm(d,dup,dup_med,dup_mean) #,ddown,ddown_med,ddown_mean)
}

combined_data$Dpauseincrease_med<-numeric(nrow(combined_data))
combined_data$Dpausedecrease_med<-numeric(nrow(combined_data))
combined_data$Dpauseincrease_mean<-numeric(nrow(combined_data))
combined_data$Dpausedecrease_mean<-numeric(nrow(combined_data))

for(i in 1:nrow(combined_data)){
  d<-as.numeric(strsplit(combined_data[i,'delta_CDS_counts'],split=',')[[1]])
  dup<-d[d>0]
  
  dup_med<-median(dup)
  dup_mean<-mean(dup)
  combined_data[i,'Dpauseincrease_med']<-length(dup[which(dup>(10*dup_med))])
  combined_data[i,'Dpauseincrease_mean']<-length(dup[which(dup>(10*dup_mean))])
  
  ddown<-d[d<0]
  ddown_med<-median(ddown)
  ddown_mean<-mean(ddown)
  combined_data[i,'Dpausedecrease_med']<-length(ddown[which(ddown<(10*ddown_med))])
  combined_data[i,'Dpausedecrease_mean']<-length(ddown[which(ddown<(10*ddown_mean))])
  
  rm(d,dup,ddown,dup_med,dup_mean,ddown_med,ddown_mean)
}



#focus on +45 -15nt -> therefore can be described as translation ..elongation.. rates
#currently only using 28,29,30 - use all read lengths with different offsets? incase of slight digestion differences cause loss of peaks? Also include all frames?
#to identify peaks that are resolved with CNOT1 knockdown
#criteria 1: have 10x greater than the average peak on an mRNA in control conditions
#criteria 2: 10x > than delta change <- directional average? or absolute average?



#to identify peaks that are induced by CNOT1 knockdown
#criteria 1: have 10x greater than the average peak on an mRNA in CNOT1 KD conditions
#criteria 2: 10x > than delta change <- directional average? or absolute average?


#need to show that these two criteria are met at the exact same position
#additional column with exact positions of those matching the criteria which can be used in downstream filtering



##########################################################################################################################################
#########################################################################################################################################
#########################################################################################################################################
names(control_combined)[2]<-'con_CDS_counts'
names(not_combined)[2]<-'not_CDS_counts'

combined_data<-merge(control_combined,not_combined,by='Gene_ID')
combined_data<-merge(combined_data,delta_combined[,c(1,2)],by='Gene_ID')

combined_data$Cpauseincrease_med<-numeric(nrow(combined_data))
combined_data$Cpauseincrease_mean<-numeric(nrow(combined_data))
combined_data$Cpauseincrease_mean_pos<-numeric(nrow(combined_data))


for(i in 1:nrow(combined_data)){
  d<-as.numeric(strsplit(combined_data[i,'con_CDS_counts'],split=',')[[1]])
  dup<-d[d>0]
  
  dup_med<-median(dup)
  dup_mean<-mean(dup)
  combined_data[i,'Cpauseincrease_med']<-length(dup[which(dup>(10*dup_med))])
  combined_data[i,'Cpauseincrease_mean']<-length(dup[which(dup>(10*dup_mean))])
  combined_data[i,'Cpauseincrease_mean_pos']<-paste0(which(d>(10*dup_mean)),collapse=',')
  
  
  rm(d,dup,dup_med,dup_mean)#,ddown,ddown_med,ddown_mean)
}



combined_data$Npauseincrease_med<-numeric(nrow(combined_data))
combined_data$Npauseincrease_mean<-numeric(nrow(combined_data))
combined_data$Npauseincrease_mean_pos<-numeric(nrow(combined_data))

for(i in 1:nrow(combined_data)){
  d<-as.numeric(strsplit(combined_data[i,'not_CDS_counts'],split=',')[[1]])
  dup<-d[d>0]
  
  dup_med<-median(dup)
  dup_mean<-mean(dup)
  combined_data[i,'Npauseincrease_med']<-length(dup[which(dup>(10*dup_med))])
  combined_data[i,'Npauseincrease_mean']<-length(dup[which(dup>(10*dup_mean))])
  combined_data[i,'Npauseincrease_mean_pos']<-paste0(which(d>(10*dup_mean)),collapse=',')
  
  rm(d,dup,dup_med,dup_mean) #,ddown,ddown_med,ddown_mean)
}

combined_data$Dpauseincrease_med<-numeric(nrow(combined_data))
combined_data$Dpausedecrease_med<-numeric(nrow(combined_data))
combined_data$Dpauseincrease_mean<-numeric(nrow(combined_data))
combined_data$Dpausedecrease_mean<-numeric(nrow(combined_data))
combined_data$Dpauseincrease_mean_pos<-numeric(nrow(combined_data))
combined_data$Dpausedecrease_mean_pos<-numeric(nrow(combined_data))

for(i in 1:nrow(combined_data)){
  d<-as.numeric(strsplit(combined_data[i,'delta_CDS_counts'],split=',')[[1]])
  dup<-d[d>0]
  
  dup_med<-median(dup)
  dup_mean<-mean(dup)
  combined_data[i,'Dpauseincrease_med']<-length(dup[which(dup>(10*dup_med))])
  combined_data[i,'Dpauseincrease_mean']<-length(dup[which(dup>(10*dup_mean))])
  combined_data[i,'Dpauseincrease_mean_pos']<-paste0(which(d>(10*dup_mean)),collapse=',')
  
  ddown<-d[d<0]
  ddown_med<-median(ddown)
  ddown_mean<-mean(ddown)
  combined_data[i,'Dpausedecrease_med']<-length(ddown[which(ddown<(10*ddown_med))])
  combined_data[i,'Dpausedecrease_mean']<-length(ddown[which(ddown<(10*ddown_mean))])
  combined_data[i,'Dpausedecrease_mean_pos']<-paste0(which(d<(10*ddown_mean)),collapse=',')
  
  rm(d,dup,ddown,dup_med,dup_mean,ddown_med,ddown_mean)
}

######################################################################################################


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

nrow(sustainedpause)
nrow(inducedpause)
nrow(resolvedpause)

length(unique(sustainedpause$gene_name))
length(unique(inducedpause$gene_name))
length(unique(resolvedpause$gene_name))

write.table(inducedpause$gene_name,'inducedpauses_282930.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(sustainedpause$gene_name,'sustainedpauses_282930.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(resolvedpause$gene_name,'resolvedpauses_282930.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)

resolvedonly<-subset(resolvedpause,(gene_name %in% inducedpause$gene_name)==FALSE)
resolvedonly<-subset(resolvedonly,(gene_name %in% sustainedpause$gene_name)==FALSE)

inducedonly<-subset(inducedpause,(gene_name %in% resolvedpause$gene_name)==FALSE)
inducedonly<-subset(inducedonly,(gene_name %in% sustainedpause$gene_name)==FALSE)

sustainedonly<-subset(sustainedpause,(gene_name %in% resolvedpause$gene_name)==FALSE)
sustainedonly<-subset(sustainedonly,(gene_name %in% inducedpause$gene_name)==FALSE)

nrow(inducedonly)
nrow(resolvedonly)
nrow(sustainedonly)

length(unique(sustainedonly$gene_name))
length(unique(inducedonly$gene_name))
length(unique(resolvedonly$gene_name))

length(unique(subset(pauses,(gene_name %in% sustainedpause$gene_name)==TRUE & (gene_name %in% inducedpause$gene_name)==TRUE & (gene_name %in% resolvedpause$gene_name)==TRUE)$gene_name))
length(unique(subset(pauses,(gene_name %in% inducedpause$gene_name)==TRUE & (gene_name %in% resolvedpause$gene_name)==TRUE)$gene_name))
length(unique(subset(pauses,(gene_name %in% sustainedpause$gene_name)==TRUE & (gene_name %in% resolvedpause$gene_name)==TRUE)$gene_name))
length(unique(subset(pauses,(gene_name %in% inducedpause$gene_name)==TRUE & (gene_name %in% sustainedpause$gene_name)==TRUE)$gene_name))



#look at these in RPF v pSILAC
pS<-merge(RPF27to31,M[,c(1,10)],by.x='Gene_ID',by.y='gene_name')
pS$ind<-rep('other',nrow(pS))
pS$sust<-rep('other',nrow(pS))
pS$res<-rep('other',nrow(pS))
for(i in 1:nrow(pS)){
  if((pS[i,'Gene_ID'] %in% inducedonly$gene_name)==TRUE){
    pS[i,'ind']<-'induced'
  }
  if((pS[i,'Gene_ID'] %in% sustainedonly$gene_name)==TRUE){
    pS[i,'sust']<-'sustained'
  }
  if((pS[i,'Gene_ID'] %in% resolvedonly$gene_name)==TRUE){
    pS[i,'res']<-'resolved'
  }
}

pS$ind<-factor(pS$ind,levels=c('other','induced'))
pS<-pS[order(pS$ind),]
colourset<-c('grey68','darkorange')
f1<-ggplot(pS,aes(x=RPF_log2FC,y=mean_FR,col=ind))+
  geom_point()+theme_bw()+coord_trans(ylim=c(-1.5,1.5),xlim=c(-1.5,1.5))+xlab('log2FC RPF (siCNOT1/siControl)')+ylab('log2FC protein production (siCNOT1/siControl)')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),plot.title=element_text(size=16),legend.title=element_blank())+
  scale_colour_manual(values=colourset)+scale_size_manual(values=c(1,1.5))+
  geom_abline(slope=1,intercept=0,col='black')+geom_abline(slope=1,intercept=0.25,col='black',linetype='dashed')+geom_abline(slope=1,intercept=-0.25,col='black',linetype='dashed')
f1



pS$res<-factor(pS$res,levels=c('other','resolved'))
pS<-pS[order(pS$res),]
colourset<-c('grey68','darkmagenta')
f1<-ggplot(pS,aes(x=RPF_log2FC,y=mean_FR,col=res))+
  geom_point()+theme_bw()+coord_trans(ylim=c(-1.5,1.5),xlim=c(-1.5,1.5))+xlab('log2FC RPF (siCNOT1/siControl)')+ylab('log2FC protein production (siCNOT1/siControl)')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),plot.title=element_text(size=16),legend.title=element_blank())+
  scale_colour_manual(values=colourset)+scale_size_manual(values=c(1,1.5))+
  geom_abline(slope=1,intercept=0,col='black')+geom_abline(slope=1,intercept=0.25,col='black',linetype='dashed')+geom_abline(slope=1,intercept=-0.25,col='black',linetype='dashed')
f1

pS$sust<-factor(pS$sust,levels=c('other','sustained'))
pS<-pS[order(pS$sust),]
colourset<-c('grey68','darkorange')
f1<-ggplot(pS,aes(x=RPF_log2FC,y=mean_FR,col=sust))+
  geom_point()+theme_bw()+coord_trans(ylim=c(-1.5,1.5),xlim=c(-1.5,1.5))+xlab('log2FC RPF (siCNOT1/siControl)')+ylab('log2FC protein production (siCNOT1/siControl)')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),plot.title=element_text(size=16),legend.title=element_blank())+
  scale_colour_manual(values=colourset)+scale_size_manual(values=c(1,1.5))+
  geom_abline(slope=1,intercept=0,col='black')+geom_abline(slope=1,intercept=0.25,col='black',linetype='dashed')+geom_abline(slope=1,intercept=-0.25,col='black',linetype='dashed')
f1


#see presence in the halflifeTE groups = from halflife_and_TE_120720.R
for(i in 1:length(Groups)){
  print(Groups[i])
  print(nrow(subset(halflife_TE,newGroup==Groups[i])))
  print(length(which(subset(halflife_TE,newGroup==Groups[i])$Gene_ID %in% inducedpause$gene_name)==TRUE))
  print(length(which(subset(halflife_TE,newGroup==Groups[i])$Gene_ID %in% sustainedpause$gene_name)==TRUE))
  print(length(which(subset(halflife_TE,newGroup==Groups[i])$Gene_ID %in% resolvedpause$gene_name)==TRUE))
}


########################################
pS<-merge(RPF27to31,M[,c(1,10)],by.x='Gene_ID',by.y='gene_name')
pS$pausetype<-rep('all other mRNAs',nrow(pS))

for(i in 1:nrow(pS)){
  if((pS[i,'Gene_ID'] %in% inducedonly$gene_name)==TRUE){
    pS[i,'pausetype']<-'induced'
  }
  if((pS[i,'Gene_ID'] %in% sustainedonly$gene_name)==TRUE){
    pS[i,'pausetype']<-'sustained'
  }
  if((pS[i,'Gene_ID'] %in% resolvedonly$gene_name)==TRUE){
    pS[i,'pausetype']<-'resolved'
  }
}

pS$pStoRPF<-pS$mean_FR-pS$RPF_log2FC

colourset<-c('darkorange','grey68','darkmagenta','green4')
f1<-ggplot(pS,aes(x=pausetype,y=pStoRPF,fill=pausetype))+
  geom_violin(draw_quantiles=c(0.25,0.5,0.75))+theme_bw()+ylab('protein synthesis to RPF (siCNOT1/siControl)')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),plot.title=element_text(size=16),legend.title=element_blank())+
  scale_fill_manual(values=colourset)
f1

#stats
DT<-dunnTest(RPF_log2FC~newGroup,data=halflife_TE,method='bh')
write.table(DT$res,file='halflifeTEgroups_RPFchange_padj.csv',col.names=TRUE,row.names=FALSE,sep=',')


colourset<-c('darkorange','grey68','darkmagenta','green4')
f1<-ggplot(pS,aes(x=pausetype,y=mean_FR,fill=pausetype))+
  geom_violin(draw_quantiles=c(0.25,0.5,0.75))+theme_bw()+ylab('protein synthesis (siCNOT1/siControl)')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),plot.title=element_text(size=16),legend.title=element_blank())+
  scale_fill_manual(values=colourset)
f1




pS<-merge(RPF27to31,M[,c(1,10)],by.x='Gene_ID',by.y='gene_name')
pS$ind<-rep('other',nrow(pS))
pS$sust<-rep('other',nrow(pS))
pS$res<-rep('other',nrow(pS))
for(i in 1:nrow(pS)){
  if((pS[i,'Gene_ID'] %in% inducedpause$gene_name)==TRUE){
    pS[i,'ind']<-'induced'
  }
  if((pS[i,'Gene_ID'] %in% sustainedpause$gene_name)==TRUE){
    pS[i,'sust']<-'sustained'
  }
  if((pS[i,'Gene_ID'] %in% resolvedpause$gene_name)==TRUE){
    pS[i,'res']<-'resolved'
  }
}

pS$pStoRPF<-pS$mean_FR-pS$RPF_log2FC

colourset<-c('darkorange','grey68')
f1<-ggplot(pS,aes(x=ind,y=pStoRPF,fill=ind))+
  geom_violin(draw_quantiles=c(0.25,0.5,0.75))+theme_bw()+ylab('protein synthesis to RPF (siCNOT1/siControl)')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),plot.title=element_text(size=16),legend.title=element_blank())+
  scale_fill_manual(values=colourset)
f1

colourset<-c('grey68','green4')
f1<-ggplot(pS,aes(x=sust,y=pStoRPF,fill=sust))+
  geom_violin(draw_quantiles=c(0.25,0.5,0.75))+theme_bw()+ylab('protein synthesis to RPF (siCNOT1/siControl)')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),plot.title=element_text(size=16),legend.title=element_blank())+
  scale_fill_manual(values=colourset)
f1

colourset<-c('grey68','darkmagenta')
f1<-ggplot(pS,aes(x=res,y=pStoRPF,fill=res))+
  geom_violin(draw_quantiles=c(0.25,0.5,0.75))+theme_bw()+ylab('protein synthesis to RPF (siCNOT1/siControl)')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),plot.title=element_text(size=16),legend.title=element_blank())+
  scale_fill_manual(values=colourset)
f1

pS$pausetype<-factor(pS$pausetype,levels=c('all other mRNAs','sustained','induced','resolved'))
pS<-pS[order(pS$pausetype),]
colourset<-c('grey68','green4','darkmagenta','darkorange2')

f1<-ggplot(pS,aes(x=pausetype,y=pStoRPF,fill=pausetype))+
  geom_boxplot(outlier.size=0.5)+theme_bw()+ylab('log2FC protein production / log2FC RPF')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.ticks.x=element_blank(),axis.title = element_text(size = 12),axis.text.x=element_blank(),axis.title.x=element_blank(),legend.title=element_blank())+
  scale_fill_manual(values=colourset)+coord_trans(ylim=c(-1,1))
f1
ggsave(plot=f1,filename='pause_PPvRPF.png',height=4,width=4)

#stats
DT<-dunnTest(pStoRPF~pausetype,data=pS,method='bh')
write.table(DT$res,file='pausegroups_PPtoRPF_padj.csv',col.names=TRUE,row.names=FALSE,sep=',')


############################

pS<-halflife_TE
pS<-RPF27to31

pS$ind<-rep('other',nrow(pS))
pS$sust<-rep('other',nrow(pS))
pS$res<-rep('other',nrow(pS))
for(i in 1:nrow(pS)){
  if((pS[i,'Gene_ID'] %in% inducedpause$gene_name)==TRUE){
    pS[i,'ind']<-'induced'
  }
  if((pS[i,'Gene_ID'] %in% sustainedpause$gene_name)==TRUE){
    pS[i,'sust']<-'sustained'
  }
  if((pS[i,'Gene_ID'] %in% resolvedpause$gene_name)==TRUE){
    pS[i,'res']<-'resolved'
  }
}

pS$pausetype<-rep('all other mRNAs',nrow(pS))
for(i in 1:nrow(pS)){
  if((pS[i,'Gene_ID'] %in% inducedonly$gene_name)==TRUE){
    pS[i,'pausetype']<-'induced'
  }
  if((pS[i,'Gene_ID'] %in% sustainedonly$gene_name)==TRUE){
    pS[i,'pausetype']<-'sustained'
  }
  if((pS[i,'Gene_ID'] %in% resolvedonly$gene_name)==TRUE){
    pS[i,'pausetype']<-'resolved'
  }
}

pS$pausetype<-factor(pS$pausetype,levels=c('all other mRNAs','sustained','induced','resolved'))
pS<-pS[order(pS$pausetype),]
colourset<-c('grey68','green4','darkmagenta','darkorange2')

f1<-ggplot(pS,aes(x=pausetype,y=log2FC_halflife,fill=pausetype))+
  geom_violin(draw_quantiles=c(0.25,0.5,0.75))+theme_bw()+ylab('log2FC mRNA half-life (siCNOT1/siControl)')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),plot.title=element_text(size=16),legend.title=element_blank())+
  scale_fill_manual(values=colourset)
f1


f1<-ggplot(pS,aes(x=pausetype,y=log2FC_halflife,fill=pausetype))+
  geom_boxplot(outlier.size=0.5)+theme_bw()+ylab('log2FC mRNA half-life (siCNOT1/siControl)')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.ticks.x=element_blank(),axis.title = element_text(size = 12),axis.text.x=element_blank(),axis.title.x=element_blank(),legend.title=element_blank())+
  scale_fill_manual(values=colourset)+coord_trans(ylim=c(-2,6))
f1
ggsave(plot=f1,filename='pause_log2FChalflife.png',height=4,width=4)

#stats
DT<-dunnTest(log2FC_halflife~pausetype,data=pS,method='bh')
write.table(DT$res,file='pausegroups_loghalflife_padj.csv',col.names=TRUE,row.names=FALSE,sep=',')

f1<-ggplot(pS,aes(x=pausetype,y=siControl_half_life,fill=pausetype))+
  geom_boxplot(outlier.size=0.5)+theme_bw()+ylab('log2FC mRNA half-life (siCNOT1/siControl)')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.ticks.x=element_blank(),axis.title = element_text(size = 12),axis.text.x=element_blank(),axis.title.x=element_blank(),legend.title=element_blank())+
  scale_fill_manual(values=colourset)+coord_trans(ylim=c(0.1,100))+scale_y_log10()
f1

f1<-ggplot(pS,aes(x=pausetype,y=siCNOT1_half_life,fill=pausetype))+
  geom_boxplot(outlier.size=0.5)+theme_bw()+ylab('log2FC mRNA half-life (siCNOT1/siControl)')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.ticks.x=element_blank(),axis.title = element_text(size = 12),axis.text.x=element_blank(),axis.title.x=element_blank(),legend.title=element_blank())+
  scale_fill_manual(values=colourset)+coord_trans(ylim=c(0.1,100))+scale_y_log10()
f1

colourset<-c('darkorange','grey68','darkmagenta','green4')
f1<-ggplot(pS,aes(x=pausetype,y=RPF_TOT_DIFF,fill=pausetype))+
  geom_violin(draw_quantiles=c(0.25,0.5,0.75))+theme_bw()+ylab('log2FC TE (siCNOT1/siControl)')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),plot.title=element_text(size=16),legend.title=element_blank())+
  scale_fill_manual(values=colourset)
f1


colourset<-c('darkorange','grey68','darkmagenta','green4')
f1<-ggplot(pS,aes(x=pausetype,y=RPF_log2FC,fill=pausetype))+
  geom_violin(draw_quantiles=c(0.25,0.5,0.75))+theme_bw()+ylab('log2FC RPF (siCNOT1/siControl)')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),plot.title=element_text(size=16),legend.title=element_blank())+
  scale_fill_manual(values=colourset)
f1

colourset<-c('darkorange','grey68','darkmagenta','green4')
f1<-ggplot(pS,aes(x=pausetype,y=RNA_log2FC,fill=pausetype))+
  geom_violin(draw_quantiles=c(0.25,0.5,0.75))+theme_bw()+ylab('log2FC RNA (siCNOT1/siControl)')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),plot.title=element_text(size=16),legend.title=element_blank())+
  scale_fill_manual(values=colourset)
f1


pS$pausetype<-factor(pS$pausetype,levels=c('other','sustained','induced','resolved'))
pS<-pS[order(pS$pausetype),]
colourset<-c('grey68','green4','darkmagenta','darkorange')
f1<-ggplot(pS,aes(y=RPF_TOT_DIFF,x=log2FC_halflife,col=pausetype))+
  geom_point()+theme_bw()+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),plot.title=element_text(size=16),legend.title=element_blank())+
  scale_colour_manual(values=colourset)
f1



pS$pausetype<-factor(pS$pausetype,levels=c('other','sustained','induced','resolved'))
pS<-pS[order(pS$pausetype),]
colourset<-c('grey68','grey68','grey68','darkmagenta')
f1<-ggplot(pS,aes(y=RPF_TOT_DIFF,x=log2FC_halflife,col=pausetype))+
  geom_point()+theme_bw()+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),plot.title=element_text(size=16),legend.title=element_blank())+
  scale_colour_manual(values=colourset)
f1

pS$pausetype<-factor(pS$pausetype,levels=c('other','induced','resolved','sustained'))
pS<-pS[order(pS$pausetype),]
colourset<-c('grey68','grey68','grey68','green4')
f1<-ggplot(pS,aes(y=RPF_TOT_DIFF,x=log2FC_halflife,col=pausetype))+
  geom_point()+theme_bw()+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),plot.title=element_text(size=16),legend.title=element_blank())+
  scale_colour_manual(values=colourset)
f1

pS$pausetype<-factor(pS$pausetype,levels=c('other','sustained','resolved','induced'))
pS<-pS[order(pS$pausetype),]
colourset<-c('grey68','grey68','grey68','darkorange')
f1<-ggplot(pS,aes(y=RPF_TOT_DIFF,x=log2FC_halflife,col=pausetype))+
  geom_point()+theme_bw()+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), axis.title = element_text(size = 12),plot.title=element_text(size=16),legend.title=element_blank())+
  scale_colour_manual(values=colourset)
f1




################# make supplemental table of pauses ##################
nrow(inducedonly)
nrow(resolvedonly)
nrow(sustainedonly)

inducedonly$pausetype<-rep('induced',nrow(inducedonly))
resolvedonly$pausetype<-rep('resolved',nrow(resolvedonly))
sustainedonly$pausetype<-rep('sustained',nrow(sustainedonly))
pausesupptable<-rbind(inducedonly[,c(1,13,19)],resolvedonly[,c(1,13,19)],sustainedonly[,c(1,13,19)])

names(pausesupptable)<-c('gene_name','nt_position_of_pause_in_CDS','pause_type')
write.csv(pausesupptable,file='Supplemental_Table_6.csv',quote=FALSE,sep=',',col.names = TRUE,row.names=FALSE)









