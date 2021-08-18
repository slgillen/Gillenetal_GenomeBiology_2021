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

########## RPF data - normalise for library size ############
sample_names<-c('siControl_rep1','siControl_rep2','siControl_rep3','siCNOT1_rep1','siCNOT1_rep2','siCNOT1_rep3')
RPF_counts<-NULL
for(i in 1:length(sample_names)){
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
    CDSsection_p2<-CDSsection_p1[(UTR5pos-11):(CDSpos-11)]# includes start and stop
    RPF_Total[[i]][j,'CDS_section_pos_counts']<-paste0(CDSsection_p2,collapse=',')
    CDSsection_p2_InF<-CDSsection_p2[seq(1,length(CDSsection_p2),3)]
    RPF_Total[[i]][j,'CDS_section_InF_pos_counts']<-paste0(CDSsection_p2_InF,collapse=',')
  }
}

control_data<-RPF_Total[1:3]
unique_genes<-unique(control_data[[1]]$Gene_ID)
length(unique_genes)
unique_genes<-unique_genes[(unique_genes %in% control_data[[2]]$Gene_ID)==TRUE]
length(unique_genes)
unique_genes<-unique_genes[(unique_genes %in% control_data[[3]]$Gene_ID)==TRUE]
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
delta_combined<-data.frame(matrix(nrow=length(unique_genes),ncol=3,0),stringsAsFactors = FALSE)
names(delta_combined)<-c('Gene_ID','delta_CDS_counts','ratio_CDS_counts')

for(i in 1:length(unique_genes)){
  delta_combined[i,'Gene_ID']<-unique_genes[i]
  
  r1C<-as.numeric(strsplit(subset(control_data[[1]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r1N<-as.numeric(strsplit(subset(not_data[[1]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r1D<-r1N-r1C
  r1Cx<-r1C+(1/subset(control_data[[1]],Gene_ID==genes[j])$TPM)
  r1Nx<-r1N+(1/subset(not_data[[1]],Gene_ID==genes[j])$TPM)
  r1R<-log(r1Nx/r1Cx,2)
  
  r2C<-as.numeric(strsplit(subset(control_data[[2]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r2N<-as.numeric(strsplit(subset(not_data[[2]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r2D<-r2N-r2C
  r2Cx<-r2C+(1/subset(control_data[[2]],Gene_ID==genes[j])$TPM)
  r2Nx<-r2N+(1/subset(not_data[[2]],Gene_ID==genes[j])$TPM)
  r2R<-log(r2Nx/r2Cx,2)
  
  r3C<-as.numeric(strsplit(subset(control_data[[3]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r3N<-as.numeric(strsplit(subset(not_data[[3]], Gene_ID==unique_genes[i])$CDS_section_InF_pos_counts,split=',')[[1]])
  r3D<-r3N-r3C
  r3Cx<-r3C+(1/subset(control_data[[3]],Gene_ID==genes[j])$TPM)
  r3Nx<-r3N+(1/subset(not_data[[3]],Gene_ID==genes[j])$TPM)
  r3R<-log(r3Nx/r3Cx,2)
  
  
  delta_combined[i,'delta_CDS_counts']<-paste(colMeans(rbind(r1D,r2D,r3D)),collapse=',')
  delta_combined[i,'ratio_CDS_counts']<-paste(colMeans(rbind(r1R,r2R,r3R)),collapse=',')
  rm(r1C,r1N,r1D,r1Cx,r1Nx,r1R)
  rm(r2C,r2N,r2D,r2Cx,r2Nx,r2R)
  rm(r3C,r3N,r3D,r3Cx,r3Nx,r3R)
}



##########################################################

#get dataframe of positions of RPFs
get_plotdf_individual<-function(countsdf,readtype){
  #readcounts<-c(strsplit(countsdf[1,'UTR5_counts'],split=',')[[1]],strsplit(countsdf[1,'CDS_counts'],split=',')[[1]])
  readcounts<-strsplit(countsdf[1,readtype],split=',')[[1]]
  readcounts<-as.numeric(readcounts)
  
  output<-data.frame(matrix(ncol=2,nrow=length(readcounts),0))
  #count<-count
  output[,2]<-readcounts
  output[,1]<-seq(1,length(readcounts),1)
  
  names(output)<-c('position','value')
  output$position<-as.numeric(output$position)
  output$value<-as.numeric(output$value)
  return(output)
}


#need to add counts for 5'UTR, CDS - split plots?
setwd("//john-doe/gw/Systems/Sarah/CNOTpaper/cotranslationalassembly/individual_pauses")

genelist<-bCATgenes
colours<-c('grey48','darkmagenta')

genelist<-c('ATF4', 'CNTRL', 'CIR1', 'CEP126', 'CEP295', 'MKKS', 'CENPJ', 'CLIP1', 'CEP72', 'HOXB4', 'SLF1', 'CCDC66', 'SDCCAG8', 'SSNA1', 'CEP152', 'PLK3', 'SMAD7', 'ID1', 'CETN3', 'SCLT1', 'IFT74', 'BOD1L1', 'SAC3D1', 'CSPP1', 'IFT27', 'WDR13', 'OFD1', 'CCDC77', 'PDE4DIP', 'AXIN2', 'PLK2', 'CEP162', 'POC5', 'C2CD5', 'PLEKHA7', 'LRRC45', 'TUBGCP5', 'MPLKIP', 'ZFYVE26', 'PRKACB', 'KIF20B', 'NUP62', 'KLHL22', 'LATS2', 'GNAI2', 'CCDC61', 'TSEN2', 'PLK4', 'TTC19', 'RAB11FIP3', 'KATNA1', 'IFT22', 'RGS14', 'NDC80', 'GEN1', 'ASPM', 'DLGAP5', 'KIZ', 'PROCR', 'ERC1', 'IFT52', 'HAUS6', 'KCTD1', 'CCDC14', 'TAF1D', 'IFT20', 'KATNAL1', 'DCTN5', 'TCP1', 'NEK1', 'ANKRD26', 'TRAF5', 'STIL', 'SPAST', 'CALM3', 'TUBG2', 'SPICE1', 'HMMR', 'TRIP4', 'DDAH2', 'NPHP4', 'XRCC2', 'VPS37A', 'TMEM67', 'LEO1', 'KIF3B', 'SNAP29', 'CCDC8', 'MZT2A', 'DCLRE1B', 'BBS9', 'FNIP2', 'TUBE1', 'IFT46', 'OLA1', 'CDKL2')
colours<-c('grey48','blue')

genelist<-subset(plotdata,proteingroup=='protein > RPF' & (gene_name %in% halflife_TE$Gene_ID)==FALSE)$gene_name #green = elongation block removed or 
length(genelist)
colours<-c('grey48','green4')

genelist<-subset(plotdata,proteingroup=='protein < RPF' & (gene_name %in% halflife_TE$Gene_ID)==FALSE)$gene_name  #red = induced elongation block
length(genelist)
colours<-c('grey48','red3')

genelist<-c('KIF23','UBE4B','SLC7A5','NEK7')
colours<-c('grey58','darkorange2')

genelist<-c('CHMP5','NAMPT','MAVS')
colours<-c('grey58','darkmagenta')

genelist<-c('PUS1','GAB1','N4BP1','USP3') 
colours<-c('grey58','springgreen4')

genelist<-c('ILVBL','USP15','WDR73','NVL')
colours<-c('grey58','violetred')

genelist<-c('ILVBL','PUS1')
colours<-c('grey58','violetred')

genelist<-c('PSMC2')
colours<-c('grey58','red3')



for(k in 1:length(genelist)){ #might be good to try with total mRNA TPM normalisation
  gene<-genelist[k]
  
  if((gene %in% delta_combined$Gene_ID)==TRUE){
    z<-subset(control_combined,Gene_ID==gene)
    print(nrow(z))
    RPFdata<-get_plotdf_individual(z,'CDS_counts') 
    RPFdata$group<-rep('mRNA',nrow(RPFdata))  
    
    
    gene<-genelist[k]
    z<-subset(not_combined,Gene_ID==gene)
    print(nrow(z))
    RPFdata2<-get_plotdf_individual(z,'CDS_counts') 
    RPFdata2$group<-rep('mRNA',nrow(RPFdata2))  
    
    RPFdata$condition<-rep('siControl',nrow(RPFdata))
    RPFdata2$condition<-rep('siCNOT1',nrow(RPFdata2))
    findmax<-max(RPFdata$value,RPFdata2$value)
    
    
    f1<-ggplot(data=RPFdata, aes(x=position,y=value,fill=group,col=group))+geom_bar(position = 'dodge2',stat='identity')+scale_x_continuous(expand=c(0.025,0))+ylab('RPF read count')+ggtitle(gene)+xlab('position in CDS (codons)')+ylab('ribosome occupancy')
    f1<-f1+scale_colour_manual(values=colours[1])+theme_bw()+theme(legend.text=element_text(size=10),axis.line = element_line(colour = "black"),axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8), axis.title = element_text(size = 8), panel.grid.major = element_blank(),panel.grid.minor=element_blank(),plot.title=element_text(size=10),legend.title=element_blank(),legend.position='none',panel.border=element_blank())+coord_trans(limy=c(0,findmax))
    
    filenamex<-paste(gene,'_siControl_TPMnorm.png',sep='',collapse='')
    ggsave(filename=filenamex,plot=f1,width=3.4, height=1.4)
    
    
    ##################################
    
    
    f1<-ggplot(data=RPFdata2, aes(x=position,y=value,fill=group,col=group))+geom_bar(position = 'dodge2',stat='identity')+scale_x_continuous(expand=c(0.025,0))+ylab('RPF read count')+ggtitle(gene)+xlab('position in CDS (codons)')+ylab('ribosome occupancy')
    f1<-f1+scale_colour_manual(values=colours[2])+theme_bw()+theme(legend.text=element_text(size=10),axis.line = element_line(colour = "black"),axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8), axis.title = element_text(size = 8), panel.grid.major = element_blank(),panel.grid.minor=element_blank(),plot.title=element_text(size=10),legend.title=element_blank(),legend.position='none',panel.border=element_blank())+coord_trans(limy=c(0,findmax))
    
    
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
    
    f1<-ggplot(data=RPFdata3, aes(x=position,y=value,fill=group,col=group))+geom_bar(position = 'dodge2',stat='identity')+scale_x_continuous(expand=c(0.025,0))+ylab('delta')+ggtitle(gene)+xlab('position in CDS (codons)')+ylab('delta')
    f1<-f1+scale_colour_manual(values=colours[2])+theme_bw()+theme(legend.text=element_text(size=10),axis.line = element_line(colour = "black"),axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8), axis.title = element_text(size = 8), panel.grid.major = element_blank(),panel.grid.minor=element_blank(),plot.title=element_text(size=10),legend.title=element_blank(),legend.position='none',panel.border=element_blank())+coord_trans(limy=c(findmin,findmax))
    
    
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
