#load libraries
library(seqinr)
library(Biostrings)
library(ggplot2)

#other functions
searchFeatures<-function(x,feature){
    section<-strsplit(x,split='|',fixed=TRUE)[[1]][[grep(feature,strsplit(x,split='|',fixed=TRUE)[[1]])]]
    return(section)
  }


###################### get data on frame of reads #####################
#input info
directory<-'/CountsOutput/'
allsamplenames<-c('con1, con2, con3')
lengths<-seq(27,31,1) #if have run the counts python script separately for multiple lengths

#based on format of alignment with bowtie to gencode protein-coding transcriptome
#get frame of read start positions in relation to translation start codon nts
for(x in 1:length(allsamplenames)){
  samplename<-allsamplenames[x]

  #read in data
  samples<-NULL
  for(i in 1:length(lengths)){
    samples[[i]]<-read.delim(paste0(directory,'RPFcountsTable_',samplename,'_',lengths[i],'.txt'),stringsAsFactors = FALSE,header=TRUE)
  } 
  lapply(samples, function(x) nrow(x))

  #sort input format       
  for(k in 1:length(samples)){
    s<-samples[[k]]
    s$ENST<-sapply(strsplit(as.character(s[,'Transcript_ID']), "|",fixed=TRUE), "[[", 1)
    s$ENSG<-sapply(strsplit(as.character(s[,'Transcript_ID']), "|",fixed=TRUE), "[[", 2)
    s$Gene_ID<-sapply(strsplit(as.character(s[,'Transcript_ID']), "|",fixed=TRUE), "[[", 6)
    s$length<-sapply(strsplit(as.character(s[,'Transcript_ID']), "|",fixed=TRUE), "[[", 7)
    
    #get ends of each region (5'UTR, CDS and 3'UTR) position across the mRNA
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
    samples[[k]]<-s
  }
  lapply(samples, function(x) nrow(x))
  
  #check frame of CDS aligned reads       
  for(m in 1:length(samples)){
    s<-samples[[m]]
    s$UTR5pos<-as.numeric(s$UTR5pos)
    s$CDSpos<-as.numeric(s$CDSpos)
    s$UTR3pos<-as.numeric(s$UTR3pos)    
    s<-subset(s,UTR5pos>30)
    s<-subset(s,(UTR3pos-CDSpos)>30)
    s$frame0<-numeric(nrow(s))
    s$frame1<-numeric(nrow(s))
    s$frame2<-numeric(nrow(s))
    for(n in 1:nrow(s)){
      region<-as.numeric(strsplit(s[n,'pos_RPF_counts'],' ')[[1]])
      regionx<-region[c((s[n,'UTR5pos']-30):(s[n,'CDSpos']+30))] #just look at CDS reads
      s[n,'UTR5_sum']<-sum(region[1:s[n,'UTR5pos']])
      s[n,'CDS_sum']<-sum(region[(s[n,'UTR5pos']+1):(s[n,'CDSpos'])])
      s[n,'UTR3_sum']<-sum(region[(s[n,'CDSpos']+1):(s[n,'UTR3pos'])])
      f0<-regionx[seq(1,length(regionx),3)]
      f1<-regionx[seq(2,length(regionx),3)]
      f2<-regionx[seq(3,length(regionx),3)]
      s[n,'f0']<-sum(f0)
      s[n,'f1']<-sum(f1)
      s[n,'f2']<-sum(f2)
    }
    samples[[m]]<-s
    print(lengths[m])
    print(sum(as.numeric(s$Total_RPF_counts)))
    print(c(sum(s$UTR5_sum),sum(s$CDS_sum),sum(s$UTR3_sum)))
    frameline<-paste0(sum(s$f0,na.rm=TRUE),'\t',sum(s$f1,na.rm=TRUE),'\t',sum(s$f2,na.rm=TRUE))
    write(x=frameline,file=paste0(samplename,'_',lengths[m],'_frameinfo.txt'),sep='\n',append=TRUE)

  }
}  



####################### P-site offset plots ######################

#input info
directory<-'/CountsOutput/'
allsamplenames<-c('con1, con2, con3')
lengths<-seq(27,31,1) #if have run the counts python script separately for multiple lengths

#based on format of alignment with bowtie to gencode protein-coding transcriptome
#get frame of read start positions in relation to translation start codon nts
for(x in 1:length(allsamplenames)){
  samplename<-allsamplenames[x]

  #read in data
  samples<-NULL
  for(i in 1:length(lengths)){
    samples[[i]]<-read.delim(paste0(directory,'RPFcountsTable_',samplename,'_',lengths[i],'.txt'),stringsAsFactors = FALSE,header=TRUE)
  } 
  lapply(samples, function(x) nrow(x))

  #sort input format       
  for(k in 1:length(samples)){
    s<-samples[[k]]
    s$ENST<-sapply(strsplit(as.character(s[,'Transcript_ID']), "|",fixed=TRUE), "[[", 1)
    s$ENSG<-sapply(strsplit(as.character(s[,'Transcript_ID']), "|",fixed=TRUE), "[[", 2)
    s$Gene_ID<-sapply(strsplit(as.character(s[,'Transcript_ID']), "|",fixed=TRUE), "[[", 6)
    s$length<-sapply(strsplit(as.character(s[,'Transcript_ID']), "|",fixed=TRUE), "[[", 7)
    
    #get ends of each region (5'UTR, CDS and 3'UTR) position across the mRNA
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
    samples[[k]]<-s
  }
  lapply(samples, function(x) nrow(x))
  
  #reads around the AUG to determine offset length
  offset_tables<-NULL
  for(k in 1:length(samples)){
    ss<-samples[[k]]
    ss$UTR5pos<-as.numeric(ss$UTR5pos)
    ss$CDSpos<-as.numeric(ss$CDSpos)
    ss$UTR3pos<-as.numeric(ss$UTR3pos)
    print(nrow(ss))
    ss<-subset(ss,UTR5pos>32)
    print(nrow(ss))
    
    #make plots of RPF read starts around AUG
    offsetplot<-data.frame(matrix(ncol=60,nrow=nrow(ss),0),stringsAsFactors = FALSE)
    for(n in 1:nrow(ss)){
      region<-as.numeric(strsplit(ss[n,'pos_RPF_counts'],' ')[[1]])
      regionx<-region[(ss[n,'UTR5pos']-29):(ss[n,'UTR5pos']+30)]
      offsetplot[n,]<-regionx
    }
    offset_tables[[k]]<-offsetplot
    offsetplotdata<-data.frame(matrix(ncol=2,nrow=60,0),stringsAsFactors = FALSE)
    names(offsetplotdata)<-c('position','counts')
    offsetplotdata$position<-seq(-30,29,1)
    offsetplotdata$counts<-unlist(apply(offsetplot,2, function(x) mean(x,na.rm=TRUE)))
    print(lengths(k))
    print(offsetplotdata[which.max(offsetplotdata$counts),'position'])
    f1<-ggplot(offsetplotdata,aes(y=counts,x=position))+theme_classic()+geom_line(size=1)+ylab('no. RPFs x 10^6')+xlab('nucleotides from AUG')
    ggsave(paste0(samplename,'_',lengths[k],'.png'),width=4,height=2.5)
  }

}  

