#load libraries
library(Biostrings)
library(seqinr)


#Sort ID conversion to gene name (can get from fasta file details) #####################
fastaFile <- readDNAStringSet("~/Desktop/CNOTdataAnalysis/CountsFilesToUse/gencode.v28.transcripts.fa")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
humanDF_alltranscripts <- data.frame(seq_name, sequence,stringsAsFactors = FALSE)
humanDF_alltranscripts$ENST<-sapply(strsplit(as.character(humanDF_alltranscripts[,'seq_name']), " ",fixed=TRUE), "[[", 2)

fastaFile <- readDNAStringSet("~/Desktop/CNOTdataAnalysis/CountsFilesToUse/gencode.v28.pc_transcripts.fa")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
humanDF <- data.frame(seq_name, sequence,stringsAsFactors = FALSE)

#convert IDs and split the sequence 5UTR/CDS/3UTR
humanDF$ENST<-sapply(strsplit(as.character(humanDF[,'seq_name']), "|",fixed=TRUE), "[[", 1)
humanDF$ENSG<-sapply(strsplit(as.character(humanDF[,'seq_name']), "|",fixed=TRUE), "[[", 2)
humanDF$Gene_ID<-sapply(strsplit(as.character(humanDF[,'seq_name']), "|",fixed=TRUE), "[[", 6)
humanDF$length<-sapply(strsplit(as.character(humanDF[,'seq_name']), "|",fixed=TRUE), "[[", 7)

searchFeatures<-function(x,feature){
  section<-strsplit(x,split='|',fixed=TRUE)[[1]][[grep(feature,strsplit(x,split='|',fixed=TRUE)[[1]])]]
  return(section)
}

humanDF$UTR5pos<-rep(NA,nrow(humanDF))
humanDF$CDSpos<-rep(NA,nrow(humanDF))
humanDF$UTR3pos<-rep(NA,nrow(humanDF))
for(i in 1:nrow(humanDF)){
  info<-humanDF[i,1]
  
  if(grepl('CDS:',info)==TRUE){
    CDSp<-searchFeatures(info,'CDS:') #add : to deal with problem in Gene_IDs
    humanDF[i,'CDSpos']<-strsplit(CDSp,split='-',fixed=TRUE)[[1]][2]
  }
  if(grepl('UTR5',info)==TRUE){
    UTR5p<-searchFeatures(info,'UTR5:')
    humanDF[i,'UTR5pos']<-strsplit(UTR5p,split='-',fixed=TRUE)[[1]][2]
  }
  if(grepl('UTR3',info)==TRUE){
    UTR3p<-searchFeatures(info,'UTR3:')
    humanDF[i,'UTR3pos']<-strsplit(UTR3p,split='-',fixed=TRUE)[[1]][2]
  }
}

humanDF$UTR5_sequence<-rep('',nrow(humanDF))
humanDF$CDS_sequence<-rep('',nrow(humanDF))
humanDF$UTR3_sequence<-rep('',nrow(humanDF))

for(i in 1:nrow(humanDF)){
  UTR5p<-as.numeric(humanDF[i,'UTR5pos'])
  CDSp<-as.numeric(humanDF[i,'CDSpos'])
  UTR3p<-as.numeric(humanDF[i,'UTR3pos'])
  allseq<-humanDF[i,'sequence']
  
  if((is.na(UTR5p)==FALSE) & (is.na(UTR3p)==FALSE)){
    humanDF[i,'UTR5_sequence']<-substring(allseq,1,UTR5p)
    humanDF[i,'CDS_sequence']<-substring(allseq,(UTR5p+1),CDSp)
    humanDF[i,'UTR3_sequence']<-substring(allseq,(CDSp+1),UTR3p)
  }else{
    if((is.na(UTR5p)==TRUE) & (is.na(UTR3p)==FALSE)){
      humanDF[i,'UTR5_sequence']<-NA
      humanDF[i,'CDS_sequence']<-substring(allseq,1,CDSp)
      humanDF[i,'UTR3_sequence']<-substring(allseq,(CDSp+1),UTR3p)
    }else{
      if((is.na(UTR5p)==FALSE) & (is.na(UTR3p)==TRUE)){
        humanDF[i,'UTR5_sequence']<-substring(allseq,1,UTR5p)
        humanDF[i,'CDS_sequence']<-substring(allseq,(UTR5p+1),CDSp)
        humanDF[i,'UTR3_sequence']<-NA
      }else{2
        if((is.na(UTR5p)==TRUE) & (is.na(UTR3p)==TRUE)){
          humanDF[i,'UTR5_sequence']<-NA
          humanDF[i,'CDS_sequence']<-substring(allseq,1,CDSp)
          humanDF[i,'UTR3_sequence']<-NA
        }
      }
    }  
  }
}

nrow(humanDF_alltranscripts)
nrow(humanDF) #protein coding transcripts
humanDFx<-merge(humanDF,humanDF_alltranscripts[,c(1,3)],by=c('ENST'))
nrow(humanDFx)
humanDFx<-subset(humanDFx,is.na(humanDFx$UTR5_sequence)==FALSE)
nrow(humanDFx)
humanDFx<-subset(humanDFx,is.na(humanDFx$UTR3_sequence)==FALSE) #most mRNAs lost in this step
nrow(humanDFx)
humanDFx<-subset(humanDFx,(backports::startsWith(humanDFx$CDS_sequence,'ATG')==TRUE))
nrow(humanDFx)
length(unique(humanDFx$Gene_ID)) #19039 unique genes



#####
#Sort one transcript per gene ###############################
#based on Controls? #also look if different transcripts change? And corresponding RPFs?
#picking the most abundant transcript per gene -> counts / length (then pick the greatest for gene X)

filenames<-c('Lib-1','Lib-3','Lib-5')
geneIDbased<-list()
TIDbased<-list()
for(i in 1:length(filenames)){
  geneIDbased[[i]]<-read.table(paste0('~/Desktop/CNOTdataAnalysis/CountsFilesToUse/Total_geneID_unique_transcriptomeonly_FS_',filenames[i]),stringsAsFactors = FALSE,header=FALSE)
  names(geneIDbased[[i]])<-c('ENSG','ENSG_counts')
  geneIDbased[[i]]<-merge(humanDFx,geneIDbased[[i]],by='ENSG')
  TIDbased[[i]]<-read.table(paste0('~/Desktop/CNOTdataAnalysis/CountsFilesToUse/Total_TID_nonunique_transcriptomeonly_FS_',filenames[i]),stringsAsFactors = FALSE,header=FALSE)
  names(TIDbased[[i]])<-c('ENST','ENST_counts')
  TIDbased[[i]]<-merge(humanDFx,TIDbased[[i]],by='ENST')
}

#mean of all 3 replicates
allControls<-merge(TIDbased[[1]],TIDbased[[2]][,c(1,14)],by='ENST')
allControls<-merge(allControls,TIDbased[[3]][c(1,14)],by='ENST')
#allControls$averageCountsENSG<-unlist(lapply(allControl[,c('ENSG_counts','ENSG_counts.x','ENSG_counts.y')], function(x) mean(x)))
allControls$averageCountsENST<-unlist(apply(allControls[,c('ENST_counts','ENST_counts.x','ENST_counts.y')], 1, function(x) mean(x)))
#allControls$averageCountsENSG_lenadj<-allControls$averageCountsENSG/nchar(allControls$sequence)
allControls$averageCountsENST_lenadj<-allControls$averageCountsENST/nchar(allControls$sequence)


keep_most_abundant<-function(dataset){
  new_dataset<-NULL
  allgenes<-unique(dataset$Gene_ID)
  for(k in 1:length(allgenes)){
    geneset<-subset(dataset,Gene_ID==allgenes[k])
    geneset<-geneset[order('averageCountsENST_lenadj',asc=FALSE),]
    new_dataset<-rbind(new_dataset,geneset[1,])
  }
  return(new_dataset)
} # this is not the most efficient way of doing it but it works and in theory only needs to be done once

allControls_onetranscript<-keep_most_abundant(allControls)

humanDFx<-subset(humanDFx,(humanDFx$ENST %in% allControls_onetranscript$ENST)==TRUE)



