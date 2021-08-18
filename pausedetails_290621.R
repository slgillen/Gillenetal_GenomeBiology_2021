
#corresponds with data from pausedata_newanalysis_100621.R
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










########################################################################
########################################################################
########################################################################

#check sequences
checkseq<-RPF_counts[[1]][,c(7,8:10,2)]
checkseq$UTR5pos<-as.numeric(checkseq$UTR5pos)
checkseq$CDSpos<-as.numeric(checkseq$CDSpos)
checkseq$CDSlen<-as.numeric(checkseq$CDSpos)-as.numeric(checkseq$UTR5pos)
summary(checkseq$CDSlen)
checkseq$CDSseq<-rep('',nrow(checkseq))
for(i in 1:nrow(checkseq)){
  checkseq[i,'CDSseq']<-substring(checkseq[i,'Transcript_Sequence'],checkseq[i,'UTR5pos']+1, checkseq[i,'CDSpos'])
}



#used this one
pauses<-merge(pauses,checkseq,by.x='gene_name',by.y='Gene_ID')
pauses$Psitepos<-((pauses$pos*3)+45)-2
pauses$CDSlen<-nchar(pauses$CDSseq)
pauses$percCDS<-(pauses$Psitepos/pauses$CDSlen)*100
summary(pauses$percCDS)

#
pauses$pausetype<-rep('all other mRNAs',nrow(pauses))
for(i in 1:nrow(pauses)){
  if((pauses[i,'gene_name'] %in% inducedonly$gene_name)==TRUE){
    pauses[i,'pausetype']<-'induced'
  }
  if((pauses[i,'gene_name'] %in% sustainedonly$gene_name)==TRUE){
    pauses[i,'pausetype']<-'sustained'
  }
  if((pauses[i,'gene_name'] %in% resolvedonly$gene_name)==TRUE){
    pauses[i,'pausetype']<-'resolved'
  }
}

pauses$pausetype<-factor(pauses$pausetype,levels=c('all other mRNAs','sustained','induced','resolved'))
pauses<-pauses[order(pauses$pausetype),]
colourset<-c('grey68','green4','darkmagenta','darkorange2')

f1<-ggplot(pauses,aes(x=percCDS,col=pausetype))+
  geom_density()+theme_bw()+xlab('% CDS')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), 
                        axis.title = element_text(size = 12),legend.title=element_blank())+
  scale_colour_manual(values=colourset)+coord_trans(xlim=c(0,100))
f1
ggsave(plot=f1,filename='pause_distribution.png',height=3,width=4.5)


########################should check have actually got codons and haven't accidentally taken all out-of-frame
pauses$EPA<-rep('',nrow(pauses))
pauses$usEPAds<-rep('',nrow(pauses))

pauses$EPA_aminoacid<-rep('',nrow(pauses))
pauses$usEPAds_aminoacid<-rep('',nrow(pauses))

for(i in 1:nrow(pauses)){
  pauses[i,'EPA']<-substring(pauses[i,'CDSseq'],pauses[i,'Psitepos']-3,pauses[i,'Psitepos']+5)
  pauses[i,'EPA_aminoacid']<-paste0(translate(s2c(substring(pauses[i,'CDSseq'],pauses[i,'Psitepos']-3,pauses[i,'Psitepos']+5)),frame=0,sens = "F", numcode = 1, NAstring = "X", ambiguous = FALSE),collapse='')
  
}

for(i in 1:nrow(pauses)){
  pauses[i,'usEPAds']<-substring(pauses[i,'CDSseq'],pauses[i,'Psitepos']-12,pauses[i,'Psitepos']+14) #3 codons either side
  pauses[i,'usEPAds_aminoacid']<-paste0(translate(s2c(substring(pauses[i,'CDSseq'],pauses[i,'Psitepos']-12,pauses[i,'Psitepos']+14)),frame=0,sens = "F", numcode = 1, NAstring = "X", ambiguous = FALSE),collapse='')
  
}


sustainedpause<-subset(pauses,con_peak=='Y' & not_peak=='Y' & delta_decrease=='N' & delta_increase=='N')
inducedpause<-subset(pauses,con_peak=='N' & not_peak=='Y' & delta_decrease=='N' & delta_increase=='Y')
resolvedpause<-subset(pauses,con_peak=='Y' & not_peak=='N' & delta_decrease=='Y' & delta_increase=='N')

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

#write fasta files for seq2logo web app use
library(seqinr)
write.fasta(sequences=as.list(inducedonly$EPA), names=seq(1,nrow(inducedonly),by=1), file.out='inducedonly_pause_EPA.fa', open = "w", nbchar = 60, as.string = FALSE)
write.fasta(sequences=as.list(inducedonly$EPA_aminoacid), names=seq(1,nrow(inducedonly),by=1), file.out='inducedonly_pause_EPA_aminoacid.fa', open = "w", nbchar = 60, as.string = FALSE)
write.fasta(sequences=as.list(inducedonly$usEPAds), names=seq(1,nrow(inducedonly),by=1), file.out='inducedonly_pause_usEPAds.fa', open = "w", nbchar = 60, as.string = FALSE)
write.fasta(sequences=as.list(inducedonly$usEPAds_aminoacid), names=seq(1,nrow(inducedonly),by=1), file.out='inducedonly_pause_usEPAds_aminoacid.fa', open = "w", nbchar = 60, as.string = FALSE)

write.fasta(sequences=as.list(sustainedonly$EPA), names=seq(1,nrow(sustainedonly),by=1), file.out='sustainedonly_pause_EPA.fa', open = "w", nbchar = 60, as.string = FALSE)
write.fasta(sequences=as.list(sustainedonly$EPA_aminoacid), names=seq(1,nrow(sustainedonly),by=1), file.out='sustainedonly_pause_EPA_aminoacid.fa', open = "w", nbchar = 60, as.string = FALSE)
write.fasta(sequences=as.list(sustainedonly$usEPAds), names=seq(1,nrow(sustainedonly),by=1), file.out='sustainedonly_pause_usEPAds.fa', open = "w", nbchar = 60, as.string = FALSE)
write.fasta(sequences=as.list(sustainedonly$usEPAds_aminoacid), names=seq(1,nrow(sustainedonly),by=1), file.out='sustainedonly_pause_usEPAds_aminoacid.fa', open = "w", nbchar = 60, as.string = FALSE)

write.fasta(sequences=as.list(resolvedonly$EPA), names=seq(1,nrow(resolvedonly),by=1), file.out='resolvedonly_pause_EPA.fa', open = "w", nbchar = 60, as.string = FALSE)
write.fasta(sequences=as.list(resolvedonly$EPA_aminoacid), names=seq(1,nrow(resolvedonly),by=1), file.out='resolvedonly_pause_EPA_aminoacid.fa', open = "w", nbchar = 60, as.string = FALSE)
write.fasta(sequences=as.list(resolvedonly$usEPAds), names=seq(1,nrow(resolvedonly),by=1), file.out='resolvedonly_pause_usEPAds.fa', open = "w", nbchar = 60, as.string = FALSE)
write.fasta(sequences=as.list(resolvedonly$usEPAds_aminoacid), names=seq(1,nrow(resolvedonly),by=1), file.out='resolvedonly_pause_usEPAds_aminoacid.fa', open = "w", nbchar = 60, as.string = FALSE)




#########################################################################################################

pauses$ind<-rep('all other mRNAs',nrow(pauses))
pauses$sust<-rep('all other mRNAs',nrow(pauses))
pauses$res<-rep('all other mRNAs',nrow(pauses))
for(i in 1:nrow(pauses)){
  if((pauses[i,'gene_name'] %in% inducedpause$gene_name)==TRUE){
    pauses[i,'ind']<-'induced'
  }
  if((pauses[i,'gene_name'] %in% sustainedpause$gene_name)==TRUE){
    pauses[i,'sust']<-'sustained'
  }
  if((pauses[i,'gene_name'] %in% resolvedpause$gene_name)==TRUE){
    pauses[i,'res']<-'resolved'
  }
}

pauses$ind<-factor(pauses$ind,levels=c('all other mRNAs','induced'))
pauses<-pauses[order(pauses$ind),]
colourset<-c('grey68','darkmagenta')

f1<-ggplot(pauses,aes(x=percCDS,col=ind))+
  geom_density()+theme_bw()+xlab('% CDS')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), 
                        axis.title = element_text(size = 12),legend.title=element_blank())+
  scale_colour_manual(values=colourset)+coord_trans(xlim=c(0,100))
f1
ggsave(plot=f1,filename='pause_distribution_ind.png',height=3,width=4.5)



pauses$res<-factor(pauses$res,levels=c('all other mRNAs','resolved'))
pauses<-pauses[order(pauses$res),]
colourset<-c('grey68','darkorange')

f1<-ggplot(pauses,aes(x=percCDS,col=res))+
  geom_density()+theme_bw()+xlab('% CDS')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), 
                        axis.title = element_text(size = 12),legend.title=element_blank())+
  scale_colour_manual(values=colourset)+coord_trans(xlim=c(0,100))
f1
ggsave(plot=f1,filename='pause_distribution_res.png',height=3,width=4.5)


pauses$sust<-factor(pauses$sust,levels=c('all other mRNAs','sustained'))
pauses<-pauses[order(pauses$sust),]
colourset<-c('grey68','green4')

f1<-ggplot(pauses,aes(x=percCDS,col=sust))+
  geom_density()+theme_bw()+xlab('% CDS')+
  theme_classic()+theme(legend.text=element_text(size=12),axis.text = element_text(size = 12), 
                        axis.title = element_text(size = 12),legend.title=element_blank())+
  scale_colour_manual(values=colourset)+coord_trans(xlim=c(0,100))
f1
ggsave(plot=f1,filename='pause_distribution_sust.png',height=3,width=4.5)




########################### codons and amino acids motifs at  E/P/A in more detail ######################
library(ggplot2)


##
ind_EPA<-inducedonly[,c(1,13,15:18)]
for(i in 1:nrow(ind_EPA)){
  ind_EPA[i,'E']<-substring(ind_EPA[i,'EPA_aminoacid'],1,1)
  ind_EPA[i,'P']<-substring(ind_EPA[i,'EPA_aminoacid'],2,2)
  ind_EPA[i,'A']<-substring(ind_EPA[i,'EPA_aminoacid'],3,3)
}


singlets<-read.table('\\\\john-doe/gw/Systems/Sarah/CNOTpaper/pausesites_revisions/codon_AAs.csv',sep=',',stringsAsFactors = FALSE,header=TRUE)
#singlets$AA<-a(singlets$AminoAcid)

uniqueAAs<-singlets[!duplicated(singlets$AA_1letter),c(2,3)]

for(i in 1:nrow(uniqueAAs)){
  uniqueAAs[i,'Esite']<-length(which(ind_EPA$E==uniqueAAs[i,'AA_1letter']))
  uniqueAAs[i,'Psite']<-length(which(ind_EPA$P==uniqueAAs[i,'AA_1letter']))
  uniqueAAs[i,'Asite']<-length(which(ind_EPA$A==uniqueAAs[i,'AA_1letter']))
}

#restructuredf
Edf<-uniqueAAs[,c(2,3)]
Edf$site<-rep('E',nrow(Edf))
names(Edf)[2]<-'AA_count'

Pdf<-uniqueAAs[,c(2,4)]
Pdf$site<-rep('P',nrow(Pdf))
names(Pdf)[2]<-'AA_count'

Adf<-uniqueAAs[,c(2,5)]
Adf$site<-rep('A',nrow(Adf))
names(Adf)[2]<-'AA_count'

plotdf<-rbind(Edf,Pdf,Adf)
rm(Edf,Pdf,Adf)

plotdf$site<-factor(plotdf$site,levels=c('E','P','A'))
plotdf<-plotdf[order(plotdf$site),]

# Stacked
ggplot(plotdf, aes(fill=AA_1letter, y=AA_count, x=site)) + geom_bar(position="fill", stat="identity")+
  theme_classic()+scale_y_continuous(expand=c(0,0))+xlab('')+ylab('Proportion')+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12),axis.title = element_text(size = 14))





##################################################
res_EPA<-resolvedonly[,c(1,13,15:18)]
for(i in 1:nrow(res_EPA)){
  res_EPA[i,'E']<-substring(res_EPA[i,'EPA_aminoacid'],1,1)
  res_EPA[i,'P']<-substring(res_EPA[i,'EPA_aminoacid'],2,2)
  res_EPA[i,'A']<-substring(res_EPA[i,'EPA_aminoacid'],3,3)
}


singlets<-read.table('\\\\john-doe/gw/Systems/Sarah/CNOTpaper/pausesites_revisions/codon_AAs.csv',sep=',',stringsAsFactors = FALSE,header=TRUE)
#singlets$AA<-a(singlets$AminoAcid)

uniqueAAs<-singlets[!duplicated(singlets$AA_1letter),c(2,3)]

for(i in 1:nrow(uniqueAAs)){
  uniqueAAs[i,'Esite']<-length(which(res_EPA$E==uniqueAAs[i,'AA_1letter']))
  uniqueAAs[i,'Psite']<-length(which(res_EPA$P==uniqueAAs[i,'AA_1letter']))
  uniqueAAs[i,'Asite']<-length(which(res_EPA$A==uniqueAAs[i,'AA_1letter']))
}

#restructuredf
Edf<-uniqueAAs[,c(2,3)]
Edf$site<-rep('E',nrow(Edf))
names(Edf)[2]<-'AA_count'

Pdf<-uniqueAAs[,c(2,4)]
Pdf$site<-rep('P',nrow(Pdf))
names(Pdf)[2]<-'AA_count'

Adf<-uniqueAAs[,c(2,5)]
Adf$site<-rep('A',nrow(Adf))
names(Adf)[2]<-'AA_count'

plotdf<-rbind(Edf,Pdf,Adf)
rm(Edf,Pdf,Adf)

plotdf$site<-factor(plotdf$site,levels=c('E','P','A'))
plotdf<-plotdf[order(plotdf$site),]

# Stacked
ggplot(plotdf, aes(fill=AA_1letter, y=AA_count, x=site)) + geom_bar(position="fill", stat="identity")+
  theme_classic()+scale_y_continuous(expand=c(0,0))+xlab('')+ylab('Proportion')+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12),axis.title = element_text(size = 14))



##################################################
sust_EPA<-sustainedonly[,c(1,13,15:18)]
for(i in 1:nrow(sust_EPA)){
  sust_EPA[i,'E']<-substring(sust_EPA[i,'EPA_aminoacid'],1,1)
  sust_EPA[i,'P']<-substring(sust_EPA[i,'EPA_aminoacid'],2,2)
  sust_EPA[i,'A']<-substring(sust_EPA[i,'EPA_aminoacid'],3,3)
}


singlets<-read.table('\\\\john-doe/gw/Systems/Sarah/CNOTpaper/pausesites_revisions/codon_AAs.csv',sep=',',stringsAsFactors = FALSE,header=TRUE)
#singlets$AA<-a(singlets$AminoAcid)

uniqueAAs<-singlets[!duplicated(singlets$AA_1letter),c(2,3)]

for(i in 1:nrow(uniqueAAs)){
  uniqueAAs[i,'Esite']<-length(which(sust_EPA$E==uniqueAAs[i,'AA_1letter']))
  uniqueAAs[i,'Psite']<-length(which(sust_EPA$P==uniqueAAs[i,'AA_1letter']))
  uniqueAAs[i,'Asite']<-length(which(sust_EPA$A==uniqueAAs[i,'AA_1letter']))
}

#susttructuredf
Edf<-uniqueAAs[,c(2,3)]
Edf$site<-rep('E',nrow(Edf))
names(Edf)[2]<-'AA_count'

Pdf<-uniqueAAs[,c(2,4)]
Pdf$site<-rep('P',nrow(Pdf))
names(Pdf)[2]<-'AA_count'

Adf<-uniqueAAs[,c(2,5)]
Adf$site<-rep('A',nrow(Adf))
names(Adf)[2]<-'AA_count'

plotdf<-rbind(Edf,Pdf,Adf)
rm(Edf,Pdf,Adf)

plotdf$site<-factor(plotdf$site,levels=c('E','P','A'))
plotdf<-plotdf[order(plotdf$site),]

# Stacked
ggplot(plotdf, aes(fill=AA_1letter, y=AA_count, x=site)) + geom_bar(position="fill", stat="identity")+
  theme_classic()+scale_y_continuous(expand=c(0,0))+xlab('')+ylab('Proportion')+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12),axis.title = element_text(size = 14))


