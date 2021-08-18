library(seqinr)
library(Biostrings)
library(zoo)
library(ggplot2)
library(dplyr)
library(stringr)
library(gplots)
library(matrixStats)




############################## get data ####################################
notRPGs<-read.delim('\\\\john-doe/gw/Systems/Sarah/CNOTdata/finaldraftdata_171019/mRNAs_noRPGsnoHistones.txt',header=TRUE,stringsAsFactors = FALSE,sep='\t')
RPF27to31groups<-read.delim('\\\\john-doe/gw/Systems/Sarah/CNOTdata/finaldraftdata_171019/RPFCDSfilt_27to31_withgroups.txt',header=TRUE,stringsAsFactors = FALSE)
nrow(notRPGs)
nrow(RPF27to31groups)
RPF27to31groups<-subset(RPF27to31groups, (Gene_ID %in% notRPGs$Gene_ID)==TRUE)
nrow(RPF27to31groups)

test<-merge(RPF27to31groups,humanDF,by='Gene_ID')
nrow(humanDF)
nrow(RPF27to31groups)
nrow(test) #8613
test<-subset(test,startsWith(Gene_ID, 'RP11-')==FALSE)
nrow(test) #8542

test<-subset(test, Group!='Totaldown RPFup')
test<-subset(test, Group!='Totalup RPFdown')
test<-subset(test, Group!='not significant')

test$Group<-sub('very ','',test$Group)
#test$Group<-factor(test$Group,levels=c('not significant','RNA up','RNA down','Both up','Both down','RPF up','RPF down'))
#test<-test[order(test$Group),]

test<-subset(test, nchar(CDS_sequence)>300) 
nrow(test) #5966

for(i in 1:nrow(test)){
  midpoint<-nchar(test[i,'CDS_sequence'])/2
  if((midpoint%%3)!=0){
    midpoint<-midpoint-(midpoint%%3)
  }
  test[i,'CDS_firsthalf']<-substring(test[i,'CDS_sequence'],1,midpoint)
  test[i,'CDS_lasthalf']<-substring(test[i,'CDS_sequence'],midpoint+1,nchar(test[i,'CDS_sequence']))
  test[i,'CDS_first150']<-substring(test[i,'CDS_sequence'],1,150)
  test[i,'CDS_last150']<-substring(test[i,'CDS_sequence'],nchar(test[i,'CDS_sequence'])-149,nchar(test[i,'CDS_sequence']))
}




All_mRNAs<-test
All_mRNAs$Group<-rep('All mRNAs',nrow(All_mRNAs))
RNA_up<-subset(test,Group=='RNA up')
RNA_down<-subset(test,Group=='RNA down')
Both_up<-subset(test,Group=='Both up')
Both_down<-subset(test,Group=='Both down')
RPF_up<-subset(test,Group=='RPF up')
RPF_down<-subset(test,Group=='RPF down')
non_sig<-subset(test,Group=='not significant')



################################### try codon frequency plots ######################################


codon_plot<-function(plotname,datasets,samplenames,colours,region){
  codon_AA<-read.table('\\\\john-doe/gw/Systems/Sarah/CNOTdataAnalysis/Codons/codon_AA.csv',sep=',',stringsAsFactors = FALSE,header=TRUE)
  codon_AA$Codon<-tolower(codon_AA$Codon)
  codonset<-as.vector(data.frame(uco(s2c('ATGATCGAACAATGGCCCGGGGCCTGA'),frame=0,as.data.frame=FALSE,index='eff'),stringsAsFactors = FALSE)$Var1)
  excludelist<-c('atg','tag','taa','tga')
  codonset<-codonset[(codonset %in% excludelist)==FALSE]
  #per mRNA codon usage
  eff_data<-list()
  for(k in 1:length(datasets)){
    
    UCO_alltypes<-lapply(datasets[[k]][,region], function(x) uco(s2c(x), frame = 0, as.data.frame = FALSE, index='eff', NA.rscu = NA))
    UCOalltypesDF_eff<-do.call(rbind,UCO_alltypes) 
    UCOalltypesDF_eff<-data.frame(UCOalltypesDF_eff,stringsAsFactors = FALSE)
    row.names(UCOalltypesDF_eff)<-datasets[[k]]$Gene_ID 
    UCOalltypesDF_eff<-UCOalltypesDF_eff[,(names(UCOalltypesDF_eff) %in% excludelist)==FALSE]
    print(dim(UCOalltypesDF_eff))
    
    eff_data[[k]]<-as.data.frame(UCOalltypesDF_eff)
    
  }
  print('eff done')
  
  freq_data<-eff_data
  for(z in 1:length(freq_data)){
    for(zz in 1:nrow(freq_data[[z]])){
      freq_data[[z]][zz,]<-freq_data[[z]][zz,]/sum(freq_data[[z]][zz,])
    }
  }
  print('freq done')
  
  toplot_data<-NULL
  for(y in 1:length(freq_data)){
    toplot_df<-data.frame(matrix(nrow=60,ncol=4,0),stringsAsFactors = FALSE)
    names(toplot_df)<-c('codon','Group','codon_freq_mean','codon_freq_median')
    toplot_df$codon<-colnames(freq_data[[y]])
    toplot_df$Group<-rep(samplenames[y],nrow(toplot_df))
    toplot_df$codon_freq_mean<-unlist(apply(freq_data[[y]],2,function(x) mean(x)))
    toplot_df$codon_freq_median<-unlist(apply(freq_data[[y]],2,function(x) median(x)))
    toplot_data<-rbind(toplot_data,toplot_df)
  }
  
  
  toplot_data$third_nt
  for(m in 1:nrow(toplot_data)){
    if(substring(toplot_data[m,'codon'],3,3)=='a'){
      toplot_data[m,'third_nt']<-'A'
    }
    if(substring(toplot_data[m,'codon'],3,3)=='c'){
      toplot_data[m,'third_nt']<-'C'
    }
    if(substring(toplot_data[m,'codon'],3,3)=='g'){
      toplot_data[m,'third_nt']<-'G'
    }
    if(substring(toplot_data[m,'codon'],3,3)=='t'){
      toplot_data[m,'third_nt']<-'T'
    }
  }
  #return(toplot_data)

  f1<-ggplot(toplot_data,aes(x=codon,y=codon_freq_mean,fill=Group,colour=third_nt))+geom_bar(stat='identity',position='dodge')+geom_text(aes(label=codon),size=4)+theme_bw()+scale_colour_manual(values=c('magenta','deepskyblue','limegreen','purple'))+ylab('codon freq')+ggtitle(paste0('mean codon frequency (per mRNA)'))+scale_fill_manual(values=colours)
  f1<-f1+theme(legend.text=element_text(size=18),axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16), axis.title = element_text(size = 18), panel.grid.major = element_line(color = "gray",linetype = "dotted"),plot.title=element_blank(),legend.title=element_blank())
  ggsave(paste0(region,'_',plotname,'_mean.png'),height=6,width=24)
  
  f1<-ggplot(toplot_data,aes(x=codon,y=codon_freq_median,fill=Group,colour=third_nt))+geom_bar(stat='identity',position='dodge')+geom_text(aes(label=codon),size=4)+theme_bw()+scale_colour_manual(values=c('magenta','deepskyblue','limegreen','purple'))+ylab('codon freq')+ggtitle(paste0('median codon frequency (per mRNA)'))+scale_fill_manual(values=colours)
  f1<-f1+theme(legend.text=element_text(size=18),axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16), axis.title = element_text(size = 18), panel.grid.major = element_line(color = "gray",linetype = "dotted"),plot.title=element_blank(),legend.title=element_blank())
  ggsave(paste0(region,'_',plotname,'_median.png'),height=6,width=24)
  
}



setwd("//john-doe/gw/Systems/Sarah/CNOTdataAnalysis/Codons/170520")
regionselection<-c('CDS_sequence','CDS_firsthalf','CDS_lasthalf','CDS_first150','CDS_last150')
regiontouse<-regionselection[5]

group_datasets<-list(RPF_up,RNA_up,non_sig)
group_names<-c('RPF up','RNA up','not significant')
colourstouse<-c('grey48','darkorange','green4')
plotdata<-codon_plot('RPFupvRNAup',group_datasets,group_names,colourstouse,regiontouse)

group_datasets<-list(RPF_up,RPF_down,non_sig)
group_names<-c('RPF up','RPF down','not significant')
colourstouse<-c('grey48','red3','green4')
plotdata<-codon_plot('RPFupvRPFdown',group_datasets,group_names,colourstouse,regiontouse)

group_datasets<-list(Both_up,Both_down,non_sig)
group_names<-c('Both up','Both down','not significant')
colourstouse<-c('gold2','darkmagenta','grey48')
plotdata<-codon_plot('BothupvBothdown',group_datasets,group_names,colourstouse,regiontouse)

group_datasets<-list(Both_up,RPF_up,non_sig)
group_names<-c('Both up','RPF up','not significant')
colourstouse<-c('darkmagenta','grey48','green4')
plotdata<-codon_plot('BothupvRPFup',group_datasets,group_names,colourstouse,regiontouse)

group_datasets<-list(Both_down,RPF_down,non_sig)
group_names<-c('Both down','RPF down','not significant')
colourstouse<-c('gold2','grey48','red3')
plotdata<-codon_plot('BothdownvRPFdown',group_datasets,group_names,colourstouse,regiontouse)

group_datasets<-list(RNA_up,RNA_down,non_sig)
group_names<-c('RNA up','RNA down','not significant')
colourstouse<-c('grey48','steelblue4','darkorange')
plotdata<-codon_plot('RNAupvRNAdown',group_datasets,group_names,colourstouse,regiontouse)



############################################## examine AAA | GAA | CAA in more detail ######################################################
regexprtouse<-c('(AAA|GAA){2,}')
test$AAA_GAA_stretches<-lapply(test$CDS_sequence,function(x) length(which(((str_locate_all(x,regexprtouse)[[1]][,1])%%3)!=1)))
test$AAA_GAA_stretches<-as.numeric(test$AAA_GAA_stretches)

regexprtouse<-c('GAAAAA')
test$AAA_GAA_stretches<-lapply(test$CDS_sequence,function(x) length(which(((str_locate_all(x,regexprtouse)[[1]][,1])%%3)!=1)))
test$AAA_GAA_stretches<-as.numeric(test$AAA_GAA_stretches)

regexprtouse<-c('(GAAAAA|CAAAAA)')
test$AAA_GAA_CAA_stretches<-lapply(test$CDS_sequence,function(x) length(which(((str_locate_all(x,regexprtouse)[[1]][,1])%%3)!=1)))
test$AAA_GAA_CAA_stretches<-as.numeric(test$AAA_GAA_CAA_stretches)


regexprtouse<-c('(AAA|GAA|CAA){2,}')
test$AAA_GAA_CAA_stretches<-lapply(test$CDS_sequence,function(x) length(which(((str_locate_all(x,regexprtouse)[[1]][,1])%%3)!=1)))
test$AAA_GAA_CAA_stretches<-as.numeric(test$AAA_GAA_CAA_stretches)

regexprtouse<-c('(AAA){2,}')
test$AAA_stretches<-lapply(test$CDS_sequence,function(x) length(which(((str_locate_all(x,regexprtouse)[[1]][,1])%%3)!=1)))
test$AAA_stretches<-as.numeric(test$AAA_stretches)

regexprtouse<-c('(GAA){2,}')
test$GAA_stretches<-lapply(test$CDS_sequence,function(x) length(which(((str_locate_all(x,regexprtouse)[[1]][,1])%%3)!=1)))
test$GAA_stretches<-as.numeric(test$GAA_stretches)


groupstouse<-c('not significant','RNA up','RNA down','Both up','Both down','RPF up','RPF down')
for(i in 1:length(groupstouse)){
  print(groupstouse[i])
  #print(summary(subset(test,Group==groupstouse[i])$AAA_stretches))
  #print(summary(subset(test,Group==groupstouse[i])$GAA_stretches))
  print(summary(subset(test,Group==groupstouse[i])$AAA_GAA_stretches))
  print(summary(subset(test,Group==groupstouse[i])$AAA_GAA_CAA_stretches))
}



#plot data
groupstouse<-c('not significant','RNA up','RNA down','Both up','Both down','RPF up','RPF down')
binned_data<-NULL
for(i in 1:length(groupstouse)){
  ss<-subset(test,Group==groupstouse[i])
  print(nrow(ss))
  
  xx<-as.data.frame(table(cut(ss$AAA_GAA_CAA_stretches, breaks=c(0,1,2,5,10,20,50))))
  xx$Group<-rep(groupstouse[i],nrow(xx))
  xx$Var1<-factor(xx$Var1,levels=c('(20,50]','(10,20]','(5,10]','(2,5]','(1,2]','(0,1]'))
  xx<-xx[order(xx$Var1),]
  
  binned_data<-rbind(binned_data,xx)
  rm(xx)
  rm(ss)
}
binned_data$Group<-factor(binned_data$Group,levels=c('not significant','RNA up','RNA down','Both up','Both down','RPF up','RPF down'))
ggplot(binned_data,aes(x=Group,y=Freq,fill=Var1))+geom_bar(stat='identity',position='fill')+theme_classic()


