#codon level
##
ind_EPA<-inducedonly[,c(1,13,15:18)]
for(i in 1:nrow(ind_EPA)){
  ind_EPA[i,'E']<-substring(ind_EPA[i,'EPA'],1,3)
  ind_EPA[i,'P']<-substring(ind_EPA[i,'EPA'],4,6)
  ind_EPA[i,'A']<-substring(ind_EPA[i,'EPA'],7,9)
}

singlets<-read.table('\\\\john-doe/gw/Systems/Sarah/CNOTpaper/pausesites_revisions/codon_AAs.csv',sep=',',stringsAsFactors = FALSE,header=TRUE)
singlets$third_nt<-sapply(singlets$Codon,function(x) substring(x,3,3))

for(i in 1:nrow(singlets)){
  singlets[i,'Esite']<-length(which(ind_EPA$E==singlets[i,'Codon']))
  singlets[i,'Psite']<-length(which(ind_EPA$P==singlets[i,'Codon']))
  singlets[i,'Asite']<-length(which(ind_EPA$A==singlets[i,'Codon']))
}

#restructuredf
Edf<-singlets[,c(1,4,5,3)]
Edf$site<-rep('E',nrow(Edf))
names(Edf)[3]<-'codon_count'

Pdf<-singlets[,c(1,4,6,3)]
Pdf$site<-rep('P',nrow(Pdf))
names(Pdf)[3]<-'codon_count'

Adf<-singlets[,c(1,4,7,3)]
Adf$site<-rep('A',nrow(Adf))
names(Adf)[3]<-'codon_count'

plotdf<-rbind(Edf,Pdf,Adf)
rm(Edf,Pdf,Adf)

plotdf$site<-factor(plotdf$site,levels=c('E','P','A'))
plotdf<-plotdf[order(plotdf$site),]
#plotdf$AA_1letter<-factor(plotdf$AA_1letter,levels=unique(plotdf$AA_1letter))
#plotdf<-plotdf[order(plotdf$AA_1letter),]

# Stacked
ggplot(plotdf, aes(fill=Codon, y=codon_count, x=site)) + geom_bar(position="fill", stat="identity")+
  theme_classic()+scale_y_continuous(expand=c(0,0))+xlab('')+ylab('Proportion')+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12),axis.title = element_text(size = 14))

ggplot(plotdf, aes(fill=third_nt, y=codon_count, x=site)) + geom_bar(position="fill", stat="identity")+
  theme_classic()+scale_y_continuous(expand=c(0,0))+xlab('')+ylab('Proportion')+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12),axis.title = element_text(size = 14))




#codon level
##
res_EPA<-resolvedonly[,c(1,13,15:18)]
for(i in 1:nrow(res_EPA)){
  res_EPA[i,'E']<-substring(res_EPA[i,'EPA'],1,3)
  res_EPA[i,'P']<-substring(res_EPA[i,'EPA'],4,6)
  res_EPA[i,'A']<-substring(res_EPA[i,'EPA'],7,9)
}

singlets<-read.table('\\\\john-doe/gw/Systems/Sarah/CNOTpaper/pausesites_revisions/codon_AAs.csv',sep=',',stringsAsFactors = FALSE,header=TRUE)
singlets$third_nt<-sapply(singlets$Codon,function(x) substring(x,3,3))

for(i in 1:nrow(singlets)){
  singlets[i,'Esite']<-length(which(res_EPA$E==singlets[i,'Codon']))
  singlets[i,'Psite']<-length(which(res_EPA$P==singlets[i,'Codon']))
  singlets[i,'Asite']<-length(which(res_EPA$A==singlets[i,'Codon']))
}

#restructuredf
Edf<-singlets[,c(1,4,5,3)]
Edf$site<-rep('E',nrow(Edf))
names(Edf)[3]<-'codon_count'

Pdf<-singlets[,c(1,4,6,3)]
Pdf$site<-rep('P',nrow(Pdf))
names(Pdf)[3]<-'codon_count'

Adf<-singlets[,c(1,4,7,3)]
Adf$site<-rep('A',nrow(Adf))
names(Adf)[3]<-'codon_count'

plotdf<-rbind(Edf,Pdf,Adf)
rm(Edf,Pdf,Adf)

plotdf$site<-factor(plotdf$site,levels=c('E','P','A'))
plotdf<-plotdf[order(plotdf$site),]
#plotdf$AA_1letter<-factor(plotdf$AA_1letter,levels=unique(plotdf$AA_1letter))
#plotdf<-plotdf[order(plotdf$AA_1letter),]

# Stacked
ggplot(plotdf, aes(fill=Codon, y=codon_count, x=site)) + geom_bar(position="fill", stat="identity")+
  theme_classic()+scale_y_continuous(expand=c(0,0))+xlab('')+ylab('Proportion')+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12),axis.title = element_text(size = 14))

ggplot(plotdf, aes(fill=third_nt, y=codon_count, x=site)) + geom_bar(position="fill", stat="identity")+
  theme_classic()+scale_y_continuous(expand=c(0,0))+xlab('')+ylab('Proportion')+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12),axis.title = element_text(size = 14))




#codon level
##
sust_EPA<-sustainedonly[,c(1,13,15:18)]
for(i in 1:nrow(sust_EPA)){
  sust_EPA[i,'E']<-substring(sust_EPA[i,'EPA'],1,3)
  sust_EPA[i,'P']<-substring(sust_EPA[i,'EPA'],4,6)
  sust_EPA[i,'A']<-substring(sust_EPA[i,'EPA'],7,9)
}

singlets<-read.table('\\\\john-doe/gw/Systems/Sarah/CNOTpaper/pausesites_revisions/codon_AAs.csv',sep=',',stringsAsFactors = FALSE,header=TRUE)
singlets$third_nt<-sapply(singlets$Codon,function(x) substring(x,3,3))

for(i in 1:nrow(singlets)){
  singlets[i,'Esite']<-length(which(sust_EPA$E==singlets[i,'Codon']))
  singlets[i,'Psite']<-length(which(sust_EPA$P==singlets[i,'Codon']))
  singlets[i,'Asite']<-length(which(sust_EPA$A==singlets[i,'Codon']))
}

#susttructuredf
Edf<-singlets[,c(1,4,5,3)]
Edf$site<-rep('E',nrow(Edf))
names(Edf)[3]<-'codon_count'

Pdf<-singlets[,c(1,4,6,3)]
Pdf$site<-rep('P',nrow(Pdf))
names(Pdf)[3]<-'codon_count'

Adf<-singlets[,c(1,4,7,3)]
Adf$site<-rep('A',nrow(Adf))
names(Adf)[3]<-'codon_count'

plotdf<-rbind(Edf,Pdf,Adf)
rm(Edf,Pdf,Adf)

plotdf$site<-factor(plotdf$site,levels=c('E','P','A'))
plotdf<-plotdf[order(plotdf$site),]
#plotdf$AA_1letter<-factor(plotdf$AA_1letter,levels=unique(plotdf$AA_1letter))
#plotdf<-plotdf[order(plotdf$AA_1letter),]

# Stacked
ggplot(plotdf, aes(fill=Codon, y=codon_count, x=site)) + geom_bar(position="fill", stat="identity")+
  theme_classic()+scale_y_continuous(expand=c(0,0))+xlab('')+ylab('Proportion')+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12),axis.title = element_text(size = 14))

ggplot(plotdf, aes(fill=third_nt, y=codon_count, x=site)) + geom_bar(position="fill", stat="identity")+
  theme_classic()+scale_y_continuous(expand=c(0,0))+xlab('')+ylab('Proportion')+
  theme(legend.text=element_text(size=12),axis.text = element_text(size = 12),axis.title = element_text(size = 14))

