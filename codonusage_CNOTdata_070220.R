#to each other / All
#to mitotic
#to pattern

#############################
#source("http://bioconductor.org/biocLite.R")
#biocLite("org.Hs.eg.db")
#biocLite("GO.db")
library(org.Hs.eg.db)
library(GO.db)

go_id = GOID( GOTERM[ Term(GOTERM) == "mitotic cell cycle"])
allegs = get(go_id, org.Hs.egGO2ALLEGS)

genes = unlist(mget(allegs,org.Hs.egSYMBOL))
new2_MitoticCellCycle<-data.frame(unique(genes),stringsAsFactors = FALSE)
names(new2_MitoticCellCycle)[1]<-'Gene_ID'
#write.table(new2_MitoticCellCycle,file='new2_MitoticCellCycle_Mm.txt',quote=FALSE,row.names=FALSE)

go_id = GOID( GOTERM[ Term(GOTERM) == "pattern specification process"])
allegs = get(go_id, org.Hs.egGO2ALLEGS)

genes = unlist(mget(allegs,org.Hs.egSYMBOL))
new2_PatternSpecificationProcess<-data.frame(unique(genes),stringsAsFactors = FALSE)
names(new2_PatternSpecificationProcess)[1]<-'Gene_ID'

new2_MitoticCellCycle<-subset(humanDF[,c(5,11)],(Gene_ID %in% new2_MitoticCellCycle$Gene_ID)==TRUE)
new2_PatternSpecificationProcess<-subset(humanDF[,c(5,11)],(Gene_ID %in% new2_PatternSpecificationProcess$Gene_ID)==TRUE)

nrow(new2_MitoticCellCycle)
nrow(new2_PatternSpecificationProcess)

All_mRNAs<-test
All_mRNAs$Group<-rep('All mRNAs',nrow(All_mRNAs))
RNA_up<-subset(test,Group=='RNA up')
RNA_down<-subset(test,Group=='RNA down')
Both_up<-subset(test,Group=='Both up')
Both_down<-subset(test,Group=='Both down')
RPF_up<-subset(test,Group=='RPF up')
RPF_down<-subset(test,Group=='RPF down')
non_sig<-subset(test,Group=='very not significant')

######################################################

setwd('~/CNOTdata/finaldraftdata_171019/Codons_070220')

type<-'mostabundant'
datasets_list<-list(RNA_up,RNA_down,RPF_up,RPF_down,Both_up,Both_down,All_mRNAs,non_sig)
lapply(datasets_list,function(x) nrow(x))
samplenames<-list('RNA up','RNA down','RPF up','RPF down','Both up','Both down','All mRNAs','not significant')
for(j in 1:length(datasets_list)){
  for(k in 1:length(datasets_list)){
    print(samplenames[[j]])
    print(samplenames[[k]])
    codon_plot(paste(samplenames[[j]],'_',samplenames[[k]],collapse=''),list(datasets_list[[j]],datasets_list[[k]]),c(samplenames[[j]],samplenames[[k]]),type)
  }
}

#added 230421
new2_MitoticCellCycle<-subset(humanDF[,c(5,11)],(Gene_ID %in% new2_MitoticCellCycle$Gene_ID)==TRUE)
new2_PatternSpecificationProcess<-subset(humanDF[,c(5,11)],(Gene_ID %in% new2_PatternSpecificationProcess$Gene_ID)==TRUE)

SAup0<-subset(humanDF[,c(5,11)],(Gene_ID %in% SAup0$Gene_ID)==TRUE)
SAdown0<-subset(humanDF[,c(5,11)],(Gene_ID %in% SAdown0$Gene_ID)==TRUE)

SAup05<-subset(humanDF[,c(5,11)],(Gene_ID %in% SAup05$Gene_ID)==TRUE)
SAdown05<-subset(humanDF[,c(5,11)],(Gene_ID %in% SAdown05$Gene_ID)==TRUE)

type<-'mostabundant'
datasets_list<-list(new2_PatternSpecificationProcess,new2_MitoticCellCycle,SAup0,SAdown0,SAup05,SAdown05)
lapply(datasets_list,function(x) nrow(x))
samplenames<-list('GO: Pattern Specification Process', 'GO: Mitotic Cell Cycle')
samplenames<-list('PatternSpecificationProcess', 'MitoticCellCycle','SAup0','SAdown0','SAup05','SAdown05')
for(j in 1:length(datasets_list)){
  for(k in 1:length(datasets_list)){
    print(samplenames[[j]])
    print(samplenames[[k]])
    codon_plot(paste(samplenames[[j]],'_',samplenames[[k]],collapse=''),list(datasets_list[[j]],datasets_list[[k]]),c(samplenames[[j]],samplenames[[k]]),type)
  }
}

#########################################################################
codon_plot<-function(plotname,datasets,names,type){
  codon_AA<-read.table('//john-doe/gw/Systems/Sarah/CNOTdataAnalysis/Codons/codon_AA.csv',sep=',',stringsAsFactors = FALSE,header=TRUE)
  codon_AA$Codon<-tolower(codon_AA$Codon)
  codonset<-as.vector(data.frame(uco(s2c('ATGATCGAACAATGGCCCGGGGCCTGA'),frame=0,as.data.frame=FALSE,index='eff'),stringsAsFactors = FALSE)$Var1)
  excludelist<-c('atg','tag','taa','tga','tgg')
  codonset<-codonset[(codonset %in% excludelist)==FALSE]
  #per mRNA codon usage
  eff_data<-list()
  for(k in 1:length(datasets)){
    
    UCO_alltypes<-lapply(datasets[[k]]$CDS_sequence, function(x) uco(s2c(x), frame = 0, as.data.frame = FALSE, index='eff', NA.rscu = NA))
    UCOalltypesDF_eff<-do.call(rbind,UCO_alltypes) 
    UCOalltypesDF_eff<-data.frame(UCOalltypesDF_eff,stringsAsFactors = FALSE)
    row.names(UCOalltypesDF_eff)<-datasets[[k]]$Gene_ID 
    UCOalltypesDF_eff<-UCOalltypesDF_eff[,(names(UCOalltypesDF_eff) %in% excludelist)==FALSE]
    #names(UCOalltypesDF_eff)<-codonset
    print(dim(UCOalltypesDF_eff))
    #return(UCOalltypesDF_eff)
    #names(UCOalltypesDF_eff)<-codonset
    #print('eff')
    
    
    eff_data[[k]]<-as.data.frame(UCOalltypesDF_eff)
    
  }
  
  uniqueAA<-unique(codon_AA$AminoAcid)
  uniqueAA<-c("Ile" ,"Leu","Val","Phe","Cys","Ala","Gly","Pro","Thr","Ser","Tyr","Gln","Asn","His","Glu","Asp","Lys","Arg")
  eff_output<-list()
  for(n in 1:length(eff_data)){
    all_means<-NULL #add in group label // column as groups // rownames as codons 
    for(AAcid in uniqueAA){
      #print(AAcid)
      means<-NULL
      codons<-subset(codon_AA,AminoAcid==AAcid)$Codon
      codon_table<-eff_data[[n]][,codons]
      codon_output<-apply(codon_table,1, function(x) x/(sum(x))) 
      #codon_output2<-do.call(cbind,t(codon_output))
      codon_output<-t(codon_output)
      means<-data.frame(apply(codon_output,2,function(x) mean(x,na.rm=TRUE)),stringsAsFactors = FALSE) #remove NaN from means # need to think about zeros 
      #for(c in 1:length(codons)){
      #  means<-rbind(means,c(codons[c],mean(codon_output2[,c])))
      #}
      names(means)<-c(paste0('SCU_',n))
      means$codon<-rownames(means)
      all_means<-rbind(all_means,means)
    }
    eff_output[[n]]<-all_means
  }
  
  toplot_data<-merge(eff_output[[1]],eff_output[[2]],by='codon')
  #means gives synonymous codon usage means
  
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
      toplot_data[m,'third_nt']<-'U'
    }
  }
  toplot_data[,'codon']<-gsub('t','u',toplot_data[,'codon'])
  #put up usage on y-axis
  f1<-ggplot(toplot_data,aes(x=SCU_2,y=SCU_1,fill=codon,colour=third_nt))+geom_point(alpha=0.5,size=0.75)+geom_text(aes(label=codon),size=4)+theme_bw()+scale_colour_manual(values=c('magenta','deepskyblue','limegreen','purple'))+xlab(names[2])+ylab(names[1])+ggtitle(paste0('synonymous codon usage (per mRNA) ',type))+coord_trans(limx=c(0,1),limy=c(0,1))+theme(legend.position = "none") #need to remove legend
  f1<-f1+theme(legend.text=element_text(size=18),axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16), axis.title = element_text(size = 18), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_blank(),legend.title=element_blank())
  ggsave(paste0(plotname,'_',type,'.png'),height=6,width=6)
  
}

#############################################################################
library(gplots)
library("RColorBrewer")

codon_data_for_heatmap<-function(datasets,names){
  codon_AA<-read.table('//john-doe/gw/Systems/Sarah/CNOTdataAnalysis/Codons/codon_AA.csv',sep=',',stringsAsFactors = FALSE,header=TRUE)
  codon_AA$Codon<-tolower(codon_AA$Codon)
  codonset<-as.vector(data.frame(uco(s2c('ATGATCGAACAATGGCCCGGGGCCTGA'),frame=0,as.data.frame=FALSE,index='eff'),stringsAsFactors = FALSE)$Var1)
  excludelist<-c('atg','tag','taa','tga','tgg')
  #codonset<-codonset[(codonset %in% excludelist)==FALSE]
  #per mRNA codon usage
  eff_data<-list()
  for(k in 1:length(datasets)){
    
    UCO_alltypes<-lapply(datasets[[k]]$CDS_sequence, function(x) uco(s2c(x), frame = 0, as.data.frame = FALSE, index='eff', NA.rscu = NA))
    UCOalltypesDF_eff<-do.call(rbind,UCO_alltypes) 
    UCOalltypesDF_eff<-data.frame(UCOalltypesDF_eff,stringsAsFactors = FALSE)
    row.names(UCOalltypesDF_eff)<-datasets[[k]]$Gene_ID 
    UCOalltypesDF_eff<-UCOalltypesDF_eff[,(names(UCOalltypesDF_eff) %in% excludelist)==FALSE]
    #names(UCOalltypesDF_eff)<-codonset
    print(dim(UCOalltypesDF_eff))
    #return(UCOalltypesDF_eff)
    
    eff_data[[k]]<-as.data.frame(UCOalltypesDF_eff)
    
  }
  
  uniqueAA<-unique(codon_AA$AminoAcid)
  uniqueAA<-c("Ile" ,"Leu","Val","Phe","Cys","Ala","Gly","Pro","Thr","Ser","Tyr","Gln","Asn","His","Glu","Asp","Lys","Arg")
  eff_output<-list()
  for(n in 1:length(eff_data)){
    all_means<-NULL #add in group label // column as groups // rownames as codons 
    for(AAcid in uniqueAA){
      #print(AAcid)
      means<-NULL
      codons<-subset(codon_AA,AminoAcid==AAcid)$Codon
      codon_table<-eff_data[[n]][,codons]
      codon_output<-apply(codon_table,1, function(x) x/(sum(x))) 
      #codon_output2<-do.call(cbind,t(codon_output))
      codon_output<-t(codon_output)
      means<-data.frame(apply(codon_output,2,function(x) mean(x,na.rm=TRUE)),stringsAsFactors = FALSE) #remove NaN from means # need to think about zeros 
      #for(c in 1:length(codons)){
      #  means<-rbind(means,c(codons[c],mean(codon_output2[,c])))
      #}
      names(means)<-c(paste0('SCU_',names[n]))
      means$codon<-rownames(means)
      all_means<-rbind(all_means,means)
    }
    eff_output[[n]]<-all_means
  }
  
  toplot_data<-merge(eff_output[[1]],eff_output[[2]],by='codon')
  for(p in 3:length(eff_output)){
    toplot_data<-merge(toplot_data,eff_output[[p]],by='codon')
  }
  #means gives synonymous codon usage means
  toplot_data$codon<-gsub('t','u',toplot_data$codon)
  toplot_data$third_nt
  for(m in 1:nrow(toplot_data)){
    if(substring(toplot_data[m,'codon'],3,3)=='a'){
      toplot_data[m,'third_nt']<-'magenta'
    }
    if(substring(toplot_data[m,'codon'],3,3)=='c'){
      toplot_data[m,'third_nt']<-'deepskyblue'
    }
    if(substring(toplot_data[m,'codon'],3,3)=='g'){
      toplot_data[m,'third_nt']<-'limegreen'
    }  
    if(substring(toplot_data[m,'codon'],3,3)=='u'){
      toplot_data[m,'third_nt']<-'purple'
    }
  }
  
  return(toplot_data)
  
}


datasets_list<-list(RNA_up,RNA_down,RPF_up,RPF_down,Both_up,Both_down,All_mRNAs,non_sig)
samplenames<-list('RNA_up','RNA_down','RPF_up','RPF_down','Both_up','Both_down','All_mRNAs','non_sig')
dataforheatmap<-codon_data_for_heatmap(datasets_list,samplenames)
rownames(dataforheatmap)<-dataforheatmap$codon
heatmap_matrix<-as.matrix(dataforheatmap[,c(2:7)])

datasets_list<-list(RNA_up,RNA_down,RPF_up,RPF_down,Both_up,Both_down,non_sig)
samplenames<-list('RNA up','RNA down','RPF up','RPF down','Both up','Both down','very not significant')
dataforheatmap<-codon_data_for_heatmap(datasets_list,samplenames)
rownames(dataforheatmap)<-dataforheatmap$codon
heatmap_matrix<-as.matrix(dataforheatmap[,c(2:8)])


######

#apply z-score to rows before clustering? Keep column order as in previous plots?
datasets_list<-list(RNA_up,RNA_down,RPF_up,RPF_down,Both_up,Both_down,non_sig)
samplenames<-list('RNA up','RNA down','RPF up','RPF down','Both up','Both down','very not significant')
notthese<-c('aug','uag','uaa','uga','ugg')
dataforheatmap<-codon_data_for_heatmap(datasets_list,samplenames)
dataforheatmap<-subset(dataforheatmap,(dataforheatmap$codon %in% notthese)==FALSE)
for(k in 1:nrow(dataforheatmap)){
  dataforheatmap[k,2:8]<-scale(as.numeric(dataforheatmap[k,2:8]))
}
rownames(dataforheatmap)<-dataforheatmap$codon
heatmap_matrix<-as.matrix(dataforheatmap[,c(2:8)])


#zscore already applied
tiff(file = "zscore_preapplied_codon_centroid_noscale.tiff",width=350,height=850)
heatmap.2(as.matrix(heatmap_df),col=colorspace::diverge_hsv(15),Rowv=TRUE,Colv=FALSE,cexRow=1,dendrogram='row',scale="none", RowSideColors = dataforheatmap$third_nt,trace="none",density.info='none', hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(20,12))
dev.off()
#zscore already applied
tiff(file = "zscore_preapplied_codon_centroid_noscale_v2.tiff",width=350,height=850)
heatmap.2(as.matrix(heatmap_df),col=colorspace::diverge_hsv(15),Rowv=TRUE,Colv=TRUE,cexRow=1,dendrogram='both',RowSideColors = dataforheatmap$third_nt,scale="none", trace="none",density.info='none', hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(20,12))
dev.off()



######
tiff(file = "SCU_groups.tiff",width=600,height=1000)
heatmap.2(heatmap_matrix,col=brewer.pal(11,"RdBu"),scale="row", trace="none",density.info='none')
dev.off()

tiff(file = "SCU_groups_centroid_rowscale_attempt2.tiff",width=500,height=1200)
heatmap.2(heatmap_matrix,col=colorspace::diverge_hsv(15),cexRow=1.5,scale="row", trace="none",density.info='none',dendrogram='row',RowSideColors = dataforheatmap$third_nt, hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(18,12))
dev.off()

tiff(file = "SCU_groups_centroid_noscale.tiff",width=600,height=1200)
heatmap.2(heatmap_matrix,col=colorspace::diverge_hsv(15),cexRow=1,scale="none", trace="none",density.info='none',RowSideColors = dataforheatmap$third_nt, hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(18,12))
dev.off()


tiff(file = "SCU_groups_centroid_colscale.tiff",width=600,height=1200)
heatmap.2(heatmap_matrix,col=colorspace::diverge_hsv(15),cexRow=1,scale="col", trace="none",density.info='none',RowSideColors = dataforheatmap$third_nt, hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(18,12))
dev.off()


########################## look at on an amino acid basis ######################################

AAfreq_data_for_heatmap<-function(datasets,names){
  codon_AA<-read.table('~/CNOTdata/finaldraftdata_171019/Codons/codon_AA.csv',sep=',',stringsAsFactors = FALSE,header=TRUE)
  codon_AA$Codon<-tolower(codon_AA$Codon)
  codonset<-data.frame(uco(s2c('ATGATCGAACAATGGCCCGGGGCCTGA'),frame=0,as.data.frame=FALSE,index='eff'),stringsAsFactors = FALSE)$Var1
  
  #per mRNA codon usage
  eff_data<-list()
  for(k in 1:length(datasets)){
    
    UCO_alltypes<-lapply(datasets[[k]]$CDS_sequence, function(x) uco(s2c(x), frame = 0, as.data.frame = FALSE, index='eff', NA.rscu = NA))
    UCOalltypesDF_eff<-do.call(rbind,UCO_alltypes) 
    row.names(UCOalltypesDF_eff)<-datasets[[k]]$Gene_ID 
    names(UCOalltypesDF_eff)<-codonset
    print('eff')
    
    
    eff_data[[k]]<-as.data.frame(UCOalltypesDF_eff)
    
  }
  
  
  eff_output<-list()
  for(n in 1:length(eff_data)){
    all_freqs<-data.frame(matrix(nrow=nrow(codon_AA),ncol=2,0),stringsAsFactors = FALSE)
    names(all_freqs)<-c('codon',paste0(names[n],'_freq'))
    all_freqs$codon<-codon_AA$Codon
    for(y in 1:nrow(all_freqs)){
      codonX<-all_freqs[y,1]
      freqs<-NULL
      codon_table<-eff_data[[n]][,codonX]
      all_freqs[y,2]<-sum(codon_table)/sum(eff_data[[n]])
    }
    eff_output[[n]]<-all_freqs
    print('eff2')
  }
  
  toplot_data<-merge(eff_output[[1]],eff_output[[2]],by='codon')
  for(p in 3:length(eff_output)){
    toplot_data<-merge(toplot_data,eff_output[[p]],by='codon')
  }
  #means gives synonymous codon usage means
  
  toplot_data$third_nt
  for(m in 1:nrow(toplot_data)){
    if(substring(toplot_data[m,'codon'],3,3)=='a'){
      toplot_data[m,'third_nt']<-'magenta'
    }
    if(substring(toplot_data[m,'codon'],3,3)=='c'){
      toplot_data[m,'third_nt']<-'deepskyblue'
    }
    if(substring(toplot_data[m,'codon'],3,3)=='g'){
      toplot_data[m,'third_nt']<-'limegreen'
    }
    if(substring(toplot_data[m,'codon'],3,3)=='t'){
      toplot_data[m,'third_nt']<-'purple'
    }
  }
  
  return(toplot_data)
  
}










###################################################################################
library(colorRamps)
codonfreq_data_for_heatmap<-function(datasets,names){
  codon_AA<-read.table('~/CNOTdata/finaldraftdata_171019/Codons/codon_AA.csv',sep=',',stringsAsFactors = FALSE,header=TRUE)
  codon_AA$Codon<-tolower(codon_AA$Codon)
  codonset<-data.frame(uco(s2c('ATGATCGAACAATGGCCCGGGGCCTGA'),frame=0,as.data.frame=FALSE,index='eff'),stringsAsFactors = FALSE)$Var1
  
  #per mRNA codon usage
  eff_data<-list()
  for(k in 1:length(datasets)){
    
    UCO_alltypes<-lapply(datasets[[k]]$CDS_sequence, function(x) uco(s2c(x), frame = 0, as.data.frame = FALSE, index='eff', NA.rscu = NA))
    UCOalltypesDF_eff<-do.call(rbind,UCO_alltypes) 
    row.names(UCOalltypesDF_eff)<-datasets[[k]]$Gene_ID 
    names(UCOalltypesDF_eff)<-codonset
    print('eff')
    
    
    eff_data[[k]]<-as.data.frame(UCOalltypesDF_eff)
    
  }
  
  
  eff_output<-list()
  for(n in 1:length(eff_data)){
    all_freqs<-data.frame(matrix(nrow=nrow(codon_AA),ncol=2,0),stringsAsFactors = FALSE)
    names(all_freqs)<-c('codon',paste0(names[n],'_freq'))
    all_freqs$codon<-codon_AA$Codon
    for(y in 1:nrow(all_freqs)){
      codonX<-all_freqs[y,1]
      freqs<-NULL
      codon_table<-eff_data[[n]][,codonX]
      all_freqs[y,2]<-sum(codon_table)/sum(eff_data[[n]])
    }
    eff_output[[n]]<-all_freqs
    print('eff2')
  }

  toplot_data<-merge(eff_output[[1]],eff_output[[2]],by='codon')
  for(p in 3:length(eff_output)){
    toplot_data<-merge(toplot_data,eff_output[[p]],by='codon')
  }
  #means gives synonymous codon usage means
  
  toplot_data$third_nt
  for(m in 1:nrow(toplot_data)){
    if(substring(toplot_data[m,'codon'],3,3)=='a'){
      toplot_data[m,'third_nt']<-'magenta'
    }
    if(substring(toplot_data[m,'codon'],3,3)=='c'){
      toplot_data[m,'third_nt']<-'deepskyblue'
    }
    if(substring(toplot_data[m,'codon'],3,3)=='g'){
      toplot_data[m,'third_nt']<-'limegreen'
    }
    if(substring(toplot_data[m,'codon'],3,3)=='t'){
      toplot_data[m,'third_nt']<-'purple'
    }
  }
  
  return(toplot_data)
  
}

datasets_list<-list(RNA_up,RNA_down,RPF_up,RPF_down,Both_up,Both_down,non_sig)
samplenames<-list('RNA up','RNA down','RPF up','RPF down','Both up','Both down','not significant')
freqdataforheatmap<-codonfreq_data_for_heatmap(datasets_list,samplenames)

rownames(freqdataforheatmap)<-freqdataforheatmap$codon
heatmap_matrix<-as.matrix(freqdataforheatmap[,c(2:8)])
tiff(file = "codonfreq_groups_centroid_rowscale.tiff",width=500,height=1200)
heatmap.2(heatmap_matrix,col=colorRamps::blue2green,cexRow=1.2,scale="row", trace="none",density.info='none',RowSideColors = freqdataforheatmap$third_nt, hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(18,12))
dev.off()


tiff(file = "codonfreq_groups_centroid_noscale.tiff",width=500,height=1200)
heatmap.2(heatmap_matrix,col=colorRamps::blue2green,cexRow=1.2,scale="none", trace="none",density.info='none',RowSideColors = freqdataforheatmap$third_nt, hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(18,12))
dev.off()


tiff(file = "codonfreq_groups_centroid_colscale.tiff",width=500,height=1200)
heatmap.2(heatmap_matrix,col=colorRamps::blue2green,cexRow=1.2,scale="col", trace="none",density.info='none',RowSideColors = freqdataforheatmap$third_nt, hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(18,12))
dev.off()



