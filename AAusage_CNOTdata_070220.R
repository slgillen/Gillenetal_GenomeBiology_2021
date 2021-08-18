library(gplots)
library("RColorBrewer")
library(Biostrings)
library(seqinr)
library(ggplot2)


setwd("~/CNOTdata/finaldraftdata_171019/AAfreq_060220")

AA_data_for_heatmap<-function(datasets,names){
  codon_AA<-read.table('~/CNOTdata/finaldraftdata_171019/Codons/codon_AA.csv',sep=',',stringsAsFactors = FALSE,header=TRUE)
  codon_AA$Codon<-tolower(codon_AA$Codon)
  codonset<-as.vector(data.frame(uco(s2c('ATGATCGAACAATGGCCCGGGGCCTGA'),frame=0,as.data.frame=FALSE,index='eff'),stringsAsFactors = FALSE)$Var1)
  excludelist<-c('tag','taa','tga')
  #codonset<-codonset[(codonset %in% excludelist)==FALSE]
  #per mRNA codon usage
  eff_data<-list()
  All_AA_meanfreqs<-NULL
  for(k in 1:length(datasets)){
    datasets[[k]]$CDS_sequence_noATG<-unlist(lapply(datasets[[k]]$CDS_sequence,function(x) substring(x,4,nchar(x)))) #check this
    UCO_alltypes<-lapply(datasets[[k]]$CDS_sequence_noATG, function(x) uco(s2c(x), frame = 0, as.data.frame = FALSE, index='eff', NA.rscu = NA))
    
    UCOalltypesDF_eff<-do.call(rbind,UCO_alltypes) 
    UCOalltypesDF_eff<-data.frame(UCOalltypesDF_eff,stringsAsFactors = FALSE)
    row.names(UCOalltypesDF_eff)<-datasets[[k]]$Gene_ID 
    UCOalltypesDF_eff<-UCOalltypesDF_eff[,(names(UCOalltypesDF_eff) %in% excludelist)==FALSE]
    #names(UCOalltypesDF_eff)<-codonset
    #print(dim(UCOalltypesDF_eff))
    
    eff_data[[k]]<-as.data.frame(UCOalltypesDF_eff)
    #print(names(UCOalltypesDF_eff))
  }

  #return(eff_data)
  uniqueAA<-unique(codon_AA$AminoAcid)
  uniqueAA<-c("Ile" ,"Leu","Val","Phe","Cys","Ala","Gly","Pro","Thr","Ser","Tyr","Gln","Asn","His","Glu","Asp","Lys","Arg","Trp","Met")
  eff_output<-list()
  for(n in 1:length(eff_data)){
    all_means<-NULL #add in group label // column as groups // rownames as codons 
    for(AAcid in uniqueAA){
      #print(AAcid)
      #for each amino acid sum the data for each codon
      #then get a frequency along the mRNA - may want to consider the pool of codons #may want to consider when taking into account mRNA abundance
      codons<-subset(codon_AA,AminoAcid==AAcid)$Codon
      #print(codons)
      codon_table<-eff_data[[n]][,codons]
      #print(dim(codon_table))
      #print(length(codon_table))
      #print(AAcid)
      if(length(codons)>1){
        codon_table[,AAcid]<-rowSums(codon_table)
        all_means<-cbind(all_means,codon_table[,AAcid])
        #names(all_means)[ncol(all_means)]<-AAcid
      }else{
        
        all_means<-cbind(all_means,codon_table)
        #names(all_means)[ncol(all_means)]<-AAcid
      }
      
      #return(codon_table)
    }
    #return(all_means)
    all_freqs<-data.frame(t(as.matrix(apply(all_means,1,function(x) x/sum(x)),stringsAsFactors = FALSE)))

    names(all_freqs)<-uniqueAA
    #return(all_freqs)
    meanfreqs<-data.frame(apply(all_freqs,2,function(x) mean(x,na.rm=TRUE)),stringsAsFactors = FALSE) #remove NaN from means # need to think about zeros 
    #for(c in 1:length(codons)){
    #  means<-rbind(means,c(codons[c],mean(codon_output2[,c])))
    #}
    #return(meanfreqs)
    names(meanfreqs)<-c(paste0('AA_freq_',names[n]))
    meanfreqs$AA<-rownames(meanfreqs)
    
    if(n==1){
      All_AA_meanfreqs<-meanfreqs
    }else{
      All_AA_meanfreqs<-merge(All_AA_meanfreqs,meanfreqs,by='AA')
    }
  }
  return(All_AA_meanfreqs)
}

datasets_list<-list(subset(GC_table,Group=='not significant'),new2_MitoticCellCycle,new2_PatternSpecificationProcess)
samplenames<-list('not significant','Mitotic','Pattern')



datasets_list<-list(subset(GC_table,Group=='not significant'),subset(GC_table,Group=='RNA up'),subset(GC_table,Group=='RNA down'),subset(GC_table,Group=='Both up'),subset(GC_table,Group=='Both down'),subset(GC_table,Group=='RPF up'),subset(GC_table,Group=='RPF down'))
samplenames<-list('not significant','RNA up','RNA down','Both up','Both down','RPF up','RPF down')
dataforheatmap<-AA_data_for_heatmap(datasets_list,samplenames)
rownames(dataforheatmap)<-dataforheatmap$AA
heatmap_matrix<-as.matrix(dataforheatmap[,c(2:8)])


tiff(file = "AAfreq_groups_centroid_rowscale.tiff",width=500,height=650)
heatmap.2(heatmap_matrix,col=colorspace::diverge_hsv(15),cexRow=1,scale="row", trace="none",density.info='none',dendrogram='row', hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(18,12))
dev.off()

tiff(file = "AAfreq_groups_centroid_noscale.tiff",width=600,height=650)
heatmap.2(heatmap_matrix,col=colorspace::diverge_hsv(15),cexRow=1,scale="none", trace="none",density.info='none', hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(18,12))
dev.off()


tiff(file = "AAfreq_groups_centroid_colscale.tiff",width=600,height=650)
heatmap.2(heatmap_matrix,col=colorspace::diverge_hsv(15),cexRow=1,scale="col", trace="none",density.info='none', hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(18,12))
dev.off()





#apply z-score to rows before clustering? Keep column order as in previous plots?
datasets_list<-list(subset(GC_table,Group=='not significant'),subset(GC_table,Group=='RNA up'),subset(GC_table,Group=='RNA down'),subset(GC_table,Group=='Both up'),subset(GC_table,Group=='Both down'),subset(GC_table,Group=='RPF up'),subset(GC_table,Group=='RPF down'))
samplenames<-list('not significant','RNA up','RNA down','Both up','Both down','RPF up','RPF down')
dataforheatmap<-AA_data_for_heatmap(datasets_list,samplenames)
rownames(dataforheatmap)<-dataforheatmap$AA
#heatmap_df<-as.data.frame(heatmap_matrix)
for(k in 1:nrow(dataforheatmap)){
  dataforheatmap[k,c(2:8)]<-scale(as.numeric(dataforheatmap[k,c(2:8)]))
}
AA_neg<-c('Glu','Asp')
AA_pos<-c('Lys','Arg','His')
AA_polar<-c('Cys','Gln','Ser','Asn','Tyr','Thr')

dataforheatmap$AA_type
for(m in 1:nrow(dataforheatmap)){
  if(dataforheatmap[m,'AA'] %in% AA_polar){
    dataforheatmap[m,'AA_type']<-'darkorange2'
  }else{
    if(dataforheatmap[m,'AA'] %in% AA_pos){
      dataforheatmap[m,'AA_type']<-'forestgreen'
    }else{
      if(dataforheatmap[m,'AA'] %in% AA_neg){
        dataforheatmap[m,'AA_type']<-'red3'
      }else{
        dataforheatmap[m,'AA_type']<-'grey66'
      }
    }
  }
}


#zscore already applied
tiff(file = "zscore_preapplied_AAfreq_groups_centroid_noscale.tiff",width=500,height=650)
heatmap.2(as.matrix(dataforheatmap[,c(2:8)]),col=colorspace::diverge_hsv(15),Rowv=TRUE,Colv=FALSE,cexRow=1,dendrogram='row',scale="none",RowSideColors = dataforheatmap$AA_type, trace="none",density.info='none', hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(20,12))
dev.off()

tiff(file = "zscore_preapplied_AAfreq_groups_centroid_noscale_v2.tiff",width=500,height=650)
heatmap.2(as.matrix(dataforheatmap[,c(2:8)]),col=colorspace::diverge_hsv(15),Rowv=TRUE,Colv=TRUE,cexRow=1,dendrogram='both',scale="none", RowSideColors = dataforheatmap$AA_type,trace="none",density.info='none', hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(20,12))
dev.off()

tiff(file = "zscore_preapplied_AAfreq_groups_centroid_noscale_v3.tiff",width=400,height=650)
heatmap.2(as.matrix(dataforheatmap[,c(2:8)]),col=colorspace::diverge_hsv(15),Rowv=TRUE,Colv=TRUE,cexRow=1,dendrogram='none',scale="none", RowSideColors = dataforheatmap$AA_type,trace="none",density.info='none', hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(16,10))
dev.off()





tiff(file = "mvp.tiff",width=500,height=650)
heatmap.2(as.matrix(dataforheatmap[,2:4]),col=colorspace::diverge_hsv(15),cexRow=1,scale="row", trace="none",density.info='none',dendrogram='row', hclustfun = function(x) hclust(x,method = 'centroid'),margins = c(18,12))
dev.off()



