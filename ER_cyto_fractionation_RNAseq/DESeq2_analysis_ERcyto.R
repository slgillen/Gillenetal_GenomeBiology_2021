#load libraries
library(DESeq2)
library(limma)
library(apeglm)


#read in samples
sample_names<-c('Rep1_siControl_cyto','Rep2_siControl_cyto','Rep3_siControl_cyto','Rep1_siCNOT1_cyto','Rep2_siCNOT1_cyto','Rep3_siCNOT1_cyto','Rep1_siControl_ER','Rep2_siControl_ER','Rep3_siControl_ER','Rep1_siCNOT1_ER','Rep2_siCNOT1_ER','Rep3_siCNOT1_ER')

samples<-list()
for(i in 1:length(sample_names)){
 xx<-read.table(paste0('fcounts_',sample_names[i],'_matrix.txt'),stringsAsFactors = FALSE,header=TRUE) #inputs are formatted outputs from featureCounts
 names(xx)<-c('ENSG',sample_names[i])
 samples[[i]]<-xx
 rm(xx)
}

ERcyto_readcounts<-samples[[1]]
for(i in 2:length(samples)){
  ERcyto_readcounts<-merge(ERcyto_readcounts,samples[[i]],by='ENSG')
}
dim(ERcyto_readcounts)
ERcyto_readcounts<-subset(ERcyto_readcounts,Rep1_siControl_cyto>5 & Rep2_siControl_cyto>5 & Rep3_siControl_cyto>5) #pre-filtering to remove genes with no/very low reads
nrow(ERcyto_readcounts)
names(ERcyto_readcounts)

#normalise library sizes for MDS plots #but for the DESeq2 use the unnormalised raw counts
ERcyto_readcounts2<-ERcyto_readcounts
apply(ERcyto_readcounts2[,2:13],2,function(x) sum(x))
ERcyto_readcounts2[,2:13]<-apply(ERcyto_readcounts2[,2:13],2,function(x) (x/sum(x))*1000000)
apply(ERcyto_readcounts2[,2:13],2,function(x) sum(x))

#MDS plot of data - sanity check of the data
plotMDS(as.matrix(ERcyto_readcounts2[,c(2:13)])) #all samples
plotMDS(as.matrix(ERcyto_readcounts2[,c(2:7)])) #cyto
plotMDS(as.matrix(ERcyto_readcounts2[,c(8:13)])) #ER
plotMDS(as.matrix(ERcyto_readcounts2[,c(2,3,4,8,9,10)])) #siControl
plotMDS(as.matrix(ERcyto_readcounts2[,c(5,6,7,11,12,13)])) #siCNOT1

#look at CNOT1 in the data table to check for the knockdown - i.e. that all samples look like they are labelled correctly and no mix up has occurred
subset(ERcyto_readcounts2,startsWith(ENSG,'ENSG00000125107')==TRUE)



################################## CYTO samples ########################################

#DESeq2 requires information about the data #assay condition batch
sampleInfo<-read.table('ERcyto_DESeq2_sampleinfo.txt',stringsAsFactors = FALSE,header=TRUE)
sampleInfo$assay<-factor(sampleInfo$assay,levels=c('cyto','ER'))
sampleInfo$condition<-factor(sampleInfo$condition,levels=c('siControl','siCNOT1'))
sampleInfo$batch<-factor(sampleInfo$batch,levels=c('R1','R2','R3'))
sampleInfo

#get sample info needed for this comparison
sampleInfo<-sampleInfo[1:6,]
sampleInfo

#Select data
countsMatrix<-ReadCountsTable[,c(1,2:7)] 
names(countsMatrix)
rownames(countsMatrix)<-countsMatrix[,1]
countsMatrix<-as.matrix(countsMatrix[,c(-1)])
head(countsMatrix)

#do the analysis
DESeq2data<-DESeqDataSetFromMatrix(countData = countsMatrix,colData = sampleInfo,design= ~ batch + condition)
DESeq2output <- DESeq(DESeq2data)
resultsNames(DESeq2output) 
resLFC <- lfcShrink(DESeq2output, coef='condition_siCNOT1_vs_siControl', type="apeglm") 

#Sort results
resOrdered <- resLFC[order(resLFC$pvalue),]
summary(resLFC)
sum(resLFC$padj < 0.05, na.rm=TRUE) 

#look at results
plotMA(resLFC, ylim=c(-2,2))

#results table
cyto_results<-as.data.frame(resOrdered)
nrow(cyto_results)
cyto_results$ENSG<-rownames(cyto_results)
cyto_results<-merge(GeneIDs,cyto_results,by='ENSG')
nrow(cyto_results)
write.csv(cyto_results, file="Cyto_siCNOT1vsiControl_DESeq2output.csv",row.names=FALSE,quote=FALSE)



################################## ER samples ########################################

#DESeq2 requires information about the data #assay condition batch
sampleInfo<-read.table('ERcyto_DESeq2_sampleinfo.txt',stringsAsFactors = FALSE,header=TRUE)
sampleInfo$assay<-factor(sampleInfo$assay,levels=c('Cyto','ER'))
sampleInfo$condition<-factor(sampleInfo$condition,levels=c('siControl','siCNOT1'))
sampleInfo$batch<-factor(sampleInfo$batch,levels=c('R1','R2','R3'))
sampleInfo

#get sample info needed for this comparison
sampleInfo<-sampleInfo[7:12,]
sampleInfo

#select data
countsMatrix<-ReadCountsTable[,c(1,8:13)] 
names(countsMatrix)
rownames(countsMatrix)<-countsMatrix[,1]
countsMatrix<-as.matrix(countsMatrix[,c(-1)])
head(countsMatrix)

#do the analysis
DESeq2data<-DESeqDataSetFromMatrix(countData = countsMatrix,colData = sampleInfo,design= ~ batch + condition)
DESeq2output <- DESeq(DESeq2data)
resultsNames(DESeq2output) 
resLFC <- lfcShrink(DESeq2output, coef='condition_siCNOT1_vs_siControl', type="apeglm") 

#Sort results
resOrdered <- resLFC[order(resLFC$pvalue),]
summary(resLFC)
sum(resLFC$padj < 0.05, na.rm=TRUE) 

#look at results
plotMA(resLFC, ylim=c(-2,2))

#results table
ER_results<-as.data.frame(resOrdered)
nrow(ER_results)
ER_results$ENSG<-rownames(ER_results)
ER_results<-merge(GeneIDs,ER_results,by='ENSG')
nrow(ER_results)
write.csv(ER_results, file="ER_siCNOT1vsiControl_DESeq2output.csv",row.names=FALSE,quote=FALSE)




