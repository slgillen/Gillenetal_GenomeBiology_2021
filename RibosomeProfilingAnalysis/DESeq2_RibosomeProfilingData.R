#load libraries
library(DESeq2)
library(limma)
library(apeglm)


#MDS plot of data - sanity check of the data
plotMDS(as.matrix(ReadCounts[,c(2:13)])) #all samples
plotMDS(as.matrix(ReadCounts[,c(2:7)])) #Total RNA
plotMDS(as.matrix(ReadCounts[,c(8:13)])) #RPF

ReadCounts<-subset(ReadCounts,Rep1_siControl_TotalRNA>5 & Rep2_siControl_TotalRNA>5 & Rep3_siControl_TotalRNA>5) #pre-filtering to remove genes with no/very low reads


################################## Total RNA samples ########################################

#DESeq2 requires information about the data #assay condition batch
sampleInfo<-read.table('TotalRNA_DESeq2_sampleinfo.txt',stringsAsFactors = FALSE,header=TRUE)
sampleInfo$condition<-factor(sampleInfo$condition,levels=c('siControl','siCNOT1'))
sampleInfo$batch<-factor(sampleInfo$batch,levels=c('R1','R2','R3'))
sampleInfo

#Select data
countsMatrix<-ReadCounts[,c(1,2:7)] 
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
TotalRNA_results<-as.data.frame(resOrdered)
nrow(TotalRNA_results)
TotalRNA_results$ENSG<-rownames(TotalRNA_results)
TotalRNA_results<-merge(GeneIDs,TotalRNA_results,by='ENSG')
nrow(TotalRNA_results)
write.csv(TotalRNA_results, file="TotalRNA_siCNOT1vsiControl_DESeq2output.csv",row.names=FALSE,quote=FALSE)



################################## RPF samples ########################################

#DESeq2 requires information about the data #assay condition batch
sampleInfo<-read.table('RPF_DESeq2_sampleinfo.txt',stringsAsFactors = FALSE,header=TRUE)
sampleInfo$condition<-factor(sampleInfo$condition,levels=c('siControl','siCNOT1'))
sampleInfo$batch<-factor(sampleInfo$batch,levels=c('R1','R2','R3'))
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
RPF_results<-as.data.frame(resOrdered)
nrow(RPF_results)
RPF_results$ENSG<-rownames(RPF_results)
RPF_results<-merge(GeneIDs,RPF_results,by='ENSG')
nrow(RPF_results)
write.csv(RPF_results, file="RPF_siCNOT1vsiControl_DESeq2output.csv",row.names=FALSE,quote=FALSE)
