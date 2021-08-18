#useful parallelising info
library("BiocParallel")
register(MulticoreParam(4))
#then just set parallel=TRUE in the functions when called


# input format ############################
#files	assay	condition	batch
#count_con1_RPF	RPF Control	1

#don't need assay if separating out (although will compare DESeq2 results to RiboDiff)

# DESeq2 #need to update R to >=3.5 ##############
library(DESeq2)

DESeq2output <- DESeqDataSetFromMatrix(countData = countsMatrix,
                              colData = sampleInfo,
                              design= ~ batch + condition) #batch is the replicate #design = ~ var1 + var2 .... #put variable of interest last
DESeq2output$condition <- factor(DESeq2output$condition, levels = c("siControl","siNOT1"))

#DESeqDataSetFromHTSeq could be used
#ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
#                                       directory = directory,
#                                       design= ~ condition)


DESeq2output <- DESeq(DESeq2output)
resultsNames(DESeq2output) # lists the coefficients
res <- results(DESeq2output, name="Totals_NOTvCON") #or can use numbers for the coefficients
res <- results(DESeq2output, contrast=c("condition","Control","NOT")) #same as above just uses contrasts and doesn't require previously set factors

# or to shrink log fold changes association with condition:
resLFC <- lfcShrink(DESeq2output, coef="Totals_NOTvCON", type="apeglm") #can try different shrinkage options in the case that sometimes shrinkage can be too strong


#order results by smallest p-value
resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE) #see number of genes that are significant

res05 <- results(dds, alpha=0.05) #alpha = q-value???
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

plotMA(res, ylim=c(-2,2)) #points colours red if p.adj<0.1 #points outside if the specified window are as triangles
plotMA(resLFC, ylim=c(-2,2))

#After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. 
#One can then recover the gene identifiers by saving the resulting indices
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]


#can plot the counts for a single gene
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE) #if would prefer to have the df to plot in ggplot



#exporting results
write.csv(as.data.frame(resOrdered), file="condition_treated_results.csv")




########################   DESeq2 Totals   ###############################
countsMatrix<-as.matrix(countsMatrix)
rownames(sampleInfo)<-sampleInfo[,1]
sampleInfo<-sampleInfo[,c(-1)]
sampleInfo[,3]<-c('R1','R2','R3','R1','R2','R3')
DESeq2output <- DESeqDataSetFromMatrix(countData = countsMatrix,colData = sampleInfo,design= ~ batch + condition) #batch is the replicate #design = ~ var1 + var2 .... #put variable of interest last
DESeq2output$condition <- factor(DESeq2output$condition, levels = c("Control","NOT"))
DESeq2output <- DESeq(DESeq2output)
resultsNames(DESeq2output) # lists the coefficients
#head(DESeq2output)
#res <- results(DESeq2output, name="conditionNOT") #or can use numbers for the coefficients
res <- results(DESeq2output, contrast=c("condition","NOT","Control")) #same as above just uses contrasts and doesn't require previously set factors
head(res)
# or to shrink log fold changes association with condition:
#resLFC <- lfcShrink(DESeq2output, coef="Totals_NOTvCON", type="apeglm") #can try different shrinkage options in the case that sometimes shrinkage can be too strong
#resLFC <- lfcShrink(DESeq2output, contrast=c("condition","Control","NOT"), type="apeglm") #can try different shrinkage options in the case that sometimes shrinkage can be too strong

#head(resLFC)  #need updated DESeq2 for lfcShrink
#order results by smallest p-value
resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE) #see number of genes that are significant

res05 <- results(DESeq2output, alpha=0.05) #alpha = q-value???
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

plotMA(res, ylim=c(-2,2)) #points colours red if p.adj<0.1 #points outside if the specified window are as triangles
#plotMA(resLFC, ylim=c(-2,2))

#exporting results
write.csv(as.data.frame(resOrdered), file="~/Desktop/CNOTdataAnalysis/DESeq2/RPFsCDSfilt_DESeq2output.csv")



