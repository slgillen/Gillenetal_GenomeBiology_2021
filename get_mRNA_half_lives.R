#load libraries
library(stats)
library(ggplot2)
library(tidyr)
library(dplyr)
library(magrittr)
library(broom)
library('robustbase')
library('aomisc')


#import counts data
#repeat for each replicate
siControl<-read.table('siControl_counts.txt',stringsAsFactors=FALSE)
names(siControl)[3:9]<-c('con_0','con_05','con_1','con_2','con_4','con_8','con_16')

apply(siControl[,3:9], 2,function(x) sum(x))

#normalise library size
for(i in 1:7){
  siControl[,i+2]<-(siControl[,i+2]/sum(siControl[,i+2]))*1000000
}

siCNOT1<-read.table('siCNOT1_counts.txt',stringsAsFactors=FALSE)
names(siCNOT1)[3:9]<-c('not_0','not_05','not_1','not_2','not_4','not_8','not_16')

apply(siCNOT1[,3:9], 2,function(x) sum(x))

#normalise library size
for(i in 1:7){
  siCNOT1[,i+2]<-(siCNOT1[,i+2]/sum(siCNOT1[,i+2]))*1000000
}

#normalise to the nanodrop concentrations for the RNA samples
data_tables<-list(siControl_rep1,siCNOT1_rep1,siControl_rep2,siCNOT1_rep2,siControl_rep3,siCNOT1_rep3)

nanodrop_concs<-read.table('~/CNOTdata/mRNAhalflives/nanodrop_RNA_concentrations.csv',sep=',',stringsAsFactors = FALSE,header=TRUE)
nanodrop_concs<-nanodrop_concs[,c(2:7)]
for(j in 1:length(data_tables)){
  dfx<-data_tables[[j]]
  for(m in 1:7){
    dfx[,m+2]<-dfx[,m+2]*nanodrop_concs[m,j]
  }
  dfx<-subset(dfx,dfx[,3]>20) #Q.C remove those with less than 20 read counts at the 0hr time point
  data_tables[[j]]<-dfx 
  print(nrow(data_tables[[j]]))
  print(apply(data_tables[[j]][,3:9],2, function(x) sum(x)))
}


######## modelling to get decay rate and mRNA half life ############

#sort data format
long_format_data<-list()
for(n in 1:6){
  prefix<-c('con_','not_','con_','not_','con_','not_')
  tomodel<-data_tables[[n]]
  
  #normalise each gene to 0hr values (0hr=100%) 
  tomodelnorm[,3:9]<-t(data.frame(apply(tomodel[,3:9],1,function(x) as.vector(as.vector(x)/as.vector(x)[1])*100)))
  
  #reformat the data
  tomodel_long <- gather(tomodel[,c(2:9)], condition, measurement, colnames(tomodel)[3]:colnames(tomodel)[9], factor_key=TRUE)
  tomodel_long$condition<-sub(prefix[n],'',tomodel_long$condition) 
  tomodel_long$condition<-sub('05','0.5',tomodel_long$condition)
  tomodel_long$condition<-as.numeric(tomodel_long$condition)
  tomodel_long<-tomodel_long[,c(2,3,1)]
  names(tomodel_long)<-c('time','y','gene_name') 
  long_format_data[[n]]<-tomodel_long
  rm(tomodel)
  rm(tomodel_long)  
}



################ model decay rates to obtain half-life estimates ###############

#only use genes detected in all three replicates
control_genes<-unique(subset(long_format_data[[1]],(unique(long_format_data[[1]]$gene_name) %in% unique(long_format_data[[3]]$gene_name)) == TRUE)$gene_name)
control_genes<-control_genes[(control_genes %in% unique(long_format_data[[5]]$gene_name)) == TRUE]
siControl_tomodel<-rbind(subset(long_format_data[[1]],(gene_name %in% control_genes)==TRUE),subset(long_format_data[[3]],(gene_name %in% control_genes)==TRUE),subset(long_format_data[[5]],(gene_name %in% control_genes)==TRUE))

not_genes<-unique(subset(long_format_data[[2]],(unique(long_format_data[[2]]$gene_name) %in% unique(long_format_data[[4]]$gene_name)) == TRUE)$gene_name)
not_genes<-not_genes[(not_genes %in% unique(long_format_data[[6]]$gene_name)) == TRUE]
siCNOT1_tomodel<-rbind(subset(long_format_data[[2]],(gene_name %in% not_genes)==TRUE),subset(long_format_data[[4]],(gene_name %in% not_genes)==TRUE),subset(long_format_data[[6]],(gene_name %in% not_genes)==TRUE))


#data from all three biological replicates used together to get most accurate model of mRNA half-life 

######## model mRNA half-life in siControl condition ##########
siControl_tomodel$y<-siControl_tomodel$y+0.01 #pseudo adjustment because 0 values (present at 8/16hr for ~10 genes) causes error in modelling

#self start function for start parameter estimation
expoDecay.Init <- function(RHS, LHS, data) {
  xy <- sortedXyData(RHS, LHS, data)
  lmFit <- lm(log(xy[, "y"]) ~ xy[, "x"])
  coefs <- coef(lmFit)
  C0 <- exp(coefs[1])
  k <- - coefs[2]
  value <- c(C0, k)
  names(value) <- RHS[c("C0", "k")]
  value
}

#detect problem genes that are unable to be modelled
genes<-unique(siControl_tomodel$gene_name)

for(i in 1:length(genes)){
  tryCatch({
    gene_data <- subset(siControl_tomodel, gene_name==genes[i])
    fits[[i]] <- nlrob(formula=y~y0*exp(-k*time), data=gene_data,start=list(y0=expoDecay.Init(gene_data$time,gene_data$y,gene_data)[1],k= expoDecay.Init(gene_data$time,gene_data$y,gene_data)[2]),maxit=50)
  }, error=function(e){cat("ERROR :",conditionMessage(e),genes[i], "\n")} ) #retain gene names of error genes
}

siControl_tomodel<-subset(siControl_tomodel, (gene_name %in% problem_genes)==FALSE)


#robust model for outlier removal
#parameters
K<-2 #number of parameters modelled -> y0 and k
N<-21
dof<-19 #(N-K)


outputdf<-siControl_tomodel
genes<-unique(siControl_tomodel$gene_name)
length(genes)
outliers_rownames<-c()

for(g in 1:length(genes)){
  gene_data<-subset(outputdf,gene_name==genes[g])
  robmodel<-nlrob(formula=y~y0*exp(-k*ttime), data=gene_data,start=list(y0=expoDecay.Init(gene_data$time,gene_data$y,gene_data)[1],k= expoDecay.Init(gene_data$time,gene_data$y,gene_data)[2]),maxit=50)
  gene_data<-as.data.frame(gene_data)
  gene_data$weighted.resid<-abs(robmodel$residuals)/robmodel$fitted.values
  gene_data<-gene_data[order(gene_data$weighted.resid),]
  RSDR<-gene_data[,'weighted.resid'][round(0.6827*N)] 
  RSDR<-RSDR*(N/(N-K))  
  outmin<-ceiling(0.8*N) #mset max outliers to 20%
  pos_outliers<-as.vector(gene_data[,'weighted.resid'][outmin:N])
  
  for(i in 1:length(pos_outliers)){
    ai<-(0.05*(N-((outmin+(i-1))-1)))/N 
    tr<-pos_outliers[i]/RSDR
    tdistpval<-pt(-tr,df=dof)
    if((2*tdistpval)<ai){
      outliers_rownames<-c(outliers_rownames,rownames(tt)[(outmin+(i-1)):N])
      break #stop if condition met because then all other values are outliers
    }
  }
}

length(outliers_rownames)/nrow(siControl_tomodel) #fraction of all data points that are outliers


#apply non-linear model after outlier removal to get decay rate estimates
fit_output<-siControl_tomodel %>% 
  group_by(gene_name) %>% 
  do(fit=nls(formula=y~y0*exp(-k*ttime), data=.,start=list(y0=expoDecay.Init(.$ttime,.$y,.)[1],k= expoDecay.Init(.$ttime,.$y,.)[2]))) %>% 
  tidy(fit) %>% 
  select(gene_name, term, estimate) %>% 
  spread(term, estimate) 

View(fit_output)


# convert decay rates to mRNA half-lives
siControl<-fit_output
siControl$half_life<-unlist(sapply(siControl$k, function(x) log(2)/x))
names(siControl)[2:4]<-paste0('siControl_',names(siControl)[2:4])



######## model mRNA half-life in siCNOT1 condition ##########
siCNOT1_tomodel$y<-siCNOT1_tomodel$y+0.01 #pseudo adjustment because 0 values (present at 8/16hr for ~10 genes) causes error in modelling

#self start function for start parameter estimation
expoDecay.Init <- function(RHS, LHS, data) {
  xy <- sortedXyData(RHS, LHS, data)
  lmFit <- lm(log(xy[, "y"]) ~ xy[, "x"])
  coefs <- coef(lmFit)
  C0 <- exp(coefs[1])
  k <- - coefs[2]
  value <- c(C0, k)
  names(value) <- RHS[c("C0", "k")]
  value
}

#detect problem genes that are unable to be modelled
genes<-unique(siCNOT1_tomodel$gene_name)

for(i in 1:length(genes)){
  tryCatch({
    gene_data <- subset(siCNOT1_tomodel, gene_name==genes[i])
    fits[[i]] <- nlrob(formula=y~y0*exp(-k*time), data=gene_data,start=list(y0=expoDecay.Init(gene_data$time,gene_data$y,gene_data)[1],k= expoDecay.Init(gene_data$time,gene_data$y,gene_data)[2]),maxit=50)
  }, error=function(e){cat("ERROR :",conditionMessage(e),genes[i], "\n")} ) #retain gene names of error genes
}

siCNOT1_tomodel<-subset(siCNOT1_tomodel, (gene_name %in% problem_genes)==FALSE)


#robust model for outlier removal
#parameters
K<-2 #number of parameters modelled -> y0 and k
N<-21
dof<-19 #(N-K)


outputdf<-siControl_tomodel
genes<-unique(siCNOT1_tomodel$gene_name)
length(genes)
outliers_rownames<-c()

for(g in 1:length(genes)){
  gene_data<-subset(outputdf,gene_name==genes[g])
  robmodel<-nlrob(formula=y~y0*exp(-k*ttime), data=gene_data,start=list(y0=expoDecay.Init(gene_data$time,gene_data$y,gene_data)[1],k= expoDecay.Init(gene_data$time,gene_data$y,gene_data)[2]),maxit=50)
  gene_data<-as.data.frame(gene_data)
  gene_data$weighted.resid<-abs(robmodel$residuals)/robmodel$fitted.values
  gene_data<-gene_data[order(gene_data$weighted.resid),]
  RSDR<-gene_data[,'weighted.resid'][round(0.6827*N)] 
  RSDR<-RSDR*(N/(N-K))  
  outmin<-ceiling(0.8*N) #mset max outliers to 20%
  pos_outliers<-as.vector(gene_data[,'weighted.resid'][outmin:N])
  
  for(i in 1:length(pos_outliers)){
    ai<-(0.05*(N-((outmin+(i-1))-1)))/N 
    tr<-pos_outliers[i]/RSDR
    tdistpval<-pt(-tr,df=dof)
    if((2*tdistpval)<ai){
      outliers_rownames<-c(outliers_rownames,rownames(tt)[(outmin+(i-1)):N])
      break #stop if condition met because then all other values are outliers
    }
  }
}

length(outliers_rownames)/nrow(siCNOT1_tomodel) #fraction of all data points that are outliers


#apply non-linear model after outlier removal to get decay rate estimates
fit_output<-siCNOT1_tomodel %>% 
  group_by(gene_name) %>% 
  do(fit=nls(formula=y~y0*exp(-k*ttime), data=.,start=list(y0=expoDecay.Init(.$ttime,.$y,.)[1],k= expoDecay.Init(.$ttime,.$y,.)[2]))) %>% 
  tidy(fit) %>% 
  select(gene_name, term, estimate) %>% 
  spread(term, estimate) 

View(fit_output)

# convert decay rates to mRNA half-lives
siCNOT1<-fit_output
siCNOT1$half_life<-unlist(sapply(siCNOT1$k, function(x) log(2)/x))
names(siCNOT1)[2:4]<-paste0('siCNOT1_',names(siCNOT1)[2:4])


combined_data<-merge(siControl,siCNOT1,by='gene_name')

