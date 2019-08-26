
library(dplyr)
library(DESeq2)
library(bigWig)
library(ggplot2)

setwd("/media/wa3j/Seagate2/Documents/PRO/adipogenesis/package/SM/degpca")
setwd(dir)

############################################################
## load data objects and functions
## see pro_norm.R for basic normalization 
## and TUmap_analysis.R for mapping TUs
############################################################

# import pro data
normfun = paste0("/run/user/1001/gvfs/smb-share:server=home1.virginia.edu,",
                 "share=cphgdesk/users/wa3j/MGlab/Analysis/Adipogenesis",
                 "/spike_norm/pro_normalization/pro_norm_functions.R")
source(normfun)

# TU annotation 
bed0 = read.table("TU_20190731.bed",stringsAsFactors=F,header=F)
names(bed0) = c("chr","start","end","gene","xy","strand")

############################################################
## import pro-seq data and set directory for data output
# focus only on preadipogenisis time points (omit 6d)
############################################################

# import pro data
bigWig.path = paste0("/media/wa3j/Seagate2/Documents/PRO/adipogenesis",
                     "/July2018/bigWig_20181001/all_pro_sample_bigWig")
file.prefix = "3T3_"
file.suffix = "_sample_plus.bigWig"
min.length = 1
reads0 = get.counts.interval(bed=bed0, bigWig.path=bigWig.path,
                             file.prefix=file.prefix,
                             file.suffix=file.suffix,
                             min.length=min.length)

# remove 6d data
reads0 = reads0[,-grep("t6d",colnames(reads0))]

############################################################
## normalize raw counts to size factors, generate DESeq object
############################################################

# estimate size factors
size_factors = estimateSizeFactorsForMatrix(reads0)

## generate DESeqDataSet 
conditions = sapply(names(reads0),function(x){strsplit(x,"_")[[1]][1]})
colData = as.data.frame(conditions)
deseq_obj_sizefac = DESeqDataSetFromMatrix(countData=reads0, 
                                           colData=colData, design=~conditions)
sizeFactors(deseq_obj_sizefac) = size_factors

# note that each column of the count matrix was divided 
# by the respective size factor
head(sapply(c(1:21),function(x){reads0[,x]/size_factors[x]})[,1:5])
head(counts(deseq_obj_sizefac, normalized=TRUE)[,1:5])

############################################################
## adjust data organization
############################################################

# basic counts, size factors, and annotation
counts = counts(deseq_obj_sizefac)
sizefac = sizeFactors(deseq_obj_sizefac)
times = colData(deseq_obj_sizefac)$conditions
t.order = cbind(c("t0","t20min","t40min","t60min","t2h","t3h","t4h"),
                c(0,0.33,0.67,1,2,3,4)) %>% 
  as.data.frame(stringsAsFactors=FALSE)
names(t.order) = c("condition","time")

# specify new conditions as times (factor)
conditions = sapply(as.character(times),function(x){
  t.order$time[which(t.order$condition==x)]})
conditions = factor(conditions, levels=c(0,0.33,0.67,1,2,3,4))

# set new deseq object
colData = as.data.frame(conditions)
deseq_obj = DESeqDataSetFromMatrix(countData=counts, 
                                   colData=colData, design=~conditions)
sizeFactors(deseq_obj) = estimateSizeFactorsForMatrix(counts)


############################################################
## implement differential expression analysis - LRT
############################################################

# basic analysis
deg = DESeq(deseq_obj, test="LRT", full=~conditions,  reduced=~1)
res = results(deg)

# out = list(deg=deg, deseq_obj=deseq_obj)
# save(out, file="LRTres_pro_TU.RData")




