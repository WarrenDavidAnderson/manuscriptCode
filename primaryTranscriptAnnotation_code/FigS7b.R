
library(dplyr)
library(DESeq2)
library(bigWig)
library(ggplot2)

setwd("/media/wa3j/Seagate2/Documents/PRO/adipogenesis/package/SM/degpca")

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
bed0 = read.table("largest.interval.bed",stringsAsFactors=F,header=F)
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
## PCA 
############################################################

# log transfor data
log2_size = rlogTransformation(deseq_obj, blind=F)
norm.expr.size = assay(log2_size)

# perform PCA
pca.size = prcomp(t(norm.expr.size), center=TRUE, scale.=FALSE)
scores.size = pca.size$x
loadings.size = pca.size$rotation
percent.size = pca.size$sdev^2 / sum( pca.size$sdev^2 )

# plot 2d scores
varcols = c(rep("gray",3), # t0
            rep("blue",3), # t20min
            rep("red",3), # t2h
            rep("cyan4",3), # t3h
            rep("darkgreen",3), # t40min
            rep("magenta",3), # t4h
            rep("purple",3), # t60min
            rep("black",3)) # t6d

# 2d pca plot function
pca.plot2d.single = function(scores,var,cols,fname){
  pdf(fname)
  par(mfrow=c(2,2))
  xlab = paste0("PC1: ",round(var[1] * 100),"% variance")
  ylab = paste0("PC2: ",round(var[2] * 100),"% variance")
  plot(scores[,1], scores[,2], xlab=xlab, ylab=ylab)
  points(scores[,1], scores[,2],col=cols,pch=19,cex=2)
  ylab = paste0("PC3: ",round(var[3] * 100),"% variance")
  plot(scores[,1], scores[,3], xlab=xlab, ylab=ylab)
  points(scores[,1], scores[,3],col=cols,pch=19,cex=2)
  dev.off()
}
pca.plot2d.single(scores.size,percent.size,varcols,"pca_ge.pdf")



