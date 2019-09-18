
library(dplyr)
library(DESeq2)
library(ggplot2)
library(bigWig)
library(prodlim)

setwd("/run/user/1001/gvfs/smb-share:server=home1.virginia.edu,share=cphgdesk/users/wa3j/MGlab/Analysis/Adipogenesis/ATAC_time_DEG")
source("atac_norm_functions.R")

dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/atac_time_deg")
setwd(dir)


############################################################
## key analysis parameters
############################################################

# set number of base pairs around peaks for motif analysis
pk.dist = 50

# parameters for designating dynamic peaks
sig.thresh = 0.001
fc.thresh = 1

# parameters for designating non-dynamic peaks
sig.un = 0.5
fc.un = 0.25


############################################################
## import macs2 peak data
############################################################

# peak and summit coordinates (summits.bed, from empiricalSummits.R)
load("ATACsummits_20190914.RData")
bed0 = summits.bed

# chrom sizes
chrm.size = read.table("mm10.chrom.sizes",header=F,stringsAsFactors=F,sep="\t")
chrm.size = cbind(chrm.size[,1], 1, chrm.size[,2])
chrm.size = as.data.frame(chrm.size, stringsAsFactors=F)
names(chrm.size) = c("chr","start","end")
chrm.size[,2:3] = apply(chrm.size[,2:3],2,function(x){data.matrix(x) %>% as.numeric})

# write peaks
# write.table(bed0[,1:3],"atacPeaks.bed0",col.names=F,row.names=F,sep="\t",quote=F)

############################################################
## process peak/summit data 
############################################################

# filter chromosomes
unique(bed0$chr)
dim(bed0) # 86345
chr.keep = paste0("chr",c(1:19))
bed0 = bed0[bed0$chr %in% chr.keep,]
dim(bed0) # 83716

# look at peak distances
d = bed0$end - bed0$start
median(d)
quantile(d)

# save coords
save(bed0, file="bed.map20190827.RData")


############################################################
## import atac-seq data and set directory for data output
############################################################

# atac data mapped to peak regions
bigWig.path = paste0("/media/wa3j/Seagate2/Documents/PRO/adipogenesis",
                     "/July2018/atac_bw_20181127/atac_all")
file.prefix = "3T3_"
file.suffix = ".bigWig"
min.length = 1
reads0 = get.counts.interval(bed=bed0[,1:3], bigWig.path=bigWig.path, 
                             file.prefix=file.prefix,
                             file.suffix=file.suffix,
                             min.length=min.length)


############################################################
## normalize raw counts to size factors, generate DESeq object
############################################################

# estimate size factors
size_factors = estimateSizeFactorsForMatrix(reads0)

## generate DESeqDataSet object and transform count data to log2 scale
conditions = sapply(names(reads0),function(x){strsplit(x,"_")[[1]][1]})
colData = as.data.frame(conditions)
deseq_obj_sizefac = DESeqDataSetFromMatrix(countData=reads0, 
                                           colData=colData, design=~conditions)
sizeFactors(deseq_obj_sizefac) = size_factors

# note that each column of the count matrix was divided 
# by the respective size factor
head(sapply(c(1:24),function(x){reads0[,x]/size_factors[x]})[,1:5])
head(counts(deseq_obj_sizefac, normalized=TRUE)[,1:5])

# save.image("atac.RData")
# load("atac.RData")

############################################################
## adjust data organization
############################################################

# basic counts, size factors, and annotation
counts = counts(deseq_obj_sizefac)
sizefac = sizeFactors(deseq_obj_sizefac)
times = colData(deseq_obj_sizefac)$conditions
t.order = cbind(c("t0","20min","40min","60min","2hr","3hr","4hr","6d"),
                c(0,0.33,0.67,1,2,3,4,144)) %>% 
  as.data.frame(stringsAsFactors=FALSE)
names(t.order) = c("condition","time")

# specify new conditions as times (factor)
conditions = sapply(as.character(times),function(x){
  t.order$time[which(t.order$condition==x)]})
conditions = factor(conditions, levels=c(0,0.33,0.67,1,2,3,4,144))

# set new deseq object
colData = as.data.frame(conditions)
deseq_obj = DESeqDataSetFromMatrix(countData=counts, 
                                   colData=colData, design=~conditions)
sizeFactors(deseq_obj) = estimateSizeFactorsForMatrix(counts)


############################################################
## exclude 6day data
############################################################

# remove the 6day data
counts = counts(deseq_obj) %>% as.data.frame
counts_preadip = counts[,-grep("6d",names(counts))]
sizefac_preadip = sizeFactors(deseq_obj)[-grep("6d",names(counts))]
conditions_preadip = colData(deseq_obj)$conditions[-grep("6d",names(counts))]
condit.map = cbind(as.character(conditions_preadip), names(counts_preadip)) %>% 
  as.data.frame(stringsAsFactors=F)
names(condit.map) = c("time","rep")

# set new deseq object
colData_preadip = as.data.frame(conditions_preadip)
deseq_obj_preadip = DESeqDataSetFromMatrix(countData=counts_preadip, 
                                           colData=colData_preadip, 
                                           design=~conditions_preadip)
sizeFactors(deseq_obj_preadip) = sizefac_preadip

############################################################
## generate all pairwise comparisons
############################################################

# loop through all pairwise combinations and run DESeq 
condits = c(0,0.33,0.67,1,2,3,4) %>% as.character
res.pairs = list()
for(ii in 1:(length(condits)-1)){
  for(jj in (ii+1):length(condits)){
    
    print(paste0("compare ",condits[ii]," to ",condits[jj]))
    
    # basic data annotation
    ind_ii = which(conditions_preadip == condits[ii])
    ind_jj = which(conditions_preadip == condits[jj])
    
    # set new deseq object
    des = conditions_preadip[c(ind_ii,ind_jj)]
    colData_ij = as.data.frame(des)
    cnt_ij = counts_preadip[,c(ind_ii,ind_jj)]
    deseq_obj_preadip = DESeqDataSetFromMatrix(countData=cnt_ij, 
                                               colData=colData_ij, 
                                               design=~des)
    sizeFactors(deseq_obj_preadip) = sizefac_preadip[c(ind_ii,ind_jj)]
    colData(deseq_obj_preadip)$condition = des
    
    # pairwise DEG analysis
    dds = DESeq(deseq_obj_preadip)
    res = results(dds)[order(results(dds)$padj),]
    res.pairs[[paste0(condits[ii],"_",condits[jj])]] = res
    
  } # jj
} # ii

# save(res.pairs, file="pairwise.deg.RData")
# load("pairwise.deg.RData") 

############################################################
## generate output for meme
############################################################

# function to find peak info from rownames of the results frame
get.map.from.res = function(res,map){
  namen = rownames(res)
  chrs = sapply(namen,function(x){strsplit(x,":")[[1]][1]})
  ends = sapply(namen,function(x){strsplit(x,"-")[[1]][2]})
  strs = sapply(namen,function(x){y=strsplit(x,"-")[[1]][1]; 
        return(strsplit(y,":")[[1]][2]) })
  coords = cbind(chrs,strs,ends) %>% as.data.frame(stringsAsFactors=F)
  out = row.match(coords, map[,1:3])
  return( cbind(res, map[out,]) )
} # get.map.from.res

# set coordinates around each summit based on pk.dist
# generate output format for meme
pk.coords.for.meme = function(res=NULL, pk.dist=NULL){
  out = c()
  for(ii in 1:nrow(res)){
    pks = res$summits[ii]
    new = c(res$chr[ii], pks-pk.dist, pks+pk.dist, rownames(res)[ii],
              signif(res$log2FoldChange[ii],3), "+", signif(res$padj[ii],3) )
    out = rbind(out,new) 
  } # ii
  out = out %>% as.data.frame(stringsAsFactors=FALSE)
  names(out) = c("chr","start","end","peak_name","log2fc","str","fdr")
  rownames(out) = c(1:nrow(out))
  return(out)
} # pk.coords.for.meme

# loop through all pairwise data and write output for meme
for(ii in 1:length(res.pairs)){
  
  # basic annotation
  comp_ii = names(res.pairs)[[ii]]
  res_ii = res.pairs[[ii]]
  indsig = which(res_ii$padj<sig.thresh & abs(res_ii$log2FoldChange)>fc.thresh)
  ind.up = which(res_ii$padj<sig.thresh & res_ii$log2FoldChange>fc.thresh)
  ind.dn = which(res_ii$padj<sig.thresh & res_ii$log2FoldChange<(-1)*fc.thresh)
  ind.un = which(res_ii$padj>sig.un & abs(res_ii$log2FoldChange)<fc.un)
  
  # if there are more insignificant peaks, use the number of sig peaks
  # select those with the highest FDRs
  if(length(indsig)==0){next}
  if(length(ind.un) > length(indsig)){ind.un = rev(ind.un)[1:length(indsig)]}
  
  # get macs2 data for sig and unsig peaks
  sig.pks.up = get.map.from.res(res=res_ii[ind.up,], map=bed0)
  sig.pks.dn = get.map.from.res(res=res_ii[ind.dn,], map=bed0)
  uns.pks = get.map.from.res(res=res_ii[ind.un,], map=bed0)
  out.sig.up = pk.coords.for.meme(res=sig.pks.up, pk.dist=pk.dist)
  out.sig.dn = pk.coords.for.meme(res=sig.pks.dn, pk.dist=pk.dist)
  out.uns = pk.coords.for.meme(res=uns.pks, pk.dist=pk.dist)
  
  # set directory and output data for meme analysis
  dir_ii = paste0("meme_",comp_ii)
  system(paste0("mkdir ",dir_ii))
  setwd(dir_ii)
  fname = paste0("upsig_",comp_ii,".bed")
  write.table(out.sig.up,fname,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  fname = paste0("downsig_",comp_ii,".bed")
  write.table(out.sig.dn,fname,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  fname = paste0("unsig_",comp_ii,".bed")
  write.table(out.uns,fname,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  setwd(dir)
  
} # ii



