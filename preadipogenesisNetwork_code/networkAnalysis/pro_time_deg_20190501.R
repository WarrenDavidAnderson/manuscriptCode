
library(dplyr)
library(DESeq2)
library(bigWig)
library(ggplot2)

dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/network_initial/round2")
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
bed0 = read.table("TU_20181107.bed",stringsAsFactors=F,header=F)
names(bed0) = c("chr","start","end","gene","xy","strand")

# TF community genes - tfclass.lists
# see communityTFs.R
load("TFcommunity.gene.lists.RData")

# create TF/comm annotation frame
tf.comm.ann = c()
for(ii in 1:length(tfclass.lists)){
  comm = names(tfclass.lists)[ii]
  tfs = tfclass.lists[[ii]]
  new = cbind(comm, tfs)
  tf.comm.ann = rbind(tf.comm.ann, new)
}
tf.comm.ann = as.data.frame(tf.comm.ann, stringsAsFactors=F)
rownames(tf.comm.ann) = 1:nrow(tf.comm.ann)

############################################################
## key analysis parameters
############################################################

# sig thresh for LRT
sig.thresh = 0.001
fc.thresh = 1

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
# out = list(deg=deg, deseq_obj=deseq_obj)
# save(out, file="LRTres_pro_20190501.RData")
# load("LRTres_pro_20190501.RData")
deg = DESeq(deseq_obj, test="LRT", full=~conditions,  reduced=~1)
res = results(deg)

############################################################
## isolate TUs with early and late dynamics
############################################################

# function to get specific contrast results
get.contrast.res = function(deg=NULL,contrast=NULL,sig.thresh=NULL,fc.thresh=NULL,reads=NULL){
  ctrst = results(deg, contrast=contrast)
  ind.up = which(ctrst$padj < sig.thresh & ctrst$log2FoldChange > fc.thresh)
  ind.down = which(ctrst$padj < sig.thresh & ctrst$log2FoldChange < -fc.thresh)
  genes.up = rownames(ctrst)[ind.up]
  genes.dn = rownames(ctrst)[ind.down]
  nsig = length(which(ctrst$padj < sig.thresh & abs(ctrst$log2FoldChange) > fc.thresh))
  ind.ref=which(conditions==contrast[3])
  ind.cpr=which(conditions==contrast[2])
  ind.pks = sapply(genes.up,function(x)which(rownames(deseq_obj)==x))
  print(head(reads[ind.pks,c(ind.ref,ind.cpr)]))
  mean.ref = apply(reads[ind.pks,c(ind.ref)],1,mean)
  mean.cpr = apply(reads[ind.pks,c(ind.cpr)],1,mean)
  print(head(cbind(mean.ref, mean.cpr),10))
  out = list(contrast=ctrst, genes.up=genes.up, genes.dn=genes.dn, nsig=nsig)
  return(out)
}

# early response contrasts: t0 vs 20 min
ctrst_0_20_pro = get.contrast.res(deg=deg,contrast=c("conditions","0.33","0"),
                              sig.thresh=sig.thresh,fc.thresh=fc.thresh,
                              reads=counts(deseq_obj) %>% as.data.frame)

ctrst_0_40_pro = get.contrast.res(deg=deg,contrast=c("conditions","0.67","0"),
                              sig.thresh=sig.thresh,fc.thresh=fc.thresh,
                              reads=counts(deseq_obj) %>% as.data.frame)

# late response contrasts: t40 vs 3h
ctrst_40_3_pro = get.contrast.res(deg=deg,contrast=c("conditions","3","0.67"),
                              sig.thresh=sig.thresh,fc.thresh=fc.thresh,
                              reads=counts(deseq_obj) %>% as.data.frame)

# late response contrasts: t40 vs 4h
ctrst_40_4_pro = get.contrast.res(deg=deg,contrast=c("conditions","4","0.67"),
                              sig.thresh=sig.thresh,fc.thresh=fc.thresh,
                              reads=counts(deseq_obj) %>% as.data.frame)


# early response increased TUs
tu.early.up.pro = intersect(ctrst_0_20_pro$genes.up, ctrst_0_40_pro$genes.up)
length(tu.early.up.pro)

# early response decreased TUs
tu.early.dn.pro = intersect(ctrst_0_20_pro$genes.dn, ctrst_0_40_pro$genes.dn)
length(tu.early.dn.pro)

# late response increased TUs
tu.late.up.pro = intersect(ctrst_40_3_pro$genes.up, ctrst_40_4_pro$genes.up)
length(tu.late.up.pro)

# late response decreased TUs
tu.late.dn.pro = intersect(ctrst_40_3_pro$genes.dn, ctrst_40_4_pro$genes.dn)
length(tu.late.dn.pro)

############################################################
## isolate early response TFs
############################################################

# function to isolate significant TFs for a particular comparison
# input: sig.genes = char vector of genes with coordinates (e.g., chr1:60434582-60566695_Raph1)
# input: tf.comm.ann = frame with "comm" and "tfs" annotations
get.sig.tfs = function(sig.genes=NULL, tf.comm.ann=NULL){
  sig.genes = sig.genes[-grep("tu_class",sig.genes)]
  sig.genes = sapply(sig.genes,function(x){strsplit(x,"_")[[1]][2]})
  sig.genes = cbind(names(sig.genes),sig.genes) %>% as.data.frame(stringsAsFactors=F)
  names(sig.genes) = c("coord","gene")
  sig.tf = tf.comm.ann$tfs[ tf.comm.ann$tfs %in% sig.genes$gene ]
  ind.tf = sapply(sig.tf,function(x){which(tf.comm.ann$tfs==x)[1]}) %>% unlist
  ind.tu = sapply(sig.tf,function(x){which(sig.genes$gene==x)[1]}) %>% unlist
  sig.tfs = cbind(tf.comm.ann[ind.tf,], sig.genes[ind.tu,])
  return(sig.tfs)
} # get.sig.tfs

### early response TFs -- increased
sig.tfs.early.up = get.sig.tfs(sig.genes=tu.early.up.pro, tf.comm.ann=tf.comm.ann)

### early response TFs -- decreased
sig.tfs.early.dn = get.sig.tfs(sig.genes=tu.early.dn.pro, tf.comm.ann=tf.comm.ann)

### late response TFs -- increased
sig.tfs.late.up = get.sig.tfs(sig.genes=tu.late.up.pro, tf.comm.ann=tf.comm.ann)

### late response TFs -- decreased
sig.tfs.late.dn = get.sig.tfs(sig.genes=tu.late.dn.pro, tf.comm.ann=tf.comm.ann)


############################################################
## MA plots
############################################################

# MA plot function - with subsampling over insignificant peaks
ma.plot.fn = function(ma.dat=NULL,sig.peaks=NULL,col.ref=NULL,col.sig=NULL,nsamp=10000){
  subsamp = sample(c(1:nrow(ma.dat)), size=nsamp, replace =F)
  ind.sig = sapply(sig.peaks,function(x){which(rownames(ma.dat)==x)})
  xlab = "log average"
  ylab = "log fold change"
  xaxis = log2( ma.dat$baseMean )
  yaxis = ma.dat$log2FoldChange
  yy = yaxis[is.na(yaxis)==FALSE]
  xx = xaxis[is.na(xaxis)==FALSE & abs(xaxis)!=Inf]
  ymin = min(yy) - 0.1 * (max(yy) - min(yy))
  ymax = max(yy) + 0.1 * (max(yy) - min(yy))
  xmin = min(xx) - 0.1 * (max(xx) - min(xx))
  xmax = max(xx) + 0.1 * (max(xx) - min(xx))
  limy = max(abs(ymin), ymax)
  plot(xaxis[subsamp],yaxis[subsamp],pch=19,col=col.ref,
       xlab=xlab,ylab=ylab,ylim=c(-limy,limy),xlim=c(xmin,xmax))
  points(xaxis[ind.sig],yaxis[ind.sig],pch=19,col=col.sig)
}

########################
# plot pro dynamic peaks - increased
svg("pro_up_ma.svg"); par(mfrow=c(2,2))
# plot early peaks
ma.dat = ctrst_0_40_pro$contrast
sig.peaks = tu.early.up.pro 
col.ref = "gray"
col.sig = "darkgreen"
ma.plot.fn(ma.dat=ma.dat,sig.peaks=sig.peaks,col.ref=col.ref,col.sig=col.sig)
# plot late peaks
ma.dat = ctrst_40_4_pro$contrast
sig.peaks = tu.late.up.pro
col.ref = "gray"
col.sig = "magenta"
ma.plot.fn(ma.dat=ma.dat,sig.peaks=sig.peaks,col.ref=col.ref,col.sig=col.sig)
dev.off()


############################################################
## get sig TF bed coords and remove redundant clusters
############################################################

# function to apply bed coords and remove unwanted communities
format.fun = function(sig=NULL, bed0=NULL, rem=c(61,267)){
  
  # remove redundant communities
  for(ii in 1:length(rem)){
    sig = sig[-which(sig$comm==paste0("community_",rem[ii])),]
  }
  out = sig
  
  # apply bed coords
  chr = start = end = strand = c()
  for(ii in 1:nrow(sig)){
    gene = sig$gene[ii]
    coords = bed0[bed0$gene==gene,]
    chr = c(chr, coords$chr)
    start = c(start, coords$start)
    end = c(end, coords$end)
    strand = c(strand, coords$strand)
  }
  
  # output
  out = out %>% mutate(chr=chr, start=start, end=end, strand=strand)
  
} # format.fun


# function to apply bed coords to coordinates
coord.bed = function(sig=NULL, bed0=NULL){
  
  # isolate tu ids
  tus = c()
  for(ii in 1:length(sig)){
    id0 = strsplit(sig[ii],"_")[[1]]
    if(length(id0)==4){
      new = paste0(id0[2],"_",id0[3],"_",id0[4])
    }
    if(length(id0)==3){
      new = paste0(id0[2],"_",id0[3])
    }
    if(length(id0)==2){
      new = id0[2]
    }
    tus = c(tus,new)
  }
  
  
  # bed coords
  chr = sapply(sig,function(x){strsplit(x,":")[[1]][1]})
  start = sapply(sig,function(x){strsplit(x,"-")[[1]][1]})
  start = sapply(start,function(x){strsplit(x,":")[[1]][2]})
  end = sapply(sig,function(x){strsplit(x,"_")[[1]][1]})
  end = sapply(end,function(x){strsplit(x,"-")[[1]][2]})
  
  # find strand
  strand = c()
  for(ii in 1:length(sig)){
    gene = tus[ii]
    coords = bed0[bed0$gene==gene,]
    strand = c(strand, coords$strand)
  }
  
  # output
  out = data.frame(chr=chr, start=start, end=end, gene=tus, xy="xy", strand=strand)
  return(out)
} # coord.bed


# add bed data - tu
tu.early.up.pro = coord.bed(sig=tu.early.up.pro, bed0=bed0)
tu.early.dn.pro = coord.bed(sig=tu.early.dn.pro, bed0=bed0)
tu.late.up.pro = coord.bed(sig=tu.late.up.pro, bed0=bed0)
tu.late.dn.pro = coord.bed(sig=tu.late.dn.pro, bed0=bed0)

# add bed data - tf
sig.tfs.early.up = format.fun(sig=sig.tfs.early.up, bed0=bed0)
sig.tfs.late.up = format.fun(sig=sig.tfs.late.up, bed0=bed0)
sig.tfs.early.dn = format.fun(sig=sig.tfs.early.dn, bed0=bed0)
sig.tfs.late.dn = format.fun(sig=sig.tfs.late.dn, bed0=bed0)


############################################################
## add annotation and output
############################################################

# add annotation for class of changes
sig.tfs.early.up = sig.tfs.early.up %>% mutate(class = "early_up")
sig.tfs.late.up = sig.tfs.late.up %>% mutate(class = "late_up")
sig.tfs.early.dn = sig.tfs.early.dn %>% mutate(class = "early_dn")
sig.tfs.late.dn = sig.tfs.late.dn %>% mutate(class = "late_dn")

# output 
fname = "deg_proTF.txt"
out = rbind(sig.tfs.early.up, sig.tfs.late.up, sig.tfs.early.dn, sig.tfs.late.dn)
write.table(out,fname,col.names=T,row.names=F,sep="\t",quote=F)

# save data in a list
out = list(tu.early.up.pro=tu.early.up.pro,
           tu.early.dn.pro=tu.early.dn.pro,
           tu.late.up.pro=tu.late.up.pro,
           tu.late.dn.pro=tu.late.dn.pro)

save(out, file="deg_proTU.RData")

############################################################
## adjust log normalized data for time series plotting
## z-score scale and then average replicates
############################################################

# log transfor data
log2_expr = rlogTransformation(deseq_obj, blind=F)
norm.expr = assay(log2_expr)

# z-score the data for each gene
tf.z.scores = scale(t(norm.expr)) %>% t

# get mean z-scores for each time point
times = c("t0","t20min","t40min","t60min","t2h","t3h","t4h")
t.profiles = sapply(times,function(x){ 
  apply(tf.z.scores[,grep(x,colnames(tf.z.scores))],1,mean) 
})
colnames(t.profiles) = c(0,0.33,0.67,1,2,3,4)

# remove unannotated TUs and re-name genes
t.profiles = t.profiles[-grep("tu_class",rownames(t.profiles)),]
rownames(t.profiles) = sapply(rownames(t.profiles),function(x){strsplit(x,"_")[[1]][2]})


############################################################
## plot temporal profiles for significant TFs
############################################################

# time vector for plotting
tdat = colnames(t.profiles) %>% as.numeric

# plot function
time.plot = function(deg=NULL,res=NULL,condits=NULL,cond.name=NULL,
                     span=0.5,nreps=3,ind=NULL,type=NULL){
  
  tlabs = as.character(condits) %>% as.numeric
  if(is.null(ind)==TRUE){
    ind = which.min(res$padj)
  }
  
  if(type=="pro"){
    gene = strsplit(rownames(res)[ind],"_")[[1]][2]
  }
  
  if(type=="atac"){
    gene = rownames(res)[ind]
  }
  
  plt <- plotCounts(deg, ind, 
                    intgroup=cond.name, returnData=T)
  p1 = ggplot(plt, aes(x=as.numeric(condits), y = count)) + 
    geom_point() + 
    geom_smooth(se=T, method="loess", colour="black", span=span) +
    scale_y_log10() + ylab(paste0(gene, " reads")) +
    scale_x_continuous("adipogenesis time (hrs)",
                       labels=tlabs,
                       breaks=as.numeric(condits))
  return(p1)
} # time.plot

# plot atac data
condits = conditions
ind = grep("Cebpa",rownames(res))
p1 = time.plot(deg=deg,res=res,condits=condits,cond.name="conditions",type="pro",ind=ind)
print(p1)

library(gridExtra)
pdf("Maff.pdf", onefile = TRUE, height=11, width=8.5)
marrangeGrob(grobs=list(p1=p1), nrow=3, ncol=2, top=NULL)
dev.off()

gg = c("Klf4","Cebpb","Egr2","Cebpa")
plts = list()
for(ii in gg){
  condits = conditions
  ind = grep(ii,rownames(res))
  plts[[ii]] = time.plot(deg=deg,res=res,condits=condits,cond.name="conditions",type="pro",ind=ind)
}

library(gridExtra)
pdf("EarlyGenes.pdf", onefile = TRUE, height=11, width=8.5)
marrangeGrob(grobs=plts, nrow=3, ncol=2, top=NULL)
dev.off()

