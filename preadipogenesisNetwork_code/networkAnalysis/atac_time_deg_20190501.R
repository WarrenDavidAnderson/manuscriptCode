
library(dplyr)
library(DESeq2)
library(ggplot2)
library(bigWig)
library(prodlim)


setwd("/run/user/1001/gvfs/smb-share:server=home1.virginia.edu,share=cphgdesk/users/wa3j/MGlab/Analysis/Adipogenesis/ATAC_time_DEG")
source("atac_norm_functions.R")

dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/network_initial/round2")
setwd(dir)

############################################################
## import data
############################################################

# load atac peak coordinates - bed.map
load("bed.map.RData")


############################################################
## key analysis parameters
############################################################

# parameters for designating dynamic peaks
sig.thresh = 0.001
fc.thresh = 1


############################################################
## import atac-seq data and set directory for data output
############################################################

# atac data mapped to peak regions
bigWig.path = paste0("/media/wa3j/Seagate2/Documents/PRO/adipogenesis",
                     "/July2018/atac_bw_20181127/atac_all")
file.prefix = "3T3_"
file.suffix = ".bigWig"
min.length = 1
reads0 = get.counts.interval(bed=bed.map[,1:3], bigWig.path=bigWig.path,
                             file.prefix=file.prefix,
                             file.suffix=file.suffix,
                              min.length=min.length)

# remove 6d data
reads0 = reads0[,-grep("6d",colnames(reads0))]

############################################################
## normalize raw counts to size factors, generate DESeq object
############################################################

# estimate size factors
size_factors = estimateSizeFactorsForMatrix(reads0)

## generate DESeqDataSet object 
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
## exclude 6day data
############################################################

# basic counts, size factors, and annotation
counts = counts(deseq_obj_sizefac)
sizefac = sizeFactors(deseq_obj_sizefac)
times = colData(deseq_obj_sizefac)$conditions
t.order = cbind(c("t0","20min","40min","60min","2hr","3hr","4hr"),
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

### OUTPUT THE SIZEFAC
atac.sizefac = sizeFactors(deseq_obj)
save(atac.sizefac, file="atac.sizefac.RData")

###############
# check
cnt = counts(deseq_obj)
sf = sizeFactors(deseq_obj)
ind = which(rownames(cnt) == "chr4:55555384-55555643")
cnt[46791:46792,]

############################################################
## implement differential expression analysis with LRT
############################################################

# basic analysis
# out = list(deg=deg, deseq_obj=deseq_obj)
# save(out, file="LRTres_atac_20190501.RData")
# load("LRTres_atac_20190501.RData")
deg = DESeq(deseq_obj, test="LRT", full=~conditions,  reduced=~1)
res = results(deg)
ind.sig = which(res$padj<sig.thresh & abs(res$log2FoldChange)>fc.thresh)
n.sig.pks = length(ind.sig)
frac.sig.pks = length(ind.sig) / nrow(res)

############################################################
## isolate peaks with early and late increases
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
  ind.pks = sapply(genes.up,function(x)which(rownames(deg)==x))
  print(head(reads[ind.pks,c(ind.ref,ind.cpr)]))
  mean.ref = apply(reads[ind.pks,c(ind.ref)],1,mean)
  mean.cpr = apply(reads[ind.pks,c(ind.cpr)],1,mean)
  print(head(cbind(mean.ref, mean.cpr),10))
  out = list(contrast=ctrst, genes.up=genes.up, genes.dn=genes.dn, nsig=nsig)
  return(out)
}

# early response contrasts: t0 vs 20 min
ctrst_0_20_atac = get.contrast.res(deg=deg,contrast=c("conditions","0.33","0"),
                 sig.thresh=sig.thresh,fc.thresh=fc.thresh,
                 reads=counts(deseq_obj) %>% as.data.frame)

ctrst_0_40_atac = get.contrast.res(deg=deg,contrast=c("conditions","0.67","0"),
                              sig.thresh=sig.thresh,fc.thresh=fc.thresh,
                              reads=counts(deseq_obj) %>% as.data.frame)

# late response contrasts: t40 vs 3h
ctrst_40_3_atac = get.contrast.res(deg=deg,contrast=c("conditions","3","0.67"),
                              sig.thresh=sig.thresh,fc.thresh=fc.thresh,
                              reads=counts(deseq_obj) %>% as.data.frame)

# late response contrasts: t40 vs 4h
ctrst_40_4_atac = get.contrast.res(deg=deg,contrast=c("conditions","4","0.67"),
                              sig.thresh=sig.thresh,fc.thresh=fc.thresh,
                              reads=counts(deseq_obj) %>% as.data.frame)


# early response increased REs
re.early.up.atac = intersect(ctrst_0_20_atac$genes.up, ctrst_0_40_atac$genes.up)
length(re.early.up.atac)

# late response increased REs
re.late.up.atac = intersect(ctrst_40_3_atac$genes.up, ctrst_40_4_atac$genes.up)
length(re.late.up.atac)

# early response decreased REs
re.early.dn.atac = intersect(ctrst_0_20_atac$genes.dn, ctrst_0_40_atac$genes.dn)
length(re.early.dn.atac)

# late response decreased REs
re.late.dn.atac = intersect(ctrst_40_3_atac$genes.dn, ctrst_40_4_atac$genes.dn)
length(re.late.dn.atac)


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

###########################
# plot atac dynamic peaks
svg("atac_ma.svg"); par(mfrow=c(2,2))

# plot early peaks
ma.dat = ctrst_0_40_atac$contrast
sig.peaks = re.early.up.atac
col.ref = "gray"
col.sig = "darkgreen"
ma.plot.fn(ma.dat=ma.dat,sig.peaks=sig.peaks,col.ref=col.ref,col.sig=col.sig)

# plot late peaks
ma.dat = ctrst_40_4_atac$contrast
sig.peaks = re.late.up.atac
col.ref = "gray"
col.sig = "magenta"
ma.plot.fn(ma.dat=ma.dat,sig.peaks=sig.peaks,col.ref=col.ref,col.sig=col.sig)
dev.off()

############################################################
## output data
############################################################

# save data in a list
out = list(re.early.up.atac=re.early.up.atac,
           re.late.up.atac=re.late.up.atac,
           re.early.dn.atac=re.early.dn.atac,
           re.late.dn.atac=re.late.dn.atac)

save(out, file="deg_atac.RData")


############################################################
## adjust log normalized data for time series plotting
## z-score scale and then average replicates
############################################################

# log transfor data
log2_expr = rlogTransformation(deseq_obj, blind=F)
norm.expr = assay(log2_expr)

# z-score the data for each gene
z.scores = scale(t(norm.expr)) %>% t

# get mean z-scores for each time point
times = c("t0","20min","40min","60min","2h","3h","4h")
t.profiles = sapply(times,function(x){ 
  apply(z.scores[,grep(x,colnames(z.scores))],1,mean) 
})
colnames(t.profiles) = c(0,0.33,0.67,1,2,3,4)


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

load("net.filt.RData")
txMaf = net.filt %>% filter(geneTU=="Maff",tt=="0_0.33,0.67")

# plot atac data
condits = conditions
ind = grep(txMaf$idRE[2],rownames(res))
p1 = time.plot(deg=deg,res=res,condits=condits,cond.name="conditions",type="atac",ind=ind)
print(p1)

library(gridExtra)
pdf("Maffre.pdf", onefile = TRUE, height=11, width=8.5)
marrangeGrob(grobs=list(p1=p1), nrow=3, ncol=2, top=NULL)
dev.off()
