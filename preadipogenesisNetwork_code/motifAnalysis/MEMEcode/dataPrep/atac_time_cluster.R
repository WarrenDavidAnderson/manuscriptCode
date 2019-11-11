
library(dplyr)
library(DESeq2)
library(ggplot2)
library(bigWig)
library(prodlim)

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

# parameters for designating non-dynamic peaks
sig.un = 0.5


############################################################
## import macs2 peak data
############################################################

# load bed0 - filtered ATAC peak coords, see atac_time_deg_meme.R
load("bed.map20191014.RData")

############################################################
## import atac-seq data and set directory for data output
############################################################

# atac data mapped to peak regions
bigWig.path = paste0("/media/wa3j/Seagate2/Documents/PRO/adipogenesis",
                     "/July2018/atac_bw_20191014/atac_all")
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
## note. enter size factors before implementing log transform
conditions = sapply(names(reads0),function(x){strsplit(x,"_")[[1]][1]})
colData = as.data.frame(conditions)
deseq_obj_sizefac = DESeqDataSetFromMatrix(countData=reads0, 
                                           colData=colData, design=~conditions)
sizeFactors(deseq_obj_sizefac) = size_factors


############################################################
## adjust data organization
## exclude 6day data
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

# remove the 6day data
counts = counts(deseq_obj) %>% as.data.frame
counts_preadip = counts[,-grep("6d",names(counts))]
sizefac_preadip = sizeFactors(deseq_obj)[-grep("6d",names(counts))]
conditions_preadip = colData(deseq_obj)$conditions[-grep("6d",names(counts))]

# set new deseq object
colData_preadip = as.data.frame(conditions_preadip)
deseq_obj_preadip = DESeqDataSetFromMatrix(countData=counts_preadip, 
                                           colData=colData_preadip, 
                                           design=~conditions_preadip)
sizeFactors(deseq_obj_preadip) = sizefac_preadip

############################################################
## implement differential peak analysis
############################################################

# differential expreassion analysis - temporal dynamics based on LRT
deg_preadip = DESeq(deseq_obj_preadip, test="LRT", 
                    full=~conditions_preadip,  reduced=~1)
res_preadip = results(deg_preadip)

############################################################
## isolate dynamic peaks for meme
############################################################

# set significant and unsignificant peaks
ind.sig = which(res_preadip$padj<sig.thresh)
ind.un = which(res_preadip$padj>sig.un)
length(ind.sig)
length(ind.un)

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

# get data for sig and unsig peaks
sig.pks = get.map.from.res(res=res_preadip[ind.sig,], map=bed0)
uns.pks = get.map.from.res(res=res_preadip[ind.un,], map=bed0)

# set coordinates around each summit based on pk.dist
pk.coords.for.meme = function(res=NULL, pk.dist=NULL){
  out = c()
  for(ii in 1:nrow(res)){
    pks = res$summits[ii]
    for(jj in 1:length(pks)){
      new = c(res$chr[ii], pks[jj]-pk.dist, pks[jj]+pk.dist, rownames(res)[ii],
              signif(res$log2FoldChange[ii],3), signif(res$padj[ii],3) )
      out = rbind(out,new) 
    } # jj
  } # ii
  out = out %>% as.data.frame(stringsAsFactors=FALSE)
  names(out) = c("chr","start","end","peakID","log2fc","fdr")
  rownames(out) = c(1:nrow(out))
  return(out)
} # pk.coords.for.meme

out.sig = pk.coords.for.meme(res=sig.pks, pk.dist=pk.dist)
out.uns = pk.coords.for.meme(res=uns.pks, pk.dist=pk.dist)

# save insignificant peaks for enrichment analysis
# see Motif_in_peak_diffDyn.R and enrichSummary_diffDyn.sh
save(out.uns, file="out.uns.RData")

############################################################
## data transformations for STEM
## /media/wa3j/Seagate2/Documents/software/stem
## java -mx1024M -jar /media/wa3j/Seagate2/Documents/software/stem/stem.jar
############################################################

# significance filter
ind.sig = sapply(out.sig$peakID,function(x){which(rownames(deseq_obj_preadip)==x)})

# log transformation
log2_exp = rlogTransformation(deseq_obj_preadip[ind.sig,], blind=TRUE)
log.pk.cnt = assay(log2_exp)

# z-score the data for each peak
pk.z.scores = scale(t(log.pk.cnt)) %>% t

# get mean z-scores for each time point
times = c("t0","20min","40min","60min","2h","3h","4h")
t.profiles = sapply(times,function(x){ 
  apply(pk.z.scores[,grep(x,colnames(pk.z.scores))],1,mean) 
})
colnames(t.profiles) = c(0,0.33,0.67,1,2,3,4)

# revise names for timeseries clustering
gene_symbol = rownames(t.profiles)
out = cbind(gene_symbol, t.profiles)
write.table(out,"pkdata_20191014.txt",sep="\t",quote=F,col.names=T,row.names=F)

############################################################
## process/graph STEM clustering results
############################################################

# function to capitalize first letter
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

# import STEM results, adjust gene names
stem.res0 = read.table("time_cluster_results_20191014.txt",sep="\t",
                       header=T,fill=T,stringsAsFactors=F,skip=1)
stem.res = stem.res0 %>% select(gene_symbol, Profile)
stem.res$gene_symbol = tolower(stem.res$gene_symbol)

# write out coordinate information for enrichment analysis
# see Motif_in_peak_diffDyn.R and enrichSummary_diffDyn.sh
save(stem.res, file="stem.res.RData")

# threshold number of profiles for plotting
n.profiles = 0

# generate plots
library(graphics)
tt = c(0,0.33,0.67,1,2,3,4)
pdf("ATACpeak_timeseries_clusters.pdf"); par(mfrow=c(3,3))
for(ii in 1:length(unique(stem.res$Profile))){
  pr = unique(stem.res$Profile)[ii]
  pk.pr = stem.res$gene_symbol[which(stem.res$Profile==pr)]
  if(length(pk.pr) < n.profiles){next}
  ex.pr = t.profiles[rownames(t.profiles) %in% pk.pr,]
  av.pr = apply(ex.pr,2,mean)
  tp = rep("l",nrow(ex.pr))
  ly = rep(1,nrow(ex.pr))
  matplot(tt,t(ex.pr),type=tp,lty=ly,col="black",
          ylab=paste0("profile ",pr),main="",xlab="time (hrs)")
  lines(tt,av.pr,col="red",cex=4)
}
dev.off()


############################################################
## output profiles of interest for meme analysis
############################################################

# clusters
profiles = unique(stem.res$Profile)

# loop through profiles and write data for meme
for(ii in profiles){
  pk.pr = stem.res$gene_symbol[which(stem.res$Profile==ii)] %>%
    unlist
  pk.un = out.uns$peakID[sample.int( nrow(out.uns), min(nrow(out.uns), length(pk.pr)) )]
  fname = paste0("preadip_sigDynamics_profile",ii,".bed")
  out = out.sig[out.sig$peakID %in% pk.pr,]
  write.table(out,fname,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  fname = paste0("preadip_unDynamics_profile",ii,".bed")
  out = out.uns[out.uns$peakID %in% pk.un,]
  write.table(out,fname,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
} ## ii, loop profiles



