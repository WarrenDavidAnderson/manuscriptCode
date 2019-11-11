

dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
      "adipogenesis/July2018/dREG/results_20181222")
setwd(dir)

options(digits = 10)

library(dplyr)
library(DESeq2)
library(bigWig)

############################################################
## key analysis parameters
############################################################

# set number of base pairs around peaks for motif analysis
pk.dist = 50

# parameters for designating dynamic peaks
sig.thresh = 0.1
fc.thresh = 1

# parameters for designating non-dynamic peaks
sig.un = 0.5
fc.un = 0.25

###################################################################
## read in dREG data and perform some necessary processing
###################################################################

# read in data
fname = "out.dREG.peak.full.bed"
dreg0 = read.table(fname,header=F,stringsAsFactors=F,sep="\t")
names(dreg0) = c("chr","start","end","score","prob","center")

# note that dREG gives some summits outside of 
# the coordinate range (0.2 %)
hist(dreg0$end - dreg0$center)
100*length(which(dreg0$end < dreg0$center)) / nrow(dreg0)

# shift the summits to avert outside summits
dreg = dreg0
ind.hi = which(dreg0$end <= dreg0$center)
ind.lo = which(dreg0$start >= dreg0$center)
dreg$center[ind.hi] = dreg0$center[ind.hi] - (dreg0$center-dreg0$end)[ind.hi] - 1
dreg$center[ind.lo] = dreg0$center[ind.lo] + (dreg0$start-dreg0$center)[ind.lo] + 1
hist(dreg$end - dreg$center)
hist(dreg$center - dreg$start)

# bed.map = dreg
# save(bed.map, file="bed.mapDreg.RData")

# separate positive and negative regions
# negative: start to summit
# positive: summit to end
dreg.minus = dreg %>% select(chr,start,center)
dreg.plus = dreg %>% select(chr,center,end)
min(dreg.minus[,3] - dreg.minus[,2])
min(dreg.plus[,3] - dreg.plus[,2])

###################################################################
## map PRO reads to dREG elements
###################################################################

# function to import proseq read data maped to specific coordinates
# input: bed = bed file for mapping reads
# input: bigWig.path = path to bigWig files
# input: file.prefix = file prefix
# input: file.suffix = file suffix
# input: min.length = min gene length
# output: data frame with reads mapped for each condition
get.counts.interval <- function(bed.plus=NULL, bed.minus=NULL,
                                file.prefix=NULL, file.suffix=NULL,
                                bigWig.path=NULL) {
  
  # output files
  vec.names = c()
  out.plus = data.frame(matrix(ncol = 0, nrow = nrow(bed.plus)))
  out.minus = data.frame(matrix(ncol = 0, nrow = nrow(bed.plus)))
  output = data.frame(matrix(ncol = 0, nrow = nrow(bed.plus)))
  
  # loop through each condition and load reads
  for(bigWig in Sys.glob(file.path(bigWig.path, 
                                   paste0(file.prefix, paste0("*",file.suffix))))) {
    factor.name = strsplit(bigWig, "/")[[1]]
    factor.name = strsplit(factor.name[length(factor.name)],file.suffix)[[1]][1]
    factor.name = strsplit(factor.name,file.prefix)[[1]][2]
    strand = strsplit(factor.name,"sample_")[[1]][2]
    factor.name = strsplit(factor.name,"_sample")[[1]][1]
    loaded.bw = load.bigWig(bigWig)
    if(strand=="plus"){
      vec.names = c(vec.names, factor.name)
      reads = bed.region.bpQuery.bigWig(loaded.bw, bed.plus)
      out.plus = cbind(out.plus, reads)
      print(factor.name)
    }
    if(strand=="minus"){
      reads = bed.region.bpQuery.bigWig(loaded.bw, bed.minus)
      out.minus = cbind(out.minus, reads)
    }
  } # bigWig
  
  # annotation and output
  output = as.data.frame(out.plus + out.minus)
  names(output) = vec.names
  r.names = paste0(bed.plus[,1],':',bed.minus[,3])
  row.names(output) = r.names
  return(output)
} # get.counts.interval

# map PRO data to dreg regions
# note that feature names are chr:summit
bigWig.path = paste0("/media/wa3j/Seagate2/Documents/PRO/adipogenesis",
                     "/July2018/bigWig_20181001/pro_preadip_bigWig")
file.prefix = "3T3_"
file.suffix = ".bigWig"
reads0 = get.counts.interval(bed.plus=dreg.plus,
                             bed.minus=dreg.minus,
                             bigWig.path=bigWig.path, 
                             file.prefix=file.prefix,
                             file.suffix=file.suffix)

# apply this row name ajustment to generate output for motif enrichment analysis
# do not apply this for generating coordinates for motif identification 
# see get.map.from.res() below
# rownamen = paste0(dreg$chr,":",dreg$start,"-",dreg$end)
# rownames(reads0) = rownamen

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
############################################################

# basic counts, size factors, and annotation
counts = counts(deseq_obj_sizefac)
sizefac = sizeFactors(deseq_obj_sizefac)
times = colData(deseq_obj_sizefac)$conditions
t.order = cbind(c("t0","t20min","t40min","t60min","t2h","t3h","t4h","t6d"),
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

# deseq_obj_preadip = deseq_obj
# save(deseq_obj_preadip, file="dregDeseq.RData")

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
    ind_ii = which(conditions == condits[ii])
    ind_jj = which(conditions == condits[jj])
    
    # set new deseq object
    des = conditions[c(ind_ii,ind_jj)]
    colData_ij = as.data.frame(des)
    cnt_ij = counts[,c(ind_ii,ind_jj)]
    deseq_obj_preadip = DESeqDataSetFromMatrix(countData=cnt_ij, 
                                               colData=colData_ij, 
                                               design=~des)
    sizeFactors(deseq_obj_preadip) = sizefac[c(ind_ii,ind_jj)]
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
get.map.from.res = function(res=NULL, pk.dist=NULL){
  namen = rownames(res)
  chrs = sapply(namen,function(x){strsplit(x,":")[[1]][1]})
  mids = sapply(namen,function(x){strsplit(x,":")[[1]][2]}) %>% as.numeric
  strs = mids - pk.dist
  ends = mids + pk.dist
  out = cbind(chrs,strs,ends,rownames(res),res$log2FoldChange,"+",res$padj) %>% 
    as.data.frame(stringsAsFactors=F)
  names(out) = c("chr","start","end","peak_name","log2fc","str","fdr")
  return(out)
} # get.map.from.res

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
  
  # generate output for meme
  out.sig.up = get.map.from.res(res=res_ii[ind.up,], pk.dist=pk.dist)
  out.sig.dn = get.map.from.res(res=res_ii[ind.dn,], pk.dist=pk.dist)
  out.uns = get.map.from.res(res=res_ii[ind.un,], pk.dist=pk.dist)

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


