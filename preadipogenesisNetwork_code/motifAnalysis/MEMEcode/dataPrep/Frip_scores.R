
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

# unmodified peak data
fname = paste0(dir,"/atac4h_0.1_5_50/3T3_atac_peaks.narrowPeak")
macs0 = read.table(fname,stringsAsFactors=F,header=F)
names(macs0) = c("chr","start","end","name","score",
                 "na","foldchange","logp","logq","rel")

# peak coordinates
# https://github.com/taoliu/MACS
fname = "3T3_atac_4hpeaks_sorted.bed"
bed0 = read.table(fname,stringsAsFactors=F,header=F)
names(bed0) = c("chr","start","end","rel","summit")
bed.filtered.raw = bed0

# chrom sizes
chrm.size = read.table("mm10.chrom.sizes",header=F,stringsAsFactors=F,sep="\t")
chrm.size = cbind(chrm.size[,1], 1, chrm.size[,2])
chrm.size = as.data.frame(chrm.size, stringsAsFactors=F)
names(chrm.size) = c("chr","start","end")
chrm.size[,2:3] = apply(chrm.size[,2:3],2,function(x){data.matrix(x) %>% as.numeric})

# write peaks
# write.table(bed0[,1:3],"atacPeaks.bed0",col.names=F,row.names=F,sep="\t",quote=F)

############################################################
## basic characteristics of MACS2 output
############################################################

d = macs0$end - macs0$start
min(d)
median(d)
quantile(d)

dim(macs0)
min(macs0$start)

unique(macs0$chr)
100* length(which(macs0$chr == "chrM")) / nrow(macs0)

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

# check for the max number of summits in a combined peak
num.summits = sapply(bed0$summit,function(x){
  strsplit(x,",")[[1]] %>% length}) 
max.summits = max(num.summits)
cnt.summits = cbind(c(1:max.summits),sapply(c(1:max.summits),
              function(x){which(num.summits==x) %>% length})) 
cnt.summits = as.data.frame(cnt.summits)
names(cnt.summits) = c("num_summit","count")
hist(num.summits,breaks=c(1:max.summits))

# loop through each peak and keep peaks with summits inside
bed.sing = c() # peaks with one summit
bed.mult = c() # peaks with multiple summits
aa = 0
for(ii in 1:nrow(bed0)){
  
  # get summit(s) and peak coordinates
  pk = strsplit(bed0$summit[ii],",")[[1]] %>% as.numeric
  st = bed0$start[ii]
  ed = bed0$end[ii]
    
  # if there is one summit inside the peak, keep the peak
  if(length(pk)==1 & st<pk[1] & ed>pk[1]){
    bed.sing = rbind(bed.sing,bed0[ii,])
    next
  }
   
  # if there is more than one summit, 
  # keep peaks with summits inside the peak coords
  # remove summits < pk.dist from previous summit
  # use peak center if sumits are outside of the peak coords
  sumit = sapply(pk,function(x){
    if(x>=st & x<=ed){return(x)}
    else{return(round(st+(ed-st)/2))}
  }) %>% unlist
  sumit = sumit[order(sumit)]
  for(jj in 2:length(sumit)){
    if(jj > length(sumit)){break}
    if(sumit[jj]-sumit[jj-1] < pk.dist){sumit = sumit[-c(jj)]; aa=aa+1}
  }
  sumit = paste(sumit, collapse=",")
  out = bed0[ii,]
  out$summit = sumit
  bed.mult = rbind(bed.mult,out) 
    
} # ii, peak loop

# combine results
nrow(bed.sing) + nrow(bed.mult)
dim(bed0)
bed.map = rbind(bed.sing, bed.mult)
save(bed.map, file="bed.map20190827.RData")


############################################################
## calculate FRiP
############################################################

# atac data mapped to peak regions and chrom.sizes
bigWig.path = paste0("/media/wa3j/Seagate2/Documents/PRO/adipogenesis",
                     "/July2018/atac_bw_20181127/atac_all")
file.prefix = "3T3_"
file.suffix = ".bigWig"
min.length = 1
reads.pks = get.counts.interval(bed=bed.filtered.raw[,1:3], 
                                bigWig.path=bigWig.path, 
                                file.prefix=file.prefix,
                                file.suffix=file.suffix,
                                min.length=min.length)
reads.tot = get.counts.interval(bed=chrm.size[,1:3], 
                                bigWig.path=bigWig.path, 
                                file.prefix=file.prefix,
                                file.suffix=file.suffix,
                                min.length=min.length)
all(names(reads.pks) == names(reads.tot))



# calculate FRiP
frip0 = colSums(reads.pks) / colSums(reads.tot)
tmap1 = c("t0", "20min", "40min", "60min", "2hr", "3hr", "4hr", "6d")
tmap2 = c(0, 20/60, 40/60, 1, 2, 3, 4, 6*24)
tmap = cbind(tmap1, tmap2) %>% as.data.frame(stringsAsFactors=F)
names(tmap) = c("condit", "hr")
frip = cbind(names(frip0), frip0) %>% as.data.frame(stringsAsFactors=F)
names(frip) = c("crep","FRiP")
condit = sapply(frip$crep,function(x)strsplit(x,"_")[[1]][1])
time = sapply(condit,function(x)tmap$hr[tmap$condit==x])
frip = frip %>% mutate(condit=condit, hr=time)
frip[,c(2,4)] = apply(frip[,c(2,4)],2,function(x){data.matrix(x)%>%as.numeric})

# summarize the FRiP data
summary <- frip %>% # the names of the new data frame and the data frame to be summarised
  group_by(hr) %>%   # the grouping variable
  summarise(mean = mean(FRiP),  # calculates the mean of each group
            sd = sd(FRiP), # calculates the standard deviation of each group
            n = n(),  # calculates the sample size per group
            SE = sd(FRiP)/sqrt(n()))
summary$hr = signif(summary$hr,2)
summary$hr = as.factor(summary$hr)

# plot FRiP data
ggplot(summary, aes(hr, mean)) + 
  geom_col() + ylab("mean FRiP") +
  geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width=0.2)