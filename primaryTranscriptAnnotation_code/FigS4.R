
# package libraries
library(NMF)
library(dplyr)
library(bigWig)
library(pracma)
library(RColorBrewer)
library(primaryTranscriptAnnotation)

setwd("/media/wa3j/Seagate2/Documents/PRO/adipogenesis/package/SM/shortfigs")
load("coordsFigS3.RData")


#######################################################################################
## get the largest interval and inferred annotations, and exon 1 data
#######################################################################################

bed.long = largest.interval.bed
bed.tss.tts = coords

# note. this annotation set is part of the primaryTranscriptAnnotation package
head(gencode.firstExon)


# check how many genes have a single first exon
# number of exon1s for dTSS=0 cases
nexon1all = sapply(bed.tss.tts$gene,function(x){
  inds = which(gencode.firstExon$gene==x)
  strd = gencode.firstExon$strand[inds][1]
  if(strd=="+"){e1 = gencode.firstExon$start[inds]}
  if(strd=="-"){e1 = gencode.firstExon$end[inds]}
  out = length(unique(e1))
  return(out)
})
length(which(nexon1all>1))
length(which(nexon1all==1))
dim(bed.tss.tts)

#######################################################################################
## compare with inferred coords
#######################################################################################

# separate by strand
long.plus = bed.long %>% filter(strand=="+")
long.minus = bed.long %>% filter(strand=="-")
infr.plus = bed.tss.tts %>% filter(strand=="+")
infr.minus = bed.tss.tts %>% filter(strand=="-")

# organize data
ind = sapply(infr.plus$gene,function(x){which(long.plus$gene==x)}) %>% unlist
long.plus = long.plus[ind,]
ind = sapply(infr.minus$gene,function(x){which(long.minus$gene==x)}) %>% unlist
long.minus = long.minus[ind,]
all(infr.minus$gene == long.minus$gene)
all(infr.plus$gene == long.plus$gene)

# 
# #######################################################################################
# ## compare with inferred coords - TSS - panel c
# #######################################################################################
# 
# # TSS dists (pos = downstream for id)
# d1 = infr.plus$start - long.plus$start
# d2 = long.minus$end - infr.minus$end
# dTSS = c(d1,d2)
# names(dTSS) = c(infr.plus$gene, infr.minus$gene)
# 
# # get number of exon1 instances for each dTSS
# nex1 = sapply(names(dTSS),function(x){length(which(gencode.firstExon$gene==x))})
# frac.rean = rep(0, length(1:max(nex1)))
# for(ii in 1:max(nex1)){
#   ind.nex1 = which(nex1 == ii)
#   ge.nex1 = names(nex1)[ind.nex1]
#   ind.tss = sapply(ge.nex1,function(x){which(names(dTSS)==x)}) %>% unlist()
#   if(length(ind.tss)>0){
#     frac.rean[ii] = length(which(dTSS[ind.tss]!=0)) / length(ind.nex1) 
#   }
# }
# frac.rean = data.frame(nex1=1:max(nex1),frac.rean=frac.rean)
# hypo.rean = data.frame(nex1=1:max(nex1),frac.rean=1/(1:max(nex1)))
# df.bar = barplot(frac.rean$frac.rean, names.arg=frac.rean$nex1)
# lines(x=df.bar, y=hypo.rean$frac.rean,lwd=2)


#######################################################################################
## compare with inferred coords - TSS
#######################################################################################

pdf("annDiff.pdf")
par(mfrow=c(2,2))

# TSS dists (pos = downstream for id)
d1 = infr.plus$start - long.plus$start
d2 = long.minus$end - infr.minus$end
dTSS = c(d1,d2)
names(dTSS) = c(infr.plus$gene, infr.minus$gene)

# assess extreme case
ii = which(abs(d2)==max(abs(d2)))
long.minus[ii,]
infr.minus[ii,]

# genes for dTSS=0
ind0 = which(dTSS==0)
genes0 = names(dTSS)[ind0]

# number of exon1s for dTSS=0 cases
nexon1 = sapply(genes0,function(x){
  inds = which(gencode.firstExon$gene==x)
  strd = gencode.firstExon$strand[inds][1]
  if(strd=="+"){e1 = gencode.firstExon$start[inds]}
  if(strd=="-"){e1 = gencode.firstExon$end[inds]}
  out = length(unique(e1))
  return(out)
})

# remove cases with dTSS=0 and one exon
ind.rem = which(nexon1==1)
dTSSnoe1 = dTSS[-ind0[ind.rem]]
dTSS1e1 = dTSS[-ind0[-ind.rem]]

# hist with single bin for dTSS=0
plt = log10(dTSS+1)
del = min(plt[plt>0])
bk = seq(min(plt), max(plt)+del, by=del)
hist(plt, col="black",main="",xlab="dTSS, identified - largest (bp)",
     ylab="frequency",breaks=bk, ylim=c(0,8000))

# hist with single exon for dTSS=0 in first bin (49%)
# gray bar for panel a
plt = log10(dTSS1e1+1)
del = min(plt[plt>0])
bk = seq(min(plt), max(plt)+del, by=del)
hist(plt, col="black",main="",xlab="dTSS, identified - largest (bp)",
     ylab="frequency",breaks=bk, ylim=c(0,8000))


# hist with only re-annotated TSSs
nexon1 = sapply(names(dTSS),function(x){
  inds = which(gencode.firstExon$gene==x)
  strd = gencode.firstExon$strand[inds][1]
  if(strd=="+"){e1 = gencode.firstExon$start[inds]}
  if(strd=="-"){e1 = gencode.firstExon$end[inds]}
  out = length(unique(e1))
  return(out)
})
ind1 = which(nexon1 > 1)
ind2 = which(dTSS > 0)
plt = dTSS[intersect(ind1, ind2)]
del = 80
qt = quantile(plt)
plt = plt[plt < qt[4]]
bk = seq(min(plt), max(plt)+del, by=del)
hist(plt, col="black",main="",xlab="dTSS, identified - largest (bp)",ylab="frequency")

# get number of exon1 instances for each dTSS
nex1 = sapply(names(dTSS),function(x){length(which(gencode.firstExon$gene==x))})
frac.rean = rep(0, length(1:max(nex1)))
for(ii in 1:max(nex1)){
  ind.nex1 = which(nex1 == ii)
  ge.nex1 = names(nex1)[ind.nex1]
  ind.tss = sapply(ge.nex1,function(x){which(names(dTSS)==x)}) %>% unlist()
  if(length(ind.tss)>0){
    frac.rean[ii] = length(which(dTSS[ind.tss]!=0)) / length(ind.nex1) 
  }
}
frac.rean = data.frame(nex1=1:max(nex1),frac.rean=frac.rean)
hypo.rean = data.frame(nex1=1:max(nex1),frac.rean=1/(1:max(nex1)))
frac.rean = frac.rean[-1,]
hypo.rean = hypo.rean[-1,]
hypo.rean$frac.rean = 1 - hypo.rean$frac.rean
df.bar = barplot(frac.rean$frac.rean, names.arg=frac.rean$nex1)
lines(x=df.bar, y=hypo.rean$frac.rean,lwd=2)

#######################################################################################
## Look at the number of exon1 for cases of dTSS=0
#######################################################################################

# genes for dTSS=0
ind0 = which(dTSS==0)
genes0 = names(dTSS)[ind0]

# number of exon1s for dTSS=0 cases
nexon1 = sapply(genes0,function(x){
  inds = which(gencode.firstExon$gene==x)
  strd = gencode.firstExon$strand[inds][1]
  if(strd=="+"){e1 = gencode.firstExon$start[inds]}
  if(strd=="-"){e1 = gencode.firstExon$end[inds]}
  out = length(unique(e1))
  return(out)
})

# hist(nexon1,col="black",breaks=c(0:26))

length(which(nexon1==1))
length(which(nexon1>1))
length(nexon1)

#######################################################################################
## compare with inferred coords - TTS
#######################################################################################

# TTS dists (pos = downstream for id)
d1 = infr.plus$end - long.plus$end
d2 = long.minus$start - infr.minus$start
dTTS = c(d1,d2)
qt = quantile(abs(dTTS), probs = 0.75)
length(which(abs(dTTS)>100000))

# assess extreme case
ii = which(abs(d2)==max(abs(d2)))
long.minus[ii,]
infr.minus[ii,]

ii = which(abs(d1)==max(abs(d1)))
long.plus[ii,]
infr.plus[ii,]

# all data hist
hist(dTTS/1000, col="black",main="",xlab="dTTS, identified - largest (kb)",ylab="frequency")

# subset hist
hist(dTTS[abs(dTTS)<qt]/1000, col="black",main="",xlab="dTTS, identified - largest (kb)",ylab="frequency")

dev.off()


