
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

# hist with all dTSS, single bin for dTSS=0 (Fig S4a)
plt = log10(dTSS+1)
del = min(plt[plt>0])
bk = seq(min(plt), max(plt)+del, by=del)
hist(plt, col="black",main="",xlab="dTSS, identified - largest (bp)",ylab="frequency",breaks=bk)

# hist with only re-annotated TSSs at 75% quantile (Fig S4b)
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

# all data hist (Fig S4c)
hist(dTTS/1000, col="black",main="",xlab="dTTS, identified - largest (kb)",ylab="frequency")

# subset hist at 75% quantile (Fig S4d)
hist(dTTS[abs(dTTS)<qt]/1000, col="black",main="",xlab="dTTS, identified - largest (kb)",ylab="frequency")

dev.off()


