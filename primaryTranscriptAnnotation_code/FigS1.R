
# fresh package install
library(devtools)
install_github("WarrenDavidAnderson/genomicsRpackage/primaryTranscriptAnnotation")

# package libraries
library(NMF)
library(dplyr)
library(bigWig)
library(pracma)
library(RColorBrewer)
library(primaryTranscriptAnnotation)

setwd("/media/wa3j/Seagate2/Documents/PRO/adipogenesis/package/SM/shortfigs")

#######################################################################################
## get merged read data from bigWigs
#######################################################################################

plus.file = "preadip_plus_merged.bigWig"
minus.file = "preadip_minus_merged.bigWig"
bw.plus = load.bigWig(plus.file)
bw.minus = load.bigWig(minus.file)


#######################################################################################
## get the largest interval for each gene
#######################################################################################

# get intervals for furthest TSS and TTS +/- interval
largest.interval.bed = get.largest.interval(bed=gencode.transcript)


#######################################################################################
## filter annotations to remove unexpressed genes
#######################################################################################

# get read counts and densities in annotated gene
transcript.reads = read.count.transcript(bed=gencode.transcript, 
                                         bw.plus=bw.plus, bw.minus=bw.minus)

# specify which genes to cut based on low expression, visualize cutoffs
den.cut = -6
cnt.cut = 3
ind.cut.den = which(log(transcript.reads$density) < den.cut)
ind.cut.cnt = which(log(transcript.reads$counts) < cnt.cut)
ind.cut = union(ind.cut.den, ind.cut.cnt)

# remove "unexpressed" genes
unexp = names(transcript.reads$counts)[ind.cut]
largest.interval.expr.bed = largest.interval.bed[
  !(largest.interval.bed$gene %in% unexp),]

#######################################################################################
## get a single TSS for each gene
#######################################################################################

# select the TSS for each gene and incorporate these TSSs into the largest interval coordinates
bp.range = c(20,120)
cnt.thresh = 5
bed.out = largest.interval.expr.bed
bed.in = gencode.firstExon[gencode.firstExon$gene %in% bed.out$gene,]
TSS.gene = get.TSS(bed.in=bed.in, bed.out=bed.out,
                   bw.plus=bw.plus, bw.minus=bw.minus,
                   bp.range=bp.range, cnt.thresh=cnt.thresh)

TSS.gene = TSS.gene$bed

#######################################################################################
## TSS performance analysis 
#######################################################################################

# parameters and analysis
window = 1000
bp.bin = 50
bed1 = TSS.gene
bed2 = largest.interval.expr.bed
bed2 = bed2[bed2$gene %in% bed1$gene,]
tss.eval = eval.tss(bed1=bed1, bed2=bed2, bin.thresh=2, 
                    bw.plus=bw.plus, bw.minus=bw.minus,
                    window=window, bp.bin=bp.bin, fname="TSSres.pdf")


infr = tss.eval$tss.dists.inf$dist
larg = tss.eval$tss.dists.lng$dist

# save.image("TSSeval.RData")

ylim = c(0,3000)
bk = seq(-window/2, window/2, bp.bin/2)
pdf("figs1.pdf"); par(mfrow=c(2,2))
hist(infr,main="inferred",xlab="dist from TSS to read max (bp)",
     breaks=bk,col="black",ylim=ylim)
hist(larg,main="largest",xlab="dist from TSS to read max (bp)",
     breaks=bk,col="black",ylim=ylim)
dev.off()

# run levene's test for variance differences
library(car)
var(infr)
var(larg)
labs = c(rep("infr",length(infr)),rep("larg",length(larg)))
Ldat = data.frame(vals=c(infr,larg),labs=as.factor(labs))
res = leveneTest(vals~labs,data=Ldat)
res$`Pr(>F)`
