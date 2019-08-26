
# package libraries
library(NMF)
library(dplyr)
library(bigWig)
library(pracma)
library(RColorBrewer)
library(primaryTranscriptAnnotation)

setwd("/media/wa3j/Seagate2/Documents/PRO/adipogenesis/package/SM/shortfigs")

#######################################################################################
## ppar fig
#######################################################################################

# see FigS3ab.R
# load("coordsFigS3.RData")

plus.file = "preadip_plus_merged.bigWig"
minus.file = "preadip_minus_merged.bigWig"
bw.plus = load.bigWig(plus.file)
bw.minus = load.bigWig(minus.file)

# plot new coord id
pdf("pparplt.pdf")
par(mfrow=c(2,2))
gene.end = bed.for.tts.eval
long.gene = largest.interval.bed
tts.plot(coords=coords, gene.end=gene.end, long.gene=long.gene,
         gene = "Pparg", xper=0.1, yper=0.2,
         bw.plus=bw.plus, bw.minus=bw.minus, bp.bin=5,
         frac.min=frac.min, frac.max=frac.max,
         add.to.end=add.to.end, tau.dist=tau.dist)
dev.off()
