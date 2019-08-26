
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
## check for genes with identical annotations, manually curate
#######################################################################################

# remove genes based on manual analysis
genes.remove1 = c(remove.genes.id$remove, fix.genes.id$upstream, fix.genes.id$downstream)
TSS.gene.filtered1 = TSS.gene[!(TSS.gene$gene %in% genes.remove1),]


#######################################################################################
## address non-identical overlaps
#######################################################################################

# run overlap analysis
overlap.data = gene.overlaps( bed = TSS.gene.filtered1 )
has.start.inside = overlap.data$has.start.inside
is.a.start.inside = overlap.data$is.a.start.inside

# document numbers and filter overlap data to remove specific gene annotations
overlap.genes = overlap.data$cases$gene %>% unique
genes.remove2 = overlap.genes[grep("AC",overlap.genes)]
genes.remove2 = c(genes.remove2, overlap.genes[grep("AL",overlap.genes)])
genes.remove2 = c(genes.remove2, overlap.genes[grep("Gm",overlap.genes)])
genes.remove2 = c(genes.remove2, overlap.genes[grep("Rik",overlap.genes)])

# identify genes with multiple starts inside (i.e., 'big' genes)
# set these genes for removal
mult.inside.starts = inside.starts(vec = is.a.start.inside$xy)
genes.remove2 = c(genes.remove2, mult.inside.starts) %>% unique

# remove filtered genes and re-run overlap analysis
in.dat = TSS.gene.filtered1[!(TSS.gene.filtered1$gene %in% genes.remove2),]
overlap.data = gene.overlaps( bed = in.dat )
has.start.inside = overlap.data$has.start.inside
is.a.start.inside = overlap.data$is.a.start.inside
case.dat = overlap.data$cases


# remove genes based on manual curation
genes.remove2 = c(genes.remove2, remove.genes.ov$remove, 
                  fix.genes.ov$upstream, fix.genes.ov$downstream)
TSS.gene.filtered2 = TSS.gene.filtered1[
  !(TSS.gene.filtered1$gene %in% genes.remove2),]


#######################################################################################
## address adjacent gene pairs that were manually identified
#######################################################################################

# get the coordinates for adjacent gene pairs
fix.genes = rbind(fix.genes.id, fix.genes.ov)
bp.bin = 5
knot.div = 40
shift.up = 100
delta.tss = 50
diff.tss = 1000
dist.from.start = 50
adjacent.coords = adjacent.gene.coords(fix.genes=fix.genes, bed.long=TSS.gene,
                                 exon1=gencode.firstExon,
                                 bw.plus=bw.plus, bw.minus=bw.minus,
                                 knot.div=knot.div, bp.bin=bp.bin,
                                 shift.up=shift.up, delta.tss=delta.tss,
                                 dist.from.start=dist.from.start,
                                 diff.tss=diff.tss, fname="adjacentSplines.pdf")

# aggregate downstream adjacent genes with main data
TSS.gene.filtered3 = rbind(TSS.gene.filtered2, adjacent.coords)

# verify the absence of overlaps
overlap.data = gene.overlaps( bed = TSS.gene.filtered3 )
overlap.data$cases # NULL

#######################################################################################
## isolate gene ends for TTS estimation
#######################################################################################

# get intervals for TTS evaluation
add.to.end = 100000
fraction.end = 0.2
dist.from.start = 50
bed.for.tts.eval = get.end.intervals(bed=TSS.gene.filtered3,
                                     add.to.end=add.to.end,
                                     fraction.end=fraction.end,
                                     dist.from.start=dist.from.start)

# distribution of clip distances
pdf("s2hist.pdf")
par(mfrow=c(2,2))
hist(bed.for.tts.eval$xy,xlab="clip distance (bp)",col="black",main="")
dev.off()






