

# directory
dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/atac_time_deg/motifEnrich")
setwd(dir)

# key packages
library(dplyr)
library(DESeq2)
library(bigWig)
library(stringr)

############################################################
## key data
############################################################

# load motif identification for enriched factors - tf.community.ann
# see communityTFs.R
load("TFcommunity.annotations.RData")

# load community/motif annotation information
fname = "TFclusters750k.txt"
tf.ids = read.table(fname,header=T,stringsAsFactors=F)
fname = "TFclusterAnnotation750k.txt"
motif.ctrs = read.table(fname,header=T,stringsAsFactors=F)

# load peak coordinates - bed.map
load("bed.map.RData")

# load enrichment data (Motif_in_peak_20190401.R)
load("enrichdat_600bp_step20.RData")

# load motif mapping counts from plotMotifDynamics_20190418.R
# pwm.peaks
load("pwm.peaks.RData")

# upload mapping between TF motif and database names
pathToFile <- paste0(getwd(),"/motif.id.key.txt")
f <- file(pathToFile, "rb")  
rawText <- readLines(f)
close(f)
motifs0 <- sapply(rawText, str_split, "\t")
motifs.raw = unlist(motifs0)[-c(1)] %>% unname
ids.new = sapply(motifs.raw,function(x){strsplit(x," ",fixed=TRUE)[[1]][1]}) %>% unname
ids.orig = sapply(motifs.raw,function(x){
  id = strsplit(x," ",fixed=TRUE)[[1]]
  out = paste(id[2:length(id)], collapse=" ")
}) %>% unlist
id.map = as.data.frame(cbind(ids.new,ids.orig),stringsAsFactors=F)
rownames(id.map) = 1:nrow(id.map)

# go to the fimo directory
fimo.bigWig.path = "/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/motifEnrich/fimo_bigWig750k"

############################################################
## aggregate annotations for community groupings
############################################################

# atac peak coordinates for read mapping
atac.bed = bed.map[,1:3]
rownames(atac.bed) = paste0(atac.bed$chr,":",atac.bed$start,"-",atac.bed$end)

# data storage
motif.mapping = as.data.frame(matrix(0,nrow(atac.bed),length(unique(pwm.peaks$maincom))))
rownames(motif.mapping) = rownames(atac.bed)

# loop through every main community
for( ii in 1:length(unique(pwm.peaks$maincom)) ){
  
  # main community
  maincom = unique(pwm.peaks$maincom)[ii]
  
  # find subcommunity with max fimo matches
  comm.pks = pwm.peaks[pwm.peaks$maincom==maincom,1:3]
  ind.max = which(comm.pks$npeaks == max(comm.pks$npeaks))
  maxcom = comm.pks$comm[ind.max]
  
  # find the motif corresponding to the max subcommunity
  tf.id = motif.ctrs$ids.orig[motif.ctrs$comm == maxcom]
  motif.id = motif.ctrs$motif[motif.ctrs$comm == maxcom]
  
  # load fimo data
  fimo.dat = paste0(fimo.bigWig.path,"/fimo_",motif.id,".bigWig")
  bw.fimo = load.bigWig(fimo.dat)
  
  # load bigwig factor mapping data into peak coordinates
  # number of TFBS in each peak region
  motif.count = bed.region.bpQuery.bigWig(bw.fimo, atac.bed)
  names(motif.count) = rownames(atac.bed)
  
  # store motif counts for each factor
  motif.mapping[,ii] = motif.count
  names(motif.mapping)[ii] = tf.id
  
} # ii, main community loop

# save(motif.mapping, file="peak.annotation.20180423.RData")
# load("peak.annotation.20180423.RData")

############################################################
## plot peak summary data
############################################################

library(NMF)

# remove peaks with out motifs
max(rowSums(motif.mapping))
min(rowSums(motif.mapping))
ind.rem = which(rowSums(motif.mapping)==0)
plt = motif.mapping[-ind.rem,]

# log scale and perform PCA
plt.log0 = log10( plt + 1 )

# perform PCA for TFs
pca = prcomp(t(plt.log0), center=TRUE, scale.=FALSE)
scores = pca$x
loadings = pca$rotation
percent = pca$sdev^2 / sum( pca$sdev^2 )
plot(scores[,1], scores[,2])

# perform PCA for ATAC peaks
pca = prcomp(plt.log0, center=TRUE, scale.=FALSE)
scores = pca$x
loadings = pca$rotation
percent = pca$sdev^2 / sum( pca$sdev^2 )
plot(scores[,1], scores[,2])

library(heatmap3)
pdf("PeakHeatmap1.pdf")
heatmap3(plt.log,useRaster=TRUE)
dev.off()

# pdf("PeakHeatmap1.pdf")
# aheatmap(plt.log)
# dev.off()




