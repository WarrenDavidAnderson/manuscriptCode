
############################################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Figure 4f
############################################################################

library(dplyr)
library(ggplot2)
library(gridExtra)

setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5")

# eqtl & deg
deg.eqtl.genes = c("CCDC3","CLIC6","FADS1","GLDN","HSPA12A","MAP1B","MLPH","MMD","MYOT","NDRG4","NEO1","PDZD2","TBC1D9")

############################################################################
# integrate with ATAC data, see shell script fig4_atac.sh
############################################################################

setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5")

# peer associations with ld data and hg19 coords
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/CIdata.ldsnps.hg19.RData"
load(fname)

# output data for command line analysis
out = CIdata.ldsnps.hg19[,c(11:13,1:10,14:15)]
out = out %>% filter(R2 > 0.8)
fname = "espn.ld.hg19.bed"
write.table(out,fname,col.names=F,row.names=F,sep="\t",quote=F)

# directory for intersect results
datdir = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5/atac/"

# load intersect results
fname = paste0(datdir,"overlap_Inter.bed")
ovr.int = read.table(fname,header=F,stringsAsFactors=F)[,-c(1:3)]
fname = paste0(datdir,"overlap_preadip.bed")
ovr.pre = read.table(fname,header=F,stringsAsFactors=F)[,-c(1:3)]
fname = paste0(datdir,"overlap_adip.bed")
ovr.adp = read.table(fname,header=F,stringsAsFactors=F)[,-c(1:3)]
fname = paste0(datdir,"overlap_All.bed")
ovr.all = read.table(fname,header=F,stringsAsFactors=F)[,-c(1:3)]
names(ovr.int)=names(ovr.pre)=names(ovr.adp)=names(ovr.all)=names(out)

# document count of association (cases) loc
tot = length(unique(out$gene))
tot = 2408

# document fraction of overlapping loci for each atac peak overlap class
frac.pre = 100 * length(unique(ovr.pre$gene)) / tot
frac.adp = 100 * length(unique(ovr.adp$gene)) / tot
frac.all = 100 * length(unique(ovr.all$gene)) / tot
frac.int = 100 * length(unique(ovr.int$gene)) / tot


############################################################################
# examine enrichment of TFBSs overlapping esnps in open chromatin
############################################################################

# load data for all tfbs overlapping ld snps (R2 > 0.8)
load("ld.tfbs.RData")

# preadip peaks
nqtl = length(unique(out$gene))
inTFinPK = merge(ld.tfbs, ovr.pre, by=c("chr","start","end","gene"))
inTF.inPK = length(unique(inTFinPK$gene))
inTF.outPK = out[!(out$gene %in% inTFinPK$gene),] 
inTF.outPK = merge(inTF.outPK, ld.tfbs, by=c("chr","start","end","gene"))
inTF.outPK = length(unique(inTF.outPK$gene))
outTF.inPK = merge(out, ovr.pre, by=c("chr","start","end","gene"))
outTF.inPK = outTF.inPK[!(outTF.inPK$gene %in% ld.tfbs$gene),] 
outTF.inPK = length(unique(outTF.inPK$gene))
outTF.outPK = nqtl - (inTF.inPK + inTF.outPK + outTF.inPK)
FETtable = matrix(c(inTF.inPK, outTF.inPK, inTF.outPK, outTF.outPK),
                     nrow = 2, 
                     dimnames = list(TF = c("in", "out"),
                                     PK = c("in", "out")))
FETres.pre = fisher.test(FETtable, alternative = "two.sided")
100 * inTF.inPK / nqtl

# adipocyte/adipose peaks
nqtl = length(unique(out$gene))
inTFinPK = merge(ld.tfbs, ovr.adp, by=c("chr","start","end","gene"))
inTF.inPK = length(unique(inTFinPK$gene))
inTF.outPK = out[!(out$gene %in% inTFinPK$gene),] 
inTF.outPK = merge(inTF.outPK, ld.tfbs, by=c("chr","start","end","gene"))
inTF.outPK = length(unique(inTF.outPK$gene))
outTF.inPK = merge(out, ovr.adp, by=c("chr","start","end","gene"))
outTF.inPK = outTF.inPK[!(outTF.inPK$gene %in% ld.tfbs$gene),] 
outTF.inPK = length(unique(outTF.inPK$gene))
outTF.outPK = nqtl - (inTF.inPK + inTF.outPK + outTF.inPK)
FETtable = matrix(c(inTF.inPK, outTF.inPK, inTF.outPK, outTF.outPK),
                  nrow = 2, 
                  dimnames = list(TF = c("in", "out"),
                                  PK = c("in", "out")))
FETres.adp = fisher.test(FETtable, alternative = "two.sided")
100 * inTF.inPK / nqtl

# all peaks
nqtl = length(unique(out$gene))
inTFinPK = merge(ld.tfbs, ovr.all, by=c("chr","start","end","gene"))
inTF.inPK = length(unique(inTFinPK$gene))
inTF.outPK = out[!(out$gene %in% inTFinPK$gene),] 
inTF.outPK = merge(inTF.outPK, ld.tfbs, by=c("chr","start","end","gene"))
inTF.outPK = length(unique(inTF.outPK$gene))
outTF.inPK = merge(out, ovr.all, by=c("chr","start","end","gene"))
outTF.inPK = outTF.inPK[!(outTF.inPK$gene %in% ld.tfbs$gene),] 
outTF.inPK = length(unique(outTF.inPK$gene))
outTF.outPK = nqtl - (inTF.inPK + inTF.outPK + outTF.inPK)
FETtable = matrix(c(inTF.inPK, outTF.inPK, inTF.outPK, outTF.outPK),
                  nrow = 2, 
                  dimnames = list(TF = c("in", "out"),
                                  PK = c("in", "out")))
FETres.all = fisher.test(FETtable, alternative = "two.sided")
100 * inTF.inPK / nqtl

# intersect peaks
nqtl = length(unique(out$gene))
inTFinPK = merge(ld.tfbs, ovr.int, by=c("chr","start","end","gene"))
inTF.inPK = length(unique(inTFinPK$gene))
inTF.outPK = out[!(out$gene %in% inTFinPK$gene),] 
inTF.outPK = merge(inTF.outPK, ld.tfbs, by=c("chr","start","end","gene"))
inTF.outPK = length(unique(inTF.outPK$gene))
outTF.inPK = merge(out, ovr.int, by=c("chr","start","end","gene"))
outTF.inPK = outTF.inPK[!(outTF.inPK$gene %in% ld.tfbs$gene),] 
outTF.inPK = length(unique(outTF.inPK$gene))
outTF.outPK = nqtl - (inTF.inPK + inTF.outPK + outTF.inPK)
FETtable = matrix(c(inTF.inPK, outTF.inPK, inTF.outPK, outTF.outPK),
                  nrow = 2, 
                  dimnames = list(TF = c("in", "out"),
                                  PK = c("in", "out")))
FETres.int = fisher.test(FETtable, alternative = "two.sided")
100 * inTF.inPK / nqtl

# adjust for multiple testing
pval = c(FETres.pre$p.value, FETres.adp$p.value, FETres.all$p.value, FETres.int$p.value)
fdr = p.adjust(pval, method="BH")

############################################################################
# peak selection for plotting in the browser
############################################################################

# data
inTFinPK = merge(ld.tfbs, ovr.all, by=c("chr","start","end","gene", "genename", "case", "ref", "alt"))
inTFinPK = inTFinPK %>% filter(R2.x > 0.8)

# list all overlapping TFs
unique(inTFinPK$tfs)

# check key genes in atac loci, for text
ovr.adp[which(ovr.adp$genename==deg.eqtl.genes[13]),]


# exploratory loci
ind = grep("FADS1",inTFinPK$genename)
inTFinPK[ind,]
ind = grep("NEO1",inTFinPK$genename)
inTFinPK[ind,]
ind = grep("MYOT",inTFinPK$genename)
inTFinPK[ind,]

# exploratory TFs
ind = grep("PPARG",inTFinPK$tfs)
loc.ppg = inTFinPK[ind,]
length(unique(loc.ppg$gene))
ind = grep("KLF5",inTFinPK$tfs)
loc.klf = inTFinPK[ind,]
length(unique(loc.klf$gene))

# load and merge case 1 data
load("ci.deg.case1.filtered.RData")
case1.atac = merge(ci.deg.case1.filtered, inTFinPK, by=c("tfs","ntfs","gene","leadSNP","genename","ref","alt"))
length(unique(case1.atac$gene))

# single case 1 in atac peaks
case1.atac


# EGR1 - example case 1 mechanism
ci.deg.case1.filtered[12,]

# Table S12
ts12 = inTFinPK[,-c(6,11,13,20:26)]
names(ts12)[c(4,5,8,10:12)] = c("ensid","gene","ldSNP","pvalue","beta","R2")
write.table(ts12,"TableS12.txt",col.names=T,row.names=F,quote=F,sep="\t")

# check key genes in the intersect
inTFinPK = merge(ld.tfbs, ovr.all, by=c("chr","start","end","gene", "genename", "case", "ref", "alt"))
inTFinPK = inTFinPK %>% filter(R2.x > 0.8)

# list all overlapping TFs
unique(inTFinPK$tfs)

# check key genes in atac loci
ovr.adp[which(ovr.adp$genename=="FADS1"),]
ovr.adp[which(ovr.adp$genename=="NEO1"),]
ovr.adp[which(ovr.adp$genename=="MYOT"),]


############################################################################
# plot motif signatures
############################################################################

require(ggseqlogo)
library(gridExtra)
library(dplyr)
library(seqLogo)

# get pwd
pwm0 = read.table("PWM",header=F,sep=" ")[c(2,4,6,8)]
pwm0 = read.table("PWM2",header=F,sep=" ")[c(2,4,6,8)]
pwm0 = read.table("PWM3",header=F,sep=" ")[c(2,4,6,8)]

# function to get the reverse complement for a matrix
# this function assumes rows ACGT as input
revcomp = function(mat=NULL){
  targetRC = matrix(0,nrow(mat), ncol(mat))
  targetRC[1,] = rev(mat[4,])
  targetRC[2,] = rev(mat[3,])
  targetRC[3,] = rev(mat[2,])
  targetRC[4,] = rev(mat[1,])
  rownames(targetRC) = rownames(mat)
  return(targetRC)
}

# function to plot motif signature
logo.plt = function(pwm=NULL, title=NULL){
  ggplot() + geom_logo(pwm) + ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ylim(c(0,2.1))
} # logo.plt

# make plot, pparg
mat = t(pwm0)
rownames(mat) = c("A","C","G","T")
colnames(mat) = c(1:ncol(mat))
pdf("ppargmotif.pdf")
plt = logo.plt(mat,"PPARG_RXRA")
print(plt)
dev.off()

# make plot, nrf1
mat = revcomp( t(pwm0) )
rownames(mat) = c("A","C","G","T")
colnames(mat) = c(1:ncol(mat))
pdf("nrf1motif.pdf")
plt = logo.plt(mat,"NRF1")
print(plt)
dev.off()



