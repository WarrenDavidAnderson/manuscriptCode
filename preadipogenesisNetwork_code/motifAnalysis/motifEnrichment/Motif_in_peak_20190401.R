
library(dplyr)
library(DESeq2)
library(ggplot2)
library(bigWig)
library(prodlim)
library(lattice)
library(stringr)


dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/atac_time_deg/motifEnrich")
setwd(dir)

############################################################
# load data from differential peak analysis
# load data from de novo motif analysis
############################################################

# load annotation for all peaks (summits.bed)
# see empiricalSummits_20190325.R
load("ATACsummits_20190325.RData")
all.peaks = paste0(summits.bed[,1],":",summits.bed[,2],"-",summits.bed[,3])
all.summits = paste0(summits.bed[,1],":",summits.bed[,2],
                     "-",summits.bed[,3],",",summits.bed[,4])

# load differential peak data
# see atac_time_deg_20190214.R
load("pairwise.deg.RData")

# load annotation for all peaks (bed.map)
# see atac_time_deg_20190214.R
load("bed.map.RData")
all.peaks = paste0(bed.map[,1],":",bed.map[,2],"-",bed.map[,3])

# functions to load and plot atac data
code.dir = paste0("/run/user/1001/gvfs/smb-share:server=home1.virginia.edu,share=cphgdesk/",
                  "users/wa3j/MGlab/Analysis/Adipogenesis/ATAC_time_DEG/peak_enrich")
source(paste0(code.dir,"/motifEnrich_functions.R"))

# load atac deseq2 data object with sizefactors, deseq_obj_preadip
# see atac_time_deg_20190214.R
load("atacDeseq.RData")
sizefac = sizeFactors(deseq_obj_preadip)

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

############################################################
# analysis parameters
############################################################

# preadipogenesis times
t.order = cbind(c("t0","20min","40min","60min","2hr","3hr","4hr"),
                c(0,0.33,0.67,1,2,3,4)) %>% 
  as.data.frame(stringsAsFactors=FALSE)
names(t.order) = c("condition","time")


# bigwig import params
file.prefix = "3T3_"
file.suffix = ".bigWig"
min.length = 1

# fdr threshold for filtering all fimo results per factor
fdr.thresh = 0.001
maxdiff.thresh = 10

# peak count threshold for considering optimal chi-sq result
pk.cnt.thresh = 100

# parameters for designating dynamic peaks
sig.thresh = 0.001
fc.thresh = 1

# parameters for designating non-dynamic peaks
sig.un = 0.5
fc.un = 0.25

# motif enrichment params 
step = 20
half.win = 600
delbp = 60

# smoothing param
loess.span = 0.3

# save.image("enrichdat_600bp_step20.RData")
# see enrichSummary.sh for analysis 

############################################################
# perform enrichment analyses for each factor
# see enrichSummary.sh for the implementation of this analysis on Rivanna
# enrich1-8.R does this analysis, too
# note that p-values are BH adjusted within factor in the .sh script in Rivanna
############################################################


dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/atac_time_deg/motifEnrich")
setwd(dir)

# data location for fimo bigWigs
# bigwigs produced in TFfimo216.sh
fimo.bigWig.path = paste0("/media/wa3j/Seagate2/Documents/PRO/adipogenesis",
                          "/July2018/atac_time_deg/motifEnrich/fimo_bigWig750k")

# load fac.chisq
# load("enrichdatProcessed500k.RData")


############################################################
# format, sort, and filter enrichment data
############################################################

# decide which time comparisons to consider based on peak counts
# note that nmotifUN is the number of unchanged peaks for the given time comparison
fac1 = fac.chisq[[1]]
fac1 = fac1[,c(3,16:18)]
indlo1 = which(fac1$npeakINC < 200)
indlo2 = which(fac1$npeakDEC < 200)
indlo = union(indlo1, indlo2)
fac1[indlo,]
fac1[-indlo,]

# look at relationship between variability and motif percentages
plt.dat = c()
for(ii in 1:length(fac.chisq)){
  new = fac.chisq[[ii]]
  plt.dat = rbind(plt.dat, new[-indlo,])
}
plot(plt.dat$INCminusUN, plt.dat$varINC)
plot(plt.dat$DECminusUN, plt.dat$varDEC)
hist(plt.dat$varINC %>% log)
hist(plt.dat$varDEC %>% log)

# separate into enrichment direction classes
enrich.thresh = 5
pkind.thresh = 0.5
enrich.inc = plt.dat %>% filter( (INCminusUN>enrich.thresh & (pkindINC-pkindUN)>pkind.thresh) | 
                                   (INCminusDEC>enrich.thresh & (pkindINC-pkindDEC)>pkind.thresh))
enrich.dec = plt.dat %>% filter((DECminusUN>enrich.thresh & (pkindDEC-pkindUN)>pkind.thresh) | 
                                  (DECminusINC>enrich.thresh & (pkindDEC-pkindINC)>pkind.thresh))
enrich.un = plt.dat %>% filter((DECminusUN < -enrich.thresh & (pkindUN-pkindDEC)>pkind.thresh) | 
                                 (INCminusUN < -enrich.thresh & (pkindUN-pkindINC)>pkind.thresh))

# filter based on enrichment peak index
enrich.pk.thresh = 1.8
enrich.inc = enrich.inc %>% filter(pkindINC > enrich.pk.thresh) %>% mutate(class="inc")
enrich.dec = enrich.dec %>% filter(pkindDEC > enrich.pk.thresh) %>% mutate(class="dec")
enrich.un = enrich.un %>% filter(pkindUN > enrich.pk.thresh) %>% mutate(class="un")

# check fdrs
max(enrich.inc$pval) # 2.762924e-11
max(enrich.dec$pval) # 0.002106126
max(enrich.un$pval)  # 0.001951505

############################################################
# aggregate and summarize enrichment data
# filter for plotting enrichment data for each factor
############################################################

# aggregate all data
processed.enrich.dat = rbind(enrich.inc, enrich.dec, enrich.un)
length(unique(processed.enrich.dat$id))
unique(processed.enrich.dat$TF)

# summarize data for each factor
enrich.summary = c()
for(ii in unique(processed.enrich.dat$id)){
  motif = ii
  ind = which(processed.enrich.dat$id == motif)
  TF = processed.enrich.dat$TF[ind][1]
  dat = processed.enrich.dat[ind,]
  minp = min(dat$pval)
  ind.inc = grep("inc",dat$class)
  ind.dec = grep("dec",dat$class)
  ind.unc = grep("un",dat$class)
  ninc = length(ind.inc)
  ndec = length(ind.dec)
  nunc = length(ind.unc)
  comp.inc = paste(dat$tcomp[ind.inc],collapse=",")
  comp.dec = paste(dat$tcomp[ind.dec],collapse=",")
  comp.unc = paste(dat$tcomp[ind.unc],collapse=",")
  new = c(TF,motif,minp,ninc,ndec,nunc,comp.inc,comp.dec,comp.unc)
  enrich.summary = rbind(enrich.summary, new)
}
enrich.summary = as.data.frame(enrich.summary,stringsAsFactors=F)
rownames(enrich.summary) = 1:nrow(enrich.summary)
names(enrich.summary) = c("factor","id","minFDR","ninc","ndec","nunc",
                          "compsINC","compsDEC","compsUNC")
#write.table(enrich.summary,"enrichSummary750k.txt",col.names=T,row.names=F,quote=F,sep="\t")

# prioritize data for plotting
# select the time comparison with the largest peak index
enrich.inc.plt = enrich.plot.select(enrich.dat=enrich.inc, comp.mode="inc")
enrich.dec.plt = enrich.plot.select(enrich.dat=enrich.dec, comp.mode="dec")
enrich.unc.plt = enrich.plot.select(enrich.dat=enrich.un, comp.mode="un")


############################################################
# plot enrichment data
############################################################

plot.enrich(comp.mode="inc",fname="enrichINCfiltered500k.pdf",mfrow=c(3,3),
            enrich.dat=enrich.inc.plt[order(-enrich.inc.plt$pkindINC),], 
            res.pairs=res.pairs, delbp=delbp, loess.span=loess.span,
            sig.thresh=sig.thresh, fc.thresh=fc.thresh,
            sig.un=sig.un, fc.un=fc.un,
            half.win=half.win, step=step,
            fimo.bigWig.path=fimo.bigWig.path,plt.trace=F)

plot.enrich(comp.mode="dec",fname="enrichDECfiltered500k.pdf",mfrow=c(3,3),
            enrich.dat=enrich.dec.plt[order(-enrich.dec.plt$pkindDEC),], 
            res.pairs=res.pairs, delbp=delbp, loess.span=loess.span,
            sig.thresh=sig.thresh, fc.thresh=fc.thresh,
            sig.un=sig.un, fc.un=fc.un,
            half.win=half.win, step=step,
            fimo.bigWig.path=fimo.bigWig.path,plt.trace=F)

plot.enrich(comp.mode="un",fname="enrichUNCfiltered500k.pdf",mfrow=c(3,3),
            enrich.dat=enrich.unc.plt[order(-enrich.unc.plt$pkindUN),], 
            res.pairs=res.pairs, delbp=delbp, loess.span=loess.span,
            sig.thresh=sig.thresh, fc.thresh=fc.thresh,
            sig.un=sig.un, fc.un=fc.un,
            half.win=half.win, step=step,
            fimo.bigWig.path=fimo.bigWig.path,plt.trace=F)


############################################################
# plot motifs and summarize motif clusters
############################################################

dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/atac_time_deg/motifEnrich")
setwd(dir)

# load data from motif_communities_TF_20190317.R
# load("TFanalysis_20190320.RData")

require(ggplot2)
require(ggseqlogo)
library(gridExtra)
library(dplyr)
library(seqLogo)
library(prodlim)
library(lattice)

# load motif plotting functions
code.dir = paste0("/run/user/1001/gvfs/smb-share:server=home1.virginia.edu,share=cphgdesk",
          "/users/wa3j/MGlab/Analysis/Adipogenesis/ATAC_time_DEG/TFcluster/motif_clust_functions.R")
source(code.dir)

# plot the motif signatures
final.plts = list()
p=1
for(id in enrich.summary$id){
  fac = enrich.summary$factor[enrich.summary$id==id]
  mat = flt.pwm.clust[[id]]
  final.plts[[p]] = logo.plt(pwm=mat, title=fac)
  p=p+1
}
pdf("TFs_filtered500k.pdf", onefile = TRUE, height=9, width=4.5)
marrangeGrob(grobs=final.plts, nrow=6, ncol=1, top=NULL)
dev.off()

# isolate filtered factors
motif.community5 = motif.community4 %>% mutate(comm = paste0("community_",community))
motif.community5 = motif.community5[motif.community5$motif %in% enrich.summary$id,]
cols = sapply(motif.community5$comm,function(x){which(names(datID)==x)})
TFspreadsheet = datTF[,cols]
# write.table(TFspreadsheet,"TFclusters750k.txt",col.names=T,row.names=F,sep="\t")
# write.table(motif.community5,"TFclusterAnnotation750k.txt",col.names=T,row.names=F,sep="\t")

# generate lists of pwm ids for each community
TF.list = list()
for(ii in 1:ncol(datTF)){
  mtfs = datID[,ii]
  TF.list[[names(datTF)[ii]]] = mtfs[mtfs != 0]
}

############################################################
# plot selected pwms
############################################################

# specify communities of interest
plt.comms = c(5,94,122,137,378)
plt.comms = c(25)
plt.comms = c(61,83)
plt.comms = paste0("community_",plt.comms)

# directory for all pwms
pwm.dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
                 "adipogenesis/July2018/atac_time_deg/communities/TF/TF679")

# function to get pwm
get.pwm = function(pwm.dir=NULL,motif=NULL){
  fname = paste0(pwm.dir,"/",motif,".txt")
  mat = read.table(fname,sep="\t",header=F,stringsAsFactors=F,fill=T) %>% t
  nas = apply(mat,1,function(x){all(is.na(x))})
  ind = which(nas==TRUE)
  if(length(ind)>0){mat = mat[-ind,]}
  rownames(mat) = c("A","C","G","T")
  colnames(mat) = c(1:ncol(mat))
  return(mat)
}

# plot the PWMs for all community members
comm.plts = list()
p=1
for(ii in 1:length(plt.comms)){
  comm = plt.comms[ii]
  mtfs = TF.list[[comm]]
  for(id in mtfs){
    print(paste(comm,", ",id))
    mat = get.pwm(pwm.dir=pwm.dir, motif=id)
    comm.plts[[p]] = logo.plt(pwm=mat, title=id)
    p=p+1
  }
}
pdf("pwm_mixed3.pdf", onefile = TRUE, height=9, width=4.5)
marrangeGrob(grobs=comm.plts, nrow=6, ncol=1, top=NULL)
dev.off()  



