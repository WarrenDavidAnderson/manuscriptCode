
library(dplyr)
library(DESeq2)
library(ggplot2)
library(bigWig)
library(prodlim)
library(lattice)
library(stringr)

# functions to load and plot atac data
source("motifEnrich_functions.R")

dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/atac_time_deg/motifEnrich/dreg")
setwd(dir)

############################################################
# load data from differential dreg peak analysis
# load data from de novo motif analysis
############################################################

# load differential peak data
# see dREGanalysis.R (res.pairs list)
load("pairwise.deg.RData")

# load annotation for all peaks (bed.map)
# see dREGanalysis.R
load("bed.mapDreg.RData")

# load atac deseq2 data object with sizefactors, deseq_obj_preadip
# see dREGanalysis.R
load("dregDeseq.RData")
sizefac = sizeFactors(deseq_obj_preadip)

# upload mapping between TF motif and database names
# from allTFDB.sh
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

# parameters for designating dynamic dreg peaks
sig.thresh = 0.1
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
# see enrichSummary_diffDREG.sh for analysis 


############################################################
# perform enrichment analyses for each factor
# see enrichSummary_diffDREG.sh for the implementation of this analysis on Rivanna
# note that p-values are BH adjusted within factor in the .sh script in Rivanna
############################################################

library(dplyr)
library(DESeq2)
library(pracma)
library(ggplot2)
library(bigWig)
library(prodlim)
library(lattice)
library(stringr)

dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/atac_time_deg/motifEnrich/dreg")
setwd(dir)

# data location for fimo bigWigs
# bigwigs produced in TFfimo.sh
# from=wa3j@interactive.hpc.virginia.edu:/scratch/wa3j/TFfimo/f750000/*.bigWig*
# to=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/motifEnrich/fimo_bigWig750k
# scp -r $from $to
fimo.bigWig.path = paste0("/media/wa3j/Seagate2/Documents/PRO/adipogenesis",
                          "/July2018/atac_time_deg/motifEnrich/fimo_bigWig750k")

# load fac.chisq
load("enrichdat_600bp_step20.RData")
load("enrichdatProcessed750k_dreg.RData")

# functions to load and plot atac data
source("/run/user/1001/gvfs/smb-share:server=home1.virginia.edu,share=cphgdesk/users/wa3j/MS/adipogenesis/vignettes/ATAC_analysis/motifEnrich/code/motifEnrich_functions.R")


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

# look at relationship between variability and motif percentage dfferences
plt.dat = c()
for(ii in 1:length(fac.chisq)){
  new = fac.chisq[[ii]]
  plt.dat = rbind(plt.dat, new[-indlo,])
}
pdf("var_cor_diff.pdf");par(mfrow=c(2,2))
plot(plt.dat$INCminusUN, plt.dat$varINC)
plot(plt.dat$DECminusUN, plt.dat$varDEC)
hist(plt.dat$varINC %>% log)
hist(plt.dat$varDEC %>% log)
dev.off()

# separate into enrichment direction classes
enrich.thresh = 7 # threshold for percentage differences
enrich.inc = plt.dat %>% filter( INCminusUN>enrich.thresh | INCminusDEC>enrich.thresh )
enrich.dec = plt.dat %>% filter( DECminusUN>enrich.thresh | DECminusINC>enrich.thresh )
enrich.un = plt.dat %>% filter( DECminusUN < -enrich.thresh | INCminusUN < -enrich.thresh )

# filter based on enrichment peak index
enrich.pk.thresh = 1.2
enrich.inc = enrich.inc %>% filter(pkindINC > enrich.pk.thresh) %>% mutate(class="inc")
enrich.dec = enrich.dec %>% filter(pkindDEC > enrich.pk.thresh) %>% mutate(class="dec")
enrich.un = enrich.un %>% filter(pkindUN > enrich.pk.thresh) %>% mutate(class="un")

# filter based on peak density difference (fraction of range)
# e.g., max(up-dn, up-un) > pk.den.diff * (up-peak - up-flank)
pk.den.diff = 0.3
enrich.inc = enrich.inc %>% filter(dratUP > pk.den.diff)
enrich.dec = enrich.dec %>% filter(dratDN > pk.den.diff)
enrich.un = enrich.un %>% filter(dratUN > pk.den.diff)

# filter based on the baseline shifted integral
integ.thresh = 13
enrich.inc = enrich.inc %>% filter(integUP > integ.thresh)
enrich.dec = enrich.dec %>% filter(integDN > integ.thresh)
enrich.un = enrich.un %>% filter(integUN > integ.thresh)

# filter based on number of peaks with motifs
npk.thresh = 300
enrich.inc = enrich.inc %>% filter(nmotINC > npk.thresh)
enrich.dec = enrich.dec %>% filter(nmotDEC > npk.thresh)
enrich.un = enrich.un %>% filter(nmotUN > npk.thresh)

# > 300 peaks with motifs for the class with the max diff percentage
enrich.inc = filt.maxdiff.pks(dat=enrich.inc, npk.thresh=npk.thresh)
enrich.dec = filt.maxdiff.pks(dat=enrich.dec, npk.thresh=npk.thresh)
enrich.un = filt.maxdiff.pks(dat=enrich.un, npk.thresh=npk.thresh)

# check fdrs
max(enrich.inc$pval) # 1.075155e-10
max(enrich.dec$pval) # 9.608468e-26
max(enrich.un$pval)  # 1.405724e-08

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
write.table(enrich.summary,"enrichSummary750k_dreg.txt",col.names=T,row.names=F,quote=F,sep="\t")

# prioritize data for plotting
# select the time comparison with the largest peak index
enrich.inc.plt = enrich.plot.select(enrich.dat=enrich.inc, comp.mode="inc")
enrich.dec.plt = enrich.plot.select(enrich.dat=enrich.dec, comp.mode="dec")
enrich.unc.plt = enrich.plot.select(enrich.dat=enrich.un, comp.mode="un")


############################################################
# plot enrichment data
############################################################

# add coordinate ids for dreg regions to the bed.map
coord.id = paste0(bed.map$chr,":",bed.map$start,"-",bed.map$end)
bed.map = bed.map %>% mutate(id = coord.id)
chrs = paste0("chr",1:19)
bed.map = bed.map[bed.map$chr %in% chrs,]

plot.enrich(comp.mode="inc",fname="enrichINCfiltered750k_dreg.pdf",mfrow=c(3,3),
            enrich.dat=enrich.inc.plt[order(-enrich.inc.plt$pkindINC),], 
            res.pairs=res.pairs, delbp=delbp, loess.span=loess.span,
            sig.thresh=sig.thresh, fc.thresh=fc.thresh,
            sig.un=sig.un, fc.un=fc.un, type="dreg",
            half.win=half.win, step=step,
            fimo.bigWig.path=fimo.bigWig.path,plt.trace=F)

plot.enrich(comp.mode="dec",fname="enrichDECfiltered750k.pdf",mfrow=c(3,3),
            enrich.dat=enrich.dec.plt[order(-enrich.dec.plt$pkindDEC),], 
            res.pairs=res.pairs, delbp=delbp, loess.span=loess.span,
            sig.thresh=sig.thresh, fc.thresh=fc.thresh,
            sig.un=sig.un, fc.un=fc.un, type="dreg",
            half.win=half.win, step=step,
            fimo.bigWig.path=fimo.bigWig.path,plt.trace=F)

plot.enrich(comp.mode="un",fname="enrichUNCfiltered750k.pdf",mfrow=c(3,3),
            enrich.dat=enrich.unc.plt[order(-enrich.unc.plt$pkindUN),], 
            res.pairs=res.pairs, delbp=delbp, loess.span=loess.span,
            sig.thresh=sig.thresh, fc.thresh=fc.thresh,
            sig.un=sig.un, fc.un=fc.un, type="dreg",
            half.win=half.win, step=step,
            fimo.bigWig.path=fimo.bigWig.path,plt.trace=F)


############################################################
# plot motifs and summarize motif clusters
############################################################

dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/atac_time_deg/motifEnrich/dreg")
setwd(dir)

# load data from motif_communities_TF.R
load("TFanalysis_20191021.RData")

require(ggplot2)
require(ggseqlogo)
library(gridExtra)
library(dplyr)
library(seqLogo)
library(prodlim)
library(lattice)

# load motif plotting functions
code.dir = paste0("/run/user/1001/gvfs/smb-share:server=home1.virginia.edu,share=cphgdesk",
                  "/users/wa3j/MS/adipogenesis/vignettes/ATAC_analysis/motifEnrich/code/motif_clust_functions.R")
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
pdf("TFs_filtered750k_dreg.pdf", onefile = TRUE, height=9, width=4.5)
marrangeGrob(grobs=final.plts, nrow=6, ncol=1, top=NULL)
dev.off()

# isolate filtered factors
# generate spreadsheets for manual analysis of TF families
motif.community5 = motif.community4 %>% mutate(comm = paste0("community_",community))
motif.community5 = motif.community5[motif.community5$motif %in% enrich.summary$id,]
cols = sapply(motif.community5$comm,function(x){which(names(datID)==x)})
TFspreadsheet = datTF[,cols]
nfac = apply(TFspreadsheet,2,function(x){length(which(x !=0))})
TFspreadsheet = TFspreadsheet[,order(nfac,decreasing=T)]
for(jj in 1:ncol(TFspreadsheet)){
  comm = names(TFspreadsheet)[jj]
  ind = which(motif.community5$comm == comm)
  ctr = motif.community5$ids.orig[ind]
  top = TFspreadsheet[1,jj]
  ind.ctr = which(TFspreadsheet[,jj] == ctr)
  TFspreadsheet[1,jj] = ctr
  TFspreadsheet[ind.ctr,jj] = top
}
write.table(TFspreadsheet,"TFclusters750k_dreg.txt",col.names=T,row.names=F,sep="\t")
write.table(motif.community5,"TFclusterAnnotation750k_dreg.txt",col.names=T,row.names=F,sep="\t")

# generate lists of pwm ids for each community
TF.list = list()
for(ii in 1:ncol(datTF)){
  mtfs = datID[,ii]
  TF.list[[names(datTF)[ii]]] = mtfs[mtfs != 0]
}


############################################################
# generate comprehensive summary table
############################################################

# comps for each TF motif
head(motif.community.summary)

# filter for enriched TFs
summary.atacDEG = c()
summary1 = motif.community.summary[motif.community.summary$id %in% motif.community5$motif,]
for(ii in 1:nrow(summary1)){
  left = summary1[ii,1:3]
  right = summary1[ii,4]
  cmps = strsplit(right,",") %>% unlist
  cnds = sapply(cmps,function(x){strsplit(x,"_")[[1]][1]})
  ldeg = length(which(cnds == "deg"))
  ldrg = length(which(cnds == "dreg"))
  ldyn = length(which(cnds == "dyn"))
  new = c(left, ldeg, ldrg, ldyn) %>% unlist
  summary.atacDEG = rbind(summary.atacDEG, new)
}
summary.atacDEG = as.data.frame(summary.atacDEG,stringsAsFactors=F)
rownames(summary.atacDEG) = 1:nrow(summary.atacDEG)
names(summary.atacDEG)[4:6] = c("nDEG", "nDREG", "nDYN")
summary.atacDEG[,1:2] = apply(summary.atacDEG[,1:2],2,function(x){as.character(x)})
summary.atacDEG[,3:6] = apply(summary.atacDEG[,3:6],2,function(x){data.matrix(x)%>%as.numeric})
summary.atacDEG = summary.atacDEG[order(summary.atacDEG$ncomps,decreasing=T),]

# combine with enrichment analysis
summary.atacDEG = summary.atacDEG %>% mutate(enrichUP=0, enrichDN=0, enrichUN=0)
for(ii in 1:nrow(summary.atacDEG)){
  id = summary.atacDEG$id[ii]
  inc.ind = which(enrich.inc$id == id)
  dec.ind = which(enrich.dec$id == id)
  unc.ind = which(enrich.un$id == id)
  if(length(inc.ind)>0){
    tcomps = paste(enrich.inc$tcomp[inc.ind], collapse=", ")
    summary.atacDEG$enrichUP[ii] = tcomps 
  }
  if(length(dec.ind)>0){
    tcomps = paste(enrich.dec$tcomp[dec.ind], collapse=", ")
    summary.atacDEG$enrichDN[ii] = tcomps 
  }
  if(length(unc.ind)>0){
    tcomps = paste(enrich.un$tcomp[unc.ind], collapse=", ")
    summary.atacDEG$enrichUN[ii] = tcomps 
  }
}

fname = "summary.atacDREG.txt"
write.table(summary.atacDEG,fname,row.names=F,col.names=T,sep="\t")


############################################################
# plot selected pwms
############################################################

# specify communities of interest
plt.comms = c(77,79)
plt.comms = paste0("community_",plt.comms)

# directory for all pwms
pwm.dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
                 "adipogenesis/July2018/atac_time_deg/communities/TF/TF735")

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
pdf("pwm_select.pdf", onefile = TRUE, height=9, width=4.5)
marrangeGrob(grobs=comm.plts, nrow=6, ncol=1, top=NULL)
dev.off()  


############################################################
# aggregate information regarding all motif clusters
############################################################

############# info from motif_communities_TF.R

# comps for each TF motif
head(motif.community.summary)

# matches with de novo motifs
head(filtered.denovo)



############################################################
# plot de novo and cluster motifs
############################################################ 

############# de novo PWMs from generate_denovoPWM_db.sh
# from=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motifdb/alldenovo
# to=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/motifEnrich
# scp -r wa3j@interactive.hpc.virginia.edu:$from $to
denovo.dir="/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/motifEnrich/alldenovo"

############# cluster average PWMs
pwm.dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
                 "adipogenesis/July2018/atac_time_deg/communities/TF/TF211")

### seqlogo plots for a given TF motif
head(motif.community5)
comm = 77 # GR
comm = 92 # EGR1

# get TF ids & PWM
tf.id1 = motif.community5$motif[motif.community5$community == comm]
tf.id2 = motif.community5$ids.orig[motif.community5$community == comm]
tfmat = get.pwm(pwm.dir=pwm.dir, motif=tf.id1)

# get de novo id
de.id = filtered.denovo$Query_ID[filtered.denovo$Target_ID == tf.id1][1]
demat = get.pwm(pwm.dir=denovo.dir, motif=de.id)

# plot the logos
pdf("EGR1.pdf")
seqLogo(tfmat)
seqLogo(demat)
dev.off()
