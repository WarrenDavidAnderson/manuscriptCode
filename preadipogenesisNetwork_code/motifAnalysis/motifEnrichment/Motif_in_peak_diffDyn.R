
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
             "adipogenesis/July2018/atac_time_deg/motifEnrich/dyn")
setwd(dir)

############################################################
# load data from differential dreg peak analysis
# load data from de novo motif analysis
############################################################

# load annotation for all peaks (bed0)
# see atac_time_deg_meme.R (see also empiricalSummits.R)
load("bed.map20191014.RData")
bed.map = bed0
all.peaks = paste0(bed0[,1],":",bed0[,2],"-",bed0[,3])
all.summits = paste0(bed0[,1],":",bed0[,2],"-",bed0[,3],",",bed0[,4])

# load coordinates for dynamic peak clusters
# analoges to res.pairs from deg/dreg
# see atac_time_cluster.R (object stem.res)
load("stem.res.RData")

# insignificant peaks for enrichment analysis
# see atac_time_cluster.R (object out.uns)
load("out.uns.RData")

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

# fdr threshold for filtering all fimo results per factor
fdr.thresh = 0.001
maxdiff.thresh = 10

# peak count threshold for considering optimal chi-sq result
pk.cnt.thresh = 100

# motif enrichment params 
step = 20
half.win = 600
delbp = 60

# smoothing param
loess.span = 0.3

# save.image("enrichdat_600bp_step20.RData")
# see enrichSummary_diffDyn.sh for analysis 



############################################################
# perform enrichment analyses for each factor
# see enrichSummary_diffDyn.sh for the implementation of this analysis on Rivanna
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
             "adipogenesis/July2018/atac_time_deg/motifEnrich/dyn")
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
load("enrichdatProcessed750k_dyn.RData")

# functions to load and plot atac data
source("/run/user/1001/gvfs/smb-share:server=home1.virginia.edu,share=cphgdesk/users/wa3j/MS/adipogenesis/vignettes/ATAC_analysis/motifEnrich/code/motifEnrich_functions.R")


############################################################
# format, sort, and filter enrichment data
############################################################

# look at relationship between variability and motif percentage dfferences
plt.dat = c()
for(ii in 1:length(fac.chisq)){
  new = fac.chisq[[ii]]
  plt.dat = rbind(plt.dat, new)
}
pdf("var_cor_diff.pdf");par(mfrow=c(2,2))
plot(plt.dat$DYNminusUNC, plt.dat$varDYN)
plot(plt.dat$UNCminusDYN, plt.dat$varUNC)
hist(plt.dat$varDYN %>% log)
hist(plt.dat$varUNC %>% log)
dev.off()

# separate into enrichment direction classes
enrich.thresh = 7 # threshold for percentage differences
enrich.dyn = plt.dat %>% filter( DYNminusUNC>enrich.thresh )
enrich.unc = plt.dat %>% filter( UNCminusDYN>enrich.thresh )

# filter based on enrichment peak index
enrich.pk.thresh = 1.2
enrich.dyn = enrich.dyn %>% filter(pkindDYN > enrich.pk.thresh) %>% mutate(class="dyn")
enrich.unc = enrich.unc %>% filter(pkindUNC > enrich.pk.thresh) %>% mutate(class="un")

# filter based on peak density difference (fraction of range)
# e.g., max(up-dn, up-un) > pk.den.diff * (up-peak - up-flank)
pk.den.diff = 0.3
enrich.dyn = enrich.dyn %>% filter(dratDYN > pk.den.diff)
enrich.unc = enrich.unc %>% filter(dratUNC > pk.den.diff)

# filter based on the baseline shifted integral
integ.thresh = 13
enrich.dyn = enrich.dyn %>% filter(integDYN > integ.thresh)
enrich.unc = enrich.unc %>% filter(integUNC > integ.thresh)

# filter based on number of peaks with motifs
npk.thresh = 300
enrich.dyn = enrich.dyn %>% filter(nmotDYN > npk.thresh, nmotUNC > npk.thresh)
enrich.unc = enrich.unc %>% filter(nmotUNC > npk.thresh, nmotDYN > npk.thresh)

# check fdrs
max(enrich.dyn$pval) # 4.508767e-15
max(enrich.unc$pval)  # 1.077544e-14


############################################################
# aggregate and summarize enrichment data
# filter for plotting enrichment data for each factor
############################################################

# aggregate all data
processed.enrich.dat = rbind(enrich.dyn, enrich.unc)
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
  ind.dyn = grep("dyn",dat$class)
  ind.unc = grep("un",dat$class)
  ndyn = length(ind.dyn)
  nunc = length(ind.unc)
  comp.dyn = paste(dat$Profile[ind.dyn],collapse=",")
  comp.unc = paste(dat$Profile[ind.unc],collapse=",")
  new = c(TF,motif,minp,ndyn,nunc,comp.dyn,comp.unc)
  enrich.summary = rbind(enrich.summary, new)
}
enrich.summary = as.data.frame(enrich.summary,stringsAsFactors=F)
rownames(enrich.summary) = 1:nrow(enrich.summary)
names(enrich.summary) = c("factor","id","minFDR","ndyn","nunc",
                          "compsDYN","compsUNC")
write.table(enrich.summary,"enrichSummary750k_dyn.txt",col.names=T,row.names=F,quote=F,sep="\t")

# prioritize data for plotting
# select the time comparison with the largest peak index
enrich.dyn.plt = enrich.plot.select(enrich.dat=enrich.dyn, comp.mode="dyn")
enrich.unc.plt = enrich.plot.select(enrich.dat=enrich.unc, comp.mode="unc")


############################################################
# plot enrichment data
############################################################

plot.enrich.dyn(comp.mode="dyn",fname="enrichINCfiltered750k_dyn.pdf",mfrow=c(3,3),
            enrich.dat=enrich.dyn.plt[order(-enrich.dyn.plt$pkindDYN),], 
            stem.res=stem.res,
            out.uns=out.uns,
            all.peaks=all.peaks,
            all.summits=all.summits,
            type="dyn",
            half.win=half.win, step=step,
            fimo.bigWig.path=fimo.bigWig.path)

plot.enrich.dyn(comp.mode="unc",fname="enrichINCfiltered750k_unc.pdf",mfrow=c(3,3),
                enrich.dat=enrich.unc.plt[order(-enrich.unc.plt$pkindUNC),], 
                stem.res=stem.res,
                out.uns=out.uns,
                all.peaks=all.peaks,
                all.summits=all.summits,
                type="unc",
                half.win=half.win, step=step,
                fimo.bigWig.path=fimo.bigWig.path)


############################################################
# plot motifs and summarize motif clusters
############################################################

dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/atac_time_deg/motifEnrich/dyn")
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
pdf("TFs_filtered750k_dyn.pdf", onefile = TRUE, height=9, width=4.5)
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
write.table(TFspreadsheet,"TFclusters750k_dyn.txt",col.names=T,row.names=F,sep="\t")
write.table(motif.community5,"TFclusterAnnotation750k_dyn.txt",col.names=T,row.names=F,sep="\t")

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
summary.atacDEG = summary.atacDEG %>% mutate(enrichDYN=0, enrichUNC=0)
for(ii in 1:nrow(summary.atacDEG)){
  id = summary.atacDEG$id[ii]
  dyn.ind = which(enrich.dyn$id == id)
  unc.ind = which(enrich.unc$id == id)
  if(length(dyn.ind)>0){
    tcomps = paste(enrich.dyn$Profile[dyn.ind], collapse=", ")
    summary.atacDEG$enrichDYN[ii] = tcomps 
  }
  if(length(unc.ind)>0){
    tcomps = paste(enrich.unc$Profile[unc.ind], collapse=", ")
    summary.atacDEG$enrichUNC[ii] = tcomps 
  }
}

fname = "summary.atacDyn.txt"
write.table(summary.atacDEG,fname,row.names=F,col.names=T,sep="\t")

###########################################################
# plot selected pwms
############################################################

# specify communities of interest
plt.comms = c(8)
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

