

setwd("/media/wa3j/Seagate2/Documents/PRO/adipogenesis/package/SM/degpca")

library(dplyr)
library(ggplot2)
library(DESeq2)
library(bigWig)


#######################################################################################
## data and annotations
#######################################################################################

# TU annotation 
tu.ann = read.table("TU_20190731.bed",stringsAsFactors=F,header=F)
names(tu.ann) = c("chr","start","end","gene","xy","strand")

# gene annotation 
ge.ann = read.table("largest.interval.bed",stringsAsFactors=F,header=F)
names(ge.ann) = c("chr","start","end","gene","xy","strand")

# pro read data path
bigWig.path = paste0("/media/wa3j/Seagate2/Documents/PRO/adipogenesis",
                     "/July2018/bigWig_20181001/all_pro_sample_bigWig")

# out = list(deg=deg, deseq_obj=deseq_obj)
load("LRTres_pro_TU.RData")
deg.tu = out$deg
deseq.tu = out$deseq_obj
rm(out)

load("LRTres_pro_GE.RData")
deg.ge = out$deg
deseq.ge = out$deseq_obj
rm(out)

#######################################################################################
## functions
#######################################################################################

# gene ids
get.id = function(input=NULL){
  out = sapply(input,function(x){
    strsplit(x,"_")[[1]][2]
  })
  return(out)
} # get.id


# matching genes, indices, deg res
match.inds = function(match.ids=NULL, all.ids=NULL){
  inds = sapply(match.ids,function(x){
    which(all.ids == x)
  })
  return(inds)
} # match.inds


# organize bed data identically
bed.org = function(bed1=NULL, bed2=NULL){
  ids.tu = bed1$gene 
  ids.ge = bed2$gene 
  match.ids = ids.tu[ids.tu %in% ids.ge]
  ind.tu = match.inds(match.ids=match.ids, all.ids=ids.tu)
  ind.ge = match.inds(match.ids=match.ids, all.ids=ids.ge)
  res.tu.match = res.tu[ind.tu,]
  res.ge.match = res.ge[ind.ge,]
  if(all(ids.tu[ind.tu] == ids.ge[ind.ge])==TRUE){
    out1 = bed1[ind.tu,]
    out2 = bed2[ind.ge,]
    out = list(bed1=out1, bed2=out2)
  } else(out=c(0))
  return(out)
} # bed.org


# function to read in data and average over replicates
# set step = 1 to sum reads over entire intervals
cnt.reads = function(bw.plus=NULL, bw.minus=NULL, step=NULL,
                     bed.plus=NULL, bed.minus, sf=NULL,
                     bigWig.path=NULL, bed.both=NULL){
  dat = list()
  for(ii in 1:length(bw.plus)){
    namen = paste( strsplit(bw.plus[ii],"_")[[1]][2:3], collapse="_")
    bwplus = load.bigWig(paste0(bigWig.path,"/",bw.plus[ii]))
    bwminus = load.bigWig(paste0(bigWig.path,"/",bw.minus[ii]))
    if(step>1){
      dat0 = bed6.step.bpQuery.bigWig(bwplus, bwminus, bed.plus, step, op="sum", follow.strand=TRUE)
      plus.dat = do.call(rbind, dat0)
      dat0 = bed6.step.bpQuery.bigWig(bwplus, bwminus, bed.minus, step, op="sum", follow.strand=TRUE)
      minus.dat = do.call(rbind, dat0)
      dati = rbind(plus.dat,minus.dat)
    }
    if(step==1){
      dati = bed6.region.bpQuery.bigWig(bwplus, bwminus, bed.both, op="sum")
    }
    norm = sf[grep(namen,names(sf))]
    dat[[namen]] = dati / norm
  }
  avrep = Reduce("+",dat) / length(bw.plus)
  return(avrep)
} # cnt.reads


# function to get composite data and pause indices
composites.pauseIndex <- function(sf=NULL, bigWig.path=NULL, bed=NULL,
                               time.id=NULL, halfWindow=NULL, step=NULL,
                               pause=NULL, body=NULL) {
  
  # omit genes that are too short based on the body parameter
  min.body = pause[2] - pause[1]
  bed = bed[bed$end - bed$start > body + min.body,]
  
  # relative coords
  crel = seq((-1 * halfWindow) + 0.5*step, halfWindow - 0.5*step, by=step)
  
  # bigwig files
  bw.all = list.files(bigWig.path)
  bw.plus.all = bw.all[grep("plus",bw.all)]
  bw.minus.all = bw.all[grep("minus",bw.all)]
  
  # conditions, note - hardcoded removal of the 6d data
  conds = sapply(bw.plus.all,function(x){strsplit(x,"_")[[1]][2]})
  conds = conds[-grep("t6d",conds)]
  
  ############################################
  ## get the composite profiles
  ############################################
  
  # adjust coordinates based on window
  bed.plus0 = bed %>% filter(strand=="+")
  bed.minus0 = bed %>% filter(strand=="-")
  bed.plus = bed.plus0; bed.minus = bed.minus0
  bed.plus[,2] = bed.plus0[,2] - halfWindow
  bed.plus[,3] = bed.plus0[,2] + halfWindow
  bed.minus[,2] = bed.minus0[,3] - halfWindow
  bed.minus[,3] = bed.minus0[,3] + halfWindow
  
  # get read density data, normalize to sizefactor, average across reps
  avrep = list()
  avcomp = list()
  for(ii in conds){
    bw.plus = bw.plus.all[grep(ii,bw.plus.all)]
    bw.minus = bw.minus.all[grep(ii,bw.minus.all)]
    avrep[[ii]] = cnt.reads(bw.plus=bw.plus, bw.minus=bw.minus, step=step,
                      bed.plus=bed.plus, bed.minus=bed.minus, sf=sf,
                      bigWig.path=bigWig.path)
    avcomp[[ii]] = colMeans(avrep[[ii]]) / step
  } # loop through conditions (times)
  
  
  ############################################
  ## get the pause indices
  ############################################
  
  avrep.pause = list()
  pause.density = list()
  avrep.body = list()
  pause.index = list()
  body.density = list()
  for(ii in conds){
    
    # read density for the pause region
    bed.plus0 = bed %>% filter(strand=="+")
    bed.minus0 = bed %>% filter(strand=="-")
    bed.plus = bed.plus0; bed.minus = bed.minus0;
    bed.plus[,2] = bed.plus0[,2] + pause[1]
    bed.plus[,3] = bed.plus0[,2] + pause[2]
    bed.minus[,3] = bed.minus0[,3] - pause[1]
    bed.minus[,2] = bed.minus0[,3] - pause[2]
    bed.both = rbind(bed.plus, bed.minus)
        bw.plus = bw.plus.all[grep(ii,bw.plus.all)]
    bw.minus = bw.minus.all[grep(ii,bw.minus.all)]
    avrep.pause[[ii]] = cnt.reads(bw.plus=bw.plus, bw.minus=bw.minus, step=1, 
                            bed.plus=bed.plus, bed.minus=bed.minus, sf=sf,
                            bed.both=bed.both, bigWig.path=bigWig.path)
    pause.density[[ii]] = avrep.pause[[ii]]/(pause[2] - pause[1])
    
    # read density for the body region
    bed.plus0 = bed %>% filter(strand=="+")
    bed.minus0 = bed %>% filter(strand=="-")
    bed.plus = bed.plus0; bed.minus = bed.minus0;
    bed.plus[,2] = bed.plus0[,2] + body
    bed.minus[,3] = bed.minus0[,3] - body
    bed.both = rbind(bed.plus, bed.minus)
    avrep.body[[ii]] = cnt.reads(bw.plus=bw.plus, bw.minus=bw.minus, step=1, 
                           bed.plus=bed.plus, bed.minus=bed.minus, sf=sf,
                           bed.both=bed.both, bigWig.path=bigWig.path)
    body.density[[ii]] = avrep.body[[ii]] / (bed.both[,3] - bed.both[,2])
    
    # pause index
    pause.index[[ii]] = pause.density[[ii]] / body.density[[ii]]
    pause.index[[ii]] = pause.index[[ii]][!is.na(pause.index[[ii]])]
    pause.index[[ii]] = pause.index[[ii]][!(pause.index[[ii]]==Inf)]
    pause.index[[ii]] = pause.index[[ii]][!(pause.index[[ii]]==-Inf)]
  } # loop through conditions (times)
  
  # output data
  out1 = list()
  for(ii in conds){
    out1[[ii]] = data.frame(bp=crel, comp=avcomp[[ii]])
  }
  out = list(df=out1, pause.index=pause.index)
  return(out)
  
} # composites.pauseIndex


#######################################################################################
## select bed coords for all
#######################################################################################

# genes of interest
all.genes = tu.ann$gene[tu.ann$gene %in% ge.ann$gene]

# match ids -- intersect for differential expression
ind.tu = match.inds(match.ids=all.genes, all.ids=tu.ann$gene)
ind.ge = match.inds(match.ids=all.genes, all.ids=ge.ann$gene)
all(tu.ann$gene[ind.tu] == ge.ann$gene[ind.ge])

#######################################################################################
## generate composites and pause indices
#######################################################################################

# size factors
sf.tu = sizeFactors(deseq.tu)
sf.ge = sizeFactors(deseq.ge)

# params
halfWindow = 1000
step = 10
pause = c(20,80)
body = 500

# get composite and pause index data
comp.tu = composites.pauseIndex(sf=sf.tu, bed=tu.ann[ind.tu,], bigWig.path=bigWig.path, 
                              halfWindow=halfWindow, step=step, pause=pause, body=body)
comp.ge = composites.pauseIndex(sf=sf.ge, bed=ge.ann[ind.ge,], bigWig.path=bigWig.path, 
                              halfWindow=halfWindow, step=step, pause=pause, body=body)


# plot composite data
xlim = c(-50, 600)
ylim1 = c(0, 1)
ylim2 = c(0, 0.5)
span = 0.04
lwd = 3
bpfit = seq(xlim[1], xlim[2], length.out=1000)
condits = names(comp.tu$df)

pdf("compositePlots.pdf")
par(mfrow=c(3,3))

for(ii in condits){
  
  # composites
  plot(comp.tu$df[[ii]]$bp, comp.tu$df[[ii]]$comp, type="l", lwd=2, 
       xlim=xlim, ylim=ylim1, col="white",xlab="",ylab="")
  points(comp.ge$df$bp, comp.ge$df$comp, type="l", lwd=2, col="white")
  
  mod1 = suppressWarnings(loess(comp.tu$df[[ii]]$comp~comp.tu$df[[ii]]$bp, span=span))
  fit1 = predict(mod1, bpfit)
  mod2 = suppressWarnings(loess(comp.ge$df[[ii]]$comp~comp.ge$df[[ii]]$bp, span=span))
  fit2 = predict(mod2, bpfit)
  lines(bpfit[order(bpfit)], fit1[order(bpfit)], col="darkgray", lwd=lwd)
  lines(bpfit[order(bpfit)], fit2[order(bpfit)], col="black", lwd=lwd)
  
  # check the fit quality
  # plot(comp.tu$df$bp, comp.tu$df$comp, type="l", lwd=2, xlim=xlim, ylim=ylim1, col="black")
  # lines(bpfit[order(bpfit)], fit1[order(bpfit)], col="red", lwd=lwd)
  # 
  # plot(comp.ge$df$bp, comp.ge$df$comp, type="l", lwd=2, xlim=xlim, ylim=ylim1, col="black")
  # lines(bpfit[order(bpfit)], fit2[order(bpfit)], col="red", lwd=lwd)
  
}

dev.off()


pdf("scaledComposite.pdf")
par(mfrow=c(3,3))

for(ii in condits){
  
  # composite scaled
  mod1 = suppressWarnings(loess(comp.tu$df[[ii]]$comp~comp.tu$df[[ii]]$bp, span=span))
  fit1 = predict(mod1, bpfit)
  mod2 = suppressWarnings(loess(comp.ge$df[[ii]]$comp~comp.ge$df[[ii]]$bp, span=span))
  fit2 = predict(mod2, bpfit)
  fit1 = (fit1-min(fit1)) / (max(fit1) - min(fit1))
  fit2 = (fit2-min(fit2)) / (max(fit2) - min(fit2))
  plot(bpfit[order(bpfit)], fit1[order(bpfit)], col="white", lwd=lwd,xlab="",ylab="")
  lines(bpfit[order(bpfit)], fit1[order(bpfit)], col="darkgray", lwd=lwd)
  lines(bpfit[order(bpfit)], fit2[order(bpfit)], col="black", lwd=lwd)
  
}

dev.off()



pdf("pauseIndex.pdf")
par(mfrow=c(3,3))
stat.dat = c()
for(ii in condits){
  
  # pause index
  pind.tu = data.frame(pind=comp.tu$pause.index[[ii]],type="tu")
  pind.ge = data.frame(pind=comp.ge$pause.index[[ii]],type="ge")
  pind = rbind(pind.tu, pind.ge)
  pind$pind = log2(pind$pind)
  boxplot(pind~type, data=pind, outline=FALSE)
  wt = wilcox.test(x=pind.tu$pind, y=pind.ge$pind, alternative="two.sided") # p-value < 2.2e-16
  tu.med = median(pind.tu$pind) # 5.83
  ge.med = median(pind.ge$pind) # 2.82 
  new = c(tu.med, ge.med, wt$p.value)
  stat.dat = rbind(stat.dat, new)
}

dev.off()


# look at statistics
stat.dat = as.data.frame(stat.dat,stringsAsFactors=F)
names(stat.dat) = c("tu.pi","ge.pi","pval")
rownames(stat.dat) = 1:nrow(stat.dat)
stat.dat = stat.dat %>% mutate(fdr = p.adjust(stat.dat$pval,method="BH"))
max(stat.dat$fdr)
