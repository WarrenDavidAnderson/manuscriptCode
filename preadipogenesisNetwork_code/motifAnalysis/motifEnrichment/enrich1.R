
# run params
run = 8
tot = 215
per = floor(tot/8)
ind1 = per * (run - 1) + 1
ind2 = ind1 + per - 1
if(run == 8){ind2 = tot}

# directory
dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/atac_time_deg/motifEnrich")
setwd(dir)

# key packages
library(dplyr)
library(DESeq2)
library(bigWig)

# load data
load("enrichdat.RData")

# go to the fimo directory
setwd(fimo.bigWig.path)

# factor files
files = list.files()[ind1:ind2]

# factor chi-sq data - list with a frame for each factor
fac.chisq = list()

# loop through each factor
ii=1
for(fac in files){
  
  # get factor id and load fimo data
  mot = strsplit(fac,"fimo_")[[1]][2]
  mot = strsplit(mot,".bigWig")[[1]][1]
  id = id.map$ids.orig[id.map$ids.new == mot]
  bw.fimo = load.bigWig(fac)
  
  # loop through each pairwise comparison
  fimo.frame = c()
  for(comp in names(res.pairs)){
    
    # differential peak analysis significant indices
    res_ii = res.pairs[[comp]]
    ind.up = which(res_ii$padj<sig.thresh & res_ii$log2FoldChange>fc.thresh)
    ind.dn = which(res_ii$padj<sig.thresh & res_ii$log2FoldChange<(-1)*fc.thresh)
    ind.un = which(res_ii$padj>sig.un & abs(res_ii$log2FoldChange)<fc.un)
    if(length(ind.up)==0){next}
    
    # get peak data in bed format and document peak counts
    up.bed = get.map.from.res( rownames(res_ii[ind.up,]) )
    dn.bed = get.map.from.res( rownames(res_ii[ind.dn,]) )
    un.bed = get.map.from.res( rownames(res_ii[ind.un,]) )
    len.up = length(ind.up)
    len.dn = length(ind.dn)
    len.un = length(ind.un)
    
    # load bigwig factor mapping data into peak coordinates
    # number of TFBS in each peak region
    up = bed.region.bpQuery.bigWig(bw.fimo, up.bed)
    dn = bed.region.bpQuery.bigWig(bw.fimo, dn.bed)
    un = bed.region.bpQuery.bigWig(bw.fimo, un.bed)
    names(up) = rownames(up.bed)
    names(dn) = rownames(dn.bed)
    names(un) = rownames(un.bed)
    
    # match atac peaks with motif coordinates
    # get numbers of peaks with >0 motifs
    # chi square test result
    result = get.xsquare.table(up=up, dn=dn, un=un)
    resultn = apply(result,2,function(x){100*x/(sum(x))})
    chi = chisq.test(result)
    
    # get read mapping coords with half.win around each summit
    up.ind = sapply(rownames(res_ii)[ind.up],function(x){which(all.peaks==x)})
    dn.ind = sapply(rownames(res_ii)[ind.dn],function(x){which(all.peaks==x)})
    un.ind = sapply(rownames(res_ii)[ind.un],function(x){which(all.peaks==x)})
    bed.up = bed.window(summits=all.summits[up.ind], half.win=half.win)
    bed.dn = bed.window(summits=all.summits[dn.ind], half.win=half.win)
    bed.un = bed.window(summits=all.summits[un.ind], half.win=half.win)
    
    # matching peak data with the motif and computing features of the trace
    updat = motif.map(bigwig=bw.fimo,coords=bed.up,step=step,bps=bps)
    dndat = motif.map(bigwig=bw.fimo,coords=bed.dn,step=step,bps=bps)
    undat = motif.map(bigwig=bw.fimo,coords=bed.un,step=step,bps=bps)
    var.up = var(updat$mean)
    var.dn = var(dndat$mean)
    var.un = var(undat$mean)
    cv.up = var(updat$mean) / mean(updat$mean)
    cv.dn = var(dndat$mean) / mean(dndat$mean)
    cv.un = var(undat$mean) / mean(undat$mean)
    
    # summary metrics
    IU = resultn[1,2] - resultn[1,1] # inc with motif - un with motif
    DU = resultn[1,3] - resultn[1,1] # dec with motif - un with motif
    ID = resultn[1,2] - resultn[1,3] # inc with motif - dec with motif
    DI = resultn[1,3] - resultn[1,2] # dec with motif - inc with motif
    new = c(id, mot, comp, chi$p.value, chi$statistic, IU, DU, ID, DI,
            var.up, var.dn, var.un, cv.up, cv.dn, cv.un, len.up, len.dn, len.un)
    names(new) = c("TF", "id", "tcomp","pval","xsq", "INCmunisUN", "DECminusUN", 
                   "INCmunisDEC", "DECminusINC", "varINC", "varDEC", "varUN",
                   "cvINC", "cvDEC", "cvUN", "nmotifINC", "nmotifDEC", "nmotifUN")
    fimo.frame = rbind(fimo.frame, new)
    
  } # pairwise comparison
  
  fac.chisq[[mot]] = fimo.frame
  
  print(ii / length(files))
  ii=ii+1
  
} # factor loop

# set to main directory and output
dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/atac_time_deg/motifEnrich")
setwd(dir)
fout = paste0("motifChi_20190328_",run,".RData")
save(fac.chisq, file=fout)
