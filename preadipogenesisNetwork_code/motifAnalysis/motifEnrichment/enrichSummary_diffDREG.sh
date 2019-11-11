


####################################################
## enrichment data processing
## see Motif_in_peak_diffPeak.R
## see TFfimo.sh for generation of fimo bigWigs
####################################################

dir=/scratch/wa3j/enrich_diffDREG
cd ${dir}

# move data - differential dreg peak data
from=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/motifEnrich/dreg/enrichdat_600bp_step20.RData
to=/scratch/wa3j/enrich_diffDREG
scp -r ${from} wa3j@interactive.hpc.virginia.edu:${to}


####################################################
## main analysis scripts
####################################################

###########################
# set file names for an array job
# see TFfimo.sh
from=/scratch/wa3j/TFfimo/f750000/*array*
to=/scratch/wa3j/enrich_diffDREG
cp ${from} ${to} 


####### org_main.sh #######
#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem=128G
#SBATCH --partition=standard

dir=/scratch/wa3j/enrich_diffDREG
cd ${dir}

# sbatch --array=1-211 org_main.sh
cnt=${SLURM_ARRAY_TASK_ID}
${dir}/org_run.sh ${cnt}


##########################
####### org_run.sh #######
#!/bin/bash

# chmod u+x *.sh

cnt="$1"

dir=/scratch/wa3j/enrich_diffDREG
cd ${dir}

# modules
module purge
module load gcc/7.1.0 R/3.5.1

mkdir run_${cnt}
cp -t run_${cnt} enrichdat_600bp_step20.RData getenrichdat.R array${cnt} 
cd run_${cnt}

# call the R script
Rscript --vanilla getenrichdat.R ${cnt}

#########################################
####### R script - getenrichdat.R #######
# note set the directory for the fimo bigwigs
# based on a set number of peaks (e.g., 50k or 75k)
# see TFfimo.sh for details

library(dplyr)
library(DESeq2)
library(bigWig)
library(pracma)

# input
args = commandArgs(trailingOnly=TRUE)

# load data
load("enrichdat_600bp_step20.RData")

# region size for evaluating the peak index
delbp = 60

# add coordinate ids for dreg regions to the bed.map
coord.id = paste0(bed.map$chr,":",bed.map$start,"-",bed.map$end)
bed.map = bed.map %>% mutate(id = coord.id)
chrs = paste0("chr",1:19)
bed.map = bed.map[bed.map$chr %in% chrs,]

# fimo data dir
fname = paste0("array",args[1])
id = read.table(fname,stringsAsFactors=F)[1,1]
id = gsub(".txt", ".bigWig", id)
fimo.dir = "/scratch/wa3j/TFfimo/f750000/"
fimo.dir = paste0(fimo.dir,id)

# get factor id and load fimo data
fac = id
mot = strsplit(fac,"fimo_")[[1]][2]
mot = strsplit(mot,".bigWig")[[1]][1]
id = id.map$ids.orig[id.map$ids.new == mot]
bw.fimo = load.bigWig(fimo.dir)

# revised bed window function
bed.window2 = function(summits=NULL, half.win=NULL, chr=NULL){
	out = cbind(chr, summits-half.win, summits+half.win)
	out = as.data.frame(out,stringsAsFactors=FALSE)
  	names(out)[2:3] = c("start", "end")
 	rownames(out) = 1:nrow(out)
  	out[,2:3] = apply(out[,2:3],2,function(x){data.matrix(x) %>% as.numeric})
	out$start[out$start<0] = 0
  	return(out)
}

# loop through each pairwise comparison
fimo.frame = c()
for(comp in names(res.pairs)){
    
    # differential peak analysis significant indices
    res_ii = as.data.frame(res.pairs[[comp]])
	rownamen = rownames(res_ii)
	chrs = sapply(rownames(res_ii),function(x)strsplit(x,":")[[1]][1])
	res_ii = res_ii %>% mutate(chrs=chrs)
	rownames(res_ii) = rownamen
	keep = paste0("chr",1:19)
	res_ii = res_ii[res_ii$chrs %in% keep,]
    ind.up = which(res_ii$padj<sig.thresh & res_ii$log2FoldChange>fc.thresh)
    ind.dn = which(res_ii$padj<sig.thresh & res_ii$log2FoldChange<(-1)*fc.thresh)
    ind.un = which(res_ii$padj>sig.un & abs(res_ii$log2FoldChange)<fc.un)
    if(length(ind.up)==0){next}
	if(length(ind.dn)==0){next}
	if(length(ind.un)==0){next}
    
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
    up.ind = sapply(rownames(res_ii)[ind.up],function(x){which(bed.map$id==x)}) %>% unlist()
    dn.ind = sapply(rownames(res_ii)[ind.dn],function(x){which(bed.map$id==x)}) %>% unlist()
    un.ind = sapply(rownames(res_ii)[ind.un],function(x){which(bed.map$id==x)}) %>% unlist()
    bed.up = bed.window2(summits=bed.map$center[up.ind], half.win=half.win, chr=bed.map$chr[up.ind])
    bed.dn = bed.window2(summits=bed.map$center[dn.ind], half.win=half.win, chr=bed.map$chr[dn.ind])
    bed.un = bed.window2(summits=bed.map$center[un.ind], half.win=half.win, chr=bed.map$chr[un.ind])
    
    # matching peak data with the motif and computing features of the trace
	bps = seq(-half.win, half.win, length.out=2*half.win/step)
    updat = motif.map(bigwig=bw.fimo,coords=bed.up,step=step,bps=bps)
    dndat = motif.map(bigwig=bw.fimo,coords=bed.dn,step=step,bps=bps)
    undat = motif.map(bigwig=bw.fimo,coords=bed.un,step=step,bps=bps)
    var.up = var(updat$mean)
    var.dn = var(dndat$mean)
    var.un = var(undat$mean)
    cv.up = var(updat$mean) / mean(updat$mean)
    cv.dn = var(dndat$mean) / mean(dndat$mean)
    cv.un = var(undat$mean) / mean(undat$mean)

	# compute peak index 
    indL = which(updat$bp < -(half.win-2*delbp))
    indM = which(updat$bp >= -delbp/2 & updat$bp <= delbp/2)
    indR = which(updat$bp > (half.win-2*delbp))
    pk.index.up = mean(updat$mean[indM]) / mean(c(updat$mean[indL],updat$mean[indR]))
    pk.index.dn = mean(dndat$mean[indM]) / mean(c(dndat$mean[indL],dndat$mean[indR]))
    pk.index.un = mean(undat$mean[indM]) / mean(c(undat$mean[indL],undat$mean[indR]))

	# compute peak density differences
	up.range = mean(updat$mean[indM]) - mean(c(updat$mean[indL],updat$mean[indR]))
    dn.range = mean(dndat$mean[indM]) - mean(c(dndat$mean[indL],dndat$mean[indR]))
    un.range = mean(undat$mean[indM]) - mean(c(undat$mean[indL],undat$mean[indR]))
	upmdn = mean(updat$mean[indM]) - mean(dndat$mean[indM])
    upmun = mean(updat$mean[indM]) - mean(undat$mean[indM])
    pkdiff.up = max(upmdn, upmun)
	dnmup = mean(dndat$mean[indM]) - mean(updat$mean[indM])
    dnmun = mean(dndat$mean[indM]) - mean(undat$mean[indM])
    pkdiff.dn = max(dnmup, dnmun)
	unmup = mean(undat$mean[indM]) - mean(updat$mean[indM])
    unmdn = mean(undat$mean[indM]) - mean(dndat$mean[indM])
    pkdiff.un = max(unmup, unmdn)
	den.diff.rat.up = pkdiff.up / up.range
	den.diff.rat.dn = pkdiff.dn / dn.range
	den.diff.rat.un = pkdiff.un / un.range

	# compute the baseline adjusted integrals
    upnorm = updat$mean - mean(c(updat$mean[indL],updat$mean[indR]))
    dnnorm = dndat$mean - mean(c(dndat$mean[indL],dndat$mean[indR]))
    unnorm = undat$mean - mean(c(undat$mean[indL],undat$mean[indR]))
    x = seq(0,step*length(upnorm)-step,by=step)
    integup = trapz(x,upnorm)
    integdn = trapz(x,dnnorm)
    integun = trapz(x,unnorm) 

    # summary metrics
    IU = resultn[1,2] - resultn[1,1] # inc with motif - un with motif
    DU = resultn[1,3] - resultn[1,1] # dec with motif - un with motif
    ID = resultn[1,2] - resultn[1,3] # inc with motif - dec with motif
    DI = resultn[1,3] - resultn[1,2] # dec with motif - inc with motif
    new <- c(id, mot, comp, chi$p.value, chi$statistic, IU, DU, ID, DI, var.up, var.dn, var.un, cv.up, cv.dn, cv.un, len.up, len.dn, len.un, pk.index.up, pk.index.dn, pk.index.un, den.diff.rat.up, den.diff.rat.dn, den.diff.rat.un, integup, integdn, integun, result[1,1],result[1,2], result[1,3])
    names(new) <- c("TF", "id", "tcomp","pval","xsq", "INCminusUN", "DECminusUN", "INCminusDEC", "DECminusINC", "varINC", "varDEC", "varUN",    "cvINC", "cvDEC", "cvUN", "npeakINC", "npeakDEC", "npeakUN", "pkindINC", "pkindDEC", "pkindUN", "dratUP", "dratDN", "dratUN", "integUP", "integDN", "integUN","nmotUN","nmotINC","nmotDEC")
    fimo.frame = rbind(fimo.frame, new)
    
} # pairwise comparison

# data processing
dat = as.data.frame(fimo.frame,stringsAsFactors=F)
dat[,4:30] = apply(dat[,4:30],2,function(x){data.matrix(x) %>% as.numeric})
rownames(dat) = 1:nrow(dat)
dat$pval = p.adjust(dat$pval, method="BH")

# write data out
fname = paste0("fimoSummary_",mot,".RData")
save(dat,file=fname)

system("rm enrichdat_600bp_step20.RData")


####################################################
## aggregation of the results
####################################################

dir=/scratch/wa3j/enrich_diffDREG
cd ${dir}

module purge
module load gcc/7.1.0 R/3.5.1

# R script, aggregate.R
files = list.files()
dat.file = files[grep("fimoSummary_",files)]
motif = strsplit(dat.file,"fimoSummary_")[[1]][2]
motif = strsplit(motif,".RData")[[1]][1]
load(dat.file)
fname = paste0("fimoSummary_",motif,".txt")
write.table(dat,fname,sep="\t",quote=F,col.names=T,row.names=F)

# loop through all data folders and collect R data in .txt format
folders=$(ls -d *run_*)
for fold in ${folders}
do
cp aggregate.R ${fold}
cd ${fold}
Rscript --vanilla aggregate.R
cp *.txt* ..
cd ..
done

# isolate and aggregate text files into a list
mkdir Rtxt
mv *fimoSummary_motif* Rtxt
cd Rtxt

# R code
files = list.files()
fac.chisq = list()
for(ii in files){
	motif = strsplit(ii,"fimoSummary_")[[1]][2]
	motif = strsplit(motif,".txt")[[1]][1]
	dat = read.table(ii,stringsAsFactors=F,header=T,sep="\t")
	fac.chisq[[motif]] = dat
}
save(fac.chisq, file="enrichdatProcessed750k_dreg.RData")

# move data for local analysis
to=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/motifEnrich/dreg
from=/scratch/wa3j/enrich_diffDREG/Rtxt/enrichdatProcessed750k_dreg.RData
scp -r wa3j@interactive.hpc.virginia.edu:${from} ${to}

# analysis
Motif_in_peak_diffDREG.R

