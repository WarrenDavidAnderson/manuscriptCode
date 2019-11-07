


####################################################
## enrichment data processing
## see Motif_in_peak_diffDyn.R
## see TFfimo.sh for generation of fimo bigWigs
####################################################

dir=/scratch/wa3j/enrich_diffDyn
cd ${dir}

# move data - differential dreg peak data
from=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/motifEnrich/dyn/enrichdat_600bp_step20.RData
to=/scratch/wa3j/enrich_diffDyn
scp -r ${from} wa3j@interactive.hpc.virginia.edu:${to}


####################################################
## main analysis scripts
####################################################

###########################
# set file names for an array job
# see TFfimo.sh
from=/scratch/wa3j/TFfimo/f750000/*array*
to=/scratch/wa3j/enrich_diffDyn
cp ${from} ${to} 


####### org_main.sh #######
#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem=128G
#SBATCH --partition=standard

dir=/scratch/wa3j/enrich_diffDyn
cd ${dir}

# sbatch --array=1-211 org_main.sh
cnt=${SLURM_ARRAY_TASK_ID}
${dir}/org_run.sh ${cnt}


##########################
####### org_run.sh #######
#!/bin/bash

# chmod u+x *.sh

cnt="$1"

dir=/scratch/wa3j/enrich_diffDyn
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

# revised function to get the contingency table
get.xsquare.table2 = function(dyn=NULL, unc=NULL){
  dyn.with = length(which(dyn > 0))
  dyn.without = length(which(dyn == 0))
  nodif.with = length(which(unc > 0))
  nodif.without = length(which(unc == 0))
  result = cbind(c(nodif.with, nodif.without),
                 c(dyn.with, dyn.without))
  colnames(result) <- c("Unchanged", "Dynamic")
  rownames(result) <- c("with motif", "without motif")
  return(result)
}


# loop through each dynamic cluster
fimo.frame = c()
for(comp in unique(stem.res$Profile)){
    
	# get bed data for dynamic cluster peaks
	ind = which(stem.res$Profile == comp) 
	pkdat = stem.res$gene_symbol[ind]
	dyn.bed = get.map.from.res( pkdat )
	len.dyn = length(ind)

	# get bed data for non-dynamic peaks
	pk.un = out.uns$peakID[sample.int( nrow(out.uns), min(nrow(out.uns), length(pkdat)) )]
	unc.bed = get.map.from.res( pk.un )
	len.unc = length(pk.un)

    # load bigwig factor mapping data into peak coordinates
    # number of TFBS in each peak region
    dyn = bed.region.bpQuery.bigWig(bw.fimo, dyn.bed)
	unc = bed.region.bpQuery.bigWig(bw.fimo, unc.bed)
    names(dyn) = rownames(dyn.bed)
	names(unc) = rownames(unc.bed)

    # match atac peaks with motif coordinates
    # get numbers of peaks with >0 motifs
    # chi square test result
    result = get.xsquare.table2(dyn=dyn, unc=unc)
    resultn = apply(result,2,function(x){100*x/(sum(x))})
    chi = chisq.test(result)
    
    # get read mapping coords with half.win around each summit
	dyn.ind = sapply(pkdat,function(x){which(all.peaks==x)}) %>% unlist()
    unc.ind = sapply(pk.un,function(x){which(all.peaks==x)}) %>% unlist()
    bed.dyn = bed.window(summits=all.summits[dyn.ind], half.win=half.win)
    bed.unc = bed.window(summits=all.summits[unc.ind], half.win=half.win)
    
    # matching peak data with the motif and computing features of the trace
	bps = seq(-half.win, half.win, length.out=2*half.win/step)
    dyndat = motif.map(bigwig=bw.fimo,coords=bed.dyn,step=step,bps=bps)
    uncdat = motif.map(bigwig=bw.fimo,coords=bed.unc,step=step,bps=bps)
    var.dyn = var(dyndat$mean)
    var.unc = var(uncdat$mean)
    cv.dyn = var(dyndat$mean) / mean(dyndat$mean)
    cv.unc = var(uncdat$mean) / mean(uncdat$mean)

	# compute peak index 
    indL = which(dyndat$bp < -(half.win-2*delbp))
    indM = which(dyndat$bp >= -delbp/2 & dyndat$bp <= delbp/2)
    indR = which(dyndat$bp > (half.win-2*delbp))
    pk.index.dyn = mean(dyndat$mean[indM]) / mean(c(dyndat$mean[indL],dyndat$mean[indR]))
    pk.index.unc = mean(uncdat$mean[indM]) / mean(c(uncdat$mean[indL],uncdat$mean[indR]))

	# compute peak density differences
    dyn.range = mean(dyndat$mean[indM]) - mean(c(dyndat$mean[indL],dyndat$mean[indR]))
    unc.range = mean(uncdat$mean[indM]) - mean(c(uncdat$mean[indL],uncdat$mean[indR]))
    pkdiff.dyn = mean(dyndat$mean[indM]) - mean(uncdat$mean[indM])
    pkdiff.unc = mean(uncdat$mean[indM]) - mean(dyndat$mean[indM])
	den.diff.rat.dyn = pkdiff.dyn / dyn.range
	den.diff.rat.unc = pkdiff.unc / unc.range

	# compute the baseline adjusted integrals
    dynnorm = dyndat$mean - mean(c(dyndat$mean[indL],dyndat$mean[indR]))
    uncnorm = uncdat$mean - mean(c(uncdat$mean[indL],uncdat$mean[indR]))
    x = seq(0,step*length(dynnorm)-step,by=step)
    integdyn = trapz(x,dynnorm)
    integunc = trapz(x,uncnorm) 

    # summary metrics
    DU = resultn[1,2] - resultn[1,1] # dyn with motif - unc with motif
    UD = resultn[1,1] - resultn[1,2] # unc with motif - dyn with motif

    new <- c(id, mot, comp, chi$p.value, chi$statistic, DU, UD, var.dyn, var.unc, cv.dyn, cv.unc, len.dyn, len.unc, pk.index.dyn, pk.index.unc,  den.diff.rat.dyn, den.diff.rat.unc, integdyn, integunc, result[1,1], result[1,2])
    names(new) <- c("TF", "id", "Profile","pval","xsq", "DYNminusUNC", "UNCminusDYN", "varDYN", "varUNC", "cvDYN", "cvUNC", "npeakDYN", "npeakUN", "pkindDYN", "pkindUNC", "dratDYN", "dratUNC", "integDYN", "integUNC", "nmotUNC","nmotDYN")
    fimo.frame = rbind(fimo.frame, new)
    
} # pairwise comparison

# data processing
dat = as.data.frame(fimo.frame,stringsAsFactors=F)
dat[,4:21] = apply(dat[,4:21],2,function(x){data.matrix(x) %>% as.numeric})
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
save(fac.chisq, file="enrichdatProcessed750k_dyn.RData")

# move data for local analysis
to=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/motifEnrich/dyn
from=/scratch/wa3j/enrich_diffDyn/Rtxt/enrichdatProcessed750k_dyn.RData
scp -r wa3j@interactive.hpc.virginia.edu:${from} ${to}

# analysis
Motif_in_peak_diffDyn.R

