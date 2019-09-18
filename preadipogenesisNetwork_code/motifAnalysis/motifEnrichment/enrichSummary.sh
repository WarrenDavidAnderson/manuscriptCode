


####################################################
## enrichment data processing
## see Motif_in_peak_20190327.R
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/summary/Renrich
cd ${dir}

# move data
from=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/motifEnrich/enrichdat_600bp_step20.RData
to=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/summary/Renrich
scp -r ${from} wa3j@interactive.hpc.virginia.edu:${to}


####################################################
## main analysis scripts
####################################################

####### org_main.sh #######
#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --partition=standard

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/summary/Renrich
cd ${dir}


# sbatch --array=1-215 org_main.sh
cnt=${SLURM_ARRAY_TASK_ID}
${dir}/org_run.sh ${cnt}


##########################
####### org_run.sh #######
#!/bin/bash

# chmod u+x *.sh

cnt="$1"

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/summary/Renrich
cd ${dir}

# modules
module purge
module load gcc R/3.5.1_test

mkdir run_${cnt}
cp -t run_${cnt} enrichdat_600bp_step20.RData getenrichdat.R array${cnt} 
cd run_${cnt}

# call the R script
Rscript --vanilla getenrichdat.R ${cnt}

#########################################
####### R script - getenrichdat.R #######

library(dplyr)
library(DESeq2)
library(bigWig)

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
fimo.dir = "/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/res/summary750k/fimo_bigWig750k/"
fimo.dir = paste0(fimo.dir,id)

# get factor id and load fimo data
fac = id
mot = strsplit(fac,"fimo_")[[1]][2]
mot = strsplit(mot,".bigWig")[[1]][1]
id = id.map$ids.orig[id.map$ids.new == mot]
bw.fimo = load.bigWig(fimo.dir)

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
    indL = which(updat$bp < -(half.win-delbp))
    indM = which(updat$bp >= -delbp/2 & updat$bp <= delbp/2)
    indR = which(updat$bp > (half.win-delbp))
    pk.index.up = mean(updat$mean[indM]) / mean(c(updat$mean[indL],updat$mean[indR]))
    pk.index.dn = mean(dndat$mean[indM]) / mean(c(dndat$mean[indL],dndat$mean[indR]))
    pk.index.un = mean(undat$mean[indM]) / mean(c(undat$mean[indL],undat$mean[indR]))
    
    # summary metrics
    IU = resultn[1,2] - resultn[1,1] # inc with motif - un with motif
    DU = resultn[1,3] - resultn[1,1] # dec with motif - un with motif
    ID = resultn[1,2] - resultn[1,3] # inc with motif - dec with motif
    DI = resultn[1,3] - resultn[1,2] # dec with motif - inc with motif
    new = c(id, mot, comp, chi$p.value, chi$statistic, IU, DU, ID, DI,
            var.up, var.dn, var.un, cv.up, cv.dn, cv.un, len.up, len.dn, len.un,
			pk.index.up, pk.index.dn, pk.index.un)
    names(new) = c("TF", "id", "tcomp","pval","xsq", "INCminusUN", "DECminusUN", 
                   "INCminusDEC", "DECminusINC", "varINC", "varDEC", "varUN",
                   "cvINC", "cvDEC", "cvUN", "npeakINC", "npeakDEC", "npeakUN",
					"pkindINC", "pkindDEC", "pkindUN")
    fimo.frame = rbind(fimo.frame, new)
    
} # pairwise comparison

# data processing
dat = as.data.frame(fimo.frame,stringsAsFactors=F)
dat[,4:18] = apply(dat[,4:18],2,function(x){data.matrix(x) %>% as.numeric})
rownames(dat) = 1:nrow(dat)
dat$pval = p.adjust(dat$pval, method="BH")

# write data out
fname = paste0("fimoSummary_",mot,".RData")
save(dat,file=fname)

####################################################
## aggregation of the results
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/summary/Renrich/fimo1M
cd ${dir}

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
module purge
module load gcc R/3.5.1_test
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
save(fac.chisq, file="enrichdatProcessed2M.RData")

# move data for local analysis
to=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/motifEnrich
from=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/summary/Renrich/enrichdatProcessed750k.RData
scp -r wa3j@interactive.hpc.virginia.edu:${from} ${to}


# analysis
Motif_in_peak_20190401.R

