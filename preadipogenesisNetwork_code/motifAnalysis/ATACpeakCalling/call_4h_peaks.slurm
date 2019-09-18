
dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/preadip_bams
cd ${dir}



#############################
## generate params in R
#############################

module load gcc/7.1.0
module load R/3.5.1
library(dplyr)
fdr = c(0.1, 0.05, 0.01)
m1 = c(1,3,5,10,20)
m2 = c(10,50,100,200)
vars = expand.grid(fdr,m1,m2)
names(vars) = c("fdr","m1","m2")
vars = vars %>% filter(m2 > m1)
for(ii in 1:nrow(vars)){
	fname = paste0("vars_",ii,".txt")
	out = vars[ii,]
	write.table(out,fname,col.names=F,row.names=F,quote=F,sep="\t")
}


#############################
## main array script, main_macs2.sh
#############################


#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --partition=standard

# load macs2
module load macs2

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/preadip_bams
cd ${dir}

# sbatch --array=1-54 main_macs2.sh
cnt=${SLURM_ARRAY_TASK_ID}
${dir}/macs2_run.sh ${cnt}


##############################################################
## macs2 run script, macs2_run.sh
##############################################################

#!/bin/bash

# chmod u+x *.sh

cnt="$1"

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/preadip_bams
cd ${dir}

# subdirectory for specific analyses
mkdir params${cnt}
cp 3T3*.bam* params${cnt}
cp vars_${cnt}.txt params${cnt}
cd params${cnt}

# isolate all motif identifiers
params=$(cat vars_${cnt}.txt)

# load macs2
module load macs2

# params
files=$(ls 3T3*.bam*)
name=3T3_atac
species=mm
fdr=$(echo ${params} | awk -F" " '{print $1}')
m1=$(echo ${params} | awk -F" " '{print $2}')
m2=$(echo ${params} | awk -F" " '{print $3}')
ndups=50
output=atac4h_${fdr}_${m1}_${m2}

# initialize log file
exec &> log_${fdr}_${m1}_${m2}_3t3atac.txt
echo ${name}
echo fdr ${fdr}
echo m1 ${m1}
echo m2 ${m2}
echo ""

# run for all 3T3 files
macs2 callpeak -t ${files} -f BAMPE -n ${name} --outdir ${output} -g ${species} -B --call-summits --keep-dup ${ndups} -q ${fdr} -m ${m1} ${m2}

rm *.bam

##############################################################
## process the results
##############################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/preadip_bams
cd ${dir}

# get peak counts for each param set
# nohup sh getNpeaks.sh &
folds=$(ls -d *params*)
exec &> Npeaks.txt
for ii in ${folds}
do
cd ${dir}
cd ${ii}
d=$(ls -d */)
cd ${d}
echo ${d}
wc -l *.narrowPeak
echo " "
done

# folder - param mapping
module load gcc/7.1.0  openmpi/2.1.5 R/3.6.0
library(dplyr)
fdr = c(0.1, 0.05, 0.01)
m1 = c(1,3,5,10,20)
m2 = c(10,50,100,200)
vars = expand.grid(fdr,m1,m2)
names(vars) = c("fdr","m1","m2")
vars = vars %>% filter(m2 > m1)
vars = vars %>% mutate(fold=paste0("params",1:nrow(vars)))
vars = vars %>% mutate(cond=paste0(vars$fdr,"_",vars$m1,"_",vars$m2))

# R code to parse the results
module load gcc/7.1.0  openmpi/2.1.5 R/3.6.0
library(dplyr)
dat0 = read.table("Npeaks.txt",header=F,stringsAsFactors=F,fill=T)
ind.param = grep("atac6d_",dat0[,1])
ind.peaks = grep("3T3_atac_peaks",dat0[,2])
out = matrix(NA,length(ind.param),4)
for(ii in 1:nrow(out)){
p = strsplit(dat0[ind.param[ii],1],"_")[[1]]
c = strsplit(dat0[ind.peaks[ii],1],"_")[[1]]
out[ii,1:2] = p[2:3]
out[ii,3] = strsplit(p[4],"/")[[1]][1]
out[ii,4] = c
}
out = apply(out,2,function(x){data.matrix(x) %>% as.numeric})
res = data.frame(fdr=out[,1],m1=out[,2],m2=out[,3],pkcnt=out[,4],stringsAsFactors=F)
min(res$pkcnt) # 76804
max(res$pkcnt) # 101743

# aggregate information, output
pkdat = merge(res, vars, by=c("fdr","m1","m2"))
fname = "peaks4h_paramvar.txt"
write.table(pkdat,fname,col.names=T,row.names=F,sep="\t",quote=F)



#### move data 20190706
from=wa3j@interactive.hpc.virginia.edu:/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/preadip_bams/params16/atac4h_0.1_5_50
to=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg
scp -r $from $to


