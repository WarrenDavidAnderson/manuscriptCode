

# based on TF cluster center averages identified in TFcluster/motif_communities_TF.R
# pwms files: /media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/communities/TF/TF211


# analysis directory
dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo
cd ${dir}

# move motif pwms over
from=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/communities/TF/TF211
to=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo
scp -r ${from} wa3j@interactive.hpc.virginia.edu:${to}


####################################################
## generate motif lists for array analysis
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo
cd ${dir}
cd TF211

# get ids
ls *.txt | awk -F".txt" '{print $1}' > motif.ids

# read motif key into R
module load gcc/7.1.0 openmpi/3.1.4 R/3.5.3
library(dplyr)
motifs = read.table("motif.ids",header=F,stringsAsFactors=F)[,1]

# R code to generate a set of motif lists
# each list will be processed with an individual job
narray = 211
nmotif = length(motifs)
motif.per = floor( nmotif / narray )
for(ii in 1:narray){
	ind1 = (ii-1) * motif.per + 1
	ind2 = ind1 + motif.per - 1
	ind = ind1:ind2
	out = motifs[ind]
	fname = paste0("motifs_",ii,".txt")
	write.table(out,fname,col.names=F,row.names=F,quote=F,sep="\t")
	if(ii == narray & ind2 != nmotif){
		ind = (ind2+1):nmotif
		ii = ii + 1
		out = motifs[ind]
		fname = paste0("motifs_",ii,".txt")
		write.table(out,fname,col.names=F,row.names=F,quote=F,sep="\t")
	}
}

mv *motifs_* ..
rm motif.ids
cd ..

####################################################
## generate meme format pwm files
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo
cd ${dir}/TF211

pwmdir=${dir}/TF211
motifs=$(ls | awk -F".txt" '{print $1}')

# create header file (header)
printf '%s\n' 'MEME version 5' '' 'ALPHABET= ACGT' '' 'strands: + -'  '' > header1
printf '%s\n' 'Background letter frequencies (from uniform background):' > header2
printf '%s\n' 'A 0.29000 C 0.21000 G 0.21000 T 0.29000' '' > header3
cat header1 header2 header3 > header
rm header1 header2 header3

# loop through all motifs and create meme files
for ii in ${motifs}
do
cat ${pwmdir}/${ii}.txt > pwm
echo MOTIF ${ii} > line1
echo letter-probability matrix: > line2
cat header line1 line2 pwm > ${ii}.txt
rm pwm line1 line2
done

# meme format pwm file dir
memedir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/TF211

####################################################
## implement fimo for all TF motifs
####################################################

###### slurm script (fimo_slurm.sh)
cd /scratch/wa3j/TFfimo

#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH -p largemem
#SBATCH --time=96:00:00
#SBATCH --mem=128G
#SBATCH --partition=standard

# sbatch --array=1-211 fimo_slurm.sh
cnt=${SLURM_ARRAY_TASK_ID}
${dir}/slurm_fimorun.sh ${cnt}


###### analysis script (slurm_fimorun.sh)

#!/bin/bash
# chmod u+x *.sh

# threshold
thresh=0.01

# count indicator
cnt="$1"

# navigate to a new directory for each factor
mkdir condit_id_${cnt}
cp -t condit_id_${cnt} motifs_${cnt}.txt 
cd condit_id_${cnt}

# loop through each unmatched motif
motif=$(cat motifs_${cnt}.txt)

# get the appropriate pwm file
memedir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/TF211
for ii in ${motif}
do
cp ${memedir}/${ii}.txt ${dir}/condit_id_${cnt}
done 

# modules, paths
module load gcc/7.1.0 bedtools
module load ucsc-tools
fimo=/home/wa3j/meme-5.0.3/bin/fimo

# fimo params
mm10=/nv/vol192/civeleklab/warren/MGlab/genomes/mm10/mm10.fa
max=100000000
outdir=out

for ii in ${motif}
do

# fimo header - initiate results file
outfile=fimo_${ii}.txt
printf '%s\n' 'motif_id chr start end' | tr ' ' '\t' > ${outfile}

# impliment fimo 
${fimo} -thresh ${thresh} --max-stored-scores ${max} -motif ${ii} -o ${outdir} ${ii}.txt ${mm10}

# process fimo results
cd ${outdir}
cp fimo.tsv ..
cd ..
cat fimo.tsv | grep -v sequence_name | cut -f1,3-5 | \
sort -k1,1 -k2,2n | awk '!seen[$0]++' | \
tail -n +5 > fimo0.bed
mv ${outfile} tmp
cat tmp fimo0.bed > ${outfile}
rm fimo0.bed 
#rm -r ${outdir}

done


####################################################
## for each factor, get fimo data for a set number of matches
## nmatch is varied in multiple analyses in Rivanna
####################################################

####### org_main.sh #######
#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --mem=128G
#SBATCH --partition=standard

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo
dir=/scratch/wa3j/TFfimo
cd ${dir}

# sbatch --array=1-211 org_main.sh
cnt=${SLURM_ARRAY_TASK_ID}
${dir}/org_run.sh ${cnt}

####### org_run.sh #######
#!/bin/bash

# chmod u+x *.sh

cnt="$1"

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo
dir=/scratch/wa3j/TFfimo
cd ${dir}

# param
nmatch=750000

# data directory
datdir=condit_id_${cnt}

# filter the data
cd ${datdir}
motif=$(cat motifs_*)
echo ${motif}
fname=fimo_${motif}.txt
pval=$(head -${nmatch} fimo.tsv | tail | tail -n +10 | cut -f8)
pval=$(expr ${pval})
nline=$(cut -f8 fimo.tsv | LC_ALL=C fgrep -n ${pval} | tail | tail -n +10 | awk -F":" '{print $1}')
head -${nline} fimo.tsv > ${nmatch}_${pval}_${fname}
mv ${nmatch}_${pval}_${fname} ..
cd ${dir}


####################################################
## get the number of sites and highest pval for each motif
####################################################

dir=/scratch/wa3j/TFfimo/f750000
cd ${dir}

# param
nmatch=750000
outfile=fimo.summary.${nmatch}

# set space separated data file
echo motif nsites highp > ${outfile}

# loop through and collect information
for ii in *fimo_motif*
do
motif=$(echo ${ii} | awk -F"fimo_" '{print $2}' | awk -F".txt" '{print $1}')
nsites=$(wc -l ${ii} | awk -F" " '{print $1}')
highp=$(tail ${ii} | tail -n +10 | cut -f8)
echo ${motif} ${nsites} ${highp} > new
mv ${outfile} tmp
cat tmp new > ${outfile}
rm new tmp
done


####################################################
## isolate sites of interest based on fimo summary
####################################################

dir=/scratch/wa3j/TFfimo/f750000
cd ${dir}

module load gcc/7.1.0 openmpi/3.1.4 R/3.5.3 rstudio

# import data into R
fname = list.files(pattern="fimo.summary")
dat = read.table(fname,header=T,sep=" ",stringsAsFactors=F)

# plot
pdf("pval75k.pdf"); par(mfrow=c(2,2))
plot(-log10(dat$highp), log10(dat$nsites)); abline(h=6.3)
plot(dat$highp, dat$nsites)
dev.off()

# identify the outliers in pval for motif mapping
# consider these motifs for removal
dat[which(log10(dat$nsites)>6.3),]

# 75k
        motif  nsites    highp
85  motif7252 4874797 1.05e-03
97  motif8781 3790284 2.61e-05
114 motif8448 2592191 5.03e-06
121 motif8629 3219029 5.48e-09
132 motif2312 5358393 9.95e-05
175 motif8633 2018510 7.90e-08


####################################################
## generate bw and bed files for fimo data
####################################################

cd /scratch/wa3j/TFfimo/f750000

# set id files for the array job
cnt=1
for ii in *fimo_motif*
do
echo ${ii} > array${cnt}
cnt=$(expr ${cnt} + 1)
done

# get mouse chromosome annotation
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

####### org_main.sh #######
#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem=128G
#SBATCH --partition=standard

dir=/scratch/wa3j/TFfimo/f750000
cd ${dir}

# sbatch --array=1-211 org_main.sh
cnt=${SLURM_ARRAY_TASK_ID}
${dir}/org_run.sh ${cnt}

##########################
####### org_run.sh #######
#!/bin/bash

# chmod u+x *.sh

cnt="$1"

# modules
module load bedtools ucsc-tools/3.7.4

dir=/scratch/wa3j/TFfimo/f750000
cd ${dir}

# data directory
fimodat=$(cat array${cnt})
motifid=$(echo ${fimodat} | awk -F".txt" '{print $1}')

# data processing
tail -n +2 ${fimodat} | cut -f3-5 | sort -k1,1 -k2,2n | awk '!seen[$0]++' > fimo0_${motifid}.bed
bedtools merge -i fimo0_${motifid}.bed | awk '{OFS="\t";} {print $1, $2, $3, 1}' > ${motifid}.bed
bedGraphToBigWig ${motifid}.bed mm10.chrom.sizes ${motifid}.bigWig
rm fimo0_${motifid}.bed 


####################################################
## move data for local analysis
## note. 750k was used for all analyses 
####################################################


# next analysis
Motif_in_peak.R



