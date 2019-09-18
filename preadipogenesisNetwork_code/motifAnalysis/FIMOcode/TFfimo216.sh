

# based on TF cluster center averages identified in TFcluster/motif_communities_TF_20190317.R
# pwms files: /media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/communities/TF/TF216


# analysis directory
dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo
cd ${dir}

# move motif pwms over
from=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/communities/TF/TF216
to=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo
scp -r ${from} wa3j@interactive.hpc.virginia.edu:${to}


####################################################
## generate motif lists for array analysis
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo
cd ${dir}
cd TF216

# get ids
ls *.txt | awk -F".txt" '{print $1}' > motif.ids

# read motif key into R
library(dplyr)
motifs = read.table("motif.ids",header=F,stringsAsFactors=F)[,1]

# R code to generate a set of motif lists
# each list will be processed with an individual job
narray = 216
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
cd ${dir}/TF216

pwmdir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/TF216
motifs=$(ls | awk -F".txt" '{print $1}')

# create header file (header)
printf '%s\n' 'MEME version 5' '' 'ALPHABET= ACGT' '' 'strands: + -'  '' > header1
printf '%s\n' 'Background letter frequencies (from uniform background):' > header2
printf '%s\n' 'A 0.25000 C 0.25000 G 0.25000 T 0.25000' '' > header3
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
memedir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/TF216

####################################################
## implement fimo for all TF motifs
####################################################

###### slurm script (fimo_slurm.sh)

#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH -p largemem
#SBATCH --time=96:00:00
#SBATCH --mem=128G
#SBATCH --partition=standard

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo
cd ${dir}

# sbatch --array=1-216 fimo_slurm.sh
cnt=${SLURM_ARRAY_TASK_ID}
${dir}/slurm_fimorun.sh ${cnt}


###### analysis script (slurm_fimorun.sh)

#!/bin/bash
# chmod u+x *.sh

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo
cd ${dir}

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
memedir=${dir}/TF216
for ii in ${motif}
do
cp ${memedir}/${ii}.txt ${dir}/condit_id_${cnt}
done 

# modules, paths
module load bedtools
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
#SBATCH --partition=standard

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/res
cd ${dir}

# sbatch --array=1-216 org_main.sh
cnt=${SLURM_ARRAY_TASK_ID}
${dir}/org_run.sh ${cnt}

####### org_run.sh #######
#!/bin/bash

# chmod u+x *.sh

cnt="$1"

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/res
cd ${dir}

# param
nmatch=500000

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
head -${nline} fimo.tsv > ${fname}
mv ${fname} ..
cd ${dir}


####################################################
## get the number of sites and highest pval for each motif
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/res/summary750k
cd ${dir}

# set space separated data file
echo motif nsites highp > fimo.summary

# loop through and collect information
for ii in *fimo_motif*
do
motif=$(echo ${ii} | awk -F"fimo_" '{print $2}' | awk -F".txt" '{print $1}')
nsites=$(wc -l ${ii} | awk -F" " '{print $1}')
highp=$(tail ${ii} | tail -n +10 | cut -f8)
echo ${motif} ${nsites} ${highp} > new
mv fimo.summary tmp
cat tmp new > fimo.summary
rm new tmp
done


####################################################
## isolate sites of interest based on fimo summary
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/res/summary750k
cd ${dir}

# import data into R
dat = read.table("fimo.summary",header=T,sep=" ",stringsAsFactors=F)

# plot
plot(-log10(dat$highp), log10(dat$nsites)); abline(h=6.5)
plot(dat$highp, dat$nsites)

# identify the outlier - motif6814
dat[which(dat$nsites>6500000),]

# isolate this fimo data set
mkdir skip
mv fimo_motif6814.txt skip

####################################################
## generate bw and bed files for fimo data
####################################################

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
#SBATCH --partition=standard

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/res/summary3M
cd ${dir}

# sbatch --array=1-215 org_main.sh
cnt=${SLURM_ARRAY_TASK_ID}
${dir}/org_run.sh ${cnt}

##########################
####### org_run.sh #######
#!/bin/bash

# chmod u+x *.sh

cnt="$1"

# modules
module load bedtools
module load ucsc-tools/3.7.4

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/res/summary3M
cd ${dir}

# data directory
fimodat=$(cat array${cnt})
motifid=$(echo ${fimodat} | awk -F"fimo_" '{print $2}' | awk -F".txt" '{print $1}')

# data processing
tail -n +2 ${fimodat} | cut -f3-5 | sort -k1,1 -k2,2n | awk '!seen[$0]++' > fimo0_${motifid}.bed
bedtools merge -i fimo0_${motifid}.bed | awk '{OFS="\t";} {print $1, $2, $3, 1}' > fimo_${motifid}.bed
bedGraphToBigWig fimo_${motifid}.bed mm10.chrom.sizes fimo_${motifid}.bigWig
rm fimo0_${motifid}.bed 


####################################################
## move data for local analysis
## note. 750k was used for all analyses - update above
####################################################


mkdir fimo_bigWig1p5M
mv *.bigWig* fimo_bigWig1p5M

mkdir fimo_bigWig2M
mv *.bigWig* fimo_bigWig2M

from=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/res/summary750k/fimo_bed750k
to=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/motifEnrich
scp -r wa3j@interactive.hpc.virginia.edu:${from} ${to}

# data loc
/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/res/summary1M/fimo_bigWig1M
/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/res/summary750k/fimo_bigWig750k
/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/res/summary1p5M/fimo_bigWig1p5M
/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/res/summary2M/fimo_bigWig2M

# copy data for local analysis
to=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/motifEnrich
from=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/TFfimo/res/summary1p5M/fimo_bigWig1p5M
scp -r wa3j@interactive.hpc.virginia.edu:${from} ${to}

# analysis
Motif_in_peak_20190401.R



