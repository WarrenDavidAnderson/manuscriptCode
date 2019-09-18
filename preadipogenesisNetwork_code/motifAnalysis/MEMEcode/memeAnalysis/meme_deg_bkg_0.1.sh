
# /scratch/wa3j/deseq_res

##############################################################
## slurm script, bkg_de_slurm.sh
## data for this analysis were generated based on reformat_deg.sh
## data for this analysis were organized based on atac_time_deg_meme.R
##############################################################

#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --partition=standard

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/deseq_res/bkg_0.1_de
dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/dreg_for_meme/bkg_0.1_de
dir=/scratch/wa3j/deseq_res/bkg_0.1_de
cd ${dir}

# sbatch --array=1-20 bkg_de_slurm.sh
cnt=${SLURM_ARRAY_TASK_ID}
${dir}/bkg_de_slurm_memerun.sh ${cnt}

##############################################################
## analysis script, bkg_de_slurm_memerun.sh
##############################################################

#!/bin/bash

# chmod u+x *.sh

cnt="$1"

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/deseq_res/bkg_0.1_de
dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/dreg_for_meme/bkg_0.1_de
dir=/scratch/wa3j/deseq_res/bkg_0.1_de
cd ${dir}

# subdirectory for specific analyses
mkdir condit_id_${cnt}
mv -t condit_id_${cnt} *_${cnt}.bed
cd condit_id_${cnt}

# modules, paths, & files
module load bedtools
module load ucsc-tools
PATH=$PATH:/home/wa3j/meme-5.0.3/bin
PATH=$PATH:/home/wa3j/meme-5.0.3/libexec/meme-5.0.3
mm10=/nv/vol192/civeleklab/warren/MGlab/genomes/mm10/mm10.fa

# meme parameters
nmotifs=50
minw=5
maxw=15
maxsize=10000000
thresh=0.1
bkg=3
objfun=de

# condition identifiers
fsigup=upsig_${cnt}
fsigdn=downsig_${cnt}
funs=unsig_${cnt}
cond=${cnt}

# convert bed to fasta
fastaFromBed -fi ${mm10} -bed ${fsigup}.bed -fo ${fsigup}.fasta
fastaFromBed -fi ${mm10} -bed ${fsigdn}.bed -fo ${fsigdn}.fasta
fastaFromBed -fi ${mm10} -bed ${funs}.bed -fo ${funs}.fasta

# generate background models (order 3)
fasta-get-markov -m ${bkg} ${fsigup}.fasta up.bkg.fasta
fasta-get-markov -m ${bkg} ${fsigdn}.fasta dn.bkg.fasta

############################
# increased versus unchanged
comp=upvun

# generate log
exec &> log_${comp}_${cond}.txt
echo comparison ${cond}

# run meme 
meme -o ${comp}_${cond} \
-dna -revcomp -objfun ${objfun} \
-nmotifs ${nmotifs} \
-minw ${minw} -maxw ${maxw} \
-evt ${thresh} -maxsize ${maxsize} \
-bfile up.bkg.fasta \
-neg ${funs}.fasta \
${fsigup}.fasta

############################
# increased versus decreased
comp=upvdn

# generate log
exec &> log_${comp}_${cond}.txt
echo comparison ${cond}

# run meme 
meme -o ${comp}_${cond} \
-dna -revcomp -objfun ${objfun} \
-nmotifs ${nmotifs} \
-minw ${minw} -maxw ${maxw} \
-evt ${thresh} -maxsize ${maxsize} \
-bfile up.bkg.fasta \
-neg ${fsigdn}.fasta \
${fsigup}.fasta

############################
# decreased versus unchanged
comp=dnvun

# generate log
exec &> log_${comp}_${cond}.txt
echo comparison ${cond}

# run meme 
meme -o ${comp}_${cond} \
-dna -revcomp -objfun ${objfun} \
-nmotifs ${nmotifs} \
-minw ${minw} -maxw ${maxw} \
-evt ${thresh} -maxsize ${maxsize} \
-bfile dn.bkg.fasta \
-neg ${funs}.fasta \
${fsigdn}.fasta

############################
# decreased versus increased
comp=dnvup

# generate log
exec &> log_${comp}_${cond}.txt
echo comparison ${cond}

# run meme 
meme -o ${comp}_${cond} \
-dna -revcomp -objfun ${objfun} \
-nmotifs ${nmotifs} \
-minw ${minw} -maxw ${maxw} \
-evt ${thresh} -maxsize ${maxsize} \
-bfile dn.bkg.fasta \
-neg ${fsigup}.fasta \
${fsigdn}.fasta


