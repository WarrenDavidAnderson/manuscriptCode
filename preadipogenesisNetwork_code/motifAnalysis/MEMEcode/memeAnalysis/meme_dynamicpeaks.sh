
##############################################################
## slurm script, dyn_slurm.sh
##############################################################

#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --partition=standard

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/meme_dynDat
cd ${dir}

# sbatch --array=1-18 dyn_slurm.sh
cnt=${SLURM_ARRAY_TASK_ID}
${dir}/dyn_slurm_memerun.sh ${cnt}

##############################################################
## analysis script, dyn_slurm_memerun.sh
##############################################################

#!/bin/bash

# chmod u+x *.sh

cnt="$1"

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/meme_dynDat
cd ${dir}

# subdirectory for specific analyses
mkdir profile_id_${cnt}
mv -t profile_id_${cnt} *profile${cnt}.bed
cd profile_id_${cnt}

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
fsig=dyn_${cnt}
funs=unsig_${cnt}
cond=${cnt}

# convert bed to fasta
fastaFromBed -fi ${mm10} -bed preadip_sigDynamics_profile${cnt}.bed -fo ${fsigup}.fasta
fastaFromBed -fi ${mm10} -bed preadip_unDynamics_profile${cnt}.bed -fo ${funs}.fasta

# generate background models (order 3)
fasta-get-markov -m ${bkg} ${fsigup}.fasta up.bkg.fasta

############################
# dynamic versus unchanged
comp=upvun_bkg_de

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
# dynamic versus unchanged
comp=upvun_nobkg_de

# meme parameters
nmotifs=50
minw=5
maxw=15
maxsize=10000000
thresh=0.1
objfun=de

# generate log
exec &> log_${comp}_${cond}.txt
echo comparison ${cond}

# run meme 
meme -o ${comp}_${cond} \
-dna -revcomp -objfun ${objfun} \
-nmotifs ${nmotifs} \
-minw ${minw} -maxw ${maxw} \
-evt ${thresh} -maxsize ${maxsize} \
-neg ${funs}.fasta \
${fsigup}.fasta

############################
# increased 
comp=upvun_bkg_classic

# meme parameters
nmotifs=50
minw=5
maxw=15
maxsize=10000000
thresh=0.1
objfun=classic

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
${fsigup}.fasta


############################
# increased 
comp=upvun_nobkg_classic

# meme parameters
nmotifs=50
minw=5
maxw=15
maxsize=10000000
thresh=0.1
objfun=classic

# generate log
exec &> log_${comp}_${cond}.txt
echo comparison ${cond}

# run meme 
meme -o ${comp}_${cond} \
-dna -revcomp -objfun ${objfun} \
-nmotifs ${nmotifs} \
-minw ${minw} -maxw ${maxw} \
-evt ${thresh} -maxsize ${maxsize} \
${fsigup}.fasta









