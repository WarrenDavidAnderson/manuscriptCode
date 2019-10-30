
# /nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/deseq_res

##############################################################
## slurm script, nobkg_classic_slurm.sh
##############################################################

#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH --time=16:00:00
#SBATCH --mem 128000
#SBATCH --partition=standard

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/deseq_res/nobkg_0.1_classic
dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/dreg_for_meme/nobkg_0.1_classic
cd ${dir}

# sbatch --array=1-20 nobkg_classic_slurm.sh
cnt=${SLURM_ARRAY_TASK_ID}
${dir}/nobkg_classic_slurm_memerun.sh ${cnt}


##############################################################
## analysis script, nobkg_classic_slurm_memerun.sh
##############################################################

#!/bin/bash

# chmod u+x *.sh

cnt="$1"

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/deseq_res/nobkg_0.1_classic
dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/dreg_for_meme/nobkg_0.1_classic
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
objfun=classic

# condition identifiers
fsigup=upsig_${cnt}
fsigdn=downsig_${cnt}
funs=unsig_${cnt}
cond=${cnt}

# convert bed to fasta
fastaFromBed -fi ${mm10} -bed ${fsigup}.bed -fo ${fsigup}.fasta
fastaFromBed -fi ${mm10} -bed ${fsigdn}.bed -fo ${fsigdn}.fasta
fastaFromBed -fi ${mm10} -bed ${funs}.bed -fo ${funs}.fasta

############################
# increased 
comp=up

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

############################
# decreased 
comp=dn

# generate log
exec &> log_${comp}_${cond}.txt
echo comparison ${cond}

# run meme 
meme -o ${comp}_${cond} \
-dna -revcomp -objfun ${objfun} \
-nmotifs ${nmotifs} \
-minw ${minw} -maxw ${maxw} \
-evt ${thresh} -maxsize ${maxsize} \
${fsigdn}.fasta

