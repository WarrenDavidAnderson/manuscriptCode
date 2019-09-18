

####################################################
## overview 
####################################################

# after running all de novo motifs versus select TF dbs 
# tomtom_denovo_vs_db_array.sh
# and generating a TF motif db (allTFDB.sh),
# this code takes all TF identifiers and runs tomtom for all of the
# corresponding PWMs - each versus all - to identify central factor classes

# single file with all pwms - TFmotifs.txt (allTFDB.sh)
# key with TF motif identifiers - motif.id.key.txt


####################################################
####################################################
## Run tomtom for each TF PWM against all others
####################################################
####################################################

# set directory for tomtom and go there
dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomTFvsTF
dir=/scratch/wa3j/tomtomTFvsTF
cd ${dir}

# copy data
from=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomDENOVOvsTF
from=/scratch/wa3j/tomtomDENOVOvsTF
cp ${from}/TFfromDENOVO.txt ${dir}


####################################################
## generate motif text files for array implementation
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomTFvsTF
dir=/scratch/wa3j/tomtomTFvsTF
cd ${dir}

module load gcc/7.1.0 openmpi/2.1.5 R/3.5.3

# read motif key into R
library(dplyr)
motifs = read.table("TFfromDENOVO.txt",header=F,stringsAsFactors=F)[,1]

# R code to generate a set of motif lists
# each list will be processed with an individual job
narray = 49
nmotif = length(motifs)
motif.per = floor( nmotif / narray )
for(ii in 1:narray){
	ind1 = (ii-1) * motif.per + 1
	ind2 = ind1 + motif.per - 1
	ind = ind1:ind2
	out = motifs[ind]
	fname = paste0("motifs_",ii,".txt")
	write.table(out,fname,col.names=F,row.names=F,quote=F,sep="\t")
	if(ii == narray){
		ind = (ind2+1):nmotif
		ii = ii + 1
		out = motifs[ind]
		fname = paste0("motifs_",ii,".txt")
		write.table(out,fname,col.names=F,row.names=F,quote=F,sep="\t")
	}
}


##############################################################
## slurm script, tomtom_main.sh
##############################################################

#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --partition=standard

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomTFvsTF
dir=/scratch/wa3j/tomtomTFvsTF
cd ${dir}

# sbatch --array=1-50 tomtom_main.sh
cnt=${SLURM_ARRAY_TASK_ID}
${dir}/tomtom_run.sh ${cnt}


##############################################################
## tomtom run script, tomtom_run.sh
## Run tomtom for de novo motifs against database TFs
##############################################################

#!/bin/bash

# chmod u+x *.sh

cnt="$1"

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomTFvsTF
dir=/scratch/wa3j/tomtomTFvsTF
cd ${dir}

# subdirectory for specific analyses
mkdir condit_id_${cnt}
cp motifs_${cnt}.txt condit_id_${cnt} 
cd condit_id_${cnt}

# isolate all motif identifiers
motifs=$(cat motifs_${cnt}.txt)

# output
outfile=tomtom${cnt}.txt

# create header file (header)
printf '%s\n' 'MEME version 5' '' 'ALPHABET= ACGT' '' 'strands: + -'  '' > header1
printf '%s\n' 'Background letter frequencies (from uniform background):' > header2
printf '%s\n' 'A 0.25000 C 0.25000 G 0.25000 T 0.25000' '' > header3
cat header1 header2 header3 > header
rm header1 header2 header3

# tomtom header
printf '%s\n' 'Query_ID Target_ID Optimal_offset p-value E-value q-value Overlap Query_consensus Target_consensus Orientation' | tr ' ' '\t' > ${outfile}

# tomtom information
motiffile=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motif_databases/TFreduced/TFmotifs.txt
motifdb=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motif_databases/TFreduced/TFmotifs.txt
tomtom=/home/wa3j/meme-5.0.3/bin/tomtom
eval=100 
minoverlap=5

# loop through motifs and implement tomtom
for ii in ${motifs}
do

# get meme format for the motif
linenum=$(cat ${motiffile} | grep -n -w ${ii})
linenumID=$(echo ${linenum} | awk -F":" '{print $1}')
motif=$(echo ${linenum} | awk -F":" '{print $2}' | awk -F"MOTIF " '{print $2}')
motif=$(echo ${motif} | awk -F" " '{print $1}')
linenumWM=$(expr ${linenumID} + 2)
tail -n +${linenumWM} ${motiffile} > part
lineend=$(grep -E --line-number --with-filename '^$' part | head -1 | awk -F":" '{print $2}')
lineend=$(expr ${lineend} + 1)
tail -n +${linenumID} ${motiffile} | head -${lineend} > PWM
cat header PWM > meme.txt
rm part PWM

# implement tomtom
${tomtom} -dist pearson -no-ssc -evalue -thresh ${eval} -m ${motif} \
-min-overlap ${minoverlap} meme.txt ${motifdb}
rm meme.txt

# aggregate tomtom outputs
cd tomtom_out
mv tomtom.tsv .. 
cd ..
rm -r tomtom_out
tail -n +2 tomtom.tsv | head -n -4 > new
mv ${outfile} tmp
cat tmp new > ${outfile}
rm tomtom.tsv new tmp
cp ${outfile} ..

done

####################################################
## generate summary file with all tomtom results
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomTFvsTF/res
dir=/scratch/wa3j/tomtomTFvsTF/res
cd ${dir}

# tomtom data and header
outfile=tomtomTFvsTF.txt
printf '%s\n' 'Query_ID Target_ID Optimal_offset p-value E-value q-value Overlap Query_consensus Target_consensus Orientation' | tr ' ' '\t' > ${outfile}

# loop through all results and combine data
for ii in *tomtom*
do
tail -n +2 ${ii} > new
mv ${outfile} tmp
cat tmp new > ${outfile}
rm tmp new
done


####################################################
## aggregate tomtom outputs (motif_communities_TF_20190314.R)
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomTFvsTF/res
dir=/scratch/wa3j/tomtomTFvsTF/res
cd ${dir}

# load data
library(dplyr)
fname = "tomtomTFvsTF.txt"
dat0 = read.table(fname,header=T,stringsAsFactors=F,sep="\t")

# isolate data for the factors that matched de novo motifs
facs = unique(dat0$Query_ID)
length(facs) # 732
dat0 = dat0[dat0$Target_ID %in% facs,]
length(unique(dat0$Target_ID)) # 732

# filter out poorly matched motifs
# only one motif did not match
thresh = 0.001
dat1 = dat0 %>% filter(p.value < 0.001)

fname = "tomtomTFvsTF001.txt"
write.table(dat1,fname,sep="\t",quote=F,col.names=T,row.names=F)

####################################################
## isolate PWMs for TF motifs
####################################################

cd /nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomTFvsTF
cd /scratch/wa3j/tomtomTFvsTF

# set motifs for analysis
motifs=$(cat TFfromDENOVO.txt)
motiffile=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motif_databases/TFreduced/TFmotifs.txt

# loop through each motif and create a PWM file
for ii in ${motifs}
do
linenum=$(cat ${motiffile} | grep -n -w ${ii} | awk -F":" '{print $1}')
linenumWM=$(expr ${linenum} + 2)
tail -n +${linenumWM} ${motiffile} > part
lineend=$(grep -n MOTIF part | head -1 | awk -F":" '{print $1}')
if [ -z "$lineend" ]
then
    lineend=$(wc -l part | awk -F" " '{print $1}')
	lineend=$(expr ${lineend} - 1)
else
   	lineend=$(expr ${lineend} - 2)
fi
head -${lineend} part > PWM_${ii}.txt
tr ' ' '\t' < PWM_${ii}.txt 1<> PWM_${ii}.txt
rm part
done

# move data
mkdir TF732
mv *PWM_* TF732

# revise file names 
for ii in *.txt
do
fnew=$(echo ${ii} | awk -F"PWM_" '{print $2}')
mv ${ii} ${fnew}
done


####################################################
## analysis of tomtom output
####################################################

/scratch/wa3j/tomtomTFvsTF/TF732

# copy to local
from=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo2/tomtomTFvsTF/res/tomtomTFvsTF.txt
from=/scratch/wa3j/tomtomTFvsTF/res/tomtomTFvsTF.txt
to=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/communities/TF
scp -r wa3j@interactive.hpc.virginia.edu:${from} ${to}

from=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo2/tomtomTFvsTF/TF679
from=/scratch/wa3j/tomtomTFvsTF/TF732
to=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/communities/TF
scp -r wa3j@interactive.hpc.virginia.edu:${from} ${to}

from=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo2/TFdb/TFreduced/motif.id.key.txt
from=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motif_databases/TFreduced/motif.id.key.txt
to=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/communities/TF
scp -r wa3j@interactive.hpc.virginia.edu:${from} ${to}

# use r script: motif_communities_TF.R

####################################################
## Classification of Transcription Factors in Mammalia
## http://tfclass.bioinf.med.uni-goettingen.de/
####################################################



