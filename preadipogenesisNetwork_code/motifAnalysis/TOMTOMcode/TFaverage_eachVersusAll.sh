

####################################################
## overview 
####################################################

# TF match data was generated using eachVersusAllTF.sh
# TF cluster averages  were found with motif_communities_TF.R

# now match all averages against each other to reduce the factor number

####################################################
####################################################
## Run tomtom for each TF PWM against all others
####################################################
####################################################

# set directory for tomtom and go there
dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomTFvsTF/TFvTFaverage
cd ${dir}

# copy data
from=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/communities/TF/center_av_pwm
to=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomTFvsTF/TFvTFaverage
scp -r ${from} wa3j@interactive.hpc.virginia.edu:${to}


####################################################
## generate motif text files for array implementation
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomTFvsTF/TFvTFaverage
cd ${dir}/center_av_pwm

# get ids
ls *.txt | awk -F".txt" '{print $1}' > motif.ids

# read motif key into R
module load gcc/7.1.0  openmpi/3.1.4 R/3.5.3
library(dplyr)
motifs = read.table("motif.ids",header=F,stringsAsFactors=F)[,1]

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

##############################################################
## generate meme database formated file
##############################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomTFvsTF/TFvTFaverage
cd ${dir}/center_av_pwm

# create file header  
outfile=TFavDB.txt
printf '%s\n' 'MEME version 5' '' 'ALPHABET= ACGT' '' 'strands: + -'  '' > header1
printf '%s\n' 'Background letter frequencies (from uniform background):' > header2
printf '%s\n' 'A 0.25000 C 0.25000 G 0.25000 T 0.25000' '' > header3
cat header1 header2 header3 > ${outfile}
rm header1 header2 header3

# loop through motifs and generate database file
for ii in motif*
do
motif=$(echo ${ii} | awk -F".txt" '{print $1}') 
echo MOTIF ${motif} > line1
echo letter-probability matrix: > line2
echo > linen
cat line1 line2 ${ii} linen > PWM
mv ${outfile} tmp
cat tmp PWM > ${outfile}
rm line1 line2 linen PWM tmp
done

mv TFavDB.txt ..


##############################################################
## slurm script, tomtom_main.sh
##############################################################

#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH --time=16:00:00
#SBATCH --partition=standard

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomTFvsTF/TFvTFaverage
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

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomTFvsTF/TFvTFaverage
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
motiffile=${dir}/TFavDB.txt
motifdb=${dir}/TFavDB.txt
tomtom=/home/wa3j/meme-5.0.3/bin/tomtom
eval=100 # p < 0.001 for orig run, filter in R
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

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomTFvsTF/TFvTFaverage/res
cd ${dir}

# tomtom data and header
outfile=tomtomTFAVvsTFAV.txt
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
## analysis of tomtom output
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomTFvsTF/TFvTFaverage/res
cd ${dir}

# copy to local
from=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomTFvsTF/TFvTFaverage/tomtomTFAVvsTFAV.txt
to=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/communities/TF
scp -r wa3j@interactive.hpc.virginia.edu:${from} ${to}

# use r script: motif_communities_TF.R

####################################################
## Classification of Transcription Factors in Mammalia
## http://tfclass.bioinf.med.uni-goettingen.de/
####################################################



