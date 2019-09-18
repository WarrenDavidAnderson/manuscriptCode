
####################################################
## generate de novo motif list
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motifdb
cd ${dir}

cat denovoPWM.txt | grep MOTIF | wc -l
cat denovoPWM.txt | grep MOTIF | awk -F" " '{print $2}' > denovo2743.txt

# move motif list
from=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motifdb/denovo2743.txt
to=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomDENOVOvsTF
to=/scratch/wa3j/tomtomDENOVOvsTF
scp -r ${from} ${to}




####################################################
## generate motif text files for array implementation
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomDENOVOvsTF
dir=/scratch/wa3j/tomtomDENOVOvsTF
cd ${dir}

# R code to generate a set of motif lists
# each list will be processed with an individual job
module load gcc/7.1.0 openmpi/2.1.5 R/3.5.3

narray = 59
library(dplyr)
motifs = read.table("denovo2743.txt",stringsAsFactors=F,header=F)
nmotif = nrow(motifs)
motif.per = floor( nmotif / narray )
for(ii in 1:narray){
	ind1 = (ii-1) * motif.per + 1
	ind2 = ind1 + motif.per - 1
	ind = ind1:ind2
	out = motifs[ind,1]
	fname = paste0("motifs_",ii,".txt")
	write.table(out,fname,col.names=F,row.names=F,quote=F,sep="\t")
	if(ii==narray & ind[length(ind)]!=nmotif){
		ii = ii + 1
		ind1 = ind2 + 1
		ind2 = nmotif
		ind = ind1:ind2
		out = motifs[ind,1]
		fname = paste0("motifs_",ii,".txt")
		write.table(out,fname,col.names=F,row.names=F,quote=F,sep="\t")
	}
}


#### analysis in scratch
from=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomDENOVOvsTF
to=/scratch/wa3j
cp -r $from $to

cd /scratch/wa3j/tomtomDENOVOvsTF


##############################################################
## slurm script, tomtom_main.sh
##############################################################

#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --partition=standard

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomDENOVOvsTF
dir=/scratch/wa3j/tomtomDENOVOvsTF
cd ${dir}

# sbatch --array=1-60 tomtom_main.sh
cnt=${SLURM_ARRAY_TASK_ID}
${dir}/tomtom_run.sh ${cnt}


##############################################################
## tomtom run script, tomtom_run.sh
## Run tomtom for de novo motifs against database TFs
##############################################################

#!/bin/bash

# chmod u+x *.sh

cnt="$1"

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomDENOVOvsTF
dir=/scratch/wa3j/tomtomDENOVOvsTF
cd ${dir}

# subdirectory for specific analyses
mkdir condit_id_${cnt}
cp motifs_${cnt}.txt condit_id_${cnt} 
cd condit_id_${cnt}

# specify output file name
outfile=tomtom_${cnt}_denovo.txt

# directory with de novo motif pwms
mdir=/nv/vol192/civeleklab/warren/mglab/atac_wafd/3t3_atac1-3/motifs/meme/tomtom_denovo/motifdb/alldenovo

# specify motifs files
motifs=$(cat motifs_${cnt}.txt)

# create header file (header)
printf '%s\n' 'MEME version 5' '' 'ALPHABET= ACGT' '' 'strands: + -'  '' > header1
printf '%s\n' 'Background letter frequencies (from uniform background):' > header2
printf '%s\n' 'A 0.25000 C 0.25000 G 0.25000 T 0.25000' '' > header3
cat header1 header2 header3 > header
rm header1 header2 header3

# tomtom header
printf '%s\n' 'Query_ID Target_ID Optimal_offset p-value E-value q-value Overlap Query_consensus Target_consensus Orientation' | tr ' ' '\t' > ${outfile}

# tomtom information
tomtom=/home/wa3j/meme-5.0.3/bin/tomtom
motifdb=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motif_databases/TFreduced/TFmotifs.txt
eval=100
minoverlap=5

# loop analysis
for ii in ${motifs}
do

# get the motif file
echo MOTIF ${ii} > line1
echo letter-probability matrix: > line2
cat header line1 line2 ${mdir}/${ii}.txt > meme.txt

# implement tomtom
${tomtom} -dist pearson -no-ssc -evalue -thresh ${eval} -m ${ii} \
-min-overlap ${minoverlap} meme.txt ${motifdb}

# store data
cd tomtom_out
mv tomtom.tsv .. 
cd ..
rm -r tomtom_out
tail -n +2 tomtom.tsv | head -n -4 > new
mv ${outfile} tmp
cat tmp new > ${outfile}

# remove extras
rm tomtom.tsv new tmp
rm line1 line2 meme.txt
done
cp ${outfile} ..


####################################################
## generate summary file with all tomtom results
####################################################

# tomtom data and header
outfile=tomtomDenovoTFall.txt
printf '%s\n' 'Query_ID Target_ID Optimal_offset p-value E-value q-value Overlap Query_consensus Target_consensus Orientation' | tr ' ' '\t' > ${outfile}

# loop through all results and combine data
for ii in *_denovo.txt
do
tail -n +2 ${ii} > new
mv ${outfile} tmp
cat tmp new > ${outfile}
rm tmp new
done

####################################################
## get the top database hit for each motif
####################################################

cd /nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomDENOVOvsTF/res
cd /scratch/wa3j/tomtomDENOVOvsTF/res

module load gcc/7.1.0 openmpi/2.1.5 R/3.5.3

library(dplyr)
fname = "tomtomDenovoTFall.txt"
dat0 = read.table(fname,header=T,stringsAsFactors=F,sep="\t")

# filter out poorly matched de novo motifs
# only one motif did not match
thresh = 0.001
dat1 = dat0 %>% filter(p.value < 0.001)
denovo.all = unique(dat0$Query_ID) # 2743
denovo.flt = unique(dat1$Query_ID) # 2466
unmatched = denovo.all[!(denovo.all %in% denovo.flt)] # 277

# loop through each de novo motif and get the top TF match
filtered = c()
for(ii in 1:length(denovo.flt)){
	mot = denovo.flt[ii]
	ttdat = dat1[dat1$Query_ID == mot,] 
	ttdat = ttdat[with(ttdat, order(p.value,Target_ID)),]
	filtered = rbind(filtered, ttdat[1,])
}

# 732 factor PWMs for 2466 motifs
length(unique(filtered$Target_ID))
hist(-log10(filtered$p.value))
max(filtered$p.value) # 0.000999632

# output information for downstream analysis
fname = "TFfromDENOVO.txt"
out = unique(filtered$Target_ID)
write.table(out,fname,row.names=F,col.names=F,quote=F,sep="\t")

fname = "DENOVOunmatchedtoTF001.txt"
out = unmatched
write.table(out,fname,row.names=F,col.names=F,quote=F,sep="\t")

fname = "TFfromDENOVOfiltered.txt"
out = filtered
write.table(out,fname,row.names=F,col.names=F,quote=F,sep="\t")

# copy data to main dir
cp *.txt ..


# main data dir
cd /nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomDENOVOvsTF
cd /scratch/wa3j/tomtomDENOVOvsTF

####################################################
## next steps 
####################################################

# match 751 TFs (TFfromDENOVO.txt) against each other and cluster
# eachVersusAllTF.sh

# match 288 unmatched de novo motifs against all databases 
# newDB_all.sh
# tomtom_denovo_vs_db_array_ALL.sh









