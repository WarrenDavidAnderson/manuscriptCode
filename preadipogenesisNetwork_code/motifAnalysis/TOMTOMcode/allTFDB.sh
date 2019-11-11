
####################################################
## download and isolate key database TF pwms



####################################################
## download databases 3/12/2018 
####################################################

cd /nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motif_databases

# download link
http://meme-suite.org/doc/download.html 
Motif Databases (updated 28 Aug 2018)
mv motif_databases.12.18\(1\).tgz memeTF.tgz
gunzip memeTF.tgz

# move to rivanna
from=/home/wa3j/Downloads/memeTF.tar
to=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motif_databases
scp -r ${from} wa3j@interactive.hpc.virginia.edu:${to}

# unpack databases
tar -xvf memeTF.tar


####################################################
## select specific databases for further analysis
####################################################

cd /nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motif_databases

# cisbp 
cd CIS-BP
cp Mus_musculus.meme Rattus_norvegicus.meme Homo_sapiens.meme ..
cd ..

# all eukaryote
cd EUKARYOTE
cp *.meme ..
cd ..

# all mouse
cd MOUSE
cp *.meme ..
cd ..

# jaspar
cd JASPAR
cp *vertebrates* ..
cd ..

# TFBS shape
cd TFBSshape
cp *.meme ..
cd ..

# move all .meme files to a new location
mkdir TFreduced
mv *.meme TFreduced

# move HOMER PWMs to the group of TFdb pwms
from=/home/wa3j/meme-5.0.3/motif_databases/WAdbs/homer/HOMER.meme
to=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motif_databases/TFreduced
cp ${from} ${to}


####################################################
## isolate PWMs and rename database TFs
####################################################

cd /nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motif_databases/TFreduced

# aggregate all meme files
cat *.meme > ALLmeme.txt

# get list of all motif ids
cat ALLmeme.txt | grep MOTIF | awk -F"MOTIF " '{print $2}' > all.TFmotif.ids.txt

# get TFid array and number of TFs
readarray tfs < all.TFmotif.ids.txt
ntf=$(wc -l all.TFmotif.ids.txt | awk -F" " '{print $1}')

# raw data file with all meme motifs
motiffile=ALLmeme.txt

# motif identifier key
echo "idnum orig" > motif.id.key.txt

# create motif database file
printf '%s\n' 'MEME version 5' '' 'ALPHABET= ACGT' '' 'strands: + -'  '' > header1
printf '%s\n' 'Background letter frequencies (from uniform background):' > header2
printf '%s\n' 'A 0.25000 C 0.25000 G 0.25000 T 0.25000' '' > header3
cat header1 header2 header3 > TFmotifs.txt
rm header1 header2 header3

# loop through all TFs and collect PWMs into a database file
# retain the motif id key
for ii in `seq 0 $(expr ${ntf} - 1)` 
do

tfid=$(echo ${tfs[ii]})
linenum=$(cat ${motiffile} | grep -n -w "${tfid}" | hexdump -C | \
head -1 | awk -F":" '{print $1}' | awk -F"|" '{print $2}')
tail -n +${linenum} ${motiffile} > part1
linestart=$(cat part1 | grep -n "letter-probability matrix:" | head -1 | cut -d':' -f1)
tail -n +${linestart} part1 > part2
lineend=$(cat part2 | grep -n -e '^[[:space:]]*$' | head -1 | cut -d':' -f1) 
cat part2 | head -${lineend} > PWM

echo "MOTIF motif${ii}" > newid
mv TFmotifs.txt tmp
cat tmp newid PWM > TFmotifs.txt
mv motif.id.key.txt tmp 
echo "motif${ii} ${tfid}" > newid
cat tmp newid > motif.id.key.txt
rm part1 part2 PWM tmp newid 

done

# remove URLs
grep -vwE "(URL)" TFmotifs.txt > tmp
rm TFmotifs.txt
mv tmp TFmotifs.txt

#######################
## next run tomtom for de novo motifs against db motifs
## see motif_communities_memeMotif_20190310.R for the list of de novo motifs
## /media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/communities/denovo2901.txt
## next analysis: tomtom_denovo_vs_db_array.sh







