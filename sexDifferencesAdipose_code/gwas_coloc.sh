


############################################################
## overlap ld snps with gwas variants (Pulit 2019, Table S8)
############################################################

# data dir
cd /media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4

# rename file
cat WHRadjBMI_dimorphic_snps_merged_bfp_association_statistics.txt > tables8.txt


############################################################
## download cardiometabolic gwas 
## https://www.ebi.ac.uk/gwas/downloads 
############################################################

# data dir
cd /media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4

# download all association data
# column annotation: https://www.ebi.ac.uk/gwas/docs/fileheaders
wget https://www.ebi.ac.uk/gwas/api/search/downloads/full
cat full > gwascat.txt
cat gwascat.txt | awk 'BEGIN {OFS="\t"}; {print $1,$2,$3,$4,$5,$6,$7,$8,$12,$13,$15,$16,$17,$18}' > gwascat.reduced.txt
cat gwascat.txt | cut -f6,8,21 | awk 'BEGIN {FS="\t"}; {print $1,$2,$3}' | head > gwascat.reduced.txt

cat gwascat.txt | cut -f6,8,21 -d$'\t' > gwascat.reduced.txt

cat gwascat.txt | cut -f6,8,21,28 -d$'\t' > gwascat.reduced2.txt

cat gwascat.txt | cut -f6,8,21,27,28 -d$'\t' > gwascat.reduced3.txt
 
