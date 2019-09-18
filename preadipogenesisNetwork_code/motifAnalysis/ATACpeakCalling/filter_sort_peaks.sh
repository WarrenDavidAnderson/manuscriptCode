
# redo 20190706
# cd /media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/atac_all3T3_0.05_5_50
cd /media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/atac4h_0.1_5_50


# get blacklisted regions
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz
gunzip mm10.blacklist.bed.gz

# bedtools
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/bedtools2/bin/

# document summit coordinate
awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$2+$10}' 3T3_atac_peaks.narrowPeak > allpeaks.txt

# sort and filter, pre-adipogenesis
sort -k1,1 -k2,2n allpeaks.txt | \
bedtools merge -i stdin -c 10,11 -o collapse | \
bedtools subtract -a stdin -b mm10.blacklist.bed | \
awk ' $2 >= 0 ' > 3T3_atac_4hpeaks_sorted.bed


# sort and filter, adipocyte data
awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$2+$10}' 3T3_atac_peaks.narrowPeak > allpeaks.txt
sort -k1,1 -k2,2n allpeaks.txt | \
bedtools merge -i stdin -c 10,11 -o collapse | \
bedtools subtract -a stdin -b mm10.blacklist.bed | \
awk ' $2 >= 0 ' > 3T3_atac_6dpeaks_sorted.bed
