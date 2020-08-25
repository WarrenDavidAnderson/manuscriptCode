

#########################################################################
## download TFBS2SNP data
#########################################################################

# analysis directory
cd /media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5

# snp2tfbs data
ftp://ccg.vital-it.ch/snp2tfbs/mapped_files/

# download data
wget ftp://ccg.vital-it.ch/snp2tfbs/mapped_files/snp2tfbs_JASPAR_CORE_2014_vert.bed.gz
wget ftp://ccg.vital-it.ch/snp2tfbs/mapped_files/snp2tfbs_JASPAR_CORE_2014_vert.txt.gz
wget ftp://ccg.vital-it.ch/snp2tfbs/mapped_files/snp2tfbs_rsIDmatch.bed.gz
gunzip *.gz

#########################################################################
## download accessible chromatin data
#########################################################################

# data dir
cd /media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5/atac

# Mohlke atac data (SGBS cells, male)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110734
gunzip *.gz
cat GSE110734_SGBS_preadipocyte_representative_peaks.narrowPeak > sgbs.preadip.pks
cat GSE110734_SGBS_adipocyte_representative_peaks.narrowPeak > sgbs.adipo.pks
tar -xvf GSE110734_RAW.tar

# Encode atac data
# Homo sapiens subcutaneous adipose tissue female adult (53 years)
# https://www.encodeproject.org/experiments/ENCSR540BML/
cat ENCFF159RKV.bam > encode.adipo.bam

# gtex data
http://hgdownload.soe.ucsc.edu/hubs/gtex/gtexHub.html
https://www.biostars.org/p/394439/


#########################################################################
## peak calling and generating a bigwig for the encode atac data (hg38)
# https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html
#########################################################################

# Encode atac data
cd /scratch/wa3j/adipo

# check if paired or single end reads
module load samtools
samtools view -c -f 1 encode.adipo.bam
samtools view -H encode.adipo.bam

# call peaks using macs2
module load macs2
file=encode.adipo.bam
out=encode.adipo.peaks
species=hs
fdr=0.05
ndups=1
macs2 callpeak -t ${file} -f BAMPE -n ${out} --outdir peaks \
-g ${species} --keep-dup ${ndups} -q ${fdr} --nomodel \
--shift -100 --extsize 200 

# get blacklisted regions
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz
gunzip hg38.blacklist.bed.gz

# processing peaks
module load gcc/7.1.0 bedtools
file=encode.adipo.peaks_peaks.narrowPeak
awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$2+$10}' \
${file} > allpeaks.txt

# sort and filter, pre-adipogenesis
sort -k1,1 -k2,2n allpeaks.txt | \
bedtools merge -i stdin -c 10,11 -o collapse | \
bedtools subtract -a stdin -b hg38.blacklist.bed | \
awk ' $2 >= 0 ' > encode_peaks_sorted.bed

# convert to hg19
wget http://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/hg38ToHg19.over.chain.gz
ml anaconda/5.2.0-py3.6
source activate crossmap
in=encode_peaks_sorted.bed
out=encodepeaks.hg19.bed
python CrossMap.py bed hg38ToHg19.over.chain.gz $in $out
source deactivate

# clean up hg19 version
sort -k1,1 -k2,2n encodepeaks.hg19.bed | \
bedtools merge -i stdin -c 2,3 -o collapse | \
cut -f1-4 | awk '{print $1, $2, $3, "1"}' > sortedpeaks.hg19.bed 

# generate bigwig for visualization
module load ucsc-tools
wget https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
echo "track type=bedGraph name=ENCODE color=0,0,0 altColor=0,0,0 alwaysZero=on visibility=full" > hdr
cat hdr sortedpeaks.hg19.bed > encodepeaks.hg19.bg
bedGraphToBigWig encodepeaks.hg19.bg hg19.chrom.sizes encodepeaks.hg19.bigWig

# main peaks files
sortedpeaks.hg19.bed
encodepeaks.hg19.bigWig



# convert the read data to bedgraph for visualization
file=encode.adipo.bam
bedtools genomecov -ibam $file -bg > encode.vis.bg

# convert visualization file to hg19
source activate crossmap
in=encode.vis.bg
out=encode.vis.hg19.bg
python CrossMap.py bed hg38ToHg19.over.chain.gz $in $out
source deactivate

# clean up hg19 visualization version
sort -k1,1 -k2,2n encode.vis.hg19.bg > encode.vis.hg19.sorted.bg
bedRemoveOverlap encode.vis.hg19.sorted.bg encode.vis.hg19.removr.bg

# convert hg19 visualization file to bedgraph
cat hdr encode.vis.hg19.removr.bg > encode.vis.bedGraph

# convert the read data to bigwig for visualization
bedGraphToBigWig encode.vis.bedGraph hg19.chrom.sizes encode.vis.bigWig

# main visualization files
encode.vis.bedGraph
encode.vis.bigWig

# move files for local analysis
to=/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5/atac/encode
from=wa3j@interactive.hpc.virginia.edu:/scratch/wa3j/adipo
file=encode.vis.bedGraph
scp -r $from/$file $to
file=encode.vis.bigWig
scp -r $from/$file $to
file=sortedpeaks.hg19.bed
scp -r $from/$file $to

# change file name
cat sortedpeaks.hg19.bed > endocdePeaks.bed

#########################################################################
## processing Mohlke peak data
#########################################################################

# data directory and bedtools
cd /media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5/atac/Mohlke/MohlkeAll
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/bedtools2/bin

# merging peaks for specific data
cat *_preadipocyte*.narrowPeak | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > preadip.bed
cat *_adipocyte*.narrowPeak | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > adipocyte.bed
cat *_Adipose*.narrowPeak | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > tissue.bed

# merge adipose cells with adipose tissue
cat adipocyte.bed tissue.bed | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > adipAll.bed

# merge pre-adipocytes with adipose
cat adipAll.bed preadip.bed | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > mohlkeAll.bed

# bedgraph for visualization
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/kentUtils/bin/linux.x86_64
bigWigMerge *_preadipocyte*.bw preadip.vis.bg
bigWigMerge *_adipocyte*.bw adip.vis.bg
bigWigMerge *_Adipose*.bw adipose.vis.bg

# downsample to ~1.8G for simple upload to cyverse
shuf -n 67578345 preadip.vis.bg | sort -k1,1 -k2,2n > preadip.vis.dn.bg
shuf -n 68334978 adip.vis.bg | sort -k1,1 -k2,2n > adip.vis.dn.bg
shuf -n 68459163 adipose.vis.bg | sort -k1,1 -k2,2n > adipose.vis.dn.bg

# add header for visualization
echo "track type=bedGraph name=PREADIP color=0,0,0 altColor=0,0,0 alwaysZero=on visibility=full" > hdr
cat hdr preadip.vis.bg > preadip.bedGraph
echo "track type=bedGraph name=ADIPOCYTE color=0,0,0 altColor=0,0,0 alwaysZero=on visibility=full" > hdr
cat hdr adip.vis.bg > adipocyte.bedGraph
echo "track type=bedGraph name=ADIPOSE color=0,0,0 altColor=0,0,0 alwaysZero=on visibility=full" > hdr
cat hdr adipose.vis.bg > adipose.bedGraph

# convert to bigwig
wget https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
bedGraphToBigWig preadip.bedGraph hg19.chrom.sizes preadip.bigWig
bedGraphToBigWig adipocyte.bedGraph hg19.chrom.sizes adipocyte.bigWig
bedGraphToBigWig adipose.bedGraph hg19.chrom.sizes adipose.bigWig




#########################################################################
## integrating ENCODE and Mohlke peak data
#########################################################################

# data directory and bedtools
cd /media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5/atac
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/bedtools2/bin

# merge all data - encode and all Mohlke
cat mohlkeAll.bed | awk -v OFS='\t' '{print $1, $2, $3}' > a
cat endocdePeaks.bed | awk -v OFS='\t' '{print $1, $2, $3}' > b
cat a b | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > main_All.bed
rm a b

# merge all adipose data - encode and adipose/adipocyte Mohlke
cat adipAll.bed | awk -v OFS='\t' '{print $1, $2, $3}' > a
cat endocdePeaks.bed | awk -v OFS='\t' '{print $1, $2, $3}' > b
cat a b | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > main_adip.bed
rm a b

# intersect all data - encode and all Mohlke
cat tissue.bed | awk -v OFS='\t' '{print $1, $2, $3}' | sort -k1,1 -k2,2n > a
cat adipocyte.bed | awk -v OFS='\t' '{print $1, $2, $3}' | sort -k1,1 -k2,2n > b
cat preadip.bed | awk -v OFS='\t' '{print $1, $2, $3}' | sort -k1,1 -k2,2n > c
cat endocdePeaks.bed | awk -v OFS='\t' '{print $1, $2, $3}' | sort -k1,1 -k2,2n > d
bedtools intersect -a a -b b > i1
bedtools intersect -a c -b d > i2
bedtools intersect -a i1 -b i2 | sort -k1,1 -k2,2n | bedtools merge -i - > main_Inter.bed
rm a b c d i1 i2


#########################################################################
## intersect eQTL coords with ENCODE and Mohlke peak data
#########################################################################

# data directory and bedtools
cd /media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5/atac
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/bedtools2/bin

# eqtl data
edat=/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5/espn.ld.hg19.bed
cat ${edat} | sort -k1,1 -k2,2n > edat.bed

# intersect eqtl data with all atac regions
bedtools intersect -wb -a main_All.bed -b edat.bed > overlap_All.bed

# intersect eqtl data with adipose atac regions
bedtools intersect -wb -a main_adip.bed -b edat.bed > overlap_adip.bed

# intersect eqtl data with preadip atac regions
bedtools intersect -wb -a preadip.bed -b edat.bed > overlap_preadip.bed

# intersect eqtl data with the intersect of all atac regions
bedtools intersect -wb -a main_Inter.bed -b edat.bed > overlap_Inter.bed


#########################################################################
## isolate TF motif PWMs
#########################################################################

##### pparg motif

# dir
cd /nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motif_databases/JASPAR

# id the motif
# MOTIF MA0065.1 PPARG::RXRA
grep RXR JASPAR2018_CORE_vertebrates_redundant.meme

# isolate the pwm
tfid=PPARG::RXRA
motiffile=JASPAR2018_CORE_vertebrates_redundant.meme
linenum=$(cat ${motiffile} | grep -n -w ${tfid} | hexdump -C | \
          head -1 | awk -F":" '{print $1}' | awk -F"|" '{print $2}')
tail -n +${linenum} ${motiffile} > part1
linestart=$(cat part1 | grep -n "letter-probability matrix:" | head -1 | cut -d':' -f1)
tail -n +${linestart} part1 > part2
lineend=$(cat part2 | grep -n -e '^[[:space:]]*$' | head -1 | cut -d':' -f1) 
cat part2 | head -${lineend} > PWM

# move the data to local
from=wa3j@interactive.hpc.virginia.edu:/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motif_databases/JASPAR/PWM
to=/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5
scp -r $from $to


##### egr1 motif

# dir
cd /nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motif_databases/JASPAR

# id the motif
grep EGR1 JASPAR2018_CORE_vertebrates_redundant.meme

# isolate the pwm
tfid=MA0162.2
motiffile=JASPAR2018_CORE_vertebrates_redundant.meme
linenum=$(cat ${motiffile} | grep -n -w ${tfid} | hexdump -C | \
          head -1 | awk -F":" '{print $1}' | awk -F"|" '{print $2}')
tail -n +${linenum} ${motiffile} > part1
linestart=$(cat part1 | grep -n "letter-probability matrix:" | head -1 | cut -d':' -f1)
tail -n +${linestart} part1 > part2
lineend=$(cat part2 | grep -n -e '^[[:space:]]*$' | head -1 | cut -d':' -f1) 
cat part2 | head -${lineend} > PWM2

# move the data to local
from=wa3j@interactive.hpc.virginia.edu:/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motif_databases/JASPAR/PWM2
to=/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5
scp -r $from $to


##### nrf1 motif

# dir
cd /nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motif_databases/JASPAR

# id the motif
grep NRF1 JASPAR2018_CORE_vertebrates_redundant.meme

# isolate the pwm
tfid=NRF1
motiffile=JASPAR2018_CORE_vertebrates_redundant.meme
linenum=$(cat ${motiffile} | grep -n -w ${tfid} | hexdump -C | \
          head -1 | awk -F":" '{print $1}' | awk -F"|" '{print $2}')
tail -n +${linenum} ${motiffile} > part1
linestart=$(cat part1 | grep -n "letter-probability matrix:" | head -1 | cut -d':' -f1)
tail -n +${linestart} part1 > part2
lineend=$(cat part2 | grep -n -e '^[[:space:]]*$' | head -1 | cut -d':' -f1) 
cat part2 | head -${lineend} > PWM3

# move the data to local
from=wa3j@interactive.hpc.virginia.edu:/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motif_databases/JASPAR/PWM3
to=/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5
scp -r $from $to






