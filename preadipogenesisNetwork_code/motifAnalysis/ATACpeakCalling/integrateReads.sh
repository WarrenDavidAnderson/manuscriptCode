

####################################################
## integrated reads for all preadip time points
## used in empiricalSummits.R to get summits
####################################################

cd /nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/preadip_bams

# morege all preadipogenesis bams
module load samtools
samtools merge preadipMerged.bam *.bam

# isolate only concordant alignments
mv preadipMerged.bam pre.bam
samtools view -b -f 0x2 pre.bam -o preadipMerged.bam
samtools sort -n preadipMerged.bam -o sorted.bam
samtools fixmate sorted.bam fixed.bam

# produce a bed file
module load bedtools
bedtools bamtobed -i fixed.bam -bedpe > preadip0.bed
rm fixed.bam sorted.bam 

# filter and sort the bed file
awk '{OFS="\t";} {print $1,$2,$6,$7,$8,$9}' preadip0.bed > bed
sort -k 1,1 bed > sorted.bed
rm bed preadip0.bed

# from bed to bedgraph
bedtools genomecov -bg -trackline -trackopts name=preadip -i sorted.bed -g mm10.chrom.sizes > preadip.bedGraph

# generate the bigwig
module load ucsc-tools
bedGraphToBigWig preadip.bedGraph mm10.chrom.sizes preadip.bigWig




