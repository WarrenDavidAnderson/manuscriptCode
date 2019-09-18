

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

# select just chromosome 1
zcat preadip.bedGraph.gz | grep -w "chr1" > chr1.bg
echo "track type=bedGraph name=atac_integrate_chr1 color=0,0,255 altColor=0,0,255 alwaysZero=on visibility=full" >> header
cat header chr1.bg > chr1.bedGraph
gzip chr1.bedGraph

# downsample bg for visualization
gzip preadip.bedGraph
touch header
echo "track type=bedGraph name=preadip_minus color=0,0,255 altColor=0,0,255 alwaysZero=on visibility=full" >> header
zcat preadip.bedGraph.gz | tail -n +2 > body
m=10000000
cat body | shuf > tmp
head -n $m tmp | sort -k1,1 -k2,2n > temp.bg
cat header temp.bg > preadip_downsamp.bedGraph
rm header body tmp temp.bg
gzip preadip_downsamp.bedGraph

# move data to local
scp -r wa3j@interactive.hpc.virginia.edu:/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/preadip_bams/preadip.bigWig /media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/integratedBW
scp -r wa3j@interactive.hpc.virginia.edu:/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/preadip_bams/preadip_downsamp.bedGraph.gz /media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/integratedBW

scp -r wa3j@interactive.hpc.virginia.edu:/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/preadip_bams/chr1.bedGraph.gz /media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/integratedBW


####################################################
## integrated reads for all preadip time points
## individual bigwigs
## used in fimo peak enrichment analyses
####################################################

cd /nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/preadip_bams/integrateIndiv/

# set id files for the array job (R)
cnt=1
for ii in *.bam*
do
echo ${ii} > array${cnt}
cnt=$(expr ${cnt} + 1)
done

####### org_main.sh #######
#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --partition=standard

#dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/preadip_bams/integrateIndiv
#cd ${dir}
dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/bams2_atac_3T3
cd ${dir}

# sbatch --array=1-21 org_main.sh
cnt=${SLURM_ARRAY_TASK_ID}
${dir}/org_run.sh ${cnt}

##########################
####### org_run.sh #######
#!/bin/bash

# chmod u+x *.sh

#dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/preadip_bams/integrateIndiv
#cd ${dir}

# try2
cd /nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/bams2_atac_3T3

module load bedtools
module load samtools
module load ucsc-tools

cnt="$1"

# data 
dat=$(cat array${cnt})
condit=$(echo ${dat} | awk -F"_sample_atac.bam" '{print $1}')

# isolate only concordant alignments
samtools view -b -f 0x2 ${dat} -o ${condit}.bam
samtools sort -n ${condit}.bam -o sorted_${condit}.bam
samtools fixmate sorted_${condit}.bam fixed_${condit}.bam

# produce a bed file
bedtools bamtobed -i fixed_${condit}.bam -bedpe > ${condit}.bed
rm ${condit}.bam sorted_${condit}.bam fixed_${condit}.bam

# filter and sort the bed file
awk '{OFS="\t";} {print $1,$2,$6,$7,$8,$9}' ${condit}.bed | sort -k 1,1 > sorted_${condit}.bed
rm ${condit}.bed

# from bed to bedgraph
bedtools genomecov -bg -trackline -trackopts name=preadip -i sorted_${condit}.bed -g mm10.chrom.sizes > ${condit}.bedGraph

# generate the bigwig
bedGraphToBigWig ${condit}.bedGraph mm10.chrom.sizes ${condit}.bigWig
rm sorted_${condit}.bed

# move data2
from=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/bams2_atac_3T3/*.bigWig
to=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/integratedBW/indiv2
scp -r wa3j@interactive.hpc.virginia.edu:${from} ${to}

# move data
from=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/preadip_bams/integrateIndiv/*.bigWig
to=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/integratedBW
scp -r wa3j@interactive.hpc.virginia.edu:${from} ${to}



