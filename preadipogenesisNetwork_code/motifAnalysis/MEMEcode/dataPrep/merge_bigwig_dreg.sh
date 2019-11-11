

############################################################
# convert bigwig to bedgraph
############################################################

# data directories
bigwigdir=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/bigWig_20181001/pro_preadip_bigWig/
dir=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/dREG
cd ${dir}

# path for UCSC tools
PATH=$PATH:/media/wa3j/Seagate2/Documents/software/kentUtils/bin/linux.x86_64

# strand separate primary aligned bams
# process plus and minus aligned reads separately
for i in *.bam*
do
name=$(echo $i | awk -F".bam" '{print $1}')
cat > proScript${name}.sh <<EOF
#!/bin/bash
samtools view -bh -F 20 ${name}.bam > ${name}_pro_plus.bam
samtools view -bh -f 0x10 ${name}.bam > ${name}_pro_minus.bam
EOF

# call a script for each core
echo calling proScript${name}.sh
chmod 700 proScript${name}.sh
nohup ./proScript${name}.sh &
done

# aggregate bigwigs across all preadipogenesis time points
minusfiles=$(ls ${bigwigdir}/*minus.bigWig)
plusfiles=$(ls ${bigwigdir}/*plus.bigWig)
bigWigMerge ${minusfiles} merged_pro_preadip_minus.bg
bigWigMerge ${plusfiles} merged_pro_preadip_plus.bg

# multiple minus values by -1
awk '{ print $1, $2, $3, $4*(-1) }' merged_pro_preadip_minus.bg > adipogen_minus_scaled.bg
cat merged_pro_preadip_plus.bg > adipogen_plus_scaled.bg

# generate standard bedgraphs with headers
touch minus.txt
touch plus.txt
echo "track type=bedGraph name=adipogen_minus_combined color=0,255,0 altColor=0,255,0 alwaysZero=on visibility=full" >> minus.txt
echo "track type=bedGraph name=adipogen_plus_combined color=255,0,0 altColor=255,0,0 alwaysZero=on visibility=full" >> plus.txt
cat minus.txt adipogen_minus_scaled.bg > adipogen_minus_combined.bedGraph
cat plus.txt adipogen_plus_scaled.bg > adipogen_plus_combined.bedGraph
rm minus.txt plus.txt *merged_pro* *_scaled.bg

# get chromosome sizes
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

# convert back to bigwig
mm10=mm10.chrom.sizes
bedGraphToBigWig adipogen_minus_combined.bedGraph ${mm10} adipogen_minus_combined.bigWig
bedGraphToBigWig adipogen_plus_combined.bedGraph ${mm10} adipogen_plus_combined.bigWig

# run dreg
https://dreg.dnasequence.org

# unpack results
cd results_20181222/
tar xvzf out.tar.gz
gunzip out.dREG.peak.full.bed.gz
