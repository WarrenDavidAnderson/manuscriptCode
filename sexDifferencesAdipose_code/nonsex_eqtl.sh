

cd /media/wa3j/DATAPART1/Documents/adipose_sex_reviews/eQTL


##########################################################
## parse general eQTL data, output SNP bed coords
##########################################################

library(dplyr)

# read in eqtl data
fname = "Adipose_Subcutaneous.v8.signif_variant_gene_pairs.txt"
edat = read.table(fname,header=T,stringsAsFactors=F)

# filter associations
pthresh = 1e-4
eqtl.pfilt = edat %>% filter(pval_nominal < pthresh)

# reduce data to one association per gene
one.assoc = function(dat=NULL){
  out0 = dat[order(dat$pval_nominal),]
  out = out0[!duplicated(out0$gene_id),]
  return(out)
}
eqtl.pfilt.uniq = one.assoc(dat=eqtl.pfilt)

# bed coords 
chr = sapply(eqtl.pfilt.uniq$variant_id, function(x){
	strsplit(x,"_")[[1]][1]
})
end = sapply(eqtl.pfilt.uniq$variant_id, function(x){
	strsplit(x,"_")[[1]][2] %>% as.numeric
})
start = end - 1
bed = data.frame(chr=chr, start=start, end=end, stringsAsFactors=F)

# output
fname = "bedhg38.bed"
write.table(bed,fname,col.names=F,row.names=F,sep="\t",quote=F)


##########################################################
## convert SNP coords to hg19
##########################################################

# move data to rivanna
from=/media/wa3j/DATAPART1/Documents/adipose_sex_reviews/eQTL/*hg38*
to=interactive.hpc.virginia.edu:/scratch/wa3j/paper
scp -r $from $to

# run crossmap
ml anaconda/5.2.0-py3.6
source activate crossmap
python CrossMap.py bed hg38ToHg19.over.chain.gz bedhg38.bed bedhg19.bed
source deactivate

# move data back
from=interactive.hpc.virginia.edu:/scratch/wa3j/paper/bedhg19.bed
to=/media/wa3j/DATAPART1/Documents/adipose_sex_reviews/eQTL
scp -r $from $to


##########################################################
## overlap eQTL and ATAC data
##########################################################

cd /media/wa3j/DATAPART1/Documents/adipose_sex_reviews/eQTL

# sort eQTL SNP coords, 14906 SNPs at P<1e-4
sort bedhg19.bed | uniq | sort -k1,1 -k2,2n > edat.bed

# intersect eqtl data with all atac regions
file=/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5/atac/main_All.bed
bedtools intersect -wb -a ${file} -b edat.bed > overlap_All.bed

# intersect eqtl data with adipose atac regions
file=/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5/atac/main_adip.bed
bedtools intersect -wb -a ${file} -b edat.bed > overlap_adip.bed

# intersect eqtl data with preadip atac regions
file=/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5/atac/preadip.bed
bedtools intersect -wb -a ${file} -b edat.bed > overlap_preadip.bed

# intersect eqtl data with the intersect of all atac regions
file=/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5/atac/main_Inter.bed
bedtools intersect -wb -a ${file} -b edat.bed > overlap_Inter.bed


##########################################################
## process results
##########################################################

cd /media/wa3j/DATAPART1/Documents/adipose_sex_reviews/eQTL
library(dplyr)

# load intersect results
fname = "overlap_Inter.bed"
ovr.int = read.table(fname,header=F,stringsAsFactors=F)[,-c(1:3)]
fname = "overlap_preadip.bed"
ovr.pre = read.table(fname,header=F,stringsAsFactors=F)[,-c(1:3)]
fname = "overlap_adip.bed"
ovr.adp = read.table(fname,header=F,stringsAsFactors=F)[,-c(1:3)]
fname = "overlap_All.bed"
ovr.all = read.table(fname,header=F,stringsAsFactors=F)[,-c(1:3)]


# document count of association (cases) loc
tot = 14906

# document fraction of overlapping loci for each atac peak overlap class
frac.pre = 100 * nrow(ovr.pre) / tot # 13.60526
frac.adp = 100 * nrow(ovr.adp) / tot # 17.36884
frac.all = 100 * nrow(ovr.all) / tot # 18.42211
frac.int = 100 * nrow(ovr.int) / tot # 6.413525















