

dir=/media/wa3j/DATAPART1/Documents/adipose_sex_reviews/HiC
cd ${dir}


###############################################################
## Download data
###############################################################

# download and extract promoter-capture HiC
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110619
tar -xvf GSE110619_RAW.tar
gunzip GSM3004355_BSA1_Down_hg19_cisInt.txt.gz

# download and process GENCODE34 gene annotation data
# Comprehensive gene annotation, GRCh37 (hg19)
# https://www.gencodegenes.org/human/
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz
gunzip gencode.v34lift37.annotation.gtf.gz
cat gencode.v34lift37.annotation.gtf > gencode.gtf

###############################################################
## Subset the HiC data 
###############################################################

# remove annotation line from the HiC coords
tail -n +2 GSM3004355_BSA1_Down_hg19_cisInt.txt > HiC.all.bed

# separate the bait coordinates
cut -f1-3,7 HiC.all.bed > bait.bed

###############################################################
## Process gene/TSS annotation data
###############################################################

# get the first exons for protein coding genes
# this consists of the most distal TTSs
grep 'transcript_type "protein_coding"' gencode.gtf | \
awk '{if($3=="exon"){print $0}}' | \
grep -w "exon_number 1" | \
cut -f1,4,5,7,9 | tr ";" "\t" | \
awk '{for(i=5;i<=NF;i++){if($i~/^gene_name/){a=$(i+1)}} print $1,$2,$3,a,"na",$4}' | \
tr " " "\t" | tr -d '"' > gencode.firstExon.bed

########################
# isolate distal TSSs based on the strand for each first exon

library(dplyr)

# import data (R) for first exons, annotate, and remove duplicate transcripts
fname = "gencode.firstExon.bed"
dat0 = read.table(fname,header=F,stringsAsFactors=F)
names(dat0) = c('chr', 'start', 'end', 'gene', 'xy', 'strand')
dat0 = unique(dat0)
gencode.firstExon = dat0

# separate by strand and sort to identify TSSs
plus = gencode.firstExon %>% filter(strand=="+")
minus = gencode.firstExon %>% filter(strand=="-")
plusU = plus[order(plus$start,decreasing=FALSE),] 
minusU = minus[order(minus$end,decreasing=TRUE),] 
plusD = plus[order(plus$start,decreasing=TRUE),] 
minusD = minus[order(minus$end,decreasing=FALSE),]
plusU = plusU[!duplicated(plusU$gene),]
minusU = minusU[!duplicated(minusU$gene),]
plusD = plusD[!duplicated(plusD$gene),]
minusD = minusD[!duplicated(minusD$gene),]

# get TSS ranges for all first exons and combine strands
plus = merge(plusU,plusD,by=c('gene','strand'))
plus = plus[,c(3,4,8,1,6,2)]
minus = merge(minusU,minusD,by=c('gene','strand'))
minus = minus[,c(3,9,5,1,6,2)]
names(plus) = names(minus) = names(dat0)
TSS = rbind(plus,minus)
TSS$start = TSS$start - 1
TSS$end = TSS$end + 1

# output TSS data for further analysis
fname = "TSS.bed"
write.table(TSS,fname,col.names=F,row.names=F,quote=F,sep="\t")


###############################################################
## Annotate HiC data gene/TSS in bait coordinates
###############################################################

dir=/media/wa3j/DATAPART1/Documents/adipose_sex_reviews/HiC
cd ${dir}

PATH=$PATH:/media/wa3j/Seagate2/Documents/software/bedtools2/bin

# sort the bed files
sort -k1,1 -k2,2n HiC.all.bed > HiC.sorted.bed
sort -k1,1 -k2,2n TSS.bed > TSS.sorted.bed

# intersect the TSS regions with the bait regions
bedtools intersect -wo -a HiC.sorted.bed -b TSS.sorted.bed > HiC.TSS.bed


###############################################################
## Main Analysis, in R
## Integrate HiC and eQTL data - filter for baits for eGenes
## Select LD ranges for loci corresponding to eGenes in baits
## Output filtered HiC data and filter LD ranges
## Intersect capture regions with LD ranges
###############################################################

dir=/media/wa3j/DATAPART1/Documents/adipose_sex_reviews/HiC
cd ${dir}

library(dplyr)
library(iterators)
library(foreach)
library(doParallel)

###################
## Import eQTL and HiC data
###################

# peer associations with ld data and hg19 coords (CIdata.ldsnps.hg19)
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/CIdata.ldsnps.hg19.RData"
load(fname)

# use this as input to generate null probability
eqtl.dat = CIdata.ldsnps.hg19

# HiC data
fname = "HiC.TSS.bed"
hic0 = read.table(fname,header=FALSE,stringsAsFactors=FALSE,sep="\t")
hic0 = hic0[,c(1:7,11,14)]
names(hic0) = c("chrBait","startBait","endBait","chrCap","startCap","endCap","ChicScore","gene","overlap")
hic0 = unique(hic0)

###################
## Determine ranges for SNPs in LD
## with the lead SNPs for each eGene
###################

# isolate lead LD SNP ranges for each eGene
snp.data = eqtl.dat[with(eqtl.dat, order(genename,pvalue)),]
loci <- isplit(snp.data,snp.data$genename)
coords = foreach(locus=loci) %dopar% {
	gene = locus$key[[1]]
	chr = locus$value$chr[1]
	start = min(locus$value$start)[1]
	end = max(locus$value$end)[1]
	out = c(chr,start,end,gene)
	return(out)
}
coords <- data.frame(matrix(unlist(coords), nrow=length(coords), byrow=T), stringsAsFactors=F)
coords[,2:3] = apply(coords[,2:3],2,function(x){data.matrix(x)%>%as.numeric})
names(coords) = c("chr","start","end","gene")

###################
## filter HiC and LD range data
## to include baits/loci for 
## the same set of eGenes
###################

# filter HiC data and eQTL data to contain the same set of eGenes
hic.egenes = hic0[hic0$gene %in% coords$gene,]
qtl.egenes = coords[coords$gene %in% hic.egenes$gene,]

# metric checks
length(unique(coords$gene)) # 2394 total eGenes
length(unique(hic.egenes$gene)) # 688 eGenes in chromatin contacts
length(unique(qtl.egenes$gene)) # 688 eGenes in chromatin contacts

# output eQTL loci coordinates for eGenes in bait fragments
fname = "eSNP.bed"
write.table(qtl.egenes,fname,col.names=F,row.names=F,quote=F,sep="\t")

# output capture coordinates with eGenes in bait
out = hic.egenes[,c(4:6,1:3,7:9)]
fname = "eGene.bed"
write.table(out,fname,col.names=F,row.names=F,quote=F,sep="\t")

###################
## Intersect capture coordinates, 
## for baits encompassing eGene TSS ranges, 
## with eSNP LD ranges 
###################

# sort the bed files
comm1 = "sort -k1,1 -k2,2n eSNP.bed > eSNP.sorted.bed"
comm2 = "sort -k1,1 -k2,2n eGene.bed > eGene.sorted.bed"
system(comm1)
system(comm2)

# intersect the filtered capture regions with the eGene loci ranges
comm = "bedtools intersect -wo -a eGene.sorted.bed -b eSNP.sorted.bed > eGene.eSNP.bed"
system(comm)

###################
## remove intermediate files
###################
system("rm eSNP.bed eGene.bed eSNP.sorted.bed eGene.sorted.bed")



###################
## process final intersection results
###################

# input data
fname = "eGene.eSNP.bed"
dat0 = read.table(fname,header=FALSE,stringsAsFactors=FALSE,sep="\t")
dat0 = dat0[dat0$V8 == dat0$V13,]

# number of eGenes identified
egenes = unique(dat0$V13)
length(egenes) # 39



###############################################################
## Overlap HiC loci - TSS pairs with the 162 DEGs
###############################################################

# HiC data
fname = "eGene.eSNP.bed"
dat0 = read.table(fname,header=FALSE,stringsAsFactors=FALSE,sep="\t")
dat0 = dat0[dat0$V8 == dat0$V13,]

# DEG data
fname = "TableS1.txt"
deg0 = read.table(fname,header=T,stringsAsFactors=F)

hic.deg = intersect(dat0$V13, deg0$gene) # 0



###############################################################
## Estimate background probability of overlapping random loci 
## with HiC coordinates
###############################################################

dir=/media/wa3j/DATAPART1/Documents/adipose_sex_reviews/HiC
cd ${dir}

library(dplyr)
library(iterators)
library(foreach)
library(doParallel)

###################
## annotated eQTL data set
###################

# load the eqtl results based on the PEER analysis
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/eqtl_data.RData"
load(fname)
head(eqtldat)

# add gene ids
names(eqtldat)[2] = "gene_id"
edat = merge(eqtldat, ann_gene0, by="gene_id")

# save data set
save(edat, file="edat.RData")

###################
## load data and basic documentation
###################

# PEER eQTL data
load("edat.RData")

# LD data
fname="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/LD08.txt"
ld.dat = read.table(fname,header=T,stringsAsFactors=F)

# HiC data
fname = "HiC.TSS.bed"
hic0 = read.table(fname,header=FALSE,stringsAsFactors=FALSE,sep="\t")
hic0 = hic0[,c(1:7,11,14)]
names(hic0) = c("chrBait","startBait","endBait","chrCap","startCap","endCap","ChicScore","gene","overlap")
hic0 = unique(hic0)

# genes
genes = unique(edat$gene_name) %>% as.character

# number of eQTLs, number of permutations, ld range
# https://journals.physiology.org/doi/full/10.1152/physiolgenomics.00178.2002
nloci = 2408
nperm = 1000
haplo = 60000 / 2

###################
## loop to identify counts of overlapping loci
###################

# loop analysis
Nloci = rep(0,nperm)
for(ii in 1:nperm){

	# select loci
	ensid = sample(genes, nloci, replace=FALSE)
	assoc = edat[edat$gene_name %in% ensid,]
	assoc = assoc[order(assoc$pvalue),]
	assoc = assoc[!duplicated(assoc$gene_name),]

	# filter HiC data and eQTL data to contain the same set of eGenes
	hic.egenes = hic0[hic0$gene %in% assoc$gene_name,]
	qtl.egenes = assoc[assoc$gene_name %in% hic.egenes$gene,]
	if(nrow(qtl.egenes)==0){Nloci[ii] = 0; next}
	
	# output eQTL loci coordinates for eGenes in bait fragments
	qtl.egenes$snps = as.character(qtl.egenes$snps)
	chr = sapply(qtl.egenes$snps,function(x){strsplit(x,"_")[[1]][1]}) 
	snp = sapply(qtl.egenes$snps,function(x){strsplit(x,"_")[[1]][2]})
	start = as.numeric(snp) - haplo
	end = as.numeric(snp) + haplo
	#start = as.numeric(snp) - 1
	#end = as.numeric(snp) 
	ind = which(start < 0)
	if(length(ind)>0){start[ind] = 0} 
	out = data.frame(chr=chr, start=start, end=end, gene=qtl.egenes$gene_name, stringsAsFactors=F)
	fname = "eSNP.bed"
	write.table(out,fname,col.names=F,row.names=F,quote=F,sep="\t")

	# output capture coordinates with eGenes in bait
	out = hic.egenes[,c(4:6,1:3,7:9)]
	fname = "eGene.bed"
	write.table(out,fname,col.names=F,row.names=F,quote=F,sep="\t")

	# sort the bed files
	comm1 = "sort -k1,1 -k2,2n eSNP.bed > eSNP.sorted.bed"
	comm2 = "sort -k1,1 -k2,2n eGene.bed > eGene.sorted.bed"
	system(comm1)
	system(comm2)

	# intersect the filtered capture regions with the eGene loci ranges
	comm = "bedtools intersect -wo -a eGene.sorted.bed -b eSNP.sorted.bed > eGene.eSNP.bed"
	system(comm)

	# process final intersection results
	fname = "eGene.eSNP.bed"
	dat0 = read.table(fname,header=FALSE,stringsAsFactors=FALSE,sep="\t")
	dat0 = dat0[dat0$V8 == dat0$V13,]
	Nloci[ii] = length(unique(dat0$V13))

	# remove files
	system("rm eSNP.bed eGene.bed eSNP.sorted.bed eGene.sorted.bed eGene.eSNP.bed")

}

save(Nloci, file="Nloci.1bp.RData")
save(Nloci, file="Nloci.60kb.RData")
save(Nloci, file="Nloci.100kb.RData")
save(Nloci, file="Nloci.150kb.RData")



###################
## evaluate results
###################

load("Nloci.60kb.RData")
frac.loci = 100*Nloci/2408
mean(frac.loci)
min(frac.loci)
max(frac.loci)

Pperm = length(which(frac.loci > 1.6)) / 1000




