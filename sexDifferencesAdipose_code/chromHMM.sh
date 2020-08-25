
dir=/media/wa3j/DATAPART1/Documents/adipose_sex_reviews/chromhmm
cd ${dir}

##########################################################
## roadmap information
##########################################################

# sample annotation, metadata
# https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15

# state description
# https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#core_15state

# samples
# 25: Adipose Derived Mesenchymal Stem Cell Cultured Cells
# 23: Mesenchymal Stem Cell Derived Adipocyte Cultured Cells

# download files
# https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/
E023_15_coreMarks_hg38lift_mnemonics.bed.gz
E025_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip *.gz

# rename files
cat E023_15_coreMarks_hg38lift_mnemonics.bed > adipo.bed
cat E025_15_coreMarks_hg38lift_mnemonics.bed > pread.bed


##########################################################
## generate SNP bed file for the eQTL data in R
##########################################################

library(dplyr)

# load eQTL association loci
fname="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/CIdata.cases.RData"
load(fname)

# generate bed file
chr = sapply(CIdata.cases$snps, function(x){
	strsplit(x,"_")[[1]][1]
})
end = sapply(CIdata.cases$snps, function(x){
	strsplit(x,"_")[[1]][2] %>% as.numeric
})
start = end - 1
gene = CIdata.cases$genename
snp = CIdata.cases$snps
assoc = paste0(snp,"_",gene)
bed = data.frame(chr=chr,start=start,end=end,
				snp=snp,gene=gene,assoc=assoc,
				stringsAsFactors=F)

# output data
fname = "eqtl.bed"
write.table(bed,fname,col.names=F,row.names=F,sep="\t",quote=F)

# sort the eqtl data
comm = "sort -k1,1 -k2,2n eqtl.bed > eqtl.sorted.bed"
system(comm)


##########################################################
## intersect eQTL data with chromHMM data in R
##########################################################

# intersect with preadipocyte data
comm = "bedtools intersect -a pread.bed -b eqtl.sorted.bed -wa -wb > eqtl.pread.bed"
system(comm)

# intersect with adipocyte data
comm = "bedtools intersect -a adipo.bed -b eqtl.sorted.bed -wa -wb > eqtl.adipo.bed"
system(comm)

# import results into R
fname = "eqtl.pread.bed"
eqtl.pread = read.table(fname,header=F,stringsAsFactors=F)
names(eqtl.pread) = c("chr_hmm","start_hmm","end_hmm","hmm_ann",
					"chr_snp","start_snp","end_snp",
					"snp_id","gene","assoc")
fname = "eqtl.adipo.bed"
eqtl.adipo = read.table(fname,header=F,stringsAsFactors=F)
names(eqtl.adipo) = c("chr_hmm","start_hmm","end_hmm","hmm_ann",
					"chr_snp","start_snp","end_snp",
					"snp_id","gene","assoc")

# summary
dim(eqtl.pread) # 2379
length(unique(eqtl.pread$assoc)) # 2379
dim(eqtl.adipo) # 2379
length(unique(eqtl.adipo$assoc)) # 2379



##########################################################
## summary plots
##########################################################

# convert hmm ids to single numbers
num = sapply(eqtl.pread$hmm_ann,function(x){strsplit(x,"_")[[1]][1]})
eqtl.pread = eqtl.pread %>% mutate(hmmnum = num)
num = sapply(eqtl.adipo$hmm_ann,function(x){strsplit(x,"_")[[1]][1]})
eqtl.adipo = eqtl.adipo %>% mutate(hmmnum = num)

# generate annotation key
TSS = c(1,2,10)
transcription = c(3,4,5)
enhancer = c(6,7,12)
repressed = c(9,13,14)
quiescent = c(15)
other = c(8,11)

# function to get percentages of each annotation
get.percent.ann = function(dat=NULL){
	per.TSS = 100 * length(dat[dat %in% TSS]) / length(dat)
	per.trs = 100 * length(dat[dat %in% transcription]) / length(dat)
	per.enh = 100 * length(dat[dat %in% enhancer]) / length(dat)
	per.rep = 100 * length(dat[dat %in% repressed]) / length(dat)
	per.sil = 100 * length(dat[dat %in% quiescent]) / length(dat)
	per.oth = 100 * length(dat[dat %in% other]) / length(dat)
	percent = c(per.TSS,per.trs,per.enh,per.rep,per.sil,per.oth)
	annotation = factor(c("tss","transcr","enhancer","repress","quiescent","other"),
					levels=c("tss","transcr","enhancer","repress","quiescent","other"))
	out = data.frame(percent=percent,annotation=annotation) 
	return(out)
}

hmm.adipo = get.percent.ann(dat=eqtl.adipo$hmmnum)
hmm.pread = get.percent.ann(dat=eqtl.pread$hmmnum)

# generate pie charts
library(ggplot2)
library(gridExtra)
values = c(	rgb(255,0,0, maxColorValue = 255),
			rgb(0,128,0, maxColorValue = 255),
			rgb(255,255,0, maxColorValue = 255),
			rgb(128,128,128, maxColorValue = 255),
			rgb(255,255,255, maxColorValue = 255),
			rgb(0,0,0, maxColorValue = 255))
plt.gen = function(dat=NULL){
	plt = ggplot(dat, aes(x="", y=percent, fill=annotation)) +
	geom_bar(width = 1, stat = "identity") +
	coord_polar("y", start=0) +
	scale_fill_manual(values=values)
}
plt.pread = plt.gen(hmm.pread)
plt.adipo = plt.gen(hmm.adipo)
plt = list(plt.pread=plt.pread, plt.adipo=plt.adipo)

pdf("pies.pdf", onefile=TRUE, height=12, width=14)
marrangeGrob(grobs=plt, nrow=3, ncol=2, top=NULL)
dev.off()







