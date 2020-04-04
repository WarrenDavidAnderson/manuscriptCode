

################################################################
# test run heritability calculation
################################################################

# dir, local
cd /media/wa3j/DATAPART1/Documents/software/gcta_1.93.0beta

# dir, rivanna
cd /nv/vol192/civeleklab/warren/software/gcta_1.93.0beta

# calculate the genetic relationships between pairwise individuals
# data: test.bed, test.bim and test.fam (output: test.grm.bin)
./gcta64 --bfile test --autosome --maf 0.01 --make-grm --out test --thread-num 4

# estimating the variance explained by the SNPs (output: test.hsq)
./gcta64 --grm test --pheno test.phen --reml --out test --thread-num 4


################################################################
# Calculate GRM
################################################################

# filter the genotype data and obtain plink format
# cd /nv/vol192/civeleklab/warren/adipose_project/eqtlLD
# cd /m/civeleklab/civeleklab/warren/adipose_project/eqtlLD
PATH=$PATH:/h4/t1/apps/statgen/plink1.90b6.7
PATH=$PATH:/h4/t1/apps/statgen/bcftools-1.9
PATH=$PATH:/h4/t1/apps/statgen/gcta_1.93.0beta

# genotype file (from gtex eqtl analysis)
dir=/m/civeleklab/civeleklab/mete/GTEx/GTEx_v8/geno/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1
fname=GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz
file=${dir}/${fname}

# convert vcf to plink format
plink --vcf ${file}

# fix the bim
awk '{gsub(/^chr/,""); print}' plink.bim > plink0.bim
rm plink.bim 
mv plink0.bim plink.bim

# calculate the genetic relationships between pairwise individuals
# plink.grm.bin
gcta=/nv/vol192/civeleklab/warren/software/gcta_1.93.0beta/gcta64
${gcta} --bfile plink --autosome --maf 0.01 --make-grm --out plink --thread-num 10


################################################################
# parse sva-normalized gtex expression files into phenotype files 
# formatted for gcta
################################################################

# dir
cd /scratch/wa3j/v8/6

# load R
module load gcc/7.1.0 openmpi/3.1.4 R/3.5.3

# load and format expression data
fname = "gencode_gene_map.txt"
ann_gene0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F,quote = "")
fname = "genes_XYM.txt"
genes_XYM = read.table(fname,header=F,sep="\t",stringsAsFactors=F,quote = "")
fname = "gtex_subq_covars.txt"
demo0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F,quote = "")
fname = "gtex_subq_expr_sva.txt"
expr0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F,check.names=F,quote = "")

# get gene annotation for expressed transcripts
inds0 = match(rownames(expr0), ann_gene0$gene_id)
inds = inds0[!is.na(inds0)]
gene.map = (ann_gene0[inds,])
all(gene.map$gene_id == rownames(expr0))

# omit transcripts on sex chromosomes
ind.sex0 = match(genes_XYM[,1],gene.map$gene_name)
ind.sex = ind.sex0[!is.na(ind.sex0)]
expr1 = expr0[-ind.sex,]
gene.map = gene.map[-ind.sex,]

# specify groups for comparison
ind_F = which(demo0$GENDER == "Female")
ind_M = which(demo0$GENDER == "Male")
expr_dat = expr1
namen = sapply(names(expr_dat),function(x){
	sp = strsplit(x,"-")[[1]]
	out = paste0(sp[1:2],collapse="-")
})
names(expr_dat) = namen
dat_Null = t(expr_dat[,ind_F]) # Null model - female
dat_Alt = t(expr_dat[,ind_M]) # Alt model - male

# generate test file for heritability analysis
tst = cbind(rownames(dat_Null), rownames(dat_Null), dat_Null[,1])
write.table(tst,"plink.phen",col.names=F,row.names=F,quote=F,sep="\t")

# write expression data outputs
write.table(dat_Null,"h_exprF.txt",col.names=F,row.names=F,quote=F,sep="\t")
write.table(dat_Alt,"h_exprM.txt",col.names=F,row.names=F,quote=F,sep="\t")
out = cbind(rownames(dat_Null), rownames(dat_Null))
write.table(out,"h_idsF.txt",col.names=F,row.names=F,quote=F,sep="\t")
out = cbind(rownames(dat_Alt), rownames(dat_Alt))
write.table(out,"h_idsM.txt",col.names=F,row.names=F,quote=F,sep="\t")
all(colnames(dat_Alt) == colnames(dat_Null))
out = colnames(dat_Alt)
write.table(out,"h_genes.txt",col.names=F,row.names=F,quote=F,sep="\t")

# move data
comm = "cp *h_*.txt /nv/vol192/civeleklab/warren/adipose_project/eqtlLD"
system(comm)

# generate female expr files
comm="mkdir exprF"
system(comm)
setwd("exprF")

ids = cbind(rownames(dat_Null), rownames(dat_Null))
for(ii in 1:ncol(dat_Null)){
	namen = colnames(dat_Null)[ii]
	fname = paste0(namen,"_F.txt")
	out = cbind(ids, dat_Null[,ii])
	write.table(out,fname,col.names=F,row.names=F,quote=F,sep="\t")
	if(ii %% 100 == 0){cat(ii/ncol(dat_Null), '\n')}
}

# generate male expr files
setwd("/scratch/wa3j/v8/6")
comm="mkdir exprM"
system(comm)
setwd("exprM")

ids = cbind(rownames(dat_Alt), rownames(dat_Alt))
for(ii in 1:ncol(dat_Alt)){
	namen = colnames(dat_Alt)[ii]
	fname = paste0(namen,"_M.txt")
	out = cbind(ids, dat_Alt[,ii])
	write.table(out,fname,col.names=F,row.names=F,quote=F,sep="\t")
	if(ii %% 100 == 0){cat(ii/ncol(dat_Alt), '\n')}
}


################################################################
# trial run of heritability
################################################################

cd /nv/vol192/civeleklab/warren/adipose_project/eqtlLD
pheno=/scratch/wa3j/v8/6/plink.phen
gcta=/nv/vol192/civeleklab/warren/software/gcta_1.93.0beta/gcta64
${gcta} --grm plink --pheno ${pheno} --reml --out test1 --thread-num 4


################################################################
# run heritability across all genes
################################################################

# directory for analysis results
resdir=/nv/vol192/civeleklab/warren/adipose_project/eqtlLD/heritability
cd ${resdir}

# key directories and information for the analysis
dat=/nv/vol192/civeleklab/warren/adipose_project/eqtlLD
gcta=/nv/vol192/civeleklab/warren/software/gcta_1.93.0beta/gcta64
grm=${dat}/plink

# loop through all ens ids
# compute heritibility for females 
# ii=ENSG00000000419.12_F.txt
cd ${resdir}/resF
datsex=${dat}/exprF
gg=$(ls ${datsex})
for ii in ${gg}
do
pheno=${datsex}/${ii}
out0=$(echo ${ii} | awk -F".txt" '{print $1}')
out=${out0}_hres.txt
${gcta} --grm ${grm} --pheno ${pheno} --reml --out ${out} --thread-num 4
rm *.log*
done

# loop through all ens ids
# compute heritibility for males 
cd ${resdir}/resM
datsex=${dat}/exprM
gg=$(ls ${datsex})
for ii in ${gg}
do
pheno=${datsex}/${ii}
out0=$(echo ${ii} | awk -F".txt" '{print $1}')
out=${out0}_hres.txt
${gcta} --grm ${grm} --pheno ${pheno} --reml --out ${out} --thread-num 4
rm *.log*
done


################################################################
# summarize heritability results
################################################################

cd /nv/vol192/civeleklab/warren/adipose_project/eqtlLD/heritability


# specify sex
sex=F

# directory and data
datdir=/nv/vol192/civeleklab/warren/adipose_project/eqtlLD/heritability/res${sex}
files=$(ls ${datdir})

# results file
outfile=heritability${sex}.txt
printf '%s\n' 'id Vg Ve Vp VgVp logL logL0 LRT df Pval n' | tr ' ' '\t' > ${outfile}

for ii in ${files}
do
fname=${datdir}/${ii}
ensid=$(echo ${ii} | awk -F"_${sex}" '{print $1}')
vars=$(tail -n +2 ${fname} | cut -f2 | awk 'BEGIN { OFS="\t"; } { print }')
echo ${ensid} ${vars} > res
mv ${outfile} tmp
cat tmp res > ${outfile}
rm tmp res
done

################
## move for local analysis
to=/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4
from=wa3j@interactive.hpc.virginia.edu:/nv/vol192/civeleklab/warren/adipose_project/eqtlLD/heritability/*heritability*.txt
scp -r $from $to

