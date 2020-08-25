

library(coloc)
library(snpStats)

##########################################################
## import eqtl snp data and get MAF
##########################################################

library(MatrixEQTL)
library(dplyr)

# main directory
dir = "/media/wa3j/DATAPART1/Documents/adipose_sex_reviews/coloc"
setwd(dir)

# data directory
datdir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl"
setwd(datdir)

# eQTL analysis parameters/annotation
SNP_file_name = "gtex_snp_mat.txt"

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";
snps$fileOmitCharacters = "NA";
snps$fileSkipRows = 1;
snps$fileSkipColumns = 1;
snps$fileSliceSize = 2000;
snps$LoadFile(SNP_file_name);

# iterate through slices and get MAF
nSlice = snps$nSlices()
snp.ids = rownames(snps)
slicelen = 2000
freqs = data.frame(snp.ids=snp.ids,freq1=0,freq2=0,maf=0)
for(ii in 1:nSlice){
  inds = ((ii-1)*slicelen+1):(ii*slicelen)
  if(ii == nSlice){inds = ((ii-1)*slicelen+1):nrow(freqs)}
  X = as.data.frame(snps$getSlice(ii))
  fr = rowMeans(X)/2
  freqs$freq1[inds] = fr
}
freqs$freq2 = 1-freqs$freq1
freqs$maf = apply(freqs[,2:3],1,min)
freqs = freqs[-which(is.na(freqs$freq1)==TRUE),]

# main directory
dir = "/media/wa3j/DATAPART1/Documents/adipose_sex_reviews/coloc"
setwd(dir)
save(freqs, file="freqs.RData")


##########################################################
## get LD data (R>0.6) for PEER eQTLs
##########################################################

library(MatrixEQTL)
library(dplyr)

# main directory
dir = "/media/wa3j/DATAPART1/Documents/adipose_sex_reviews/coloc"
setwd(dir)

# load eQTL association loci
fname="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/CIdata.cases.RData"
load(fname)

# load LD data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/LD06.txt"
lddat0 = read.table(fname,stringsAsFactors=F,header=T)

# remove non-snps from the ld data
snp.chk = function(snpids){
  nnuc = sapply(snpids,function(x){
    n1 = strsplit(x,"_")[[1]][3]
    n2 = strsplit(x,"_")[[1]][4]
    if(nchar(n1)==1 & nchar(n2)==1){
      return("snp")
    } else {return("indel")}
  })
  return(nnuc)
}
checkA = snp.chk(lddat0$snpA)
checkB = snp.chk(lddat0$snpB)
indA = which(checkA == "snp")
indB = which(checkB == "snp")
lddat = lddat0[intersect(indA,indB),]
length(unique(c(lddat$snpA,lddat$snpB)))

# get LD snps for every association snp
# note that the index snps are included in ldSNP
CIdata.ldsnps = c()
for(ii in 1:nrow(CIdata.cases)){
  dat0 = CIdata.cases[ii,]
  snp0 = dat0$snps
  inds1 = grep(snp0, lddat[,1])
  inds2 = grep(snp0, lddat[,2])
  ld1 = lddat[inds1,c(2,3)]
  ld2 = lddat[inds2,c(1,3)]
  new0 = dat0 %>% mutate(ldSNP=snp0, R2=1)
  new1 = dat0[rep(seq_len(nrow(dat0)), each = nrow(ld1)), ]
  new2 = dat0[rep(seq_len(nrow(dat0)), each = nrow(ld2)), ]
  new1 = cbind(new1, ld1)
  new2 = cbind(new2, ld2)
  names(new1)[9] = names(new2)[9] = "ldSNP"
  new = rbind(new0, new1, new2)
  new = new[order(new$R2,decreasing=T),]
  CIdata.ldsnps = rbind(CIdata.ldsnps, new)
}

# save(CIdata.ldsnps, file="CIdata.ldsnps.RData")


##########################################################
## identify traits with overlaps at LD R2>0.8
##########################################################

# GWAS catalog data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4/gwascat.reduced2.txt"
gwascat0 = read.table(fname,header=T,stringsAsFactors=F,sep="\t",quote="")
gwascat0 = gwascat0 %>% filter(P.VALUE < 5 * 10^(-8))

# eQTL data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/CIdata.ldsnps.hg19.RData"
load(fname)
CIdata.ldsnps.hg19 = CIdata.ldsnps.hg19 %>% filter(R2>0.8)

# genotype annotation
# cat GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt > gtex.geno.ann.txt
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/gtex.geno.ann.txt"
geno.ann = read.table(fname,stringsAsFactors=F,header=T)

# get rsids for LD SNPS
names(geno.ann)[1] = "ldSNP"
eqtl.ldsnps = merge(CIdata.ldsnps.hg19, geno.ann[,c(1,7)], by="ldSNP")
names(eqtl.ldsnps)[ncol(eqtl.ldsnps)] = "ldrsid"

# format gwas
gwascat = gwascat0
names(gwascat) = c("pmid","trait","rsid","pval")
ldrsid = sapply(gwascat$rsid,function(x){strsplit(x,"-")[[1]][1]})
gwascat = gwascat %>% mutate(ldrsid=ldrsid)

# look at the intersection with gwas
# save(all.gwas.snps, file="all.gwas.snps.RData")
all.gwas.snps = merge(eqtl.ldsnps, gwascat, by="ldrsid")

# output traits for manual analysis
all.gwas.traits = all.gwas.snps$trait %>% unique
fname = "all.gwas.traits.txt"
write.table(all.gwas.traits,fname,col.names=F,row.names=F,quote=F)
trait.assoc = all.gwas.snps[,c(1,3,6:14)]
fname = "TableS9.txt"
write.table(trait.assoc,fname,col.names=T,row.names=F,quote=F,sep="\t")

##########################################################
## import eqtl and gwas data
##########################################################

library(dplyr)

# main directory
dir = "/media/wa3j/DATAPART1/Documents/adipose_sex_reviews/coloc"
setwd(dir)

# load MAF data
load("freqs.RData")

# GWAS catalog data (p >= 9e-06)
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4/gwascat.reduced2.txt"
gwascat0 = read.table(fname,header=T,stringsAsFactors=F,sep="\t",quote="")

# eQTL data (CIdata.ldsnps)
file="CIdata.ldsnps.RData"
load(file)

# get rsids
names(geno.ann)[1] = "ldSNP"
eqtl.ldsnps = merge(CIdata.ldsnps, geno.ann[,c(1,7)], by="ldSNP")
names(eqtl.ldsnps)[ncol(eqtl.ldsnps)] = "ldrsid"

# get maf
names(freqs)[1] = "ldSNP"
eqtl.ldsnps = merge(eqtl.ldsnps, freqs[,c(1,4)], by="ldSNP")

# format gwas
gwascat = gwascat0
names(gwascat) = c("pmid","trait","rsid","pval")
ldrsid = sapply(gwascat$rsid,function(x){strsplit(x,"-")[[1]][1]})
gwascat = gwascat %>% mutate(ldrsid=ldrsid)

# integrate gwas and eqtl data and filter for overlapping traits
eqtl.gwas.data = merge(eqtl.ldsnps, gwascat, by="ldrsid")


##########################################################
## overlap associations and implement coloc
##########################################################

# coloc loop for all gwas traits/associations
# max(PP4) = 0.646491
coloc.res = c()
ngene.gwas = length(unique(eqtl.gwas.data$genename))
for(ii in 1:ngene.gwas){
  
  gene = unique(eqtl.gwas.data$genename)[ii]
  gene.dat = eqtl.gwas.data[eqtl.gwas.data$genename %in% gene,]
  traits = unique(gene.dat$trait)
  
  for(jj in 1:length(traits)){
    trait = traits[jj]
    datjj = gene.dat %>% filter(trait==(!!trait))
    p.eqtl = datjj$pvalue
    p.gwas = datjj$pval
    maf = datjj$maf
    if(length(maf)==1){next}
    
    res <- coloc.abf(dataset1=list(pvalues=p.eqtl,N=nrow(datjj),type="quant"),
                     dataset2=list(pvalues=p.gwas,N=nrow(datjj),type="quant"),
                     MAF=maf)
    PP4 = res$summary[6]
    summary = data.frame(gene=gene,trait=trait,snp=datjj$snps[1],PP4=PP4)
    coloc.res = rbind(coloc.res, summary)
  } # jj, trait loop

} # ii, gene loop
