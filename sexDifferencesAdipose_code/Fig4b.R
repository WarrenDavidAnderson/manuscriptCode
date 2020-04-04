
############################################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Figure 4b
############################################################################

library(dplyr)
library(gridExtra)
library(beanplot)

setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4")

############################################################################
# load data
############################################################################

# intersect eqtl data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4/intersect.eqtl.RData"
load(fname)

# union eqtl data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4/union.eqtl.RData"
load(fname)
union.eqtl$snp = as.character(union.eqtl$snp)
union.eqtl$gene = as.character(union.eqtl$gene)

# subq adipose DEG data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/all_fc_v8.txt"
deg0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F)

# METSIM correlation data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig3/metsimdegcor.RData"
load(fname)

# Pulit GWAS data
gwas0 = read.table("tables8.txt",header=T,stringsAsFactors=F,sep="\t")
unique(gwas0$Male.or.Female.specific)

# genotype annotation
# cat GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt > gtex.geno.ann.txt
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/gtex.geno.ann.txt"
geno.ann = read.table(fname,stringsAsFactors=F,header=T)

# LD data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/LD08.txt"
lddat0 = read.table(fname,stringsAsFactors=F,header=T)

# load eQTL classification data (CIdata.cases)
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/CIdata.cases.RData"
load(fname)

# download nat comm data
natcom = read.table("Table1.csv",header=T,stringsAsFactors=F,sep=",")

# GWAS catalog data
fname = "gwascat.reduced.txt"
gwascat0 = read.table(fname,header=T,stringsAsFactors=F,sep="\t",quote="")

# heritability data
fname = "heritabilityF.txt"
h2f = read.table(fname,header=F,stringsAsFactors=F,sep=" ",skip=1)
fname = "heritabilityM.txt"
h2m = read.table(fname,header=F,stringsAsFactors=F,sep=" ",skip=1)
names(h2f) = names(h2m) = c("id","Vg","Ve","Vp","h2","logL","logL0","LRT","df","Pval","n")

# data from sex-stratified eqtl analysis
load("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/CIdata.RData")

############################################################################
# Fig 4b - effect sizes and heritability
############################################################################

pdf("fig4b.pdf")
par(mfrow=c(2,4))

# eqtl effect sizes
bean1 = data.frame(beta=CIdata$coefF,type="eqtl f",stringsAsFactors=T)
bean2 = data.frame(beta=CIdata$coefM,type="eqtl m",stringsAsFactors=T)
beandatEffect = rbind(bean1,bean2)
col = list("red", "blue")
ylim = c(-max(abs(beandatEffect$beta)), max(abs(beandatEffect$beta)))
beanplot(beta~type, data=beandatEffect, side="both", ll=0, 
         col=col, beanlinewd=0, innerborder=NA, 
         ylim=ylim,
         what = c(FALSE, TRUE, FALSE, TRUE))

# heritability
h2 = merge(h2f[,c(1,5)], h2m[,c(1,5)], by="id")
mean(h2$h2.x)
mean(h2$h2.y)
median(h2$h2.x)
median(h2$h2.y)
cut = 0.5
h2p = h2[which(h2$h2.x>cut | h2$h2.y>cut),]
bean1 = data.frame(beta=h2p$h2.x,type="eqtl f",stringsAsFactors=T)
bean2 = data.frame(beta=h2p$h2.y,type="eqtl m",stringsAsFactors=T)
beandatEffect = rbind(bean1,bean2)
col = list("red", "blue")
ylim = c(-max(abs(beandatEffect$beta)), max(abs(beandatEffect$beta)))
beanplot(beta~type, data=beandatEffect, side="both", ll=0, 
         col=col, beanlinewd=0, innerborder=NA, 
         #ylim=ylim,
         what = c(FALSE, TRUE, FALSE, TRUE))

dev.off()

# supp heritability fig
library(LSD)
pdf("heritHeat.pdf")
par(mfrow=c(2,2))
cut=0.1
h2p = h2[which(h2$h2.x>cut | h2$h2.y>cut),]
x = log10(h2p$h2.y)
y = log10(h2p$h2.x)
xlab = "male log10 h2"
ylab = "female log10 h2"
heatscatter(x, y, add.contour=F, nlevels=10, xlab=xlab, ylab=ylab)
abline(0,1,lwd=3)
dev.off()

############################################################################
# Supplementary - norm method effect size correlations
############################################################################

pdf("supfig_eqtlCor.pdf")
par(mfrow=c(2,2))

# peer versus sva
x = union.eqtl$b_peer
y = union.eqtl$b_sva
cor(x,y) # 0.9916561
plot(x,y, ylab="SVA effect size", xlab="PEER effect size")

# peer versus no
x = union.eqtl$b_peer
y = union.eqtl$b_no
cor(x,y) # 0.8427661
plot(x,y, ylab="covariate effect size", xlab="PEER effect size")


# no versus sva
x = union.eqtl$b_no
y = union.eqtl$b_sva
cor(x,y) # 0.8327513
plot(x,y, ylab="SVA effect size", xlab="covariate effect size")

dev.off()

# table S6 (ref then alt)
tables6 = union.eqtl[,c(2,3,10,4:9)]
names(tables6)[1:3] = c("variant_id","ensid","gene")
tables6 = merge(tables6,geno.ann[,c(1,7)],by="variant_id")
tables6 = tables6[,c(10,1:9)]
names(tables6)[1] = "rsid"
write.table(tables6,"TableS6.txt",col.names=T,row.names=F,quote=F,sep="\t")

# table S7 (ref then alt)
tables7 = intersect.eqtl[,c(2,3,10,4:9)]
names(tables7)[1:3] = c("variant_id","ensid","gene")
tables7 = merge(tables7,geno.ann[,c(1,7)],by="variant_id")
tables7 = tables7[,c(10,1:9)]
names(tables7)[1] = "rsid"
write.table(tables7,"TableS7.txt",col.names=T,row.names=F,quote=F,sep="\t")

############################################################################
# basic processing of LD data
############################################################################

# incorporate rsids
g0 = geno.ann[geno.ann$variant_id %in% CIdata.cases$snps,c(1,7)]
names(g0)[1] = "snps"
CIdata.cases = merge(CIdata.cases, g0, by="snps")
names(CIdata.cases)[9] = "rsid"

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

# get LD snps for every association snp
# note that the index snps are included in ldSNP
# verify that ld snps have the same sign
upbnd = sort(table(lddat$snpB),decreasing=TRUE)[1]
upbound = nrow(union.eqtl) * upbnd 
union.eqtl.ldsnps = as.data.frame(matrix(NA, nrow=upbound, ncol=12))
i1 = 0
for(ii in 1:nrow(union.eqtl)){

    
  dat0 = union.eqtl[ii,]
  snp0 = dat0$snp
  inds1 = grep(snp0, lddat[,1])
  inds2 = grep(snp0, lddat[,2])
  ld1 = lddat[inds1,c(2,3)]
  ld2 = lddat[inds2,c(1,3)]
  new0 = dat0 %>% mutate(ldSNP=snp0, R2=1)
  new1 = dat0[rep(seq_len(nrow(dat0)), each = nrow(ld1)), ]
  new2 = dat0[rep(seq_len(nrow(dat0)), each = nrow(ld2)), ]
  new1 = cbind(new1, ld1)
  new2 = cbind(new2, ld2)
  names(new1)[11] = names(new2)[11] = "ldSNP"
  r1 = nrow(new1)
  r2 = nrow(new2)
  new = as.data.frame(matrix(0,nrow=(1+r1+r2),ncol=12))
  new[1,] = new0
  if(r1>0){new[2:(r1+1),] = new1}
  if(r2>0){new[(r1+2):nrow(new),] = new2}
  names(new) = names(new0)
  new = new[order(new$R2,decreasing=T),]
  inds = c( (i1 + 1) : (i1 + nrow(new)) )
  i1 = i1 + nrow(new)
  union.eqtl.ldsnps[inds,] = new
  
}
names(union.eqtl.ldsnps) = names(new)
ind.rem = which(is.na(union.eqtl.ldsnps$assoc)==T)
union.eqtl.ldsnps = union.eqtl.ldsnps[-ind.rem,]
# save(union.eqtl.ldsnps, file="union.eqtl.ldsnps.RData")
# load("union.eqtl.ldsnps.RData")


############################################################################
# intersect LD SNPs with Pulit and Rask-Andersen GWAS SNPs (all hg19)
############################################################################

# get rsids for LD SNPS
names(geno.ann)[1] = "ldSNP"
eqtl.ldsnps = merge(union.eqtl.ldsnps, geno.ann[,c(1,7)], by="ldSNP")
names(eqtl.ldsnps)[13] = "ldrsid"

# look at the intersect
intersect(eqtl.ldsnps$ldrsid, gwas0$SNP)
intersect(eqtl.ldsnps$ldrsid, natcom$Lead.SNP)

# output rsids for liftover
rsHMG = sapply(gwas0$SNP,function(x){strsplit(x,"rs")[[1]][2]})
rsNCM = sapply(natcom$Lead.SNP,function(x){strsplit(x,"rs")[[1]][2]})
write.table(rsHMG,"rsHMG.txt",col.names=F,row.names=F,quote=F)
write.table(rsNCM,"rsNCM.txt",col.names=F,row.names=F,quote=F)

# unix liftover code - no changes in rsids
# https://genome.sph.umich.edu/wiki/LiftOver#Lift_dbSNP_rs_numbers
cd /media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4
wget /net/fantasia/home/zhanxw/amd/analyze/verifyBamID/liftRsNumber.py
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/database/organism_data/RsMergeArch.bcp.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/database/organism_data/SNPHistory.bcp.gz
python liftRsNumber.py rsHMG.txt > output1.rs
python liftRsNumber.py rsNCM.txt > output2.rs
cut -f1 output1.rs | sort -u
cut -f1 output2.rs | sort -u


############################################################################
# check interaction eqtl ld snps against the entire gwas catalog (hg38)
############################################################################

library(dplyr)
setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4") 

# GWAS catalog data
fname = "gwascat.reduced2.txt"
gwascat0 = read.table(fname,header=T,stringsAsFactors=F,sep="\t",quote="")
gwascat0 = gwascat0 %>% filter(P.VALUE < 5 * 10^(-8))

# genotype annotation
# cat GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt > gtex.geno.ann.txt
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/gtex.geno.ann.txt"
geno.ann = read.table(fname,stringsAsFactors=F,header=T)

# associations with LD
load("union.eqtl.ldsnps.RData")

# get rsids for LD SNPS
names(geno.ann)[1] = "ldSNP"
eqtl.ldsnps = merge(union.eqtl.ldsnps, geno.ann[,c(1,7)], by="ldSNP")
names(eqtl.ldsnps)[13] = "ldrsid"

# format gwas
gwascat = gwascat0
names(gwascat) = c("pmid","trait","rsid","pval")
ldrsid = sapply(gwascat$rsid,function(x){strsplit(x,"-")[[1]][1]})
gwascat = gwascat %>% mutate(ldrsid=ldrsid)

# look at the intersection with gwas
# save(all.gwas.snps, file="all.gwas.snps.RData")
all.gwas.snps = merge(eqtl.ldsnps, gwascat, by="ldrsid")
rsid2 = sapply(all.gwas.snps$rsid,function(x){strsplit(x,"-")[[1]][1]})
all.gwas.snps = all.gwas.snps %>% mutate(rsid2=rsid2)
all.gwas.snps = unique(all.gwas.snps[,c(1:13,15,17)])

# output traits for manual analysis
all.gwas.traits = all.gwas.snps$trait %>% unique
fname = "all.gwas.traits.txt"
write.table(all.gwas.traits,fname,col.names=F,row.names=F,quote=F)
trait.assoc = all.gwas.snps[,c(1,3,6:14)]
fname = "TableS9.txt"
write.table(trait.assoc,fname,col.names=T,row.names=F,quote=F,sep="\t")


############################################################################
# process Pulit GWAS data
############################################################################

library(data.table)
library(dplyr)

# Pulit GWAS data
datdir = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4/pulitGWAS/"
fname = "bmi.giant-ukbb.meta-analysis.females.23May2018.txt"
Pulit.gwas.bmi.f = fread( paste0(datdir,fname) )
fname = "bmi.giant-ukbb.meta-analysis.males.23May2018.txt"
Pulit.gwas.bmi.m = fread( paste0(datdir,fname) )
fname = "whradjbmi.giant-ukbb.meta-analysis.females.23May2018.txt"
Pulit.gwas.whradjbmi.f = fread( paste0(datdir,fname) )
fname = "whradjbmi.giant-ukbb.meta-analysis.males.23May2018.txt"
Pulit.gwas.whradjbmi.m = fread( paste0(datdir,fname) )
fname = "whr.giant-ukbb.meta-analysis.females.23May2018.txt"
Pulit.gwas.whr.f = fread( paste0(datdir,fname) )
fname = "whr.giant-ukbb.meta-analysis.males.23May2018.txt"
Pulit.gwas.whr.m = fread( paste0(datdir,fname) )

# combine the data
d1 = Pulit.gwas.bmi.f[,c(1:5,7,9)] %>% as.data.frame(stringsAsFactors=F)
d2 = Pulit.gwas.bmi.m[,c(1:5,7,9)] %>% as.data.frame(stringsAsFactors=F)
d3 = Pulit.gwas.whradjbmi.f[,c(1:5,7,9)] %>% as.data.frame(stringsAsFactors=F)
d4 = Pulit.gwas.whradjbmi.m[,c(1:5,7,9)] %>% as.data.frame(stringsAsFactors=F)
d5 = Pulit.gwas.whr.f[,c(1:5,7,9)] %>% as.data.frame(stringsAsFactors=F)
d6 = Pulit.gwas.whr.m[,c(1:5,7,9)] %>% as.data.frame(stringsAsFactors=F)
rm(Pulit.gwas.bmi.f,Pulit.gwas.bmi.m,
   Pulit.gwas.whradjbmi.f,Pulit.gwas.whradjbmi.m,
   Pulit.gwas.whr.f,Pulit.gwas.whr.m)
names(d1)[6:7] = c("beta-f-bmi", "pval-f-bmi")
names(d2)[6:7] = c("beta-m-bmi", "pval-m-bmi")
names(d3)[6:7] = c("beta-f-whradjbmi", "pval-f-whradjbmi")
names(d4)[6:7] = c("beta-m-whradjbmi", "pval-m-whradjbmi")
names(d5)[6:7] = c("beta-f-whr", "pval-f-whr")
names(d6)[6:7] = c("beta-m-whr", "pval-m-whr")
datlist = list(d1,d2,d3,d4,d5,d6)
rm(d1,d2,d3,d4,d5,d6)
Pulit.gwas <- Reduce(
  function(x, y) merge(x, y, all = TRUE, 
    by=c("CHR","POS","SNP","Tested_Allele","Other_Allele")),
  datlist
)
rm(datlist)

# revise SNP annotation
get.rs = function(rs0=NULL){
  lapply(rs0,function(x){strsplit(x,":")[[1]][1]}) %>% unlist
}
ldrsid = get.rs(Pulit.gwas$SNP)
Pulit.gwas = Pulit.gwas %>% mutate(ldrsid = ldrsid)
rm(ldrsid)
# save(Pulit.gwas, file="Pulit.gwas.RData")
# load("Pulit.gwas.RData")

############################################################################
# Match Pulit GWAS data with GTEx eqtls and filter
############################################################################

load("union.eqtl.ldsnps.RData")
load("Pulit.gwas.RData")

# gtex ld snp allele info
ldref = sapply(union.eqtl.ldsnps$ldSNP,function(x){strsplit(x,"_")[[1]][3]})
ldalt = sapply(union.eqtl.ldsnps$ldSNP,function(x){strsplit(x,"_")[[1]][4]})
union.eqtl.ldsnps = union.eqtl.ldsnps %>% mutate(ldref=ldref, ldalt=ldalt)
names(geno.ann)[7] = "ldrsid"
names(geno.ann)[1] = "ldSNP"
eqtl.ldsnps = merge(union.eqtl.ldsnps, geno.ann[,c(1,7)], by="ldSNP")

# merge gtex with pulit (consider reverse allele coding)
names(Pulit.gwas)[c(4:5)] = c("ldref", "ldalt")
gtex.pulit.merge1 = merge(eqtl.ldsnps, Pulit.gwas, by=c("ldref", "ldalt", "ldrsid"))
names(Pulit.gwas)[c(5,4)] = c("ldref", "ldalt")
gtex.pulit.merge2 = merge(eqtl.ldsnps, Pulit.gwas, by=c("ldref", "ldalt", "ldrsid"))
gtex.pulit.merge2$`beta-f-bmi` = (-1) * gtex.pulit.merge2$`beta-f-bmi`
gtex.pulit.merge2$`beta-m-bmi` = (-1) * gtex.pulit.merge2$`beta-m-bmi`
gtex.pulit.merge2$`beta-f-whradjbmi` = (-1) * gtex.pulit.merge2$`beta-f-whradjbmi`
gtex.pulit.merge2$`beta-m-whradjbmi` = (-1) * gtex.pulit.merge2$`beta-m-whradjbmi`
gtex.pulit.merge2$`beta-f-whr` = (-1) * gtex.pulit.merge2$`beta-f-whr`
gtex.pulit.merge2$`beta-m-whr` = (-1) * gtex.pulit.merge2$`beta-m-whr`
gtex.pulit.merge = rbind(gtex.pulit.merge1, gtex.pulit.merge2)
rm(gtex.pulit.merge1, gtex.pulit.merge2)
# save(gtex.pulit.merge, file="gtex.pulit.mergeAll.RData")
# load("gtex.pulit.mergeAll.RData")

# filter associations at 0.05 for the min gwas assoc at each snp
pmin = apply(gtex.pulit.merge[,c(20,22,24,26,28,30)],1,function(x){min(x)})
gtex.pulit.merge = gtex.pulit.merge %>% mutate(pmin = pmin)
gtex.pulit.merge.filt = gtex.pulit.merge %>% filter(pmin < 0.05)

# consider one ld snp per association, min of all p-values
nassoc = length(unique(gtex.pulit.merge.filt$snp))
gtex.pulit.assoc = data.frame(matrix(0,nrow=nassoc,ncol=ncol(gtex.pulit.merge.filt)))
names(gtex.pulit.assoc) = names(gtex.pulit.merge.filt)
for(ii in 1:nassoc){
  snp.ii = unique(gtex.pulit.merge.filt$snp)[ii]
  ind.snp = which(gtex.pulit.merge.filt$snp == snp.ii)
  pvs = gtex.pulit.merge.filt$pmin[ind.snp]
  ind.min = which(pvs == min(pvs))[1]
  ind = ind.snp[ind.min]
  gtex.pulit.assoc[ii,] = gtex.pulit.merge.filt[ind,]
}

# save(gtex.pulit.assoc, file="gtex.pulit.assoc.RData")
# load("gtex.pulit.assoc.RData")

############################################################################
# generate sex-stratified eqtl data for all associations (PEER)
############################################################################

setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4")

library(MatrixEQTL)
library(dplyr)
library(stats)

# see gtex_eqtl
load("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4/gtexeqtlv2.RData")
nSlice = snps$nSlices()
sliceSize = snps$fileSliceSize
snp.ids = rownames(snps)
rm(resids_all,resids_pr,PC,pc_loadings,pc_scores,fit,fitdat)

# get covariate data for sex
covar = as.matrix(cvrt)
model = t(covar) %>% as.data.frame
sex = as.factor(model$GENDER)
names(sex) = rownames(model)

# unfiltered associations from PEER analysis
load("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/eqtl_data.RData")
head(eqtldat)
eqtldatfilt = eqtldat
rm(eqtls,eqtldat,fit,fitdat,inds)
rm(expr,expr0,expr1,plts,geneNames)

# nominal filter
eqtldatfilt = eqtldatfilt %>% filter(pvalue<0.001)
eqtldatfilt = unique(eqtldatfilt)

# retain only associations for SNPs QCd by GTEx
eqtldatfilt = eqtldatfilt[eqtldatfilt$snps %in% geno.ann$variant_id,]

# implement linear regression
vars = c("snp","sex","snp:sex")
ind.m = which(covar[nrow(covar),]==0)
ind.f = which(covar[nrow(covar),]==1)
namen0 = c("snp","sex","snp:sex")
namen1 = c("beta","pval","95ci")
namen2 = sapply(namen1,function(x){paste0(x,"_",namen0)}) %>% as.vector
namen = sapply(c("ld","lead"),function(x){paste0(x,"_",namen2)}) %>% as.vector

CIdata.unfilt = data.frame(matrix(0,nrow=nrow(eqtldatfilt),ncol=12))

for(ii in 1:nrow(eqtldatfilt)){

  # get expression data (resides based on above regression)
  trs.ii = eqtldatfilt$gene[ii]
  edat = resids[which(rownames(resids)==trs.ii),] %>% t 
  edat.ii = edat %>% t %>% as.vector
  names(edat.ii) = colnames(edat)
  
  # get genotype data
  snp.tgt = eqtldatfilt$snps[ii]
  snp.ind = which(snp.ids == snp.tgt)
  sl = floor(snp.ind / sliceSize)
  if(snp.ind / sliceSize == floor(snp.ind / sliceSize)){
    sl = sl-1
  }
  snp.gen = as.data.frame(snps$getSlice(sl+1))
  rw = c((1+sliceSize*sl) : (sliceSize*(sl+1)))
  if(sl == (nSlice-1)){
    rw = c((1+sliceSize*sl) : (sliceSize*sl+nrow(snp.gen)))
  }
  names(snp.gen) = snp.sms
  rownames(snp.gen) = snp.ids[rw]
  geno = snp.gen[rownames(snp.gen)==snp.tgt,] 
  geno.ii = geno %>% t %>% as.vector
  names(geno.ii) = names(geno)
  
  # all(colnames(edat.ii) == names(geno.ii))
  # all(colnames(edat.ii) == names(sex))
  
  # implement full regression analysis
  regdat = data.frame(expr=edat.ii, snp=geno.ii, sex=sex) 
  res.lm = lm(expr ~. + sex*snp, data=regdat)
  coefs = res.lm$coefficients
  pvals = summary(res.lm)$coefficients[,4]
  conf95 = confint(res.lm)
  
  # implement sex-specific regression analysis
  indM = which(sex==0)
  indF = which(sex==1)
  regdatM = data.frame(expr=edat.ii[indM], snp=geno.ii[indM]) 
  regdatF = data.frame(expr=edat.ii[indF], snp=geno.ii[indF]) 
  resM = lm(expr ~ snp, data=regdatM)
  resF = lm(expr ~ snp, data=regdatF)
  coefsM = resM$coefficients[2]
  confM = confint(resM)[2,]
  coefsF = resF$coefficients[2]
  confF = confint(resF)[2,]
  
  # pack data into the new frame
  new = data.frame(eqtldatfilt[ii,],coefM=coefsM,coefF=coefsF,
                   ci1M=confM[1],ci2M=confM[2],ci1F=confF[1],ci2F=confF[2])
  CIdata.unfilt[ii,] = new
  
  # track time
  report = 100*ii/nrow(eqtldatfilt)
  if(ii %% 10 == 0){
    cat("\r", round(report,1),"%")
  }
  
} # ii, eqtl regression loop
Sys.time() 
names(CIdata.unfilt) = names(new)
CIdata.unfilt$snps = eqtldatfilt$snps
CIdata.unfilt$gene = eqtldatfilt$gene
# save(CIdata.unfilt,file="CIdata.unfilt.RData")


############################################################################
# overlap sex-stratified eqtls and gwas 
# fig 4c
############################################################################

library(NMF)
library(dplyr)
library(beanplot)
setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4")

# load gtex/gwas data
load("gtex.pulit.assoc.RData")
load("CIdata.unfilt.RData")

#########################################
# classify unfiltered eQTLs
CIdata = CIdata.unfilt

# check if the CIs bracket zero
b0M = rep("no",nrow(CIdata))
b0F = rep("no",nrow(CIdata))
for(ii in 1:nrow(CIdata)){
  if(CIdata$ci1M[ii]<=0 & CIdata$ci2M[ii]>=0){
    b0M[ii] = "yes"
  }
  if(CIdata$ci1F[ii]<=0 & CIdata$ci2F[ii]>=0){
    b0F[ii] = "yes"
  }
}

# document cases with one zero and one non-zero slope (case 1)
ind1 = which(b0M=="no" & b0F=="yes")
ind2 = which(b0M=="yes" & b0F=="no")
case1 = rep("no",nrow(CIdata))
case1[c(ind1,ind2)] = "yes"

# document cases with opposing slopes (case 2)
case2 = rep("no",nrow(CIdata))
for(ii in 1:length(case2)){
  s1 = sign(CIdata$ci1M[ii]) * sign(CIdata$ci2M[ii])
  s2 = sign(CIdata$ci1F[ii]) * sign(CIdata$ci2F[ii])
  if(case1[ii]=="yes"){next}
  if(s1==1 & s2==1 & sign(CIdata$ci1M[ii]) != sign(CIdata$ci1F[ii])){
    case2[ii] = "yes"
  }
}

# document cases with non-overlapping same direction slopes (case 3)
case3 = rep("no",nrow(CIdata))
for(ii in 1:length(case3)){
  s1a = CIdata$ci1M[ii]>CIdata$ci1F[ii] & CIdata$ci1M[ii]>CIdata$ci2F[ii] 
  s1b = CIdata$ci2M[ii]>CIdata$ci1F[ii] & CIdata$ci2M[ii]>CIdata$ci2F[ii]
  s2a = CIdata$ci1F[ii]>CIdata$ci1M[ii] & CIdata$ci1F[ii]>CIdata$ci2M[ii] 
  s2b = CIdata$ci2F[ii]>CIdata$ci1M[ii] & CIdata$ci2F[ii]>CIdata$ci2M[ii]  
  if(case1[ii]=="yes" | case2[ii]=="yes"){next}
  if((s1a & s1b) | (s2a & s2b)){
    case3[ii] = "yes"
  }
}

# summarize cases
cases = rep(0,nrow(CIdata))
for(ii in 1:length(case3)){
  c1 = case1[ii]
  c2 = case2[ii]
  c3 = case3[ii]
  if(c1=="yes" & c2=="no" & c3=="no"){cases[ii] = 1}
  if(c1=="no" & c2=="yes" & c3=="no"){cases[ii] = 2}
  if(c1=="no" & c2=="no" & c3=="yes"){cases[ii] = 3}
}

# counts of each case
length(which(cases==1)) # 43821
length(which(cases==2)) # 36256
length(which(cases==3)) # 1301
length(which(cases==1))+length(which(cases==2))+length(which(cases==3)) # 81378

CIdata.unfilt = CIdata.unfilt %>% mutate(cases=cases)
#########################################


# subset data before merging
gwas = gtex.pulit.assoc[,c(1:4,6:7,14,15,19:30)]
gtex = CIdata.unfilt[,c(1:2,7,8,13)]
names(gtex)[1] = "snp"

# merge the data
# save(gwas.eqtl, file="gwas.eqtl.RData")
# load("gwas.eqtl.RData")
gwas.eqtl = merge(gwas, gtex, by=c("snp", "gene"))

# look at beta correlations
x = gwas.eqtl$coefF
y = gwas.eqtl$`beta-f-bmi`
plot(x,y)
cor(x,y) # 0.01319762
x = gwas.eqtl$coefM
y = gwas.eqtl$`beta-m-bmi`
plot(x,y)
cor(x,y) # -0.03615713
x = gwas.eqtl$coefF
y = gwas.eqtl$`beta-f-whradjbmi`
plot(x,y)
cor(x,y) # 0.01392124
x = gwas.eqtl$coefM
y = gwas.eqtl$`beta-m-whradjbmi`
plot(x,y)
cor(x,y) # -0.02956909
x = gwas.eqtl$coefF
y = gwas.eqtl$`beta-f-whr`
plot(x,y)
cor(x,y) # 0.01836171
x = gwas.eqtl$coefM
y = gwas.eqtl$`beta-m-whr`
plot(x,y)
cor(x,y) # -0.0495762

# identify instances of sign match
get.match = function(dat=NULL){
  m1 = sign(dat[,1]) * sign(dat[,3])
  m2 = sign(dat[,2]) * sign(dat[,4])
  ind1 = which(m1 == 1)
  ind2 = which(m2 == 1)
  ind.match = intersect(ind1, ind2)
  return(ind.match)
}

# scale the entire matrix by a specified cutoff
scl.max = function(plt=NULL, cut=NULL){
  sc = apply(plt,2,function(x){
    out = x
    ind.max = which(x > cut)
    ind.min = which(x < -1*cut)
    if(length(ind.max)>0){out[ind.max] = cut}
    if(length(ind.min)>0){out[ind.min] = -1*cut}
    return(out)
  })
  return(sc)
}

# function to generate a set of plots
# plot eqtl/gwas effect size for matches
make.plots = function(dat=NULL,fname=NULL, cw=80, ch=1){
  
  # basics
  dat0 = dat
  inds = get.match(dat)
  dat = dat[inds,]
  inds = order(dat[,1])
  overlap.percent = 100*length(inds)/nrow(dat0)
  pltE = scl.max(dat[inds,3:4], cut = 0.4) # F then M
  pltG = scl.max(dat[inds,1:2], cut = 0.005)
  var = c("darkgreen","purple","black")
  names(var) = c(1,2,3)
  ann_colors = list(var=var)
  annotation = data.frame(var = factor(dat$cases[inds], labels = c("1","2","3")))
  
  pdf(fname); #par(mfrow=c(2,2))
  aheatmap(pltE, breaks=0, Colv=NA, Rowv=NA, annRow=annotation, annColors=ann_colors, cellwidth=cw, cellheight=ch)
  aheatmap(pltG, breaks=0, Colv=NA, Rowv=NA, annRow=annotation, annColors=ann_colors, cellwidth=cw, cellheight=ch)
  
  # bean1 = data.frame(beta=dat[,1],type="gwas f",stringsAsFactors=T)
  # bean2 = data.frame(beta=dat[,2],type="gwas m",stringsAsFactors=T)
  # bean3 = data.frame(beta=dat[,3],type="eqtl f",stringsAsFactors=T)
  # bean4 = data.frame(beta=dat[,4],type="eqtl m",stringsAsFactors=T)
  # beandatG = rbind(bean1,bean2)
  # beandatE = rbind(bean3,bean4)
  # col = list("red", "blue")
  # beanplot(beta~type, data=beandatG, side="both", ll=0, 
  #          col=col, beanlinewd=0, innerborder=NA, 
  #          ylim=c(-max(abs(beandatG$beta)), max(abs(beandatG$beta))),
  #          what = c(FALSE, TRUE, TRUE, TRUE))
  # beanplot(beta~type, data=beandatE, side="both", ll=0,
  #          ylim=c(-max(abs(beandatE$beta)), max(abs(beandatE$beta))),
  #          col=col, beanlinewd=0, innerborder=NA,
  #          what = c(FALSE, TRUE, TRUE, TRUE))
  
  f = dat[,c(1,3)]
  m = dat[,c(2,4)]
  plot(f[,1], f[,2], ylab="eqtl f", xlab="gwas f")
  plot(m[,1], m[,2], ylab="eqtl m", xlab="gwas m")
  corf = cor.test(f[,1], f[,2])
  corm = cor.test(m[,1], m[,2])
  
  dev.off()
  
  return(list(overlap.percent=overlap.percent, corf=corf, corm=corm))
  
} # make.plots

# check total number of associations
assoc = paste0(gwas.eqtl$snp,"_",gwas.eqtl$geneid)
length(assoc) # 1350
length(unique(assoc)) # 1350

# plots and correlations for sign identical matches
dat = gwas.eqtl %>% select('beta-f-bmi', 'beta-m-bmi', coefF, coefM, cases)
cor.bmi = make.plots(dat,fname="BMI.pdf")
dat = gwas.eqtl %>% select('beta-f-whradjbmi', 'beta-m-whradjbmi', coefF, coefM, cases)
cor.whradjbmi = make.plots(dat,fname="WHRadjBMI.pdf")
dat = gwas.eqtl %>% select('beta-f-whr', 'beta-m-whr', coefF, coefM, cases)
cor.whr = make.plots(dat,fname="WHR.pdf")

############################################################################
# Fig 4c
############################################################################

library(dplyr)
library(ggplot2)
library(gridExtra)

setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4")

# multi color bar for all association cases
pdf("fig4d.pdf")
par(mfrow=c(2,3))
col = c("darkgreen","purple","black")
counts = data.frame(table(CIdata.cases$case))
counts = as.matrix(counts[,2])
counts = 100*counts/sum(counts)
barplot(counts,col=col)
dev.off()
