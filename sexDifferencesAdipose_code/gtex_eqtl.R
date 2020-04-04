
############################################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Figure 4,5
# eQTL data processing
############################################################################

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl"
setwd(dir)


############################################################################
# load peer-based eqtl data
############################################################################

library(MatrixEQTL)
library(dplyr)

# genotype annotation
# cat GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt > gtex.geno.ann.txt
geno.ann = read.table("gtex.geno.ann.txt",stringsAsFactors=F,header=T)

# load the eqtl results based on the PEER analysis
load("eqtl_data.RData")
head(eqtldat)

# eQTL analysis parameters/annotation
SNP_file_name = "gtex_snp_mat.txt"
expression_file_name = "expr_eqtl.txt"
covariates_file_name = "covar_eqtl.txt"
errorCovariance = numeric()
cisDist = 1e6

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";
snps$fileOmitCharacters = "NA";
snps$fileSkipRows = 1;
snps$fileSkipColumns = 1;
snps$fileSliceSize = 2000;
snps$LoadFile(SNP_file_name);

## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";
gene$fileOmitCharacters = "NA";
gene$fileSkipRows = 1;
gene$fileSkipColumns = 1;
gene$fileSliceSize = 2000;
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      
cvrt$fileOmitCharacters = "NA"; 
cvrt$fileSkipRows = 1;          
cvrt$fileSkipColumns = 1;       
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4")
# save.image("gtexeqtlv2.RData")
# load("gtexeqtl.RData")

############################################################################
# basic plots, metrics and genes
############################################################################

library(ggplot2)

## Plot the Q-Q plot of local and distant p-values
pdf("ciseqtlpeer_qqhist.pdf")
par(mfrow=c(2,2))
plot(eqtls,main="")
hist(eqtldat$pvalue,main="",xlab="p-value")
dev.off()

# association counts
min(eqtldat$pvalue)
pcut = 1e-4
length(which(eqtldat$pvalue < pcut))

# combine expression with genotype and sex annotation
combine.expr.geno.sex = function(expression=NULL, genotype=NULL, sex=NULL, male=NULL, female=NULL){
  sex[sex == male] = "male"
  sex[sex == female] = "female"
  out = data.frame(expression, genotype, sex) 
  out[,1:2] = apply(out[,1:2],2,function(x){data.matrix(x) %>% as.numeric})
  out$sex = as.factor(out$sex)
  return(out)
}

# generate a specific scatter plot
single.plot = function(dat=NULL,gene=NULL,rs=NULL,geno=NULL){
  ref0 = strsplit(geno,"_")[[1]][3]
  alt0 = strsplit(geno,"_")[[1]][4]
  ref = paste0(ref0,ref0)
  het = paste0(ref0,alt0)
  alt = paste0(alt0,alt0)
  dat = dat[is.na(dat$genotype)==F,]
  plt = ggplot(dat, aes(as.factor(genotype), expression)) + 
    geom_boxplot(aes(as.factor(genotype), expression,fill=sex)) + 
    geom_point() + 
    geom_jitter(position=position_jitter(0.2)) + 
    facet_grid(. ~ sex) +
    theme(strip.background = element_blank(),
          strip.text.y = element_blank()) +
    geom_smooth(data=dat,aes(genotype+1,expression),method='lm',col="white") +
    scale_fill_manual(values=c("red", "blue")) +
    xlab(paste0("genotype (",rs,")")) +
    ylab(paste0(gene," expression")) +
    scale_x_discrete(labels=c("0" = ref, "1" = het, "2" = alt))
  return(plt)
}

# get expression residuals, do not regress out sex
expr = as.matrix(gene)
covar = as.matrix(cvrt)
model = t(covar) %>% as.data.frame
fitdat = t(expr) %>% as.data.frame
model$Platform = as.factor(model$Platform)
model = model[,-ncol(model)]
fit = apply(fitdat,2,function(x){
  dat = as.data.frame(cbind(x,model))
  fit_regdat = lm(x~.,data=dat)
  return(fit_regdat$residuals)
})
resids = t(fit)
#save(resids, file="resids.RData")

# specify associations for plotting
assoc.ind = 1
snp.tgt = eqtldat$snps[assoc.ind]
snp.ids = rownames(snps)
snp.sms = colnames(snps)
nSlice = snps$nSlices()
sliceSize = snps$fileSliceSize
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
gt = snp.gen[rownames(snp.gen)==snp.tgt,] 
  
# get expression data
gene.tgt = eqtldat$gene[assoc.ind]
ind = which(rownames(resids) == gene.tgt)
ge = resids[ind,]
ge.ind = sapply(names(gt),function(x){which(names(ge)==x)})
ge = ge[ge.ind]
model = model[ge.ind,]
all(rownames(model) == names(ge))
all(names(ge) == names(gt))

# plot the data
rsid = geno.ann$rs_id_dbSNP151_GRCh38p7[which(geno.ann$variant_id == snp.tgt)]
geid = ann_gene0$gene_name[which(ann_gene0$gene_id == gene.tgt)]
sex = as.character(model$GENDER)
pltdat = combine.expr.geno.sex(expression=ge, genotype=t(gt)[,1], sex=sex, male=0, female=1)
plt = single.plot(dat=pltdat, gene=geid, rs=rsid, geno=as.character(snp.tgt))
print(plt)

plts[[3]] = plt

# print multiple
library(gridExtra)
pdf("eQTL1.pdf", onefile = FALSE)
marrangeGrob(grobs=plts, nrow=2, ncol=2, top=NULL)
dev.off()

############################################################################
# identify associations based on different normalization methods
############################################################################

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl"
setwd(dir)

library(MatrixEQTL)
library(greenbrown)
library(dplyr)
library(stats)

# load the eqtl results based on the PEER correction
load("eqtl_data.RData")
head(eqtldat)
eqtl.peer0 = eqtldat
remove(eqtldat)

# load the eqtl results based on the SVA correction
load("eqtl_data_S.RData")
head(eqtldat)
eqtl.sva0 = eqtldat
remove(eqtldat)

# load the eqtl results without correction
load("eqtl_data_NS.RData")
head(eqtldat)
eqtl.nonorm0 = eqtldat
remove(eqtldat)

setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4")


# filter associations
pthresh = 1e-4
eqtl.peer.pfilt = eqtl.peer0 %>% filter(pvalue < pthresh)
eqtl.sva.pfilt = eqtl.sva0 %>% filter(pvalue < pthresh)
eqtl.nonorm.pfilt = eqtl.nonorm0 %>% filter(pvalue < pthresh)

# retain only associations for SNPs QCd by GTEx
eqtl.peer.pfilt.geno = eqtl.peer.pfilt[eqtl.peer.pfilt$snps %in% geno.ann$variant_id,]
eqtl.sva.pfilt.geno = eqtl.sva.pfilt[eqtl.sva.pfilt$snps %in% geno.ann$variant_id,]
eqtl.nonorm.pfilt.geno = eqtl.nonorm.pfilt[eqtl.nonorm.pfilt$snps %in% geno.ann$variant_id,]


# reduce data to one association per gene
one.assoc = function(dat=NULL){
  out0 = dat[order(dat$pvalue),]
  out = out0[!duplicated(out0$gene),]
  return(out)
}
eqtl.peer.pfilt.uniq = one.assoc(dat=eqtl.peer.pfilt.geno)
eqtl.sva.pfilt.uniq = one.assoc(dat=eqtl.sva.pfilt.geno)
eqtl.nonorm.pfilt.uniq = one.assoc(dat=eqtl.nonorm.pfilt.geno)
save(eqtl.peer.pfilt.uniq, file="eqtl.peer.pfilt.uniq.RData")
save(eqtl.sva.pfilt.uniq, file="eqtl.sva.pfilt.uniq.RData")
save(eqtl.nonorm.pfilt.uniq, file="eqtl.nonorm.pfilt.uniq.RData")

# load("eqtl.peer.pfilt.uniq.RData")
# load("eqtl.sva.pfilt.uniq.RData")
# load("eqtl.nonorm.pfilt.uniq.RData")

# document associations
assoc = paste0(eqtl.peer.pfilt.uniq$snps,"_",eqtl.peer.pfilt.uniq$gene)
eqtl.peer.pfilt.uniq = eqtl.peer.pfilt.uniq %>% mutate(assoc = assoc)
assoc = paste0(eqtl.sva.pfilt.uniq$snps,"_",eqtl.sva.pfilt.uniq$gene)
eqtl.sva.pfilt.uniq = eqtl.sva.pfilt.uniq %>% mutate(assoc = assoc)
assoc = paste0(eqtl.nonorm.pfilt.uniq$snps,"_",eqtl.nonorm.pfilt.uniq$gene)
eqtl.nonorm.pfilt.uniq = eqtl.nonorm.pfilt.uniq %>% mutate(assoc = assoc)

# identify matching associations from all three normalization methods (intersect)
# note that identical effect direction is determined
a = eqtl.peer.pfilt.uniq$assoc
b = eqtl.sva.pfilt.uniq$assoc
c = eqtl.nonorm.pfilt.uniq$assoc
intersect.assoc = Reduce(intersect, list(a,b,c)) # 56
eqtl.peer.pfilt.uniq.intersect = eqtl.peer.pfilt.uniq[eqtl.peer.pfilt.uniq$assoc %in% intersect.assoc,c(1,2,4,6,7)]
eqtl.sva.pfilt.uniq.intersect = eqtl.sva.pfilt.uniq[eqtl.sva.pfilt.uniq$assoc  %in% intersect.assoc,c(1,2,4,6,7)]
eqtl.nonorm.pfilt.uniq.intersect = eqtl.nonorm.pfilt.uniq[eqtl.nonorm.pfilt.uniq$assoc %in% intersect.assoc,c(1,2,4,6,7)]
intersect.eqtl <- Reduce(
  function(x, y) merge(x, y, all=TRUE, by="assoc"),
  list(eqtl.peer.pfilt.uniq.intersect,eqtl.sva.pfilt.uniq.intersect,eqtl.nonorm.pfilt.uniq.intersect)
)[,c(1:5,8,9,12,13)]
names(intersect.eqtl)[2:9] = c("snp","gene","p_peer","b_peer","p_sva","b_sva","p_no","b_no")
intersect.eqtl[,4:9] = apply(intersect.eqtl[,4:9],2,function(x){
  data.matrix(x) %>% as.numeric
})
eq = apply(intersect.eqtl,1,function(x){
  s = sign( as.numeric(x[c(5,7,9)]) ) 
  return( AllEqual(s) )
})
ind = which(eq == TRUE)
intersect.eqtl = intersect.eqtl[ind,] # 56 / 56

# identify matching associations from all three normalization methods (union)
# note that identical effect direction is determined
union.assoc = unique(c(a,b,c)) # 4927
assoc = paste0(eqtl.peer0$snps,"_",eqtl.peer0$gene)
eqtl.peer.union = eqtl.peer0 %>% mutate(assoc = assoc)
assoc = paste0(eqtl.sva0$snps,"_",eqtl.sva0$gene)
eqtl.sva.union = eqtl.sva0 %>% mutate(assoc = assoc)
assoc = paste0(eqtl.nonorm0$snps,"_",eqtl.nonorm0$gene)
eqtl.nonorm.union = eqtl.nonorm0 %>% mutate(assoc = assoc)
eqtl.peer.union = eqtl.peer.union[eqtl.peer.union$assoc %in% union.assoc,]
eqtl.sva.union = eqtl.sva.union[eqtl.sva.union$assoc %in% union.assoc,]
eqtl.nonorm.union = eqtl.nonorm.union[eqtl.nonorm.union$assoc %in% union.assoc,]
union.eqtl <- Reduce(
  function(x, y) merge(x, y, all=TRUE, by="assoc"),
  list(eqtl.peer.union,eqtl.sva.union,eqtl.nonorm.union)
)[,c(1:3,5,7,11,13,17,19)]
names(union.eqtl)[2:9] = c("snp","gene","p_peer","b_peer","p_sva","b_sva","p_no","b_no")
union.eqtl[,4:9] = apply(union.eqtl[,4:9],2,function(x){
  data.matrix(x) %>% as.numeric
})
eq = apply(union.eqtl,1,function(x){
  s = sign( as.numeric(x[c(5,7,9)]) ) 
  return( AllEqual(s) )
})
ind = which(eq == TRUE)
union.eqtl = union.eqtl[ind,] # 4811 / 4927

# add gene ids 
getid = function(ens){
  gid = sapply(ens,function(x){
    ann_gene0$gene_name[which(ann_gene0$gene_id == x)]
  })
  return(gid)
}
geid.intersect = getid(intersect.eqtl$gene)
geid.union = getid(union.eqtl$gene)
intersect.eqtl = intersect.eqtl %>% mutate(geneid = geid.intersect)
union.eqtl = union.eqtl %>% mutate(geneid = geid.union)

# enforce a single association per gene and save
intersect.eqtl = intersect.eqtl[order(intersect.eqtl$p_peer),]
intersect.eqtl = intersect.eqtl[!duplicated(intersect.eqtl$assoc),]
union.eqtl = union.eqtl[order(union.eqtl$p_peer),]
union.eqtl = union.eqtl[!duplicated(union.eqtl$assoc),]
save(intersect.eqtl, file="intersect.eqtl.RData") # 4811
save(union.eqtl, file="union.eqtl.RData") # 56


############################################################################
# run regression for individual associations
############################################################################

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl"
setwd(dir)

library(MatrixEQTL)
library(dplyr)
library(stats)

# load resids
load("resids.RData")

# load("gtexeqtl.RData")
nSlice = snps$nSlices()
sliceSize = snps$fileSliceSize
snp.ids = rownames(snps)
snp.sms = colnames(snps)

# get covariate data for sex
covar = as.matrix(cvrt)
model = t(covar) %>% as.data.frame
sex = as.factor(model$GENDER)
names(sex) = rownames(model)

# isolate significant associations from PEER analysis
# one association per gene
load("eqtl_data.RData")
head(eqtldat)
eqtldatfilt = eqtldat %>% filter(pvalue < 1e-4)
eqtldatfilt = eqtldatfilt[order(eqtldatfilt$pvalue),]
eqtldatfilt = eqtldatfilt[!duplicated(eqtldatfilt$gene),]

# implement linear regression
CIdata = c()
vars = c("snp","sex","snp:sex")
ind.m = which(covar[nrow(covar),]==0)
ind.f = which(covar[nrow(covar),]==1)
namen0 = c("snp","sex","snp:sex")
namen1 = c("beta","pval","95ci")
namen2 = sapply(namen1,function(x){paste0(x,"_",namen0)}) %>% as.vector
namen = sapply(c("ld","lead"),function(x){paste0(x,"_",namen2)}) %>% as.vector
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
  CIdata = rbind(CIdata, new)
  
  # track time
  report = 100*ii/nrow(eqtldatfilt)
  if(ii %% 10 == 0){
    cat("\r", round(report,1),"%")
  }
  
} # ii, eqtl regression loop
# save(CIdata,file="CIdata.RData")

############################################################################
# classify the associations
############################################################################

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl"
setwd(dir)

library(MatrixEQTL)
library(dplyr)
library(stats)

# load("CIdata.RData")

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
CIdata.cases = CIdata[,c(1:6)] %>% mutate(case=cases)

# counts of each case
length(which(cases==1)) # 803
length(which(cases==2)) # 1584
length(which(cases==3)) # 21
length(which(cases==1))+length(which(cases==2))+length(which(cases==3)) # 2408

# incorporate gene ids
load("eqtl_data.RData")
ids = sapply(CIdata.cases$gene,function(x){ann_gene0$gene_name[ann_gene0$gene_id==x]})
CIdata.cases = CIdata.cases %>% mutate(genename = ids)
CIdata.cases[,1] = as.character(CIdata.cases[,1])
CIdata.cases[,2] = as.character(CIdata.cases[,2])


# save(CIdata.cases, file="CIdata.cases.RData")


############################################################################
# combine peer eqtl data with LD data
############################################################################

library(dplyr)

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl"
setwd(dir)

# genotype annotation
# cat GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt > gtex.geno.ann.txt
geno.ann = read.table("gtex.geno.ann.txt",stringsAsFactors=F,header=T)

# load LD data
fname = "LD06.txt"
lddat0 = read.table(fname,stringsAsFactors=F,header=T)

# load eQTL classification data (CIdata.cases)
load("CIdata.cases.RData")

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


############################################################################
# revise coordinates from hg38 to hg19 for compatibility with crossmap
# for compatibility with snp2tfbs data
############################################################################

# generate bed coords for the ldSNP column of CIdata.ldsnps
# note id format: chromosome_position_ref_alt_build
ref = sapply(CIdata.ldsnps$ldSNP,function(x){strsplit(x,"_")[[1]][3]})
alt = sapply(CIdata.ldsnps$ldSNP,function(x){strsplit(x,"_")[[1]][4]})
chr = sapply(CIdata.ldsnps$ldSNP,function(x){strsplit(x,"_")[[1]][1]})
end = sapply(CIdata.ldsnps$ldSNP,function(x){strsplit(x,"_")[[1]][2]})
end = as.numeric(end)
start = end - 1
bedhg38 = data.frame(chr=chr,start=start,end=end,ref=ref,alt=alt,id=CIdata.ldsnps$ldSNP,stringsAsFactors=F)

#########################################
## unix code
#########################################

# output data
fname = "bedhg38.bed"
write.table(bedhg38,fname,col.names=F,row.names=F,sep="\t",quote=F)

# move data
from=/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/bedhg38.bed
to=wa3j@interactive.hpc.virginia.edu:/nv/vol192/civeleklab/warren/adipose_project/eqtlLD
scp -r $from $to

# get the mapping file
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

# run crossmap
cd /nv/vol192/civeleklab/warren/adipose_project/eqtlLD
ml anaconda/5.2.0-py3.6
source activate crossmap
python CrossMap.py bed hg38ToHg19.over.chain.gz bedhg38.bed bedhg19.bed
source deactivate

# move data
from=wa3j@interactive.hpc.virginia.edu:/nv/vol192/civeleklab/warren/adipose_project/eqtlLD/bedhg19.bed
to=/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl
scp -r $from $to


#########################################
## load crossmap results back into R and integrate
#########################################

# load hg19 coordinates
bedhg19 = read.table("bedhg19.bed",header=F,stringsAsFactors=F)
names(bedhg19) = c("chr","start","end","ref","alt","ldSNP")

# integrate hg19 information with the ld data
CIdata.ldsnps.hg19 = merge(CIdata.ldsnps, bedhg19, by="ldSNP")

# save results for preparing fig 5
save(CIdata.ldsnps.hg19, file="CIdata.ldsnps.hg19.RData")
