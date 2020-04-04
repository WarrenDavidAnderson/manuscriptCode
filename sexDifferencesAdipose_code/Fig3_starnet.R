
############################################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Figure 3
# STARNET analysis
############################################################################

library(NMF)
library(dplyr)
library(cocor)
library(limma)
library(impute)

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig3"
setwd(dir)


######################################################################################################
## import data and basic formatting
######################################################################################################

# basic function
head2<-function(obj,n=6){head(obj[,1:n])}
tail2<-function(obj,n=6){tail(obj[,1:n])}

# subq adipose DEG data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/all_fc_v8.txt"
deg0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F)

# import starnet phenotype data 
fname = "clinicaltraits.txt"
phen0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F)
phen0 = phen0[,-c(3,17)]
phen0 = phen0[,c(1,15,2:14)]
namen = c("ID","Sex","Triglyceride","Blood glucose","Diastolic blood pressure","Systolic blood pressure",
          "Total cholesterol","LDL cholesterol","Waist to hip ratio","Waist circumference",
          "Body mass index","C reactive protein","Glycated haemoglobin (HbA1c)","Hip circumference","HDL cholesterol")
cbind(names(phen0),namen)
names(phen0) = namen

# sample sizes
nF = length(which(phen0$Sex=="female"))
nM = length(which(phen0$Sex=="male"))

# import starnet expression data 
fname = "162genes_matrix.txt"
expr0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F)
expr0[,2:ncol(expr0)] = apply(expr0[,2:ncol(expr0)],2,as.numeric)

# checks for expression data
length(expr0$gene)
all(expr0$gene %in% deg0$gene)
hist(expr0[123,2:ncol(expr0)] %>% t)

# log2 normalize the data
expr.norm = expr0
expr.norm[,2:ncol(expr.norm)] = apply(expr.norm[,2:ncol(expr.norm)],2,function(x){log2(x+1)})
hist(expr.norm[123,2:ncol(expr.norm)] %>% t)

# average for duplicated genes
ave.dups = function(dat=NULL){
  genes = dat$gene
  out = matrix(0,nrow=length(unique(genes)),ncol=(ncol(dat)-1))
  for(ii in 1:length(unique(genes))){
    ind = which(dat$gene == unique(genes)[ii])
    if(length(ind) > 1){
      av = apply(dat[ind,2:ncol(dat)],2,mean)
    } else {
      av = dat[ind,2:ncol(dat)]
    }
    out[ii,] = data.matrix(av)
  }
  out = data.frame(out)
  names(out) = names(dat)[2:ncol(dat)]
  rownames(out) = unique(genes)
  return(out)
} # ave.dups
expr = ave.dups(dat=expr.norm)


# filter phenotype samples for NAs (remove samples with more than 50% NAs)
row.nas = apply(phen0,1,function(x){length(which(is.na(x)==TRUE))})
row.nas = row.nas / ncol(phen0)
ind.rem = which(row.nas > 0.5)
phen0 = phen0[-ind.rem,]

# filter phenotypes for NAs (remove phenotypes with more than 50% NAs)
col.nas = apply(phen0,2,function(x){length(which(is.na(x)==TRUE))})
col.nas = col.nas / nrow(phen0)
which(col.nas > 0.5) # not an issue

# organize pheno and expr data based on the deg data set
ind.deg = sapply(deg0$gene,function(x){which(rownames(expr)==x)})
expr = expr[ind.deg,]
all(rownames(expr) == deg0$gene)
ind.expr = sapply(phen0$ID,function(x){which(names(expr)==x)}) %>% unlist
expr = expr[,ind.expr]
all(names(expr) == phen0$ID)

# impute missing phenotype data
pheno.imputed = phen0
pheno.imputed[,3:15] = data.matrix(pheno.imputed[,3:15])
pheno.imputed[,3:15] = impute.knn(t(pheno.imputed[,3:15]), k=10, rowmax=0.5, colmax=0.8)$data %>% t

######################################################################################################
## adjust data for BMI
######################################################################################################

# adjust expression data for BMI
BMI = pheno.imputed$`Body mass index`
expr.bmi.adj = data.frame(matrix(0,nrow=nrow(expr),ncol=ncol(expr)))
names(expr.bmi.adj) = names(expr)
rownames(expr.bmi.adj) = rownames(expr)
for(ii in 1:nrow(expr)){
  y = expr[ii,] %>% as.numeric
  regdat = data.frame(y=y,x=BMI)
  res = lm(y~x,regdat)
  expr.bmi.adj[ii,] = res$residuals
}

# adjust phenotype data for BMI
BMI = pheno.imputed$`Body mass index`
phno.bmi.adj = data.frame(matrix(0,nrow=nrow(pheno.imputed),ncol=ncol(pheno.imputed)-3))
names(phno.bmi.adj) = names(pheno.imputed)[c(3:10,12:15)]
rownames(phno.bmi.adj) = pheno.imputed$ID
ind = 1
for(ii in c(3:10,12:15)){
  y = pheno.imputed[,ii] %>% as.numeric
  regdat = data.frame(y=y,x=BMI)
  res = lm(y~x,regdat)
  phno.bmi.adj[,ind] = res$residuals
  ind = ind + 1
}

# save.image("starnetCor.RData")
# load("starnetCor.RData")

######################################################################################################
## gene expression test set analysis
######################################################################################################

# function for DEG analysis
deg.analysis = function(dat_Null=NULL,dat_Alt=NULL,null=NULL,alt=NULL){
  
  # implement linear model analysis with eBayes and BH adjustments
  subQall = rbind(dat_Null, dat_Alt)
  subQsex = c(rep(null,nrow(dat_Null)), rep(alt,nrow(dat_Alt)))
  design <- model.matrix(~0+subQsex)
  colnames(design) <- c(null,alt)
  contrast <- makeContrasts(Female - Male, levels = design)
  fit <- lmFit(t(subQall), design)
  fit <- contrasts.fit(fit, contrast) %>% eBayes
  output0 <- topTable(fit,number=ncol(subQall),adjust.method="BH")
  
  # convert log2FC to foldchange increase/decrease
  output = output0 %>% mutate(FC = 2^logFC)
  fc = rep(1,nrow(output))
  for (ii in 1:nrow(output)){
    if(output$FC[ii] > 1){fc[ii] = output$FC[ii]}
    if(output$FC[ii] < 1){fc[ii] = -1/output$FC[ii]}
  }
  output = output %>% mutate(ratioFC = fc) %>% mutate(absFC = abs(ratioFC))
  rownames(output) = rownames(output0)
  
  return(output)
} # deg.analysis

# sex annotation
sex = pheno.imputed$Sex 

# specify groups for comparison
ind_F = which(sex == "female")
ind_M = which(sex == "male")
null <- "Female"
alt <- "Male"
expr_dat = expr.bmi.adj
dat_Null = t(expr_dat[,ind_F]) # Null model - female
dat_Alt = t(expr_dat[,ind_M]) # Alt model - male

# implement deg analysis for starnet
output = deg.analysis(dat_Null=dat_Null,dat_Alt=dat_Alt,null=null,alt=alt)
starnet.deg = output %>% mutate(gene = rownames(output)) %>% select(gene,logFC,P.Value,ratioFC)

# merge with 3data degs
starnet.merged = merge(starnet.deg, deg0, by="gene")

# evaluate similarity of results
sign.prod = sign(starnet.merged$ratioFC * starnet.merged$gtex_fc)
ind.mismatch.sign = which(sign.prod == -1)
ind.mismatch.pval = which(starnet.merged$P.Value >= 0.05)
starnet.merged$gene[ind.mismatch.pval]
ind.not = union(ind.mismatch.sign,ind.mismatch.pval)
replicated = starnet.merged[-ind.not,]

# lookup FADS1
starnet.merged[47,]

######################################################################################################
## evaluate phenotype correlations
## partial correlation
######################################################################################################

# partial correlation code
# http://www.yilab.gatech.edu/pcor.html
source("pcor.R")

# sex annotation
sex = pheno.imputed$Sex

# female phenotype correlations
ind = which(sex == "female")
dat = phno.bmi.adj[ind,]
cor.mat.f = cor(dat)
hm.f = aheatmap(cor.mat.f, breaks=0)

# male phenotype correlations
ind = which(sex == "male")
dat = phno.bmi.adj[ind,]
cor.mat.m = cor(dat)
hm.m = aheatmap(cor.mat.m, breaks=0)

# all phenotype correlations
dat = phno.bmi.adj
cor.mat.a = cor(dat)
hm.a = aheatmap(cor.mat.m, breaks=0)

# plot heatmaps according to female organization
rr = hm.f$rowInd
cc = hm.f$colInd
pdf("FigS8a.pdf")
aheatmap(cor.mat.f[rr,cc], Colv=NA, Rowv=NA, breaks=0)
aheatmap(cor.mat.m[rr,cc], Colv=NA, Rowv=NA, breaks=0)
aheatmap(cor.mat.a[rr,cc], Colv=NA, Rowv=NA, breaks=0)
dev.off()

######################################################################################################
## evaluate gene-phenotype partial correlations
######################################################################################################

# data storage
cors.coef.m = matrix(0,nrow=nrow(expr.bmi.adj),ncol=ncol(phno.bmi.adj)) %>% as.data.frame
cors.pval.m = matrix(0,nrow=nrow(expr.bmi.adj),ncol=ncol(phno.bmi.adj)) %>% as.data.frame
cors.coef.f = matrix(0,nrow=nrow(expr.bmi.adj),ncol=ncol(phno.bmi.adj)) %>% as.data.frame
cors.pval.f = matrix(0,nrow=nrow(expr.bmi.adj),ncol=ncol(phno.bmi.adj)) %>% as.data.frame
cors.coef.a = matrix(0,nrow=nrow(expr.bmi.adj),ncol=ncol(phno.bmi.adj)) %>% as.data.frame
cors.pval.a = matrix(0,nrow=nrow(expr.bmi.adj),ncol=ncol(phno.bmi.adj)) %>% as.data.frame
rownames(cors.coef.m) = rownames(cors.pval.m) = rownames(expr.bmi.adj)
rownames(cors.coef.f) = rownames(cors.pval.f) = rownames(expr.bmi.adj)
rownames(cors.coef.a) = rownames(cors.pval.a) = rownames(expr.bmi.adj)
names(cors.coef.m) = names(cors.pval.m) = names(phno.bmi.adj)
names(cors.coef.f) = names(cors.pval.f) = names(phno.bmi.adj)
names(cors.coef.a) = names(cors.pval.a) = names(phno.bmi.adj)

# loop through data and compute correlations
for(ii in 1:nrow(expr.bmi.adj)){
  for(jj in 1:ncol(phno.bmi.adj)){
    
    # annotation
    gene = rownames(expr.bmi.adj)[ii]
    trait = names(phno.bmi.adj)[jj] 
    
    # male correlations
    ind = which(sex == "male")
    x = expr.bmi.adj[ii,ind] %>% as.numeric
    y = phno.bmi.adj[ind,jj]
    z = phno.bmi.adj[ind,-jj]
    res = pcor.test(x,y,z,method="pearson")
    cors.coef.m[ii,jj] = res$estimate
    cors.pval.m[ii,jj] = res$p.value
    
    # female correlations
    ind = which(sex == "female")
    x = expr.bmi.adj[ii,ind] %>% as.numeric
    y = phno.bmi.adj[ind,jj]
    z = phno.bmi.adj[ind,-jj]
    res = pcor.test(x,y,z,method="pearson")
    cors.coef.f[ii,jj] = res$estimate
    cors.pval.f[ii,jj] = res$p.value
    
    # all correlations
    x = expr.bmi.adj[ii,] %>% as.numeric
    y = phno.bmi.adj[,jj]
    z = phno.bmi.adj[,-jj]
    res = pcor.test(x,y,z,method="pearson")
    cors.coef.a[ii,jj] = res$estimate
    cors.pval.a[ii,jj] = res$p.value
    
  } # phenotype loop
  
} # gene loop

# function to clip matrix
matrixClip = function(mat,cut=0.2){
  mat1 = apply(mat,2,function(x){
    xout = x
    ind_hi = which(x > cut)
    ind_lo = which(x < -1*cut)
    if(length(ind_hi)>0){xout[ind_hi] = cut}
    if(length(ind_lo)>0){xout[ind_lo] = -1*cut}
    return(xout)
  })
  return(mat1)
}

# function to scale matrix by col max (max for for each phenotype)
scalemat <- function(mat,n=2){ 
  ranges = apply(mat,2,function(x){max(abs(x))})
  newmat = mat
  for(ii in 1:length(ranges)){
    newmat[,ii] = mat[,ii] / ranges[ii]
  }
  return(newmat)
}

# initial heatmaps
plt = matrixClip( scalemat(cors.coef.f), cut=0.7)
hm.gpcor.f = aheatmap(plt, breaks=0)
plt = matrixClip( scalemat(cors.coef.m), cut=0.7)
hm.gpcor.m = aheatmap(plt, breaks=0)
plt = matrixClip( scalemat(cors.coef.a), cut=0.7)
hm.gpcor.a = aheatmap(plt, breaks=0)

# heatmaps according to female organization
rr = hm.gpcor.f$rowInd
cc = hm.gpcor.f$colInd
pdf("FigS8b.pdf")
plt = matrixClip( scalemat(cors.coef.f[rr,cc]), cut=0.7)
aheatmap(plt, breaks=0, Colv=NA, Rowv=NA)
plt = matrixClip( scalemat(cors.coef.m[rr,cc]), cut=0.7)
aheatmap(plt, breaks=0, Colv=NA, Rowv=NA)
plt = matrixClip( scalemat(cors.coef.a[rr,cc]), cut=0.7)
aheatmap(plt, breaks=0, Colv=NA, Rowv=NA)
dev.off()

######################################################################################################
## filter correlations based on FDRs
######################################################################################################

# get column-wise fdrs (computed for each phenotype)
get.fdr = function(dat=NULL){
  out = data.frame(matrix(0,nrow(dat),ncol(dat)))
  names(out) = names(dat)
  rownames(out) = rownames(dat)
  for(ii in 1:ncol(out)){
    out[,ii] = p.adjust(dat[,ii], method="BH")
  }
  return(out)
} # get.fdr
cors.fdr.f = get.fdr(dat = cors.pval.f)
cors.fdr.m = get.fdr(dat = cors.pval.m)
cors.fdr.a = get.fdr(dat = cors.pval.a)

# select genes with at least 1 phenotype cor at fdr < 0.05
min.fdr.f = apply(cors.fdr.f, 1, min)
min.fdr.m = apply(cors.fdr.m, 1, min)
indf = which(min.fdr.f < 0.05)
indm = which(min.fdr.m < 0.05)
ind = union(indf, indm)
filt.cors.coef.f = cors.coef.f[ind,]
filt.cors.coef.m = cors.coef.m[ind,]
filt.cors.coef.a = cors.coef.a[ind,]

100* length(ind) / nrow(cors.coef.f)

# data for final heatmaps
plt.f = matrixClip( scalemat(filt.cors.coef.f), cut=0.7)
plt.m = matrixClip( scalemat(filt.cors.coef.m), cut=0.7)

# output supplementary table 5
out = filt.cors.coef.a %>% mutate(gene = rownames(filt.cors.coef.a))
write.table(out,"TableS5.txt",col.names=T,row.names=F,quote=F,sep="\t")

# data for final heatmaps
plt.f = matrixClip( scalemat(filt.cors.coef.f), cut=0.7)
plt.m = matrixClip( scalemat(filt.cors.coef.m), cut=0.7)

# check specific genes
genes = c("FADS1","MAP1B","HSPA12A","CLIC6","MMD","PDZD2")
grep(genes[6], rownames(plt.f))

######################################################################################################
## get overlap between STARNET correlation data and subq DEG data
######################################################################################################

# isolate genes with 3-way differential expression and a significant correlation
corgenes0 = rownames(filt.cors.coef.a)
deg_cor_genes = intersect(deg0$gene, corgenes0)
length(deg_cor_genes) # 149

# binom.test for more than 50% correlated out of 162
n = nrow(deg0)
x = length(deg_cor_genes) 
res = binom.test(x, n, p=0.5, alternative="two.sided")
res$p.value


######################################################################################################
## fig 3 heatmaps
######################################################################################################

# deg direction annotation (row annotation)
# http://nmf.r-forge.r-project.org/demo-aheatmap.html

# annotation for differential expression direction
pltann = data.frame(gene=rownames(plt.f))
degann = merge(pltann, deg0, by="gene") 
degann = degann %>% mutate(dir = sign(gtex_fc))
ind = sapply(rownames(plt.f),function(x){which(degann$gene==x)})  
degann = degann[ind,]
direction = degann$dir
direction = gsub("-1","M",direction)
direction = gsub("1","F",direction)
ann = data.frame(dir = direction)
dir = c("red","blue")
names(dir) = c("F", "M")
ann_colors = list(dir=dir)

# set to zero differences less than a specified value
clip2 = function(mat,cut=0.2){
  mat1 = apply(mat,2,function(x){
    xout = x
    ind = which(abs(x) < cut)
    if(length(ind)>0){xout[ind] = 0}
    return(xout)
  })
  return(mat1)
} # clip2

# generate organization
hh = aheatmap(plt.f)
rr = hh$rowInd
cc = hh$colInd
diff = data.matrix(filt.cors.coef.f - filt.cors.coef.m)
diff = matrixClip( clip2(diff,cut=0.15), cut=0.3)
ann = data.frame(dir=ann[rr,])

# heatmap plot
pdf("Fig3a.pdf")
hm.f = aheatmap(plt.f[rr,cc], Colv=NA, Rowv=NA, annRow=ann, annColors=ann_colors, breaks=0)
hm.m = aheatmap(plt.m[rr,cc], Colv=NA, Rowv=NA, annRow=ann, annColors=ann_colors, breaks=0)
hm.d = aheatmap(diff[rr,cc], Colv=NA, Rowv=NA, annRow=ann, annColors=ann_colors, breaks=0)
dev.off()

######################################################################################################
## fig 8
######################################################################################################

# basic clustering, phenotype organization
hh = aheatmap(plt.f)
cc = hh$colInd

# get data for genes of interest
genes = c("FADS1","MAP1B","HSPA12A","CLIC6","MMD","PDZD2")
inds.gene = sapply(genes,function(x){which(rownames(cors.coef.m)==x)})
diff = cors.coef.f[inds.gene,cc] - cors.coef.m[inds.gene,cc]
diff = matrixClip(diff, cut=0.2)
plt.f = matrixClip(cors.coef.f[inds.gene,cc], cut=0.2)
plt.m = matrixClip(cors.coef.m[inds.gene,cc], cut=0.2)
fdr.f = cors.fdr.f[inds.gene,cc]
fdr.m = cors.fdr.m[inds.gene,cc]

# annotation for differential expression direction
pltann = data.frame(gene=rownames(plt.f))
degann = merge(pltann, deg0, by="gene") 
degann = degann %>% mutate(dir = sign(gtex_fc))
ind = sapply(rownames(plt.f),function(x){which(degann$gene==x)})  
degann = degann[ind,]
direction = degann$dir
direction = gsub("-1","M",direction)
direction = gsub("1","F",direction)
ann = data.frame(dir = direction)
dir = c("red","blue")
names(dir) = c("F", "M")
ann_colors = list(dir=dir)

# heatmap plot
pdf("Fig8a.pdf")
hm.f = aheatmap(plt.f, Colv=NA, Rowv=NA, annRow=ann, annColors=ann_colors, breaks=0)
hm.m = aheatmap(plt.m, Colv=NA, Rowv=NA, annRow=ann, annColors=ann_colors, breaks=0)
hm.d = aheatmap(diff, Colv=NA, Rowv=NA, annRow=ann, annColors=ann_colors, breaks=0)
dev.off()


# get select info for paper
gene = "MAP1B"
fdat = cors.coef.f[inds.gene,cc]
mdat = cors.coef.m[inds.gene,cc]
ind = which(rownames(fdat) == gene)

fdat[ind,]
fdr.f[ind,]

mdat[ind,]
fdr.m[ind,]




