
############################################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Figure 3
# AAGMEx analysis
############################################################################
library(NMF)
library(dplyr)
library(cocor)
library(limma)
library(impute)

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig3"
setwd(dir)

# load aagmex data
fname = "female.csv"
corF0 = read.table(fname,header=T,sep=",",stringsAsFactors=F)
fname = "male.csv"
corM0 = read.table(fname,header=T,sep=",",stringsAsFactors=F)

# subq adipose DEG data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/all_fc_v8.txt"
deg0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F)


######################################################################################################
## basic formatting
######################################################################################################

# basic function
head2<-function(obj,n=6){head(obj[,1:n])}
tail2<-function(obj,n=6){tail(obj[,1:n])}

# female data
vars = names(corF0)
ind.cor = grep("Pearson_partial_",vars)
ind.fdr = grep("corr_fdr",vars)
cors.coef.f0 = corF0[,ind.cor]
cors.fdr.f0 = corF0[,ind.fdr]

# male data
vars = names(corM0)
ind.cor = grep("Pearson_partial_",vars)
ind.fdr = grep("corr_fdr",vars)
cors.coef.m0 = corM0[,ind.cor]
cors.fdr.m0 = corM0[,ind.fdr]

# check organization and set rownames
all(corM0$Probe_Id == corF0$Probe_Id)
pheno.m.cor = sapply(names(cors.coef.m0),function(x){strsplit(x,"[.]")[[1]][1]})
pheno.m.fdr = sapply(names(cors.fdr.m0),function(x){strsplit(x,"[.]")[[1]][1]})
pheno.f.cor = sapply(names(cors.coef.f0),function(x){strsplit(x,"[.]")[[1]][1]})
pheno.f.fdr = sapply(names(cors.fdr.f0),function(x){strsplit(x,"[.]")[[1]][1]})
all(pheno.m.cor == pheno.m.fdr)
all(pheno.m.cor == pheno.f.cor)
all(pheno.m.cor == pheno.f.fdr)
rownames(cors.coef.m0) = rownames(cors.coef.f0) = corM0$Probe_Id

# select the probe with the min FDR for each trait
get.probe = function(cors=NULL,fdr=NULL,full=NULL){
  min.fdr = apply(fdr,1,min)
  genes = unique(full$Symbol)
  inds = rep(0, length(genes))
  for(ii in 1:length(genes)){
    ind = which(full$Symbol == genes[ii])
    inds[ii] = ind[which(min.fdr[ind] == min(min.fdr[ind]))[1]]
  }
  cor.filt = cors[inds,]
  fdr.filt = fdr[inds,]
  full.filt = full[inds,]
  rownames(cor.filt) = rownames(fdr.filt) = full.filt$Symbol
  out = list(cor.filt=cor.filt, fdr.filt=fdr.filt, full.filt=full.filt)
  return(out)
} # get.probe
corM = get.probe(cors=cors.coef.m0, fdr=cors.fdr.m0, full=corM0)
cors.coef.m = corM$cor.filt
cors.fdr.m = corM$fdr.filt
corF = get.probe(cors=cors.coef.f0, fdr=cors.fdr.f0, full=corF0)
cors.coef.f = corF$cor.filt
cors.fdr.f = corF$fdr.filt


######################################################################################################
## initial heatmaps
######################################################################################################

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

# heatmaps according to female organization
rr = hm.gpcor.f$rowInd
cc = hm.gpcor.f$colInd
plt = matrixClip( scalemat(cors.coef.f[rr,cc]), cut=0.7)
aheatmap(plt, breaks=0, Colv=NA, Rowv=NA)
plt = matrixClip( scalemat(cors.coef.m[rr,cc]), cut=0.7)
aheatmap(plt, breaks=0, Colv=NA, Rowv=NA)

# heatmaps according to female organization
rr = hm.gpcor.f$rowInd
cc = hm.gpcor.f$colInd
pdf("FigS8c.pdf")
plt = matrixClip( scalemat(cors.coef.f[rr,cc]), cut=0.7)
aheatmap(plt, breaks=0, Colv=NA, Rowv=NA)
plt = matrixClip( scalemat(cors.coef.m[rr,cc]), cut=0.7)
aheatmap(plt, breaks=0, Colv=NA, Rowv=NA)
dev.off()


######################################################################################################
## filter correlations based on FDRs
######################################################################################################

# select genes with at least 1 phenotype cor at fdr < 0.05
min.fdr.f = apply(cors.fdr.f, 1, min)
min.fdr.m = apply(cors.fdr.m, 1, min)
indf = which(min.fdr.f < 0.05)
indm = which(min.fdr.m < 0.05)
ind = union(indf, indm)
filt.cors.coef.f = cors.coef.f[ind,]
filt.cors.coef.m = cors.coef.m[ind,]

100* length(ind) / nrow(cors.coef.f)

# data for final heatmaps
plt.f = matrixClip( scalemat(filt.cors.coef.f), cut=0.7)
plt.m = matrixClip( scalemat(filt.cors.coef.m), cut=0.7)

# check specific genes
genes = c("FADS1","MAP1B","HSPA12A","CLIC6","MMD","PDZD2")
grep(genes[6], rownames(plt.f))

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
pdf("Fig3b.pdf")
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
pdf("Fig8b.pdf")
hm.f = aheatmap(plt.f, Colv=NA, Rowv=NA, annRow=ann, annColors=ann_colors, breaks=0)
hm.m = aheatmap(plt.m, Colv=NA, Rowv=NA, annRow=ann, annColors=ann_colors, breaks=0)
hm.d = aheatmap(diff, Colv=NA, Rowv=NA, annRow=ann, annColors=ann_colors, breaks=0)
dev.off()


# get select info for paper
gene = "FADS1"
fdat = cors.coef.f[inds.gene,cc]
mdat = cors.coef.m[inds.gene,cc]
ind = which(rownames(fdat) == gene)

fdat[ind,]
fdr.f[ind,]

mdat[ind,]
fdr.m[ind,]





