
############################################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Figure 6
############################################################################

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig6"
setwd(dir)

library(Rcpp)
library(reshape)
library(dplyr)
library(ggplot2)
library(NMF) 
library(limma)
library(VennDiagram)
library(WGCNA)
library(gridExtra)
library(stats)
library(data.table)

###########################################################
## read and format data
############################################################

fname = "HMDP_Adipose_Gene_Expression.csv"
expr0 = read.table(fname,sep=",",header=T,stringsAsFactors=F)

fname = "HMDP_Traits.csv"
pheno0 = read.table(fname,sep=",",header=T,stringsAsFactors=F)

# initial annotation for the phenotype set
annotation = pheno0 %>% select(mouse_number,Strain,Sex)

# function for finding a given annotation associated with a given vector
var.annot = function(vecdat=NULL, annot=NULL, invar=NULL, outvar=NULL){
  ann.in = ( annot %>% select(invar) )[,1]
  ann.out = ( annot %>% select(outvar) )[,1]
  out = sapply(vecdat,function(x){ann.out[which(ann.in==x)][1]}) %>% unlist %>% unname
  return(out)
}

# set exression matrix and annotation
# main objects: 
##### expr2: gene expression matrix
##### expr2.samp.annotation: column/sample annotation for expr2; mouse, strain, sex
##### expr2.gene.annotation: row/probe annotation for expr2; probe, gene name
##### note reshape warning: "multiple rows match for mouse_number=339: first taken"
expr1 = reshape(expr0[,c(1,2,5)], idvar="probesetID", timevar="mouse_number", direction="wide")
expr.mouse.nums = sapply(names(expr1)[2:ncol(expr1)],function(x){strsplit(x,"expression_value.")[[1]][2]})
expr.mouse.strains = var.annot(vecdat=expr.mouse.nums, annot=annotation, invar="mouse_number", outvar="Strain")
expr.mouse.sex = var.annot(vecdat=expr.mouse.nums, annot=annotation, invar="mouse_number", outvar="Sex")
expr2 = expr1[,2:ncol(expr1)]
rownames(expr2) = expr1$probesetID
names(expr2) = expr.mouse.nums
expr2.samp.annotation = cbind(expr.mouse.nums, expr.mouse.strains, expr.mouse.sex) %>% as.data.frame(stringsAsFactors=F)
names(expr2.samp.annotation) = c("mouse_number","strain","sex")
expr.genes = var.annot(vecdat=rownames(expr2), annot=expr0, invar="probesetID", outvar="gene_symbol")
expr2.gene.annotation = cbind(rownames(expr2), expr.genes) %>% as.data.frame(stringsAsFactors=F)
names(expr2.gene.annotation) = c("probe","gene")

# set phenotype matrix and annotation
# main objects: 
##### pheno2: phenotype matrix
##### pheno2.samp.annotation: column/sample annotation for pheno2; mouse, strain, sex
##### note reshape warning: "multiple rows match for mouse_number=339: first taken"
pheno1 = reshape(pheno0[,c(1,2,5)], idvar="trait_name", timevar="mouse_number", direction="wide")
pheno.mouse.nums = sapply(names(pheno1)[2:ncol(pheno1)],function(x){strsplit(x,"value.")[[1]][2]})
pheno.mouse.strains = var.annot(vecdat=pheno.mouse.nums, annot=annotation, invar="mouse_number", outvar="Strain")
pheno.mouse.sex = var.annot(vecdat=pheno.mouse.nums, annot=annotation, invar="mouse_number", outvar="Sex")
pheno2 = pheno1[,2:ncol(pheno1)]
rownames(pheno2) = pheno1$trait_name
names(pheno2) = pheno.mouse.nums
pheno2.samp.annotation = cbind(pheno.mouse.nums, pheno.mouse.strains, pheno.mouse.sex) %>% as.data.frame(stringsAsFactors=F)
names(pheno2.samp.annotation) = c("mouse_number","strain","sex")


############################################################
## select strains for analysis
############################################################

# identify strains with both M and F - expression and phenotype
strainsM.expr = expr2.samp.annotation %>% filter(sex == "M") %>% select(strain) %>% unique
strainsM.pheno = pheno2.samp.annotation %>% filter(sex == "M") %>% select(strain) %>% unique
strainsF.expr = expr2.samp.annotation %>% filter(sex == "F") %>% select(strain) %>% unique
strainsF.pheno = pheno2.samp.annotation %>% filter(sex == "F") %>% select(strain) %>% unique
strainsMF = Reduce(intersect, list(strainsM.expr$strain,strainsM.pheno$strain,strainsF.expr$strain,strainsF.pheno$strain))

# filter expression data for intersect strains
# main objects: 
##### expr3: gene expression matrix
##### expr3.samp.annotation: column/sample annotation for expr3; mouse, strain, sex
##### expr2.gene.annotation: row/probe annotation for expr3; probe, gene name
expr.strain.inds = expr2.samp.annotation$strain %in% strainsMF
expr3.samp.annotation = expr2.samp.annotation[expr.strain.inds,]
expr3 = expr2[,expr.strain.inds]

# filter expression data for intersect strains
# main objects: 
##### pheno3: phenotype matrix
##### pheno3.samp.annotation: column/sample annotation for pheno3; mouse, strain, sex
pheno.strain.inds = pheno2.samp.annotation$strain %in% strainsMF
pheno3.samp.annotation = pheno2.samp.annotation[pheno.strain.inds,]
pheno3 = pheno2[,pheno.strain.inds]

# modify organization
# expr and pheno data are organized identically by mouse number
# main objects: 
##### pheno4: phenotype matrix
##### pheno4.samp.annotation: column/sample annotation for pheno4; mouse, strain, sex
expr.nums = expr3.samp.annotation$mouse_number
pheno.inds = sapply(expr.nums,function(x){which(pheno3.samp.annotation$mouse_number==x)}) %>% unlist
pheno4.samp.annotation = pheno3.samp.annotation[pheno.inds,]
pheno4 = pheno3[,pheno.inds]

# check data organization
# note that the annotation files are identical now
all(expr3.samp.annotation$mouse_number == pheno4.samp.annotation$mouse_number)
all(names(expr3) == names(pheno4))
all(names(pheno4) == expr3.samp.annotation$mouse_number)
samp.annotation = pheno4.samp.annotation

# save.image("procdata.RData")
# load("procdata.RData")


############################################################
## perform normality checks
## data are approximately normal and will not be further modified
############################################################

# exploratory analysis
pheno.mins = apply(pheno4,1,function(x){min(x[is.na(x)==F])})
hist(pheno4[142,] %>% data.matrix %>% as.numeric)
hist(pheno4[which(pheno.mins==min(pheno.mins)),] %>% data.matrix %>% as.numeric)
head(pheno4[which(pheno.mins==min(pheno.mins)),1:5])
hist(pheno4[which(rownames(pheno4)=="Total cholesterol (liver)"),] %>% data.matrix %>% as.numeric)
hist(pheno4[which(rownames(pheno4)=="log10(Total cholesterol (liver))"),] %>% data.matrix %>% as.numeric)
head(pheno4[75,1:5])


############################################################
## separate by sex
## take means of multiple samples in a given strain
############################################################

# separate data by sex
ind.m = which(samp.annotation$sex=="M")
ind.f = which(samp.annotation$sex=="F")
anno.M = samp.annotation[ind.m,] # annotation, male
anno.F = samp.annotation[ind.f,] # annotation, female
expr.M = expr3[,ind.m] # expression, male
expr.F = expr3[,ind.f] # expression, female
pheno.M = pheno4[,ind.m] # phenotype, male
pheno.F = pheno4[,ind.f] # phenotype, female

# function to get means and revised annotations
get.means = function(expr=NULL,anno=NULL,strains=NULL){
  mat = matrix(0, nrow(expr), length(strains)) %>% as.data.frame
  nreps = rep(0, length(strains))
  for(ii in 1:length(strains)){
    ind.str = which(anno$strain == strains[ii])
    nreps[ii] = length(ind.str)
    if(length(ind.str)>1){
      mat[,ii] = apply(expr[,ind.str],1,mean)
    } else {mat[,ii] = expr[,ind.str]}
  } # ii loop
  names(mat) = strains
  rownames(mat) = rownames(expr)
  ann.out = cbind(strains,nreps) %>% as.data.frame(stringsAsFactors=F)
  return(list(expr=mat, anno=ann.out))
}

# get the mean expression data and strain annotation for males and females
expr.M.mean = get.means(expr=expr.M, anno=anno.M, strains=unique(anno.M$strain))
expr.F.mean = get.means(expr=expr.F, anno=anno.F, strains=unique(anno.F$strain))

# get the mean phenotype data and strain annotation for males and females
pheno.M.mean = get.means(expr=pheno.M, anno=anno.M, strains=unique(anno.M$strain))
pheno.F.mean = get.means(expr=pheno.F, anno=anno.F, strains=unique(anno.F$strain))

# adjust organization for mean data
ind.f.new = sapply(names(expr.M.mean$expr),function(x){which(names(expr.F.mean$expr)==x)}) %>% unlist
expr.F.mean$expr = expr.F.mean$expr[,ind.f.new]
ind.f.new = sapply(names(pheno.M.mean$expr),function(x){which(names(pheno.F.mean$expr)==x)}) %>% unlist
pheno.F.mean$expr = pheno.F.mean$expr[,ind.f.new]

# check data organization
all(names(expr.F.mean$expr) == names(expr.M.mean$expr))
all(names(pheno.F.mean$expr) == names(pheno.M.mean$expr))
all(names(pheno.F.mean$expr) == names(expr.F.mean$expr))

# data for analysis
expr.mean.M = expr.M.mean$expr
expr.mean.F = expr.F.mean$expr
pheno.mean.M = pheno.M.mean$expr
pheno.mean.F = pheno.F.mean$expr


############################################################
## for each gene, take the average across probes
############################################################

# get probe averages for each gene
get.probe.avg = function(dat=NULL, ann=NULL){
  pb = data.frame(probe=rownames(dat),stringsAsFactors=F)
  ann0 = merge(pb,ann,by="probe")
  genes = unique(ann0$gene)
  new = data.frame(matrix(0,nrow=length(genes),ncol=ncol(dat)))
  rownames(new) = genes
  names(new) = names(dat)
  for(ii in 1:nrow(new)){
    gn = genes[ii]
    prb = ann0$probe[which(ann0$gene == gn)]
    datprb = dat[pb$probe %in% prb,]
    new[ii,] = apply(datprb,2,mean)
  }
  return(new)
} # get.probe.avg
exprM = get.probe.avg(dat=expr.mean.M, ann=expr2.gene.annotation)
exprF = get.probe.avg(dat=expr.mean.F, ann=expr2.gene.annotation)


############################################################
## Select phenotypes of interest based on assessment 
## of the phenotype clusters
############################################################

# NMR_BF%_8wks
ind.bf.m = which(rownames(pheno.mean.M) == "NMR_BF%_8wks")
ind.bf.f = which(rownames(pheno.mean.F)=="NMR_BF%_8wks")
fat8.M = pheno.mean.M[ind.bf.m,]
fat8.F = pheno.mean.F[ind.bf.f,]

# HOMA-IR
ind.bf.m = which(rownames(pheno.mean.M) == "HOMA-IR")
ind.bf.f = which(rownames(pheno.mean.F)=="HOMA-IR")
homair.M = pheno.mean.M[ind.bf.m,]
homair.F = pheno.mean.F[ind.bf.f,]

# BF_Percent_Growth_0to8wks
ind.bf.m = which(rownames(pheno.mean.M) == "BF_Percent_Growth_0to8wks")
ind.bf.f = which(rownames(pheno.mean.F)=="BF_Percent_Growth_0to8wks")
growth.M = pheno.mean.M[ind.bf.m,]
growth.F = pheno.mean.F[ind.bf.f,]


############################################################
## implement correlation analysis
############################################################

# cor.fn: function to compute the correlations of interest
# correlations with respect to a single phenotype
cor.fn = function(de=NULL, dp=NULL, thresh=0.05, na.prop=0.4){
  
  # data storage
  res = data.frame(matrix(NA,nrow(de),4))
  
  # x is the the current row/gene of the dExpr  
  # y is the the phenotype 
  y = dp %>% data.matrix %>% as.numeric
  pheno = rownames(dp)
  for(ii in 1:nrow(de)){

    # specify expression vector
    x = de[ii,] %>% data.matrix %>% as.numeric
    
    # if there is a certain proportion of missing values, skip to next row
    if( length(which(is.na(x)==TRUE)) > na.prop*ncol(de) | length(which(is.na(y)==TRUE)) > na.prop*ncol(de) ){next}
    
    # compute the correlation
    rr = bicorAndPvalue(x,y,use="pairwise.complete.obs")
    if( is.na(rr$bicor) == TRUE ){next}
    
    # record the data
    res[ii,] = c(rownames(de)[ii], pheno, rr$bicor, rr$p)

  } # ii

  # after it has gone through all of the rows, converts the results to a data frame
  res[,3:4] = apply(res[,3:4],2,function(x){data.matrix(x)%>%as.numeric})
  names(res) = c("gene","pheno","bicor","p-value")
  na.ind = which(is.na(res$gene)==TRUE)
  if(length(na.ind)>0){res = res[-na.ind,]}
  res = res %>% mutate(fdr = p.adjust(res$`p-value`,method="BH"))
  
  # return the results
  return(res)
  
} # cor.fn

# identify correlations for body fat
de = exprF - exprM
dp = fat8.F - fat8.M
fat.res = cor.fn(de=de, dp=dp, thresh=0.05, na.prop=0.4)

# identify correlations for insulin resistance
de = exprF - exprM
dp = homair.F - homair.M
insulin.res = cor.fn(de=de, dp=dp, thresh=0.05, na.prop=0.4)

# identify correlations for percent growth
de = exprF - exprM
dp = growth.F - growth.M
growth.res = cor.fn(de=de, dp=dp, thresh=0.05, na.prop=0.4)


############################################################
## implement differential expression analysis
############################################################

# specify data
null <- "Female"
alt <- "Male"
expr.F.Null = t(exprF) 
expr.M.Alt = t(exprM)

# function to calculate differential expression
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
}

# evaluate differential expression
degs0 = deg.analysis(dat_Null=expr.F.Null, dat_Alt=expr.M.Alt, null=null, alt=alt)

# check
ind = which(colnames(dat_Null) == "Ddx3y")
fav = mean(dat_Null[,ind])
mav = mean(dat_Alt[,ind])
fc = fav - mav

# output the results
output = degs0 %>% mutate(gene=rownames(degs0))
fname = "hmdpDEG.txt"
write.table(output,fname,col.names=T,row.names=F,sep="\t",quote=F)

# save.image("cordata.RData")
# load("cordata.RData")

############################################################
## filter deg and cor results, assess overlap (Venn)
############################################################

# parameters
fdr.thresh = 0.05
fc.thresh = log2(1.05)
cor.thresh = 0.2

# filter the correlation data
fat.res.filt = fat.res %>% filter(abs(bicor)>cor.thresh & fdr<fdr.thresh) # 3602 genes
insulin.res.filt = insulin.res %>% filter(abs(bicor)>cor.thresh & fdr<fdr.thresh) # 1081 genes
growth.res.filt = growth.res %>% filter(abs(bicor)>cor.thresh & fdr<fdr.thresh) # 894 genes

# filter the deg data
degs0 = degs0 %>% mutate(gene = rownames(degs0))
deg.filt = degs0 %>% filter(abs(logFC)>fc.thresh & adj.P.Val<fdr.thresh) # 5352 genes

# generate information for supplementary table
ts13.fat = merge(fat.res,degs0[,c(1,4,5,10)],by="gene")
ts13.grw = merge(growth.res,degs0[,c(1,4,5,10)],by="gene")
ts13.ins = merge(insulin.res,degs0[,c(1,4,5,10)],by="gene")
tables13 = rbind(ts13.fat,ts13.grw,ts13.ins)
write.table(tables13,"TableS13.txt",col.names=T,row.names=F,quote=F,sep="\t")

###############################
# intersection venn plots

inter = Reduce(intersect, list(fat.res.filt$gene,growth.res.filt$gene,insulin.res.filt$gene)) %>% length
fat.grow = intersect(fat.res.filt$gene,growth.res.filt$gene) %>% length
fat.ir = intersect(fat.res.filt$gene,insulin.res.filt$gene) %>% length
grow.ir = intersect(growth.res.filt$gene,insulin.res.filt$gene) %>% length

pdf("degcorVennFig6.pdf")

# 3 pheno
venn.plot <- draw.triple.venn(area1 = nrow(fat.res.filt), 
                              area2 = nrow(growth.res.filt), 
                              area3 = nrow(insulin.res.filt),
                              n12 = fat.grow, 
                              n23 = grow.ir, 
                              n13 = fat.ir,
                              n123 = inter, 
                              category = c("fat", "grow", "ir"),
                              euler.d=F, scaled=F, overrideTriple=T)

# fat intersection
a1 = nrow(fat.res.filt)
a2 = nrow(deg.filt)
a3 = length(intersect(fat.res.filt$gene,deg.filt$gene))
grid.newpage()
draw.pairwise.venn(area1=a1, area2=a2, cross.area=a3, category=c("fat","deg"))

# insulin intersection
a1 = nrow(insulin.res.filt)
a2 = nrow(deg.filt)
a3 = length(intersect(insulin.res.filt$gene,deg.filt$gene))
grid.newpage()
draw.pairwise.venn(area1=a1, area2=a2, cross.area=a3, category=c("ins","deg"))

# growth intersection
a1 = nrow(growth.res.filt)
a2 = nrow(deg.filt)
a3 = length(intersect(growth.res.filt$gene,deg.filt$gene))
grid.newpage()
draw.pairwise.venn(area1=a1, area2=a2, cross.area=a3, category=c("growth","deg"))

dev.off()


############################################################
## generate plots for the fold change correlation analysis
############################################################

# functions to identify an optimal plotting example
make.barplot = function(pltdat=NULL,main=NULL){
  indpos = which(pltdat > 0)
  indneg = which(pltdat < 0)
  col = rep("red",length(c(indpos,indneg)))
  col[indneg] = "blue"
  barplot(pltdat, col=col, xlab="Strains", ylab="Fold-Change", 
          main=main, cex.axis = 1.5, cex.lab = 1.5) 
} # make.barplot
fcg.plot.gen = function(pheno=NULL,degF=NULL,degM=NULL,fname=NULL,
                        class=NULL,fc.gene=NULL,dg.gene=NULL){
  
  pdf(fname)
  par(mfrow=c(2,2))
  
  # order the phenotype data
  pheno = pheno[order(pheno,decreasing=TRUE)]
  pheno = pheno[which(is.na(pheno)==FALSE)]
  strains = names(pheno)
  
  # bar graph of the phenotype for every strain
  make.barplot(as.numeric(pheno[1,]), main=class)
  
  # deg bar graph
  ddeg = degF - degM
  deg.strain = rownames(degM)
  ind.str = sapply(strains,function(x){which(deg.strain==x)})
  ind.gene = which(colnames(degM) == dg.gene)
  make.barplot(as.numeric(ddeg[ind.str,ind.gene]), main=dg.gene)
  
  # fc cor bar graph
  ind.gene = which(colnames(degM) == fc.gene)
  make.barplot(as.numeric(ddeg[ind.str,ind.gene]), main=fc.gene)
  
  # deg scatter
  x = as.numeric(ddeg[ind.str,which(colnames(degM) == dg.gene)])
  y = as.numeric(pheno[1,])
  plot(x,y,pch=19,ylab=class,xlab=dg.gene)
  abline(lm(y ~ x))
  
  # fc scatter
  x = as.numeric(ddeg[ind.str,which(colnames(degM) == fc.gene)])
  y = as.numeric(pheno[1,])
  plot(x,y,pch=19,ylab=class,xlab=fc.gene)
  abline(lm(y ~ x))

  dev.off()
  
} # fcg.plot.gen

# select genes and plot for fat
deg = degs0[,c(1,4,5,10)] %>% mutate(sig = -log10(adj.P.Val))
fccor = fat.res %>% mutate(sig = -log10(fdr))
deg.fccor = merge(deg,fccor,by="gene")
deg.fccor = deg.fccor %>% mutate(ratio = sig.x/sig.y)
deg.fccor = deg.fccor[order(deg.fccor$ratio),] 
fc.gene = "Xbp1" 
dg.gene = "Uty"
deg.fccor[deg.fccor$gene %in% fc.gene,]
deg.fccor[deg.fccor$gene %in% dg.gene,]
pheno = fat8.F - fat8.M
class = "fat"
degF = expr.F.Null
degM = expr.M.Alt
fcg.plot.gen(pheno,degF,degM,fname="fatFig6.pdf",
            class,fc.gene,dg.gene)

# select genes and plot for insulin
deg = degs0[,c(1,4,5,10)] %>% mutate(sig = -log10(adj.P.Val))
fccor = insulin.res %>% mutate(sig = -log10(fdr))
deg.fccor = merge(deg,fccor,by="gene")
deg.fccor = deg.fccor %>% mutate(ratio = sig.x/sig.y)
deg.fccor = deg.fccor[order(deg.fccor$ratio),]
fc.gene = "Rgs3"
dg.gene = "Fabp3"
deg.fccor[deg.fccor$gene %in% fc.gene,]
deg.fccor[deg.fccor$gene %in% dg.gene,]
pheno = homair.F - homair.M
class = "HOMA-IR"
degF = expr.F.Null
degM = expr.M.Alt
fcg.plot.gen(pheno,degF,degM,fname="irFig6.pdf",
             class,fc.gene,dg.gene)

# select genes and plot for growth
deg = degs0[,c(1,4,5,10)] %>% mutate(sig = -log10(adj.P.Val))
fccor = growth.res %>% mutate(sig = -log10(fdr))
deg.fccor = merge(deg,fccor,by="gene")
deg.fccor = deg.fccor %>% mutate(ratio = sig.x/sig.y)
deg.fccor = deg.fccor[order(deg.fccor$ratio),]
fc.gene = "Plekhg6" 
dg.gene = "Hoxa10" 
deg.fccor[deg.fccor$gene %in% fc.gene,]
deg.fccor[deg.fccor$gene %in% dg.gene,]
pheno = growth.F - growth.M
class = "growth"
degF = expr.F.Null
degM = expr.M.Alt
fcg.plot.gen(pheno,degF,degM,fname="growFig6.pdf",
             class,fc.gene,dg.gene)


de = exprF - exprM
dp = fat8.F - fat8.M
fat.res = cor.fn(de=de, dp=dp, thresh=0.05, na.prop=0.4)

# identify correlations for insulin resistance
de = exprF - exprM
dp = homair.F - homair.M
insulin.res = cor.fn(de=de, dp=dp, thresh=0.05, na.prop=0.4)

# identify correlations for percent growth
de = exprF - exprM
dp = growth.F - growth.M
growth.res = cor.fn(de=de, dp=dp, thresh=0.05, na.prop=0.4)


############################################################
## overlap fc genes with interaction eqtls
## for fig 5
############################################################

library(dplyr)
library(biomaRt)

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig6"
setwd(dir)

# union eqtl data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4/union.eqtl.RData"
load(fname)

# combine data for fc cor analysis
all(fat.res$gene == insulin.res$gene)
all(fat.res$gene == growth.res$gene)
all(insulin.res$gene == growth.res$gene)
fc.cor = cbind(fat.res[,c(1,3:5)], insulin.res[,3:5], growth.res[,3:5])
names(fc.cor)[2:4] = c("bicor.fat", "pval.fat", "fdr.fat")
names(fc.cor)[5:7] = c("bicor.ir", "pval.ir", "fdr.ir")
names(fc.cor)[8:10] = c("bicor.growth", "pval.growth", "fdr.growth")

# filter for significant fc genes
fc.genes = Reduce(union, list(fat.res.filt$gene,insulin.res.filt$gene,growth.res.filt$gene))
fc.cor.sig = fc.cor[fc.cor$gene %in% fc.genes,]
length(unique(fc.cor.sig$gene))

## get human orthologs of mouse genes
## get overlapping lists of genes

union.eqtl = union.eqtl[order(union.eqtl$p_peer),]
union.eqtl = union.eqtl[!duplicated(union.eqtl$gene),]
eqtl = union.eqtl[,-c(2,3)]
names(eqtl)[8] = "gene"

genesH = eqtl$gene
genesM = unique(fc.cor.sig$gene)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesHM = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                 values = genesH, mart = mouse, attributesL = c("hgnc_symbol"), 
                 martL = human, uniqueRows=T)
orthologs = sapply(genesH,function(x){genesHM$MGI.symbol[which(genesHM$HGNC.symbol==x)]}) %>% unlist


# filter orthologs to include only human genes with murine orthologs
orthos = c()
for(ii in 1:length(orthologs)){
  orth = names(orthologs)[ii]
  orths = orthologs[[ii]]
  if(length(orths)>0){
    orths = paste(orths,collapse=",")
    new = c(orth, orths)
    orthos = rbind(orthos,new)
  }
}
orthos = data.frame(human=orthos[,1],mouse=orthos[,2],stringsAsFactors=F)
dim(orthos)
length(unique(orthos$human))
length(unique(orthos$mouse))

# merge fc cor data with interaction eqtl data
# save(eqtl.fccor, file="eqtl.fccor.RData")

names(fc.cor.sig)[1] = "mouse"
fc.cor.ortho = merge(fc.cor.sig,orthos,by="mouse")
names(fc.cor.ortho)[11] = "gene"
names(eqtl)[8] = "gene"
eqtl.fccor = merge(eqtl,fc.cor.ortho,by="gene")

# functional enrichment of eqtl.fccor genes
library(enrichR)
100 * nrow(eqtl.fccor) / nrow(eqtl) 
100 * nrow(eqtl.fccor) / nrow(fc.cor.sig) 

# select databases
db1 = c("Reactome_2016","KEGG_2019_Human","WikiPathways_2019_Human","Panther_2016")
db2 = c("GO_Biological_Process_2018","GO_Cellular_Component_2018","GO_Molecular_Function_2018")
db = c(db1,db2)

# get enrichment results
res = enrichr(genes=eqtl.fccor$gene, databases=db)

# filter each set of results
fdrthresh = 0.05
pthresh = 0.001
ngenes = 1
res_filtered = c()
for(ii in 1:length(res)){
  dat = res[[ii]]
  dat = dat %>% filter(Adjusted.P.value < fdrthresh)
  #dat = dat %>% filter(P.value < pthresh)
  dat = dat %>% mutate(source = names(res)[ii])
  nge = dat$Genes
  nge = sapply(nge,function(x){
    ge = strsplit(x,";")[[1]]
    return(length(ge))
  })
  dat = dat %>% mutate(nge=nge) 
  dat = dat %>% filter(nge > ngenes)
  res_filtered = rbind(res_filtered, dat)
}
res_filtered = res_filtered[order(res_filtered$Odds.Ratio,decreasing=T),]
res_filtered = res_filtered[,c(1,3,4,7,9,10,11)]
save(res_filtered, file="res_filtered.RData")

# table S14
out = res_filtered
write.table(out,"TableS14.txt",col.names=T,row.names=F,quote=F,sep="\t")

