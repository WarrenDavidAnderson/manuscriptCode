
############################################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Figure 7
############################################################################

library(dplyr)
library(ggplot2)
library(gridExtra)
library(Biobase)
library(GEOquery)
library(limma)
library(lmtest)

setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig7")


############################################################################
# download and annotate SGBS timeseries data
############################################################################

# load series and platform data from GEO
gset <- getGEO("GSE76131", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# log2 transform (not applied)
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { 
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) 
  }
expr_data = exprs(gset)

# get gene list
gene_info = fData(gset) %>% select(ID, Gene.symbol)
ind_missing = which(gene_info$Gene.symbol == "")
gene_info = apply(gene_info,2,as.character) %>% as.data.frame(stringsAsFactors=FALSE)
gene_info$Gene.symbol[ind_missing] = gene_info$ID[ind_missing]

# sample annotation
pdat0 = pData(gset) %>% select(title,source_name_ch1,organism_ch1,characteristics_ch1.1)
names(pdat0) = c("sample","cell","org","time")
unique(pdat0$cell)
unique(pdat0$org)
unique(pdat0$time)
pdat1 = cbind(rownames(pdat0), pdat0)
names(pdat1)[1] = "sample_id"
pdat1[] = apply(pdat1,2,as.character) 

# check organization of gene and sample IDs
all.equal(colnames(expr_data), pdat1$sample_id)
all.equal(rownames(expr_data), gene_info$ID)

# revise sample annotaton
pdat1$time = sapply(pdat1$time,function(x){
  out1 = strsplit(x,": ")[[1]][2]
  return(out1)
})

# SGBS times
unique(pdat1$time)

# average probes for individual genes
avg.probe = function(expr=NULL,gene=NULL){
  gene = gene[gene$ID %in% rownames(expr),]
  out = matrix(0,length(unique(gene$Gene.symbol)), ncol(expr))
  rownames(out) = unique(gene$Gene.symbol)
  colnames(out) = colnames(expr)
  for(ii in 1:nrow(out)){
    gg = unique(gene$Gene.symbol)[ii]
    pb = gene$ID[which(gene$Gene.symbol==gg)]
    dat = expr[rownames(expr) %in% pb,]
    if(class(dat) == "numeric"){
      out[ii,] = dat
    } else {
      out[ii,] = apply(dat,2,mean)
    }
  }
  out = data.frame(out)
}
sgbs = avg.probe(expr=expr_data, gene=gene_info)

# save(sgbs, file="sgbs.RData")

###############################################################
###### check SGBS expression quantiles and timeseries dymanics
## https://www.ncbi.nlm.nih.gov/pubmed/27385551
###############################################################

library(lmtest)

setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig7")
load("sgbs.RData")

# quantiles for differentiated samples
ind = which(pdat1$time=="384h")
av.adip = apply(sgbs[,ind],1,mean)
quantile(av.adip)  

# function to compute a specific quantile value
get.quant = function(dat=NULL, gene=NULL){
  percentile <- ecdf(dat)
  gg = dat[which(names(dat) == gene)]
  out = percentile(gg)
  return(out)
}

# print the percentiles
genes = c("CCDC3","CLIC6","FADS1","GLDN","HSPA12A","MAP1B","MLPH","MMD","MYOT","NDRG4","NEO1","PDZD2","TBC1D9")
quants = sapply(genes,function(x){get.quant(dat=av.adip, gene=x)})

# function for evaluating timeseries dynamics using linear model LRT
tdyn.fun = function(dat=NULL, times=NULL, gene=NULL){
  gg = dat[which(rownames(dat) == gene),] %>% as.numeric
  df = data.frame(expr=gg, time=times)
  hA = lm(expr ~ time, data=df)
  h0 = lm(expr ~ 1, data=df)
  lrt = lrtest(h0, hA)
  return(lrt$`Pr(>Chisq)`[2])
}

# compute LRTs
times = sapply(pdat1$time,function(x){strsplit(x,"h")[[1]][1] %>% as.numeric})
LRTs = sapply(genes,function(x){
  tdyn.fun(dat=sgbs, times=pdat1$time, gene=x)
})
fdr = p.adjust(LRTs,method="BH")

# timeseries plot function
time.plt = function(dat=NULL, times=NULL, genes=NULL, span=0.5){
  plts = list()
  for(gene in genes){
    gg = dat[which(rownames(dat) == gene),] %>% as.numeric
    df = data.frame(expr=gg, time=c(times/24))
    p1 = ggplot(df, aes(x=time, y=expr)) + 
      geom_point() + 
      geom_smooth(se=T, method="loess", colour="black", span=span) +
      ylab(paste0(gene, " expression")) + 
      xlab("adipogenesis time (days)")
    plts[[gene]] = p1
  }
  return(plts)
}

# get the timeseries plots
times = sapply(pdat1$time,function(x){strsplit(x,"h")[[1]][1] %>% as.numeric})
plts = time.plt(dat=sgbs, times=times, genes=genes)
library(gridExtra)
pdf("SGBSadipogen2.pdf", onefile = TRUE, height=11, width=8.5)
marrangeGrob(grobs=plts, nrow=6, ncol=4, top=NULL)
dev.off()

# save.image("SGBS.RData")
# load("SGBS.RData")


###############################################################
###### evaluate human adipocyte ASC data
## [1] tissue: abdominal subcutaneous adipose tissue
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77532
###############################################################

setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig7")

library(Biobase)
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)

# download and manipulate gene identifiers
grep -n -w "platform_table_begin" GSE77532_family.soft # 114
tail -n +115 GSE77532_family.soft | head -n -1 > anno0.txt
anno0 = read.table("anno0.txt",fill=T,header=T,stringsAsFactors=F)

# process gene ids
ind.rm = grep("chr",anno0$GB_ACC)
anno1 = anno0[-ind.rm,]
ind.rm = which(anno1$GB_ACC=="")
anno1 = anno1[-ind.rm,]

# get gene ids biomart (86-91)
# save(anno.adip, file="anno.adip.RData")
# load("anno.adip.RData")
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
refseq <- anno1$GB_ACC
ensids0 = getBM(filters="refseq_mrna", attributes=c("refseq_mrna","external_gene_name"), values=refseq, mart=mart)
names(anno1)[6] = "refseq_mrna" 
anno.adip = merge(anno1, ensids0, by="refseq_mrna")

# download expression data
gset <- getGEO("GSE77532", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL16686", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# log2 transform - not applied
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { 
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) 
}

# expression data set
expr_data = data.frame(exprs(gset))

# gene annotation, in order of expression data frame
# remove non-annotated probes
# load("anno.adip.RData")
names(anno.adip)[2] = "probe"
anno.adip = anno.adip[anno.adip$probe %in% rownames(expr_data),]
anno.adip = anno.adip[!duplicated(anno.adip$probe),]
expr_data = expr_data[rownames(expr_data) %in% anno.adip$probe,]
ydat = data.frame(ind=c(1:nrow(expr_data)), probe=rownames(expr_data))
gene_info = merge(anno.adip[,c(2,9)], ydat, by="probe")
gene_info = gene_info[order(gene_info$ind,decreasing=FALSE),]
all(gene_info$probe == rownames(expr_data))

# sample annotation, (rnaseq = GPL16791)
pdat0 = pData(gset) %>% dplyr::select(geo_accession,
                               `gender:ch1`,
                               `culture time:ch1`,
                               `induction:ch1`,
                               `age (years):ch1`)
names(pdat0) = c("id","sex","day","induction","age")

# explore and organize phenotype data
# same organization as expression object
pdat0$day = sapply(pdat0$day,function(x){
  a = as.character(x)
  strsplit(a,"day ")[[1]][2]
})
ind = sapply(names(expr_data),function(x){which(pdat0$id==x)})
pdat = pdat0[ind,]
all(pdat$id == names(expr_data))

# filter and organize pheno/expr data
# save.image("adipocyte.RData")
ind1 = which(pdat$induction == "post-adipogenic induction")
ind2 = which(pdat$day == "1" & pdat$induction == "control")
pdat = pdat[c(ind1,ind2),]
pdat$day[which(pdat$day == "1" & pdat$induction == "control")] = 0
expr_data = expr_data[,c(ind1,ind2)]

# plot dynamics
ind = which(gene_info$external_gene_name == "FADS1")
y = expr_data[ind,] %>% as.numeric
plot(y)

###############################################################
###### evaluate human adipocyte expression data
## Statistical analysis and plotting
###############################################################

library(dplyr)
library(ggplot2)
library(lmtest)

# gene list
genes = c("CCDC3","CLIC6","FADS1","GLDN","HSPA12A","MAP1B","MLPH","MMD","MYOT","NDRG4","NEO1","PDZD2","TBC1D9")

# function for evaluating timeseries dynamics using linear model LRT
tdyn.fun = function(dat=NULL, times=NULL, gene=NULL, gene_info=NULL){
  ind = which(gene_info$external_gene_name == gene)
  gg = dat[ind,] 
  gg = apply(gg,2,mean) %>% as.numeric
  if(all(is.na(gg) == TRUE)){return(NA);break}
  df = data.frame(expr=gg, time=times)
  hA = lm(expr ~ time, data=df)
  h0 = lm(expr ~ 1, data=df)
  lrt = lrtest(h0, hA)
  out = lrt$`Pr(>Chisq)`[2]
  return(out)
}

# compute LRTs
times = pdat$day
LRTs = sapply(genes,function(x){
  tdyn.fun(dat=expr_data, times=times, gene=x, gene_info=gene_info)
})
LRT.fdr = p.adjust(LRTs,method="BH")

# timeseries plot function
time.plt = function(dat=NULL, times=NULL, genes=NULL, gene_info=NULL, span=0.5){
  plts = list()
  for(gene in genes){
    ind = which(gene_info$external_gene_name == gene)
    gg = dat[ind,] 
    gg = apply(gg,2,mean) %>% as.numeric
    if(all(is.na(gg) == TRUE)){next}
    df = data.frame(expr=gg, time=times)
    p1 = ggplot(df, aes(x=time, y=expr)) + 
      geom_point() + 
      geom_smooth(se=T, method="loess", colour="black", span=span) +
      ylab(paste0(gene, " expression")) + 
      xlab("adipogenesis time (days)")
    plts[[gene]] = p1
  }
  return(plts)
}

# get the timeseries plots
times = as.numeric(pdat$day)
plts = time.plt(dat=expr_data, times=times, genes=genes, gene_info=gene_info)
library(gridExtra)
pdf("hADIPadipogen2.pdf", onefile = TRUE, height=11, width=8.5)
marrangeGrob(grobs=plts, nrow=6, ncol=4, top=NULL)
dev.off()

# quantiles for differentiated samples
ind = which(times==21)
av.adip = apply(expr_data[,ind],1,mean)
quantile(av.adip)  

# function to compute a specific quantile value
get.quant = function(dat=NULL, gene=NULL, gene_info=NULL){
  percentile <- ecdf(dat)
  ind = which(gene_info$external_gene_name == gene)
  gg = dat[ind] %>% as.numeric
  out = percentile(gg)
  return(out)
}

# print the percentiles
genes = c("CCDC3","CLIC6","FADS1","GLDN","HSPA12A","MAP1B","MLPH","MMD","MYOT","NDRG4","NEO1","PDZD2","TBC1D9")
quants = sapply(genes,function(x){get.quant(dat=av.adip, gene=x, gene_info=gene_info)}) %>% unlist

# save.image("HADIP.RData")
# load("HADIP.RData")


###############################################################
###### check 3T3 expression quantiles and timeseries dymanics
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
###############################################################

setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig7")
setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig7/GSE95533")

library(Biobase)
library(GEOquery)
library(limma)
library(dplyr)
library(DESeq2)
library(ggplot2)

# get data
cd /media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig7/GSE95533
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE95nnn/GSE95533/matrix/GSE95533_series_matrix.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE95nnn/GSE95533/suppl/GSE95533%5FRNA%5Fcounts%2Etxt%2Egz

# import read count data
pth = "GSE95533_RNA_counts.txt"
expr0 = read.table(pth,header=T,stringsAsFactors=F)

# convert refseq ids into gene names
library(biomaRt)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
refseq <- expr0$RefSeq
geneids0 = getBM(filters="refseq_mrna", attributes=c("refseq_mrna","external_gene_name"), values=refseq, mart=mart)
names(geneids0) = c("RefSeq","gene")
expr1 = merge(geneids0, expr0, by="RefSeq")

# get orthologs
genes = c("CCDC3","CLIC6","FADS1","GLDN","HSPA12A","MAP1B","MLPH","MMD","MYOT","NDRG4","NEO1","PDZD2","TBC1D9")
genesH = genes
genesM = geneids0$gene
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesHM = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                 values = genesH, mart = mouse, attributesL = c("hgnc_symbol"), 
                 martL = human, uniqueRows=T)
orthologs = sapply(genesH,function(x){genesHM$MGI.symbol[which(genesHM$HGNC.symbol==x)]})

# gene annotation file
gene.ann = expr1[,1:6]

# sample annotation file
samp.ann = names(expr1)[7:16]
rep = c(1,1,1,1,1,2,2,2,2,2)
day = c(0,4,24,48,168, 0,4,24,48,168) / 24
samp.ann = data.frame(id=samp.ann, rep=rep, day=day, stringsAsFactors=F)

# un-normalized expression matrix
# remove small number of duplicate gene entries
expr.unnorm = expr1[,7:16]
ind.rem = which(duplicated(gene.ann$gene)==T)
expr.unnorm = expr.unnorm[-ind.rem,]
gene.ann = gene.ann[-ind.rem,]
rownames(expr.unnorm) = gene.ann$gene

#####################
## normalize raw counts to size factors, generate DESeq object
#####################

# estimate size factors
size_factors = estimateSizeFactorsForMatrix(expr.unnorm)

## generate DESeqDataSet 
conditions = factor(samp.ann$day,levels=c(0,4,24,48,168) / 24)
colData = as.data.frame(conditions)
deseq_obj = DESeqDataSetFromMatrix(countData=expr.unnorm, 
                                           colData=colData, design=~conditions)
sizeFactors(deseq_obj) = size_factors

# note that each column of the count matrix was divided 
# by the respective size factor
head(sapply(c(1:10),function(x){expr.unnorm[,x]/size_factors[x]})[,1:5])
head(counts(deseq_obj, normalized=TRUE)[,1:5])

#####################
## implement differential expression analysis - LRT
#####################

# basic analysis
deg = DESeq(deseq_obj, test="LRT", full=~conditions,  reduced=~1)
res = results(deg)

# look at data of interest
data.frame(res[rownames(res) %in% orthologs,])

# timeseries plot function
time.plt = function(dat=NULL, times=NULL, genes=NULL, span=0.5){
  plts = list()
  for(gene in genes){
    ind = which(rownames(dat) == gene)
    if(length(ind)==0){next}
    gg = dat[ind,] %>% as.numeric
    df = data.frame(expr=gg, time=times)
    p1 = ggplot(df, aes(x=time, y=expr)) + 
      geom_point() + 
      geom_smooth(se=T, method="loess", colour="black", span=span) +
      ylab(paste0(gene, " expression")) + 
      xlab("adipogenesis time (days)")
    plts[[gene]] = p1
  }
  return(plts)
}

# plot the data - log2(norm counts + 1)
normdat = counts(deseq_obj, normalized=TRUE) + 1 %>% log2
plts = time.plt(dat=normdat, genes=orthologs, times=samp.ann$day)
library(gridExtra)
pdf("3T3adipogen2.pdf", onefile = TRUE, height=11, width=8.5)
marrangeGrob(grobs=plts, nrow=6, ncol=4, top=NULL)
dev.off()

# quantiles for differentiated samples
times = samp.ann$day
ind = which(times==7)
av.adip = apply(normdat[,ind],1,mean)
quantile(av.adip)  

# function to compute a specific quantile value
get.quant = function(dat=NULL, gene=NULL){
  percentile <- ecdf(dat)
  ind = which(names(dat) == gene)
  gg = dat[ind] %>% as.numeric
  out = percentile(gg)
  return(out)
}

# print the percentiles
quants = sapply(orthologs,function(x){get.quant(dat=av.adip, gene=x)}) %>% unlist

# save.image("HADIP.RData")
# load("HADIP.RData")