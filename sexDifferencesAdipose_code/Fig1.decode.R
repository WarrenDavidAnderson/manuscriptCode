
############################################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Figure 1
# acquire and process deCODE data
############################################################################

setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/get_geo_data")

# load relevant libraries
library(Biobase)
library(GEOquery)
library(dplyr)

# load series and platform data from GEO
gset <- getGEO("GSE7965", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL3991", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))
gene_info = fData(gset) %>% select(ID, Gene.symbol)

# log2 transform if necessary
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) {
  print("log2 transform performed")
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex)
} else {print("log2 transform not performed")}
expr_data0 = exprs(gset)

# sample annotation
pdat2 = phenoData(gset)
datainfo0 = pdat2@data[,1:2]
tissueType0 = strsplit(as.character(datainfo0$title)," ")
tissueType1 = unlist(lapply(tissueType0,function(x)return(x[1])))
datainfo1 = datainfo0
datainfo1$title = tissueType1

# tissue types
unique(tissueType1)

# get gene ID list
gene_list = sapply(rownames(expr_data0),function(x){
  ind0 = which(gene_info$ID==x)
  gene_0 = gene_info$Gene.symbol[ind0]
  if(gene_0==''){out=x}
  if(gene_0!=''){out=gene_0}
  return(as.character(out))
}) %>% unname

# add gene names to the expression data frame
geneList = cbind(rownames(expr_data0),gene_list)
colnames(geneList)[1] = "ID"
nrow(expr_data0) == nrow(geneList)

# check data organization (all TRUE)
all.equal(rownames(expr_data0), gene_info$ID)
all.equal(geneList[,1], rownames(expr_data0))
all.equal(colnames(expr_data0), rownames(datainfo0))

# isolate subcutaneous adipose data
adipose_samples = rownames(datainfo1)[which(datainfo1$title=="Adipose")]
adipose_ind = getIndices(adipose_samples, colnames(expr_data0))
length(adipose_ind) == length(adipose_samples)
expr_data1 = expr_data0[,adipose_ind]

# isolate XIST data
indXIST = which(geneList[,2]=="XIST")
dat_xist = expr_data1[indXIST,]

# plot XIST data with reasonable cutoffs demarkating the thresholds for attribution
# of sex identity
pdf("decode_XIST.pdf")
par(mar=c(6,6,6,6))
hist(dat_xist,main="XIST distribution",xlab="XIST expression",
     cex.lab=2,cex.main=2,cex.axis=1.5)
cutHi = -0.2
cutLo = -0.75
abline(v=cutHi,col="red")
abline(v=cutLo,col="red")
dev.off()

# set sex indices
indF = which(expr_data1[indXIST,] > cutHi)
indM = which(expr_data1[indXIST,] < cutLo)

# set sex annotation
MF_ann = matrix(0,ncol(expr_data1),2) %>% as.data.frame
names(MF_ann) = c("sample","sex")
MF_ann$sample = colnames(expr_data1)
MF_ann$sex[indF] = "Female"
MF_ann$sex[indM] = "Male"

# double check dimensions and data organization
# note that sample/gene organization in the expression data set
# matches the organization in the annotation frames
dim(expr_data1)
dim(MF_ann)
dim(geneList)
all( names(expr_data1) == MF_ann$sample )
all( rownames(expr_data1) == geneList[,1])

# write data
write.table(expr_data1,"expr_data_decode.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(MF_ann,"phenotypes_decode.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(geneList,"genes_decode.txt",col.names=T,row.names=F,quote=F,sep="\t")


