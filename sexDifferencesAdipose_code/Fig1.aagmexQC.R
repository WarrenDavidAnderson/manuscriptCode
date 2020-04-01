
############################################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Figure 1
# quality control of the AAGMEx data
############################################################################

setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/get_geo_data")

# load relevant libraries
library(dplyr)
library(reshape2)
library(beeswarm)

# read in data previously imported from GEO
expr_data = read.table("expr_data_aagmex.txt",header=T,sep="\t",stringsAsFactors=F)
pheno_data = read.table("phenotypes_aagmex.txt",header=T,sep="\t",stringsAsFactors=F)

# correlation matrix
cor_expr = cor(expr_data)
colnames(cor_expr) = colnames(expr_data)
rownames(cor_expr) = colnames(expr_data)

# isolate vector of correlations
cor_melted = corVector(cor_expr)
max(cor_melted$cor)

# look for discontinuities in the correlation distributions
# this might indicate duplicated or mixed samples
pdf("aagmex_sample_correlations.pdf",width=5,height=4)
par(mfrow=c(1,2))
hist(cor_melted$cor,xlab="Pearson correlation",main="",cex.lab=1.2)
hist(cor_melted$cor,xlim=c(0.98,1),ylim=c(0,50),breaks=seq(-1,1,0.001),
     xlab="Pearson correlation",main="",cex.lab=1.2)
abline(v=0.995,col="red")
dev.off()

# remove questionable high correlations
ind = which(cor_melted$cor > 0.995)
samp_rem = as.vector( c(cor_melted$s1[ind], cor_melted$s2[ind]) )
ind_rem = sapply(samp_rem,function(x)which(colnames(expr_data)==x)) %>% unlist

expr_data1 = expr_data[,-ind_rem]

# remove corresponding samples from sex annotation
pheno_data1 = pheno_data[-ind_rem,]

# check
all(pheno_data1$sample == names(expr_data1))

# plot XIST
pheno_data1$sex[pheno_data1$sex=="M"] = "Male"
pheno_data1$sex[pheno_data1$sex=="F"] = "Female"
pdf("aagmex_XIST_mf.pdf",width=4,height=4)
indXIST = which(rownames(expr_data1)=="XIST")
dat_xist = expr_data1[indXIST,] %>% t
pltdat = cbind(dat_xist,pheno_data1$sex) %>% as.data.frame
names(pltdat) = c("XIST", "sex")
pltdat$XIST = pltdat$XIST %>% data.matrix %>% as.numeric
beeswarm(XIST~sex, data = pltdat, cex=0.25, pch=16, col=c("red","blue"),
         cex.lab=1.3, xlab="")
abline(h=8)
dev.off()

# remove male samples with high XIST
cut = 8
M_remove = dat_xist[which(dat_xist[pheno_data1$sex=="Male"] > cut),] %>% names
ind_remove = sapply(M_remove,function(x)which(colnames(expr_data1)==x)) %>% unlist
expr_data2 = expr_data1[,-ind_remove]
pheno_data2 = pheno_data1[-ind_remove,]

# check
all(names(expr_data2) == pheno_data2$accession)

# write data to file
write.table(expr_data2,"expr_data_aagmex_qc.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(pheno_data2,"phenotypes_aagmex_qc.txt",col.names=T,row.names=F,quote=F,sep="\t")


