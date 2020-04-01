
###########################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Figure 1
# HMDP DEG analysis
############################################################

library(dplyr)

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/HMDP"
setwd(dir)


################################################################################################
## load gene annotation data 
## use GEO microarray
################################################################################################

library(Biobase)
library(GEOquery)
library(limma)
library(dplyr)

gset <- getGEO("GSE64769", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL8759", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
fvarLabels(gset) <- make.names(fvarLabels(gset))

# gene annotation
gene_info0 = fData(gset) %>% select(ID, Gene.symbol, Chromosome.annotation)
chr = sapply(gene_info0$Chromosome.annotation,function(x){strsplit(x,",")[[1]][1]})
chr = sapply(chr,function(x){strsplit(x," ")[[1]][2]})

# remove sex chromosomes
rems = c("X","Y","NA")
gene_info = gene_info0 %>% mutate(chr = chr)
gene_info = gene_info[!(gene_info$chr %in% rems),]


############################################################
## read and format data
############################################################

fname = "HMDP_Adipose_Gene_Expression.csv"
expr0 = read.table(fname,sep=",",header=T,stringsAsFactors=F)

fname = "HMDP_Traits.csv"
pheno0 = read.table(fname,sep=",",header=T,stringsAsFactors=F)

# initial annotation for the phenotype set
annotation = pheno0 %>% select(mouse_number,Strain,Sex)

# set exression matrix 
# name expression columns by mouse number
expr = reshape(expr0[,c(1,2,5)], idvar="probesetID", timevar="mouse_number", direction="wide")
namen = sapply(names(expr)[-1],function(x){strsplit(x,"expression_value.")[[1]][2]})
names(expr)[2:ncol(expr)] = namen

# set phenotype matrix and annotation
pheno = reshape(pheno0[,c(1,2,5)], idvar="trait_name", timevar="mouse_number", direction="wide")
namen = sapply(names(pheno)[-1],function(x){strsplit(x,"value.")[[1]][2]})
names(pheno)[2:ncol(pheno)] = namen


############################################################
## select strains for expression analysis
## selarate male and female data, remove autosomes
############################################################

# identify strains with both M and F - expression and phenotype
strains.expr = annotation[names(expr)[-1] %in% annotation$mouse_number,]
strains.m = strains.expr %>% filter(Sex == "M") %>% select(Strain) %>% unique
strains.f = strains.expr %>% filter(Sex == "F") %>% select(Strain) %>% unique
strains.mf = intersect(strains.m$Strain, strains.f$Strain)

# filter expression data to contain selected strains
strains = sapply(names(expr)[-1],function(x){annotation$Strain[which(annotation$mouse_number==x)[1]]})
inds = which(strains %in% strains.mf == TRUE)
expr = expr[,c(1,inds+1)]

# remove sex chromosomes
inds.keep = which(expr$probesetID %in% gene_info$ID == TRUE)
expr = expr[inds.keep,]

# isolate male and female expression data
nums = names(expr)[-1]
sex = sapply(nums,function(x){annotation$Sex[which(annotation$mouse_number==x)[1]]})
ind.m = which(sex == "M")
ind.f = which(sex == "F")
expr.m = expr[,c(1,ind.m+1)]
expr.f = expr[,c(1,ind.f+1)]


############################################################
## differential expression analysis
############################################################

library(limma)

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
}

# get differential expression results
null = "Female"
alt = "Male"
m.expr = t(expr.m[,-1])
f.expr = t(expr.f[,-1])
colnames(m.expr) = colnames(f.expr) = expr$probesetID
output = deg.analysis(dat_Null=f.expr, dat_Alt=m.expr, null=null, alt=alt)

# add probe and gene information
output = output %>% mutate(probe = rownames(output))
genes = sapply(output$probe,function(x){gene_info$Gene.symbol[which(gene_info$ID==x)]})
genes[which(genes == "")] = "NA"
output = output %>% mutate(gene = genes)

# remove NAs and multiple annotations
output = output %>% filter(gene != "NA")
inds = grep("///",output$gene)
output = output[-inds,]

# output the results
fname = "hmdpDEG.txt"
write.table(output,fname,col.names=T,row.names=F,sep="\t",quote=F)






