

############################################################################
# Integrative analysis of sex differences in adipose tissue gene regulation
# Figure 2
# see vignette: 
# see script: 
############################################################################

############################################################################
# Run GSEA for GTEx separated by AA and non-AA subjects
# compare with AAGMEx
############################################################################


############################################################
## load data
############################################################

library(dplyr)

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/GSEA"
setwd(dir)

# Human DEGs for three data sets
dat = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/GSEA/deg.RData"
load(dat)

# GTEx subject annotation
fname = "Subject_Phenotypes.csv"
ann0 = read.delim(fname,header=T,stringsAsFactors=F,sep="\t", fill=T)
ann = ann0[,1:11]


############################################################
## separate GTEx data based on race
############################################################

# use gene ids for gtex
all(colnames(dat_F_gtex) == colnames(dat_M_gtex))
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/gencode_gene_map.txt"
ann_gene0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F,quote = "")
ggenes = sapply(colnames(dat_F_gtex),function(x){ann_gene0$gene_name[which(ann_gene0$gene_id==x)]})
colnames(dat_F_gtex) = colnames(dat_M_gtex) = ggenes

# get race annotation
subjidF = sapply(rownames(dat_F_gtex),function(x){
  res = strsplit(x,"-")[[1]]
  out = paste0(res[1],"-",res[2])
  return(out)
})
subjidM = sapply(rownames(dat_M_gtex),function(x){
  res = strsplit(x,"-")[[1]]
  out = paste0(res[1],"-",res[2])
  return(out)
})
indF = sapply(subjidF,function(x){which(ann$SUBJID==x)})
indM = sapply(subjidM,function(x){which(ann$SUBJID==x)})
raceF = ann$RACE[indF]
raceM = ann$RACE[indM]

# look and number of AAs for males and females
# https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/
# variable.cgi?study_id=phs000424.v7.p2&phv=169064&phd=3910&pha=&pht=2742&phvf=&phdf=&phaf=&phtf=&dssp=2&consent=&temp=1
nAAf = length(which(raceF==2)) # 25
nAAm = length(which(raceM==2)) # 46



