

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
library(ggplot2)
library(gridExtra)

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/GSEA/GTExAA"
setwd(dir)

# GTEx subject annotation
fname = "Subject_Phenotypes.csv"
ann0 = read.delim(fname,header=T,stringsAsFactors=F,sep="\t", fill=T)
ann = ann0[,1:11]

# load in unprocessed data
datdir = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtexCovar/"
fname = "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
fname = paste0(datdir,fname)
covar0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F,quote="")
fname = "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
fname = paste0(datdir,fname)
subj0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F,check.names=F)
fname = "subq_gtex_invNorm.txt"
fname = paste0(datdir,fname)
expr0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F,check.names=F)
fname = "Adipose_Subcutaneous.v8.covariates.txt"
fname = paste0(datdir,fname)
ecov0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F,check.names=F) %>% t
names = ecov0[1,]
ecov0 = ecov0[-1,] %>% as.data.frame
names(ecov0) = names


############################################################
## basic data organization
############################################################

# isolate subcutaneous data in covariate file
covar1 = covar0[covar0$SAMPID %in% names(expr0),]
all(covar1$SAMPID == names(expr0))

# isolate subcutaneous data in subject attribute file
subq_subs = sapply(names(expr0),function(n){
  x = strsplit(n,"-")[[1]]
  return(paste0(x[1],"-",x[2]))
})
subj1 = subj0[subj0$SUBJID %in% subq_subs,]
subj1 = subj1[match(subq_subs,subj1$SUBJID),]
all(names(expr0) == names(subq_subs))
all(subj1$SUBJID == subq_subs)

# verify data organization in eQTL covariate file
all(rownames(ecov0) == subj1$SUBJID)
all(rownames(ecov0) == subq_subs)

# combine covariates
covar = cbind(covar1,subj1)

# recode the sex indicator
covar$GENDER[covar$SEX == 1] = "Male"
covar$GENDER[covar$SEX == 2] = "Female"

# incorporate race information
subjid = sapply(covar$SAMPID,function(x){
  res = strsplit(x,"-")[[1]]
  return(paste0(res[1],"-",res[2]))
})
inds = sapply(subjid,function(x){which(ann$SUBJID==x)})
race = ann$RACE[inds]
race = gsub("1","Asian",race)
race = gsub("2","AA",race)
race = gsub("3","White",race)
race = gsub("4","Native",race)
race = gsub("98","Noreport",race)
race = gsub("99","Unknown",race)
covar = covar %>% mutate(race=race)

# use centers of the age ranges
ageCenter = function(vec=NULL){
  age_ranges = unique(vec)
  age_table = matrix(c(0),length(age_ranges),2) %>% as.data.frame
  names(age_table) = c("range","mean")
  for(ii in 1:nrow(age_table)){
    age_mean = strsplit(age_ranges[ii],"-") %>% unlist %>% as.numeric %>% mean
    age_table[ii,1] = age_ranges[ii]
    age_table[ii,2] = age_mean
  }
  ages = age_table$mean[sapply(vec,function(x)which(age_table$range==x))%>%unlist]
  return(ages)
}
covar$AGE = ageCenter(covar$AGE)

############################################################
## implement sva (independent of race and other key covariates)
############################################################

library(sva)

# specify full and null model matrices
all(rownames(ecov0) == subjid)
covar$Platform = ecov0$platform
modFull = model.matrix(~as.factor(GENDER)+AGE+SMRIN+as.factor(Platform)+as.factor(race), data=covar)
modNull = model.matrix(~1,data=covar)

# identify the surrogate variables
svadat0 = expr0 %>% data.matrix
LVs = sva(svadat0, modFull, modNull)
LV = LVs$sv
colnames(LV) = paste0("sv",c(1:ncol(LV)))

############################################################
## regress out covariates/svs
############################################################

# implement regression and keep residuals
# correct for everything other than sex and race
model0 = covar %>% select(AGE,SMRIN,Platform)
model0$Platform = as.factor(model0$Platform)
model_resid = cbind(model0,LV)
fitdat = expr0 %>% t
fit = apply(fitdat,2,function(x){
  dat = as.data.frame(cbind(x,model_resid))
  fit_regdat = lm(x~.,data=dat)
  return(fit_regdat$residuals)
})
resids_expr = t(fit)

############################################################
## implement differential expression analysis
############################################################

library(limma)

# function for DEG analysis
# implement linear model analysis with eBayes and BH adjustments
deg.analysis = function(dat_Null=NULL,dat_Alt=NULL,null=NULL,alt=NULL){
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

# function for getting indices of unique genes with min p-values
min.p.indices = function(dat=NULL){
  ind_keep = c()
  for(ii in 1:length(unique(dat$ID))){
    ind = which(dat$ID == unique(dat$ID)[ii])
    if(length(ind)==1){ind_keep = c(ind_keep, ind)}
    if(length(ind)>1){
      ind2 = which(dat$adj.P.Val[ind] == min(dat$adj.P.Val[ind]))
      ind_keep = c(ind_keep, ind[ind2])
    }
  }
  if(length(ind_keep)==length(unique(dat$ID))){return(ind_keep)}
}

# specify groups for comparison
ind_F = which(covar$GENDER == "Female")
ind_M = which(covar$GENDER == "Male")
null <- "Female"
alt <- "Male"
dat_Null = t(resids_expr[,ind_F]) # Null model - female
dat_Alt = t(resids_expr[,ind_M]) # Alt model - male

# implement DEG analysis and select gene ids with lowest FDRs
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/gencode_gene_map.txt"
ann_gene0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F,quote = "")
output = deg.analysis(dat_Null=dat_Null,dat_Alt=dat_Alt,null=null,alt=alt)
genes = sapply(rownames(output),function(x){ann_gene0$gene_name[which(ann_gene0$gene_id==x)]})
output = output %>% mutate(ID = genes, ens = rownames(output))
output = output[min.p.indices(output),]

# use gene ids for gtex
resids = resids_expr[rownames(resids_expr) %in% output$ens,]
ggenes = sapply(rownames(resids),function(x){ann_gene0$gene_name[which(ann_gene0$gene_id==x)]})
rownames(resids) = ggenes

# save.image("gtexAA.RData")
# load("gtexAA.RData")

############################################################
## run gsea for aa and non-aa gtex subjects
############################################################

library(dplyr)
library(ggplot2)
library(gridExtra)
library(GSEA.plot)

# run GSEA.R and newGSEAplots.R

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/GSEA/GTExAA"
setwd(dir)

load("gtexAA.RData")

# separate by race and sex
all(covar$SAMPID == names(resids))
indM = which(covar$GENDER == "Male")
indF = which(covar$GENDER == "Female")
indAA = which(covar$race == "AA")
indna = which(covar$race != "AA")

##########################
# set up gene set library

# hallmark
data(hallmark.gs)
d0 = hallmark.gs

# KLF14 targets
data(transf)
data(transm)
tx = t(c("HALLMARK_KLF14","http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_KLF14",transf[,1],transm[,1]))
tx = t(c("HALLMARK_KLF14","source",transf[,1],transm[,1]))
d1 = add_to_database(database=d0, addition=tx)

# receptors/ligands
data(Kadoki_ligands.db)
data(Kadoki_receptors.db)
d2 = c(d1, Kadoki_ligands.db, Kadoki_receptors.db)

# transcription factors
data(ENCODE.db)
d3 = c(d2, ENCODE.db)

################################
## African American GSEA

# set the data frames
dat_F_gtex = resids[,intersect(indF,indAA)]
dat_M_gtex = resids[,intersect(indM,indAA)]

# specify phenotypes
pheno.input = list()
pheno.input[["phen"]] = c("Female","Male")
pheno.input[["class.v"]] = c(rep(0,ncol(dat_F_gtex)), rep(1,ncol(dat_M_gtex)))

# format expression data
expr.input = cbind(dat_F_gtex, dat_M_gtex) %>% as.data.frame(stringsAsFactors=F)

# run GSEA
AA.gtex = GSEAplots(input.ds.name=expr.input,
                    input.cls.name=pheno.input, gene.set.input=d3,
                    doc.string="gtex_AA", nperm=1000,
                    fdr.q.val.threshold = 10,abs.val=F,gs.size.threshold.max=1e50, bar_percent=0.1, gap_percent=0.1,
                    under_percent=0.05,upper_percent=0.05,color_line="black",
                    color_tick="black")

# generate plots
plot.ES(list.of.plots=AA.gtex$plots, plotname="gtex_AA")


################################
## non African American GSEA

# set the data frames
dat_F_gtex = resids[,intersect(indF,indna)]
dat_M_gtex = resids[,intersect(indM,indna)]

# specify phenotypes
pheno.input = list()
pheno.input[["phen"]] = c("Female","Male")
pheno.input[["class.v"]] = c(rep(0,ncol(dat_F_gtex)), rep(1,ncol(dat_M_gtex)))

# format expression data
expr.input = cbind(dat_F_gtex, dat_M_gtex) %>% as.data.frame(stringsAsFactors=F)

# run GSEA
AA.gtex = GSEAplots(input.ds.name=expr.input,
                    input.cls.name=pheno.input, gene.set.input=d3,
                    doc.string="gtex_nonAA", nperm=1000,
                    fdr.q.val.threshold = 10,abs.val=F,gs.size.threshold.max=1e50, bar_percent=0.1, gap_percent=0.1,
                    under_percent=0.05,upper_percent=0.05,color_line="black",
                    color_tick="black")

# generate plots
plot.ES(list.of.plots=AA.gtex$plots, plotname="gtex_nonAA")



############################################################
## process and plot AA gsea results
############################################################

library(dplyr)
library(ggplot2)
library(gridExtra)

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/GSEA/GTExAA"
setwd(dir)

# initial GSEA results
load("gsea.res.filt.RData")

# race-separated GSEA results
baseF = ".SUMMARY.RESULTS.REPORT.Female.txt"
baseM = ".SUMMARY.RESULTS.REPORT.Male.txt"
datasets = c("gtex_AA","gtex_nonAA")
gsea.res0 = list()
for(ii in datasets){
  fname = paste0(ii,baseF)
  datF = read.table(fname,header=T,stringsAsFactors=F,sep="\t")
  fname = paste0(ii,baseM)
  datM = read.table(fname,header=T,stringsAsFactors=F,sep="\t")
  dat = rbind(datF, datM)
  dat = dat[,c(1,2,4:8)]
  gsea.res0[[ii]] = dat
}

# merge the frames
namen = c("size","ES","NES","P","FDR","FWER")
n1 = paste0(namen,"_",datasets[1])
n2 = paste0(namen,"_",datasets[2])
gsea.res = Reduce(function(x,y) merge(x = x, y = y, by = "GS"), gsea.res0)
names(gsea.res)[2:ncol(gsea.res)] = c(n1, n2)
klf = grep("KLF14",gsea.res$GS)
gsea.res$GS[klf] = "KLF14targets"
gtex_aagmex0 = merge(gsea.res.filt, gsea.res, by="GS")

# subset the data
gtex_aagmex = gtex_aagmex0[,c(1,4,16,22,28)]
plt0 = gtex_aagmex[c(3,1,5,4,2,6),]
plt = plt0[,c(2,5,4,3)]
rownames(plt) = plt0$GS

# generate the heatmap
library(NMF)
pdf("hmAA.pdf")
aheatmap(plt,Colv=NA,Rowv=NA,breaks=0)
dev.off()


