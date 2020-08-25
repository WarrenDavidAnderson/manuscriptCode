
############################################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Figure 3a
############################################################################

library(dplyr)
library(gridExtra)

setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4")

############################################################################
# load data
############################################################################

# load data for the eqtl plots
file = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4/gtexeqtlv2.RData"
load(file)
file = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/resids.RData"
load(file)

# intersect eqtl data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4/intersect.eqtl.RData"
load(fname)

# union eqtl data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig4/union.eqtl.RData"
load(fname)

# subq adipose DEG data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/all_fc_v8.txt"
deg0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F)

# METSIM correlation data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig3/metsimdegcor.RData"
load(fname)

# PEER cases
file = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/CIdata.cases.RData"
load(file)

# genotype annotation
# cat GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt > gtex.geno.ann.txt
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/gtex.geno.ann.txt"
geno.ann = read.table(fname,stringsAsFactors=F,header=T)


############################################################################
## filter eqtl-cases data based on QCd SNPs
############################################################################

# filter eqtl data based on QCd SNPs
CIdata.cases = CIdata.cases[CIdata.cases$snps %in% geno.ann$variant_id,]


############################################################################
# load peer-based eqtl data
############################################################################

library(MatrixEQTL)
library(dplyr)

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl"
setwd(dir)


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


############################################################################
# eqtl plots
############################################################################

library(ggplot2)

# deg genes with eqtl interaction (union of 3 normalization methods)
egene.deg = intersect(deg0$gene, union.eqtl$geneid) # 33 / 162

# further intersect with metsim correlation genes
deg.cor = intersect(deg0$gene, metsimdegcor$id) # 142 / 33
egene.deg.cor = intersect(egene.deg, metsimdegcor$id) # 31 / 33

# look at PEER eqtl case data (includes cases 1,2)
# save(pcase, file="pcase.RData")
pcase = CIdata.cases[CIdata.cases$genename %in% egene.deg.cor,] # 8

# select for case 1,2,3
c1 = CIdata.cases[CIdata.cases$case==1,]
c2 = CIdata.cases[CIdata.cases$case==2,]
c3 = CIdata.cases[CIdata.cases$case==3,]


# genes for ploting cases
case.genes = c("HSPA12A","MAP1B","URAD")
case.genes = c("FADS1","MYOT","NEO1")
case.genes = c1$genename[1:20]

case.genes = "BRD4"
case.genes = "FLT3"

case.genes = c("SSC5D","TBC1D24","URAD")

# female case 1
case.genes = c("LINC00884","RIN2","RP11-55K13.1","SQRDL","PABPC4","TIAM1","KRBA2","AKR1E2","CEP170")

# male case 1
case.genes = c("ATP5G3","XXYLT1","FLT3","LINC01238","AIFM3","LINC01314","TNFSF12" )

# eqtl & deg
case.genes = c("CCDC3","CLIC6","FADS1","GLDN","HSPA12A","MAP1B","MLPH","MMD","MYOT","NDRG4","NEO1","PDZD2","TBC1D9")


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
    geom_smooth(data=dat,aes(genotype+1,expression),method='lm',col="white",size=2) +
    scale_fill_manual(values=c("red", "blue")) +
    xlab(paste0("genotype (",rs,")")) +
    ylab(paste0(gene," expression")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    scale_x_discrete(labels=c("0" = ref, "1" = het, "2" = alt))
  return(plt)
}

# get the plots for the stated genes
snp.ids = rownames(snps)
snp.sms = colnames(snps)
nSlice = snps$nSlices()
sliceSize = snps$fileSliceSize
plts = list()
for(assoc.gene in case.genes){
  
  # specify association for plotting
  assoc.ind = which(CIdata.cases$genename == assoc.gene)
  snp.tgt = CIdata.cases$snps[assoc.ind]
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
  gene.tgt = CIdata.cases$gene[assoc.ind]
  ind = which(rownames(resids) == gene.tgt)
  ge = resids[ind,]
  ge.ind = sapply(names(gt),function(x){which(names(ge)==x)})
  ge = ge[ge.ind]
  model = model[ge.ind,]
  
  # generate the plot
  rsid = geno.ann$rs_id_dbSNP151_GRCh38p7[which(geno.ann$variant_id == snp.tgt)]
  sex = as.character(model$GENDER)
  pltdat = combine.expr.geno.sex(expression=ge, genotype=t(gt)[,1], sex=sex, male=0, female=1)
  plt = single.plot(dat=pltdat, gene=assoc.gene, rs=rsid, geno=as.character(snp.tgt))
  plts[[assoc.gene]] = plt
  
}

# Fig 4 plots
pdf("eqtl&deg.pdf", onefile=TRUE, height=12, width=14)
marrangeGrob(grobs=plts, nrow=3, ncol=4, top=NULL)
dev.off()






