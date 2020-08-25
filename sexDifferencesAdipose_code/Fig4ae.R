
############################################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Figure 4a,e
############################################################################

library(dplyr)
library(ggplot2)
library(gridExtra)

setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5")

############################################################################
# load data
############################################################################

# eqtl interaction case data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/CIdata.cases.RData"
load(fname)

# peer associations with ld data and hg19 coords
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/CIdata.ldsnps.hg19.RData"
load(fname)

# snp2tfbs bed data
fname = "snp2tfbs_JASPAR_CORE_2014_vert.bed"
snp2tfbs0 = read.table(fname,header=F,stringsAsFactors=F)
names(snp2tfbs0) = c("chr","start","end","ref","alt","rsid","ntfs","tfs","scorediff")

# eqtl & deg
deg.eqtl.genes = c("CCDC3","CLIC6","FADS1","GLDN","HSPA12A","MAP1B","MLPH","MMD","MYOT","NDRG4","NEO1","PDZD2","TBC1D9")

############################################################################
# integrate ld-association data with tfbs data
############################################################################

# combine data
ld.tfbs = merge(CIdata.ldsnps.hg19, snp2tfbs0, by=c("chr","start","end","ref","alt"))
names(ld.tfbs)[c(7,16)] = c("leadSNP","ldrsid")

# filter to remove snps with no effect on the tfbs score
# a positive score implies a larger PWM score in the alternate allele
ld.tfbs = ld.tfbs %>% filter(scorediff != 0)
ld.tfbs = ld.tfbs %>% filter(scorediff != "0,0") 
ld.tfbs = ld.tfbs %>% filter(scorediff != "0,0,0")
ld.tfbs = ld.tfbs %>% filter(scorediff != "0,0,0,0")

# filter for ld = 0.8
# save(ld.tfbs, file="ld.tfbs.RData")
# load("ld.tfbs.RData")
ld.tfbs = ld.tfbs %>% filter(R2 > 0.8)

# for each association, determine the number of LD variants overlapping tfbs
ntfbs = data.frame(assoc = unique(ld.tfbs$genename), ntfbs = 0)
for(ii in 1:nrow(ntfbs)){
  ge = ntfbs$assoc[ii]
  ind = which(ld.tfbs$genename == ge)
  ntfbs$ntfbs[ii] = length(ind)
}

# basic quantification of loci: 881
length(unique(ld.tfbs$genename))

# identify associations of interest
ind = grep("PPARG",ld.tfbs$tfs)
pparg = ld.tfbs[ind,]
length(unique(pparg$gene)) # 36
ind = grep("CEBP",ld.tfbs$tfs)
cebp = ld.tfbs[ind,]
length(unique(cebp$gene)) # 16
ind = grep("FADS1",ld.tfbs$genename)
fads1 = ld.tfbs[ind,]
ind = grep("NEO1",ld.tfbs$genename)
neo1 = ld.tfbs[ind,]
ind = grep("MYOT",ld.tfbs$genename)
myot = ld.tfbs[ind,]

# for manuscript text
ind = grep(deg.eqtl.genes[12],ld.tfbs$genename)
ld.tfbs[ind,]

############################################################################
# examine enrichment of dominant TFs
############################################################################

library(scales)
library(enrichR)
library(ggplot2)

# filter the associations to include only the top TF
toptf = sapply(ld.tfbs$tfs,function(x){strsplit(x,",")[[1]][1]})

# for dimer sites, keep both
toptf = sapply(toptf,function(x){strsplit(x,"_")[[1]]}) %>% unlist
toptf = sapply(toptf,function(x){strsplit(x,"-")[[1]][1]}) %>% unlist

# convert to upper case for non-human identifiers
toptf = unique(toupper(toptf))

# manual filter
rem = c(28,96,125,163)
toptf[rem]
toptf = toptf[-rem]

# select databases
listEnrichrDbs()
db1 = c("Reactome_2016","KEGG_2019_Human","WikiPathways_2019_Human","Panther_2016")
db2 = c("GO_Biological_Process_2018","GO_Cellular_Component_2018","GO_Molecular_Function_2018")
db = c(db1,db2)

# get enrichment results
res = enrichr(genes=toptf, databases=db)

# filter each set of results
fdrthresh = 0.01
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

# filter for the number of genes
nge.thresh = 10
res_filtered = res_filtered %>% filter(nge > nge.thresh)

# select plot results
res_plt = res_filtered[,c(1,3,4,7)]
namen = sapply(res_plt$Term,function(x){strsplit(x,"[(]")[[1]][1]})
namen = sapply(namen,function(x){strsplit(x,"_")[[1]][1]})
res_plt$Term = namen

# generate plot
res_plt$Term = factor(res_plt$Term,levels=rev(res_plt$Term))
res_plt$fdr = -log10(res_plt$Adjusted.P.value)
plt = ggplot(res_plt, aes(x=Odds.Ratio, y=Term, size=nge, color=fdr)) + 
  geom_point() + xlim(c(4,45)) + 
  scale_color_gradient(low="gray", high="black", limits=c(4,30), oob=squish) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(plt)

# save plot & results table
pdf("Fig5a.pdf", onefile = FALSE, height=4, width=5)
print(plt)
dev.off()
fname = "TableS11.txt"
write.table(res_filtered,fname,col.names=T,row.names=F,quote=F,sep="\t")


############################################################################
# integrate with TF differential expression and eqtl data - case 1
############################################################################

# overlap data
load("ld.tfbs.RData")

# gtex deg data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/gtex_deg.txt"
gtex_deg = read.table(fname,header=T,stringsAsFactors=F)

# eqtl confidence interval data
file = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/CIdata.RData"
load(file)

# get TFs from DEG data
qfilt = 0.05
gtex.deg.tf = gtex_deg[gtex_deg$gene %in% toptf,]
qval = p.adjust(gtex.deg.tf$P.Value, method="BH")
gtex.deg.tf = gtex.deg.tf %>% mutate(qval = qval)
gtex.deg.tf = gtex.deg.tf %>% filter(qval < qfilt)
gtex.deg.tf = gtex.deg.tf[order(gtex.deg.tf$absFC,decreasing=T),]

# integrate TF DEG data with confidence interval data
CI = CIdata[,c(1,2,4,7:12)]
CI.case = merge(CIdata.cases, CI, by=c("snps","gene","pvalue"))

# find cases where CIs bracket zero (0, else sign of beta)
ci0 = function(dat){
  ind = apply(dat,1,function(x){
    if(sign(x[1]) != sign(x[2])){
      return(0)
    } else {
      return(sign(x[1]))
    } })
  return(ind)
}
indm = ci0(CI.case[,11:12])
indf = ci0(CI.case[,13:14])
CI.case = CI.case %>% mutate(m0 = indm, f0 = indf)
CI.case = CI.case[,c(1:3,7:10,15:16)]
# save(CI.case, file="CI.case.RData")
# load("CI.case.RData")

# integrate with DEG data
# for dimer motif, keep the tf with the most significant differential expression
toptf = sapply(ld.tfbs$tfs,function(x){strsplit(x,",")[[1]][1]})
toptf = sapply(toptf,function(x){
  tt = strsplit(x,"_")[[1]]
  if(length(tt) > 1){
    d = gtex.deg.tf[gtex.deg.tf$gene %in% tt,]
    d = d[order(d$qval,decreasing=F),]
    if(nrow(d)>0){
      return(d$gene[1])
    } else{return(tt[1])}
  } else(return(tt))
}) %>% unlist
ld.tfbs = ld.tfbs %>% mutate(toptf = toupper(toptf))
scorediff = sapply(ld.tfbs$scorediff,function(x){strsplit(x,",")[[1]][1]}) %>% as.numeric
ld.tfbs$scorediff = scorediff
ci.deg = merge(CI.case[,c(2,6:9)], ld.tfbs, by="gene")
names(gtex.deg.tf)[10] = "toptf"
ci.deg = merge(ci.deg, gtex.deg.tf[,c(8,10:11)], by="toptf") 

# process case 1 data 
# save(ci.deg.case1, file="ci.deg.case1.RData")
# save(ci.deg.case1.filtered, file="ci.deg.case1.filtered.RData")
# load("ci.deg.case1.filtered.RData")
ci.deg.case1 = ci.deg %>% filter(case == 1)
ci.deg.case1.M = ci.deg.case1 %>% filter(f0==0, ratioFC<0)
ci.deg.case1.F = ci.deg.case1 %>% filter(m0==0, ratioFC>0)
ci.deg.case1.filtered = rbind(ci.deg.case1.F, ci.deg.case1.M)

# loci accounted
length(unique(ci.deg.case1.F$genename)) # 9
length(unique(ci.deg.case1.M$genename)) # 7
unique(ci.deg.case1.F$genename)
unique(ci.deg.case1.M$genename)
unique(ci.deg.case1.F$toptf)
unique(ci.deg.case1.M$toptf)


############################################################################
# plot TF differential expression - case 1
############################################################################

library(dplyr)
library(ggplot2)
library(gridExtra)

setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig5")

fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig1/fig1a.RData"
load(fname)

# function to generate a single plot (gtex)
plt.gen = function(datF=NULL,datM=NULL,gene=NULL){
  
  # zscore scaled data
  ind.gtex = which(colnames(datF) == gene)
  gtex.dat = c(datF[,ind.gtex], datM[,ind.gtex])
  gtex.sex = c(rep("F",length(datF[,ind.gtex])), rep("M",length(datM[,ind.gtex])))
  gtex = data.frame(expr=scale(gtex.dat), sex=gtex.sex, data="gtex")
  gene.plt.dat = gtex
  
  # generate the plot
  gene.plot = ggplot(data=gene.plt.dat, aes(x=data,y=expr,fill=sex)) + geom_boxplot() +
    scale_fill_manual(values=c("red", "blue")) + ylab(paste0(gene," z-score")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  return(gene.plot)
  
} # plt.gen

# function for plotting a vector of genes
plot.genes = function(datF=NULL, datM=NULL, genes=NULL){
  plts = list()
  for(ii in genes){
    plts[[ii]] = plt.gen(datF=dat_F_gtex, datM=dat_M_gtex, gene=ii)
  }
  return(plts)
} # plot.genes


# plots for female case 1 factors
genes = c("FOXF2","FOXO1","MECOM","PBX1","STAT1","STAT5A")
plts = plot.genes(datF=dat_F_gtex, datM=dat_M_gtex, genes=genes)
pdf("fig5_case1TFs_F.pdf", onefile=FALSE, height=12, width=14)
marrangeGrob(grobs=plts, nrow=3, ncol=4, top=NULL)
dev.off()

# plots for male case 1 factors
genes = c("E2F3","EGR1","MYC","NRF1","SREBF2","TCF3","ZEB1")
plts = plot.genes(datF=dat_F_gtex, datM=dat_M_gtex, genes=genes)
pdf("fig5_case1TFs_M.pdf", onefile=FALSE, height=12, width=14)
marrangeGrob(grobs=plts, nrow=3, ncol=4, top=NULL)
dev.off()



############################################################################
# integrate with TF differential expression and eqtl data - case 2
############################################################################

# eqtl confidence intervals
load("CI.case.RData")

# ld-eqtl data
load("ld.tfbs.RData")

# gtex deg data
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/gtex_deg.txt"
gtex_deg = read.table(fname,header=T,stringsAsFactors=F)

# eqtl confidence interval data
file = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/CIdata.RData"
load(file)

fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig1/fig1a.RData"
load(fname)

# isolate case 2 data with >1 TBBS
# write data for manual analysis
case2 = ld.tfbs %>% filter(case==2, ntfs>1)
write.table(case2,"case2.txt",col.names=T,row.names=F,sep="\t")

length(unique(case2$genename))
