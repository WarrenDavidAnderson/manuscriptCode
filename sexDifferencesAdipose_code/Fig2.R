
############################################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Figure 2
############################################################################

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/GSEA"
setwd(dir)

# Human DEGs for three data sets
dat = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/GSEA/deg.RData"
load(dat)

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

############################################################################
# gtex GSEA
############################################################################

# use gene ids for gtex
all(colnames(dat_F_gtex) == colnames(dat_M_gtex))
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/gencode_gene_map.txt"
ann_gene0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F,quote = "")
ggenes = sapply(colnames(dat_F_gtex),function(x){ann_gene0$gene_name[which(ann_gene0$gene_id==x)]})
colnames(dat_F_gtex) = colnames(dat_M_gtex) = ggenes

# specify phenotypes
pheno.input = list()
pheno.input[["phen"]] = c("Female","Male")
pheno.input[["class.v"]] = c(rep(0,nrow(dat_F_gtex)), rep(1,nrow(dat_M_gtex)))

# format expression data
expr.input = cbind(t(dat_F_gtex), t(dat_M_gtex)) %>% as.data.frame(stringsAsFactors=F)

# run GSEA.R and newGSEAplots.R
library(GSEA.plot)

##########################
# set up gene set library

# hallmark
data(hallmark.gs)
d0 = hallmark.gs

# KLF14 targets
data(transf)
data(transm)
tx = t(c("KLF14targets","source",transf[,1],transm[,1]))
d1 = add_to_database(database=d0, addition=tx)

# receptors/ligands
data(Kadoki_ligands.db)
data(Kadoki_receptors.db)
d2 = c(d1, Kadoki_ligands.db, Kadoki_receptors.db)

# transcription factors
data(ENCODE.db)
d3 = c(d2, ENCODE.db)

# run GSEA
pp.gtex = GSEAplots(input.ds.name=expr.input,
              input.cls.name=pheno.input, 
              gene.set.input=d3,
              doc.string="gtex", 
              nperm=1000,
              fdr.q.val.threshold=10,
              abs.val=F,
              gs.size.threshold.max=1e50, 
              bar_percent=0.1, gap_percent=0.1,
              under_percent=0.05,
              upper_percent=0.05,
              color_line="black",
              color_tick="black")

# generate plots
plot.ES(list.of.plots=pp.gtex$plots, plotname="gtex")


############################################################################
# decode GSEA
############################################################################

# use gene ids for decode
all(colnames(dat_F_decode) == colnames(dat_M_decode))
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/genes_decode.txt"
ann_gene0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F,quote="")
ggenes = sapply(colnames(dat_F_decode),function(x){ann_gene0$gene_list[which(ann_gene0$ID==x)]})
colnames(dat_F_decode) = colnames(dat_M_decode) = ggenes

# specify phenotypes
pheno.input = list()
pheno.input[["phen"]] = c("Female","Male")
pheno.input[["class.v"]] = c(rep(0,nrow(dat_F_decode)), rep(1,nrow(dat_M_decode)))

# format expression data
expr.input = cbind(t(dat_F_decode), t(dat_M_decode)) %>% as.data.frame(stringsAsFactors=F)

# run GSEA
pp.decode = GSEAplots(input.ds.name=expr.input,
                    input.cls.name=pheno.input, 
                    gene.set.input=d3,
                    doc.string="decode", 
                    nperm=1000,
                    fdr.q.val.threshold = 10,
                    abs.val=F,
                    gs.size.threshold.max=1e50, 
                    bar_percent=0.1, gap_percent=0.1,
                    under_percent=0.05,
                    upper_percent=0.05,
                    color_line="black",
                    color_tick="black")

# generate plots
plot.ES(list.of.plots=pp.decode$plots, plotname="decode")


############################################################################
# aagmex GSEA
############################################################################

# specify phenotypes
pheno.input = list()
pheno.input[["phen"]] = c("Female","Male")
pheno.input[["class.v"]] = c(rep(0,nrow(dat_F_aagmex)), rep(1,nrow(dat_M_aagmex)))

# format expression data
expr.input = cbind(t(dat_F_aagmex), t(dat_M_aagmex)) %>% as.data.frame(stringsAsFactors=F)

# run GSEA
pp.aagmex = GSEAplots(input.ds.name=expr.input,
                      input.cls.name=pheno.input, 
                      gene.set.input=d3,
                      doc.string="aagmex", nperm=1000,
                      fdr.q.val.threshold = 10,
                      abs.val=F,
                      gs.size.threshold.max=1e50, 
                      bar_percent=0.1, gap_percent=0.1,
                      under_percent=0.05,
                      upper_percent=0.05,
                      color_line="black",
                      color_tick="black")

# generate plots
plot.ES(list.of.plots=pp.aagmex$plots, plotname="aagmex")


############################################################################
# process GSEA results
############################################################################

dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/GSEA"
setwd(dir)

library(dplyr)

# load data
baseF = ".SUMMARY.RESULTS.REPORT.Female.txt"
baseM = ".SUMMARY.RESULTS.REPORT.Male.txt"
datasets = c("gtex","decode","aagmex")
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
n3 = paste0(namen,"_",datasets[3])
gsea.res = Reduce(function(x,y) merge(x = x, y = y, by = "GS"), gsea.res0)
names(gsea.res)[2:ncol(gsea.res)] = c(n1, n2, n3)

# filter by FDR
fdr.cut = 0.1
gsea.res.filt = gsea.res %>% filter(FDR_gtex < fdr.cut & FDR_decode < fdr.cut & FDR_aagmex < fdr.cut) # zero
gsea.res.filt = gsea.res %>% filter(FDR_gtex < fdr.cut & FDR_decode < fdr.cut) # 6

# aagmes sig results
gsea.res %>% filter(FDR_aagmex < fdr.cut)

# generate heatmap
library(NMF)
pltdat = gsea.res.filt[,c(4,10,16)]
names(pltdat) = datasets
rownames(pltdat) = gsea.res.filt$GS
pdf("hm1.pdf")
aheatmap(pltdat,Colv=NA,breaks=0)
dev.off()

# gtex enrich plots
library(gridExtra)
namen = names(pp.gtex$gene.set.leading)
ind1 = which(namen == "HALLMARK_ADIPOGENESIS")
ind2 = which(namen == "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
ind3 = which(namen == "HALLMARK_FATTY_ACID_METABOLISM")
plts = pp.gtex$plots[c(ind1,ind2,ind3)]
pdf("gtex1.pdf", onefile = FALSE)
marrangeGrob(grobs=plts, nrow=2, ncol=2, top=NULL)
dev.off()


############################################################################
# examine fold change correlations for annotations of interest
############################################################################

# download full differential expression data
fname = "gtex_deg.txt"
gtex = read.table(fname,header=T,stringsAsFactors=F)
fname = "decode_deg.txt"
decode = read.table(fname,header=T,stringsAsFactors=F)
fname = "aagmex_deg.txt"
aagmex = read.table(fname,header=T,stringsAsFactors=F)

# get unique genes for gtex and decode
min.p.indices = function(dat=NULL){
  ind_keep = c()
  for(ii in 1:length(unique(dat$gene))){
    ind = which(dat$gene == unique(dat$gene)[ii])
    if(length(ind)==1){ind_keep = c(ind_keep, ind)}
    if(length(ind)>1){
      ind2 = which(dat$adj.P.Val[ind] == min(dat$adj.P.Val[ind]))
      ind_keep = c(ind_keep, ind[ind2])
    }
  }
  if(length(ind_keep)==length(unique(dat$gene))){return(ind_keep)}
}
gtex = gtex[min.p.indices(gtex),]
decode = decode[min.p.indices(decode),]

# merge the fold change data sets
gtex = gtex %>% select(gene, logFC)
decode = decode %>% select(gene, logFC)
aagmex = aagmex %>% select(gene, logFC)
names(gtex)[2] = "gtex"
names(decode)[2] = "decode"
names(aagmex)[2] = "aagmex"
FC  = Reduce(function(x,y)merge(x,y,by="gene",all=TRUE), list(gtex,decode,aagmex))
FC = FC[rowSums(is.na(FC))==0,]

# generate correlation hm plots
cw = 20
ch = 20
pdf("gseaCorHm.pdf")
for(ii in 1:nrow(gsea.res.filt)){
  namen = gsea.res.filt$GS[ii]
  ann.genes = pp.aagmex$gene.set.reference.matrix[[namen]]
  fc = FC[FC$gene %in% ann.genes, 2:4]
  cmat = cor(fc)
  aheatmap(cmat,main=namen,breaks=c(0),Colv=NA,Rowv=NA,cellwidth=cw,cellheight=ch)
}
dev.off()

# save(gsea.res.filt, file="gsea.res.filt.RData")


############################################################################
# standard fet enrichment analysis
############################################################################

res_plt = bartF.filt.fdr[1:10,]
res_plt$TF = factor(res_plt$TF,levels=rev(res_plt$TF))
res_plt$P.value = -log10(res_plt$fdr)
plt = ggplot(res_plt, aes(x=statistic, y=TF, size=max_auc, color=P.value)) + 
  geom_point() + xlim(c(5,20)) + 
  scale_color_gradient(low="gray", high="black", limits=c(6,16), oob=squish) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(plt)

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/GSEA"
setwd(dir)

# Human DEGs for three data sets
dat = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/GSEA/deg.RData"
load(dat)

library(dplyr)
library(scales)
library(ggplot2)
library(enrichR)

# select databases
listEnrichrDbs()
db1 = c("Reactome_2016","KEGG_2019_Human","WikiPathways_2019_Human","Panther_2016")
db2 = c("GO_Biological_Process_2018","GO_Cellular_Component_2018","GO_Molecular_Function_2018")
db = c(db1,db2)

# get enrichment results
# save(res, file="enrichr.res.RData")
res = enrichr(genes=all_fc$gene, databases=db)

# filter each set of results
fdrthresh = 0.1
pthresh = 0.001
ngenes = 1
res_filtered = c()
for(ii in 1:length(res)){
  dat = res[[ii]]
  #dat = dat %>% filter(Adjusted.P.value < fdrthresh)
  dat = dat %>% filter(P.value < pthresh)
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

# n enrichments
length(unique(res_filtered$Term))

# output results
fname = "TableS3.txt"
write.table(res_filtered,fname,col.names=T,row.names=F,quote=F, sep="\t")

# select plot results
# https://www.ncbi.nlm.nih.gov/pubmed/30543876
# https://www.ncbi.nlm.nih.gov/pubmed/24838110
res_plt = res_filtered[,c(1,2,4,7)]
ind = c(1:9,11,14)
res_plt = res_plt[ind,]
namen = sapply(res_plt$Term,function(x){strsplit(x,"[(]")[[1]][1]})
namen = sapply(namen,function(x){strsplit(x,"_")[[1]][1]})
res_plt$Term = namen

# generate plot
res_plt$Term = factor(res_plt$Term,levels=rev(res_plt$Term))
res_plt$P.value = -log10(res_plt$P.value)
plt = ggplot(res_plt, aes(x=Odds.Ratio, y=Term, size=nge, color=P.value)) + 
  geom_point() + xlim(c(0,45)) + 
  scale_color_gradient(low="gray", high="black", limits=c(3,4.5), oob=squish) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(plt)

# save plot
pdf("enrichFET.pdf", onefile = FALSE, height=2.9, width=6)
print(plt)
dev.off()

############################################################################
# bart analysis
############################################################################

bart.dir = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/BART"
setwd(bart.dir)

library(dplyr)

# import data
fname = paste0(bart.dir,"/Female/femaleHuman_1576773085946449___Geneset_bart_results.txt")
bartF = read.table(fname,header=T,stringsAsFactors=F)
fname = paste0(bart.dir,"/Male/maleHuman_15767735883041396___Geneset_bart_results.txt")
bartM = read.table(fname,header=T,stringsAsFactors=F)

# identify FDRs
fdr = p.adjust(bartF$pvalue,method="BH")
bartF = bartF %>% mutate(fdr=fdr)
fdr = p.adjust(bartM$pvalue,method="BH")
bartM = bartM %>% mutate(fdr=fdr)

# filter the data by p-value
pcut = 0.001
bartF.filt.pval = bartF %>% filter(irwin_hall_pvalue < pcut)
bartM.filt.pval = bartM %>% filter(irwin_hall_pvalue < pcut)

# filter the data by FDR
fdrcut = 0.05
bartF.filt.fdr = bartF %>% filter(fdr < fdrcut) %>% mutate(class="female")
bartM.filt.fdr = bartM %>% filter(fdr < fdrcut) %>% mutate(class="male")
bartF.filt.fdr = bartF.filt.fdr[order(bartF.filt.fdr$fdr),]
bartM.filt.fdr = bartM.filt.fdr[order(bartM.filt.fdr$fdr),]

# output summary table
tableS4 = rbind(bartF.filt.fdr,bartM.filt.fdr)[c(1:3,8,5,9)]
write.table(tableS4,"TableS4.txt",col.names=T,row.names=F,quote=F,sep="\t")

# generate plot, female
library(scales)
library(ggplot2)
res_plt = bartF.filt.fdr[1:10,]
res_plt$TF = factor(res_plt$TF,levels=rev(res_plt$TF))
res_plt$P.value = -log10(res_plt$fdr)
plt = ggplot(res_plt, aes(x=statistic, y=TF, size=max_auc, color=P.value)) + 
  geom_point() + xlim(c(5,20)) + 
  scale_color_gradient(low="gray", high="black", limits=c(6,16), oob=squish) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(plt)
plts=list()
plts[[1]] = plt

# save plot
pdf("bartF.pdf", onefile = FALSE, height=2.9, width=6)
print(plt)
dev.off()

# generate plot, male
res_plt = bartM.filt.fdr[1:10,]
res_plt$TF = factor(res_plt$TF,levels=rev(res_plt$TF))
res_plt$P.value = -log10(res_plt$fdr)
plt = ggplot(res_plt, aes(x=statistic, y=TF, size=max_auc, color=P.value)) + 
  geom_point() + xlim(c(5,30)) + 
  scale_color_gradient(low="gray", high="black", limits=c(6,16), oob=squish) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(plt)
plts[[2]] = plt

# save plot
pdf("bartM.pdf", onefile = FALSE, height=2.9, width=6)
print(plt)
dev.off()

# save plots together
library(gridExtra)
pdf("bartRes.pdf", onefile = FALSE, height=2.9, width=6)
marrangeGrob(grobs=plts, nrow=1, ncol=2, top=NULL)
dev.off()

