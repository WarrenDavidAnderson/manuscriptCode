


############################################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Figure 1
############################################################################

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig1"
setwd(dir)

# Human DEGs for three data sets
dat = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/deg.RData"
load(dat)

library(NMF)
library(dplyr)
library(limma)
library(VennDiagram)
library(SuperExactTest)

############################################################################
# ma plots
############################################################################

# function for plotting DEG analysis results
deg.ma.plot = function(data=NULL,fc.cut=1.05, fdr.cut=0.05, all_fc=NULL){
  outf = output %>% filter(logFC > log2(fc.cut), adj.P.Val < fdr.cut)
  outm = output %>% filter(logFC < -log2(fc.cut), adj.P.Val < fdr.cut)
  ymax = max(output$logFC)
  ymin = min(output$logFC)
  yext = max(abs(ymax), abs(ymin))
  ylim = c(-1.05*yext, 1.05*yext)
  plot(output$AveExpr, output$logFC, col="gray", ylim=ylim,
       xlab="Average expression",ylab="Log2 fold change")
  points(outf$AveExpr, outf$logFC, col="red")
  points(outm$AveExpr, outm$logFC, col="blue")
  indsig = sapply(all_fc$gene,function(x){which(output$ID==x)}) %>% unlist
  points(output$AveExpr[indsig], output$logFC[indsig], col="cyan", pch=19)
  abline(h=0,lty=2)
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

# use gene ids for gtex
all(colnames(dat_F_gtex) == colnames(dat_M_gtex))
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/gencode_gene_map.txt"
ann_gene0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F,quote = "")
ggenes = sapply(colnames(dat_F_gtex),function(x){ann_gene0$gene_name[which(ann_gene0$gene_id==x)]})
colnames(dat_F_gtex) = colnames(dat_M_gtex) = ggenes

# use gene ids for decode
all(colnames(dat_F_decode) == colnames(dat_M_decode))
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/genes_decode.txt"
ann_gene0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F,quote="")
ggenes = sapply(colnames(dat_F_decode),function(x){ann_gene0$gene_list[which(ann_gene0$ID==x)]})
colnames(dat_F_decode) = colnames(dat_M_decode) = ggenes


# put the ma plots in a pdf
fname = "degMA.pdf"
pdf(fname)
par(mfrow=c(3,3))
output = deg.analysis(dat_Null=dat_F_gtex,dat_Alt=dat_M_gtex,null=null,alt=alt)
output = output[min.p.indices(output),]
deg.ma.plot(data=output, all_fc=all_fc)

output = deg.analysis(dat_Null=dat_F_decode,dat_Alt=dat_M_decode,null=null,alt=alt)
output = output[min.p.indices(output),]
deg.ma.plot(data=output, all_fc=all_fc)

output = deg.analysis(dat_Null=dat_F_aagmex,dat_Alt=dat_M_aagmex,null=null,alt=alt)
output = output %>% mutate(ID = rownames(output))
deg.ma.plot(data=output, all_fc=all_fc)

dev.off()

# save.image("fig1a.RData")
# load("fig1a.RData")

############################################################################
# significance of 3-way overlap --> 162
############################################################################

# elements of the Venn
match = paste0("match",1:162)
gtex = paste0("gtex",1:2692)
decode = paste0("decode",1:503)
aagmex = paste0("aagmex",1:1749)
gtex.decode = paste0("gtexdecode",1:150)
gtex.aagmex = paste0("gtexaagmex",1:287)
decode.aagmex = paste0("decodeaagmex",1:132)

# individual data sets
d1 = c(match,gtex,gtex.decode,gtex.aagmex)
d2 = c(match,decode,gtex.decode,decode.aagmex)
d3 = c(match,aagmex,gtex.aagmex,decode.aagmex)

# statistical test
x = list(d1, d2, d3)
n = max(ncol(dat_F_gtex), ncol(dat_F_decode), ncol(dat_F_aagmex))
res = MSET(x,n,lower.tail=FALSE,log.p=FALSE)

############################################################################
# generate histograms for 3 m/f deg analyses
############################################################################

fname = "degHIST.pdf"
pdf(fname)
par(mfrow=c(3,3)) 

histdat = all_fc$gtex_fc
bk = seq(floor(min(histdat)),ceiling(max(histdat)),by=0.05)
hist(histdat,breaks=bk, col="black",xlab="gtex")

histdat = all_fc$decode_fc
bk = seq(floor(min(histdat)),ceiling(max(histdat)),by=0.05)
hist(histdat,breaks=bk, col="black",xlab="decode")

histdat = all_fc$aagmex_fc
bk = seq(floor(min(histdat)),ceiling(max(histdat)),by=0.05)
hist(histdat,breaks=bk, col="black",xlab="aagmex")

dev.off()

############################################################################
# generate heatmap and 3d scatter for 3 m/f deg analyses
############################################################################

fname = "degHEAT2.pdf"
pdf(fname)

cw = 12
ch = 2

hmdat = all_fc %>% select(gtex_fc,decode_fc,aagmex_fc)
aheatmap(hmdat,breaks=c(0),Colv=NA,cellwidth=cw,cellheight=ch)

dev.off()

############################################################################
# generate venn diagram for 3 m/f deg analyses
############################################################################

gtex = nrow(gtex1 %>% filter(adj.P.Val<0.05,abs(ratioFC)>1.05))
decode = nrow(decode1 %>% filter(adj.P.Val<0.05,abs(ratioFC)>1.05))
aagmex = nrow(aagmex0 %>% filter(adj.P.Val<0.05,abs(ratioFC)>1.05))
gt_dc = nrow(gtex_decode_fc)
dc_aa = nrow(decode_aagmex_fc)
gt_aa = nrow(gtex_aagmex_fc)
gt_cd_aa = nrow(all_fc)

pdf("degVENN.pdf", height = 6, width = 6)
venn.plot <- draw.triple.venn(area1 = gtex, area2 = decode, area3 = aagmex,
                              n12 = gt_dc, n23 = dc_aa, n13 = gt_aa,
                              n123 = gt_cd_aa, category = c("GTEx", "deCODE", "AAGMEx"),
                              euler.d=F, scaled=F, overrideTriple=T)
dev.off()

############################################################################
# example gene expression
############################################################################


# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig1"
setwd(dir)

# Human DEGs for three data sets
dat = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/deg.RData"
load(dat)

library(ggplot2)
library(gridExtra)

# use gene ids for gtex
all(colnames(dat_F_gtex) == colnames(dat_M_gtex))
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/gencode_gene_map.txt"
ann_gene0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F,quote = "")
ggenes = sapply(colnames(dat_F_gtex),function(x){ann_gene0$gene_name[which(ann_gene0$gene_id==x)]})
colnames(dat_F_gtex) = colnames(dat_M_gtex) = ggenes

# use gene ids for decode
all(colnames(dat_F_decode) == colnames(dat_M_decode))
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/genes_decode.txt"
ann_gene0 = read.table(fname,header=T,sep="\t",stringsAsFactors=F,quote="")
ggenes = sapply(colnames(dat_F_decode),function(x){ann_gene0$gene_list[which(ann_gene0$ID==x)]})
colnames(dat_F_decode) = colnames(dat_M_decode) = ggenes

# all_fc
genes = c("FADS1","MYOT","NEO1" ,"DEGS1")

# female case 1 TFs
genes = c("FOXF2","FOXO1","MECOM","PBX1","STAT1","STAT5A")

# male case 1 TFs
genes = c("E2F3","EGR1","MYC","NRF1","SREBF2","TCF3","ZEB1" )

# eqtl deg overlap
genes = c("CCDC3","CLIC6","FADS1","GLDN","HSPA12A","MAP1B","MLPH","MMD","MYOT","NDRG4","NEO1","PDZD2","TBC1D9")

# plot for all genes
plts = list()
for(gene in genes){
  
  # zscore scaled data
  ind.gtex = which(colnames(dat_F_gtex) == gene)
  ind.dcde = which(colnames(dat_F_decode) == gene)
  ind.agmx = which(colnames(dat_F_aagmex) == gene)
  gtex.dat = c(dat_F_gtex[,ind.gtex], dat_M_gtex[,ind.gtex])
  dcde.dat = c(dat_F_decode[,ind.dcde], dat_M_decode[,ind.dcde])
  agmx.dat = c(dat_F_aagmex[,ind.agmx], dat_M_aagmex[,ind.agmx])
  gtex.sex = c(rep("F",length(dat_F_gtex[,ind.gtex])), rep("M",length(dat_M_gtex[,ind.gtex])))
  dcde.sex = c(rep("F",length(dat_F_decode[,ind.dcde])), rep("M",length(dat_M_decode[,ind.dcde])))
  agmx.sex = c(rep("F",length(dat_F_aagmex[,ind.agmx])), rep("M",length(dat_M_aagmex[,ind.agmx])))
  gtex = data.frame(expr=scale(gtex.dat), sex=gtex.sex, data="gtex")
  dcde = data.frame(expr=scale(dcde.dat), sex=dcde.sex, data="decode")
  agmx = data.frame(expr=scale(agmx.dat), sex=agmx.sex, data="aagmex")
  gene.plt.dat = rbind(gtex, dcde, agmx)
  
  # generate the plot
  plts[[gene]] = ggplot(data=gene.plt.dat, aes(x=data,y=expr,fill=sex)) + geom_boxplot() +
    scale_fill_manual(values=c("red", "blue")) + ylab(paste0(gene," z-score")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
}

pdf("fig7_13genes.pdf", onefile=TRUE, height=12, width=14)
marrangeGrob(grobs=plts, nrow=3, ncol=2, top=NULL)
dev.off()


############################################################
## integrate murine and human DEGs - generate Venn diagram
# see fig6.R
############################################################

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig1"
setwd(dir)

# Human DEGs for three data sets
dat = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/deg.RData"
load(dat)

load("fig1a.RData")

library(dplyr)
library(NMF)
library(VennDiagram)
library(biomaRt)

# HMDP DEG analysis results
fname = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig6/hmdpDEG.txt"
hmdp1 = read.table(fname,header=T,sep="\t",stringsAsFactors=F,quote = "")

## get human orthologs of mouse genes
## get overlapping lists of genes
genesH = all_fc$gene
genesM = hmdp1$gene
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesHM = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                 values = genesH, mart = mouse, attributesL = c("hgnc_symbol"), 
                 martL = human, uniqueRows=T)
orthologs = sapply(genesH,function(x){genesHM$MGI.symbol[which(genesHM$HGNC.symbol==x)]})

# filter orthologs to include only human genes with murine orthologs
# e.g., remove CPAMD8 which has no murine ortholog
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

# add human (gtex) and murine fold change information for orthologs
# this will filter the data based on expressed probes
mouse.ind = sapply(orthos$mouse,function(x){which(hmdp1$gene==x)}) %>% unlist
orthos = orthos[orthos$mouse %in% hmdp1$gene[mouse.ind],]
mouse.ind = sapply(orthos$mouse,function(x){which(hmdp1$gene==x)}) %>% unlist
human.ind = sapply(orthos$human,function(x){which(all_fc$gene==x)}) %>% unlist
orthofc = cbind(orthos, all_fc[human.ind,])
orthofc = orthofc %>% mutate(mouse_fc = hmdp1$ratioFC[mouse.ind])
orthofc = orthofc %>% mutate(mouse_fdr = hmdp1$adj.P.Val[mouse.ind])

# identify cases where the mouse fold changes match those for human data
# and |fc| > 1.05 and FDR < 0.05
# save(orthofc_hm, file="orthofc_hm.RData")
prodhm = orthofc$gtex_fc * orthofc$mouse_fc
mfdr = orthofc$mouse_fdr
ind1 = which(sign(prodhm) == 1 & abs(orthofc$mouse_fc) > 1.05)
ind2 = which(mfdr < 0.05)
inds = intersect(ind1, ind2)
orthofc_hm = orthofc[inds,]

# draw the Venn diagram
# a1 = intersect genes from degs of 3 data sets
# a2 = mouse genes with orthologs from a1
library(VennDiagram)
a1 = nrow(all_fc)
a2 = nrow(orthos)
ax = nrow(orthofc_hm)
pdf("degVENN_hm.pdf", height = 6, width = 6)
venn.plot <- draw.pairwise.venn(area1 = a1, area2 = a2, 
                              cross.area = ax,
                              category = c("human", "mouse"),
                              scaled=T)
dev.off()


############################################################
## integrate healthy human DEGs and obese DEGs - generate Venn diagram
## see vignette: obesity.pdf
## see R code: obesity.R
############################################################

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig1"
setwd(dir)

# Human DEGs for three data sets
dat = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/deg.RData"
load(dat)

load("fig1a.RData")
load("ob.RData")
ob0 = output

library(dplyr)
library(VennDiagram)

# retain the probe with the lowest p-value
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
ob0 = ob0[min.p.indices(ob0),]

# reduce the obesity set to match the 3data set 
ob_overlap = ob0[ob0$gene %in% all_fc$gene,]

# compare obese and healthy data
ind.health = sapply(ob_overlap$gene,function(x){which(all_fc$gene==x)})
combined = cbind(all_fc[ind.health,], ob_overlap)
fcprod = sign(combined$gtex_fc * combined$ratioFC)
ind1 = which(fcprod == 1 & abs(combined$ratioFC) > 1.05)
ind2 = which(combined$adj.P.Val < 0.05)
inds = intersect(ind1, ind2)
combined_filtered = combined[inds,]

# draw the Venn diagram
# a1 = intersect genes from degs of 3 data sets
# a2 = obese genes with matches from a1 (here a1 = a2 = 162)
library(VennDiagram)
a1 = nrow(all_fc)
a2 = nrow(combined)
ax = nrow(combined_filtered)
pdf("degVENN_ob.pdf", height = 6, width = 6)
venn.plot <- draw.pairwise.venn(area1 = a1, area2 = a2, 
                                cross.area = ax,
                                category = c("human", "mouse"),
                                scaled=T)
dev.off()


############################################################
## integrate healthy human DEGs and multi-tissue DEGs - generate Venn diagram
## see vignette: 
## see R code: tissue_deg.R
############################################################

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig1"
setwd(dir)

# Human DEGs for three data sets
dat = "/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/DEG/deg.RData"
load(dat)

# multi-tissue DEGs
load("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig1/tissues/tissueFOLD.RData")
load("/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig1/tissues/tissuePVAL.RData")

library(dplyr)
library(VennDiagram)

# params
fcthresh = 1.05
pvthresh = 0.05

# compare with all other tissues
ntissues = rep(0,nrow(all_fc))
names(ntissues) = all_fc$gene
for(ii in 1:nrow(all_fc)){
  fc1 = sign(all_fc$gtex_fc[ii])
  fc2 = tissueFOLD[ii,-1]
  pv2 = tissuePVAL[ii,-1]
  if(fc1 == 1){ind1 = which(fc2 > fcthresh)}
  if(fc1 == -1){ind1 = which(fc2 < -1*fcthresh)}
  ind2 = which(pv2 < pvthresh)
  ind = intersect(ind1,ind2)
  ntissues[ii] = length(ind)
}

# draw the Venn diagram
# a1 = intersect genes from degs of 3 data sets
# a2 = obese genes with matches from a1 (here a1 = a2 = 162)
library(VennDiagram)
a1 = nrow(all_fc)
a2 = nrow(all_fc)
ax = length(which(ntissues>0))
pdf("degVENN_tissues.pdf", height = 6, width = 6)
venn.plot <- draw.pairwise.venn(area1 = a1, area2 = a2, 
                                cross.area = ax,
                                category = c("subq", "tissues"),
                                scaled=T)
dev.off()


######################
# check just for visceral tissue

# compare with all other tissues
nvisc = rep(0,nrow(all_fc))
names(nvisc) = all_fc$gene
for(ii in 1:nrow(all_fc)){
  fc1 = sign(all_fc$gtex_fc[ii])
  fc2 = tissueFOLD[ii,2]
  pv2 = tissuePVAL[ii,2]
  if(fc1 == 1){ind1 = which(fc2 > fcthresh)}
  if(fc1 == -1){ind1 = which(fc2 < -1*fcthresh)}
  ind2 = which(pv2 < pvthresh)
  ind = intersect(ind1,ind2)
  nvisc[ii] = length(ind)
}

# draw the Venn diagram
# a1 = intersect genes from degs of 3 data sets
# a2 = obese genes with matches from a1 (here a1 = a2 = 162)
library(VennDiagram)
a1 = nrow(all_fc)
a2 = nrow(all_fc)
ax = length(which(nvisc>0))
pdf("degVENN_visc.pdf", height = 6, width = 6)
venn.plot <- draw.pairwise.venn(area1 = a1, area2 = a2, 
                                cross.area = ax,
                                category = c("subq", "visc"),
                                scaled=T)
dev.off()


############################################################
## aggregate differential expression information
## see vignette: 
## see R code: tissue_deg.R
############################################################

# orthofc
# ob_overlap

# main data
main = all_fc

# add obesity
ind = sapply(main$gene,function(x){which(ob_overlap$gene==x)})
main = main %>% mutate(mgh_fc = ob_overlap$ratioFC[ind], 
                       mgh_fdr = ob_overlap$adj.P.Val[ind])

# add gtex visc
ind = sapply(main$gene,function(x){which(tissueFOLD$gene==x)})
main = main %>% mutate(visc_fc = tissueFOLD$visceralAdipose[ind], 
                       visc_fdr = tissuePVAL$visceralAdipose[ind])

# add hmdp
miss = main$gene[!(main$gene %in% orthofc$human)]
fchmdp = c(orthofc$mouse_fc, rep(NA,length(miss)))
pvhmdp = c(orthofc$mouse_fdr, rep(NA,length(miss)))
names(fchmdp) = names(pvhmdp) = c(orthofc$human,miss)
ind = sapply(main$gene,function(x){which(names(fchmdp)==x)})
main = main %>% mutate(hmdp_fc = fchmdp[ind], 
                       hmdp_fdr = pvhmdp[ind])

# add the number of tissues
ind = sapply(main$gene,function(x){which(names(ntissues)==x)})
main = main %>% mutate(ntissues = ntissues[ind])

# save the data table
fname = "supplementaryTable1.txt"
write.table(main,fname,col.names=T,row.names=F,quote=F,sep="\t")


# hist of other tissue counts
fname = "ntissuesHIST.pdf"
pdf(fname)
par(mfrow=c(3,3)) 
hist(ntissues,col="black",xlab="n other tissues")
dev.off()


################################333
# checks
# save.image("fig1dat.RData")
# load("fig1dat.RData")
length(which(main$ntissues>0))
mainH = main %>% filter(hmdp_fc != "NA")
mainH = mainH %>% filter(abs(hmdp_fc)>1.05, hmdp_fdr<0.05)
sH = sign(mainH$gtex_fc * mainH$hmdp_fc)
length(which(sH==1))
mainD = main %>% filter(abs(mgh_fc) > 1.05, mgh_fdr < 0.05)
sD = sign(mainD$gtex_fc * mainD$mgh_fc)
length(which(sD==1))
mainV = main %>% filter(abs(visc_fc) > 1.05, visc_fdr < 0.05)
sV = sign(mainV$gtex_fc * mainV$visc_fc)
length(which(sV==1))
length(which(main$ntissues>0))


############################################################
## analyze suplementary information
############################################################

# save the data table
fname = "supplementaryTable1.txt"
table0 = read.table(fname,header=T,stringsAsFactors=F,sep="\t")

# identify healthy subcutaneous genes
table0_exclusive = table0 %>% filter(ntissues==0, abs(mgh_fc)<1.05) %>% mutate(specific="yes")
table0_inclusive = table0[!(table0$gene %in% table0_exclusive$gene),] %>% mutate(specific="no")
tableS1 = rbind(table0_exclusive,table0_inclusive)
write.table(tableS1,"TableS1.txt",col.names=T,row.names=F,quote=F,sep="\t")

# identify multiple tissue and obesity genes
table0_multi = table0 %>% filter(ntissues!=0, abs(mgh_fc)>1.05, mgh_fdr<0.05) 

# identify murine matches
table0_murine = table0 %>% filter(sign(decode_fc)==sign(hmdp_fc), abs(hmdp_fc)>1.05, hmdp_fdr<0.05)
table0_exclusive_mice = table0_exclusive %>% filter(sign(decode_fc)==sign(hmdp_fc), abs(hmdp_fc)>1.05, hmdp_fdr<0.05)
write.table(table0_exclusive_mice,"TableS2.txt",col.names=T,row.names=F,quote=F,sep="\t")




