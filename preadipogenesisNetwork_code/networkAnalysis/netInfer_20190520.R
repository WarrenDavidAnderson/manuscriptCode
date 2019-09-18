

library(dplyr)
library(DESeq2)
library(ggplot2)
library(bigWig)
library(stringr)
library(prodlim)

dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/network_initial/round3")
setwd(dir)


############################################################
## import data, specify data paths
############################################################

# load differential expression data, deg, see round2
load("LRTres_pro_20190501.RData")
deg.pro = out[["deg"]]
deseq.obj.pro = out[["deseq_obj"]]
load("LRTres_atac_20190501.RData")
deg.atac = out[["deg"]]
deseq.obj.atac = out[["deseq_obj"]]

# load atact peak annotation (motif.mapping, peak.ann.key)
# see peak_annotation.R
load("peak.annotation.20180423.RData")
load("peak.annotation.key.20180423.RData")
motif.mapping = motif.mapping[-which(rowSums(motif.mapping)==0),]

# load TF community gene information (tfclass.lists)
# see communityTFs.R
load("TFcommunity.gene.lists.RData")  

# TU annotation 
bed0 = read.table("TU_20181107.bed",stringsAsFactors=F,header=F)
names(bed0) = c("chr","start","end","gene","xy","strand")

# load atac peak coords and empirically defined summits
# see /ATAC_peaks/integratedReads.sh, empiricalSummits_20190325.R
# summits.bed = find.summits(peaks=bed.map[,1:3], window=50, bin=5, bw.file=loaded.bw)
# save(summits.bed, file="ATACsummits_20190325.RData")
load("ATACsummits_20190325.RData")
all.peaks = paste0(summits.bed[,1],":",summits.bed[,2],"-",summits.bed[,3])
all.summits = paste0(summits.bed[,1],":",summits.bed[,2],
                     "-",summits.bed[,3],",",summits.bed[,4])

# upload mapping between TF motif and database names
pathToFile <- paste0(getwd(),"/motif.id.key.txt")
f <- file(pathToFile, "rb")  
rawText <- readLines(f)
close(f)
motifs0 <- sapply(rawText, str_split, "\t")
motifs.raw = unlist(motifs0)[-c(1)] %>% unname
ids.new = sapply(motifs.raw,function(x){strsplit(x," ",fixed=TRUE)[[1]][1]}) %>% unname
ids.orig = sapply(motifs.raw,function(x){
  id = strsplit(x," ",fixed=TRUE)[[1]]
  out = paste(id[2:length(id)], collapse=" ")
}) %>% unlist
id.map = as.data.frame(cbind(ids.new,ids.orig),stringsAsFactors=F)
rownames(id.map) = 1:nrow(id.map)

# atac bigWig directory
bigWig.path.atac = paste0("/media/wa3j/Seagate2/Documents/PRO/adipogenesis",
                          "/July2018/atac_bw_20181127/atac_all")

# motif coordinates, data location for fimo bigWigs
# bigwigs produced in TFfimo216.sh
fimo.bed.path = paste0("/media/wa3j/Seagate2/Documents/PRO/adipogenesis",
                       "/July2018/atac_time_deg/motifEnrich/fimo_bed750k/")

# bedtools
bed.dir = "/media/wa3j/Seagate2/Documents/software/bedtools2/bin/"

############################################################
## key analysis parameters
############################################################

# sig thresh for LRT
sig.thresh = 0.05
fc.thresh = 1

# search limit for re -> tu edge weights (kb)
atac.lim = 100

# exponential decay at half max for re -> tu eage weights (kb)
atac.half = 10

# rate constant for exponential decay model of the regulatory weight (bp)
xhalf = 10 # kb
tau = -1 / (log(1/3) / (xhalf*1000))

# network times and corresponding deseq condition factor ids
tnet = c("t0", "60m", "2h", "3h", "4h")
tcnd = c(0, 1, 2, 3, 4)

############################################################
## TF annotation
############################################################

# community - gene annotation, select the motif with....
tfclass.ann = c()
for(ii in 1:length(tfclass.lists)){
  new = cbind(names(tfclass.lists)[ii], tfclass.lists[[ii]])
  tfclass.ann = rbind(tfclass.ann, new)
}
tfclass.ann = as.data.frame(tfclass.ann, stringsAsFactors=F)
names(tfclass.ann) = c("comm","tf")

# for every main community, select the one with the most fimo sites
# see peak.ann.key in peak_annotation.R
# alternatively, remove communities 61 and 267
ind.rem1 = which(tfclass.ann$comm == "community_61")
ind.rem2 = which(tfclass.ann$comm == "community_267")
tfclass.ann = tfclass.ann[-c(ind.rem1,ind.rem2),]

# add representative community TF identifier
tfclass.ann = tfclass.ann %>% mutate(TFid=NA, motifid=NA)
for(ii in 1:nrow(tfclass.ann)){
  ind = which(peak.ann.key$maincom == tfclass.ann$comm[ii])
  tfclass.ann$TFid[ii] = peak.ann.key$TF[ind]
  tfclass.ann$motifid[ii] = peak.ann.key$motif[ind]
}


############################################################
## implement network analysis for a list with multiple times
## initial exploratory analysis
############################################################

# specify the comparison lists tnet and tcnd
t0 = c("t0")
t1 = c("20m", "40m")
t2 = c("3h", "4h")
tnet = list(t0=t0, t1=t1, t2=t2)
t0 = c(0)
t1 = c(0.33, 0.67)
t2 = c(3, 4)
tcnd = list(t0=t0, t1=t1, t2=t2)

net2 = infer.net.list(deg.pro=deg.pro, deseq.obj.pro=deseq.obj.pro,
                      deg.atac=deg.atac, deseq.obj.atac=deseq.obj.atac,
                      sig.thresh=sig.thresh, fc.thresh=fc.thresh, 
                      bigWig.path.atac=bigWig.path.atac, 
                      fimo.bed.path=fimo.bed.path, bed.dir=bed.dir,
                      tfclass.ann=tfclass.ann, bed0=bed0, 
                      tnet=tnet, tcnd=tcnd, all.summits=all.summits,
                      atac.lim=100, atac.half=10, full=FALSE,
                      motif.mapping=motif.mapping, dist1=50, dist2=300)

# save(net2, file="net2.RData")
# load("net2.RData")

net2$conn
net2$contin
net2$cnts
net = net2$net

check = force.connect(network=net, lim=10000)
dim(check$net)
dim(net)

# filter edges and check for connectivity
hist(net$weight, col="black")
quantile(abs(net$weight))
mean( abs(net$weight) ) + sd( abs(net$weight) )

# filter
filt0 = net %>% filter(abs(weight) > 0.9)
res.filt0 = force.connect(filt0)
res.filt0$cnts
hist(res.filt0$net$weight, col="black", xlab="edge weight")

# check2
ind1 = which(net$tt=="0_0.33,0.67" & net$type=="REtoTF")
ind2 = which(net$tt=="0.33,0.67_3,4" & net$type=="TFtoRE")
dat1 = net[ind1,]
dat2 = net[ind2,]
all(dat1$geneTU %in% dat2$geneTU)
all(dat2$geneTU %in% dat1$geneTU)
ind1 = which(net$tt=="0.33,0.67_3,4" & net$type=="REtoTF")
ind2 = which(net$tt=="0.33,0.67_3,4" & net$type=="TFtoRE")
dat1 = net[ind1,]
dat2 = net[ind2,]
all(dat1$idRE %in% dat2$idRE)
all(dat2$idRE %in% dat1$idRE)



############################################################
## format data for visualization in cytoscape
############################################################

# filter
q = quantile( abs(net$weight) )
filt0 = net %>% filter(abs(weight) > q[4])
res.filt0 = force.connect(filt0)
res.filt0$cnts
net.filt = res.filt0$net

save(net.filt, file="net.filt.RData")

# check
re = "chr5:32135323-32136049"
ind = grep(re, net.filt$idRE)
net.filt[ind,]

nodes[which(nodes$node==re),]

#####################
# cytoscape edges

# re -> tf
ind = which(net.filt$type=="REtoTF")
edges1 = net.filt[ind,] %>% select(idRE,geneTU,weight,tt,type,tu)

# tf -> re
ind = which(net.filt$type=="TFtoRE")
edges2 = net.filt[ind,] %>% select(geneTU,idRE,weight,tt,type,tu)

# all edges
names(edges1) = names(edges2) = c("source","target","weight","time","type","tu")
edges = rbind(edges1, edges2)
ind1 = grep("tu_class",edges$target)
ind2 = grep("unann_",edges$target)
ind3 = grep("Rik",edges$target)
ind = Reduce(union, list(ind1,ind2,ind3))
edges$target[ind] = paste0("unann",c(1:length(ind)))
edges = edges %>% mutate(sign = sign(weight))

#####################
# cytoscape nodes

# classify nodes
unique(edges$time)
ind = which(edges$time=="0_0.33,0.67" & edges$type=="REtoTF")
node1 = edges[ind,] %>% select(source,time,type,tu) %>% mutate(class="earlyRE")
node2 = edges[ind,] %>% select(target,time,type,tu) %>% mutate(class="earlyTF")
ind = which(edges$time=="0.33,0.67_3,4" & edges$type=="REtoTF")
node3 = edges[ind,] %>% select(source,time,type,tu) %>% mutate(class="lateRE")
node4 = edges[ind,] %>% select(target,time,type,tu) %>% mutate(class="lateTF")
names(node1)=names(node2)=names(node3)=names(node4) = c("node","time","type","tu","class")
node = rbind(node1,node2,node3,node4)

# check
unique(node$class)
allnode = unique(node$node)
alledge = unique(c(edges$source,edges$target))
all(allnode %in% alledge)
all(alledge %in% allnode)

# tus
nodes1 = node %>% filter(class=="earlyTF") %>% unique
nodes2 = node %>% filter(class=="lateTF") %>% unique
nodes0 = rbind(nodes1,nodes2) %>% unique
dups = nodes0$node[ which(duplicated(nodes0$node)==T) ] %>% unique
for(ii in 1:length(dups)){
  ind = which(nodes0$node == dups[ii])
  dat = nodes0[ind,]
  tt = paste0(dat$class,collapse="--")
  dat$class = tt
  nodes0[nodes0$node==dat$node[1],] = dat[1,]
}
nodesTF = unique(nodes0) %>% mutate(ann="tftu")

# res
nodes1 = node %>% filter(class=="earlyRE") %>% unique
nodes2 = node %>% filter(class=="lateRE") %>% unique
nodes0 = rbind(nodes1,nodes2) %>% unique
dups = nodes0$node[ which(duplicated(nodes0$node)==T) ] %>% unique
# for(ii in 1:length(dups)){
#   ind = which(nodes0$node == dups[ii])
#   dat = nodes0[ind,]
#   tt = paste0(dat$time,collapse="--")
#   dat$time = tt
#   nodes0[nodes0$node==dat$node[1],] = dat[1,]
# }
nodesRE = unique(nodes0) %>% mutate(ann="re")
nodesRE$tu = "na"

# all nodes
nodes = rbind(nodesTF, nodesRE)

# output for cytoscape
write.table(edges,"cyto_edges.txt",sep="\t",quote=F,col.names=T,row.names=F)
write.table(nodes,"cyto_nodes.txt",sep="\t",quote=F,col.names=T,row.names=F)



############################################################
## plot re/tf dynamics
############################################################

library(graphics)
# load("net2.RData")

# basic data
t1 = c("t0","20m","40m","60m","2h","3h","4h")
t2 = c(0,0.33,0.67,1,2,3,4)
tt = data.frame(cond=t1, time=t2, stringsAsFactors=F)
atac.reads = counts(deseq.obj.atac)
atac.sf = sizeFactors(deseq.obj.atac)
atac.norm = sf.norm(cnts=atac.reads, sf=atac.sf, ann=tt)
pro.reads = counts(deseq.obj.pro)
pro.sf = sizeFactors(deseq.obj.pro)
pro.norm = sf.norm(cnts=pro.reads, sf=pro.sf, ann=tt)

# earlt TFs
TFi = nodes %>% filter(class=="earlyTF" | class=="earlyTF--lateTF") %>% select(node)


pdf("tfre.pdf"); par(mfrow=c(3,3))
for(ii in 1:nrow(TFi)){
  
  # targets of early TFs
  tf = TFi[ii,1]
  net = edges %>% filter(source==tf)
  net.up = net %>% filter(sign==1)
  net.dn = net %>% filter(sign==-1)
  rds.up = atac.norm[rownames(atac.norm) %in% net.up$target,]
  rds.dn = atac.norm[rownames(atac.norm) %in% net.dn$target,]
  
  rds.up.log = log2(rds.up)
  rds.dn.log = log2(rds.dn)
  
  curves.up = read.curve.smooth(dat=rds.up.log)
  curves.dn = read.curve.smooth(dat=rds.dn.log)
  
  tf.dat = pro.norm[grep(tf,rownames(pro.norm)),] 
  tf.dat = log2(tf.dat)
  tf.dat = read.curve.smooth(dat=tf.dat)
  tf.dat = plot.curves(curves=tf.dat, id=tf, plt=FALSE)
  
  # plot up/down traces
  up.av = plot.curves(curves=curves.up, id=tf, plt=TRUE, col="red")
  dn.av = plot.curves(curves=curves.dn, id=tf, plt=TRUE, col="blue")
  
  # plot means with TF data
  tax = up.av$time 
  tf.curve = scale01(tf.dat$mean)
  up.curve = scale01(up.av$mean)
  dn.curve = scale01(dn.av$mean)
  plot(tax, tf.curve, type="l", xlab="time (hr)", ylab="scaled signal", main=tf, lwd=2)
  lines(tax, up.curve, type="l", col="red", lwd=2)
  lines(tax, dn.curve, type="l", col="blue", lwd=2)
  
}
dev.off()




# targets of early TFs
tf = "Fosl2"
net = edges %>% filter(source==tf)
net.up = net %>% filter(sign==1)
net.dn = net %>% filter(sign==-1)
rds.up = atac.norm[rownames(atac.norm) %in% net.up$target,]
rds.dn = atac.norm[rownames(atac.norm) %in% net.dn$target,]

rds.up.log = log2(rds.up)
rds.dn.log = log2(rds.dn)

curves.up = read.curve.smooth(dat=rds.up.log)
curves.dn = read.curve.smooth(dat=rds.dn.log)

tf.dat = pro.norm[grep(tf,rownames(pro.norm)),] 
tf.dat = log2(tf.dat)
tf.dat = read.curve.smooth(dat=tf.dat)
tf.dat = plot.curves(curves=tf.dat, id=tf, plt=FALSE)

# plot up/down traces
up.av = plot.curves(curves=curves.up, id=tf, plt=TRUE, col="red")
dn.av = plot.curves(curves=curves.dn, id=tf, plt=TRUE, col="blue")

# plot means with TF data
tax = up.av$time 
tf.curve = scale01(tf.dat$mean)
up.curve = scale01(up.av$mean)
dn.curve = scale01(dn.av$mean)
plot(tax, tf.curve, type="l", xlab="time (hr)", ylab="scaled signal", main=tf, lwd=2)
lines(tax, up.curve, type="l", col="red", lwd=2)
lines(tax, dn.curve, type="l", col="blue", lwd=2)































