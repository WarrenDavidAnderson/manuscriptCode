
fndir = paste0("/run/user/1001/gvfs/smb-share:server=home1.virginia.edu,",
  "share=cphgdesk/users/wa3j/MS/adipogenesis/vignettes/ATAC_analysis/motifMAP/",
  "motif_clust_functions.R")
source(fndir)

dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
  "adipogenesis/July2018/atac_time_deg/communities/TF")
setwd(dir)

library(NMF)
library(dplyr)
library(igraph)
library(stringr)
library(randomcoloR)
library(RColorBrewer)

##############################################################
## import data and generate three column graph format
##############################################################

# basic tomtom data and pval filter
fname = "tomtomTFvsTF.txt"
dat0 = read.table(fname,header=T,stringsAsFactors=F)
facs = unique(dat0$Query_ID)
length(facs) # 735
dat0 = dat0[dat0$Target_ID %in% facs,]
length(unique(dat0$Target_ID)) # 735
dat0 = dat0 %>% filter(p.value<0.001)

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

# directory for database TF pwms
pwm.dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
                 "adipogenesis/July2018/atac_time_deg/communities/TF/TF735")

##############################################################
## generate three column graph format
##############################################################

# remove self to self edges
ind.rem = which(dat0$Query_ID == dat0$Target_ID)
dat0 = dat0[-ind.rem,]

# generate edge data with -log10(pval)
threecol = dat0 %>% select(Query_ID, Target_ID, p.value)
names(threecol) = c('from','to','pval')
minp = threecol %>% filter(pval>0) %>% select(pval) %>% min
threecol$pval[threecol$pval == 0] = minp
weight = -log10(threecol$pval)
threecol$weight = weight

# filter weights
threecol = threecol %>% filter(weight>6)

# counts
tomtom.uniq = unique(c(threecol$from, threecol$to))
length(tomtom.uniq)
matched = tomtom.uniq
unmatched = facs[!(facs %in% matched)]

# weight distribution
hist(threecol$weight)

##############################################################
## create the graph format and detect communities
##############################################################

# create the graph
g = graph.data.frame(threecol, directed=F)
g = simplify(g)

# commnity detection
comm.g = fastgreedy.community(g)

# community annotation
ann.comm = sapply(unique(comm.g$membership),function(x){length(which(comm.g$membership==x))})
ann.comm = cbind(unique(comm.g$membership),ann.comm) %>% as.data.frame(stringsAsFactors=F)
names(ann.comm) = c("comm","count")
ann.comm = ann.comm[order(ann.comm$count,decreasing=T),]


##############################################################
## subcluster large communities 
##############################################################

# generate subclusters
subclust.dat = subcluster(clusters=comm.g, graph=threecol,
                  max.clust=10, max.iter=10, change.thresh=3)


##############################################################
## plot graph of all communities (original clustering)
##############################################################

# plot with motif colors
cols = distinctColorPalette(length(unique(ann.comm$comm)))

# comunity cluster plot function
comm.plt = function(g=NULL,layout=NULL,vertex.color=NULL,
                    edge.width=NULL,vertex.size=NULL){
  #pdf("clust.pdf")
  plot(g,layout=layout, rescale=F, vertex.label.cex=.5,
       xlim=range(layout[,1]), ylim=range(layout[,2]),
       edge.width=edge.width, vertex.size=vertex.size,
       edge.curved=T, vertex.label=NA, 
       vertex.color=vertex.color, edge.color="black",
       margin=0, asp=0)
  #dev.off()
}

# annotate community colors
vertex.color = cols[comm.g$membership]
names(vertex.color) = comm.g$membership

# layout.fruchterman.reingold
layout = layout.fruchterman.reingold(g)
layout = layout.norm(layout,-1,1,-1,1)
vertex.size = (degree(g,mode='out') + 10)/3
edge.width = E(g)$weight/4

pdf("clusters.pdf", width=6); par(mfrow=c(2,1))
comm.plt(g, layout, vertex.color, edge.width, vertex.size)
dev.off()

# annotate communities 
length(comm.g$membership)
length(comm.g$names)
comm.ann = list()
for(comm in ann.comm$comm){
  ind.comm = which(comm.g$membership == comm)
  tfs = comm.g$names[ind.comm]
  color = vertex.color[ind.comm[1]]
  comm.ann[[paste0("community_",comm)]] = list(color=color, tfs=tfs)
}


##############################################################
## output community annotation
## incorporate subclusters and unmatched motifs
##############################################################

# include un-matched motifs as additional clusters
motif.community = subclust.dat
ncom = max(subclust.dat$clust)
newcomms = c( (ncom+1) : (ncom+length(unmatched)) )
new = cbind(newcomms, unmatched) %>% as.data.frame(stringsAsFactors=F)
names(new) = names(motif.community)
motif.community = rbind(motif.community, new)
motif.community$clust = data.matrix(motif.community$clust) %>% as.numeric

# add original TF ids
inds = sapply(motif.community$motif,function(x){which(id.map$ids.new==x)})
ids.orig = id.map$ids.orig[inds]
motif.community = motif.community %>% mutate(ids.orig = ids.orig)

# checks
length(unique(motif.community$clust))
length(unique(motif.community$motif))

# generate motif list
motif.clusters = list()
for(ii in 1:length(unique(motif.community$clust))){
  clust = unique(motif.community$clust)[ii]
  ind = which(motif.community$clust == clust)
  namen = paste0("community_",clust)
  motif.clusters[[namen]] = motif.community$motif[ind]
}

# generate cluster TF 'spreadsheet' format
maxlen = lapply(motif.clusters,function(x){length(x)}) %>% unlist %>% max
datTF = as.data.frame(matrix(0,maxlen,length(motif.clusters)))
dat = as.data.frame(matrix(0,maxlen,length(motif.clusters)))
for(ii in 1:ncol(dat)){
  len = length(motif.clusters[[ii]])
  dat[1:len,ii] = motif.clusters[[ii]]
  inds = sapply(motif.clusters[[ii]],function(x){which(id.map$ids.new==x)})
  datTF[1:len,ii] = id.map$ids.orig[inds]
}
names(dat) = names(datTF) = names(motif.clusters)
write.table(dat,"TFclusters.txt",col.names=T,row.names=F,sep="\t")


##############################################################
## find avaerage motif PWM for each cluster
##############################################################

require(ggplot2)
require(ggseqlogo)
library(gridExtra)
library(dplyr)
library(seqLogo)
library(prodlim)
library(lattice)

# find the most connected motif for each community based on connection sum
names(motif.community)[1] = "community"
max.connect.motifs = get.max.connect.motifs2(motif.community=motif.community, graph=threecol)

# average motif analysis
site.thr = 0.1 # threshold for keeping a margin site in the pwm
pad.len = 60
max.offset = 10
bkg = c(0.25, 0.25, 0.25, 0.25)
pltdat = denovo.mean.plot.TF(dat=dat, pwm.dir=pwm.dir, tomtom=dat0, 
                max.connect.motifs=max.connect.motifs,bkg=bkg, 
                pad.len=pad.len,max.offset=max.offset)

# plot the motif information
pdf("TFsubclust10Filt1e6.pdf", onefile = TRUE, height=9, width=4.5)
marrangeGrob(grobs=pltdat$plt, nrow=6, ncol=1, top=NULL)
dev.off()


##############################################################
## filter for low quality TFs
##############################################################

# params
info.thresh = 1
n.info.match.sites = 3
nsite.in.row = 4
n.sites.info = 1.6
n.sites.info.thresh = 12

# loop through all clusters/motifs
bad.plt = list()
bad.nme = list()
flt.plt = list()
flt.pwm = list()
bad.facs = c()
m = n = p = 1
for(ii in 1:ncol(dat)){
  
  # get max connected motif for each community
  comm = strsplit(colnames(dat)[ii], "community_")[[1]][2]
  center = max.connect.motifs$motif[max.connect.motifs$community == comm]
  id = motif.community$ids.orig[which(motif.community$motif == center)]
  
  # get the center pwm
  fname = paste0(pwm.dir,"/",center,".txt")
  mat = read.table(fname,sep="\t",header=F,stringsAsFactors=F,fill=T) %>% t
  nas = apply(mat,1,function(x){all(is.na(x))})
  ind = which(nas==TRUE)
  if(length(ind)>0){mat = mat[-ind,]}
  rownames(mat) = c("A","C","G","T")
  colnames(mat) = c(1:ncol(mat))
  
  # get the corresponding average and align to the center
  av0 = pltdat$pwm.av[[center]]
  offset = find.offset(query=mat, target=av0, max.offset=max.offset, bkg=bkg)
  orient = offset$orient
  offset = offset$offset.dist
  if(orient == "rc"){av0 = revcomp(av0)}
  padded = get.padded.pwms(off=offset, query=mat, target=av0, bkg=bkg)
  center.pwm = padded$query.pad
  mean.pwm = padded$target.pad
  
  ###### implement filters #####
  # summary data and filter information
  mean.pwm = add.pseudocount(PWM=mean.pwm, cnt=1e-6)
  center.pwm = add.pseudocount(PWM=center.pwm, cnt=1e-6)
  site.info.tf = apply(center.pwm,2,function(x){x %*% (log2(x)-log2(bkg))})
  site.info.av = apply(mean.pwm,2,function(x){x %*% (log2(x)-log2(bkg))})
  ind.more.tf = which(site.info.tf > info.thresh) %>% unname
  ind.more.av = which(site.info.av > info.thresh) %>% unname
  ind.more.tf2 = which(site.info.tf > n.sites.info) %>% unname
  
  # filter by matching high info sites
  hi = length(intersect(ind.more.tf,ind.more.av))

  # n letters in a row in the denovo motif (output = 0)
  samen.cls = samen(mat=mean.pwm, n=nsite.in.row)
  
  # use TDdb motif name
  motifname = motif.community$ids.orig[motif.community$motif==center]
  
  if(hi < n.info.match.sites){
    reason = paste0("less than ",n.info.match.sites," matching sites with info >",info.thresh)
    bad.plt[[n]] = logo.plt(pwm=center.pwm, title=reason)
    bad.nme[[m]] = logo.plt(pwm=center.pwm, title=motifname)
    bad.facs = c(bad.facs, center); n=n+1; m=m+1
    next
  }
  
  if(samen.cls == 0){
    reason = paste0(nsite.in.row, " identical bases in a row")
    bad.plt[[n]] = logo.plt(pwm=center.pwm, title=reason)
    bad.nme[[m]] = logo.plt(pwm=center.pwm, title=motifname)
    bad.facs = c(bad.facs, center); n=n+1; m=m+1
    next
  }
  
  if(length(ind.more.tf2) > n.sites.info.thresh){
    reason = paste0("more than ",n.sites.info.thresh," sites with info >",n.sites.info)
    bad.plt[[n]] = logo.plt(pwm=center.pwm, title=reason)
    bad.nme[[m]] = logo.plt(pwm=center.pwm, title=motifname)
    bad.facs = c(bad.facs, center); n=n+1; m=m+1
    next
  }
  
  # plot and store good pwm averages
  out = remove.lowinfo.margin(pwm.plt=mean.pwm, site.thr=site.thr, bkg=bkg, cnt=1e-6)
  flt.plt[[p]] = logo.plt(pwm=out, title=motifname); p=p+1
  flt.pwm[[center]] = out
  
}

# plot the motif information
pdf("TF_filteredin.pdf", onefile = TRUE, height=9, width=4.5)
marrangeGrob(grobs=flt.plt, nrow=6, ncol=1, top=NULL)
dev.off()

# plot the motif information
pdf("TF_filteredout.pdf", onefile = TRUE, height=9, width=4.5)
marrangeGrob(grobs=bad.plt, nrow=6, ncol=1, top=NULL)
dev.off()

# plot the motif information
pdf("TF_filteredoutNames.pdf", onefile = TRUE, height=9, width=4.5)
marrangeGrob(grobs=bad.nme, nrow=6, ncol=1, top=NULL)
dev.off()

##############################################################
## output PWMs for tomtom 
##############################################################

length(flt.pwm) # 378

pwm.out = flt.pwm
cdir = getwd()
pwm.dir = paste0(cdir,"/center_av_pwm")
setwd(pwm.dir)

for(ii in 1:length(pwm.out)){
  tf = names(pwm.out)[ii]
  mat = t( pwm.out[[ii]] )
  fname = paste0(tf,".txt")
  write.table(mat,fname,sep="\t",quote=F,col.names=F,row.names=F)
}

setwd(cdir)

# setwd("/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/communities/TF")
# save.image("TFclusteravs.RData")

# analyze using tomtom
# each against all
# /TFcluster/TFaverage_eachVersusAll.sh


##############################################################
## cluster average TFs based on tomtom results
##############################################################

setwd("/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/communities/TF")
# load("TFclusteravs.RData")

library(NMF)
library(dplyr)
library(igraph)
library(stringr)
library(randomcoloR)
library(RColorBrewer)

# import tomtom data and pval filter
fname = "tomtomTFAVvsTFAV.txt"
tomtomAv0 = read.table(fname,header=T,stringsAsFactors=F)
tomtomAv0 = tomtomAv0[duplicated(tomtomAv0)==F,]
tomtomAv0 = tomtomAv0 %>% filter(p.value<0.001)
length(unique(tomtomAv0$Query_ID)) # 374
length(unique(tomtomAv0$Target_ID)) # 374

# remove self to self edges
ind.rem = which(tomtomAv0$Query_ID == tomtomAv0$Target_ID)
tomtomAv0 = tomtomAv0[-ind.rem,]

# generate edge data with -log10(pval)
threecolAv = tomtomAv0 %>% select(Query_ID, Target_ID, p.value)
names(threecolAv) = c('from','to','pval')
minp = threecolAv %>% filter(pval>0) %>% select(pval) %>% min
threecolAv$pval[threecolAv$pval == 0] = minp
weight = -log10(threecolAv$pval)
threecolAv$weight = weight

# filter weights
threecolAv = threecolAv %>% filter(weight>6)

# create the graph
gAv = graph.data.frame(threecolAv, directed=F)
gAv = simplify(gAv)

# commnity detection
comm.gAv = fastgreedy.community(gAv)

# annotate community colors
cols = distinctColorPalette(length(unique(comm.gAv$membership)))
vertex.color = cols[comm.gAv$membership]
names(vertex.color) = comm.gAv$membership

# layout.fruchterman.reingold
layout = layout.fruchterman.reingold(gAv)
layout = layout.norm(layout,-1,1,-1,1)
vertex.size = (degree(gAv,mode='out') + 10)/3
edge.width = E(gAv)$weight/4
comm.plt(g=gAv, layout, vertex.color, edge.width, vertex.size)


##############################################################
## for each cluster, pick the average with the highest information
## integrate with unclustered TFs to get the final set of TFs
##############################################################

require(ggplot2)
require(ggseqlogo)
library(gridExtra)
library(dplyr)
library(seqLogo)
library(prodlim)
library(lattice)

aa.flt.pwm = names(flt.pwm)
clustered = comm.gAv$names # 17
unclusted = aa.flt.pwm[!(aa.flt.pwm %in% clustered)] # 357
length(aa.flt.pwm) == length(clustered) + length(unclusted)
newclust = comm.gAv$membership %>% unique # 7
nclust = length(newclust)

# loop through clusters, pick the element with highest info
clustTFs = c()
for(ii in newclust){
  
  # get cluster motifs
  ind = which(comm.gAv$membership == ii)
  motifs = comm.gAv$names[ind]
  
  # get the cluster pwms
  clust.pwm = list()
  for(jj in motifs){
    fname = paste0(pwm.dir,"/",jj,".txt")
    mat = read.table(fname,sep="\t",header=F,stringsAsFactors=F,fill=T) %>% t
    nas = apply(mat,1,function(x){all(is.na(x))})
    ind = which(nas==TRUE)
    if(length(ind)>0){mat = mat[-ind,]}
    rownames(mat) = c("A","C","G","T")
    colnames(mat) = c(1:ncol(mat))
    clust.pwm[[jj]] = mat
  }
  
  # get summed info for each pwm
  infos = lapply(clust.pwm,function(x){
    pwm = add.pseudocount(PWM=x, cnt=1e-6)
    info.tf = apply(pwm,2,function(x){x %*% (log2(x)-log2(bkg))})
    return(sum(info.tf)) 
  }) %>% unlist
  
  out = names(infos)[infos == max(infos)]
  clustTFs = c(clustTFs, out)
}

# aggregate filtered/clustered TF pwms for plotting
length(flt.pwm) # 378
length(clustered) # 22
length(unclusted) # 356
length(clustTFs) # 9
flt.pwm.clust = list()
for(ii in unclusted){
  flt.pwm.clust[[ii]]=flt.pwm[[ii]]
}
for(ii in clustTFs){
  flt.pwm.clust[[ii]]=flt.pwm[[ii]]
}
length(flt.pwm.clust) # 365

# plot filtered/clustered TFs
filt.clust.plts = list()
p=1
for(ii in names(flt.pwm.clust)){
  id = motif.community$ids.orig[motif.community$motif==ii]
  mat = flt.pwm.clust[[ii]]
  filt.clust.plts[[p]] = logo.plt(pwm=mat, title=id)
  p=p+1
}
pdf("TF_filteredClust365.pdf", onefile = TRUE, height=9, width=4.5)
marrangeGrob(grobs=filt.clust.plts, nrow=6, ncol=1, top=NULL)
dev.off()

##############################################################
## filter TF pwm set based on number of de novo motifs matching 
##############################################################

# import de novo motif ids (E < 0.1)
# /nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motifdb
# from=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motifdb/denovo2775.txt
# to=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/communities/TF
# scp -r wa3j@interactive.hpc.virginia.edu:${from} ${to}
fname = "denovo2775.txt"
denovo0 = read.table(fname,header=F,stringsAsFactors=F,sep="\t")
names(denovo0) = c("denovoid")

# get the number of comparisons for each de novo motif
shortid = sapply(denovo0$denovoid,function(x){strsplit(x,"_")[[1]][1]})
denovo1 = denovo0 %>% mutate(shortid=shortid)
denovo2 = c()
for(id in unique(shortid)){
  ind = which(denovo1$shortid==id)
  cond = sapply(denovo1$denovoid[ind],function(x){strsplit(x,"_")})
  cond = lapply(cond,function(x){paste(x[c(2,5,6)],collapse="_")}) %>% unlist
  ncomps = length(cond)
  comps = paste(cond,collapse=",")
  new = data.frame(motif=id, comps=comps, ncomps=ncomps, stringsAsFactors=F)
  denovo2 = rbind(denovo2, new)
}

# get mapping between de novo motifs and database TFs
# see tomtom_denovo_vs_db_array.sh
# from=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/tomtomDENOVOvsTF/TFfromDENOVOfiltered.txt
# to=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/communities/TF
# scp -r wa3j@interactive.hpc.virginia.edu:${from} ${to}
fname = "TFfromDENOVOfiltered.txt"
tomtom_DETF = read.table(fname,header=F,stringsAsFactors=F,sep="\t")
names(tomtom_DETF) = c("Query_ID",	"Target_ID",	"Optimal_offset",	"pvalue",	"Evalue",
            "qvalue",	"Overlap",	"Query_consensus",	"Target_consensus",	"Orientation")
tomtom_DETF = tomtom_DETF %>% filter(pvalue < 0.001)
consensus = sapply(tomtom_DETF$Query_ID,function(x){strsplit(x,"_")[[1]][1]})
tomtom_DETF = tomtom_DETF %>% mutate(consensus=consensus)
denovo.flt = unique(tomtom_DETF$Query_ID) # 2515

# loop through each de novo motif and get the top TF match
filtered.denovo = c()
for(ii in 1:length(denovo.flt)){
  mot = denovo.flt[ii]
  ttdat = tomtom_DETF[tomtom_DETF$Query_ID == mot,] 
  ttdat = ttdat[with(ttdat, order(pvalue,Target_ID)),]
  filtered.denovo = rbind(filtered.denovo, ttdat[1,])
}
length(unique(filtered.denovo$Query_ID)) # 2515
length(unique(filtered.denovo$Target_ID)) # 735

# initial annotation for TF clusters
# motif.community
motif.community0 = motif.community

# remove filtered TF clusters
motif.community1 = motif.community0
for(ii in bad.facs){
  comm = motif.community1$community[motif.community1$motif==ii]
  ind = which(motif.community1$community == comm)
  motif.community1 = motif.community1[-ind,]
}
length(motif.community1$community %>% unique) # 378

# remove all of the unclustered clusters 
motif.community2 = motif.community1
clustclust = c()
for(ii in clustered){
  comm = motif.community2$community[which(motif.community2$motif==ii)]
  ind = which(motif.community2$community == comm)
  motif.community2 = motif.community2[-ind,]
}
length(motif.community2$community %>% unique) # 356

# add back cluster representatives
motif.community3 = motif.community2
for(comm in newclust){
  
  # get new cluster motifs
  ind = which(comm.gAv$membership == comm)
  motifs = comm.gAv$names[ind]
  
  # aggregate old clusters
  oldcoms = sapply(motifs,function(x){motif.community1$community[motif.community1$motif==x]})
  indcoms = sapply(oldcoms,function(x){which(motif.community1$community==x)}) %>% unlist
  newcom = motif.community1[indcoms,]
  newcom$community = min(oldcoms)
  motif.community3 = rbind(motif.community3, newcom)
}
length(motif.community3$community %>% unique) # 365

#### get the num of de novo hits for each cluster ###
motif.community.summary = c()
for(ii in 1:length(unique(motif.community3$community))){
  comm = unique(motif.community3$community)[ii]
  tf.motifs = motif.community3$motif[motif.community3$community == comm]
  ind.de = sapply(tf.motifs,function(x){which(filtered.denovo$Target_ID==x)}) %>% unlist
  de.motifs = filtered.denovo$consensus[ind.de] %>% unique
  comps = sapply(de.motifs,function(x){
    all = denovo2$comps[denovo2$motif==x]
    out = strsplit(all,",")[[1]]
    return(out)
  }) %>% unlist %>% unique
  compout = paste(comps,collapse=",")
  center = names(flt.pwm.clust)[names(flt.pwm.clust) %in% tf.motifs]
  id = motif.community3$ids.orig[motif.community3$motif == center]
  new = c(center, id, length(comps), compout)
  motif.community.summary = rbind(motif.community.summary, new)
}
motif.community.summary = as.data.frame(motif.community.summary, stringsAsFactors=F)
names(motif.community.summary) = c("id","tf","ncomps","comps")
rownames(motif.community.summary) = 1:nrow(motif.community.summary)

# filter TFs based on those identified by multiple comparisons
onecomp = motif.community.summary %>% filter(ncomps==1)
motif.community.summary = motif.community.summary %>% filter(ncomps>1)
flt.pwm.clust.final = list()
for(ii in 1:nrow(motif.community.summary)){
  motif = motif.community.summary$id[ii]
  ind = which(names(flt.pwm.clust) == motif)
  flt.pwm.clust.final[[motif]] = flt.pwm.clust[[ind]]
}

# filter annotations of all cluster TFs
motif.community4 = c()
for(ii in 1:nrow(motif.community.summary)){
  center = motif.community.summary$id[ii]
  comm = motif.community3$community[motif.community3$motif == center]
  ind = which(motif.community3$community == comm)
  new = motif.community3[ind,]
  motif.community4 = rbind(motif.community4, new)
}
dim(motif.community4) # 523
length(unique(motif.community4$community)) # 211

# plot filtered/clustered TFs
final.plts = list()
p=1
for(ii in names(flt.pwm.clust.final)){
  id = motif.community4$ids.orig[motif.community4$motif==ii]
  mat = flt.pwm.clust[[ii]]
  final.plts[[p]] = logo.plt(pwm=mat, title=id)
  p=p+1
}
pdf("TF_final211.pdf", onefile = TRUE, height=9, width=4.5)
marrangeGrob(grobs=final.plts, nrow=6, ncol=1, top=NULL)
dev.off()

# generate speadsheet format for clusters
motif.clusters.id = list()
motif.clusters.tf = list()
for(ii in 1:length(unique(motif.community4$community))){
  clust = unique(motif.community4$community)[ii]
  ind = which(motif.community4$community == clust)
  namen = paste0("community_",clust)
  motif.clusters.id[[namen]] = motif.community4$motif[ind]
  motif.clusters.tf[[namen]] = motif.community4$ids.orig[ind]
}
maxlen = lapply(motif.clusters.id,function(x){length(x)}) %>% unlist %>% max
datTF = as.data.frame(matrix(0,maxlen,length(motif.clusters.tf)))
datID = as.data.frame(matrix(0,maxlen,length(motif.clusters.id)))
for(ii in 1:ncol(datID)){
  len = length(motif.clusters.id[[ii]])
  datID[1:len,ii] = motif.clusters.id[[ii]]
  datTF[1:len,ii] = motif.clusters.tf[[ii]]
}
names(datID) = names(datTF) = names(motif.clusters.tf)
datID = datID[,order(apply(datID,2,function(x)length(which(x!=0))),decreasing=T)]
datTF = datTF[,order(apply(datTF,2,function(x)length(which(x!=0))),decreasing=T)]
write.table(dat,"TFclusters.txt",col.names=T,row.names=F,sep="\t")

# generate pwm files for the 211 TF motifs of interest
pwm.out = flt.pwm.clust.final
cdir = getwd()
pwm.dir = paste0(cdir,"/TF211")
setwd(pwm.dir)
for(ii in 1:length(pwm.out)){
  tf = names(pwm.out)[ii]
  mat = t( pwm.out[[ii]] )
  fname = paste0(tf,".txt")
  write.table(mat,fname,sep="\t",quote=F,col.names=F,row.names=F)
}
setwd(cdir)

# save workspace
# save.image("TFanalysis_20191021.RData")

# next: TFfimo.sh



