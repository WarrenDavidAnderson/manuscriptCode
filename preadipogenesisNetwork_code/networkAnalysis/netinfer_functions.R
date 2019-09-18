
############################################################
## general functions
############################################################

# get bed with tss for each gene
get.tss.bed = function(gene.file=NULL,gene.list=NULL){
  inds = unlist( sapply(gene.list,function(x){which(gene.file[,4]==x)}) )
  mat = gene.file[inds,]
  ind.plus = which(mat[,6]=="+")
  ind.minus = which(mat[,6]=="-")
  mat[ind.plus,3] = mat[ind.plus,2] + 1
  mat[ind.minus,2] = mat[ind.minus,3] - 1
  return(mat)
}

# get PRO bed coords
get.pro.bed.coords = function(gene.list=NULL,bed=NULL){
  sig.tf.genes = sapply(gene.list,function(x){strsplit(x,"_")[[1]][2]})
  chr = sapply(gene.list,function(x){strsplit(x,":")[[1]][1]})
  start = sapply(gene.list,function(x){strsplit(x,":")[[1]][2]})
  start = sapply(start,function(x){strsplit(x,"-")[[1]][1]})
  end = sapply(gene.list,function(x){strsplit(x,"-")[[1]][2]})
  end = sapply(end,function(x){strsplit(x,"_")[[1]][1]})
  strand = sapply(sig.tf.genes,function(x){bed$strand[bed0$gene==x]})
  tf = cbind(gene.list, chr, start, end) %>% as.data.frame(stringsAsFactors=F)
  tf.bed = tf %>% select(chr,start,end)
  tf.bed = tf.bed %>% mutate(gene=sig.tf.genes, xy="xy", strand=strand)
  tf.bed[,2:3] = apply(tf.bed[,2:3],2,function(x){data.matrix(x)%>%as.numeric})
  return(tf.bed)
}

# get bed coords for early REs
get.atac.bed.coords = function(gene.list=NULL){
  chr = sapply(gene.list,function(x){strsplit(x,":")[[1]][1]})
  start = sapply(gene.list,function(x){strsplit(x,":")[[1]][2]})
  start = sapply(start,function(x){strsplit(x,"-")[[1]][1]})
  end = sapply(gene.list,function(x){strsplit(x,"-")[[1]][2]})
  re.bed = as.data.frame(cbind(chr,start,end),stringsAsFactors=FALSE)
  re.bed = re.bed %>% mutate(gene=gene.list, xy="xy",strand=NA)
  re.bed[,2:3] = apply(re.bed[,2:3],2,function(x){data.matrix(x)%>%as.numeric})
  return(re.bed)
}

# function to get specific contrast results
get.contrast.res = function(deg=NULL,contrast=NULL, conditions=NULL,
                            sig.thresh=NULL,fc.thresh=NULL,deseq=NULL){
  reads=counts(deseq) %>% as.data.frame
  ctrst = results(deg, contrast=contrast)
  ind.up = which(ctrst$padj < sig.thresh & ctrst$log2FoldChange > fc.thresh)
  ind.down = which(ctrst$padj < sig.thresh & ctrst$log2FoldChange < -fc.thresh)
  genes.up = rownames(ctrst)[ind.up]
  genes.dn = rownames(ctrst)[ind.down]
  nsig = length(which(ctrst$padj < sig.thresh & abs(ctrst$log2FoldChange) > fc.thresh))
  ind.ref=which(conditions==contrast[3])
  ind.cpr=which(conditions==contrast[2])
  ind.pks = sapply(genes.up,function(x)which(rownames(deseq)==x))
  #print(head(reads[ind.pks,c(ind.ref,ind.cpr)]))
  mean.ref = apply(reads[ind.pks,c(ind.ref)],1,mean)
  mean.cpr = apply(reads[ind.pks,c(ind.cpr)],1,mean)
  #print(head(cbind(mean.ref, mean.cpr),10))
  out = list(contrast=ctrst, genes.up=genes.up, genes.dn=genes.dn, nsig=nsig)
  return(out)
} # get.contrast.res

# function to get gene id
get.gene.id = function(vec=NULL){
  out = sapply(vec,function(x){
    g = strsplit(x,"_")[[1]]
    g = g[2:length(g)]
    if(length(g)>1){g = paste0(g, collapse="_")}
    return(g)
  })
  return(out)
} # get.gene.id

# function to isolate significant TFs for a particular comparison
# input: sig.genes = char vector of genes with coordinates (e.g., chr1:60434582-60566695_Raph1)
# input: tf.comm.ann = frame with "comm" and "tfs" annotations
get.sig.tfs = function(sig.genes=NULL, tf.comm.ann=NULL){
  sig.genes = sig.genes[-grep("tu_class",sig.genes)]
  sig.genes = sapply(sig.genes,function(x){strsplit(x,"_")[[1]][2]})
  sig.genes = cbind(names(sig.genes),sig.genes) %>% as.data.frame(stringsAsFactors=F)
  names(sig.genes) = c("coord","gene")
  sig.tf = tf.comm.ann$tf[ tf.comm.ann$tf %in% sig.genes$gene ]
  ind.tf = sapply(sig.tf,function(x){which(tf.comm.ann$tf==x)[1]}) %>% unlist
  ind.tu = sapply(sig.tf,function(x){which(sig.genes$gene==x)[1]}) %>% unlist
  sig.tfs = cbind(tf.comm.ann[ind.tf,], sig.genes[ind.tu,])
  return(sig.tfs)
} # get.sig.tfs

# function to apply bed coords to tu coordinates
coord.bed.tu = function(sig=NULL, bed0=NULL, ctrst=NULL){
  deg.genes = sapply(rownames(ctrst),function(x){strsplit(x,"_")[[1]]})
  deg.genes = lapply(deg.genes,function(x){paste(x[2:length(x)],collapse="_")})
  ind.bed = sapply(sig,function(x){which(bed0$gene==x)})
  ind.fc = sapply(sig,function(x){which(deg.genes==x)})
  bed.dat = data.frame(bed0[ind.bed,], log2fc=ctrst$log2FoldChange[ind.fc], stringsAsFactors=F)
  bed.dat[,c(2,3,7)] = apply(bed.dat[,c(2,3,7)],2,function(x){data.matrix(x) %>% as.numeric})
  rownames(bed.dat) = 1:nrow(bed.dat)
  return(bed.dat)
} # coord.bed.tu

# function to get bed coords for atac peaks
coord.bed.pk = function(pks=NULL){
  chr = sapply(pks,function(x){strsplit(x,":")[[1]][1]})
  start = sapply(pks,function(x){strsplit(x,":")[[1]][2]})
  start = sapply(start,function(x){strsplit(x,"-")[[1]][1]})
  end = sapply(pks,function(x){strsplit(x,"-")[[1]][2]})
  bed.dat = data.frame(chr=chr, start=start, end=end, gene=pks, xy=NA, strand=NA, stringsAsFactors=F)
  bed.dat[,c(2,3)] = apply(bed.dat[,c(2,3)],2,function(x){data.matrix(x) %>% as.numeric})
  rownames(bed.dat) = 1:nrow(bed.dat)
  return(bed.dat)
} # coord.bed.pk

############################################################
## TF -> RE edges 
############################################################

# function for computing w for TF --> RE edges
# this function is called inside pro.edge.weight()
# input: dist = a set of distances between atac summits and TF motifs
# inputs = dist1 and dist2
# w will be 1 for dist < dist1
# w will be m * dist + b for dist1 < dist < dist2
tfre.compute.w = function(dist=NULL, dist1=50, dist2=300){
  m = -1 / (dist2 - dist1)
  b = -dist2 * m
  w = sapply(dist,function(x){
    if(x <= dist1){w=1}
    if(x>dist1 & x<=dist2){w=m*x+b}
    if(x>dist2){w=0}
    return(w)
  })
  return(w)
} # tfre.compute.w

# get tf --> re edge weights
# input: fc.dat = fold change data, deseq lrt contrast data frame
# input: summits = atac peaks with summits (e.g., "chr1:1-10,5")
# input: peak.ann = annotation for peak motif counts
# input: gen.ann = general annotation for TF communities and motifs
# input: bed.dir = bedtools directory, e.g., "bed/directory/"
# input: fimo.bed.path = path for fimo bed files, e.g., "fimo/path/"
# input: dist1 = distance up to which w=1
# input: dist2 = linear interpolation of w betweed dist1 and dist2
pro.edge.weight = function(fc.dat=NULL, summits=NULL, peak.ann=NULL, gen.ann=NULL,
                           bed.dir="/media/wa3j/Seagate2/Documents/software/bedtools2/bin/",
                           fimo.bed.path=NULL, dist1=50, dist2=300,
                           bigWig.path.atac=NULL, sizefac=NULL,
                           cond.atac1=NULL, cond.atac2=NULL){
  
  # atac read files and size factors
  atac.files = list.files(bigWig.path.atac)
  atac.files1 = atac.files[sapply(cond.atac1,function(x){grep(x,atac.files)}) %>% as.vector]
  atac.files2 = atac.files[sapply(cond.atac2,function(x){grep(x,atac.files)}) %>% as.vector]
  atac.sf1 = sizefac[sapply(cond.atac1,function(x){grep(x,names(sizefac))}) %>% as.vector]
  atac.sf2 = sizefac[sapply(cond.atac2,function(x){grep(x,names(sizefac))}) %>% as.vector]
  
  # loop through each tu, get all peaks within interval, compute w * s
  edge.weights = c()
  for(ii in 1:nrow(fc.dat)){
    
    # basic tf info
    namen = rownames(fc.dat)[ii]
    tf.chr = strsplit(namen,":")[[1]][1]
    tf.start = strsplit(namen,"-")[[1]][1]
    tf.start = strsplit(tf.start,":")[[1]][2]
    tf.end = strsplit(namen,"_")[[1]][1]
    tf.end = strsplit(tf.end,"-")[[1]][2]
    tf.id = strsplit(rownames(fc.dat)[ii],"_")[[1]][2]
    tf.fc = fc.dat$log2FoldChange[ii]
    motif.id = gen.ann$motifid[gen.ann$tf == tf.id]
    
    # isolate all atac peaks for the tf
    ind = which(names(peak.ann) == gen.ann$TFid[gen.ann$tf==tf.id])
    pks = peak.ann[,ind]
    names(pks) = rownames(peak.ann)
    pks = pks[pks>0]
    pks = names(pks)
    
    # isolate the peak summits
    pk.summits = sapply(pks,function(c){summits[grep(c,summits)]})
    chr.summits = sapply(pk.summits,function(x){strsplit(x,":")[[1]][1]})
    pk.summits = sapply(pk.summits,function(x){strsplit(x,",")[[1]][2]}) 
    pk.summits = as.numeric(pk.summits)
    
    # isolate the closest motif to each summit
    fname = paste0(fimo.bed.path,"fimo_",motif.id,".bed")
    bed.fimo = read.table(fname,header=F,stringsAsFactors=F)
    bed.com = paste0(bed.dir,"closestBed ")
    
    # get the closest motif distance for every peak
    bed.summit = data.frame(chr=chr.summits, start=pk.summits, end=pk.summits+1, pks=pks)
    write.table(bed.summit[,1:4],"a1.bed",col.names=F,row.names=F,sep="\t",quote=F)
    write.table(bed.fimo[,1:3],"b1.bed",col.names=F,row.names=F,sep="\t",quote=F)
    system("sort -k1,1 -k2,2n a1.bed > a.bed")
    system("sort -k1,1 -k2,2n b1.bed > b.bed")
    system(paste0(bed.com,"-d -a a.bed -b b.bed > dat.bed"))
    info = file.info("dat.bed")
    if(info$size == 0){next}
    dat = read.table("dat.bed",header=F,stringsAsFactors=F)
    system("rm a.bed b.bed a1.bed b1.bed dat.bed") # SU-summit, MT-motif
    names(dat) = c("chrSU","startSU","endSU","coordSU","chrMT","startMT","endMT","dist")
    
    # remove any duplicates and compute closest distances
    ind.dup = which(duplicated(dat$coordSU)==T)
    if(length(ind.dup)>0){
      dups = dat$coordSU[ind.dup]
      for(jj in 1:length(unique(dups))){
        ind.all = which(dat$coordSU == unique(dups)[jj])
        ind.min = which(dat$dist[ind.all] == min(dat$dist[ind.all]))[1]
        dat = dat[-ind.all[-ind.min],]
      }
    }
    ind = sapply(pk.summits,function(x){which(dat$startSU==x)}) %>% unlist
    dat = dat[ind,]
    d = sapply(c(1:nrow(dat)),function(jj){
      d = min( abs(pk.summits[jj] - dat$startMT), abs(pk.summits[jj] - dat$endMT) )
      return(d)
    })
    
    # calculate w
    w = tfre.compute.w(dist=d, dist1=dist1, dist2=dist2)
    
    # calculate s2 - atac read density difference
    # for every peak, calculate the read difference over the interval
    atac.chr = sapply(dat$coordSU,function(x){strsplit(x,":")[[1]][1]})
    atac.str = sapply(dat$coordSU,function(x){strsplit(x,":")[[1]][2]})
    atac.str = sapply(atac.str,function(x){strsplit(x,"-")[[1]][1]})
    atac.end = sapply(dat$coordSU,function(x){strsplit(x,"-")[[1]][2]})
    bed.atac = data.frame(chr=atac.chr, start=atac.str, end=atac.end)
    bed.atac[,2:3] = apply(bed.atac[,2:3],2,function(x){data.matrix(x) %>% as.numeric})
    s2 = rep(0,nrow(dat))
    for(jj in 1:nrow(dat)){
      atac.reads1 = sapply(atac.files1,function(x){
        bw = load.bigWig(paste0(bigWig.path.atac,"/",x))
        cnt = bed.region.bpQuery.bigWig(bw=bw, bed=bed.atac[jj,])
        unload.bigWig(bw)
        return(cnt)
      })
      atac.reads2 = sapply(atac.files2,function(x){
        bw = load.bigWig(paste0(bigWig.path.atac,"/",x))
        cnt = bed.region.bpQuery.bigWig(bw=bw, bed=bed.atac[jj,])
        unload.bigWig(bw)
        return(cnt)
      })
      atac.reads1 = atac.reads1 / atac.sf1
      atac.reads2 = atac.reads2 / atac.sf2
      s2[jj] = (mean(atac.reads2)-mean(atac.reads1)) / (bed.atac[jj,3] - bed.atac[jj,2])
    } # jj, calculate s2
    
    # generate output
    re.chr = sapply(pks,function(x){strsplit(x,":")[[1]][1]})
    re.start = sapply(pks,function(x){strsplit(x,"-")[[1]][1]})
    re.start = sapply(re.start,function(x){strsplit(x,":")[[1]][2]})
    re.end = sapply(pks,function(x){strsplit(x,"-")[[1]][2]})
    new = data.frame(chrTU=tf.chr, startTU=tf.start, endTU=tf.end, geneTU=tf.id,
                     chrRE=re.chr, startRE=re.start, endRE=re.end, idRE=pks,
                     dist=d, w=w, s1=tf.fc, s2=s2, weight = w * tf.fc * s2)
    edge.weights = rbind(edge.weights, new)
    print(ii / nrow(fc.dat))
    
  } # ii, tu loop
  
  # output results
  rownames(edge.weights) = 1:nrow(edge.weights)
  edge.weights[,c(2,3,6,7,9:13)] = apply(edge.weights[,c(2,3,6,7,9:13)],2,function(x){
    data.matrix(x) %>% as.numeric
  })
  out = edge.weights[edge.weights$weight != 0,]
  out = out[order(abs(out$weight),decreasing=TRUE),]
  return(out)

} # pro.edge.weight

############################################################
## RE -> TU edges 
############################################################

# get re --> tu edge weights
# input: bed.dat = bed data for dynamic TUs
# input: peak.bed = atac peak coordinates in bed6 format
# input: summits = atac peaks with summits (e.g., "chr1:1-10,5")
# input: atac.lim = search limit for re -> tu eage weights
# input: atac.half = exponential decay at half max for re -> tu eage weights (kb)
# input: tau = rate constant for exponential decay model of the regulatory weight (bp) 
# input: cond.atac = condition(s) for computing average peak density (e.g., c("20m","40m"))
# input: bed.dat = TU bed data (e.g., dynTF0)
# input: class = comparison class (e.g., "early_up")
# input: bed.dir = bedtools directory, e.g., "bed/directory/"
atac.edge.weight = function(bed.dat=NULL, peak.bed=NULL, summits=NULL,
                            atac.lim=100, tau=NULL,
                            cond.atac1=NULL, cond.atac2=NULL, sizefac=NULL,
                            bigWig.path.atac=NULL,
                            bed.dir="/media/wa3j/Seagate2/Documents/software/bedtools2/bin/"){
  
  # get bed coords and fold change data for dynamic TUs
  bed.dat = bed.dat[duplicated(bed.dat)==FALSE,]
  bed.data = data.frame(chr=bed.dat$chr, start=bed.dat$start, 
                       end=bed.dat$end, gene=bed.dat$gene, 
                       xy="na", strand=bed.dat$strand)
  tu.fc = bed.dat$log2fc
  
  # get bed with tss for each gene
  tss.bed = get.tss.bed(gene.file=bed.data, gene.list=bed.data$gene)
  
  # atac read files
  atac.files = list.files(bigWig.path.atac)
  atac.files1 = atac.files[sapply(cond.atac1,function(x){grep(x,atac.files)}) %>% as.vector]
  atac.files2 = atac.files[sapply(cond.atac2,function(x){grep(x,atac.files)}) %>% as.vector]
  
  # atac size factors
  atac.sf1 = sizefac[sapply(cond.atac1,function(x){grep(x,names(sizefac))}) %>% as.vector]
  atac.sf2 = sizefac[sapply(cond.atac2,function(x){grep(x,names(sizefac))}) %>% as.vector]
  
  # loop through each tss, get all peaks within interval, compute w * s
  edge.weights = c()
  for(ii in 1:nrow(tss.bed)){
    
    # expand tss bed coords to include the search region
    tss.bed.ii = tss.bed[ii,]
    tss.bed.ii$start = tss.bed.ii$start - atac.lim*1000
    tss.bed.ii$end = tss.bed.ii$end + atac.lim*1000
    
    # find all sig peaks near the tss
    bed.com = paste0(bed.dir,"intersectBed ")
    write.table(tss.bed.ii[,1:4],"a.bed",col.names=F,row.names=F,sep="\t",quote=F)
    write.table(peak.bed[,1:4],"b1.bed",col.names=F,row.names=F,sep="\t",quote=F)
    system("sort -k1,1 -k2,2n b1.bed > b.bed")
    system(paste0(bed.com,"-wo -a a.bed -b b.bed > dat.bed"))
    info = file.info("dat.bed")
    if(info$size == 0){next}
    dat = read.table("dat.bed",header=F,stringsAsFactors=F)
    system("rm a.bed b1.bed b.bed dat.bed")
    names(dat) = c("chrTU","startTU","endTU","geneTU","chrRE","startRE","endRE","idRE","weight")
    
    # get weight parameter for each overlapping atac peak
    w = rep(0,nrow(dat))
    d = rep(0,nrow(dat))
    for(jj in 1:nrow(dat)){
      summit = summits[grep(dat$idRE[jj],summits)]
      summit = strsplit(summit,",")[[1]][2] %>% as.numeric
      if(tss.bed.ii$strand=="+"){tss = tss.bed$start[ii]}
      if(tss.bed.ii$strand=="-"){tss = tss.bed$end[ii]}
      dd = abs(tss - summit)
      d[jj] = dd
      w[jj] = 2*exp(-dd/tau) / (1 + exp(-dd/tau))
    } # jj, calculate w
    
    # get the difference in read density for each overlapping atac peak
    s2 = rep(0,nrow(dat))
    for(jj in 1:nrow(dat)){
      atac.reads1 = sapply(atac.files1,function(x){
        bw = load.bigWig(paste0(bigWig.path.atac,"/",x))
        cnt = bed.region.bpQuery.bigWig(bw=bw, bed=dat[jj,5:7])
        unload.bigWig(bw)
        return(cnt)
      })
      atac.reads2 = sapply(atac.files2,function(x){
        bw = load.bigWig(paste0(bigWig.path.atac,"/",x))
        cnt = bed.region.bpQuery.bigWig(bw=bw, bed=dat[jj,5:7])
        unload.bigWig(bw)
        return(cnt)
      })
      atac.reads1 = atac.reads1 / atac.sf1
      atac.reads2 = atac.reads2 / atac.sf2
      s2[jj] = (mean(atac.reads2)-mean(atac.reads1)) / (dat$endRE[jj] - dat$startRE[jj])
    } # jj, calculate s
    
    # incorporate edge weights
    new = data.frame(dat[,1:8], dist=d, w=w, s1=tu.fc[ii], s2=s2, weight=w*tu.fc[ii]*s2, stringsAsFactors=F)
    edge.weights = rbind(edge.weights, new)
    
  } # ii, tss loop
  
  # output results
  if(is.null(edge.weights)==TRUE){
    return("nodata")
  } else {
    edge.weights[,c(2,3,6,7,9:13)] = apply(edge.weights[,c(2,3,6,7,9:13)],2,function(x){
      data.matrix(x) %>% as.numeric
    })
    out = edge.weights[edge.weights$weight != 0,]
    out = out[order(abs(out$weight),decreasing=TRUE),]
    return(out) 
  }
  
} # atac.edge.weight


############################################################
## functions for building the entire network 
## and verifying connectedness 
############################################################

# establish the network for a vector of time points
infer.net.vec = function(deg.pro=NULL, deseq.obj.pro=NULL,
                         deg.atac=NULL, deseq.obj.atac=NULL,
                         sig.thresh=0.001, fc.thresh=1, 
                         bigWig.path.atac=NULL, fimo.bed.path=NULL,
                         tfclass.ann=NULL, bed0=NULL, 
                         tnet=NULL, tcnd=NULL, all.summits=NULL,
                         atac.lim=100, atac.half=10,
                         motif.mapping=NULL, dist1=50, dist2=300,
                         bed.dir=NULL, full=FALSE){
  
  # space constand for the atach edge weighting term w
  tau = -1 / (log(1/3) / (atac.half*1000))
  
  # loop through each time point, starting at second time point
  # even: re --> tu, odd: tf --> re
  net.full = c()
  Eretf = c()
  Etfre = c()
  for(ii in 2:length(tnet)){
    
    ################################
    # get tf --> re edges 
    if(ii >= 3){
      
      # pro and atac contrasts, distinct time point comparisons
      contrast = c("conditions",tcnd[ii-1],tcnd[ii-2])
      conditions = colData(deseq.obj.pro)$conditions
      ctrst.pro = get.contrast.res(deg=deg.pro,contrast=contrast,
                                   sig.thresh=sig.thresh,fc.thresh=fc.thresh,
                                   conditions=conditions, deseq=deseq.obj.pro)
      contrast = c("conditions",tcnd[ii],tcnd[ii-1])
      conditions = colData(deseq.obj.atac)$conditions
      ctrst.atac = get.contrast.res(deg=deg.atac,contrast=contrast,
                                    sig.thresh=sig.thresh,fc.thresh=fc.thresh,
                                    conditions=conditions, deseq=deseq.obj.atac)
      
      # isolate dynamic peaks 
      res.atac.up0 = ctrst.atac$genes.up
      res.atac.dn0 = ctrst.atac$genes.dn
      
      # subset peak/TF matrix to include dynamic peaks
      peak.matrix.up = motif.mapping[rownames(motif.mapping) %in% res.atac.up0,]
      peak.matrix.dn = motif.mapping[rownames(motif.mapping) %in% res.atac.dn0,]
      
      # get fold change data for dynamic tfs, filter for regulated tfs
      pro.up = ctrst.pro$genes.up
      pro.dn = ctrst.pro$genes.dn
      res.pro.up = ctrst.pro$contrast[rownames(ctrst.pro$contrast) %in% pro.up,]
      res.pro.dn = ctrst.pro$contrast[rownames(ctrst.pro$contrast) %in% pro.dn,]
      res.pro.up = as.data.frame(res.pro.up, stringsAsFactors=FALSE)
      res.pro.dn = as.data.frame(res.pro.dn, stringsAsFactors=FALSE)
      up.genes = sapply(rownames(res.pro.up),function(x){strsplit(x,"_")[[1]][2]})
      dn.genes = sapply(rownames(res.pro.dn),function(x){strsplit(x,"_")[[1]][2]})
      
      # filter new tfs to include only previous re targets
      if(full==FALSE){
        res.pro.up = res.pro.up[up.genes %in% Eretf$geneTU,]
        res.pro.dn = res.pro.dn[dn.genes %in% Eretf$geneTU,]
      }
      
      # get tf --> re edges
      cond.atac1 = tnet[ii-1]
      cond.atac2 = tnet[ii]
      sizefac = sizeFactors(deseq.obj.atac)
      if(nrow(res.pro.up) != 0){
        Eupup = pro.edge.weight(fc.dat=res.pro.up, peak.ann=peak.matrix.up, summits=all.summits, 
                                gen.ann=tfclass.ann, dist1=dist1, dist2=dist2,
                                bed.dir=bed.dir, fimo.bed.path=fimo.bed.path, 
                                bigWig.path.atac=bigWig.path.atac, sizefac=sizefac, 
                                cond.atac1=cond.atac1, cond.atac2=cond.atac2)
        Eupdn = pro.edge.weight(fc.dat=res.pro.up, peak.ann=peak.matrix.dn, summits=all.summits, 
                                gen.ann=tfclass.ann, dist1=dist1, dist2=dist2,
                                bed.dir=bed.dir, fimo.bed.path=fimo.bed.path, 
                                bigWig.path.atac=bigWig.path.atac, sizefac=sizefac, 
                                cond.atac1=cond.atac1, cond.atac2=cond.atac2)
      }
      if(nrow(res.pro.dn) != 0){
        Ednup = pro.edge.weight(fc.dat=res.pro.dn, peak.ann=peak.matrix.up, summits=all.summits, 
                                gen.ann=tfclass.ann, dist1=dist1, dist2=dist2,
                                bed.dir=bed.dir, fimo.bed.path=fimo.bed.path, 
                                bigWig.path.atac=bigWig.path.atac, sizefac=sizefac, 
                                cond.atac1=cond.atac1, cond.atac2=cond.atac2)
        Edndn = pro.edge.weight(fc.dat=res.pro.dn, peak.ann=peak.matrix.dn, summits=all.summits, 
                                gen.ann=tfclass.ann, dist1=dist1, dist2=dist2,
                                bed.dir=bed.dir, fimo.bed.path=fimo.bed.path, 
                                bigWig.path.atac=bigWig.path.atac, sizefac=sizefac, 
                                cond.atac1=cond.atac1, cond.atac2=cond.atac2)
      }
      
      # aggregate results and annotate the net
      Etfre0 = rbind(Eupup, Eupdn, Ednup, Edndn) %>% mutate(type="TFtoRE",t1=tcnd[ii-1], t2=tcnd[ii])
      Etfre = Etfre0 %>% mutate(tu="tf") %>% mutate(tt=paste0(tcnd[ii-1],",",tcnd[ii]))
      net.full = rbind(net.full, Etfre)
      Eupup=Eupdn=Ednup=Edndn=c()
      Eretf = c()
      
    } # get tf --> re edges 
    ################################
    
    
    ################################
    # get re --> tu/tf edges 
    if(ii >= 2){
      
      # pro and atac contrasts, same time point comparison
      contrast = c("conditions",tcnd[ii],tcnd[ii-1])
      conditions = colData(deseq.obj.pro)$conditions
      ctrst.pro = get.contrast.res(deg=deg.pro,contrast=contrast,
                                   sig.thresh=sig.thresh,fc.thresh=fc.thresh,
                                   conditions=conditions, deseq=deseq.obj.pro)
      conditions = colData(deseq.obj.atac)$conditions
      ctrst.atac = get.contrast.res(deg=deg.atac,contrast=contrast,
                                    sig.thresh=sig.thresh,fc.thresh=fc.thresh,
                                    conditions=conditions, deseq=deseq.obj.atac)
      
      # separate significant tus and tfs
      genes.up = get.gene.id(ctrst.pro[["genes.up"]])
      genes.dn = get.gene.id(ctrst.pro[["genes.dn"]])
      tfs.up = get.sig.tfs(sig.genes=ctrst.pro[["genes.up"]], tf.comm.ann=tfclass.ann)
      tfs.dn = get.sig.tfs(sig.genes=ctrst.pro[["genes.dn"]], tf.comm.ann=tfclass.ann)
      tus.up = genes.up[ !(genes.up %in% tfs.up$tf) ]
      tus.dn = genes.dn[ !(genes.dn %in% tfs.dn$tf) ]
      
      # get tf/tu bed coords and fold changes
      if(ii != length(tnet)){
        tu.bed.up = coord.bed.tu(sig=tfs.up$tf, bed0=bed0, ctrst=ctrst.pro[["contrast"]])
        tu.bed.dn = coord.bed.tu(sig=tfs.dn$tf, bed0=bed0, ctrst=ctrst.pro[["contrast"]])
      } 
      if (ii == length(tnet)){ # use tus only for the last time comparison
        tu.bed.up = coord.bed.tu(sig=genes.up, bed0=bed0, ctrst=ctrst.pro[["contrast"]])
        tu.bed.dn = coord.bed.tu(sig=genes.dn, bed0=bed0, ctrst=ctrst.pro[["contrast"]])
      }
      
      # get atac peak bed coords
      # filter for peaks associated with dynamic tfs
      pk.bed.up = coord.bed.pk(pks=ctrst.atac$genes.up)
      pk.bed.dn = coord.bed.pk(pks=ctrst.atac$genes.dn)
      if( full==FALSE & ii!=2 ){ # use atac peaks assoc with tfs
        pk.bed.up = pk.bed.up[pk.bed.up$gene %in% Etfre$idRE,]
        pk.bed.dn = pk.bed.dn[pk.bed.dn$gene %in% Etfre$idRE,]
      }
      
      # get re --> tu edges
      cond.atac1 = tnet[ii-1]
      cond.atac2 = tnet[ii]
      sizefac = sizeFactors(deseq.obj.atac)
      Eup = atac.edge.weight(bed.dat=tu.bed.up, peak.bed=pk.bed.up, summits=all.summits,
                             bigWig.path.atac=bigWig.path.atac, bed.dir=bed.dir,
                             atac.lim=atac.lim, tau=tau, sizefac=sizefac,
                             cond.atac1=cond.atac1, cond.atac2=cond.atac2)
      Edn = atac.edge.weight(bed.dat=tu.bed.dn, peak.bed=pk.bed.dn, summits=all.summits,
                             bigWig.path.atac=bigWig.path.atac, bed.dir=bed.dir,
                             atac.lim=atac.lim, tau=tau, sizefac=sizefac,
                             cond.atac1=cond.atac1, cond.atac2=cond.atac2)
      
      # aggregate results and annotate the net
      Eretf0 = rbind(Eup, Edn) %>% mutate(type="REtoTF",t1=tcnd[ii-1], t2=tcnd[ii])
      Etf = Eretf0[Eretf0$geneTU %in% tfclass.ann$tf,] %>% mutate(tu="tf")
      Etu = Eretf0[!(Eretf0$geneTU %in% tfclass.ann$tf),] %>% mutate(tu="tu")
      Eretf = Etf %>% mutate(tt=paste0(tcnd[ii-1],",",tcnd[ii]))
      Eall = rbind(Etf, Etu) %>% mutate(tt=paste0(tcnd[ii-1],",",tcnd[ii]))
      net.full = rbind(net.full, Eall)
      Eup=Edn=c()
      Etfre = c()
      
    } # get re --> tu/tf edges 
    ################################
    
  } # ii, net inference loop
  
  # filter output
  ind.rm = grep("nodata",net.full[,1])
  if(length(ind.rm)>0){
    net = net.full[-ind.rm,]
  } else {
    net = net.full
  }
  
  # generate a strictly connected network and output results
  net[,c(2,3,6,7,9:13)] = apply(net[,c(2,3,6,7,9:13)],2,function(x){
    data.matrix(x) %>% as.numeric
  })
  res = force.connect(net)
  return(res)
  
} # infer.net.vec


# function to enforce conectivity from the start to end time
force.connect = function(network=NULL, lim=10000){
  
  # loop foward, but repeat until there are no changes
  df1 = 1
  df2 = 2
  net = network
  all.time = unique(net$tt)
  types = c("REtoTF", "TFtoRE")
  iter = 0
  contin = "yes"
  connections = c()
  while(identical(df1,df2) != TRUE){
    
    # loop through time frames, starting at the second
    new = c()
    for(ii in 2:length(all.time)){
      
      # get indices for subsetting data
      comp.ii = all.time[ii]
      comp.de = all.time[ii-1]
      ind.ii = which(net$tt == comp.ii)
      ind.de = which(net$tt == comp.de)
      ind.re = which(net$type == "REtoTF")
      ind.tf = which(net$type == "TFtoRE")
      
      # limit the tfs from the initial re-->tf edges
      dat.tfre.ii = net[intersect(ind.ii,ind.tf),]
      if(ii == 2){
        dat.retf.de = net[intersect(ind.de,ind.re),]
        dat.retf1 = dat.retf.de[dat.retf.de$geneTU %in% dat.tfre.ii$geneTU, ]
      }
      if(ii > 2){
        ind.de = which(new$tt == comp.de & new$type == "REtoTF")
        dat.retf.de = new[ind.de,]
        new = new[-ind.de,]
        dat.retf1 = dat.retf.de[dat.retf.de$geneTU %in% dat.tfre.ii$geneTU, ]
      }
      
      # limit current tfs based on previous tfs
      dat.tfre = dat.tfre.ii[dat.tfre.ii$geneTU %in% dat.retf1$geneTU,]
      #all(dat.retf.de$geneTU %in% dat.tfre.ii$geneTU)
      #all(dat.tfre.ii$geneTU %in% dat.retf.de$geneTU)
      
      # limit current res
      dat.retf.ii = net[intersect(ind.ii,ind.re),]
      dat.retf2 = dat.retf.ii[dat.retf.ii$idRE %in% dat.tfre$idRE,]
      dat.tfre = dat.tfre[dat.tfre$idRE %in% dat.retf2$idRE,]
      
      # update data output
      new = rbind(new, dat.retf1, dat.tfre, dat.retf2)
      
    } # time loop
    
    # update the network after each loop iteration
    net = new
    iter = iter + 1
    connections = c(connections, nrow(new))
    if(iter == lim){contin = "no"; break}
    
    # keep the current and previous frames
    if(iter==1){
      df1 = new
    } 
    if(iter==2){
      df2 = new
    }
    if(iter>=3){
      df1 = df2
      df2 = new
    }
    
  } # while loop
  
  # check for connectivity
  conn = "yes"
  cnts = data.frame( matrix(0,length(all.time),3) )
  names(cnts) = c("time","nTFtoRE","nREtoTF")
  cnts$time = all.time
  for(ii in 1:length(all.time)){
    ind.re = which(new$type == "REtoTF" & new$tt==all.time[ii])
    ind.tf = which(new$type == "TFtoRE" & new$tt==all.time[ii])
    cnts$nREtoTF[ii] = length(ind.re)
    cnts$nTFtoRE[ii] = length(ind.tf)
    if(length(ind.re)==0 | length(ind.tf)==0){
      conn = "no"
    }
    if(length(ind.re)!=0 & length(ind.tf)==0 & ii==1){
      conn = "yes"
    }
  }
  
  # output results
  if(connections[length(connections)]==0){contin = "no"}
  out = list(net=new, cnts=cnts, conn=conn, contin=contin, ncon=connections)
  return(out)
  
} # force.connect



# establish the network for a list of time points
infer.net.list = function(deg.pro=NULL, deseq.obj.pro=NULL,
                         deg.atac=NULL, deseq.obj.atac=NULL,
                         sig.thresh=0.001, fc.thresh=1, 
                         bigWig.path.atac=NULL, fimo.bed.path=NULL,
                         tfclass.ann=NULL, bed0=NULL, 
                         tnet=NULL, tcnd=NULL, 
                         all.summits=NULL,
                         atac.lim=100, atac.half=10,
                         motif.mapping=NULL, dist1=50, dist2=300,
                         bed.dir=NULL, full=FALSE){
  
  # space constand for the atach edge weighting term w
  tau = -1 / (log(1/3) / (atac.half*1000))
  
  # loop through each time point, starting at second time point
  # even: re --> tu, odd: tf --> re
  net.full = c()
  Eretf = c()
  Etfre = c()
  for(ii in 2:length(tnet)){
    
    ################################
    # get tf --> re edges 
    if(ii >= 3){
      
      # pro contrasts, intersect of timepoints
      tci = tcnd[[ii-2]]
      tcf = tcnd[[ii-1]]
      cnt = 1
      ctrst.pro = list()
      for(xx in 1:length(tci)){
        for(yy in 1:length(tcf)){
          contrast = c("conditions",tcf[yy],tci[xx])
          conditions = colData(deseq.obj.pro)$conditions
          ctrst.pro[[cnt]] = get.contrast.res(deg=deg.pro,contrast=contrast,
                                       sig.thresh=sig.thresh,fc.thresh=fc.thresh,
                                       conditions=conditions, deseq=deseq.obj.pro)
          cnt = cnt + 1
        }
      }
      ind.up = list()
      ind.dn = list()
      for(jj in 1:length(ctrst.pro)){
        dat = ctrst.pro[[jj]]$contrast
        ind.up[[jj]] = which(dat$log2FoldChange > fc.thresh & dat$padj < sig.thresh)
        ind.dn[[jj]] = which(dat$log2FoldChange < (-1)*fc.thresh & dat$padj < sig.thresh)
      }
      indup.pro = Reduce(intersect, ind.up)
      inddn.pro = Reduce(intersect, ind.dn)
      
      # atac contrasts, intersect of timepoints
      tci = tcnd[[ii-1]]
      tcf = tcnd[[ii]]
      cnt = 1
      ctrst.atac = list()
      for(xx in 1:length(tci)){
        for(yy in 1:length(tcf)){
          contrast = c("conditions",tcf[yy],tci[xx])
          conditions = colData(deseq.obj.atac)$conditions
          ctrst.atac[[cnt]] = get.contrast.res(deg=deg.atac,contrast=contrast,
                                        sig.thresh=sig.thresh,fc.thresh=fc.thresh,
                                        conditions=conditions, deseq=deseq.obj.atac)
          cnt = cnt + 1
        }
      }
      ind.up = list()
      ind.dn = list()
      for(jj in 1:length(ctrst.atac)){
        dat = ctrst.atac[[jj]]$contrast
        ind.up[[jj]] = which(dat$log2FoldChange > fc.thresh & dat$padj < sig.thresh)
        ind.dn[[jj]] = which(dat$log2FoldChange < (-1)*fc.thresh & dat$padj < sig.thresh)
      }
      indup.atac = Reduce(intersect, ind.up)
      inddn.atac = Reduce(intersect, ind.dn)
      
      # isolate dynamic peaks 
      res.atac.up0 = rownames(ctrst.atac[[1]]$contrast)[indup.atac]
      res.atac.dn0 = rownames(ctrst.atac[[1]]$contrast)[inddn.atac]
      
      # subset peak/TF matrix to include dynamic peaks
      peak.matrix.up = motif.mapping[rownames(motif.mapping) %in% res.atac.up0,]
      peak.matrix.dn = motif.mapping[rownames(motif.mapping) %in% res.atac.dn0,]
      
      # get fold change data for dynamic tfs, filter for regulated tfs
      # compute mean fold change
      fc = list()
      for(jj in 1:length(ctrst.pro)){
        dat = ctrst.pro[[jj]]$contrast
        fc[[jj]] = dat$log2FoldChange
      }
      meanfc = colMeans( do.call(rbind, fc) )
      ctrst = ctrst.pro[[1]]$contrast
      ctrst$log2FoldChange = meanfc
      pro.up = rownames(ctrst)[indup.pro]
      pro.dn = rownames(ctrst)[inddn.pro]
      res.pro.up = ctrst[rownames(ctrst) %in% pro.up,]
      res.pro.dn = ctrst[rownames(ctrst) %in% pro.dn,]
      res.pro.up = as.data.frame(res.pro.up, stringsAsFactors=FALSE)
      res.pro.dn = as.data.frame(res.pro.dn, stringsAsFactors=FALSE)
      up.genes = sapply(rownames(res.pro.up),function(x){strsplit(x,"_")[[1]][2]})
      dn.genes = sapply(rownames(res.pro.dn),function(x){strsplit(x,"_")[[1]][2]})
      
      # filter new tfs to include only previous re targets
      if(full==FALSE){
        res.pro.up = res.pro.up[up.genes %in% Eretf$geneTU,]
        res.pro.dn = res.pro.dn[dn.genes %in% Eretf$geneTU,]
      }
      
      # get tf --> re edges
      cond.atac1 = tnet[[ii-1]]
      cond.atac2 = tnet[[ii]]
      sizefac = sizeFactors(deseq.obj.atac)
      Eupup=Eupdn=Ednup=Edndn=c()
      if(nrow(res.pro.up) != 0){
        Eupup = pro.edge.weight(fc.dat=res.pro.up, peak.ann=peak.matrix.up, summits=all.summits, 
                                gen.ann=tfclass.ann, dist1=dist1, dist2=dist2,
                                bed.dir=bed.dir, fimo.bed.path=fimo.bed.path, 
                                bigWig.path.atac=bigWig.path.atac, sizefac=sizefac, 
                                cond.atac1=cond.atac1, cond.atac2=cond.atac2)
        Eupdn = pro.edge.weight(fc.dat=res.pro.up, peak.ann=peak.matrix.dn, summits=all.summits, 
                                gen.ann=tfclass.ann, dist1=dist1, dist2=dist2,
                                bed.dir=bed.dir, fimo.bed.path=fimo.bed.path, 
                                bigWig.path.atac=bigWig.path.atac, sizefac=sizefac, 
                                cond.atac1=cond.atac1, cond.atac2=cond.atac2)
      }
      if(nrow(res.pro.dn) != 0){
        Ednup = pro.edge.weight(fc.dat=res.pro.dn, peak.ann=peak.matrix.up, summits=all.summits, 
                                gen.ann=tfclass.ann, dist1=dist1, dist2=dist2,
                                bed.dir=bed.dir, fimo.bed.path=fimo.bed.path, 
                                bigWig.path.atac=bigWig.path.atac, sizefac=sizefac, 
                                cond.atac1=cond.atac1, cond.atac2=cond.atac2)
        Edndn = pro.edge.weight(fc.dat=res.pro.dn, peak.ann=peak.matrix.dn, summits=all.summits, 
                                gen.ann=tfclass.ann, dist1=dist1, dist2=dist2,
                                bed.dir=bed.dir, fimo.bed.path=fimo.bed.path, 
                                bigWig.path.atac=bigWig.path.atac, sizefac=sizefac, 
                                cond.atac1=cond.atac1, cond.atac2=cond.atac2)
      }
      
      # aggregate results and annotate the net
      Etfre0 = rbind(Eupup, Eupdn, Ednup, Edndn) %>% mutate(type="TFtoRE",
                      t1=paste0(tcnd[[ii-1]],collapse=","), t2=paste0(tcnd[[ii]],collapse=","))
      Etfre = Etfre0 %>% mutate(tu="tf") %>% mutate(tt=paste0(t1,"_",t2))
      net.full = rbind(net.full, Etfre)
      Eupup=Eupdn=Ednup=Edndn=c()
      Eretf = c()
      
    } # get tf --> re edges 
    ################################
    
    
    ################################
    # get re --> tu/tf edges 
    if(ii >= 2){
      
      # pro contrasts, intersect of timepoints
      tci = tcnd[[ii-1]]
      tcf = tcnd[[ii]]
      cnt = 1
      ctrst.pro = list()
      for(xx in 1:length(tci)){
        for(yy in 1:length(tcf)){
          contrast = c("conditions",tcf[yy],tci[xx])
          conditions = colData(deseq.obj.pro)$conditions
          ctrst.pro[[cnt]] = get.contrast.res(deg=deg.pro,contrast=contrast,
                                              sig.thresh=sig.thresh,fc.thresh=fc.thresh,
                                              conditions=conditions, deseq=deseq.obj.pro)
          cnt = cnt + 1
        }
      }
      ind.up = list()
      ind.dn = list()
      for(jj in 1:length(ctrst.pro)){
        dat = ctrst.pro[[jj]]$contrast
        ind.up[[jj]] = which(dat$log2FoldChange > fc.thresh & dat$padj < sig.thresh)
        ind.dn[[jj]] = which(dat$log2FoldChange < (-1)*fc.thresh & dat$padj < sig.thresh)
      }
      indup.pro = Reduce(intersect, ind.up)
      inddn.pro = Reduce(intersect, ind.dn)
      
      # atac contrasts, intersect of timepoints
      tci = tcnd[[ii-1]]
      tcf = tcnd[[ii]]
      cnt = 1
      ctrst.atac = list()
      for(xx in 1:length(tci)){
        for(yy in 1:length(tcf)){
          contrast = c("conditions",tcf[yy],tci[xx])
          conditions = colData(deseq.obj.atac)$conditions
          ctrst.atac[[cnt]] = get.contrast.res(deg=deg.atac,contrast=contrast,
                                               sig.thresh=sig.thresh,fc.thresh=fc.thresh,
                                               conditions=conditions, deseq=deseq.obj.atac)
          cnt = cnt + 1
        }
      }
      ind.up = list()
      ind.dn = list()
      for(jj in 1:length(ctrst.atac)){
        dat = ctrst.atac[[jj]]$contrast
        ind.up[[jj]] = which(dat$log2FoldChange > fc.thresh & dat$padj < sig.thresh)
        ind.dn[[jj]] = which(dat$log2FoldChange < (-1)*fc.thresh & dat$padj < sig.thresh)
      }
      indup.atac = Reduce(intersect, ind.up)
      inddn.atac = Reduce(intersect, ind.dn)
    
      # separate significant tus and tfs
      upgene = rownames(ctrst.pro[[1]]$contrast)[indup.pro]
      dngene = rownames(ctrst.pro[[1]]$contrast)[inddn.pro]
      genes.up = get.gene.id( upgene )
      genes.dn = get.gene.id( dngene )
      tfs.up = get.sig.tfs(sig.genes=upgene, tf.comm.ann=tfclass.ann)
      tfs.dn = get.sig.tfs(sig.genes=dngene, tf.comm.ann=tfclass.ann)
      tus.up = genes.up[ !(genes.up %in% tfs.up$tf) ]
      tus.dn = genes.dn[ !(genes.dn %in% tfs.dn$tf) ]
      
      # get tf/tu bed coords and mean fold changes
      fc = list()
      for(jj in 1:length(ctrst.pro)){
        dat = ctrst.pro[[jj]]$contrast
        fc[[jj]] = dat$log2FoldChange
      }
      meanfc = colMeans( do.call(rbind, fc) )
      ctrst = ctrst.pro[[1]]$contrast
      ctrst$log2FoldChange = meanfc
      if(ii != length(tnet)){
        tu.bed.up = coord.bed.tu(sig=tfs.up$tf, bed0=bed0, ctrst=ctrst)
        tu.bed.dn = coord.bed.tu(sig=tfs.dn$tf, bed0=bed0, ctrst=ctrst)
      } 
      if (ii == length(tnet)){ # use tus only for the last time comparison
        tu.bed.up = coord.bed.tu(sig=genes.up, bed0=bed0, ctrst=ctrst)
        tu.bed.dn = coord.bed.tu(sig=genes.dn, bed0=bed0, ctrst=ctrst)
      }
      
      # get atac peak bed coords
      # filter for peaks associated with dynamic tfs
      upgene = rownames(ctrst.atac[[1]]$contrast)[indup.atac]
      dngene = rownames(ctrst.atac[[1]]$contrast)[inddn.atac]
      pk.bed.up = coord.bed.pk(pks=upgene)
      pk.bed.dn = coord.bed.pk(pks=dngene)
      if( full==FALSE & ii!=2 ){ # use atac peaks assoc with tfs
        pk.bed.up = pk.bed.up[pk.bed.up$gene %in% Etfre$idRE,]
        pk.bed.dn = pk.bed.dn[pk.bed.dn$gene %in% Etfre$idRE,]
      }
      
      # get re --> tu edges
      cond.atac1 = tnet[[ii-1]]
      cond.atac2 = tnet[[ii]]
      sizefac = sizeFactors(deseq.obj.atac)
      Eup=Edn=c()
      Eup = atac.edge.weight(bed.dat=tu.bed.up, peak.bed=pk.bed.up, summits=all.summits,
                             bigWig.path.atac=bigWig.path.atac, bed.dir=bed.dir,
                             atac.lim=atac.lim, tau=tau, sizefac=sizefac,
                             cond.atac1=cond.atac1, cond.atac2=cond.atac2)
      Edn = atac.edge.weight(bed.dat=tu.bed.dn, peak.bed=pk.bed.dn, summits=all.summits,
                             bigWig.path.atac=bigWig.path.atac, bed.dir=bed.dir,
                             atac.lim=atac.lim, tau=tau, sizefac=sizefac,
                             cond.atac1=cond.atac1, cond.atac2=cond.atac2)
      
      # aggregate results and annotate the net
      Eretf0 = rbind(Eup, Edn) %>% mutate(type="REtoTF",t1=paste0(tcnd[[ii-1]],collapse=","), t2=paste0(tcnd[[ii]],collapse=","))
      Etf = Eretf0[Eretf0$geneTU %in% tfclass.ann$tf,] %>% mutate(tu="tf")
      Etu = Eretf0[!(Eretf0$geneTU %in% tfclass.ann$tf),] %>% mutate(tu="tu")
      Eretf = Etf %>% mutate(tt=paste0(tcnd[ii-1],",",tcnd[ii]))
      Eall = rbind(Etf, Etu) %>% mutate(tt=paste0(t1,"_",t2))
      net.full = rbind(net.full, Eall)
      Eup=Edn=c()
      Etfre = c()
      
    } # get re --> tu/tf edges 
    ################################
    
  } # ii, net inference loop
  
  # filter output if needed
  ind.rm = grep("nodata",net.full[,1])
  if(length(ind.rm)>0){
    net = net.full[-ind.rm,]
  } else {
    net = net.full
  }
  
  # generate a strictly connected network and output results
  net[,c(2,3,6,7,9:13)] = apply(net[,c(2,3,6,7,9:13)],2,function(x){
    data.matrix(x) %>% as.numeric
  })
  res = force.connect(net)
  res[["fullnet"]] = net
  return(res)
  
} # infer.net.list


############################################################
## plot functions
############################################################


# normalize reads by size factor and average across replicates
sf.norm = function(cnts=NULL, sf=NULL, ann=NULL){
  
  # normalization
  out1 = cnts
  for(ii in 1:ncol(cnts)){
    ind = which(names(sf) == colnames(cnts)[ii])
    out1[,ii] = cnts[,ii] / sf[ind]
  }
  
  # average reps with time organization
  out2 = c()
  for(ii in 1:nrow(ann)){
    ind = grep(ann[ii,1], colnames(out1))
    av = apply(out1[,ind],1,mean)
    out2 = cbind(out2, av)
  }
  out2 = as.data.frame(out2)
  names(out2) = ann[,2]
  
  return(out2)
} # sf.norm


# function to generate smooth curves for each peak
read.curve.smooth = function(dat=NULL, len=200){
  x = colnames(dat) %>% as.numeric
  taxis = seq(min(x), max(x), length.out=len)
  paxis = matrix(0,nrow(dat),length(taxis))
  for(jj in 1:nrow(paxis)){
    datjj = cbind(x, t(dat[jj,])) %>% as.data.frame(stringsAsFactors=T)
    names(datjj) = c("time","reads")
    fit = loess(reads~time, datjj)
    paxis[jj,] = predict(fit, taxis)
  }
  rownames(paxis) = rownames(dat)
  colnames(paxis) = taxis
  return(paxis)
} # read.curve.smooth


# function to plot a set of curves
plot.curves = function(curves=NULL, id=NULL, plt=NULL, col="black"){
  tt = colnames(curves)
  tp = rep("l",nrow(curves))
  ly = rep(1,nrow(curves))
  av.pr = apply(curves,2,mean)
  if(plt==TRUE){
    matplot(tt,t(curves),type=tp,lty=ly,col=col,
            ylab=paste0(id),main="",xlab="time (hrs)")
    lines(tt,av.pr,col="darkgray",cex=4, lwd=3) 
  }
  out = data.frame(time=tt, mean=av.pr, stringsAsFactors=F)
  return(out)
} # plot.curves


# scale a curve to 0,1
scale01 = function(vec=NULL){
  out = (vec - min(vec)) / (max(vec) - min(vec))
  return(out)
} # scale01
