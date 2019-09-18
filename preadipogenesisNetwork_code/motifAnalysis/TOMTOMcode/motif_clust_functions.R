
##############################################################
## functions for subclustering
##############################################################

# function to iteratively subcluster data 
# until all clusters are below a specified level
# the input graph must contain four columns: from, to, pval, weight
subcluster = function(clusters=NULL, graph=NULL, max.clust=40, max.iter=100, change.thresh=3){
  
  # first reset data structure and community ids
  clustdat0 = cbind(clusters$membership, clusters$names)
  clustdat0 = as.data.frame(clustdat0, stringsAsFactors=F)
  names(clustdat0) = c("clust", "motif")
  clustdat0$clust = data.matrix(clustdat0$clust) %>% as.numeric
  
  # get motif counts per cluster
  cnts0 = get.motif.cnt(clustdat = clustdat0)
  cluster.data = cnts0$cluster.data
  cluster.counts = cnts0$cluster.counts
  maxcnt = max( cluster.counts$count )
  
  # loop through iterations of subclustering
  iter = 1
  while(maxcnt > max.clust){
    
    # isolate 'small' clusters and reset cluster identifiers
    ind = which(cluster.counts$count <= max.clust)
    cluster.data.small = cluster.data[cluster.data$clust %in% cluster.counts$clust[ind],]
    cnts.small = get.motif.cnt(clustdat = cluster.data.small)
    data.small = cnts.small$cluster.data
    counts.small = cnts.small$cluster.counts
    n.clust.small = max(counts.small$clust)
    
    # iteratively re-cluster 'large' clusters
    cluster.data.loop = data.small
    ind = which(cluster.counts$count > max.clust)
    for(ii in 1:length(ind)){
      
      # generate new clusters
      comm = cluster.counts$clust[ ind[ii] ]
      mots = cluster.data$motif[cluster.data$clust == comm]
      nmots = length(mots)
      ind.o = sapply(mots,function(x){which(graph$from==x)}) %>% unlist
      ind.i = sapply(mots,function(x){which(graph$to==x)}) %>% unlist
      threecol_filtered = graph[intersect(ind.o,ind.i),]
      g.filt = graph.data.frame(threecol_filtered, directed=F)
      g.filt = simplify(g.filt)
      comm.g.filt = fastgreedy.community(g.filt)
      data.big0 = cbind(clust=comm.g.filt$membership, motif=comm.g.filt$names)
      data.big0 = as.data.frame(data.big0, stringsAsFactors=F)
      data.big0$clust = data.matrix(data.big0$clust) %>% as.numeric
      
      # check if new clusters should be used or not
      sizes = sapply(unique(data.big0$clust),function(x){length(which(data.big0$clust==x))})
      nmots.new = max(sizes)
      if(nmots-nmots.new < change.thresh){
        data.big0 = cluster.data[cluster.data$clust == comm,]
      }
      
      # adjust cluster identifiers
      n.clust.big = length(unique(data.big0$clust))
      newids = (n.clust.small + 1) : (n.clust.small + n.clust.big)
      newkey = cbind(oldid=unique(data.big0$clust), newid=newids) %>% as.data.frame
      data.big0$clust = sapply(data.big0$clust,function(x){newkey$newid[which(newkey$oldid==x)]})
      
      # update the data
      n.clust.small = max(data.big0$clust)
      cluster.data.loop = rbind(cluster.data.loop, data.big0)
    } # ii, loop through large cluster
    
    # update the data for all re-clusterings
    cluster.counts.old = cluster.counts
    cnts0 = get.motif.cnt(clustdat = cluster.data.loop)
    cluster.data = cnts0$cluster.data
    cluster.counts = cnts0$cluster.counts
    dcount = maxcnt - max( cluster.counts$count )
    maxcnt = max( cluster.counts$count )
    if(dcount<change.thresh){
      cluster.counts = cluster.counts.old
      break
    }
    cat(paste0(iter, " round(s) of re-clustering completed, max = ", maxcnt,"\n"))
    if(iter==max.iter){break}
    if(maxcnt<max.clust){break}
    iter = iter + 1
    
  } # while loop - successive rounds of re-clustering
  
  return(cluster.data)
} # subcluster

# function to get motif counts for each cluster and re-name cluster with an ordinal set of integers
# assumes input dataframe with two columns: clust (numeric) and motif (char)
get.motif.cnt = function(clustdat){
  clustdat.new = clustdat
  clust0 = (unique(clustdat$clust))
  newid = c(1:length(clust0))
  key = cbind(clust0, newid) %>% as.data.frame
  clustdat.new$clust = sapply(clustdat.new$clust,function(x){key$newid[which(key$clust0==x)]})
  cnt = sapply(newid,function(x){length(which(clustdat.new$clust==x))})
  cnt.mat = cbind(clust=newid, count=cnt) %>% as.data.frame
  out = list(cluster.data=clustdat.new, cluster.counts=cnt.mat)
  return(out)
} # get.motif.cnt

##############################################################
## functions for finding the average motif
##############################################################

# function to find the offset between two PWMs
# input: query = query PWM that will serve as a reference
# input: target = target PWM
# input: max.offset = maximal ofset to consider (default 10)
# input: bkg = backround nucleotide frequencies, A,C,G,T (default all 0.25)
# output: the output in bp of the query relative to the target
#         a positive value indicates the query is shifted right
find.offset = function(query=NULL, target=NULL, max.offset=10, bkg=rep(0.25,4)){
  
  # reverse complement the target
  targetRC = revcomp(target)
  
  # run the loop to compute offsets for regular and reverse complement
  offset.data = offset.loop(query=query, target=target, max.offset=max.offset, bkg=bkg)
  offset.dataRC = offset.loop(query=query, target=targetRC, max.offset=max.offset, bkg=bkg)
  
  # combine data
  offset.data = offset.data %>% mutate(or = "reg")
  offset.dataRC = offset.dataRC %>% mutate(or = "rc")
  offset.data = rbind(offset.data, offset.dataRC)
  
  # determine the best offset based on max info and min distance
  info = offset.data$offset[ which(offset.data$meanInfo == max(offset.data$meanInfo))[1]  ] 
  dist = offset.data$offset[ which(offset.data$sumDist == min(offset.data$sumDist))[1]  ]
  orient = offset.data$or[ which(offset.data$sumDist == min(offset.data$sumDist))[1]  ]
  
  # output results
  out = list(offset.info=info, offset.dist=dist, orient=orient)
  return(out)
  
} # find.offset

# function to find the padded pwms for a given offset
# input: off = offset in number of base pairs
# input: target = target matrix, oriented with respect to the query matrix
# input: query = query matrix as a reference
# input: bkg = genomic background (ACGT, e.g., rep(0.25,4))
# output: list with query.pad matrix and target.pad matrix
get.padded.pwms = function(off=NULL, query=NULL, target=NULL, bkg=NULL){
  
  # negative offsets
  if(off < 0){
    padL.target = matrix(bkg, 4, abs(off))
    Ltx = ncol(target) - off
    if(Ltx > ncol(query)){
      padR.query = matrix(bkg, 4, (Ltx-ncol(query)))
      target.pad = cbind(padL.target, target)
      query.pad = cbind(query, padR.query)
    }
    if(Ltx < ncol(query)){
      padR.target = matrix(bkg, 4, (ncol(query)-Ltx))
      target.pad = cbind(padL.target, target, padR.target)
      query.pad = query
    }
    if(Ltx == ncol(query)){
      target.pad = cbind(padL.target, target)
      query.pad = query
    }
  } # off < 0
  
  # positive offsets
  if(off > 0){
    padL.query = matrix(bkg, 4, abs(off))
    Rtx = ncol(query) + off
    if(Rtx > ncol(target)){
      padR.target = matrix(bkg, 4, (Rtx - ncol(target)))
      target.pad = cbind(target, padR.target)
      query.pad = cbind(padL.query, query)
    }
    if(Rtx < ncol(target)){
      padR.query = matrix(bkg, 4, (ncol(target)-Rtx))
      target.pad = target
      query.pad = cbind(padL.query, query, padR.query)
    }
    if(Rtx == ncol(target)){
      target.pad = target
      query.pad = cbind(padL.query, query)
    }
  } # (off > 0)
  
  # zero offset
  if(off == 0){
    if(ncol(target) > ncol(query)){
      padR.query = matrix(bkg, 4, (ncol(target)-ncol(query)))
      target.pad = target
      query.pad = cbind(query, padR.query)
    }
    if(ncol(target) < ncol(query)){
      padR.target = matrix(bkg, 4, (ncol(query)-ncol(target)))
      target.pad = cbind(target, padR.target)
      query.pad = query
    }
    if(ncol(target) == ncol(query)){
      target.pad = target
      query.pad = query
    }
  } # (off == 0)
  colnames(query.pad) = colnames(target.pad) = 1:ncol(target.pad)
  rownames(target.pad) = rownames(query.pad) = rownames(query)
  return(list(query.pad=query.pad, target.pad=target.pad))
  
} # get.padded.pwms


# offset loop, called by find.offset()
offset.loop = function(query=NULL, target=NULL, max.offset=10, bkg=rep(0.25,4)){
  # specify the offsets
  offset = seq(-max.offset, max.offset, 1)
  
  # data storage: track the info or the mean matrix 
  # and summed absolute differences
  offset.data = matrix(0, length(offset), 3) %>% as.data.frame
  names(offset.data) = c("offset", "meanInfo", "sumDist")
  offset.data$offset = offset
  
  # loop through all offsets, perform calculations
  for(off in offset){
    
    ########################################
    ## generate padded query and target matrices for comparison
    padded = get.padded.pwms(off=off, query=query, target=target, bkg=bkg)
    query.pad = padded$query.pad
    target.pad = padded$target.pad
    
    ########################################
    ## perform comparisons for a given offset
    
    # method 1: total information of the mean matrix
    # consider only sites with info > 1 bit in the query
    # seqLogo(target.pad)
    av = 0.5 * (target.pad + query.pad)
    query.pad = add.pseudocount(PWM=query.pad, cnt=1e-6)
    av = add.pseudocount(PWM=av, cnt=1e-6)
    site.info.query = apply(query.pad,2,function(x){x %*% (log2(x)-log2(bkg)) })
    site.info.mean = apply(av,2,function(x){ x %*% (log2(x)-log2(bkg)) })
    ind.more1 = which(site.info.query > 1)
    if(length(ind.more1) > 1){
      av.info.sum = sum( site.info.mean[ind.more1]  )
    } else {
      av.info.sum = max( site.info.mean  )
      ind.more1 = which(site.info.query == max(site.info.query))
    }
    
    # method 2: total absolute difference between the matrices
    # consider only sites with info > 1 bit in the query
    diff = colSums( abs( target.pad - query.pad ) )
    diff.sum = sum( diff[ind.more1]  )
    
    # keep track of the data
    offset.data$meanInfo[offset.data$offset==off] = av.info.sum
    offset.data$sumDist[offset.data$offset==off] = diff.sum
    
  } # offset loop
  
  return(offset.data)
} # offset.loop


# function to get the reverse complement for a matrix
# this function assumes rows ACGT as input
revcomp = function(mat=NULL){
  targetRC = matrix(0,nrow(mat), ncol(mat))
  targetRC[1,] = rev(mat[4,])
  targetRC[2,] = rev(mat[3,])
  targetRC[3,] = rev(mat[2,])
  targetRC[4,] = rev(mat[1,])
  rownames(targetRC) = rownames(mat)
  return(targetRC)
}

# function to add a pseudo count to a PWM
add.pseudocount = function(PWM=NULL, cnt=1e-6){
  mat = PWM + cnt
  rr = nrow(PWM)
  mat = apply(mat,2,function(x){
    new = x
    ind = which(x==max(x))[1]
    new[ind] = x[ind] - rr*cnt
    return(new)
  })
  return(mat)
}

# function to remove low info sites at the margins of a pwm
# input: pwm.plt = probability matrix to be adjusted
# input: site.thr = site infor below which marginal sites are removed
# input: cnt = pseudocount number (recommend ~ 1e-6) 
#        for calling add.pseudocount{}
# input: bkg = genomic background (ACGT)
# output: pwm with low info sides removed
remove.lowinfo.margin = function(pwm.plt=NULL, site.thr=NULL, bkg=NULL, cnt=NULL){
  pwm.plt = add.pseudocount(PWM=pwm.plt, cnt=cnt)
  site.info = apply(pwm.plt,2,function(x){x %*% (log2(x)-log2(bkg))})
  ind.lo = which(site.info < site.thr)
  remL = 999
  remR = 999
  if(length(ind.lo) > 0){
    for(jj in 2:length(site.info)){
      if(site.info[jj-1] > site.thr){break}
      if(site.info[jj-1] < site.thr & site.info[jj] < site.thr){
        remL = unique(c(remL, jj-1, jj))
      }
    } # left loop
    for(jj in (length(site.info)-1):1){
      if(site.info[jj+1] > site.thr){break}
      if(site.info[jj+1] < site.thr & site.info[jj] < site.thr){
        remR = unique(c(remR, jj+1, jj))
      }
    } # right loop
  }
  col.rem = unique(c(remL,remR))
  pwm.plt = pwm.plt[,-col.rem]
  return(pwm.plt)
} # remove.lowinfo.margin


# get degree sum for each motif within each community
# use -log10(Eval) edge weights
get.degree = function(motif.community=NULL, threecol=NULL){
  motif.degree = c()
  for(ii in 1:length(unique(motif.community$community))){
    ind = which(motif.community$community == 
                  unique(motif.community$community)[ii])
    motifs = motif.community$motif[ind]
    deg = sapply(motifs,function(x){
      ind1 = which(threecol$from == x)
      ind2 = which(threecol$to == x)
      ind3 = which(threecol$from %in% motifs)
      ind4 = which(threecol$to %in% motifs)
      ind.deg = intersect( union(ind1,ind2), union(ind3,ind4) )
      new = c(x, unique(motif.community$community)[ii], 
              sum(threecol$weight[ind.deg]))
      return(new)
    }) %>% t
    motif.degree = rbind(motif.degree, deg)
  }
  motif.degree = motif.degree %>% as.data.frame(stringsAsFactors=F)
  rownames(motif.degree) = 1:nrow(motif.degree)
  names(motif.degree) = c("motif", "community", "degsum")
  motif.degree$degsum = as.numeric(motif.degree$degsum)
  return(motif.degree)
} # get.degree


# find the most connected motif for each community based on degree sum
# use motif length to decide for ties
# input: motif.degree = three col frame with names() = motif, community, degsum
#        motif is the motif identifier, community is an integer
#        degsum is the degree sum for each motif, see get.degree()
# output: a frame with a column of community numbers and the max connected motif
get.max.connect.motifs = function(motif.degree=NULL){
  max.connect.motifs = sapply(unique(motif.degree$community),function(x){
    ind = which(motif.degree$community == x)
    ind.max = which(motif.degree$degsum[ind] == max(motif.degree$degsum[ind]))
    motifs = motif.degree$motif[ind][ind.max]
    max.char = which(nchar(motifs) == max(nchar(motifs)))[1]
    return(motifs[max.char])
  })
  max.connect.motifs = cbind(unique(motif.degree$community), max.connect.motifs)
  max.connect.motifs = as.data.frame(max.connect.motifs, stringsAsFactors=F)
  names(max.connect.motifs) = c("community", "motif")
  max.connect.motifs$community = as.numeric(max.connect.motifs$community)
  max.connect.motifs = max.connect.motifs[order(max.connect.motifs$community),]
  return(max.connect.motifs)
} # get.max.connect.motifs

# find the most connected motif for each community based on number of connections
# use motif length to decide for ties
# output: a frame with a column of community numbers and the max connected motif
get.max.connect.motifs2 = function(motif.community=NULL, graph=NULL){
  
  max.connect.motifs = sapply(unique(motif.community$community),function(x){
    ind = which(motif.community$community == x)
    motifs = motif.community$motif[ind]
    nedge = c()
    for(ii in motifs){
      ind.i = which(graph$from == ii)
      ind.o = sapply(motifs,function(x){which(graph$to==x)}) %>% unlist
      ind = intersect(ind.i, ind.o)
      nedge = c(nedge, length(ind))
    } # ii
    names(nedge) = motifs
    maxcon = names(nedge)[which(nedge==max(nedge))[1]]
    return(maxcon)
  }) # sapply
  community = unique(motif.community$community)
  
  max.connect.motifs = cbind(community, max.connect.motifs)
  max.connect.motifs = as.data.frame(max.connect.motifs, stringsAsFactors=F)
  names(max.connect.motifs) = c("community", "motif")
  max.connect.motifs$community = as.numeric(max.connect.motifs$community)
  max.connect.motifs = max.connect.motifs[order(max.connect.motifs$community),]
  return(max.connect.motifs)
} # get.max.connect.motifs2


# function to get the de novo pwm
get.denovo.pwm = function(id=NULL){
  fname = paste0(pwm.dir,"/",id,".txt")
  mat = read.table(fname,sep="\t",header=F,stringsAsFactors=F) %>% t
  nas = apply(mat,1,function(x){all(is.na(x))})
  ind = which(nas==TRUE)
  if(length(ind)>0){mat = mat[-ind,]}
  rownames(mat) = c("A","C","G","T")
  colnames(mat) = c(1:ncol(mat))
  return(mat)
} # get.denovo.pwm


##############################################################
## functions for filtering motifs
##############################################################

# function to plot motif signature
logo.plt = function(pwm=NULL, title=NULL){
  ggplot() + geom_logo(pwm) + ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ylim(c(0,2.1))
} # logo.plt


# function to evaluate whether a matrix is dinucleotide
dinuc.metric = function(mat=NULL){
  out = c()
  if((ncol(mat) %% 2) == 0){mat = mat[,1:(ncol(mat)-1)]}
  sums = apply(mat,1,function(x){sum(x)})
  adjs = apply(mat,1,function(x){
    sum.adj = c()
    for(ii in 2:length(x)){
      sum.adj = max(sum.adj, c(x[ii]+x[ii-1]))
    }
    return(sum.adj)
  })
  for(ii in 1:(nrow(mat)-1)){
    for(jj in (ii+1):nrow(mat)){
      m1 = mat[ii,]
      m2 = rev( mat[jj,] )
      dot = m1 %*% m2
      out = sum(out, dot)
    } # jj
  } # ii
  if(max(adjs)>1){out=1}
  if(length(which(sums>0))>2){out=1}
  return(out)
} # dinuc.metric


# function to evaluate whether a matrix is ~ all GC
# higher values, more GC
gc.metric = function(mat=NULL){
  p1 = mat[1,] %*% mat[2,] # A dot C
  p2 = mat[3,] %*% mat[4,] # G dot T
  p3 = mat[3,] %*% mat[2,] # G dot C
  if( p3 > 1 & p1+p2 > 0 ){
    out = p3 / (p1 + p2)
  } else {out = 0}
  
  # gc = apply(mat,2,function(x){
  #   rownames(mat)[which(x==max(x))]
  # })
  # gclen = length(which(gc=="G" | gc=="C"))
  # if(gclen == ncol(mat)){out=0}
  
  return(out)
} # gc.metric


# function to evaluate whether a matrix is ~ all CT
# higher values, more CT
ct.metric = function(mat=NULL){
  p1 = mat[1,] %*% mat[2,] # A dot C
  p2 = mat[3,] %*% mat[4,] # G dot T
  p3 = mat[4,] %*% mat[2,] # T dot C
  if( p3 > 0.5 & p1+p2 > 0 ){
    out = p3 / (p1 + p2)
  } else {out = 0}
  return(out)
} # ct.metric


# function to evaluate whether a matrix is ~ all TG
# higher values, more T
tg.metric = function(mat=NULL){
  p1 = sum( mat[3,] + mat[4,] )
  p2 = sum( mat[1,] + mat[2,] )
  if( p1 > 1 & p2 > 0 ){
    out = p1 / p2
  } else {out = 0}
  return(out)
} # tg.metric

# function to evaluate info corraltions between de novo and TF
info.cor = function(mat=NULL, ref=NULL, bkg=NULL){
  mat = add.pseudocount(PWM=mat)
  ref = add.pseudocount(PWM=ref)
  info.mat = apply(mat,2,function(x){x %*% (log2(x)-log2(bkg))})
  info.ref = apply(ref,2,function(x){x %*% (log2(x)-log2(bkg))})
  ind.mat = which(info.mat > 0.15)
  ind.ref = which(info.ref > 0.15)
  ind = intersect(ind.mat, ind.ref)
  if(length(ind)==0){return(1); break}
  base.cor = cor(as.vector(mat[,ind]), as.vector(ref[,ind]))
  info.cor = cor(info.mat[ind], info.ref[ind])
  if(is.na(base.cor)==T){base.cor=1}
  if(is.na(info.cor)==T){info.cor=1}
  return(min(base.cor,info.cor))
} # info.cor

# function to evaluate the fraction of de novo bases high info in the TF
frac.denovo = function(mat=NULL, ref=NULL, bkg=NULL){
  mat = add.pseudocount(PWM=mat)
  ref = add.pseudocount(PWM=ref)
  info.mat = apply(mat,2,function(x){x %*% (log2(x)-log2(bkg))})
  info.ref = apply(ref,2,function(x){x %*% (log2(x)-log2(bkg))})
  ind.mat = which(info.mat > 0.5)
  ind.ref = which(info.ref[ind.mat] > 0.5)
  frac = length(ind.mat) / length(ind.ref)
  frac2 = length(which(info.mat>0.5 & info.ref>0.5))
  if(length(ind.mat)==0){frac=0}
  if(length(ind.mat)!=0 & length(ind.ref)==0){frac=10}
  if(length(ind.ref)>0 & length(frac2)>0 & frac2 / length(ind.ref) < 2){frac=0}
  return(frac)
} # frac.denovo

# function to determine if the motif is ~ all Gs (>0.9)
gg.metric = function(mat=NULL, bkg=NULL){
  mat = add.pseudocount(PWM=mat)
  info.mat = apply(mat,2,function(x){x %*% (log2(x)-log2(bkg))})
  top.prob = top.prob(mat)
  lenG = length(which(names(top.prob)=="G"))
  lenC = length(which(names(top.prob)=="C"))
  lenA = length(which(names(top.prob)=="A"))
  out = lenG / length(top.prob)
  if(lenG==length(top.prob)-1 & lenC==1){out=0} # sp1 filter
  if(lenG==length(top.prob)-1 & lenA==1){out=0} # HNRNPH2 filter
  return(out)
} # gg.metric

# function to determine if the motif is ~ all Cs (>0.9)
cc.metric = function(mat=NULL, bkg=NULL){
  mat = add.pseudocount(PWM=mat)
  info.mat = apply(mat,2,function(x){x %*% (log2(x)-log2(bkg))})
  top.prob = top.prob(mat)
  out = length(which(names(top.prob)=="C")) / length(top.prob)
  return(out)
} # cc.metric


# function to determine of a de novo motif has n identical bases in a row
# zero is returned if there are four bases in a row
samen = function(mat=NULL, n=4){
  out=1
  spann = apply(mat,1,function(x){
    res = 0
    for(ii in n:length(x)){
      new = sum ( x[(ii-n+1):ii] )
      res = max(res,new)
    }
    return(res)
  })
  ind.max1 = which(spann == max(spann))
  if( n-max(spann) < 0.03 & sum(spann[-ind.max1]) < n-1){out=0}
  if(spann[2]==1 & spann[3]>=n){out=1} # sp1 filter
  return(out)
}

# function to determine if only two bases are in the de novo motif
# zero is returned in this case
only2 = function(mat=NULL){
  out = 1
  bp = apply(mat,2,function(x){
    rownames(mat)[which(x!=0)]
  }) %>% unlist
  if(length(unique(bp))==2){out=0}
  return(out)
}

# function to evaluate the letter matches at single sites between two matrices
# the input matrices must already be aligned
site.prob = function(mat=NULL, ref=NULL, bkg=NULL){
  mat = add.pseudocount(PWM=mat)
  ref = add.pseudocount(PWM=ref)
  info.mat = apply(mat,2,function(x){x %*% (log2(x)-log2(bkg))})
  info.ref = apply(ref,2,function(x){x %*% (log2(x)-log2(bkg))})
  ind.mat = which(info.mat > 1)
  ind.ref = which(info.ref > 1)
  
  prop.mat = apply(mat,2,function(x){x[which(x==max(x))[1]]})
  prop.ref = apply(ref,2,function(x){x[which(x==max(x))[1]]})
  ind.pr.mat = which(prop.mat>0.8)
  ind.pr.ref = which(prop.ref>0.8)
  ind.mat = intersect(ind.mat, ind.pr.mat)
  ind.ref = intersect(ind.ref, ind.pr.ref)
  
  mat = mat[,intersect(ind.mat,ind.ref)] %>% as.data.frame
  ref = ref[,intersect(ind.mat,ind.ref)] %>% as.data.frame
  if(ncol(mat)==0 | ncol(ref)==0){return(0); break}
  ref.inds = which(top.prob(ref) > 0.7)
  top.prob.mat = top.prob(mat)[ref.inds]
  top.prob.ref = top.prob(ref)[ref.inds]
  matches = names(top.prob.mat) == names(top.prob.ref)
  #if(length(matches)<2){return(10); break}
  out = length(which(matches==FALSE)) / ncol(mat)
  return(out)
} # site.prob

# function to get the base with the highest probability at every position
# this returns a numeric vector of probabilities
# the names refer to the base with the highest probability
top.prob = function(ref=NULL){
  out1 = c()
  namen = rownames(ref)
  for(ii in 1:ncol(ref)){
    x = ref[,ii]
    ind = which(x==max(x))
    new = cbind(namen[ind], x[ind])
    out1 = rbind(out1, new)
  }
  out = data.matrix(out1[,2]) %>% as.numeric
  names(out) = out1[,1]
  return(out)
}


##############################################################
## functions for plotting motifs
##############################################################

# function for getting average denovo motifs and generating plots
# colnames(dat) = "community_3"  "community_11" "community_10" ...
# pwm.dir = directory for de novo motif pwms
# tomtom = tomtom data matrix
# names(max.connect.motifs) = "community" "motif"
# denovo.mean.plot = function(dat=NULL, pwm.dir=NULL, tomtom=NULL, 
#                             max.connect.motifs=NULL,
#                             bkg=NULL, pad.len=NULL,
#                             max.offset=NULL){
#   
#   plt = list()
#   pwm.av = list()
#   p=1
#   
#   dat0 = tomtom
#   
#   # loop through every motif and find the average motif
#   for(ii in 1:ncol(dat)){
#     
#     # get max connected motif for each community
#     comm = strsplit(colnames(dat)[ii], "community_")[[1]][2]
#     center = max.connect.motifs$motif[max.connect.motifs$community == comm]
#     
#     #if(center=="TSACC_deg6_1174"){break}
#     
#     # get all TFs from each cluster
#     TFs = dat[,ii]
#     TFs = TFs[which(TFs != 0)]
#     
#     # get a list of PWMs for cluster TFs
#     pwms = list()
#     for(jj in TFs){
#       fname = paste0(pwm.dir,"/",jj,".txt")
#       mat = read.table(fname,sep="\t",header=F,stringsAsFactors=F,fill=T) %>% t
#       nas = apply(mat,1,function(x){all(is.na(x))})
#       ind = which(nas==TRUE)
#       if(length(ind)>0){mat = mat[-ind,]}
#       rownames(mat) = c("A","C","G","T")
#       colnames(mat) = c(1:ncol(mat))
#       pwms[[jj]] = mat
#     }
#     
#     # get tomtom mapping information
#     ind1 = sapply(TFs,function(x){which(dat0$Query_ID == x)}) %>% unlist
#     ind2 = sapply(TFs,function(x){which(dat0$Target_ID == x)}) %>% unlist
#     tomtom = dat0[intersect(ind1, ind2),]
#     
#     # isolate the direct neighbors of the cluster center
#     query = center
#     tomtom0 = tomtom
#     tomtom = tomtom[tomtom$Query_ID==query,]
#     tomtom = tomtom[with(tomtom, order(Target_ID,p.value)),]
#     tomtom = tomtom[!duplicated(tomtom[,1:2]),]
#     if(nrow(tomtom)==0){
#       pwm.av[[center]] = pwms[[center]]
#       plt[[p]] = ggplot() + geom_logo(pwms[[center]]) + ggtitle("center alone") + ylim(0,2.1) +
#         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#               panel.background = element_blank(), axis.line = element_line(colour = "black"))
#       p=p+1
#       plt[[p]] = ggplot() + geom_logo(pwms[[center]]) + ggtitle("average") + ylim(0,2.1) +
#         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#               panel.background = element_blank(), axis.line = element_line(colour = "black"))
#       p=p+1
#       next
#     }
#     
#     # central pwm with padding
#     center.pwm = pwms[[query]]
#     pad = matrix(bkg, 4, pad.len)
#     center.pad = cbind(pad, center.pwm, pad)
#     colnames(center.pad) = 1:ncol(center.pad)
#     center.start = pad.len + 1
#     
#     # store neighbor PWMs relative to center
#     pwm.deg2 = list()
#     offsets0 = c()
#     for(jj in 1:nrow(tomtom)){
#       tf = tomtom$Target_ID[jj]
#       pwm.ii = pwms[[tf]]
#       offset = find.offset(query=center.pwm, target=pwm.ii, max.offset=max.offset, bkg=bkg)
#       orient = offset$orient
#       if(orient == "rc"){pwm.ii = revcomp(pwm.ii)}
#       offset = offset$offset.dist
#       new = c(tf, offset)
#       offsets0 = rbind(offsets0, new)
#       padL = matrix(bkg, 4, pad.len-offset)
#       padR = matrix(bkg, 4, ncol(center.pad)-(ncol(padL)+ncol(pwm.ii)))
#       pwm.pad = cbind(padL, pwm.ii, padR)
#       colnames(pwm.pad) = 1:ncol(pwm.pad)
#       pwm.deg2[[tf]] = pwm.pad
#     }
#     offsets0 = as.data.frame(offsets0,stringsAsFactors=F)
#     names(offsets0) = c("motif","offset")
#     offsets0$offset = as.numeric(offsets0$offset)
#     
#     # get pwms for all other motifs in the cluster
#     higher = TFs[!(TFs %in% names(pwm.deg2))]
#     while(length(higher) > 0){
#       ind.un = sapply(higher,function(x){which(tomtom0$Query_ID==x | tomtom0$Target_ID==x)}) %>% unlist
#       ind.pw = sapply(names(pwm.deg2),function(x){which(tomtom0$Query_ID==x | tomtom0$Target_ID==x)}) %>% unlist
#       ind = intersect(ind.un, ind.pw)[1]
#       facs = c( tomtom0$Query_ID[ind], tomtom0$Target_ID[ind] )
#       unmatched = higher[higher %in% facs]
#       matched = names(pwm.deg2)[names(pwm.deg2) %in% facs]
#       #if(unmatched=="CCHCWCCYGCTGSR_dyn6_1500"){break}
#       pwm.un = pwms[[unmatched]]
#       pwm.mt = pwms[[matched]]
#       offset = find.offset(query=pwm.mt, target=pwm.un, max.offset=max.offset, bkg=bkg)
#       orient = offset$orient
#       if(orient == "rc"){pwm.un = revcomp(pwm.un)}
#       offset = offset$offset.dist + offsets0$offset[offsets0$motif==matched]
#       new = c(unmatched, offset)
#       offsets0 = rbind(offsets0, new)
#       offsets0$offset = as.numeric(offsets0$offset)
#       padL = matrix(bkg, 4, pad.len-offset)
#       padR = matrix(bkg, 4, ncol(center.pad)-(ncol(padL)+ncol(pwm.un)))
#       pwm.pad = cbind(padL, pwm.un, padR)
#       colnames(pwm.pad) = 1:ncol(pwm.pad)
#       pwm.deg2[[unmatched]] = pwm.pad
#       higher = higher[-which(higher==unmatched)]
#     }
#     pwm.deg2[[center]] = center.pad
#     
#     # average
#     sum = matrix(0, nrow(center.pad), ncol(center.pad))
#     for(jj in 1:length(pwm.deg2)){
#       sum = sum + pwm.deg2[[jj]]
#     }
#     av = (1/length(pwm.deg2)) * sum
#     
#     # remove pad
#     pad = apply(av,2,function(x){x==bkg})
#     pad = apply(pad,2,function(x){all(x==TRUE)})
#     ind.keep = which(pad == FALSE)
#     pwm.plt = av[,ind.keep]
#     
#     # remove low info sites from the average pwm
#     pwm.plt = remove.lowinfo.margin(pwm.plt=pwm.plt, site.thr=site.thr, bkg=bkg, cnt=1e-6)
#     
#     # align the center motif to the average pwm
#     res = find.offset(query=center.pwm, target=pwm.plt, max.offset=10, bkg=bkg)
#     offset = res$offset.dist
#     orient = res$orient
#     if(orient == "rc"){pwm.plt = revcomp(pwm.plt)}
#     padded = get.padded.pwms(off=offset, query=center.pwm, target=pwm.plt, bkg=bkg)
#     center.plot = padded$query.pad
#     mean.plot = padded$target.pad
#     
#     # save plots 
#     plt[[p]] = ggplot() + geom_logo(center.plot) + ggtitle("center") + ylim(0,2.1) +
#       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#             panel.background = element_blank(), axis.line = element_line(colour = "black"))
#     p=p+1
#     plt[[p]] = ggplot() + geom_logo(mean.plot) + ggtitle("average") + ylim(0,2.1) +
#       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#             panel.background = element_blank(), axis.line = element_line(colour = "black"))
#     p=p+1
#     
#     # save the matrix for downstream analysis
#     pwm.av[[center]] = mean.plot
#     
#   }
#   out = list(plt=plt, pwm.av=pwm.av)
#   return(out)
# } # denovo.mean.plot


# function for getting average denovo motifs and generating plots
# colnames(dat) = "community_3"  "community_11" "community_10" ...
# pwm.dir = directory for de novo motif pwms
# tomtom = tomtom data matrix
# names(max.connect.motifs) = "community" "motif"
denovo.mean.plot.TF = function(dat=NULL, pwm.dir=NULL, tomtom=NULL, 
                            max.connect.motifs=NULL,
                            bkg=NULL, pad.len=NULL,
                            max.offset=NULL){
  
  dat0 = tomtom

  # loop through every TF and find the average motif
  plt = list()
  pwm.av = list()
  p=1
  for(ii in 1:ncol(dat)){
    
    # get max connected motif for each community
    comm = strsplit(colnames(dat)[ii], "community_")[[1]][2]
    center = max.connect.motifs$motif[max.connect.motifs$community == comm]
    id = motif.community$ids.orig[which(motif.community$motif == center)]
    #if(id=="FOXO1_MOUSE.H10MO.C"){break}
    
    # get all TFs from each cluster
    TFs = dat[,ii]
    TFs = TFs[which(TFs != 0)]
    
    # get a list of PWMs for cluster TFs
    pwms = list()
    for(jj in TFs){
      #fname = paste0(pwm.dir,"/","PWM_",jj,".txt")
      fname = paste0(pwm.dir,"/",jj,".txt")
      mat = read.table(fname,sep="\t",header=F,stringsAsFactors=F,fill=T) %>% t
      nas = apply(mat,1,function(x){all(is.na(x))})
      ind = which(nas==TRUE)
      if(length(ind)>0){mat = mat[-ind,]}
      rownames(mat) = c("A","C","G","T")
      colnames(mat) = c(1:ncol(mat))
      pwms[[jj]] = mat
    }
    
    # get tomtom mapping information
    ind1 = sapply(TFs,function(x){which(dat0$Query_ID == x)}) %>% unlist
    ind2 = sapply(TFs,function(x){which(dat0$Target_ID == x)}) %>% unlist
    tomtom = dat0[intersect(ind1, ind2),]
    
    # isolate the direct neighbors of the cluster center
    query = center
    tomtom0 = tomtom
    tomtom = tomtom[tomtom$Query_ID==query,]
    tomtom = tomtom[with(tomtom, order(Target_ID,p.value)),]
    tomtom = tomtom[!duplicated(tomtom[,1:2]),]
    if(nrow(tomtom)==0){
      pwm.av[[center]] = pwms[[center]]
      plt[[p]] = ggplot() + geom_logo(pwms[[center]]) + ggtitle(id) + ylim(0,2.1) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
      p=p+1
      plt[[p]] = ggplot() + geom_logo(pwms[[center]]) + ggtitle("loner") + ylim(0,2.1) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
      p=p+1
      next
    }
    
    # central pwm with padding
    center.pwm = pwms[[query]]
    pad = matrix(bkg, 4, pad.len)
    center.pad = cbind(pad, center.pwm, pad)
    colnames(center.pad) = 1:ncol(center.pad)
    center.start = pad.len + 1
    
    # neighbors to the center
    pwm.deg2 = list()
    offsets0 = c()
    rc0 = c()
    for(jj in 1:nrow(tomtom)){
      tf = tomtom$Target_ID[jj]
      pwm.ii = pwms[[tf]]
      offset = find.offset(query=center.pwm, target=pwm.ii, max.offset=max.offset, bkg=bkg)
      orient = offset$orient
      rc0 = c(rc0, "reg")
      if(orient == "rc"){pwm.ii = revcomp(pwm.ii); rc0[length(rc0)]="rc"}
      offset = offset$offset.dist
      new = c(tf, offset)
      offsets0 = rbind(offsets0, new)
      padL = matrix(bkg, 4, pad.len-offset)
      padR = matrix(bkg, 4, ncol(center.pad)-(ncol(padL)+ncol(pwm.ii)))
      pwm.pad = cbind(padL, pwm.ii, padR)
      colnames(pwm.pad) = 1:ncol(pwm.pad)
      pwm.deg2[[tf]] = pwm.pad
    }
    offsets0 = as.data.frame(offsets0,stringsAsFactors=F)
    names(offsets0) = c("motif","offset")
    offsets0$offset = as.numeric(offsets0$offset)
    rc0 = cbind(offsets0$motif,rc0)
    rc0 = as.data.frame(rc0,stringsAsFactors=F)
    names(rc0) = c("motif","rc")
    
    # get pwms for all other motifs in the cluster
    higher = TFs[!(TFs %in% names(pwm.deg2))]
    xx=1
    while(length(higher) > 0){
      ind.un = sapply(higher,function(x){which(tomtom0$Query_ID==x | tomtom0$Target_ID==x)}) %>% unlist
      ind.pw = sapply(names(pwm.deg2),function(x){which(tomtom0$Query_ID==x | tomtom0$Target_ID==x)}) %>% unlist
      ind = intersect(ind.un, ind.pw)[1]
      facs = c( tomtom0$Query_ID[ind], tomtom0$Target_ID[ind] )
      unmatched = higher[higher %in% facs]
      matched = names(pwm.deg2)[names(pwm.deg2) %in% facs]
      #if(unmatched == "motif6814"){break}
      pwm.un = pwms[[unmatched]]
      pwm.mt = pwms[[matched]]
      rc.mat = rc0$rc[rc0$motif==matched]
      if(rc.mat=="rc"){
        pwm.mt = revcomp(pwm.mt)
      }
      offset = find.offset(query=pwm.mt, target=pwm.un, max.offset=max.offset, bkg=bkg)
      orient = offset$orient
      rc1 = c(unmatched, "reg")
      if(orient == "rc"){
        pwm.un = revcomp(pwm.un)
      }
      # if(rc.mat=="rc" & orient=="reg"){
      #   rc1 = c(unmatched, "rc")
      #   xx=xx+1
      # }
      if(rc.mat=="reg" & orient=="rc"){
        rc1 = c(unmatched, "rc")
      }
      rc0 = rbind(rc0, rc1)
      offset = offset$offset.dist + offsets0$offset[offsets0$motif==matched]
      new = c(unmatched, offset)
      offsets0 = rbind(offsets0, new)
      offsets0$offset = as.numeric(offsets0$offset)
      padL = matrix(bkg, 4, pad.len-offset)
      padR = matrix(bkg, 4, ncol(center.pad)-(ncol(padL)+ncol(pwm.un)))
      pwm.pad = cbind(padL, pwm.un, padR)
      colnames(pwm.pad) = 1:ncol(pwm.pad)
      pwm.deg2[[unmatched]] = pwm.pad
      higher = higher[-which(higher==unmatched)]
    }
    pwm.deg2[[center]] = center.pad
    
    #print(c(id,xx))
    #if(xx != 1){print(c(id,xx))}

    # average
    sum = matrix(0, nrow(center.pad), ncol(center.pad))
    for(jj in 1:length(pwm.deg2)){
      sum = sum + pwm.deg2[[jj]]
    }
    av = (1/length(pwm.deg2)) * sum
    
    # remove pad
    pad = apply(av,2,function(x){x==bkg})
    pad = apply(pad,2,function(x){all(x==TRUE)})
    ind.keep = which(pad == FALSE)
    pwm.plt = av[,ind.keep]
    
    # remove low info sites from the average pwm
    pwm.plt = remove.lowinfo.margin(pwm.plt=pwm.plt, site.thr=site.thr, bkg=bkg, cnt=1e-6)
    
    # align the center motif to the average pwm
    res = find.offset(query=center.pwm, target=pwm.plt, max.offset=10, bkg=bkg)
    offset = res$offset.dist
    orient = res$orient
    if(orient == "rc"){pwm.plt = revcomp(pwm.plt)}
    padded = get.padded.pwms(off=offset, query=center.pwm, target=pwm.plt, bkg=bkg)
    center.plot = padded$query.pad
    mean.plot = padded$target.pad
    
    # save plots 
    plt[[p]] = ggplot() + geom_logo(center.plot) + ggtitle(id) + ylim(0,2.1) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    p=p+1
    plt[[p]] = ggplot() + geom_logo(mean.plot) + ggtitle("average") + ylim(0,2.1) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    p=p+1
    
    # save the matrix for downstream analysis
    out = remove.lowinfo.margin(pwm.plt=mean.plot, site.thr=site.thr, bkg=bkg, cnt=1e-6)
    pwm.av[[center]] = out
    
  }
  out = list(plt=plt, pwm.av=pwm.av)
  return(out)
} # denovo.mean.plot.TF




# test plot
# p=1
# for(ii in 1:length(pwm.deg2)){
#   pwm = pwm.deg2[[ii]][,58:72]
#   id = names(pwm.deg2)[ii]
#   plt[[p]] = ggplot() + geom_logo(pwm) + ggtitle(id) + ylim(0,2.1) +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           panel.background = element_blank(), axis.line = element_line(colour = "black"))
#   p=p+1
# }
# pdf("TFexample3.pdf", onefile = TRUE, height=9, width=4.5)
# marrangeGrob(grobs=plt, nrow=6, ncol=1, top=NULL)
# dev.off()

# # ex info plt
# cnt=1e-6
# bkg = rep(0.25,4)
# topfr = seq(0.25,1,0.01)
# info = c()
# for(ii in topfr){
#   other = (1 - ii) / 3
#   x = c(ii, other, other, other) + cnt
#   imax = which(x==max(x))[1]
#   x[imax] = x[imax] - length(x)*cnt
#   new = x %*% (log2(x)-log2(bkg))
#   info = c(info, new)
# }
# plot(topfr, info)


