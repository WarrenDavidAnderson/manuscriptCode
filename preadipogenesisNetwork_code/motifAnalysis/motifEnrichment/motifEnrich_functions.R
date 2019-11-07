
############################################################
# motif enrichment plotting functions
############################################################

# function to prioritize data for plotting for each factor
# select the time comparison with the largest peak index
# input: enrich.dat = data frame with enrichment data (e.g., enrich.inc)
# input: comp.mode = c("inc","dec","un")
enrich.plot.select = function(enrich.dat=NULL, comp.mode=NULL){
  out = c()
  for(ii in unique(enrich.dat$id)){
    ind = which(enrich.dat$id == ii)
    dat = enrich.dat[ind,]
    if(comp.mode=="inc"){
      dat = dat[with(dat, order(-pkindINC, -varINC)),] 
    }
    if(comp.mode=="dec"){
      dat = dat[with(dat, order(-pkindDEC, -varDEC)),] 
    }
    if(comp.mode=="un"){
      dat = dat[with(dat, order(-pkindUN, -varUN)),] 
    }
    out = rbind(out, dat[1,])
  }
  return(out)
} # enrich.plot.select

# function for plotting enrichments
# input: comp.mode = c("inc","dec","un")
# input: enrich.dat = data frame with enrichment data (e.g., enrich.inc)
# input: res.pairs = differential peak data
# input: sig.thresh = FDR threshold for the differential peak analysis
# input: fc.thresh = fold change threshold for dynamics peaks
# input: sig.un = significance threshold for designating unchanged peaks
# input: fc.un = fold change threshold for designating unchanged peaks
# input: half.win = window around the summits for enrichment analysis
# input: step = step size for density enrichment plots (bp)
# input: fimo.bigWig.path = path to fimo bigWig files
# input: pwm.list = list of PWMs
plot.enrich = function(comp.mode=NULL,fname=NULL,mfrow=NULL,
                       enrich.dat=NULL, res.pairs=NULL, 
                       sig.thresh=NULL, fc.thresh=NULL,
                       sig.un=NULL, fc.un=NULL, type=NULL,
                       half.win=NULL, step=NULL, delbp = 60,
                       fimo.bigWig.path=NULL, loess.span = 0.3,
                       plt.trace=F, plt.motif=F, pwm.list=NULL){
  
  diff.filtered = c()
  
  # set plot file
  pdf(fname)
  par(mfrow=mfrow)
  
  # loop through all entries
  for(ii in 1:nrow(enrich.dat)){
    
    # specify data
    dat = enrich.dat[ii,]
    tfid = dat$id
    
    # text for plotting
    if(comp.mode=="inc"){
      per.dat1 = signif(dat$INCminusUN,2)
      per.dat2 = signif(dat$INCminusDEC,2)
      var.dat1 = signif(dat$varINC,3)
      var.dat2 = signif(dat$cvINC,3)
      per.diff = max(dat$INCminusUN, dat$INCminusDEC)
    }
    if(comp.mode=="dec"){
      per.dat1 = signif(dat$DECminusUN,2)
      per.dat2 = signif(dat$DECminusINC,2)
      var.dat1 = signif(dat$varDEC,3)
      var.dat2 = signif(dat$cvDEC,3)
      per.diff = max(dat$DECminusUN, dat$DECminusINC)
    }
    if(comp.mode=="un"){
      per.dat1 = signif(-dat$INCminusUN,2)
      per.dat2 = signif(-dat$DECminusUN,2)
      var.dat1 = signif(dat$varUN,3)
      var.dat2 = signif(dat$cvUN,3)
      per.diff = max(-dat$INCminusUN, -dat$DECminusUN)
    }
    
    # comparison
    comp = dat$tcomp

    # differential atac peaks
    if(type == "deg"){
      
      # peak data indices
      res_ii = res.pairs[[comp]]
      ind.up = which(res_ii$padj<sig.thresh & res_ii$log2FoldChange>fc.thresh)
      ind.dn = which(res_ii$padj<sig.thresh & res_ii$log2FoldChange<(-1)*fc.thresh)
      ind.un = which(res_ii$padj>sig.un & abs(res_ii$log2FoldChange)<fc.un)
      if(length(ind.up)==0){next}
      
      # get read mapping coords with half.win around each summit
      up.ind = sapply(rownames(res_ii)[ind.up],function(x){which(all.peaks==x)})
      dn.ind = sapply(rownames(res_ii)[ind.dn],function(x){which(all.peaks==x)})
      un.ind = sapply(rownames(res_ii)[ind.un],function(x){which(all.peaks==x)})
      bed.up = bed.window(summits=all.summits[up.ind], half.win=half.win)
      bed.dn = bed.window(summits=all.summits[dn.ind], half.win=half.win)
      bed.un = bed.window(summits=all.summits[un.ind], half.win=half.win)
    }
    
    # dreg peaks
    if(type == "dreg"){
      
      # peak data indices
      res_ii = as.data.frame(res.pairs[[comp]])
      rownamen = rownames(res_ii)
      chrs = sapply(rownames(res_ii),function(x)strsplit(x,":")[[1]][1])
      res_ii = res_ii %>% mutate(chrs=chrs)
      rownames(res_ii) = rownamen
      keep = paste0("chr",1:19)
      res_ii = res_ii[res_ii$chrs %in% keep,]
      ind.up = which(res_ii$padj<sig.thresh & res_ii$log2FoldChange>fc.thresh)
      ind.dn = which(res_ii$padj<sig.thresh & res_ii$log2FoldChange<(-1)*fc.thresh)
      ind.un = which(res_ii$padj>sig.un & abs(res_ii$log2FoldChange)<fc.un)
      if(length(ind.up)==0){next}
      
      # get read mapping coords with half.win around each summit
      up.ind = sapply(rownames(res_ii)[ind.up],function(x){which(bed.map$id==x)}) %>% unlist()
      dn.ind = sapply(rownames(res_ii)[ind.dn],function(x){which(bed.map$id==x)}) %>% unlist()
      un.ind = sapply(rownames(res_ii)[ind.un],function(x){which(bed.map$id==x)}) %>% unlist()
      bed.up = bed.window2(summits=bed.map$center[up.ind], half.win=half.win, chr=bed.map$chr[up.ind])
      bed.dn = bed.window2(summits=bed.map$center[dn.ind], half.win=half.win, chr=bed.map$chr[dn.ind])
      bed.un = bed.window2(summits=bed.map$center[un.ind], half.win=half.win, chr=bed.map$chr[un.ind])
    }
    
    # factor bigWig file
    fname = list.files(path=fimo.bigWig.path, pattern=paste0(dat$id,".bigWig"))
    bw.fimo = paste0(fimo.bigWig.path,"/",fname)
    bw.fimo = load.bigWig(bw.fimo)
    
    # bp range centered on zero
    bps = seq(-half.win, half.win, length.out=2*half.win/step)
    
    # matching peak data with the motif 
    updat = motif.map(bigwig=bw.fimo,coords=bed.up,step=step,bps=bps)
    dndat = motif.map(bigwig=bw.fimo,coords=bed.dn,step=step,bps=bps)
    undat = motif.map(bigwig=bw.fimo,coords=bed.un,step=step,bps=bps)
    
    # compute the baseline adjusted integrals
    indL = which(updat$bp < -(half.win-2*delbp))
    indM = which(updat$bp >= -delbp/2 & updat$bp <= delbp/2)
    indR = which(updat$bp > (half.win-2*delbp))
    upnorm = updat$mean - mean(c(updat$mean[indL],updat$mean[indR]))
    dnnorm = dndat$mean - mean(c(dndat$mean[indL],dndat$mean[indR]))
    unnorm = undat$mean - mean(c(undat$mean[indL],undat$mean[indR]))
    x = seq(0,step*length(upnorm)-step,by=step)
    integup = trapz(x,upnorm)
    integdn = trapz(x,dnnorm)
    integun = trapz(x,unnorm)
    
    # compute peak index and pk density difference
    pk.up = mean(updat$mean[indM]) / mean(c(updat$mean[indL],updat$mean[indR]))
    pk.dn = mean(dndat$mean[indM]) / mean(c(dndat$mean[indL],dndat$mean[indR]))
    pk.un = mean(undat$mean[indM]) / mean(c(undat$mean[indL],undat$mean[indR]))
    if(comp.mode=="inc"){
      pk.index = pk.up
      integ = integup
    }
    if(comp.mode=="dec"){
      pk.index = pk.dn
      integ = integdn
    }
    if(comp.mode=="un"){
      pk.index = pk.un
      integ = integun
    }
    pk.up = format(pk.up,digits=3)
    pk.dn = format(pk.dn,digits=3)
    pk.un = format(pk.un,digits=3)
    pk.index = format(pk.index,digits=3)
    
    # plot the peak traces
    fnt=1
    motif = dat$TF
    all = c(updat$mean, dndat$mean, undat$mean)
    dd = max(all) - min(all)
    ylim = c( min(all)-0.1*dd, max(all)+0.1*dd )
    ylab = paste0(motif, " density")
    tt = strsplit(comp,"_")[[1]]
    main = paste0(tt[1], "hr vs ",tt[2], "hr, pk.index = ",pk.index, "\n percent diff: ",
                  per.dat1,"%, ",per.dat2,"% \n var: ",var.dat1,", ",var.dat2,
                  "\n npeaks: ",dat$npeakINC,", ",dat$npeakDEC,", ",dat$npeakUN)
    plot(updat$bp, updat$mean, type="l", ylim=ylim, 
         ylab=ylab, xlab="distance from summit (bp)",
         lwd=2, col="red", main=main,cex.lab=fnt,cex.main=fnt)
    lines(dndat$bp, dndat$mean,lwd=2,col="blue")
    lines(undat$bp, undat$mean,lwd=2,col="darkgray")
    
    # plot smoothed curves
    main = paste0("peak indices = ",pk.up,", ",pk.dn,", ",pk.un,
                  "\n max % diff = ",format(per.diff,digits=3),
                  "\n integral = ", format(integ,digits=3))
    plot(updat$bp, updat$mean, type="l", ylim=ylim, 
         ylab="lowess smoothed", xlab="distance from summit (bp)",
         lwd=2, col="white", main=main,cex.lab=fnt,cex.main=fnt)
    fit.up = loess(mean~bp, updat, span=loess.span)
    fit.dn = loess(mean~bp, dndat, span=loess.span)
    fit.un = loess(mean~bp, undat, span=loess.span)
    xaxis = seq(fit.up$x[1], fit.up$x[length(fit.up$x)], length.out=200)
    plt.up = predict(fit.up, xaxis)
    plt.dn = predict(fit.dn, xaxis)
    plt.un = predict(fit.un, xaxis)
    lines(xaxis, plt.up, lwd=2, col="red")
    lines(xaxis, plt.dn, lwd=2, col="blue")
    lines(xaxis, plt.un, lwd=2, col="darkgray")

    # get peak data in bed format
    up.bed = get.map.from.res( rownames(res_ii[ind.up,]) )
    dn.bed = get.map.from.res( rownames(res_ii[ind.dn,]) )
    un.bed = get.map.from.res( rownames(res_ii[ind.un,]) )
    up.bed[,2:3] = apply(up.bed[,2:3],2,function(x){data.matrix(x) %>% as.numeric})
    dn.bed[,2:3] = apply(dn.bed[,2:3],2,function(x){data.matrix(x) %>% as.numeric})
    un.bed[,2:3] = apply(un.bed[,2:3],2,function(x){data.matrix(x) %>% as.numeric})
    
    # load bigwig factor mapping data into peak coordinates
    # number of TFBS in each peak region
    up = bed.region.bpQuery.bigWig(bw.fimo, up.bed)
    dn = bed.region.bpQuery.bigWig(bw.fimo, dn.bed)
    un = bed.region.bpQuery.bigWig(bw.fimo, un.bed)
    names(up) = rownames(up.bed)
    names(dn) = rownames(dn.bed)
    names(un) = rownames(un.bed)
    
    # match atac peaks with motif coordinates
    # get numbers of peaks with >0 motifs
    inc.with = length(which(up > 0))
    inc.without = length(which(up == 0))
    dec.with = length(which(dn > 0))
    dec.without = length(which(dn == 0))
    nodif.with = length(which(un > 0))
    nodif.without = length(which(un == 0))
    result = cbind(c(nodif.with, nodif.without),
                   c(inc.with, inc.without),
                   c(dec.with, dec.without))
    colnames(result) <- c("Un", "Inc", "Dec")
    rownames(result) <- c("with motif", "without motif")
    resultn = apply(result,2,function(x){100*x/(sum(x))})
    
    # chi square test result
    chi = chisq.test(result)
    
    # plot data
    pv = format(chi$p.value,digits=3)
    main = paste0(motif,"\n p = ",pv,
          "\n countsUID: ",result[1,1],", ",result[1,2],", ",result[1,3],
          "\n %: ",resultn[1,1]%>%format(digits=3),", ",
          resultn[1,2]%>%format(digits=3),", ",resultn[1,3]%>%format(digits=3))
    barplot(resultn[1,], col = c("gray","red","blue"), 
            main=main, cex.main=fnt, ylim=c(0,100),
            ylab="percent", xlab="peak class",
            cex.lab=fnt,cex.names=fnt)
    barplot(resultn[2,], col = c("white"), 
            axes=F, axisnames=F,
            offset = resultn[1,], add=T)
    
    print(ii / nrow(enrich.dat))
    
  } # for loop
  dev.off()
} # plot.enrich


# function for plotting enrichments for dynamic cluster peaks
plot.enrich.dyn = function(comp.mode=NULL,fname=NULL,mfrow=c(3,3),
                           enrich.dat=NULL, 
                           stem.res=NULL,
                           out.uns=NULL,
                           all.peaks=NULL,
                           all.summits=NULL,
                           type=NULL,
                           half.win=NULL, step=NULL,
                           fimo.bigWig.path=NULL){
  
  diff.filtered = c()
  
  # set plot file
  pdf(fname)
  par(mfrow=mfrow)
  
  # loop through all entries
  for(ii in 1:nrow(enrich.dat)){
    
    # specify data
    dat = enrich.dat[ii,]
    tfid = dat$id
    
    # text for plotting

    if(comp.mode=="dyn"){
      per.dat1 = signif(dat$DYNminusUNC,2)
      var.dat1 = signif(dat$varDYN,3)
      per.diff = dat$DYNminusUNC
    }
    
    if(comp.mode=="unc"){
      per.dat1 = signif(dat$UNCminusDYN,2)
      var.dat1 = signif(dat$varUNC,3)
      per.diff = dat$UNCminusDYN
    }
    
    # comparison
    comp = dat$Profile
    
    # differential atac peaks
    # get bed data for dynamic cluster peaks
    # get bed data for non-dynamic peaks
    ind = which(stem.res$Profile == comp) 
    pkdat = stem.res$gene_symbol[ind]
    dyn.bed = get.map.from.res( pkdat )
    len.dyn = length(ind)
    pk.un = out.uns$peakID[sample.int( nrow(out.uns), min(nrow(out.uns), length(pkdat)) )]
    unc.bed = get.map.from.res( pk.un )
    len.unc = length(pk.un)
    
    # get read mapping coords with half.win around each summit
    dyn.ind = sapply(pkdat,function(x){which(all.peaks==x)}) %>% unlist()
    unc.ind = sapply(pk.un,function(x){which(all.peaks==x)}) %>% unlist()
    bed.dyn = bed.window(summits=all.summits[dyn.ind], half.win=half.win)
    bed.unc = bed.window(summits=all.summits[unc.ind], half.win=half.win)
    
    # factor bigWig file
    fname = list.files(path=fimo.bigWig.path, pattern=paste0(dat$id,".bigWig"))
    bw.fimo = paste0(fimo.bigWig.path,"/",fname)
    bw.fimo = load.bigWig(bw.fimo)
    
    # bp range centered on zero
    bps = seq(-half.win, half.win, length.out=2*half.win/step)
    
    # matching peak data with the motif 
    dyndat = motif.map(bigwig=bw.fimo,coords=bed.dyn,step=step,bps=bps)
    uncdat = motif.map(bigwig=bw.fimo,coords=bed.unc,step=step,bps=bps)
    
    # compute the baseline adjusted integrals
    indL = which(dyndat$bp < -(half.win-2*delbp))
    indM = which(dyndat$bp >= -delbp/2 & dyndat$bp <= delbp/2)
    indR = which(dyndat$bp > (half.win-2*delbp))
    dynnorm = dyndat$mean - mean(c(dyndat$mean[indL],dyndat$mean[indR]))
    uncnorm = uncdat$mean - mean(c(uncdat$mean[indL],uncdat$mean[indR]))
    x = seq(0,step*length(dynnorm)-step,by=step)
    integdyn = trapz(x,dynnorm)
    integunc = trapz(x,uncnorm)
    
    # compute peak index and pk density difference
    pk.dyn = mean(dyndat$mean[indM]) / mean(c(dyndat$mean[indL],dyndat$mean[indR]))
    pk.unc = mean(uncdat$mean[indM]) / mean(c(uncdat$mean[indL],uncdat$mean[indR]))
    if(comp.mode=="dyn"){
      pk.index = pk.dyn
      integ = integdyn
    }
    if(comp.mode=="unc"){
      pk.index = pk.unc
      integ = integunc
    }
    pk.dyn = format(pk.dyn,digits=3)
    pk.unc = format(pk.unc,digits=3)
    pk.index = format(pk.index,digits=3)
    
    # plot the peak traces
    fnt=1
    motif = dat$TF
    all = c(dyndat$mean, uncdat$mean)
    dd = max(all) - min(all)
    ylim = c( min(all)-0.1*dd, max(all)+0.1*dd )
    ylab = paste0(motif, " density")
    tt = comp
    main = paste0("Profile: ",tt[1], ", pk.index = ",pk.index, "\n percent diff: ",
                  per.dat1,"%, ", "% \n var: ",var.dat1,", ",
                  "\n npeaks: ",dat$npeakDYN,", ",dat$npeakUNC)
    plot(dyndat$bp, dyndat$mean, type="l", ylim=ylim, 
         ylab=ylab, xlab="distance from summit (bp)",
         lwd=2, col="red", main=main,cex.lab=fnt,cex.main=fnt)
    lines(uncdat$bp, uncdat$mean,lwd=2,col="darkgray")
    
    # plot smoothed curves
    main = paste0("peak indices = ",pk.dyn,", ",pk.unc,
                  "\n max % diff = ",format(per.diff,digits=3),
                  "\n integral = ", format(integ,digits=3))
    plot(dyndat$bp, dyndat$mean, type="l", ylim=ylim, 
         ylab="lowess smoothed", xlab="distance from summit (bp)",
         lwd=2, col="white", main=main,cex.lab=fnt,cex.main=fnt)
    fit.dyn = loess(mean~bp, dyndat, span=loess.span)
    fit.unc = loess(mean~bp, uncdat, span=loess.span)
    xaxis = seq(fit.dyn$x[1], fit.dyn$x[length(fit.dyn$x)], length.out=200)
    plt.dyn = predict(fit.dyn, xaxis)
    plt.unc = predict(fit.unc, xaxis)
    lines(xaxis, plt.dyn, lwd=2, col="red")
    lines(xaxis, plt.unc, lwd=2, col="darkgray")
    
    # load bigwig factor mapping data into peak coordinates
    # number of TFBS in each peak region
    dyn = bed.region.bpQuery.bigWig(bw.fimo, dyn.bed)
    unc = bed.region.bpQuery.bigWig(bw.fimo, unc.bed)
    names(dyn) = rownames(dyn.bed)
    names(unc) = rownames(unc.bed)
    
    # match atac peaks with motif coordinates
    # get numbers of peaks with >0 motifs
    dyn.with = length(which(dyn > 0))
    dyn.without = length(which(dyn == 0))
    nodif.with = length(which(unc > 0))
    nodif.without = length(which(unc == 0))
    result = cbind(c(nodif.with, nodif.without),
                   c(dyn.with, dyn.without))
    colnames(result) <- c("Unc", "Dyn")
    rownames(result) <- c("with motif", "without motif")
    resultn = apply(result,2,function(x){100*x/(sum(x))})
    
    # chi square test result
    chi = chisq.test(result)
    
    # plot data
    pv = format(chi$p.value,digits=3)
    main = paste0(motif,"\n p = ",pv,
                  "\n countsUnDy: ",result[1,1],", ",result[1,2],
                  "\n %: ",resultn[1,1]%>%format(digits=3),", ",
                  resultn[1,2]%>%format(digits=3))
    barplot(resultn[1,], col = c("gray","red"), 
            main=main, cex.main=fnt, ylim=c(0,100),
            ylab="percent", xlab="peak class",
            cex.lab=fnt,cex.names=fnt)
    barplot(resultn[2,], col = c("white"), 
            axes=F, axisnames=F,
            offset = resultn[1,], add=T)
    
    print(ii / nrow(enrich.dat))
    
  } # for loop
  dev.off()
} # plot.enrich.dyn

# > 300 peaks with motifs for the class with the max diff percentage
filt.maxdiff.pks = function(dat=NULL, npk.thresh=NULL){
  
  out = dat
  
  # unchanged class
  if(dat$class[1] == "un"){
    
    diff.per = dat %>% select(INCminusUN, DECminusUN)
    diff.per = apply(diff.per,2,function(x){-1*x})
    ind.rem = c()
    cl.max.diff = apply(diff.per,1,function(x){
      mm = max(x)
      ind = which(x == max(x))
      if(ind==1){return("inc")}
      if(ind==2){return("dec")}
    })
    
    for(ii in 1:length(cl.max.diff)){
      if(cl.max.diff[ii] == "inc"){
        if(dat$nmotINC[ii] < npk.thresh){ind.rem = c(ind.rem,ii)}
      }
      if(cl.max.diff[ii] == "dec"){
        if(dat$nmotDEC[ii] < npk.thresh){ind.rem = c(ind.rem,ii)}
      }
    }
    if(length(ind.rem)>0){out = dat[-ind.rem,]}
    
  } # unchanged class
  
  # increased class
  if(dat$class[1] == "inc"){
    
    diff.per = dat %>% select(INCminusUN, INCminusDEC)
    ind.rem = c()
    cl.max.diff = apply(diff.per,1,function(x){
      mm = max(x)
      ind = which(x == max(x))
      if(ind==1){return("un")}
      if(ind==2){return("dec")}
    })
    
    for(ii in 1:length(cl.max.diff)){
      if(cl.max.diff[ii] == "un"){
        if(dat$nmotUN[ii] < npk.thresh){ind.rem = c(ind.rem,ii)}
      }
      if(cl.max.diff[ii] == "dec"){
        if(dat$nmotDEC[ii] < npk.thresh ){ind.rem = c(ind.rem,ii)}
      }
    }
    if(length(ind.rem)>0){out = dat[-ind.rem,]}
    
  } # increased class
  
  # decreased class
  if(dat$class[1] == "dec"){
    
    diff.per = dat %>% select(DECminusUN, DECminusINC)
    ind.rem = c()
    cl.max.diff = apply(diff.per,1,function(x){
      mm = max(x)
      ind = which(x == max(x))
      if(ind==1){return("un")}
      if(ind==2){return("inc")}
    })
    
    for(ii in 1:length(cl.max.diff)){
      if(cl.max.diff[ii] == "un"){
        if(dat$nmotUN[ii] < npk.thresh){ind.rem = c(ind.rem,ii)}
      }
      if(cl.max.diff[ii] == "inc"){
        if(dat$nmotINC[ii] < npk.thresh){ind.rem = c(ind.rem,ii)}
      }
    }
    if(length(ind.rem)>0){out = dat[-ind.rem,]}
    
  } # decreased class
  
  return(out)
} # filt.maxdiff.pks


############################################################
# analysis functions
############################################################

# revised bed window function
bed.window2 = function(summits=NULL, half.win=NULL, chr=NULL){
  out = cbind(chr, summits-half.win, summits+half.win)
  out = as.data.frame(out,stringsAsFactors=FALSE)
  names(out)[2:3] = c("start", "end")
  rownames(out) = 1:nrow(out)
  out[,2:3] = apply(out[,2:3],2,function(x){data.matrix(x) %>% as.numeric})
  out$start[out$start<0] = 0
  return(out)
}

# reset bed coords based on window width
# input: summits =  vector of form "chr:start-end,summit"
# input: half.win = half window for plotting composites
# output: bed format coordinates around the peaks
bed.window = function(summits=NULL, half.win=NULL){
  chr = sapply(summits,function(x){strsplit(x,":")[[1]][1]})
  smt = sapply(summits,function(x){strsplit(x,",")[[1]][2]}) %>% as.numeric
  out = cbind(chr, smt-half.win, smt+half.win)
  out = as.data.frame(out,stringsAsFactors=FALSE)
  names(out)[2:3] = c("start", "end")
  rownames(out) = 1:nrow(out)
  out[,2:3] = apply(out[,2:3],2,function(x){data.matrix(x) %>% as.numeric})
  return(out)
} # bed.window

# function to find peak info from rownames of the results frame
get.map.from.res = function(namen=NULL){
  chrs = sapply(namen,function(x){strsplit(x,":")[[1]][1]})
  ends = sapply(namen,function(x){strsplit(x,"-")[[1]][2]})
  strs = sapply(namen,function(x){y=strsplit(x,"-")[[1]][1]; 
  return(strsplit(y,":")[[1]][2]) })
  coords = cbind(chrs,strs,ends) %>% as.data.frame(stringsAsFactors=F)
  coords[,2:3] = apply(coords[,2:3],2,function(x){data.matrix(x) %>% as.numeric})
  return( coords )
} # get.map.from.res

# matching peak data to the motif
motif.map = function(bigwig=NULL,coords=NULL,step=NULL,bps=NULL){
  dat = bed.step.bpQuery.bigWig(bw=bigwig, bed=coords, gap.value=0, step=step, as.matrix=T)
  dat = data.frame(colMeans(dat), bps)
  names(dat) = c("mean","bp")
  return(dat)
} # motif.map

# function to get average normalized read profiles
# input: reads = read count data frame
# input: times = condition identifier and numeric time columns
#        dataframe with names "condition" and "time"
# output: frame with means with respect to time for plotting
organize.reads = function(reads=NULL, times=NULL){
  
  # organized raw reads
  inds = sapply(times$condition,function(x){grep(x,names(reads))}) %>% as.vector
  read.org = reads[,inds] 
  
  # average at every time point
  t.profiles = sapply(times$condition,function(x){ 
    apply(read.org[,grep(x,colnames(read.org))],1,mean) 
  })
  colnames(t.profiles) = times$time
  return(t.profiles)
} # organize.reads


# function to get chi-square analysis table
# inputs are counts of motifs in peaks (integer counts)
# input: up = counts of motifs in increased peaks
# input: up = counts of motifs in decreased peaks
# input: up = counts of motifs in unchanged peaks
# the names of up/dn/un are peak identifiers, 
# e.g., chr8:11381934-11382229
get.xsquare.table = function(up=NULL, dn=NULL, un=NULL){
  inc.with = length(which(up > 0))
  inc.without = length(which(up == 0))
  dec.with = length(which(dn > 0))
  dec.without = length(which(dn == 0))
  nodif.with = length(which(un > 0))
  nodif.without = length(which(un == 0))
  result = cbind(c(nodif.with, nodif.without),
                 c(inc.with, inc.without),
                 c(dec.with, dec.without))
  colnames(result) <- c("Unchanged", "Increased", "Decreased")
  rownames(result) <- c("with motif", "without motif")
  return(result)
}

############################################################
# data import functions
############################################################

# function to import proseq read data maped to specific coordinates
# input: bed = bed file for mapping reads
# input: bigWig.path = path to bigWig files
# input: file.prefix = file prefix
# input: file.suffix = file suffix
# input: min.length = min gene length
# input: half.win = half of window size for plotting around the summit
# input: step = increment size within the plotting window (bp)
# input: sizefac = deseq2 size factors from the mapping to peaks
# input: vector of condition identifiers for omission
# output: data frame with averaged reads mapped for each condition
get.counts.interval <- function(bed=NULL, bigWig.path=NULL, 
                                file.prefix=NULL, 
                                file.suffix=NULL,
                                min.length=NULL,
                                half.win=NULL,
                                step=NULL,
                                sizefac=NULL,
                                omit=NULL) {
  
  # bp range centered on zero
  bps = seq(-half.win, half.win, length.out=2*half.win/step)
  
  # filter bed based on min gene length
  len = bed$end - bed$start
  indrem = which(len < min.length)
  if(length(indrem)>0){bed = bed[-indrem,]}
  
  # output files
  vec.names = c()
  output = data.frame(matrix(ncol = 0, nrow = length(bps)))
  
  # loop through each condition and load reads
  for(bigWig in Sys.glob(file.path(bigWig.path, 
      paste0(file.prefix, paste0("*",file.suffix))))) {
    
    # timepoint annotation
    factor.name = strsplit(bigWig, "/")[[1]]
    factor.name = strsplit(factor.name[length(factor.name)],
                           file.suffix)[[1]][1]
    factor.name = strsplit(factor.name,file.prefix)[[1]][2]
    ind = sapply(omit,function(x){grep(x,factor.name)}) %>% unlist
    if(length(ind)>0){next}
    vec.names = c(vec.names, factor.name)
    
    # get the sizefactor
    sf = sizefac[grep(factor.name,names(sizefac))]
    
    # load the data and normalize
    loaded.bw = load.bigWig(bigWig)
    reads = bed.step.bpQuery.bigWig(bw=loaded.bw, bed=bed, gap.value=0, step=step, as.matrix=T)
    reads = (1/sf) * reads
    
    # generate output
    output = cbind(output, colMeans(reads))
  } # bigWig
  
  # annotation and output
  colnames(output) = vec.names
  rownames(output) = bps
  return(output)
} # get.counts.interval


# function to import proseq read data maped to specific coordinates
# input: bed = bed file for mapping reads
# input: bigWig.path = path to bigWig files
# input: file.prefix = file prefix
# input: file.suffix = file suffix
# input: min.length = min gene length
# input: sizefac = deseq2 size factors from the mapping to peaks
# input: vector of condition identifiers for omission
# output: data frame with averaged reads mapped for each condition
get.counts.interval2 <- function(bed=NULL, bigWig.path=NULL, 
                                file.prefix=NULL, 
                                file.suffix=NULL,
                                min.length=NULL,
                                sizefac=NULL,
                                omit=NULL) {
  
  # filter bed based on min gene length
  len = bed$end - bed$start
  indrem = which(len < min.length)
  if(length(indrem)>0){bed = bed[-indrem,]}
  
  # output files
  vec.names = c()
  output = c()
  
  # loop through each condition and load reads
  for(bigWig in Sys.glob(file.path(bigWig.path, 
                                   paste0(file.prefix, paste0("*",file.suffix))))) {
    
    # timepoint annotation
    factor.name = strsplit(bigWig, "/")[[1]]
    factor.name = strsplit(factor.name[length(factor.name)],
                           file.suffix)[[1]][1]
    factor.name = strsplit(factor.name,file.prefix)[[1]][2]
    ind = sapply(omit,function(x){grep(x,factor.name)}) %>% unlist
    if(length(ind)>0){next}
    vec.names = c(vec.names, factor.name)
    
    # get the sizefactor
    sf = sizefac[grep(factor.name,names(sizefac))]
    
    # load the data and normalize
    loaded.bw = load.bigWig(bigWig)
    reads = bed.region.bpQuery.bigWig(bw=loaded.bw, bed=bed)
    reads = (1/sf) * reads
    
    # generate output
    output = cbind(output, reads)
  } # bigWig
  
  # annotation and output
  colnames(output) = vec.names
  rownames(output) = paste0(bed[,1],":",bed[,2],"-",bed[,3])
  return(output)
} # get.counts.interval2


# function to import proseq read data maped to specific coordinates
# input: bed = bed file for mapping reads
# input: bigWig.path = path to bigWig files
# input: file.prefix = file prefix
# input: file.suffix = file suffix
# input: min.length = min gene length
# input: half.win = half of window size for plotting around the summit
# input: step = increment size within the plotting window (bp)
# input: sizefac = deseq2 size factors from the mapping to peaks
# input: vector of condition identifiers for omission
# output: data frame with averaged reads mapped for each condition
get.counts.interval.log <- function(bed=NULL, bigWig.path=NULL, 
                                file.prefix=NULL, 
                                file.suffix=NULL,
                                min.length=NULL,
                                half.win=NULL,
                                step=NULL,
                                sizefac=NULL,
                                omit=NULL) {
  
  # bp range centered on zero
  bps = seq(-half.win, half.win, length.out=2*half.win/step)
  
  # filter bed based on min gene length
  len = bed$end - bed$start
  indrem = which(len < min.length)
  if(length(indrem)>0){bed = bed[-indrem,]}
  
  # output files
  vec.names = c()
  output = data.frame(matrix(ncol = 0, nrow = length(bps)))
  
  # loop through each condition and load reads
  for(bigWig in Sys.glob(file.path(bigWig.path, 
                                   paste0(file.prefix, paste0("*",file.suffix))))) {
    
    # timepoint annotation
    factor.name = strsplit(bigWig, "/")[[1]]
    factor.name = strsplit(factor.name[length(factor.name)],
                           file.suffix)[[1]][1]
    factor.name = strsplit(factor.name,file.prefix)[[1]][2]
    ind = sapply(omit,function(x){grep(x,factor.name)}) %>% unlist
    if(length(ind)>0){next}
    vec.names = c(vec.names, factor.name)
    
    # get the sizefactor
    sf = sizefac[grep(factor.name,names(sizefac))]
    
    # load the data and normalize
    loaded.bw = load.bigWig(bigWig)
    reads = bed.step.bpQuery.bigWig(bw=loaded.bw, bed=bed, gap.value=0, step=step, as.matrix=T)
    reads = log2( (1/sf) * reads + 0.01)
    
    # generate output
    output = cbind(output, colMeans(reads))
  } # bigWig
  
  # annotation and output
  colnames(output) = vec.names
  rownames(output) = bps
  return(output)
} # get.counts.interval.log
