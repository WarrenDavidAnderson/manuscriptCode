
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
                                omit=NULL,
                                mode=NULL) {
  
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
    
    # load the data and normalize (or not)
    loaded.bw = load.bigWig(bigWig)
    reads = bed.step.bpQuery.bigWig(bw=loaded.bw, bed=bed, gap.value=0, step=step, as.matrix=T)
    if(mode == "log2"){
      reads = log2( (1/sf) * reads + 0.01)
    }
    if(mode == "raw"){
      reads = (1/sf) * reads 
    }
    
    # generate output
    output = cbind(output, colMeans(reads))
  } # bigWig
  
  # annotation and output
  colnames(output) = vec.names
  rownames(output) = bps
  return(output)
} # get.counts.interval.log


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
counts.interval.heatmap <- function(bed=NULL, bigWig.path=NULL, 
                                    file.prefix=NULL, 
                                    file.suffix=NULL,
                                    min.length=NULL,
                                    half.win=NULL,
                                    step=NULL,
                                    sizefac=NULL,
                                    omit=NULL,
                                    mode=NULL) {
  

  # filter bed based on min gene length
  len = bed$end - bed$start
  indrem = which(len < min.length)
  if(length(indrem)>0){bed = bed[-indrem,]}
  
  # output files
  vec.names = c()
  output = data.frame(matrix(ncol = 0, nrow = nrow(bed)))
  
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
    
    # load the data and normalize (or not)
    loaded.bw = load.bigWig(bigWig)
    reads = bed.region.bpQuery.bigWig(bw=loaded.bw, bed=bed)
    if(mode == "log2"){
      reads = log2( (1/sf) * reads + 0.01)
    }
    if(mode == "raw"){
      reads = (1/sf) * reads 
    }
    
    # generate output
    output = cbind(output, reads)
  } # bigWig
  
  # annotation and output
  colnames(output) = vec.names
  rownames(output) = paste0(bed[,1],":",bed[,2],"-",bed[,3])
  return(output)
} # get.counts.interval.log
