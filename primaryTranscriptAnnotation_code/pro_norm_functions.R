
# function to import proseq read data maped to specific coordinates
# input: bed = bed6 file for mapping reads
# input: bigWig.path = path to bigWig files
# input: file.prefix = file prefix
# input: file.suffix = file suffix
# input: min.length = min gene length
# output: data frame with reads mapped for each condition
get.counts.interval <- function(bed=NULL, bigWig.path=NULL, 
                                file.prefix=NULL, file.suffix=NULL,
                                min.length=NULL) {
  
  # filter bed based on min gene length
  len = bed$end - bed$start
  ind = which(len < min.length)
  if(length(ind)>0){
    bed = bed[-ind,] 
  }
  
  # output files
  vec.names = c()
  output = data.frame(matrix(ncol = 0, nrow = nrow(bed)))
  
  # loop through each condition and load reads
  for(bigWig in Sys.glob(file.path(bigWig.path, 
      paste0(file.prefix, paste0("*",file.suffix))))) {
    factor.name = strsplit(bigWig, "/")[[1]]
    factor.name = strsplit(factor.name[length(factor.name)],file.suffix)[[1]][1]
    factor.name = strsplit(factor.name,file.prefix)[[1]][2]
    print(factor.name)
    vec.names = c(vec.names, factor.name)
    loaded.bw.plus = load.bigWig(bigWig)
    bg.minus = paste0(bigWig.path,'/',file.prefix,factor.name,
                      gsub("plus","minus",file.suffix))
    loaded.bw.minus = load.bigWig(bg.minus)
    reads = bed6.region.bpQuery.bigWig(loaded.bw.plus, loaded.bw.minus, bed)
    output = cbind(output, reads)
  } # bigWig
  
  # annotation and output
  colnames(output) = vec.names
  r.names = paste0(bed[,1],':',bed[,2],'-',bed[,3],'_',bed[,4])
  row.names(output) = r.names
  return(output)
} # get.counts.interval




