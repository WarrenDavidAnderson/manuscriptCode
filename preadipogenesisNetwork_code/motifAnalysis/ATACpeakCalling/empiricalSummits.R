

library(dplyr)
library(DESeq2)
library(ggplot2)
library(bigWig)
library(prodlim)

setwd("/run/user/1001/gvfs/smb-share:server=home1.virginia.edu,share=cphgdesk/users/wa3j/MGlab/Analysis/Adipogenesis/ATAC_time_DEG")
source("atac_norm_functions.R")

dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/atac_time_deg")
setwd(dir)


############################################################
## key analysis parameters
############################################################

window = 50 # window = sliding window size (bp)
bin = 5 # bin = increment by which the sliding window will move (bp)

############################################################
## import macs2 peak data and ATAC bigWigs
############################################################

# peak coordinates
# https://github.com/taoliu/MACS
fname = "3T3_atac_4hpeaks_sorted.bed"
bed0 = read.table(fname,stringsAsFactors=F,header=F)
names(bed0) = c("chr","start","end","rel","summit")
all.peaks = paste0(bed0[,1],":",bed0[,2],"-",bed0[,3])

# aggregated/integrated preadipogenic ATAC bigWig
bw.file = paste0("/media/wa3j/Seagate2/Documents/PRO/",
                 "adipogenesis/July2018/integratedBW/preadip.bigWig")
loaded.bw = load.bigWig(bw.file)


############################################################
# function to get coordinates for all peaks
############################################################

# function to find peak info from rownames of the results frame
get.map.from.res = function(res){
  namen = res
  chrs = sapply(namen,function(x){strsplit(x,":")[[1]][1]})
  ends = sapply(namen,function(x){strsplit(x,"-")[[1]][2]})
  strs = sapply(namen,function(x){y=strsplit(x,"-")[[1]][1]; 
  return(strsplit(y,":")[[1]][2]) })
  coords = cbind(chrs,strs,ends) %>% as.data.frame(stringsAsFactors=F)
  return( coords )
} # get.map.from.res

coords0 = get.map.from.res(all.peaks)
coords0[,2:3] = apply(coords0[,2:3],2,function(x){data.matrix(x) %>% as.numeric})

############################################################
#  functions to loop through peaks and empirically find summits
############################################################

# function to get coordinates for mapping reads to identify a peak window
# input: coords = coordinates in a single line with bed format
# input: win = window size (bp)
# input: del = delta for sliding the window over the coordinates (bp)
# output: mat = dataframe with bed coordinates
get.coords = function(coords=NULL, win=NULL, del=NULL){
  chr = coords[1] %>% as.character
  range = (coords[3] - coords[2]) %>% data.matrix %>% as.numeric
  n.win = ceiling( (range-win)/del )
  len = n.win * del + win
  diff = len - range
  addL = floor(diff/2)
  addR = ceiling(diff/2)
  str = (coords[2] - addL) %>% data.matrix %>% as.numeric
  end = (coords[3] + addR) %>% data.matrix %>% as.numeric
  mat = c()
  for(ii in 1:(n.win+1)){
    st = str + (ii-1) * del
    ed = st + win
    mat = rbind(mat, c(chr, st, ed))
  }
  mat = as.data.frame(mat,stringsAsFactors=FALSE)
  names(mat) = c("chr","start","end")
  mat[,2:3] = apply(mat[,2:3],2,function(x){data.matrix(x)%>%as.numeric})
  return(mat)
} # get.coords

# function to empirically find summits for a set of peaks
# we slide a window along the peak and identify the center
# of the window with the highest reads - that is the summit
# input: peaks = bed file with ATAC peaks
# input: window = sliding window size (bp)
# input: bin = increment by which the sliding window will move (bp)
# input: bw.file = bigWig file with path included
# output = bed with peaks and summits
find.summits = function(peaks=NULL, window=NULL, bin=NULL, bw.file=NULL){
  
  # loop through all peaks
  summits = rep(0,nrow(peaks))
  for(ii in 1:nrow(peaks)){
    
    # set the summit to the middle for peak < window
    dpeak = peaks[ii,3] - peaks[ii,2]
    if(dpeak < window){
      summits[ii] = floor(peaks[ii,2] + dpeak/2)
      next
    }
    
    # identify the summit as the center of the max peak
    windows = get.coords(coords=peaks[ii,], win=window, del=bin)
    counts = bed.region.bpQuery.bigWig(bw.file, windows)
    ind = which(counts == max(counts))[1]
    summits[ii] = floor(windows[ind,2] + window/2)
    
    # progress tracking
    if((100*ii/nrow(peaks)) %% 1 == 0){
      Sys.sleep(0.5)
      print(paste0(100*ii/nrow(peaks),"% completed"))
      flush.console()
    }
    
  } ## ii, peak loop
  
  # return bed with peaks and summits
  return(cbind(peaks,summits))
  
} # find.summits

############################################################
#  get all summits empirically
############################################################

summits.bed = find.summits(peaks=bed0[,1:3], window=window, bin=bin, bw.file=loaded.bw)
save(summits.bed, file="ATACsummits_20190914.RData")

