

library(dplyr)

############################################################
# integrate enrichment data from 3 peak classes
############################################################

dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/atac_time_deg/motifEnrich/combined")
setwd(dir)

############################################################
# load data
############################################################

######################################
# summary tables for motif enrichment

fname = "enrichSummary750k_deg.txt"
enrich.summary.deg = read.table(fname,header=T,sep="\t",stringsAsFactors=F)

fname = "enrichSummary750k_dreg.txt"
enrich.summary.dreg = read.table(fname,header=T,sep="\t",stringsAsFactors=F)

fname = "enrichSummary750k_dyn.txt"
enrich.summary.dyn = read.table(fname,header=T,sep="\t",stringsAsFactors=F)


######################################
# summary tables for combined de novo and motif enrichment

fname = "summary.atacDEG.txt"
comb.summary.deg = read.table(fname,header=T,sep="\t",stringsAsFactors=F)

fname = "summary.atacDREG.txt"
comb.summary.dreg = read.table(fname,header=T,sep="\t",stringsAsFactors=F)

fname = "summary.atacDyn.txt"
comb.summary.dyn = read.table(fname,header=T,sep="\t",stringsAsFactors=F)


######################################
# TF cluster annotations

fname = "TFclusterAnnotation750k_deg.txt"
ann.clust.deg = read.table(fname,header=T,sep="\t",stringsAsFactors=F)

fname = "TFclusterAnnotation750k_dreg.txt"
ann.clust.dreg = read.table(fname,header=T,sep="\t",stringsAsFactors=F)

fname = "TFclusterAnnotation750k_dyn.txt"
ann.clust.dyn = read.table(fname,header=T,sep="\t",stringsAsFactors=F)

############################################################
# integrate information from the three analysis
############################################################

# TFids for each analysis
degTF = enrich.summary.deg$factor
dynTF = enrich.summary.dyn$factor
dregTF = enrich.summary.dreg$factor

# TF motifs identified by all three analyses
TF.all3 = Reduce(intersect, list(degTF,dynTF,dregTF))

# TF motifs identified by two analyses
TF.deg.dreg = intersect(degTF, dregTF)
TF.dyn.dreg = intersect(dynTF, dregTF)
TF.deg.dyn = intersect(degTF, dynTF)
TF.deg.dreg = TF.deg.dreg[!(TF.deg.dreg %in% TF.all3)]
TF.dyn.dreg = TF.dyn.dreg[!(TF.dyn.dreg %in% TF.all3)]
TF.deg.dyn = TF.deg.dyn[!(TF.deg.dyn %in% TF.all3)]

# TF motifs unique to a given analysis
TF.mult = unique(c(TF.all3,TF.deg.dreg,TF.dyn.dreg,TF.deg.dyn))
TF.deg1 = degTF[!(degTF %in% TF.mult)]
TF.dyn1 = dynTF[!(dynTF %in% TF.mult)]
TF.dreg1 = dregTF[!(dregTF %in% TF.mult)]

# checks
length(degTF) == length(TF.all3) + length(unique(c(TF.deg.dreg,TF.deg.dyn))) + length(TF.deg1)
length(dynTF) == length(TF.all3) + length(unique(c(TF.dyn.dreg,TF.deg.dyn))) + length(TF.dyn1)
length(dregTF) == length(TF.all3) + length(unique(c(TF.deg.dreg,TF.dyn.dreg))) + length(TF.dreg1)
length(unique(c(degTF, dynTF, dregTF))) == length(TF.mult) + length(TF.deg1) + length(TF.dyn1) + length(TF.dreg1)

############################################################
# get community annotations for the analyses 
############################################################

