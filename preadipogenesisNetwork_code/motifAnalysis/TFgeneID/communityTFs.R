

dir = paste0("/media/wa3j/Seagate2/Documents/PRO/",
             "adipogenesis/July2018/atac_time_deg/TFlist")
setwd(dir)

library(dplyr)

############################################################
# import data
############################################################

# import genencode gene annotations
fname = "gencode.vM18.annotation.bed"
dat0 = read.table(fname,header=F,stringsAsFactors=F)
names(dat0) = c('chr', 'start', 'end', 'gene', 'xy', 'strand')
dat0 = unique(dat0)
gencode = unique(dat0$gene)

# import manually aggregated gene lists based on TFclass
tfclass.lists0 = list()
all.files = list.files()
tf.files = all.files[grep("_genes.txt",all.files)]
for(ii in tf.files){
  comm = strsplit(ii,"_genes")[[1]][1]
  comm = strsplit(comm,"comm")[[1]][2]
  comm = paste0("community_",comm)
  tfclass.lists0[[comm]] = read.table(ii,stringsAsFactors=F)[,1]
}

############################################################
# document aggregated communities and TFclass information
############################################################

community_42.comms = c("community_42")
community_42.annos = c("1.1 Basic leucine zipper factors (bZIP)","1.1.1 Jun-related")

community_67.comms = c("community_67","community_69","community_157","community_268")
community_67.annos = c("1.1 Basic leucine zipper factors (bZIP)","1.1.2 Fos-related")

community_58.comms = c("community_58","community_90","community_96","community_159","community_166")
community_58.annos = c("1.1 Basic leucine zipper factors (bZIP)","1.1.3 Maf-related")

community_23.comms = c("community_23")
community_23.annos = c("1.1 Basic leucine zipper factors (bZIP)","1.1.8 CEBP-related") 

community_47.comms = c("community_47","community_62","community_63","community_64","community_66",
                       "community_76","community_88","community_140","community_364")
community_47.annos = c("1.2 Basic helix-loop-helix factors (bHLH)","1.2.1 E2A-related factors",
                       "1.2.2 MyoD/ASC-related factors","1.2.3 Tal-related")

community_101.comms = c("community_101")
community_101.annos = c("1.2 Basic helix-loop-helix factors (bHLH)","1.2.4 Hairy-related")

community_158.comms = c("community_158")
community_158.annos = c("1.2 Basic helix-loop-helix factors (bHLH)","1.2.6 bHLH-ZIP","1.2.6.1 TFE3")

community_44.comms = c("community_44")
community_44.annos = c("1.2 Basic helix-loop-helix factors (bHLH)","1.2.6 bHLH-ZIP","1.2.6.5 MYC")

community_4.comms = c("community_4")
community_4.annos = c("1.3 Basic helix-span-helix factors (bHSH)","1.3.1 AP2")

community_77.comms = c("community_77")
community_77.annos = c("2.1 Nuclear receptors with C4 zinc fingers","2.1.1 Steroid hormone receptors",
                       "2.1.1.1 GR-like(NR3C)")

community_356.comms = c("community_356")
community_356.annos = c("2.1 Nuclear receptors with C4 zinc fingers",
                        "2.1.2 Thyroid hormone receptor-related factors","2.1.2.1 RAR(NR1B)")

community_338.comms = c("community_338")
community_338.annos = c("2.1 Nuclear receptors with C4 zinc fingers",
                        "2.1.2 Thyroid hormone receptor-related factors","2.1.2.2 T3R(NR1A)")

community_46.comms = c("community_46")
community_46.annos = c("2.1 Nuclear receptors with C4 zinc fingers",
                       "2.1.2 Thyroid hormone receptor-related factors","2.1.2.5 PPAR(NR1C)")

community_384.comms = c("community_384")
community_384.annos = c("2.1 Nuclear receptors with C4 zinc fingers",
                        "2.1.3 RXR-related receptors","2.1.3.3 TLX(NR2E1)")

community_16.comms = c("community_16")
community_16.annos = c("2.1 Nuclear receptors with C4 zinc fingers","2.1.4 NGFI(NR4A)")

community_79.comms = c("community_79")
community_79.annos = c("2.3 C2H2 zinc finger factors",
                       "2.3.1 Three-zinc finger Krüppel-related","2.3.1.2 Kr-like")

community_244.comms = c("community_244")
community_244.annos = c("2.3 C2H2 zinc finger factors","2.3.3 More than 3 adjacent zinc fingers","2.3.3.0.39")

community_382.comms = c("community_382")
community_382.annos = c("2.3 C2H2 zinc finger factors",
                        "2.3.3 More than 3 adjacent zinc fingers","2.3.3.0.127 ZNF416")

community_134.comms = c("community_134")
community_134.annos = c("2.3 C2H2 zinc finger factors",
                        "2.3.3 More than 3 adjacent zinc fingers","2.3.3.0.197 ZNF189")

community_156.comms = c("community_156")
community_156.annos = c("2.3 C2H2 zinc finger factors",
                        "2.3.3 More than 3 adjacent zinc fingers","2.3.3.11 ZBTB6-like")

community_126.comms = c("community_126")
community_126.annos = c("2.3 C2H2 zinc finger factors",
                        "2.3.3 More than 3 adjacent zinc fingers",
                        "2.3.3.13 ZNF148-like")

community_75.comms = c("community_75")
community_75.annos = c("2.3 C2H2 zinc finger factors",
                       "2.3.3 More than 3 adjacent zinc fingers","2.3.3.16 ZNF238-like")

community_6.comms = c("community_6")
community_6.annos = c("2.3 C2H2 zinc finger factors",
                      "2.3.3 More than 3 adjacent zinc fingers","2.3.3.22 BCL6")

community_20.comms = c("community_20")
community_20.annos = c("2.3 C2H2 zinc finger factors",
                       "2.3.3 More than 3 adjacent zinc fingers","2.3.3.50 CTCF-like")

community_259.comms = c("community_259","community_121")
community_259.annos = c("2.3 C2H2 zinc finger factors",
                        "2.3.3 More than 3 adjacent zinc fingers","2.3.3.65 ZFX-ZFY")

community_141.comms = c("community_141","community_178")
community_141.annos = c("2.3 C2H2 zinc finger factors",
                        "2.3.3 More than 3 adjacent zinc fingers","2.3.3.80 ZNF99-like")

community_302.comms = c("community_302")
community_302.annos = c("2.3 C2H2 zinc finger factors",
                        "2.3.4 Factors with multiple dispersed zinc fingers","2.3.4.0.68 ZNF519")

community_311.comms = c("community_311")
community_311.annos = c("2.3 C2H2 zinc finger factors",
                        "2.3.4 Factors with multiple dispersed zinc fingers","2.3.4.17 HIC")

community_61.comms = c("community_61","community_83")
community_61.annos = c("2.3 C2H2 zinc finger factors",
                       "2.3.1 Three-zinc finger Krüppel-related factors","2.3.1.1 Sp1-like",
                       "2.3.1.3 EGR","2.3.3 More than 3 adjacent zinc finger factors",
                       "2.3.4 Factors with multiple dispersed zink fingers")

community_376.comms = c("community_376")
community_376.annos = c("3.1 Homeo domain factors","3.1.4 TALE-type HD","3.1.4.2 MEIS")

community_36.comms = c("community_36")
community_36.annos = c("3.1 Homeo domain factors","3.1.4.4 PBX")

community_265.comms = c("community_265")
community_265.annos = c("3.2 Paired box factors","3.2.2 PD","3.2.2.2 PAX2-like")

community_45.comms = c("community_45")
community_45.annos = c("3.3 Fork head / winged helix factors","3.3.1 FOX")

community_13.comms = c("community_13","community_235") 
community_13.annos = c("3.3 Fork head / winged helix factors","3.3.2 E2F")

community_74.comms = c("community_74","community_131","community_153","community_217")
community_74.annos = c("3.5 Tryptophan cluster factors","3.5.2 Ets-related")

community_2.comms = c("community_2")
community_2.annos = c("3.6 TEA domain factors","3.6.1 TEF1-related")

community_39.comms = c("community_39")
community_39.annos = c("4.2 Heteromeric CCAAT-binding factors",
                       "4.2.1 Heteromeric CCAAT-binding factors")

community_114.comms = c("community_114")
community_114.annos = c("5.1 MADS box factors","5.1.2 Responders to external signals")

community_11.comms = c("community_11","community_143")
community_11.annos = c("6.2 STAT domain factors","6.2.1 STAT")

community_12.comms = c("community_12")
community_12.annos = c("6.4 Runt domain factors","6.4.1 Runt-related")

community_41.comms = c("community_41")
community_41.annos = c("6.4 Runt domain factors","6.4.1 Runt-related")

community_146.comms = c("community_146")
community_146.annos = c("6.7 Grainyhead","6.7.2 CP2-related")

community_1.comms = c("community_1","community_273","community_369")
community_1.annos = c("7.1 SMAD/NF-1 DNA-binding domain factors","7.1.1 SMAD")

community_5.comms = c("community_5","community_94","community_122","community_137","community_378")
community_5.annos = c("2.3 C2H2 zinc finger factors","2.3.3 More than 3 adjacent zinc fingers",
                      "2.3.3.52 ZNF322-like","3.1 Homeo domain factors","3.1.2.21 TLX",
                      "7.1 SMAD/NF-1 DNA-binding domain factors","7.1.2 NF-1")

community_25.comms = c("community_25")
community_25.annos = c("2.3 C2H2 zinc finger factors","2.3.4 Factors with multiple dispersed zinc fingers",
                       "2.3.4.14 EVI1-like","6.1 Rel homology region (RHR) factors",
                       "6.1.5 Early B-Cell Factor-related factors")

community_34.comms = c("community_34")
community_34.annos = c("0 Yet undefined DNA-binding domains")

community_136.comms = c("community_136")
community_136.annos = c("no DNA binding domain found")

##### dimer
community_267.comms = c("community_267")
community_267.annos = c("3.2 Paired box factors","3.2.1 PD+HD",
                        "3.3 Fork head / winged helix factors",
                        "3.3.1 FOX","3.3.1.15 FOXO")

############################################################
# get gene lists for each PWM group (TF class)
############################################################

# function to make first letter uppercase
# https://stackoverflow.com/questions/6364783/capitalize-the-first-letter-of-both-words-in-a-two-word-string
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

# function to make letters lower case
make.lower = function(namen=NULL){
  out = c()
  for(ii in 1:length(namen)){
    raw = namen[ii]
    lowr = tolower(raw)
    out = c(out, simpleCap(lowr))
  }
  return(out)
}

# function to include genes that match gencode and identify
# classes with mismatches
get.tf.list = function(tf.class=NULL, gencode=NULL){
  
  # new gene name list and document mismatches
  tfclass.lists = list()
  mismatches = c()
  
  # loop through each set of genes
  for(ii in 1:length(tf.class)){
    comm = names(tf.class)[ii]
    gen0 = tf.class[[ii]]
    gen1 = make.lower(gen0)
    gncd = gen1[gen1 %in% gencode]
    if( length(gen1) == length(gncd) ){
      tfclass.lists[[comm]] = gncd
    } else {mismatches = c(mismatches, comm)}
  } #for loop
  out = list(newlist=tfclass.lists, mismatches=mismatches)
} # get.tf.list

# get new lists and mismatches
tfdat = get.tf.list(tf.class=tfclass.lists0, gencode=gencode)
tfclass.lists = tfdat$newlist 
mismatches = tfdat$mismatches

########## manually address mismatches ########## 
# https://www.ncbi.nlm.nih.gov/gene
# https://www.genecards.org --> for identifying orthologs
# e.g., ZNF vs Zfp genes
#################################################

# community_101
ii=1
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Bhlhe40","Bhlhe41")
new = c(gncd, new)
length(gen1) == length(new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_126
ii=2
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Zfp148","Zfp281")
new = c(gncd, new)
length(gen1) == length(new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_13
ii=3
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Tfdp1")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_134
ii=4
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Zfp189")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_141
ii=5
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Zbtb14")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_146
ii=6
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Tfcp2")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_16
ii=7
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Nr4a1","Nr4a3","Nr4a2")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_25
ii=8
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Mecom")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_302
ii=9
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Zfp386","Zfp493")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_338
ii=10
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Thra","Thrb")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_382
ii=11
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Zfp418")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new


# community_384
ii=12
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Nr2e3","Nr2e1")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_4
ii=13
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Tfap2c","Tfap2d","Tfap2a","Tfap2b","Tfap2e")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_42
ii=14
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Jun")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_44
ii=15
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Mycl","Myc","Mycn")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_46
ii=16
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Ppara","Ppard","Pparg")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_47
ii=17
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Tcf3","Figla","Nhlh1","Nhlh2","Tcf12","Myod1","Neurog1","Neurog2","Neurog3","Tcf4")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_5
ii=18
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Nfia","Nfib","Nfic","Nfix","Zfp322a")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_58
ii=19
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Maf")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_61
ii=20
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Zbtb14","Zfp131","Zkscan7","Zfp184","Zscan26","Zfp189","Zkscan6","Zkscan8","Zfp202","Zfp13",
        "Zfp422","Zfp239","Zfp647","Zfp263","Zfp266","Zfp846","Zfp110","Zfp369","Zfp275","Zfp287","Zfp317","Zfp329",
        "Zfp354c","Zfp358","Zfp105","Zkscan14","Zfp408","Zfp444","Zfp449","Zfp454","Zkscan16","Zkscan17",
        "Zscan25","Zfp74","Zfp575","Zfp623","Gm3055","Zfp641","Zfp648","Zfp655","Zfp664",
        "Zfp667","Zfp708","Zfp691","Zfp707","Zfp746","Zfp768","Gm3055","Zfp790","Zfp940","Zfp791","Zfp7",
        "Zfp879","Zfp148","Zfp281","Zbtb7b","Zfp142","Zfp236","Zfp251","Zfp277","Zbtb21","Zfp316",
        "Zfp335","Zfp341","Zfp407","Zfp438","Zfp445","Zfp451","Zfp462","Zfp467","Zfp469","Zfp553",
        "Zfp507","Zfp513","Zfp579","Zfp597","Zfp606","Zfp618","Zfp628","Zfp629","Zfp644",
        "Zfp646","Zfp668","Zfp697","Zfp770","Zfp775","Zfp784","Zfp787","Zfp800","Vezf1")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new[new%in%gencode]

# community_67
ii=21
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Fos","Fosl1","Fosl2")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_74
ii=22
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Ets1")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_75
ii=23
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Zbtb18")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

# community_77
ii=24
comm = mismatches[ii]
gen0 = tfclass.lists0[[comm]] %>% unique
gen1 = make.lower(gen0)
gncd = gen1[gen1 %in% gencode]
not = gen1[!(gen1 %in% gencode)]
length(gen1) == length(gncd)
new = c("Nr3c1","Nr3c2","Pgr")
new = c(gncd, new)
length(new) == length(new[new%in%gencode])
tfclass.lists[[comm]] = new

#### run a check with the original function
tfdat = get.tf.list(tf.class=tfclass.lists, gencode=gencode)
length(tfdat$newlist)


############################################################
# generate annotation list for community groupings
############################################################

tf.community.ann = list()
for(ii in names(tfclass.lists)){
  comms = get( paste0(ii,".comms") )
  anns = get( paste0(ii,".annos") )
  new = list(comms = comms, anns = anns)
  tf.community.ann[[ii]] = new
}

# check counts
ncnt = lapply(tf.community.ann,function(x){length(x$comms)}) %>% unlist
length(ncnt)
sum(ncnt)

############################################################
# save key information:
# get.tf.list - community gene lists
# tf.community.ann - community annotation
############################################################

save(tfclass.lists, file="TFcommunity.gene.lists.RData")
save(tf.community.ann, file="TFcommunity.annotations.RData")





