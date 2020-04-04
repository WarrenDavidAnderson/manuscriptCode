

cd /media/wa3j/Seagate2/Documents/adipose_sex_ms/final/BART

##########################################
# download data table and parse for BART analyses
##########################################

library(dplyr)

# get data
fname = "supplementaryTable1.txt"
dat0 = read.table(fname,header=T,stringsAsFactors=F)

# data for all human genes
out = dat0$gene
fname = "allHuman.txt"
write.table(out,fname,col.names=F,row.names=F,quote=F)

# data for F>M human genes
out = dat0 %>% filter(gtex_fc > 0)
out = out$gene
fname = "femaleHuman.txt"
write.table(out,fname,col.names=F,row.names=F,quote=F)

# data for M>F human genes
out = dat0 %>% filter(gtex_fc < 0)
out = out$gene
fname = "maleHuman.txt"
write.table(out,fname,col.names=F,row.names=F,quote=F)

# data for mouse ids
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
out = dat0 %>% filter(is.na(hmdp_fc)==F)
out = out[which(abs(out$hmdp_fc)>1.05),]
out = out %>% filter(hmdp_fdr<0.05)
ss = sign(out$gtex_fc*out$hmdp_fc)
out = out$gene[which(ss==1)]
out = tolower(out)
out = firstup(out)
fname = "mouse.txt"
write.table(out,fname,col.names=F,row.names=F,quote=F)
