############################################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Figure 1
# acquire and process AAGMEx data
############################################################################

setwd("/media/wa3j/Seagate2/Documents/adipose_sex_ms/get_geo_data")

# load relevant libraries
library(Biobase)
library(GEOquery)
library(dplyr)

# load series and platform data from GEO
gset <- getGEO("GSE95674", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL10904", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# log2 transform if necessary
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) {
  print("log2 transform performed")
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex)
} else {print("log2 transform not performed")}
expr_data = exprs(gset)

# get gene list
gene_list = rownames(expr_data)

# sample annotation
pdat = pData(gset) %>% select(geo_accession,source_name_ch1,organism_ch1,
                              characteristics_ch1,characteristics_ch1.1,characteristics_ch1.2)
names(pdat) = c("accession","pheno","organism","sex","age","tissue")
pdat[] = apply(pdat,2,as.character)

# process sample data
diabetic_status = sapply(pdat$pheno,function(x){strsplit(x," ")[[1]][1]})
ethnicity = sapply(pdat$pheno,function(x){
  split1 = strsplit(x," ") %>% unlist
  out = paste0(split1[2],"_",split1[3])
  return(out)})
age = sapply(pdat$age,function(x){strsplit(x," ")[[1]][5] %>% as.numeric})
sex = sapply(pdat$sex,function(x){strsplit(x," ")[[1]][2]})
unique(diabetic_status)
unique(sex)
unique(pdat$tissue)
all.equal(nrow(pdat), length(diabetic_status), length(ethnicity),
          length(age), length(sex))
phenotypes0 = pdat %>% select(accession,tissue) %>%
  mutate(diabetic_status=diabetic_status,age=age,sex=sex)

# check data organization
# note that gene names are the row identifiers for this data set
all(phenotypes0$accession == names(expr_data))

# write data
write.table(expr_data,"expr_data_aagmex.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(phenotypes0,"phenotypes_aagmex.txt",col.names=T,row.names=F,
            quote=F,sep="\t")
