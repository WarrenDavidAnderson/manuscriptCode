############################################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Figure S5
############################################################################

# data directory
dir="/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/fig6"
setwd(dir)

library(reshape)
library(dplyr)
library(NMF)


###########################################################
## read and format data
############################################################

fname = "HMDP_Adipose_Gene_Expression.csv"
expr0 = read.table(fname,sep=",",header=T,stringsAsFactors=F)

fname = "HMDP_Traits.csv"
pheno0 = read.table(fname,sep=",",header=T,stringsAsFactors=F)

# initial annotation for the phenotype set
annotation = pheno0 %>% select(mouse_number,Strain,Sex)

# function for finding a given annotation associated with a given vector
var.annot = function(vecdat=NULL, annot=NULL, invar=NULL, outvar=NULL){
  ann.in = ( annot %>% select(invar) )[,1]
  ann.out = ( annot %>% select(outvar) )[,1]
  out = sapply(vecdat,function(x){ann.out[which(ann.in==x)][1]}) %>% unlist %>% unname
  return(out)
}

# set exression matrix and annotation
# main objects: 
##### expr2: gene expression matrix
##### expr2.samp.annotation: column/sample annotation for expr2; mouse, strain, sex
##### expr2.gene.annotation: row/probe annotation for expr2; probe, gene name
##### note reshape warning: "multiple rows match for mouse_number=339: first taken"

## reorganizes data so that there is a column for expression_value.mouse_number for each mouse
expr1 = reshape(expr0[,c(1,2,5)], idvar="probesetID", timevar="mouse_number", direction="wide")

## applies function to column names that gets rid of expression_value. and leaves just the number
expr.mouse.nums = sapply(names(expr1)[2:ncol(expr1)],function(x){strsplit(x,"expression_value.")[[1]][2]})

## get annotation for strain values
expr.mouse.strains = var.annot(vecdat=expr.mouse.nums, annot=annotation, invar="mouse_number", outvar="Strain")
expr.mouse.sex = var.annot(vecdat=expr.mouse.nums, annot=annotation, invar="mouse_number", outvar="Sex")

## just expression values
expr2 = expr1[,2:ncol(expr1)]

# adds row names that are the probesetIDs
rownames(expr2) = expr1$probesetID

#rename columns with just the mouse number
names(expr2) = expr.mouse.nums

#group together the annotations created with the function previously
#there is one annotation that is mouse number, strain, and sex
expr2.samp.annotation = cbind(expr.mouse.nums, expr.mouse.strains, expr.mouse.sex) %>% as.data.frame(stringsAsFactors=F)
names(expr2.samp.annotation) = c("mouse_number","strain","sex")

#there is another annotation that is probesetID with gene symbol
expr.genes = var.annot(vecdat=rownames(expr2), annot=expr0, invar="probesetID", outvar="gene_symbol")
expr2.gene.annotation = cbind(rownames(expr2), expr.genes) %>% as.data.frame(stringsAsFactors=F)
names(expr2.gene.annotation) = c("probe","gene")

# set phenotype matrix and annotation
# main objects: 
##### pheno2: phenotype matrix
##### pheno2.samp.annotation: column/sample annotation for pheno2; mouse, strain, sex
##### note reshape warning: "multiple rows match for mouse_number=339: first taken"
#reshape pheno0 using the mouse number as column heads, the keys, and the values are the the expression value, row names are trait name
pheno1 = reshape(pheno0[,c(1,2,5)], idvar="trait_name", timevar="mouse_number", direction="wide")

#make column names just the mouse number
pheno.mouse.nums = sapply(names(pheno1)[2:ncol(pheno1)],function(x){strsplit(x,"value.")[[1]][2]})
pheno.mouse.strains = var.annot(vecdat=pheno.mouse.nums, annot=annotation, invar="mouse_number", outvar="Strain")
pheno.mouse.sex = var.annot(vecdat=pheno.mouse.nums, annot=annotation, invar="mouse_number", outvar="Sex")

# phenotype expressions
pheno2 = pheno1[,2:ncol(pheno1)]

#reset rownmaes into the trait_name
rownames(pheno2) = pheno1$trait_name

#column names are the split mouse numbers
names(pheno2) = pheno.mouse.nums

#create an annotation matrix with mouse number, strain, and sex
pheno2.samp.annotation = cbind(pheno.mouse.nums, pheno.mouse.strains, pheno.mouse.sex) %>% as.data.frame(stringsAsFactors=F)
names(pheno2.samp.annotation) = c("mouse_number","strain","sex")



############################################################
## select strains for analysis
############################################################

# identify strains with both M and F - expression and phenotype
#for the phenotype and expression data, select the strains corresponding to male and females
# remove duplicates
#using intersection: find a strain that has a phenotype associated with it and has male and female
strainsM.expr = expr2.samp.annotation %>% filter(sex == "M") %>% select(strain) %>% unique
strainsM.pheno = pheno2.samp.annotation %>% filter(sex == "M") %>% select(strain) %>% unique
strainsF.expr = expr2.samp.annotation %>% filter(sex == "F") %>% select(strain) %>% unique
strainsF.pheno = pheno2.samp.annotation %>% filter(sex == "F") %>% select(strain) %>% unique
strainsMF = Reduce(intersect, list(strainsM.expr$strain,strainsM.pheno$strain,strainsF.expr$strain,strainsF.pheno$strain))

# filter expression data for intersect strains
# main objects: 
##### expr3: gene expression matrix
##### expr3.samp.annotation: column/sample annotation for expr3; mouse, strain, sex
##### expr2.gene.annotation: row/probe annotation for expr3; probe, gene name
#returns true or false for each annotation saying whether or not the strain has male and females
expr.strain.inds = expr2.samp.annotation$strain %in% strainsMF

#expr 3 annotations only include strains with both male and female
expr3.samp.annotation = expr2.samp.annotation[expr.strain.inds,]

#exressions for strains with M and F and unique
expr3 = expr2[,expr.strain.inds]

# do same thing you did for expr. selecting the ones that are in MF intersect
# main objects: 
##### pheno3: phenotype matrix
##### pheno3.samp.annotation: column/sample annotation for pheno3; mouse, strain, sex
pheno.strain.inds = pheno2.samp.annotation$strain %in% strainsMF
pheno3.samp.annotation = pheno2.samp.annotation[pheno.strain.inds,]
pheno3 = pheno2[,pheno.strain.inds]

# modify organization
# expr and pheno data are organized identically by mouse number
# main objects: 
##### pheno4: phenotype matrix
##### pheno4.samp.annotation: column/sample annotation for pheno4; mouse, strain, sex
expr.nums = expr3.samp.annotation$mouse_number

#gives the indices where the mouse number of pheno is equal to the mouse number from expr
pheno.inds = sapply(expr.nums,function(x){which(pheno3.samp.annotation$mouse_number==x)}) %>% unlist

# edits pheno3 annotation to only include annotations for the mouse numbers we have expressions for
pheno4.samp.annotation = pheno3.samp.annotation[pheno.inds,]

#same thing as above but for phenotype values
pheno4 = pheno3[,pheno.inds]

# check data organization
# note that the annotation files are identical now
all(expr3.samp.annotation$mouse_number == pheno4.samp.annotation$mouse_number)

# expr 3 has row names as the probe, col names as mouse number and then expression values
#pheno 4 has row names= trait/diet, col names=mouse number and then weight?
all(names(expr3) == names(pheno4))
all(names(pheno4) == expr3.samp.annotation$mouse_number)
samp.annotation = pheno4.samp.annotation

#save.image("SIprocdata.RData")
#load("SIprocdata.RData")


#########################################################
##filtering phenotype data
########################################################

### filter NAs
pheno5 <- t(pheno4)
remove.na.func= function(pheno=NULL,na.prop=0.3){
  res=c()
  column.names=c()
  for (ii in 1:ncol(pheno)){
    x=pheno[,ii] %>% data.matrix() %>% as.numeric()
    # if(ii==3){break}
    if(length(which(is.na(x)==TRUE)) < na.prop*nrow(pheno)){
      res=cbind(res, pheno[,ii])
      new.column.name=colnames(pheno)[ii]
      column.names=c(column.names, new.column.name)
      column.names=t(column.names)
    }
  }
  #pheno6 <- t(res)
  colnames(res)=column.names
  pheno6 <- res 
  return(pheno6)
}
pheno6=remove.na.func(pheno=pheno5, na.prop=0.3)

# impute missing data
library(impute)
pheno7=impute.knn(pheno6) 



# phenotype correlation
library(NMF)
phenocor = cor(pheno7,method="pearson", use="pairwise.complete.obs")
pdf("heat_20181108_Innis.pdf", width=10, height=10)
aheatmap(phenocor,breaks=0)
dev.off()

hm = aheatmap(phenocor,breaks=0)
clusters = colnames(phenocor)[rev(hm$colInd)]
write.table(clusters,"clusters.txt",sep="\t",quote=F,col.names=F,row.names=F)
