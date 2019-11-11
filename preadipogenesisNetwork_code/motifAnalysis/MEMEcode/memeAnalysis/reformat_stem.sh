

######################################################
## R script to change file names and save name key
# revise profile identifiers in R (rivanna)
######################################################

cd /nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/meme_dynDat/dyn_orig 
mkdir newFiles

module load gcc/7.1.0 openmpi/3.1.4 R/3.5.3

# directory info
main.dir = paste0("/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/meme_dynDat/dyn_orig")
new.dir = paste0(main.dir,"/newFiles")
setwd(main.dir)

# sort original profile ids
files = list.files(pattern = "sigDyn")
prof = sapply(files,function(x){strsplit(x,"_profile")[[1]][2]})
prof = sapply(prof,function(x){strsplit(x,".bed")[[1]][1]})
prof.sort = sort( as.numeric(prof) )

# loop through files, change file names, save key
setwd(new.dir)
cnt = 1
key = c()
for(ii in prof.sort){

	sg.old = paste0(main.dir,"/preadip_sigDynamics_profile",ii,".bed")
	sg.new = paste0("preadip_sigDynamics_profile",cnt,".bed")
	un.old = paste0(main.dir,"/preadip_unDynamics_profile",ii,".bed")
	un.new = paste0("preadip_unDynamics_profile",cnt,".bed")
	
	comm1 = paste0("cat ",sg.old," > ",sg.new)
	comm2 = paste0("cat ",un.old," > ",un.new)
	system(comm1)
	system(comm2)
	
	new = c(ii, cnt)
	key = rbind(key, new)
	cnt = cnt + 1

}
colnames(key) = c("origprofile","id")
rownames(key) = key[,2]
key = as.data.frame(key,stringsAsFactors=F)

# write out the key
fname = "bedKey_STEM.txt"
write.table(key,fname,sep="\t",quote=F,col.names=T,row.names=F)

mkdir oldFiles
mv *.bed oldFiles


