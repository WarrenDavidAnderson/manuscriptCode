

######################################################
## move deseq data
######################################################

from=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/meme_degDat/*meme*
to=wa3j@interactive.hpc.virginia.edu:/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/deseq_res
scp -r $from $to


######################################################
## get all individual bed files in the same place
######################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/deseq_res
dir=/scratch/wa3j/deseq_res/meme_degDat
cd ${dir}

folders=$(ls -d */)
for ii in ${folders}
do
cd ${ii}
cp *.bed ..
cd ..
done

######################################################
## R script to change file names and save name key
# generates mapping between 1-20 and time comparisons
######################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/deseq_res
dir=/scratch/wa3j/deseq_res/meme_degDat
cd ${dir}

# load modules
module load gcc/7.1.0 openmpi/2.1.5 R/3.5.3

# files with condition ids
f.dn = list.files(pattern="downsig")

# loop through files, change file names, save key
cnt = 1
key = c()
for(ii in f.dn){

	condit = strsplit(ii,"sig_")[[1]][2]
	condit = strsplit(condit,".bed")[[1]][1]
	
	comm1 = paste0("cat downsig_",condit,".bed > downsig_",cnt,".bed")
	comm2 = paste0("cat upsig_",condit,".bed > upsig_",cnt,".bed")
	comm3 = paste0("cat unsig_",condit,".bed > unsig_",cnt,".bed")
	system(comm1)
	system(comm2)
	system(comm3)
	
	new = c(condit, cnt)
	key = rbind(key, new)
	cnt = cnt + 1

}
colnames(key) = c("comparison","id")
rownames(key) = key[,2]
key = as.data.frame(key,stringsAsFactors=F)

# output key
fname = "memeTimeKey.txt"
write.table(key,fname,sep="\t",col.names=T,row.names=F,quote=F)

# separate files based on identifier
system("mkdir array") # ordinal identifier for slurm array
for(ii in 1:nrow(key)){
	fname = list.files(pattern=paste0("sig_",ii,".bed"))
	comm = paste0("mv -t array ",paste(fname,collapse=" "))
	system(comm)
}









