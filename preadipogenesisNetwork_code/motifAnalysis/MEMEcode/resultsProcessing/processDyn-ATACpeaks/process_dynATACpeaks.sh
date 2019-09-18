

# directory and main files
dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/meme_dynDat
cd ${dir}


####################################################
## after the meme loop has run, collect all meme results
####################################################

# main directory
dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/meme_dynDat
cd ${dir}

# directory for meme data
datdir=meme_pairwise_memeResults
mkdir ${datdir}

# folders with meme data
folders=$(ls -d *id_*/)

# loop through all data and aggregate meme.txt results
for ii in ${folders}
do
	
cd ${ii}

# loop through each subfolder and copy the Evalue data
folds=$(ls -d */)
for jj in ${folds}
do
cd ${jj}
id=$(echo ${jj} | awk -F"/" '{print $1}')	
cat meme.txt > ${id}_meme.txt
cp ${id}_meme.txt ${dir}/${datdir}
cd .. 
done

cd ..
done



####################################################
## isolate all unique motifs with E-values
## that are less than a set threshold
####################################################

subdir=meme_pairwise_memeResults
dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/meme_dynDat
cd ${dir}/${subdir}

# document meme run conditions
comp=dynamic_meme

# isolate E-value data for each motif
for ii in *.txt
do
cond=$(echo ${ii} | awk -F"_meme.txt" '{print $1}')
cat ${ii} | grep E-value > ${comp}_${cond}_Eval.txt
done

# separate motif Evalue data and total meme output files
mkdir memeEval
mkdir memeAll
mv *Eval.txt* memeEval
mv *meme.txt* memeAll

cd ${dir}/${subdir}/memeEval

# remove all "Stopped because motif E-value >" lines
for ii in *.txt
do
sed -i "/\b\Stopped\b/d" ${ii}
sed -i "/\b\product\b/d" ${ii}
sed -i "/\b\objective\b/d" ${ii}
done

# first use R to loop through data and set identity all motifs of interest
module load gcc/7.1.0 openmpi/2.1.5 R/3.5.3

Eval.thresh = 0.1 # E-value threshold

# loop to get all significant motifs
files = list.files(getwd()) 
info = file.info(files)
nonempty = rownames(info[info$size != 0, ]) 
Edat = matrix(NA,nrow=0,ncol=6)
for(ii in nonempty){
	cond = strsplit(ii,"_Eval.txt")[[1]][1]
	datii = read.table(ii,stringsAsFactors=FALSE)
	indsig = which(datii[,15] < Eval.thresh)
	if(length(indsig)==0){next}
	new = cbind(cond, datii[indsig,2], datii[indsig,6], datii[indsig,9], 
			datii[indsig,12], datii[indsig,15])
	Edat = rbind(Edat, matrix(new,length(indsig),6))
}
Edat = as.data.frame(Edat,stringsAsFactors=FALSE)
names(Edat) = c("condition","motif","width","nsites","llr","Eval")
Edat[,c(3:6)] = apply(Edat[,c(3:6)],2,function(x){as.numeric(data.matrix(x))})

# isolate unique motifs
# document lowest eval and highest llr
motifs.all = Edat
motifs.unq = matrix(0,0,7)
while( nrow(motifs.all)>0 ){
	motif = motifs.all$motif[1]
	ind = which(motifs.all$motif == motif)
	condits = unique(motifs.all$condition[ind])
	con = paste(condits, collapse=",")
	ind.min = which(motifs.all$Eval[ind] == min(motifs.all$Eval[ind]))
	ind.max = which(motifs.all$llr[ind][ind.min] == max(motifs.all$llr[ind][ind.min]))
 	new = cbind(motifs.all[ind[ind.min][ind.max],c(2:6)], length(ind), con)
	motifs.unq = rbind(motifs.unq, new)
	motifs.all = motifs.all[-ind,]
}
motifs.unq = as.data.frame(motifs.unq,stringsAsFactors=FALSE)
names(motifs.unq) = c("motif","width","nsites","llr","Eval","ncond","condits")
nrow(Edat) # 307 motifs
nrow(motifs.unq) # 263 unique motifs
motifs.unq$condits = as.character(motifs.unq$condits)

# save.image("meme.res.for.tomtom.RData")
# load("meme.res.for.tomtom.RData")

# generate output file for generating a db of meme motifs
fname = "all_sig_motifs_dynamics.txt"
write.table(Edat,fname,quote=F,sep="\t",col.names=T,row.names=F)

####################################################
## isolate key results for downstream analysis
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/meme_dynDat/meme_pairwise_memeResults
cd ${dir}
mkdir pswm

# directory for meme pswms to be used with tomtom
# mkdir pswm
cp ${dir}/memeEval/*dynamics.txt ${dir}/pswm
cp ${dir}/memeEval/*RData ${dir}/pswm
cd ${dir}/pswm


-bash-4.2$pwd
/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/meme_dynDat/meme_pairwise_memeResults/pswm













