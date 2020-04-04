
#########################################
## unzip the file
#########################################

PATH=$PATH:/h4/t1/users/wa3j/software/pigz-2.4
cd /scratch/wa3j/v8/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1
tar -xvf phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1.GRU.tar
cd phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1
fname=GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_AllVar_QC_metrics.vcf.gz
pigz -dk ${fname}


cd /scratch/wa3j/v8/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1

#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH -p largemem
#SBATCH --time=96:00:00
#SBATCH --mem=975G
#SBATCH --cpus-per-task=16

gunzip GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_AllVar_QC_metrics.vcf.gz

# move data for analysis using powhatan
dir=/scratch/wa3j/v8/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1
file=GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_AllVar_QC_metrics.vcf
to=/nv/vol192/civeleklab/warren/adipose_project/eqtlLD
mv $dir/${file} ${to}

from=/nv/vol192/civeleklab/mete/GTEx/GTEx_v8/geno/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_AllVar_QC_metrics.vcf.gz


#########################################
## perform LD calculations
#########################################

# use the unfiltered vcf file to evaluate LD
cd /m/civeleklab/civeleklab/warren/adipose_project/eqtlLD
PATH=$PATH:/h4/t1/apps/statgen/plink1.90b6.7

# evaluate LD
dir=/m/civeleklab/civeleklab/warren/adipose_project/eqtlLD
fname=GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_AllVar_QC_metrics.vcf.gz
file=${dir}/${fname}
plink --noweb --vcf ${file} --ld-window-kb 2000 --ld-window 999999 --r2 --ld-window-r2 0.8 --out gtex_geno_LD08

# subset the LD file
awk -v OFS='\t' '{print $3, $6, $7}' gtex_geno_LD08.ld > gtex_geno_LD08.txt
awk -v OFS='\t' '{print $3, $6, $7}' gtex_geno_LD06.ld > gtex_geno_LD06.txt

#########################################
## get LD SNPs for eQTL SNPs
#########################################

# move the association data
from=/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl/CIdata.cases.RData
to=wa3j@interactive.hpc.virginia.edu:/nv/vol192/civeleklab/warren/adipose_project/eqtlLD
scp -r $from $to

# set up files for slurm array analysis
module load gcc/7.1.0 openmpi/3.1.4 R/3.5.3
cd /nv/vol192/civeleklab/warren/adipose_project/eqtlLD
load("CIdata.cases.RData")
narray = 200
cntper = floor(nrow(CIdata.cases) / narray)
for(ii in 1:narray){
	start = (ii-1)*cntper + 1
	end = start + cntper - 1
	if(ii == narray){
		end = nrow(CIdata.cases)
	}
	inds = start:end
	out = CIdata.cases$snps[inds]
	fname = paste0("array",ii)
	write.table(out,fname,col.names=F,row.names=F,quote=F)
}

############################
## main slurm script

#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --mem 96000
#SBATCH --partition=standard

module load gcc/7.1.0 openmpi/3.1.4 R/3.5.3
cd /nv/vol192/civeleklab/warren/adipose_project/eqtlLD

# sbatch --array=1-200 main.sh
cnt=${SLURM_ARRAY_TASK_ID}
sh run.sh ${cnt}

############################
## run slurm script

#!/bin/bash

# chmod u+x *.sh *.R

cnt="$1"
newdir=res${cnt}
mkdir ${newdir}
cp -t ${newdir} run.R array${cnt} 
cd ${newdir}
Rscript --vanilla run.R array${cnt}

############################
## run R script

# ld file name
fld = "/nv/vol192/civeleklab/warren/adipose_project/eqtlLD/gtex_geno_LD08.txt"

# import data
args = commandArgs(trailingOnly=TRUE)
namen = args
esnps = read.table(namen, header=F, stringsAsFactors=F)
esnps = esnps[,1]

# loop through esnps and isolate ld data
out = c()
for(ii in 1:length(esnps)){

	print(ii)

	snp = esnps[ii]
	comm = paste0("cat ",fld," | LC_ALL=C fgrep ",snp," > ldii.txt")
	comm = paste0("LC_ALL=C fgrep ",snp," ",fld," > ldii.txt")
  	system(comm)

	info = file.info("ldii.txt")
	empty = info$size

	if(empty == 0){
		system("rm ldii.txt")
		next
	}

	new = read.table("ldii.txt",sep="\t",header=F,stringsAsFactors=F)
	out = rbind(out, new)
	system("rm ldii.txt")

}

# output the results
cnt = strsplit(namen,"array")[[1]][2]
fname = paste0("ld_",cnt,".txt")
write.table(out,fname,col.names=F,row.names=F,quote=F,sep="\t")


############################
## aggregate results

# data directory
cd /nv/vol192/civeleklab/warren/adipose_project/eqtlLD

# loop through all folders with data
folds=$(ls -d *res*/)
for ii in ${folds}
do
cd ${ii}
cp *.txt ..
cd ..
done

# combine all of the data
fname=LD08.txt
echo snpA snpB R2 | awk -v OFS='\t' '{print $1, $2, $3}' > header
cat header *ld_*.txt > ${fname}
rm *ld_*.txt header


############################
## move data for local analysis

# data directory
cd /nv/vol192/civeleklab/warren/adipose_project/eqtlLD

from=wa3j@interactive.hpc.virginia.edu:/nv/vol192/civeleklab/warren/adipose_project/eqtlLD/LD06.txt
to=/media/wa3j/Seagate2/Documents/adipose_sex_ms/final/gtex_eqtl
scp -r $from $to






