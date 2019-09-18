
####################################################
## generate pswm file with all sig meme motifs for tomtom
####################################################


####################################################
## get pswms for the differential peak motifs
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motifdb
cd ${dir}

# locations for differential ATAC peak meme motif data
f1=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/deseq_res/bkg_0.1_classic/bkg_0.1_classic_pairwise_memeResults/memeAll
f2=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/deseq_res/bkg_0.1_de/bkg_0.1_de_pairwise_memeResults/memeAll
f3=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/deseq_res/nobkg_0.1_classic/nobkg_0.1_classic_pairwise_memeResults/memeAll
f4=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/deseq_res/nobkg_0.1_de/nobkg_0.1_de_pairwise_memeResults/memeAll

# motif db file
outfile=denovoPWM.txt

# create header file (header)
printf '%s\n' 'MEME version 5' '' 'ALPHABET= ACGT' '' 'strands: + -'  '' > header1
printf '%s\n' 'Background letter frequencies (from uniform background):' > header2
printf '%s\n' 'A 0.25000 C 0.25000 G 0.25000 T 0.25000' '' > header3
cat header1 header2 header3 > ${outfile}
rm header1 header2 header3

#################
# loop through background classic meme motifs (f1)
subdir=${f1}
condit=deg_bkg_classic_

# loop through all motif ids and add the PWM if it is there
ind=1
for ii in ${subdir}/*.txt
do
file=$(echo ${ii} | awk -F"${subdir}/" '{print $2}')
cond=$(echo ${file} | awk 'BEGIN{FS="_"; OFS="_"} {print $1,$2}')
motifs=$(cat ${ii} | grep MOTIF | grep E-value | awk -F"MOTIF " '{print $2}' | awk -F" MEME" '{print $1}')
for jj in ${motifs}
do
linenum=$(cat ${ii} | grep -n -w ${jj} | grep probability | awk -F":" '{print $1}')
for kk in ${linenum}
do
if [ ! -z "${kk}" ]; then
	linenum=$(expr ${kk} + 2)
	tail -n +${linenum} ${ii} > topWM
	lineend=$(cat topWM | grep -n "\------" | head -1 | awk -F":" '{print $1}')
	head -$(expr ${lineend} - 1) topWM > PWM
	printf '%s\n' '' 'MOTIF '${jj}_${condit}${cond}_${ind} > headerWM
	mv ${outfile} tmp
	cat tmp headerWM PWM > ${outfile}
	rm tmp headerWM PWM topWM
	ind=$(expr ${ind} + 1)
fi
done
done
done

#################
# loop through background differential meme motifs (f2)
subdir=${f2}
condit=deg_bkg_de_
for ii in ${subdir}/*.txt
do
file=$(echo ${ii} | awk -F"${subdir}/" '{print $2}')
cond=$(echo ${file} | awk 'BEGIN{FS="_"; OFS="_"} {print $1,$2}')
motifs=$(cat ${ii} | grep MOTIF | grep E-value | awk -F"MOTIF " '{print $2}' | awk -F" MEME" '{print $1}')
for jj in ${motifs}
do
linenum=$(cat ${ii} | grep -n -w ${jj} | grep probability | awk -F":" '{print $1}')
for kk in ${linenum}
do
if [ ! -z "${kk}" ]; then
	linenum=$(expr ${kk} + 2)
	tail -n +${linenum} ${ii} > topWM
	lineend=$(cat topWM | grep -n "\------" | head -1 | awk -F":" '{print $1}')
	head -$(expr ${lineend} - 1) topWM > PWM
	printf '%s\n' '' 'MOTIF '${jj}_${condit}${cond}_${ind} > headerWM
	mv ${outfile} tmp
	cat tmp headerWM PWM > ${outfile}
	rm tmp headerWM PWM topWM
	ind=$(expr ${ind} + 1)
fi
done
done
done

#################
# loop through no-background classic meme motifs (f3)
subdir=${f3}
condit=deg_nobkg_classic_
for ii in ${subdir}/*.txt
do
file=$(echo ${ii} | awk -F"${subdir}/" '{print $2}')
cond=$(echo ${file} | awk 'BEGIN{FS="_"; OFS="_"} {print $1,$2}')
motifs=$(cat ${ii} | grep MOTIF | grep E-value | awk -F"MOTIF " '{print $2}' | awk -F" MEME" '{print $1}')
for jj in ${motifs}
do
linenum=$(cat ${ii} | grep -n -w ${jj} | grep probability | awk -F":" '{print $1}')
for kk in ${linenum}
do
if [ ! -z "${kk}" ]; then
	linenum=$(expr ${kk} + 2)
	tail -n +${linenum} ${ii} > topWM
	lineend=$(cat topWM | grep -n "\------" | head -1 | awk -F":" '{print $1}')
	head -$(expr ${lineend} - 1) topWM > PWM
	printf '%s\n' '' 'MOTIF '${jj}_${condit}${cond}_${ind} > headerWM
	mv ${outfile} tmp
	cat tmp headerWM PWM > ${outfile}
	rm tmp headerWM PWM topWM
	ind=$(expr ${ind} + 1)
fi
done
done
done


#################
# loop through no-background differential meme motifs (f4)
subdir=${f4}
condit=deg_nobkg_de_
for ii in ${subdir}/*.txt
do
file=$(echo ${ii} | awk -F"${subdir}/" '{print $2}')
cond=$(echo ${file} | awk 'BEGIN{FS="_"; OFS="_"} {print $1,$2}')
motifs=$(cat ${ii} | grep MOTIF | grep E-value | awk -F"MOTIF " '{print $2}' | awk -F" MEME" '{print $1}')
for jj in ${motifs}
do
linenum=$(cat ${ii} | grep -n -w ${jj} | grep probability | awk -F":" '{print $1}')
for kk in ${linenum}
do
if [ ! -z "${kk}" ]; then
	linenum=$(expr ${kk} + 2)
	tail -n +${linenum} ${ii} > topWM
	lineend=$(cat topWM | grep -n "\------" | head -1 | awk -F":" '{print $1}')
	head -$(expr ${lineend} - 1) topWM > PWM
	printf '%s\n' '' 'MOTIF '${jj}_${condit}${cond}_${ind} > headerWM
	mv ${outfile} tmp
	cat tmp headerWM PWM > ${outfile}
	rm tmp headerWM PWM topWM
	ind=$(expr ${ind} + 1)
fi
done
done
done


####################################################
## get pswms for the differential dynamics motifs
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motifdb
cd ${dir}

# locations for all meme motif comparison data
f1=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/meme_dynDat/meme_pairwise_memeResults/memeAll

#################
# loop through dynamic
subdir=${f1}
condit=dyn_
for ii in ${subdir}/*.txt
do
file=$(echo ${ii} | awk -F"${subdir}/" '{print $2}')
cond=$(echo ${file} | awk 'BEGIN{FS="_"; OFS="_"} {print $2,$3,$4}')
motifs=$(cat ${ii} | grep MOTIF | grep E-value | awk -F"MOTIF " '{print $2}' | awk -F" MEME" '{print $1}')
for jj in ${motifs}
do
linenum=$(cat ${ii} | grep -n -w ${jj} | grep probability | awk -F":" '{print $1}')
for kk in ${linenum}
do
if [ ! -z "${kk}" ]; then
	linenum=$(expr ${kk} + 2)
	tail -n +${linenum} ${ii} > topWM
	lineend=$(cat topWM | grep -n "\------" | head -1 | awk -F":" '{print $1}')
	head -$(expr ${lineend} - 1) topWM > PWM
	printf '%s\n' '' 'MOTIF '${jj}_${condit}${cond}_${ind} > headerWM
	mv ${outfile} tmp
	cat tmp headerWM PWM > ${outfile}
	rm tmp headerWM PWM topWM
	ind=$(expr ${ind} + 1)
fi
done
done
done


####################################################
## get pswms for the differential dreg motifs
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motifdb
cd ${dir}

# locations for all meme motif comparison data
f1=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/dreg_for_meme/bkg_0.1_de/bkg_0.1_de_pairwise_memeResults/memeAll
f2=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/dreg_for_meme/bkg_0.1_classic/bkg_0.1_classic_pairwise_memeResults/memeAll
f3=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/dreg_for_meme/nobkg_0.1_classic/nobkg_0.1_classic_pairwise_memeResults/memeAll
f4=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/dreg_for_meme/nobkg_0.1_de/nobkg_0.1_de_pairwise_memeResults/memeAll

#################
# loop through bkg_de motifs (f1)
subdir=${f1}
condit=dreg_bkg_de_
for ii in ${subdir}/*.txt
do
file=$(echo ${ii} | awk -F"${subdir}/" '{print $2}')
cond=$(echo ${file} | awk 'BEGIN{FS="_"; OFS="_"} {print $1,$2}')
motifs=$(cat ${ii} | grep MOTIF | grep E-value | awk -F"MOTIF " '{print $2}' | awk -F" MEME" '{print $1}')
for jj in ${motifs}
do
linenum=$(cat ${ii} | grep -n -w ${jj} | grep probability | awk -F":" '{print $1}')
for kk in ${linenum}
do
if [ ! -z "${kk}" ]; then
	linenum=$(expr ${kk} + 2)
	tail -n +${linenum} ${ii} > topWM
	lineend=$(cat topWM | grep -n "\------" | head -1 | awk -F":" '{print $1}')
	head -$(expr ${lineend} - 1) topWM > PWM
	printf '%s\n' '' 'MOTIF '${jj}_${condit}${cond}_${ind} > headerWM
	mv ${outfile} tmp
	cat tmp headerWM PWM > ${outfile}
	rm tmp headerWM PWM topWM
	ind=$(expr ${ind} + 1)
fi
done
done
done

#################
# loop through bkg_cl motifs (f2)
subdir=${f2}
condit=dreg_bkg_classic_
for ii in ${subdir}/*.txt
do
file=$(echo ${ii} | awk -F"${subdir}/" '{print $2}')
cond=$(echo ${file} | awk 'BEGIN{FS="_"; OFS="_"} {print $1,$2}')
motifs=$(cat ${ii} | grep MOTIF | grep E-value | awk -F"MOTIF " '{print $2}' | awk -F" MEME" '{print $1}')
for jj in ${motifs}
do
linenum=$(cat ${ii} | grep -n -w ${jj} | grep probability | awk -F":" '{print $1}')
for kk in ${linenum}
do
if [ ! -z "${kk}" ]; then
	linenum=$(expr ${kk} + 2)
	tail -n +${linenum} ${ii} > topWM
	lineend=$(cat topWM | grep -n "\------" | head -1 | awk -F":" '{print $1}')
	head -$(expr ${lineend} - 1) topWM > PWM
	printf '%s\n' '' 'MOTIF '${jj}_${condit}${cond}_${ind} > headerWM
	mv ${outfile} tmp
	cat tmp headerWM PWM > ${outfile}
	rm tmp headerWM PWM topWM
	ind=$(expr ${ind} + 1)
fi
done
done
done


#################
# loop through nobkg_cl motifs (f3)
subdir=${f3}
condit=dreg_nobkg_classic_
for ii in ${subdir}/*.txt
do
file=$(echo ${ii} | awk -F"${subdir}/" '{print $2}')
cond=$(echo ${file} | awk 'BEGIN{FS="_"; OFS="_"} {print $1,$2}')
motifs=$(cat ${ii} | grep MOTIF | grep E-value | awk -F"MOTIF " '{print $2}' | awk -F" MEME" '{print $1}')
for jj in ${motifs}
do
linenum=$(cat ${ii} | grep -n -w ${jj} | grep probability | awk -F":" '{print $1}')
for kk in ${linenum}
do
if [ ! -z "${kk}" ]; then
	linenum=$(expr ${kk} + 2)
	tail -n +${linenum} ${ii} > topWM
	lineend=$(cat topWM | grep -n "\------" | head -1 | awk -F":" '{print $1}')
	head -$(expr ${lineend} - 1) topWM > PWM
	printf '%s\n' '' 'MOTIF '${jj}_${condit}${cond}_${ind} > headerWM
	mv ${outfile} tmp
	cat tmp headerWM PWM > ${outfile}
	rm tmp headerWM PWM topWM
	ind=$(expr ${ind} + 1)
fi
done
done
done

#################
# loop through nobkg_de motifs (f4)
subdir=${f4}
condit=dreg_nobkg_de_
for ii in ${subdir}/*.txt
do
file=$(echo ${ii} | awk -F"${subdir}/" '{print $2}')
cond=$(echo ${file} | awk 'BEGIN{FS="_"; OFS="_"} {print $1,$2}')
motifs=$(cat ${ii} | grep MOTIF | grep E-value | awk -F"MOTIF " '{print $2}' | awk -F" MEME" '{print $1}')
for jj in ${motifs}
do
linenum=$(cat ${ii} | grep -n -w ${jj} | grep probability | awk -F":" '{print $1}')
for kk in ${linenum}
do
if [ ! -z "${kk}" ]; then
	linenum=$(expr ${kk} + 2)
	tail -n +${linenum} ${ii} > topWM
	lineend=$(cat topWM | grep -n "\------" | head -1 | awk -F":" '{print $1}')
	head -$(expr ${lineend} - 1) topWM > PWM
	printf '%s\n' '' 'MOTIF '${jj}_${condit}${cond}_${ind} > headerWM
	mv ${outfile} tmp
	cat tmp headerWM PWM > ${outfile}
	rm tmp headerWM PWM topWM
	ind=$(expr ${ind} + 1)
fi
done
done
done



####################################################
## generate pwm files for each de novo motif
####################################################

dir=/nv/vol192/civeleklab/warren/MGlab/ATAC_WAFD/3T3_ATAC1-3/motifs/meme/tomtom_denovo/motifdb
cd ${dir}

# location for the new PWM files
mkdir alldenovo
cp denovoPWM.txt alldenovo
cd alldenovo

# isolate de novo motif ids and specify the original motif db file
motiffile=denovoPWM.txt
motifs=$(cat ${motiffile} | grep MOTIF | awk -F"MOTIF " '{print $2}')

# loop through de novo motifs and create pwm files
for ii in ${motifs}
do
linenum=$(cat ${motiffile} | grep -n -w ${ii} | awk -F":" '{print $1}')
line=$(expr ${linenum} + 2)
tail -n +${line} ${motiffile} > part
lineend=$(grep -E --line-number --with-filename '^$' part | head -1 | awk -F":" '{print $2}')
if [ -z "${lineend}" ]; then
endline=$(wc -l part | awk -F" " '{print $1}')
else
endline=$(expr ${lineend} - 1)
fi
head -${endline} part > tmp
awk '{OFS="\t";} {print $1,$2,$3,$4}' tmp > ${ii}.txt
rm part tmp
done

rm denovoPWM.txt run.sh *.out


## next analysis: tomtom_denovo_vs_db_array.sh

