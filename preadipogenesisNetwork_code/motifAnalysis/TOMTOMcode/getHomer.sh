

dir=/home/wa3j/meme-5.0.3/motif_databases/WAdbs
cd ${dir}

# get the homer motifs
cd /home/wa3j/meme-5.0.3/motif_databases/WAdbs/homer
wget http://homer.ucsd.edu/homer/custom.motifs

# specify name for home motif file
outfile=HOMER.meme
infile=custom.motifs

# specify motifs
motifsids=( $(cat custom.motifs | grep ">" | awk -F" " '{print $2}') )
motifsseq=( $(cat custom.motifs | grep ">" | awk -F" " '{print $1}' | awk -F">" '{print $2}') )
nmotifs=$(cat custom.motifs | grep ">" | wc -l)

# create header file (header)
printf '%s\n' 'MEME version 5' '' 'ALPHABET= ACGT' '' 'strands: + -'  '' > header1
printf '%s\n' 'Background letter frequencies (from uniform background):' > header2
printf '%s\n' 'A 0.25000 C 0.25000 G 0.25000 T 0.25000' '' > header3
cat header1 header2 header3 > ${outfile}
rm header1 header2 header3

# loop through motifs and generate PWM file
for ii in `seq 0 $(expr ${nmotifs} - 1)`
do
tf=$(echo ${motifsids[ii]})
id=$(echo ${motifsseq[ii]})
linestart=$(grep -w -n ${tf} ${infile} | grep -w ${id} | awk -F":" '{print $1}') 
linestart=$(expr ${linestart} + 1)
tail -n +${linestart} ${infile} > topWM
lineend=$(cat topWM | grep -n ">" | head -1 | awk -F":" '{print $1}')
if [ -z "$lineend" ] 
then 
lineend=$(wc -l topWM | awk -F" " '{print $1}')
lineend=$(expr ${lineend} + 1)
fi
head -$(expr ${lineend} - 1) topWM > PWM 
echo MOTIF ${tf}_${id} > line1
echo > line3
echo letter-probability matrix: > line2
mv ${outfile} tmp
cat tmp line1 line2 PWM line3 > ${outfile}
rm topWM PWM line1 line2 line3 tmp
done


#################
# homer motif file
fi=/home/wa3j/meme-5.0.3/motif_databases/WAdbs/homer/HOMER.meme
cp ${fi} /home/wa3j/meme-5.0.3/motif_databases/WAdbs
chmod u+x *meme

# copy over to the total set of all meme motifs
db=/scratch/wa3j/adipogenesis/orphan/alldb
cp ${fi} ${db}



