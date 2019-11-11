
####################################################
## Classification of Transcription Factors in Mammalia
## http://tfclass.bioinf.med.uni-goettingen.de/
## http://hocomoco10.autosome.ru/motif
## http://cisbp.ccbr.utoronto.ca
####################################################

##################################################################
# community_42
# 1.1 Basic leucine zipper factors (bZIP)
# 1.1.1 Jun-related
cd comm42
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.1.1_mammalia_dbd_fasta.fasta
cat 1.1.1_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-annot" '{print $1}' > comm42_genes.txt

##################################################################
# community_67
# community_69
# community_157
# community_268
# 1.1 Basic leucine zipper factors (bZIP)
# 1.1.2 Fos-related
cd comm67
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.1.2_mammalia_dbd_fasta.fasta
cat 1.1.2_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-annot" '{print $1}' > comm67_genes.txt

##################################################################
# community_58
# community_90
# community_96
# community_159
# community_166
# 1.1 Basic leucine zipper factors (bZIP)
# 1.1.3 Maf-related
cd comm58
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.1.3_mammalia_dbd_fasta.fasta
cat 1.1.3_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-annot" '{print $1}' > comm58_genes.txt

##################################################################
# community_23
# 1.1 Basic leucine zipper factors (bZIP)
# 1.1.8 CEBP-related
# CEBP, etc (convert lower case)
cd comm23
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.1.8_mammalia_dbd_fasta.fasta
cat 1.1.8_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm23_genes.txt

##################################################################
# community_47
# community_62
# community_63
# community_64
# community_66
# community_76
# community_88
# community_140
# community_364
# 1.2 Basic helix-loop-helix factors (bHLH)
# 1.2.1 E2A-related factors
# 1.2.2 MyoD/ASC-related factors
# 1.2.3 Tal-related
cd comm47
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.2.1_mammalia_dbd_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.2.2_mammalia_dbd_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.2.3_mammalia_dbd_fasta.fasta
cat *.fasta > Basic.fasta
cat Basic.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-annot" '{print $1}' |\
awk -F"annot" '{print $1}' | awk -F"_" '{print $1}' |\
awk -F"-" '{print $1}' | sort -u > comm47_genes.txt

##################################################################
# community_101
# 1.2 Basic helix-loop-helix factors (bHLH)
# 1.2.4 Hairy-related
cd comm101
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.2.4_mammalia_dbd_fasta.fasta
cat 1.2.4_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-annot" '{print $1}' |\
awk -F"annot" '{print $1}' | awk -F"_" '{print $1}' > comm101_genes.txt

##################################################################
# community_158
# 1.2 Basic helix-loop-helix factors (bHLH)
# 1.2.6 bHLH-ZIP
# 1.2.6.1 TFE3
cd comm158
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.2.6.1_mammalia_dbd_fasta.fasta
cat 1.2.6.1_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-annot" '{print $1}' > comm158_genes.txt

##################################################################
# community_44
# 1.2 Basic helix-loop-helix factors (bHLH)
# 1.2.6 bHLH-ZIP
# 1.2.6.5 MYC
cd comm44
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.2.6.5_mammalia_dbd_fasta.fasta
cat 1.2.6.5_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-annot" '{print $1}' > comm44_genes.txt

##################################################################
# community_4
# 1.3 Basic helix-span-helix factors (bHSH)
# 1.3.1 AP2
# Tfap2a,b,c,d,e (manually modify identifiers)
cd comm4
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.3.1_mammalia_dbd_fasta.fasta
cat 1.3.1_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm4_genes.txt

##################################################################
# community_77
# 2.1 Nuclear receptors with C4 zinc fingers
# 2.1.1 Steroid hormone receptors
# 2.1.1.1 GR-like(NR3C)
cd comm77
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.1.1.1_mammalia_dbd_fasta.fasta
cat 2.1.1.1_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm77_genes.txt

##################################################################
# community_356
# 2.1 Nuclear receptors with C4 zinc fingers
# 2.1.2 Thyroid hormone receptor-related factors
# 2.1.2.1 RAR(NR1B)
cd comm356
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.1.2.1_mammalia_dbd_fasta.fasta
cat 2.1.2.1_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' |\
awk -F"_ma" '{print $1}' > comm356_genes.txt

##################################################################
# community_338
# 2.1 Nuclear receptors with C4 zinc fingers
# 2.1.2 Thyroid hormone receptor-related factors
# 2.1.2.2 T3R(NR1A)
cd comm338
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.1.2.2_mammalia_dbd_fasta.fasta
cat 2.1.2.2_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm338_genes.txt

##################################################################
# community_46
# 2.1 Nuclear receptors with C4 zinc fingers
# 2.1.2 Thyroid hormone receptor-related factors
# 2.1.2.5 PPAR(NR1C)
cd comm46
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.1.2.5_mammalia_dbd_fasta.fasta
cat 2.1.2.5_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm46_genes.txt

##################################################################
# community_384
# 2.1 Nuclear receptors with C4 zinc fingers
# 2.1.3 RXR-related receptors
# 2.1.3.3 TLX(NR2E1)
cd comm384
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.1.3.3_mammalia_dbd_fasta.fasta
cat 2.1.3.3_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm384_genes.txt

##################################################################
# community_16
# 2.1 Nuclear receptors with C4 zinc fingers
# 2.1.4 NGFI(NR4A)
cd comm16
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.1.4_mammalia_dbd_fasta.fasta
cat 2.1.4_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm16_genes.txt

##################################################################
# community_79
# 2.3 C2H2 zinc finger factors
# 2.3.1 Three-zinc finger Krüppel-related
# 2.3.1.2 Kr-like
cd comm79
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.1.2_mammalia_dbd_fasta.fasta
cat 2.3.1.2_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm79_genes.txt

##################################################################
*** note see community_61 (2.3.3.x and others)
# community_244
# 2.3 C2H2 zinc finger factors
# 2.3.3 More than 3 adjacent zinc fingers
# 2.3.3.0.39
cd comm244
echo ZFP691 > comm244_genes.txt

##################################################################
*** note see community_61 (2.3.3.x and others)
# community_382
# 2.3 C2H2 zinc finger factors
# 2.3.3 More than 3 adjacent zinc fingers
# 2.3.3.0.127 ZNF416
cd comm382
echo ZNF416 > comm382_genes.txt

##################################################################
*** note see community_61 (2.3.3.x and others)
# community_134
# 2.3 C2H2 zinc finger factors
# 2.3.3 More than 3 adjacent zinc fingers
# 2.3.3.0.197 ZNF189
cd comm134
echo ZNF189 > comm134_genes.txt

##################################################################
*** note see community_61 (2.3.3.x and others)
# community_156
# 2.3 C2H2 zinc finger factors
# 2.3.3 More than 3 adjacent zinc fingers
# 2.3.3.11 ZBTB6-like
cd comm156
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.3.11_mammalia_dbd_fasta.fasta
cat 2.3.3.11_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm156_genes.txt

##################################################################
*** note see community_61 (2.3.3.x and others)
# community_126
# 2.3 C2H2 zinc finger factors
# 2.3.3 More than 3 adjacent zinc fingers
# 2.3.3.13 ZNF148-like
cd comm126
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.3.13_mammalia_dbd_fasta.fasta
cat 2.3.3.13_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm126_genes.txt

##################################################################
*** note see community_61 (2.3.3.x and others)
# community_75
# 2.3 C2H2 zinc finger factors
# 2.3.3 More than 3 adjacent zinc fingers
# 2.3.3.16 ZNF238-like
cd comm75
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.3.16_mammalia_dbd_fasta.fasta
cat 2.3.3.16_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm75_genes.txt

##################################################################
*** note see community_61 (2.3.3.x and others)
# community_6
# 2.3 C2H2 zinc finger factors
# 2.3.3 More than 3 adjacent zinc fingers
# 2.3.3.22 BCL6
cd comm6
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.3.22_mammalia_prot_fasta.fasta
cat 2.3.3.22_mammalia_prot_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm6_genes.txt

##################################################################
*** note see community_61 (2.3.3.x and others)
# community_20
# 2.3 C2H2 zinc finger factors
# 2.3.3 More than 3 adjacent zinc fingers
# 2.3.3.50 CTCF-like
# Ctcf (convert lower case)
cd comm20
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.3.50_mammalia_prot_fasta.fasta
cat 2.3.3.50_mammalia_prot_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' |\
sort -u > comm20_genes.txt

##################################################################
*** note see community_61 (2.3.3.x and others)
# community_259
# community_121
# 2.3 C2H2 zinc finger factors
# 2.3.3 More than 3 adjacent zinc fingers
# 2.3.3.65 ZFX-ZFY
cd comm259
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.3.65_mammalia_prot_fasta.fasta
cat 2.3.3.65_mammalia_prot_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm259_genes.txt

##################################################################
*** note see community_61 (2.3.3.x and others)
# community_141
# community_178
# 2.3 C2H2 zinc finger factors
# 2.3.3 More than 3 adjacent zinc fingers
# 2.3.3.80 ZNF99-like
cd comm141
echo Zfp161 > comm141_genes.txt

##################################################################
*** note see community_61 (2.3.3.x and others)
# community_302
# 2.3 C2H2 zinc finger factors
# 2.3.4 Factors with multiple dispersed zinc fingers
# 2.3.4.0.68 ZNF519
cd comm302
echo ZNF519 > comm302_genes.txt

##################################################################
*** note see community_61 (2.3.3.x and others)
# community_311
# 2.3 C2H2 zinc finger factors
# 2.3.4 Factors with multiple dispersed zinc fingers
# 2.3.4.17 HIC
cd comm311
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.4.17_mammalia_prot_fasta.fasta
cat 2.3.4.17_mammalia_prot_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm311_genes.txt

##################################################################
# community_376
# 3.1 Homeo domain factors
# 3.1.4 TALE-type HD
# 3.1.4.2 MEIS
cd comm376
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/3.1.4.2_mammalia_dbd_fasta.fasta
cat 3.1.4.2_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm376_genes.txt

##################################################################
# community_36
# 3.1 Homeo domain factors
# 3.1.4.4 PBX
cd comm36
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/3.1.4.4_mammalia_dbd_fasta.fasta
cat 3.1.4.4_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm36_genes.txt

##################################################################
# community_265
# 3.2 Paired box factors
# 3.2.2 PD
# 3.2.2.2 PAX2-like
cd comm265
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/3.2.2.2_mammalia_dbd_fasta.fasta
cat 3.2.2.2_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm265_genes.txt

##################################################################
# community_45 *** note dimer, comm267
# 3.3 Fork head / winged helix factors
# 3.3.1 FOX
cd comm45
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/3.3.1_mammalia_dbd_fasta.fasta
cat 3.3.1_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm45_genes.txt

##################################################################
# community_13
# community_235
# 3.3 Fork head / winged helix factors
# 3.3.2 E2F
cd comm13
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/3.3.2_mammalia_dbd_fasta.fasta
cat 3.3.2_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' |\
awk -F"(" '{print $1}' | sort -u > comm13_genes.txt

##################################################################
# community_74
# community_131
# community_153
# community_217
# 3.5 Tryptophan cluster factors
# 3.5.2 Ets-related
cd comm74
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/3.5.2_mammalia_dbd_fasta.fasta
cat 3.5.2_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' |\
sort -u > comm74_genes.txt

##################################################################
# community_2
# 3.6 TEA domain factors
# 3.6.1 TEF1-related
# Tead1-4 (convert lower case)
cd comm2
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/3.5.2.5_mammalia_dbd_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/3.6.1_mammalia_dbd_fasta.fasta
cat 3.6.1_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm2_genes.txt

##################################################################
# community_39
# 4.2 Heteromeric CCAAT-binding factors
# 4.2.1 Heteromeric CCAAT-binding factors
# Add CEBPZ
cd comm39
echo CEBPZ > CEBPZ.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/4.2.1_mammalia_prot_fasta.fasta
cat 4.2.1_mammalia_prot_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' |\
sort -u > comm39_genes
cat CEBPZ.fasta comm39_genes > comm39_genes.txt
rm comm39_genes

##################################################################
# community_114
# 5.1 MADS box factors
# 5.1.2 Responders to external signals
cd comm114
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/5.1.2_mammalia_dbd_fasta.fasta
cat 5.1.2_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' |\
sort -u > comm114_genes.txt
 
##################################################################
# community_11
# community_143
# 6.2 STAT domain factors
# 6.2.1 STAT
# Stat1-4,6 and Stat5a,b (convert lower case)
cd comm11
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/6.2.1_mammalia_dbd_fasta.fasta
cat 6.2.1_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"_" '{print $1}' > comm11_genes.txt

##################################################################
# community_12
# 6.4 Runt domain factors
# 6.4.1 Runt-related
cd comm12
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/6.4.1_mammalia_dbd_fasta.fasta
cat 6.4.1_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' |\
sort -u > comm12_genes.txt

##################################################################
# community_41
# 6.5 T-Box factors
# 6.5.3 TBX1-related
cd comm41
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/6.5.3_mammalia_dbd_fasta.fasta
cat 6.5.3_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' |\
sort -u > comm41_genes.txt

##################################################################
# community_146
# 6.7 Grainyhead
# 6.7.2 CP2-related
cd comm146
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/6.7.2_mammalia_prot_fasta.fasta
cat 6.7.2_mammalia_prot_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' |\
sort -u > comm146_genes.txt

##################################################################
# community_1
# community_273
# community_369
# 7.1 SMAD/NF-1 DNA-binding domain factors
# 7.1.1 SMAD
cd comm1
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/7.1.1_mammalia_dbd_fasta.fasta
cat 7.1.1_mammalia_dbd_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' |\
sort -u > comm1_genes.txt


##################################################################
# community_5
# community_94
# community_122
# community_137
# community_378
# 2.3 C2H2 zinc finger factors
# 2.3.3 More than 3 adjacent zinc fingers
# 2.3.3.52 ZNF322-like
# 3.1 Homeo domain factors
# 3.1.2.21 TLX
# 7.1 SMAD/NF-1 DNA-binding domain factors
# 7.1.2 NF-1
cd comm5
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/7.1.2_mammalia_dbd_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.3.52_mammalia_prot_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/3.1.2.21_mammalia_dbd_fasta.fasta
cat *.fasta > C2H2.fasta
cat C2H2.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' |\
sort -u > comm5_genes.txt

##################################################################
# community_25
# 2.3 C2H2 zinc finger factors
# 2.3.4 Factors with multiple dispersed zinc fingers
# 2.3.4.14 EVI1-like
# 6.1 Rel homology region (RHR) factors
# 6.1.5 Early B-Cell Factor-related factors
cd comm25
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.4.14_mammalia_prot_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/6.1.5_mammalia_dbd_fasta.fasta
cat *.fasta > EVI1.fasta
cat EVI1.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' |\
sort -u > comm25_genes.txt


##################################################################
##
## undefined classes
##
##################################################################


##################################################################
# community_34
# 0 Yet undefined DNA-binding domains
# Nrf1 (convert lower case)
cd comm34
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/0.0.6_mammalia_prot_fasta.fasta
cat 0.0.6_mammalia_prot_fasta.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm34_genes.txt


##################################################################
# community_136
# no DNA binding domain found
cd comm136
echo EP300 > comm136_genes.txt

##################################################################
# no homolog
# community_164


##################################################################
##
## dimer and redundant
##
##################################################################

##################################################################
# community_267 ***dimer
# 3.2 Paired box factors
# 3.2.1 PD+HD
# 3.3 Fork head / winged helix factors
# 3.3.1 FOX
cd comm267
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/3.2.1_mammalia_dbd_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/3.3.1.15_mammalia_dbd_fasta.fasta
cat *.fasta > dimer.fasta
cat dimer.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm267_genes.txt

##################################################################
*** note this group has factor above from 2.3.3.x
# community_61
# community_83
# 2.3 C2H2 zinc finger factors
# 2.3.1 Three-zinc finger Krüppel-related factors
# 2.3.1.1 Sp1-like
# 2.3.1.3 EGR
# 2.3.3 More than 3 adjacent zinc finger factors
# 2.3.4 Factors with multiple dispersed zink fingers
cd comm61
printf "Znf263\nRreb1\n" > ZnfZNF148263.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.1.1_mammalia_dbd_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.1.3_mammalia_dbd_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.3.22_mammalia_prot_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.3.0_mammalia_prot_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.3.13_mammalia_dbd_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.3.8_mammalia_dbd_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.4.8_mammalia_prot_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.4.3_mammalia_prot_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.4.0_mammalia_prot_fasta.fasta
cat *.fasta* > C2H2.fasta
cat C2H2.fasta | grep Mus_musculus |\
awk -F"Mus_musculus_" '{print $2}' | awk -F"-" '{print $1}' > comm61_genes.txt



##################################################################
##
## isolate all lists from TFclass
##
##################################################################

dir=/media/wa3j/Seagate2/Documents/PRO/adipogenesis/July2018/atac_time_deg/TFlist/comdat
cd ${dir}

dirs=$(ls -d *comm*)
for dir in ${dirs}
do
cd ${dir}
cp *genes.txt* ..
cd ..
done



#####################################3
## see communityTFs.R






