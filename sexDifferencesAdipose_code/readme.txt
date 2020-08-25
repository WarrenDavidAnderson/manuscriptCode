
############################################################################
# Sex differences in human adipose tissue gene expression and genetic regulation involve adipogenesis
# Genome Research, 2020
# Information on the files in this repository
# Contact: warrena@virginia.edu, mete@virginia.edu
############################################################################


############################################################################
# Data acquisition and basic processing/analysis
# DEG and eQTL analysis
############################################################################

###################################################
## Fig 1 - GEO data

# Download and initial processing of GEO data (deCODE and AAGMEx)
Fig1.GEOdata.pdf
Fig1.decode.R
Fig1.aagmex.R

# Quality control analysis of GEO data (deCODE and AAGMEx)
Fig1.GEOqc.pdf
Fig1.decodeQC.R
Fig1.aagmexQC.R

###################################################
## Fig 1 - Annotation of autosomal genes
Fig1.autosomal.pdf

###################################################
## Fig 1 - GTEx data acquisition, normalization, qc
Fig1.GTEx.pdf

###################################################
## Fig S3 - GTEx covariate correction and PEER analysis
FigS3.GTEx.pdf

###################################################
## Fig 1 - Differential expression analysis
Fig1_deg.R
Fig1_deg.pdf

###################################################
## Fig 1 - DEG analysis for obesity
Fig1_obesity.R
Fig1_obesity.pdf

###################################################
## Fig 3, Fig S4 - GTEx PEER/eQTL analysis
gtex_peer.pdf



############################################################################
# Downstream analysis and figure generation
############################################################################

###################################################
## Fig 1
Fig1_human.R 	# integration of human data
Fig1_hmdp.R		# HMDP DEG analysis

###################################################
## Fig 2
BART.sh			# data formatting for BART analysis
GSEA.R			# function for implementing GSEA
newGSEAplots.R 	# function for plotting GSEA results
Fig2.R			# Fig 2 analysis and plots

###################################################
## Fig 3
getLDsnps.sh 	# code for identifying LD SNPS
heritability.sh # code for the heritability analysis
crossmap.sh 	# code for converting between genome builds
gtex_eqtl.R 	# GTEx eQTL data processing (for Figs 4,5)
coloc.R 		# Final GWAS overlap and coloc analysis
Fig3.sh			# Parsing GWAS catalog results
Fig3a.R			# Fig 3a analysis and plots
Fig3b.R			# Fig 3b analysis and plots
chromHMM.sh		# Fig 3d analysis and plots
HiCanalysis.sh	# promoter-capture HiC analysis

###################################################
## Fig 4
Fig4_atac.sh 	# ATACseq analysis and command line code
Fig4ae.R 		# Fig 4 panels a and e
Fig4f.R 		# Fig 4 panel f

###################################################
## Fig 5
Fig5.R 			# Fig 5 analysis and figure panels

###################################################
## Fig 6
Fig1_human.R 	# expression plots
Fig3a.R			# eQTL plots

###################################################
## Fig 7
Fig7.R 			# Fig 7 analysis and figure panels

###################################################
## Fig S1
## see Fig 1 vignettes (Data acquisition and basic processing/analysis)
FigS1.R			# Fig S1D analysis

###################################################
## Fig S2
FigS2.R 		# Fig S2 simulations

###################################################
## Fig S5
FigS5.R			# Fig S5 code

###################################################
## Fig S6
## See Fig 1 code

###################################################
## Fig S7
gtexAAcounts.R 	# parse the GTEx African American data
FigS7.R 		# Fig S7 analysis and plot

###################################################
## Other
gwas_coloc.sh	# coloc analysis
coloc.R			# coloc analysis





