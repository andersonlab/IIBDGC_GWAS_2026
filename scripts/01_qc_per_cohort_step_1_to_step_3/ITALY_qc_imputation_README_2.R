# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
#########
# italy #
#########

# text files downloaded to /path/to/project

# FIVE SEPARATED FILES
cd /path/to/project
  
md5sum -c italy_gsa.txt.gz.md5
# italy_gsa.txt.gz: OK


zcat /path/to/project | head -2
# Sample ID	Sample Group	SNP Name	Chr	Position	Allele1 - Top	Allele2 - Top	X	Y'
# sample_id		1:103380393	1	103380393	G	G	0.044	

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/italy/


############################################################################## 
# make a temporary sorted fcr file:

# farm4:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh italy basement gsa
# Job <8600207> is submitted to queue <basement>.


# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/sort_fcr* /path/to/user/scripts/IIBDGC/

##############################################################################  
# generate map, fam and lgen files from this:

#############
# MAP FILES #
#############

# MAP file describes a single marker and must contain exactly 4 columns:
# chromosome (1-22, X, Y or 0 if unplaced)
# rs# or snp identifier
# Genetic distance (morgans)
# Base-pair position (bp units)

# farm4:

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_map_from_fcr.sh italy long gsa
#Job <8626745> is submitted to queue <long>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_map_from_fcr* /path/to/user/scripts/IIBDGC/

#############
# FAM FILES #
#############

# FAM (6 first columns of PED) PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:
# Family ID
# Individual ID
# Paternal ID
# Maternal ID
# Sex (1=male; 2=female; other=unknown)
# Phenotype

# farm4:

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_fam_from_fcr.sh italy normal gsa
# Job <8626747> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_fam_from_fcr* /path/to/user/scripts/IIBDGC/

##############
# LGEN FILES #
##############

# LGEN file, test.lgen
# family ID
# individual ID
# snp ID
# allele 1 of this genotype
# allele 2 of this genotype

# farm4:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_lgen_from_fcr.sh italy normal gsa
# Job <8626748> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/


wc -l italy_gsa.fam 
# 1016 italy_gsa.fam # OK

wc -l italy_gsa.map 
# 712189 italy_gsa.map

zcat /path/to/project | awk 'END { print NR }'
# 723584025
# 723584025/712189
# [1] 1016 # OK


##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/italy_gsa/
mkdir ${path_gwas}pre_imputation/QC/italy_gsa/logs/
  
MEM=1000
path_gwas=/path/to/ibdgwas/IIBDGC/
bsub -J"plink0" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/italy_gsa/logs/stderr_plink_0_italy_gsa \
-o ${path_gwas}pre_imputation/QC/italy_gsa/logs/stdout_plink_0_italy_gsa \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/italy/italy_gsa \
--missing-genotype - \
--make-bed --out ${path_gwas}pre_imputation/QC/italy_gsa/italy_gsa_hg19"
# Job <912882> is submitted to queue <normal>.

# 712189 variants loaded from .bim file.
# 1016 people (0 males, 0 females, 1016 ambiguous) loaded from .fam.
# 712189 variants and 1016 people pass filters and QC.


#######################################################################################################################################

#####################
# 1.- REMOVE INDELS #
#####################


### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

bim<-read.table(paste(path,"pre_imputation/QC/italy_gsa/italy_gsa_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        -      A      C      D      G      I      T
# -  18780  34311  13615    950  52669   3114   1056
# A      0      0  59414      0 273288      0   1240
# C      0  49429      0      0   1779      0      0
# D      0      0      0      0      0    608      0
# G      0 199104   1639      0      0      0      0
# I      0      0      0    249      0      0      0
# T      0    944      0      0      0      0      0

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 125352      6

table(indels$V1)
# 0     1     2     3     4     5     6     7     8     9    10    11    12 
# 6 10358 10957  7604  5079  5956  5895  5990  4610  4484  4837  7322  5705 
# 13    14    15    16    17    18    19    20    21    22    23    24    25 
# 4907  3457  4467  4622  6752  2260  5971  2363  1205  2324  4327  2840   355 
# 26 
# 699 


# X    X chromosome                    -> 23
# Y    Y chromosome                    -> 24
# XY   Pseudo-autosomal region of X    -> 25
# MT   Mitochondrial                   -> 26

mt<-bim[which(bim$V1 %in% c(0,26)),]
table(mt$V1)
# 0   26 
# 8 1180  

chry<-bim[which(bim$V1 %in% c(24)),]
table(chry$V1)
# 24 
# 4974

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])
# 0     1     2     3     4     5     6     7     8     9    10    11    12 
# 8 10358 10957  7604  5079  5956  5895  5990  4610  4484  4837  7322  5705 
# 13    14    15    16    17    18    19    20    21    22    23    24    25 
# 4907  3457  4467  4622  6752  2260  5971  2363  1205  2324  4327  2840   355 
# 26 
# 1180 

dim(all_remove)
#[1] 125835      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/italy_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

##########################################################
# KEEP ONLY THE CASES AND CONTROLS AS DONE IN PREVIOUS QC


pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/italy_gsa/italy_gsa_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# italy 
# 1016 

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#           Crohn's Disease Ulcerative Colitis
# italy 396             310                310

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#       Affected Unaffected
# italy      620        396

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#         0   1
# italy 620 396

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                    Affected Unaffected
#                           0        396
# Crohn's Disease         310          0
# Ulcerative Colitis      310          0

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#       Female Male
# italy    465  551


################################################
# Affection status, by default, should be coded:
# -9 missing 
# 0 missing
# 1 unaffected
# 2 affected

fam[,"V6"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$control=="1")]),"V6"]<-"1"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$affection=="Affected")]),"V6"]<-"2"

table(fam$V6)
# 1   2 
# 396 620   

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# 1   2 
# 551 465   

fam_test<-read.table(paste(path,"pre_imputation/QC/italy_gsa/italy_gsa_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 1016

write.table(fam,paste(path,"pre_imputation/QC/italy_gsa/italy_gsa_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 0   2

write.table(samples_remove,paste(path,"pre_imputation/QC/italy_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")


######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/italy_gsa/italy_gsa_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/italy_gsa/italy_gsa_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/italy_gsa/italy_gsa_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/italy_gsa/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/italy_gsa/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/italy_gsa/italy_gsa_hg19_noind
# 586354 variants and 1016 people pass filters and QC.
# Among remaining phenotypes, 620 are cases and 396 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37


