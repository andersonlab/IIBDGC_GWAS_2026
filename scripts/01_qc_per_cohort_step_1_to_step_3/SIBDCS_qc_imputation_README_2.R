# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
#########################
# kiel, austria, sibdcs #
#########################

# text files downloaded to /path/to/project

# FIVE SEPARATED FILES
cd /path/to/project
  
md5sum -c sibdcs_gsa.txt.gz.md5
# sibdcs_gsa.txt.gz: OK

zcat /path/to/project | head -2

# create folders for the study:

mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/sibdcs/
  
############################################################################## 
# make a temporary sorted fcr file:

bash ${path_gwas}scripts/sort_fcr.sh sibdcs basement gsa
# Job <8639813> is submitted to queue <basement>.

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
  bash ${path_gwas}scripts/create_map_from_fcr.sh sibdcs normal gsa
# Job <8882207> is submitted to queue <long>.

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
  bash ${path_gwas}scripts/create_fam_from_fcr.sh sibdcs normal gsa
# Job <8657082> is submitted to queue <normal>.

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
  bash ${path_gwas}scripts/create_lgen_from_fcr.sh sibdcs normal gsa
# Job <8657083> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/

wc -l sibdcs_gsa.fam
# 2322

wc -l sibdcs_gsa.map
# 712189

zcat /path/to/project | awk 'END { print NR }'
# 1653702859

1653702859/2322
# [1] 712189 OK


##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/sibdcs_gsa/
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/sibdcs/sibdcs_gsa \
--missing-genotype - \
--output-missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/sibdcs_gsa/sibdcs_gsa_hg19
# 712189 variants and 2322 people pass filters and QC.


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

bim<-read.table(paste(path,"pre_imputation/QC/sibdcs_gsa/sibdcs_gsa_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        0      A      C      D      G      I      T
# 0  18748  27304  13111   1099  48138   3475   1096
# A      0      0  59905      0 277696      0   1198
# C      0  50650      0      0   1816      0      0
# D      0      0      0      0      0    248      0
# G      0 204940   1692      0      0      0      0
# I      0      0      0     99      0      0      0
# T      0    974      0      0      0      0      0

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
# [1] 113318      6

table(indels$V1)
# 0     1     2     3     4     5     6     7     8     9    10    11    12 
# 6  9298 10023  6757  4499  5319  5111  5563  4020  3880  4223  6632  5032 
# 13    14    15    16    17    18    19    20    21    22    23    24    25 
# 4775  3032  4051  4123  6272  1936  5505  2009  1036  2100  4478  2761   331 
# 26 
# 546 

# X    X chromosome                    -> 23
# Y    Y chromosome                    -> 24
# XY   Pseudo-autosomal region of X    -> 25
# MT   Mitochondrial                   -> 26

mt<-bim[which(bim$V1 %in% c(0,26)),]
table(mt$V1)
#  0   26 
#  8 1180 

chry<-bim[which(bim$V1 %in% c(24)),]
table(chry$V1)
# 24 
# 4974

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])
# 0     1     2     3     4     5     6     7     8     9    10    11    12 
# 8  9298 10023  6757  4499  5319  5111  5563  4020  3880  4223  6632  5032 
# 13    14    15    16    17    18    19    20    21    22    23    24    25 
# 4775  3032  4051  4123  6272  1936  5505  2009  1036  2100  4478  2761   331 
# 26 
# 1180 

dim(all_remove)
#[1] 113954     6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/sibdcs_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

##########################################################
# KEEP ONLY THE CASES AND CONTROLS AS DONE IN PREVIOUS QC


pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/sibdcs_gsa/sibdcs_gsa_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# sibdcs 
# 2322


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#             Crohn's Disease Indeterminate Ulcerative Colitis Unknown
# sibdcs    1            1296            87                918      20


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#             Affected
# sibdcs    1     2321

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#           0
# sibdcs 2321

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                         Affected
#                       1        0
# Crohn's Disease       0     1296
# Indeterminate         0       87
# Ulcerative Colitis    0      918
# Unknown               0       20

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#             Female Male Unknown
# sibdcs    1   1123 1197       1


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
# 2   -9 
# 2321    1 

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# 1    2   -9 
# 1197 1123    2 

fam_test<-read.table(paste(path,"pre_imputation/QC/sibdcs_gsa/sibdcs_gsa_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 2322

write.table(fam,paste(path,"pre_imputation/QC/sibdcs_gsa/sibdcs_gsa_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 1  2

write.table(samples_remove,paste(path,"pre_imputation/QC/sibdcs_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")


######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/sibdcs_gsa/sibdcs_gsa_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/sibdcs_gsa/sibdcs_gsa_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/sibdcs_gsa/sibdcs_gsa_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/sibdcs_gsa/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/sibdcs_gsa/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/sibdcs_gsa/sibdcs_gsa_hg19_noind
# 598235 variants and 2321 people pass filters and QC.
# Among remaining phenotypes, 2321 are cases and 0 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37

