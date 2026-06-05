# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
###########
# austria #
###########

# text files downloaded to /path/to/project

# FIVE SEPARATED FILES
cd /path/to/project

md5sum -c austria_gsa.txt.gz.md5
# austria_gsa.txt.gz: OK

zcat /path/to/project | head -2

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/austria/

  
############################################################################## 
# make a temporary sorted fcr file:

# farm4:
bash ${path_gwas}scripts/sort_fcr.sh austria basement gsa
# Job <8600195> is submitted to queue <basement>.

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
  bash ${path_gwas}scripts/create_map_from_fcr.sh austria long gsa
# Job <8625705> is submitted to queue <long>.

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
  bash ${path_gwas}scripts/create_fam_from_fcr.sh austria normal gsa
# Job <8625706> is submitted to queue <normal>.

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
  bash ${path_gwas}scripts/create_lgen_from_fcr.sh austria normal gsa
# Job <8625708> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/

# (base) user-server:/path/to/ibdgwas/IIBDGC/pre_imputation/raw$ wc -l kiel_*/*.fam
# 1946 kiel_bc/kiel_bc_gsa.fam
# 358 kiel_eze/kiel_eze_gsa.fam
# 1098 kiel_foc/kiel_foc_gsa.fam
# 3582 kiel_hlitm/kiel_hlitm_gsa.fam
# 3523 kiel_ibd/kiel_ibd_gsa.fam
# 10507 total
# 
# (base) user-server:/path/to/ibdgwas/IIBDGC/pre_imputation/raw$ wc -l kiel_*/*.map
# 700078 kiel_bc/kiel_bc_gsa.map
# 700078 kiel_eze/kiel_eze_gsa.map
# 700078 kiel_foc/kiel_foc_gsa.map
# 700078 kiel_hlitm/kiel_hlitm_gsa.map
# 700078 kiel_ibd/kiel_ibd_gsa.map
# 3500390 total
# 

# wc -l austria_gsa.fam
# # 1832 austria_gsa.fam
# 
# wc -l austria_gsa.map
# # 712189 austria_gsa.map
# 
# zcat /path/to/project | awk 'END { print NR }'
# # 1304730249

# 1304730249/712189
# # [1] 1832 OK


##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/austria_gsa/
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/austria/austria_gsa \
--missing-genotype - \
--output-missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/austria_gsa/austria_gsa_hg19
# 712189 variants and 1832 people pass filters and QC.


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

bim<-read.table(paste(path,"pre_imputation/QC/austria_gsa/austria_gsa_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        0      A      C      D      G      I      T
# 0  18748  31640  14595   1105  54826   3484   1151
# A      0      0  58602      0 271136      0   1139
# C      0  49847      0      0   1666      0      0
# D      0      0      0      0      0    239      0
# G      0 201493   1511      0      0      0      0
# I      0      0      0     93      0      0      0
# T      0    914      0      0      0      0      0

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
# [1] 125881      6

table(indels$V1)
# 0     1     2     3     4     5     6     7     8     9    10    11    12 
# 6 10402 11090  7505  5166  6005  5919  6173  4674  4456  4772  7300  5656 
# 13    14    15    16    17    18    19    20    21    22    23    24    25 
# 5117  3439  4431  4518  6717  2184  5972  2274  1193  2290  4853  2801   346 
# 26 
# 622 

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
#

dim(all_remove)
#[1] 126441      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/austria_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

##########################################################
# KEEP ONLY THE CASES AND CONTROLS AS DONE IN PREVIOUS QC


pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/austria_gsa/austria_gsa_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
#austria 
#   1832


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#         Crohn's Disease Indeterminate Ulcerative Colitis Unknown
# austria            1126            35                375     296

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#         Affected
# austria     1832

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#            0    1
# austria 1832    0


table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                    Affected
# Crohn's Disease        1126
# Indeterminate            35
# Ulcerative Colitis      375
# Unknown                 296

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#              Female Male Unknown
# austria    0    943  888       1



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
#    2 
# 1832

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# 1   2  -9 
# 888 943   1 

fam_test<-read.table(paste(path,"pre_imputation/QC/austria_gsa/austria_gsa_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 1832

write.table(fam,paste(path,"pre_imputation/QC/austria_gsa/austria_gsa_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 0  2

write.table(samples_remove,paste(path,"pre_imputation/QC/austria_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")


######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/austria_gsa/austria_gsa_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/austria_gsa/austria_gsa_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/austria_gsa/austria_gsa_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/austria_gsa/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/austria_gsa/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/austria_gsa/austria_gsa_hg19_noind
# 585748 variants and 1832 people pass filters and QC.
# Among remaining phenotypes, 1832 are cases and 0 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37

