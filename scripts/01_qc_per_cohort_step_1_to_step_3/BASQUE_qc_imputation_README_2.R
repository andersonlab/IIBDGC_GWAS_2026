# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
##########
# basque #
##########

# text files downloaded to /path/to/project

# FIVE SEPARATED FILES
cd /path/to/project
  
md5sum -c basque.txt.gz.md5
# basque.txt.gz: OK


zcat /path/to/project | head -2
#  Sample ID	Sample Group	   SNP Name	  Chr	  Position	  Allele1 - Top	Allele2 - Top	X Raw	Y Raw	X	Y	GC Score	GT Score
#  sample_id		            1:103380393	    1	 103380393	              G	            G	594	5329	0.040	0.804	0.8999	0.8634


# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/basque/
  
mv /path/to/project \
/path/to/project

############################################################################## 
# make a temporary sorted fcr file:

# farm4:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr_vs1.sh basque basement gsa
# Job <492973> is submitted to queue <basement>.


# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/sort_fcr* /path/to/user/scripts/IIBDGC/

zcat /path/to/project | wc -l
# 1084663848
zcat /path/to/ibdgwas/IIBDGC/pre_imputation/raw/basque/basque_gsa_sorted.txt.gz | wc -l
# 1084663848

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
bash ${path_gwas}scripts/create_map_from_fcr.sh basque long gsa
#Job <678700> is submitted to queue <long>.

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
bash ${path_gwas}scripts/create_fam_from_fcr.sh basque normal gsa
# Job <678702> is submitted to queue <normal>.

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
bash ${path_gwas}scripts/create_lgen_from_fcr.sh basque normal gsa
# Job <678703> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/


wc -l basque_gsa.fam 
# 1523  basque_gsa.fam # OK

wc -l basque_gsa.map 
# 712189 basque_gsa.map

zcat /path/to/project | awk 'END { print NR }'
# 1084663848
# 1084663848/712189
# [1] 1523 # OK

##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/basque_gsa/

/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/basque/basque_gsa \
--missing-genotype - \
--make-bed --out ${path_gwas}pre_imputation/QC/basque_gsa/basque_gsa_hg19

# 712189 variants and 1523 people pass filters and QC.
# Note: No phenotypes present.

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

bim<-read.table(paste(path,"pre_imputation/QC/basque_gsa/basque_gsa_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        -      A      C      D      G      I      T
# -  18748  35638  14768    992  57147   3232   1110
# A      0      0  58305      0 268872      0   1177
# C      0  49129      0      0   1721      0      0
# D      0      0      0      0      0    489      0
# G      0 198106   1610      0      0      0      0
# I      0      0      0    208      0      0      0
# T      0    937      0      0      0      0      0

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 132332      6

table(indels$V1)
#


# X    X chromosome                    -> 23
# Y    Y chromosome                    -> 24
# XY   Pseudo-autosomal region of X    -> 25
# MT   Mitochondrial                   -> 26

mt<-bim[which(bim$V1 %in% c(0,26)),]
table(mt$V1)
#    0   26 
# 8 1180 

chry<-bim[which(bim$V1 %in% c(24)),]
table(chry$V1)
# 24 
# 4974

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])
# 0     1     2     3     4     5     6     7     8     9    10    11    12 
# 8 10940 11655  8075  5481  6338  6312  6485  5026  4780  5162  7646  6051 
# 13    14    15    16    17    18    19    20    21    22    23    24    25 
# 5228  3642  4645  4789  6969  2406  6069  2496  1266  2436  4564  2823   357 
# 26 
# 1180 

dim(all_remove)
#[1] 132829      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/basque_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/basque_gsa/basque_gsa_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# basque 
#   1523 


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#            Crohn's Disease Ulcerative Colitis 
# basque 994             301                228

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#            Affected Unaffected
# basque   7      529        987

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#          0   1
# basque 529 987

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                        Affected Unaffected
#                      7        0        987
# Crohn's Disease      0      301          0
# Ulcerative Colitis   0      228          0

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#        Female Male Unknown
# basque    585  933       5


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
# 1   2  -9 
# 987 529   7 

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# 1   2  -9 
# 933 585   5  

fam_test<-read.table(paste(path,"pre_imputation/QC/basque_gsa/basque_gsa_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 1523

write.table(fam,paste(path,"pre_imputation/QC/basque_gsa/basque_gsa_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 7   2

write.table(samples_remove,paste(path,"pre_imputation/QC/basque_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")


######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/basque_gsa/basque_gsa_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/basque_gsa/basque_gsa_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/basque_gsa/basque_gsa_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/basque_gsa/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/basque_gsa/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/basque_gsa/basque_gsa_hg19_noind
# 1523 people (933 males, 585 females, 5 ambiguous) loaded from .fam.
# 579360 variants and 1516 people pass filters and QC.
# Among remaining phenotypes, 529 are cases and 987 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37


