# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
###########
# finland #
############

# text files downloaded to /path/to/project

md5sum -c finland.txt.gz.md5
# finland.txt.gz.md5: OK

zcat /path/to/project | head -2
# Sample ID	   SNP Name	Chr	  Position	REF	  ALT	   GT	    X Raw	  Y Raw	     X      Y  GC Score
# sample_id	rs3131972	  1	    752721 	  A	   G	  0/1	     1977	   1480	 0.567	0.635	0.8260000

# check GT codes
zcat finland_illugwas.txt.gz | awk -F '\t' '{print $7}'| sort | uniq -c | sort -nr


# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/finland/
  
############################################################################## 
# make a temporary sorted fcr file:

# to keep consistency in naming
# mv /path/to/project \
# /path/to/project

# rename file to reflect gwas

mv /path/to/project \
/path/to/project

# farm4:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh finland long illugwas
# Job <366165> is submitted to default queue <long>.


# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/sort_fcr* /path/to/user/scripts/IIBDGC/

zcat /path/to/ibdgwas/IIBDGC/pre_imputation/raw/finland/finland_illugwas_sorted.txt.gz | wc -l
# 149323105

zcat /path/to/project | wc -l
# 149323105

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

# farm5:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_map_from_fcr.sh finland normal illugwas
#job <401525> is submitted to queue <normal>.

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
bash ${path_gwas}scripts/create_fam_from_fcr.sh finland normal illugwas
#Job <401433> is submitted to queue <normal>.

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

# farm5:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_lgen_from_fcr.sh finland normal illugwas
# Job <401573> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/
  



wc -l /path/to/ibdgwas/IIBDGC/pre_imputation/raw/finland/*.fam
# 612 /path/to/ibdgwas/IIBDGC/pre_imputation/raw/finland/finland_illugwas.fam

wc -l /path/to/ibdgwas/IIBDGC/pre_imputation/raw/finland/*.map
# 243992 /path/to/ibdgwas/IIBDGC/pre_imputation/raw/finland/finland_illugwas.map

zcat /path/to/project | awk 'END { print NR }'
# 149323105

# 149323105/243992
# [1] 612 OK

############################
# remove sorted file:
# rm /path/to/ibdgwas/IIBDGC/pre_imputation/raw/finland/finland_illugwas_sorted.txt.gz

##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/finland_illugwas/
mkdir ${path_gwas}pre_imputation/QC/finland_illugwas/logs/
  
MEM=1000
path_gwas=/path/to/ibdgwas/IIBDGC/
bsub -J"plink0" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/finland_illugwas/logs/stderr_plink_0_finland_illugwas \
-o ${path_gwas}pre_imputation/QC/finland_illugwas/logs/stdout_plink_0_finland_illugwas \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/finland/finland_illugwas \
--missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/finland_illugwas/finland_illugwas_hg19"

# 243992 variants and 612 people pass filters and QC.
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

bim<-read.table(paste(path,"pre_imputation/QC/finland_illugwas/finland_illugwas_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#       0     A     C     G     T
# 0     8  1297  1118  1089  1236
# A     0     0 11732 50840    78
# C     0 10881     0   119 45683
# G     0 46019   112     0 10729
# T     0    70 51124 11857     0

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 4748    6

table(indels$V1)
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
# 445 408 345 275 240 223 173 215 150 253 236 216 171 165 121 156 154 168  83 114 
# 21  22  23 
# 54  68 315 

# X    X chromosome                    -> 23
# Y    Y chromosome                    -> 24
# XY   Pseudo-autosomal region of X    -> 25
# MT   Mitochondrial                   -> 26

mt<-bim[which(bim$V1 %in% c(0,26)),]
table(mt$V1)
# 0  

chry<-bim[which(bim$V1 %in% c(24)),]
table(chry$V1)
# 0

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])
#  

dim(all_remove)
#[1] 4748      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/finland_illugwas/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/finland_illugwas/finland_illugwas_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# finland 
# 612


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#             Crohn's Disease Indeterminate Ulcerative Colitis
# finland 157              99            15                341

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#             Affected Unaffected
# finland  12      457        143

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#           0
# finland 264

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                        Affected Unaffected
#                     12        2        143
# Crohn's Disease      0       99          0
# Indeterminate        0       15          0
# Ulcerative Colitis   0      341          0

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#             Female Male
# finland   6    240  366


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
#   2  -9 
# 457 155 



################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# 1   2  -9 
# 366 243   3 

fam_test<-read.table(paste(path,"pre_imputation/QC/finland_illugwas/finland_illugwas_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 612

write.table(fam,paste(path,"pre_imputation/QC/finland_illugwas/finland_illugwas_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 155   2

write.table(samples_remove,paste(path,"pre_imputation/QC/finland_illugwas/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")


######

mkdir ${path_gwas}pre_imputation/QC/finland_illugwas/
mkdir ${path_gwas}pre_imputation/QC/finland_illugwas/logs/
  
path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

bsub -J"pl1" -M"$MEM" -n4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G "ibdgwas" -q normal \
-e ${path_gwas}pre_imputation/QC/finland_illugwas/logs/stderr_plink_finland_illugwas_1 \
-o ${path_gwas}pre_imputation/QC/finland_illugwas/logs/stdout_plink_finland_illugwas_1 \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/finland_illugwas/finland_illugwas_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/finland_illugwas/finland_illugwas_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/finland_illugwas/finland_illugwas_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/finland_illugwas/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/finland_illugwas/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/finland_illugwas/finland_illugwas_hg19_noind"
# Job <426053> is submitted to queue <normal>.

less ${path_gwas}pre_imputation/QC/finland_illugwas/logs/stdout_plink_finland_illugwas_1
# 239244 variants and 457 people pass filters and QC.
# Among remaining phenotypes, 457 are cases and 0 are controls.




#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37


