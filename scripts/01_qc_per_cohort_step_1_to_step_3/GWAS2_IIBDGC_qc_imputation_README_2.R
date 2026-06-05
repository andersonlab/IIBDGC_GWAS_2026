# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# UK-IIBGBD:

# created new folder to store the files, with similar structure as previous analysis:

/path/to/ibdgwas/IIBDGC/
  
path_gwas=/path/to/ibdgwas/IIBDGC/
    
#########
# gwas2 #
#########

# For previous steps see: /path/to/ibdgwas/pre_imputation/raw/wtccc-hrc/liftOver/sample_id/liftOver_steps.txt

# start with files that have been updated to build 37:

/path/to/software/plink_linux_x86_64_20181202/./plink --file /path/to/ibdgwas/pre_imputation/raw/wtccc-hrc/liftOver/sample_id/wtccc2_hg19 \
--make-bed --out ${path_gwas}pre_imputation/raw/gwas2/wtccc2_hg19

/path/to/software/plink_linux_x86_64_20181202/./plink \
--file /path/to/ibdgwas/pre_imputation/raw/wtccc-hrc/liftOver/sample_id/wtccc2_hg19 \
--make-bed --out ${path_gwas}pre_imputation/QC/gwas2/gwas2_hg19
#801371 variants and 7778 people pass filters and QC

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/gwas2/gwas2_hg19 \
--freq --out ${path_gwas}pre_imputation/QC/gwas2/gwas2_hg19

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

bim<-read.table(paste(path,"pre_imputation/QC/gwas2/gwas2_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        A      C      G      T
# 0    355    599    647    314
# A      0  37156 145077  24859
# C  33706      0  36649 124216
# G 124454  36067      0  24097
# T  26183 156340  30652      0

table(bim[which(bim$V5==0),"V1"])
# 23  24  26 
# 968 694 253

# all maf=0, ok to exclude
frq<-read.table(paste(path,"pre_imputation/QC/gwas2/gwas2_hg19.frq",sep=""),head=T)
summary(frq$MAF[which(frq$A1==0)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       0       0       0       0       0 
length(frq$MAF[which(frq$A1==0)])
# [1] 1915

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 12635      6

table(indels$V1)

# X    X chromosome                    -> 23
# Y    Y chromosome                    -> 24
# XY   Pseudo-autosomal region of X    -> 25
# MT   Mitochondrial                   -> 26

mt<-bim[which(bim$V1 %in% c(0,26)),]
table(mt$V1)
# 26 
# 395 

chry<-bim[which(bim$V1 %in% c(24)),]
table(chry$V1)
# 24 
# 901 

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])

dim(all_remove)
#[1] 2107      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/gwas2/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

##########################################################
# KEEP ONLY THE CASES AND CONTROLS AS DONE IN PREVIOUS QC

## R

fam<-read.table(paste(path,"pre_imputation/QC/gwas2/gwas2_hg19.fam",sep=""),head=F)

# SEX:
table(fam$V5,useNA="ifany")
# 1    2 
# 3908 3870 

# pheno:
table(fam$V6,useNA="ifany")
# 1    2 
# 5417 2361

############################################
# double check assigned pheno from phenoDB, plus sex. There is an issue reported by Loukas see "/path/to/ibdgwas/pre_imputation/qc/gwas2/qc_steps_gwas2.txt":

# see script "manifest_approach.R" to see the matches

# fam_file_ids<-read.table(paste(path,"/manifests/match_sanger_ids_supplier_ids_fam_manifest.tsv",sep=""),sep="\t",head=T)
# for (i in 1:nrow(fam_file_ids)){
#   fam$id[which(fam$id %in% fam_file_ids$sanger_sample_name_fam_file[i])]<-as.character(fam_file_ids$sanger_sample_name_manifest[i])
# }


# AFTER MEETING WITH CARL AND ALEX, KEEP PHENO IN FAM, THEY RELY MORE IN KATIE'S FILES THAN IN PHENODB

sample<-read.table("/path/to/ibdgwas/post_imputation/GWAS2/GWAS2.sample",head=T)
sample<-sample[-1,] # remove subheader
dim(sample)
#[1] 4697   15
dim(fam)
#[1] 7778    6

dim(fam[which(fam$V1 %in% sample$ID_1),])
#[1] 4697    6

fam_subset<-fam[which(fam$V1 %in% sample$ID_1),]
fam_subset<-fam_subset[match(sample$ID_1,fam_subset$V1),]
table(as.character(fam_subset$V1)==as.character(sample$ID_1))
#TRUE 
#4697 

table(as.character(fam_subset$V6),as.character(sample$bin1))
#      0    1
# 1 2776    0
# 2    0 1921


table(as.character(fam_subset$V5),as.character(sample$sex))
#      1    2
# 1 2380    0
# 2    0 2317


# same phenotypes and sex as in postQC sample file

################################################
# Affection status, by default, should be coded:
# -9 missing 
# 0 missing
# 1 unaffected
# 2 affected

table(fam$V6)
# 1    2 
# 5417 2361  

################################################
# Sex (1=male; 2=female; other=unknown)

table(fam$V5)
# 1    2 
# 3908 3870 

fam_test<-read.table(paste(path,"/pre_imputation/QC/gwas2/gwas2_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 7778

write.table(fam,paste(path,"pre_imputation/QC/gwas2/gwas2_hg19_tmp_edited.fam",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 0   2

write.table(samples_remove,paste(path,"pre_imputation/QC/gwas2/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/gwas2/gwas2_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/gwas2/gwas2_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/gwas2/gwas2_hg19_tmp_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/gwas2/list_indel_var_exclude \
--make-bed --out ${path_gwas}pre_imputation/QC/gwas2/gwas2_hg19_noind
# 7778 people (3908 males, 3870 females) loaded from .fam.
# 799264 variants and 7778 people pass filters and QC.
# Among remaining phenotypes, 2361 are cases and 5417 are controls.


#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37


