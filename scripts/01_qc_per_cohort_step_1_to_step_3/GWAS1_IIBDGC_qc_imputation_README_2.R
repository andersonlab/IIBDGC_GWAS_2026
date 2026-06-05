# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# UK-IIBGBD:

# created new folder to store the files, with similar structure as previous analysis:

cd /path/to/ibdgwas/IIBDGC/
  
path_gwas=/path/to/ibdgwas/IIBDGC/
    
#########
# gwas1 #
#########

# For previous steps see: /path/to/ibdgwas/pre_imputation/raw/wtccc-hrc/liftOver/gwas1/liftOver_steps_gwas1.txt

# start with files that have been updated to build 37:

/path/to/software/plink_linux_x86_64_20181202/./plink --file /path/to/ibdgwas/pre_imputation/raw/wtccc-hrc/liftOver/sample_id/wtccc1_hg19 \
--make-bed --out ${path_gwas}/pre_imputation/raw/gwas1/wtccc1_hg19


/path/to/software/plink_linux_x86_64_20181202/./plink --bfile ${path_gwas}/pre_imputation/raw/gwas1/wtccc1_hg19 \
--keep ${path_gwas}/pre_imputation/raw/gwas1/list_cd_ctrl_samples.txt --make-bed \
--out ${path_gwas}/pre_imputation/QC/gwas1/gwas1_hg19
# 469281 variants and 4684 people pass filters and QC. 

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

bim<-read.table(paste(path,"pre_imputation/QC/gwas1/gwas1_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#       A     C     G     T
# 0  2096  2877  2807  2075
# A     0 19087 84978 14880
# C 16463     0 22502 71301
# G 71627 22552     0 16562
# T 15147 85316 19011     0

table(bim[which(bim$V5==0),"V1"])
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
# 937 902 613 757 588 469 465 577 347 498 551 497 407 288 280 359 223 359 131 235 
# 21  22  23 
# 149 112 111 

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 9855      6

table(indels$V1)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
# 937 902 613 757 588 469 465 577 347 498 551 497 407 288 280 359 223 359 131 235 
# 21  22  23 
# 149 112 111 


# X    X chromosome                    -> 23
# Y    Y chromosome                    -> 24
# XY   Pseudo-autosomal region of X    -> 25
# MT   Mitochondrial                   -> 26

mt<-bim[which(bim$V1 %in% c(0,26)),]
table(mt$V1)
# 0

chry<-bim[which(bim$V1 %in% c(24)),]
table(chry$V1)
# 24 
# 1480

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
# 937 902 613 757 588 469 465 577 347 498 551 497 407 288 280 359 223 359 131 235 
# 21  22  23 
# 149 112 111 

dim(all_remove)
#[1] 9855      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/gwas1/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

##########################################################
# KEEP ONLY THE CASES AND CONTROLS AS DONE IN PREVIOUS QC

## R

keep<-read.table("/path/to/ibdgwas/pre_imputation/raw/wtccc-hrc/liftOver/sample_id/CD+CTRL.inds",head=F)
keep$V1<-gsub(".*_","",keep$V1)
keep$V2<-keep$V1
# write.table(keep,paste(path,"pre_imputation/raw/gwas1/list_cd_ctrl_samples.txt"),col.names=F,row.names=F,quote=F,sep="\t")

# NOTE THERE ARE 54 UC AMONG THE SAMPLES TO KEEP, SEE BELOW

############################################
# edit fam files to add CD and CTR info:

### R

fam<-read.table(paste(path,"pre_imputation/QC/gwas1/gwas1_hg19.fam",sep=""),head=F)

cd<-read.table("/nfs/wtccc/wtccc_analysis/data/CD.inds",head=F)
ctr<-read.table("/nfs/wtccc/wtccc_analysis/data/CTRL.inds",head=F)

dim(keep[which(keep$V1 %in% fam$V1),])
#[1] 4684    2
dim(keep)
#[1] 5009    2
dim(keep[which(!keep$V1 %in% fam$V1),])
#[1] 325   2 # not available data in FAM files, excluded in a previous step?

dim(fam[which(fam$V1 %in% cd$V1),])
#[1] 1748    6
dim(fam[which(fam$V1 %in% ctr$V1),])
#[1] 2938    6

# SEX:
table(fam$V5,useNA="ifany")
#   1    2 
#2126 2558

############################################
# double check assigned pheno from phenoDB:

pheno_gwas1<-read.table(paste(path,"pheno/gwas1.pheno",sep=""),head=T,sep="\t")

dim(fam[which(fam$V1 %in% pheno_gwas1$FID),])
#[1] 1747    6   ### THERE SHOULD BE 1748!

tmp<-fam[which(fam$V1 %in% cd$V1),]

tmp[which(!tmp$V1 %in% pheno_gwas1$FID),]
#             V1         V2 V3 V4 V5 V6
#9438 sample_id sample_id  0  0  2 -9 # 

study<-fread(paste(path,"pheno/phenodb_phenotypes_20180221.tsv",sep=""),head=T,sep="\t")
study<-as.data.frame(study)

study[which(study$sanger_sample_id=="sample_id"),]
#NA

# however it is among the post imputation samples:
# grep sample_id /path/to/ibdgwas/post_imputation/gwas1/gwas1.sample 
# sample_id sample_id 0 2 1 -0.00410341 -0.000808381 0.00129952 -0.00521764 -0.0150869 -0.00447408 0.00254699 0.0136704 0.0325098 -0.00244892

# when the search is performed in the website, for samples incluced in gwas1, it appears:
# pdb4783	sample_id, sample_id, sample_id, sample_id, sample_id, sample_id

# but in the data we retrieve from the flat files it is only:
study[which(study$phenodb_id=="4783"),]
#       phenodb_id study_name sanger_sample_id Age_at_diagnosis
#12927       4783     gwas1      sample_id               NA
# Crohn's Disease


# HOWEVER SAMPLE IS NOT INCLUDED IN FAM FILE:
fam[which(fam$V1=="sample_id"),]
#[1] V1 V2 V3 V4 V5 V6
#<0 rows> (or 0-length row.names)

table(fam$V6)
#   -9 
# 4684 

fam$V6[which(fam$V1 %in% ctr$V1)]<-1
fam$V6[which(fam$V1 %in% pheno_gwas1$FID[which(pheno_gwas1$pheno_ibd==2)])]<-2
fam$V6[which(fam$V1 %in% c("sample_id"))]<-2

table(fam$V6)
#  -9    1    2 
#   6 2936 1742 

pheno_gwas1[which(pheno_gwas1$FID %in% fam$V1[which(fam$V6=="-9")]),]
#             FID         IID pheno_ibd pheno_uc pheno_cd
#777  sample_id sample_id        -9       -9       -9
#804  sample_id sample_id        -9       -9       -9
#1157 sample_id sample_id        -9       -9       -9
#1995  sample_id  sample_id       -9       -9       -9
#2015 sample_id sample_id        -9       -9       -9
#2027 sample_id sample_id        -9       -9       -9

study[which(study$sanger %in% pheno_gwas1[which(pheno_gwas1$FID %in% fam$V1[which(fam$V6=="-9")]),c("FID")]),]

# all defined as cases (1) in last analysis, bin1 is the pheno, 
# see /path/to/project

#grep -E 'sample_id|sample_id|sample_id|sample_id|sample_id|sample_id' /path/to/ibdgwas/post_imputation/gwas1/gwas1.sample 
#ID_1 ID_2 missing sex bin1 cov1 cov2 cov3 cov4 cov5 cov6 cov7 cov8 cov9 cov10
#sample_id sample_id 0 1 1 -0.00311595 0.0040771 -0.00571062 0.00452708 -0.0176381 -0.0159592 0.0192715 0.0125079 0.00598984 0.00396433
#sample_id sample_id 0 2 1 0.000692544 0.0104836 -0.00352755 -0.00407356 -0.0141501 -0.00759309 0.00399558 0.00120655 -0.00658289 -0.00187177
#sample_id sample_id 0 2 1 0.0104388 0.00048112 0.00196847 -0.0122741 0.0110871 0.0146428 -0.0128835 0.00587426 -0.00365975 -0.00552439
#sample_id sample_id 0 2 1 0.00119136 0.00131514 -0.000576971 0.00394638 0.00400345 0.00552555 0.00394038 -0.0135627 0.0131808 0.0167976
#sample_id sample_id 0 1 1 0.00669121 0.000910273 0.00572968 -0.00174094 0.0148642 0.026408 0.00581771 -0.00453025 0.0123796 -0.018945
#sample_id sample_id 0 2 1 -0.00046699 -0.00320205 -0.0076374 0.0293367 0.0137688 -0.0142926 -0.0118935 0.0103862 0.0141927 0.0146828

# keep them as cases, two with data for oral steroids (1), other 3 as IBD_affection_status =1 

# AFTER TALKING WITH CARL AND ALEX, THEY RELY ON HOW STATUS WAS ASIGNED BEFORE, so use that info to update fam files

sample<-read.table("/path/to/ibdgwas/post_imputation/GWAS1/GWAS1.sample",head=T)
sample<-sample[-1,] # remove subheader
dim(sample)
#[1] 4126   15
dim(fam)
#[1] 4684    6

dim(fam[which(fam$V1 %in% sample$ID_1),])
#[1] 4126    6

fam_subset<-fam[which(fam$V1 %in% sample$ID_1),]
fam_subset<-fam_subset[match(sample$ID_1,fam_subset$V1),]
table(as.character(fam_subset$V1)==as.character(sample$ID_1))
#TRUE 
#4126 

table(as.character(fam_subset$V6),as.character(sample$bin1))
#      0    1
#-9    0    6
#1  2920    0
#2     0 1200

table(as.character(fam_subset$V5)==as.character(sample$sex))
#TRUE 
#4126 

# update according to sample file used in analyiss
fam$V6[which(fam$V6=="-9")]<-2

################################################
# Affection status, by default, should be coded:
# -9 missing 
# 0 missing
# 1 unaffected
# 2 affected

table(fam$V6)
# 1    2 
# 2936 1748 

################################################
# Sex (1=male; 2=female; other=unknown)

table(fam$V5)
# 1    2 
# 2126 2558 

fam_test<-read.table(paste(path,"/pre_imputation/QC/gwas1/gwas1_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 4684

write.table(fam,paste(path,"pre_imputation/QC/gwas1/gwas1_hg19_edited.fam",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 0   2

write.table(samples_remove,paste(path,"pre_imputation/QC/gwas1/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/gwas1/gwas1_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/gwas1/gwas1_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/gwas1/gwas1_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/gwas1/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/gwas1/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/gwas1/gwas1_hg19_noind
# 4684 people (2126 males, 2558 females) loaded from .fam.
# 459426 variants and 4684 people pass filters and QC.
# Among remaining phenotypes, 1748 are cases and 2936 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37



#######################################################################################################################################
