# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
############
# OLD GWAS #
############

# drwxrwsr-x 2 username team152 4096 Mar 28  2011 belgian
# drwxrwsr-x 2 username team152 4096 Apr  5  2011 chop
# -rw-rw---- 1 username team152  320 Mar  7  2011 DATA_USE_REQUIREMENTS.txt
# drwxrwsr-x 2 username team152 4096 Mar 28  2011 german
# drwxrwsr-x 2 username team152 4096 Jul  6  2012 hapmap3-imputation
# drwxrwsr-x 2 username team152 4096 Dec 19  2011 mds
# drwxrwsr-x 2 username team152 4096 Mar 28  2011 newids
# drwxrwsr-x 2 username team152 4096 Mar 25  2011 niddk_cd
# drwxrwsr-x 2 username team152 4096 Mar 25  2011 niddk_uc
# drwxrwsr-x 2 username team152 4096 Jan 24  2013 pkgs
# drwxrwsr-x 2 username team152 4096 Jul 14  2011 qc
# drwxrwsr-x 2 username team152 4096 Mar 25  2011 swedish
# drwxrwsr-x 2 username team152 4096 Mar 25  2011 wtccc_cd
# drwxrwsr-x 2 username team152 4096 Mar 25  2011 wtccc_uc

# Use of these datasets is restricted to projects approved by the Management
# Committee of the International IBD Genetics Consortium.  In addition, any use
# of the WTCCC data requires submitting an application form to the WTCCC
# Consortium Data Access Committee (CDAC) and signing the corresponding Data
# Use Agreement (DUA).
# DATA_USE_REQUIREMENTS.txt 

# /path/to/project


##############
# swedish uc #
##############

# These files were uploaded by Mauro D'Amato.

md5sum -c db126-2356samples-297031snps-uc1-uc2-eira1_ctl_wo_UIBD_wo_RAW.ped.gz.md5
# db126-2356samples-297031snps-uc1-uc2-eira1_ctl_wo_UIBD_wo_RAW.ped.gz: OK

gunzip db126-2356samples-297031snps-uc1-uc2-eira1_ctl_wo_UIBD_wo_RAW.ped.gz

gunzip Swedish_UC_GWAS.map.gz
gunzip Swedish_UC_GWAS.ped.gz

## explore files:
wc -l Swedish_UC_GWAS.map
# 297031 Swedish_UC_GWAS.map

head Swedish_UC_GWAS.map
# 1	rs3934834	0	995669
# 1	rs3737728	0	1011278
# in hg18


head Swedish_UC_GWAS.ped 
# FAM715800346 715800346 0 0 1 2
# FAM715800348 715800348 0 0 2 2
# FAM715800349 715800349 0 0 1 2
# FAM FILE INSTEAD OF PED!!

# RENAME TWO FILES:
mv Swedish_UC_GWAS.ped Swedish_UC_GWAS.fam
mv db126-2356samples-297031snps-uc1-uc2-eira1_ctl_wo_UIBD_wo_RAW.ped Swedish_UC_GWAS.ped

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/swedish_uc_old_gwas/
  
cp Swedish_UC_GWAS.* /path/to/ibdgwas/IIBDGC/pre_imputation/raw/swedish_uc_old_gwas/

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--ped /path/to/ibdgwas/IIBDGC/pre_imputation/raw/swedish_uc_old_gwas/Swedish_UC_GWAS.ped \
--map /path/to/ibdgwas/IIBDGC/pre_imputation/raw/swedish_uc_old_gwas/Swedish_UC_GWAS.map \
--freq \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18
# 297031 variants and 1264 people pass filters and QC.
# Among remaining phenotypes, 923 are cases and 341 are controls.



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

bim<-read.table(paste(path,"pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        A      C      G
# A      0  28414 128532
# C  26614      0      0
# G 113471      0      0

table(bim$V1)
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 22640 24409 20856 18248 18453 20120 16154 17475 15289 15049 14108 14476 11146 
# 14    15    16    17    18    19    20    21    22 
# 9533  8509  8682  8075 10054  5727  7480  5296  5252 

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 0  6

table(indels$V1)
# 0 
# 0

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


dim(all_remove)
# [1] 0  6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/swedish_uc_old_gwas/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

# pheno<-read.csv("/path/to/project",head=T)
# no addtional pheno file for this cohort

fam<-read.table(paste(path,"pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18.fam",sep=""),head=F)
# table(nrow(pheno)==nrow(fam))
# # TRUE 
# # 1 

# test<-merge(fam,pheno,by.y="Inf2_IndividualID",by.x="V2",sort=F)

################################################
# Sex (1=male; 2=female; other=unknown)

# table(test$Inf2_Gender,test$V5,useNA="ifany")
table(fam$V5,useNA="ifany")
# 1   2 
# 697 567 


################################################
# Affection status, by default, should be coded:
# -9 missing 
# 0 missing
# 1 unaffected
# 2 affected

# table(test$Inf2_Phenotype,test$V6,useNA="ifany")
table(fam$V6,useNA="ifany")
# 1   2 
# 341 923 

# edit FID field, to reflect a sample per family:
# table(fam$V1)
# 0 
# 550

fam$V2<-fam$V1
write.table(fam,paste(path,"pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 0   2

write.table(samples_remove,paste(path,"pre_imputation/QC/swedish_uc_old_gwas/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18.bed \
--bim ${path_gwas}/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18.bim \
--fam ${path_gwas}/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/swedish_uc_old_gwas/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/swedish_uc_old_gwas/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_noind
# 297031 variants and 1264 people pass filters and QC.
# Among remaining phenotypes, 923 are cases and 341 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# NOTE: data in build 18, needs to be liftover first:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_noind \
--recode tab --out ${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_noind

# ${path_gwas}previous_qced_b38/liftover/liftOverx
# usage:
# liftOver oldFile map.chain newFile unMapped

#### R convert map into bed:

map<-read.table(paste(path,"pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_noind.map",sep=""),head=F)

map$chr<-paste("chr",map$V1,sep="")
map$chr[which(map$V1=="23")]<-"chrX"
map$chr[which(map$V1=="25")]<-"chrX"
map$chr[which(map$V1=="24")]<-"chrY"
map$chromStart<-map$V4
map$chromEnd<-map$V4+1

map$chromStart<-format(map$chromStart, scientific=F)
map$chromEnd<-format(map$chromEnd, scientific=F)


write.table(map[,c(5:7,2)],paste(path,"pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_noind.bed",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#### lift positions

${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_noind.bed \
${path_gwas}previous_qced_b38/liftover/hg18ToHg19.over.chain.gz \
${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_noind_lifted_hg19 \
${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_noind_no_lifted_hg19

cut -f 4 ${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_noind_no_lifted_hg19 | sed "/^#/d" \
> ${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/nonlifted_variants_to_exclude.dat

wc -l ${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/nonlifted_variants_to_exclude.dat
# 27

### exclude non lifted:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--file ${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_noind \
--exclude ${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/nonlifted_variants_to_exclude.dat \
--recode tab --out ${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_noind_liftedvariants
# 297004 variants and 1264 people pass filters and QC.

## R

bed_lifted<-read.table(paste(path,"pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_noind_lifted_hg19",sep=""),head=F)
map<-read.table(paste(path,"pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_noind_liftedvariants.map",sep=""),head=F)

table(map$V2==bed_lifted$V4)
# TRUE 
# 297004 

map$pos<-bed_lifted$V2
table(map$V2==bed_lifted$V4)
table(map$pos==bed_lifted$V2)
# TRUE 
# 297004 

write.table(map[,c(1:3,5)],paste(path,"pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_noind_liftedvariants_edited.map",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


########

# put together updated .ped plus lifted.map
/path/to/software/plink_linux_x86_64_20181202/./plink \
--ped ${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_noind_liftedvariants.ped \
--map ${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_noind_liftedvariants_edited.map \
--make-bed --out ${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_noind
# 297004 variants and 1264 people pass filters and QC.
# Among remaining phenotypes, 923 are cases and 341 are controls.

# remove intermediate files:

rm ${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18_noind_*
  
  
# alleles are recoded, so now there are variants (monomorphic) with 0 alleles that need to be excluded:
  
#### R
  
bim<-read.table(paste(path,"pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_noind.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        A      C      G
# A      0  28411 128519
# C  26614      0      0
# G 113460      0      0

table(bim[which(bim$V5=="0"),"V1"])
#

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 0   6

write.table(indels[,"V2",drop=F],paste(path,"pre_imputation/QC/swedish_uc_old_gwas/list_indel_var_exclude_2",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_noind \
--exclude ${path_gwas}/pre_imputation/QC/swedish_uc_old_gwas/list_indel_var_exclude_2 \
--remove ${path_gwas}/pre_imputation/QC/swedish_uc_old_gwas/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_noind
# 297004 variants and 1264 people pass filters and QC.
# Among remaining phenotypes, 923 are cases and 341 are controls.
