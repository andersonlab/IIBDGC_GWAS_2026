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

# 204 samples contributed by Cedars-Sinai are not included here, but will be
# supplied separately by that center.


# NOTE: this file includes QC steps followed with NIDDK CD and UC all togegher, both files merged after liftover step. In the past both cohorts were 
# merged and analysed together


# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/niddk_old_gwas/
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/
  
  
############
# niddk cd #
############

md5sum -c cdgwa_no_cedars.tgz.md5
# cdgwa_no_cedars.tgz: OK

gunzip cdgwa_no_cedars.tgz
tar -xvf cdgwa_no_cedars.tar

# cdgwa_no_cedars.bed
# cdgwa_no_cedars.bim
# cdgwa_no_cedars.fam


## explore files:
wc -l cdgwa_no_cedars.bim
# 317503 cdgwa_no_cedars.bim

head -5000 cdgwa_no_cedars.bim | tail -10
# 1	rs923028	0	48766900	A	G
# 1	rs319960	0	48771635	A	G

# in hg17

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/niddk_cd_old_gwas/
  
cp cdgwa_no_cedars.* /path/to/ibdgwas/IIBDGC/pre_imputation/raw/niddk_cd_old_gwas/
  
path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/niddk_cd_old_gwas/cdgwa_no_cedars \
--freq \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17
# 317503 variants and 1792 people pass filters and QC.
# Among remaining phenotypes, 833 are cases and 959 are controls.
# 1792 people (863 males, 929 females) loaded from .fam.


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

bim<-read.table(paste(path,"pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        A      C      G
# 0      1      0      0
# A      0  30618 137203
# C  28591      0      0
# G 121090      0      0

table(bim$V1)
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 23275 25351 21580 19113 19272 20811 16675 18274 15835 15592 14660 15032 11526 
# 14    15    16    17    18    19    20    21    22    23    25 
# 9829  8900  9006  8343 10495  5927  7836  5493  5505  9171     2 

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 1  6

table(indels$V1)
# 2 
# 1

# X    X chromosome                    -> 23
# Y    Y chromosome                    -> 24
# XY   Pseudo-autosomal region of X    -> 25
# MT   Mitochondrial                   -> 26

mt<-bim[which(bim$V1 %in% c(0,26)),]
table(mt$V1)



chry<-bim[which(bim$V1 %in% c(24)),]
table(chry$V1)

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])


dim(all_remove)
# [1] 1  6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/niddk_cd_old_gwas/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

# pheno<-read.csv("/path/to/project",head=T)
# no addtional pheno file for this cohort

fam<-read.table(paste(path,"pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17.fam",sep=""),head=F)
# table(nrow(pheno)==nrow(fam))
# # TRUE 
# # 1 

# test<-merge(fam,pheno,by.y="Inf2_IndividualID",by.x="V2",sort=F)

################################################
# Sex (1=male; 2=female; other=unknown)

# table(test$Inf2_Gender,test$V5,useNA="ifany")
table(fam$V5,useNA="ifany")
# 1   2 
# 863 929


################################################
# Affection status, by default, should be coded:
# -9 missing 
# 0 missing
# 1 unaffected
# 2 affected

# table(test$Inf2_Phenotype,test$V6,useNA="ifany")
table(fam$V6,useNA="ifany")
# 1   2 
# 282 268

# edit FID field, to reflect a sample per family:
# table(fam$V1)
# 0 
# 550

# fam$V1<-fam$V2
write.table(fam,paste(path,"pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 0   2

write.table(samples_remove,paste(path,"pre_imputation/QC/niddk_cd_old_gwas/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17.bed \
--bim ${path_gwas}/pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17.bim \
--fam ${path_gwas}/pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/niddk_cd_old_gwas/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/niddk_cd_old_gwas/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_noind
# 1792 people (863 males, 929 females) loaded from .fam.
# 317502 variants and 1792 people pass filters and QC.
# Among remaining phenotypes, 833 are cases and 959 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# NOTE: data in build 18, needs to be liftover first:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_noind \
--recode tab --out ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_noind

# ${path_gwas}previous_qced_b38/liftover/liftOverx
# usage:
# liftOver oldFile map.chain newFile unMapped

#### R convert map into bed:

map<-read.table(paste(path,"pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_noind.map",sep=""),head=F)

map$chr<-paste("chr",map$V1,sep="")
map$chr[which(map$V1=="23")]<-"chrX"
map$chr[which(map$V1=="25")]<-"chrX"
map$chr[which(map$V1=="24")]<-"chrY"
map$chromStart<-map$V4
map$chromEnd<-map$V4+1

map$chromStart<-format(map$chromStart, scientific=F)
map$chromEnd<-format(map$chromEnd, scientific=F)


write.table(map[,c(5:7,2)],paste(path,"pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_noind.bed",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#### lift positions

${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_noind.bed \
${path_gwas}previous_qced_b38/liftover/hg17ToHg19.over.chain.gz \
${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_noind_lifted_hg19 \
${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_noind_no_lifted_hg19

cut -f 4 ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_noind_no_lifted_hg19 | sed "/^#/d" \
> ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/nonlifted_variants_to_exclude.dat

wc -l ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/nonlifted_variants_to_exclude.dat
# 31

### exclude non lifted:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--file ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_noind \
--exclude ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/nonlifted_variants_to_exclude.dat \
--recode tab --out ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_noind_liftedvariants
# 317471 variants and 1792 people pass filters and QC.

## R

bed_lifted<-read.table(paste(path,"pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_noind_lifted_hg19",sep=""),head=F)
map<-read.table(paste(path,"pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_noind_liftedvariants.map",sep=""),head=F)

table(map$V2==bed_lifted$V4)
# TRUE 
# 317471 

map$pos<-bed_lifted$V2
table(map$V2==bed_lifted$V4)
table(map$pos==bed_lifted$V2)
# TRUE 
# 317471 

write.table(map[,c(1:3,5)],paste(path,"pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_noind_liftedvariants_edited.map",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


########

# put together updated .ped plus lifted.map
/path/to/software/plink_linux_x86_64_20181202/./plink \
--ped ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_noind_liftedvariants.ped \
--map ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_noind_liftedvariants_edited.map \
--make-bed --out ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg19_noind
# 317471 variants and 1792 people pass filters and QC.
# Among remaining phenotypes, 833 are cases and 959 are controls.

# remove intermediate files:

rm ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17_noind_*
  
  
# alleles are recoded, so now there are variants (monomorphic) with 0 alleles that need to be excluded:
  
#### R
  
bim<-read.table(paste(path,"pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg19_noind.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        A      C      G
# 0      1      0      1
# A      0  30614 137180
# C  28591      0      0
# G 121084      0      0

table(bim[which(bim$V5=="0"),"V1"])
# 10 18 
# 1  1 

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 2   6

write.table(indels[,"V2",drop=F],paste(path,"pre_imputation/QC/niddk_cd_old_gwas/list_indel_var_exclude_2",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg19_noind \
--exclude ${path_gwas}/pre_imputation/QC/niddk_cd_old_gwas/list_indel_var_exclude_2 \
--remove ${path_gwas}/pre_imputation/QC/niddk_cd_old_gwas/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg19_noind_2
# 317469 variants and 1792 people pass filters and QC.
# Among remaining phenotypes, 833 are cases and 959 are controls.


