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


###########
# BELGIAN #
###########

# /path/to/project

# -rw-r-----  1 username team152 14961652 Mar 25  2011  Inf2_20110228.tgz
# -rw-r-----  1 username team152       59 Mar 25  2011  Inf2_20110228.tgz.md5
# -rw-r-----  1 username team152    68201 Mar  1  2011 'InfosAllIndividualsInfinium1&2.xlsx'
# -rw-r--r--  1 username team152       45 Mar 28  2011  README.txt

# These files were uploaded by Emilie Theatre.

md5sum -c Inf2_20110228.tgz.md5
# Inf2_20110228.tgz: OK

gunzip Inf2_20110228.tgz 
tar -xvf Inf2_20110228.tar

# -rw-r-----  1 username team152  19778075 Mar  1  2011  Inf2_20110228.bed
# -rw-r-----  1 username team152   8340562 Mar  1  2011  Inf2_20110228.bim
# -rw-r-----  1 username team152      6454 Mar  1  2011  Inf2_20110228.fam
# -rw-r-----  1 username team152  28129280 Mar 25  2011  Inf2_20110228.tar
# -rw-r-----  1 username team152        59 Mar 25  2011  Inf2_20110228.tgz.md5
# -rw-r-----  1 username team152     68201 Mar  1  2011 'InfosAllIndividualsInfinium1&2.xlsx'
# -rw-r--r--  1 username team152        45 Mar 28  2011  README.txt

## explore files:
wc -l inf2_20110228.bim
# 302789 inf2_20110228.bim


# excel file includes two spreadsheets with info for each array, saved as two separated csv files:
# -rw-r--r-- 1 username team152 44415 Feb 20 11:20 InfosAllIndividualsInfinium1.csv
# -rw-r--r-- 1 username team152  8734 Feb 20 11:20 InfosAllIndividualsInfinium2.csv


# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/belgium_inf2_old_gwas/

cp Inf2_20110228.* /path/to/ibdgwas/IIBDGC/pre_imputation/raw/belgium_inf2_old_gwas/  

  
path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/
  
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/belgium_inf2_old_gwas/Inf2_20110228 \
--freq \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17
# 290854 variants loaded from .bim file.
# 272 people (59 males, 98 females, 115 ambiguous) loaded from .fam.




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

bim<-read.table(paste(path,"pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#      A     C     G     T
# 0    21    30    26    17
# A     0 14165 63479     0
# C 12774     0     0 54968
# G 55221     0     0 12757
# T     0 63485 13911     0

table(bim$V1)
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 21527 23281 19735 17191 17569 18996 15145 16776 14639 14404 13459 13790 10460 
# 14    15    16    17    18    19    20    21    22    23 
# 9097  8250  8359  7761  9619  5405  7353  5029  5078  7931 

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 94      6

table(indels$V1)


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
# [1] 94  6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/belgium_inf2_old_gwas/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

pheno<-read.csv("/path/to/project",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17.fam",sep=""),head=F)
table(nrow(pheno)==nrow(fam))
# TRUE 
# 1 

test<-merge(fam,pheno,by.y="Inf2_IndividualID",by.x="V2",sort=F)

################################################
# Sex (1=male; 2=female; other=unknown)

table(test$Inf2_Gender,test$V5,useNA="ifany")
#     0   1   2
# - 115   0   0
# F   0   0  98
# M   0  59   0


################################################
# Affection status, by default, should be coded:
# -9 missing 
# 0 missing
# 1 unaffected
# 2 affected

table(test$Inf2_Phenotype,test$V6,useNA="ifany")
#           1   2   3
# CD        0 159   0
# healthy 111   0   0
# UC        0   0   2

# same data; edit fam to be ibd
fam$V6[which(fam$V6=="3")]<-2
# 1   2 
# 111 161 

write.table(fam,paste(path,"pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 0   2

write.table(samples_remove,paste(path,"pre_imputation/QC/belgium_inf2_old_gwas/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17.bed \
--bim ${path_gwas}/pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17.bim \
--fam ${path_gwas}/pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/belgium_inf2_old_gwas/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/belgium_inf2_old_gwas/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_noind
# 290760 variants and 272 people pass filters and QC.
# Among remaining phenotypes, 161 are cases and 111 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# NOTE: data in build 35, needs to be liftover first:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_noind \
--recode tab --out ${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_noind

# ${path_gwas}previous_qced_b38/liftover/liftOverx
# usage:
# liftOver oldFile map.chain newFile unMapped

#### R convert map into bed:

map<-read.table(paste(path,"pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_noind.map",sep=""),head=F)

map$chr<-paste("chr",map$V1,sep="")
map$chr[which(map$V1=="23")]<-"chrX"
map$chr[which(map$V1=="25")]<-"chrX"
map$chr[which(map$V1=="24")]<-"chrY"
map$chromStart<-map$V4
map$chromEnd<-map$V4+1

write.table(map[,c(5:7,2)],paste(path,"pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_noind.bed",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#### lift positions

${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_noind.bed \
${path_gwas}previous_qced_b38/liftover/hg17ToHg19.over.chain.gz \
${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_noind_lifted_hg19 \
${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_noind_no_lifted_hg19

cut -f 4 ${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_noind_no_lifted_hg19 | sed "/^#/d" \
> ${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/nonlifted_variants_to_exclude.dat

wc -l ${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/nonlifted_variants_to_exclude.dat
# 30

### exclude non lifted:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--file ${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_noind \
--exclude ${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/nonlifted_variants_to_exclude.dat \
--recode tab --out ${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_noind_liftedvariants
# 290730 variants and 272 people pass filters and QC.

## R

bed_lifted<-read.table(paste(path,"pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_noind_lifted_hg19",sep=""),head=F)
map<-read.table(paste(path,"pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_noind_liftedvariants.map",sep=""),head=F)

table(map$V2==bed_lifted$V4)
# TRUE 
# 290730 

map$pos<-bed_lifted$V2
table(map$V2==bed_lifted$V4)
table(map$pos==bed_lifted$V2)
# TRUE 
# 290730 

write.table(map[,c(1:3,5)],paste(path,"pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_noind_liftedvariants_edited.map",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


########

# put together updated .ped plus lifted.map
/path/to/software/plink_linux_x86_64_20181202/./plink \
--ped ${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_noind_liftedvariants.ped \
--map ${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_noind_liftedvariants_edited.map \
--make-bed --out ${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg19_noind
# 290730 variants and 272 people pass filters and QC.

# remove intermediate files:

rm ${path_gwas}pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg17_noind_*
  
#######################################################################################################################################
