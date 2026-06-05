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

# -rw-r-----  1 username team152 67904111 Mar 25  2011  Inf1_20110228.tgz
# -rw-r-----  1 username team152       59 Mar 25  2011  Inf1_20110228.tgz.md5
# -rw-r-----  1 username team152 14961652 Mar 25  2011  Inf2_20110228.tgz
# -rw-r-----  1 username team152       59 Mar 25  2011  Inf2_20110228.tgz.md5
# -rw-r-----  1 username team152    68201 Mar  1  2011 'InfosAllIndividualsInfinium1&2.xlsx'
# -rw-r--r--  1 username team152       45 Mar 28  2011  README.txt

# These files were uploaded by Emilie Theatre.

md5sum -c Inf1_20110228.tgz.md5
# Inf1_20110228.tgz: OK 

gunzip Inf1_20110228.tgz 
tar -xvf Inf1_20110228.tar

# -rw-r-----  1 username team152 107490098 Mar  1  2011  Inf1_20110228.bed
# -rw-r-----  1 username team152   8682845 Mar  1  2011  Inf1_20110228.bim
# -rw-r-----  1 username team152     30901 Mar  1  2011  Inf1_20110228.fam
# -rw-r-----  1 username team152 116213760 Mar 25  2011  Inf1_20110228.tar
# -rw-r-----  1 username team152        59 Mar 25  2011  Inf1_20110228.tgz.md5
# -rw-r-----  1 username team152  19778075 Mar  1  2011  Inf2_20110228.bed
# -rw-r-----  1 username team152   8340562 Mar  1  2011  Inf2_20110228.bim
# -rw-r-----  1 username team152      6454 Mar  1  2011  Inf2_20110228.fam
# -rw-r-----  1 username team152  28129280 Mar 25  2011  Inf2_20110228.tar
# -rw-r-----  1 username team152        59 Mar 25  2011  Inf2_20110228.tgz.md5
# -rw-r-----  1 username team152     68201 Mar  1  2011 'InfosAllIndividualsInfinium1&2.xlsx'
# -rw-r--r--  1 username team152        45 Mar 28  2011  README.txt

# ## explore files:
# wc -l Inf2_20110228.bim
# # 290854 Inf2_20110228.bim
# wc -l Inf1_20110228.bim
# # 302789 Inf1_20110228.bim

# explore overlap between both:

# ### /software/R-4.3.1/bin/R
# 
# bim1<-read.table("/path/to/project",head=F)
# bim2<-read.table("/path/to/project",head=F)
# 
# dim(bim1)
# # [1] 302789      6
# 
# dim(bim2)
# # [1] 290854      6
# 
# dim(bim1[which(bim1$V2 %in% bim2$V2),])
# # [1] 285294      6
# 
# dim(bim2[which(!bim2$V2 %in% bim1$V2),])
# # [1] 5560    6
# dim(bim1[which(!bim1$V2 %in% bim2$V2),])
# # [1] 17495     6
# 
# # explore build:
# head(bim1)
# V1         V2 V3      V4 V5 V6
# 1  1  rs3934834  0 1045729  T  C
# 2  1  rs3737728  0 1061338  T  C
# 3  1  rs6687776  0 1070488  T  C

# keep separated



# excel file includes two spreadsheets with info for each array, saved as two separated csv files:
# -rw-r--r-- 1 username team152 44415 Feb 20 11:20 InfosAllIndividualsInfinium1.csv
# -rw-r--r-- 1 username team152  8734 Feb 20 11:20 InfosAllIndividualsInfinium2.csv


# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/belgium_inf1_old_gwas/
  
cp Inf1_20110228.* /path/to/ibdgwas/IIBDGC/pre_imputation/raw/belgium_inf1_old_gwas/

  
path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/
  

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/belgium_inf1_old_gwas/Inf1_20110228 \
--freq \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17
# 302789 variants and 1417 people pass filters and QC.
# Among remaining phenotypes, 517 are cases and 900 are controls.




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

bim<-read.table(paste(path,"pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#       A     C     G     T
# 0    11    19    14    15
# A     0 14365 65900     0
# C 13154     0     0 57864
# G 58083     0     0 13135
# T     0 66038 14191     0

table(bim$V1)
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 22414 24261 20634 18074 18407 19922 15882 17505 15220 14963 14018 14424 10961 
# 14    15    16    17    18    19    20    21    22    23 
# 9431  8554  8667  8033 10043  5664  7575  5264  5304  7569 

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 59      6

table(indels$V1)
# 1  2  3  4  5  6  7  8  9 10 11 12 14 15 16 22 
# 2  4  2  3  3  8  1  4  3  3  5 10  1  5  4  1 

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
# [1] 59  6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/belgium_inf1_old_gwas/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

pheno<-read.csv("/path/to/project",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17.fam",sep=""),head=F)
table(nrow(pheno)==nrow(fam))
# TRUE 
# 1 

test<-merge(fam,pheno,by.y="Inf1_IndividualID",by.x="V2",sort=F)

################################################
# Sex (1=male; 2=female; other=unknown)

table(test$Inf1_Gender,test$V5,useNA="ifany")
#     0   1   2
# -  89   0   0
# F   0   0 367
# M   0 961   0

################################################
# Affection status, by default, should be coded:
# -9 missing 
# 0 missing
# 1 unaffected
# 2 affected

table(test$Inf1_Phenotype,test$V6,useNA="ifany")
#           1   2
# CD        0 517
# healthy 900   0

# same data

# FID and IID do not match, edit:
fam$V1<-fam$V2

write.table(fam,paste(path,"pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 0   2

write.table(samples_remove,paste(path,"pre_imputation/QC/belgium_inf1_old_gwas/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17.bed \
--bim ${path_gwas}/pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17.bim \
--fam ${path_gwas}/pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/belgium_inf1_old_gwas/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/belgium_inf1_old_gwas/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_noind
# 302730 variants and 1417 people pass filters and QC.
# Among remaining phenotypes, 517 are cases and 900 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# NOTE: data in build 35, needs to be liftover first:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_noind \
--recode tab --out ${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_noind

# ${path_gwas}previous_qced_b38/liftover/liftOverx
# usage:
# liftOver oldFile map.chain newFile unMapped

#### R convert map into bed:

map<-read.table(paste(path,"pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_noind.map",sep=""),head=F)

map$chr<-paste("chr",map$V1,sep="")
map$chr[which(map$V1=="23")]<-"chrX"
map$chr[which(map$V1=="25")]<-"chrX"
map$chr[which(map$V1=="24")]<-"chrY"
map$chromStart<-map$V4
map$chromEnd<-map$V4+1

write.table(map[,c(5:7,2)],paste(path,"pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_noind.bed",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#### lift positions

${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_noind.bed \
${path_gwas}previous_qced_b38/liftover/hg17ToHg19.over.chain.gz \
${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_noind_lifted_hg19 \
${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_noind_no_lifted_hg19

cut -f 4 ${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_noind_no_lifted_hg19 | sed "/^#/d" \
> ${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/nonlifted_variants_to_exclude.dat

wc -l ${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/nonlifted_variants_to_exclude.dat
# 30

### exclude non lifted:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--file ${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_noind \
--exclude ${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/nonlifted_variants_to_exclude.dat \
--recode tab --out ${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_noind_liftedvariants
# 302700 variants and 1417 people pass filters and QC.

## R

bed_lifted<-read.table(paste(path,"pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_noind_lifted_hg19",sep=""),head=F)
map<-read.table(paste(path,"pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_noind_liftedvariants.map",sep=""),head=F)

table(map$V2==bed_lifted$V4)
# TRUE 
# 302700 

map$pos<-bed_lifted$V2
table(map$V2==bed_lifted$V4)
table(map$pos==bed_lifted$V2)
# TRUE 
# 302700 

write.table(map[,c(1:3,5)],paste(path,"pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_noind_liftedvariants_edited.map",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


########

# put together updated .ped plus lifted.map
/path/to/software/plink_linux_x86_64_20181202/./plink \
--ped ${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_noind_liftedvariants.ped \
--map ${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_noind_liftedvariants_edited.map \
--make-bed --out ${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg19_noind
# 302700 variants and 1417 people pass filters and QC

# remove intermediate files:

rm ${path_gwas}pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg17_noind_*

#######################################################################################################################################
