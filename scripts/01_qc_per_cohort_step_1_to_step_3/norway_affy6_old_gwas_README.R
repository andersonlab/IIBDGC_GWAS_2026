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

#################
# NORWAY Affy 6 #
#################

# These files were uploaded by David Ellinghaus.

md5sum -c Uc_Nor_Affy6.tgz.md5
# md5sum -c Uc_Nor_Affy6.tgz.md5: OK

gunzip Uc_Nor_Affy6.tgz
tar -xvf Uc_Nor_Affy6.tar

Uc_Nor_Affy6/
# -rw-r----- 1 username team152       134 Feb 15  2011 README
# -rw-r----- 1 username team152 128612829 Feb 15  2011 Uc_Nor_Affy6.bed
# -rw-r----- 1 username team152  26023094 Feb 15  2011 Uc_Nor_Affy6.bim
# -rw-r----- 1 username team152     13893 Feb 15  2011 Uc_Nor_Affy6.fam
# -rw-r----- 1 username team152   1049857 Feb 15  2011 Uc_Nor_Affy6.hh
# -rw-r----- 1 username team152      3257 Feb 15  2011 Uc_Nor_Affy6.log
# -rw-r----- 1 username team152    263234 Feb 15  2011 Uc_Nor_Affy6.nof
  
# Uc_Nor_Affy6: Ulcerative colitis GWAS from Norway, Affymetrix6.0 snp array.
# Alleles in Uc_Nor_Affy6.bim are annotated on plus strand.

## explore files:
wc -l Uc_Nor_Affy6.bim
# 931977 Uc_Nor_Affy6.bim

head -5000 Uc_Nor_Affy6.bim | tail -10
# 1	rs938249	0	15209109	C	T
# 1	rs6429703	0	15212547	T	C

# in hg18

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/norway_affy6_old_gwas/
  
cp Uc_Nor_Affy6.* /path/to/ibdgwas/IIBDGC/pre_imputation/raw/norway_affy6_old_gwas/
  
path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/norway_affy6_old_gwas/Uc_Nor_Affy6 \
--freq \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18
# 931977 variants and 550 people pass filters and QC.
# Among remaining phenotypes, 268 are cases and 282 are controls.



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

bim<-read.table(paste(path,"pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        A      C      G      T
# A      0  40072 177193  29267
# C  33241      0  42680 143381
# G 143393  42715      0  33418
# T  29551 177368  39698      0

table(bim$V1)
# 0     1     2     3     4     5     6     7     8     9    10    11    12 
# 792 73234 75706 62116 57363 57764 57659 48244 49849 42582 49458 45784 43671 
# 13    14    15    16    17    18    19    20    21    22    23    24    26 
# 35044 28830 26822 28569 21272 27117 12380 23420 12906 12027 38001   902   465 

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
# 0  26 
# 792 465


chry<-bim[which(bim$V1 %in% c(24)),]
table(chry$V1)
# 24 
# 902 

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])


dim(all_remove)
# [1] 1257  6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/norway_affy6_old_gwas/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

# pheno<-read.csv("/path/to/project",head=T)
# no addtional pheno file for this cohort

fam<-read.table(paste(path,"pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18.fam",sep=""),head=F)
# table(nrow(pheno)==nrow(fam))
# # TRUE 
# # 1 

# test<-merge(fam,pheno,by.y="Inf2_IndividualID",by.x="V2",sort=F)

################################################
# Sex (1=male; 2=female; other=unknown)

# table(test$Inf2_Gender,test$V5,useNA="ifany")
table(fam$V5,useNA="ifany")
# 1   2 
# 305 245 


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
table(fam$V1)
# 0 
# 550

fam$V1<-fam$V2
write.table(fam,paste(path,"pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 0   2

write.table(samples_remove,paste(path,"pre_imputation/QC/norway_affy6_old_gwas/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18.bed \
--bim ${path_gwas}/pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18.bim \
--fam ${path_gwas}/pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/norway_affy6_old_gwas/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/norway_affy6_old_gwas/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_noind
# 930720 variants and 550 people pass filters and QC.
# Among remaining phenotypes, 268 are cases and 282 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# NOTE: data in build 18, needs to be liftover first:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_noind \
--recode tab --out ${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_noind

# ${path_gwas}previous_qced_b38/liftover/liftOverx
# usage:
# liftOver oldFile map.chain newFile unMapped

#### R convert map into bed:

map<-read.table(paste(path,"pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_noind.map",sep=""),head=F)

map$chr<-paste("chr",map$V1,sep="")
map$chr[which(map$V1=="23")]<-"chrX"
map$chr[which(map$V1=="25")]<-"chrX"
map$chr[which(map$V1=="24")]<-"chrY"
map$chromStart<-map$V4
map$chromEnd<-map$V4+1

map$chromStart<-format(map$chromStart, scientific=F)
map$chromEnd<-format(map$chromEnd, scientific=F)


write.table(map[,c(5:7,2)],paste(path,"pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_noind.bed",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#### lift positions

${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_noind.bed \
${path_gwas}previous_qced_b38/liftover/hg18ToHg19.over.chain.gz \
${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_noind_lifted_hg19 \
${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_noind_no_lifted_hg19

cut -f 4 ${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_noind_no_lifted_hg19 | sed "/^#/d" \
> ${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/nonlifted_variants_to_exclude.dat

wc -l ${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/nonlifted_variants_to_exclude.dat
# 114

### exclude non lifted:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--file ${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_noind \
--exclude ${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/nonlifted_variants_to_exclude.dat \
--recode tab --out ${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_noind_liftedvariants
# 930606 variants and 550 people pass filters and QC.

## R

bed_lifted<-read.table(paste(path,"pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_noind_lifted_hg19",sep=""),head=F)
map<-read.table(paste(path,"pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_noind_liftedvariants.map",sep=""),head=F)

table(map$V2==bed_lifted$V4)
# TRUE 
# 930606 

map$pos<-bed_lifted$V2
table(map$V2==bed_lifted$V4)
table(map$pos==bed_lifted$V2)
# TRUE 
# 930606 

write.table(map[,c(1:3,5)],paste(path,"pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_noind_liftedvariants_edited.map",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


########

# put together updated .ped plus lifted.map
/path/to/software/plink_linux_x86_64_20181202/./plink \
--ped ${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_noind_liftedvariants.ped \
--map ${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_noind_liftedvariants_edited.map \
--make-bed --out ${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg19_noind
# 930606 variants and 550 people pass filters and QC.
# Among remaining phenotypes, 268 are cases and 282 are controls.

# remove intermediate files:

rm ${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18_noind_*
  
  
# alleles are recoded, so now there are variants (monomorphic) with 0 alleles that need to be excluded:
  
#### R
  
bim<-read.table(paste(path,"pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg19_noind.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        0      A      C      G      T
# 0  24328  13265  16911  16722  13092
# A      0      0  36285 160961  26718
# C      0  30324      0  38443 129946
# G      0 130056  38670      0  30440
# T      0  27045 161207  36193      0

table(bim[which(bim$V5=="0"),"V1"])
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 6869 6788 5286 5146 4715 4506 3837 4363 3453 4313 4127 3728 2961 2481 2525 2861 
# 17   18   19   20   21   22   23   24 
# 2060 2428 1181 2147 1055 1230 5540  718 

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 84318   6

write.table(indels[,"V2",drop=F],paste(path,"pre_imputation/QC/norway_affy6_old_gwas/list_indel_var_exclude_2",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg19_noind \
--exclude ${path_gwas}/pre_imputation/QC/norway_affy6_old_gwas/list_indel_var_exclude_2 \
--remove ${path_gwas}/pre_imputation/QC/norway_affy6_old_gwas/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg19_noind
# 846288 variants and 550 people pass filters and QC.
# Among remaining phenotypes, 268 are cases and 282 are controls.
