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

# -rw-rw----  1 username team152  43722352 Mar 28  2011 cedars_370k_newids.tgz
# -rw-rw----  1 username team152        64 Mar 28  2011 cedars_370k_newids.tgz.md5
# -rw-rw----  1 username team152  81686907 Mar 28  2011 cedars_610k_newids.tgz
# -rw-rw----  1 username team152        64 Mar 28  2011 cedars_610k_newids.tgz.md5
# -rw-rw----  1 username team152 112294186 Mar 28  2011 cedars_omni_exp_newids.tgz
# -rw-rw----  1 username team152        68 Mar 28  2011 cedars_omni_exp_newids.tgz.md5

###############
# cedars omni #
###############

md5sum -c cedars_omni_exp_newids.tgz.md5
# cedars_omni_exp_newids.tgz: OK

gunzip cedars_omni_exp_newids.tgz
tar -xvf cedars_omni_exp_newids.tar

# -rw-r--r--  1 username team152 225067843 Mar 28  2011 cedars_omni_exp_newids.bed
# -rw-r--r--  1 username team152  20419536 Mar 28  2011 cedars_omni_exp_newids.bim
# -rw-r--r--  1 username team152     39296 Mar 28  2011 cedars_omni_exp_newids.fam


## explore files:
wc -l cedars_omni_exp_newids.bim
# 733120 cedars_omni_exp_newids.bim

head -5000 cedars_omni_exp_newids.bim | tail -10
# 1	rs2297916	0	17837781	A	G
# 1	rs6665682	0	17841575	A	G

# in hg18

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/cedars_omni_old_gwas/
  
cp cedars_omni_exp_newids.* /path/to/ibdgwas/IIBDGC/pre_imputation/raw/cedars_omni_old_gwas/
  
path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/cedars_omni_old_gwas/cedars_omni_exp_newids \
--freq \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18
# 733120 variants and 1228 people pass filters and QC.
# Among remaining phenotypes, 1228 are cases and 0 are controls.
# 1228 people (0 males, 0 females, 1228 ambiguous) loaded from .fam.



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

bim<-read.table(paste(path,"pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        A      C      D      G      I      T
# 0   3505   1596      1   4807      5    109
# A      0  71035      0 316269      0    555
# C  64027      0      0    870      0      0
# D      0      0      0      0      3      0
# G 269281    671      0      0      0      0
# T    386      0      0      0      0      0

table(bim$V1)
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 59782 58184 47601 40778 42441 48781 38552 37356 33142 39445 36978 35890 28059 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 23532 21876 22997 20504 21874 15289 18574 10340 10676 18236  1693   540 

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 10026  6

table(indels$V1)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 1053  535  441  293  361  812  578  291  330  376  468  482  198  165  264  447 
# 17   18   19   20   21   22   23   24   25 
# 418  170  393  261   70  192  594  799   35 


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
# 1693 

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])


dim(all_remove)
# [1] 10026    6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/cedars_omni_old_gwas/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

# pheno file provided by Talin
pheno<-fread(paste(path,"pheno/list_cedars_old_gwas_ids_TH.txt",sep=""),head=T,sep="\t")
pheno<-as.data.frame(pheno)


fam<-read.table(paste(path,"pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18.fam",sep=""),head=F)
table(pheno$study[which(pheno$sample_id %in% fam$V1)])
# cedars_omni 
# 1228

table(pheno$Sex[which(pheno$sample_id %in% fam$V1)])
#   Female   Male 
# 1    631    596 

table(pheno$IBD_Affection_status[which(pheno$sample_id %in% fam$V1)],pheno$Unrelated_control[which(pheno$sample_id %in% fam$V1)],useNA="ifany")
#               0    1 <NA>
#               0    0    1
# Affected      0    0 1222
# Unaffected    1    4    0

table(pheno$IBD_Affection_status[which(pheno$sample_id %in% fam$V1)],pheno$Disease_type[which(pheno$sample_id %in% fam$V1)],useNA="ifany")
#                Crohn's Disease Indeterminate Ulcerative Colitis
#              1               0             0                  0
# Affected     0             735            70                417
# Unaffected   5               0             0                  0


################################################
# Sex (1=male; 2=female; other=unknown)

table(fam$V5,useNA="ifany")
# 0 
# 1228 

fam$V5[which(fam$V1 %in% pheno$sample_id[which(pheno$Sex=="Male")])]<-1
fam$V5[which(fam$V1 %in% pheno$sample_id[which(pheno$Sex=="Female")])]<-2

table(fam$V5,useNA="ifany")
# 0   1   2 
# 1 596 631 

################################################
# Affection status, by default, should be coded:
# -9 missing 
# 0 missing
# 1 unaffected
# 2 affected

table(fam$V6,useNA="ifany")
# 2 
# 1228

fam[,"V6"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$Unrelated_control=="1")]),"V6"]<-"1"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$IBD_Affection_status=="Affected")]),"V6"]<-"2"

table(fam$V6)
# 1    2   -9 
# 4 1222    2

write.table(fam,paste(path,"pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

# exclude also from analysis sample sample_id 
pheno[which(pheno$sample_id=="sample_id"),]
# sample_id       study sample_comment Genetic_ID Genetic_ID_2 Sex
# 1480 sample_id cedars_omni     withdrawn?    sample_id                 
# IBD_Affection_status Disease_type Unrelated_control self_race
# 1480                                                  NA          
# self_hispanic self_jewish
# 1480     

samples_remove<-fam[which((fam$V6=="-9") | (fam$V1=="sample_id")),c(1,2)]
dim(samples_remove)
# [1] 2   2


write.table(samples_remove,paste(path,"pre_imputation/QC/cedars_omni_old_gwas/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18.bed \
--bim ${path_gwas}/pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18.bim \
--fam ${path_gwas}/pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/cedars_omni_old_gwas/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/cedars_omni_old_gwas/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_noind
# 723094 variants and 1226 people pass filters and QC.
# Among remaining phenotypes, 1222 are cases and 4 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# NOTE: data in build 18, needs to be liftover first:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_noind \
--recode tab --out ${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_noind

# ${path_gwas}previous_qced_b38/liftover/liftOverx
# usage:
# liftOver oldFile map.chain newFile unMapped

#### R convert map into bed:

map<-read.table(paste(path,"pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_noind.map",sep=""),head=F)

map$chr<-paste("chr",map$V1,sep="")
map$chr[which(map$V1=="23")]<-"chrX"
map$chr[which(map$V1=="25")]<-"chrX"
map$chr[which(map$V1=="24")]<-"chrY"
map$chromStart<-map$V4
map$chromEnd<-map$V4+1

map$chromStart<-format(map$chromStart, scientific=F)
map$chromEnd<-format(map$chromEnd, scientific=F)


write.table(map[,c(5:7,2)],paste(path,"pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_noind.bed",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#### lift positions

${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_noind.bed \
${path_gwas}previous_qced_b38/liftover/hg18ToHg19.over.chain.gz \
${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_noind_lifted_hg19 \
${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_noind_no_lifted_hg19

cut -f 4 ${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_noind_no_lifted_hg19 | sed "/^#/d" \
> ${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/nonlifted_variants_to_exclude.dat

wc -l ${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/nonlifted_variants_to_exclude.dat
# 87

### exclude non lifted:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--file ${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_noind \
--exclude ${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/nonlifted_variants_to_exclude.dat \
--recode tab --out ${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_noind_liftedvariants
# 723007 variants and 1226 people pass filters and QC.

## R

bed_lifted<-read.table(paste(path,"pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_noind_lifted_hg19",sep=""),head=F)
map<-read.table(paste(path,"pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_noind_liftedvariants.map",sep=""),head=F)

table(map$V2==bed_lifted$V4)
# TRUE 
# 723007 

map$pos<-bed_lifted$V2
table(map$V2==bed_lifted$V4)
table(map$pos==bed_lifted$V2)
# TRUE 
# 723007 

write.table(map[,c(1:3,5)],paste(path,"pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_noind_liftedvariants_edited.map",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


########

# put together updated .ped plus lifted.map
/path/to/software/plink_linux_x86_64_20181202/./plink \
--ped ${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_noind_liftedvariants.ped \
--map ${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_noind_liftedvariants_edited.map \
--make-bed --out ${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg19_noind
# 723007 variants and 1226 people pass filters and QC.
# Among remaining phenotypes, 1222 are cases and 4 are controls.

# remove intermediate files:

rm ${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18_noind_*
  
  
# alleles are recoded, so now there are variants (monomorphic) with 0 alleles that need to be excluded:
  
#### R
  
bim<-read.table(paste(path,"pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg19_noind.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        A      C      G      T
# 0      2      0      1      0
# A      0  71042 316257    556
# C  64002      0    872      0
# G 269223    668      0      0
# T    384      0      0      0

table(bim[which(bim$V5=="0"),"V1"])
# 16 19 24 
# 1  1  1 

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 3   6

write.table(indels[,"V2",drop=F],paste(path,"pre_imputation/QC/cedars_omni_old_gwas/list_indel_var_exclude_2",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg19_noind \
--exclude ${path_gwas}/pre_imputation/QC/cedars_omni_old_gwas/list_indel_var_exclude_2 \
--remove ${path_gwas}/pre_imputation/QC/cedars_omni_old_gwas/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg19_noind
# 723004 variants and 1226 people pass filters and QC.
# Among remaining phenotypes, 1222 are cases and 4 are controls.

#######################################################################################################################################
