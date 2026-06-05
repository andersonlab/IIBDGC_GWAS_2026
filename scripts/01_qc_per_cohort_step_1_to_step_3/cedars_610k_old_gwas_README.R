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
# cedars 610k #
###############

md5sum -c cedars_610k_newids.tgz.md5
# md5sum -c cedars_610k_newids.tgz.md5: OK

gunzip cedars_610k_newids.tgz
tar -xvf cedars_610k_newids.tar


# cedars_610k_newids.bed
# cedars_610k_newids.bim
# cedars_610k_newids.fam


## explore files:
wc -l cedars_610k_newids.bim
# 592532 cedars_610k_newids.bim

head -5000 cedars_610k_newids.bim | tail -10
# 1	rs10903034	0	24353492	G	A
# 1	rs12028945	0	24354591	A	G

# in hg18

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/cedars_610k_old_gwas/
  
cp cedars_610k_newids.* /path/to/ibdgwas/IIBDGC/pre_imputation/raw/cedars_610k_old_gwas/
  
path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/cedars_610k_old_gwas/cedars_610k_newids \
--freq \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18
# 592532 variants and 1033 people pass filters and QC.
# Among remaining phenotypes, 1033 are cases and 0 are controls.
# 1033 people (0 males, 0 females, 1033 ambiguous) loaded from .fam.


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

bim<-read.table(paste(path,"pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        0      A      C      G      T
# 0     45   2132    808   2439    138
# A     35      0  57004 254003    655
# C      8  52103      0    995      0
# G     35 220599    936      0      0
# T      3    594      0      0      0

table(bim$V1)
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 45183 47741 39707 35541 35961 40724 32265 32618 27721 30725 28679 28528 22089 
# 14    15    16    17    18    19    20    21    22    23    24    25    26 
# 19213 17581 17596 15379 17337 10307 14686  8535  8883 15157   158    80   138 

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 5643  6

table(indels$V1)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
# 484 338 259 220 235 681 358 177 204 235 264 242 161 117 178 196 218  82 138 152 
# 21  22  23  24  25  26 
# 67 114 431  81   7   4 

# X    X chromosome                    -> 23
# Y    Y chromosome                    -> 24
# XY   Pseudo-autosomal region of X    -> 25
# MT   Mitochondrial                   -> 26

mt<-bim[which(bim$V1 %in% c(0,26)),]
table(mt$V1)
# 26 
# 138 


chry<-bim[which(bim$V1 %in% c(24)),]
table(chry$V1)
# 24 
# 158 

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])


dim(all_remove)
# [1] 5777    6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/cedars_610k_old_gwas/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

# pheno file provided by Talin
pheno<-fread(paste(path,"pheno/list_cedars_old_gwas_ids_TH.txt",sep=""),head=T,sep="\t")
pheno<-as.data.frame(pheno)

fam<-read.table(paste(path,"pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18.fam",sep=""),head=F)
table(pheno$study[which(pheno$sample_id %in% fam$V1)])
# cedars_610k 
# 826 

table(pheno$Sex[which(pheno$sample_id %in% fam$V1)])
#     Female   Male 
# 138    429    466 

table(pheno$IBD_Affection_status[which(pheno$sample_id %in% fam$V1)],pheno$Unrelated_control[which(pheno$sample_id %in% fam$V1)],useNA="ifany")
#              0 <NA>
#              0  138
# Affected     0  893
# Unaffected   2    0

table(pheno$IBD_Affection_status[which(pheno$sample_id %in% fam$V1)],pheno$Disease_type[which(pheno$sample_id %in% fam$V1)],useNA="ifany")
#                Crohn's Disease Indeterminate Ulcerative Colitis
#            138               0             0                  0
# Affected     0             823             2                 68
# Unaffected   2               0             0                  0



################################################
# Sex (1=male; 2=female; other=unknown)

table(fam$V5,useNA="ifany")
# 0 
# 1033 

fam$V5[which(fam$V1 %in% pheno$sample_id[which(pheno$Sex=="Male")])]<-1
fam$V5[which(fam$V1 %in% pheno$sample_id[which(pheno$Sex=="Female")])]<-2

table(fam$V5,useNA="ifany")
#   0   1   2 
# 138 466 429 

################################################
# Affection status, by default, should be coded:
# -9 missing 
# 0 missing
# 1 unaffected
# 2 affected

table(fam$V6,useNA="ifany")
# 2 
# 1033

fam[,"V6"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$Unrelated_control=="1")]),"V6"]<-"1"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$IBD_Affection_status=="Affected")]),"V6"]<-"2"

table(fam$V6)
# 2  -9 
# 893 140 

write.table(fam,paste(path,"pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

# exclude also from analysis sample sample_id 
pheno[which(pheno$sample_id=="sample_id"),]
#        sample_id       study sample_comment Genetic_ID Genetic_ID_2 Sex
# 1480 sample_id cedars_omni     withdrawn?    sample_id                 
#      IBD_Affection_status Disease_type Unrelated_control self_race
# 1480                                                  NA          
#         self_hispanic self_jewish
# 1480     

samples_remove<-fam[which((fam$V6=="-9") | (fam$V1=="sample_id")),c(1,2)]
dim(samples_remove)
# [1] 140   2


write.table(samples_remove,paste(path,"pre_imputation/QC/cedars_610k_old_gwas/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18.bed \
--bim ${path_gwas}/pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18.bim \
--fam ${path_gwas}/pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/cedars_610k_old_gwas/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/cedars_610k_old_gwas/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_noind
# 586755 variants and 893 people pass filters and QC.
# Among remaining phenotypes, 893 are cases and 0 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# NOTE: data in build 18, needs to be liftover first:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_noind \
--recode tab --out ${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_noind

# ${path_gwas}previous_qced_b38/liftover/liftOverx
# usage:
# liftOver oldFile map.chain newFile unMapped

#### R convert map into bed:

map<-read.table(paste(path,"pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_noind.map",sep=""),head=F)

map$chr<-paste("chr",map$V1,sep="")
map$chr[which(map$V1=="23")]<-"chrX"
map$chr[which(map$V1=="25")]<-"chrX"
map$chr[which(map$V1=="24")]<-"chrY"
map$chromStart<-map$V4
map$chromEnd<-map$V4+1

map$chromStart<-format(map$chromStart, scientific=F)
map$chromEnd<-format(map$chromEnd, scientific=F)


write.table(map[,c(5:7,2)],paste(path,"pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_noind.bed",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#### lift positions

${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_noind.bed \
${path_gwas}previous_qced_b38/liftover/hg18ToHg19.over.chain.gz \
${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_noind_lifted_hg19 \
${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_noind_no_lifted_hg19

cut -f 4 ${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_noind_no_lifted_hg19 | sed "/^#/d" \
> ${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/nonlifted_variants_to_exclude.dat

wc -l ${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/nonlifted_variants_to_exclude.dat
# 82

### exclude non lifted:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--file ${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_noind \
--exclude ${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/nonlifted_variants_to_exclude.dat \
--recode tab --out ${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_noind_liftedvariants
# 586673 variants and 893 people pass filters and QC.


## R

bed_lifted<-read.table(paste(path,"pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_noind_lifted_hg19",sep=""),head=F)
map<-read.table(paste(path,"pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_noind_liftedvariants.map",sep=""),head=F)

table(map$V2==bed_lifted$V4)
# TRUE 
# 586673 

map$pos<-bed_lifted$V2
table(map$V2==bed_lifted$V4)
table(map$pos==bed_lifted$V2)
# TRUE 
# 586673 

write.table(map[,c(1:3,5)],paste(path,"pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_noind_liftedvariants_edited.map",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


########

# put together updated .ped plus lifted.map
/path/to/software/plink_linux_x86_64_20181202/./plink \
--ped ${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_noind_liftedvariants.ped \
--map ${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_noind_liftedvariants_edited.map \
--make-bed --out ${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg19_noind
# 586673 variants and 893 people pass filters and QC.
# Among remaining phenotypes, 893 are cases and 0 are controls.

# remove intermediate files:

rm ${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18_noind_*
  
  
# alleles are recoded, so now there are variants (monomorphic) with 0 alleles that need to be excluded:
  
#### R
  
bim<-read.table(paste(path,"pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg19_noind.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        0      A      C      G      T
# 0      2    123     25    143      0
# A      0      0  56898 253767    650
# C      0  52148      0    987      0
# G      0 220403    932      0      0
# T      0    595      0      0      0

table(bim[which(bim$V5=="0"),"V1"])
# 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 
# 40 26 16 14 16 28 11 16  8 14  6 19  9  3  7 11 10  7  5  3  1  4 19 

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 293   6

write.table(indels[,"V2",drop=F],paste(path,"pre_imputation/QC/cedars_610k_old_gwas/list_indel_var_exclude_2",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg19_noind \
--exclude ${path_gwas}/pre_imputation/QC/cedars_610k_old_gwas/list_indel_var_exclude_2 \
--remove ${path_gwas}/pre_imputation/QC/cedars_610k_old_gwas/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg19_noind
# 586380 variants and 893 people pass filters and QC.
# Among remaining phenotypes, 893 are cases and 0 are controls.

#######################################################################################################################################
