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
# cedars 370k #
###############

md5sum -c cedars_370k_newids.tgz.md5
# md5sum -c cedars_370k_newids.tgz.md5: OK

gunzip cedars_370k_newids.tgz
tar -xvf cedars_370k_newids.tar


# -rw-r--r--  1 username team152  71437980 Mar 28  2011 cedars_370k_newids.bed
# -rw-r--r--  1 username team152   9560101 Mar 28  2011 cedars_370k_newids.bim
# -rw-r--r--  1 username team152     26432 Mar 28  2011 cedars_370k_newids.fam

  
## explore files:
wc -l cedars_370k_newids.bim
# 345111 cedars_370k_newids.bim

head -5000 cedars_370k_newids.bim | tail -10
# 1	rs11210896	0	43875862	A	G
# 1	rs669446	0	43902366	A	G

# in hg18

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/cedars_370k_old_gwas/
  
cp cedars_370k_newids.* /path/to/ibdgwas/IIBDGC/pre_imputation/raw/cedars_370k_old_gwas/
  
path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/cedars_370k_old_gwas/cedars_370k_newids \
--freq \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18
# 345111 variants and 826 people pass filters and QC.
# Among remaining phenotypes, 826 are cases and 0 are controls.


# check phenotyped data within the mds folder:
cd  /path/to/project
  
### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

pheno_ibd<-read.table(paste(path,"raw_data_from_data_commons/UK_IIBDGC/ibdgc_other_cohorts/gwas/mds/IBD15b.mepr.mds_cov",sep=""),head=T)

dim(fam[which(fam$V1 %in% pheno_ibd$IID),])
# [1] 0 6

pheno_ibd$cohort<-gsub("\\*.*","",pheno_ibd$FID)
table(pheno_ibd$cohort)
# case_cd_bel1_eur_I317    case_cd_bel2_eur_I317    case_cd_cedr_eur_I317 
# 513                      153                      835 
# case_cd_chop_eur_I550    case_cd_germ_eur_I550    case_cd_nidk_eur_I317 
# 1482                      479                      756 
# case_cd_wtcc_eur_A500    case_uc_cedr_eur_I317    case_uc_chop_eur_I550 
# 1719                      835                      671 
# case_uc_germ_eur_A6.0    case_uc_nid3_eur_I550    case_uc_nid5_eur_I550 
# 988                      486                      447 
# case_uc_norw_eur_A6.0    case_uc_swed_eur_I317    case_uc_wtcc_eur_A6.0 
# 258                      918                     2342 
# control_cd_bel1_eur_I317 control_cd_bel2_eur_I317 control_cd_cedr_eur_I317 
# 884                       94                     1364 
# control_cd_chop_eur_I550 control_cd_germ_eur_I550 control_cd_nidk_eur_I317 
# 3054                      573                      462 
# control_cd_wtcc_eur_A500 control_uc_cedr_eur_I317 control_uc_chop_eur_I550 
# 1612                     1566                     3038 
# control_uc_germ_eur_A6.0 control_uc_nid3_eur_I550 control_uc_nid5_eur_I550 
# 2383                      624                     1420 
# control_uc_norw_eur_A6.0 control_uc_swed_eur_I317 control_uc_wtcc_eur_A6.0 
# 279                      341                     4076 
# nocc_cd_bel2_eur_I317    nocc_cd_cedr_eur_I317 
# 2                        5 


##

pheno_uc<-read.table("/path/to/project",head=T)

dim(fam[which(fam$V1 %in% pheno_uc$IID),])
# [1] 0 6

pheno_uc$cohort<-gsub("\\*.*","",pheno_uc$FID)
table(pheno_uc$cohort)
# case_uc_cedr_eur_I317    case_uc_chop_eur_I550    case_uc_germ_eur_A6.0 
# 836                      664                      990 
# case_uc_nid3_eur_I550    case_uc_nid5_eur_I550    case_uc_norw_eur_A6.0 
# 498                      451                      258 
# case_uc_swed_eur_I317    case_uc_wtcc_eur_A6.0 control_uc_cedr_eur_I317 
# 918                     2353                     2928 
# control_uc_chop_eur_I550 control_uc_germ_eur_A6.0 control_uc_nid3_eur_I550 
# 6091                     2915                     1070 
# control_uc_nid5_eur_I550 control_uc_norw_eur_A6.0 control_uc_swed_eur_I317 
# 1428                      279                      341 
# control_uc_wtcc_eur_A6.0 
# 5412 




### check overlaps between cedars samples:

fam1<-read.table(paste(path,"pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18.fam",sep=""),head=F)
fam2<-read.table(paste(path,"pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg18.fam",sep=""),head=F)
fam3<-read.table(paste(path,"pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg18.fam",sep=""),head=F)

dim(fam1)
# [1] 826   6
dim(fam2)
# [1] 1033    6
dim(fam3)
# [1] 1228    6

dim(fam1[which(fam1$V1 %in% fam2$V1),])
# [1] 0 6
dim(fam1[which(fam1$V1 %in% fam3$V1),])
# [1] 0 6
dim(fam2[which(fam2$V1 %in% fam3$V1),])
# [1] 0 6


# create a list of IDs for these samples to send to Talin

fam1$study<-"cedars_370k"
fam2$study<-"cedars_610k"
fam3$study<-"cedars_omni"
all<-rbind(fam1,fam2,fam3)
table(all$V1==all$V2)
# TRUE 
# 3087
colnames(all)[1]<-"sample_id"
all<-all[,c(1,7)]

write.table(all,"~/files_iibdgc/list_cedars_old_gwas_ids.tsv",col.names=T,row.names=F,quote=F,sep="\t")

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

bim<-read.table(paste(path,"pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        0      A      C      G      T
# 0      3   1956    728   2179    142
# A      0      0  32752 146229    339
# C      0  30549      0    447      0
# G      0 129041    439      0      0
# T      0    307      0      0      0

table(bim$V1)
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 25597 27276 23283 20811 20694 24908 18549 19229 16633 16953 15959 16117 12573 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 10574  9670  9572  9083 10863  6574  8178  5667  5793 10338   150    67 

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 5008  6

table(indels$V1)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
# 405 286 218 197 198 626 301 170 173 212 226 206 122 100 165 170 200  72 121 140 
# 21  22  23  24  25 
# 66 103 439  84   8 

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
# 150 

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])


dim(all_remove)
# [1] 5008    6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/cedars_370k_old_gwas/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

# pheno file provided by Talin
pheno<-fread(paste(path,"pheno/list_cedars_old_gwas_ids_TH.txt",sep=""),head=T,sep="\t")
pheno<-as.data.frame(pheno)

table(pheno$Sex)
#     Female   Male 
# 356   1344   1387

table(pheno$IBD_Affection_status[which(pheno$sample_id %in% fam$V1)],pheno$Unrelated_control[which(pheno$sample_id %in% fam$V1)],useNA="ifany")
#              0 <NA>
#              0  217
# Affected     0  608
# Unaffected   1    0

table(pheno$IBD_Affection_status[which(pheno$sample_id %in% fam$V1)],pheno$Disease_type[which(pheno$sample_id %in% fam$V1)],useNA="ifany")
#                Crohn's Disease Indeterminate Ulcerative Colitis
#            217               0             0                  0
# Affected     0              68            21                519
# Unaffected   1               0             0                  0


fam<-read.table(paste(path,"pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18.fam",sep=""),head=F)
table(pheno$study[which(pheno$sample_id %in% fam$V1)])
# cedars_370k 
# 826 

################################################
# Sex (1=male; 2=female; other=unknown)

table(fam$V5,useNA="ifany")
# 0 
# 826 

fam$V5[which(fam$V1 %in% pheno$sample_id[which(pheno$Sex=="Male")])]<-1
fam$V5[which(fam$V1 %in% pheno$sample_id[which(pheno$Sex=="Female")])]<-2

table(fam$V5,useNA="ifany")
#   0   1   2 
# 217 325 284

################################################
# Affection status, by default, should be coded:
# -9 missing 
# 0 missing
# 1 unaffected
# 2 affected

table(fam$V6,useNA="ifany")
# 2 
# 826

fam[,"V6"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$Unrelated_control=="1")]),"V6"]<-"1"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$IBD_Affection_status=="Affected")]),"V6"]<-"2"

table(fam$V6)
# 2  -9 
# 608 218

write.table(fam,paste(path,"pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

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
# [1] 218   2

            
write.table(samples_remove,paste(path,"pre_imputation/QC/cedars_370k_old_gwas/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18.bed \
--bim ${path_gwas}/pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18.bim \
--fam ${path_gwas}/pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/cedars_370k_old_gwas/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/cedars_370k_old_gwas/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_noind
# 340103 variants and 608 people pass filters and QC.
# Among remaining phenotypes, 608 are cases and 0 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# NOTE: data in build 18, needs to be liftover first:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_noind \
--recode tab --out ${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_noind

# ${path_gwas}previous_qced_b38/liftover/liftOverx
# usage:
# liftOver oldFile map.chain newFile unMapped

#### R convert map into bed:

map<-read.table(paste(path,"pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_noind.map",sep=""),head=F)

map$chr<-paste("chr",map$V1,sep="")
map$chr[which(map$V1=="23")]<-"chrX"
map$chr[which(map$V1=="25")]<-"chrX"
map$chr[which(map$V1=="24")]<-"chrY"
map$chromStart<-map$V4
map$chromEnd<-map$V4+1

map$chromStart<-format(map$chromStart, scientific=F)
map$chromEnd<-format(map$chromEnd, scientific=F)


write.table(map[,c(5:7,2)],paste(path,"pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_noind.bed",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#### lift positions

${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_noind.bed \
${path_gwas}previous_qced_b38/liftover/hg18ToHg19.over.chain.gz \
${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_noind_lifted_hg19 \
${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_noind_no_lifted_hg19

cut -f 4 ${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_noind_no_lifted_hg19 | sed "/^#/d" \
> ${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/nonlifted_variants_to_exclude.dat

wc -l ${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/nonlifted_variants_to_exclude.dat
# 41

### exclude non lifted:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--file ${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_noind \
--exclude ${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/nonlifted_variants_to_exclude.dat \
--recode tab --out ${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_noind_liftedvariants
# 340062 variants and 608 people pass filters and QC.

## R

bed_lifted<-read.table(paste(path,"pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_noind_lifted_hg19",sep=""),head=F)
map<-read.table(paste(path,"pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_noind_liftedvariants.map",sep=""),head=F)

table(map$V2==bed_lifted$V4)
# TRUE 
# 340062 

map$pos<-bed_lifted$V2
table(map$V2==bed_lifted$V4)
table(map$pos==bed_lifted$V2)
# TRUE 
# 340062 

write.table(map[,c(1:3,5)],paste(path,"pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_noind_liftedvariants_edited.map",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


########

# put together updated .ped plus lifted.map
/path/to/software/plink_linux_x86_64_20181202/./plink \
--ped ${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_noind_liftedvariants.ped \
--map ${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_noind_liftedvariants_edited.map \
--make-bed --out ${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg19_noind
# 340062 variants and 608 people pass filters and QC.
# Among remaining phenotypes, 608 are cases and 0 are controls.

# remove intermediate files:

rm ${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg18_noind_*
  
  
# alleles are recoded, so now there are variants (monomorphic) with 0 alleles that need to be excluded:
  
#### R
  
bim<-read.table(paste(path,"pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg19_noind.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#         A      C      G      T
# 0     33     11     39      2
# A      0  32735 146166    336
# C  30546      0    445      0
# G 129007    436      0      0
# T    306      0      0      0

table(bim[which(bim$V5=="0"),"V1"])
# 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 22 23 
# 7  5  1  3  5 19  6  3  1  4  1  5  3  2  2  4  3  1  2  2  1  5 

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 85   6

write.table(indels[,"V2",drop=F],paste(path,"pre_imputation/QC/cedars_370k_old_gwas/list_indel_var_exclude_2",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg19_noind \
--exclude ${path_gwas}/pre_imputation/QC/cedars_370k_old_gwas/list_indel_var_exclude_2 \
--remove ${path_gwas}/pre_imputation/QC/cedars_370k_old_gwas/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg19_noind
# 339977 variants and 608 people pass filters and QC.
# Among remaining phenotypes, 608 are cases and 0 are controls.

#######################################################################################################################################
