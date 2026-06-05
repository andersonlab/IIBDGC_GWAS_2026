# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# text files downloaded to /path/to/project

cd /path/to/project
  
md5sum -c australia_2012_omniexome.txt.gz.md5
# australia_2012_omniexome.txt.gz: OK
md5sum -c australia_2014_omniexome.txt.gz.md5
# australia_2014_omniexome.txt.gz: OK

zcat /path/to/project | head -2
# Sample ID	Sample Group	SNP Name	Chr	Position	Allele1 - Top	Allele2 - Top	X Raw	Y Raw	GC Score	GT Score
# sample_id		200610-104	MT	212	A	A	21455	533	2.300	0.058	0.4097	0.8207

zcat /path/to/project | head -2
# Sample ID	Sample Group	SNP Name	Chr	Position	Allele1 - Top	Allele2 - Top	X Raw	Y Raw	GC Score	GT Score
# sample_id		exm-IND1-200449980	1	202183358	I	I	623	7971	0.026	0.471	0.5392	0.6815


zcat /path/to/project | awk 'END { print NR }'
# 886441045

zcat /path/to/project | awk 'END { print NR }'
# 364464955
# 364464955/964193
# [1] 378

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/australia_2012/
  
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/australia_2014/

############################################################################## 
# make a temporary sorted fcr file, updated file to indicate queue:

# farm3:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh australia_2012 basement omniexome
#Job <5759709> is submitted to queue <basement>.

bash ${path_gwas}scripts/sort_fcr.sh australia_2014 long omniexome
# Job <5754215> is submitted to queue <long>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/sort_fcr* /path/to/user/scripts/IIBDGC/

##############################################################################  
# generate map, fam and lgen files from this:

#############
# MAP FILES #
#############

# MAP file describes a single marker and must contain exactly 4 columns:
# chromosome (1-22, X, Y or 0 if unplaced)
# rs# or snp identifier
# Genetic distance (morgans)
# Base-pair position (bp units)

# farm3:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_map_from_fcr.sh australia_2014 normal omniexome
# Job <5757129> is submitted to default queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_map_from_fcr.sh australia_2012 long omniexome
# Job <6268982> is submitted to queue <long>.

# keep updated back-up copy:
cp /path/to/ibdgwas/IIBDGC/scripts/create_map_from_fcr* /path/to/user/scripts/IIBDGC/
  
#############
# FAM FILES #
#############

# FAM (6 first columns of PED) PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:
# Family ID
# Individual ID
# Paternal ID
# Maternal ID
# Sex (1=male; 2=female; other=unknown)
# Phenotype

# farm3:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_fam_from_fcr.sh australia_2014 normal omniexome
# Job <5757130> is submitted to default queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_fam_from_fcr.sh australia_2012 normal omniexome
# Job <6268991> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_fam_from_fcr* /path/to/user/scripts/IIBDGC/


##############
# LGEN FILES #
##############

# LGEN file, test.lgen
# family ID
# individual ID
# snp ID
# allele 1 of this genotype
# allele 2 of this genotype

# farm3:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_lgen_from_fcr.sh australia_2014 normal omniexome
#Job <5757135> is submitted to default queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_lgen_from_fcr.sh australia_2012 normal omniexome
# Job <6268996> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/


##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/australia_omniexome/
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/australia_2012/australia_2012_omniexome \
--missing-genotype - \
--output-missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/australia_omniexome/australia_2012_omniexome_hg19
# 951117 variants and 932 people pass filters and QC.


/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/australia_2014/australia_2014_omniexome \
--missing-genotype - \
--output-missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/australia_omniexome/australia_2014_omniexome_hg19
# 964193 variants and 378 people pass filters and QC.


##############################################
# /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

a<-read.table(paste(path,"pre_imputation/QC/australia_omniexome/australia_2012_omniexome_hg19.bim",sep=""),head=F)
b<-read.table(paste(path,"pre_imputation/QC/australia_omniexome/australia_2014_omniexome_hg19.bim",sep=""),head=F)

dim(a)
# [1] 951117      6
dim(b)
# [1] 964193      6

dim(a[which(a$V2 %in% b$V2),])
# [1] 932068      6

# not same set of IDs

a[which(a$V2=="rs1044145"),]
# V1        V2 V3 V4 V5 V6  XX
# 186  0 rs1044145  0  0  A  G 0_0
b[which(b$V2=="rs1044145"),]
# V1        V2 V3        V4 V5 V6          XX
# 68404  1 rs1044145  0 207217359  A  G 1_207217359

##########################

# MERGE THESE TWO SETS:

# using merge mode 3
# The default (mode 1) behaviour is to call the merged genotype as missing if the original and new files contain different, non-missing calls; otherwise: i.e.
# Merge mode
# data1.ped ,  data2.ped  ->  1    2    3    4    5    
# ---------    ---------      -----------------------
# 0/0      ,   0/0       ->  0/0  0/0  0/0  0/0  0/0
# 0/0      ,   A/A       ->  A/A  A/A  A/A  0/0  A/A
# A/A      ,   0/0       ->  A/A  A/A  A/A  A/A  0/0
# A/A      ,   A/T       ->  0/0  A/A  A/T  A/A  A/T

# although there should not be any need to overwrite, no sample overlap


/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/australia_2014_omniexome_hg19 \
--bmerge /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/australia_2012_omniexome_hg19.bed \
/path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/australia_2012_omniexome_hg19.bim \
/path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/australia_2012_omniexome_hg19.fam \
--merge-mode 3 \
--make-bed --out ${path_gwas}pre_imputation/QC/australia_omniexome/australia_omniexome_hg19

# Warning: Multiple chromosomes seen for variant 'rs1044145'.
# Warning: Multiple chromosomes seen for variant 'rs10830214'.  # it will keep the one from the first file
# Warning: Multiple chromosomes seen for variant 'rs10843891'.

# Of these, 32125 are new, while 932068 are present in the base dataset.
# Warning: Variants 'rs3748597' and 'exm340' have the same position.
# Warning: Variants 'rs13303010' and 'exm2264981' have the same position.
# Warning: Variants 'exm596' and 'exm2253593' have the same position.
# 20599 more same-position warnings: see log file.
# Performing single-pass merge (1310 people, 983242 variants).
# 983242 variants and 1310 people pass filters and QC.

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/australia_omniexome/australia_omniexome_hg19 \
--freq \
--out ${path_gwas}pre_imputation/QC/australia_omniexome/australia_omniexome_hg19



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

bim<-read.table(paste(path,"pre_imputation/QC/australia_omniexome/australia_omniexome_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        0      A      C      D      G      I      T
# 0     10  50288  22735     15  98937     49   3422
# A      0      0  77303      0 360467      0   2637
# C      0  68070      0      0   5246      0      0
# D      0      0      0      0      0     49      0
# G      0 287143   4827      0      0      0      0
# I      0      0      0     28      0      0      0
# T      0   2016      0      0      0      0      0


indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 175533      6

table(indels$V1)
# 0     1     2     3     4     5     6     7     8     9    10    11    12 
# 377 17762 12740 10678  7215  8236  9459  8166  6531  7142  7027 11045  8938 
# 13    14    15    16    17    18    19    20    21    22    23    24    25 
# 3584  5438  5717  7363  9252  3137 10211  4744  1961  3548  4509   574   110 
# 26 
# 69 

# X    X chromosome                    -> 23
# Y    Y chromosome                    -> 24
# XY   Pseudo-autosomal region of X    -> 25
# MT   Mitochondrial                   -> 26

mt<-bim[which(bim$V1 %in% c(0,26)),]
table(mt$V1)
# 0  26 
# 825 222

chry<-bim[which(bim$V1 %in% c(24)),]
table(chry$V1)
# 24 
# 1575

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])
#

dim(all_remove)
#[1] 176135      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/australia_omniexome/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

##########################################################
# KEEP ONLY THE CASES AND CONTROLS AS DONE IN PREVIOUS QC

pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/australia_omniexome/australia_omniexome_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# anz 
# 1310


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#         Ulcerative Colitis
# anz 617                693

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#     Affected Unaffected
# anz      693        617

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#       0   1
# anz 693 617

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                    Affected Unaffected
#                           0        617
# Ulcerative Colitis      693          0

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#     Female Male
# anz    706  604


################################################
# Affection status, by default, should be coded:
# -9 missing 
# 0 missing
# 1 unaffected
# 2 affected

fam[,"V6"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$control=="1")]),"V6"]<-"1"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$affection=="Affected")]),"V6"]<-"2"

table(fam$V6)
# 1   2 
# 617 693 

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# 1   2 
# 604 706   

fam_test<-read.table(paste(path,"pre_imputation/QC/australia_omniexome/australia_omniexome_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 1310

write.table(fam,paste(path,"pre_imputation/QC/australia_omniexome/australia_omniexome_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 0  2

write.table(samples_remove,paste(path,"pre_imputation/QC/australia_omniexome/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/australia_omniexome/australia_omniexome_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/australia_omniexome/australia_omniexome_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/australia_omniexome/australia_omniexome_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/australia_omniexome/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/australia_omniexome/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/australia_omniexome/australia_omniexome_hg19_noind
# 807107 variants and 1310 people pass filters and QC.
# Among remaining phenotypes, 693 are cases and 617 are controls.


#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37
