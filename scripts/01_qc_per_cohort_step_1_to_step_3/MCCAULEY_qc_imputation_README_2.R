# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# text files downloaded to /path/to/project

cd /path/to/project
cd /path/to/user/preprocessed/
    
get mccauley.*
  
md5sum -c mccauley.txt.gz.md5
# mccauley.txt.gz: OK

zcat /path/to/project | head -10
# Sample ID	SNP Name	Chr	Position	REF	ALT	GT	X Raw	Y Raw	X	Y	GC Score
# sample_id	GSA-rs114420996	1	58814	G	A	0/0	576	3193	0.066	1.087	0.4790000
# sample_id	GSA-rs114420996	1	58814	G	A	0/0	342	2911	0.025	0.958	0.4790000

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/mccauley/


############################################################################## 
# make a temporary sorted fcr file:

# to keep consistency in naming
mv /path/to/project \
/path/to/project

# farm4:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh mccauley long gsa
# Job <193881> is submitted to queue <long>.


# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/sort_fcr* /path/to/user/scripts/IIBDGC/

zcat /path/to/ibdgwas/IIBDGC/pre_imputation/raw/mccauley/mccauley_gsa_sorted.txt.gz | wc -l
# 542169217

zcat /path/to/project | wc -l
# 542169217

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

# farm5:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_map_from_fcr.sh mccauley long gsa
#Job <293432> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_map_from_fcr* /path/to/user/scripts/IIBDGC/

  
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

# farm4:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_fam_from_fcr.sh mccauley normal gsa
#Job <293445> is submitted to queue <normal>.

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

# farm4:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_lgen_from_fcr.sh mccauley normal gsa
#Job <293449> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/
  

###########################

wc -l /path/to/ibdgwas/IIBDGC/pre_imputation/raw/mccauley/*.fam
# 788 /path/to/ibdgwas/IIBDGC/pre_imputation/raw/mccauley/mccauley_gsa.fam

wc -l /path/to/ibdgwas/IIBDGC/pre_imputation/raw/mccauley/*.map
# 688032 /path/to/ibdgwas/IIBDGC/pre_imputation/raw/mccauley/mccauley_gsa.map

zcat /path/to/project | awk 'END { print NR }'
# 542169217

# 542169217/688032
# [1] 788 OK

##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/mccauley_gsa/
mkdir ${path_gwas}pre_imputation/QC/mccauley_gsa/logs/

MEM=1000
path_gwas=/path/to/ibdgwas/IIBDGC/
bsub -J"plink0" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/mccauley_gsa/logs/stderr_plink_0_mccauley_gsa \
-o ${path_gwas}pre_imputation/QC/mccauley_gsa/logs/stdout_plink_0_mccauley_gsa \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/mccauley_old/mccauley_gsa \
--missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19"

# 688032 variants and 788 people pass filters and QC.
# Note: No phenotypes present.



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

bim<-read.table(paste(path,"pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 116771     6

table(indels$V1)
# 

# X    X chromosome                    -> 23
# Y    Y chromosome                    -> 24
# XY   Pseudo-autosomal region of X    -> 25
# MT   Mitochondrial                   -> 26

mt<-bim[which(bim$V1 %in% c(0,26)),]
table(mt$V1)
# 26 
# 133 

chry<-bim[which(bim$V1 %in% c(24)),]
table(chry$V1)
# 24 
# 1473 

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 10510 10033  7884  6063  6435  6485  6364  5091  4675  4891  6495  5640  3099 
# 14    15    16    17    18    19    20    21    22    23    24    26 
# 3560  3914  4009  4938  2398  4308  2527  1459  1895  3278   775   133 

dim(all_remove)
#[1] 116859      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/mccauley_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# mccauley 
# 788

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#          Crohn's Disease Indeterminate Ulcerative Colitis
# mccauley             453            17                318


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#          Affected
# mccauley      788

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#            0
# mccauley 788

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                    Affected
# Crohn's Disease         453
# Indeterminate            17
# Ulcerative Colitis      318

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#          Female Male
# mccauley    384  404


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
# 2 
# 788

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# 1   2 
# 404 384 

fam_test<-read.table(paste(path,"pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 788

write.table(fam,paste(path,"pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 0  2

write.table(samples_remove,paste(path,"pre_imputation/QC/mccauley_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")


######

mkdir ${path_gwas}pre_imputation/QC/mccauley_gsa/
mkdir ${path_gwas}pre_imputation/QC/mccauley_gsa/logs/
  
path_gwas=/path/to/ibdgwas/IIBDGC/
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/mccauley_gsa/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/mccauley_gsa/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19_noind
# 571173 variants and 788 people pass filters and QC.
# Among remaining phenotypes, 788 are cases and 0 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37


