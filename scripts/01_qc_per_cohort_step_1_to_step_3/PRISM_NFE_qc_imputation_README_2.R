# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
#############
# prism_nfe #
#############

# text files downloaded to /path/to/project
  
md5sum -c prism_nfe.txt.gz.md5
# prism_nfe.txt.gz.md5: OK

zcat /path/to/project | head -2

# Sample ID	SNP Name	Chr	Position	Allele1 - Top	Allele2 - Top	X Raw	Y Raw	X	Y	GC Score
# sample_id	GSA-rs114420996	1	58814	G	A	462	2679	0.084	1.112	0.4790000



zcat /path/to/project | awk -F '\t' '{print $7}'| sort | uniq -c | sort -nr 


# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/prism_nfe/
  
############################################################################## 
# make a temporary sorted fcr file:

# to keep consistency in naming
mv /path/to/project \
/path/to/project

# farm4:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh prism_nfe long gsa
# Job <547508> is submitted to default queue <long>.


# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/sort_fcr* /path/to/user/scripts/IIBDGC/

zcat /path/to/project | wc -l
# 331631425

zcat /path/to/ibdgwas/IIBDGC/pre_imputation/raw/prism_nfe/prism_nfe_gsa_sorted.txt.gz | wc -l
# 331631425


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

# farm4:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_map_from_fcr.sh prism_nfe normal gsa
#Job <678688> is submitted to queue <normal>.

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
bash ${path_gwas}scripts/create_fam_from_fcr.sh prism_nfe normal gsa
#Job <678689> is submitted to queue <normal>.

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
bash ${path_gwas}scripts/create_lgen_from_fcr.sh prism_nfe normal gsa
#Job <678690> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/
  

wc -l /path/to/ibdgwas/IIBDGC/pre_imputation/raw/prism_nfe/*.fam
# 482 /path/to/ibdgwas/IIBDGC/pre_imputation/raw/prism_nfe/prism_nfe_gsa.fam

wc -l /path/to/ibdgwas/IIBDGC/pre_imputation/raw/prism_nfe/*.map
# 688032 /path/to/ibdgwas/IIBDGC/pre_imputation/raw/prism_nfe/prism_nfe_gsa.map

zcat /path/to/project | awk 'END { print NR }'
# 331631425

# 331631425/688032
# [1] 482 OK

############################
# remove sorted file:
# rm /path/to/ibdgwas/IIBDGC/pre_imputation/raw/prism_nfe_gsa/prism_nfe_gsa_illugwas_sorted.txt.gz

##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/prism_nfe_gsa/
mkdir ${path_gwas}pre_imputation/QC/prism_nfe_gsa/logs/
  
MEM=1000
path_gwas=/path/to/ibdgwas/IIBDGC/
bsub -J"plink0" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/prism_nfe_gsa/logs/stderr_plink_0_prism_nfe_gsa \
-o ${path_gwas}pre_imputation/QC/prism_nfe_gsa/logs/stdout_plink_0_prism_nfe_gsa \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/prism_nfe/prism_nfe_gsa \
--missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19"

# 688032 variants and 482 people pass filters and QC.
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

bim<-read.table(paste(path,"pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 114185     6

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
# 

dim(all_remove)
#[1] 114272      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/prism_nfe_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# prism_nfe 
# 482


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#               Crohn's Disease Ulcerative Colitis Unknown
# prism_nfe  49             278            17                138

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#               Affected Unaffected Unknown
# prism_nfe  13      433         33       3

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#             0   1
# prism_nfe 433  33

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                        Affected Unaffected Unknown
#                     13        0         33       3
# Crohn's Disease      0      278          0       0
# Indeterminate        0       17          0       0
# Ulcerative Colitis   0      138          0       0

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#               Female Male
# prism_nfe  13    251  218


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
# 1   2  -9 
# 33 433  16  

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# 1   2  -9 
# 218 251  13 

fam_test<-read.table(paste(path,"pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 482

write.table(fam,paste(path,"pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 16  2

write.table(samples_remove,paste(path,"pre_imputation/QC/prism_nfe_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")


######

mkdir ${path_gwas}pre_imputation/QC/prism_nfe_gsa/
mkdir ${path_gwas}pre_imputation/QC/prism_nfe_gsa/logs/
  
path_gwas=/path/to/ibdgwas/IIBDGC/
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/prism_nfe_gsa/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/prism_nfe_gsa/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19_noind
# 573760 variants and 466 people pass filters and QC.
# Among remaining phenotypes, 433 are cases and 33 are controls.




#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37


