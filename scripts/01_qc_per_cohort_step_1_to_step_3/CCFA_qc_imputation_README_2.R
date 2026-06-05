# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# text files downloaded to /path/to/project

cd /path/to/project
cd /path/to/user/preprocessed/
  
get ccfa.*
  
md5sum -c ccfa.txt.gz.md5
# ccfa.txt.gz: OK

zcat /path/to/project | head -10
# Sample ID	SNP Name	Chr	Position	REF	ALT	GT	X Raw	Y Raw	X	Y	GC Score
# sample_id	GSA-rs114420996	1	58814	G	A	0/0	750	3711	0.037	0.976	0.4790000
# sample_id	GSA-rs114420996	1	58814	G	A	0/0	1165	3466	0.07	0.854	0.4790000
# sample_id	GSA-rs114420996	1	58814	G	A	0/0	975	3842	0.045	0.992	0.4790000

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/ccfa/
  
  
############################################################################## 
# make a temporary sorted fcr file:

# to keep consistency in naming
mv /path/to/project \
/path/to/project

# farm5:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh ccfa basement gsa
# Job <62129> is submitted to queue <basement>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_map_from_fcr* /path/to/user/scripts/IIBDGC/

zcat /path/to/ibdgwas/IIBDGC/pre_imputation/raw/ccfa/ccfa_gsa_sorted.txt.gz | wc -l
# 1505414017

zcat /path/to/project | wc -l
# 1505414017


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
bash ${path_gwas}scripts/create_map_from_fcr.sh ccfa long gsa
#Job <196016> is submitted to queue <long>.

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
bash ${path_gwas}scripts/create_fam_from_fcr.sh ccfa normal gsa
#Job <196094> is submitted to queue <normal>.

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
bash ${path_gwas}scripts/create_lgen_from_fcr.sh ccfa normal gsa
# job <196098> is submitted to queue <normal>.


# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/
  
###########################

wc -l /path/to/ibdgwas/IIBDGC/pre_imputation/raw/ccfa/*.fam
# 2188 /path/to/ibdgwas/IIBDGC/pre_imputation/raw/ccfa/ccfa_gsa.fam

wc -l /path/to/ibdgwas/IIBDGC/pre_imputation/raw/ccfa/*.map
# 688032 /path/to/ibdgwas/IIBDGC/pre_imputation/raw/ccfa/ccfa_gsa.map

zcat /path/to/project | awk 'END { print NR }'
# 1505414017

# 1505414017/688032
# [1] 2188 OK

##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/ccfa_gsa/
mkdir ${path_gwas}pre_imputation/QC/ccfa_gsa/logs/
  
MEM=1000
path_gwas=/path/to/ibdgwas/IIBDGC/
bsub -J"plink0" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/ccfa_gsa/logs/stderr_plink_0_ccfa_gsa \
-o ${path_gwas}pre_imputation/QC/ccfa_gsa/logs/stdout_plink_0_ccfa_gsa \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/ccfa/ccfa_gsa \
--missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19"

# 688032 variants and 2188 people pass filters and QC.
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

bim<-read.table(paste(path,"pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 76562     6

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
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 6814 6620 5081 3736 4155 4158 4114 3167 3032 3170 4318 3646 2027 2306 2637 2712 
# 17   18   19   20   21   22   23   24   26 
# 3419 1464 2852 1619  995 1276 2507  720  133 

dim(all_remove)
#[1] 76678      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/ccfa_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# ccfa 
# 2188

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#      Crohn's Disease Indeterminate Ulcerative Colitis
# ccfa            1449            29                710

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#      Affected
# ccfa     2188

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#         0
# ccfa 2188

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                    Affected
# Crohn's Disease        1449
# Indeterminate            29
# Ulcerative Colitis      710

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#        Female Male
#   ccfa   1183 1005


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
#    2
# 2188

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# 1    2 
# 1005 1183

fam_test<-read.table(paste(path,"pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 2188

write.table(fam,paste(path,"pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 0  2

write.table(samples_remove,paste(path,"pre_imputation/QC/ccfa_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")


######

mkdir ${path_gwas}pre_imputation/QC/ccfa_gsa/
mkdir ${path_gwas}pre_imputation/QC/ccfa_gsa/logs/
  
path_gwas=/path/to/ibdgwas/IIBDGC/
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/ccfa_gsa/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/ccfa_gsa/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19_noind
# 611354 variants and 2188 people pass filters and QC.
# Among remaining phenotypes, 2188 are cases and 0 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37


