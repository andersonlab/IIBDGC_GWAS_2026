# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# text files downloaded to /path/to/project

md5sum -c niddk_broad_gsa.txt.gz.md5


zcat /path/to/project | head -10
# Sample ID	Sample Group	SNP Name	Chr	Position	Allele1 - Top	Allele2 - Top	X Raw	Y Raw	X	Y	GC Score	GT Score
# sample_id		1:100292476	1	100292476	A	A	2155	485	0.241	0.103	0.2593	0.6671
# sample_id		1:101064936	1	101064936	A	A	25059	3019	2.930	0.615	0.2559	0.6628

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/niddk_broad/
  
  
############################################################################## 
# make a temporary sorted fcr file:

# farm4:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh niddk_broad basement gsa
# Job <7389490> is submitted to queue <basement>.

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

# farm4:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_map_from_fcr.sh niddk_broad long gsa 
# Job <7631579> is submitted to default queue <normal>.

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
#Job <7631641> is submitted to default queue <normal>.

bash ${path_gwas}scripts/create_fam_from_fcr.sh niddk_broad normal gsa 
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
  bash ${path_gwas}scripts/create_lgen_from_fcr.sh niddk_broad normal gsa 
#Job <7631663> is submitted to default queue <normal>.

# keep updated back-up copy:
cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/
  
  
############################
# remove sorted file:

rm /path/to/ibdgwas/IIBDGC/pre_imputation/raw/niddk_broad/niddk_broad_gsa_sorted.txt.gz


##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/niddk_broad_gsa/
mkdir ${path_gwas}pre_imputation/QC/niddk_broad_gsa/logs/
MEM=18000

bsub -J"bcft_gwas3" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group \
-e ${path_gwas}pre_imputation/QC/niddk_broad_gsa/logs/stderr_lfile_to_bed \
-o ${path_gwas}pre_imputation/QC/niddk_broad_gsa/logs/stdout_lfile_to_bed \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/niddk_broad/niddk_broad_gsa \
--missing-genotype - \
--make-bed --out ${path_gwas}pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19"
# Job <703650> is submitted to default queue <normal>.

# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.



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

bim<-read.table(paste(path,"pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        -      A      C      D      G      I      T
# -   3055  23148  10107    911  35627   2238    791
# A      0      0  63372      0 290336      0   1152
# C      0  53120      0      0   1686      0      0
# D      0      0      0      0      0    585      0
# G      0 211159   1624      0      0      0      0
# I      0      0      0    208      0      0      0
# T      0    959      0      0      0      0      0

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 76670      6

table(indels$V1)

# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 6252 6235 6300 4674 3099 3683 3599 3795 2666 2771 2806 4073 3335 1953 2142 2502 
# 16   17   18   19   20   21   22   23   24   25   26 
# 2488 3394 1293 2693 1415  895 1155 2617  792   29   14 


# X    X chromosome                    -> 23
# Y    Y chromosome                    -> 24
# XY   Pseudo-autosomal region of X    -> 25
# MT   Mitochondrial                   -> 26

mt<-bim[which(bim$V1 %in% c(0,26)),]
table(mt$V1)
# 0   26 
# 7387  140  

chry<-bim[which(bim$V1 %in% c(24)),]
table(chry$V1)
# 24 
# 1480

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 7387 6235 6300 4674 3099 3683 3599 3795 2666 2771 2806 4073 3335 1953 2142 2502 
# 16   17   18   19   20   21   22   23   24   25   26 
# 2488 3394 1293 2693 1415  895 1155 2617  792   29  140 


dim(all_remove)
#[1] 77931      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/niddk_broad_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# niddk_broad 
# 5531 


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#                   Crohn's Disease Indeterminate Ulcerative Colitis
# niddk_broad  1043            2747             7               1734

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                  Affected Unaffected
# niddk_broad   16     4488       1027

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#                0    1
# niddk_broad 4488 1027

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                         Affected Unaffected
#                      16        0       1027
# Crohn's Disease       0     2747          0
# Indeterminate         0        7          0
# Ulcerative Colitis    0     1734          0


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#                  Female Male
# niddk_broad   16   2657 2858


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
# -9    1    2 
# 16 1027 4488

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# -9    1    2 
# 16 2858 2657   

fam_test<-read.table(paste(path,"pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 5531

write.table(fam,paste(path,"pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 16  2

write.table(samples_remove,paste(path,"pre_imputation/QC/niddk_broad_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")


######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/niddk_broad_gsa/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/niddk_broad_gsa/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19_noind
# 622147 variants and 5515 people pass filters and QC.
# Among remaining phenotypes, 4488 are cases and 1027 are controls.


#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37


