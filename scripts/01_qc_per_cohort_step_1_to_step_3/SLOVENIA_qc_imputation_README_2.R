# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
############
# slovenia #
############

# text files downloaded to /path/to/project

md5sum -c slovenia_gsa.txt.gz.md5
# slovenia_gsa.txt.gz: OK

zcat /path/to/project | head -2
# Sample ID	Sample Group	SNP Name	Chr	Position	Allele1 - Top	Allele2 - Top	X Raw	Y Raw	X	Y	GC Score	GT Score
# sample_id		1:100292476	1	100292476	A	A	7034	1118	0.745	0.288	0.2505	0.6560

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/slovenia/

############################################################################## 
# make a temporary sorted fcr file:

# farm4:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh slovenia 
# Job <1070109> is submitted to default queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/sort_fcr* /path/to/user/scripts/IIBDGC/
  
# many variants with chr = 0 and position = 0, investigate later and ask Phil

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
bash ${path_gwas}scripts/create_map_from_fcr.sh slovenia 
# Job <1074395> is submitted to default queue <normal>. # NEW

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

# farm4:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_fam_from_fcr.sh slovenia 
#Job <1099148> is submitted to default queue <normal>. # NEW SCRIPT

# keep updated back-up copy:
cp /path/to/ibdgwas/IIBDGC/scripts/create_fam_from_fcr* /path/to/user/scripts/IIBDGC/
  

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
bash ${path_gwas}scripts/create_lgen_from_fcr.sh slovenia 
#Job <1099162> is submitted to default queue <normal>.

# keep updated back-up copy:
cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/

  
############################
# remove sorted file:
rm /path/to/ibdgwas/IIBDGC/pre_imputation/raw/slovenia/slovenia_gsa_sorted.txt.gz
   
##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/slovenia_gsa/
mkdir ${path_gwas}pre_imputation/QC/slovenia_gsa/logs/

MEM=1000
path_gwas=/path/to/ibdgwas/IIBDGC/
bsub -J"plink0" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/slovenia_gsa/logs/stderr_plink_0_slovenia_gsa \
-o ${path_gwas}pre_imputation/QC/slovenia_gsa/logs/stdout_plink_0_slovenia_gsa \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/slovenia/slovenia_gsa \
--missing-genotype - \
--make-bed --out ${path_gwas}pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19"
# 700078 variants loaded from .bim file.
# 270 people (0 males, 0 females, 270 ambiguous) loaded from .fam.
# 700078 variants and 270 people pass filters and QC.

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

bim<-read.table(paste(path,"pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        -      A      C      D      G      I      T
# -   2530  48365  17477    969  68028   2535   1105
# A      0      0  56561      0 258475      0    846
# C      0  48008      0      0   1095      0      0
# D      0      0      0      0      0    289      0
# G      0 191855   1095      0      0      0      0
# I      0      0      0    148      0      0      0
# T      0    697      0      0      0      0      0

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 141446      6

table(indels$V1)

# 0     1     2     3     4     5     6     7     8     9    10    11    12 
# 6512 12009 11650  9143  6881  7273  7497  7257  5761  5477  5734  7687  6470 
# 13    14    15    16    17    18    19    20    21    22    23    24    25 
# 3731  4089  4471  4663  5781  2782  4859  2867  1698  2199  4045   708   129 
# 26 
# 73 


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
# 0     1     2     3     4     5     6     7     8     9    10    11    12 
# 7387 12009 11650  9143  6881  7273  7497  7257  5761  5477  5734  7687  6470 
# 13    14    15    16    17    18    19    20    21    22    23    24    25 
# 3731  4089  4471  4663  5781  2782  4859  2867  1698  2199  4045   708   129 
# 26 
# 140 

dim(all_remove)
#[1] 142388      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/slovenia_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# slovenia 
#      270 


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#              Crohn's Disease
# slovenia 181              89

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#          Affected Unaffected
# slovenia       89        181

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#            0   1
# slovenia  89 181

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                 Affected Unaffected
#                        0        181
# Crohn's Disease       89          0

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#              Female Male
# slovenia   1    210   59


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
# 181  89   

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# -9   1   2 
#  1  59 210   

fam_test<-read.table(paste(path,"pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 270

write.table(fam,paste(path,"pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 0   2

write.table(samples_remove,paste(path,"pre_imputation/QC/slovenia_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")


######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/slovenia_gsa/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/slovenia_gsa/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19_noind
# 270 people (59 males, 210 females, 1 ambiguous) loaded from .fam.
# Among remaining phenotypes, 89 are cases and 181 are controls.
# 556469 variants and 270 people pass filters and QC.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37

