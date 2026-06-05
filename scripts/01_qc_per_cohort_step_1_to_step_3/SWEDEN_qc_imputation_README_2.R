# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
##########
# sweden #
##########

# text files downloaded to /path/to/project

# TWO SEPARATED FILES, sweden_controls sweden_sicbio

md5sum -c sweden_controls_gsa.txt.gz.md5
# sweden_controls_gsa.txt.gz: OK

md5sum -c sweden_sicbio_gsa.txt.gz.md5
# sweden_sicbio_gsa.txt.gz: OK

zcat /path/to/project | head -2
# Sample ID	Sample Group	SNP Name	Chr	Position	Allele1 - Top	Allele2 - Top	X Raw	Y Raw	X	Y	GC Score	GT Score
# sample_id		1:100292476	1	100292476	A	A	6879	1768	0.693	0.309	0.2505	0.6560

zcat /path/to/project | head -2
# Sample ID	Sample Group	SNP Name	Chr	Position	Allele1 - Top	Allele2 - Top	X Raw	Y Raw	X	Y	GC Score	GT Score
# sample_id		1:100292476	1	100292476	A	A	9634	1761	0.724	0.289	0.2505	0.6560


# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/sweden_sicbio/
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/sweden_controls/
  
############################################################################## 
# make a temporary sorted fcr file:

# farm4:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh sicbio_gsa
# Job <1070113> is submitted to default queue <normal>.

bash ${path_gwas}scripts/sort_fcr.sh sweden_controls 
# Job <1074042> is submitted to default queue <long>.

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
bash ${path_gwas}scripts/create_map_from_fcr.sh sweden_sicbio 
# Job <1074806> is submitted to default queue <normal>.
# bash ${path_gwas}scripts/create_map_from_fcr.sh sweden_controls 

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_map_from_fcr.sh sweden_controls 
# Job <1098598> is submitted to default queue <normal>.

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
bash ${path_gwas}scripts/create_fam_from_fcr.sh sweden_sicbio 
# Job <1099165> is submitted to default queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_fam_from_fcr.sh sweden_controls 
# Job <1099166> is submitted to default queue <normal>.

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
bash ${path_gwas}scripts/create_lgen_from_fcr.sh sweden_sicbio 
#Job <1074809> is submitted to default queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_lgen_from_fcr.sh sweden_controls 
# Job <1098600> is submitted to default queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/
  

############################
# remove sorted file:
rm /path/to/ibdgwas/IIBDGC/pre_imputation/raw/sweden_controls/sweden_controls_gsa_sorted.txt.gz

##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/sweden_gsa/
mkdir ${path_gwas}pre_imputation/QC/sweden_gsa/logs/
 
MEM=1000
path_gwas=/path/to/ibdgwas/IIBDGC/
bsub -J"plink0" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/sweden_gsa/logs/stderr_plink_0_sweden_controls_gsa \
-o ${path_gwas}pre_imputation/QC/sweden_gsa/logs/stdout_plink_0_sweden_controls_gsa \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/sweden_controls/sweden_controls_gsa \
--missing-genotype - \
--output-missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/sweden_gsa/sweden_controls_gsa_hg19"
# 700078 variants and 949 people pass filters and QC.


MEM=1000
path_gwas=/path/to/ibdgwas/IIBDGC/
  bsub -J"plink0" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/sweden_gsa/logs/stderr_plink_0_sweden_sicbio_gsa \
-o ${path_gwas}pre_imputation/QC/sweden_gsa/logs/stdout_plink_0_sweden_sicbio_gsa \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/sweden_sicbio/sweden_sicbio_gsa \
--missing-genotype - \
--output-missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/sweden_gsa/sweden_sicbio_gsa_hg19"
# 700078 variants and 574 people pass filters and QC.

# MERGE THESE TWO SETS:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/sweden_controls_gsa_hg19 \
--bmerge /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/sweden_sicbio_gsa_hg19.bed \
/path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/sweden_sicbio_gsa_hg19.bim \
/path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/sweden_sicbio_gsa_hg19.fam \
--merge-mode 3 \
--make-bed --out ${path_gwas}pre_imputation/QC/sweden_gsa/sweden_gsa_hg19
# Warning: Variants 'rs121434297' and '1:11855218' have the same position.
# 1600 more same-position warnings: see log file.
# 700078 variants and 1523 people pass filters and QC.

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/sweden_gsa/sweden_gsa_hg19 \
--freq \
--out ${path_gwas}pre_imputation/QC/sweden_gsa/sweden_gsa_hg19



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

bim<-read.table(paste(path,"pre_imputation/QC/sweden_gsa/sweden_gsa_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        0      A      C      D      G      I      T
# 0   2528  40293  15920   1034  60267   2626   1035
# A      0      0  58037      0 266298      0    912
# C      0  49524      0      0   1237      0      0
# D      0      0      0      0      0    196      0
# G      0 198105   1204      0      0      0      0
# I      0      0      0     85      0      0      0
# T      0    777      0      0      0      0      0

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 123984      6

table(indels$V1)
# 0     1     2     3     4     5     6     7     8     9    10    11    12 
# 6587 10488 10198  7950  5851  6405  6453  6375  4852  4687  4875  6726  5655 
# 13    14    15    16    17    18    19    20    21    22    23    24    25 
# 3266  3556  3924  4052  5210  2287  4392  2431  1463  1921  3588   668    86 
# 26 
# 38 

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
# 

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])
# 0     1     2     3     4     5     6     7     8     9    10    11    12 
# 7387 10488 10198  7950  5851  6405  6453  6375  4852  4687  4875  6726  5655 
# 13    14    15    16    17    18    19    20    21    22    23    24    25 
# 3266  3556  3924  4052  5210  2287  4392  2431  1463  1921  3588   668    86 
# 26 
# 140 

dim(all_remove)
#[1] 124886      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/sweden_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/sweden_gsa/sweden_gsa_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# sweden_controls   sweden_sicbio 
# 949             574 


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#                     Crohn's Disease Indeterminate Ulcerative Colitis
# sweden_controls 949               0             0                  0
# sweden_sicbio   157             192             6                219

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                     Affected Unaffected
# sweden_controls        0        949       0
# sweden_sicbio        417         48     109

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#                   0   1
# sweden_controls   0 949
# sweden_sicbio   526  48

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                        Affected Unaffected
  #                           0        997     109
  # Crohn's Disease         192          0       0
  # Indeterminate             6          0       0
  # Ulcerative Colitis      219          0       0

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#                 Female Male Unknown
# sweden_controls    374  575       0
# sweden_sicbio      286  286       2

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
# 997 417 109    


################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# 1   2  -9 
# 861 660   2 

fam_test<-read.table(paste(path,"pre_imputation/QC/sweden_gsa/sweden_gsa_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 1523

write.table(fam,paste(path,"pre_imputation/QC/sweden_gsa/sweden_gsa_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 109   2

write.table(samples_remove,paste(path,"pre_imputation/QC/sweden_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/sweden_gsa/sweden_gsa_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/sweden_gsa/sweden_gsa_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/sweden_gsa/sweden_gsa_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/sweden_gsa/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/sweden_gsa/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/sweden_gsa/sweden_gsa_hg19_noind
# 575192 variants and 1414 people pass filters and QC.
# Among remaining phenotypes, 417 are cases and 997 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37


