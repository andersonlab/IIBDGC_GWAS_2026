# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
###############
# netherlands #
###############

# text files downloaded to /path/to/project

# TWO SEPARATED FILES, netherlands_controls netherlands_cases

md5sum -c netherlands_controls_gsa.txt.gz.md5
# netherlands_controls_gsa.txt.gz: OK

md5sum -c netherlands_cases_gsa.txt.gz.md5
# netherlands_cases_gsa.txt.gz: OK


zcat /path/to/project | head -2
# Sample ID	Sample Group	SNP Name	Chr	Position	Allele1 - Top	Allele2 - Top	X RaY Raw	X	Y	GC Score	GT Score
# sample_id		1:100292476	1	100292476	A	A	7272	1022	0.645	0.200.2669	0.6764

zcat /path/to/project | head -2
# Sample ID	Sample Group	SNP Name	Chr	Position	Allele1 - Top	Allele2 - Top	X Raw	Y Raw	X	Y	GC Score	GT Score
# sample_id		1:100292476	1	100292476	A	A	5667	1077	0.659	0.163	0.2607	0.6688

zcat /path/to/project | awk 'END { print NR }'
# 2908824091
# 2908824091/700078
# [1] 4155

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/netherlands_controls/
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/netherlands_cases/
  
  
############################################################################## 
# make a temporary sorted fcr file:

# farm4:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh netherlands_controls 
# Job <1246828> is submitted to queue <long>. 

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh netherlands_cases basement gsa
# Job <6251974> is submitted to queue <basement>.

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
bash ${path_gwas}scripts/create_map_from_fcr.sh netherlands_controls 
# Job <1282275> is submitted to default queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_map_from_fcr.sh netherlands_cases long gsa
# Job <6956532> is submitted to queue <long>.

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
bash ${path_gwas}scripts/create_fam_from_fcr.sh netherlands_controls 
# Job <1282276> is submitted to default queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_fam_from_fcr.sh netherlands_cases long gsa
# Job <6956528> is submitted to queue <long>.

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
bash ${path_gwas}scripts/create_lgen_from_fcr.sh netherlands_controls 
#Job <1282277> is submitted to default queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_lgen_from_fcr.sh netherlands_cases normal gsa
# Job <6956530> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/



############################
# remove sorted file:

rm /path/to/ibdgwas/IIBDGC/pre_imputation/raw/netherlands_controls/netherlands_controls_gsa_sorted.txt.gz
rm /path/to/ibdgwas/IIBDGC/pre_imputation/raw/netherlands_cases/netherlands_cases_gsa_sorted.txt.gz


##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/netherlands_gsa/
mkdir ${path_gwas}pre_imputation/QC/netherlands_gsa/logs
  

MEM=1000
path_gwas=/path/to/ibdgwas/IIBDGC/
bsub -J"plink0" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/netherlands_gsa/logs/stderr_plink_0_netherlands_controls_gsa \
-o ${path_gwas}pre_imputation/QC/netherlands_gsa/logs/stdout_plink_0_netherlands_controls_gsa \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/netherlands_controls/netherlands_controls_gsa \
--missing-genotype - \
--output-missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/netherlands_gsa/netherlands_controls_gsa_hg19"
#700078 variants and 550 people pass filters and QC.


bsub -J"plink0" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/netherlands_gsa/logs/stderr_plink_0_netherlands_cases_gsa \
-o ${path_gwas}pre_imputation/QC/netherlands_gsa/logs/stdout_plink_0_netherlands_cases_gsa \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/netherlands_cases/netherlands_cases_gsa \
--missing-genotype - \
--output-missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/netherlands_gsa/netherlands_cases_gsa_hg19"
# 700078 variants and 4155 people pass filters and QC.


# MERGE THESE TWO SETS:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/netherlands_controls_gsa_hg19 \
--bmerge /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/netherlands_cases_gsa_hg19.bed \
/path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/netherlands_cases_gsa_hg19.bim \
/path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/netherlands_cases_gsa_hg19.fam \
--merge-mode 3 \
--make-bed --out ${path_gwas}pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19
# 700078 variants and 4705 people pass filters and QC.

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19 \
--freq \
--out ${path_gwas}pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19


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

bim<-read.table(paste(path,"pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        0      A      C      D      G      I      T
# 0   1483   2963    305      0    785      8      6
# A      0      0  72118      0 324308      0   1948
# C      0  58229      0      0   2872      0      0
# D      0      0      0      0      0   2817      0
# G      0 226801   2837      0      0      0      0
# I      0      0      0   1123      0      0      0
# T      0   1475      0      0      0      0      0

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 9490      6

table(indels$V1)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 1751  533  846  572  278  443  315  379  222  299  247  489  280  372  164  298 
# 16   17   18   19   20   21   22   23   24   26 
# 258  467  113  193   97  103  132  496  129   14 

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
# 7387  533  846  572  278  443  315  379  222  299  247  489  280  372  164  298 
# 16   17   18   19   20   21   22   23   24   26 
# 258  467  113  193   97  103  132  496  129  140

dim(all_remove)
#[1] 15252      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/netherlands_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

##########################################################
# KEEP ONLY THE CASES AND CONTROLS AS DONE IN PREVIOUS QC

pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# netherlands_cases netherlands_controls 
# 4155                  550 


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#                           Crohn's Disease Indeterminate Ulcerative Colitis
# netherlands_cases       1            2440           106               1548
# netherlands_controls  550               0             0                  0
# 
#                      Unknown
# netherlands_cases         60
# netherlands_controls       0

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                      Affected Unaffected
# netherlands_cases        4155          0
# netherlands_controls        0        550

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#                         0    1
# netherlands_cases    4155    0
# netherlands_controls    0  550  

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                    Affected Unaffected
#                           1        550
# Crohn's Disease        2440          0
# Indeterminate           106          0
# Ulcerative Colitis     1548          0
# Unknown                  60          0

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#                      Female Male
# netherlands_cases      2441 1714
# netherlands_controls    317  233


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
# 1    2 
# 550 4155 

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# 1    2 
# 1947 2758 

fam_test<-read.table(paste(path,"pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 4705

write.table(fam,paste(path,"pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 0  2

write.table(samples_remove,paste(path,"pre_imputation/QC/netherlands_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/netherlands_gsa/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/netherlands_gsa/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19_noind
# 684826 variants and 4705 people pass filters and QC.
# Among remaining phenotypes, 4155 are cases and 550 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37


