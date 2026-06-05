# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# text files downloaded to /path/to/project

md5sum -c niddk_feinstein_gsa.txt.gz.md5
# niddk_feinstein_gsa.txt.gz: OK

zcat /path/to/project | head -10

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/niddk_feinstein/

  
############################################################################## 
# make a temporary sorted fcr file:

# farm3:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh niddk_feinstein basement
# Job <4402770> is submitted to queue <basement>.

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

# farm3, updated script to indicate which queue:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_map_from_fcr.sh niddk_feinstein long
# Job <5610101> is submitted to queue <long>.

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
bash ${path_gwas}scripts/create_fam_from_fcr.sh niddk_feinstein 
#Job <5440684> is submitted to default queue <normal>.

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
bash ${path_gwas}scripts/create_lgen_from_fcr.sh niddk_feinstein 
#Job <5440802> is submitted to default queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/


############################
# remove sorted file:
rm /path/to/ibdgwas/IIBDGC/pre_imputation/raw/niddk_feinstein/niddk_feinstein_gsa_sorted.txt.gz

##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/niddk_feinstein_gsa/
mkdir ${path_gwas}pre_imputation/QC/niddk_feinstein_gsa/logs/
MEM=25000

bsub -J"bcft_gwas3" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}pre_imputation/QC/niddk_feinstein_gsa/logs/stderr_lfile_to_bed \
-o ${path_gwas}pre_imputation/QC/niddk_feinstein_gsa/logs/stdout_lfile_to_bed \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/niddk_feinstein/niddk_feinstein_gsa \
--missing-genotype - \
--make-bed --out ${path_gwas}pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19"
# Job <706876> is submitted to default queue <normal>.

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

bim<-read.table(paste(path,"pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        -      A      C      D      G      I      T
# -   2322  19719   8373    906  27879   2289    757
# A      0      0  65941      0 303113      0   1457
# C      0  54516      0      0   2358      0      0
# D      0      0      0      0      0    976      0
# G      0 216116   2224      0      0      0      0
# I      0      0      0    379      0      0      0
# T      0   1143      0      0      0      0      0

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 63600      6

table(indels$V1)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 6188 5058 5150 3775 2481 2974 2848 3086 2002 2292 2240 3393 2702 1636 1715 2056 
# 16   17   18   19   20   21   22   23   24   25   26 
# 2006 2916  966 2280 1155  727  970 2473  470   32    9 

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
# 7387 5058 5150 3775 2481 2974 2848 3086 2002 2292 2240 3393 2702 1636 1715 2056 
# 16   17   18   19   20   21   22   23   24   25   26 
# 2006 2916  966 2280 1155  727  970 2473  470   32  140 

dim(all_remove)
#[1] 64930     6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/niddk_feinstein_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# niddk_feinstein 
# 8687 


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#                      Crohn's Disease Indeterminate Ulcerative Colitis Unknown
# niddk_feinstein 2824            4145            98               1611       9

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                      Affected Unaffected Unknown
# niddk_feinstein  378     5855       2446       8

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#                    0    1
# niddk_feinstein 6241 2446

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)])
      ,as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                         Affected Unaffected Unknown
#                     378        0       2446       0
# Crohn's Disease       0     4145          0       0
# Indeterminate         0       98          0       0
# Ulcerative Colitis    0     1611          0       0
# Unknown               0        1          0       8

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#                      Female Male Unknown
# niddk_feinstein   16   4300 4365       6


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
# 1    2   -9 
# 2446 5855  386 

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# 1    2   -9 
# 4365 4300   22 

fam_test<-read.table(paste(path,"pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 8687

write.table(fam,paste(path,"pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 386  2

write.table(samples_remove,paste(path,"pre_imputation/QC/niddk_feinstein_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")


######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/niddk_feinstein_gsa/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/niddk_feinstein_gsa/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19_noind
# 645538 variants and 8301 people pass filters and QC.
# Among remaining phenotypes, 5855 are cases and 2446 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37


