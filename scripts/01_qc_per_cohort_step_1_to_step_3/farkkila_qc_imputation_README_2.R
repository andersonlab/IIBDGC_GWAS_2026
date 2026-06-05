# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
hiesp_broad_ibd_daly_farkkila_gsa.txt.gz

# new batch of text files downloaded to /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/

cd /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/
md5sum hiesp_broad_ibd_daly_farkkila_gsa.txt.gz
# da6944e10173634dd299f7f12bd8c982  hiesp_broad_ibd_daly_farkkila_gsa.txt.gz

less /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/file-manifest.json 
# "file_name": "hiesp_broad_ibd_daly_farkkila_gsa.txt.gz"  - OK


zcat /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/hiesp_broad_ibd_daly_farkkila_gsa.txt.gz \
| head -10

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/farkkila/
  
# create a tmp copy with final file ID, # to keep consistency in naming - only copy while files are still downloading (hpc-server client checks out 
# downloaded files to evaluate whether those have been downloaded), remove afterwards
cp /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/hiesp_broad_ibd_daly_farkkila_gsa.txt.gz \
/path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/farkkila_gsa.txt.gz

############################################################################## 
# make a temporary sorted fcr file:

# farm5:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh farkkila long gsa
# Job <798704> is submitted to queue <long>.

study=farkkila
less ${path_gwas}scripts/logs/fcrsort_stderr_${study} 

zcat /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/farkkila_gsa.txt.gz | wc -l
# 74995489
zcat /path/to/ibdgwas/IIBDGC/pre_imputation/raw/farkkila/farkkila_gsa_sorted.txt.gz | wc -l
# 74995489

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
bash ${path_gwas}scripts/create_map_from_fcr.sh farkkila long gsa
# Job <840474> is submitted to queue <long>.

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

# farm5:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_fam_from_fcr.sh farkkila normal gsa
#Job <840486> is submitted to queue <normal>.

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

# farm5:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_lgen_from_fcr.sh farkkila normal gsa
# Job <840493> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/


###########################

wc -l /path/to/ibdgwas/IIBDGC/pre_imputation/raw/farkkila/*.fam
# 109 /path/to/ibdgwas/IIBDGC/pre_imputation/raw/farkkila/farkkila_gsa.fam

wc -l /path/to/ibdgwas/IIBDGC/pre_imputation/raw/farkkila/*.map
# 688032 /path/to/ibdgwas/IIBDGC/pre_imputation/raw/farkkila/farkkila_gsa.map

zcat /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/farkkila_gsa.txt.gz \
| awk 'END { print NR }'

# 74995489/688032
# [1] 109 OK

# remove initial copy:
rm /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/farkkila_gsa.txt.gz

##########################################################################################################################################
##########################################################################################################################################


##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/farkkila_gsa/
mkdir ${path_gwas}pre_imputation/QC/farkkila_gsa/logs/
  
MEM=800

bsub -J"pl" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/farkkila_gsa/logs/stderr_plink_0_farkkila_gsa \
-o ${path_gwas}pre_imputation/QC/farkkila_gsa/logs/stdout_plink_0_farkkila_gsa \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/farkkila/farkkila_gsa \
--missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19"
# Job <699758> is submitted to queue <normal>.

tail -200 ${path_gwas}pre_imputation/QC/farkkila_gsa/logs/stdout_plink_0_farkkila_gsa
# 688032 variants and 109 people pass filters and QC.
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

bim<-read.table(paste(path,"pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 159675     6

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
#[1] 159721      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/farkkila_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

pheno<-fread("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes-4313226.csv.gz")
pheno<-as.data.frame(pheno)
fam<-read.table(paste(path,"pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# hiesp_broad_ibd_daly_farkkila_gsa 
# 109 

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
  #                                      Crohn's Disease Ulcerative Colitis
  # hiesp_broad_ibd_daly_farkkila_gsa 41              10                 56
  #                                  
  #                                   Unknown
  # hiesp_broad_ibd_daly_farkkila_gsa       2


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                                      Affected
# hiesp_broad_ibd_daly_farkkila_gsa 41       68

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
# 0
# hiesp_broad_ibd_daly_farkkila_gsa 104

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
  #                       Affected
  #                    41        0
  #  Crohn's Disease    0       10
  # Ulcerative Colitis  0       56
  # Unknown             0        2

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
# Female Male
# hiesp_broad_ibd_daly_farkkila_gsa  5     49   55


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
# 2 -9 
# 68 41

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# 1   2  -9 
# 55 49  5

fam_test<-read.table(paste(path,"pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 109

write.table(fam,paste(path,"pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 41 2

write.table(samples_remove,paste(path,"pre_imputation/QC/farkkila_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")


######

path_gwas=/path/to/ibdgwas/IIBDGC/
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/farkkila_gsa/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/farkkila_gsa/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19_noind
# 528311 variants and 68 people pass filters and QC.
# Among remaining phenotypes, 68 are cases and 0 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37

