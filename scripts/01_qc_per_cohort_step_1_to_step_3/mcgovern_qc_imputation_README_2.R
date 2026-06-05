# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# {
#   "object_id": "dg.EA80/4f6b3478-22ad-4467-87d9-0e9da2a8cd3f",
#   "md5sum": "1d25d5f964c6c51fa4bb773ac100e910",
#   "file_name": "anvil_ccdg_broad_ai_ibd_daly_mcgovern_gsa.txt.gz",
#   "file_size": 56552775388
# },


# new batch of text files downloaded to /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/

cd /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/
md5sum anvil_ccdg_broad_ai_ibd_daly_mcgovern_gsa.txt.gz
# 1d25d5f964c6c51fa4bb773ac100e910  anvil_ccdg_broad_ai_ibd_daly_mcgovern_gsa.txt.gz

zcat /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/anvil_ccdg_broad_ai_ibd_daly_mcgovern_gsa.txt.gz \
| head -10

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/mcgovern/
  
# create a tmp copy with final file ID, # to keep consistency in naming - only copy while files are still downloading (hpc-server client checks out 
# downloaded files to evaluate whether those have been downloaded), remove afterwards
cp /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/anvil_ccdg_broad_ai_ibd_daly_mcgovern_gsa.txt.gz \
/path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/mcgovern_gsa.txt.gz

############################################################################## 
# make a temporary sorted fcr file:

# farm5 - mem=2500 
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh mcgovern basement gsa
# Job <877743> is submitted to queue <long>.

study=mcgovern
less ${path_gwas}scripts/logs/fcrsort_stderr_${study} 

zcat /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/mcgovern_gsa.txt.gz | wc -l
# 4188738817
zcat /path/to/ibdgwas/IIBDGC/pre_imputation/raw/mcgovern/mcgovern_gsa_sorted.txt.gz | wc -l
# 4188738817

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
  bash ${path_gwas}scripts/create_map_from_fcr.sh mcgovern long gsa
# Job <541391> is submitted to queue <long>.

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
  bash ${path_gwas}scripts/create_fam_from_fcr.sh mcgovern normal gsa
#Job <541397> is submitted to queue <normal>.

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
  bash ${path_gwas}scripts/create_lgen_from_fcr.sh mcgovern normal gsa
# Job <541409> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/


###########################

wc -l /path/to/ibdgwas/IIBDGC/pre_imputation/raw/mcgovern/*.fam
# 6088 /path/to/ibdgwas/IIBDGC/pre_imputation/raw/mcgovern/mcgovern_gsa.fam

wc -l /path/to/ibdgwas/IIBDGC/pre_imputation/raw/mcgovern/*.map
# 688032 /path/to/ibdgwas/IIBDGC/pre_imputation/raw/mcgovern/mcgovern_gsa.map

zcat /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/mcgovern_gsa.txt.gz \
| awk 'END { print NR }'
# 4188738817

# 4188738817/688032
# [1] 6088 OK

# remove initial copy:
rm /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/mcgovern_gsa.txt.gz

##########################################################################################################################################
##########################################################################################################################################


##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/mcgovern_gsa/
mkdir ${path_gwas}pre_imputation/QC/mcgovern_gsa/logs/
  
MEM=4000

bsub -J"pl" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/mcgovern_gsa/logs/stderr_plink_0_mcgovern_gsa \
-o ${path_gwas}pre_imputation/QC/mcgovern_gsa/logs/stdout_plink_0_mcgovern_gsa \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/mcgovern/mcgovern_gsa \
--missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19"
# Job <706215> is submitted to queue <normal>.

tail -200 ${path_gwas}pre_imputation/QC/mcgovern_gsa/logs/stdout_plink_0_mcgovern_gsa
# 688032 variants and 6088 people pass filters and QC.
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

bim<-read.table(paste(path,"pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 35857     6

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
#[1] 35986      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/mcgovern_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

pheno<-fread("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes-4313226.csv.gz")
pheno<-as.data.frame(pheno)
fam<-read.table(paste(path,"pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# anvil_ccdg_broad_ai_ibd_daly_mcgovern_gsa 
# 6088 

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#                                                Crohn's Disease Indeterminate
# anvil_ccdg_broad_ai_ibd_daly_mcgovern_gsa 1245            3005            13
#                                            
#                                           Ulcerative Colitis Unknown
# anvil_ccdg_broad_ai_ibd_daly_mcgovern_gsa               1562     263

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                                                Affected Unaffected Unknown
# anvil_ccdg_broad_ai_ibd_daly_mcgovern_gsa   35     4840       1210       3

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#                                              0    1
# anvil_ccdg_broad_ai_ibd_daly_mcgovern_gsa 4840 1210

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
  #                         Affected Unaffected Unknown
  #                      35        0       1210       0
  # Crohn's Disease       0     3005          0       0
  # Indeterminate         0       13          0       0
  # Ulcerative Colitis    0     1562          0       0
  # Unknown               0      260          0       3

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#                                                Female Male Unknown
# anvil_ccdg_broad_ai_ibd_daly_mcgovern_gsa   35   3068 2967      18


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
#   1    2   -9 
# 1210 4840   38

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
#    1    2   -9 
# 2967 3068   53

fam_test<-read.table(paste(path,"pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 6088


write.table(fam,paste(path,"pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 38 2

write.table(samples_remove,paste(path,"pre_imputation/QC/mcgovern_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")


######

path_gwas=/path/to/ibdgwas/IIBDGC/
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/mcgovern_gsa/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/mcgovern_gsa/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19_noind
# 652046 variants and 6050 people pass filters and QC.
# Among remaining phenotypes, 4840 are cases and 1210 are controls. 



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37


