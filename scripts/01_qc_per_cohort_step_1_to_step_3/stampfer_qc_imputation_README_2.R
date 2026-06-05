# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
{
  "object_id": "dg.EA80/5f2c43cb-a63f-4c51-aa2c-8912d0c64871",
  "md5sum": "1dc85fcf60d547fd92c2baaf5703200e",
  "file_name": "ccdg_ibd_daly_stampfer_ibd_gsa-array_v1_4.txt.gz",
  "file_size": 18240303117
},

# new batch of text files downloaded to /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/

cd /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/
md5sum ccdg_ibd_daly_stampfer_ibd_gsa-array_v1_4.txt.gz
#  1dc85fcf60d547fd92c2baaf5703200e  ccdg_ibd_daly_stampfer_ibd_gsa-array_v1_4.txt.gz

zcat /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/ccdg_ibd_daly_stampfer_ibd_gsa-array_v1_4.txt.gz \
| head -10

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/stampfer/
  
# create a tmp copy with final file ID, # to keep consistency in naming - only copy while files are still downloading (hpc-server client checks out 
# downloaded files to evaluate whether those have been downloaded), remove afterwards
cp /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/ccdg_ibd_daly_stampfer_ibd_gsa-array_v1_4.txt.gz \
/path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/stampfer_gsa.txt.gz

############################################################################## 
# make a temporary sorted fcr file:

# farm5:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh stampfer long gsa
# Job <877983> is submitted to queue <long>.

study=stampfer
less ${path_gwas}scripts/logs/fcrsort_stderr_${study} 

zcat /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/stampfer_gsa.txt.gz | wc -l
# 1321709473
zcat /path/to/ibdgwas/IIBDGC/pre_imputation/raw/stampfer/stampfer_gsa_sorted.txt.gz | wc -l
# 1321709473

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
  bash ${path_gwas}scripts/create_map_from_fcr.sh stampfer long gsa
# Job <924310> is submitted to queue <long>.

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
  bash ${path_gwas}scripts/create_fam_from_fcr.sh stampfer normal gsa
#Job <924311> is submitted to queue <normal>.

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
  bash ${path_gwas}scripts/create_lgen_from_fcr.sh stampfer normal gsa
# Job <924312> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/


###########################

wc -l /path/to/ibdgwas/IIBDGC/pre_imputation/raw/stampfer/*.fam
# 1921 /path/to/ibdgwas/IIBDGC/pre_imputation/raw/stampfer/stampfer_gsa.fam

wc -l /path/to/ibdgwas/IIBDGC/pre_imputation/raw/stampfer/*.map
# 688032 /path/to/ibdgwas/IIBDGC/pre_imputation/raw/stampfer/stampfer_gsa.map

zcat /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/stampfer_gsa.txt.gz \
| awk 'END { print NR }'


# 1321709473/688032
# [1] 1921 OK

# remove initial copy:
rm /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/stampfer_gsa.txt.gz

##########################################################################################################################################
##########################################################################################################################################



##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
  mkdir ${path_gwas}pre_imputation/QC/stampfer_gsa/
  mkdir ${path_gwas}pre_imputation/QC/stampfer_gsa/logs/
  
  MEM=1000

bsub -J"pl" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/stampfer_gsa/logs/stderr_plink_0_stampfer_gsa \
-o ${path_gwas}pre_imputation/QC/stampfer_gsa/logs/stdout_plink_0_stampfer_gsa \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/stampfer/stampfer_gsa \
--missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19"
# Job <701014> is submitted to queue <normal>.

tail -200 ${path_gwas}pre_imputation/QC/stampfer_gsa/logs/stdout_plink_0_stampfer_gsa
# 688032 variants and 1921 people pass filters and QC.
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

bim<-read.table(paste(path,"pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 15267     6

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
#[1] 95324      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/stampfer_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

pheno<-fread("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes-4313226.csv.gz")
pheno<-as.data.frame(pheno)
fam<-read.table(paste(path,"pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# ccdg_ibd_daly_stampfer_ibd_gsa-array_v1_4 
# 1921

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#                                                Crohn's Disease
# ccdg_ibd_daly_stampfer_ibd_gsa-array_v1_4 1483             197
#                                            
#                                             Ulcerative Colitis
# ccdg_ibd_daly_stampfer_ibd_gsa-array_v1_4                241


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                                                Affected Unaffected
# ccdg_ibd_daly_stampfer_ibd_gsa-array_v1_4  444      438       1039

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#                                              0    1
# ccdg_ibd_daly_stampfer_ibd_gsa-array_v1_4  438 1039

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                         Affected Unaffected
#                     444        0       1039
# Crohn's Disease       0      197          0
# Ulcerative Colitis    0      241          0

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#                                                Female Male
# ccdg_ibd_daly_stampfer_ibd_gsa-array_v1_4  444   1238  239 


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
# 1039  438  444 

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# 1    2   -9 
# 239 1238  444 

fam_test<-read.table(paste(path,"pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 1921

write.table(fam,paste(path,"pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 444 2

write.table(samples_remove,paste(path,"pre_imputation/QC/stampfer_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")


######

path_gwas=/path/to/ibdgwas/IIBDGC/
  
  /path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/stampfer_gsa/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/stampfer_gsa/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19_noind
# 672634 variants and 1477 people pass filters and QC.
# Among remaining phenotypes, 438 are cases and 1039 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37


