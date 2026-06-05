# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

#########################
# kiel, austria, sibdcs #
#########################

# text files downloaded to /path/to/project

# FIVE SEPARATED FILES
cd /path/to/project
  
md5sum -c kiel_foc_gsa.txt.gz.md5
# kiel_foc_gsa.txt.gz: OK
md5sum -c kiel_eze_gsa.txt.gz.md5
# kiel_eze_gsa.txt.gz: OK
md5sum -c kiel_bc_gsa.txt.gz.md5
# kiel_bc_gsa.txt.gz: OK
md5sum -c kiel_ibd_gsa.txt.gz.md5
# kiel_ibd_gsa.txt.gz: OK
md5sum -c kiel_hlitm_gsa.txt.gz.md5
# kiel_hlitm_gsa.txt.gz: OK
md5sum -c austria_gsa.txt.gz.md5
# austria_gsa.txt.gz: OK
md5sum -c sibdcs_gsa.txt.gz.md5
# sibdcs_gsa.txt.gz: OK

zcat /path/to/project | head -2
zcat /path/to/project | head -2
zcat /path/to/project | head -2
zcat /path/to/project | head -2
zcat /path/to/project | head -2
zcat /path/to/project | head -2
zcat /path/to/project | head -2

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/kiel_foc/
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/kiel_eze/
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/kiel_bc/
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/kiel_ibd/
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/kiel_hlitm/
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/austria/
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/sibdcs/

############################################################################## 
# make a temporary sorted fcr file:

# farm4:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh kiel_foc long gsa
# Job <1241894> is submitted to queue <long>.

bash ${path_gwas}scripts/sort_fcr.sh kiel_eze long gsa
# Job <1241895> is submitted to queue <long>.

bash ${path_gwas}scripts/sort_fcr.sh kiel_bc basement gsa
# Job <6015864> is submitted to queue <basement>.

bash ${path_gwas}scripts/sort_fcr.sh kiel_ibd basement gsa
# Job <5805371> is submitted to queue <basement>.

bash ${path_gwas}scripts/sort_fcr.sh kiel_hlitm basement gsa
# Job <5805373> is submitted to queue <basement>.

bash ${path_gwas}scripts/sort_fcr.sh austria basement gsa
# Job <8600195> is submitted to queue <basement>.

bash ${path_gwas}scripts/sort_fcr.sh sibdcs basement gsa
# Job <8639813> is submitted to queue <basement>.

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
bash ${path_gwas}scripts/create_map_from_fcr.sh kiel_foc normal gsa
# Job <1282259> is submitted to default queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_map_from_fcr.sh kiel_eze normal gsa
# Job <1282260> is submitted to default queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_map_from_fcr.sh kiel_bc normal gsa
# Job <6044622> is submitted to default queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_map_from_fcr.sh kiel_ibd normal gsa
# Job <6004936> is submitted to queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_map_from_fcr.sh kiel_hlitm normal gsa
# Job <6004937> is submitted to queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_map_from_fcr.sh austria long gsa
# Job <8625705> is submitted to queue <long>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_map_from_fcr.sh sibdcs normal gsa
# Job <8882207> is submitted to queue <long>.

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
bash ${path_gwas}scripts/create_fam_from_fcr.sh kiel_foc normal gsa 
# Job <1282261> is submitted to default queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_fam_from_fcr.sh kiel_eze normal gsa
# Job <1282262> is submitted to default queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_fam_from_fcr.sh kiel_bc normal gsa
#Job <6044624> is submitted to default queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_fam_from_fcr.sh kiel_ibd normal gsa
# Job <6004948> is submitted to queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_fam_from_fcr.sh kiel_hlitm normal gsa
# Job <6004949> is submitted to queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_fam_from_fcr.sh austria normal gsa
# Job <8625706> is submitted to queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_fam_from_fcr.sh sibdcs normal gsa
# Job <8657082> is submitted to queue <normal>.

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
bash ${path_gwas}scripts/create_lgen_from_fcr.sh kiel_foc normal gsa 
#Job <1282264> is submitted to default queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_lgen_from_fcr.sh kiel_eze normal gsa 
# Job <1282265> is submitted to default queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_lgen_from_fcr.sh kiel_bc normal gsa 
# Job <6044625> is submitted to default queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_lgen_from_fcr.sh kiel_ibd normal gsa
# Job <6004951> is submitted to queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_lgen_from_fcr.sh kiel_hlitm normal gsa
# Job <6004952> is submitted to queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_lgen_from_fcr.sh austria normal gsa
# Job <8625708> is submitted to queue <normal>.

path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_lgen_from_fcr.sh sibdcs normal gsa
# Job <8657083> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/

# (base) user-server:/path/to/ibdgwas/IIBDGC/pre_imputation/raw$ wc -l kiel_*/*.fam
# 1946 kiel_bc/kiel_bc_gsa.fam
# 358 kiel_eze/kiel_eze_gsa.fam
# 1098 kiel_foc/kiel_foc_gsa.fam
# 3582 kiel_hlitm/kiel_hlitm_gsa.fam
# 3523 kiel_ibd/kiel_ibd_gsa.fam
# 10507 total
# 
# (base) user-server:/path/to/ibdgwas/IIBDGC/pre_imputation/raw$ wc -l kiel_*/*.map
# 700078 kiel_bc/kiel_bc_gsa.map
# 700078 kiel_eze/kiel_eze_gsa.map
# 700078 kiel_foc/kiel_foc_gsa.map
# 700078 kiel_hlitm/kiel_hlitm_gsa.map
# 700078 kiel_ibd/kiel_ibd_gsa.map
# 3500390 total
# 
# 
# zcat /path/to/project | awk 'END { print NR }'
# 768685645
# 768685645/700078
# [1] 1098 # OK
# 
# zcat /path/to/project | awk 'END { print NR }'
# 250627925
# 250627925/700078
# [1] 358 # OK
# 
# zcat /path/to/project | awk 'END { print NR }'
# 1362351789
# 1362351789/700078
# [1] 1946 # OK
# 
# zcat /path/to/project | awk 'END { print NR }'
# 2466374795
# 2466374795/700078
# [1] 3523 # OK
# 
# zcat /path/to/project | awk 'END { print NR }'
# 2507679397
# 2507679397/700078
# [1] 3582 # OK
# 
# zcat /path/to/project | awk 'END { print NR }'



wc -l austria_gsa.fam
# 1832 austria_gsa.fam

wc -l austria_gsa.map
# 712189 austria_gsa.map

zcat /path/to/project | awk 'END { print NR }'
# 1304730249

1304730249/712189
# [1] 1832 OK



wc -l sibdcs_gsa.fam
# 2322

wc -l sibdcs_gsa.map
# 712189

zcat /path/to/project | awk 'END { print NR }'
# 1653702859

1653702859/2322
# [1] 712189 OK


##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/kiel_austria_sibdcs_gsa
mkdir ${path_gwas}pre_imputation/QC/kiel_austria_sibdcs_gsa/logs

/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/kiel_foc/kiel_foc_gsa \
--missing-genotype - \
--output-missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_foc_gsa_hg19
# 700078 variants and 1098 people pass filters and QC.

/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/kiel_eze/kiel_eze_gsa \
--missing-genotype - \
--output-missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_eze_gsa_hg19
# 700078 variants and 358 people pass filters and QC.

/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/kiel_bc/kiel_bc_gsa \
--missing-genotype - \
--output-missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_bc_gsa_hg19
# 700078 variants and 1946 people pass filters and QC.

/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/kiel_ibd/kiel_ibd_gsa \
--missing-genotype - \
--output-missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_ibd_gsa_hg19
# 700078 variants and 3523 people pass filters and QC.

/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/kiel_hlitm/kiel_hlitm_gsa \
--missing-genotype - \
--output-missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_hlitm_gsa_hg19
# 700078 variants and 3582 people pass filters and QC.

/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/austria/austria_gsa \
--missing-genotype - \
--output-missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/kiel_austria_sibdcs_gsa/austria_gsa_hg19
# 712189 variants and 1832 people pass filters and QC.

/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/sibdcs/sibdcs_gsa \
--missing-genotype - \
--output-missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/kiel_austria_sibdcs_gsa/sibdcs_gsa_hg19
# 712189 variants and  people pass filters and QC.


####################################################################################################################
# double check variants shared:

### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

a<-read.table(paste(path,"pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_hlitm_gsa_hg19.bim",sep=""),head=F)
b<-read.table(paste(path,"pre_imputation/QC/kiel_austria_sibdcs_gsa/sibdcs_gsa_hg19.bim",sep=""),head=F)

dim(a)
# [1] 700078      6
dim(b)
# [1] 712189      6

dim(a[which(a$V2 %in% b$V2),])
# [1] 662507      

dim(b[which(b$V2 %in% a$V2),])
# [1] 662507      6


####################################################################################################################

# MERGE THESE TWO SETS:

### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohort<-c("austria","kiel_foc","kiel_eze","kiel_bc","kiel_ibd","kiel_hlitm")

dat<-matrix(ncol=3,nrow=length(cohort))
dat<-as.data.frame(dat)

for (i in 1:length(cohort)){
  dat[i,1]<-paste(path,"pre_imputation/QC/kiel_austria_sibdcs_gsa/",cohort[i],"_gsa_hg19.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/kiel_austria_sibdcs_gsa/",cohort[i],"_gsa_hg19.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/kiel_austria_sibdcs_gsa/",cohort[i],"_gsa_hg19.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/kiel_austria_sibdcs_gsa/list_kiel_files_tomerge.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/kiel_austria_sibdcs_gsa/sibdcs_gsa_hg19 \
--merge-list ${path_gwas}pre_imputation/QC/kiel_austria_sibdcs_gsa/list_kiel_files_tomerge.txt \
--merge-mode 3 \
--make-bed --out ${path_gwas}pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19
# 749760 variants and 14661 people pass filters and QC.
# Note: No phenotypes present.


# Warning: Variants 'seq-rs387907306.1' and 'seq-rs387907306' have the same
# position.
# Warning: Variants 'seq-rs143229412.1' and 'seq-rs143229412' have the same
# position.
# Warning: Variants 'rs121434297' and '1:11855218' have the same position.
# 1600 more same-position warnings: see log file.
# Performing single-pass merge (10507 people, 700078 variants).


/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19 \
--freq \
--out ${path_gwas}pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19

#######################################################################################################################################

#####################
# 1.- REMOVE INDELS #
#####################

### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)

path<-"/path/to/ibdgwas/IIBDGC/"

bim<-read.table(paste(path,"pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        0      A      C      D      G      I      T
# 0  20716  25746   8053    773  23429   2211    831
# A      0      0  69218      0 315749      0   1682
# C      0  54716      0      0   2605      0      0
# D      0      0      0      0      0   1620      0
# G      0 218233   2218      0      0      0      0
# I      0      0      0    732      0      0      0
# T      0   1228      0      0      0      0      0

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)

table(indels$V1)


# X    X chromosome                    -> 23
# Y    Y chromosome                    -> 24
# XY   Pseudo-autosomal region of X    -> 25
# MT   Mitochondrial                   -> 26

mt<-bim[which(bim$V1 %in% c(0,26)),]
table(mt$V1)
# 0   26 
# 6954 1195 

chry<-bim[which(bim$V1 %in% c(24)),]
table(chry$V1)
# 24 
# 5262

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])
# 

dim(all_remove)
#[1] 86055     6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/kiel_austria_sibdcs_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

##########################################################
# KEEP ONLY THE CASES AND CONTROLS AS DONE IN PREVIOUS QC


pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# austria    kiel  sibdcs 
#    1832   10507    2322 


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#              Crohn's Disease Indeterminate Ulcerative Colitis Unknown
# austria    0            1126            35                375     296
# kiel    4680            4551           397                879       0
# sibdcs     1            1296            87                918      20

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#              Affected Unaffected
# austria    0     1832          0
# kiel       0     5827       4680
# sibdcs     1     2321          0

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#            0    1
# austria 1832    0
# kiel    5827 4680
# sibdcs  2321    0

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                         Affected Unaffected
#                       1        0       4680
# Crohn's Disease       0     6973          0
# Indeterminate         0      519          0
# Ulcerative Colitis    0     2172          0
# Unknown               0      316          0

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#              Female Male Unknown
# austria    0    943  888       1
# kiel       0   5342 5165       0
# sibdcs     1   1123 1197       1


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
#  -9    1    2 
# 1 4680 9980  

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# -9    1    2 
#  3 7250 7408

fam_test<-read.table(paste(path,"pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 14661

write.table(fam,paste(path,"pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 1  2

write.table(samples_remove,paste(path,"pre_imputation/QC/kiel_austria_sibdcs_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")


######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/kiel_austria_sibdcs_gsa/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/kiel_austria_sibdcs_gsa/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19_noind
# 663705 variants and 14660 people pass filters and QC.
# Among remaining phenotypes, 9980 are cases and 4680 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37


