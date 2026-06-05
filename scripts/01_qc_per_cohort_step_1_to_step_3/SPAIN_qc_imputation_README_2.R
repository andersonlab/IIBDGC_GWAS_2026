# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# text files downloaded to /path/to/project

cd /path/to/project
md5sum -c spain_gsa.txt.gz.md5
# spain_gsa.txt.gz: OK

zcat /path/to/project | head -10
# Sample ID	Sample Group	SNP Name	Chr	Position	Allele1 - Top	Allele2 - Top	X Raw	Y Raw	XGC Score	GT Score
# sample_id		200003	9	139026180	A	A	2724	165	0.9670094	0.04792601	0.890285	
# sample_id		200003	9	139026180	A	A	5991	284	0.8249298	0.03843912	0.890285

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/spain/
  
############################################################################## 
# make a temporary sorted fcr file, updated file to indicate queue:

# farm3:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh spain basement
# Job <4669565> is submitted to queue <basement>.

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
  
bash ${path_gwas}scripts/create_map_from_fcr.sh spain normal
# Job <5643505> is submitted to queue <normal>.

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
bash ${path_gwas}scripts/create_fam_from_fcr.sh spain 
#Job <5643509> is submitted to default queue <normal>.

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
bash ${path_gwas}scripts/create_lgen_from_fcr.sh spain 
#Job <5643515> is submitted to default queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/


############################
# remove sorted file:

rm /path/to/ibdgwas/IIBDGC/pre_imputation/raw/spain/spain_gsa_sorted.txt.gz


##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/spain_gsa/
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/spain/spain_gsa \
--missing-genotype - \
--make-bed --out ${path_gwas}pre_imputation/QC/spain_gsa/spain_gsa_hg18
# 620901 variants and 3443 people pass filters and QC.
# Note: No phenotypes present.

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

bim<-read.table(paste(path,"pre_imputation/QC/spain_gsa/spain_gsa_hg18.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#        -      A      C      G      T
# -  28824   2405    868   2706    155
# A      0      0  56935 253319    638
# C      0  51990      0    981      0
# G      0 220577    922      0      0
# T      0    581      0      0      0

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 34958     6

table(indels$V1)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 2489 2191 1476 1557 1732 2666 2255 1300 1144 1642 1338 1105  805  872 1278 1434 
# 17   18   19   20   21   22   23   24   25   26 
# 1250  577  910  574  398  535 2974 2067  378   11


# X    X chromosome                    -> 23
# Y    Y chromosome                    -> 24
# XY   Pseudo-autosomal region of X    -> 25
# MT   Mitochondrial                   -> 26

mt<-bim[which(bim$V1 %in% c(0,26)),]
table(mt$V1)
# 26 
# 138


all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 2489 2191 1476 1557 1732 2666 2255 1300 1144 1642 1338 1105  805  872 1278 1434 
# 17   18   19   20   21   22   23   24   25   26 
# 1250  577  910  574  398  535 2974 2067  378  138 

dim(all_remove)
#[1] 35085      6
write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/spain_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

##########################################################
# KEEP ONLY THE CASES AND CONTROLS AS DONE IN PREVIOUS QC


pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
fam<-read.table(paste(path,"pre_imputation/QC/spain_gsa/spain_gsa_hg18.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# spain 
# 3443 


table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
#            Crohn's Disease Ulcerative Colitis
# spain 1482            1164                797

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#       Affected Unaffected
# spain     1961       1482

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#          0    1
# spain 1961 1482

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                    Affected Unaffected
#                           0       1482
# Crohn's Disease        1164          0
# Ulcerative Colitis      797          0

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#       Female Male
# spain   1531 1912


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
# 1482 1961 

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# 1    2 
# 1912 1531 

fam_test<-read.table(paste(path,"pre_imputation/QC/spain_gsa/spain_gsa_hg18.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 3443

write.table(fam,paste(path,"pre_imputation/QC/spain_gsa/spain_gsa_hg18_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

######

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 0  2

write.table(samples_remove,paste(path,"pre_imputation/QC/spain_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")


################################################

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/spain_gsa/spain_gsa_hg18.bed \
--bim ${path_gwas}/pre_imputation/QC/spain_gsa/spain_gsa_hg18.bim \
--fam ${path_gwas}/pre_imputation/QC/spain_gsa/spain_gsa_hg18_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/spain_gsa/list_indel_var_exclude \
--make-bed --out ${path_gwas}pre_imputation/QC/spain_gsa/spain_gsa_hg18_noind
# 585816 variants and 3443 people pass filters and QC.
# Among remaining phenotypes, 1961 are cases and 1482 are controls.


#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# NOTE: data in build 36, needs to be liftover first:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/spain_gsa/spain_gsa_hg18_noind \
--recode tab --out ${path_gwas}pre_imputation/QC/spain_gsa/spain_gsa_hg18_noind

# ${path_gwas}previous_qced_b38/liftover/liftOverx
# usage:
# liftOver oldFile map.chain newFile unMapped

#### R convert map into bed:

map<-read.table(paste(path,"pre_imputation/QC/spain_gsa/spain_gsa_hg18_noind.map",sep=""),head=F)

map$chr<-paste("chr",map$V1,sep="")
map$chr[which(map$V1=="23")]<-"chrX"
map$chr[which(map$V1=="25")]<-"chrX"
map$chr[which(map$V1=="24")]<-"chrY"
map$chromStart<-map$V4
map$chromEnd<-map$V4+1

write.table(map[,c(5:7,2)],paste(path,"pre_imputation/QC/spain_gsa/spain_gsa_hg18_noind.bed",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#### lift positions

${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}pre_imputation/QC/spain_gsa/spain_gsa_hg18_noind.bed \
${path_gwas}previous_qced_b38/liftover/hg18ToHg19.over.chain.gz \
${path_gwas}pre_imputation/QC/spain_gsa/spain_gsa_hg18_noind_lifted_hg19 \
${path_gwas}pre_imputation/QC/spain_gsa/spain_gsa_hg18_noind_no_lifted_hg19

cut -f 4 ${path_gwas}pre_imputation/QC/spain_gsa/spain_gsa_hg18_noind_no_lifted_hg19 | sed "/^#/d" \
> ${path_gwas}pre_imputation/QC/spain_gsa/nonlifted_variants_to_exclude.dat

wc -l ${path_gwas}pre_imputation/QC/spain_gsa/nonlifted_variants_to_exclude.dat
# 82

### exclude non lifted:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--file ${path_gwas}pre_imputation/QC/spain_gsa/spain_gsa_hg18_noind \
--exclude ${path_gwas}pre_imputation/QC/spain_gsa/nonlifted_variants_to_exclude.dat \
--recode tab --out ${path_gwas}pre_imputation/QC/spain_gsa/spain_gsa_hg18_noind_liftedvariants
# 585734 variants and 3443 people pass filters and QC.

####

bed_lifted<-read.table(paste(path,"pre_imputation/QC/spain_gsa/spain_gsa_hg18_noind_lifted_hg19",sep=""),head=F)
map<-read.table(paste(path,"pre_imputation/QC/spain_gsa/spain_gsa_hg18_noind_liftedvariants.map",sep=""),head=F)
                
table(map$V2==bed_lifted$V4)
# TRUE 
# 585734 

map$pos<-bed_lifted$V2
table(map$V2==bed_lifted$V4)
table(map$pos==bed_lifted$V2)
# TRUE 
# 585734 

write.table(map[,c(1:3,5)],paste(path,"pre_imputation/QC/spain_gsa/spain_gsa_hg18_noind_liftedvariants_edited.map",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


########

# put together updated .ped plus lifted.map
/path/to/software/plink_linux_x86_64_20181202/./plink \
--ped ${path_gwas}pre_imputation/QC/spain_gsa/spain_gsa_hg18_noind_liftedvariants.ped \
--map ${path_gwas}pre_imputation/QC/spain_gsa/spain_gsa_hg18_noind_liftedvariants_edited.map \
--make-bed --out ${path_gwas}pre_imputation/QC/spain_gsa/spain_gsa_hg19_noind
# 585734 variants and 3443 people pass filters and QC.

# remove intermediate files:

rm ${path_gwas}pre_imputation/QC/spain_gsa/spain_gsa_hg18_noind_*
  

    