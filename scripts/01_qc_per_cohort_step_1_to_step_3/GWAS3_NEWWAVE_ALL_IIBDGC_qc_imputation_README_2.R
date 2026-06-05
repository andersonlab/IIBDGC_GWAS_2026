# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
###############################
# GWAS3, NewWave all together #
###############################

path_gwas=/path/to/ibdgwas/IIBDGC/

mkdir ${path_gwas}pre_imputation/QC/all_hce/

#############################################################################################################
# /path/to/project less analysis_steps.txt 

# Create a merged new_wave call set that combines:
#       29 GWAS3 samples released after March 2015 (Batch 5, August 2015, not included in initial analysis)
#       274 GWAS3 samples removed from initial analysis due to overlap with sequencing data
#       1516 "New wave" IBD GWAS samples, released June 2016
#       4201 Dupytren's disease samples, to act as controls


#############################################################################################################
# /path/to/project less analysis_steps.txt 
  
# Create merged "new wave" dataset including Katie's previously processed "new wave" samples plus extra 554 (final "new wave" samples).
# previously processed "new wave" samples: /path/to/project
# 1818 are cases and 4201 are controls
# latest "new wave" batch: /path/to/project
# 554 are cases

#1) Merge both datasets
# bsub -J "merge" -R"select[mem>12000] rusage[mem=12000]" -M12000 -G crohns -o /path/to/project 
# "/path/to/software/plink2 --noweb --allow-no-sex 
# --bfile /path/to/project 
# --bmerge /path/to/project
# /path/to/project 
# /path/to/project --make-bed --out merged_callset_2.gencall"

# exclude any variant NewWave array-specific variant that is not included in GWAS3

##############################################
# /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(dplyr)

path<-"/path/to/ibdgwas/IIBDGC/"

fam<-read.table(paste(path,"pre_imputation/QC/new_wave/new_wave_hg19.fam",sep=""),head=F)
table(fam$V6)
# 1    2 
# 4201 2372

table(fam$V6,fam$V5)
#      0    1    2
# 1    0 3312  889
# 2  811  766  795

fam_remove<-fam[which(fam$V6==1),]
fam_remove<-fam_remove[,1:2]
colnames(fam_remove)<-c("FID","IID")
dim(fam_remove)
# [1] 4201    2

fam_remove2<-read.table(paste(path,"pre_imputation/QC/gwas3/list_gwas3_ids_duplicated_in_newwave",sep=""),head=T)
dim(fam_remove2)
# [1] 273   2

fam_remove2$FID<-paste("urn:wtsi:",fam_remove2$FID,sep="")
fam_remove2$IID<-paste("urn:wtsi:",fam_remove2$IID,sep="")

fam_remove<-rbind(fam_remove,fam_remove2)
dim(fam_remove)
# [1] 4474    2

fam_remove<-fam_remove[!duplicated(fam_remove$FID),]
dim(fam_remove)
# [1] 4474    2

dim(fam[which(fam$V1 %in% fam_remove$FID),])
# [1] 4474    6 # OK

write.table(fam_remove,paste(path,"pre_imputation/QC/all_hce/list_new_wave_controls_gwas3_samples_to_exclude",sep=""),col.names=T,row.names=F,quote=F)

################

g3<-read.table(paste(path,"pre_imputation/QC/gwas3/gwas3_hg19.bim",sep=""),head=F)
nw<-read.table(paste(path,"pre_imputation/QC/new_wave/new_wave_hg19.bim",sep=""),head=F)

dim(g3)
# [1] 535434      6
dim(nw)
# [1] 557662      6

dim(nw[which(nw$V2 %in% g3$V2),])
# [1] 535434      6

var<-nw[which(!nw$V2 %in% g3$V2),]
dim(var)
# [1] 22228     6
write.table(var,paste(path,"pre_imputation/QC/all_hce/list_new_wave_variants_to_exclude_before_merging_with_gwas3",sep=""),col.names=F,row.names=F,quote=F)

##########################

# exclude and remove
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}pre_imputation/QC/new_wave/new_wave_hg19.bed \
--bim ${path_gwas}pre_imputation/QC/new_wave/new_wave_hg19.bim \
--fam ${path_gwas}pre_imputation/QC/new_wave/new_wave_hg19_tmp.fam \
--exclude ${path_gwas}/pre_imputation/QC/all_hce/list_new_wave_variants_to_exclude_before_merging_with_gwas3 \
--remove ${path_gwas}/pre_imputation/QC/all_hce/list_new_wave_controls_gwas3_samples_to_exclude \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_tmp
# --exclude: 535434 variants remaining.
# --remove: 2099 people remaining.
# 535434 variants and 2372 people pass filters and QC.
# Among remaining phenotypes, 2372 are cases and 0 are controls.



# # merge all files
# /path/to/software/plink_linux_x86_64_20181202/./plink \
# --bed ${path_gwas}pre_imputation/QC/gwas3/gwas3_hg19_tmp.bed \
# --bim ${path_gwas}pre_imputation/QC/gwas3/gwas3_hg19_tmp.bim \
# --fam ${path_gwas}pre_imputation/QC/gwas3/gwas3_hg19_tmp_edited.fam \
# --bmerge ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_tmp.bed ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_tmp.bim \
# ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_tmp.fam \
# --allow-no-sex \
# --make-bed --out ${path_gwas}pre_imputation/QC/all_hce/all_hce_hg19

# Error: 246077 variants with 3+ alleles present.
# * If you believe this is due to strand inconsistency, try --flip with
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/all_hce/all_hce_hg19-merge.missnp.
# (Warning: if this seems to work, strand errors involving SNPs with A/T or C/G
#   alleles probably remain in your data.  If LD between nearby SNPs is high,
#   --flip-scan should detect them.)
# * If you are dealing with genuine multiallelic variants, we recommend exporting
# that subset of the data to VCF (via e.g. '--recode vcf'), merging with
# another tool/script, and then importing the result; PLINK is not yet suited
# to handling them.

# issue when attempted, forgot both were not + strand

#### run intermediate steps to (1) align to plus strand, and (2) use 1000GP to update A/T G/C




#######################################################################################################################################

#####################
# 1.- REMOVE INDELS #
#####################

### /software/R-4.3.1/bin/R

g3<-read.table(paste(path,"pre_imputation/QC/gwas3/gwas3_hg19.bim",sep=""),head=F)
nw<-read.table(paste(path,"pre_imputation/QC/all_hce/new_wave_hg19_tmp.bim",sep=""),head=F)

table(g3$V5,g3$V6)
#        A      C      D      G      I      T
# 0      0      3      0      0      0      0
# A      0  25295      0 127574      0   4271
# C  19105      0      0   9618      0  76397
# D      0      0      0      0   8482      0
# G  76608  10046      0      0      0  18936
# I      0      0   4150      0      0      0
# T   4207 126113      0  24629      0      0

table(nw$V5,nw$V6)
#        A      C      D      G      I      T
# A      0  49904      0 253808      0   4931
# C  38061      0      0  10388      0      0
# D      0      0      0      0   8483      0
# G 152886   9277      0      0      0      0
# I      0      0   4149      0      0      0
# T   3547      0      0      0      0      0

table(g3$V1)
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 48011 41021 34517 28492 28740 33486 26908 23772 22563 23957 29741 26472 14919 
# 14    15    16    17    18    19    20    21    22    23    24    26 
# 17037 16925 19268 21070 11954 21696 13290  6873  9116 13113  2098   395 

table(nw$V1)
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 47998 41021 34518 28493 28738 33487 26910 23772 22561 23955 29740 26471 14917 
# 14    15    16    17    18    19    20    21    22    23    24    25    26 
# 17037 16926 19269 21069 11953 21696 13288  6873  9116 13095  2115     4   412 

g3i<-g3[which(!g3$V5 %in% c("A","G","C","T") | !g3$V6 %in% c("A","G","C","T")),]
dim(g3i)
# [1] 12635     6

nwi<-nw[which(!nw$V5 %in% c("A","G","C","T") | !nw$V6 %in% c("A","G","C","T")),]
dim(nwi)
# [1] 12632     6

dim(nwi[which(nwi$V2 %in% g3i$V2),])
# [1] 12632     6

g3mt<-g3[which(g3$V1==26),]
dim(g3mt)
# [1] 395   6

nwmt<-nw[which(nw$V1==26),]
dim(nwmt)
# [1] 412   6


dim(nwmt[which(nwmt$V2 %in% g3mt$V2),])
# [1] 395   6

dim(nwmt[which(nwmt$V2 %in% g3$V2),])
# [1] 412   6

table(g3[which(g3$V2 %in% nwmt$V2),"V1"])
# 1  11  13  26 
# 14   1   2 395 


g3[which(g3$V2 %in% nwmt$V2),]
# V1             V2 V3        V4 V5 V6
# 1       1     exm2216283  0    564766  C  T
# 2       1     200610-117  0    564768  C  T
# 3       1     exm2216284  0    564862  T  C
# 4       1     exm2216292  0    565111  C  T
# 5       1     MitoG4821A  0    565370  A  G
# 6       1     MitoA4825G  0    565374  G  A
# 7       1     MitoA5657G  0    566206  G  A
# 8       1     MitoA5952G  0    566501  G  A
# 9       1      200610-82  0    567284  A  G
# 10      1      200610-10  0    567302  G  A
# 11      1     MitoT6777C  0    567327  C  T
# 12      1     MitoA7056G  0    567606  G  A
# 13      1     MitoG7522A  0    568072  A  G
# 14      1     MitoA8870G  0    569418  G  A
# 316109 11     MitoT1191C  0  10531217  G  A
# 381422 13     200610-109  0 110076590  G  A
# 381423 13     200610-108  0 110076602  G  A

# annotated to wrong chromosome!!!! exclude all


# keep common lists:
indels<-g3[which(!g3$V5 %in% c("A","G","C","T") | !g3$V6 %in% c("A","G","C","T")),]
mt<-nw[which(nw$V1==26),]

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

dim(all_remove)
#[1] 13047      6

### some variants are in different chr and position, remove as well:
g3<-g3[match(nw$V2,g3$V2),]
table(g3$V2==nw$V2)
# TRUE 
# 535434


colnames(g3)[c(1,3:6)]<-paste(colnames(g3)[c(1,3:6)],"_g3",sep="")
colnames(nw)[c(1,3:6)]<-paste(colnames(nw)[c(1,3:6)],"_nw",sep="")

all<-merge(g3,nw,by="V2")

dim(all[which(all$V1_g3!=all$V1_nw),])
# [1] 62 11
head(all[which(all$V1_g3!=all$V1_nw),])
#               V2 V1_g3 V3_g3     V4_g3 V5_g3 V6_g3 V1_nw V3_nw    V4_nw V5_nw
# 1062   200610-10     1     0    567302     G     A    26     0     6753     G
# 1067  200610-108    13     0 110076602     G     A    26     0     1109     G
# 1068  200610-109    13     0 110076590     G     A    26     0     1121     G
# 1075  200610-117     1     0    564768     C     T    26     0     4219     G
# 1243   200610-82     1     0    567284     A     G    26     0     6735     A
# 14983 exm1097142     9     0 108285772     A     G    14     0 36143784     A
#       V6_nw
# 1062      A
# 1067      A
# 1068      A
# 1075      A
# 1243      G
# 14983     G

dim(all[which(all$V4_g3!=all$V4_nw),])
# [1] 2571   11


dim(all[which((all$V4_g3!=all$V4_nw) | (all$V1_g3!=all$V1_nw)),])
# [1] 2575   11

all_remove3<-all[which((all$V4_g3!=all$V4_nw) | (all$V1_g3!=all$V1_nw)),c("V1_g3","V2","V3_g3","V4_g3","V5_g3","V6_g3")]
colnames(all_remove3)<-colnames(all_remove)

dim(all_remove3[which(all_remove3$V2 %in% all_remove$V2),])
# [1] 1277    6

all_remove<-rbind(all_remove,all_remove3)
all_remove<-all_remove[!duplicated(all_remove$V2),]

dim(all_remove)
# [1] 14345     6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/all_hce/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

#########
# GWAS3 #
#########

##########################################################
# KEEP ONLY THE CASES AND CONTROLS AS DONE IN PREVIOUS QC

## R

fam<-read.table(paste(path,"pre_imputation/QC/gwas3/keep/gwas3_hg19.fam",sep=""),head=F)

# SEX:
table(fam$V5,useNA="ifany")
# 0     1     2 
# 2883  8917 10452

# pheno:
table(fam$V6,useNA="ifany")
# 1     2 
# 10484 11768 

############################################
# double check assigned pheno from phenoDB, plus sex. There is an issue reported by Loukas see "/path/to/ibdgwas/pre_imputation/qc/gwas3/qc_steps_gwas3.txt":

# see script "manifest_approach.R" to see the matches

# fam_file_ids<-read.table(paste(path,"/manifests/match_sanger_ids_supplier_ids_fam_manifest.tsv",sep=""),sep="\t",head=T)
# for (i in 1:nrow(fam_file_ids)){
#   fam$id[which(fam$id %in% fam_file_ids$sanger_sample_name_fam_file[i])]<-as.character(fam_file_ids$sanger_sample_name_manifest[i])
# }


# AFTER MEETING WITH CARL AND ALEX, KEEP PHENO IN FAM, THEY RELY MORE IN KATIE'S FILES THAN IN PHENODB

sample<-read.table("/path/to/ibdgwas/post_imputation/GWAS3/GWAS3.sample",head=T)
sample<-sample[-1,] # remove subheader
dim(sample)
#[1] 18355   15
dim(fam)
#[1] 21979    6

dim(fam[which(fam$V1 %in% sample$ID_1),])
#[1] 18355    6

fam_subset<-fam[which(fam$V1 %in% sample$ID_1),]
fam_subset<-fam_subset[match(sample$ID_1,fam_subset$V1),]
table(as.character(fam_subset$V1)==as.character(sample$ID_1))
#TRUE 
#18355 

table(as.character(fam_subset$V6),as.character(sample$bin1))
#      0    1
# 1 9495    0
# 2    0 8860


table(as.character(fam_subset$V5),as.character(sample$sex))
#      0    1    2
# 0 1757    0    0
# 1    0 7550    0
# 2    0    0 9048


# same phenotypes and sex as in postQC sample file

################################################
# Affection status, by default, should be coded:
# -9 missing 
# 0 missing
# 1 unaffected
# 2 affected

table(fam$V6)
# 1     2 
# 10484 11768 

################################################
# Sex (1=male; 2=female; other=unknown)

table(fam$V5)
# 0     1     2 
# 2833  8818 10328 

# summary of phenotypes:
pheg3<-fread("/path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/UK_IIBDGC/coreex_gaibdc_phenotypes.tsv",head=T)

dim(pheg3[which(pheg3$orig_sample_id %in% fam$V1),])

dim(pheg3[which(pheg3$orig_sample_id %in% fam$V1),])
# [1] 11768    45
dim(pheg3)
# [1] 11930    45

table(fam$V6)
# 1     2 
# 10484 11768 

pheg3<-pheg3[which(pheg3$orig_sample_id %in% fam$V1),]
table(pheg3$diag)
# Crohn's Disease      Indeterminate Ulcerative Colitis            Unknown 
#            5923                 19               5452                374 

table(pheg3$sex)
# Female    Male Unknown 
#  2942    2561    6265

table(fam$V5)
# 0     1     2 
# 2833  8818 10328 

# worse thatn in fam file

fam_test<-read.table(paste(path,"/pre_imputation/QC/gwas3/gwas3_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 22252

write.table(fam,paste(path,"pre_imputation/QC/gwas3/gwas3_hg19_edited.fam",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



### NEW WAVE:
fam<-read.table(paste(path,"pre_imputation/QC/new_wave/new_wave_hg19.fam",sep=""),head=F)

# SEX:
table(fam$V5,useNA="ifany")
#   0    1    2 
# 811 4078 1684 

# pheno:
table(fam$V6,useNA="ifany")
#   1    2 
#4201 2372 

fam$id<-gsub(".*_","",fam$V1)

# AFTER MEETING WITH CARL AND ALEX, KEEP PHENO IN FAM, THEY RELY MORE IN KATIE'S FILES THAN IN PHENODB

sample<-read.table(paste(path,"pre_imputation/raw/new_wave/new_wave.ibd.sample",sep=""),head=T)
sample<-sample[-1,] # remove subheader
dim(sample)
#[1] 5384   15
dim(fam)
#[1] 6573    6

dim(fam[which(fam$V1 %in% sample$ID_1),])
#[1] 5384    6

fam_subset<-fam[which(fam$V1 %in% sample$ID_1),]
fam_subset<-fam_subset[match(sample$ID_1,fam_subset$V1),]
table(as.character(fam_subset$V1)==as.character(sample$ID_1))
#TRUE 
#5384 

table(as.character(fam_subset$V6),as.character(sample$bin1),useNA="ifany")
#      0    1 <NA>
# 1 3905    0    0
# 2    0 1477    2

table(as.character(fam_subset$V5)==as.character(sample$sex))
#TRUE 
#5384

# same phenotypes and sex as in postQC sample file

# edit famility IDs to 0:
table(fam$V3)
#-9 
#6573 

table(fam$V4)
#-9 
#6573 

# force the samples to be founders by recoding the father and mother IDs
fam$V3<-0
fam$V4<-0

fam_test<-read.table(paste(path,"pre_imputation/raw/new_wave/new_wave_hg19.fam",sep=""),head=F)
table(fam_test$V1==fam$V1)
#TRUE 
#6573

write.table(fam,paste(path,"pre_imputation/QC/new_wave/new_wave_hg19_tmp.fam",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


#########################################################################
# Identify new_wave samples, and final pheno of remaining samples:

fam<-read.table(paste(path,"pre_imputation/QC/all_hce/new_wave_hg19_tmp.fam",sep=""))

phen_new_wave<-fread("/path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/UK_IIBDGC/coreex_gaibdc_phenotypes.tsv",head=T,sep="\t")
phen_nw<-fread("/path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/UK_IIBDGC/coreex_idbgwas_Final_phenotypes.tsv",head=T,sep="\t")

phen<-rbind(phen_new_wave,phen_nw)

table(fam$V6)
# 2 
# 2099

fam$id<-gsub("urn:wtsi:","",fam$V1)

dim(phen_nw[which(phen_nw$orig_sample_id %in% fam$id),])
# [1] 2070   45

dim(phen[which(phen$orig_sample_id %in% fam$id),])
# [1] 2099   45

phen<-phen[which(phen$orig_sample_id %in% fam$id),]
table(phen$diag,useNA="ifany")
# Crohn's Disease      Indeterminate Ulcerative Colitis            Unknown 
#             281                 93                303               1422  



##############


/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}pre_imputation/QC/gwas3/gwas3_hg19.bed \
--bim ${path_gwas}pre_imputation/QC/gwas3/gwas3_hg19.bim \
--fam ${path_gwas}pre_imputation/QC/gwas3/gwas3_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/all_hce/list_indel_var_exclude \
--make-bed --out ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_noind
# 521089 variants and 22252 people pass filters and QC.
# Among remaining phenotypes, 11768 are cases and 10484 are controls.

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_tmp \
--exclude ${path_gwas}/pre_imputation/QC/all_hce/list_indel_var_exclude \
--make-bed --out ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind
# 521089 variants and 2099 people pass filters and QC.
# Among remaining phenotypes, 2099 are cases and 0 are controls.


#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37


#######################################################################################################################################

#########################
# 4.- ALIGN TO + STRAND #
#########################


/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_noind \
--allow-no-sex \
--output-chr M --recode vcf-iid --out ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19
# 521089 variants and 22252 people pass filters and QC.
# Among remaining phenotypes, 11768 are cases and 10484 are controls.


/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind \
--allow-no-sex \
--output-chr M --recode vcf-iid --out ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19
# 521089 variants and 2099 people pass filters and QC.
# Among remaining phenotypes, 2099 are cases and 0 are controls.

###########################################################
# 4.2 FLIP ALLELES TO + STRAND:

export BCFTOOLS_PLUGINS=/path/to/software/bcftools-1.9/plugins

/path/to/software/bcftools-1.9/./bcftools +fixref ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19.vcf \
-Oz -o ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_posstrandaligned.vcf.gz \
-- -f /path/to/project -m top
# # SC, guessed strand convention
# SC	TOP-compatible	0
# SC	BOT-compatible	0
# # ST, substitution types
# ST	A>C	19049	3.7%
# ST	A>G	76269	14.6%
# ST	A>T	4177	0.8%
# ST	C>A	25216	4.8%
# ST	C>G	10024	1.9%
# ST	C>T	125807	24.1%
# ST	G>A	127239	24.4%
# ST	G>C	9592	1.8%
# ST	G>T	24554	4.7%
# ST	T>A	4243	0.8%
# ST	T>C	76055	14.6%
# ST	T>G	18864	3.6%
# # NS, Number of sites:
# NS	total        	521089
# NS	ref match    	431279	82.8%
# NS	ref mismatch 	89810	17.2%
# NS	flipped      	509	0.1%
# NS	swapped      	89300	17.1%
# NS	flip+swap    	13477	2.6%
# NS	unresolved   	1	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0


/path/to/software/bcftools-1.9/./bcftools +fixref ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19.vcf \
-Oz -o ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_posstrandaligned.vcf.gz \
-- -f /path/to/project -m top
# # SC, guessed strand convention
# SC	TOP-compatible	0
# SC	BOT-compatible	0
# # ST, substitution types
# ST	A>C	37930	7.3%
# ST	A>G	152203	29.2%
# ST	A>T	3524	0.7%
# ST	C>A	49753	9.5%
# ST	C>G	9250	1.8%
# ST	C>T	0	0.0%
# ST	G>A	253167	48.6%
# ST	G>C	10366	2.0%
# ST	G>T	0	0.0%
# ST	T>A	4896	0.9%
# ST	T>C	0	0.0%
# ST	T>G	0	0.0%
# # NS, Number of sites:
# NS	total        	521089
# NS	ref match    	217439	41.7%
# NS	ref mismatch 	303650	58.3%
# NS	flipped      	214408	41.1%
# NS	swapped      	44894	8.6%
# NS	flip+swap    	44826	8.6%
# NS	unresolved   	1	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0


# REMOVE MISSMATCH ALLELES, THEY  IMPUTATION SERVER TO CRASH:
/path/to/software/bcftools-1.9/./bcftools norm \
--check-ref x \
-f /path/to/project \
${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_posstrandaligned.vcf.gz \
-Oz -o ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_posstrandaligned_2.vcf.gz
# Lines   total/split/realigned/skipped:	521089/0/0/1

/path/to/software/bcftools-1.9/./bcftools norm \
--check-ref x \
-f /path/to/project \
${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_posstrandaligned.vcf.gz \
-Oz -o ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_posstrandaligned_2.vcf.gz
# Lines   total/split/realigned/skipped:	521089/0/0/1

###########################################################
# 4.3 DOUBLE-CHECK THE CONVERSION

/path/to/software/bcftools-1.9/./bcftools +fixref ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_posstrandaligned_2.vcf.gz \
-- -f /path/to/project


/path/to/software/bcftools-1.9/./bcftools +fixref ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_posstrandaligned_2.vcf.gz \
-- -f /path/to/project


###########################################################
# 4.4 VCF to BED

/path/to/software/plink_linux_x86_64_20181202/./plink \
--vcf ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_posstrandaligned_2.vcf.gz \
--double-id \
--make-bed --out ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr \
--missing \
--out ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr



/path/to/software/plink_linux_x86_64_20181202/./plink \
--vcf ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_posstrandaligned_2.vcf.gz \
--double-id \
--make-bed --out ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr \
--missing \
--out ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr


########################################################################################################################################

################################################################################
# 5.- UPDATE NAME VARIANTS TO CHR:POSITION_REF_ALT; REMOVE DUPLICATED VARIANTS #
################################################################################

zcat ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_posstrandaligned_2.vcf.gz | cut -f '1-5' | awk '{print $3,$1":"$2"_"$4"_"$5}' \
> ${path_gwas}pre_imputation/QC/all_hce/list_variants_gwas3_hg19_posstrandaligned


###### R

## GWAS3

bim<-read.table(paste(path,"pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr.bim",sep=""),sep="\t",head=F)
ids<-read.table(paste(path,"pre_imputation/QC/all_hce/list_variants_gwas3_hg19_posstrandaligned",sep=""),sep=" ",head=F,skip=35)

varmiss<-read.table(paste(path,"pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr.lmiss",sep=""),head=T)

colnames(ids)[2]<-"ids"

table(bim$V2==ids$V1)
# TRUE 
# 521088 


#the major allele is set to A2 by default by Plink, keep ids with real ref/alt as in vcf using the ids file
bim.1<-cbind(bim,ids[,"ids",drop=F])

bim.1<-merge(bim.1,varmiss[,c("SNP","F_MISS")],by.x="V2",by.y="SNP",sort=F)
table(bim.1$V2==bim$V2)
# TRUE 
# 521088
bim.1$ids<-as.character(bim.1$ids)

# identify duplicated variants (same chr position ref and alt)
dups<-bim.1[which(duplicated(bim.1$ids)),"ids"]

length(dups)
#[1] 7013

bim.1[which(duplicated(bim.1$V2)),]
#<0 rows> (or 0-length row.names) OK

for (i in 1:length(dups)){
  
  tmp<-bim.1[which(bim.1$ids %in% dups[i]),]
  
  keep<-tmp[which(tmp$F_MISS==min(tmp$F_MISS)),]
  
  if (nrow(keep)>1){
    keep<-keep[1,]
  }
  
  exclude<-tmp[which(!tmp$V2 %in% keep$V2),]
  
  bim.1$ids[which(bim.1$V2 %in% exclude$V2)]<-paste(bim.1$ids[which(bim.1$V2 %in% exclude$V2)],"_rm",sep="")
  
}

dim(bim.1[which(duplicated(bim.1$ids)),])
#[1] 16  8
 

bim.1[which(duplicated(bim.1$ids)),]
# 

duplicated_variants<-bim.1[grep("_rm",bim.1$ids),"ids"]
length(duplicated_variants)
# [1] 7014


write.table(duplicated_variants,paste(path,"pre_imputation/QC/all_hce/list_duplicated_gwas3_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

# to exclude for NW
duplicated_variants_nw<-bim.1[grep("_rm",bim.1$ids),]

table(bim.1$V2==bim$V2)
# TRUE 
# 521088 

table(bim.1$V1)
# 


# NOTE chrx pseudoautosomal now in chr23, recode this

bim_old<-read.table(paste(path,"pre_imputation/QC/all_hce/gwas3_hg19_noind.bim",sep=""),head=F)
table(bim_old$V1)
#    1     2     3     4     5     6     7     8     9    10    11    12    13 
# 46630 39935 33635 27792 27973 32665 26252 23310 22089 23438 28892 25737 14640 
# 14    15    16    17    18    19    20    21    22    23    24 
# 16537 16485 18894 20458 11710 20851 13052  6737  8917 12899  1561 

table(bim.1$V1[which(bim.1$V2 %in% bim_old$V2[which(bim_old$V1==25)])])
# 0


# https://www.ncbi.nlm.nih.gov/grc/human
# Name Chr Start Stop
# PAR#1	X	60,001	2,699,520
# PAR#2	X	154,931,044	155,260,560
# PAR#1	Y	10,001	2,649,520
# PAR#2	Y	59,034,050	59,363,566

dim(bim.1[which( (bim.1$V1==23) & ((bim.1$V4>=60001 & bim$V4<=2699520) | (bim.1$V4>=154931044 & bim$V4<=155260560)) ),])
#[1] 149   8

dim(bim.1[which( (bim.1$V1==23) & ((bim.1$V4>=60001 & bim$V4<=2699520) | (bim.1$V4>=154931044 & bim$V4<=155260560)) & (!bim.1$ids %in% duplicated_variants) ),])
#[1] 149   8

dim(bim.1[which( (bim.1$V1==24) & ((bim.1$V4>=10001 & bim$V4<=2649520) | (bim.1$V4>=59034050 & bim$V4<=59363566)) ),])
# 0

write.table(bim.1[,c(2,7,3:6)],paste(path,"pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr_edited.bim",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



## NEW WAVE

bim<-read.table(paste(path,"pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr.bim",sep=""),sep="\t",head=F)
ids<-read.table(paste(path,"pre_imputation/QC/all_hce/list_variants_new_wave_hg19_posstrandaligned",sep=""),sep=" ",head=F,skip=35)

varmiss<-read.table(paste(path,"pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr.lmiss",sep=""),head=T)

colnames(ids)[2]<-"ids"

table(bim$V2==ids$V1)
# TRUE 
# 521088 


#the major allele is set to A2 by default by Plink, keep ids with real ref/alt as in vcf using the ids file
bim.1<-cbind(bim,ids[,"ids",drop=F])

bim.1<-merge(bim.1,varmiss[,c("SNP","F_MISS")],by.x="V2",by.y="SNP",sort=F)
table(bim.1$V2==bim$V2)
# TRUE 
# 521088
bim.1$ids<-as.character(bim.1$ids)

# identify duplicated variants (same chr position ref and alt)
dups<-bim.1[which(duplicated(bim.1$ids)),"ids"]

length(dups)
#[1] 7013


dim(bim.1[which(bim.1$V2 %in% duplicated_variants_nw$V2),])

bim.1$ids[which(bim.1$V2 %in% duplicated_variants_nw$V2)]<-paste(bim.1$ids[which(bim.1$V2 %in% duplicated_variants_nw$V2)],"_rm",sep="")
duplicated_variants<-bim.1[grep("_rm",bim.1$ids),"ids"]
length(duplicated_variants)
# [1] 7014

write.table(duplicated_variants,paste(path,"pre_imputation/QC/all_hce/list_duplicated_new_wave_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

table(bim.1$V2==bim$V2)
# TRUE 
# 521088 

table(bim.1$V1)
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 46630 39935 33635 27792 27973 32665 26252 23310 22089 23438 28892 25737 14640 
# 14    15    16    17    18    19    20    21    22    23    24 
# 16537 16485 18894 20458 11710 20851 13052  6736  8917 12899  1561 


# NOTE chrx pseudoautosomal now in chr23, recode this

bim_old<-read.table(paste(path,"pre_imputation/QC/all_hce/new_wave_hg19_noind.bim",sep=""),head=F)
table(bim_old$V1)
#     1     2     3     4     5     6     7     8     9    10    11    12    13 
# 46630 39935 33635 27792 27973 32665 26252 23310 22089 23438 28892 25737 14640 
# 14    15    16    17    18    19    20    21    22    23    24 
# 16537 16485 18894 20458 11710 20851 13052  6737  8917 12899  1561 

table(bim.1$V1[which(bim.1$V2 %in% bim_old$V2[which(bim_old$V1==25)])])
# 0


# https://www.ncbi.nlm.nih.gov/grc/human
# Name Chr Start Stop
# PAR#1	X	60,001	2,699,520
# PAR#2	X	154,931,044	155,260,560
# PAR#1	Y	10,001	2,649,520
# PAR#2	Y	59,034,050	59,363,566

dim(bim.1[which( (bim.1$V1==23) & ((bim.1$V4>=60001 & bim$V4<=2699520) | (bim.1$V4>=154931044 & bim$V4<=155260560)) ),])
# [1] 149   8

dim(bim.1[which( (bim.1$V1==23) & ((bim.1$V4>=60001 & bim$V4<=2699520) | (bim.1$V4>=154931044 & bim$V4<=155260560)) & (!bim.1$ids %in% duplicated_variants) ),])
#[1] 149   8

dim(bim.1[which( (bim.1$V1==24) & ((bim.1$V4>=10001 & bim$V4<=2649520) | (bim.1$V4>=59034050 & bim$V4<=59363566)) ),])
# 0

write.table(bim.1[,c(2,7,3:6)],paste(path,"pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr_edited.bim",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


##########

# NOTE: use pre-vcf fam file which keeps the ca/ctr info:


/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr.bed \
--bim ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr_edited.bim \
--fam ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_noind.fam \
--exclude ${path_gwas}pre_imputation/QC/all_hce/list_duplicated_gwas3_var_exclude \
--split-x 'b37' \
--make-bed --out ${path_gwas}/pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr_nodup
# 514074 variants and 22252 people pass filters and QC.
# Among remaining phenotypes, 11768 are cases and 10484 are controls.
# --split-x: 149 chromosome codes changed.


/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr.bed \
--bim ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr_edited.bim \
--fam ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind.fam \
--exclude ${path_gwas}pre_imputation/QC/all_hce/list_duplicated_new_wave_var_exclude \
--split-x 'b37' \
--make-bed --out ${path_gwas}/pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr_nodup
# 514074 variants and 2099 people pass filters and QC.
# Among remaining phenotypes, 2099 are cases and 0 are controls.
# --split-x: 149 chromosome codes changed.

##############################################################################################################################################

##############################################
# 6.- COMPARE ALLELE FREQUENCIES WITH 1000GP #
##############################################

#########################################################################
# 6.1 COMPARE ALLELE FREQUENCIES WITH 1000GP, GENERATE LIST OF VARIANTS:

cat ${path_gwas}/pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr_nodup.bim | cut -f 2 > ${path_gwas}/pre_imputation/QC/all_hce/list_variants_all_hce_hg19_noind_posstr_nodup


##############################################
# 6.2 EXTRACT VARIANTS FROM 1000GP

# see script 1000gp_data_to_compare.R about how 1000GP data was generated

for chr in {1..24}; do /path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/1000gp/1000GP_EUR_chr${chr}_b37 \
--extract ${path_gwas}/pre_imputation/QC/all_hce/list_variants_all_hce_hg19_noind_posstr_nodup \
--make-bed --out ${path_gwas}pre_imputation/QC/all_hce/1000GP_EUR_chr${chr}_b37_all_hce_variants;done


##############################################
# 6.3 COMBINE FILES

#### R

dat<-matrix(ncol=3,nrow=23)
dat<-as.data.frame(dat)
for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/all_hce/1000GP_EUR_chr",i+1,"_b37_all_hce_variants.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/all_hce/1000GP_EUR_chr",i+1,"_b37_all_hce_variants.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/all_hce/1000GP_EUR_chr",i+1,"_b37_all_hce_variants.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/all_hce/list_1000GP_files_all_hce_variants_tomerge.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/all_hce/1000GP_EUR_chr1_b37_all_hce_variants \
--merge-list ${path_gwas}pre_imputation/QC/all_hce/list_1000GP_files_all_hce_variants_tomerge.txt \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/all_hce/1000GP_EUR_b37_all_hce_variants
#422805 variants and 503 people pass filters and QC.

# remove intermediate files:
rm ${path_gwas}pre_imputation/QC/all_hce/1000GP_EUR_chr*_b37_all_hce_variants*
  
  
##############################################
# 6.4 ESTIMATE FREQUENCIES:

# 1000GP

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/all_hce/1000GP_EUR_b37_all_hce_variants \
--freq \
--out ${path_gwas}pre_imputation/QC/all_hce/1000GP_EUR_b37_all_hce_variants

# gwas3 

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr_nodup \
--freq \
--out ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr_nodup

# new_wave 

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr_nodup \
--allow-no-sex \
--freq \
--out ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr_nodup

#### R

gp<-read.table(paste(path,"pre_imputation/QC/all_hce/1000GP_EUR_b37_all_hce_variants.frq",sep=""),head=T)

# gwas3

g1<-read.table(paste(path,"pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr_nodup.frq",sep=""),head=T)

dim(g1)
#[1] 514074      6
dim(gp)
#[1] 422805      6

dim(g1[which(!g1$SNP %in% gp$SNP),])
#[1] 91269     6 # LARGE NUMBER OF VARIANTS

colnames(gp)[3:6]<-paste(colnames(gp)[3:6],"_gp",sep="")
colnames(g1)[3:6]<-paste(colnames(g1)[3:6],"_g1",sep="")

all<-merge(g1,gp[,2:6],by="SNP",all=T)

check<-all[which(all$A1_g1!=all$A1_gp),]
dim(check)
#[1] 15207   10

# keep only A/T C/G
check<-check[which( (check$A1_g1=="G" & check$A2_g1=="C") | (check$A1_g1=="C" & check$A2_g1=="G") | (check$A1_g1=="A" & check$A2_g1=="T") | (check$A1_g1=="T" & check$A2_g1=="A")),]
dim(check)
#[1] 9691    10

# LIST OF VARIANTS TO REMOVE, WE CANNOT REALLY BE SURE WHETHER THERE IS STRAND ISSUE OR NOT
remove<-check[which(check$MAF_g1>=0.45),]
dim(remove)
#[1] 108 10 - remove


# LIST OF VARIANTS TO FLIP:
flip<-check[which(check$MAF_g1<0.45),]
dim(flip)
#[1] 9473 10

flip<-flip[order(flip$MAF_g1,decreasing=T),]

pdf(paste(path,"pre_imputation/QC/all_hce/plot_maf_AT_CG_all_hce_1000gp.pdf",sep=""))
plot(flip$MAF_g1,flip$MAF_gp)
dev.off()

remove_2<-flip[which(flip$MAF_g1>0.2 & flip$MAF_gp<0.1),]
remove_3<-flip[which(flip$MAF_g1<0.1 & flip$MAF_gp>0.2),]

remove<-rbind(remove,remove_2,remove_3)
dim(remove)
#[1] 138 10

write.table(remove[,"SNP"],paste(path,"pre_imputation/QC/all_hce/list_variants_to_remove_AT_CG_gwas3",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

flip<-flip[which(!flip$SNP %in% remove$SNP),]
dim(flip)
#[1] 9553    10

write.table(flip[,"SNP"],paste(path,"pre_imputation/QC/all_hce/list_variants_to_flip_AT_CG_gwas3",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

# new_wave
g1<-read.table(paste(path,"pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr_nodup.frq",sep=""),head=T)

dim(g1)
#[1] 514074      6
dim(gp)
#[1] 422805      6

dim(g1[which(!g1$SNP %in% gp$SNP),])
#[1] 91269     6 # LARGE NUMBER OF VARIANTS

colnames(gp)[3:6]<-paste(colnames(gp)[3:6],"_gp",sep="")
colnames(g1)[3:6]<-paste(colnames(g1)[3:6],"_g1",sep="")

all<-merge(g1,gp[,2:6],by="SNP",all=T)

check<-all[which(all$A1_g1!=all$A1_gp),]
dim(check)
#[1] 6088   10

# keep only A/T C/G
check<-check[which( (check$A1_g1=="G" & check$A2_g1=="C") | (check$A1_g1=="C" & check$A2_g1=="G") | (check$A1_g1=="A" & check$A2_g1=="T") | (check$A1_g1=="T" & check$A2_g1=="A")),]
dim(check)
#[1] 56    10

# LIST OF VARIANTS TO REMOVE, WE CANNOT REALLY BE SURE WHETHER THERE IS STRAND ISSUE OR NOT
remove<-check[which(check$MAF_g1>=0.45),]
dim(remove)
#[1] 43 10 - remove


# LIST OF VARIANTS TO FLIP:
flip<-check[which(check$MAF_g1<0.45),]
dim(flip)
#[1] 13 10

flip<-flip[order(flip$MAF_g1,decreasing=T),]

pdf(paste(path,"pre_imputation/QC/all_hce/plot_maf_AT_CG_all_hce_1000gp.pdf",sep=""))
plot(flip$MAF_g1,flip$MAF_gp)
dev.off()

remove_2<-flip[which(flip$MAF_g1>0.2 & flip$MAF_gp<0.1),]
remove_3<-flip[which(flip$MAF_g1<0.1 & flip$MAF_gp>0.2),]

remove<-rbind(remove,remove_2,remove_3)
dim(remove)
#[1] 44 10

write.table(remove[,"SNP"],paste(path,"pre_imputation/QC/all_hce/list_variants_to_remove_AT_CG_new_wave",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

flip<-flip[which(!flip$SNP %in% remove$SNP),]
dim(flip)
#[1] 12 10

write.table(flip[,"SNP"],paste(path,"pre_imputation/QC/all_hce/list_variants_to_flip_AT_CG_new_wave",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

# create an unique to_remove file:

nw<-read.table(paste(path,"pre_imputation/QC/all_hce/list_variants_to_remove_AT_CG_new_wave",sep=""),head=F)
g3<-read.table(paste(path,"pre_imputation/QC/all_hce/list_variants_to_remove_AT_CG_gwas3",sep=""),head=F)

dim(nw)
# [1] 44  1
dim(g3)
# [1] 138   1

remove<-rbind(nw,g3)
dim(remove)
# [1] 182   1


remove<-remove[!duplicated(remove$V1),,drop=F]
dim(remove)
# [1] 160   1

write.table(remove[,"V1"],paste(path,"pre_imputation/QC/all_hce/list_variants_to_remove_AT_CG",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


##########################################################################################
# 6.5 REMOVE VARIANTS WE CANNOT BE SURE ARE IN THE RIGHT STRAND, AND FLIP THE OTHERS

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr_nodup \
--exclude ${path_gwas}pre_imputation/QC/all_hce/list_variants_to_remove_AT_CG \
--flip ${path_gwas}pre_imputation/QC/all_hce/list_variants_to_flip_AT_CG_gwas3 \
--make-bed --out ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr_nodup_flip
# --flip: 9553 SNPs flipped.
# 513914 variants and 22252 people pass filters and QC.
# Among remaining phenotypes, 11768 are cases and 10484 are controls.


/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr_nodup \
--exclude ${path_gwas}pre_imputation/QC/all_hce/list_variants_to_remove_AT_CG \
--flip ${path_gwas}pre_imputation/QC/all_hce/list_variants_to_flip_AT_CG_new_wave \
--make-bed --out ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr_nodup_flip
# --flip: 12 SNPs flipped.
# 513914 variants and 2099 people pass filters and QC.
# Among remaining phenotypes, 2099 are cases and 0 are controls.



# DOUBLE CHECK MAF IS SIMILAR:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr_nodup_flip \
--freq \
--out ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr_nodup_flip

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr_nodup_flip \
--freq \
--out ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr_nodup_flip

###### R

nw<-read.table(paste(path,"pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr_nodup_flip.frq",sep=""),head=T)
g3<-read.table(paste(path,"pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr_nodup_flip.frq",sep=""),head=T)
                              
table(nw$SNP==g3$SNP)
# TRUE 
# 513914            
table(nw$CHR==g3$CHR)
# TRUE 
# 513914 

colnames(nw)[c(3:6)]<-paste(colnames(nw)[c(3:6)],"_nw",sep="")
colnames(g3)[c(3:6)]<-paste(colnames(g3)[c(3:6)],"_g3",sep="")

all<-merge(g3[,2:5],nw[,2:5],by="SNP")

check<-(all[which(all$A1_g3!=all$A1_nw),])
dim(check)
# [1] 7104    7


table(check$A1_g3,check$A2_g3)
#      A    C    G    T
# A    0  145  520  755
# C  144    0 1614  665
# G  654 1287    0  167
# T  540  486  125    0

table(check$A1_nw,check$A2_nw)
#      A    C    G    T
# A    0  144  654  540
# C  145    0 1287  486
# G  520 1614    0  125
# T  755  665  167    0


table(check$A1_g3,check$A2_nw)
#      A    C    G    T
# A 1420    0    0    0
# C    0 2423    0    0
# G    0    0 2108    0
# T    0    0    0 1151

table(check$A2_g3,check$A1_nw)
#      A    C    G    T
# A 1338    0    0    0
# C    0 1918    0    0
# G    0    0 2259    0
# T    0    0    0 1587


# variants in 100GP:
check<-check[which(!check$SNP %in% gp$SNP),]
# [1] 2816    7

check<-check[which( (check$A1_g3=="G" & check$A2_g3=="C") | (check$A1_g3=="C" & check$A2_g3=="G") | (check$A1_g3=="A" & check$A2_g3=="T") | (check$A1_g3=="T" & check$A2_g3=="A")),]
dim(check)
# [1] 4196    7

write.table(check[,"SNP",drop=F],paste(path,"pre_imputation/QC/all_hce/list_variants_to_remove_AT_CG_diff_maf_nw_gws3",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

####################

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr_nodup_flip \
--exclude ${path_gwas}pre_imputation/QC/all_hce/list_variants_to_remove_AT_CG_diff_maf_nw_gws3 \
--make-bed --out ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr_nodup_flip_2
# 509718 variants and 22252 people pass filters and QC.
# Among remaining phenotypes, 11768 are cases and 10484 are controls.


/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr_nodup_flip \
--exclude ${path_gwas}pre_imputation/QC/all_hce/list_variants_to_remove_AT_CG_diff_maf_nw_gws3 \
--make-bed --out ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr_nodup_flip_2
# 509718 variants and 2099 people pass filters and QC.
# Among remaining phenotypes, 2099 are cases and 0 are controls.


### merge:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/all_hce/gwas3_hg19_noind_posstr_nodup_flip_2 \
--bmerge ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr_nodup_flip_2.bed ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr_nodup_flip_2.bim ${path_gwas}pre_imputation/QC/all_hce/new_wave_hg19_noind_posstr_nodup_flip_2.fam \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/all_hce/all_hce_hg19_noind_posstr_nodup_flip
# 509718 variants and 24351 people pass filters and QC.
# Among remaining phenotypes, 13867 are cases and 10484 are controls.

########################################################################################################################################
