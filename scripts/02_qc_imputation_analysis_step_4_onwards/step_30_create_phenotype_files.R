# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# #####################################################
# #  2.- CREATE PHENOTYPE/COVARIATE FILES FOR REGENIE #
# #####################################################
# 

# Phenotype file format
# Line 1 : Header with FID, IID and P phenotypes names.
#
# Followed by lines of P + 2 values. Space/tab separated. Each line contains individual FID and IID followed by P phenotype values
# (for binary traits, must be coded as 0=control, 1=case, NA=missing unless using --1).
#
# Samples listed in this file that are not in bgen/bed/pgen file are ignored. Genotyped samples that are not in this file are removed from the analysis.
#
# Missing values must be coded as NA.
#
# With QTs, missing values are mean-imputed in Step 1 and they are dropped when testing each phenotype in Step 2 (unless using --force-impute).
#
# With BTs, missing values are mean-imputed in Step 1 when fitting the level 0 linear ridge regression and they are dropped when fitting the
# level 1 logistic ridge regression for each trait . In Step 2, missing values are dropped when testing each trait.
#
# To remove all samples that have missing values at any of the P phenotypes, use option --strict in Step 1 and 2. This is also useful
# when analyzing a single trait to avoid making a new bed/pgen/bgen file just for the complete data set of individuals (so setting the
# phenotype values of individuals to remove to NA), although --remove can also be used in that situation.


#  1.4.1 CREATE HARMONISED PHENOTYPE FILES FOR REGENIE

#### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

array<-c("illumina370","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

illumina370<-c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","swedish_uc_old_gwas","niddk_old_gwas")
affymetrix6<-c("german_affy6_old_gwas","norway_affy6_old_gwas","gwas2")
humanomniexpress<-c("australia_omniexome")
affymetrix500<-c("gwas1")
humancoreexome<-c("all_hce")
humanomni1<-c("pittsburgh_gsa")
quad610<-c("spain_gsa")
gsa<-c("italy_gsa","kiel_austria_sibdcs_gsa"
       ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa",
       "lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa",
       "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
       "franke_gsa","hyams_protect_gsa","lewis_sparc_gsa",
       "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
       "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
       "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
       "xavier_share_gsa")
illuminaexome<-c("prism_nfe_gwas","helmsley_prism_gsa","helmsley_xavier_prism_gsa","finland_illugwas")


studies<-c(illumina370,affymetrix6,humanomniexpress,affymetrix500,humancoreexome,humanomni1,quad610,gsa,illuminaexome)

ancestry<-c("eur","noneur")


############################################################################################################################

# illumina370<-c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","swedish_uc_old_gwas","niddk_old_gwas")

# belgium_inf1_old_gwas

ph_tmp<-read.csv("/path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/UK_IIBDGC/ibdgc_other_cohorts/gwas/belgian/InfosAllIndividualsInfinium1.csv",head=T)
ph_tmp<-ph_tmp[,c("Inf1_IndividualID","Inf1_Phenotype")]
ph_tmp<-ph_tmp[,c(1,1,2)]
table(ph_tmp$Inf1_Phenotype)
# CD healthy
# 517     900

ph_tmp$ibd<-NA
ph_tmp$ibd[which(ph_tmp$Inf1_Phenotype=="healthy")]<-"0"
ph_tmp$ibd[which(ph_tmp$Inf1_Phenotype=="CD")]<-"1"

ph_tmp$cd<-NA
ph_tmp$cd[which(ph_tmp$Inf1_Phenotype=="healthy")]<-"0"
ph_tmp$cd[which(ph_tmp$Inf1_Phenotype=="CD")]<-"1"

ph_tmp$uc<-NA

ph_tmp<-ph_tmp[,c(1:2,4:6)]
colnames(ph_tmp)<-c("FID","IID","ibd","cd","uc")

ph_tmp$study<-"belgium_inf1_old_gwas"
ph_tmp$array<-"illumina370"

pheno<-ph_tmp


# belgium_inf2_old_gwas

ph_tmp<-read.csv("/path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/UK_IIBDGC/ibdgc_other_cohorts/gwas/belgian/InfosAllIndividualsInfinium2.csv",head=T)
ph_tmp<-ph_tmp[,c("Inf2_IndividualID","Inf2_Phenotype")]
ph_tmp<-ph_tmp[,c(1,1,2)]
table(ph_tmp$Inf2_Phenotype)
# CD healthy      UC
# 159     111       2

ph_tmp$ibd<-NA
ph_tmp$ibd[which(ph_tmp$Inf2_Phenotype=="healthy")]<-"0"
ph_tmp$ibd[which(ph_tmp$Inf2_Phenotype %in% c("CD","UC"))]<-"1"

ph_tmp$cd<-NA
ph_tmp$cd[which(ph_tmp$Inf2_Phenotype=="healthy")]<-"0"
ph_tmp$cd[which(ph_tmp$Inf2_Phenotype %in% c("CD"))]<-"1"

ph_tmp$uc<-NA
ph_tmp$uc[which(ph_tmp$Inf2_Phenotype=="healthy")]<-"0"
ph_tmp$uc[which(ph_tmp$Inf2_Phenotype %in% c("UC"))]<-"1"

table(ph_tmp$ibd)
table(ph_tmp$cd)
table(ph_tmp$uc)

ph_tmp<-ph_tmp[,c(1:2,4:6)]
colnames(ph_tmp)<-c("FID","IID","ibd","cd","uc")

ph_tmp$study<-"belgium_inf2_old_gwas"
ph_tmp$array<-"illumina370"

pheno<-rbind(pheno,ph_tmp)


# swedish_uc_old_gwas

ph_tmp<-read.table(paste(path,"pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18.fam",sep=""),head=F)
ph_tmp<-ph_tmp[,c(1,1,6)] # force same FID as IID
table(ph_tmp$V6)
# 1   2
# 341 923

ph_tmp$ibd<-NA
ph_tmp$ibd[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$ibd[which(ph_tmp$V6==2)]<-"1"

ph_tmp$cd<-NA

ph_tmp$uc<-NA
ph_tmp$uc[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$uc[which(ph_tmp$V6=="2")]<-"1"

table(ph_tmp$ibd)
table(ph_tmp$cd)
table(ph_tmp$uc)

ph_tmp<-ph_tmp[,c(1:2,4:6)]
colnames(ph_tmp)<-c("FID","IID","ibd","cd","uc")

ph_tmp$study<-"swedish_uc_old_gwas"
ph_tmp$array<-"illumina370"

pheno<-rbind(pheno,ph_tmp)


# niddk_old_gwas

ph_tmp<-read.table(paste(path,"pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17.fam",sep=""),head=F)
ph_tmp<-ph_tmp[,c(1,1,6)] # force same FID as IID
table(ph_tmp$V6)
# 1   2
# 959 833

ph_tmp$ibd<-NA
ph_tmp$ibd[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$ibd[which(ph_tmp$V6=="2")]<-"1"

ph_tmp$cd<-NA
ph_tmp$cd[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$cd[which(ph_tmp$V6=="2")]<-"1"

ph_tmp$uc<-NA
ph_tmp$uc[which(ph_tmp$V6=="1")]<-"0"
# ph_tmp$uc[which(ph_tmp$V6=="2")]<-"1"

table(ph_tmp$ibd)
table(ph_tmp$cd)
table(ph_tmp$uc)

ph_tmp<-ph_tmp[,c(1:2,4:6)]
colnames(ph_tmp)<-c("FID","IID","ibd","cd","uc")

ph_tmp$study<-"niddk_old_gwas"
ph_tmp$array<-"illumina370"

pheno<-rbind(pheno,ph_tmp)

# confirmed by Phil that all are cases:
ph_tmp<-read.table(paste(path,"pre_imputation/QC/niddk_uc_old_gwas/niddk_uc_old_gwas_hg18.fam",sep=""),head=F)
ph_tmp<-ph_tmp[,c(1,1,6)] # force same FID as IID
table(ph_tmp$V6)
# 1   2
# 552 500

ph_tmp$ibd<-NA
ph_tmp$ibd[which(ph_tmp$V6=="1")]<-"1"
ph_tmp$ibd[which(ph_tmp$V6=="2")]<-"1"

ph_tmp$cd<-NA
# ph_tmp$cd[which(ph_tmp$V6=="1")]<-"0"
# ph_tmp$cd[which(ph_tmp$V6=="2")]<-"1"

ph_tmp$uc<-NA
ph_tmp$uc[which(ph_tmp$V6=="1")]<-"1"
ph_tmp$uc[which(ph_tmp$V6=="2")]<-"1"

table(ph_tmp$ibd)
table(ph_tmp$cd)
table(ph_tmp$uc)

ph_tmp<-ph_tmp[,c(1:2,4:6)]
colnames(ph_tmp)<-c("FID","IID","ibd","cd","uc")

ph_tmp$study<-"niddk_old_gwas"
ph_tmp$array<-"illumina370"

pheno<-rbind(pheno,ph_tmp)

##############################################
# affymetrix6<-c("german_affy6_old_gwas","norway_affy6_old_gwas","gwas2")

# german_affy6_old_gwas

ph_tmp<-read.table(paste(path,"pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_hg18.fam",sep=""),head=F)
ph_tmp<-ph_tmp[,c(2,2,6)] # force same FID as IID
table(ph_tmp$V6)
# 1    2
# 1770 1057

ph_tmp$ibd<-NA
ph_tmp$ibd[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$ibd[which(ph_tmp$V6=="2")]<-"1"

ph_tmp$cd<-NA
# ph_tmp$cd[which(ph_tmp$V6=="1")]<-"0"
# ph_tmp$cd[which(ph_tmp$V6=="2")]<-"1"

ph_tmp$uc<-NA
ph_tmp$uc[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$uc[which(ph_tmp$V6=="2")]<-"1"

table(ph_tmp$ibd)
table(ph_tmp$cd)
table(ph_tmp$uc)

ph_tmp<-ph_tmp[,c(1:2,4:6)]
colnames(ph_tmp)<-c("FID","IID","ibd","cd","uc")

ph_tmp$study<-"german_affy6_old_gwas"
ph_tmp$array<-"affymetrix6"

pheno<-rbind(pheno,ph_tmp)


# norway_affy6_old_gwas

ph_tmp<-read.table(paste(path,"pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg18.fam",sep=""),head=F)
ph_tmp<-ph_tmp[,c(2,2,6)] # force same FID as IID
table(ph_tmp$V6)
# 1   2
# 282 268

ph_tmp$ibd<-NA
ph_tmp$ibd[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$ibd[which(ph_tmp$V6=="2")]<-"1"

ph_tmp$cd<-NA
# ph_tmp$cd[which(ph_tmp$V6=="1")]<-"0"
# ph_tmp$cd[which(ph_tmp$V6=="2")]<-"1"

ph_tmp$uc<-NA
ph_tmp$uc[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$uc[which(ph_tmp$V6=="2")]<-"1"

table(ph_tmp$ibd)
table(ph_tmp$cd)
table(ph_tmp$uc)

ph_tmp<-ph_tmp[,c(1:2,4:6)]
colnames(ph_tmp)<-c("FID","IID","ibd","cd","uc")

ph_tmp$study<-"norway_affy6_old_gwas"
ph_tmp$array<-"affymetrix6"

pheno<-rbind(pheno,ph_tmp)



# gwas2 - update new data with latest phenotype data freeze - use already the data in unix folder ukibdgwas:

ph_ukibdgc_update<-fread("/path/to/ibdgwas_ukibdgc/phenotype_data/phenotype_data_freezes/20221101/phenotype_data_UKIBDGCGWAS_20221101.tsv.gz",head=T)

ph_tmp<-read.table(paste(path,"/pre_imputation/QC/gwas2/gwas2_hg19.fam",sep=""),head=F)
ph_tmp<-ph_tmp[,c(2,2,6)] # force same FID as IID
table(ph_tmp$V6)
# 1    2
# 5417 2361

ph_tmp$ibd<-NA
ph_tmp$ibd[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$ibd[which(ph_tmp$V6=="2")]<-"1"

ph_tmp$cd<-NA
# ph_tmp$cd[which(ph_tmp$V6=="1")]<-"0"
# ph_tmp$cd[which(ph_tmp$V6=="2")]<-"1"

ph_tmp$uc<-NA
ph_tmp$uc[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$uc[which(ph_tmp$V6=="2")]<-"1"

dim(ph_ukibdgc_update[which(ph_ukibdgc_update$orig_sample_id %in% ph_tmp$V2),])
# [1] 7778   46
dim(ph_tmp)
# [1] 7778    6

table(ph_tmp$ibd)
# 0    1 
# 5417 2361 
table(ph_tmp$cd)
# < table of extent 0 >
table(ph_tmp$uc)
# 0    1 
# 5417 2361


# re-do with the updated phenotypes

ph_tmp$ibd<-NA
ph_tmp$ibd[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$ibd[which(ph_tmp$V2 %in% 
                   ph_ukibdgc_update$orig_sample_id[which(ph_ukibdgc_update$diag %in% c("Crohn's Disease","Indeterminate","Ulcerative Colitis"))])]<-"1"
table(ph_tmp$ibd)
# 0    1 
# 5417 2353

ph_tmp$cd<-NA
ph_tmp$cd[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$cd[which(ph_tmp$V2 %in% 
                   ph_ukibdgc_update$orig_sample_id[which(ph_ukibdgc_update$diag %in% c("Crohn's Disease"))])]<-"1"
table(ph_tmp$cd)
# 0    1 
# 5417   56 

ph_tmp$uc<-NA
ph_tmp$uc[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$uc[which(ph_tmp$V2 %in% 
                  ph_ukibdgc_update$orig_sample_id[which(ph_ukibdgc_update$diag %in% c("Ulcerative Colitis"))])]<-"1"
table(ph_tmp$uc)
# 0    1 
# 5417 2294 


ph_tmp<-ph_tmp[,c(1:2,4:6)]
colnames(ph_tmp)<-c("FID","IID","ibd","cd","uc")

ph_tmp$study<-"gwas2"
ph_tmp$array<-"affymetrix6"

pheno<-rbind(pheno,ph_tmp)

##############################################

affymetrix500<-c("gwas1")

ph_tmp<-read.table(paste(path,"pre_imputation/QC/gwas1/gwas1_hg19_edited.fam",sep=""),head=F)
ph_tmp<-ph_tmp[,c(2,2,6)] # force same FID as IID
table(ph_tmp$V6)
# 1    2
# 2936 1748

ph_tmp$ibd<-NA
ph_tmp$ibd[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$ibd[which(ph_tmp$V6=="2")]<-"1"

ph_tmp$cd<-NA
ph_tmp$cd[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$cd[which(ph_tmp$V6=="2")]<-"1"

ph_tmp$uc<-NA
# ph_tmp$uc[which(ph_tmp$V6=="1")]<-"0"
# ph_tmp$uc[which(ph_tmp$V6=="2")]<-"1"

table(ph_tmp$ibd)
# 0    1 
# 2936 1748 
table(ph_tmp$cd)
# 0    1 
# 2936 1748 
table(ph_tmp$uc)
# < table of extent 0 >

# re-do with the updated phenotypes

ph_tmp$ibd<-NA
ph_tmp$ibd[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$ibd[which(ph_tmp$V2 %in% 
                   ph_ukibdgc_update$orig_sample_id[which(ph_ukibdgc_update$diag %in% c("Crohn's Disease","Indeterminate","Ulcerative Colitis"))])]<-"1"
table(ph_tmp$ibd)
# 0    1 
# 2936 1744 

ph_tmp$cd<-NA
ph_tmp$cd[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$cd[which(ph_tmp$V2 %in% 
                  ph_ukibdgc_update$orig_sample_id[which(ph_ukibdgc_update$diag %in% c("Crohn's Disease"))])]<-"1"
table(ph_tmp$cd)
# 0    1 
# 2936 1692 

ph_tmp$uc<-NA
ph_tmp$uc[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$uc[which(ph_tmp$V2 %in% 
                  ph_ukibdgc_update$orig_sample_id[which(ph_ukibdgc_update$diag %in% c("Ulcerative Colitis"))])]<-"1"
table(ph_tmp$uc)
# 0    1 
# 2936   51 


ph_tmp<-ph_tmp[,c(1:2,4:6)]
colnames(ph_tmp)<-c("FID","IID","ibd","cd","uc")

ph_tmp$study<-"gwas1"
ph_tmp$array<-"affymetrix500"

pheno<-rbind(pheno,ph_tmp)


##############################################

humancoreexome<-c("all_hce")

ph_tmp_g3<-read.table(paste(path,"pre_imputation/QC/gwas3/gwas3_hg19_edited.fam",sep=""),head=F)
pheg3<-fread("/path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/UK_IIBDGC/coreex_gaibdc_phenotypes.tsv",head=T)

dim(pheg3[which(pheg3$orig_sample_id %in% ph_tmp_g3$V1),])
# [1] 11768    45
dim(pheg3)
# [1] 11930    45

table(ph_tmp_g3$V6)
# 1     2
# 10484 11768

pheg3<-pheg3[which(pheg3$orig_sample_id %in% ph_tmp_g3$V1),]
table(pheg3$diag)
# Crohn's Disease      Indeterminate Ulcerative Colitis            Unknown
#            5923                 19               5452                374

ph_tmp_g3<-merge(ph_tmp_g3,pheg3,by.x="V1",by.y="orig_sample_id",all.x=T)

ph_tmp_g3$ibd<-NA
ph_tmp_g3$ibd[which(ph_tmp_g3$V6=="1")]<-"0"
ph_tmp_g3$ibd[which(ph_tmp_g3$V6=="2")]<-"1"

ph_tmp_g3$cd<-NA
ph_tmp_g3$cd[which(ph_tmp_g3$V6=="1")]<-"0"
ph_tmp_g3$cd[which(ph_tmp_g3$diag=="Crohn's Disease")]<-"1"

ph_tmp_g3$uc<-NA
ph_tmp_g3$uc[which(ph_tmp_g3$V6=="1")]<-"0"
ph_tmp_g3$uc[which(ph_tmp_g3$diag=="Ulcerative Colitis")]<-"1"

ph_tmp_g3<-ph_tmp_g3[,c("V1","V2","V6","ibd","cd","uc")]

####

ph_tmp_nw<-read.table(paste(path,"pre_imputation/QC/all_hce/new_wave_hg19_noind.fam",sep=""),head=F)

sample1<-read.table("/path/to/project",head=T)
sample1<-sample1[-1,] # remove subheader
table(sample1$bin1)
sample2<-read.table("/path/to/project",head=T)
sample2<-sample2[-1,] # remove subheader
table(sample2$bin1)
sample3<-read.table("/path/to/project",head=T)
sample3<-sample3[-1,] # remove subheader

ph_tmp_nw$ibd<-NA
ph_tmp_nw$ibd[which(ph_tmp_nw$V6=="1")]<-"0"
ph_tmp_nw$ibd[which(ph_tmp_nw$V6=="2")]<-"1"

ph_tmp_nw$cd<-NA
ph_tmp_nw$cd[which(ph_tmp_nw$V6=="1")]<-"0"
ph_tmp_nw$cd[which(ph_tmp_nw$V1 %in% sample2$ID_1[which(sample2$bin1==1)])]<-"1"

ph_tmp_nw$uc<-NA
ph_tmp_nw$uc[which(ph_tmp_nw$V6=="1")]<-"0"
ph_tmp_nw$uc[which(ph_tmp_nw$V1 %in% sample2$ID_1[which(sample3$bin1==1)])]<-"1"

table(ph_tmp_nw$ibd)
table(ph_tmp_nw$cd)
table(ph_tmp_nw$uc)

# update missing ones with:

phen_new_wave<-fread("/path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/UK_IIBDGC/coreex_gaibdc_phenotypes.tsv",head=T,sep="\t")
phen_nw<-fread("/path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/UK_IIBDGC/coreex_idbgwas_Final_phenotypes.tsv",head=T,sep="\t")
phen<-rbind(phen_new_wave,phen_nw)

table(ph_tmp_nw$V6)
# 2
# 2099

ph_tmp_nw$id<-gsub("urn:wtsi:","",ph_tmp_nw$V1)
ph_tmp_nw<-merge(ph_tmp_nw,phen,by.x="id",by.y="orig_sample_id",all.x=T)
table(ph_tmp_nw$V6,ph_tmp_nw$diag,useNA="ifany")
#   Crohn's Disease Indeterminate Ulcerative Colitis Unknown
# 2             281            93                303    1422

table(ph_tmp_nw$cd)
# 1
# 744

table(ph_tmp_nw$cd[which(ph_tmp_nw$diag=="Crohn's Disease")],useNA="ifany")
# 1 <NA>
# 175  106

ph_tmp_nw$cd[which(ph_tmp_nw$diag=="Crohn's Disease")]<-"1"
table(ph_tmp_nw$cd)
# 1
# 850


table(ph_tmp_nw$uc)
# 1
# 372

table(ph_tmp_nw$uc[which(ph_tmp_nw$diag=="Ulcerative Colitis")],useNA="ifany")
# 1 <NA>
# 251   52

ph_tmp_nw$uc[which(ph_tmp_nw$diag=="Ulcerative Colitis")]<-"1"
table(ph_tmp_nw$uc,useNA="ifany")
# 1 <NA>
# 424 1675

ph_tmp_nw<-ph_tmp_nw[,c("V1","V2","V6","ibd","cd","uc")]

ph_tmp<-rbind(ph_tmp_nw,ph_tmp_g3)

table(ph_tmp$ibd)
# 0     1 
# 10484 13867 
table(ph_tmp$cd)
# 0     1 
# 10484  6773 
table(ph_tmp$uc)
# 0     1 
# 10484  5876 


# re-do with the updated phenotypes

ph_tmp$ibd<-NA
ph_tmp$ibd[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$ibd[which(ph_tmp$V2 %in% 
                   ph_ukibdgc_update$orig_sample_id[which(ph_ukibdgc_update$diag %in% c("Crohn's Disease","Indeterminate","Ulcerative Colitis"))])]<-"1"
table(ph_tmp$ibd)
# 0     1 
# 10484 13106


ph_tmp$cd<-NA
ph_tmp$cd[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$cd[which(ph_tmp$V2 %in% 
                  ph_ukibdgc_update$orig_sample_id[which(ph_ukibdgc_update$diag %in% c("Crohn's Disease"))])]<-"1"
table(ph_tmp$cd)
# 0     1 
# 10484  7071 


ph_tmp$uc<-NA
ph_tmp$uc[which(ph_tmp$V6=="1")]<-"0"
ph_tmp$uc[which(ph_tmp$V2 %in% 
                  ph_ukibdgc_update$orig_sample_id[which(ph_ukibdgc_update$diag %in% c("Ulcerative Colitis"))])]<-"1"
table(ph_tmp$uc)
# 0     1 
# 10484  5765 

ph_tmp<-ph_tmp[,c(1:2,4:6)]
colnames(ph_tmp)<-c("FID","IID","ibd","cd","uc")

ph_tmp$study<-"all_hce"
ph_tmp$array<-"humancoreexome"

pheno<-rbind(pheno,ph_tmp)

table(pheno$study,pheno$ibd)
#                           0     1
# all_hce               10484 13106
# belgium_inf1_old_gwas   900   517
# belgium_inf2_old_gwas   111   161
# german_affy6_old_gwas  1770  1057
# gwas1                  2936  1744
# gwas2                  5417  2353
# niddk_old_gwas          959  1885
# norway_affy6_old_gwas   282   268
# swedish_uc_old_gwas     341   923

table(pheno$study,pheno$cd)
#                           0     1
# all_hce               10484  7071
# belgium_inf1_old_gwas   900   517
# belgium_inf2_old_gwas   111   159
# german_affy6_old_gwas     0     0
# gwas1                  2936  1692
# gwas2                  5417    56
# niddk_old_gwas          959   833
# norway_affy6_old_gwas     0     0
# swedish_uc_old_gwas       0     0


table(pheno$study,pheno$uc)
#                           0     1
# all_hce               10484  5765
# belgium_inf1_old_gwas     0     0
# belgium_inf2_old_gwas   111     2
# german_affy6_old_gwas  1770  1057
# gwas1                  2936    51
# gwas2                  5417  2294
# niddk_old_gwas          959  1052
# norway_affy6_old_gwas   282   268
# swedish_uc_old_gwas     341   923


table(pheno$array,pheno$ibd)
#                    0     1
# affymetrix500   2936  1744
# affymetrix6     7469  3678
# humancoreexome 10484 13106
# illumina370     2311  3486

table(pheno$array,pheno$cd)
#                    0     1
# affymetrix500   2936  1692
# affymetrix6     5417    56
# humancoreexome 10484  7071
# illumina370     1970  1509

table(pheno$array,pheno$uc)
#                    0     1
# affymetrix500   2936    51
# affymetrix6     7469  3619
# humancoreexome 10484  5765
# illumina370     1411  1977


##############################################

# pheno data for studies below in core pheno file:

# set to missing data from finish studies - due to overlap with finnGen, plus old studies from cedars - for which we do not have appropiate Ctr genotyped with same array

studies_2<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
             "kiel_austria_sibdcs_gsa","netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa",
             "niddk_feinstein_gsa","basque_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa",
             "belgium_vermeire_gsa","prism_nfe_gsa","prism_nfe_gwas","finland_illugwas","german_affy6_old_gwas",
             "norway_affy6_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
             "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas",
             "swedish_uc_old_gwas",
             "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
             "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
             "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
             "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
             "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
             "xavier_share_gsa")


studies_2[which(!studies_2 %in% studies)]
# "cedars_370k_old_gwas"
# "cedars_610k_old_gwas"     
# "cedars_omni_old_gwas"


# exclude pheno from samples recruited at JH
jh<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/hopkins_gsa_exclude.csv",head=T)



## latest udpate phenotype file - sent by Phil on 18/Jan/2022

ph<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes-827da7e_updated.csv.gz",head=T)


ph$ibd<-NA
ph$ibd[which(ph$control=="1")]<-"0"
ph$ibd[which(ph$affection=="Affected")]<-"1"

ph$cd<-NA
ph$cd[which(ph$control=="1")]<-"0"
ph$cd[which(ph$diag=="Crohn's Disease")]<-"1"

ph$uc<-NA
ph$uc[which(ph$control=="1")]<-"0"
ph$uc[which(ph$diag=="Ulcerative Colitis")]<-"1"

table(as.character(ph$batch[which(is.na(ph$ibd))]))
# anvil_ccdg_broad_ai_ibd_daly_bernstein_gsa 
# 23 
# anvil_ccdg_broad_ai_ibd_daly_brant_niddk_gsa 
# 23 
# anvil_ccdg_broad_ai_ibd_daly_cho_niddk_gsa 
# 10 
# anvil_ccdg_broad_ai_ibd_daly_franchimont_gsa 
# 37 
# anvil_ccdg_broad_ai_ibd_daly_lewis_sparc_gsa 
# 4 
# anvil_ccdg_broad_ai_ibd_daly_mccauley_gsa 
# 14 
# anvil_ccdg_broad_ai_ibd_daly_mcgovern_gsa 
# 126 
# anvil_ccdg_broad_ai_ibd_daly_moayyedi_imagine_gsa 
# 1422 
# anvil_ccdg_broad_ai_ibd_daly_newberry_share_gsa 
# 4 
# anvil_ccdg_broad_ai_ibd_daly_pekow_share_gsa 
# 72 
# anvil_ccdg_broad_ai_ibd_daly_rioux_igenomed_gsa 
# 4 
# anvil_ccdg_broad_ai_ibd_daly_sands_msccr_gsa 
# 89 
# anvil_ccdg_broad_ai_ibd_daly_vermeire_gsa 
# 138 
# anvil_ccdg_broad_ai_ibd_daly_xavier_share_gsa 
# 13 
# basque 
# 7 
# belgium_franchimont 
# 721 
# belgium_louis 
# 17 
# belgium_vermeire 
# 97 
# ccdg_broad_ai_ibd_daly_palotie_hus_gsa 
# 1 
# ccdg_broad_ai_ibd_daly_xavier_prism_gsa 
# 82 
# ccdg_hiesp_ibd_daly_weersma_ibd_gsa-array_v1_3 
# 289 
# ccdg_ibd_daly_stampfer_ibd_gsa-array_v1_4 
# 444 
# cedars 
# 27 
# finland 
# 155 
# Helmsley_PRISM_GWAS-Chip_B1 
# 94 
# Helmsley_Xavier_IBD-PRISM_WES_B2 
# 38 
# hiesp_broad_ibd_daly_farkkila_gsa 
# 41 
# lithuania 
# 24 
# niddk_broad 
# 16 
# niddk_feinstein 
# 410 
# prism_nfe 
# 16 
# prism_nfe_gwas 
# 12 
# sibdcs 
# 1 
# sweden_sicbio 
# 109 
# uk_gaibdc 
# 31 

ph<-ph[,c("sample_id","sample_id","ibd","cd","uc")]
colnames(ph)<-c("FID","IID","ibd","cd","uc")

ph$study<-NA
ph$array<-NA

array_tmp<-c("humanomniexpress","humanomni1","quad610","gsa","illuminaexome")


for (i in 1:length(array_tmp)) {

  print(array_tmp[i])
  study<-get(array_tmp[i])

  for (j in 1:length(study)){

    print(study[j])

    file.tmp<-paste(path,"pre_imputation/QC/",study[j],"/",study[j],"_hg19_edited.fam",sep="")

    if (file.exists(file.tmp)){
      tmp<-read.table(file.tmp,head=F)
    }else{
      tmp<-read.table(paste(path,"pre_imputation/QC/",study[j],"/",study[j],"_hg18_edited.fam",sep=""),head=F)
    }
    print(length(ph$study[which(ph$FID %in% tmp$V1)]))
    ph$study[which(ph$FID %in% tmp$V1)]<-study[j]
    ph$array[which(ph$FID %in% tmp$V1)]<-array_tmp[i]

  }

}

#############
# [1] "humanomniexpress"
# [1] "australia_omniexome"
# [1] 1310
# [1] "humanomni1"
# [1] "pittsburgh_gsa"
# [1] 2781
# [1] "quad610"
# [1] "spain_gsa"
# [1] 3443
# [1] "gsa"
# [1] "italy_gsa"
# [1] 1016
# [1] "kiel_austria_sibdcs_gsa"
# [1] 14661
# [1] "netherlands_gsa"
# [1] 4705
# [1] "slovenia_gsa"
# [1] 270
# [1] "sweden_gsa"
# [1] 1523
# [1] "niddk_broad_gsa"
# [1] 5531
# [1] "niddk_feinstein_gsa"
# [1] 8687
# [1] "basque_gsa"
# [1] 1523
# [1] "prism_nfe_gsa"
# [1] 482
# [1] "lithuania_gsa"
# [1] 2301
# [1] "belgium_louis_gsa"
# [1] 1548
# [1] "belgium_franchimont_gsa"
# [1] 2253
# [1] "belgium_vermeire_gsa"
# [1] 4111
# [1] "mccauley_gsa"
# [1] 788
# [1] "ccfa_gsa"
# [1] 2188
# [1] "cedars_gsa"
# [1] 3123
# [1] "bernstein_gsa"
# [1] 537
# [1] "farkkila_gsa"
# [1] 109
# [1] "franchimont_gsa"
# [1] 2826
# [1] "franke_gsa"
# [1] 885
# [1] "hyams_protect_gsa"
# [1] 418
# [1] "lewis_sparc_gsa"
# [1] 2863
# [1] "mccauley_new_gsa"
# [1] 1642
# [1] "mcgovern_gsa"
# [1] 6088
# [1] "moayyedi_imagine_gsa"
# [1] 2550
# [1] "newberry_share_gsa"
# [1] 869
# [1] "niddk_cho_gsa"
# [1] 1786
# [1] "niddk_duerr_gsa"
# [1] 1948
# [1] "niddk_rioux_gsa"
# [1] 919
# [1] "niddk_silverberg_gsa"
# [1] 2440
# [1] "palotie_hus_gsa"
# [1] 879
# [1] "pekow_share_gsa"
# [1] 706
# [1] "rioux_igenomed_gsa"
# [1] 186
# [1] "sands_msccr_gsa"
# [1] 1464
# [1] "stampfer_gsa"
# [1] 1921
# [1] "vermeire_gsa"
# [1] 4851
# [1] "weersma_gsa"
# [1] 998
# [1] "xavier_prism_gsa"
# [1] 774
# [1] "xavier_share_gsa"
# [1] 709
# [1] "illuminaexome"
# [1] "prism_nfe_gwas"
# [1] 874
# [1] "helmsley_prism_gsa"
# [1] 874
# [1] "helmsley_xavier_prism_gsa"
# [1] 1335
# [1] "finland_illugwas"
# [1] 612

#############

pheno<-rbind(pheno,ph)

table(pheno$array)
# affymetrix500      affymetrix6              gsa   humancoreexome 
#          4684            11155            93078            24351 
# humanomni1 humanomniexpress      illumina370    illuminaexome 
#       2781             1310             5797             3695 
# quad610 
#    3443 


write.table(pheno,paste(path,"pheno/phenotype_data_for_all_cohorts_Dec2022_analysis.tsv",sep=""),col.names=T,row.names=F,sep="\t")


ii<-1

# keep only nonDuplicated samples - see step 25

# ccfa_gsa - all samples already included in lewis_sparc

for (j in 1:length(studies)) {
  
  if(studies[j]!="ccfa_gsa") {
    nodup_tmp<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam",sep=""),head=F)
    
    nodup_tmp$study<-studies[j]
    if(j==1) {
      nodup<-nodup_tmp
    }else{
      nodup<-rbind(nodup,nodup_tmp)
    }
    rm(nodup_tmp) 
  }
  
}
dim(nodup)
# [1] 111638      7


# get first degree related samples  - see above

rel<-read.table(paste(path,"post_imputation/2022/analysis/list_related_samples_across_arrays_toexclude_",ancestry[ii],sep=""),head=T)
table(rel$array)
# affymetrix500    affymetrix6 humancoreexome     humanomni1    illumina370 
# 54             38             20             10            101 
# illuminaexome        quad610 
# 30              1 


# get updated sex - see step 24

sex<-fread(paste(path,"pheno/gwas-mega-2-core_updated_sex_gender_Dec22.txt.gz",sep=""),head=T)
sex$sex<-sex$SNPSEX2
sex$FID<-sex$IID

sex<-sex[,c("FID","sex")]
# for swedish_uc_old_gwas - keep sex as in fam file - no chrX or chrY data to infer


tmp<-read.table(paste(path,"pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg18.fam",sep=""),head=F)
tmp<-tmp[,c("V1","V5")]
colnames(tmp)<-colnames(sex)
sex<-rbind(sex,tmp)


# set to missing samples recruited at JH:
table(pheno$study[which(pheno$FID %in% jh$V1)])
# niddk_broad_gsa niddk_feinstein_gsa      niddk_old_gwas 
# 154                 114                  81


for (i in 1:length(array)) {

  print(array[i])

  # see step 27
  tmp<-read.table(paste(path,"post_imputation/2022/",array[i],"/genotyped_data/",array[i],"_all_studies_merged_",ancestry[ii],"_pruned.fam",sep=""),head=F)
  
  # Number EUR samples
  print(paste("EUR samples: ",nrow(tmp),sep=""))
  
  phen_tmp<-pheno[which(pheno$FID %in% tmp$V1),]
  phen_tmp<-merge(phen_tmp,sex[,c("FID","sex")],by="FID",sort=F,all.x=T,order=F)

  phen_tmp<-phen_tmp[match(tmp$V1,phen_tmp$FID),]

  print(table(as.character(phen_tmp$FID)==as.character(tmp$V1)))
  print(table(nrow(phen_tmp)==nrow(tmp)))

  # print(table(phen_tmp$ibd))
  # print(table(phen_tmp$cd))
  # print(table(phen_tmp$uc))
  # print(table(phen_tmp$sex))
  
  # samples with no pheno data:
  # if(array[i] %in% array_tmp) {
    print(paste("N samples with no pheno: ",nrow(phen_tmp[which(is.na(phen_tmp$ibd) | phen_tmp$ibd==""),]),sep=""))
  # }
  
  # set to missing pheno from samples from finish studies - likely to be included in FinnGen
  print(paste("N Finish Samples: ",nrow(phen_tmp[which(phen_tmp$study %in% c("finland_illugwas","farkkila_gsa","palotie_hus_gsa")),]),sep=""))
  
  phen_tmp$ibd[which(phen_tmp$study %in% c("finland_illugwas","farkkila_gsa","palotie_hus_gsa"))]<-NA
  phen_tmp$cd[which(phen_tmp$study %in% c("finland_illugwas","farkkila_gsa","palotie_hus_gsa"))]<-NA
  phen_tmp$uc[which(phen_tmp$study %in% c("finland_illugwas","farkkila_gsa","palotie_hus_gsa"))]<-NA
  phen_tmp$sex[which(phen_tmp$study %in% c("finland_illugwas","farkkila_gsa","palotie_hus_gsa"))]<-NA
  

  # duplicated samples:
  print(paste("N Duplicated Samples: ",nrow(phen_tmp[which(!phen_tmp$IID %in% nodup$V2),]),sep=""))
  
  # set to missing pheno from samples not in no-duplicates EUR samples file"
  phen_tmp$ibd[which(!phen_tmp$IID %in% nodup$V2)]<-NA
  phen_tmp$cd[which(!phen_tmp$IID %in% nodup$V2)]<-NA
  phen_tmp$uc[which(!phen_tmp$IID %in% nodup$V2)]<-NA
  phen_tmp$sex[which(!phen_tmp$IID %in% nodup$V2)]<-NA
  
  # set to missing pheno from samples collected at John Hopkins"
  print(paste("N JH samples: ",nrow(phen_tmp[which(phen_tmp$IID %in% jh$V1),]),sep=""))
  print(table(phen_tmp$study[which(phen_tmp$IID %in% jh$V1)]))
  
  phen_tmp$ibd[which(phen_tmp$IID %in% jh$V1)]<-NA
  phen_tmp$cd[which(phen_tmp$IID %in% jh$V1)]<-NA
  phen_tmp$uc[which(phen_tmp$IID %in% jh$V1)]<-NA
  phen_tmp$sex[which(phen_tmp$IID %in% jh$V1)]<-NA
  

  # set to missing pheno from samples in rel file:
  print(paste("N relatives to exclude ",array[i],": ",nrow(phen_tmp[which(phen_tmp$IID %in% rel$V1),]),sep=""))

  phen_tmp$ibd[which(phen_tmp$IID %in% rel$V1)]<-NA
  phen_tmp$cd[which(phen_tmp$IID %in% rel$V1)]<-NA
  phen_tmp$uc[which(phen_tmp$IID %in% rel$V1)]<-NA
  phen_tmp$sex[which(phen_tmp$IID %in% rel$V1)]<-NA
  
  # print("SEX study")
  # print(table(phen_tmp$sex,phen_tmp$study))

  # set to missing pheno from samples where sex = 0:
  print(paste("N sex=0 to exclude ",array[i],": ",nrow(phen_tmp[which(phen_tmp$sex==0),]),sep=""))
  
  print("N sex=0 per sex ")
  print(table(phen_tmp$ibd[which(phen_tmp$sex==0)]))

  phen_tmp$ibd[which(phen_tmp$sex==0)]<-NA
  phen_tmp$cd[which(phen_tmp$sex==0)]<-NA
  phen_tmp$uc[which(phen_tmp$sex==0)]<-NA
  phen_tmp$sex[which(phen_tmp$sex==0)]<-NA
  
  phen_tmp$sex[which(is.na(phen_tmp$ibd))]<-NA

  print(paste("N samples to include in analysis ",array[i],": ",nrow(phen_tmp[which(!is.na(phen_tmp$ibd)),]),sep=""))

  print("IBD")
  print(table(phen_tmp$ibd))
  print("CD")
  print(table(phen_tmp$cd))
  print("UC")
  print(table(phen_tmp$uc))
  print("SEX")
  print(table(phen_tmp$sex))

  print("SEX IBD")
  print(table(phen_tmp$sex,phen_tmp$ibd,useNA="ifany"))

  # print(table(phen_tmp$array[which(!is.na(phen_tmp$ibd))]))
  print(table(phen_tmp$study[which(!is.na(phen_tmp$ibd))]))

  print(table(as.character(phen_tmp$FID)==as.character(tmp$V1)))

  # save phenotype data:
  write.table(phen_tmp[,c("FID","IID","ibd")],paste(path,"post_imputation/2022/",array[i],"/phenotype_data/",array[i],"_all_studies_merged_",ancestry[ii],"_phenotype_ibd",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  write.table(phen_tmp[,c("FID","IID","cd")],paste(path,"post_imputation/2022/",array[i],"/phenotype_data/",array[i],"_all_studies_merged_",ancestry[ii],"_phenotype_cd",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  write.table(phen_tmp[,c("FID","IID","uc")],paste(path,"post_imputation/2022/",array[i],"/phenotype_data/",array[i],"_all_studies_merged_",ancestry[ii],"_phenotype_uc",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

  if(i==1) {
    final_phen<-phen_tmp
  }else{
    final_phen<-rbind(final_phen,phen_tmp)
  }
  
  # create covariate file with sex to adjust for, no NA allowed:
  phen_tmp<-phen_tmp[which(!is.na(phen_tmp$ibd)),]
  write.table(phen_tmp[,c("FID","IID","sex")],paste(path,"post_imputation/2022/",array[i],"/phenotype_data/",array[i],"_all_studies_merged_",ancestry[ii],"_covariate_sex",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

}

###################
# 
# [1] "illumina370"
# [1] "EUR samples: 5633"
# 
# TRUE 
# 5633 
# 
# TRUE 
# 1 
# [1] "N samples with no pheno: 0"
# [1] "N Finish Samples: 0"
# [1] "N Duplicated Samples: 2181"
# [1] "N JH samples: 80"
# 
# niddk_old_gwas 
# 80 
# [1] "N relatives to exclude illumina370: 101"
# [1] "N sex=0 to exclude illumina370: 9"
# [1] "N sex=0 per sex "
# 
# 0 1 
# 8 1 
# [1] "N samples to include in analysis illumina370: 3263"
# [1] "IBD"
# 0    1 
# 1852 1411 
# [1] "CD"
# 0    1 
# 1536  425 
# [1] "UC"
# 0    1 
# 1040  986 
# [1] "SEX"
# 1    2 
# 2046 1217 
# [1] "SEX IBD"
#         0    1 <NA>
# 1    1296  750    0
# 2     556  661    0
# <NA>    0    0 2370
# 
# belgium_inf1_old_gwas belgium_inf2_old_gwas        niddk_old_gwas 
# 986                   152                   907 
# swedish_uc_old_gwas 
# 1218 
# 
# TRUE 
# 5633 


# [1] "affymetrix6"
# [1] "EUR samples: 11087"
# 
# TRUE 
# 11087 
# 
# TRUE 
# 1 
# [1] "N samples with no pheno: 8"
# [1] "N Finish Samples: 0"
# [1] "N Duplicated Samples: 2841"
# [1] "N JH samples: 0"
# < table of extent 0 >
#   [1] "N relatives to exclude affymetrix6: 38"
# [1] "N sex=0 to exclude affymetrix6: 21"
# [1] "N sex=0 per sex "
# 0  1 
# 10 11 
# [1] "N samples to include in analysis affymetrix6: 8180"
# [1] "IBD"
# 0    1 
# 4643 3537 
# [1] "CD"
# 0    1 
# 2762   56 
# [1] "UC"
# 0    1 
# 4643 3478 
# [1] "SEX"
# 
# 1    2 
# 4093 4087 
# [1] "SEX IBD"
# 
# 0    1 <NA>
#   1    2443 1650    0
# 2    2200 1887    0
# <NA>    0    0 2907
# 
# german_affy6_old_gwas                 gwas2 norway_affy6_old_gwas 
#                  2546                  5093                   541 
# 
# TRUE 
# 11087 


# [1] "humanomniexpress"
# [1] "EUR samples: 1245"
# 
# TRUE 
# 1245 
# 
# TRUE 
# 1 
# [1] "N samples with no pheno: 0"
# [1] "N Finish Samples: 0"
# [1] "N Duplicated Samples: 7"
# [1] "N JH samples: 0"
# < table of extent 0 >
#   [1] "N relatives to exclude humanomniexpress: 0"
# [1] "N sex=0 to exclude humanomniexpress: 3"
# [1] "N sex=0 per sex "
# 
# 1 
# 3 
# [1] "N samples to include in analysis humanomniexpress: 1235"
# [1] "IBD"
# 0   1 
# 589 646 
# [1] "CD"
# 0 
# 589 
# [1] "UC"
# 
# 0   1 
# 589 646 
# [1] "SEX"
# 1   2 
# 578 657 
# [1] "SEX IBD"
# 
# 0   1 <NA>
#   1    251 327    0
# 2    338 319    0
# <NA>   0   0   10
# 
# australia_omniexome 
# 1235 
# 
# TRUE 
# 1245 


# [1] "affymetrix500"
# [1] "EUR samples: 4652"
# 
# TRUE 
# 4652 
# 
# TRUE 
# 1 
# [1] "N samples with no pheno: 3"
# [1] "N Finish Samples: 0"
# [1] "N Duplicated Samples: 4"
# [1] "N JH samples: 0"
# < table of extent 0 >
#   [1] "N relatives to exclude affymetrix500: 54"
# [1] "N sex=0 to exclude affymetrix500: 0"
# [1] "N sex=0 per sex "
# < table of extent 0 >
#   [1] "N samples to include in analysis affymetrix500: 4591"
# [1] "IBD"
# 0    1 
# 2886 1705 
# [1] "CD"
# 0    1 
# 2886 1656 
# [1] "UC"
# 0    1 
# 2886   48 
# [1] "SEX"
# 1    2 
# 2088 2503 
# [1] "SEX IBD"
# 0    1 <NA>
#   1    1425  663    0
# 2    1461 1042    0
# <NA>    0    0   61
# 
# gwas1 
# 4591 
# 
# TRUE 
# 4652 


# [1] "humancoreexome"
# [1] "EUR samples: 22588"
# 
# TRUE 
# 22588 
# 
# TRUE 
# 1 
# [1] "N samples with no pheno: 694"
# [1] "N Finish Samples: 0"
# [1] "N Duplicated Samples: 1345"
# [1] "N JH samples: 0"
# < table of extent 0 >
#   [1] "N relatives to exclude humancoreexome: 20"
# [1] "N sex=0 to exclude humancoreexome: 62"
# [1] "N sex=0 per sex "
# 
# 0  1 
# 21 39 
# [1] "N samples to include in analysis humancoreexome: 20563"
# [1] "IBD"
# 0     1 
# 10313 10250 
# [1] "CD"
# 0     1 
# 10313  5489 
# [1] "UC"
# 0     1 
# 10313  4553 
# [1] "SEX"
# 1     2 
# 9417 11146 
# [1] "SEX IBD"
# 
# 0    1 <NA>
#   1    4551 4866    0
# 2    5762 5384    0
# <NA>    0    0 2025
# 
# all_hce 
# 20563 
# 
# TRUE 
# 22588 


# [1] "humanomni1"
# [1] "EUR samples: 2701"
# 
# TRUE 
# 2701 
# 
# TRUE 
# 1 
# [1] "N samples with no pheno: 0"
# [1] "N Finish Samples: 0"
# [1] "N Duplicated Samples: 648"
# [1] "N JH samples: 0"
# < table of extent 0 >
#   [1] "N relatives to exclude humanomni1: 10"
# [1] "N sex=0 to exclude humanomni1: 9"
# [1] "N sex=0 per sex "
# 
# 0 1 
# 8 1 
# [1] "N samples to include in analysis humanomni1: 2034"
# [1] "IBD"
# 0    1 
# 1273  761 
# [1] "CD"
# 0    1 
# 1273  386 
# [1] "UC"
# 0    1 
# 1273  369 
# [1] "SEX"
# 
# 1    2 
# 867 1167 
# [1] "SEX IBD"
# 
# 0   1 <NA>
#   1    500 367    0
# 2    773 394    0
# <NA>   0   0  667
# 
# pittsburgh_gsa 
# 2034 
# 
# TRUE 
# 2701 

# [1] "quad610"
# [1] "EUR samples: 3396"
# 
# TRUE 
# 3396 
# 
# TRUE 
# 1 
# [1] "N samples with no pheno: 0"
# [1] "N Finish Samples: 0"
# [1] "N Duplicated Samples: 0"
# [1] "N JH samples: 0"
# < table of extent 0 >
#   [1] "N relatives to exclude quad610: 1"
# [1] "N sex=0 to exclude quad610: 31"
# [1] "N sex=0 per sex "
# 
# 0  1 
# 7 24 
# [1] "N samples to include in analysis quad610: 3364"
# [1] "IBD"
# 0    1 
# 1453 1911 
# [1] "CD"
# 0    1 
# 1453 1130 
# [1] "UC"
# 0    1 
# 1453  781 
# [1] "SEX"
# 
# 1    2 
# 1859 1505 
# [1] "SEX IBD"
# 
# 0   1 <NA>
#   1    871 988    0
# 2    582 923    0
# <NA>   0   0   32
# 
# spain_gsa 
# 3364 
# 
# TRUE 
# 3396 

# [1] "gsa"
# [1] "EUR samples: 81334"
# 
# TRUE 
# 81334 
# 
# TRUE 
# 1 
# [1] "N samples with no pheno: 205"
# [1] "N Finish Samples: 921"
# [1] "N Duplicated Samples: 21825"
# [1] "N JH samples: 257"
# 
# niddk_broad_gsa niddk_feinstein_gsa 
# 144                 113 
# [1] "N relatives to exclude gsa: 0"
# [1] "N sex=0 to exclude gsa: 142"
# [1] "N sex=0 per sex "
# 0   1 
# 33 109 
# [1] "N samples to include in analysis gsa: 57998"
# [1] "IBD"
# 
# 0     1 
# 16730 41268 
# [1] "CD"
# 
# 0     1 
# 16730 25869 
# [1] "UC"
# 
# 0     1 
# 16730 13847 
# [1] "SEX"
# 
# 1     2 
# 28078 29920 
# [1] "SEX IBD"
# 
# 0     1  <NA>
#   1     8467 19611     0
# 2     8263 21657     0
# <NA>     0     0 23336
# 
# basque_gsa belgium_franchimont_gsa       belgium_louis_gsa 
# 1468                       3                    1473 
# belgium_vermeire_gsa           bernstein_gsa              cedars_gsa 
# 120                     439                      25 
# franchimont_gsa              franke_gsa       hyams_protect_gsa 
# 2194                     755                     318 
# italy_gsa kiel_austria_sibdcs_gsa         lewis_sparc_gsa 
# 929                   13541                    2416 
# lithuania_gsa            mccauley_gsa        mccauley_new_gsa 
# 2188                       4                    1293 
# mcgovern_gsa    moayyedi_imagine_gsa         netherlands_gsa 
# 4928                     984                    4266 
# newberry_share_gsa         niddk_broad_gsa           niddk_cho_gsa 
# 733                      29                    1166 
# niddk_duerr_gsa     niddk_feinstein_gsa         niddk_rioux_gsa 
# 1149                    5016                     590 
# niddk_silverberg_gsa         pekow_share_gsa           prism_nfe_gsa 
# 1680                     524                     387 
# rioux_igenomed_gsa         sands_msccr_gsa            slovenia_gsa 
# 151                     837                     261 
# stampfer_gsa              sweden_gsa            vermeire_gsa 
# 1438                    1340                    4502 
# weersma_gsa        xavier_prism_gsa        xavier_share_gsa 
# 10                     248                     593 
# 
# TRUE 
# 81334 

# [1] "illuminaexome"
# [1] "EUR samples: 2997"
# 
# TRUE 
# 2997 
# 
# TRUE 
# 1 
# [1] "N samples with no pheno: 0"
# [1] "N Finish Samples: 440"
# [1] "N Duplicated Samples: 1215"
# [1] "N JH samples: 0"
# < table of extent 0 >
#   [1] "N relatives to exclude illuminaexome: 30"
# [1] "N sex=0 to exclude illuminaexome: 0"
# [1] "N sex=0 per sex "
# < table of extent 0 >
#   [1] "N samples to include in analysis illuminaexome: 1333"
# [1] "IBD"
# 
# 0    1 
# 321 1012 
# [1] "CD"
# 
# 0   1 
# 321 549 
# [1] "UC"
# 
# 0   1 
# 321 424 
# [1] "SEX"
# 
# 1   2 
# 624 709 
# [1] "SEX IBD"
# 
# 0    1 <NA>
#   1     136  488    0
# 2     185  524    0
# <NA>    0    0 1664
# 
# helmsley_prism_gsa helmsley_xavier_prism_gsa            prism_nfe_gwas 
# 78                       854                       401 
# 
# TRUE 
# 2997 


###################

write.table(final_phen,paste(path,"pheno/phenotype_data_for_all_cohorts_Dec2022_analysis.tsv",sep=""),col.names=T,row.names=F,sep="\t")


dim(final_phen)
# [1] 135633      8



# CREATE A PHENO FILE PER ANCESTRY, AND ADD PCS TO THE COVARIATE FILE


ancestry<-c("eur_all","eur_nonjewish","eur_jewish")


# create files for eur_all as well as for eur_jewish, and eur_nonjewish - and add the population PCs 

# see step 13 for how this file was generated
pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/pca_1000gp_iibdgc_pca_pcaEur_withDuplicates.tsv",sep=""),head=T,sep="\t")
dim(pca)

for (jj in c(1:length(ancestry))) {
  
  print(ancestry[jj])
  
  if(ancestry[jj]=="eur_jewish") {
    array<-c("illumina370","gsa")
  }else{
    array<-c("illumina370","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")
  }
  
  for (ii in 1:length(array)) {
    
    print(array[ii])
  
    eur<-final_phen[which(final_phen$array==array[ii]),]
    
    if(ancestry[jj]=="eur_all") {
      group<-c("Jewish","Non-Jewish")
    } else if(ancestry[jj]=="eur_jewish") {
      group<-c("Jewish")
    } else if(ancestry[jj]=="eur_nonjewish") {
      group<-c("Non-Jewish")
    }
    
    # set to missing pheno anything not in ancestry group
    eur$ibd[which(!eur$FID %in% pca$FID[which(pca$pca_jewish %in% group)])]<-NA
    eur$cd[which(!eur$FID %in% pca$FID[which(pca$pca_jewish %in% group)])]<-NA
    eur$uc[which(!eur$FID %in% pca$FID[which(pca$pca_jewish %in% group)])]<-NA
    
    print("IBD")
    print(table(eur$ibd))

    print("CD")
    print(table(eur$cd))

    print("UC")
    print(table(eur$uc))

    # save phenotype data:
    write.table(eur[,c("FID","IID","ibd")],paste(path,"post_imputation/2022/",array[ii],"/phenotype_data/",array[ii],"_all_studies_merged_",ancestry[jj],"_phenotype_ibd",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    write.table(eur[,c("FID","IID","cd")],paste(path,"post_imputation/2022/",array[ii],"/phenotype_data/",array[ii],"_all_studies_merged_",ancestry[jj],"_phenotype_cd",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    write.table(eur[,c("FID","IID","uc")],paste(path,"post_imputation/2022/",array[ii],"/phenotype_data/",array[ii],"_all_studies_merged_",ancestry[jj],"_phenotype_uc",sep=""),col.names=T,row.names=F,quote=F,sep="\t")


    # create covariate file with sex to adjust for, no NA allowed:
    eur<-eur[which(!is.na(eur$ibd)),c("FID","IID","sex")]
    
    # add pcs:
    # see step 29 for how this file was generated
    pca_2<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/",array[ii],"_all_studies_merged_",ancestry[jj],"_pruned_withDuplicates_ref_pcs_new_projection.sscore",sep=""),head=F)
    colnames(pca_2)<-c("FID","IID","PHENO1","ALLELE_CT","NAMED_ALLELE_DOSAGE_SUM",paste("PC",seq(1:40),sep=""))
    
    
    # see step 29 for selection of number of PCs to include in each analysis
    if (ancestry[jj]=="eur_all") {
      if(array[ii]=="illumina370") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1:2),sep=""))],by="FID",all.x=T,sort=F)
      } else if(array[ii]=="affymetrix6") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1),sep=""))],by="FID",all.x=T,sort=F)
      } else if(array[ii]=="humanomniexpress") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1),sep=""))],by="FID",all.x=T,sort=F)
      } else if(array[ii]=="affymetrix500") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1:3),sep=""))],by="FID",all.x=T,sort=F)
      } else if(array[ii]=="humancoreexome") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1:3),sep=""))],by="FID",all.x=T,sort=F)
      } else if(array[ii]=="humanomni1") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1:3),sep=""))],by="FID",all.x=T,sort=F)
      } else if(array[ii]=="quad610") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1),sep=""))],by="FID",all.x=T,sort=F)
      } else if(array[ii]=="gsa") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1:6),sep=""))],by="FID",all.x=T,sort=F)
      } else if(array[ii]=="illuminaexome") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1),sep=""))],by="FID",all.x=T,sort=F)
      }
    }


    if (ancestry[jj]=="eur_nonjewish") {
      if(array[ii]=="illumina370") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1:5),sep=""))],by="FID",all.x=T,sort=F)
      } else if(array[ii]=="affymetrix6") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1),sep=""))],by="FID",all.x=T,sort=F)
      } else if(array[ii]=="humanomniexpress") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1),sep=""))],by="FID",all.x=T,sort=F)
      } else if(array[ii]=="affymetrix500") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1:3),sep=""))],by="FID",all.x=T,sort=F)
      } else if(array[ii]=="humancoreexome") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1:2),sep=""))],by="FID",all.x=T,sort=F)
      } else if(array[ii]=="humanomni1") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1:2),sep=""))],by="FID",all.x=T,sort=F)
      } else if(array[ii]=="quad610") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1),sep=""))],by="FID",all.x=T,sort=F)
      } else if(array[ii]=="gsa") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1:7),sep=""))],by="FID",all.x=T,sort=F)
      } else if(array[ii]=="illuminaexome") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1:4),sep=""))],by="FID",all.x=T,sort=F)
      }
    }

    if (ancestry[jj]=="eur_jewish") {
      if(array[ii]=="illumina370") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1:2),sep=""))],by="FID",all.x=T,sort=F)
      } else if(array[ii]=="gsa") {
        eur<-merge(eur,pca_2[,c("FID",paste("PC",seq(1:2),sep=""))],by="FID",all.x=T,sort=F)
      }
    }

    print(head(eur))
    # print(nrow(eur[which(is.na(eur$PC1)),]))
    
    write.table(eur,paste(path,"post_imputation/2022/",array[ii],"/phenotype_data/",array[ii],"_all_studies_merged_",ancestry[jj],"_covariate_sex_PCs",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    
  }
}

############
# [1] "eur_all"
# [1] "illumina370"
# [1] "IBD"
# 
# 0    1 
# 1852 1411 
# [1] "CD"
# 
# 0    1 
# 1536  425 
# [1] "UC"
# 
# 0    1 
# 1040  986 
# FID IID sex       PC1       PC2
# 1 111 111   1 0.0794304 0.0196129
# 2 116 116   1 0.0788237 0.0161966
# 3 119 119   2 0.0750553 0.0133555
# 4 131 131   1 0.0800506 0.0133744
# 5 156 156   2 0.0809209 0.0157113
# 6 169 169   1 0.0745985 0.0115431
# [1] "affymetrix6"
# [1] "IBD"
# 
# 0    1 
# 4643 3537 
# [1] "CD"
# 
# 0    1 
# 2762   56 
# [1] "UC"
# 
# 0    1 
# 4643 3478 
# FID         IID sex         PC1
# 1 sample_id sample_id   1 -0.01358670
# 2 sample_id sample_id   1 -0.00890804
# 3 sample_id sample_id   1 -0.00993558
# 4 sample_id sample_id   1  0.02113140
# 5 sample_id sample_id   2 -0.00795526
# 6 sample_id sample_id   2 -0.00490796
# [1] "humanomniexpress"
# [1] "IBD"
# 
# 0   1 
# 589 646 
# [1] "CD"
# 
# 0 
# 589 
# [1] "UC"
# 
# 0   1 
# 589 646 
# FID        IID sex         PC1
# 1 sample_id sample_id   2  0.00690944
# 2 sample_id sample_id   2 -0.00495356
# 3 sample_id sample_id   1  0.02929900
# 4 sample_id sample_id   1  0.09854200
# 5 sample_id sample_id   1  0.00931770
# 6 sample_id sample_id   2 -0.00861860
# [1] "affymetrix500"
# [1] "IBD"
# 
# 0    1 
# 2886 1705 
# [1] "CD"
# 
# 0    1 
# 2886 1656 
# [1] "UC"
# 
# 0    1 
# 2886   48 
# FID        IID sex          PC1          PC2         PC3
# 1 sample_id sample_id   2  1.01712e-03 -0.006480040 -0.00112107
# 2 sample_id sample_id   2 -1.19223e-02  0.019985000 -0.01068530
# 3 sample_id sample_id   2 -1.11992e-03 -0.014368300  0.00322521
# 4 sample_id sample_id   2 -5.19528e-04  0.005245670 -0.00559486
# 5 sample_id sample_id   2 -7.46729e-05 -0.004905040  0.00434556
# 6 sample_id sample_id   2  2.55296e-03 -0.000270089  0.00256630
# [1] "humancoreexome"
# [1] "IBD"
# 
# 0     1 
# 10313 10250 
# [1] "CD"
# 
# 0     1 
# 10313  5489 
# [1] "UC"
# 
# 0     1 
# 10313  4553 
# FID                      IID sex          PC1
# 1 295064_A01_usgwas5491570 295064_A01_usgwas5491570   1  0.002109500
# 2 295064_A02_usgwas5491581 295064_A02_usgwas5491581   1  0.002047930
# 3 295064_A03_usgwas5491593 295064_A03_usgwas5491593   1  0.010051600
# 4 295064_A04_usgwas5491602 295064_A04_usgwas5491602   1  0.000806616
# 5 295064_A05_usgwas5491612 295064_A05_usgwas5491612   2  0.000532126
# 6 295064_A06_usgwas5491620 295064_A06_usgwas5491620   2 -0.000666420
# PC2         PC3
# 1  0.01010100  0.00200736
# 2 -0.00258677  0.00254523
# 3 -0.01688130  0.00125684
# 4  0.01019960  0.00774624
# 5 -0.00208568 -0.00322202
# 6  0.00686395  0.00257825
# [1] "humanomni1"
# [1] "IBD"
# 
# 0    1 
# 1273  761 
# [1] "CD"
# 
# 0    1 
# 1273  386 
# [1] "UC"
# 
# 0    1 
# 1273  369 
# FID        IID sex         PC1          PC2          PC3
# 1 sample_id sample_id   2  0.03357550  0.014560500  0.015089700
# 2 sample_id sample_id   2 -0.01506530 -0.000221472 -0.003470120
# 3 sample_id sample_id   2 -0.01233630 -0.004350460 -0.000353052
# 4 sample_id sample_id   1 -0.00674724 -0.009100180  0.005353800
# 5 sample_id sample_id   1 -0.00832581 -0.011532800 -0.003757100
# 6 sample_id sample_id   2 -0.01685040 -0.012507700 -0.003507890
# [1] "quad610"
# [1] "IBD"
# 
# 0    1 
# 1453 1911 
# [1] "CD"
# 
# 0    1 
# 1453 1130 
# [1] "UC"
# 
# 0    1 
# 1453  781 
# FID        IID sex         PC1
# 1 sample_id sample_id   1  0.00419006
# 2 sample_id sample_id   2  0.00700849
# 3 sample_id sample_id   2  0.00498251
# 4 sample_id sample_id   2 -0.00780600
# 5 sample_id sample_id   2 -0.01056520
# 6 sample_id sample_id   1 -0.00402292
# [1] "gsa"
# [1] "IBD"
# 
# 0     1 
# 16730 41268 
# [1] "CD"
# 
# 0     1 
# 16730 25869 
# [1] "UC"
# 
# 0     1 
# 16730 13847 
# FID        IID sex         PC1         PC2         PC3         PC4
# 1 sample_id sample_id   1 -0.01525920  0.00788711 -0.00956283  0.00456275
# 2 sample_id sample_id   1 -0.02064310 -0.00129437 -0.00372167  0.00517153
# 3 sample_id sample_id   1 -0.01985500  0.00229600 -0.00592546  0.00732909
# 4 sample_id sample_id   2 -0.00328851  0.00947697  0.00164722  0.00283761
# 5 sample_id sample_id   2  0.05545220 -0.00705540 -0.01833060  0.00455338
# 6 sample_id sample_id   1 -0.00643564  0.01162640 -0.00193812 -0.00551566
# PC5          PC6
# 1 -0.008559680 -0.008299030
# 2 -0.004630660  0.000883184
# 3 -0.001461150 -0.001080030
# 4 -0.005378930 -0.017579400
# 5  0.000720782 -0.001315150
# 6 -0.001047890 -0.006126110
# [1] "illuminaexome"
# [1] "IBD"
# 
# 0    1 
# 321 1012 
# [1] "CD"
# 
# 0   1 
# 321 549 
# [1] "UC"
# 
# 0   1 
# 321 424 
# FID        IID sex        PC1
# 1 sample_id sample_id   2  0.0270847
# 2 sample_id sample_id   2  0.0339443
# 3 sample_id sample_id   1 -0.0761287
# 4 sample_id sample_id   2 -0.0101471
# 5 sample_id sample_id   1  0.0122968
# 6 sample_id sample_id   2 -0.0767798
# [1] "eur_nonjewish"
# [1] "illumina370"
# [1] "IBD"
# 
# 0    1 
# 1473 1276 
# [1] "CD"
# 
# 0    1 
# 1157  299 
# [1] "UC"
# 
# 0   1 
# 669 977 
# FID  IID sex       PC1       PC2          PC3       PC4          PC5
# 1  156  156   2 0.0623434 0.0665440  0.029520100 0.0328773  0.006979220
# 2  312  312   2 0.0236954 0.0121878 -0.000802426 0.0301351  0.008869870
# 3  416  416   1 0.0592152 0.0679680  0.027684600 0.0392271 -0.000904147
# 4  926  926   2 0.0249931 0.0240979  0.010105400 0.0277673 -0.011680700
# 5  986  986   2 0.0556714 0.0595461  0.029491900 0.0244286  0.005037230
# 6 1137 1137   1 0.0263075 0.0194252  0.022778900 0.0206466  0.008663090
# [1] "affymetrix6"
# [1] "IBD"
# 
# 0    1 
# 4642 3537 
# [1] "CD"
# 
# 0    1 
# 2762   56 
# [1] "UC"
# 
# 0    1 
# 4642 3478 
# FID         IID sex         PC1
# 1 sample_id sample_id   1 -0.01357660
# 2 sample_id sample_id   1 -0.00892828
# 3 sample_id sample_id   1 -0.00995059
# 4 sample_id sample_id   1  0.02109300
# 5 sample_id sample_id   2 -0.00794775
# 6 sample_id sample_id   2 -0.00491365
# [1] "humanomniexpress"
# [1] "IBD"
# 
# 0   1 
# 586 638 
# [1] "CD"
# 
# 0 
# 586 
# [1] "UC"
# 
# 0   1 
# 586 638 
# FID        IID sex         PC1
# 1 sample_id sample_id   2 -0.00880226
# 2 sample_id sample_id   2  0.00444183
# 3 sample_id sample_id   1 -0.03241790
# 4 sample_id sample_id   1 -0.10680500
# 5 sample_id sample_id   1 -0.01108450
# 6 sample_id sample_id   2  0.00905540
# [1] "affymetrix500"
# [1] "IBD"
# 
# 0    1 
# 2873 1666 
# [1] "CD"
# 
# 0    1 
# 2873 1617 
# [1] "UC"
# 
# 0    1 
# 2873   48 
# FID        IID sex          PC1         PC2          PC3
# 1 sample_id sample_id   2  0.003970850  0.00514613  0.000158041
# 2 sample_id sample_id   2 -0.020151400 -0.01051480  0.011575700
# 3 sample_id sample_id   2  0.004217240  0.01385730 -0.005378000
# 4 sample_id sample_id   2 -0.003546270 -0.00193568  0.005293670
# 5 sample_id sample_id   2  0.003318710  0.00255693 -0.004432890
# 6 sample_id sample_id   2  0.000681187  0.00185000 -0.003507210
# [1] "humancoreexome"
# [1] "IBD"
# 
# 0     1 
# 10272 10143 
# [1] "CD"
# 
# 0     1 
# 10272  5421 
# [1] "UC"
# 
# 0     1 
# 10272  4516 
# FID                      IID sex          PC1
# 1 295064_A01_usgwas5491570 295064_A01_usgwas5491570   1  0.001066890
# 2 295064_A02_usgwas5491581 295064_A02_usgwas5491581   1  0.001001610
# 3 295064_A03_usgwas5491593 295064_A03_usgwas5491593   1  0.009693110
# 4 295064_A04_usgwas5491602 295064_A04_usgwas5491602   1 -0.002135250
# 5 295064_A05_usgwas5491612 295064_A05_usgwas5491612   2  0.000768437
# 6 295064_A06_usgwas5491620 295064_A06_usgwas5491620   2 -0.001336030
# PC2
# 1 -0.01027780
# 2  0.00269282
# 3  0.01670220
# 4 -0.00981058
# 5  0.00189857
# 6 -0.00685974
# [1] "humanomni1"
# [1] "IBD"
# 
# 0    1 
# 1200  748 
# [1] "CD"
# 
# 0    1 
# 1200  377 
# [1] "UC"
# 
# 0    1 
# 1200  365 
# FID        IID sex         PC1         PC2
# 1 sample_id sample_id   2  0.03637040  0.01873660
# 2 sample_id sample_id   2 -0.01529650 -0.00169088
# 3 sample_id sample_id   2 -0.01069930 -0.00531111
# 4 sample_id sample_id   1 -0.00396598 -0.00917549
# 5 sample_id sample_id   1 -0.00773294 -0.01234630
# 6 sample_id sample_id   2 -0.01598140 -0.01410270
# [1] "quad610"
# [1] "IBD"
# 
# 0    1 
# 1433 1890 
# [1] "CD"
# 
# 0    1 
# 1433 1119 
# [1] "UC"
# 
# 0    1 
# 1433  771 
# FID        IID sex         PC1
# 1 sample_id sample_id   1 -0.00407599
# 2 sample_id sample_id   2 -0.00706029
# 3 sample_id sample_id   2 -0.00468999
# 4 sample_id sample_id   2  0.00766555
# 5 sample_id sample_id   2  0.01061370
# 6 sample_id sample_id   1  0.00384590
# [1] "gsa"
# [1] "IBD"
# 
# 0     1 
# 16038 37085 
# [1] "CD"
# 
# 0     1 
# 16038 23221 
# [1] "UC"
# 
# 0     1 
# 16038 12423 
# FID        IID sex          PC1        PC2         PC3         PC4
# 1 sample_id sample_id   1  0.011295400 0.01544890 -0.00361827 -0.00724263
# 2 sample_id sample_id   1  0.018318600 0.00656943 -0.00419269 -0.00633817
# 3 sample_id sample_id   1  0.016695800 0.01044090 -0.00636128 -0.00322785
# 4 sample_id sample_id   2 -0.002665610 0.00822851 -0.00246717 -0.00739365
# 5 sample_id sample_id   1 -0.000714269 0.01208700  0.00549392 -0.00321770
# 6 sample_id sample_id   2  0.017491800 0.01259480 -0.00426107 -0.00586299
# PC5          PC6          PC7
# 1  0.005257700 -0.008150460  0.000758192
# 2 -0.000487347  0.000572827 -0.001789670
# 3 -0.001585880 -0.001405360  0.000565754
# 4 -0.000213147 -0.017686200  0.004633680
# 5 -0.002631470 -0.006965670 -0.004658850
# 6 -0.001554880 -0.003095670  0.003331350
# [1] "illuminaexome"
# [1] "IBD"
# 
# 0   1 
# 292 877 
# [1] "CD"
# 
# 0   1 
# 292 478 
# [1] "UC"
# 
# 0   1 
# 292 370 
# FID        IID sex         PC1          PC2         PC3         PC4
# 1 sample_id sample_id   2  0.02759600  0.004130410 -0.00282040  0.00586041
# 2 sample_id sample_id   2  0.03450110  0.005509510  0.00524613  0.00582955
# 3 sample_id sample_id   2 -0.01401730  0.000238889 -0.00188893 -0.00166550
# 4 sample_id sample_id   1  0.00891672 -0.001618070 -0.02050200  0.00205244
# 5 sample_id sample_id   1  0.01013850 -0.004158590 -0.00723442  0.01684100
# 6 sample_id sample_id   1  0.02710170  0.003444860  0.00855846 -0.00510755
# [1] "eur_jewish"
# [1] "illumina370"
# [1] "IBD"
# 
# 0   1 
# 379 135 
# [1] "CD"
# 
# 0   1 
# 379 126 
# [1] "UC"
# 
# 0   1 
# 371   9 
# FID IID sex       PC1         PC2
# 1 111 111   1 0.0369600 -0.03176020
# 2 116 116   1 0.0316661 -0.02865750
# 3 119 119   2 0.0213008  0.00662362
# 4 131 131   1 0.0329191 -0.01613690
# 5 169 169   1 0.0249993 -0.04573150
# 6 181 181   2 0.0271946 -0.02646450
# [1] "gsa"
# [1] "IBD"
# 
# 0    1 
# 692 4183 
# [1] "CD"
# 
# 0    1 
# 692 2648 
# [1] "UC"
# 
# 0    1 
# 692 1424 
# FID        IID sex       PC1        PC2
# 1 sample_id sample_id   2 0.0272260 -0.0194518
# 2 sample_id sample_id   2 0.0338263 -0.0213178
# 3 sample_id sample_id   1 0.0338376 -0.0212192
# 4 sample_id sample_id   1 0.0343479 -0.0223625
# 5 sample_id sample_id   1 0.0335777 -0.0214412
# 6 sample_id sample_id   1 0.0283227 -0.0191977
  

############

# summary for phil - to identify samples from John Hophkings

dim(nodup)
# [1] 114977      7

final<-nodup[,c("V1","study")]
fwrite(final,paste(path,"pheno/list_noduplicated_samples_2022Jan31.tsv.gz",sep=""),col.names=T,row.names=F,sep="\t")



# subset UKIBDGC data to share with rest of the team:

final_phen<-fread(paste(path,"pheno/phenotype_data_for_all_cohorts_Dec2022_analysis.tsv",sep=""),head=T)
final_phen<-final_phen[which(final_phen$study %in% c("all_hce","gwas1","gwas2")),]


# see step 13 for how this file was generated
pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/pca_1000gp_iibdgc_pca_pcaEur_withDuplicates.tsv",sep=""),head=T,sep="\t")
dim(pca)

final_phen<-merge(final_phen,pca[,c("FID","inferred_population","pca_jewish")],by="FID",all.x=T)
    
array<-names(table(final_phen$array))

for (i in 1:length(array)) {
  
  pca_tmp<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/",array[i],"_all_studies_merged_eur_all_pruned_withDuplicates_ref_pcs_new_projection.sscore",sep=""),head=F)
  colnames(pca_tmp)<-c("FID","IID","PHENO1","ALLELE_CT","NAMED_ALLELE_DOSAGE_SUM",paste("PC",seq(1:40),sep=""))
  
  if (i==1) {
    pca_2<-pca_tmp
  }else{
    pca_2<-rbind(pca_2,pca_tmp)
  }
}


final_phen<-merge(final_phen,pca_2[,c("FID","PC1","PC2","PC3")],by="FID",all.x=T)

final_phen$PC2[which(final_phen$array=="affymetrix6")]<-NA
final_phen$PC3[which(final_phen$array=="affymetrix6")]<-NA

fwrite(final_phen,"/path/to/project",col.names=T,row.names=F,
       quote=F,sep="\t",na="NA")



