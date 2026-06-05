# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# Evaluate in which studies / arrays it would be feasible to analyze chrX pseudoautosomal region


## /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"
studies<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
           "kiel_austria_sibdcs_gsa","netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa",
           "niddk_feinstein_gsa","basque_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa",
           "belgium_vermeire_gsa","prism_nfe_gsa","prism_nfe_gwas","finland_illugwas","german_affy6_old_gwas",
           "norway_affy6_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
           "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas",
           # "swedish_uc_old_gwas",
           "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

dat<-as.data.frame(matrix(nrow=length(studies),ncol=5))
colnames(dat)<-c("studies","N_variants_chrX_PAR1","N_variants_chrX_PAR2","N_variants_chrX_PAR1_info_0.4","N_variants_chrX_PAR2_info_0.4")

for (j in 1:length(studies)) {
  
  print(studies[j])
  
  x<-fread(paste("/path/to/ibdgwas/IIBDGC/imputed/",studies[j],"/2022/eur/chrX.info.gz",sep=""),head=T)
  x$position<-gsub("chrX:","",x$SNP)
  x$position<-gsub(":.*","",x$position)
  x$position<-as.numeric(x$position)
  
  dat$studies[j]<-studies[j]
  dat$N_variants_chrX_PAR1[j]<-nrow(x[which(x$position< 2781479),])
  dat$N_variants_chrX_PAR2[j]<-nrow(x[which(x$position>155701383),])
  
  dat$N_variants_chrX_PAR1_info_0.4[j]<-nrow(x[which(x$position< 2781479 & x$Rsq>=0.4),])
  dat$N_variants_chrX_PAR2_info_0.4[j]<-nrow(x[which(x$position>155701383 & x$Rsq>=0.4),])
  
  rm(x)
  
}

# studies N_variants_chrX_PAR1 N_variants_chrX_PAR2
# 1                    all_hce                    0                    0
# 2             niddk_old_gwas                    0                    0
# 3        australia_omniexome                33495                 2296
# 4                      gwas1                    0                    0
# 5                      gwas2                64956                  436
# 6             pittsburgh_gsa                    0                    0
# 7                  spain_gsa                  251                    0
# 8                  italy_gsa                34592                 2740
# 9    kiel_austria_sibdcs_gsa                49771                 4398
# 10           netherlands_gsa                32690                 2922
# 11              slovenia_gsa                10866                 1533
# 12                sweden_gsa                24624                 2157
# 13           niddk_broad_gsa                54877                    0 # WHY 0 IN PAR2
# 14       niddk_feinstein_gsa                69980                 4102
# 15                basque_gsa                46737                 2781
# 16             lithuania_gsa                22762                 2172
# 17         belgium_louis_gsa                19393                 2602
# 18   belgium_franchimont_gsa                21196                 2691
# 19      belgium_vermeire_gsa                29669                 3067
# 20             prism_nfe_gsa                18012                 1969
# 21            prism_nfe_gwas                    0                    0
# 22          finland_illugwas                    0                    0
# 23     german_affy6_old_gwas                32085                  303
# 24     norway_affy6_old_gwas                18676                  447
# 25     belgium_inf1_old_gwas                    0                    0
# 26     belgium_inf2_old_gwas                    0                    0
# 27      cedars_370k_old_gwas                  161                    0
# 28      cedars_610k_old_gwas                  245                    0
# 29      cedars_omni_old_gwas                45573                 2247
# 30       swedish_uc_old_gwas                   NA                   NA
# 31              mccauley_gsa                34043                 2268
# 32                  ccfa_gsa                36416                 3018
# 33                cedars_gsa                50638                 3414
# 34             bernstein_gsa                12618                 1584
# 35              farkkila_gsa                10847                 1085
# 36           franchimont_gsa                24360                 3089
# 37                franke_gsa                17696                 2016
# 38        helmsley_prism_gsa                    0                    0
# 39 helmsley_xavier_prism_gsa                    0                    0
# 40         hyams_protect_gsa                21066                 1810
# 41           lewis_sparc_gsa                39998                 3111
# 42          mccauley_new_gsa                41379                 3028
# 43              mcgovern_gsa                57332                 3856
# 44      moayyedi_imagine_gsa                25450                 2254
# 45        newberry_share_gsa                26761                 1927
# 46             niddk_cho_gsa                42317                 2772
# 47           niddk_duerr_gsa                31444                 2877
# 48           niddk_rioux_gsa                24312                 2486
# 49      niddk_silverberg_gsa                38981                 3144
# 50           palotie_hus_gsa                18950                 1522
# 51           pekow_share_gsa                25723                 2074
# 52        rioux_igenomed_gsa                16569                 1618
# 53           sands_msccr_gsa                40392                 2586
# 54              stampfer_gsa                33946                 2377
# 55              vermeire_gsa                30832                 3209
# 56               weersma_gsa                17252                 1820
# 57          xavier_prism_gsa                18188                 2164
# 58          xavier_share_gsa                27915                 2166

  
  
  
## PRE IMPUTATION

dat<-as.data.frame(matrix(nrow=length(studies),ncol=3))
colnames(dat)<-c("studies","N_variants_chrX_PAR1","N_variants_chrX_PAR2")

for (j in 1:length(studies)) {
  
  print(studies[j])
  x<-fread(paste("/path/to/ibdgwas/IIBDGC/pre_imputation/QC/",studies[j],"/",studies[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim",sep=""),head=F)
  colnames(x)<-c("chr","SNP","cm","position","Ref","Alt")
  
  x<-x[which(x$chr==23),]
  
  file<-paste("/path/to/ibdgwas/IIBDGC/imputed/",studies[j],"/qc/2022/list_varias_to_exclude_3",sep="")
 
  x$position<-as.numeric(x$position)
  
  if(file.exists(file)) {
    exclude<-read.table(file,head=F)
  }else{
    exclude<-read.table(paste("/path/to/ibdgwas/IIBDGC/imputed/",studies[j],"/qc/2022/list_variants_to_exclude_2",sep=""),head=F)
  }
  
  exclude$chr<-gsub(":.*","",exclude$V1)
  exclude<-exclude[which(exclude$chr==23),]
  exclude$position<-gsub("23:","",exclude$V1)
  exclude$position<-as.numeric(gsub("_.*","",exclude$position))
  
  x<-x[which(!x$SNP %in% exclude$V1),]
  
  dat$studies[j]<-studies[j]
  dat$N_variants_chrX_PAR1[j]<-nrow(x[which(x$position< 2781479),])
  dat$N_variants_chrX_PAR2[j]<-nrow(x[which(x$position>155701383),])
  
}


# studies N_variants_chrX_PAR1 N_variants_chrX_PAR2
# 1                    all_hce                   22                    3
# 2             niddk_old_gwas                    2                    0
# 3        australia_omniexome                  233                   42
# 4                      gwas1                    0                    0
# 5                      gwas2                  276                    4
# 6             pittsburgh_gsa                    0                    0
# 7                  spain_gsa                    5                    0
# 8                  italy_gsa                  378                   68
# 9    kiel_austria_sibdcs_gsa                  246                   57
# 10           netherlands_gsa                  245                   58
# 11              slovenia_gsa                  224                   55
# 12                sweden_gsa                  235                   56
# 13           niddk_broad_gsa                  268                   67
# 14       niddk_feinstein_gsa                  288                   58
# 15                basque_gsa                  411                   63
# 16             lithuania_gsa                  232                   49
# 17         belgium_louis_gsa                  230                   56
# 18   belgium_franchimont_gsa                  230                   58
# 19      belgium_vermeire_gsa                  235                   58
# 20             prism_nfe_gsa                  221                   54
# 21            prism_nfe_gwas                    0                    0
# 22          finland_illugwas                    0                    0
# 23     german_affy6_old_gwas                  218                    5
# 24     norway_affy6_old_gwas                  217                    4
# 25     belgium_inf1_old_gwas                    0                    0
# 26     belgium_inf2_old_gwas                    0                    0
# 27      cedars_370k_old_gwas                    3                    0
# 28      cedars_610k_old_gwas                    6                    0
# 29      cedars_omni_old_gwas                  274                   49
# 30              mccauley_gsa                  281                   56
# 31                  ccfa_gsa                  261                   59
# 32                cedars_gsa                  281                   60
# 33             bernstein_gsa                  147                   31
# 34              farkkila_gsa                  178                   44
# 35           franchimont_gsa                  241                   60
# 36                franke_gsa                  222                   53
# 37        helmsley_prism_gsa                    0                    0
# 38 helmsley_xavier_prism_gsa                    0                    0
# 39         hyams_protect_gsa                  252                   54
# 40           lewis_sparc_gsa                  262                   60
# 41          mccauley_new_gsa                  286                   61
# 42              mcgovern_gsa                  278                   61
# 43      moayyedi_imagine_gsa                  240                   57
# 44        newberry_share_gsa                  256                   54
# 45             niddk_cho_gsa                  281                   59
# 46           niddk_duerr_gsa                  251                   60
# 47           niddk_rioux_gsa                  252                   57
# 48      niddk_silverberg_gsa                  261                   62
# 49           palotie_hus_gsa                  198                   45
# 50           pekow_share_gsa                  244                   53
# 51        rioux_igenomed_gsa                  242                   55
# 52           sands_msccr_gsa                  294                   56
# 53              stampfer_gsa                  267                   54
# 54              vermeire_gsa                  234                   56
# 55               weersma_gsa                  230                   55
# 56          xavier_prism_gsa                  232                   54
# 57          xavier_share_gsa                  256                   58


## see amount of sex=0:

${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${j}_nomonom

###################################################################

j=13
fam<-fread(paste("/path/to/ibdgwas/IIBDGC/pre_imputation/QC/",studies[j],"/",studies[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam",sep=""),head=F)
bim<-fread(paste("/path/to/ibdgwas/IIBDGC/pre_imputation/QC/",studies[j],"/",studies[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim",sep=""),head=F)

# look at frequency:

#############

i=niddk_broad_gsa
j=eur
path_gwas=/path/to/ibdgwas/IIBDGC/


emacs ${path_gwas}pre_imputation/QC/${i}/test_range_par2 
# chr23      155701383       158000000       1



/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${j}_nomonom \
--allow-no-sex \
--keep-allele-order \
--flip ${path_gwas}imputed/${i}/qc/list_variants_negative_EmpR_toflip \
--exclude ${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude_2 \
--freq counts \
--extract range ${path_gwas}pre_imputation/QC/${i}/test_range_par2 \
--missing \
--out ${path_gwas}pre_imputation/QC/${i}/test

#############

freq<-fread(paste("/path/to/ibdgwas/IIBDGC/pre_imputation/QC/",studies[j],"/test.frq.counts",sep=""),head=T)
freq$position<-gsub("23:","",freq$SNP)
freq$position<-as.numeric(gsub("_.*","",freq$position))

freq[which(freq$position>155701383),]

varmiss<-fread(paste("/path/to/ibdgwas/IIBDGC/pre_imputation/QC/",studies[j],"/test.lmiss",sep=""),head=T)
varmiss$position<-gsub("23:","",varmiss$SNP)
varmiss$position<-as.numeric(gsub("_.*","",varmiss$position))

varmiss[which(varmiss$position>155701383),]

summary(varmiss$F_MISS)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.001365 0.002534 0.006446 0.006823 0.051070 

samplemiss<-fread(paste("/path/to/ibdgwas/IIBDGC/pre_imputation/QC/",studies[j],"/test.imiss",sep=""),head=T)
summary(samplemiss$F_MISS)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000000 0.000000 0.003399 0.000000 0.552200 

samplemiss[which(samplemiss$F_MISS==max(samplemiss$F_MISS)),]
#           FID        IID MISS_PHENO N_MISS N_GENO F_MISS position
# 1: sample_id sample_id          N     37     67 0.5522       NA


x<-fread(paste("/path/to/ibdgwas/IIBDGC/imputed/",studies[j],"/2022/noneur/chrX.info.gz",sep=""),head=T)
x$position<-gsub("chrX:","",x$SNP)
x$position<-gsub(":.*","",x$position)
x$position<-as.numeric(x$position)

nrow(x[which(x$position>155701383),])
# [1] 4029 # non eur data has par2 data, issue with run - or with samples! try again with exact same file

# resubmitted, error message 1 Chunk(s) excluded: at least one sample has a call rate < 50.0% - so exclude samples with larger missingess in PAR2
samplemiss$F_MISS<-as.numeric(samplemiss$F_MISS)
samplemiss<-as.data.frame(samplemiss)
samplemiss[which(samplemiss$F_MISS==max(samplemiss$F_MISS)),]
# FID        IID MISS_PHENO N_MISS N_GENO F_MISS position
# 3642 sample_id sample_id          N   3353  12218 0.2744       NA

samplemiss[which(samplemiss$F_MISS>0.20),]









