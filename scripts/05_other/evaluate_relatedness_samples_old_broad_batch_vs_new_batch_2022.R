# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
### compare which sample IDs belong to the same individual - old batches vs new batches

#### PART 1 - TO EVALUATE ROUGHLY WHICH SAMPLES WERE INCLUDED - SEE BELOW PART 2 FOR FINAL SELECTION

###########################
# 1.- SELECT THE VARIANTS #
###########################

# SEE SCRIPT:
# create_list_common_variants_among_cohorts - USE SAME VARIANTS AS IN STEP 12.3 create_list_common_variants_among_cohorts 
path=/path/to/ibdgwas/IIBDGC/
wc -l ${path}pre_imputation/QC/relatedness/list_common_variants_old_new_broad_batch
# 118825 

########################
# 2.- EXTRACT VARIANTS #
########################

studies=(basque_gsa belgium_franchimont_gsa belgium_louis_gsa belgium_vermeire_gsa ccfa_gsa cedars_gsa italy_gsa kiel_austria_sibdcs_gsa lithuania_gsa mccauley_gsa netherlands_gsa niddk_broad_gsa niddk_feinstein_gsa prism_nfe_gsa slovenia_gsa sweden_gsa prism_nfe_gwas finland_illugwas bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip \
--allow-no-sex \
--extract ${path}pre_imputation/QC/relatedness/list_common_variants_old_new_broad_batch \
--make-bed \
--out ${path}pre_imputation/QC/relatedness/${i}_subset_old_new_broad_batch
done

wc -l ${path}pre_imputation/QC/relatedness/*_subset_old_new_broad_batch.fam
###########
# 1516 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/basque_gsa_subset_old_new_broad_batch.fam
# 1532 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_franchimont_gsa_subset_old_new_broad_batch.fam
# 1531 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_louis_gsa_subset_old_new_broad_batch.fam
# 4014 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_vermeire_gsa_subset_old_new_broad_batch.fam
# 514 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/bernstein_gsa_subset_old_new_broad_batch.fam
# 2188 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/ccfa_gsa_subset_old_new_broad_batch.fam
# 3096 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_gsa_subset_old_new_broad_batch.fam
# 68 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/farkkila_gsa_subset_old_new_broad_batch.fam
# 457 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/finland_illugwas_subset_old_new_broad_batch.fam
# 2789 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/franchimont_gsa_subset_old_new_broad_batch.fam
# 885 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/franke_gsa_subset_old_new_broad_batch.fam
# 780 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/helmsley_prism_gsa_subset_old_new_broad_batch.fam
# 1297 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/helmsley_xavier_prism_gsa_subset_old_new_broad_batch.fam
# 418 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/hyams_protect_gsa_subset_old_new_broad_batch.fam
# 1016 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/italy_gsa_subset_old_new_broad_batch.fam
# 14660 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/kiel_austria_sibdcs_gsa_subset_old_new_broad_batch.fam
# 2859 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/lewis_sparc_gsa_subset_old_new_broad_batch.fam
# 2277 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/lithuania_gsa_subset_old_new_broad_batch.fam
# 788 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mccauley_gsa_subset_old_new_broad_batch.fam
# 1628 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mccauley_new_gsa_subset_old_new_broad_batch.fam
# 6050 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mcgovern_gsa_subset_old_new_broad_batch.fam
# 1128 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/moayyedi_imagine_gsa_subset_old_new_broad_batch.fam
# 4705 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/netherlands_gsa_subset_old_new_broad_batch.fam
# 865 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/newberry_share_gsa_subset_old_new_broad_batch.fam
# 5515 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_broad_gsa_subset_old_new_broad_batch.fam
# 1764 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_cho_gsa_subset_old_new_broad_batch.fam
# 1944 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_duerr_gsa_subset_old_new_broad_batch.fam
# 8301 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_feinstein_gsa_subset_old_new_broad_batch.fam
# 919 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_rioux_gsa_subset_old_new_broad_batch.fam
# 2371 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_silverberg_gsa_subset_old_new_broad_batch.fam
# 878 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/palotie_hus_gsa_subset_old_new_broad_batch.fam
# 634 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/pekow_share_gsa_subset_old_new_broad_batch.fam
# 466 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/prism_nfe_gsa_subset_old_new_broad_batch.fam
# 862 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/prism_nfe_gwas_subset_old_new_broad_batch.fam
# 182 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/rioux_igenomed_gsa_subset_old_new_broad_batch.fam
# 1430 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/sands_msccr_gsa_subset_old_new_broad_batch.fam
# 270 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/slovenia_gsa_subset_old_new_broad_batch.fam
# 1477 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/stampfer_gsa_subset_old_new_broad_batch.fam
# 1414 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/sweden_gsa_subset_old_new_broad_batch.fam
# 4713 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/vermeire_gsa_subset_old_new_broad_batch.fam
# 709 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/weersma_gsa_subset_old_new_broad_batch.fam
# 692 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/xavier_prism_gsa_subset_old_new_broad_batch.fam
# 696 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/xavier_share_gsa_subset_old_new_broad_batch.fam
# 92298 total
###########

wc -l ${path}pre_imputation/QC/relatedness/*_subset_old_new_broad_batch.bim
# 118825  OK

##############################################################################################################################################

#######################
# 3.- MERGE ALL FILES #
#######################

#### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("basque_gsa","belgium_franchimont_gsa","belgium_louis_gsa","belgium_vermeire_gsa","ccfa_gsa","cedars_gsa",
           "italy_gsa","kiel_austria_sibdcs_gsa","lithuania_gsa","mccauley_gsa","netherlands_gsa",
           "niddk_broad_gsa","niddk_feinstein_gsa","prism_nfe_gsa","slovenia_gsa","sweden_gsa",
           "prism_nfe_gwas","finland_illugwas",
           "helmsley_prism_gsa","helmsley_xavier_prism_gsa",
           "bernstein_gsa","farkkila_gsa","franchimont_gsa","franke_gsa",
           "hyams_protect_gsa","lewis_sparc_gsa","mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa",
           "newberry_share_gsa","niddk_cho_gsa","niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa",
           "palotie_hus_gsa","pekow_share_gsa","rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa",
           "vermeire_gsa","weersma_gsa","xavier_prism_gsa","xavier_share_gsa")


length(cohorts)
# [1] 43

#####

cohorts<-cohorts[-1]

dat<-matrix(ncol=3,nrow=length(cohorts))
dat<-as.data.frame(dat)

for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_old_new_broad_batch.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_old_new_broad_batch.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_old_new_broad_batch.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/relatedness/list_cohorts_tomerge_old_new_data_from_broad_gsa.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#############################

# hpc-server

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/relatedness/basque_gsa_subset_old_new_broad_batch \
--allow-no-sex \
--merge-list ${path}pre_imputation/QC/relatedness/list_cohorts_tomerge_old_new_data_from_broad_gsa.txt \
--make-bed --out ${path}pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch

# 3111 variants and 90143 people pass filters and QC.
# Among remaining phenotypes, 65491 are cases and 24652 are controls.


##############################################################################################################################################

############################
# 3.- ESTIMATE RELATEDNESS #
############################

# The current release is KING version 2.2.4 (released on October 11, 2019). 

# --related provides integrative, fast, and accurate inference for close relationships. This option is highly recommended, 
# especially when dealing with biobank-level datasets. Integration of the IBD segment inference furthur improves the inference accuracy. 
# In our tests, --related has been successfully applied to datasets consisting of ~10 million samples. When "--rplot" is specified, 
# several relationship plots are generated automatically. --related --degree 2 specifies that only related pairs
# (up to the 2nd-degree in this case) between families are included in the output. Specifically all pairs across families with 
# a Kinship coefficient less than 0.0884 will be excluded from the output. More details of 
# --related analysis are available in the INTEGRATED RELATIONSHIP INFERENCE section later in this tutorial.

# --duplicate implements the fastest (and accurate) algorithm to identify duplicates or MZ twins. The running time is in seconds, 
# unless the number of samples is > 1,000,000 in which case a few minutes may be needed. In our tests, --duplicate has been successfully 
# applied to datasets consisting of ~10 million samples. One potential application of the duplicate analysis is to identify duplicates 
# accross different studies, in which case multiple datasets can be read in conveniently as shown in GENERAL INPUT FILES section.

# The amount of computer memory required by KING analysis is modest, at ~N ✕ M / 4 (where N is the number of samples and M is the number of SNPs)
# plus a small percentage of overhead cost. E.g., for a dataset consisting of 100,000 samples each genotyped at 1,000,000 SNPs, 
# the required memory size is ~25GB.

# (3110*107598/4)/1E+6
# [1] 83.65744


#################


# The current release is KING version 2.2.4 (released on October 11, 2019). 

# --related provides integrative, fast, and accurate inference for close relationships. This option is highly recommended, 
# especially when dealing with biobank-level datasets. Integration of the IBD segment inference furthur improves the inference accuracy. 
# In our tests, --related has been successfully applied to datasets consisting of ~10 million samples. When "--rplot" is specified, 
# several relationship plots are generated automatically. --related --degree 2 specifies that only related pairs
# (up to the 2nd-degree in this case) between families are included in the output. Specifically all pairs across families with 
# a Kinship coefficient less than 0.0884 will be excluded from the output. More details of 
# --related analysis are available in the INTEGRATED RELATIONSHIP INFERENCE section later in this tutorial.

# --duplicate implements the fastest (and accurate) algorithm to identify duplicates or MZ twins. The running time is in seconds, 
# unless the number of samples is > 1,000,000 in which case a few minutes may be needed. In our tests, --duplicate has been successfully 
# applied to datasets consisting of ~10 million samples. One potential application of the duplicate analysis is to identify duplicates 
# accross different studies, in which case multiple datasets can be read in conveniently as shown in GENERAL INPUT FILES section.

# The amount of computer memory required by KING analysis is modest, at ~N ✕ M / 4 (where N is the number of samples and M is the number of SNPs)
# plus a small percentage of overhead cost. E.g., for a dataset consisting of 100,000 samples each genotyped at 1,000,000 SNPs, 
# the required memory size is ~25GB.


# (92298*118825/4)/1E+6
# [1] 2741.827

path_ukbb=/path/to/project
MEM=5000

bsub -J"kg" -M"$MEM" -n8 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path}pre_imputation/QC/relatedness/log/stderr_king_iibdgc_merged_old_new_broad_batch \
-o ${path}pre_imputation/QC/relatedness/log/stdout_king_iibdgc_merged_old_new_broad_batch \
"/path/to/software/./king \
-b ${path}pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch.bed \
--related --cpus 8 --prefix ${path}pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch"
# Job <23607> is submitted to queue <normal>.

# Relationship summary (total relatives: 0 by pedigree, 32846 by inference)
#                 MZ      PO      FS      2nd
# =====================================================
#   Inference     27204   3303    1778    455
# 
# Between-family relatives (Kinship >= 0.17678) saved in file /path/to/ibdgwas/IIBDGC/pre_
# imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch.kin0
# 
# Note only duplicates and 1st-degree relatives are included in the inference.
# Specifying '--degree 2' if a higher degree relationship inference is needed.


# previous attempt using 3100 markers - shared between all old iibdgc (2020)
# ${path}pre_imputation/QC/relatedness/iibdgc_merged_old_new_data_from_broad_gsa_plink


###################################
# label samples with low call rate - all samples regardless QC are kept


/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch \
--missing --out ${path}pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch
# Sample missing data report written to
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch.imiss,
# and variant-based missing data report written to
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch.lmiss.


################################

#### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("basque_gsa","belgium_franchimont_gsa","belgium_louis_gsa","belgium_vermeire_gsa","ccfa_gsa","cedars_gsa",
           "italy_gsa","kiel_austria_sibdcs_gsa","lithuania_gsa","mccauley_gsa","netherlands_gsa",
           "niddk_broad_gsa","niddk_feinstein_gsa","prism_nfe_gsa","slovenia_gsa","sweden_gsa",
           "prism_nfe_gwas","finland_illugwas",
           "helmsley_prism_gsa","helmsley_xavier_prism_gsa",
           "bernstein_gsa","farkkila_gsa","franchimont_gsa","franke_gsa",
           "hyams_protect_gsa","lewis_sparc_gsa","mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa",
           "newberry_share_gsa","niddk_cho_gsa","niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa",
           "palotie_hus_gsa","pekow_share_gsa","rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa",
           "vermeire_gsa","weersma_gsa","xavier_prism_gsa","xavier_share_gsa")

new<-c("helmsley_prism_gsa","helmsley_xavier_prism_gsa",
       "bernstein_gsa","farkkila_gsa","franchimont_gsa","franke_gsa",
       "hyams_protect_gsa","lewis_sparc_gsa","mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa",
       "newberry_share_gsa","niddk_cho_gsa","niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa",
       "palotie_hus_gsa","pekow_share_gsa","rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa",
       "vermeire_gsa","weersma_gsa","xavier_prism_gsa","xavier_share_gsa")

length(cohorts)
# [1] 43

for (i in 1:length(cohorts)) {
  
  tmp<-read.table(paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_old_new_broad_batch.fam",sep=""),head=F)

  tmp$cohort<-as.character(cohorts[i])
  tmp$V1<-as.character(tmp$V1)
  tmp$V2<-as.character(tmp$V2)
  
  if(i==1){
    dat<-tmp
  }else{
    dat<-rbind(dat,tmp)
  }
}

table(dat$cohort)
# basque_gsa   belgium_franchimont_gsa         belgium_louis_gsa 
#       1516                      1532                      1531 
# belgium_vermeire_gsa             bernstein_gsa                  ccfa_gsa 
#                 4014                       514                      2188 
# cedars_gsa              farkkila_gsa          finland_illugwas 
# 3096                        68                       457 
# franchimont_gsa                franke_gsa        helmsley_prism_gsa 
#            2789                       885                       780 
# helmsley_xavier_prism_gsa         hyams_protect_gsa                 italy_gsa 
# 1297                       418                      1016 
# kiel_austria_sibdcs_gsa           lewis_sparc_gsa             lithuania_gsa 
# 14660                      2859                      2277 
# mccauley_gsa          mccauley_new_gsa              mcgovern_gsa 
# 788                      1628                      6050 
# moayyedi_imagine_gsa           netherlands_gsa        newberry_share_gsa 
#                 1128                      4705                       865 
# niddk_broad_gsa             niddk_cho_gsa           niddk_duerr_gsa 
#            5515                      1764                      1944 
# niddk_feinstein_gsa           niddk_rioux_gsa      niddk_silverberg_gsa 
# 8301                       919                      2371 
# palotie_hus_gsa           pekow_share_gsa             prism_nfe_gsa 
# 878                       634                       466 
# prism_nfe_gwas        rioux_igenomed_gsa           sands_msccr_gsa 
# 862                       182                      1430 
# slovenia_gsa              stampfer_gsa                sweden_gsa 
# 270                      1477                      1414 
# vermeire_gsa               weersma_gsa          xavier_prism_gsa 
# 4713                       709                       692 
# xavier_share_gsa 
# 696 

# double check all samples included in exercise
famall<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch.fam",sep=""),head=F)

dim(famall)
# 92298

dim(dat)
# 92298

dim(dat[which(dat$V1 %in% famall$V1),])
# [1] 92298     7 # all OK


colnames(dat)[6]<-"pheno"

aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
#                      Group.1  x.1  x.2
# 1                 basque_gsa  987  529
# 2    belgium_franchimont_gsa  611  921
# 3          belgium_louis_gsa  606  925
# 4       belgium_vermeire_gsa  822 3192
# 5              bernstein_gsa   10  504
# 6                   ccfa_gsa    0 2188
# 7                 cedars_gsa  951 2145
# 8               farkkila_gsa    0   68
# 9           finland_illugwas    0  457
# 10           franchimont_gsa 1489 1300
# 11                franke_gsa  448  437
# 12        helmsley_prism_gsa  261  519
# 13 helmsley_xavier_prism_gsa  160 1137
# 14         hyams_protect_gsa    0  418
# 15                 italy_gsa  396  620
# 16   kiel_austria_sibdcs_gsa 4680 9980
# 17           lewis_sparc_gsa    0 2859
# 18             lithuania_gsa 1173 1104
# 19              mccauley_gsa    0  788
# 20          mccauley_new_gsa  238 1390
# 21              mcgovern_gsa 1210 4840
# 22      moayyedi_imagine_gsa  251  877
# 23           netherlands_gsa  550 4155
# 24        newberry_share_gsa    0  865
# 25           niddk_broad_gsa 1027 4488
# 26             niddk_cho_gsa  740 1024
# 27           niddk_duerr_gsa  918 1026
# 28       niddk_feinstein_gsa 2446 5855
# 29           niddk_rioux_gsa  583  336
# 30      niddk_silverberg_gsa  544 1827
# 31           palotie_hus_gsa    0  878
# 32           pekow_share_gsa    0  634
# 33             prism_nfe_gsa   33  433
# 34            prism_nfe_gwas  263  599
# 35        rioux_igenomed_gsa    0  182
# 36           sands_msccr_gsa  309 1121
# 37              slovenia_gsa  181   89
# 38              stampfer_gsa 1039  438
# 39                sweden_gsa  997  417
# 40              vermeire_gsa  817 3896
# 41               weersma_gsa  385  324
# 42          xavier_prism_gsa   63  629
# 43          xavier_share_gsa    0  696


colnames(dat)[5]<-"sex"
aggregate(as.factor(dat$sex), by=list(dat$cohort), FUN=summary)
#                      Group.1  x.0  x.1  x.2
# 1                 basque_gsa    5  931  580
# 2    belgium_franchimont_gsa    1  666  865
# 3          belgium_louis_gsa   21  662  848
# 4       belgium_vermeire_gsa    0 1934 2080
# 5              bernstein_gsa  217  124  173
# 6                   ccfa_gsa    0 1005 1183
# 7                 cedars_gsa    0 1496 1600
# 8               farkkila_gsa    0   34   34
# 9           finland_illugwas    0  301  156
# 10           franchimont_gsa   71 1274 1444
# 11                franke_gsa    0  508  377
# 12        helmsley_prism_gsa    2  382  396
# 13 helmsley_xavier_prism_gsa   22  580  695
# 14         hyams_protect_gsa    0  212  206
# 15                 italy_gsa    0  551  465
# 16   kiel_austria_sibdcs_gsa    2 7250 7408
# 17           lewis_sparc_gsa    0 1291 1568
# 18             lithuania_gsa    0 1206 1071
# 19              mccauley_gsa    0  404  384
# 20          mccauley_new_gsa    0  845  783
# 21              mcgovern_gsa   18 2965 3067
# 22      moayyedi_imagine_gsa    0  490  638
# 23           netherlands_gsa    0 1947 2758
# 24        newberry_share_gsa    0  369  496
# 25           niddk_broad_gsa    0 2858 2657
# 26             niddk_cho_gsa    4  905  855
# 27           niddk_duerr_gsa    3  969  972
# 28       niddk_feinstein_gsa    8 4161 4132
# 29           niddk_rioux_gsa    2  420  497
# 30      niddk_silverberg_gsa    5 1212 1154
# 31           palotie_hus_gsa    0  420  458
# 32           pekow_share_gsa    0  317  317
# 33             prism_nfe_gsa    0  218  248
# 34            prism_nfe_gwas    0  403  459
# 35        rioux_igenomed_gsa    0   87   95
# 36           sands_msccr_gsa    1  765  664
# 37              slovenia_gsa    1   59  210
# 38              stampfer_gsa    0  239 1238
# 39                sweden_gsa    0  821  593
# 40              vermeire_gsa    0 2280 2433
# 41               weersma_gsa    0  355  354
# 42          xavier_prism_gsa    0  343  349
# 43          xavier_share_gsa    0  348  348

table(dat$V1==dat$V2)
# TRUE 
# 92298
kin<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch.kin0",sep=""),head=T)
table(kin$InfType)
# 2nd    3rd    4th Dup/MZ     FS     PO     UN 
# 455    106     42  27204   1778   3303     27 

summary(kin$Kinship)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.8063  0.4994  0.5000  0.4554  0.5000  0.5000 

# several <0 Kinship - low qual?
dim(kin[which(kin$Kinship<0),])
# [1] 29 14

smissing<-read.table("/path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch.imiss",head=T)
summary(smissing$F_MISS)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000000 0.0003114 0.0006228 0.0031901 0.0016750 0.7933000


dim(smissing[which(smissing$F_MISS>=0.05),])
# [1] 459   6

dim(kin[which(kin$Kinship<0 & 
                ((kin$FID1 %in% smissing$FID[which(smissing$F_MISS>=0.05)]) | 
                (kin$FID2 %in% smissing$FID[which(smissing$F_MISS>=0.05)])) ),])
# [1] 11 14

dim(kin[which(kin$Kinship<0 & 
                ((kin$FID1 %in% smissing$FID[which(smissing$F_MISS>=0.10)]) | 
                   (kin$FID2 %in% smissing$FID[which(smissing$F_MISS>=0.10)])) ),])
# [1] 11 14

dim(kin[which(((kin$FID1 %in% smissing$FID[which(smissing$F_MISS>=0.10)]) | 
                   (kin$FID2 %in% smissing$FID[which(smissing$F_MISS>=0.10)])) ),])
# [1] 760  14
summary(kin$Kinship[which(((kin$FID1 %in% smissing$FID[which(smissing$F_MISS>=0.10)]) | 
                 (kin$FID2 %in% smissing$FID[which(smissing$F_MISS>=0.10)])) )])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.5795  0.2123  0.2691  0.2678  0.3284  0.5000 
table(kin$InfType[which(((kin$FID1 %in% smissing$FID[which(smissing$F_MISS>=0.10)]) | 
                             (kin$FID2 %in% smissing$FID[which(smissing$F_MISS>=0.10)])) )])
# 2nd    3rd    4th Dup/MZ     FS     PO     UN 
# 423    106     42     45      7    110     27 

dim(kin[which(((kin$FID1 %in% smissing$FID[which(smissing$F_MISS>=0.10)]) | 
                 (kin$FID2 %in% smissing$FID[which(smissing$F_MISS>=0.10)])) ),])
# [1] 760  14
summary(kin$Kinship[which(((kin$FID1 %in% smissing$FID[which(smissing$F_MISS>=0.10)]) | 
                             (kin$FID2 %in% smissing$FID[which(smissing$F_MISS>=0.10)])) )])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.5795  0.2123  0.2691  0.2678  0.3284  0.5000 
table(kin$InfType[which(((kin$FID1 %in% smissing$FID[which(smissing$F_MISS>=0.10)]) | 
                           (kin$FID2 %in% smissing$FID[which(smissing$F_MISS>=0.10)])) )])
# 2nd    3rd    4th Dup/MZ     FS     PO     UN 
# 423    106     42     45      7    110     27 


lowmiss<-smissing$FID[which(smissing$F_MISS>=0.10)]


table(cut(kin$Kinship,breaks=c(1,0.354,0.177)))
# (0.177,0.354]     (0.354,1] 
#        5559         27293 

h <- hist(kin$Kinship, breaks = "FD", plot = FALSE)
pdf(paste(path,"pre_imputation/QC/relatedness/histogram_iibdgc_merged_old_new_data_from_broad_gsa_plink.pdf",sep=""),height = 7,width = 14)
ggplot(kin, aes(x=Kinship)) + geom_histogram(bins=90) + xlim(0,0.55)
dev.off()

pdf(paste(path,"pre_imputation/QC/relatedness/histogram_iibdgc_merged_old_new_data_from_broad_gsa_plink_zoomin.pdf",sep=""),height = 7,width = 14)
ggplot(kin, aes(x=Kinship)) + geom_histogram(bins=90) + xlim(0,0.55) + ylim(0,5000) + 
  geom_vline(xintercept=0.345, linetype="dashed", color = "red")
dev.off()

kin_tmp<-kin[which( !(kin$FID1 %in% lowmiss) & !(kin$FID2 %in% lowmiss)),]

summary(kin_tmp$InfType)
# 2nd    3rd    4th Dup/MZ     FS     PO     UN 
# 32      0      0  27159   1771   3193      0 
summary(kin$InfType)
# 2nd    3rd    4th Dup/MZ     FS     PO     UN 
# 455    106     42  27204   1778   3303     27 

h <- hist(kin$Kinship, breaks = "FD", plot = FALSE)
pdf(paste(path,"pre_imputation/QC/relatedness/histogram_iibdgc_merged_old_new_data_from_broad_gsa_plink_nolowmissingsamples.pdf",sep=""),height = 7,width = 14)
ggplot(kin_tmp, aes(x=Kinship)) + geom_histogram(bins=90) + xlim(0,0.55)
dev.off()

pdf(paste(path,"pre_imputation/QC/relatedness/histogram_iibdgc_merged_old_new_data_from_broad_gsa_plink_nolowmissingsamples_zoomin.pdf",sep=""),height = 7,width = 14)
ggplot(kin_tmp, aes(x=Kinship)) + geom_histogram(bins=90) + xlim(0,0.55) + ylim(0,5000) + 
  geom_vline(xintercept=0.345, linetype="dashed", color = "red")
dev.off()

# some of the grey area samples are out when high missingness samples are excluded
dim(kin[which(kin$Kinship>0.345 & kin$Kinship<=0.40),])
# [1] 86 14
dim(kin[which(kin_tmp$Kinship>0.345 & kin_tmp$Kinship<=0.40),])
# [1]  6 14

table(kin_tmp$InfType)
# 2nd    3rd    4th Dup/MZ     FS     PO     UN 
#  32      0      0  27159   1771   3193      0 

table(kin_tmp$InfType[which(kin_tmp$Kinship>0.345)])
# 2nd    3rd    4th Dup/MZ     FS     PO     UN 
#   1      0      0  27159      0     12      0 



h <- hist(kin$Kinship, breaks = "FD", plot = FALSE)
pdf(paste(path,"pre_imputation/QC/relatedness/histogram_iibdgc_merged_old_new_data_from_broad_gsa_plink.pdf",sep=""),height = 7,width = 14)
ggplot(kin, aes(x=Kinship)) + geom_histogram(bins=90) + xlim(0,0.55)
dev.off()

pdf(paste(path,"pre_imputation/QC/relatedness/histogram_iibdgc_merged_old_new_data_from_broad_gsa_plink_zoomin.pdf",sep=""),height = 7,width = 14)
ggplot(kin, aes(x=Kinship)) + geom_histogram(bins=90) + xlim(0.2,0.55) + ylim(0,2000) + 
  geom_vline(xintercept=0.345, linetype="dashed", color = "red")
dev.off()



# continue just with the non higmissingness
dup<-kin_tmp[which(kin_tmp$Kinship>0.345 & kin_tmp$InfType=="Dup/MZ"),c("FID1","FID2","Kinship")]
summary(dup$Kinship)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4428  0.4999  0.5000  0.4998  0.5000  0.5000 

dup<-merge(dup,dat[,c("V1","cohort","sex","pheno")],by.x="FID2",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-2):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-2):ncol(dup)],"FID2",sep="_")

dup<-merge(dup,dat[,c("V1","cohort","sex","pheno")],by.x="FID1",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-2):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-2):ncol(dup)],"FID1",sep="_")

x<-table(dup$cohort_FID1,dup$cohort_FID2)
x

# which samples in new cohorts are part of existing ones - current match also includes old-intra/inter cohort

dup$cohort_type_FID1<-"old"
dup$cohort_type_FID1[which(dup$cohort_FID1 %in% new)]<-"new"
table(dup$cohort_FID1,dup$cohort_type_FID1)

dup$cohort_type_FID2<-"old"
dup$cohort_type_FID2[which(dup$cohort_FID2 %in% new)]<-"new"
table(dup$cohort_FID2,dup$cohort_type_FID2)


#################################################################
### MATCHES OLD-OLD TO EXCLUDE

old<-dup[which(dup$cohort_type_FID1=="old" & dup$cohort_type_FID2=="old"),]
dim(old)
# [1] 2900   11
table(old$cohort_FID1,old$cohort_FID2)
# 
summary(old$Kinship)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.4428  0.4995  0.4998  0.4994  0.5000  0.5000 


#################################################################
### MATCHES NEW-NEW TO INSPECT:

new<-dup[which(dup$cohort_type_FID1=="new" & dup$cohort_type_FID2=="new"),]
dim(new)
# [1] 1613   11

table(new$cohort_FID1,new$cohort_FID2)
#####################################
#                      bernstein_gsa franchimont_gsa franke_gsa
# bernstein_gsa                   33               0          0
# farkkila_gsa                     0               0          0
# franchimont_gsa                  0              61          0
# franke_gsa                       0               0          7
# helmsley_prism_gsa               0               0          0
# hyams_protect_gsa                0               0          0
# lewis_sparc_gsa                  0               0          0
# mccauley_new_gsa                 0               0          0
# mcgovern_gsa                     0               0          0
# newberry_share_gsa               0               0          0
# niddk_cho_gsa                    0               0          0
# niddk_duerr_gsa                  0               0          0
# niddk_rioux_gsa                  0               0          0
# niddk_silverberg_gsa             0               0          0
# pekow_share_gsa                  0               0          0
# rioux_igenomed_gsa               0               0          0
# sands_msccr_gsa                  0               0          0
# vermeire_gsa                     0               0          0
# xavier_prism_gsa                 0               0          0
# xavier_share_gsa                 0               0          0
# 
#                      helmsley_prism_gsa helmsley_xavier_prism_gsa
# bernstein_gsa                         0                         0
# farkkila_gsa                          0                         0
# franchimont_gsa                       0                         0
# franke_gsa                            0                         0
# helmsley_prism_gsa                    1                         4
# hyams_protect_gsa                     0                         0
# lewis_sparc_gsa                      12                        29
# mccauley_new_gsa                      2                         1
# mcgovern_gsa                          5                         3
# newberry_share_gsa                    0                         0
# niddk_cho_gsa                         1                         6
# niddk_duerr_gsa                       0                         0
# niddk_rioux_gsa                       0                         0
# niddk_silverberg_gsa                  0                         1
# pekow_share_gsa                       0                         1
# rioux_igenomed_gsa                    0                         0
# sands_msccr_gsa                       0                         3
# vermeire_gsa                          0                         0
# xavier_prism_gsa                     23                         0
# xavier_share_gsa                    163                       322
# 
#                      hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa
# bernstein_gsa                        0               0                0
# farkkila_gsa                         0               0                0
# franchimont_gsa                      0               0                0
# franke_gsa                           0               0                0
# helmsley_prism_gsa                   0               0                0
# hyams_protect_gsa                    0               0                0
# lewis_sparc_gsa                      0               2                0
# mccauley_new_gsa                     0               1                9
# mcgovern_gsa                         0               0                0
# newberry_share_gsa                   0               9                0
# niddk_cho_gsa                        0               7                0
# niddk_duerr_gsa                      1             119                0
# niddk_rioux_gsa                      0               0                0
# niddk_silverberg_gsa                 0               0                0
# pekow_share_gsa                      0              23                0
# rioux_igenomed_gsa                   0               0                0
# sands_msccr_gsa                      0               0                2
# vermeire_gsa                         0               0                0
# xavier_prism_gsa                     0               0                0
# xavier_share_gsa                     0               6                0
# 
#                      mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa
# bernstein_gsa                   0                   20                  0
# farkkila_gsa                    0                    0                  0
# franchimont_gsa                 0                    0                  0
# franke_gsa                      0                    0                  0
# helmsley_prism_gsa              0                    0                  0
# hyams_protect_gsa               1                    1                  0
# lewis_sparc_gsa                 5                    0                  0
# mccauley_new_gsa               12                    0                  0
# mcgovern_gsa                   46                    4                  0
# newberry_share_gsa              3                    0                  2
# niddk_cho_gsa                   9                    0                  0
# niddk_duerr_gsa                 2                    0                  0
# niddk_rioux_gsa                 0                   13                  0
# niddk_silverberg_gsa            0                    7                  0
# pekow_share_gsa                 5                    0                  0
# rioux_igenomed_gsa              0                    6                  0
# sands_msccr_gsa                10                    0                  0
# vermeire_gsa                    1                    0                  0
# xavier_prism_gsa                0                    0                  0
# xavier_share_gsa                6                    0                  0
# 
#                      niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa
# bernstein_gsa                    0               0               0
# farkkila_gsa                     0               0               0
# franchimont_gsa                  0               0               0
# franke_gsa                       0               0               0
# helmsley_prism_gsa               0               0               0
# hyams_protect_gsa                0               0               0
# lewis_sparc_gsa                  0               0               0
# mccauley_new_gsa                 0               0               0
# mcgovern_gsa                     0               0               1
# newberry_share_gsa               0               0               0
# niddk_cho_gsa                   16               2               1
# niddk_duerr_gsa                  0             121               0
# niddk_rioux_gsa                  0               0               5
# niddk_silverberg_gsa             0               0               0
# pekow_share_gsa                  0               0               0
# rioux_igenomed_gsa               0               0               8
# sands_msccr_gsa                  0               0               0
# vermeire_gsa                     0               0               0
# xavier_prism_gsa                 0               0               0
# xavier_share_gsa                 0               0               0
# 
#                      niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa
# bernstein_gsa                           2               0               0
# farkkila_gsa                            0               1               0
# franchimont_gsa                         0               0               0
# franke_gsa                              0               0               0
# helmsley_prism_gsa                      0               0               0
# hyams_protect_gsa                       0               0               0
# lewis_sparc_gsa                         0               0               0
# mccauley_new_gsa                        0               0               0
# mcgovern_gsa                            2               0               0
# newberry_share_gsa                      0               0               0
# niddk_cho_gsa                           2               0              24
# niddk_duerr_gsa                         7               0               0
# niddk_rioux_gsa                         1               0               0
# niddk_silverberg_gsa                   40               0               0
# pekow_share_gsa                         0               0               1
# rioux_igenomed_gsa                      0               0               0
# sands_msccr_gsa                         0               0               0
# vermeire_gsa                            0               0               0
# xavier_prism_gsa                        0               0               0
# xavier_share_gsa                        0               0               0
# 
#                      rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa
# bernstein_gsa                         0               0            0
# farkkila_gsa                          0               0            0
# franchimont_gsa                       0               0            0
# franke_gsa                            0               0            0
# helmsley_prism_gsa                    0               0            0
# hyams_protect_gsa                     0               0            0
# lewis_sparc_gsa                       0               0            0
# mccauley_new_gsa                      0               0            0
# mcgovern_gsa                          0               0            0
# newberry_share_gsa                    0               1            0
# niddk_cho_gsa                         0             256            1
# niddk_duerr_gsa                       0               0            1
# niddk_rioux_gsa                       0               0            0
# niddk_silverberg_gsa                  0               0            0
# pekow_share_gsa                       0               1            0
# rioux_igenomed_gsa                    6               0            0
# sands_msccr_gsa                       0              10            0
# vermeire_gsa                          0               0            0
# xavier_prism_gsa                      0               0            0
# xavier_share_gsa                      0               0            0
# 
#                      vermeire_gsa xavier_prism_gsa xavier_share_gsa
# bernstein_gsa                   0                0                0
# farkkila_gsa                    0                0                0
# franchimont_gsa                 9                0                0
# franke_gsa                      0                0                0
# helmsley_prism_gsa              0                0                0
# hyams_protect_gsa               0                0                0
# lewis_sparc_gsa                 0                8                0
# mccauley_new_gsa                0                1                1
# mcgovern_gsa                    0                0                0
# newberry_share_gsa              0                0                0
# niddk_cho_gsa                   1                2                3
# niddk_duerr_gsa                 0                1                0
# niddk_rioux_gsa                 0                0                0
# niddk_silverberg_gsa            0                0                0
# pekow_share_gsa                 0                0                1
# rioux_igenomed_gsa              0                0                0
# sands_msccr_gsa                 0                0                1
# vermeire_gsa                   17                0                0
# xavier_prism_gsa                0                2                0
# xavier_share_gsa                0               43                3
#####################################

x<-table(new$cohort_FID1,new$cohort_FID2)
write.table(x,paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_new_new_broad_batch_table_studies.txt",sep=""),
            col.names=T,row.names=T,quote=F,sep="\t")



summary(new$Kinship)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4765  0.4991  0.4995  0.4991  0.4996  0.4999 

dim(new[which(new$Kinship<=0.49),])
# [1] 14 11

table(new$sex_FID1,new$sex_FID2)
#     0   1   2
# 0   3   6  11
# 1  18 680  96
# 2  27  89 683


#################################################################
### MATCHES OLD-NEW TO KEEP

dim(dup)
# [1] 27159    11

nrow(dup)-nrow(old)
# [1] 24259

dup<-dup[which(! (dup$cohort_type_FID1=="old" & dup$cohort_type_FID2=="old")),]
dim(dup)
# [1] 24259    11

dup2<-dup[which(dup$cohort_type_FID2=="old"),]

dup1<-dup[which(dup$cohort_type_FID1=="old"),]
colnames(dup1)
# [1] "FID1"             "FID2"             "Kinship"          "cohort_FID2"     
# [5] "sex_FID2"         "pheno_FID2"       "cohort_FID1"      "sex_FID1"        
# [9] "pheno_FID1"       "cohort_type_FID1" "cohort_type_FID2"

colnames(dup1)<-c("FID2","FID1","Kinship",
              "cohort_FID1","sex_FID1","pheno_FID1",
              "cohort_FID2","sex_FID2","pheno_FID2",
              "cohort_type_FID2","cohort_type_FID1")
dup1<-dup1[,colnames(dup2)]

dup<-rbind(dup1,dup2)
dim(dup)
# [1] 22646    11

24259-nrow(new)
# [1] 22646 OK

table(dup$cohort_type_FID1,dup$cohort_type_FID2)
#       old
# new 22646

dup$cohort_FID1<-as.character(dup$cohort_FID1)
dup$cohort_FID2<-as.character(dup$cohort_FID2)
table(dup$cohort_FID1,dup$cohort_FID2)

summary(dup$Kinship)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4428  0.5000  0.5000  0.4999  0.5000  0.5000

dim(dup[which(dup$Kinship<=0.49),])
# [1] 33 11

table(dup$cohort_FID1[which(dup$Kinship<=0.49)],
      dup$cohort_FID2[which(dup$Kinship<=0.49)])

## what replaces what - if that possible

x<-table(dup$cohort_FID1,dup$cohort_FID2)
write.table(x,paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch_table_studies.txt",sep=""),
            col.names=T,row.names=T,quote=F,sep="\t")

# for each of the old studies that can be potentially replaced by a new study, create a id match (old vs new) map:

old_studies<-names(table(dup$cohort_FID2))
old_studies<-old_studies[which(!old_studies %in% c("belgium_louis_gsa","finland_illugwas","kiel_austria_sibdcs_gsa"))]


for (i in 1:length(old_studies)) {
  
  fam<-dat[which(dat$cohort==old_studies[i]),]
  fam<-fam[,c("V1","sex","pheno","cohort")]
  colnames(fam)[1]<-"Sample_id"
  
  print(paste("N samples in old",old_studies[i],":",nrow(fam)))
  
  tmp<-dup[which(dup$cohort_FID2==old_studies[i]),]
  
  fam$new_cohort<-NA
  fam$new_cohort_sample_id<-NA
  fam$new_cohort_sex<-NA
  fam$new_cohort_pheno<-NA

  for (j in 1:nrow(fam)) {
    fam$new_cohort[j]<-paste(tmp$cohort_FID1[which(tmp$FID2==fam$Sample_id[j])],collapse="|")
    fam$new_cohort_sample_id[j]<-paste(tmp$FID1[which(tmp$FID2==fam$Sample_id[j])],collapse="|")
    fam$new_cohort_sex[j]<-paste(tmp$sex_FID1[which(tmp$FID2==fam$Sample_id[j])],collapse="|")
    fam$new_cohort_pheno[j]<-paste(tmp$pheno_FID1[which(tmp$FID2==fam$Sample_id[j])],collapse="|")
  }
  
  # print("New cohort:")
  # print(table(fam$new_cohort))
  # print("Concordance sex old new:")
  # print(table(fam$sex,fam$new_cohort_sex))
  # 
  # print("Concordance phenotype old new:")
  # print(table(fam$pheno,fam$new_cohort_pheno))
  
  print(paste("Samples with no new id :",nrow(fam[which(fam$new_cohort==""),])))
  
  if (!old_studies[i] %in% c("niddk_feinstein_gsa","netherlands_gsa")) {
    
    if(!exists("tocheck")) {
      tocheck<-fam[which(fam$new_cohort==""),]
    }else{
      tocheck<-rbind(tocheck,fam[which(fam$new_cohort==""),])
    }
  }
  fwrite(fam,paste(path,"pre_imputation/QC/relatedness/",old_studies[i],"_map_old_plus_new_ids.tsv.gz",sep=""),
                        col.names=T,row.names=T,quote=F,sep="\t")
  
}



#######################################################
# [1] "N samples in old belgium_franchimont_gsa : 1532"
# [1] "Samples with no new id : 2"

# [1] "N samples in old belgium_vermeire_gsa : 4014"
# [1] "Samples with no new id : 118"

# [1] "N samples in old ccfa_gsa : 2188"
# [1] "Samples with no new id : 1"

# [1] "N samples in old cedars_gsa : 3096"
# [1] "Samples with no new id : 20"

# [1] "N samples in old mccauley_gsa : 788"
# [1] "Samples with no new id : 6"

# [1] "N samples in old netherlands_gsa : 4705"
# [1] "Samples with no new id : 3996"


# [1] "N samples in old niddk_broad_gsa : 5515"
# [1] "Samples with no new id : 177"

# [1] "N samples in old niddk_feinstein_gsa : 8301"
# [1] "Samples with no new id : 6059"

# [1] "N samples in old prism_nfe_gsa : 466"
# [1] "Samples with no new id : 42"

# [1] "N samples in old prism_nfe_gwas : 862"
# [1] "Samples with no new id : 57"

#######################################################

# cd /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/
# mkdir sample_maps_old_new_based_on_kinship_20220818/
# mv *_map_old_plus_new_ids.tsv.gz sample_maps_old_new_based_on_kinship_20220818/
# 
# tar -zcvf sample_maps_old_new_based_on_kinship_20220818.tar.gz sample_maps_old_new_based_on_kinship_20220818/
# 
# ~/Desktop/

table(tocheck$cohort)
# belgium_franchimont_gsa    belgium_vermeire_gsa                ccfa_gsa 
#                       2                     118                       1 
# cedars_gsa            mccauley_gsa         niddk_broad_gsa 
#         20                       6                     177 
# prism_nfe_gsa          prism_nfe_gwas 
#            42                      57

tocheck<-merge(tocheck,smissing[,c("FID","N_MISS","N_GENO","F_MISS")],by.x="Sample_id",by.y="FID",all.x=T)
summary(tocheck$F_MISS)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000084 0.0006354 0.0015910 0.0230185 0.0045400 0.7302000 

# not entirely caused by variant misingess


############# look at the samples with old-new IIBDGC ids and the old new maps to their supplier IDs:


########################################################################################################
########################################################################################################

# PART 2

###########################
# 1.- SELECT THE VARIANTS #
###########################

# SEE SCRIPT:
# create_list_common_variants_among_cohorts - USE SAME VARIANTS AS IN STEP 12.3 create_list_common_variants_among_cohorts 
path=/path/to/ibdgwas/IIBDGC/
  wc -l ${path}pre_imputation/QC/relatedness/list_common_variants_old_new_broad_batch_2
# 108295 

########################
# 2.- EXTRACT VARIANTS #
########################

studies=(basque_gsa belgium_franchimont_gsa belgium_louis_gsa belgium_vermeire_gsa ccfa_gsa cedars_gsa italy_gsa kiel_austria_sibdcs_gsa lithuania_gsa mccauley_gsa netherlands_gsa niddk_broad_gsa niddk_feinstein_gsa prism_nfe_gsa slovenia_gsa sweden_gsa prism_nfe_gwas finland_illugwas bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip \
--allow-no-sex \
--extract ${path}pre_imputation/QC/relatedness/list_common_variants_old_new_broad_batch_2 \
--make-bed \
--out ${path}pre_imputation/QC/relatedness/${i}_subset_old_new_broad_batch
done

wc -l ${path}pre_imputation/QC/relatedness/*_subset_old_new_broad_batch.fam
###########
# 1491 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/basque_gsa_subset_old_new_broad_batch.fam
# 1488 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_franchimont_gsa_subset_old_new_broad_batch.fam
# 1506 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_louis_gsa_subset_old_new_broad_batch.fam
# 3979 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_vermeire_gsa_subset_old_new_broad_batch.fam
# 481 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/bernstein_gsa_subset_old_new_broad_batch.fam
# 2064 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/ccfa_gsa_subset_old_new_broad_batch.fam
# 1803 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_gsa_subset_old_new_broad_batch.fam
# 68 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/farkkila_gsa_subset_old_new_broad_batch.fam
# 444 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/finland_illugwas_subset_old_new_broad_batch.fam
# 2727 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/franchimont_gsa_subset_old_new_broad_batch.fam
# 860 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/franke_gsa_subset_old_new_broad_batch.fam
# 760 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/helmsley_prism_gsa_subset_old_new_broad_batch.fam
# 1293 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/helmsley_xavier_prism_gsa_subset_old_new_broad_batch.fam
# 418 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/hyams_protect_gsa_subset_old_new_broad_batch.fam
# 949 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/italy_gsa_subset_old_new_broad_batch.fam
# 13870 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/kiel_austria_sibdcs_gsa_subset_old_new_broad_batch.fam
# 2857 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/lewis_sparc_gsa_subset_old_new_broad_batch.fam
# 2210 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/lithuania_gsa_subset_old_new_broad_batch.fam
# 769 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mccauley_gsa_subset_old_new_broad_batch.fam
# 1612 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mccauley_new_gsa_subset_old_new_broad_batch.fam
# 5960 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mcgovern_gsa_subset_old_new_broad_batch.fam
# 1128 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/moayyedi_imagine_gsa_subset_old_new_broad_batch.fam
# 4555 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/netherlands_gsa_subset_old_new_broad_batch.fam
# 861 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/newberry_share_gsa_subset_old_new_broad_batch.fam
# 5122 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_broad_gsa_subset_old_new_broad_batch.fam
# 1746 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_cho_gsa_subset_old_new_broad_batch.fam
# 1817 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_duerr_gsa_subset_old_new_broad_batch.fam
# 7960 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_feinstein_gsa_subset_old_new_broad_batch.fam
# 913 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_rioux_gsa_subset_old_new_broad_batch.fam
# 2330 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_silverberg_gsa_subset_old_new_broad_batch.fam
# 878 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/palotie_hus_gsa_subset_old_new_broad_batch.fam
# 633 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/pekow_share_gsa_subset_old_new_broad_batch.fam
# 464 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/prism_nfe_gsa_subset_old_new_broad_batch.fam
# 735 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/prism_nfe_gwas_subset_old_new_broad_batch.fam
# 176 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/rioux_igenomed_gsa_subset_old_new_broad_batch.fam
# 1422 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/sands_msccr_gsa_subset_old_new_broad_batch.fam
# 262 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/slovenia_gsa_subset_old_new_broad_batch.fam
# 1468 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/stampfer_gsa_subset_old_new_broad_batch.fam
# 1379 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/sweden_gsa_subset_old_new_broad_batch.fam
# 4683 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/vermeire_gsa_subset_old_new_broad_batch.fam
# 709 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/weersma_gsa_subset_old_new_broad_batch.fam
# 687 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/xavier_prism_gsa_subset_old_new_broad_batch.fam
# 693 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/xavier_share_gsa_subset_old_new_broad_batch.fam
# 88230 total

###########

wc -l ${path}pre_imputation/QC/relatedness/*_subset_old_new_broad_batch.bim
# 108295  OK

##############################################################################################################################################

#######################
# 3.- MERGE ALL FILES #
#######################

#### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("basque_gsa","belgium_franchimont_gsa","belgium_louis_gsa","belgium_vermeire_gsa","ccfa_gsa","cedars_gsa",
           "italy_gsa","kiel_austria_sibdcs_gsa","lithuania_gsa","mccauley_gsa","netherlands_gsa",
           "niddk_broad_gsa","niddk_feinstein_gsa","prism_nfe_gsa","slovenia_gsa","sweden_gsa",
           "prism_nfe_gwas","finland_illugwas",
           "helmsley_prism_gsa","helmsley_xavier_prism_gsa",
           "bernstein_gsa","farkkila_gsa","franchimont_gsa","franke_gsa",
           "hyams_protect_gsa","lewis_sparc_gsa","mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa",
           "newberry_share_gsa","niddk_cho_gsa","niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa",
           "palotie_hus_gsa","pekow_share_gsa","rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa",
           "vermeire_gsa","weersma_gsa","xavier_prism_gsa","xavier_share_gsa")


length(cohorts)
# [1] 43

#####

cohorts<-cohorts[-1]

dat<-matrix(ncol=3,nrow=length(cohorts))
dat<-as.data.frame(dat)

for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_old_new_broad_batch.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_old_new_broad_batch.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_old_new_broad_batch.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/relatedness/list_cohorts_tomerge_old_new_data_from_broad_gsa.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#############################

# hpc-server

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/relatedness/basque_gsa_subset_old_new_broad_batch \
--allow-no-sex \
--merge-list ${path}pre_imputation/QC/relatedness/list_cohorts_tomerge_old_new_data_from_broad_gsa.txt \
--make-bed --out ${path}pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch

# 108295 variants and 88230 people pass filters and QC.



##############################################################################################################################################

############################
# 3.- ESTIMATE RELATEDNESS #
############################

# The current release is KING version 2.2.4 (released on October 11, 2019). 

# --related provides integrative, fast, and accurate inference for close relationships. This option is highly recommended, 
# especially when dealing with biobank-level datasets. Integration of the IBD segment inference furthur improves the inference accuracy. 
# In our tests, --related has been successfully applied to datasets consisting of ~10 million samples. When "--rplot" is specified, 
# several relationship plots are generated automatically. --related --degree 2 specifies that only related pairs
# (up to the 2nd-degree in this case) between families are included in the output. Specifically all pairs across families with 
# a Kinship coefficient less than 0.0884 will be excluded from the output. More details of 
# --related analysis are available in the INTEGRATED RELATIONSHIP INFERENCE section later in this tutorial.

# --duplicate implements the fastest (and accurate) algorithm to identify duplicates or MZ twins. The running time is in seconds, 
# unless the number of samples is > 1,000,000 in which case a few minutes may be needed. In our tests, --duplicate has been successfully 
# applied to datasets consisting of ~10 million samples. One potential application of the duplicate analysis is to identify duplicates 
# accross different studies, in which case multiple datasets can be read in conveniently as shown in GENERAL INPUT FILES section.

# The amount of computer memory required by KING analysis is modest, at ~N ✕ M / 4 (where N is the number of samples and M is the number of SNPs)
# plus a small percentage of overhead cost. E.g., for a dataset consisting of 100,000 samples each genotyped at 1,000,000 SNPs, 
# the required memory size is ~25GB.

# (3110*107598/4)/1E+6
# [1] 83.65744


#################


# The current release is KING version 2.2.4 (released on October 11, 2019). 

# --related provides integrative, fast, and accurate inference for close relationships. This option is highly recommended, 
# especially when dealing with biobank-level datasets. Integration of the IBD segment inference furthur improves the inference accuracy. 
# In our tests, --related has been successfully applied to datasets consisting of ~10 million samples. When "--rplot" is specified, 
# several relationship plots are generated automatically. --related --degree 2 specifies that only related pairs
# (up to the 2nd-degree in this case) between families are included in the output. Specifically all pairs across families with 
# a Kinship coefficient less than 0.0884 will be excluded from the output. More details of 
# --related analysis are available in the INTEGRATED RELATIONSHIP INFERENCE section later in this tutorial.

# --duplicate implements the fastest (and accurate) algorithm to identify duplicates or MZ twins. The running time is in seconds, 
# unless the number of samples is > 1,000,000 in which case a few minutes may be needed. In our tests, --duplicate has been successfully 
# applied to datasets consisting of ~10 million samples. One potential application of the duplicate analysis is to identify duplicates 
# accross different studies, in which case multiple datasets can be read in conveniently as shown in GENERAL INPUT FILES section.

# The amount of computer memory required by KING analysis is modest, at ~N ✕ M / 4 (where N is the number of samples and M is the number of SNPs)
# plus a small percentage of overhead cost. E.g., for a dataset consisting of 100,000 samples each genotyped at 1,000,000 SNPs, 
# the required memory size is ~25GB.


# (92298*118825/4)/1E+6
# [1] 2741.827

path_ukbb=/path/to/project
  MEM=5000

bsub -J"kg" -M"$MEM" -n8 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path}pre_imputation/QC/relatedness/log/stderr_king_iibdgc_merged_old_new_broad_batch \
-o ${path}pre_imputation/QC/relatedness/log/stdout_king_iibdgc_merged_old_new_broad_batch \
"/path/to/software/./king \
-b ${path}pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch.bed \
--related --cpus 8 --prefix ${path}pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch"
# Job <23607> is submitted to queue <normal>.

# Relationship summary (total relatives: 0 by pedigree, 32846 by inference)
#                 MZ      PO      FS      2nd
# =====================================================
#   Inference     20961   2765    1539    17
# 
# Between-family relatives (Kinship >= 0.17678) saved in file /path/to/ibdgwas/IIBDGC/pre_
# imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch.kin0
# 
# Note only duplicates and 1st-degree relatives are included in the inference.
# Specifying '--degree 2' if a higher degree relationship inference is needed.


# previous attempt using 3100 markers - shared between all old iibdgc (2020)
# ${path}pre_imputation/QC/relatedness/iibdgc_merged_old_new_data_from_broad_gsa_plink


###################################
# label samples with low call rate - all samples regardless QC are kept


/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch \
--missing --out ${path}pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch
# Sample missing data report written to
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch.imiss,
# and variant-based missing data report written to
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch.lmiss.


################################

#### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("basque_gsa","belgium_franchimont_gsa","belgium_louis_gsa","belgium_vermeire_gsa","ccfa_gsa","cedars_gsa",
           "italy_gsa","kiel_austria_sibdcs_gsa","lithuania_gsa","mccauley_gsa","netherlands_gsa",
           "niddk_broad_gsa","niddk_feinstein_gsa","prism_nfe_gsa","slovenia_gsa","sweden_gsa",
           "prism_nfe_gwas","finland_illugwas",
           "helmsley_prism_gsa","helmsley_xavier_prism_gsa",
           "bernstein_gsa","farkkila_gsa","franchimont_gsa","franke_gsa",
           "hyams_protect_gsa","lewis_sparc_gsa","mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa",
           "newberry_share_gsa","niddk_cho_gsa","niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa",
           "palotie_hus_gsa","pekow_share_gsa","rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa",
           "vermeire_gsa","weersma_gsa","xavier_prism_gsa","xavier_share_gsa")

new<-c("helmsley_prism_gsa","helmsley_xavier_prism_gsa",
       "bernstein_gsa","farkkila_gsa","franchimont_gsa","franke_gsa",
       "hyams_protect_gsa","lewis_sparc_gsa","mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa",
       "newberry_share_gsa","niddk_cho_gsa","niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa",
       "palotie_hus_gsa","pekow_share_gsa","rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa",
       "vermeire_gsa","weersma_gsa","xavier_prism_gsa","xavier_share_gsa")

length(cohorts)
# [1] 43

for (i in 1:length(cohorts)) {
  
  tmp<-read.table(paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_old_new_broad_batch.fam",sep=""),head=F)
  
  tmp$cohort<-as.character(cohorts[i])
  tmp$V1<-as.character(tmp$V1)
  tmp$V2<-as.character(tmp$V2)
  
  if(i==1){
    dat<-tmp
  }else{
    dat<-rbind(dat,tmp)
  }
}

sample_miss<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch.imiss",sep=""),head=T)
dat<-merge(dat,sample_miss[,c("IID","F_MISS")],by.x="V2",by.y="IID",all.x=T)

table(dat$cohort)
# basque_gsa   belgium_franchimont_gsa         belgium_louis_gsa 
# 1491                      1488                      1506 
# belgium_vermeire_gsa             bernstein_gsa                  ccfa_gsa 
# 3979                       481                      2064 
# cedars_gsa              farkkila_gsa          finland_illugwas 
# 1803                        68                       444 
# franchimont_gsa                franke_gsa        helmsley_prism_gsa 
# 2727                       860                       760 
# helmsley_xavier_prism_gsa         hyams_protect_gsa                 italy_gsa 
# 1293                       418                       949 
# kiel_austria_sibdcs_gsa           lewis_sparc_gsa             lithuania_gsa 
# 13870                      2857                      2210 
# mccauley_gsa          mccauley_new_gsa              mcgovern_gsa 
# 769                      1612                      5960 
# moayyedi_imagine_gsa           netherlands_gsa        newberry_share_gsa 
# 1128                      4555                       861 
# niddk_broad_gsa             niddk_cho_gsa           niddk_duerr_gsa 
# 5122                      1746                      1817 
# niddk_feinstein_gsa           niddk_rioux_gsa      niddk_silverberg_gsa 
# 7960                       913                      2330 
# palotie_hus_gsa           pekow_share_gsa             prism_nfe_gsa 
# 878                       633                       464 
# prism_nfe_gwas        rioux_igenomed_gsa           sands_msccr_gsa 
# 735                       176                      1422 
# slovenia_gsa              stampfer_gsa                sweden_gsa 
# 262                      1468                      1379 
# vermeire_gsa               weersma_gsa          xavier_prism_gsa 
# 4683                       709                       687 
# xavier_share_gsa 
# 693 

# double check all samples included in exercise
famall<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch.fam",sep=""),head=F)

dim(famall)
# 88230

dim(dat)
# 88230

dim(dat[which(dat$V1 %in% famall$V1),])
# [1] 88230     7 # all OK


colnames(dat)[6]<-"pheno"

aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
#                      Group.1  x.1  x.2
# 1                 basque_gsa  966  525
# 2    belgium_franchimont_gsa  592  896
# 3          belgium_louis_gsa  599  907
# 4       belgium_vermeire_gsa  813 3166
# 5              bernstein_gsa   10  471
# 6                   ccfa_gsa    0 2064
# 7                 cedars_gsa  757 1046
# 8               farkkila_gsa    0   68
# 9           finland_illugwas    0  444
# 10           franchimont_gsa 1468 1259
# 11                franke_gsa  431  429
# 12        helmsley_prism_gsa  259  501
# 13 helmsley_xavier_prism_gsa  160 1133
# 14         hyams_protect_gsa    0  418
# 15                 italy_gsa  374  575
# 16   kiel_austria_sibdcs_gsa 4455 9415
# 17           lewis_sparc_gsa    0 2857
# 18             lithuania_gsa 1139 1071
# 19              mccauley_gsa    0  769
# 20          mccauley_new_gsa  238 1374
# 21              mcgovern_gsa 1183 4777
# 22      moayyedi_imagine_gsa  251  877
# 23           netherlands_gsa  548 4007
# 24        newberry_share_gsa    0  861
# 25           niddk_broad_gsa 1012 4110
# 26             niddk_cho_gsa  737 1009
# 27           niddk_duerr_gsa  871  946
# 28       niddk_feinstein_gsa 2350 5610
# 29           niddk_rioux_gsa  582  331
# 30      niddk_silverberg_gsa  539 1791
# 31           palotie_hus_gsa    0  878
# 32           pekow_share_gsa    0  633
# 33             prism_nfe_gsa   33  431
# 34            prism_nfe_gwas  257  478
# 35        rioux_igenomed_gsa    0  176
# 36           sands_msccr_gsa  308 1114
# 37              slovenia_gsa  175   87
# 38              stampfer_gsa 1032  436
# 39                sweden_gsa  979  400
# 40              vermeire_gsa  813 3870
# 41               weersma_gsa  385  324
# 42          xavier_prism_gsa   62  625
# 43          xavier_share_gsa    0  693

colnames(dat)[5]<-"sex"
aggregate(as.factor(dat$sex), by=list(dat$cohort), FUN=summary)
# Group.1  x.0  x.1  x.2
# 1                 basque_gsa    7  912  572
# 2    belgium_franchimont_gsa    2  642  844
# 3          belgium_louis_gsa   11  650  845
# 4       belgium_vermeire_gsa   26 1903 2050
# 5              bernstein_gsa    2  188  291
# 6                   ccfa_gsa    2  941 1121
# 7                 cedars_gsa    7  875  921
# 8               farkkila_gsa   17   25   26
# 9           finland_illugwas    0  295  149
# 10           franchimont_gsa  113 1213 1401
# 11                franke_gsa   18  487  355
# 12        helmsley_prism_gsa    1  371  388
# 13 helmsley_xavier_prism_gsa    6  582  705
# 14         hyams_protect_gsa    3  210  205
# 15                 italy_gsa    3  515  431
# 16   kiel_austria_sibdcs_gsa  193 6791 6886
# 17           lewis_sparc_gsa   19 1280 1558
# 18             lithuania_gsa    0 1174 1036
# 19              mccauley_gsa    2  394  373
# 20          mccauley_new_gsa   11  827  774
# 21              mcgovern_gsa   11 2922 3027
# 22      moayyedi_imagine_gsa    7  489  632
# 23           netherlands_gsa   11 1893 2651
# 24        newberry_share_gsa    6  366  489
# 25           niddk_broad_gsa   16 2649 2457
# 26             niddk_cho_gsa    9  892  845
# 27           niddk_duerr_gsa   12  909  896
# 28       niddk_feinstein_gsa   13 4000 3947
# 29           niddk_rioux_gsa    8  417  488
# 30      niddk_silverberg_gsa   14 1194 1122
# 31           palotie_hus_gsa    9  414  455
# 32           pekow_share_gsa    2  315  316
# 33             prism_nfe_gsa    4  215  245
# 34            prism_nfe_gwas    0  359  376
# 35        rioux_igenomed_gsa    5   81   90
# 36           sands_msccr_gsa  709  409  304
# 37              slovenia_gsa    5   55  202
# 38              stampfer_gsa    4  236 1228
# 39                sweden_gsa    0  797  582
# 40              vermeire_gsa   22 2259 2402
# 41               weersma_gsa    4  353  352
# 42          xavier_prism_gsa    2  341  344
# 43          xavier_share_gsa    7  343  343

table(dat$V1==dat$V2)
# TRUE 
# 88230
kin<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch.kin0",sep=""),head=T)
table(kin$InfType)
# 2nd    4th Dup/MZ     FS     PO 
# 17      2  20961   1539   2765  

summary(kin$Kinship)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.9967  0.4999  0.5000  0.4570  0.5000  0.5000 

table(cut(kin$Kinship,breaks=c(1,0.354,0.177)))
# (0.177,0.354]     (0.354,1] 
# 4314         20966 

#  <0 Kinship
dim(kin[which(kin$Kinship<0),])
# [1] 4 14

dup<-kin[which(kin$Kinship>0.345 & kin$InfType=="Dup/MZ"),c("FID1","FID2","Kinship")]
dim(dup)
# [1] 20961     3

summary(dup$Kinship)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0 0.4444  0.5000  0.5000  0.4999  0.5000  0.5000 


dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID2",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID2",sep="_")

dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID1",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID1",sep="_")


x<-table(dup$cohort_FID1,dup$cohort_FID2)
x

# which samples in new cohorts are part of existing ones - current match also includes old-intra/inter cohort

dup$cohort_type_FID1<-"old"
dup$cohort_type_FID1[which(dup$cohort_FID1 %in% new)]<-"new"
table(dup$cohort_FID1,dup$cohort_type_FID1)
#                          new  old
# belgium_franchimont_gsa    0 1491
# belgium_louis_gsa          0    4
# belgium_vermeire_gsa       0 3864
# bernstein_gsa             22    0
# ccfa_gsa                   0 2131
# cedars_gsa                 0 1795
# farkkila_gsa               1    0
# finland_illugwas           0   19
# franchimont_gsa            9    0
# helmsley_prism_gsa         4    0
# hyams_protect_gsa          2    0
# kiel_austria_sibdcs_gsa    0   95
# lewis_sparc_gsa           54    0
# mccauley_gsa               0  771
# mccauley_new_gsa          18    0
# mcgovern_gsa              14    0
# netherlands_gsa            0  696
# newberry_share_gsa        13    0
# niddk_broad_gsa            0 5284
# niddk_cho_gsa            305    0
# niddk_duerr_gsa          117    0
# niddk_feinstein_gsa        0 2233
# niddk_rioux_gsa           14    0
# niddk_silverberg_gsa       8    0
# pekow_share_gsa           31    0
# prism_nfe_gsa              0  487
# prism_nfe_gwas             0  890
# rioux_igenomed_gsa        14    0
# sands_msccr_gsa           17    0
# vermeire_gsa               1    0
# xavier_prism_gsa          22    0
# xavier_share_gsa         535    0

dup$cohort_type_FID2<-"old"
dup$cohort_type_FID2[which(dup$cohort_FID2 %in% new)]<-"new"
table(dup$cohort_FID2,dup$cohort_type_FID2)
#                            new  old
# belgium_vermeire_gsa         0    1
# bernstein_gsa                3    0
# farkkila_gsa                 1    0
# franchimont_gsa           1492    0
# franke_gsa                  90    0
# helmsley_prism_gsa         963    0
# helmsley_xavier_prism_gsa  418    0
# hyams_protect_gsa            4    0
# lewis_sparc_gsa           2348    0
# mccauley_new_gsa           779    0
# mcgovern_gsa              3755    0
# moayyedi_imagine_gsa        71    0
# newberry_share_gsa           5    0
# niddk_cho_gsa             1240    0
# niddk_duerr_gsa           1418    0
# niddk_rioux_gsa            735    0
# niddk_silverberg_gsa      1629    0
# palotie_hus_gsa             19    0
# pekow_share_gsa             55    0
# rioux_igenomed_gsa           6    0
# sands_msccr_gsa            659    0
# stampfer_gsa                 3    0
# vermeire_gsa              3878    0
# weersma_gsa                696    0
# xavier_prism_gsa           464    0
# xavier_share_gsa           229    0

#################################################################
### MATCHES OLD-OLD TO EXCLUDE

old<-dup[which(dup$cohort_type_FID1=="old" & dup$cohort_type_FID2=="old"),]
dim(old)
# [1] 1   11
table(old$cohort_FID1,old$cohort_FID2)
#                belgium_vermeire_gsa
# prism_nfe_gwas                    1
# exclude this potential dup



#################################################################
### MATCHES NEW-NEW TO INSPECT:

new<-dup[which(dup$cohort_type_FID1=="new" & dup$cohort_type_FID2=="new"),]
dim(new)
# [1] 1201   11

table(as.character(new$cohort_FID1),as.character(new$cohort_FID2))

#####################################
#                      helmsley_prism_gsa helmsley_xavier_prism_gsa
# bernstein_gsa                         0                         0
# farkkila_gsa                          0                         0
# franchimont_gsa                       0                         0
# helmsley_prism_gsa                    0                         4
# hyams_protect_gsa                     0                         0
# lewis_sparc_gsa                      12                        29
# mccauley_new_gsa                      2                         1
# mcgovern_gsa                          4                         3
# newberry_share_gsa                    0                         0
# niddk_cho_gsa                         1                         6
# niddk_duerr_gsa                       0                         0
# niddk_rioux_gsa                       0                         0
# niddk_silverberg_gsa                  0                         1
# pekow_share_gsa                       0                         1
# rioux_igenomed_gsa                    0                         0
# sands_msccr_gsa                       0                         3
# vermeire_gsa                          0                         0
# xavier_prism_gsa                     22                         0
# xavier_share_gsa                    161                       319
# 
#                      hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa
# bernstein_gsa                        0               0                0
# farkkila_gsa                         0               0                0
# franchimont_gsa                      0               0                0
# helmsley_prism_gsa                   0               0                0
# hyams_protect_gsa                    0               0                0
# lewis_sparc_gsa                      0               0                0
# mccauley_new_gsa                     0               1                0
# mcgovern_gsa                         0               0                0
# newberry_share_gsa                   0               9                0
# niddk_cho_gsa                        0               4                0
# niddk_duerr_gsa                      1             105                0
# niddk_rioux_gsa                      0               0                0
# niddk_silverberg_gsa                 0               0                0
# pekow_share_gsa                      0              23                0
# rioux_igenomed_gsa                   0               0                0
# sands_msccr_gsa                      0               0                2
# vermeire_gsa                         0               0                0
# xavier_prism_gsa                     0               0                0
# xavier_share_gsa                     0               6                0
# 
#                      mcgovern_gsa moayyedi_imagine_gsa niddk_duerr_gsa
# bernstein_gsa                   0                   20               0
# farkkila_gsa                    0                    0               0
# franchimont_gsa                 0                    0               0
# helmsley_prism_gsa              0                    0               0
# hyams_protect_gsa               1                    1               0
# lewis_sparc_gsa                 5                    0               0
# mccauley_new_gsa               12                    0               0
# mcgovern_gsa                    0                    4               0
# newberry_share_gsa              3                    0               0
# niddk_cho_gsa                   7                    0               2
# niddk_duerr_gsa                 2                    0               0
# niddk_rioux_gsa                 0                   13               0
# niddk_silverberg_gsa            0                    7               0
# pekow_share_gsa                 5                    0               0
# rioux_igenomed_gsa              0                    6               0
# sands_msccr_gsa                10                    0               0
# vermeire_gsa                    1                    0               0
# xavier_prism_gsa                0                    0               0
# xavier_share_gsa                6                    0               0
# 
#                      niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa
# bernstein_gsa                      0                    2               0
# farkkila_gsa                       0                    0               1
# franchimont_gsa                    0                    0               0
# helmsley_prism_gsa                 0                    0               0
# hyams_protect_gsa                  0                    0               0
# lewis_sparc_gsa                    0                    0               0
# mccauley_new_gsa                   0                    0               0
# mcgovern_gsa                       1                    2               0
# newberry_share_gsa                 0                    0               0
# niddk_cho_gsa                      1                    2               0
# niddk_duerr_gsa                    0                    7               0
# niddk_rioux_gsa                    0                    1               0
# niddk_silverberg_gsa               0                    0               0
# pekow_share_gsa                    0                    0               0
# rioux_igenomed_gsa                 8                    0               0
# sands_msccr_gsa                    0                    0               0
# vermeire_gsa                       0                    0               0
# xavier_prism_gsa                   0                    0               0
# xavier_share_gsa                   0                    0               0
# 
#                      pekow_share_gsa sands_msccr_gsa stampfer_gsa
# bernstein_gsa                      0               0            0
# farkkila_gsa                       0               0            0
# franchimont_gsa                    0               0            0
# helmsley_prism_gsa                 0               0            0
# hyams_protect_gsa                  0               0            0
# lewis_sparc_gsa                    0               0            0
# mccauley_new_gsa                   0               0            0
# mcgovern_gsa                       0               0            0
# newberry_share_gsa                 0               1            0
# niddk_cho_gsa                     24             251            1
# niddk_duerr_gsa                    0               0            1
# niddk_rioux_gsa                    0               0            0
# niddk_silverberg_gsa               0               0            0
# pekow_share_gsa                    0               1            0
# rioux_igenomed_gsa                 0               0            0
# sands_msccr_gsa                    0               1            0
# vermeire_gsa                       0               0            0
# xavier_prism_gsa                   0               0            0
# xavier_share_gsa                   0               0            0
# 
#                      vermeire_gsa xavier_prism_gsa xavier_share_gsa
# bernstein_gsa                   0                0                0
# farkkila_gsa                    0                0                0
# franchimont_gsa                 9                0                0
# helmsley_prism_gsa              0                0                0
# hyams_protect_gsa               0                0                0
# lewis_sparc_gsa                 0                8                0
# mccauley_new_gsa                0                1                1
# mcgovern_gsa                    0                0                0
# newberry_share_gsa              0                0                0
# niddk_cho_gsa                   1                2                3
# niddk_duerr_gsa                 0                1                0
# niddk_rioux_gsa                 0                0                0
# niddk_silverberg_gsa            0                0                0
# pekow_share_gsa                 0                0                1
# rioux_igenomed_gsa              0                0                0
# sands_msccr_gsa                 0                0                1
# vermeire_gsa                    0                0                0
# xavier_prism_gsa                0                0                0
# xavier_share_gsa                0               43                0



#################################################################
### MATCHES OLD-NEW TO KEEP

dim(dup)
# [1] 20961    11


xx<-dup[,c("cohort_FID1","cohort_FID2")]
yy<-dup[,c("cohort_FID2","cohort_FID1")]
colnames(yy)<-colnames(xx)

dim(xx)
# [1] 20961    2

xx<-rbind(xx,yy)

x<-table(xx$cohort_FID1,xx$cohort_FID2)

write.table(x,paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_old_new_broad_batch_table_studies_2.txt",sep=""),
            col.names=T,row.names=T,quote=F,sep="\t")


#################################################################

table(dup$cohort_type_FID1,dup$cohort_type_FID2)
#       new   old
# new  1201     0
# old 19759     1

# COHORTS TO BE COMPLETELY REPLACED BY NEW ONES

# belgium_franchimont_gsa by franchimont_gsa
# belgium_vermeire_gsa by vermeire_gsa
# ccfa_gsa by lewis_sparc_gsa
# cedars_gsa by mcgovern_gsa
# mccauley_gsa by mccauley_new_gsa
# prism_nfe_gsa by xavier_prism_gsa
# prism_nfe_gwas by xavier_share_gsa

# exclude any pairs with those matches - 
dup_tmp<-dup

dup_tmp<-dup[which(! (dup$cohort_FID1=="belgium_franchimont_gsa" & 
                        dup$cohort_FID2=="franchimont_gsa")),]
dup_tmp<-dup_tmp[which(! (dup_tmp$cohort_FID1=="belgium_vermeire_gsa" & 
                            dup_tmp$cohort_FID2=="vermeire_gsa")),]
dup_tmp<-dup_tmp[which(! (dup_tmp$cohort_FID1=="ccfa_gsa" & 
                            dup_tmp$cohort_FID2=="lewis_sparc_gsa")),]
dup_tmp<-dup_tmp[which(! (dup_tmp$cohort_FID1=="cedars_gsa" & 
                            dup_tmp$cohort_FID2=="mcgovern_gsa")),]
dup_tmp<-dup_tmp[which(! (dup_tmp$cohort_FID1=="mccauley_gsa" & 
                            dup_tmp$cohort_FID2=="mccauley_new_gsa")),]
dup_tmp<-dup_tmp[which(! (dup_tmp$cohort_FID1=="prism_nfe_gsa" & 
                            dup_tmp$cohort_FID2=="xavier_prism_gsa")),]
dup_tmp<-dup_tmp[which(! (dup_tmp$cohort_FID1=="prism_nfe_gwas" & 
                            dup_tmp$cohort_FID2=="xavier_share_gsa")),]

dim(dup_tmp)
# [1] 10439    11

table(dup_tmp$cohort_type_FID1,dup_tmp$cohort_type_FID2)
#     new   old
# new 1201    0
# old 9237    1

# prioritize any cohort over sands - gender issue not resolved yet 

# new batches over old
new_gwas<-c("bernstein_gsa","farkkila_gsa","franchimont_gsa","franke_gsa",
            "hyams_protect_gsa","lewis_sparc_gsa","mccauley_new_gsa","mcgovern_gsa",
            "moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa","niddk_duerr_gsa",
            "niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
            "rioux_igenomed_gsa","stampfer_gsa",
            "vermeire_gsa","xavier_prism_gsa","xavier_share_gsa","netherlands_gsa")

# exclude "sands_msccr_gsa", from this list and include in bad array below to set lowest priority
# exclude "helmsley_prism_gsa","helmsley_xavier_prism_gsa", include below in bad arrays, non gsa
# exclude "weersma_gsa" and add "netherlands_gsa" - all should be contained in netherlands

old_gwas<-cohorts[!cohorts %in% new_gwas]

# classify old gwas into good and bad, prioritize good:

bad_array<-c("stampfer_gsa","helmsley_xavier_prism_gsa","helmsley_prism_gsa",
             "prism_nfe_gwas","finland_illugwas")
good_array<-old_gwas[!old_gwas %in% bad_array]


cohorts[which(!cohorts %in% c(old_gwas,new_gwas))] 
# 0

old_gwas[which(!old_gwas %in% c(good_array,bad_array))] 
# 0



table(dup$sex_FID1,dup$sex_FID2)
#      0    1    2
# 0    34     6    27
# 1   187 10023     0
# 2   184     0 10500


dup<-dup[,c("FID1","FID2","Kinship","cohort_FID1","pheno_FID1","F_MISS_FID1","cohort_FID2","pheno_FID2","F_MISS_FID2")]

table(dup$pheno_FID1,dup$pheno_FID2)
#      1    2
# 1  4168    76
# 2  1403 15314

dup_ids<-c(as.character(dup$FID1),as.character(dup$FID2))
length(dup_ids)
#[1] 41922
table(length(dup_ids)==(nrow(dup)*2))
# TRUE 
# 1 

dup_ids<-dup_ids[!duplicated(dup_ids)]
length(dup_ids)
#[1] 39573

dim(dup)
# [1] 20961    9

rm(data_remove,data_inconsist)

for (i in 1:length(dup_ids)) {
  
  tmp1<-dup[which(dup$FID1==dup_ids[i]),]
  tmp2<-dup[which(dup$FID2==dup_ids[i]),]
  colnames(tmp2)[1:2]<-colnames(tmp2)[2:1]
  
  tmp<-rbind(tmp1,tmp2)
  
  ids_tmp<-c(as.character(tmp$FID1),as.character(tmp$FID2))
  ids_tmp<-ids_tmp[!duplicated(ids_tmp)]
  
  #keep same order as in dup_ids and in all
  ids_tmp<-ids_tmp[match(dup_ids[which(dup_ids %in% ids_tmp)],ids_tmp)]
  
  # get number of possible combinations
  n_possible_combinations<-nrow(permutations(length(ids_tmp), 2, letters[1:length(ids_tmp)]))/2
  
  
  # DOUBLE CHEK THAT WE GET ALL EXPECTED PAIRS
  
  if (nrow(dup[which((dup$FID1 %in% ids_tmp) | (dup$FID2 %in% ids_tmp)),])==n_possible_combinations ) {
    
    data<-as.data.frame(matrix(ncol=1,nrow=length(ids_tmp)))
    data$V1<-ids_tmp
    data<-merge(data,dat[,c("V1","sex","pheno","cohort","F_MISS")],by="V1",all.x=T,sort=F)
    
    
    # OPTION A: samples have different pheno, remove all (already checked before all samples have same sex:
    
    if ( (dim(table(data$pheno))!=1) ) {
      
      # # make an exeption to deal with sex = 0
      # if ( any(data$sex==0) & any(data$sex %in% c(1,2)) & dim(table(data$pheno))==1 ) {
      #   
      #   keep_sample<-data$V1[which(data$sex!=0)]
      #   
      #   if(!exists("data_remove")) {
      #     data_remove<-data[which(!data$V1 %in% keep_sample),]
      #   } else {
      #     data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
      #   }
      #   
      # } else 
      if (any(data$cohort %in% c("sands_msccr_gsa"))){
        
        data<-data[which(data$cohort=="sands_msccr_gsa"),]
        
        if (!exists("data_remove")) {
          data_remove<-data
        } else {
          data_remove<-rbind(data_remove,data)
        }
        
      } else {
        
        if (!exists("data_remove")) {
          data_remove<-data
        } else {
          data_remove<-rbind(data_remove,data)
        }
        
        # keep this also as independent table:
        if(!exists("data_inconsist")){
          jj<-1
          data_inconsist<-data
          data_inconsist$group<-jj
        } else {
          jj<-jj+1
          data$group<-jj
          data_inconsist<-rbind(data_inconsist,data)
        }
      }
      
    } 
    
    else {
      
      # OPTION B: samples have same pheno:
      
      # B.1: rest of samples, if samples in new and old arrays, keep only samples in new arrays:
      
      if (any(data$cohort %in% old_gwas) & any(data$cohort %in% new_gwas)) {
        
        keep_sample<-data$V1[which(!data$cohort %in% old_gwas)]
        
        # if more than one new pair in new gwas, keep sample with largest call rate, make exception for CCFA and McCauley:
        if (length(keep_sample)>1) {
          
          if ( all(data$cohort %in% c("ccfa_gsa","cedars_gsa")) | all(data$cohort %in% c("ccfa_gsa","mccauley_gsa"))) {
            
            keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]# to deal with >1 with same F_miss
            
          } 
          
          else if ( any(data$cohort %in% c("ccfa_gsa","mccauley_gsa","cedars_gsa")) ) {
            
            keep_sample<-data[which(!(data$cohort %in% old_gwas) & !(data$cohort %in% c("ccfa_gsa","mccauley_gsa","cedars_gsa"))),]
            keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
            
          } 
          else {
            
            keep_sample<-data[which(!data$cohort %in% old_gwas),]
            keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
            
          }
          
        }
        
        if (!exists("data_remove")) {
          data_remove<-data[which(!data$V1 %in% keep_sample),]
        } else {
          data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
        }
        
      } 
      
      # B.2: old_gwas duplicated pairs only
      
      else if (all(data$cohort %in% old_gwas)) {
        
        # B.2.1: if in good and in bad array, keep sample in the best array:
        
        if (any(data$cohort %in% good_array) & any(data$cohort %in% bad_array)) {
          
          keep_sample<-data$V1[which(data$cohort %in% good_array)]
          
          # if more than one new pair in good array, keep sample with largest call rate:
          if (length(keep_sample)>1) {
            keep_sample<-data[which(data$cohort %in% good_array),]
            keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
          }
          
          if(!exists("data_remove")) {
            data_remove<-data[which(!data$V1 %in% keep_sample),]
          } else {
            data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
          }
          
        } 
        
        # B.2.2: if all in good OR all in bad array, keep sample in the best array:
        
        else if (all(data$cohort %in% good_array) | all(data$cohort %in% bad_array)) {
          
          keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]# to deal with same F_miss
          
          if(!exists("data_remove")) {
            data_remove<-data[which(!data$V1 %in% keep_sample),]
          } else {
            data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
          }
          
        }
      }
      
      # B.3: new_gwas duplicated pairs only, keep sample with largest call rate: 
      
      else if (all(data$cohort %in% new_gwas)) {
        
        # B.3.1: Samples in ccfa (no ctr) and other new cohort, keep other:
        
        if ( all(data$cohort %in% c("ccfa_gsa","cedars_gsa")) | all(data$cohort %in% c("ccfa_gsa","mccauley_gsa"))) {
          
          keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]# to deal with >1 with same F_miss
          
          if(!exists("data_remove")) {
            data_remove<-data[which(!data$V1 %in% keep_sample),]
          } else {
            data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
          }
          
        } else if ( any(data$cohort %in% c("ccfa_gsa","mccauley_gsa","cedars_gsa")) ) {
          
          keep_sample<-data$V1[which(!data$cohort %in% c("ccfa_gsa","mccauley_gsa","cedars_gsa"))]
          
          # if more than one new pair in new gwas, keep sample with largest call rate:
          if (length(keep_sample)>1) {
            
            keep_sample<-data[which(!data$cohort %in% c("ccfa_gsa","mccauley_gsa","cedars_gsa")),]
            keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
            
          }
          
          if(!exists("data_remove")) {
            data_remove<-data[which(!data$V1 %in% keep_sample),]
          } else {
            data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
          }
        } 
        
        # B.3.2: Samples in kiel and gwas3 OR belgium:
        
        else if ( all( any(data$cohort=="all_hce" | data$cohort=="belgium_vermeire_gsa") & any(data$cohort %in% c("kiel_austria_sibdcs_gsa")) )  ) {
          
          keep_sample<-data$V1[which(!data$cohort %in% c("kiel_austria_sibdcs_gsa"))]
          
          # if more than one new pair in new gwas, keep sample with largest call rate:
          if (length(keep_sample)>1) {
            keep_sample<-data[which(!data$cohort %in% c("kiel_austria_sibdcs_gsa")),]
            keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
          }
          
          if(!exists("data_remove")) {
            data_remove<-data[which(!data$V1 %in% keep_sample),]
          } else {
            data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
          }
        } 
        
        # B.3.3: Other keep sample with largest call rate:
        
        else {
          
          keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))]
          
          if(!exists("data_remove")) {
            data_remove<-data[which(!data$V1 %in% keep_sample),]
          } else {
            data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
          }
        }
      } 
      
    }
  }
  
  else {
    
    # not all combinations of duplicated pairs found, show list of IDs, to manually inspect issues
    print(i)
    print(paste("Number of expected combinations: ",n_possible_combinations,sep=""))
    print(paste("Number of observed combinations: ",nrow(dup[which((dup$FID1 %in% ids_tmp) | (dup$FID2 %in% ids_tmp)),]),sep=""))
    print(ids_tmp)
    
  }
  
}


dim(data_remove)
#[1] 44761  4

data_remove<-data_remove[!duplicated(data_remove$V1),]
dim(data_remove)
# [1] 21562    5
# double check which sample per pair removed in pairs:

table(data_remove$cohort)
# belgium_franchimont_gsa         belgium_louis_gsa      belgium_vermeire_gsa 
#                    1486                         4                      3860 
# bernstein_gsa                  ccfa_gsa                cedars_gsa 
#            17                      2063                      1786 
# finland_illugwas           franchimont_gsa                franke_gsa 
#               19                         5                         3 
# helmsley_prism_gsa helmsley_xavier_prism_gsa   kiel_austria_sibdcs_gsa 
#                760                       357                        93 
# lewis_sparc_gsa              mccauley_gsa          mccauley_new_gsa 
#              77                       764                        10 
# mcgovern_gsa      moayyedi_imagine_gsa           netherlands_gsa 
#           62                        21                         1 
# newberry_share_gsa           niddk_broad_gsa             niddk_cho_gsa 
#                  1                      4937                       298 
# niddk_duerr_gsa       niddk_feinstein_gsa           niddk_rioux_gsa 
#             497                      2049                       273 
# niddk_silverberg_gsa           palotie_hus_gsa           pekow_share_gsa 
#                  378                         1                        33 
# prism_nfe_gsa            prism_nfe_gwas        rioux_igenomed_gsa 
# 422                       166                         9 
# sands_msccr_gsa              stampfer_gsa              vermeire_gsa 
# 342                         1                         5 
# weersma_gsa          xavier_prism_gsa          xavier_share_gsa 
# 696                        24                        42

# vs ones in each study:
# basque_gsa   belgium_franchimont_gsa         belgium_louis_gsa 
# 1491                      1488                      1506 
# belgium_vermeire_gsa             bernstein_gsa                  ccfa_gsa 
# 3979                       481                      2064 
# cedars_gsa              farkkila_gsa          finland_illugwas 
# 1803                        68                       444 
# franchimont_gsa                franke_gsa        helmsley_prism_gsa 
# 2727                       860                       760 
# helmsley_xavier_prism_gsa         hyams_protect_gsa                 italy_gsa 
# 1293                       418                       949 
# kiel_austria_sibdcs_gsa           lewis_sparc_gsa             lithuania_gsa 
# 13870                      2857                      2210 
# mccauley_gsa          mccauley_new_gsa              mcgovern_gsa 
# 769                      1612                      5960 
# moayyedi_imagine_gsa           netherlands_gsa        newberry_share_gsa 
# 1128                      4555                       861 
# niddk_broad_gsa             niddk_cho_gsa           niddk_duerr_gsa 
# 5122                      1746                      1817 
# niddk_feinstein_gsa           niddk_rioux_gsa      niddk_silverberg_gsa 
# 7960                       913                      2330 
# palotie_hus_gsa           pekow_share_gsa             prism_nfe_gsa 
# 878                       633                       464 
# prism_nfe_gwas        rioux_igenomed_gsa           sands_msccr_gsa 
# 735                       176                      1422 
# slovenia_gsa              stampfer_gsa                sweden_gsa 
# 262                      1468                      1379 
# vermeire_gsa               weersma_gsa          xavier_prism_gsa 
# 4683                       709                       687 
# xavier_share_gsa 
# 693 

#### create a file to report inconsistent data:

dim(data_inconsist)
# [1] 5922 6
data_inconsist<-data_inconsist[!duplicated(data_inconsist$V1),]
dim(data_inconsist)
# [1] 2864  6

# explore these are report back to the suppliers

write.table(data_remove,paste(path,"pre_imputation/QC/relatedness/list_iibdgc_old_new_broad_batch_data_remove",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
write.table(data_inconsist,paste(path,"pre_imputation/QC/relatedness/list_iibdgc_old_new_broad_batch_data_inconsist",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

# data_remove<-read.table(paste(path,"pre_imputation/QC/relatedness/list_iibdgc_old_new_broad_batch_data_remove",sep=""),head=T)
# data_inconsist<-read.table(paste(path,"pre_imputation/QC/relatedness/list_iibdgc_old_new_broad_batch_data_inconsist",sep=""),head=T)

# if we exclude these, how the final datasets look like? 
dat$cohort<-as.factor(dat$cohort)
table(dat$cohort)
# basque_gsa   belgium_franchimont_gsa         belgium_louis_gsa 
#       1491                      1488                      1506 
# belgium_vermeire_gsa             bernstein_gsa                  ccfa_gsa 
#                 3979                       481                      2064 
# cedars_gsa              farkkila_gsa          finland_illugwas 
#       1803                        68                       444 
# franchimont_gsa                franke_gsa        helmsley_prism_gsa 
#            2727                       860                       760 
# helmsley_xavier_prism_gsa         hyams_protect_gsa                 italy_gsa 
#                      1293                       418                       949 
# kiel_austria_sibdcs_gsa           lewis_sparc_gsa             lithuania_gsa 
#                   13870                      2857                      2210 
# mccauley_gsa          mccauley_new_gsa              mcgovern_gsa 
#          769                      1612                      5960 
# moayyedi_imagine_gsa           netherlands_gsa        newberry_share_gsa 
#                 1128                      4555                       861 
# niddk_broad_gsa             niddk_cho_gsa           niddk_duerr_gsa 
#            5122                      1746                      1817 
# niddk_feinstein_gsa           niddk_rioux_gsa      niddk_silverberg_gsa 
#                7960                       913                      2330 
# palotie_hus_gsa           pekow_share_gsa             prism_nfe_gsa 
#             878                       633                       464 
# prism_nfe_gwas        rioux_igenomed_gsa           sands_msccr_gsa 
#            735                       176                      1422 
# slovenia_gsa              stampfer_gsa                sweden_gsa 
#          262                      1468                      1379 
# vermeire_gsa               weersma_gsa          xavier_prism_gsa 
#         4683                       709                       687 
# xavier_share_gsa 
#              693 

table(dat$cohort[which(!dat$V2 %in% data_remove$V1)])
# basque_gsa   belgium_franchimont_gsa         belgium_louis_gsa 
#       1491                         2                      1502 
# belgium_vermeire_gsa             bernstein_gsa                  ccfa_gsa 
#                  119                       464                         1 
# cedars_gsa              farkkila_gsa          finland_illugwas 
#         17                        68                       425 
# franchimont_gsa                franke_gsa        helmsley_prism_gsa 
#            2722                       857                         0 
# helmsley_xavier_prism_gsa         hyams_protect_gsa                 italy_gsa 
#                       936                       418                       949 
# kiel_austria_sibdcs_gsa           lewis_sparc_gsa             lithuania_gsa 
#                   13777                      2780                      2210 
# mccauley_gsa          mccauley_new_gsa              mcgovern_gsa 
#            5                      1602                      5898 
# moayyedi_imagine_gsa           netherlands_gsa        newberry_share_gsa 
#                 1107                      4554                       860 
# niddk_broad_gsa             niddk_cho_gsa           niddk_duerr_gsa 
#             185                      1448                      1320 
# niddk_feinstein_gsa           niddk_rioux_gsa      niddk_silverberg_gsa 
# 5911                       640                      1952 
# palotie_hus_gsa           pekow_share_gsa             prism_nfe_gsa 
# 877                       600                        42 
# prism_nfe_gwas        rioux_igenomed_gsa           sands_msccr_gsa 
# 569                       167                      1080 
# slovenia_gsa              stampfer_gsa                sweden_gsa 
# 262                      1467                      1379 
# vermeire_gsa               weersma_gsa          xavier_prism_gsa 
# 4678                        13                       663 
# xavier_share_gsa 
# 651


### save files:

data_remove<-data_remove[,c(1,1)]
colnames(data_remove)<-c("FID","IID")
dim(data_remove)
# [1] 21562     2
write.table(data_remove,paste(path,"pre_imputation/QC/relatedness/list_iibdgc_old_new_broad_batch_data_remove",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

colnames(data_inconsist)[1]<-c("IIBDGC_id")
write.table(data_inconsist[,c("IIBDGC_id","sex","pheno","cohort","group")],paste(path,"pre_imputation/QC/relatedness/list_iibdgc_old_new_broad_batch_data_inconsistent_phenotype.tsv",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")

##################################

# work out where the inconsitencies come from:

#### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"
library(data.table)

data_inconsist<-fread(paste(path,"pre_imputation/QC/relatedness/list_iibdgc_old_new_broad_batch_data_inconsistent_phenotype.tsv",sep=""),head=T,sep="\t")
data_inconsist<-as.data.frame(data_inconsist)

table(data_inconsist$cohort)[order(table(data_inconsist$cohort),decreasing=T)]
# niddk_broad_gsa           niddk_duerr_gsa      niddk_silverberg_gsa 
# 1282                       441                       371 
# niddk_cho_gsa           niddk_rioux_gsa       niddk_feinstein_gsa 
# 283                       266                       105 
# mcgovern_gsa           lewis_sparc_gsa           pekow_share_gsa 
# 38                        27                        13 
# moayyedi_imagine_gsa        rioux_igenomed_gsa                franke_gsa 
# 9                         4                         3 
# kiel_austria_sibdcs_gsa          xavier_share_gsa                cedars_gsa 
# 3                         3                         2 
# helmsley_prism_gsa helmsley_xavier_prism_gsa             prism_nfe_gsa 
# 2                         2                         2 
# prism_nfe_gwas          xavier_prism_gsa                  ccfa_gsa 
# 2                         2                         1 
# mccauley_new_gsa           netherlands_gsa               weersma_gsa 
# 1                         1                         1

cohorts<-names(table(data_inconsist$cohort))

for (i in 1:length(cohorts)) {
  
  print(cohorts[i])
  
  ids<-data_inconsist$group[which(data_inconsist$cohort %in% cohorts[i])]
  print(table(data_inconsist$cohort[which(data_inconsist$group %in% ids)]))
  
}

##################
# [1] "ccfa_gsa"
# 
# ccfa_gsa  lewis_sparc_gsa mccauley_new_gsa 
# 1                1                1 
# [1] "cedars_gsa"
# 
# cedars_gsa mcgovern_gsa 
# 2            2 
# [1] "franke_gsa"
# 
# franke_gsa kiel_austria_sibdcs_gsa 
# 3                       3 
# [1] "helmsley_prism_gsa"
# 
# helmsley_prism_gsa    niddk_broad_gsa      niddk_cho_gsa     prism_nfe_gwas 
# 2                  1                  1                  1 
# xavier_share_gsa 
# 1 
# [1] "helmsley_xavier_prism_gsa"
# 
# helmsley_xavier_prism_gsa           niddk_broad_gsa             niddk_cho_gsa 
# 2                         2                         1 
# niddk_silverberg_gsa 
# 1 
# [1] "kiel_austria_sibdcs_gsa"
# 
# franke_gsa kiel_austria_sibdcs_gsa 
# 3                       3 
# [1] "lewis_sparc_gsa"
# 
# ccfa_gsa     lewis_sparc_gsa    mccauley_new_gsa     niddk_broad_gsa 
# 1                  27                   1                  23 
# niddk_cho_gsa     niddk_duerr_gsa niddk_feinstein_gsa 
# 1                  25                   3 
# [1] "mccauley_new_gsa"
# 
# ccfa_gsa  lewis_sparc_gsa mccauley_new_gsa 
# 1                1                1 
# [1] "mcgovern_gsa"
# 
# cedars_gsa         mcgovern_gsa      niddk_broad_gsa 
# 2                   38                   28 
# niddk_cho_gsa  niddk_feinstein_gsa niddk_silverberg_gsa 
# 1                    7                    1 
# prism_nfe_gwas     xavier_share_gsa 
# 1                    1 
# [1] "moayyedi_imagine_gsa"
# 
# moayyedi_imagine_gsa      niddk_broad_gsa  niddk_feinstein_gsa 
# 9                    7                    2 
# niddk_rioux_gsa niddk_silverberg_gsa   rioux_igenomed_gsa 
# 7                    2                    1 
# [1] "netherlands_gsa"
# 
# netherlands_gsa     weersma_gsa 
# 1               1 
# [1] "niddk_broad_gsa"
# 
# helmsley_prism_gsa helmsley_xavier_prism_gsa           lewis_sparc_gsa 
# 1                         2                        23 
# mcgovern_gsa      moayyedi_imagine_gsa           niddk_broad_gsa 
# 28                         7                      1282 
# niddk_cho_gsa           niddk_duerr_gsa           niddk_rioux_gsa 
# 280                       369                       246 
# niddk_silverberg_gsa           pekow_share_gsa        rioux_igenomed_gsa 
# 364                        12                         4 
# [1] "niddk_cho_gsa"
# 
# helmsley_prism_gsa helmsley_xavier_prism_gsa           lewis_sparc_gsa 
# 1                         1                         1 
# mcgovern_gsa           niddk_broad_gsa             niddk_cho_gsa 
# 1                       280                       283 
# niddk_duerr_gsa       niddk_feinstein_gsa           pekow_share_gsa 
# 1                         2                        13 
# [1] "niddk_duerr_gsa"
# 
# lewis_sparc_gsa      niddk_broad_gsa        niddk_cho_gsa 
# 25                  369                    1 
# niddk_duerr_gsa  niddk_feinstein_gsa niddk_silverberg_gsa 
# 441                   71                    4 
# [1] "niddk_feinstein_gsa"
# 
# lewis_sparc_gsa         mcgovern_gsa moayyedi_imagine_gsa 
# 3                    7                    2 
# niddk_cho_gsa      niddk_duerr_gsa  niddk_feinstein_gsa 
# 2                   71                  105 
# niddk_rioux_gsa niddk_silverberg_gsa 
# 20                    6 
# [1] "niddk_rioux_gsa"
# 
# moayyedi_imagine_gsa      niddk_broad_gsa  niddk_feinstein_gsa 
# 7                  246                   20 
# niddk_rioux_gsa   rioux_igenomed_gsa 
# 266                    4 
# [1] "niddk_silverberg_gsa"
# 
# helmsley_xavier_prism_gsa              mcgovern_gsa      moayyedi_imagine_gsa 
# 1                         1                         2 
# niddk_broad_gsa           niddk_duerr_gsa       niddk_feinstein_gsa 
# 364                         4                         6 
# niddk_silverberg_gsa 
# 371 
# [1] "pekow_share_gsa"
# 
# niddk_broad_gsa   niddk_cho_gsa pekow_share_gsa 
# 12              13              13 
# [1] "prism_nfe_gsa"
# 
# prism_nfe_gsa xavier_prism_gsa xavier_share_gsa 
# 2                2                1 
# [1] "prism_nfe_gwas"
# 
# helmsley_prism_gsa       mcgovern_gsa     prism_nfe_gwas   xavier_share_gsa 
# 1                  1                  2                  2 
# [1] "rioux_igenomed_gsa"
# 
# moayyedi_imagine_gsa      niddk_broad_gsa      niddk_rioux_gsa 
# 1                    4                    4 
# rioux_igenomed_gsa 
# 4 
# [1] "weersma_gsa"
# 
# netherlands_gsa     weersma_gsa 
# 1               1 
# [1] "xavier_prism_gsa"
# 
# prism_nfe_gsa xavier_prism_gsa xavier_share_gsa 
# 2                2                1 
# [1] "xavier_share_gsa"
# 
# helmsley_prism_gsa       mcgovern_gsa      prism_nfe_gsa     prism_nfe_gwas 
# 1                  1                  1                  2 
# xavier_prism_gsa   xavier_share_gsa 
# 1                  3 

##################

# combine with the per study inconsistencies:

studies<-c("niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
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

for (j in 1:length(studies)) {
  
  print(studies[j])
  
  file<-paste(path,"pre_imputation/QC/",studies[j],"/_list_duplicated_samples_inconsistent_phenotype",sep="")
  if(file.exists(file)){
    
    tmp<-read.table(file,head=T)
    tmp$cohort<-studies[j]

    
    if(!exists("data_inconsist_2")) {
      tmp$group<-max(data_inconsist$group)+(tmp$group)
      data_inconsist_2<-tmp
    }else{
      tmp$group<-max(data_inconsist_2$group)+(tmp$group)
      data_inconsist_2<-rbind(data_inconsist_2,tmp)
    }
  }
}

colnames(data_inconsist_2)[1]<-"IIBDGC_id"
data_inconsist_2<-data_inconsist_2[,colnames(data_inconsist)]

data_inconsist<-rbind(data_inconsist,data_inconsist_2)

write.table(data_inconsist[,c("IIBDGC_id","sex","pheno","cohort","group")],paste(path,"pre_imputation/QC/relatedness/list_iibdgc_old_new_broad_batch_data_inconsistent_phenotype_2.tsv",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")

table(data_inconsist$cohort)[order(table(data_inconsist$cohort),decreasing=T)]

# niddk_broad_gsa           niddk_duerr_gsa      niddk_silverberg_gsa 
# 1286                       453                       379 
# niddk_cho_gsa           niddk_rioux_gsa       niddk_feinstein_gsa 
# 289                       268                       110 
# mcgovern_gsa   kiel_austria_sibdcs_gsa           lewis_sparc_gsa 
# 80                        73                        27 
# pekow_share_gsa      moayyedi_imagine_gsa           netherlands_gsa 
# 13                         9                         5 
# rioux_igenomed_gsa                franke_gsa          xavier_share_gsa 
# 4                         3                         3 
# all_hce                cedars_gsa        helmsley_prism_gsa 
# 2                         2                         2 
# helmsley_xavier_prism_gsa             prism_nfe_gsa            prism_nfe_gwas 
# 2                         2                         2 
# sweden_gsa          xavier_prism_gsa                  ccfa_gsa 
# 2                         2                         1 
# mccauley_new_gsa               weersma_gsa 
# 1                         1 

