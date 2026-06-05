# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
##############################################################################
# 10.- REMOVE VARIANTS WITH SIGNIF DISCREPANCY IN CALL RATE BETWEEN CA VS CTR #
##############################################################################

### /software/R-4.3.1/bin/R

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
           "swedish_uc_old_gwas",
           "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

for (j in 1:length(studies)) {
  
  print(studies[j])
  
  fam<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01.fam",sep=""),head=F)
  print(table(fam$V6,useNA="ifany"))

  write.table(fam[which(fam$V6==1),1:2],paste(path,"pre_imputation/QC/",studies[j],"/list_controls",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  write.table(fam[which(fam$V6==2),1:2],paste(path,"pre_imputation/QC/",studies[j],"/list_cases",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
}

#############

# [1] "all_hce"
# 1     2 
# 10455 13435 

# [1] "niddk_old_gwas"
# 1    2 
# 935 1825

# [1] "australia_omniexome"
# 1   2 
# 615 691

# [1] "gwas1"
# 1    2 
# 2933 1747 

# [1] "gwas2"
# 1    2 
# 5417 2361 

# [1] "pittsburgh_gsa"
# 1    2 
# 1468 1257 

# [1] "spain_gsa"
# 1    2 
# 1482 1961 

# [1] "italy_gsa"
# 1   2 
# 375 578 

# [1] "kiel_austria_sibdcs_gsa"
# 1    2 
# 4584 9819 

# [1] "netherlands_gsa"
# 1    2 
# 550 4031 

# [1] "slovenia_gsa"
# 1   2 
# 177  87 

# [1] "sweden_gsa"
# 1   2 
# 986 415 

# [1] "niddk_broad_gsa"
# 1    2 
# 1017 4438 

# [1] "niddk_feinstein_gsa"
# 1    2 
# 2424 5753 

# [1] "basque_gsa"
# 1   2 
# 983 525 

# [1] "lithuania_gsa"
# 1    2 
# 1143 1088 

# [1] "belgium_louis_gsa"
# 1   2 
# 600 918 

# [1] "belgium_franchimont_gsa"
# 
# 1   2 
# 595 906 
# [1] "belgium_vermeire_gsa"
# 
# 1    2 
# 814 3179 
# [1] "prism_nfe_gsa"
# 
# 1   2 
# 33 433 
# [1] "prism_nfe_gwas"
# 
# 1   2 
# 262 556 

# [1] "finland_illugwas"
# 2 
# 446 


# [1] "german_affy6_old_gwas"
# 1    2 
# 1768 1055 
# [1] "norway_affy6_old_gwas"
# 
# 1   2 
# 282 268 
# [1] "belgium_inf1_old_gwas"
# 
# 1   2 
# 900 517 
# [1] "belgium_inf2_old_gwas"
# 
# 1   2 
# 111 161 

# [1] "cedars_370k_old_gwas"
# 2 
# 605 

# [1] "cedars_610k_old_gwas"
# 2 
# 889 


# [1] "cedars_omni_old_gwas"
# 
# 1    2 
# 4 1217 
# [1] "swedish_uc_old_gwas"
# 
# 1   2 
# 341 923 

# [1] "mccauley_gsa"
# 2 
# 782 


# [1] "ccfa_gsa"
# 2 
# 2177 


# [1] "cedars_gsa"
# 1    2 
# 947 2131 

# [1] "bernstein_gsa"
# 1   2 
# 10 501 

# [1] "farkkila_gsa"
# 2 
# 68 


# [1] "franchimont_gsa"
# 
# 1    2 
# 1489 1299 
# [1] "franke_gsa"
# 
# 1   2 
# 432 435 
# [1] "helmsley_prism_gsa"
# 
# 1   2 
# 260 501 
# [1] "helmsley_xavier_prism_gsa"
# 
# 1    2 
# 160 1133

# [1] "hyams_protect_gsa"
# 2 
# 418 


# [1] "lewis_sparc_gsa"
# 2 
# 2859 


# [1] "mccauley_new_gsa"
# 
# 1    2 
# 238 1381 

# [1] "mcgovern_gsa"
# 1    2 
# 1205 4817 
# [1] "moayyedi_imagine_gsa"
# 
# 1   2 
# 251 877 

# [1] "newberry_share_gsa"
# 2 
# 863 

# [1] "niddk_cho_gsa"
# 
# 1    2 
# 740 1024 
# [1] "niddk_duerr_gsa"
# 
# 1    2 
# 918 1026 
# [1] "niddk_rioux_gsa"
# 
# 1   2 
# 583 336 
# [1] "niddk_silverberg_gsa"
# 1    2 
# 544 1827 
#

# [1] "palotie_hus_gsa"
# 2 
# 878 
# 
# [1] "pekow_share_gsa"
# 2 
# 634 
# 
# [1] "rioux_igenomed_gsa"
# 2 
# 182 
# 
# [1] "sands_msccr_gsa"
# 1    2 
# 309 1121 
# 
# [1] "stampfer_gsa"
# 1    2 
# 1032  436 
# 
# [1] "vermeire_gsa"
# 1    2 
# 814 3883 
# 
# [1] "weersma_gsa"
# 1   2 
# 385 324 
# 
# [1] "xavier_prism_gsa"
# 1   2 
# 62 627 
# 
# [1] "xavier_share_gsa"
# 2 
# 696 

#############

path_gwas=/path/to/ibdgwas/IIBDGC/
studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

MEM=500

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_29_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_29_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01 \
--allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/${i}/list_controls \
--missing \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_ctr"
done

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_30_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_30_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01 \
--allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/${i}/list_cases \
--missing \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_cases"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_29_${i} | grep -E "completed"
done

###########

# all_hce
# Successfully completed.
# niddk_old_gwas
# Successfully completed.
# australia_omniexome
# Successfully completed.
# gwas1
# Successfully completed.
# gwas2
# Successfully completed.
# pittsburgh_gsa
# Successfully completed.
# spain_gsa
# Successfully completed.
# italy_gsa
# Successfully completed.
# kiel_austria_sibdcs_gsa
# Successfully completed.
# netherlands_gsa
# Successfully completed.
# slovenia_gsa
# Successfully completed.
# sweden_gsa
# Successfully completed.
# niddk_broad_gsa
# Successfully completed.
# niddk_feinstein_gsa
# Successfully completed.
# basque_gsa
# Successfully completed.
# lithuania_gsa
# Successfully completed.
# belgium_louis_gsa
# Successfully completed.
# belgium_franchimont_gsa
# Successfully completed.
# belgium_vermeire_gsa
# Successfully completed.
# prism_nfe_gsa
# Successfully completed.
# prism_nfe_gwas
# Successfully completed.
# finland_illugwas   ### ONLY CASES
# german_affy6_old_gwas
# Successfully completed.
# norway_affy6_old_gwas
# Successfully completed.
# belgium_inf1_old_gwas
# Successfully completed.
# belgium_inf2_old_gwas
# Successfully completed.
# cedars_370k_old_gwas ### ONLY CASES
# cedars_610k_old_gwas ### ONLY CASES
# cedars_omni_old_gwas
# Successfully completed.
# swedish_uc_old_gwas
# Successfully completed.
# mccauley_gsa ### ONLY CASES
# ccfa_gsa ### ONLY CASES
# cedars_gsa 
# Successfully completed.
# bernstein_gsa
# Successfully completed.
# farkkila_gsa # ONLY CASES
# franchimont_gsa
# Successfully completed.
# franke_gsa
# Successfully completed.
# helmsley_prism_gsa
# Successfully completed.
# helmsley_xavier_prism_gsa
# Successfully completed.
# hyams_protect_gsa ### ONLY CASES
# lewis_sparc_gsa ### ONLY CASES
# mccauley_new_gsa
# Successfully completed.
# mcgovern_gsa
# Successfully completed.
# moayyedi_imagine_gsa
# Successfully completed.
# newberry_share_gsa ### ONLY CASES
# niddk_cho_gsa
# Successfully completed.
# niddk_duerr_gsa
# Successfully completed.
# niddk_rioux_gsa
# Successfully completed.
# niddk_silverberg_gsa
# Successfully completed.
# palotie_hus_gsa ### ONLY CASES
# pekow_share_gsa ### ONLY CASES
# rioux_igenomed_gsa ### ONLY CASES
# sands_msccr_gsa
# Successfully completed.
# stampfer_gsa
# Successfully completed.
# vermeire_gsa
# Successfully completed.
# weersma_gsa
# Successfully completed.
# xavier_prism_gsa
# Successfully completed.
# xavier_share_gsa ### ONLY CASES

###########


for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_30_${i} | grep -E "completed"
done

###########
# all_hce
# Successfully completed.
# niddk_old_gwas
# Successfully completed.
# australia_omniexome
# Successfully completed.
# gwas1
# Successfully completed.
# gwas2
# Successfully completed.
# pittsburgh_gsa
# Successfully completed.
# spain_gsa
# Successfully completed.
# italy_gsa
# Successfully completed.
# kiel_austria_sibdcs_gsa
# Successfully completed.
# netherlands_gsa
# Successfully completed.
# slovenia_gsa
# Successfully completed.
# sweden_gsa
# Successfully completed.
# niddk_broad_gsa
# Successfully completed.
# niddk_feinstein_gsa
# Successfully completed.
# basque_gsa
# Successfully completed.
# lithuania_gsa
# Successfully completed.
# belgium_louis_gsa
# Successfully completed.
# belgium_franchimont_gsa
# Successfully completed.
# belgium_vermeire_gsa
# Successfully completed.
# prism_nfe_gsa
# Successfully completed.
# prism_nfe_gwas
# Successfully completed.
# finland_illugwas
# Successfully completed.
# german_affy6_old_gwas
# Successfully completed.
# norway_affy6_old_gwas
# Successfully completed.
# belgium_inf1_old_gwas
# Successfully completed.
# belgium_inf2_old_gwas
# Successfully completed.
# cedars_370k_old_gwas
# Successfully completed.
# cedars_610k_old_gwas
# Successfully completed.
# cedars_omni_old_gwas
# Successfully completed.
# swedish_uc_old_gwas
# Successfully completed.
# mccauley_gsa
# Successfully completed.
# ccfa_gsa
# Successfully completed.
# cedars_gsa
# Successfully completed.
# bernstein_gsa
# Successfully completed.
# farkkila_gsa
# Successfully completed.
# franchimont_gsa
# Successfully completed.
# franke_gsa
# Successfully completed.
# helmsley_prism_gsa
# Successfully completed.
# helmsley_xavier_prism_gsa
# Successfully completed.
# hyams_protect_gsa
# Successfully completed.
# lewis_sparc_gsa
# Successfully completed.
# mccauley_new_gsa
# Successfully completed.
# mcgovern_gsa
# Successfully completed.
# moayyedi_imagine_gsa
# Successfully completed.
# newberry_share_gsa
# Successfully completed.
# niddk_cho_gsa
# Successfully completed.
# niddk_duerr_gsa
# Successfully completed.
# niddk_rioux_gsa
# Successfully completed.
# niddk_silverberg_gsa
# Successfully completed.
# palotie_hus_gsa
# Successfully completed.
# pekow_share_gsa
# Successfully completed.
# rioux_igenomed_gsa
# Successfully completed.
# sands_msccr_gsa
# Successfully completed.
# stampfer_gsa
# Successfully completed.
# vermeire_gsa
# Successfully completed.
# weersma_gsa
# Successfully completed.
# xavier_prism_gsa
# Successfully completed.
# xavier_share_gsa
# Successfully completed.

###########


###### /software/R-4.3.1/bin/R

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
           "swedish_uc_old_gwas",
           "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")


foo <- function(y){
  # include here as.numeric to be sure that your values are numeric:
  table <-  matrix(as.numeric(c(y[2], y[3], y[4], y[5])), ncol = 2, byrow = TRUE)
  if(any(is.na(table))) p <- "error" else p <- fisher.test(table, alternative="two.sided")$p.value
  p
}

for (j in 1:length(studies)) {
  
  # print(studies[j])
  
  # skip cohorts with just cases:
  file_ctr<-paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_ctr.lmiss",sep="")
  
  if(file.exists(file_ctr)) {
    
    ca_var_miss<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_cases.lmiss",sep=""),head=T)
    ctr_var_miss<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_ctr.lmiss",sep=""),head=T)

    ca_var_miss$N_NO_MISS<-ca_var_miss$N_GENO-ca_var_miss$N_MISS
    ctr_var_miss$N_NO_MISS<-ctr_var_miss$N_GENO-ctr_var_miss$N_MISS

    colnames(ca_var_miss)[c(3,6)]<-paste(colnames(ca_var_miss)[c(3,6)],"_cases",sep="")
    colnames(ctr_var_miss)[c(3,6)]<-paste(colnames(ctr_var_miss)[c(3,6)],"_ctr",sep="")

    all<-merge(ca_var_miss[,c(2,3,6)],ctr_var_miss[,c(2,3,6)],by="SNP",sort=F)

    all$fisher_pvalue <- apply(all, 1, foo)

    print(nrow(all[which(all$fisher_pvalue<1E-4),]))

    # pdf(paste(path,"pre_imputation/QC/",studies[j],"/hist_missingness_pvalue_ca_ctr.pdf",sep=""),height = 5,width = 10)
    # ggplot(all, aes(x=(-log10(fisher_pvalue)))) + geom_histogram(binwidth=1) + geom_vline(xintercept=(-log10(1E-4)), linetype="dashed",color = "red", size=1)
    # dev.off()

    write.table(all,paste(path,"pre_imputation/QC/",studies[j],"/comparison_missingness_cases_ctr.lmiss",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    write.table(all[which(all$fisher_pvalue<1E-4),1,drop=F],paste(path,"pre_imputation/QC/",studies[j],"/list_variants_exclude_by_missingness_cases_ctr.lmiss",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

  } else {
    print(studies[j])
    all<-as.data.frame(matrix(ncol=1,nrow=0))
    colnames(all)<-"SNP"
    write.table(all,paste(path,"pre_imputation/QC/",studies[j],"/list_variants_exclude_by_missingness_cases_ctr.lmiss",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  }
  
  rm(list=ls()[!ls() %in% c("j","studies","path","foo")])
}

##############
# some changes in Ca/Ctr Numbers resulted in slightly different N variants to be excluded in this 2022 round vs 2020 for the cohorts processed in 2020
# [1] "all_hce"
# [1] 63865
# [1] "niddk_old_gwas"
# [1] 7532  
# [1] "australia_omniexome"
# [1] 669
# [1] "gwas1"
# [1] 1284
# [1] "gwas2"
# [1] 5955
# [1] "pittsburgh_gsa"
# [1] 25826
# [1] "spain_gsa"
# [1] 4340
# [1] "italy_gsa"
# [1] 338
# [1] "kiel_austria_sibdcs_gsa"
# [1] 9222
# [1] "netherlands_gsa"
# [1] 2971
# [1] "slovenia_gsa"
# [1] 1
# [1] "sweden_gsa"
# [1] 4505
# [1] "niddk_broad_gsa"
# [1] 185
# [1] "niddk_feinstein_gsa"
# [1] 32714
# [1] "basque_gsa"
# [1] 3948
# [1] "lithuania_gsa"
# [1] 7
# [1] "belgium_louis_gsa"
# [1] 0
# [1] "belgium_franchimont_gsa"
# [1] 2
# [1] "belgium_vermeire_gsa"
# [1] 50
# [1] "prism_nfe_gsa"
# [1] 1
# [1] "prism_nfe_gwas"
# [1] 50
# [1] "finland_illugwas"
# [1] "german_affy6_old_gwas"
# [1] 31269
# [1] "norway_affy6_old_gwas"
# [1] 11598
# [1] "belgium_inf1_old_gwas"
# [1] 976
# [1] "belgium_inf2_old_gwas"
# [1] 532
# [1] "cedars_370k_old_gwas"
# [1] "cedars_610k_old_gwas"
# [1] "cedars_omni_old_gwas"
# [1] 0
# [1] "swedish_uc_old_gwas"
# [1] 349
# [1] "mccauley_gsa"
# [1] "ccfa_gsa"
# [1] "cedars_gsa"
# [1] 11
# [1] "bernstein_gsa"
# [1] 7
# [1] "farkkila_gsa"
# [1] "franchimont_gsa"
# [1] 100
# [1] "franke_gsa"
# [1] 2
# [1] "helmsley_prism_gsa"
# [1] 9
# [1] "helmsley_xavier_prism_gsa"
# [1] 1
# [1] "hyams_protect_gsa"
# [1] "lewis_sparc_gsa"
# [1] "mccauley_new_gsa"
# [1] 321
# [1] "mcgovern_gsa"
# [1] 494
# [1] "moayyedi_imagine_gsa"
# [1] 1
# [1] "newberry_share_gsa"
# [1] "niddk_cho_gsa"
# [1] 618
# [1] "niddk_duerr_gsa"
# [1] 1801
# [1] "niddk_rioux_gsa"
# [1] 643
# [1] "niddk_silverberg_gsa"
# [1] 435
# [1] "palotie_hus_gsa"
# [1] "pekow_share_gsa"
# [1] "rioux_igenomed_gsa"
# [1] "sands_msccr_gsa"
# [1] 38
# [1] "stampfer_gsa"
# [1] 3
# [1] "vermeire_gsa"
# [1] 46
# [1] "weersma_gsa"
# [1] 0
# [1] "xavier_prism_gsa"
# [1] 1
# [1] "xavier_share_gsa"

##############

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_31_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_31_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01 \
--allow-no-sex \
--make-bed \
--exclude ${path_gwas}pre_imputation/QC/${i}/list_variants_exclude_by_missingness_cases_ctr.lmiss \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_31_${i} | grep -E "completed"
done

###########

# all_hce
# Successfully completed.
# niddk_old_gwas
# Successfully completed.
# australia_omniexome
# Successfully completed.
# gwas1
# Successfully completed.
# gwas2
# Successfully completed.
# pittsburgh_gsa
# Successfully completed.
# spain_gsa
# Successfully completed.
# italy_gsa
# Successfully completed.
# kiel_austria_sibdcs_gsa
# Successfully completed.
# netherlands_gsa
# Successfully completed.
# slovenia_gsa
# Successfully completed.
# sweden_gsa
# Successfully completed.
# niddk_broad_gsa
# Successfully completed.
# niddk_feinstein_gsa
# Successfully completed.
# basque_gsa
# Successfully completed.
# lithuania_gsa
# Successfully completed.
# belgium_louis_gsa
# Successfully completed.
# belgium_franchimont_gsa
# Successfully completed.
# belgium_vermeire_gsa
# Successfully completed.
# prism_nfe_gsa
# Successfully completed.
# prism_nfe_gwas
# Successfully completed.
# finland_illugwas
# Successfully completed.
# german_affy6_old_gwas
# Successfully completed.
# norway_affy6_old_gwas
# Successfully completed.
# belgium_inf1_old_gwas
# Successfully completed.
# belgium_inf2_old_gwas
# Successfully completed.
# cedars_370k_old_gwas
# Successfully completed.
# cedars_610k_old_gwas
# Successfully completed.
# cedars_omni_old_gwas
# Successfully completed.
# swedish_uc_old_gwas
# Successfully completed.
# mccauley_gsa
# Successfully completed.
# ccfa_gsa
# Successfully completed.
# cedars_gsa
# Successfully completed.
# bernstein_gsa
# Successfully completed.
# farkkila_gsa
# Successfully completed.
# franchimont_gsa
# Successfully completed.
# franke_gsa
# Successfully completed.
# helmsley_prism_gsa
# Successfully completed.
# helmsley_xavier_prism_gsa
# Successfully completed.
# hyams_protect_gsa
# Successfully completed.
# lewis_sparc_gsa
# Successfully completed.
# mccauley_new_gsa
# Successfully completed.
# mcgovern_gsa
# Successfully completed.
# moayyedi_imagine_gsa
# Successfully completed.
# newberry_share_gsa
# Successfully completed.
# niddk_cho_gsa
# Successfully completed.
# niddk_duerr_gsa
# Successfully completed.
# niddk_rioux_gsa
# Successfully completed.
# niddk_silverberg_gsa
# Successfully completed.
# palotie_hus_gsa
# Successfully completed.
# pekow_share_gsa
# Successfully completed.
# rioux_igenomed_gsa
# Successfully completed.
# sands_msccr_gsa
# Successfully completed.
# stampfer_gsa
# Successfully completed.
# vermeire_gsa
# Successfully completed.
# weersma_gsa
# Successfully completed.
# xavier_prism_gsa
# Successfully completed.
# xavier_share_gsa
# Successfully completed.

###########

for i in ${studies[@]}
do echo ${i} \
&& wc -l ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam \
&& wc -l ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim
done

###########

# all_hce
# 23890 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/all_hce/all_hce_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 436292 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/all_hce/all_hce_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# niddk_old_gwas
# 2760 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 300366 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# australia_omniexome
# 1306 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/australia_omniexome_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 740026 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/australia_omniexome_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# gwas1
# 4680 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/gwas1_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 454780 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/gwas1_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# gwas2
# 7778 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas2/gwas2_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 783694 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas2/gwas2_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# pittsburgh_gsa
# 2725 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 898529 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# spain_gsa
# 3443 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/spain_gsa/spain_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 577252 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/spain_gsa/spain_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# italy_gsa
# 953 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/italy_gsa/italy_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 574072 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/italy_gsa/italy_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# kiel_austria_sibdcs_gsa
# 14403 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 566099 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# netherlands_gsa
# 4581 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 667396 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# slovenia_gsa
# 264 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 529731 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# sweden_gsa
# 1401 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/sweden_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 559371 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/sweden_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# niddk_broad_gsa
# 5455 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 591636 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# niddk_feinstein_gsa
# 8177 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 604463 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# basque_gsa
# 1508 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/basque_gsa/basque_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 567984 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/basque_gsa/basque_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# lithuania_gsa
# 2231 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lithuania_gsa/lithuania_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 556492 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lithuania_gsa/lithuania_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# belgium_louis_gsa
# 1518 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 592896 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# belgium_franchimont_gsa
# 1501 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 563237 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# belgium_vermeire_gsa
# 3993 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 655357 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# prism_nfe_gsa
# 466 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 559195 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# prism_nfe_gwas
# 818 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 243031 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# finland_illugwas
# 446 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/finland_illugwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 237095 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/finland_illugwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# german_affy6_old_gwas
# 2823 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 790342 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# norway_affy6_old_gwas
# 550 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 731823 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# belgium_inf1_old_gwas
# 1417 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 301714 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# belgium_inf2_old_gwas
# 272 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 290191 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# cedars_370k_old_gwas
# 605 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 338736 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim
# 
# cedars_610k_old_gwas
# 889 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 584005 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim
# 
# cedars_omni_old_gwas
# 1221 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 719502 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# swedish_uc_old_gwas
# 1264 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 296655 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim
# 
# mccauley_gsa
# 782 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 557000 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# ccfa_gsa
# 2177 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 601251 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# cedars_gsa
# 3078 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_gsa/cedars_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 647811 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_gsa/cedars_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# bernstein_gsa
# 511 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/bernstein_gsa/bernstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 554198 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/bernstein_gsa/bernstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# farkkila_gsa
# 68 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 511520 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# franchimont_gsa
# 2788 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franchimont_gsa/franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 650404 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franchimont_gsa/franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# franke_gsa
# 867 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franke_gsa/franke_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 558258 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franke_gsa/franke_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# helmsley_prism_gsa
# 761 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 243233 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# helmsley_xavier_prism_gsa
# 1293 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 243225 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# hyams_protect_gsa
# 418 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 569715 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# lewis_sparc_gsa
# 2859 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 609245 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim
# 
# mccauley_new_gsa
# 1619 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 579265 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# mcgovern_gsa
# 6022 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 633320 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# moayyedi_imagine_gsa
# 1128 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 608670 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# newberry_share_gsa
# 863 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 584248 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# niddk_cho_gsa
# 1764 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 594546 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# niddk_duerr_gsa
# 1944 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 582404 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# niddk_rioux_gsa
# 919 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 575873 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# niddk_silverberg_gsa
# 2371 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 612238 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# palotie_hus_gsa
# 878 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 554381 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# pekow_share_gsa
# 634 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 574453 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# rioux_igenomed_gsa
# 182 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 543416 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# sands_msccr_gsa
# 1430 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 595986 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# stampfer_gsa
# 1468 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 645334 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# vermeire_gsa
# 4697 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/vermeire_gsa/vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 655936 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/vermeire_gsa/vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# weersma_gsa
# 709 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/weersma_gsa/weersma_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 567070 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/weersma_gsa/weersma_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# xavier_prism_gsa
# 689 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 570022 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

# xavier_share_gsa
# 696 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam
# 568855 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim

###########

############################################################################################################################################


