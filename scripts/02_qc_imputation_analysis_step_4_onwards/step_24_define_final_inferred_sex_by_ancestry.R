# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
#######################################
# Define inferred (final) inferred sex

#######################################
# 24.1 CREATE LIST OF MALES AND FEMALES - UPDATED WITH SAMPLES THAT PASS QC AND BY NEW PHENOTYPE DATA

##############################
### /software/R-4.3.1/bin/R

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
           "swedish_uc_old_gwas",
           "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

studies_new<-c("bernstein_gsa","farkkila_gsa",
               "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
               "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
               "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
               "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
               "xavier_share_gsa","franchimont_gsa")



################################################################################################
# RE-do case control status for all new studies (batches from 2021 and 2022)

pheno<-fread(paste(path,"pheno/gwas-mega-2-core-phenotypes-28ceede.csv.gz",sep=""),head=T)
pheno<-as.data.frame(pheno)

colnames(pheno)[6]<-"gender"

pheno$sex<-0
pheno$sex[which(pheno$gender=="Female")]<-2
pheno$sex[which(pheno$gender=="Male")]<-1
table(pheno$gender,pheno$sex,useNA="ifany")
#             0     1     2
#          3261     0     0
# Female      0     0 55427
# Male        0 51830     0
# Unknown  8872     0     0


pheno$pheno<-NA
pheno$pheno[which(pheno$control=="1")]<-"1"
pheno$pheno[which(pheno$affection=="Affected")]<-"2"
table(pheno$pheno,pheno$affection,useNA="ifany")
#            Affected Unaffected Unknown
# 1        0        0      28970       0
# 2        0    85835          0       0
# <NA>  3967        0        143     475
table(pheno$pheno,pheno$diag,useNA="ifany")
#            Crohn's Disease Indeterminate Ulcerative Colitis Unknown
# 1    28970               0             0                  0       0
# 2       17           49394          1221              31608    3595
# <NA>  4345               0             0                  0     240

################################################################
# collect all post-QC fam files, before being split into of ancestry

for (j in 1:length(studies)) {
  
  print(studies[j])
  fam_tmp<-fread(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam",sep=""),head=F)
  fam_tmp<-as.data.frame(fam_tmp)
  fam_tmp$study<-studies[j]
  
  # for studies pre-2021, keep data in fam fike
  if (!studies[j] %in% studies_new) {
    fam_tmp$sex<-fam_tmp$V5
    fam_tmp$pheno<-fam_tmp$V6
  } else {
    fam_tmp<-merge(fam_tmp,pheno[,c("sample_id","sex","pheno")],by.x="V1",by.y="sample_id",order=F,all.x=T)
    print(table(fam_tmp$V5,fam_tmp$sex,useNA="ifany"))
    print(table(fam_tmp$V6,fam_tmp$pheno,useNA="ifany"))
  }
  
  if(j==1) {
    fam_all<-fam_tmp
  }else{
    fam_all<-rbind(fam_all,fam_tmp)
  }
  
  rm(fam_tmp)
  
}


for (j in 1:length(studies)) {
  
  print(studies[j])
  print(table(fam_all$sex[which(fam_all$study==studies[j])],useNA="ifany"))
  
  fam<-fam_all[which(fam_all$study==studies[j]),]
  
  fam.male<-fam[which(fam$sex==1),1:2] 
  colnames(fam.male)<-c("FID","IID")
  
  fam.female<-fam[which(fam$sex==2),1:2] 
  colnames(fam.female)<-c("FID","IID")
  
  if(nrow(fam.male)>0) {
    write.table(fam.male,paste(path,"pre_imputation/QC/",studies[j],"/list_male_samples_step24",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  }
  if(nrow(fam.female)>0) {
    write.table(fam.female,paste(path,"pre_imputation/QC/",studies[j],"/list_female_samples_step24",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  }
  
  rm(list=ls()[!ls() %in% c("j","studies","path","fam_all")])
}

q("no")

##############################
# [1] "all_hce"
# 0     1     2 
# 30 11022 12711 
# [1] "niddk_old_gwas"
# 
# 1    2 
# 1374 1358 
# [1] "australia_omniexome"
# 
# 0   1   2 
# 1 602 696 
# [1] "gwas1"
# 
# 0    1    2 
# 2 2112 2545 
# [1] "gwas2"
# 
# 0    1    2 
# 78 3832 3863 
# [1] "pittsburgh_gsa"
# 
# 0    1    2 
# 7 1179 1523 
# [1] "spain_gsa"
# 
# 0    1    2 
# 56 1861 1491 
# [1] "italy_gsa"
# 
# 0   1   2 
# 3 510 430 
# [1] "kiel_austria_sibdcs_gsa"
# 
# 0    1    2 
# 213 6992 7124 
# [1] "netherlands_gsa"
# 
# 0    1    2 
# 10 1883 2644 
# [1] "slovenia_gsa"
# 
# 0   1   2 
# 5  55 203 
# [1] "sweden_gsa"
# 
# 1   2 
# 812 584 
# [1] "niddk_broad_gsa"
# 
# 0    1    2 
# 18 2813 2595 
# [1] "niddk_feinstein_gsa"
# 
# 0    1    2 
# 11 4074 4039 
# [1] "basque_gsa"
# 
# 0   1   2 
# 6 919 573 
# [1] "lithuania_gsa"
# 
# 1    2 
# 1177 1048 
# [1] "belgium_louis_gsa"
# 
# 0   1   2 
# 8 656 847 
# [1] "belgium_franchimont_gsa"
# 
# 0   1   2 
# 1 640 850 
# [1] "belgium_vermeire_gsa"
# 
# 0    1    2 
# 14 1907 2053 
# [1] "prism_nfe_gsa"
# 
# 0   1   2 
# 4 213 245 
# [1] "prism_nfe_gwas"
# 
# 1   2 
# 392 418 
# [1] "finland_illugwas"
# 
# 1   2 
# 296 148 
# [1] "german_affy6_old_gwas"
# 
# 0    1    2 
# 3 1355 1439 
# [1] "norway_affy6_old_gwas"
# 
# 0   1   2 
# 1 301 242 
# [1] "belgium_inf1_old_gwas"
# 
# 0   1   2 
# 11 982 413 
# [1] "belgium_inf2_old_gwas"
# 
# 0   1   2 
# 1 136 134 
# [1] "cedars_370k_old_gwas"
# 
# 0   1   2 
# 1 320 276 
# [1] "cedars_610k_old_gwas"
# 
# 0   1   2 
# 2 460 419 
# [1] "cedars_omni_old_gwas"
# 
# 1   2 
# 590 617 
# [1] "swedish_uc_old_gwas"
# 
# 1   2 
# 693 565 
# [1] "mccauley_gsa"
# 
# 1   2 
# 399 377 
# [1] "ccfa_gsa"
# 
# 0    1    2 
# 1  993 1173 
# [1] "cedars_gsa"
# 
# 0    1    2 
# 4 1477 1578 
# [1] "bernstein_gsa"
# 
# 0   1   2 
# 217 120 172 
# [1] "farkkila_gsa"
# 
# 1  2 
# 32 36 
# [1] "franchimont_gsa"
# 
# 0    1    2 
# 71 1264 1434 
# [1] "franke_gsa"
# 
# 1   2 
# 497 367 
# [1] "helmsley_prism_gsa"
# 
# 0   1   2 
# 2 368 383 
# [1] "helmsley_xavier_prism_gsa"
# 
# 0   1   2 
# 22 568 691 
# [1] "hyams_protect_gsa"
# 
# 1   2 
# 212 204 
# [1] "lewis_sparc_gsa"
# 
# 1    2 
# 1280 1566 
# [1] "mccauley_new_gsa"
# 1   2
# 831 779 
# [1] "mcgovern_gsa"
# 0    1    2 
# 18 2930 3030 
# [1] "moayyedi_imagine_gsa"
# 1   2 
# 488 632 
# [1] "newberry_share_gsa"
# 1   2 
# 364 494 
# [1] "niddk_cho_gsa"
# 0   1   2 
# 4 901 849 
# [1] "niddk_duerr_gsa"
# 0   1   2 
# 3 967 967 
# [1] "niddk_rioux_gsa"
# 0   1   2 
# 2 417 494 
# [1] "niddk_silverberg_gsa"
# 0    1    2 
# 5 1210 1132 
# [1] "palotie_hus_gsa"
# 1   2 
# 416 453 
# [1] "pekow_share_gsa"
# 1   2 
# 316 316 
# [1] "rioux_igenomed_gsa"
# 1  2 
# 87 94 
# [1] "sands_msccr_gsa"
# 1   2 
# 759 660 
# [1] "stampfer_gsa"
# 1    2 
# 236 1228 
# [1] "vermeire_gsa"
# 1    2 
# 2258 2401 
# [1] "weersma_gsa"
# 1   2 
# 355 350 
# [1] "xavier_prism_gsa"
# 
# 1   2 
# 340 346 
# [1] "xavier_share_gsa"
# 1   2 
# 346 346 


##############################

###################################
# 7.2 KEEP ONLY FEMALES AND CHR23

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_24.1_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_24.1_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip \
--keep ${path_gwas}pre_imputation/QC/${i}/list_female_samples_step24 \
--split-x b37 'no-fail' --chr 23 --make-bed \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_females_only_step24"
done


for i in ${studies[@]}
do echo ${i} && tail -100 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_24.1_${i} | grep -E "completed"
done

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_24.2_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_24.2_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_females_only_step24 \
--missing --hardy --freq --out ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_females_only_step24"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_24.2_${i} | grep -E "completed"
done

for i in ${studies[@]}
do 
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_females_only_step24.hwe
done

# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_check_sex_females_only.hwe': No such file or directory


##############################
### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"
# exclude swedish_uc_old_gwas - no chrX or chrY data
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

length(studies)
# [1] 57

for (j in 1:length(studies)) {
  
  print(studies[j])
  
  bim<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_check_sex_females_only_step24.bim",sep=""),sep="\t",head=F)
  print(paste("N markers in chrX:",table(bim$V1)))
  
  hwe<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_check_sex_females_only_step24.hwe",sep=""),head=T)
  hwe<-hwe[which(hwe$TEST=="UNAFF"),]
  
  frq<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_check_sex_females_only_step24.frq",sep=""),head=T)
  var_miss<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_check_sex_females_only_step24.lmiss",sep=""),head=T)
  
  all<-merge(hwe[,c("SNP","P")],frq[,c("SNP","MAF")],by="SNP",sort=F)
  all<-merge(all,var_miss[,c("SNP","F_MISS")],by="SNP",sort=F)
  
  all<-all[which(all$MAF>0.1 & all$F_MISS<0.001 & all$P>1E-4),]
  print(paste("N good chrX markers:",nrow(all)))
  
  write.table(all,paste(path,"pre_imputation/QC/",studies[j],"/list_good_chrX_variants_step24",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
  rm(list=ls()[!ls() %in% c("j","studies","path")])
  
}

q("no")

##################################

# [1] "all_hce"
# [1] "N markers in chrX: 12477"
# [1] "N good chrX markers: 5002"
# [1] "niddk_old_gwas"
# [1] "N markers in chrX: 8968"
# [1] "N good chrX markers: 2454"
# [1] "australia_omniexome"
# [1] "N markers in chrX: 18586"
# [1] "N good chrX markers: 8501"
# [1] "gwas1"
# [1] "N markers in chrX: 9746"
# [1] "N good chrX markers: 1946"
# [1] "gwas2"
# [1] "N markers in chrX: 36382"
# [1] "N good chrX markers: 11557"
# [1] "pittsburgh_gsa"
# [1] "N markers in chrX: 22362"
# [1] "N good chrX markers: 8376"
# [1] "spain_gsa"
# [1] "N markers in chrX: 14735"
# [1] "N good chrX markers: 8916"
# [1] "italy_gsa"
# [1] "N markers in chrX: 24945"
# [1] "N good chrX markers: 10203"
# [1] "kiel_austria_sibdcs_gsa"
# [1] "N markers in chrX: 26756"
# [1] "N good chrX markers: 5282"
# [1] "netherlands_gsa"
# [1] "N markers in chrX: 17412"
# [1] "N good chrX markers: 4076"
# [1] "slovenia_gsa"
# [1] "N markers in chrX: 13858"
# [1] "N good chrX markers: 6404"
# [1] "sweden_gsa"
# [1] "N markers in chrX: 14315"
# [1] "N good chrX markers: 5846"
# [1] "niddk_broad_gsa"
# [1] "N markers in chrX: 15302"
# [1] "N good chrX markers: 3583"
# [1] "niddk_feinstein_gsa"
# [1] "N markers in chrX: 15451"
# [1] "N good chrX markers: 6902"
# [1] "basque_gsa"
# [1] "N markers in chrX: 24702"
# [1] "N good chrX markers: 9271"
# [1] "lithuania_gsa"
# [1] "N markers in chrX: 14423"
# [1] "N good chrX markers: 7632"
# [1] "belgium_louis_gsa"
# [1] "N markers in chrX: 15664"
# [1] "N good chrX markers: 6541"
# [1] "belgium_franchimont_gsa"
# [1] "N markers in chrX: 14798"
# [1] "N good chrX markers: 6383"
# [1] "belgium_vermeire_gsa"
# [1] "N markers in chrX: 16934"
# [1] "N good chrX markers: 7360"
# [1] "prism_nfe_gsa"
# [1] "N markers in chrX: 14522"
# [1] "N good chrX markers: 7857"
# [1] "prism_nfe_gwas"
# [1] "N markers in chrX: 6371"
# [1] "N good chrX markers: 4297"
# [1] "finland_illugwas"
# [1] "N markers in chrX: 6056"
# [1] "N good chrX markers: 4892"
# [1] "german_affy6_old_gwas"
# [1] "N markers in chrX: 34415"
# [1] "N good chrX markers: 8327"
# [1] "norway_affy6_old_gwas"
# [1] "N markers in chrX: 31902"
# [1] "N good chrX markers: 9852"
# [1] "belgium_inf1_old_gwas"
# [1] "N markers in chrX: 7568"
# [1] "N good chrX markers: 4285"
# [1] "belgium_inf2_old_gwas"
# [1] "N markers in chrX: 7922"
# [1] "N good chrX markers: 4100"
# [1] "cedars_370k_old_gwas"
# [1] "N markers in chrX: 9924"
# [1] "N good chrX markers: 6034"
# [1] "cedars_610k_old_gwas"
# [1] "N markers in chrX: 14735"
# [1] "N good chrX markers: 8169"
# [1] "cedars_omni_old_gwas"
# [1] "N markers in chrX: 17670"
# [1] "N good chrX markers: 10157"
# [1] "mccauley_gsa"
# [1] "N markers in chrX: 14417"
# [1] "N good chrX markers: 7534"
# [1] "ccfa_gsa"
# [1] "N markers in chrX: 15164"
# [1] "N good chrX markers: 8048"
# [1] "cedars_gsa"
# [1] "N markers in chrX: 16799"
# [1] "N good chrX markers: 6986"
# [1] "bernstein_gsa"
# [1] "N markers in chrX: 14244"
# [1] "N good chrX markers: 8126"
# [1] "farkkila_gsa"
# [1] "N markers in chrX: 13238"
# [1] "N good chrX markers: 8080"
# [1] "franchimont_gsa"
# [1] "N markers in chrX: 17028"
# [1] "N good chrX markers: 6900"
# [1] "franke_gsa"
# [1] "N markers in chrX: 14212"
# [1] "N good chrX markers: 7734"
# [1] "helmsley_prism_gsa"
# [1] "N markers in chrX: 6371"
# [1] "N good chrX markers: 4551"
# [1] "helmsley_xavier_prism_gsa"
# [1] "N markers in chrX: 6371"
# [1] "N good chrX markers: 4696"
# [1] "hyams_protect_gsa"
# [1] "N markers in chrX: 14714"
# [1] "N good chrX markers: 8450"
# [1] "lewis_sparc_gsa"
# [1] "N markers in chrX: 15381"
# [1] "N good chrX markers: 7976"
# [1] "mccauley_new_gsa"
# [1] "N markers in chrX: 14945"
# [1] "N good chrX markers: 7140"
# [1] "mcgovern_gsa"
# [1] "N markers in chrX: 15958"
# [1] "N good chrX markers: 7370"
# [1] "moayyedi_imagine_gsa"
# [1] "N markers in chrX: 14969"
# [1] "N good chrX markers: 7269"
# [1] "newberry_share_gsa"
# [1] "N markers in chrX: 15057"
# [1] "N good chrX markers: 7623"
# [1] "niddk_cho_gsa"
# [1] "N markers in chrX: 15106"
# [1] "N good chrX markers: 6703"
# [1] "niddk_duerr_gsa"
# [1] "N markers in chrX: 14975"
# [1] "N good chrX markers: 6605"
# [1] "niddk_rioux_gsa"
# [1] "N markers in chrX: 14769"
# [1] "N good chrX markers: 7237"
# [1] "niddk_silverberg_gsa"
# [1] "N markers in chrX: 15313"
# [1] "N good chrX markers: 7862"
# [1] "palotie_hus_gsa"
# [1] "N markers in chrX: 14252"
# [1] "N good chrX markers: 7137"
# [1] "pekow_share_gsa"
# [1] "N markers in chrX: 14678"
# [1] "N good chrX markers: 8018"
# [1] "rioux_igenomed_gsa"
# [1] "N markers in chrX: 13960"
# [1] "N good chrX markers: 9064"
# [1] "sands_msccr_gsa"
# [1] "N markers in chrX: 14937"
# [1] "N good chrX markers: 6834"
# [1] "stampfer_gsa"
# [1] "N markers in chrX: 17138"
# [1] "N good chrX markers: 7057"
# [1] "vermeire_gsa"
# [1] "N markers in chrX: 16946"
# [1] "N good chrX markers: 7308"
# [1] "weersma_gsa"
# [1] "N markers in chrX: 14441"
# [1] "N good chrX markers: 7373"
# [1] "xavier_prism_gsa"
# [1] "N markers in chrX: 14668"
# [1] "N good chrX markers: 7596"
# [1] "xavier_share_gsa"
# [1] "N markers in chrX: 14720"
# [1] "N good chrX markers: 7574"

##################################

##################################
# 7.3 KEEP ONLY MALES AND CHR24

# use data from stage before excluding chrY

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_24.3_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_24.3_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip \
--keep ${path_gwas}pre_imputation/QC/${i}/list_male_samples_step24 --chr 24 --make-bed \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_males_only_step24"
done

for i in ${studies[@]}
do 
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_males_only_step24.bim
done

# OK - no chrY
# niddk_old_gwas
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_check_sex_males_only_step24.bim': No such file or directory# gwas1
# gwas1
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/gwas1_check_sex_males_only_step24.bim': No such file or directory
# prism_nfe_gwas
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_check_sex_males_only_step24.bim': No such file or directory
# finland_illugwas
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/finland_illugwas_check_sex_males_only_step24.bim': No such file or directory
# swedish_uc_old_gwas
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_check_sex_males_only_step24.bim': No such file or directory
# belgium_inf1_old_gwas
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_check_sex_males_only_step24.bim': No such file or directory
# belgium_inf2_old_gwas
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_check_sex_males_only_step24.bim': No such file or directory
# helmsley_prism_gsa
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_check_sex_males_only_step24.bim': No such file or directory
# helmsley_xavier_prism_gsa
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_check_sex_males_only_step24.bim': No such file or directory

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_24.4_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_24.4_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_males_only_step24 \
--freq --missing \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_males_only_step24"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_24.4_${i} | grep -E "completed"
done

for i in ${studies[@]}
do 
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_males_only_step24.frq
done


#  AND FIND VARIANTS WITH NO CALLS IN MOST OF FEMALES

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_24.5_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_24.5_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip \
--keep ${path_gwas}pre_imputation/QC/${i}/list_female_samples_step24 --chr 24 --recode A \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_females_only_chr24_step24"
done

for i in ${studies[@]}
do 
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_females_only_chr24_step24.raw
done

##############################
### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

# exclude swedish_uc_old_gwas - no chrX or chrY data
# exclude swedish_uc_old_gwas - no chrY data
# gwas1","prism_nfe_gwas","finland_illugwas","swedish_uc_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas","helmsley_prism_gsa","helmsley_xavier_prism_gsa


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
to_exclude<-c("swedish_uc_old_gwas","gwas1","prism_nfe_gwas","finland_illugwas","swedish_uc_old_gwas"
              ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","helmsley_prism_gsa"
              ,"helmsley_xavier_prism_gsa","niddk_old_gwas")
studies<-studies[which(!studies %in% to_exclude)]


for (j in 1:length(studies)) {
  
  print(studies[j])
  
  frq<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_check_sex_males_only_step24.frq",sep=""),head=T)
  var_miss<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_check_sex_males_only_step24.lmiss",sep=""),head=T)
  
  all<-merge(frq[,c("SNP","MAF","NCHROBS")],var_miss[,c("SNP","F_MISS")],by="SNP",sort=F)
  
  ## females, find variants with less % of calls:
  ped<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_check_sex_females_only_chr24_step24.raw",sep=""),head=T,check.names=F)
  
  dat<-matrix(nrow=nrow(all),ncol=2)
  dat<-as.data.frame(dat)
  colnames(dat)<-c("variant","percentage_NA")
  
  for (i in 1:nrow(dat)) {
    tmp<-ped[,6+i,drop=F]
    dat$variant[i]<-gsub("_[A-Z]{1}$","",colnames(tmp))
    dat$percentage_NA[i]<-nrow(tmp[which(is.na(tmp)),,drop=F])/nrow(tmp)
  }
  
  ### handle some outliers, Sands likely to have mixed up samples
  
  # if(studies[j] %in% c("sands_msccr_gsa")) {
  #   all<-all[which(all$SNP %in% dat$variant[which(dat$percentage_NA>0.40)]),]
  # } else if(studies[j] %in% c("farkkila_gsa","rioux_igenomed_gsa","bernstein_gsa")) {
  #   all<-all[which(all$SNP %in% dat$variant[which(dat$percentage_NA>0.75)]),]
  # } else {
  #   all<-all[which(all$SNP %in% dat$variant[which(dat$percentage_NA>0.95)]),]
  # }
  # 
  # keep variants with large number of calls in males
  
  # if(studies[j] %in% c("farkkila_gsa","rioux_igenomed_gsa","sands_msccr_gsa","bernstein_gsa")) {
  #   all<-all[which(all$NCHROBS>=max(all$NCHROBS,na.rm=T)-(max(all$NCHROBS,na.rm=T)*0.1)),]
  # } else if(studies[j] %in% c("slovenia_gsa")) {
  #   all<-all[which(all$NCHROBS>=max(all$NCHROBS,na.rm=T)-(max(all$NCHROBS,na.rm=T)*0.025)),]
  # }else {
  #   all<-all[which(all$NCHROBS>=max(all$NCHROBS,na.rm=T)-(max(all$NCHROBS,na.rm=T)*0.005)),]
  # }
  
  if(studies[j] %in% c("franchimont_gsa")) {
    all<-all[which(all$SNP %in% dat$variant[which(dat$percentage_NA>0.95)]),]
  } else {
    all<-all[which(all$SNP %in% dat$variant[which(dat$percentage_NA>0.975)]),]
  }
  if(studies[j] %in% c("rioux_igenomed_gsa")) {
    all<-all[which(all$NCHROBS>=max(all$NCHROBS,na.rm=T)-(max(all$NCHROBS,na.rm=T)*0.025)),]
  } else {
    all<-all[which(all$NCHROBS>=max(all$NCHROBS,na.rm=T)-(max(all$NCHROBS,na.rm=T)*0.005)),]
  }

  
  print(paste("N good ChrY variants:",nrow(all)))
  
  write.table(all,paste(path,"pre_imputation/QC/",studies[j],"/list_good_chrY_variants_step24",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
  rm(list=ls()[!ls() %in% c("j","studies","path")])
  
}

q("no")

##############################

# [1] "all_hce"
# [1] "N good ChrY variants: 988"
# [1] "australia_omniexome"
# [1] "N good ChrY variants: 437"
# [1] "gwas2"
# [1] "N good ChrY variants: 101"
# [1] "pittsburgh_gsa"
# [1] "N good ChrY variants: 504"
# [1] "spain_gsa"
# [1] "N good ChrY variants: 37"
# [1] "italy_gsa"
# [1] "N good ChrY variants: 257"
# [1] "kiel_austria_sibdcs_gsa"
# [1] "N good ChrY variants: 159"
# [1] "netherlands_gsa"
# [1] "N good ChrY variants: 127"
# [1] "slovenia_gsa"
# [1] "N good ChrY variants: 140"
# [1] "sweden_gsa"
# [1] "N good ChrY variants: 189"
# [1] "niddk_broad_gsa"
# [1] "N good ChrY variants: 162"
# [1] "niddk_feinstein_gsa"
# [1] "N good ChrY variants: 311"
# [1] "basque_gsa"
# [1] "N good ChrY variants: 286"
# [1] "lithuania_gsa"
# [1] "N good ChrY variants: 157"
# [1] "belgium_louis_gsa"
# [1] "N good ChrY variants: 192"
# [1] "belgium_franchimont_gsa"
# [1] "N good ChrY variants: 189"
# [1] "belgium_vermeire_gsa"
# [1] "N good ChrY variants: 371"
# [1] "prism_nfe_gsa"
# [1] "N good ChrY variants: 152"
# [1] "german_affy6_old_gwas"
# [1] "N good ChrY variants: 27"
# [1] "norway_affy6_old_gwas"
# [1] "N good ChrY variants: 136"
# [1] "cedars_370k_old_gwas"
# [1] "N good ChrY variants: 55"
# [1] "cedars_610k_old_gwas"
# [1] "N good ChrY variants: 59"
# [1] "cedars_omni_old_gwas"
# [1] "N good ChrY variants: 435"
# [1] "mccauley_gsa"
# [1] "N good ChrY variants: 143"
# [1] "ccfa_gsa"
# [1] "N good ChrY variants: 202"
# [1] "cedars_gsa"
# [1] "N good ChrY variants: 298"
# [1] "bernstein_gsa"
# [1] "N good ChrY variants: 113"
# [1] "farkkila_gsa"
# [1] "N good ChrY variants: 64"
# [1] "franchimont_gsa"
# [1] "N good ChrY variants: 297"
# [1] "franke_gsa"
# [1] "N good ChrY variants: 152"
# [1] "hyams_protect_gsa"
# [1] "N good ChrY variants: 167"
# [1] "lewis_sparc_gsa"
# [1] "N good ChrY variants: 223"
# [1] "mccauley_new_gsa"
# [1] "N good ChrY variants: 184"
# [1] "mcgovern_gsa"
# [1] "N good ChrY variants: 260"
# [1] "moayyedi_imagine_gsa"
# [1] "N good ChrY variants: 223"
# [1] "newberry_share_gsa"
# [1] "N good ChrY variants: 197"
# [1] "niddk_cho_gsa"
# [1] "N good ChrY variants: 186"
# [1] "niddk_duerr_gsa"
# [1] "N good ChrY variants: 196"
# [1] "niddk_rioux_gsa"
# [1] "N good ChrY variants: 174"
# [1] "niddk_silverberg_gsa"
# [1] "N good ChrY variants: 227"
# [1] "palotie_hus_gsa"
# [1] "N good ChrY variants: 150"
# [1] "pekow_share_gsa"
# [1] "N good ChrY variants: 160"
# [1] "rioux_igenomed_gsa"
# [1] "N good ChrY variants: 79"
# [1] "sands_msccr_gsa"
# [1] "N good ChrY variants: 190"
# [1] "stampfer_gsa"
# [1] "N good ChrY variants: 261"
# [1] "vermeire_gsa"
# [1] "N good ChrY variants: 384"
# [1] "weersma_gsa"
# [1] "N good ChrY variants: 143"
# [1] "xavier_prism_gsa"
# [1] "N good ChrY variants: 138"
# [1] "xavier_share_gsa"
# [1] "N good ChrY variants: 159"
##################################


path_gwas=/path/to/ibdgwas/IIBDGC/
  
for i in ${studies[@]}
do
cat ${path_gwas}pre_imputation/QC/${i}/list_good_chrX_variants_step24 ${path_gwas}pre_imputation/QC/${i}/list_good_chrY_variants_step24 > \
${path_gwas}pre_imputation/QC/${i}/list_good_chrXY_variants_step24
done

# cat: /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/list_good_chrY_variants_step24: No such file or directory
# cat: /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/list_good_chrY_variants_step24: No such file or directory
# cat: /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/list_good_chrY_variants_step24: No such file or directory
# cat: /path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/list_good_chrY_variants_step24: No such file or directory
# cat: /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/list_good_chrY_variants_step24: No such file or directory
# cat: /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/list_good_chrY_variants_step24: No such file or directory
# cat: /path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/list_good_chrX_variants_step24: No such file or directory
# cat: /path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/list_good_chrY_variants_step24: No such file or directory
# cat: /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/list_good_chrY_variants_step24: No such file or directory
# cat: /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/list_good_chrY_variants_step24: No such file or directory

for i in ${studies[@]}
do
echo ${i} && wc -l ${path_gwas}pre_imputation/QC/${i}/list_good_chrXY_variants_step24
done

##############################
# all_hce
# 5990 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/all_hce/list_good_chrXY_variants_step24
# niddk_old_gwas
# 2454 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/list_good_chrXY_variants_step24
# australia_omniexome
# 8938 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/list_good_chrXY_variants_step24
# gwas1
# 1946 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/list_good_chrXY_variants_step24
# gwas2
# 11658 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas2/list_good_chrXY_variants_step24
# pittsburgh_gsa
# 8880 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pittsburgh_gsa/list_good_chrXY_variants_step24
# spain_gsa
# 8953 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/spain_gsa/list_good_chrXY_variants_step24
# italy_gsa
# 10460 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/italy_gsa/list_good_chrXY_variants_step24
# kiel_austria_sibdcs_gsa
# 5441 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/kiel_austria_sibdcs_gsa/list_good_chrXY_variants_step24
# netherlands_gsa
# 4203 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/list_good_chrXY_variants_step24
# slovenia_gsa
# 6544 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/slovenia_gsa/list_good_chrXY_variants_step24
# sweden_gsa
# 6035 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/list_good_chrXY_variants_step24
# niddk_broad_gsa
# 3745 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_broad_gsa/list_good_chrXY_variants_step24
# niddk_feinstein_gsa
# 7213 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_feinstein_gsa/list_good_chrXY_variants_step24
# basque_gsa
# 9557 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/basque_gsa/list_good_chrXY_variants_step24
# lithuania_gsa
# 7789 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lithuania_gsa/list_good_chrXY_variants_step24
# belgium_louis_gsa
# 6733 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_louis_gsa/list_good_chrXY_variants_step24
# belgium_franchimont_gsa
# 6572 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_franchimont_gsa/list_good_chrXY_variants_step24
# belgium_vermeire_gsa
# 7731 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_vermeire_gsa/list_good_chrXY_variants_step24
# prism_nfe_gsa
# 8009 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gsa/list_good_chrXY_variants_step24
# prism_nfe_gwas
# 4297 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/list_good_chrXY_variants_step24
# finland_illugwas
# 4892 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/list_good_chrXY_variants_step24
# german_affy6_old_gwas
# 8354 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/german_affy6_old_gwas/list_good_chrXY_variants_step24
# norway_affy6_old_gwas
# 9988 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/norway_affy6_old_gwas/list_good_chrXY_variants_step24
# belgium_inf1_old_gwas
# 4285 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/list_good_chrXY_variants_step24
# belgium_inf2_old_gwas
# 4100 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/list_good_chrXY_variants_step24
# cedars_370k_old_gwas
# 6089 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_370k_old_gwas/list_good_chrXY_variants_step24
# cedars_610k_old_gwas
# 8228 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_610k_old_gwas/list_good_chrXY_variants_step24
# cedars_omni_old_gwas
# 10592 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_omni_old_gwas/list_good_chrXY_variants_step24
# swedish_uc_old_gwas
# 0 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/list_good_chrXY_variants_step24
# mccauley_gsa
# 7677 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_gsa/list_good_chrXY_variants_step24
# ccfa_gsa
# 8250 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/ccfa_gsa/list_good_chrXY_variants_step24
# cedars_gsa
# 7284 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_gsa/list_good_chrXY_variants_step24
# bernstein_gsa
# 8239 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/bernstein_gsa/list_good_chrXY_variants_step24
# farkkila_gsa
# 8144 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/farkkila_gsa/list_good_chrXY_variants_step24
# franchimont_gsa
# 7197 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franchimont_gsa/list_good_chrXY_variants_step24
# franke_gsa
# 7886 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franke_gsa/list_good_chrXY_variants_step24
# helmsley_prism_gsa
# 4551 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/list_good_chrXY_variants_step24
# helmsley_xavier_prism_gsa
# 4696 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/list_good_chrXY_variants_step24
# hyams_protect_gsa
# 8617 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/hyams_protect_gsa/list_good_chrXY_variants_step24
# lewis_sparc_gsa
# 8199 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lewis_sparc_gsa/list_good_chrXY_variants_step24
# mccauley_new_gsa
# 7324 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_new_gsa/list_good_chrXY_variants_step24
# mcgovern_gsa
# 7630 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mcgovern_gsa/list_good_chrXY_variants_step24
# moayyedi_imagine_gsa
# 7492 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/moayyedi_imagine_gsa/list_good_chrXY_variants_step24
# newberry_share_gsa
# 7820 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/newberry_share_gsa/list_good_chrXY_variants_step24
# niddk_cho_gsa
# 6889 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_cho_gsa/list_good_chrXY_variants_step24
# niddk_duerr_gsa
# 6801 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_duerr_gsa/list_good_chrXY_variants_step24
# niddk_rioux_gsa
# 7411 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_rioux_gsa/list_good_chrXY_variants_step24
# niddk_silverberg_gsa
# 8089 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_silverberg_gsa/list_good_chrXY_variants_step24
# palotie_hus_gsa
# 7287 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/palotie_hus_gsa/list_good_chrXY_variants_step24
# pekow_share_gsa
# 8178 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pekow_share_gsa/list_good_chrXY_variants_step24
# rioux_igenomed_gsa
# 9143 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/list_good_chrXY_variants_step24
# sands_msccr_gsa
# 7024 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sands_msccr_gsa/list_good_chrXY_variants_step24
# stampfer_gsa
# 7318 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/stampfer_gsa/list_good_chrXY_variants_step24
# vermeire_gsa
# 7692 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/vermeire_gsa/list_good_chrXY_variants_step24
# weersma_gsa
# 7516 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/weersma_gsa/list_good_chrXY_variants_step24
# xavier_prism_gsa
# 7734 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_prism_gsa/list_good_chrXY_variants_step24
# xavier_share_gsa
# 7733 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_share_gsa/list_good_chrXY_variants_step24
##############################

# keep only samples that pass all QC

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_24.6_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_24.6_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip \
--extract ${path_gwas}pre_imputation/QC/${i}/list_good_chrXY_variants_step24 --make-bed \
--keep ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_check_sex_good_chrXY_step24"
done

for i in ${studies[@]}
do
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_check_sex_good_chrXY_step24.bim
done

# swedish_uc_old_gwas
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_check_sex_good_chrXY.bim': No such file or directory


for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_24.7_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_24.7_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_check_sex_good_chrXY_step24 \
--check-sex ycount --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_check_sex_good_chrXY_step24"
done

for i in ${studies[@]}
do
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_check_sex_good_chrXY_step24.sexcheck
done

# swedish_uc_old_gwas
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_check_sex_good_chrXY.sexcheck': No such file or directory

# for i in ${studies[@]}
# do 
# bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
# -e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_24.8_${i} \
# -o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_24.8_${i} \
# "/path/to/software/plink_linux_x86_64_20181202/./plink \
# --bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_check_sex_good_chrXY_step24 \
# --check-sex --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_check_sex_good_chrX_step24"
# done
# for i in ${studies[@]}
# do
# echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_check_sex_good_chrX_step24.sexcheck
# done


# In 'ycount' mode, gender is still imputed from the X chromosome, but female calls are downgraded to ambiguous whenever more than 0 
# nonmissing Y genotypes are present, and male calls are downgraded when fewer than 0 are present. (Note that these are counts, not rates.) 
# These thresholds are controllable with --check-sex ycount's optional 3rd and 4th numeric parameters.



##############################
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

studies_new<-c("bernstein_gsa","farkkila_gsa","franchimont_gsa",
               "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
               "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
               "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
               "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
               "xavier_share_gsa")

################################################################################################
# RE-do case control status for all new studies (batches from 2021 and 2022)

pheno<-fread(paste(path,"pheno/gwas-mega-2-core-phenotypes-28ceede.csv.gz",sep=""),head=T)
pheno<-as.data.frame(pheno)

colnames(pheno)[6]<-"gender"

pheno$sex<-0
pheno$sex[which(pheno$gender=="Female")]<-2
pheno$sex[which(pheno$gender=="Male")]<-1
table(pheno$gender,pheno$sex,useNA="ifany")
#             0     1     2
#          3261     0     0
# Female      0     0 55427
# Male        0 51830     0
# Unknown  8872     0     0


pheno$pheno<-NA
pheno$pheno[which(pheno$control=="1")]<-"1"
pheno$pheno[which(pheno$affection=="Affected")]<-"2"
table(pheno$pheno,pheno$affection,useNA="ifany")
#            Affected Unaffected Unknown
# 1        0        0      28970       0
# 2        0    85835          0       0
# <NA>  3967        0        143     475
table(pheno$pheno,pheno$diag,useNA="ifany")
#            Crohn's Disease Indeterminate Ulcerative Colitis Unknown
# 1    28970               0             0                  0       0
# 2       17           49394          1221              31608    3595
# <NA>  4345               0             0                  0     240

################################################################
# collect all post-QC fam files, before being split into of ancestry

for (j in 1:length(studies)) {
  
  print(studies[j])
  fam_tmp<-fread(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam",sep=""),head=F)
  fam_tmp<-as.data.frame(fam_tmp)
  fam_tmp$study<-studies[j]
  
  # for studies pre-2021, keep data in fam fike
  if (!studies[j] %in% studies_new) {
    fam_tmp$sex<-fam_tmp$V5
    fam_tmp$pheno<-fam_tmp$V6
  } else {
    fam_tmp<-merge(fam_tmp,pheno[,c("sample_id","sex","pheno")],by.x="V1",by.y="sample_id",order=F,all.x=T)
    print(table(fam_tmp$V5,fam_tmp$sex,useNA="ifany"))
    print(table(fam_tmp$V6,fam_tmp$pheno,useNA="ifany"))
  }
  
  if(j==1) {
    fam_all<-fam_tmp
  }else{
    fam_all<-rbind(fam_all,fam_tmp)
  }
  
  rm(fam_tmp)
  
}


path<-"/path/to/ibdgwas/IIBDGC/"
studies<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
           "kiel_austria_sibdcs_gsa","netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa",
           "niddk_feinstein_gsa","basque_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa",
           "belgium_vermeire_gsa","prism_nfe_gsa","prism_nfe_gwas","finland_illugwas","german_affy6_old_gwas",
           "norway_affy6_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
           "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas",
           "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")
non_chry<-c("niddk_old_gwas","swedish_uc_old_gwas","gwas1","prism_nfe_gwas","finland_illugwas","swedish_uc_old_gwas"
            ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","helmsley_prism_gsa"
            ,"helmsley_xavier_prism_gsa")
# non chrx or y
to_exclude<-c("swedish_uc_old_gwas")
studies<-studies[which(!studies %in% to_exclude)]


dat<-as.data.frame(matrix(ncol=6,nrow=length(studies)))
colnames(dat)<-c("studies","N_samples","N_samples_set_sex_0_eur","N_samples_discordant_gender_eur","N_samples_set_sex_0_noneur","N_samples_discordant_gender_noneur")
dat$studies<-studies

for (j in 1:length(studies)) {
  
  print(studies[j])
  
  sex<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_check_sex_good_chrXY_step24.sexcheck",sep=""),head=T)
  sex<-merge(sex,fam_all[,c("V2","sex")],by.x="IID",by.y="V2",all.x=T)
  print(table(sex$PEDSEX,sex$sex))
  # print(table(sex$SNPSEX,sex$sex))

  sex$PEDSEX<-sex$sex
  sex$PEDSEX<-as.factor(sex$PEDSEX)
  
  # X thresholds:
  
  fmin_male<-mean(sex[which(sex$PEDSEX==1),"F"])-(10*sd(sex[which(sex$PEDSEX==1),"F"]))
  print(fmin_male)
  
  fmax_female<-mean(sex[which(sex$PEDSEX==2),"F"])+(10*sd(sex[which(sex$PEDSEX==2),"F"]))
  print(fmax_female)
  
  # # Y thresholds: 
  
  if (studies[j] %in% non_chry) {
    
    # no chrY
    
    p1n<-ggplot(sex[which(sex$PEDSEX==1),],aes(x=F,y=YCOUNT)) + geom_point(colour = "#00AFBB") + 
      xlim(min(sex$F),max(sex$F)) + 
      ylim(0,max(sex$YCOUNT)) +  
      # geom_hline(yintercept=ymin_male, linetype="dashed", color = "#00AFBB")+
      geom_vline(xintercept=fmin_male, linetype="dashed", color = "#00AFBB")
    
    p2n<-ggplot(sex[which(sex$PEDSEX==2),],aes(x=F,y=YCOUNT)) + geom_point(colour = "#E7B800") + 
      xlim(min(sex$F),max(sex$F)) + 
      ylim(0,max(sex$YCOUNT)) +  
      # geom_hline(yintercept=ymax_female, linetype="dashed", color = "#E7B800")+
      geom_vline(xintercept=fmax_female, linetype="dashed", color = "#E7B800")
    # + geom_vline(xintercept=fmin_female, linetype="dashed", color = "#E7B800")
    
    p3n<-ggplot(sex[which(sex$PEDSEX==0),],aes(x=F,y=YCOUNT)) + geom_point(colour = "#FC4E07") + 
      xlim(min(sex$F),max(sex$F)) + 
      ylim(0,max(sex$YCOUNT)) +  
      # geom_hline(yintercept=ymin_male, linetype="dashed", color = "#00AFBB")+
      geom_vline(xintercept=fmin_male, linetype="dashed", color = "#00AFBB")+ 
      # geom_hline(yintercept=ymax_female, linetype="dashed", color = "#E7B800")+
      geom_vline(xintercept=fmax_female, linetype="dashed", color = "#E7B800")
    # + geom_vline(xintercept=fmin_female, linetype="dashed", color = "#E7B800")
    
    pn<-ggarrange(p1n,p2n,p3n,ncol=3,labels=c("Male","Female","NA"))
    
    pdf(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_histogram_homozigosity_chrXY_with_newlimits_step24.pdf",sep=""),width =21)
    print(pn)
    dev.off()
    
  } else {
    
    ymin_male<-mean(sex[which(sex$PEDSEX==1),"YCOUNT"])-(10*sd(sex[which(sex$PEDSEX==1),"YCOUNT"]))
    print(ymin_male)
    ymax_female<-mean(sex[which(sex$PEDSEX==2),"YCOUNT"])+(10*sd(sex[which(sex$PEDSEX==2),"YCOUNT"]))
    print(ymax_female)
    
    # many mismatches in Farkkila make setting limits by observed distribution is appropriate - set up new ones after inspecting plot
    if(studies[j] %in% c("franchimont_gsa","rioux_igenomed_gsa","niddk_rioux_gsa")) {
      fmin_male<-mean(sex[which(sex$PEDSEX==1),"F"])-(4*sd(sex[which(sex$PEDSEX==1),"F"]))
      print(fmin_male)
      
      fmax_female<-mean(sex[which(sex$PEDSEX==2),"F"])+(4*sd(sex[which(sex$PEDSEX==2),"F"]))
      print(fmax_female)
      
      ymin_male<-mean(sex[which(sex$PEDSEX==1),"YCOUNT"])-(4*sd(sex[which(sex$PEDSEX==1),"YCOUNT"]))
      print(ymin_male)
      ymax_female<-mean(sex[which(sex$PEDSEX==2),"YCOUNT"])+(4*sd(sex[which(sex$PEDSEX==2),"YCOUNT"]))
      print(ymax_female)
      
    }
    
    # # many mismatches in Farkkila make setting limits by observed distribution is appropriate - set up new ones after inspecting plot
    # if(studies[j]=="franchimont_gsa") {
    #   fmin_male<-0.55
    #   fmax_female<-0.52
    #   ymin_male<-200
    #   ymax_female<-10
    # }
    
    # # many mismatches in Sands make setting limits by observed distribution is appropriate - set up new ones after inspecting plot
    # if(studies[j]=="sands_msccr_gsa") {
    #   fmin_male<-0.40
    #   fmax_female<-0.25
    #   ymin_male<-75
    #   ymax_female<-25
    # }
    
    # if(studies[j] %in% c("slovenia_gsa")) {
    #   
    #   fmin_male<-mean(sex[which(sex$PEDSEX==1),"F"])-(3*sd(sex[which(sex$PEDSEX==1),"F"]))
    #   print(fmin_male)
    #   
    #   fmax_female<-mean(sex[which(sex$PEDSEX==2),"F"])+(3*sd(sex[which(sex$PEDSEX==2),"F"]))
    #   print(fmax_female)
    #   
    #   ymin_male<-mean(sex[which(sex$PEDSEX==1),"YCOUNT"])-(3*sd(sex[which(sex$PEDSEX==1),"YCOUNT"]))
    #   print(ymin_male)
    #   ymax_female<-mean(sex[which(sex$PEDSEX==2),"YCOUNT"])+(3*sd(sex[which(sex$PEDSEX==2),"YCOUNT"]))
    #   print(ymax_female)
    #   
    # }
    
    
    p1n<-ggplot(sex[which(sex$PEDSEX==1),],aes(x=F,y=YCOUNT)) + geom_point(colour = "#00AFBB") + 
      xlim(min(sex$F),max(sex$F)) + 
      ylim(0,max(sex$YCOUNT)) +  
      geom_hline(yintercept=ymin_male, linetype="dashed", color = "#00AFBB")+
      geom_vline(xintercept=fmin_male, linetype="dashed", color = "#00AFBB")
    
    p2n<-ggplot(sex[which(sex$PEDSEX==2),],aes(x=F,y=YCOUNT)) + geom_point(colour = "#E7B800") + 
      xlim(min(sex$F),max(sex$F)) + 
      ylim(0,max(sex$YCOUNT)) +  
      geom_hline(yintercept=ymax_female, linetype="dashed", color = "#E7B800")+
      geom_vline(xintercept=fmax_female, linetype="dashed", color = "#E7B800")
    # + geom_vline(xintercept=fmin_female, linetype="dashed", color = "#E7B800")
    
    p3n<-ggplot(sex[which(sex$PEDSEX==0),],aes(x=F,y=YCOUNT)) + geom_point(colour = "#FC4E07") + 
      xlim(min(sex$F),max(sex$F)) + 
      ylim(0,max(sex$YCOUNT)) +  
      geom_hline(yintercept=ymin_male, linetype="dashed", color = "#00AFBB")+
      geom_vline(xintercept=fmin_male, linetype="dashed", color = "#00AFBB")+ 
      geom_hline(yintercept=ymax_female, linetype="dashed", color = "#E7B800")+
      geom_vline(xintercept=fmax_female, linetype="dashed", color = "#E7B800")
    # + geom_vline(xintercept=fmin_female, linetype="dashed", color = "#E7B800")
    
    pn<-ggarrange(p1n,p2n,p3n,ncol=3,labels=c("Male","Female","NA"))
    
    pdf(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_histogram_homozigosity_chrXY_with_newlimits_step24.pdf",sep=""),width =21)
    print(pn)
    dev.off()
    
  }
  
  
  ##############################

  if (studies[j] %in% non_chry) {
    sex$SNPSEX2<-0
    sex$SNPSEX2[which(sex$F<=fmax_female)]<-2
    sex$SNPSEX2[which(sex$F>=fmin_male)]<-1
  } else {
    sex$SNPSEX2<-0
    sex$SNPSEX2[which(sex$F<=fmax_female & sex$YCOUNT<=ymax_female)]<-2
    sex$SNPSEX2[which(sex$F>=fmin_male & sex$YCOUNT>=ymin_male)]<-1
  }

  print("Recorded vs inferred sex:")
  print(table(sex$sex,sex$SNPSEX2))
  
  samples_exclude<-sex[which((sex$sex==2 & sex$SNPSEX2==1) | (sex$sex==1 & sex$SNPSEX2==2)),]
  print(table(samples_exclude$sex,samples_exclude$SNPSEX2))
  
  print(paste("N samples exclude - wrong sex:",nrow(samples_exclude)))
  write.table(samples_exclude[,c("FID","IID","PEDSEX","SNPSEX2","F","YCOUNT")],paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_list_samples_sex_discrepancy_step24",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  

  # # keep summary
  dat$N_samples[j]<-nrow(sex)
  
  eur<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam",sep=""),head=F)
  dat$N_samples_set_sex_0_eur[j]<-nrow(sex[which(sex$SNPSEX2==0 & sex$FID %in% eur$V1),])
  dat$N_samples_discordant_gender_eur[j]<-nrow(samples_exclude[which(samples_exclude$FID %in% eur$V1),])
  
  file_noneur<-paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam",sep="")
  
  if(file.exists(file_noneur)){
    noneur<-read.table(file_noneur,head=F)
    dat$N_samples_set_sex_0_noneur[j]<-nrow(sex[which(sex$SNPSEX2==0 & sex$FID %in% noneur$V1),])
    dat$N_samples_discordant_gender_noneur[j]<-nrow(samples_exclude[which(samples_exclude$FID %in% noneur$V1),])
  }
  
  
  # keep record of new inferred sex:
  sex$study<-studies[j]
  
  if(j==1) {
    sex_all<-sex[,c("IID","FID","PEDSEX","SNPSEX2","study")]
  } else {
    sex<-sex[,c("IID","FID","PEDSEX","SNPSEX2","study")]
    sex_all<-rbind(sex,sex_all)
  }
  
  rm(list=ls()[!ls() %in% c("j","studies","path","non_gsa","dat","non_chry","fam_all","sex_all")])
  
}

write.table(dat,"~/tmp_plots/summary_sample_gender_sex_discrepancies_step24.tsv",col.names=T,row.names=F,quote=F,sep="\t")
fwrite(sex_all,paste(path,"pheno/gwas-mega-2-core_updated_sex_gender_Dec22.txt.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

dat

# create an empty file for 


samples_exclude<-as.data.frame(matrix(ncol=5,nrow=0))
colnames(samples_exclude)<-c("FID","PEDSEX","SNPSEX2","F","YCOUNT")
write.table(samples_exclude[,c("FID","PEDSEX","SNPSEX2","F","YCOUNT")],paste(path,"pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_list_samples_sex_discrepancy_step24",sep=""),col.names=T,row.names=F,quote=F,sep="\t")



q("no")

# check plots and results!
# 
# for i in ${studies[@]}
# do cp /path/to/ibdgwas/IIBDGC/pre_imputation/QC/${i}/${i}_histogram_homozigosity_chrXY_with_newlimits_step24.pdf ~/tmp_plots/
# done
# 
# for i in ${studies[@]}
# do cp /path/to/ibdgwas/IIBDGC/pre_imputation/QC/${i}/${i}_histogram_homozigosity_chrXY_with_newlimits.pdf ~/tmp_plots/
# done
# 


MEM=500
studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
path_gwas=/path/to/ibdgwas/IIBDGC/
  

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_24.1_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_24.1_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/${i}/${i}_list_samples_sex_discrepancy_step24 \
--make-bed --out ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24"
done


for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_24.1_${i} | grep -E "completed"
done

for i in ${studies[@]}
do ls -la ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24.fam
done
