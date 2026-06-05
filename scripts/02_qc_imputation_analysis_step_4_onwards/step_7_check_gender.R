# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
####################
# 7.- CHECK GENDER #
####################

#######################################
# 7.1 CREATE LIST OF MALES AND FEMALES

##############################
### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"
studies<-c("australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
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

for (j in 1:length(studies)) {
  
  print(studies[j])
  fam<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip.fam",sep=""),head=F)
  print(table(fam$V5,useNA="ifany"))
  
  fam.male<-fam[which(fam$V5==1),1:2] 
  colnames(fam.male)<-c("FID","IID")
  
  fam.female<-fam[which(fam$V5==2),1:2] 
  colnames(fam.female)<-c("FID","IID")
  
  if(nrow(fam.male)>0) {
    write.table(fam.male,paste(path,"pre_imputation/QC/",studies[j],"/list_male_samples",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  }
  if(nrow(fam.female)>0) {
    write.table(fam.female,paste(path,"pre_imputation/QC/",studies[j],"/list_female_samples",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  }
  
  rm(list=ls()[!ls() %in% c("j","studies","path")])
}

q("no")

##############################
# [1] "australia_omniexome"
# 1   2 
# 604 706 
# [1] "gwas1"
# 
# 1    2 
# 2126 2558 
# [1] "gwas2"
# 
# 1    2 
# 3908 3870 
# [1] "pittsburgh_gsa"
# 
# 0    1    2 
# 60 1193 1528 
# [1] "spain_gsa"
# 
# 1    2 
# 1912 1531 
# [1] "italy_gsa"
# 
# 1   2 
# 551 465 
# [1] "kiel_austria_sibdcs_gsa"
# 
# 0    1    2 
# 2 7250 7408 
# [1] "netherlands_gsa"
# 
# 1    2 
# 1947 2758 
# [1] "slovenia_gsa"
# 
# 0   1   2 
# 1  59 210 
# [1] "sweden_gsa"
# 
# 1   2 
# 821 593 
# [1] "niddk_broad_gsa"
# 
# 1    2 
# 2858 2657 
# [1] "niddk_feinstein_gsa"
# 
# 0    1    2 
# 8 4161 4132 
# [1] "basque_gsa"
# 
# 0   1   2 
# 5 931 580 
# [1] "lithuania_gsa"
# 
# 1    2 
# 1206 1071 
# [1] "belgium_louis_gsa"
# 
# 0   1   2 
# 21 662 848 
# [1] "belgium_franchimont_gsa"
# 
# 0   1   2 
# 1 666 865 
# [1] "belgium_vermeire_gsa"
# 
# 1    2 
# 1934 2080 
# [1] "prism_nfe_gsa"
# 
# 1   2 
# 218 248 
# [1] "prism_nfe_gwas"
# 
# 1   2 
# 403 459 
# [1] "finland_illugwas"
# 
# 1   2 
# 301 156 
# [1] "german_affy6_old_gwas"
# 
# 1    2 
# 1373 1454 
# [1] "norway_affy6_old_gwas"
# 
# 1   2 
# 305 245 
# [1] "belgium_inf1_old_gwas"
# 
# 0   1   2 
# 89 961 367 
# [1] "belgium_inf2_old_gwas"
# 
# 0   1   2 
# 115  59  98 
# [1] "cedars_370k_old_gwas"
# 
# 1   2 
# 325 283 
# [1] "cedars_610k_old_gwas"
# 
# 1   2 
# 466 427 
# [1] "cedars_omni_old_gwas"
# 
# 1   2 
# 596 630 
# [1] "swedish_uc_old_gwas"
# 
# 1   2 
# 697 567 
# [1] "mccauley_gsa"
# 
# 1   2 
# 404 384 
# [1] "ccfa_gsa"
# 
# 1    2 
# 1005 1183 
# [1] "cedars_gsa"
# 
# 1    2 
# 1496 1600 
# [1] "bernstein_gsa"
# 
# 0   1   2 
# 217 124 173 
# [1] "farkkila_gsa"
# 
# 1  2 
# 34 34 
# [1] "franchimont_gsa"
# 
# 0    1    2 
# 71 1274 1444 
# [1] "franke_gsa"
# 
# 1   2 
# 508 377 
# [1] "helmsley_prism_gsa"
# 
# 0   1   2 
# 2 382 396 
# [1] "helmsley_xavier_prism_gsa"
# 
# 0   1   2 
# 22 580 695 
# [1] "hyams_protect_gsa"
# 
# 1   2 
# 212 206 
# [1] "lewis_sparc_gsa"
# 
# 1    2 
# 1291 1568 
# [1] "mccauley_new_gsa"
# 
# 1   2 
# 845 783 
# [1] "mcgovern_gsa"
# 
# 0    1    2 
# 18 2965 3067 
# [1] "moayyedi_imagine_gsa"
# 
# 1   2 
# 490 638 
# [1] "newberry_share_gsa"
# 
# 1   2 
# 369 496 
# [1] "niddk_cho_gsa"
# 
# 0   1   2 
# 4 905 855 
# [1] "niddk_duerr_gsa"
# 
# 0   1   2 
# 3 969 972 
# [1] "niddk_rioux_gsa"
# 
# 0   1   2 
# 2 420 497 
# [1] "niddk_silverberg_gsa"
# 
# 0    1    2 
# 5 1212 1154 
# [1] "palotie_hus_gsa"
# 
# 1   2 
# 420 458 
# [1] "pekow_share_gsa"
# 
# 1   2 
# 317 317 
# [1] "rioux_igenomed_gsa"
# 
# 1  2 
# 87 95 
# [1] "sands_msccr_gsa"
# 
# 0   1   2 
# 1 765 664 
# [1] "stampfer_gsa"
# 
# 1    2 
# 239 1238 
# [1] "vermeire_gsa"
# 
# 1    2 
# 2280 2433 
# [1] "weersma_gsa"
# 
# 1   2 
# 355 354 
# [1] "xavier_prism_gsa"
# 
# 1   2 
# 343 349 
# [1] "xavier_share_gsa"
# 
# 1   2 
# 348 348 

##############################

###################################
# 7.2 KEEP ONLY FEMALES AND CHR23


path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_11_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_11_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip \
--keep ${path_gwas}pre_imputation/QC/${i}/list_female_samples --chr 23 --make-bed \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_females_only"
done

for i in ${studies[@]}
do echo ${i} && tail -100 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_11_${i} | grep -E "completed"
done

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_12_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_12_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_females_only \
--missing --hardy --freq --out ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_females_only"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_12_${i} | grep -E "completed"
done

for i in ${studies[@]}
do 
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_females_only.hwe
done

# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_check_sex_females_only.hwe': No such file or directory
# mccauley_gsa

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


for (j in 1:length(studies)) {
  
  print(studies[j])
  
  bim<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_check_sex_females_only.bim",sep=""),sep="\t",head=F)
  print(paste("N markers in chrX:",table(bim$V1)))

  hwe<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_check_sex_females_only.hwe",sep=""),head=T)
  hwe<-hwe[which(hwe$TEST=="UNAFF"),]
  
  frq<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_check_sex_females_only.frq",sep=""),head=T)
  var_miss<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_check_sex_females_only.lmiss",sep=""),head=T)
  
  all<-merge(hwe[,c("SNP","P")],frq[,c("SNP","MAF")],by="SNP",sort=F)
  all<-merge(all,var_miss[,c("SNP","F_MISS")],by="SNP",sort=F)
  
  all<-all[which(all$MAF>0.05 & all$F_MISS<0.01 & all$P>1E-4),]
  print(paste("N good chrX markers:",nrow(all)))
  
  write.table(all,paste(path,"pre_imputation/QC/",studies[j],"/list_good_chrX_variants",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
  rm(list=ls()[!ls() %in% c("j","studies","path")])
  
}

q("no")

##################################

# [1] "all_hce"
# [1] "N markers in chrX: 12477"
# [1] "N good chrX markers: 6445"
# [1] "niddk_old_gwas"
# [1] "N markers in chrX: 8968"
# [1] "N good chrX markers: 6193"

# [1] "australia_omniexome"
# [1] "N markers in chrX: 18586"
# [1] "N good chrX markers: 12791"
# [1] "gwas1"
# [1] "N markers in chrX: 9746"
# [1] "N good chrX markers: 6377"
# [1] "gwas2"
# [1] "N markers in chrX: 36382"
# [1] "N good chrX markers: 21440"
# [1] "pittsburgh_gsa"
# [1] "N markers in chrX: 22362"
# [1] "N good chrX markers: 15153"
# [1] "spain_gsa"
# [1] "N markers in chrX: 14735"
# [1] "N good chrX markers: 12287"
# [1] "italy_gsa"
# [1] "N markers in chrX: 24945"
# [1] "N good chrX markers: 15477"
# [1] "kiel_austria_sibdcs_gsa"
# [1] "N markers in chrX: 26756"
# [1] "N good chrX markers: 8928"
# [1] "netherlands_gsa"
# [1] "N markers in chrX: 17412"
# [1] "N good chrX markers: 8537"
# [1] "slovenia_gsa"
# [1] "N markers in chrX: 13858"
# [1] "N good chrX markers: 9198"
# [1] "sweden_gsa"
# [1] "N markers in chrX: 14315"
# [1] "N good chrX markers: 10197"
# [1] "niddk_broad_gsa"
# [1] "N markers in chrX: 15302"
# [1] "N good chrX markers: 9365"
# [1] "niddk_feinstein_gsa"
# [1] "N markers in chrX: 15451"
# [1] "N good chrX markers: 9363"
# [1] "basque_gsa"
# [1] "N markers in chrX: 24702"
# [1] "N good chrX markers: 17182"
# [1] "lithuania_gsa"
# [1] "N markers in chrX: 14423"
# [1] "N good chrX markers: 10103"
# [1] "belgium_louis_gsa"
# [1] "N markers in chrX: 15664"
# [1] "N good chrX markers: 10478"
# [1] "belgium_franchimont_gsa"
# [1] "N markers in chrX: 14798"
# [1] "N good chrX markers: 10140"
# [1] "belgium_vermeire_gsa"
# [1] "N markers in chrX: 16934"
# [1] "N good chrX markers: 10254"
# [1] "prism_nfe_gsa"
# [1] "N markers in chrX: 14522"
# [1] "N good chrX markers: 10613"
# [1] "prism_nfe_gwas"
# [1] "N markers in chrX: 6371"
# [1] "N good chrX markers: 498"
# [1] "finland_illugwas"
# [1] "N markers in chrX: 6056"
# [1] "N good chrX markers: 5274"
# [1] "german_affy6_old_gwas"
# [1] "N markers in chrX: 34415"
# [1] "N good chrX markers: 18651"
# [1] "norway_affy6_old_gwas"
# [1] "N markers in chrX: 31902"
# [1] "N good chrX markers: 17610"
# [1] "belgium_inf1_old_gwas"
# [1] "N markers in chrX: 7568"
# [1] "N good chrX markers: 6848"
# [1] "belgium_inf2_old_gwas"
# [1] "N markers in chrX: 7929"
# [1] "N good chrX markers: 6028"
# [1] "cedars_370k_old_gwas"
# [1] "N markers in chrX: 9924"
# [1] "N good chrX markers: 9086"
# [1] "cedars_610k_old_gwas"
# [1] "N markers in chrX: 14735"
# [1] "N good chrX markers: 12498"
# [1] "cedars_omni_old_gwas"
# [1] "N markers in chrX: 17670"
# [1] "N good chrX markers: 13983"
# [1] "mccauley_gsa"
# [1] "N markers in chrX: 14417"
# [1] "N good chrX markers: 10492"
# [1] "ccfa_gsa"
# [1] "N markers in chrX: 15164"
# [1] "N good chrX markers: 10991"
# [1] "cedars_gsa"
# [1] "N markers in chrX: 16799"
# [1] "N good chrX markers: 10344"
# [1] "bernstein_gsa"
# [1] "N markers in chrX: 14244"
# [1] "N good chrX markers: 10337"
# [1] "farkkila_gsa"
# [1] "N markers in chrX: 13238"
# [1] "N good chrX markers: 9895"
# [1] "franchimont_gsa"
# [1] "N markers in chrX: 17028"
# [1] "N good chrX markers: 9958"
# [1] "franke_gsa"
# [1] "N markers in chrX: 14212"
# [1] "N good chrX markers: 10296"
# [1] "helmsley_prism_gsa"
# [1] "N markers in chrX: 6371"
# [1] "N good chrX markers: 5322"
# [1] "helmsley_xavier_prism_gsa"
# [1] "N markers in chrX: 6371"
# [1] "N good chrX markers: 5487"
# [1] "hyams_protect_gsa"
# [1] "N markers in chrX: 14714"
# [1] "N good chrX markers: 11177"
# [1] "lewis_sparc_gsa"
# [1] "N markers in chrX: 15381"
# [1] "N good chrX markers: 11048"
# [1] "mccauley_new_gsa"
# [1] "N markers in chrX: 14945"
# [1] "N good chrX markers: 11015"
# [1] "mcgovern_gsa"
# [1] "N markers in chrX: 15958"
# [1] "N good chrX markers: 10398"
# [1] "moayyedi_imagine_gsa"
# [1] "N markers in chrX: 14969"
# [1] "N good chrX markers: 10814"
# [1] "newberry_share_gsa"
# [1] "N markers in chrX: 15057"
# [1] "N good chrX markers: 11039"
# [1] "niddk_cho_gsa"
# [1] "N markers in chrX: 15106"
# [1] "N good chrX markers: 10802"
# [1] "niddk_duerr_gsa"
# [1] "N markers in chrX: 14975"
# [1] "N good chrX markers: 10604"
# [1] "niddk_rioux_gsa"
# [1] "N markers in chrX: 14769"
# [1] "N good chrX markers: 10605"
# [1] "niddk_silverberg_gsa"
# [1] "N markers in chrX: 15313"
# [1] "N good chrX markers: 10753"
# [1] "palotie_hus_gsa"
# [1] "N markers in chrX: 14252"
# [1] "N good chrX markers: 10449"
# [1] "pekow_share_gsa"
# [1] "N markers in chrX: 14678"
# [1] "N good chrX markers: 11240"
# [1] "rioux_igenomed_gsa"
# [1] "N markers in chrX: 13960"
# [1] "N good chrX markers: 10690"
# [1] "sands_msccr_gsa"
# [1] "N markers in chrX: 14937"
# [1] "N good chrX markers: 1531"

# [1] "N good chrX markers: 1531" # what happened here!

# > dim(all[which(all$MAF>0.05),])
# [1] 11692     4
# > dim(all[which(all$F_MISS<0.01),])
# [1] 14277     4
# > dim(all[which(all$P>1E-4),])
# [1] 4217    4 # driven by low hwe pvalues, N samples = 1,474 not particular large N females~ 660 - keep an eye on this cohort


# [1] "stampfer_gsa"
# [1] "N markers in chrX: 17138"
# [1] "N good chrX markers: 10158"
# [1] "vermeire_gsa"
# [1] "N markers in chrX: 16946"
# [1] "N good chrX markers: 10279"
# [1] "weersma_gsa"
# [1] "N markers in chrX: 14441"
# [1] "N good chrX markers: 10586"
# [1] "xavier_prism_gsa"
# [1] "N markers in chrX: 14668"
# [1] "N good chrX markers: 10678"
# [1] "xavier_share_gsa"
# [1] "N markers in chrX: 14720"
# [1] "N good chrX markers: 10644"

##################################

##################################
# 7.3 KEEP ONLY MALES AND CHR24

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_13_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_13_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip \
--keep ${path_gwas}pre_imputation/QC/${i}/list_male_samples --chr 24 --make-bed \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_males_only"
done

for i in ${studies[@]}
do 
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_males_only.bim
done

# OK - no chrY
# niddk_old_gwas
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_check_sex_males_only.bim': No such file or directory
# gwas1
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/gwas1_check_sex_males_only.bim': No such file or directory
# prism_nfe_gwas
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_check_sex_males_only.bim': No such file or directory
# finland_illugwas
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/finland_illugwas_check_sex_males_only.bim': No such file or directory
# swedish_uc_old_gwas
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_check_sex_males_only.bim': No such file or directory
# belgium_inf1_old_gwas
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_check_sex_males_only.bim': No such file or directory
# belgium_inf2_old_gwas
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_check_sex_males_only.bim': No such file or directory
# helmsley_prism_gsa
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_check_sex_males_only.bim': No such file or directory
# helmsley_xavier_prism_gsa
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_check_sex_males_only.bim': No such file or directory
# hyams_protect_gsa


for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_14_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_14_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_males_only \
--freq --missing \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_males_only"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_14_${i} | grep -E "completed"
done

for i in ${studies[@]}
do 
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_males_only.frq
done


#  AND FIND VARIANTS WITH NO CALLS IN MOST OF FEMALES

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_15_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_15_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip \
--keep ${path_gwas}pre_imputation/QC/${i}/list_female_samples --chr 24 --recode A \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_females_only_chr24"
done

for i in ${studies[@]}
do 
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_check_sex_females_only_chr24.raw
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
  
  frq<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_check_sex_males_only.frq",sep=""),head=T)
  var_miss<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_check_sex_males_only.lmiss",sep=""),head=T)
  
  all<-merge(frq[,c("SNP","MAF","NCHROBS")],var_miss[,c("SNP","F_MISS")],by="SNP",sort=F)
  
  ## females, find variants with less % of calls:
  ped<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_check_sex_females_only_chr24.raw",sep=""),head=T,check.names=F)
  
  dat<-matrix(nrow=nrow(all),ncol=2)
  dat<-as.data.frame(dat)
  colnames(dat)<-c("variant","percentage_NA")
  
  for (i in 1:nrow(dat)) {
    tmp<-ped[,6+i,drop=F]
    dat$variant[i]<-gsub("_[A-Z]{1}$","",colnames(tmp))
    dat$percentage_NA[i]<-nrow(tmp[which(is.na(tmp)),,drop=F])/nrow(tmp)
  }
  
  ### handle some outliers, Sands likely to have mixed up samples
  
  if(studies[j] %in% c("sands_msccr_gsa")) {
    all<-all[which(all$SNP %in% dat$variant[which(dat$percentage_NA>0.40)]),]
  } else if(studies[j] %in% c("farkkila_gsa","rioux_igenomed_gsa","bernstein_gsa")) {
    all<-all[which(all$SNP %in% dat$variant[which(dat$percentage_NA>0.75)]),]
  } else {
    all<-all[which(all$SNP %in% dat$variant[which(dat$percentage_NA>0.95)]),]
  }
  
  # keep variants with large number of calls in males
  
  if(studies[j] %in% c("farkkila_gsa","rioux_igenomed_gsa","sands_msccr_gsa","bernstein_gsa")) {
    all<-all[which(all$NCHROBS>=max(all$NCHROBS,na.rm=T)-(max(all$NCHROBS,na.rm=T)*0.1)),]
  } else if(studies[j] %in% c("slovenia_gsa")) {
    all<-all[which(all$NCHROBS>=max(all$NCHROBS,na.rm=T)-(max(all$NCHROBS,na.rm=T)*0.025)),]
  }else {
    all<-all[which(all$NCHROBS>=max(all$NCHROBS,na.rm=T)-(max(all$NCHROBS,na.rm=T)*0.005)),]
  }
  
  print(paste("N good ChrY variants:",nrow(all)))
  
  write.table(all,paste(path,"pre_imputation/QC/",studies[j],"/list_good_chrY_variants",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
  rm(list=ls()[!ls() %in% c("j","studies","path")])
  
}

q("no")

##############################

# [1] "all_hce"
# [1] "N good ChrY variants: 503"
# [1] "australia_omniexome"
# [1] "N good ChrY variants: 453"
# [1] "gwas2"
# [1] "N good ChrY variants: 91"
# [1] "pittsburgh_gsa"
# [1] "N good ChrY variants: 279"
# [1] "spain_gsa"
# [1] "N good ChrY variants: 49"
# [1] "italy_gsa"
# [1] "N good ChrY variants: 59"
# [1] "kiel_austria_sibdcs_gsa"
# [1] "N good ChrY variants: 161"
# [1] "netherlands_gsa"
# [1] "N good ChrY variants: 125"
# [1] "slovenia_gsa"
# [1] "N good ChrY variants: 31"
# [1] "sweden_gsa"
# [1] "N good ChrY variants: 205"
# [1] "niddk_broad_gsa"
# [1] "N good ChrY variants: 179"
# [1] "niddk_feinstein_gsa"
# [1] "N good ChrY variants: 360"
# [1] "basque_gsa"
# [1] "N good ChrY variants: 359"
# [1] "lithuania_gsa"
# [1] "N good ChrY variants: 168"
# [1] "belgium_louis_gsa"
# [1] "N good ChrY variants: 145"
# [1] "belgium_franchimont_gsa"
# [1] "N good ChrY variants: 146"
# [1] "belgium_vermeire_gsa"
# [1] "N good ChrY variants: 323"
# [1] "prism_nfe_gsa"
# [1] "N good ChrY variants: 115"
# [1] "german_affy6_old_gwas"
# [1] "N good ChrY variants: 27"
# [1] "norway_affy6_old_gwas"
# [1] "N good ChrY variants: 131"
# [1] "cedars_370k_old_gwas"
# [1] "N good ChrY variants: 55"
# [1] "cedars_610k_old_gwas"
# [1] "N good ChrY variants: 59"
# [1] "cedars_omni_old_gwas"
# [1] "N good ChrY variants: 435"
# [1] "mccauley_gsa"
# [1] "N good ChrY variants: 128"
# [1] "ccfa_gsa"
# [1] "N good ChrY variants: 223"
# [1] "cedars_gsa"
# [1] "N good ChrY variants: 311"
# [1] "bernstein_gsa"
# [1] "N good ChrY variants: 299"
# [1] "farkkila_gsa"
# [1] "N good ChrY variants: 13"
# [1] "franchimont_gsa"
# [1] "N good ChrY variants: 292"
# [1] "franke_gsa"
# [1] "N good ChrY variants: 57"
# [1] "hyams_protect_gsa"
# [1] "N good ChrY variants: 188"
# [1] "lewis_sparc_gsa"
# [1] "N good ChrY variants: 240"
# [1] "mccauley_new_gsa"
# [1] "N good ChrY variants: 188"
# [1] "mcgovern_gsa"
# [1] "N good ChrY variants: 289"
# [1] "moayyedi_imagine_gsa"
# [1] "N good ChrY variants: 255"
# [1] "newberry_share_gsa"
# [1] "N good ChrY variants: 163"
# [1] "niddk_cho_gsa"
# [1] "N good ChrY variants: 222"
# [1] "niddk_duerr_gsa"
# [1] "N good ChrY variants: 227"
# [1] "niddk_rioux_gsa"
# [1] "N good ChrY variants: 212"
# [1] "niddk_silverberg_gsa"
# [1] "N good ChrY variants: 251"
# [1] "palotie_hus_gsa"
# [1] "N good ChrY variants: 184"
# [1] "pekow_share_gsa"
# [1] "N good ChrY variants: 80"
# [1] "rioux_igenomed_gsa"
# [1] "N good ChrY variants: 248"
# [1] "sands_msccr_gsa"
# [1] "N good ChrY variants: 93"
# [1] "stampfer_gsa"
# [1] "N good ChrY variants: 294"
# [1] "vermeire_gsa"
# [1] "N good ChrY variants: 354"
# [1] "weersma_gsa"
# [1] "N good ChrY variants: 178"
# [1] "xavier_prism_gsa"
# [1] "N good ChrY variants: 90"
# [1] "xavier_share_gsa"
# [1] "N good ChrY variants: 190"

##################################


path_gwas=/path/to/ibdgwas/IIBDGC/

for i in ${studies[@]}
do
cat ${path_gwas}pre_imputation/QC/${i}/list_good_chrX_variants ${path_gwas}pre_imputation/QC/${i}/list_good_chrY_variants > \
${path_gwas}pre_imputation/QC/${i}/list_good_chrXY_variants
done


for i in ${studies[@]}
do
echo ${i} && wc -l ${path_gwas}pre_imputation/QC/${i}/list_good_chrXY_variants
done

##############################
# all_hce
# 6948 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/all_hce/list_good_chrXY_variants
# niddk_old_gwas
# 6193 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/list_good_chrXY_variants

# australia_omniexome
# 13244 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/list_good_chrXY_variants
# gwas1
# 6377 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/list_good_chrXY_variants
# gwas2
# 21531 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas2/list_good_chrXY_variants
# pittsburgh_gsa
# 15432 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pittsburgh_gsa/list_good_chrXY_variants
# spain_gsa
# 12336 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/spain_gsa/list_good_chrXY_variants
# italy_gsa
# 15536 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/italy_gsa/list_good_chrXY_variants
# kiel_austria_sibdcs_gsa
# 9089 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/kiel_austria_sibdcs_gsa/list_good_chrXY_variants
# netherlands_gsa
# 8662 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/list_good_chrXY_variants
# slovenia_gsa
# 9229 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/slovenia_gsa/list_good_chrXY_variants
# sweden_gsa
# 10402 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/list_good_chrXY_variants
# niddk_broad_gsa
# 9544 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_broad_gsa/list_good_chrXY_variants
# niddk_feinstein_gsa
# 9723 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_feinstein_gsa/list_good_chrXY_variants
# basque_gsa
# 17541 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/basque_gsa/list_good_chrXY_variants
# lithuania_gsa
# 10271 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lithuania_gsa/list_good_chrXY_variants
# belgium_louis_gsa
# 10623 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_louis_gsa/list_good_chrXY_variants
# belgium_franchimont_gsa
# 10286 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_franchimont_gsa/list_good_chrXY_variants
# belgium_vermeire_gsa
# 10577 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_vermeire_gsa/list_good_chrXY_variants
# prism_nfe_gsa
# 10728 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gsa/list_good_chrXY_variants
# prism_nfe_gwas
# 498 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/list_good_chrXY_variants
# finland_illugwas
# 5274 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/list_good_chrXY_variants
# german_affy6_old_gwas
# 18678 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/german_affy6_old_gwas/list_good_chrXY_variants
# norway_affy6_old_gwas
# 17741 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/norway_affy6_old_gwas/list_good_chrXY_variants
# belgium_inf1_old_gwas
# 6848 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/list_good_chrXY_variants
# belgium_inf2_old_gwas
# 6028 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/list_good_chrXY_variants
# cedars_370k_old_gwas
# 9141 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_370k_old_gwas/list_good_chrXY_variants
# cedars_610k_old_gwas
# 12557 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_610k_old_gwas/list_good_chrXY_variants
# cedars_omni_old_gwas
# 14418 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_omni_old_gwas/list_good_chrXY_variants
# swedish_uc_old_gwas
# 0 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/list_good_chrXY_variants
# mccauley_gsa
# 10620 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_gsa/list_good_chrXY_variants
# ccfa_gsa
# 11214 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/ccfa_gsa/list_good_chrXY_variants
# cedars_gsa
# 10655 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_gsa/list_good_chrXY_variants
# bernstein_gsa
# 10636 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/bernstein_gsa/list_good_chrXY_variants
# farkkila_gsa
# 9908 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/farkkila_gsa/list_good_chrXY_variants
# franchimont_gsa
# 10250 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franchimont_gsa/list_good_chrXY_variants
# franke_gsa
# 10353 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franke_gsa/list_good_chrXY_variants
# helmsley_prism_gsa
# 5322 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/list_good_chrXY_variants
# helmsley_xavier_prism_gsa
# 5487 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/list_good_chrXY_variants
# hyams_protect_gsa
# 11365 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/hyams_protect_gsa/list_good_chrXY_variants
# lewis_sparc_gsa
# 11288 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lewis_sparc_gsa/list_good_chrXY_variants
# mccauley_new_gsa
# 11203 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_new_gsa/list_good_chrXY_variants
# mcgovern_gsa
# 10687 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mcgovern_gsa/list_good_chrXY_variants
# moayyedi_imagine_gsa
# 11069 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/moayyedi_imagine_gsa/list_good_chrXY_variants
# newberry_share_gsa
# 11202 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/newberry_share_gsa/list_good_chrXY_variants
# niddk_cho_gsa
# 11024 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_cho_gsa/list_good_chrXY_variants
# niddk_duerr_gsa
# 10831 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_duerr_gsa/list_good_chrXY_variants
# niddk_rioux_gsa
# 10817 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_rioux_gsa/list_good_chrXY_variants
# niddk_silverberg_gsa
# 11004 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_silverberg_gsa/list_good_chrXY_variants
# palotie_hus_gsa
# 10633 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/palotie_hus_gsa/list_good_chrXY_variants
# pekow_share_gsa
# 11320 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pekow_share_gsa/list_good_chrXY_variants
# rioux_igenomed_gsa
# 10938 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/list_good_chrXY_variants
# sands_msccr_gsa
# 1624 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sands_msccr_gsa/list_good_chrXY_variants
# stampfer_gsa
# 10452 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/stampfer_gsa/list_good_chrXY_variants
# vermeire_gsa
# 10633 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/vermeire_gsa/list_good_chrXY_variants
# weersma_gsa
# 10764 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/weersma_gsa/list_good_chrXY_variants
# xavier_prism_gsa
# 10768 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_prism_gsa/list_good_chrXY_variants
# xavier_share_gsa
# 10834 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_share_gsa/list_good_chrXY_variants

##############################

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_16_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_16_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip \
--extract ${path_gwas}pre_imputation/QC/${i}/list_good_chrXY_variants --make-bed \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_check_sex_good_chrXY"
done

for i in ${studies[@]}
do
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_check_sex_good_chrXY.bim
done

# swedish_uc_old_gwas
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_check_sex_good_chrXY.bim': No such file or directory


for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_17_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_17_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_check_sex_good_chrXY \
--check-sex ycount --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_check_sex_good_chrXY"
done

for i in ${studies[@]}
do
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_check_sex_good_chrXY.sexcheck
done

# swedish_uc_old_gwas
# ls: cannot access '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_check_sex_good_chrXY.sexcheck': No such file or directory

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_18_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_18_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_check_sex_good_chrXY \
--check-sex --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_check_sex_good_chrX"
done

for i in ${studies[@]}
do
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_check_sex_good_chrX.sexcheck
done


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
to_exclude<-c("swedish_uc_old_gwas","")
studies<-studies[which(!studies %in% to_exclude)]


dat<-as.data.frame(matrix(ncol=4,nrow=length(studies)))
colnames(dat)<-c("studies","N_samples","N_samples_set_sex_0","N_samples_discordant_gender")
dat$studies<-studies

for (j in 1:length(studies)) {
  
  print(studies[j])

  sex<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_check_sex_good_chrXY.sexcheck",sep=""),head=T)
  
  sex$PEDSEX<-as.factor(sex$PEDSEX)
  
  # X thresholds:
  
  fmin_male<-mean(sex[which(sex$PEDSEX==1),"F"])-(4*sd(sex[which(sex$PEDSEX==1),"F"]))
  print(fmin_male)

  fmax_female<-mean(sex[which(sex$PEDSEX==2),"F"])+(4*sd(sex[which(sex$PEDSEX==2),"F"]))
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

    pdf(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_histogram_homozigosity_chrXY_with_newlimits.pdf",sep=""),width =21)
    print(pn)
    dev.off()
    
  } else {
    
    ymin_male<-mean(sex[which(sex$PEDSEX==1),"YCOUNT"])-(4*sd(sex[which(sex$PEDSEX==1),"YCOUNT"]))
    print(ymin_male)
    ymax_female<-mean(sex[which(sex$PEDSEX==2),"YCOUNT"])+(4*sd(sex[which(sex$PEDSEX==2),"YCOUNT"]))
    print(ymax_female)
    
    # many mismatches in Farkkila make setting limits by observed distribution is appropriate - set up new ones after inspecting plot
    if(studies[j]=="farkkila_gsa") {
      fmin_male<-0.50
      fmax_female<-0.25
      ymin_male<-8
      ymax_female<-5
    }
    
    # many mismatches in Farkkila make setting limits by observed distribution is appropriate - set up new ones after inspecting plot
    if(studies[j]=="franchimont_gsa") {
      fmin_male<-0.55
      fmax_female<-0.52
      ymin_male<-200
      ymax_female<-10
    }
    
    # many mismatches in Sands make setting limits by observed distribution is appropriate - set up new ones after inspecting plot
    if(studies[j]=="sands_msccr_gsa") {
      fmin_male<-0.40
      fmax_female<-0.25
      ymin_male<-75
      ymax_female<-25
    }
    
    if(studies[j] %in% c("slovenia_gsa")) {
      
      fmin_male<-mean(sex[which(sex$PEDSEX==1),"F"])-(3*sd(sex[which(sex$PEDSEX==1),"F"]))
      print(fmin_male)
      
      fmax_female<-mean(sex[which(sex$PEDSEX==2),"F"])+(3*sd(sex[which(sex$PEDSEX==2),"F"]))
      print(fmax_female)
      
      ymin_male<-mean(sex[which(sex$PEDSEX==1),"YCOUNT"])-(3*sd(sex[which(sex$PEDSEX==1),"YCOUNT"]))
      print(ymin_male)
      ymax_female<-mean(sex[which(sex$PEDSEX==2),"YCOUNT"])+(3*sd(sex[which(sex$PEDSEX==2),"YCOUNT"]))
      print(ymax_female)
      
    }
    
    
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
    
    pdf(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_histogram_homozigosity_chrXY_with_newlimits.pdf",sep=""),width =21)
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
  
  # sex_chrx<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_check_sex_good_chrX.sexcheck",sep=""),head=T)
  # colnames(sex_chrx)[4]<-"SNPSEX_chrx"
  # sexm<-merge(sex,sex_chrx[,c("FID","SNPSEX_chrx")],by="FID",sort=F)

  ## update data from pittsburgh -  udpate on sex according to Rich duerr updates
  if(studies[j]=="pittsburgh_gsa") {
   
    a<-read.csv(paste(path,"pre_imputation/QC/id_store/id_matches_2020/files_iibdgc/from_suppliers/pittsburgh_pheno_f41f3bb_Duerr_update_04202020_missing.csv",sep=""),head=T)
    b<-read.csv(paste(path,"pre_imputation/QC/id_store/id_matches_2020/files_iibdgc/from_suppliers/pittsburgh_pheno_f41f3bb_Duerr_update_04202020_complete.csv",sep=""),head=T)
    a<-rbind(a[,c("sample_id","sex")],b[,c("sample_id","sex")])
    a$PEDSEX<-NA
    a$PEDSEX[which(a$sex=="Female")]<-2
    a$PEDSEX[which(a$sex=="Male")]<-1
    
    sex<-merge(sex,a[,c("sample_id","PEDSEX")],by.x="FID",by.y="sample_id",all.x=T)
    sex<-sex[,c("FID","IID","PEDSEX.y","SNPSEX","STATUS","F","YCOUNT","SNPSEX2")]
    colnames(sex)<-c("FID","IID","PEDSEX","SNPSEX","STATUS","F","YCOUNT","SNPSEX2")

  }
  
  ## update data from italy -  udpate on sex according to Anna
  if(studies[j]=="italy_gsa") {
    
    path<-"/path/to/ibdgwas/IIBDGC/"
    a<-read.table(paste(path,"pre_imputation/QC/italy_gsa/italy_gsa_hg19_edited.fam",sep=""),head=F)
    
    files<-list.files(paste(path,"pre_imputation/QC/id_store/id_matches_2020/maps/",sep=""))
    files<-files[grep("italy",files)]
    b<-read.table(paste(path,"pre_imputation/QC/id_store/id_matches_2020/maps/",files[1],sep=""),head=T,sep="\t")
    colnames(b)
    colnames(b)<-c("Sample_id","Sample_name","IIBDGC_id")
    
    all<-merge(b,a,by.x="IIBDGC_id",by.y="V1")
    
    # NEW FILE FROM ANNA
    c<-read.csv(paste(path,"pre_imputation/QC/id_store/id_matches_2020/files_iibdgc/from_suppliers/LAuraobservations.csv",sep=""),head=T)
    dim(c)
    all<-merge(all,c[,c("subject_id","sex")],by.x="Sample_name",by.y="subject_id",all.x=T)
    
    all$PEDSEX<-all$sex
    
    sex<-merge(sex,all[,c("IIBDGC_id","PEDSEX")],by.x="FID",by.y="IIBDGC_id",all.x=T)
    sex<-sex[,c("FID","IID","PEDSEX.y","SNPSEX","STATUS","F","YCOUNT","SNPSEX2")]
    colnames(sex)<-c("FID","IID","PEDSEX","SNPSEX","STATUS","F","YCOUNT","SNPSEX2")
    
  }
  
  ## update data from Australia -  udpate on sex according to Graham
  if(studies[j]=="australia_omniexome") {
    
    a<-read.csv(paste(path,"pre_imputation/QC/id_store/id_matches_2020/files_iibdgc/from_suppliers/australia_ommiexome_pheno_updates.csv",sep=""),head=T)
    a$PEDSEX<-NA
    a$PEDSEX[which(a$PEDSEX..Phenotype.reported...checked.and.confirmed..1.Male.2.Female=="2 (Female)")]<-2
    a$PEDSEX[which(a$PEDSEX..Phenotype.reported...checked.and.confirmed..1.Male.2.Female=="1 (Male)")]<-1
    
    sex$PEDSEX[which(sex$FID %in% a$IIBDGC_id[which(a$PEDSEX==a$SNPSEX2)])]<-
      sex$SNPSEX2[which(sex$FID %in% a$IIBDGC_id[which(a$PEDSEX==a$SNPSEX2)])]
    
  }
  
  print("Recorded vs inferred sex:")
  print(table(sex$PEDSEX,sex$SNPSEX2))
  
  samples_exclude<-sex[which((sex$PEDSEX==2 & sex$SNPSEX2==1) | (sex$PEDSEX==1 & sex$SNPSEX2==2)),]
  # print(table(samples_exclude$PEDSEX,samples_exclude$SNPSEX2))
  
  print(paste("N samples exclude - wrong sex:",nrow(samples_exclude)))
  write.table(samples_exclude[,c("FID","PEDSEX","SNPSEX2","F","YCOUNT")],paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_list_samples_sex_discrepancy",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  
  # recode the rest:
  
  fam<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip.fam",sep=""),head=F)
  fam_ed<-merge(fam,sex[,c("FID","SNPSEX2")],by.x="V1",by.y="FID",sort=F)
  
  print(table(fam_ed$V2==fam$V2))
  
  print(paste("N samples sex set to 0:",nrow(fam_ed[which(fam_ed$SNPSEX2==0),])))
  
  # keep summary
  dat$N_samples[j]<-nrow(fam)
  dat$N_samples_set_sex_0[j]<-nrow(fam_ed[which(fam_ed$SNPSEX2==0),])
  dat$N_samples_discordant_gender[j]<-nrow(samples_exclude)
  
  fam_ed<-fam_ed[,c("V1","V2","V3","V4","SNPSEX2","V6")]
  
  # set temporarily sex = 0 the samples in the exclusion list - add new step later to exclude samples for the analysis
  fam_ed$SNPSEX2[which(fam_ed$V1 %in% samples_exclude$FID)]<-0
  
  write.table(fam_ed,paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_edited.fam",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
  # N cases and ctr excluded:
  print("N cases and ctr excluded:")
  print(table(fam$V6[which(fam$V1 %in% samples_exclude$FID)]))
  
  rm(list=ls()[!ls() %in% c("j","studies","path","non_gsa","dat","non_chry")])
  
}

write.table(dat,"~/tmp_plots/summary_sample_gender_sex_discrepancies.tsv",col.names=T,row.names=F,quote=F,sep="\t")
dat

### create one for swedish_uc_old_gwas
studies<-"swedish_uc_old_gwas"

for (j in 1:length(studies)) {
  fam<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip.fam",sep=""),head=F)
  write.table(fam,paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_edited.fam",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
}

q("no")

# check plots and results!
  
############ 

# [1] "all_hce"
# [1] 0.3699355
# [1] 0.5932043
# [1] 227.0763
# [1] 279.6589
# [1] "Recorded vs inferred sex:"
# 
# 0     1     2
# 0    29  1725  1890
# 1    75  9348   160
# 2    65   187 10872
# [1] "N samples exclude - wrong sex: 347"
# 
# TRUE 
# 24351 
# [1] "N samples sex set to 0: 169"
# [1] "N cases and ctr excluded:"
# 
# 1   2 
# 27 320 

# [1] "niddk_old_gwas"
# [1] 0.5199841
# [1] 0.5118007
# [1] "Recorded vs inferred sex:"
# 
# 1    2
# 1 1410   12
# 2   17 1405
# [1] "N samples exclude - wrong sex: 29"
# 
# TRUE 
# 2844 
# [1] "N samples sex set to 0: 0"
# [1] "N cases and ctr excluded:"
# 
# 1  2 
# 7 22 


# [1] "australia_omniexome"
# [1] 0.616271
# [1] 0.4883122
# [1] 284.2508
# [1] 210.509
# [1] "Recorded vs inferred sex:"
#     0   1   2
# 1   1 606   3
# 2   0   1 699
# [1] "N samples exclude - wrong sex: 4"
# 
# TRUE 
# 1310 
# [1] "N samples sex set to 0: 1"
# [1] "N cases and ctr excluded:"
# 
# 1  2 
# 2 12 


# [1] "gwas1"
# [1] 0.8510928
# [1] 0.2060789
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 1    0 2123    3
# 2    2    1 2555
# [1] "N samples exclude - wrong sex: 4"
# 
# TRUE 
# 4684 
# [1] "N samples sex set to 0: 2"
# [1] "N cases and ctr excluded:"
# 
# 1 2 
# 3 1 

# [1] "gwas2"
# [1] 0.9577627
# [1] 0.1793104
# [1] 86.05493
# [1] 0
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 1   76 3832    0
# 2    2    0 3868
# [1] "N samples exclude - wrong sex: 0"
# 
# TRUE 
# 7778 
# [1] "N samples sex set to 0: 78"
# [1] "N cases and ctr excluded:"
# < table of extent 0 >

#   [1] "pittsburgh_gsa"
# [1] 0.7591599
# [1] 0.280888
# [1] 207.1452
# [1] 61.48009
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 1    3 1197   14
# 2    5   16 1546
# [1] "N samples exclude - wrong sex: 30"
# 
# TRUE 
# 2781 
# [1] "N samples sex set to 0: 8"
# [1] "N cases and ctr excluded:"
# 
# 1  2 
# 25  5 

# [1] "spain_gsa"
# [1] 0.9940536
# [1] 0.1613709
# [1] 46.4446
# [1] 9.606449
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 1   29 1883    0
# 2   27    0 1504
# [1] "N samples exclude - wrong sex: 0"
# 
# TRUE 
# 3443 
# [1] "N samples sex set to 0: 56"
# [1] "N cases and ctr excluded:"
# < table of extent 0 >

#   [1] "italy_gsa"
# [1] 0.2001846
# [1] 0.2099138
# [1] 49.03787
# [1] 23.5082
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 1  10 521  20
# 2  10  16 438
# [1] "N samples exclude - wrong sex: 36"
# 
# TRUE 
# 1016 
# [1] "N samples sex set to 0: 21"
# [1] "N cases and ctr excluded:"
# 
# 1  2 
# 15 21 

# [1] "kiel_austria_sibdcs_gsa"
# [1] 0.507751
# [1] 0.523615
# [1] 80.3725
# [1] 84.49217
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 0    1    1    0
# 1    7 7140  103
# 2   11  111 7286
# [1] "N samples exclude - wrong sex: 214"
# 
# TRUE 
# 14660 
# [1] "N samples sex set to 0: 19"
# [1] "N cases and ctr excluded:"
# 
# 1   2 
# 4 210 
# [1] "netherlands_gsa"
# [1] 0.3993712
# [1] 0.4947344
# [1] 78.45216
# [1] 54.31178
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 1   15 1917   15
# 2   46   12 2700
# [1] "N samples exclude - wrong sex: 27"
# 
# TRUE 
# 4705 
# [1] "N samples sex set to 0: 61"
# [1] "N cases and ctr excluded:"
# 
# 2 
# 27 

# [1] "slovenia_gsa"
# [1] -0.05972588
# [1] 0.4709448
# [1] 9.482556
# [1] 13.0298
# [1] 0.182867
# [1] 0.3570064
# [1] 14.61192
# [1] 9.944967
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 0   0   0   1
# 1   2  55   2
# 2   2   2 206
# [1] "N samples exclude - wrong sex: 4"
# 
# TRUE 
# 270 
# [1] "N samples sex set to 0: 4"
# [1] "N cases and ctr excluded:"
# 
# 1 2 
# 3 1 
# [1] "sweden_gsa"
# [1] 0.6104244
# [1] 0.4251331
# [1] 128.5972
# [1] 85.51756
# [1] "Recorded vs inferred sex:"
# 
# 1   2
# 1 814   7
# 2   6 587
# [1] "N samples exclude - wrong sex: 13"
# 
# TRUE 
# 1414 
# [1] "N samples sex set to 0: 0"
# [1] "N cases and ctr excluded:"
# 
# 1  2 
# 11  2 

# [1] "niddk_broad_gsa"
# [1] 0.6323922
# [1] 0.5288119
# [1] 115.352
# [1] 72.13633
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 1    1 2835   22
# 2   19   24 2614
# [1] "N samples exclude - wrong sex: 46"
# 
# TRUE 
# 5515 
# [1] "N samples sex set to 0: 20"
# [1] "N cases and ctr excluded:"
# 
# 1  2 
# 9 37 

# [1] "niddk_feinstein_gsa"
# [1] 0.5171672
# [1] 0.5058649
# [1] 191.4599
# [1] 166.4619
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 0    0    4    4
# 1    3 4103   55
# 2   14   47 4071
# [1] "N samples exclude - wrong sex: 102"
# 
# TRUE 
# 8301 
# [1] "N samples sex set to 0: 17"
# [1] "N cases and ctr excluded:"
# 
# 1  2 
# 9 93 
# [1] "basque_gsa"
# [1] 0.7735553
# [1] 0.1613679
# [1] 356.4772
# [1] 65.55902
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 0   3   2   0
# 1   7 924   0
# 2   3   0 577
# [1] "N samples exclude - wrong sex: 0"
# 
# TRUE 
# 1516 
# [1] "N samples sex set to 0: 13"
# [1] "N cases and ctr excluded:"
# < table of extent 0 >
#   [1] "lithuania_gsa"
# [1] 0.4008437
# [1] 0.5782721
# [1] 69.60423
# [1] 97.20851
# [1] "Recorded vs inferred sex:"
# 
# 1    2
# 1 1181   25
# 2   21 1050
# [1] "N samples exclude - wrong sex: 46"
# 
# TRUE 
# 2277 
# [1] "N samples sex set to 0: 0"
# [1] "N cases and ctr excluded:"
# 
# 1  2 
# 30 16 
# [1] "belgium_louis_gsa"
# [1] 0.5073549
# [1] 0.3764134
# [1] 80.19207
# [1] 46.98763
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 0   5   6  10
# 1   4 651   7
# 2   3   5 840
# [1] "N samples exclude - wrong sex: 12"
# 
# TRUE 
# 1531 
# [1] "N samples sex set to 0: 12"
# [1] "N cases and ctr excluded:"
# 
# 1 2 
# 5 7 
# [1] "belgium_franchimont_gsa"
# [1] 0.2476441
# [1] 0.4753855
# [1] 41.91967
# [1] 65.32962
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 0   0   1   0
# 1   0 645  21
# 2   2  10 853
# [1] "N samples exclude - wrong sex: 31"
# 
# TRUE 
# 1532 
# [1] "N samples sex set to 0: 2"
# [1] "N cases and ctr excluded:"
# 
# 1  2 
# 16 15 
# [1] "belgium_vermeire_gsa"
# [1] 0.572042
# [1] 0.3276447
# [1] 243.4374
# [1] 96.48579
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 1   13 1915    6
# 2   13    6 2061
# [1] "N samples exclude - wrong sex: 12"
# 
# TRUE 
# 4014 
# [1] "N samples sex set to 0: 26"
# [1] "N cases and ctr excluded:"
# 
# 1 2 
# 5 7 
# [1] "prism_nfe_gsa"
# [1] 0.8000809
# [1] 0.1741934
# [1] 70.56623
# [1] 6.557341
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 1   2 216   0
# 2   2   0 246
# [1] "N samples exclude - wrong sex: 0"
# 
# TRUE 
# 466 
# [1] "N samples sex set to 0: 4"
# [1] "N cases and ctr excluded:"
# < table of extent 0 >
#   [1] "prism_nfe_gwas"
# [1] 0.6230518
# [1] 0.7525303
# [1] "Recorded vs inferred sex:"
# 
# 1   2
# 1 398   5
# 2  10 449
# [1] "N samples exclude - wrong sex: 15"
# 
# TRUE 
# 862 
# [1] "N samples sex set to 0: 0"
# [1] "N cases and ctr excluded:"
# 
# 1  2 
# 1 14 
# [1] "finland_illugwas"
# [1] 0.525892
# [1] 0.8426228
# [1] "Recorded vs inferred sex:"
# 
# 1   2
# 1 297   4
# 2   6 150
# [1] "N samples exclude - wrong sex: 10"
# 
# TRUE 
# 457 
# [1] "N samples sex set to 0: 0"
# [1] "N cases and ctr excluded:"
# 
# 2 
# 10 
# [1] "german_affy6_old_gwas"
# [1] 0.9872924
# [1] 0.1841669
# [1] 25.84039
# [1] 0
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 1    5 1368    0
# 2    2    0 1452
# [1] "N samples exclude - wrong sex: 0"
# 
# TRUE 
# 2827 
# [1] "N samples sex set to 0: 7"
# [1] "N cases and ctr excluded:"
# < table of extent 0 >
#   [1] "norway_affy6_old_gwas"
# [1] 1
# [1] 0.2123022
# [1] 127.5701
# [1] 0
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 1   1 304   0
# 2   0   0 245
# [1] "N samples exclude - wrong sex: 0"
# 
# TRUE 
# 550 
# [1] "N samples sex set to 0: 1"
# [1] "N cases and ctr excluded:"
# < table of extent 0 >
#   [1] "belgium_inf1_old_gwas"
# [1] 0.9901173
# [1] 0.1505447
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 0   0  41  48
# 1  11 950   0
# 2   0   0 367
# [1] "N samples exclude - wrong sex: 0"
# 
# TRUE 
# 1417 
# [1] "N samples sex set to 0: 11"
# [1] "N cases and ctr excluded:"
# < table of extent 0 >
#   [1] "belgium_inf2_old_gwas"
# [1] 0.9906225
# [1] 0.1436435
# [1] "Recorded vs inferred sex:"
# 
# 0  1  2
# 0  0 79 36
# 1  1 58  0
# 2  0  0 98
# [1] "N samples exclude - wrong sex: 0"
# 
# TRUE 
# 272 
# [1] "N samples sex set to 0: 1"
# [1] "N cases and ctr excluded:"
# < table of extent 0 >
#   [1] "cedars_370k_old_gwas"
# [1] 0.780269
# [1] 0.4986897
# [1] 42.59272
# [1] 18.85818
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 1   0 324   1
# 2   1   2 280
# [1] "N samples exclude - wrong sex: 3"
# 
# TRUE 
# 608 
# [1] "N samples sex set to 0: 1"
# [1] "N cases and ctr excluded:"
# 
# 2 
# 3 
# [1] "cedars_610k_old_gwas"
# [1] 0.8161607
# [1] 0.4991328
# [1] 47.67975
# [1] 20.13728
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 1   0 465   1
# 2   3   3 421
# [1] "N samples exclude - wrong sex: 4"
# 
# TRUE 
# 893 
# [1] "N samples sex set to 0: 3"
# [1] "N cases and ctr excluded:"
# 
# 2 
# 4 
# [1] "cedars_omni_old_gwas"
# [1] 0.7955312
# [1] 0.3467193
# [1] 332.6458
# [1] 121.938
# [1] "Recorded vs inferred sex:"
# 
# 1   2
# 1 594   2
# 2   3 627
# [1] "N samples exclude - wrong sex: 5"
# 
# TRUE 
# 1226 
# [1] "N samples sex set to 0: 0"
# [1] "N cases and ctr excluded:"
# 
# 2 
# 5 
# [1] "mccauley_gsa"
# [1] 0.6574553
# [1] 0.4833173
# [1] 102.4786
# [1] 60.9409
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 1   3 400   1
# 2   0   5 379
# [1] "N samples exclude - wrong sex: 6"
# 
# TRUE 
# 788 
# [1] "N samples sex set to 0: 3"
# [1] "N cases and ctr excluded:"
# 
# 2 
# 6 
# [1] "ccfa_gsa"
# [1] 0.6487664
# [1] 0.3170274
# [1] 147.8335
# [1] 58.68251
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 1    0  998    7
# 2    2    4 1177
# [1] "N samples exclude - wrong sex: 11"
# 
# TRUE 
# 2188 
# [1] "N samples sex set to 0: 2"
# [1] "N cases and ctr excluded:"
# 
# 2 
# 11 
# [1] "cedars_gsa"
# [1] 0.6564962
# [1] 0.375563
# [1] 221.6593
# [1] 98.7857
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 1    6 1483    7
# 2    4    9 1587
# [1] "N samples exclude - wrong sex: 16"
# 
# TRUE 
# 3096 
# [1] "N samples sex set to 0: 10"
# [1] "N cases and ctr excluded:"
# 
# 1  2 
# 3 13 
# [1] "bernstein_gsa"
# [1] 0.2377314
# [1] 0.3346587
# [1] 99.64907
# [1] 108.8124
# [1] "Recorded vs inferred sex:"
# 
# 1   2
# 0  81 136
# 1 120   4
# 2   1 172
# [1] "N samples exclude - wrong sex: 5"
# 
# TRUE 
# 514 
# [1] "N samples sex set to 0: 0"
# [1] "N cases and ctr excluded:"
# 
# 2 
# 5 
# [1] "farkkila_gsa"
# [1] -1.014434
# [1] 1.819207
# [1] -9.081789
# [1] 23.88197
# [1] "Recorded vs inferred sex:"
# 
# 0  1  2
# 1  0 25  9
# 2  1  7 26
# [1] "N samples exclude - wrong sex: 16"
# 
# TRUE 
# 68 
# [1] "N samples sex set to 0: 1"
# [1] "N cases and ctr excluded:"
# 
# 2 
# 16 
# [1] "franchimont_gsa"
# [1] 0.05717396
# [1] 0.8139173
# [1] 18.90075
# [1] 236.6331
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 0    0   34   37
# 1    1 1208   65
# 2    1   56 1387
# [1] "N samples exclude - wrong sex: 121"
# 
# TRUE 
# 2789 
# [1] "N samples sex set to 0: 2"
# [1] "N cases and ctr excluded:"
# 
# 1   2 
# 105  16 
# [1] "franke_gsa"
# [1] 0.2173622
# [1] 0.8847324
# [1] 13.7269
# [1] 51.35286
# [1] "Recorded vs inferred sex:"
# 
# 1   2
# 1 490  18
# 2  18 359
# [1] "N samples exclude - wrong sex: 36"
# 
# TRUE 
# 885 
# [1] "N samples sex set to 0: 0"
# [1] "N cases and ctr excluded:"
# 
# 1  2 
# 34  2 
# [1] "helmsley_prism_gsa"
# [1] 0.4674508
# [1] 0.5475899
# [1] "Recorded vs inferred sex:"
# 
# 1   2
# 0   0   2
# 1 375   7
# 2   6 390
# [1] "N samples exclude - wrong sex: 13"
# 
# TRUE 
# 780 
# [1] "N samples sex set to 0: 0"
# [1] "N cases and ctr excluded:"
# 
# 1  2 
# 2 11 
# [1] "helmsley_xavier_prism_gsa"
# [1] 0.6193141
# [1] 0.3601657
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 0   0   8  14
# 1   1 574   5
# 2   0   4 691
# [1] "N samples exclude - wrong sex: 9"
# 
# TRUE 
# 1297 
# [1] "N samples sex set to 0: 1"
# [1] "N cases and ctr excluded:"
# 
# 1 2 
# 1 8 
# [1] "hyams_protect_gsa"
# [1] 0.6200814
# [1] 0.3711209
# [1] 112.6361
# [1] 54.60204
# [1] "Recorded vs inferred sex:"
# 
# 1   2
# 1 210   2
# 2   1 205
# [1] "N samples exclude - wrong sex: 3"
# 
# TRUE 
# 418 
# [1] "N samples sex set to 0: 0"
# [1] "N cases and ctr excluded:"
# 
# 2 
# 3 
# [1] "lewis_sparc_gsa"
# [1] 0.6355259
# [1] 0.3577071
# [1] 154.1853
# [1] 71.53864
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 1    0 1281   10
# 2    1    8 1559
# [1] "N samples exclude - wrong sex: 18"
# 
# TRUE 
# 2859 
# [1] "N samples sex set to 0: 1"
# [1] "N cases and ctr excluded:"
# 
# 2 
# 18 
# [1] "mccauley_new_gsa"
# [1] 0.4810127
# [1] 0.4574291
# [1] 92.97643
# [1] 79.17084
# [1] "Recorded vs inferred sex:"
# 
# 1   2
# 1 832  13
# 2   8 775
# [1] "N samples exclude - wrong sex: 21"
# 
# TRUE 
# 1628 
# [1] "N samples sex set to 0: 0"
# [1] "N cases and ctr excluded:"
# 
# 1  2 
# 1 20 
# [1] "mcgovern_gsa"
# [1] 0.716694
# [1] 0.3634103
# [1] 212.5066
# [1] 85.18696
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 0    0   12    6
# 1    4 2949   12
# 2    8   14 3045
# [1] "N samples exclude - wrong sex: 26"
# 
# TRUE 
# 6050 
# [1] "N samples sex set to 0: 12"
# [1] "N cases and ctr excluded:"
# 
# 1  2 
# 5 21 
# [1] "moayyedi_imagine_gsa"
# [1] 0.8033076
# [1] 0.3317031
# [1] 208.5282
# [1] 43.59448
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 1   0 489   1
# 2   5   1 632
# [1] "N samples exclude - wrong sex: 2"
# 
# TRUE 
# 1128 
# [1] "N samples sex set to 0: 5"
# [1] "N cases and ctr excluded:"
# 
# 2 
# 2 
# [1] "newberry_share_gsa"
# [1] 0.7640911
# [1] 0.344738
# [1] 130.304
# [1] 8.655238
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 1   2 366   1
# 2   4   0 492
# [1] "N samples exclude - wrong sex: 1"
# 
# TRUE 
# 865 
# [1] "N samples sex set to 0: 6"
# [1] "N cases and ctr excluded:"
# 
# 2 
# 1 
# [1] "niddk_cho_gsa"
# [1] 0.7976913
# [1] 0.3691771
# [1] 180.209
# [1] 55.35015
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 0   0   2   2
# 1   1 902   2
# 2   3   3 849
# [1] "N samples exclude - wrong sex: 5"
# 
# TRUE 
# 1764 
# [1] "N samples sex set to 0: 4"
# [1] "N cases and ctr excluded:"
# 
# 1 2 
# 2 3 
# [1] "niddk_duerr_gsa"
# [1] 0.7920603
# [1] 0.3878923
# [1] 185.6476
# [1] 32.93121
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 0   0   2   1
# 1   0 967   2
# 2  10   1 961
# [1] "N samples exclude - wrong sex: 3"
# 
# TRUE 
# 1944 
# [1] "N samples sex set to 0: 10"
# [1] "N cases and ctr excluded:"
# 
# 1 2 
# 1 2 
# [1] "niddk_rioux_gsa"
# [1] 0.9345547
# [1] 0.454773
# [1] 194.0726
# [1] 78.7557
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 0   0   2   0
# 1   3 417   0
# 2   1   4 492
# [1] "N samples exclude - wrong sex: 4"
# 
# TRUE 
# 919 
# [1] "N samples sex set to 0: 4"
# [1] "N cases and ctr excluded:"
# 
# 1 2 
# 1 3 
# [1] "niddk_silverberg_gsa"
# [1] 0.8092606
# [1] 0.3969423
# [1] 208.7843
# [1] 76.02448
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 0    0    3    2
# 1    3 1207    2
# 2    5    6 1143
# [1] "N samples exclude - wrong sex: 8"
# 
# TRUE 
# 2371 
# [1] "N samples sex set to 0: 8"
# [1] "N cases and ctr excluded:"
# 
# 1 2 
# 1 7 
# [1] "palotie_hus_gsa"
# [1] 0.7745049
# [1] 0.2555769
# [1] 182.2929
# [1] 37.11109
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 1   6 414   0
# 2   2   1 455
# [1] "N samples exclude - wrong sex: 1"
# 
# TRUE 
# 878 
# [1] "N samples sex set to 0: 8"
# [1] "N cases and ctr excluded:"
# 
# 2 
# 1 
# [1] "pekow_share_gsa"
# [1] 0.6587652
# [1] 0.3093671
# [1] 54.91352
# [1] 19.20456
# [1] "Recorded vs inferred sex:"
# 
# 1   2
# 1 315   2
# 2   1 316
# [1] "N samples exclude - wrong sex: 3"
# 
# TRUE 
# 634 
# [1] "N samples sex set to 0: 0"
# [1] "N cases and ctr excluded:"
# 
# 2 
# 3 
# [1] "rioux_igenomed_gsa"
# [1] 0.2217374
# [1] 0.6033419
# [1] 68.24649
# [1] 156.7637
# [1] "Recorded vs inferred sex:"
# 
# 1  2
# 1 84  3
# 2  2 93
# [1] "N samples exclude - wrong sex: 5"
# 
# TRUE 
# 182 
# [1] "N samples sex set to 0: 0"
# [1] "N cases and ctr excluded:"
# 
# 2 
# 5 
# [1] "sands_msccr_gsa"
# [1] -1.064553
# [1] 1.693165
# [1] -121.7135
# [1] 227.0552
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 0   0   0   1
# 1   1 410 354
# 2   0 358 306
# [1] "N samples exclude - wrong sex: 712"
# 
# TRUE 
# 1430 
# [1] "N samples sex set to 0: 1"
# [1] "N cases and ctr excluded:"
# 
# 1   2 
# 143 569 
# [1] "stampfer_gsa"
# [1] 0.9363852
# [1] 0.2015263
# [1] 283.3617
# [1] 36.71733
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 1    3  236    0
# 2    2    1 1235
# [1] "N samples exclude - wrong sex: 1"
# 
# TRUE 
# 1477 
# [1] "N samples sex set to 0: 5"
# [1] "N cases and ctr excluded:"
# 
# 1 
# 1 
# [1] "vermeire_gsa"
# [1] 0.7365233
# [1] 0.3132987
# [1] 273.0142
# [1] 85.43935
# [1] "Recorded vs inferred sex:"
# 
# 0    1    2
# 1    7 2267    6
# 2   10    7 2416
# [1] "N samples exclude - wrong sex: 13"
# 
# TRUE 
# 4713 
# [1] "N samples sex set to 0: 17"
# [1] "N cases and ctr excluded:"
# 
# 1 2 
# 5 8 
# [1] "weersma_gsa"
# [1] 0.769836
# [1] 0.2734049
# [1] 175.2706
# [1] 9.552542
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 1   2 353   0
# 2   2   0 352
# [1] "N samples exclude - wrong sex: 0"
# 
# TRUE 
# 709 
# [1] "N samples sex set to 0: 4"
# [1] "N cases and ctr excluded:"
# < table of extent 0 >
#   [1] "xavier_prism_gsa"
# [1] 0.6690088
# [1] 0.4262093
# [1] 62.18069
# [1] 28.97252
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 1   0 341   2
# 2   1   2 346
# [1] "N samples exclude - wrong sex: 4"
# 
# TRUE 
# 692 
# [1] "N samples sex set to 0: 1"
# [1] "N cases and ctr excluded:"
# 
# 1 2 
# 1 3 
# [1] "xavier_share_gsa"
# [1] 0.9556248
# [1] 0.2259828
# [1] 184.1904
# [1] 9.051569
# [1] "Recorded vs inferred sex:"
# 
# 0   1   2
# 1   3 345   0
# 2   4   0 344
# [1] "N samples exclude - wrong sex: 0"
# 
# TRUE 
# 696 
# [1] "N samples sex set to 0: 7"
# [1] "N cases and ctr excluded:"
# < table of extent 0 >
#   > 
#   > write.table(dat,"~/tmp_plots/summary_sample_gender_sex_discrepancies.tsv",col.names=T,row.names=F,quote=F,sep="\t")
# > dat
#                      studies N_samples N_samples_set_sex_0
# 1        australia_omniexome      1310                   1
# 2                      gwas1      4684                   2
# 3                      gwas2      7778                  78
# 4             pittsburgh_gsa      2781                   8
# 5                  spain_gsa      3443                  56
# 6                  italy_gsa      1016                  21
# 7    kiel_austria_sibdcs_gsa     14660                  19
# 8            netherlands_gsa      4705                  61
# 9               slovenia_gsa       270                   4
# 10                sweden_gsa      1414                   0
# 11           niddk_broad_gsa      5515                  20
# 12       niddk_feinstein_gsa      8301                  17
# 13                basque_gsa      1516                  13
# 14             lithuania_gsa      2277                   0
# 15         belgium_louis_gsa      1531                  12
# 16   belgium_franchimont_gsa      1532                   2
# 17      belgium_vermeire_gsa      4014                  26
# 18             prism_nfe_gsa       466                   4
# 19            prism_nfe_gwas       862                   0
# 20          finland_illugwas       457                   0
# 21     german_affy6_old_gwas      2827                   7
# 22     norway_affy6_old_gwas       550                   1
# 23     belgium_inf1_old_gwas      1417                  11
# 24     belgium_inf2_old_gwas       272                   1
# 25      cedars_370k_old_gwas       608                   1
# 26      cedars_610k_old_gwas       893                   3
# 27      cedars_omni_old_gwas      1226                   0
# 28              mccauley_gsa       788                   3
# 29                  ccfa_gsa      2188                   2
# 30                cedars_gsa      3096                  10
# 31             bernstein_gsa       514                   0
# 32              farkkila_gsa        68                   1
# 33           franchimont_gsa      2789                   2
# 34                franke_gsa       885                   0
# 35        helmsley_prism_gsa       780                   0
# 36 helmsley_xavier_prism_gsa      1297                   1
# 37         hyams_protect_gsa       418                   0
# 38           lewis_sparc_gsa      2859                   1
# 39          mccauley_new_gsa      1628                   0
# 40              mcgovern_gsa      6050                  12
# 41      moayyedi_imagine_gsa      1128                   5
# 42        newberry_share_gsa       865                   6
# 43             niddk_cho_gsa      1764                   4
# 44           niddk_duerr_gsa      1944                  10
# 45           niddk_rioux_gsa       919                   4
# 46      niddk_silverberg_gsa      2371                   8
# 47           palotie_hus_gsa       878                   8
# 48           pekow_share_gsa       634                   0
# 49        rioux_igenomed_gsa       182                   0
# 50           sands_msccr_gsa      1430                   1
# 51              stampfer_gsa      1477                   5
# 52              vermeire_gsa      4713                  17
# 53               weersma_gsa       709                   4
# 54          xavier_prism_gsa       692                   1
# 55          xavier_share_gsa       696                   7
# N_samples_discordant_gender
# 1                           14
# 2                            4
# 3                            0
# 4                           30
# 5                            0
# 6                           36
# 7                          214
# 8                           27
# 9                            4
# 10                          13
# 11                          46
# 12                         102
# 13                           0
# 14                          46
# 15                          12
# 16                          31
# 17                          12
# 18                           0
# 19                          15
# 20                          10
# 21                           0
# 22                           0
# 23                           0
# 24                           0
# 25                           3
# 26                           4
# 27                           5
# 28                           6
# 29                          11
# 30                          16
# 31                           5
# 32                          16
# 33                         121
# 34                          36
# 35                          13
# 36                           9
# 37                           3
# 38                          18
# 39                          21
# 40                          26
# 41                           2
# 42                           1
# 43                           5
# 44                           3
# 45                           4
# 46                           8
# 47                           1
# 48                           3
# 49                           5
# 50                         712
# 51                           1
# 52                          13
# 53                           0
# 54                           4
# 55                           0
# > dat
# studies N_samples N_samples_set_sex_0
# 1        australia_omniexome      1310                   1
# 2                      gwas1      4684                   2
# 3                      gwas2      7778                  78
# 4             pittsburgh_gsa      2781                   8
# 5                  spain_gsa      3443                  56
# 6                  italy_gsa      1016                  21
# 7    kiel_austria_sibdcs_gsa     14660                  19
# 8            netherlands_gsa      4705                  61
# 9               slovenia_gsa       270                   4
# 10                sweden_gsa      1414                   0
# 11           niddk_broad_gsa      5515                  20
# 12       niddk_feinstein_gsa      8301                  17
# 13                basque_gsa      1516                  13
# 14             lithuania_gsa      2277                   0
# 15         belgium_louis_gsa      1531                  12
# 16   belgium_franchimont_gsa      1532                   2
# 17      belgium_vermeire_gsa      4014                  26
# 18             prism_nfe_gsa       466                   4
# 19            prism_nfe_gwas       862                   0
# 20          finland_illugwas       457                   0
# 21     german_affy6_old_gwas      2827                   7
# 22     norway_affy6_old_gwas       550                   1
# 23     belgium_inf1_old_gwas      1417                  11
# 24     belgium_inf2_old_gwas       272                   1
# 25      cedars_370k_old_gwas       608                   1
# 26      cedars_610k_old_gwas       893                   3
# 27      cedars_omni_old_gwas      1226                   0
# 28              mccauley_gsa       788                   3
# 29                  ccfa_gsa      2188                   2
# 30                cedars_gsa      3096                  10
# 31             bernstein_gsa       514                   0
# 32              farkkila_gsa        68                   1
# 33           franchimont_gsa      2789                   2
# 34                franke_gsa       885                   0
# 35        helmsley_prism_gsa       780                   0
# 36 helmsley_xavier_prism_gsa      1297                   1
# 37         hyams_protect_gsa       418                   0
# 38           lewis_sparc_gsa      2859                   1
# 39          mccauley_new_gsa      1628                   0
# 40              mcgovern_gsa      6050                  12
# 41      moayyedi_imagine_gsa      1128                   5
# 42        newberry_share_gsa       865                   6
# 43             niddk_cho_gsa      1764                   4
# 44           niddk_duerr_gsa      1944                  10
# 45           niddk_rioux_gsa       919                   4
# 46      niddk_silverberg_gsa      2371                   8
# 47           palotie_hus_gsa       878                   8
# 48           pekow_share_gsa       634                   0
# 49        rioux_igenomed_gsa       182                   0
# 50           sands_msccr_gsa      1430                   1
# 51              stampfer_gsa      1477                   5
# 52              vermeire_gsa      4713                  17
# 53               weersma_gsa       709                   4
# 54          xavier_prism_gsa       692                   1
# 55          xavier_share_gsa       696                   7
# N_samples_discordant_gender
# 1                           14
# 2                            4
# 3                            0
# 4                           30
# 5                            0
# 6                           36
# 7                          214
# 8                           27
# 9                            4
# 10                          13
# 11                          46
# 12                         102
# 13                           0
# 14                          46
# 15                          12
# 16                          31
# 17                          12
# 18                           0
# 19                          15
# 20                          10
# 21                           0
# 22                           0
# 23                           0
# 24                           0
# 25                           3
# 26                           4
# 27                           5
# 28                           6
# 29                          11
# 30                          16
# 31                           5
# 32                          16
# 33                         121
# 34                          36
# 35                          13
# 36                           9
# 37                           3
# 38                          18
# 39                          21
# 40                          26
# 41                           2
# 42                           1
# 43                           5
# 44                           3
# 45                           4
# 46                           8
# 47                           1
# 48                           3
# 49                           5
# 50                         712
# 51                           1
# 52                          13
# 53                           0
# 54                           4
# 55                           0





############ 

for i in ${studies[@]}
do echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/list_samples_wrong_gender
done
# update this file

############

# for the studies we have an answer for (regarding gender updates - see files_to_be_sent_to_IIBDGC_PI.R)

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

for (j in 1:length(studies)) {
  file_2020<-paste(path,"pre_imputation/QC/id_store/id_matches_2020/files_iibdgc/from_suppliers/",studies[j],"_list_samples_sex_discrepancy_with_supplierIds_updated_gender_from_suppliers.tsv",sep="")
  file_2022<-paste(path,"pre_imputation/QC/id_store/id_matches_2022/files_iibdgc/from_suppliers/",studies[j],"_list_samples_sex_discrepancy_with_supplierIds_updated_gender_from_suppliers.tsv",sep="")
  
  if(file.exists(file_2020)) {
    print(studies[j])
    tmp<-read.table(paste(path,"pre_imputation/QC/id_store/id_matches_2020/files_iibdgc/from_suppliers/",studies[j],"_list_samples_sex_discrepancy_with_supplierIds_updated_gender_from_suppliers.tsv",sep=""),head=T,sep="\t")
    print(table(tmp$action))
    tmp<-tmp[which(tmp$action=="exclude"),]
    write.table(tmp[,c(1,1)],paste(path,"pre_imputation/QC/",studies[j],"/list_samples_wrong_gender",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  } else if(file.exists(file_2022)){
    print(studies[j])
    tmp<-read.table(paste(path,"pre_imputation/QC/id_store/id_matches_2022/files_iibdgc/from_suppliers/",studies[j],"_list_samples_sex_discrepancy_with_supplierIds_updated_gender_from_suppliers.tsv",sep=""),head=T,sep="\t")
    print(table(tmp$action))
    tmp<-tmp[which(tmp$action=="exclude"),]
    write.table(tmp[,c(1,1)],paste(path,"pre_imputation/QC/",studies[j],"/list_samples_wrong_gender",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    
  }else{
    print(paste("No list samples with wrong gender for study: ",studies[j],sep=""))
    tmp<-as.data.frame(matrix(ncol=2,nrow=0))
    colnames(tmp)<-c("FID","IID")
    write.table(tmp,paste(path,"pre_imputation/QC/",studies[j],"/list_samples_wrong_gender",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  }
  rm(tmp)
}

##############

# [1] "all_hce"
# exclude 
# 347 

# [1] "niddk_old_gwas"
# exclude 
# 29 


# [1] "australia_omniexome"
# exclude 
# 4 

# [1] "gwas1"
# exclude 
# 4 

# [1] "gwas2"
# [1] "No list samples with wrong gender for study: gwas2"

# [1] "pittsburgh_gsa"
# exclude 
# 30 

# [1] "spain_gsa"
# [1] "No list samples with wrong gender for study: spain_gsa"

# [1] "italy_gsa"
# exclude 
# 36 

# [1] "kiel_austria_sibdcs_gsa"
# exclude updated_recorded_sex_keep 
# 10                       204 

# [1] "netherlands_gsa"
# exclude 
# 27 

# [1] "slovenia_gsa"
# updated_recorded_sex_keep 
# 4 

# [1] "sweden_gsa"
# exclude 
# 13 

# [1] "niddk_broad_gsa"
# exclude 
# 46 

# [1] "niddk_feinstein_gsa"
# exclude 
# 102 

# [1] "basque_gsa"
# [1] "No list samples with wrong gender for study: basque_gsa"

# [1] "lithuania_gsa"
# exclude 
# 46 

# [1] "belgium_louis_gsa"
# exclude 
# 12 

# [1] "belgium_franchimont_gsa"
# exclude 
# 31 

# [1] "belgium_vermeire_gsa"
# exclude 
# 12 

# [1] "prism_nfe_gsa"
# [1] "No list samples with wrong gender for study: prism_nfe_gsa"

# [1] "prism_nfe_gwas"
# exclude 
# 15 

# [1] "finland_illugwas"
# exclude 
# 10 

# [1] "german_affy6_old_gwas"
# [1] "No list samples with wrong gender for study: german_affy6_old_gwas"
# [1] "norway_affy6_old_gwas"
# [1] "No list samples with wrong gender for study: norway_affy6_old_gwas"
# [1] "belgium_inf1_old_gwas"
# [1] "No list samples with wrong gender for study: belgium_inf1_old_gwas"
# [1] "belgium_inf2_old_gwas"
# [1] "No list samples with wrong gender for study: belgium_inf1_old_gwas"

# [1] "cedars_370k_old_gwas"
# exclude 
# 3 

# [1] "cedars_610k_old_gwas"
# exclude 
# 4 

# [1] "cedars_omni_old_gwas"
# exclude 
# 5 

# [1] "mccauley_gsa"
# exclude 
# 6 

# [1] "ccfa_gsa"
# exclude 
# 11 

# [1] "cedars_gsa"
# exclude 
# 16 

# [1] "bernstein_gsa"
# exclude updated_recorded_sex_keep 
# 3                         2 

# [1] "farkkila_gsa"
# updated_recorded_sex_keep 
# 16 

# [1] "franchimont_gsa"
# [1] "No list samples with wrong gender for study: franchimont_gsa"

# [1] "franke_gsa"
# exclude updated_recorded_sex_keep 
# 18                        18 

# [1] "helmsley_prism_gsa"
# exclude updated_recorded_sex_keep 
# 12                         1 

# [1] "helmsley_xavier_prism_gsa"
# exclude updated_recorded_sex_keep 
# 4                         5 

# [1] "hyams_protect_gsa"
# [1] "No list samples with wrong gender for study: hyams_protect_gsa"

# [1] "lewis_sparc_gsa"
# [1] "No list samples with wrong gender for study: lewis_sparc_gsa"

# [1] "mccauley_new_gsa"
# exclude updated_recorded_sex_keep 
# 9                        12 

# [1] "mcgovern_gsa"
# exclude 
# 27 

# [1] "moayyedi_imagine_gsa"
# [1] "No list samples with wrong gender for study: moayyedi_imagine_gsa"

# [1] "newberry_share_gsa"
# [1] "No list samples with wrong gender for study: newberry_share_gsa"

# [1] "niddk_cho_gsa"
# [1] "No list samples with wrong gender for study: niddk_cho_gsa"

# [1] "niddk_duerr_gsa"
# [1] "No list samples with wrong gender for study: niddk_duerr_gsa"

# [1] "niddk_rioux_gsa"
# [1] "No list samples with wrong gender for study: niddk_rioux_gsa"

# [1] "niddk_silverberg_gsa"
# [1] "No list samples with wrong gender for study: niddk_silverberg_gsa"

# [1] "palotie_hus_gsa"
# updated_recorded_sex_keep 
# 1 
# [1] "pekow_share_gsa"
# [1] "No list samples with wrong gender for study: pekow_share_gsa"

# [1] "rioux_igenomed_gsa"
# [1] "No list samples with wrong gender for study: rioux_igenomed_gsa"

# [1] "sands_msccr_gsa"
# [1] "No list samples with wrong gender for study: sands_msccr_gsa"

# [1] "stampfer_gsa"
# exclude 
# 1 

# [1] "vermeire_gsa"
# exclude 
# 7 
# updated_recorded_sex_keep 
# 5 
# updated_recorded_sex_keep|sample_switch_should_be_ID1902_MALE 
# 1 

# [1] "weersma_gsa"
# [1] "No list samples with wrong gender for study: weersma_gsa"

# [1] "xavier_prism_gsa"
# exclude updated_recorded_sex_keep 
# 3                         1 

# [1] "xavier_share_gsa"
# [1] "No list samples with wrong gender for study: xavier_share_gsa"


##############

# do not remove (YET) any of the wrong samples - just set to 0 but keep the list of samples to exclude:
# --remove ${path_gwas}pre_imputation/QC/niddk_cho_gsa/list_samples_wrong_gender \

path_gwas=/path/to/ibdgwas/IIBDGC/
studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

MEM=500

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_19_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_19_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip.bed \
--bim ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip.bim \
--fam ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_edited.fam \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/${i}/list_samples_wrong_gender \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck"
done

# do as well for "swedish_uc_old_gwas"
i=swedish_uc_old_gwas
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_19_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_19_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip.bed \
--bim ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip.bim \
--fam ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_edited.fam \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck"



for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_19_${i} | grep -E "completed"
done

for i in ${studies[@]}
do echo ${i} \
&& wc -l ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck.fam \
&& wc -l ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck.bim
done

################

# all_hce
# 24004 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/all_hce/all_hce_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 509718 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/all_hce/all_hce_hg19_noind_posstr_nodup_flip_sexcheck.bim
# niddk_old_gwas
# 2815 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 314041 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.bim
# australia_omniexome
# 1306 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/australia_omniexome_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 788162 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/australia_omniexome_hg19_noind_posstr_nodup_flip_sexcheck.bim
# gwas1
# 4680 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/gwas1_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 456448 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/gwas1_hg19_noind_posstr_nodup_flip_sexcheck.bim
# gwas2
# 7778 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas2/gwas2_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 791554 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas2/gwas2_hg19_noind_posstr_nodup_flip_sexcheck.bim
# pittsburgh_gsa
# 2751 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 944977 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# spain_gsa
# 3443 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/spain_gsa/spain_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 585661 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/spain_gsa/spain_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# italy_gsa
# 980 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/italy_gsa/italy_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 584948 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/italy_gsa/italy_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# kiel_austria_sibdcs_gsa
# 14650 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 661984 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# netherlands_gsa
# 4678 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 684232 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# slovenia_gsa
# 270 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 557319 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# sweden_gsa
# 1401 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/sweden_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 574824 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/sweden_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# niddk_broad_gsa
# 5469 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 621749 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# niddk_feinstein_gsa
# 8199 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 644768 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# basque_gsa
# 1516 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/basque_gsa/basque_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 577964 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/basque_gsa/basque_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# lithuania_gsa
# 2231 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lithuania_gsa/lithuania_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 579645 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lithuania_gsa/lithuania_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# belgium_louis_gsa
# 1519 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 605309 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# belgium_franchimont_gsa
# 1501 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 584798 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# belgium_vermeire_gsa
# 4002 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 674160 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# prism_nfe_gsa
# 466 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 573340 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# prism_nfe_gwas
# 847 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 243968 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_hg19_noind_posstr_nodup_flip_sexcheck.bim
# finland_illugwas
# 447 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/finland_illugwas_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 239229 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/finland_illugwas_hg19_noind_posstr_nodup_flip_sexcheck.bim
# german_affy6_old_gwas
# 2827 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 885198 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.bim
# norway_affy6_old_gwas
# 550 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 841006 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.bim
# belgium_inf1_old_gwas
# 1417 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 302700 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.bim
# belgium_inf2_old_gwas
# 272 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 290730 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.bim
# cedars_370k_old_gwas
# 605 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 339945 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.bim
# cedars_610k_old_gwas
# 889 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 586293 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.bim
# cedars_omni_old_gwas
# 1221 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 721646 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.bim
# swedish_uc_old_gwas
# 1264 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 297004 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck.bim
# mccauley_gsa
# 782 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 570755 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# ccfa_gsa
# 2177 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 610890 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# cedars_gsa
# 3080 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_gsa/cedars_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 668857 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_gsa/cedars_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# bernstein_gsa
# 511 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/bernstein_gsa/bernstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 565340 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/bernstein_gsa/bernstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# farkkila_gsa
# 68 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 527922 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# franchimont_gsa
# 2789 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franchimont_gsa/franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 670502 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franchimont_gsa/franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# franke_gsa
# 867 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franke_gsa/franke_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 570698 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franke_gsa/franke_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# helmsley_prism_gsa
# 768 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 243969 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# helmsley_xavier_prism_gsa
# 1293 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 243943 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# hyams_protect_gsa
# 418 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 580493 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# lewis_sparc_gsa
# 2859 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 619419 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# mccauley_new_gsa
# 1619 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 590905 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# mcgovern_gsa
# 6023 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 651518 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# moayyedi_imagine_gsa
# 1128 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 615066 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# newberry_share_gsa
# 865 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 592265 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# niddk_cho_gsa
# 1764 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 608310 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# niddk_duerr_gsa
# 1944 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 595827 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# niddk_rioux_gsa
# 919 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 589254 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# niddk_silverberg_gsa
# 2371 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 626387 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# palotie_hus_gsa
# 878 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 565002 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# pekow_share_gsa
# 634 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 580800 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# rioux_igenomed_gsa
# 182 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 551393 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# sands_msccr_gsa
# 1430 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 603749 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# stampfer_gsa
# 1476 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 671985 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# vermeire_gsa
# 4706 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/vermeire_gsa/vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 674734 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/vermeire_gsa/vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# weersma_gsa
# 709 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/weersma_gsa/weersma_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 573475 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/weersma_gsa/weersma_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# xavier_prism_gsa
# 689 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 583917 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim
# xavier_share_gsa
# 696 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck.fam
# 584504 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck.bim

################

####################################################################################
# 7.4 SET TO MISSING MALE HET CHRX CALLS  - Finish setting hh calls to missing

MEM=500

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_20_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_20_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck \
--allow-no-sex \
--set-hh-missing \
--make-bed --out ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh"
done


for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_20_${i} | grep -E "completed"
done

########################################################################################################################################

