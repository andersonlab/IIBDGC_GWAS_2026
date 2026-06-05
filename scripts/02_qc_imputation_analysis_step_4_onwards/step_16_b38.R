# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
######################
# 16.- UPDATE TO B38 #
######################

# NOTE: data in build 37, needs to be liftover first:

### SAMPLES WITH DUPLICATES

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=1500
# MEM=3000 - for all_hce and kiel_austria_sibdcs_gsa

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_52_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_52_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het \
--recode tab --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_52_${i} | grep -E "completed"
done

# ${path_gwas}previous_qced_b38/liftover/liftOverx
# usage:
# liftOver oldFile map.chain newFile unMapped

#### convert map into bed:
###### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

length(cohorts)
# [1] 58


for (j in 1:length(cohorts)) {
  
  print(cohorts[j])
  
  map<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.map",sep=""),head=F)
  
  map$chr<-paste("chr",map$V1,sep="")
  map$chr[which(map$V1=="23")]<-"chrX"
  map$chr[which(map$V1=="25")]<-"chrX"
  map$chr[which(map$V1=="24")]<-"chrY"
  map$chromStart<-map$V4
  map$chromEnd<-map$V4+1
  
  map$chromStart<-format(map$chromStart, scientific=F)
  map$chromEnd<-format(map$chromEnd, scientific=F)
  
  print(dim(map))
  
  write.table(map[,c(5:7,2)],paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_hg19_postqc.bed",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

}

############
# [1] "australia_omniexome"
# [1] 739312      7
# [1] "gwas1"
# [1] 453416      7
# [1] "gwas2"
# [1] 773228      7
# [1] "pittsburgh_gsa"
# [1] 897842      7
# [1] "spain_gsa"
# [1] 575804      7
# [1] "italy_gsa"
# [1] 573850      7
# [1] "netherlands_gsa"
# [1] 666191      7
# [1] "slovenia_gsa"
# [1] 529608      7
# [1] "sweden_gsa"
# [1] 558607      7
# [1] "niddk_broad_gsa"
# [1] 590626      7
# [1] "niddk_feinstein_gsa"
# [1] 603214      7
# [1] "basque_gsa"
# [1] 567401      7
# [1] "prism_nfe_gsa"
# [1] 559017      7
# [1] "lithuania_gsa"
# [1] 554135      7
# [1] "belgium_louis_gsa"
# [1] 591948      7
# [1] "belgium_franchimont_gsa"
# [1] 561994      7
# [1] "belgium_vermeire_gsa"
# [1] 653147      7
# [1] "prism_nfe_gwas"
# [1] 242899      7
# [1] "finland_illugwas"
# [1] 237055      7
# [1] "norway_affy6_old_gwas"
# [1] 725899      7
# [1] "belgium_inf1_old_gwas"
# [1] 301714      7
# [1] "belgium_inf2_old_gwas"
# [1] 290191      7
# [1] "niddk_old_gwas"
# [1] 300166      7
# [1] "cedars_370k_old_gwas"
# [1] 338717      7
# [1] "cedars_610k_old_gwas"
# [1] 583907      7
# [1] "cedars_omni_old_gwas"
# [1] 719139      7
# [1] "swedish_uc_old_gwas"
# [1] 296648      7
# [1] "mccauley_gsa"
# [1] 556682      7
# [1] "ccfa_gsa"
# [1] 599981      7
# [1] "cedars_gsa"
# [1] 646302      7
# [1] "bernstein_gsa"
# [1] 553944      7
# [1] "farkkila_gsa"
# [1] 511512      7
# [1] "franchimont_gsa"
# [1] 648211      7
# [1] "franke_gsa"
# [1] 557496      7
# [1] "helmsley_prism_gsa"
# [1] 243103      7
# [1] "helmsley_xavier_prism_gsa"
# [1] 243082      7
# [1] "hyams_protect_gsa"
# [1] 569468      7
# [1] "lewis_sparc_gsa"
# [1] 607715      7
# [1] "mccauley_new_gsa"
# [1] 578685      7
# [1] "moayyedi_imagine_gsa"
# [1] 608003      7
# [1] "newberry_share_gsa"
# [1] 583694      7
# [1] "niddk_cho_gsa"
# [1] 593576      7
# [1] "niddk_duerr_gsa"
# [1] 581142      7
# [1] "niddk_rioux_gsa"
# [1] 574971      7
# [1] "niddk_silverberg_gsa"
# [1] 611240      7
# [1] "palotie_hus_gsa"
# [1] 553667      7
# [1] "pekow_share_gsa"
# [1] 574060      7
# [1] "rioux_igenomed_gsa"
# [1] 543342      7
# [1] "sands_msccr_gsa"
# [1] 595558      7
# [1] "stampfer_gsa"
# [1] 643686      7
# [1] "vermeire_gsa"
# [1] 653343      7
# [1] "weersma_gsa"
# [1] 566401      7
# [1] "xavier_prism_gsa"
# [1] 569740      7
# [1] "xavier_share_gsa"
# [1] 568510      7
# [1] "all_hce"
# [1] 434598      7
# [1] "kiel_austria_sibdcs_gsa"
# [1] 564511      7
# [1] "german_affy6_old_gwas"
# [1] 783439      7
# [1] "mcgovern_gsa"
# [1] 630955      7

############


#### lift positions

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=200

for i in ${studies[@]}
do 
bsub -J"lf" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_liftover_2_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_liftover_2_${i} \
"${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}pre_imputation/QC/${i}/${i}_hg19_postqc.bed \
${path_gwas}previous_qced_b38/liftover/hg19ToHg38.over.chain.gz \
${path_gwas}pre_imputation/QC/${i}/${i}_hg19_postqc_lifted_hg38 \
${path_gwas}pre_imputation/QC/${i}/${i}_hg19_postqc_no_lifted_hg38"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_liftover_2_${i} | grep -E "completed"
done

for i in ${studies[@]}
do 
cut -f 4 ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_postqc_no_lifted_hg38 | sed "/^#/d" \
> ${path_gwas}pre_imputation/QC/${i}/nonlifted_hg38_variants_to_exclude_tmp.dat
done

for i in ${studies[@]}
do 
grep "alt" ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_postqc_lifted_hg38 | cut -f 4 \
| cat - ${path_gwas}pre_imputation/QC/${i}/nonlifted_hg38_variants_to_exclude_tmp.dat > \
${path_gwas}pre_imputation/QC/${i}/nonlifted_hg38_variants_to_exclude.dat
done

for i in ${studies[@]}
do ls -la ${path_gwas}pre_imputation/QC/${i}/nonlifted_hg38_variants_to_exclude.dat
done

for i in ${studies[@]}
do echo ${i} && wc -l ${path_gwas}pre_imputation/QC/${i}/nonlifted_hg38_variants_to_exclude.dat
done

######
# all_hce
# 72 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/all_hce/nonlifted_hg38_variants_to_exclude.dat
# niddk_old_gwas
# 49 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/nonlifted_hg38_variants_to_exclude.dat
# australia_omniexome
# 126 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/nonlifted_hg38_variants_to_exclude.dat
# gwas1
# 85 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/nonlifted_hg38_variants_to_exclude.dat
# gwas2
# 151 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas2/nonlifted_hg38_variants_to_exclude.dat
# pittsburgh_gsa
# 178 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pittsburgh_gsa/nonlifted_hg38_variants_to_exclude.dat
# spain_gsa
# 111 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/spain_gsa/nonlifted_hg38_variants_to_exclude.dat
# italy_gsa
# 109 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/italy_gsa/nonlifted_hg38_variants_to_exclude.dat
# kiel_austria_sibdcs_gsa
# 112 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/kiel_austria_sibdcs_gsa/nonlifted_hg38_variants_to_exclude.dat
# netherlands_gsa
# 163 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/nonlifted_hg38_variants_to_exclude.dat
# slovenia_gsa
# 118 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/slovenia_gsa/nonlifted_hg38_variants_to_exclude.dat
# sweden_gsa
# 119 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/nonlifted_hg38_variants_to_exclude.dat
# niddk_broad_gsa
# 129 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_broad_gsa/nonlifted_hg38_variants_to_exclude.dat
# niddk_feinstein_gsa
# 126 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_feinstein_gsa/nonlifted_hg38_variants_to_exclude.dat
# basque_gsa
# 114 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/basque_gsa/nonlifted_hg38_variants_to_exclude.dat
# lithuania_gsa
# 124 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lithuania_gsa/nonlifted_hg38_variants_to_exclude.dat
# belgium_louis_gsa
# 143 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_louis_gsa/nonlifted_hg38_variants_to_exclude.dat
# belgium_franchimont_gsa
# 133 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_franchimont_gsa/nonlifted_hg38_variants_to_exclude.dat
# belgium_vermeire_gsa
# 155 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_vermeire_gsa/nonlifted_hg38_variants_to_exclude.dat
# prism_nfe_gsa
# 130 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gsa/nonlifted_hg38_variants_to_exclude.dat
# prism_nfe_gwas
# 46 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/nonlifted_hg38_variants_to_exclude.dat
# finland_illugwas
# 45 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/nonlifted_hg38_variants_to_exclude.dat
# german_affy6_old_gwas
# 150 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/german_affy6_old_gwas/nonlifted_hg38_variants_to_exclude.dat
# norway_affy6_old_gwas
# 132 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/norway_affy6_old_gwas/nonlifted_hg38_variants_to_exclude.dat
# belgium_inf1_old_gwas
# 49 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/nonlifted_hg38_variants_to_exclude.dat
# belgium_inf2_old_gwas
# 46 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/nonlifted_hg38_variants_to_exclude.dat
# cedars_370k_old_gwas
# 60 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_370k_old_gwas/nonlifted_hg38_variants_to_exclude.dat
# cedars_610k_old_gwas
# 115 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_610k_old_gwas/nonlifted_hg38_variants_to_exclude.dat
# cedars_omni_old_gwas
# 170 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_omni_old_gwas/nonlifted_hg38_variants_to_exclude.dat
# swedish_uc_old_gwas
# 42 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/nonlifted_hg38_variants_to_exclude.dat
# mccauley_gsa
# 132 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_gsa/nonlifted_hg38_variants_to_exclude.dat
# ccfa_gsa
# 139 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/ccfa_gsa/nonlifted_hg38_variants_to_exclude.dat
# cedars_gsa
# 151 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_gsa/nonlifted_hg38_variants_to_exclude.dat
# bernstein_gsa
# 127 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/bernstein_gsa/nonlifted_hg38_variants_to_exclude.dat
# farkkila_gsa
# 120 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/farkkila_gsa/nonlifted_hg38_variants_to_exclude.dat
# franchimont_gsa
# 154 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franchimont_gsa/nonlifted_hg38_variants_to_exclude.dat
# franke_gsa
# 127 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franke_gsa/nonlifted_hg38_variants_to_exclude.dat
# helmsley_prism_gsa
# 46 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/nonlifted_hg38_variants_to_exclude.dat
# helmsley_xavier_pbrism_gsa
# 46 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/nonlifted_hg38_variants_to_exclude.dat
# hyams_protect_gsa
# 139 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/hyams_protect_gsa/nonlifted_hg38_variants_to_exclude.dat
# lewis_sparc_gsa
# 142 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lewis_sparc_gsa/nonlifted_hg38_variants_to_exclude.dat
# mccauley_new_gsa
# 139 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_new_gsa/nonlifted_hg38_variants_to_exclude.dat
# mcgovern_gsa
# 146 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mcgovern_gsa/nonlifted_hg38_variants_to_exclude.dat
# moayyedi_imagine_gsa
# 145 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/moayyedi_imagine_gsa/nonlifted_hg38_variants_to_exclude.dat
# newberry_share_gsa
# 141 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/newberry_share_gsa/nonlifted_hg38_variants_to_exclude.dat
# niddk_cho_gsa
# 142 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_cho_gsa/nonlifted_hg38_variants_to_exclude.dat
# niddk_duerr_gsa
# 134 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_duerr_gsa/nonlifted_hg38_variants_to_exclude.dat
# niddk_rioux_gsa
# 132 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_rioux_gsa/nonlifted_hg38_variants_to_exclude.dat
# niddk_silverberg_gsa
# 144 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_silverberg_gsa/nonlifted_hg38_variants_to_exclude.dat
# palotie_hus_gsa
# 133 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/palotie_hus_gsa/nonlifted_hg38_variants_to_exclude.dat
# pekow_share_gsa
# 138 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pekow_share_gsa/nonlifted_hg38_variants_to_exclude.dat
# rioux_igenomed_gsa
# 127 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/nonlifted_hg38_variants_to_exclude.dat
# sands_msccr_gsa
# 146 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sands_msccr_gsa/nonlifted_hg38_variants_to_exclude.dat
# stampfer_gsa
# 155 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/stampfer_gsa/nonlifted_hg38_variants_to_exclude.dat
# vermeire_gsa
# 155 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/vermeire_gsa/nonlifted_hg38_variants_to_exclude.dat
# weersma_gsa
# 134 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/weersma_gsa/nonlifted_hg38_variants_to_exclude.dat
# xavier_prism_gsa
# 132 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_prism_gsa/nonlifted_hg38_variants_to_exclude.dat
# xavier_share_gsa
# 135 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_share_gsa/nonlifted_hg38_variants_to_exclude.dat
######
### exclude non lifted:


studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=800
# MEM=3000 - all_hce gwas2 kiel_austria_sibdcs_gsa niddk_broad_gsa niddk_feinstein_gsa mcgovern_gsa

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_53_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_53_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het \
--exclude ${path_gwas}pre_imputation/QC/${i}/nonlifted_hg38_variants_to_exclude.dat \
--recode tab --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_postqc_lifted_hg38_liftedvariants"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_53_${i} | grep -E "completed"
done


###### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"


cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

length(cohorts)
# [1] 58

for (j in 1:length(cohorts)) {
  
  print(cohorts[j])
  
  bed_lifted<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_hg19_postqc_lifted_hg38",sep=""),head=F)
  excluded<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/nonlifted_hg38_variants_to_exclude.dat",sep=""),head=F)
  bed_lifted<-bed_lifted[which(!bed_lifted$V4 %in% excluded$V1),]
  map<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_hg19_postqc_lifted_hg38_liftedvariants.map",sep=""),head=F)

  map$pos<-bed_lifted$V2
  
  print(table(as.character(map$V2)==as.character(bed_lifted$V4)))
  print(table(map$pos==bed_lifted$V2))
  
  write.table(map[,c(1:3,5)],paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_hg19_postqc_lifted_hg38_liftedvariants_edited.map",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
}

########
# [1] "all_hce"
# 
# TRUE 
# 434526 
# 
# TRUE 
# 434526 
# 
# [1] "australia_omniexome"
# 
# TRUE 
# 739186 
# 
# TRUE 
# 739186 
# [1] "gwas1"
# 
# TRUE 
# 453331 
# 
# TRUE 
# 453331 
# [1] "gwas2"
# 
# TRUE 
# 773077 
# 
# TRUE 
# 773077 
# [1] "pittsburgh_gsa"
# 
# TRUE 
# 897664 
# 
# TRUE 
# 897664 
# [1] "spain_gsa"
# 
# TRUE 
# 575693 
# 
# TRUE 
# 575693 
# [1] "italy_gsa"
# 
# TRUE 
# 573741 
# 
# TRUE 
# 573741 
# [1] "kiel_austria_sibdcs_gsa"
# 
# TRUE 
# 564399 
# 
# TRUE 
# 564399 
# [1] "netherlands_gsa"
# 
# TRUE 
# 666028 
# 
# TRUE 
# 666028 
# [1] "slovenia_gsa"
# 
# TRUE 
# 529490 
# 
# TRUE 
# 529490 
# [1] "sweden_gsa"
# 
# TRUE 
# 558488 
# 
# TRUE 
# 558488 
# [1] "niddk_broad_gsa"
# 
# TRUE 
# 590497 
# 
# TRUE 
# 590497 
# [1] "niddk_feinstein_gsa"
# 
# TRUE 
# 603088 
# 
# TRUE 
# 603088 
# [1] "basque_gsa"
# 
# TRUE 
# 567287 
# 
# TRUE 
# 567287 
# [1] "prism_nfe_gsa"
# 
# TRUE 
# 558887 
# 
# TRUE 
# 558887 
# [1] "lithuania_gsa"
# 
# TRUE 
# 554011 
# 
# TRUE 
# 554011 
# [1] "belgium_louis_gsa"
# 
# TRUE 
# 591805 
# 
# TRUE 
# 591805 
# [1] "belgium_franchimont_gsa"
# 
# TRUE 
# 561861 
# 
# TRUE 
# 561861 
# [1] "belgium_vermeire_gsa"
# 
# TRUE 
# 652992 
# 
# TRUE 
# 652992 
# [1] "prism_nfe_gwas"
# 
# TRUE 
# 242853 
# 
# TRUE 
# 242853 
# [1] "finland_illugwas"
# 
# TRUE 
# 237010 
# 
# TRUE 
# 237010 
# [1] "german_affy6_old_gwas"
# 
# TRUE 
# 783289 
# 
# TRUE 
# 783289 
# [1] "norway_affy6_old_gwas"
# 
# TRUE 
# 725767 
# 
# TRUE 
# 725767 
# [1] "belgium_inf1_old_gwas"
# 
# TRUE 
# 301665 
# 
# TRUE 
# 301665 
# [1] "belgium_inf2_old_gwas"
# 
# TRUE 
# 290145 
# 
# TRUE 
# 290145 
# [1] "niddk_old_gwas"
# 
# TRUE 
# 300117 
# 
# TRUE 
# 300117 
# [1] "cedars_370k_old_gwas"
# 
# TRUE 
# 338657 
# 
# TRUE 
# 338657 
# [1] "cedars_610k_old_gwas"
# 
# TRUE 
# 583792 
# 
# TRUE 
# 583792 
# [1] "cedars_omni_old_gwas"
# 
# TRUE 
# 718969 
# 
# TRUE 
# 718969 
# [1] "swedish_uc_old_gwas"
# 
# TRUE 
# 296606 
# 
# TRUE 
# 296606 
# [1] "mccauley_gsa"
# 
# TRUE 
# 556550 
# 
# TRUE 
# 556550 
# [1] "ccfa_gsa"
# 
# TRUE 
# 599842 
# 
# TRUE 
# 599842 
# [1] "cedars_gsa"
# 
# TRUE 
# 646151 
# 
# TRUE 
# 646151 
# [1] "bernstein_gsa"
# 
# TRUE 
# 553817 
# 
# TRUE 
# 553817 
# [1] "farkkila_gsa"
# 
# TRUE 
# 511392 
# 
# TRUE 
# 511392 
# [1] "franchimont_gsa"
# 
# TRUE 
# 648057 
# 
# TRUE 
# 648057 
# [1] "franke_gsa"
# 
# TRUE 
# 557369 
# 
# TRUE 
# 557369 
# [1] "helmsley_prism_gsa"
# 
# TRUE 
# 243057 
# 
# TRUE 
# 243057 
# [1] "helmsley_xavier_prism_gsa"
# 
# TRUE 
# 243036 
# 
# TRUE 
# 243036 
# [1] "hyams_protect_gsa"
# 
# TRUE 
# 569329 
# 
# TRUE 
# 569329 
# [1] "lewis_sparc_gsa"
# 
# TRUE 
# 607573 
# 
# TRUE 
# 607573 
# [1] "mccauley_new_gsa"
# 
# TRUE 
# 578546 
# 
# TRUE 
# 578546 
# [1] "mcgovern_gsa"
# 
# TRUE 
# 630809 
# 
# TRUE 
# 630809 
# [1] "moayyedi_imagine_gsa"
# 
# TRUE 
# 607858 
# 
# TRUE 
# 607858 
# [1] "newberry_share_gsa"
# 
# TRUE 
# 583553 
# 
# TRUE 
# 583553 
# [1] "niddk_cho_gsa"
# 
# TRUE 
# 593434 
# 
# TRUE 
# 593434 
# [1] "niddk_duerr_gsa"
# 
# TRUE 
# 581008 
# 
# TRUE 
# 581008 
# [1] "niddk_rioux_gsa"
# 
# TRUE 
# 574839 
# 
# TRUE 
# 574839 
# [1] "niddk_silverberg_gsa"
# 
# TRUE 
# 611096 
# 
# TRUE 
# 611096 
# [1] "palotie_hus_gsa"
# 
# TRUE 
# 553534 
# 
# TRUE 
# 553534 
# [1] "pekow_share_gsa"
# 
# TRUE 
# 573922 
# 
# TRUE 
# 573922 
# [1] "rioux_igenomed_gsa"
# 
# TRUE 
# 543215 
# 
# TRUE 
# 543215 
# [1] "sands_msccr_gsa"
# 
# TRUE 
# 595412 
# 
# TRUE 
# 595412 
# [1] "stampfer_gsa"
# 
# TRUE 
# 643531 
# 
# TRUE 
# 643531 
# [1] "vermeire_gsa"
# 
# TRUE 
# 653188 
# 
# TRUE 
# 653188 
# [1] "weersma_gsa"
# 
# TRUE 
# 566267 
# 
# TRUE 
# 566267 
# [1] "xavier_prism_gsa"
# 
# TRUE 
# 569608 
# 
# TRUE 
# 569608 
# [1] "xavier_share_gsa"
# 
# TRUE 
# 568375 
# 
# TRUE 
# 568375 
########

########

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
# kiel_austria_sibdcs_gsa run interactively

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=3000

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_54_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_54_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--ped ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_postqc_lifted_hg38_liftedvariants.ped \
--map ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_postqc_lifted_hg38_liftedvariants_edited.map \
--merge-x 'no-fail' \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_54_${i} | grep -E "completed"
done


######################################################
# 16.1 - REMOVE MONOMORPHIC VARIANTS

###### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"
cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

length(cohorts)
# [1] 58

for (j in 1:length(cohorts)) {
  
  print(cohorts[j])
  bim<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38.bim",sep=""),head=F)
  table(bim$V5,bim$V6)
  table(bim[which(bim$V5=="0"),"V1"])
  monom<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
  print(dim(monom))
  write.table(monom[,"V2",drop=F],paste(path,"pre_imputation/QC/",cohorts[j],"/list_indel_var_exclude_2",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
}

########
# [1] "kiel_austria_sibdcs_gsa"
# [1] 20147     6
# [1] "all_hce"
# [1] 12870     6
# [1] "australia_omniexome"
# [1] 936   6
# [1] "gwas1"
# [1] 7571    6
# [1] "gwas2"
# [1] 505   6
# [1] "pittsburgh_gsa"
# [1] 6653    6
# [1] "spain_gsa"
# [1] 22  6
# [1] "italy_gsa"
# [1] 13303     6
# [1] "netherlands_gsa"
# [1] 81425     6
# [1] "slovenia_gsa"
# [1] 5533    6
# [1] "sweden_gsa"
# [1] 3200    6
# [1] "niddk_broad_gsa"
# [1] 1676    6
# [1] "niddk_feinstein_gsa"
# [1] 2360    6
# [1] "basque_gsa"
# [1] 7331    6
# [1] "prism_nfe_gsa"
# [1] 410   6
# [1] "lithuania_gsa"
# [1] 614   6
# [1] "belgium_louis_gsa"
# [1] 29081     6
# [1] "belgium_franchimont_gsa"
# [1] 8612    6
# [1] "belgium_vermeire_gsa"
# [1] 16345     6
# [1] "prism_nfe_gwas"
# [1] 5 6
# [1] "finland_illugwas"
# [1] 1024    6
# [1] "german_affy6_old_gwas"
# [1] 1037    6
# [1] "norway_affy6_old_gwas"
# [1] 2501    6
# [1] "belgium_inf1_old_gwas"
# [1] 0 6
# [1] "belgium_inf2_old_gwas"
# [1] 0 6
# [1] "niddk_old_gwas"
# [1] 9 6
# [1] "cedars_370k_old_gwas"
# [1] 23  6
# [1] "cedars_610k_old_gwas"
# [1] 121   6
# [1] "cedars_omni_old_gwas"
# [1] 80  6
# [1] "swedish_uc_old_gwas"
# [1] 0 6
# [1] "mccauley_gsa"
# [1] 357   6
# [1] "ccfa_gsa"
# [1] 375   6
# [1] "cedars_gsa"
# [1] 36248     6
# [1] "bernstein_gsa"
# [1] 654   6
# [1] "farkkila_gsa"
# [1] 17048     6
# [1] "franchimont_gsa"
# [1] 75675     6
# [1] "franke_gsa"
# [1] 562   6
# [1] "helmsley_prism_gsa"
# [1] 26  6
# [1] "helmsley_xavier_prism_gsa"
# [1] 1 6
# [1] "hyams_protect_gsa"
# [1] 157   6
# [1] "lewis_sparc_gsa"
# [1] 432   6
# [1] "mccauley_new_gsa"
# [1] 539   6
# [1] "mcgovern_gsa"
# [1] 4064    6
# [1] "moayyedi_imagine_gsa"
# [1] 19041     6
# [1] "newberry_share_gsa"
# [1] 8010    6
# [1] "niddk_cho_gsa"
# [1] 534   6
# [1] "niddk_duerr_gsa"
# [1] 247   6
# [1] "niddk_rioux_gsa"
# [1] 5057    6
# [1] "niddk_silverberg_gsa"
# [1] 10397     6
# [1] "palotie_hus_gsa"
# [1] 611   6
# [1] "pekow_share_gsa"
# [1] 1301    6
# [1] "rioux_igenomed_gsa"
# [1] 304   6
# [1] "sands_msccr_gsa"
# [1] 924   6
# [1] "stampfer_gsa"
# [1] 72504     6
# [1] "vermeire_gsa"
# [1] 68923     6
# [1] "weersma_gsa"
# [1] 6007    6
# [1] "xavier_prism_gsa"
# [1] 1620    6
# [1] "xavier_share_gsa"
# [1] 484   6
########


######################################################
# 16.2 - CREATE FILE FORCING A1 ALLELE TO BE REFERENCE

# force A1 and A2 to be ref and alt alleles:
# studies=(kiel_austria_sibdcs_gsa) PENDING

studies=(australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do 
zcat ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_posstrandaligned_2.vcf.gz | cut -f '1-5' | awk '{print $1":"$2"_"$4"_"$5,$4}' | \
sed '/^##/d' > ${path_gwas}pre_imputation/QC/${i}/list_variants_${i}_hg19_posstrandaligned_with_A1
done
for i in ${studies[@]}
do 
sort ${path_gwas}pre_imputation/QC/${i}/list_variants_${i}_hg19_posstrandaligned_with_A1 | uniq \
> ${path_gwas}pre_imputation/QC/${i}/list_variants_${i}_hg19_posstrandaligned_with_A1_ed
done

# two exceptions:

i=niddk_old_gwas
zcat ${path_gwas}pre_imputation/QC/niddk_uc_old_gwas/niddk_uc_old_gwas_posstrandaligned.vcf.gz | cut -f '1-5' | awk '{print $1":"$2"_"$4"_"$5,$4}' | \
sed '/^##/d' > ${path_gwas}pre_imputation/QC/${i}/list_variants_niddk_uc_old_gwas_hg19_posstrandaligned_with_A1

zcat ${path_gwas}pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_posstrandaligned.vcf.gz | cut -f '1-5' | awk '{print $1":"$2"_"$4"_"$5,$4}' | \
sed '/^##/d' > ${path_gwas}pre_imputation/QC/${i}/list_variants_niddk_cd_old_gwas_hg19_posstrandaligned_with_A1

cat ${path_gwas}pre_imputation/QC/${i}/list_variants_niddk_uc_old_gwas_hg19_posstrandaligned_with_A1 \
${path_gwas}pre_imputation/QC/${i}/list_variants_niddk_cd_old_gwas_hg19_posstrandaligned_with_A1 \
> ${path_gwas}pre_imputation/QC/${i}/list_variants_niddk_old_gwas_hg19_posstrandaligned_with_A1

sort ${path_gwas}pre_imputation/QC/${i}/list_variants_niddk_old_gwas_hg19_posstrandaligned_with_A1 | uniq \
> ${path_gwas}pre_imputation/QC/${i}/list_variants_niddk_old_gwas_hg19_posstrandaligned_with_A1_ed


i=all_hce
zcat ${path_gwas}pre_imputation/QC/${i}/gwas3_hg19_posstrandaligned_2.vcf.gz | cut -f '1-5' | awk '{print $1":"$2"_"$4"_"$5,$4}' | \
sed '/^##/d' > ${path_gwas}pre_imputation/QC/${i}/list_variants_gwas3_hg19_posstrandaligned_with_A1

zcat ${path_gwas}pre_imputation/QC/${i}/new_wave_hg19_posstrandaligned_2.vcf.gz | cut -f '1-5' | awk '{print $1":"$2"_"$4"_"$5,$4}' | \
sed '/^##/d' > ${path_gwas}pre_imputation/QC/${i}/list_variants_new_wave_hg19_posstrandaligned_with_A1

cat ${path_gwas}pre_imputation/QC/${i}/list_variants_gwas3_hg19_posstrandaligned_with_A1 \
${path_gwas}pre_imputation/QC/${i}/list_variants_new_wave_hg19_posstrandaligned_with_A1 \
> ${path_gwas}pre_imputation/QC/${i}/list_variants_${i}_hg19_posstrandaligned_with_A1

sort ${path_gwas}pre_imputation/QC/${i}/list_variants_${i}_hg19_posstrandaligned_with_A1 | uniq > \
${path_gwas}pre_imputation/QC/${i}/list_variants_${i}_hg19_posstrandaligned_with_A1_ed


studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do 
ls -la ${path_gwas}pre_imputation/QC/${i}/list_variants_${i}_hg19_posstrandaligned_with_A1_ed
done

for i in ${studies[@]}
do 
ls -la ${path_gwas}/pre_imputation/QC/${i}/list_indel_var_exclude_2
done

##########

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=800

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_55_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_55_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38 \
--exclude ${path_gwas}/pre_imputation/QC/${i}/list_indel_var_exclude_2 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_55_${i} | grep -E "completed"
done

##########

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=800

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_56_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_56_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom \
--allow-no-sex \
--a2-allele ${path_gwas}pre_imputation/QC/${i}/list_variants_${i}_hg19_posstrandaligned_with_A1_ed \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_56_${i} | grep -E "completed"
done


##############################################################
# 16.1 DOUBLE CHECK REF ALT ALLELE WITH HG38 REF SEQUENCE:

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=800

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_57_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_57_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt \
--allow-no-sex \
--keep-allele-order \
--output-chr M --recode vcf-iid --out ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt"
done


for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_57_${i} | grep -E "completed"
done

# Note that most PLINK analyses treat the A1 (usually minor) allele as the reference allele, which makes sense when only biallelic variants are involved.
# However, since it is conventional for VCF files to set the major allele as the reference allele instead

# 4.2 FLIP ALLELES TO + STRAND:

path_gwas=/path/to/ibdgwas/IIBDGC/
export BCFTOOLS_PLUGINS=/path/to/software/bcftools-1.16/plugins
MEM=500


for i in ${studies[@]}
do
bsub -J"bcf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_bcftools_3_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_bcftools_3_${i} \
"/path/to/software/bcftools-1.16/./bcftools +fixref \
${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt.vcf \
-Oz -o ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned.vcf.gz \
-- -f ${path_gwas}resources/hg38/hg38_edited.fa -m top"
done

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
for i in ${studies[@]}
do
echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_bcftools_3_${i} | grep "Successfully completed."
done

for i in ${studies[@]}
do
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned.vcf.gz
done

for i in ${studies[@]}
do
echo ${i} && tail -24 ${path_gwas}pre_imputation/QC/${i}/logs/stderr_bcftools_3_${i}
done


##############
# all_hce
# ST	A>G	60846	14.4%
# ST	A>T	2632	0.6%
# ST	C>A	20826	4.9%
# ST	C>G	6644	1.6%
# ST	C>T	105063	24.9%
# ST	G>A	106443	25.2%
# ST	G>C	6646	1.6%
# ST	G>T	20188	4.8%
# ST	T>A	2893	0.7%
# ST	T>C	60925	14.4%
# ST	T>G	14167	3.4%
# # NS, Number of sites:
# NS	total        	421656
# NS	ref match    	420273	99.7%
# NS	ref mismatch 	1381	0.3%
# NS	flipped      	163	0.0%
# NS	swapped      	1059	0.3%
# NS	flip+swap    	7788	1.8%
# NS	unresolved   	17	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	2
# NS	non-ACGT     	2
# NS	non-SNP      	0
# NS	non-biallelic	0
# niddk_old_gwas
# ST	A>G	59547	19.8%
# ST	A>T	0	0.0%
# ST	C>A	14055	4.7%
# ST	C>G	0	0.0%
# ST	C>T	63124	21.0%
# ST	G>A	63360	21.1%
# ST	G>C	0	0.0%
# ST	G>T	13927	4.6%
# ST	T>A	0	0.0%
# ST	T>C	59481	19.8%
# ST	T>G	13310	4.4%
# # NS, Number of sites:
# NS	total        	300108
# NS	ref match    	298947	99.6%
# NS	ref mismatch 	1161	0.4%
# NS	flipped      	39	0.0%
# NS	swapped      	1067	0.4%
# NS	flip+swap    	55	0.0%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# australia_omniexome
# ST	A>G	135847	18.4%
# ST	A>T	1909	0.3%
# ST	C>A	34715	4.7%
# ST	C>G	4337	0.6%
# ST	C>T	160632	21.8%
# ST	G>A	160743	21.8%
# ST	G>C	4188	0.6%
# ST	G>T	34483	4.7%
# ST	T>A	1921	0.3%
# ST	T>C	136210	18.5%
# ST	T>G	31562	4.3%
# # NS, Number of sites:
# NS	total        	738250
# NS	ref match    	735608	99.6%
# NS	ref mismatch 	2641	0.4%
# NS	flipped      	177	0.0%
# NS	swapped      	2283	0.3%
# NS	flip+swap    	6373	0.9%
# NS	unresolved   	11	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	1
# NS	non-ACGT     	1
# NS	non-SNP      	0
# NS	non-biallelic	0
# gwas1
# ST	A>G	73294	16.4%
# ST	A>T	14006	3.1%
# ST	C>A	17807	4.0%
# ST	C>G	21205	4.8%
# ST	C>T	79890	17.9%
# ST	G>A	79720	17.9%
# ST	G>C	21106	4.7%
# ST	G>T	17913	4.0%
# ST	T>A	14080	3.2%
# ST	T>C	73006	16.4%
# ST	T>G	16838	3.8%
# # NS, Number of sites:
# NS	total        	445760
# NS	ref match    	444401	99.7%
# NS	ref mismatch 	1359	0.3%
# NS	flipped      	151	0.0%
# NS	swapped      	1135	0.3%
# NS	flip+swap    	34119	7.7%
# NS	unresolved   	22	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# gwas2
# ST	A>G	126129	16.3%
# ST	A>T	24000	3.1%
# ST	C>A	31794	4.1%
# ST	C>G	33969	4.4%
# ST	C>T	141302	18.3%
# ST	G>A	141411	18.3%
# ST	G>C	33871	4.4%
# ST	G>T	31865	4.1%
# ST	T>A	23812	3.1%
# ST	T>C	126140	16.3%
# ST	T>G	29072	3.8%
# # NS, Number of sites:
# NS	total        	772572
# NS	ref match    	769783	99.6%
# NS	ref mismatch 	2789	0.4%
# NS	flipped      	289	0.0%
# NS	swapped      	2341	0.3%
# NS	flip+swap    	57633	7.5%
# NS	unresolved   	47	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# pittsburgh_gsa
# ST	A>G	166978	18.7%
# ST	A>T	4167	0.5%
# ST	C>A	42029	4.7%
# ST	C>G	6834	0.8%
# ST	C>T	185631	20.8%
# ST	G>A	184881	20.7%
# ST	G>C	6857	0.8%
# ST	G>T	42087	4.7%
# ST	T>A	4236	0.5%
# ST	T>C	167716	18.8%
# ST	T>G	39801	4.5%
# # NS, Number of sites:
# NS	total        	891011
# NS	ref match    	887814	99.6%
# NS	ref mismatch 	3197	0.4%
# NS	flipped      	230	0.0%
# NS	swapped      	2765	0.3%
# NS	flip+swap    	11127	1.2%
# NS	unresolved   	23	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# spain_gsa
# ST	A>G	111674	19.4%
# ST	A>T	555	0.1%
# ST	C>A	27564	4.8%
# ST	C>G	875	0.2%
# ST	C>T	120920	21.0%
# ST	G>A	120919	21.0%
# ST	G>C	914	0.2%
# ST	G>T	27276	4.7%
# ST	T>A	574	0.1%
# ST	T>C	112313	19.5%
# ST	T>G	26050	4.5%
# # NS, Number of sites:
# NS	total        	575671
# NS	ref match    	573481	99.6%
# NS	ref mismatch 	2188	0.4%
# NS	flipped      	119	0.0%
# NS	swapped      	1944	0.3%
# NS	flip+swap    	1567	0.3%
# NS	unresolved   	3	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	2
# NS	non-ACGT     	2
# NS	non-SNP      	0
# NS	non-biallelic	0
# italy_gsa
# ST	A>G	98387	17.6%
# ST	A>T	891	0.2%
# ST	C>A	28001	5.0%
# ST	C>G	1386	0.2%
# ST	C>T	127847	22.8%
# ST	G>A	127814	22.8%
# ST	G>C	1318	0.2%
# ST	G>T	27514	4.9%
# ST	T>A	866	0.2%
# ST	T>C	98528	17.6%
# ST	T>G	23822	4.3%
# # NS, Number of sites:
# NS	total        	560438
# NS	ref match    	558688	99.7%
# NS	ref mismatch 	1750	0.3%
# NS	flipped      	121	0.0%
# NS	swapped      	1498	0.3%
# NS	flip+swap    	2396	0.4%
# NS	unresolved   	1	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# kiel_austria_sibdcs_gsa
# ST	A>G	90646	16.7%
# ST	A>T	776	0.1%
# ST	C>A	27623	5.1%
# ST	C>G	1385	0.3%
# ST	C>T	130019	23.9%
# ST	G>A	129660	23.8%
# ST	G>C	1358	0.2%
# ST	G>T	27275	5.0%
# ST	T>A	795	0.1%
# ST	T>C	90931	16.7%
# ST	T>G	21879	4.0%
# # NS, Number of sites:
# NS	total        	544252
# NS	ref match    	542819	99.7%
# NS	ref mismatch 	1433	0.3%
# NS	flipped      	132	0.0%
# NS	swapped      	1174	0.2%
# NS	flip+swap    	2343	0.4%
# NS	unresolved   	1	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# netherlands_gsa
# ST	A>G	102078	17.5%
# ST	A>T	802	0.1%
# ST	C>A	29495	5.0%
# ST	C>G	1309	0.2%
# ST	C>T	133771	22.9%
# ST	G>A	133398	22.8%
# ST	G>C	1312	0.2%
# ST	G>T	29092	5.0%
# ST	T>A	833	0.1%
# ST	T>C	102373	17.5%
# ST	T>G	24974	4.3%
# # NS, Number of sites:
# NS	total        	584603
# NS	ref match    	582979	99.7%
# NS	ref mismatch 	1624	0.3%
# NS	flipped      	143	0.0%
# NS	swapped      	1338	0.2%
# NS	flip+swap    	2307	0.4%
# NS	unresolved   	1	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# slovenia_gsa
# ST	A>G	91817	17.5%
# ST	A>T	651	0.1%
# ST	C>A	26377	5.0%
# ST	C>G	904	0.2%
# ST	C>T	120004	22.9%
# ST	G>A	119802	22.9%
# ST	G>C	879	0.2%
# ST	G>T	26123	5.0%
# ST	T>A	665	0.1%
# ST	T>C	91929	17.5%
# ST	T>G	22348	4.3%
# # NS, Number of sites:
# NS	total        	523957
# NS	ref match    	522432	99.7%
# NS	ref mismatch 	1525	0.3%
# NS	flipped      	120	0.0%
# NS	swapped      	1281	0.2%
# NS	flip+swap    	1679	0.3%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# sweden_gsa
# ST	A>G	96965	17.5%
# ST	A>T	739	0.1%
# ST	C>A	28029	5.0%
# ST	C>G	1109	0.2%
# ST	C>T	127301	22.9%
# ST	G>A	126994	22.9%
# ST	G>C	1054	0.2%
# ST	G>T	27647	5.0%
# ST	T>A	761	0.1%
# ST	T>C	97255	17.5%
# ST	T>G	23663	4.3%
# # NS, Number of sites:
# NS	total        	555288
# NS	ref match    	553702	99.7%
# NS	ref mismatch 	1586	0.3%
# NS	flipped      	128	0.0%
# NS	swapped      	1321	0.2%
# NS	flip+swap    	2006	0.4%
# NS	unresolved   	1	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# niddk_broad_gsa
# ST	A>G	100101	17.0%
# ST	A>T	874	0.1%
# ST	C>A	30027	5.1%
# ST	C>G	1491	0.3%
# ST	C>T	137647	23.4%
# ST	G>A	137170	23.3%
# ST	G>C	1474	0.3%
# ST	G>T	29658	5.0%
# ST	T>A	927	0.2%
# ST	T>C	100448	17.1%
# ST	T>G	24435	4.1%
# # NS, Number of sites:
# NS	total        	588821
# NS	ref match    	587209	99.7%
# NS	ref mismatch 	1612	0.3%
# NS	flipped      	154	0.0%
# NS	swapped      	1308	0.2%
# NS	flip+swap    	2552	0.4%
# NS	unresolved   	2	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# niddk_feinstein_gsa
# ST	A>G	101963	17.0%
# ST	A>T	1121	0.2%
# ST	C>A	30401	5.1%
# ST	C>G	1995	0.3%
# ST	C>T	140129	23.3%
# ST	G>A	139532	23.2%
# ST	G>C	2037	0.3%
# ST	G>T	30069	5.0%
# ST	T>A	1156	0.2%
# ST	T>C	102277	17.0%
# ST	T>G	24914	4.1%
# # NS, Number of sites:
# NS	total        	600728
# NS	ref match    	599106	99.7%
# NS	ref mismatch 	1622	0.3%
# NS	flipped      	159	0.0%
# NS	swapped      	1299	0.2%
# NS	flip+swap    	3347	0.6%
# NS	unresolved   	5	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# basque_gsa
# ST	A>G	98028	17.5%
# ST	A>T	900	0.2%
# ST	C>A	28001	5.0%
# ST	C>G	1450	0.3%
# ST	C>T	128149	22.9%
# ST	G>A	128029	22.9%
# ST	G>C	1364	0.2%
# ST	G>T	27584	4.9%
# ST	T>A	880	0.2%
# ST	T>C	97954	17.5%
# ST	T>G	23750	4.2%
# # NS, Number of sites:
# NS	total        	559956
# NS	ref match    	558200	99.7%
# NS	ref mismatch 	1756	0.3%
# NS	flipped      	120	0.0%
# NS	swapped      	1506	0.3%
# NS	flip+swap    	2472	0.4%
# NS	unresolved   	1	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# lithuania_gsa
# ST	A>G	96034	17.4%
# ST	A>T	732	0.1%
# ST	C>A	27829	5.0%
# ST	C>G	1111	0.2%
# ST	C>T	127752	23.1%
# ST	G>A	127291	23.0%
# ST	G>C	1104	0.2%
# ST	G>T	27636	5.0%
# ST	T>A	749	0.1%
# ST	T>C	96296	17.4%
# ST	T>G	23401	4.2%
# # NS, Number of sites:
# NS	total        	553397
# NS	ref match    	551758	99.7%
# NS	ref mismatch 	1637	0.3%
# NS	flipped      	149	0.0%
# NS	swapped      	1326	0.2%
# NS	flip+swap    	2029	0.4%
# NS	unresolved   	5	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	2
# NS	non-ACGT     	2
# NS	non-SNP      	0
# NS	non-biallelic	0
# belgium_louis_gsa
# ST	A>G	96941	17.2%
# ST	A>T	804	0.1%
# ST	C>A	28463	5.1%
# ST	C>G	1303	0.2%
# ST	C>T	130246	23.1%
# ST	G>A	129848	23.1%
# ST	G>C	1249	0.2%
# ST	G>T	28187	5.0%
# ST	T>A	846	0.2%
# ST	T>C	97237	17.3%
# ST	T>G	23731	4.2%
# # NS, Number of sites:
# NS	total        	562724
# NS	ref match    	561071	99.7%
# NS	ref mismatch 	1653	0.3%
# NS	flipped      	145	0.0%
# NS	swapped      	1352	0.2%
# NS	flip+swap    	2267	0.4%
# NS	unresolved   	3	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# belgium_franchimont_gsa
# ST	A>G	95834	17.3%
# ST	A>T	770	0.1%
# ST	C>A	27879	5.0%
# ST	C>G	1216	0.2%
# ST	C>T	127636	23.1%
# ST	G>A	127372	23.0%
# ST	G>C	1186	0.2%
# ST	G>T	27638	5.0%
# ST	T>A	794	0.1%
# ST	T>C	96022	17.4%
# ST	T>G	23398	4.2%
# # NS, Number of sites:
# NS	total        	553249
# NS	ref match    	551627	99.7%
# NS	ref mismatch 	1621	0.3%
# NS	flipped      	132	0.0%
# NS	swapped      	1338	0.2%
# NS	flip+swap    	2137	0.4%
# NS	unresolved   	3	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	1
# NS	non-ACGT     	1
# NS	non-SNP      	0
# NS	non-biallelic	0
# belgium_vermeire_gsa
# ST	A>G	103331	16.2%
# ST	A>T	1343	0.2%
# ST	C>A	33589	5.3%
# ST	C>G	2425	0.4%
# ST	C>T	152095	23.9%
# ST	G>A	151985	23.9%
# ST	G>C	2412	0.4%
# ST	G>T	33280	5.2%
# ST	T>A	1405	0.2%
# ST	T>C	103661	16.3%
# ST	T>G	25470	4.0%
# # NS, Number of sites:
# NS	total        	636647
# NS	ref match    	634829	99.7%
# NS	ref mismatch 	1814	0.3%
# NS	flipped      	197	0.0%
# NS	swapped      	1409	0.2%
# NS	flip+swap    	3994	0.6%
# NS	unresolved   	6	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	4
# NS	non-ACGT     	4
# NS	non-SNP      	0
# NS	non-biallelic	0
# prism_nfe_gsa
# ST	A>G	97293	17.4%
# ST	A>T	752	0.1%
# ST	C>A	28201	5.0%
# ST	C>G	1165	0.2%
# ST	C>T	128021	22.9%
# ST	G>A	127844	22.9%
# ST	G>C	1153	0.2%
# ST	G>T	27936	5.0%
# ST	T>A	787	0.1%
# ST	T>C	97538	17.5%
# ST	T>G	23828	4.3%
# # NS, Number of sites:
# NS	total        	558477
# NS	ref match    	556855	99.7%
# NS	ref mismatch 	1622	0.3%
# NS	flipped      	141	0.0%
# NS	swapped      	1332	0.2%
# NS	flip+swap    	2122	0.4%
# NS	unresolved   	2	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# prism_nfe_gwas
# ST	A>G	47859	19.7%
# ST	A>T	68	0.0%
# ST	C>A	11651	4.8%
# ST	C>G	117	0.0%
# ST	C>T	50657	20.9%
# ST	G>A	50510	20.8%
# ST	G>C	104	0.0%
# ST	G>T	11772	4.8%
# ST	T>A	76	0.0%
# ST	T>C	47623	19.6%
# ST	T>G	11119	4.6%
# # NS, Number of sites:
# NS	total        	242848
# NS	ref match    	241888	99.6%
# NS	ref mismatch 	960	0.4%
# NS	flipped      	41	0.0%
# NS	swapped      	867	0.4%
# NS	flip+swap    	230	0.1%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# finland_illugwas
# ST	A>G	46477	19.7%
# ST	A>T	65	0.0%
# ST	C>A	11340	4.8%
# ST	C>G	116	0.0%
# ST	C>T	49229	20.9%
# ST	G>A	49126	20.8%
# ST	G>C	105	0.0%
# ST	G>T	11448	4.9%
# ST	T>A	74	0.0%
# ST	T>C	46281	19.6%
# ST	T>G	10789	4.6%
# # NS, Number of sites:
# NS	total        	235986
# NS	ref match    	235040	99.6%
# NS	ref mismatch 	946	0.4%
# NS	flipped      	39	0.0%
# NS	swapped      	855	0.4%
# NS	flip+swap    	225	0.1%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# german_affy6_old_gwas
# ST	A>G	127303	16.3%
# ST	A>T	24130	3.1%
# ST	C>A	32245	4.1%
# ST	C>G	34339	4.4%
# ST	C>T	143668	18.4%
# ST	G>A	143743	18.4%
# ST	G>C	34125	4.4%
# ST	G>T	32337	4.1%
# ST	T>A	23993	3.1%
# ST	T>C	127474	16.3%
# ST	T>G	29407	3.8%
# # NS, Number of sites:
# NS	total        	782252
# NS	ref match    	779535	99.7%
# NS	ref mismatch 	2717	0.3%
# NS	flipped      	279	0.0%
# NS	swapped      	2280	0.3%
# NS	flip+swap    	56872	7.3%
# NS	unresolved   	43	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# norway_affy6_old_gwas
# ST	A>G	118789	16.4%
# ST	A>T	22226	3.1%
# ST	C>A	29592	4.1%
# ST	C>G	31759	4.4%
# ST	C>T	132063	18.3%
# ST	G>A	131810	18.2%
# ST	G>C	31534	4.4%
# ST	G>T	29577	4.1%
# ST	T>A	22123	3.1%
# ST	T>C	118718	16.4%
# ST	T>G	27542	3.8%
# # NS, Number of sites:
# NS	total        	723266
# NS	ref match    	720698	99.6%
# NS	ref mismatch 	2568	0.4%
# NS	flipped      	261	0.0%
# NS	swapped      	2164	0.3%
# NS	flip+swap    	52692	7.3%
# NS	unresolved   	43	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# belgium_inf1_old_gwas
# ST	A>G	59989	19.9%
# ST	A>T	0	0.0%
# ST	C>A	14002	4.6%
# ST	C>G	0	0.0%
# ST	C>T	63713	21.1%
# ST	G>A	63832	21.2%
# ST	G>C	0	0.0%
# ST	G>T	13902	4.6%
# ST	T>A	0	0.0%
# ST	T>C	59783	19.8%
# ST	T>G	13208	4.4%
# # NS, Number of sites:
# NS	total        	301665
# NS	ref match    	300502	99.6%
# NS	ref mismatch 	1163	0.4%
# NS	flipped      	40	0.0%
# NS	swapped      	1063	0.4%
# NS	flip+swap    	60	0.0%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# belgium_inf2_old_gwas
# ST	A>G	57358	19.8%
# ST	A>T	0	0.0%
# ST	C>A	13788	4.8%
# ST	C>G	0	0.0%
# ST	C>T	60960	21.0%
# ST	G>A	61249	21.1%
# ST	G>C	0	0.0%
# ST	G>T	13679	4.7%
# ST	T>A	0	0.0%
# ST	T>C	57101	19.7%
# ST	T>G	13022	4.5%
# # NS, Number of sites:
# NS	total        	290145
# NS	ref match    	289026	99.6%
# NS	ref mismatch 	1119	0.4%
# NS	flipped      	36	0.0%
# NS	swapped      	1027	0.4%
# NS	flip+swap    	56	0.0%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# cedars_370k_old_gwas
# ST	A>G	66443	19.6%
# ST	A>T	301	0.1%
# ST	C>A	16163	4.8%
# ST	C>G	415	0.1%
# ST	C>T	70652	20.9%
# ST	G>A	70662	20.9%
# ST	G>C	432	0.1%
# ST	G>T	16078	4.7%
# ST	T>A	306	0.1%
# ST	T>C	66401	19.6%
# ST	T>G	15324	4.5%
# # NS, Number of sites:
# NS	total        	338634
# NS	ref match    	337186	99.6%
# NS	ref mismatch 	1447	0.4%
# NS	flipped      	74	0.0%
# NS	swapped      	1295	0.4%
# NS	flip+swap    	793	0.2%
# NS	unresolved   	1	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	1
# NS	non-ACGT     	1
# NS	non-SNP      	0
# NS	non-biallelic	0
# cedars_610k_old_gwas
# ST	A>G	113047	19.4%
# ST	A>T	583	0.1%
# ST	C>A	27989	4.8%
# ST	C>G	889	0.2%
# ST	C>T	122718	21.0%
# ST	G>A	122654	21.0%
# ST	G>C	928	0.2%
# ST	G>T	27742	4.8%
# ST	T>A	590	0.1%
# ST	T>C	113708	19.5%
# ST	T>G	26386	4.5%
# # NS, Number of sites:
# NS	total        	583671
# NS	ref match    	581442	99.6%
# NS	ref mismatch 	2227	0.4%
# NS	flipped      	120	0.0%
# NS	swapped      	1975	0.3%
# NS	flip+swap    	1600	0.3%
# NS	unresolved   	4	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	2
# NS	non-ACGT     	2
# NS	non-SNP      	0
# NS	non-biallelic	0
# cedars_omni_old_gwas
# ST	A>G	144804	20.1%
# ST	A>T	255	0.0%
# ST	C>A	33944	4.7%
# ST	C>G	317	0.0%
# ST	C>T	156359	21.8%
# ST	G>A	144739	20.1%
# ST	G>C	353	0.0%
# ST	G>T	33904	4.7%
# ST	T>A	137	0.0%
# ST	T>C	137454	19.1%
# ST	T>G	33035	4.6%
# # NS, Number of sites:
# NS	total        	718889
# NS	ref match    	716484	99.7%
# NS	ref mismatch 	2405	0.3%
# NS	flipped      	130	0.0%
# NS	swapped      	2119	0.3%
# NS	flip+swap    	577	0.1%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# swedish_uc_old_gwas
# ST	A>G	58492	19.7%
# ST	A>T	0	0.0%
# ST	C>A	14123	4.8%
# ST	C>G	0	0.0%
# ST	C>T	62172	21.0%
# ST	G>A	62411	21.0%
# ST	G>C	0	0.0%
# ST	G>T	14026	4.7%
# ST	T>A	0	0.0%
# ST	T>C	58584	19.8%
# ST	T>G	13388	4.5%
# # NS, Number of sites:
# NS	total        	296606
# NS	ref match    	295668	99.7%
# NS	ref mismatch 	938	0.3%
# NS	flipped      	41	0.0%
# NS	swapped      	844	0.3%
# NS	flip+swap    	53	0.0%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# mccauley_gsa
# ST	A>G	96886	17.4%
# ST	A>T	733	0.1%
# ST	C>A	28020	5.0%
# ST	C>G	1181	0.2%
# ST	C>T	127729	23.0%
# ST	G>A	127567	22.9%
# ST	G>C	1152	0.2%
# ST	G>T	27751	5.0%
# ST	T>A	773	0.1%
# ST	T>C	96961	17.4%
# ST	T>G	23641	4.3%
# # NS, Number of sites:
# NS	total        	556193
# NS	ref match    	554556	99.7%
# NS	ref mismatch 	1637	0.3%
# NS	flipped      	139	0.0%
# NS	swapped      	1348	0.2%
# NS	flip+swap    	2092	0.4%
# NS	unresolved   	5	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# ccfa_gsa
# ST	A>G	102768	17.1%
# ST	A>T	943	0.2%
# ST	C>A	30592	5.1%
# ST	C>G	1498	0.2%
# ST	C>T	138680	23.1%
# ST	G>A	138459	23.1%
# ST	G>C	1461	0.2%
# ST	G>T	30197	5.0%
# ST	T>A	982	0.2%
# ST	T>C	103102	17.2%
# ST	T>G	25314	4.2%
# # NS, Number of sites:
# NS	total        	599467
# NS	ref match    	597734	99.7%
# NS	ref mismatch 	1731	0.3%
# NS	flipped      	174	0.0%
# NS	swapped      	1385	0.2%
# NS	flip+swap    	2605	0.4%
# NS	unresolved   	2	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	2
# NS	non-ACGT     	2
# NS	non-SNP      	0
# NS	non-biallelic	0
# cedars_gsa
# ST	A>G	104189	17.1%
# ST	A>T	943	0.2%
# ST	C>A	31158	5.1%
# ST	C>G	1529	0.3%
# ST	C>T	141601	23.2%
# ST	G>A	141231	23.2%
# ST	G>C	1525	0.3%
# ST	G>T	30868	5.1%
# ST	T>A	939	0.2%
# ST	T>C	104643	17.2%
# ST	T>G	25523	4.2%
# # NS, Number of sites:
# NS	total        	609903
# NS	ref match    	608134	99.7%
# NS	ref mismatch 	1765	0.3%
# NS	flipped      	184	0.0%
# NS	swapped      	1389	0.2%
# NS	flip+swap    	2665	0.4%
# NS	unresolved   	3	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	4
# NS	non-ACGT     	4
# NS	non-SNP      	0
# NS	non-biallelic	0
# bernstein_gsa
# ST	A>G	96786	17.5%
# ST	A>T	721	0.1%
# ST	C>A	27866	5.0%
# ST	C>G	1078	0.2%
# ST	C>T	126697	22.9%
# ST	G>A	126281	22.8%
# ST	G>C	1022	0.2%
# ST	G>T	27563	5.0%
# ST	T>A	730	0.1%
# ST	T>C	97002	17.5%
# ST	T>G	23628	4.3%
# # NS, Number of sites:
# NS	total        	553163
# NS	ref match    	551545	99.7%
# NS	ref mismatch 	1618	0.3%
# NS	flipped      	131	0.0%
# NS	swapped      	1340	0.2%
# NS	flip+swap    	1944	0.4%
# NS	unresolved   	2	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# farkkila_gsa
# ST	A>G	87681	17.7%
# ST	A>T	581	0.1%
# ST	C>A	24796	5.0%
# ST	C>G	822	0.2%
# ST	C>T	112065	22.7%
# ST	G>A	111728	22.6%
# ST	G>C	806	0.2%
# ST	G>T	24629	5.0%
# ST	T>A	586	0.1%
# ST	T>C	87994	17.8%
# ST	T>G	21251	4.3%
# # NS, Number of sites:
# NS	total        	494344
# NS	ref match    	492810	99.7%
# NS	ref mismatch 	1534	0.3%
# NS	flipped      	118	0.0%
# NS	swapped      	1290	0.3%
# NS	flip+swap    	1528	0.3%
# NS	unresolved   	4	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# franchimont_gsa
# ST	A>G	98080	17.1%
# ST	A>T	852	0.1%
# ST	C>A	29026	5.1%
# ST	C>G	1408	0.2%
# ST	C>T	133019	23.2%
# ST	G>A	132652	23.2%
# ST	G>C	1380	0.2%
# ST	G>T	28736	5.0%
# ST	T>A	878	0.2%
# ST	T>C	98282	17.2%
# ST	T>G	23953	4.2%
# # NS, Number of sites:
# NS	total        	572382
# NS	ref match    	570719	99.7%
# NS	ref mismatch 	1662	0.3%
# NS	flipped      	144	0.0%
# NS	swapped      	1350	0.2%
# NS	flip+swap    	2415	0.4%
# NS	unresolved   	5	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	1
# NS	non-ACGT     	1
# NS	non-SNP      	0
# NS	non-biallelic	0
# franke_gsa
# ST	A>G	97110	17.4%
# ST	A>T	715	0.1%
# ST	C>A	28046	5.0%
# ST	C>G	1126	0.2%
# ST	C>T	127809	23.0%
# ST	G>A	127520	22.9%
# ST	G>C	1059	0.2%
# ST	G>T	27723	5.0%
# ST	T>A	752	0.1%
# ST	T>C	97331	17.5%
# ST	T>G	23750	4.3%
# # NS, Number of sites:
# NS	total        	556807
# NS	ref match    	555194	99.7%
# NS	ref mismatch 	1613	0.3%
# NS	flipped      	139	0.0%
# NS	swapped      	1331	0.2%
# NS	flip+swap    	1967	0.4%
# NS	unresolved   	1	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# helmsley_prism_gsa
# ST	A>G	47901	19.7%
# ST	A>T	68	0.0%
# ST	C>A	11657	4.8%
# ST	C>G	117	0.0%
# ST	C>T	50671	20.8%
# ST	G>A	50543	20.8%
# ST	G>C	104	0.0%
# ST	G>T	11785	4.8%
# ST	T>A	77	0.0%
# ST	T>C	47675	19.6%
# ST	T>G	11125	4.6%
# # NS, Number of sites:
# NS	total        	243031
# NS	ref match    	242069	99.6%
# NS	ref mismatch 	962	0.4%
# NS	flipped      	41	0.0%
# NS	swapped      	869	0.4%
# NS	flip+swap    	231	0.1%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# helmsley_xavier_prism_gsa
# ST	A>G	47908	19.7%
# ST	A>T	68	0.0%
# ST	C>A	11665	4.8%
# ST	C>G	118	0.0%
# ST	C>T	50669	20.8%
# ST	G>A	50544	20.8%
# ST	G>C	109	0.0%
# ST	G>T	11780	4.8%
# ST	T>A	76	0.0%
# ST	T>C	47662	19.6%
# ST	T>G	11129	4.6%
# # NS, Number of sites:
# NS	total        	243035
# NS	ref match    	242073	99.6%
# NS	ref mismatch 	962	0.4%
# NS	flipped      	39	0.0%
# NS	swapped      	871	0.4%
# NS	flip+swap    	233	0.1%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# hyams_protect_gsa
# ST	A>G	98902	17.4%
# ST	A>T	857	0.2%
# ST	C>A	28787	5.1%
# ST	C>G	1218	0.2%
# ST	C>T	130696	23.0%
# ST	G>A	130422	22.9%
# ST	G>C	1168	0.2%
# ST	G>T	28499	5.0%
# ST	T>A	890	0.2%
# ST	T>C	99061	17.4%
# ST	T>G	24225	4.3%
# # NS, Number of sites:
# NS	total        	569172
# NS	ref match    	567515	99.7%
# NS	ref mismatch 	1657	0.3%
# NS	flipped      	143	0.0%
# NS	swapped      	1357	0.2%
# NS	flip+swap    	2218	0.4%
# NS	unresolved   	3	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# lewis_sparc_gsa
# ST	A>G	103473	17.0%
# ST	A>T	1037	0.2%
# ST	C>A	31056	5.1%
# ST	C>G	1625	0.3%
# ST	C>T	140917	23.2%
# ST	G>A	140748	23.2%
# ST	G>C	1618	0.3%
# ST	G>T	30632	5.0%
# ST	T>A	1055	0.2%
# ST	T>C	103784	17.1%
# ST	T>G	25509	4.2%
# # NS, Number of sites:
# NS	total        	607141
# NS	ref match    	605399	99.7%
# NS	ref mismatch 	1740	0.3%
# NS	flipped      	176	0.0%
# NS	swapped      	1390	0.2%
# NS	flip+swap    	2834	0.5%
# NS	unresolved   	3	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	2
# NS	non-ACGT     	2
# NS	non-SNP      	0
# NS	non-biallelic	0
# mccauley_new_gsa
# ST	A>G	99669	17.2%
# ST	A>T	853	0.1%
# ST	C>A	29209	5.1%
# ST	C>G	1425	0.2%
# ST	C>T	133414	23.1%
# ST	G>A	133394	23.1%
# ST	G>C	1413	0.2%
# ST	G>T	28986	5.0%
# ST	T>A	897	0.2%
# ST	T>C	99755	17.3%
# ST	T>G	24407	4.2%
# # NS, Number of sites:
# NS	total        	578007
# NS	ref match    	576321	99.7%
# NS	ref mismatch 	1686	0.3%
# NS	flipped      	153	0.0%
# NS	swapped      	1371	0.2%
# NS	flip+swap    	2469	0.4%
# NS	unresolved   	5	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# mcgovern_gsa
# ST	A>G	105936	16.9%
# ST	A>T	1053	0.2%
# ST	C>A	32168	5.1%
# ST	C>G	1766	0.3%
# ST	C>T	146431	23.4%
# ST	G>A	146122	23.3%
# ST	G>C	1739	0.3%
# ST	G>T	31807	5.1%
# ST	T>A	1057	0.2%
# ST	T>C	106303	17.0%
# ST	T>G	26047	4.2%
# # NS, Number of sites:
# NS	total        	626745
# NS	ref match    	624946	99.7%
# NS	ref mismatch 	1795	0.3%
# NS	flipped      	196	0.0%
# NS	swapped      	1397	0.2%
# NS	flip+swap    	2990	0.5%
# NS	unresolved   	3	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	4
# NS	non-ACGT     	4
# NS	non-SNP      	0
# NS	non-biallelic	0
# moayyedi_imagine_gsa
# ST	A>G	102335	17.4%
# ST	A>T	803	0.1%
# ST	C>A	29922	5.1%
# ST	C>G	1275	0.2%
# ST	C>T	135017	22.9%
# ST	G>A	134867	22.9%
# ST	G>C	1224	0.2%
# ST	G>T	29565	5.0%
# ST	T>A	840	0.1%
# ST	T>C	102586	17.4%
# ST	T>G	25136	4.3%
# # NS, Number of sites:
# NS	total        	588817
# NS	ref match    	587130	99.7%
# NS	ref mismatch 	1687	0.3%
# NS	flipped      	154	0.0%
# NS	swapped      	1363	0.2%
# NS	flip+swap    	2249	0.4%
# NS	unresolved   	3	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# newberry_share_gsa
# ST	A>G	99878	17.4%
# ST	A>T	829	0.1%
# ST	C>A	29136	5.1%
# ST	C>G	1300	0.2%
# ST	C>T	132226	23.0%
# ST	G>A	131916	22.9%
# ST	G>C	1271	0.2%
# ST	G>T	28868	5.0%
# ST	T>A	852	0.1%
# ST	T>C	100062	17.4%
# ST	T>G	24518	4.3%
# # NS, Number of sites:
# NS	total        	575543
# NS	ref match    	573889	99.7%
# NS	ref mismatch 	1654	0.3%
# NS	flipped      	146	0.0%
# NS	swapped      	1351	0.2%
# NS	flip+swap    	2296	0.4%
# NS	unresolved   	3	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# niddk_cho_gsa
# ST	A>G	101571	17.1%
# ST	A>T	925	0.2%
# ST	C>A	30168	5.1%
# ST	C>G	1469	0.2%
# ST	C>T	137370	23.2%
# ST	G>A	137050	23.1%
# ST	G>C	1418	0.2%
# ST	G>T	29994	5.1%
# ST	T>A	937	0.2%
# ST	T>C	101803	17.2%
# ST	T>G	25013	4.2%
# # NS, Number of sites:
# NS	total        	592900
# NS	ref match    	591209	99.7%
# NS	ref mismatch 	1689	0.3%
# NS	flipped      	163	0.0%
# NS	swapped      	1361	0.2%
# NS	flip+swap    	2543	0.4%
# NS	unresolved   	4	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	2
# NS	non-ACGT     	2
# NS	non-SNP      	0
# NS	non-biallelic	0
# niddk_duerr_gsa
# ST	A>G	99757	17.2%
# ST	A>T	912	0.2%
# ST	C>A	29479	5.1%
# ST	C>G	1444	0.2%
# ST	C>T	134357	23.1%
# ST	G>A	134098	23.1%
# ST	G>C	1390	0.2%
# ST	G>T	29189	5.0%
# ST	T>A	934	0.2%
# ST	T>C	100011	17.2%
# ST	T>G	24543	4.2%
# # NS, Number of sites:
# NS	total        	580761
# NS	ref match    	579086	99.7%
# NS	ref mismatch 	1675	0.3%
# NS	flipped      	159	0.0%
# NS	swapped      	1355	0.2%
# NS	flip+swap    	2512	0.4%
# NS	unresolved   	3	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# niddk_rioux_gsa
# ST	A>G	98615	17.3%
# ST	A>T	819	0.1%
# ST	C>A	28832	5.1%
# ST	C>G	1222	0.2%
# ST	C>T	131160	23.0%
# ST	G>A	131081	23.0%
# ST	G>C	1210	0.2%
# ST	G>T	28619	5.0%
# ST	T>A	822	0.1%
# ST	T>C	98841	17.3%
# ST	T>G	24178	4.2%
# # NS, Number of sites:
# NS	total        	569782
# NS	ref match    	568124	99.7%
# NS	ref mismatch 	1656	0.3%
# NS	flipped      	149	0.0%
# NS	swapped      	1346	0.2%
# NS	flip+swap    	2199	0.4%
# NS	unresolved   	3	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	2
# NS	non-ACGT     	2
# NS	non-SNP      	0
# NS	non-biallelic	0
# niddk_silverberg_gsa
# ST	A>G	102757	17.1%
# ST	A>T	944	0.2%
# ST	C>A	30614	5.1%
# ST	C>G	1503	0.3%
# ST	C>T	139349	23.2%
# ST	G>A	139110	23.2%
# ST	G>C	1463	0.2%
# ST	G>T	30379	5.1%
# ST	T>A	955	0.2%
# ST	T>C	102915	17.1%
# ST	T>G	25274	4.2%
# # NS, Number of sites:
# NS	total        	600699
# NS	ref match    	598974	99.7%
# NS	ref mismatch 	1722	0.3%
# NS	flipped      	172	0.0%
# NS	swapped      	1367	0.2%
# NS	flip+swap    	2610	0.4%
# NS	unresolved   	5	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	3
# NS	non-ACGT     	3
# NS	non-SNP      	0
# NS	non-biallelic	0
# palotie_hus_gsa
# ST	A>G	96616	17.5%
# ST	A>T	715	0.1%
# ST	C>A	27852	5.0%
# ST	C>G	1083	0.2%
# ST	C>T	126755	22.9%
# ST	G>A	126371	22.9%
# ST	G>C	1037	0.2%
# ST	G>T	27576	5.0%
# ST	T>A	751	0.1%
# ST	T>C	96792	17.5%
# ST	T>G	23621	4.3%
# # NS, Number of sites:
# NS	total        	552923
# NS	ref match    	551291	99.7%
# NS	ref mismatch 	1631	0.3%
# NS	flipped      	141	0.0%
# NS	swapped      	1343	0.2%
# NS	flip+swap    	1951	0.4%
# NS	unresolved   	3	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	1
# NS	non-ACGT     	1
# NS	non-SNP      	0
# NS	non-biallelic	0
# pekow_share_gsa
# ST	A>G	99480	17.4%
# ST	A>T	806	0.1%
# ST	C>A	29032	5.1%
# ST	C>G	1223	0.2%
# ST	C>T	131394	22.9%
# ST	G>A	131073	22.9%
# ST	G>C	1229	0.2%
# ST	G>T	28690	5.0%
# ST	T>A	817	0.1%
# ST	T>C	99822	17.4%
# ST	T>G	24453	4.3%
# # NS, Number of sites:
# NS	total        	572621
# NS	ref match    	570955	99.7%
# NS	ref mismatch 	1666	0.3%
# NS	flipped      	150	0.0%
# NS	swapped      	1359	0.2%
# NS	flip+swap    	2223	0.4%
# NS	unresolved   	2	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# rioux_igenomed_gsa
# ST	A>G	95483	17.6%
# ST	A>T	686	0.1%
# ST	C>A	27298	5.0%
# ST	C>G	978	0.2%
# ST	C>T	123839	22.8%
# ST	G>A	123488	22.7%
# ST	G>C	947	0.2%
# ST	G>T	27032	5.0%
# ST	T>A	703	0.1%
# ST	T>C	95618	17.6%
# ST	T>G	23341	4.3%
# # NS, Number of sites:
# NS	total        	542911
# NS	ref match    	541310	99.7%
# NS	ref mismatch 	1601	0.3%
# NS	flipped      	127	0.0%
# NS	swapped      	1334	0.2%
# NS	flip+swap    	1802	0.3%
# NS	unresolved   	2	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# sands_msccr_gsa
# ST	A>G	102644	17.3%
# ST	A>T	854	0.1%
# ST	C>A	30207	5.1%
# ST	C>G	1367	0.2%
# ST	C>T	136955	23.0%
# ST	G>A	136689	23.0%
# ST	G>C	1400	0.2%
# ST	G>T	29883	5.0%
# ST	T>A	900	0.2%
# ST	T>C	102878	17.3%
# ST	T>G	25251	4.2%
# # NS, Number of sites:
# NS	total        	594488
# NS	ref match    	592779	99.7%
# NS	ref mismatch 	1709	0.3%
# NS	flipped      	157	0.0%
# NS	swapped      	1378	0.2%
# NS	flip+swap    	2456	0.4%
# NS	unresolved   	3	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# stampfer_gsa
# ST	A>G	98037	17.2%
# ST	A>T	843	0.1%
# ST	C>A	28928	5.1%
# ST	C>G	1394	0.2%
# ST	C>T	132527	23.2%
# ST	G>A	132210	23.2%
# ST	G>C	1366	0.2%
# ST	G>T	28708	5.0%
# ST	T>A	866	0.2%
# ST	T>C	98158	17.2%
# ST	T>G	23903	4.2%
# # NS, Number of sites:
# NS	total        	571027
# NS	ref match    	569384	99.7%
# NS	ref mismatch 	1643	0.3%
# NS	flipped      	152	0.0%
# NS	swapped      	1335	0.2%
# NS	flip+swap    	2390	0.4%
# NS	unresolved   	1	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# vermeire_gsa
# ST	A>G	99559	17.0%
# ST	A>T	881	0.2%
# ST	C>A	29689	5.1%
# ST	C>G	1470	0.3%
# ST	C>T	136142	23.3%
# ST	G>A	135975	23.3%
# ST	G>C	1457	0.2%
# ST	G>T	29462	5.0%
# ST	T>A	911	0.2%
# ST	T>C	99846	17.1%
# ST	T>G	24348	4.2%
# # NS, Number of sites:
# NS	total        	584265
# NS	ref match    	582564	99.7%
# NS	ref mismatch 	1697	0.3%
# NS	flipped      	160	0.0%
# NS	swapped      	1358	0.2%
# NS	flip+swap    	2541	0.4%
# NS	unresolved   	5	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	4
# NS	non-ACGT     	4
# NS	non-SNP      	0
# NS	non-biallelic	0
# weersma_gsa
# ST	A>G	97924	17.5%
# ST	A>T	733	0.1%
# ST	C>A	28356	5.1%
# ST	C>G	1125	0.2%
# ST	C>T	128136	22.9%
# ST	G>A	127820	22.8%
# ST	G>C	1072	0.2%
# ST	G>T	27948	5.0%
# ST	T>A	749	0.1%
# ST	T>C	98179	17.5%
# ST	T>G	24021	4.3%
# # NS, Number of sites:
# NS	total        	560260
# NS	ref match    	558630	99.7%
# NS	ref mismatch 	1630	0.3%
# NS	flipped      	133	0.0%
# NS	swapped      	1348	0.2%
# NS	flip+swap    	2001	0.4%
# NS	unresolved   	3	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# xavier_prism_gsa
# ST	A>G	98708	17.4%
# ST	A>T	770	0.1%
# ST	C>A	28707	5.1%
# ST	C>G	1200	0.2%
# ST	C>T	130469	23.0%
# ST	G>A	130226	22.9%
# ST	G>C	1207	0.2%
# ST	G>T	28469	5.0%
# ST	T>A	811	0.1%
# ST	T>C	98905	17.4%
# ST	T>G	24188	4.3%
# # NS, Number of sites:
# NS	total        	567988
# NS	ref match    	566358	99.7%
# NS	ref mismatch 	1630	0.3%
# NS	flipped      	143	0.0%
# NS	swapped      	1338	0.2%
# NS	flip+swap    	2168	0.4%
# NS	unresolved   	2	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# xavier_share_gsa
# ST	A>G	98438	17.3%
# ST	A>T	793	0.1%
# ST	C>A	28683	5.1%
# ST	C>G	1225	0.2%
# ST	C>T	130763	23.0%
# ST	G>A	130412	23.0%
# ST	G>C	1265	0.2%
# ST	G>T	28491	5.0%
# ST	T>A	817	0.1%
# ST	T>C	98737	17.4%
# ST	T>G	24068	4.2%
# # NS, Number of sites:
# NS	total        	567891
# NS	ref match    	566254	99.7%
# NS	ref mismatch 	1637	0.3%
# NS	flipped      	141	0.0%
# NS	swapped      	1344	0.2%
# NS	flip+swap    	2210	0.4%
# NS	unresolved   	4	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0

##############

# REMOVE MISSMATCH ALLELES, THEY  IMPUTATION SERVER TO CRASH:

for i in ${studies[@]}
do
bsub -J"bcf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_bcftools_4_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_bcftools_4_${i} \
"/path/to/software/bcftools-1.16/./bcftools norm \
--check-ref x -f ${path_gwas}resources/hg38/hg38_edited.fa \
${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned.vcf.gz \
-Oz -o ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_2.vcf.gz"
done

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
for i in ${studies[@]}
do
echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_bcftools_4_${i} | grep "Successfully completed."
done

for i in ${studies[@]}
do
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_2.vcf.gz
done

for i in ${studies[@]}
do
echo ${i} && tail -24 ${path_gwas}pre_imputation/QC/${i}/logs/stderr_bcftools_4_${i}
done

##############
# all_hce
# Lines   total/split/realigned/skipped:	421656/0/0/19
# niddk_old_gwas
# Lines   total/split/realigned/skipped:	300108/0/0/0
# australia_omniexome
# Lines   total/split/realigned/skipped:	738250/0/0/12
# gwas1
# Lines   total/split/realigned/skipped:	445760/0/0/22
# gwas2
# Lines   total/split/realigned/skipped:	772572/0/0/46
# pittsburgh_gsa
# Lines   total/split/realigned/skipped:	891011/0/0/23
# spain_gsa
# Lines   total/split/realigned/skipped:	575671/0/0/5
# italy_gsa
# Lines   total/split/realigned/skipped:	560438/0/0/1
# kiel_austria_sibdcs_gsa
# Lines   total/split/realigned/skipped:	544252/0/0/1
# netherlands_gsa
# Lines   total/split/realigned/skipped:	584603/0/0/1
# slovenia_gsa
# Lines   total/split/realigned/skipped:	523957/0/0/0
# sweden_gsa
# Lines   total/split/realigned/skipped:	555288/0/0/1
# niddk_broad_gsa
# Lines   total/split/realigned/skipped:	588821/0/0/2
# niddk_feinstein_gsa
# Lines   total/split/realigned/skipped:	600728/0/0/5
# basque_gsa
# Lines   total/split/realigned/skipped:	559956/0/0/1
# lithuania_gsa
# Lines   total/split/realigned/skipped:	553397/0/0/7
# belgium_louis_gsa
# Lines   total/split/realigned/skipped:	562724/0/0/3
# belgium_franchimont_gsa
# Lines   total/split/realigned/skipped:	553249/0/0/4
# belgium_vermeire_gsa
# Lines   total/split/realigned/skipped:	636647/0/0/10
# prism_nfe_gsa
# Lines   total/split/realigned/skipped:	558477/0/0/2
# prism_nfe_gwas
# Lines   total/split/realigned/skipped:	242848/0/0/0
# finland_illugwas
# Lines   total/split/realigned/skipped:	235986/0/0/0
# german_affy6_old_gwas
# Lines   total/split/realigned/skipped:	782252/0/0/43
# norway_affy6_old_gwas
# Lines   total/split/realigned/skipped:	723266/0/0/43
# belgium_inf1_old_gwas
# Lines   total/split/realigned/skipped:	301665/0/0/0
# belgium_inf2_old_gwas
# Lines   total/split/realigned/skipped:	290145/0/0/0
# cedars_370k_old_gwas
# Lines   total/split/realigned/skipped:	338634/0/0/2
# cedars_610k_old_gwas
# Lines   total/split/realigned/skipped:	583671/0/0/6
# cedars_omni_old_gwas
# Lines   total/split/realigned/skipped:	718889/0/0/0
# swedish_uc_old_gwas
# Lines   total/split/realigned/skipped:	296606/0/0/0
# mccauley_gsa
# Lines   total/split/realigned/skipped:	556193/0/0/5
# ccfa_gsa
# Lines   total/split/realigned/skipped:	599467/0/0/4
# cedars_gsa
# Lines   total/split/realigned/skipped:	609903/0/0/7
# bernstein_gsa
# Lines   total/split/realigned/skipped:	553163/0/0/2
# farkkila_gsa
# Lines   total/split/realigned/skipped:	494344/0/0/4
# franchimont_gsa
# Lines   total/split/realigned/skipped:	572382/0/0/6
# franke_gsa
# Lines   total/split/realigned/skipped:	556807/0/0/1
# helmsley_prism_gsa
# Lines   total/split/realigned/skipped:	243031/0/0/0
# helmsley_xavier_prism_gsa
# Lines   total/split/realigned/skipped:	243035/0/0/0
# hyams_protect_gsa
# Lines   total/split/realigned/skipped:	569172/0/0/3
# lewis_sparc_gsa
# Lines   total/split/realigned/skipped:	607141/0/0/5
# mccauley_new_gsa
# Lines   total/split/realigned/skipped:	578007/0/0/5
# mcgovern_gsa
# Lines   total/split/realigned/skipped:	626745/0/0/7
# moayyedi_imagine_gsa
# Lines   total/split/realigned/skipped:	588817/0/0/3
# newberry_share_gsa
# Lines   total/split/realigned/skipped:	575543/0/0/3
# niddk_cho_gsa
# Lines   total/split/realigned/skipped:	592900/0/0/6
# niddk_duerr_gsa
# Lines   total/split/realigned/skipped:	580761/0/0/3
# niddk_rioux_gsa
# Lines   total/split/realigned/skipped:	569782/0/0/5
# niddk_silverberg_gsa
# Lines   total/split/realigned/skipped:	600699/0/0/8
# palotie_hus_gsa
# Lines   total/split/realigned/skipped:	552923/0/0/4
# pekow_share_gsa
# Lines   total/split/realigned/skipped:	572621/0/0/2
# rioux_igenomed_gsa
# Lines   total/split/realigned/skipped:	542911/0/0/2
# sands_msccr_gsa
# Lines   total/split/realigned/skipped:	594488/0/0/3
# stampfer_gsa
# Lines   total/split/realigned/skipped:	571027/0/0/1
# vermeire_gsa
# Lines   total/split/realigned/skipped:	584265/0/0/9
# weersma_gsa
# Lines   total/split/realigned/skipped:	560260/0/0/3
# xavier_prism_gsa
# Lines   total/split/realigned/skipped:	567988/0/0/2
# xavier_share_gsa
# Lines   total/split/realigned/skipped:	567891/0/0/4
##############

##########

for i in ${studies[@]}
do
bsub -J"bcf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_bcftools_5_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_bcftools_5_${i} \
"/path/to/software/bcftools-1.16/./bcftools +fixref \
${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_2.vcf.gz \
-- -f ${path_gwas}resources/hg38/hg38_edited.fa"
done

for i in ${studies[@]}
do
echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_bcftools_5_${i} | grep "Successfully completed."
done

for i in ${studies[@]}
do
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_2.vcf.gz
done

for i in ${studies[@]}
do
echo ${i} && tail -24 ${path_gwas}pre_imputation/QC/${i}/logs/stderr_bcftools_5_${i} | grep "ref match"
done

########
# all_hce
# NS	ref match    	421637	100.0%
# niddk_old_gwas
# NS	ref match    	300108	100.0%
# australia_omniexome
# NS	ref match    	738238	100.0%
# gwas1
# NS	ref match    	445738	100.0%
# gwas2
# NS	ref match    	772526	100.0%
# pittsburgh_gsa
# NS	ref match    	890988	100.0%
# spain_gsa
# NS	ref match    	575666	100.0%
# italy_gsa
# NS	ref match    	560437	100.0%
# kiel_austria_sibdcs_gsa
# NS	ref match    	544251	100.0%
# netherlands_gsa
# NS	ref match    	584602	100.0%
# slovenia_gsa
# NS	ref match    	523957	100.0%
# sweden_gsa
# NS	ref match    	555287	100.0%
# niddk_broad_gsa
# NS	ref match    	588819	100.0%
# niddk_feinstein_gsa
# NS	ref match    	600723	100.0%
# basque_gsa
# NS	ref match    	559955	100.0%
# lithuania_gsa
# NS	ref match    	553390	100.0%
# belgium_louis_gsa
# NS	ref match    	562721	100.0%
# belgium_franchimont_gsa
# NS	ref match    	553245	100.0%
# belgium_vermeire_gsa
# NS	ref match    	636637	100.0%
# prism_nfe_gsa
# NS	ref match    	558475	100.0%
# prism_nfe_gwas
# NS	ref match    	242848	100.0%
# finland_illugwas
# NS	ref match    	235986	100.0%
# german_affy6_old_gwas
# NS	ref match    	782209	100.0%
# norway_affy6_old_gwas
# NS	ref match    	723223	100.0%
# belgium_inf1_old_gwas
# NS	ref match    	301665	100.0%
# belgium_inf2_old_gwas
# NS	ref match    	290145	100.0%
# cedars_370k_old_gwas
# NS	ref match    	338632	100.0%
# cedars_610k_old_gwas
# NS	ref match    	583665	100.0%
# cedars_omni_old_gwas
# NS	ref match    	718889	100.0%
# swedish_uc_old_gwas
# NS	ref match    	296606	100.0%
# mccauley_gsa
# NS	ref match    	556188	100.0%
# ccfa_gsa
# NS	ref match    	599463	100.0%
# cedars_gsa
# NS	ref match    	609896	100.0%
# bernstein_gsa
# NS	ref match    	553161	100.0%
# farkkila_gsa
# NS	ref match    	494340	100.0%
# franchimont_gsa
# NS	ref match    	572376	100.0%
# franke_gsa
# NS	ref match    	556806	100.0%
# helmsley_prism_gsa
# NS	ref match    	243031	100.0%
# helmsley_xavier_prism_gsa
# NS	ref match    	243035	100.0%
# hyams_protect_gsa
# NS	ref match    	569169	100.0%
# lewis_sparc_gsa
# NS	ref match    	607136	100.0%
# mccauley_new_gsa
# NS	ref match    	578002	100.0%
# mcgovern_gsa
# NS	ref match    	626738	100.0%
# moayyedi_imagine_gsa
# NS	ref match    	588814	100.0%
# newberry_share_gsa
# NS	ref match    	575540	100.0%
# niddk_cho_gsa
# NS	ref match    	592894	100.0%
# niddk_duerr_gsa
# NS	ref match    	580758	100.0%
# niddk_rioux_gsa
# NS	ref match    	569777	100.0%
# niddk_silverberg_gsa
# NS	ref match    	600691	100.0%
# palotie_hus_gsa
# NS	ref match    	552919	100.0%
# pekow_share_gsa
# NS	ref match    	572619	100.0%
# rioux_igenomed_gsa
# NS	ref match    	542909	100.0%
# sands_msccr_gsa
# NS	ref match    	594485	100.0%
# stampfer_gsa
# NS	ref match    	571026	100.0%
# vermeire_gsa
# NS	ref match    	584256	100.0%
# weersma_gsa
# NS	ref match    	560257	100.0%
# xavier_prism_gsa
# NS	ref match    	567986	100.0%
# xavier_share_gsa
# NS	ref match    	567887	100.0%
########

#### VCF to BED

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=1000

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_58_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_58_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--vcf ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_2.vcf.gz \
--keep-allele-order --allow-no-sex \
--double-id \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_58_${i} | grep -E "completed"
done

for i in ${studies[@]}
do echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned.fam
done

##############################################################
# 16.2 UPDATE NAME VARIANTS TO CHR:POSITION_REF_ALT in b38

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

MEM=500

for i in ${studies[@]}
do 
bsub -J"b38" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_b38_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_${i} \
"zcat ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_2.vcf.gz | cut -f '1-5' | awk '{print \$3,\$1\":\"\$2\"_\"\$4\"_\"\$5}' \
> ${path_gwas}pre_imputation/QC/${i}/list_variants_all_hce_hg38_posstrandaligned"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_${i} | grep -E "completed"
done

for i in ${studies[@]}
do echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/list_variants_all_hce_hg38_posstrandaligned
done


###  /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

# all_hce
# gwas2
# kiel_austria_sibdcs_gsa
# german_affy6_old_gwas

cohorts<-c("all_hce","gwas2","kiel_austria_sibdcs_gsa","german_affy6_old_gwas")

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("australia_omniexome","gwas1"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

length(cohorts)
# [1] 58

for (j in 1:length(cohorts)) {
  
  print(cohorts[j])
  
  bim<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned.bim",sep=""),sep="\t",head=F)
  ids<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/list_variants_all_hce_hg38_posstrandaligned",sep=""),sep=" ",head=F,skip=30)
  ids<-ids[which(ids$V2!=""),]
  
  colnames(ids)[2]<-"ids"
  ids$ids<-as.character(ids$ids)
  bim$V2<-as.character(bim$V2)
  
  print(table(bim$V2==ids$V1))
  
  #the major allele is set to A2 by default by Plink, keep ids with real ref/alt as in vcf using the ids file
  bim.1<-cbind(bim,ids[,"ids",drop=F])
  
  bim.1$ids<-gsub("X:","23:",bim.1$ids)
  write.table(bim.1[,c(1,7,3:6)],paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_edited.bim",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
  rm(list=ls()[!ls() %in% c("path","j","cohorts")])
}

##############
# [1] "kiel_austria_sibdcs_gsa"
# TRUE 
# 544251 
# [1] "australia_omniexome"
# TRUE 
# 738238 
# [1] "gwas1"
# TRUE 
# 445738 
# [1] "gwas2"
# TRUE 
# 772526 
# [1] "all_hce"
# TRUE 
# 421637 
# [1] "pittsburgh_gsa"
# TRUE 
# 890988 
# [1] "spain_gsa"
# TRUE 
# 575666 
# [1] "italy_gsa"
# TRUE 
# 560437 
# [1] "netherlands_gsa"
# TRUE 
# 584602 
# [1] "slovenia_gsa"
# TRUE 
# 523957 
# [1] "sweden_gsa"
# TRUE 
# 555287 
# [1] "niddk_broad_gsa"
# TRUE 
# 588819 
# [1] "niddk_feinstein_gsa"
# TRUE 
# 600723 
# [1] "basque_gsa"
# TRUE 
# 559955 
# [1] "prism_nfe_gsa"
# TRUE 
# 558475 
# [1] "lithuania_gsa"
# TRUE 
# 553390 
# [1] "belgium_louis_gsa"
# TRUE 
# 562721 
# [1] "belgium_franchimont_gsa"
# TRUE 
# 553245 
# [1] "belgium_vermeire_gsa"
# TRUE 
# 636637 
# [1] "prism_nfe_gwas"
# TRUE 
# 242848 
# [1] "finland_illugwas"
# TRUE 
# 235986 
# [1] "german_affy6_old_gwas"
# TRUE 
# 782209 
# [1] "norway_affy6_old_gwas"
# TRUE 
# 723223 
# [1] "belgium_inf1_old_gwas"
# TRUE 
# 301665 
# [1] "belgium_inf2_old_gwas"
# TRUE 
# 290145 
# [1] "niddk_old_gwas"
# TRUE 
# 300108 
# [1] "cedars_370k_old_gwas"
# TRUE 
# 338632 
# [1] "cedars_610k_old_gwas"
# TRUE 
# 583665 
# [1] "cedars_omni_old_gwas"
# TRUE 
# 718889 
# [1] "swedish_uc_old_gwas"
# TRUE 
# 296606 
# [1] "mccauley_gsa"
# TRUE 
# 556188 
# [1] "ccfa_gsa"
# TRUE 
# 599463 
# [1] "cedars_gsa"
# TRUE 
# 609896 
# [1] "bernstein_gsa"
# TRUE 
# 553161 
# [1] "farkkila_gsa"
# TRUE 
# 494340 
# [1] "franchimont_gsa"
# TRUE 
# 572376 
# [1] "franke_gsa"
# TRUE 
# 556806 
# [1] "helmsley_prism_gsa"
# TRUE 
# 243031 
# [1] "helmsley_xavier_prism_gsa"
# TRUE 
# 243035 
# [1] "hyams_protect_gsa"
# TRUE 
# 569169 
# [1] "lewis_sparc_gsa"
# TRUE 
# 607136 
# [1] "mccauley_new_gsa"
# TRUE 
# 578002 
# [1] "mcgovern_gsa"
# TRUE 
# 626738 
# [1] "moayyedi_imagine_gsa"
# TRUE 
# 588814 
# [1] "newberry_share_gsa"
# TRUE 
# 575540 
# [1] "niddk_cho_gsa"
# TRUE 
# 592894 
# [1] "niddk_duerr_gsa"
# TRUE 
# 580758 
# [1] "niddk_rioux_gsa"
# TRUE 
# 569777 
# [1] "niddk_silverberg_gsa"
# TRUE 
# 600691 
# [1] "palotie_hus_gsa"
# TRUE 
# 552919 
# [1] "pekow_share_gsa"
# TRUE 
# 572619 
# [1] "rioux_igenomed_gsa"
# TRUE 
# 542909 
# [1] "sands_msccr_gsa"
# TRUE 
# 594485 
# [1] "stampfer_gsa"
# TRUE 
# 571026 
# [1] "vermeire_gsa"
# TRUE 
# 584256 
# [1] "weersma_gsa"
# TRUE 
# 560257 
# [1] "xavier_prism_gsa"
# TRUE 
# 567986 
# [1] "xavier_share_gsa"
# TRUE 
# 567887 


##############

##############################################################
# 16.3 RENAME VARIANTS USING B38


path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_59_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_59_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned.bed \
--bim ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_edited.bim \
--fam ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt.fam \
--keep-allele-order --allow-no-sex \
--freq counts \
--make-bed --out ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_59_${i} | grep -E "completed"
done

for i in ${studies[@]}
do echo ${i} && head -2 ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated.bim;done

###################################################################################################################################################################################
