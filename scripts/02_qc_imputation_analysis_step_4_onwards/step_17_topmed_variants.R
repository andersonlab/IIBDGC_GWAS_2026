# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
##################################
# 17.- KEEP ONLY TOPMed VARIANTS #
##################################

MEM=200
path_gwas=/path/to/ibdgwas/IIBDGC/

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do 
bsub -J"TOPMed" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_TOPMed_subset_1_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_TOPMed_subset_1_${i} \
"awk 'NR==FNR{vals[\$2];next} (\$9) in vals' ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated.bim \
<(zcat ${path_gwas}resources/TOPMed/PASS.Variantsbravo-dbsnp-all_edited_3.vcf.gz) > ${path_gwas}pre_imputation/QC/${i}/${i}_TOPMed_variants"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_TOPMed_subset_1_${i} | grep -E "completed"
done

for i in ${studies[@]}
do 
bsub -J"TOPMed" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_TOPMed_subset_2_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_TOPMed_subset_2_${i} \
"cat ${path_gwas}pre_imputation/QC/${i}/${i}_TOPMed_variants | cut -f 9 \
> ${path_gwas}pre_imputation/QC/${i}/${i}_TOPMed_variants_ids"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_TOPMed_subset_2_${i} | grep -E "completed"
done

for i in ${studies[@]}
do echo ${i} && wc -l ${path_gwas}pre_imputation/QC/${i}/${i}_TOPMed_variants_ids
done


#######

for i in ${studies[@]}
do 
bsub -J"TOPMed" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_TOPMed_subset_3_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_TOPMed_subset_3_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated \
--keep-allele-order --allow-no-sex \
--extract ${path_gwas}pre_imputation/QC/${i}/${i}_TOPMed_variants_ids \
--make-bed --out ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_TOPMed_subset_3_${i} | grep -E "completed"
done


for i in ${studies[@]}
do echo ${i} && wc -l ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
done

########
# all_hce
# 402953 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/all_hce/all_hce_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# niddk_old_gwas
# 296804 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# australia_omniexome
# 727424 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/australia_omniexome/australia_omniexome_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# gwas1
# 440085 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/gwas1/gwas1_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# gwas2
# 763089 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/gwas2/gwas2_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# pittsburgh_gsa
# 874916 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# spain_gsa
# 567733 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/spain_gsa/spain_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# italy_gsa
# 545078 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/italy_gsa/italy_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# kiel_austria_sibdcs_gsa
# 529019 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# netherlands_gsa
# 566699 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/netherlands_gsa/netherlands_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# slovenia_gsa
# 508977 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/slovenia_gsa/slovenia_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# sweden_gsa
# 538700 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/sweden_gsa/sweden_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# niddk_broad_gsa
# 570723 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# niddk_feinstein_gsa
# 582173 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# basque_gsa
# 544451 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/basque_gsa/basque_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# lithuania_gsa
# 536257 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/lithuania_gsa/lithuania_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# belgium_louis_gsa
# 544126 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# belgium_franchimont_gsa
# 536055 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# belgium_vermeire_gsa
# 604979 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# prism_nfe_gsa
# 541087 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# prism_nfe_gwas
# 239864 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# finland_illugwas
# 233117 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/finland_illugwas/finland_illugwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# german_affy6_old_gwas
# 771834 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# norway_affy6_old_gwas
# 713535 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# belgium_inf1_old_gwas
# 298330 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# belgium_inf2_old_gwas
# 286925 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# cedars_370k_old_gwas
# 333825 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# cedars_610k_old_gwas
# 575223 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# cedars_omni_old_gwas
# 618094 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# swedish_uc_old_gwas
# 293310 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# mccauley_gsa
# 538558 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mccauley_gsa/mccauley_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# ccfa_gsa
# 579393 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/ccfa_gsa/ccfa_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# cedars_gsa
# 588802 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_gsa/cedars_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# bernstein_gsa
# 536303 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/bernstein_gsa/bernstein_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# farkkila_gsa
# 479808 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/farkkila_gsa/farkkila_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# franchimont_gsa
# 553569 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/franchimont_gsa/franchimont_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# franke_gsa
# 539704 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/franke_gsa/franke_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# helmsley_prism_gsa
# 240045 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# helmsley_xavier_prism_gsa
# 240086 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# hyams_protect_gsa
# 550882 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# lewis_sparc_gsa
# 585307 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# mccauley_new_gsa
# 559060 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# mcgovern_gsa
# 603387 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# moayyedi_imagine_gsa
# 570155 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# newberry_share_gsa
# 557449 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# niddk_cho_gsa
# 572514 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# niddk_duerr_gsa
# 560943 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# niddk_rioux_gsa
# 551384 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# niddk_silverberg_gsa
# 579570 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# palotie_hus_gsa
# 535460 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# pekow_share_gsa
# 554586 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# rioux_igenomed_gsa
# 526343 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# sands_msccr_gsa
# 575463 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# stampfer_gsa
# 551132 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/stampfer_gsa/stampfer_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# vermeire_gsa
# 563873 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/vermeire_gsa/vermeire_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# weersma_gsa
# 542803 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/weersma_gsa/weersma_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# xavier_prism_gsa
# 550115 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
# xavier_share_gsa
# 549624 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim
########


## /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("gwas1","gwas2","all_hce","niddk_old_gwas","australia_omniexome"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas"
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
  
  bim<-fread(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated.bim",sep=""),head=F)
  tm<-fread(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_TOPMed_variants",sep=""),head=F)

  print(dim(bim[which(!bim$V2 %in% tm$V9),]))
  # [1] 18780     6

  rm(list=ls()[!ls() %in% c("path","j","cohorts")]) 
}

############
# [1] "gwas1"
# [1] 5653    6
# [1] "gwas2"
# [1] 9437    6
# [1] "all_hce"
# [1] 18684     6
# [1] "niddk_old_gwas"
# [1] 3304    6
# [1] "australia_omniexome"
# [1] 10814     6
# [1] "pittsburgh_gsa"
# [1] 16072     6
# [1] "spain_gsa"
# [1] 7933    6
# [1] "italy_gsa"
# [1] 15359     6
# [1] "kiel_austria_sibdcs_gsa"
# [1] 15232     6
# [1] "netherlands_gsa"
# [1] 17903     6
# [1] "slovenia_gsa"
# [1] 14980     6
# [1] "sweden_gsa"
# [1] 16587     6
# [1] "niddk_broad_gsa"
# [1] 18096     6
# [1] "niddk_feinstein_gsa"
# [1] 18550     6
# [1] "basque_gsa"
# [1] 15504     6
# [1] "prism_nfe_gsa"
# [1] 17388     6
# [1] "lithuania_gsa"
# [1] 17133     6
# [1] "belgium_louis_gsa"
# [1] 18595     6
# [1] "belgium_franchimont_gsa"
# [1] 17190     6
# [1] "belgium_vermeire_gsa"
# [1] 31658     6
# [1] "prism_nfe_gwas"
# [1] 2984    6
# [1] "finland_illugwas"
# [1] 2869    6
# [1] "german_affy6_old_gwas"
# [1] 10375     6
# [1] "norway_affy6_old_gwas"
# [1] 9688    6
# [1] "belgium_inf1_old_gwas"
# [1] 3335    6
# [1] "belgium_inf2_old_gwas"
# [1] 3220    6
# [1] "cedars_370k_old_gwas"
# [1] 4807    6
# [1] "cedars_610k_old_gwas"
# [1] 8442    6
# [1] "cedars_omni_old_gwas"
# |--------------------------------------------------|
#   |==================================================|
#   [1] 100795      6
# [1] "swedish_uc_old_gwas"
# [1] 3296    6
# [1] "mccauley_gsa"
# |--------------------------------------------------|
#   |==================================================|
#   [1] 17630     6
# [1] "ccfa_gsa"
# [1] 20070     6
# [1] "cedars_gsa"
# [1] 21094     6
# [1] "bernstein_gsa"
# [1] 16858     6
# [1] "farkkila_gsa"
# [1] 14532     6
# [1] "franchimont_gsa"
# [1] 18807     6
# [1] "franke_gsa"
# [1] 17102     6
# [1] "helmsley_prism_gsa"
# [1] 2986    6
# [1] "helmsley_xavier_prism_gsa"
# [1] 2949    6
# [1] "hyams_protect_gsa"
# [1] 18287     6
# [1] "lewis_sparc_gsa"
# [1] 21829     6
# [1] "mccauley_new_gsa"
# [1] 18942     6
# [1] "mcgovern_gsa"
# [1] 23351     6
# [1] "moayyedi_imagine_gsa"
# [1] 18659     6
# [1] "newberry_share_gsa"
# [1] 18091     6
# [1] "niddk_cho_gsa"
# [1] 20380     6
# [1] "niddk_duerr_gsa"
# [1] 19815     6
# [1] "niddk_rioux_gsa"
# [1] 18393     6
# [1] "niddk_silverberg_gsa"
# [1] 21121     6
# [1] "palotie_hus_gsa"
# [1] 17459     6
# [1] "pekow_share_gsa"
# [1] 18033     6
# [1] "rioux_igenomed_gsa"
# [1] 16566     6
# [1] "sands_msccr_gsa"
# [1] 19022     6
# [1] "stampfer_gsa"
# [1] 19894     6
# [1] "vermeire_gsa"
# [1] 20383     6
# [1] "weersma_gsa"
# [1] 17454     6
# [1] "xavier_prism_gsa"
# [1] 17871     6
# [1] "xavier_share_gsa"
# [1] 18263     6

############

