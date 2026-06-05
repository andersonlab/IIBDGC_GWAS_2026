# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# UK-IIBGBD:

# See how variants were defined in: create_list_common_variants_among_cohorts.R

######################################################################
# 13.0- CREATE FILES INCLUDING DUPLICATED SAMPLES AND NON DUPLICATES #
######################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
    
path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_37_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_37_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
  --bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy \
  --allow-no-sex \
  --extract ${path_gwas}pre_imputation/QC/relatedness/list_common_variants_among_cohorts_iibdgc \
  --make-bed \
  --out ${path_gwas}pre_imputation/QC/relatedness/${i}_subset_withDuplicates"
done


for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_38_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_38_${i} \
" /path/to/software/plink_linux_x86_64_20181202/./plink \
  --bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2 \
  --allow-no-sex \
  --extract ${path_gwas}pre_imputation/QC/relatedness/list_common_variants_among_cohorts_iibdgc \
  --make-bed \
  --out ${path_gwas}pre_imputation/QC/relatedness/${i}_subset_noDuplicates"
done



for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_37_${i} | grep -E "completed"
done

for i in ${studies[@]}
do echo ${i} && wc -l ${path_gwas}pre_imputation/QC/relatedness/${i}_subset_withDuplicates.bim \
&& wc -l ${path_gwas}pre_imputation/QC/relatedness/${i}_subset_withDuplicates.fam
done


for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_37_${i} | grep -E "completed"
done

for i in ${studies[@]}
do echo ${i} && wc -l ${path_gwas}pre_imputation/QC/relatedness/${i}_subset_noDuplicates.bim \
&& wc -l ${path_gwas}pre_imputation/QC/relatedness/${i}_subset_noDuplicates.fam
done


wc -l ${path_gwas}pre_imputation/QC/relatedness/*_subset_withDuplicates.fam
##########
# 23890 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/all_hce_subset_withDuplicates.fam
# 1306 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/australia_omniexome_subset_withDuplicates.fam
# 1508 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/basque_gsa_subset_withDuplicates.fam
# 1501 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_franchimont_gsa_subset_withDuplicates.fam
# 1417 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf1_old_gwas_subset_withDuplicates.fam
# 272 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf2_old_gwas_subset_withDuplicates.fam
# 1518 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_louis_gsa_subset_withDuplicates.fam
# 3993 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_vermeire_gsa_subset_withDuplicates.fam
# 511 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/bernstein_gsa_subset_withDuplicates.fam
# 2177 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/ccfa_gsa_subset_withDuplicates.fam
# 605 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_370k_old_gwas_subset_withDuplicates.fam
# 889 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_610k_old_gwas_subset_withDuplicates.fam
# 3078 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_gsa_subset_withDuplicates.fam
# 1221 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_omni_old_gwas_subset_withDuplicates.fam
# 68 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/farkkila_gsa_subset_withDuplicates.fam
# 446 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/finland_illugwas_subset_withDuplicates.fam
# 2788 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/franchimont_gsa_subset_withDuplicates.fam
# 867 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/franke_gsa_subset_withDuplicates.fam
# 2823 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/german_affy6_old_gwas_subset_withDuplicates.fam
# 4680 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas1_subset_withDuplicates.fam
# 7778 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas2_subset_withDuplicates.fam
# 761 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/helmsley_prism_gsa_subset_withDuplicates.fam
# 1293 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/helmsley_xavier_prism_gsa_subset_withDuplicates.fam
# 418 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/hyams_protect_gsa_subset_withDuplicates.fam
# 953 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/italy_gsa_subset_withDuplicates.fam
# 14403 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/kiel_austria_sibdcs_gsa_subset_withDuplicates.fam
# 2859 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/lewis_sparc_gsa_subset_withDuplicates.fam
# 2231 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/lithuania_gsa_subset_withDuplicates.fam
# 782 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mccauley_gsa_subset_withDuplicates.fam
# 1619 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mccauley_new_gsa_subset_withDuplicates.fam
# 6022 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mcgovern_gsa_subset_withDuplicates.fam
# 1128 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/moayyedi_imagine_gsa_subset_withDuplicates.fam
# 4581 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/netherlands_gsa_subset_withDuplicates.fam
# 863 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/newberry_share_gsa_subset_withDuplicates.fam
# 5455 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_broad_gsa_subset_withDuplicates.fam
# 1764 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_cho_gsa_subset_withDuplicates.fam
# 1944 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_duerr_gsa_subset_withDuplicates.fam
# 8177 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_feinstein_gsa_subset_withDuplicates.fam
# 2760 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_old_gwas_subset_withDuplicates.fam
# 919 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_rioux_gsa_subset_withDuplicates.fam
# 2371 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_silverberg_gsa_subset_withDuplicates.fam
# 550 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/norway_affy6_old_gwas_subset_withDuplicates.fam
# 878 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/palotie_hus_gsa_subset_withDuplicates.fam
# 634 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/pekow_share_gsa_subset_withDuplicates.fam
# 2725 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/pittsburgh_gsa_subset_withDuplicates.fam
# 466 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/prism_nfe_gsa_subset_withDuplicates.fam
# 818 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/prism_nfe_gwas_subset_withDuplicates.fam
# 182 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/rioux_igenomed_gsa_subset_withDuplicates.fam
# 1430 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/sands_msccr_gsa_subset_withDuplicates.fam
# 264 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/slovenia_gsa_subset_withDuplicates.fam
# 3443 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/spain_gsa_subset_withDuplicates.fam
# 1468 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/stampfer_gsa_subset_withDuplicates.fam
# 1401 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/sweden_gsa_subset_withDuplicates.fam
# 1264 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/swedish_uc_old_gwas_subset_withDuplicates.fam
# 4697 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/vermeire_gsa_subset_withDuplicates.fam
# 709 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/weersma_gsa_subset_withDuplicates.fam
# 689 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/xavier_prism_gsa_subset_withDuplicates.fam
# 696 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/xavier_share_gsa_subset_withDuplicates.fam
# 146953 total
##########
# 146953 total

wc -l ${path_gwas}pre_imputation/QC/relatedness/*_subset_noDuplicates.fam
##########
# 22426 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/all_hce_subset_noDuplicates.fam
# 1298 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/australia_omniexome_subset_noDuplicates.fam
# 1491 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/basque_gsa_subset_noDuplicates.fam
# 2 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_franchimont_gsa_subset_noDuplicates.fam
# 1067 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf1_old_gwas_subset_noDuplicates.fam
# 161 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf2_old_gwas_subset_noDuplicates.fam
# 1502 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_louis_gsa_subset_noDuplicates.fam
# 119 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_vermeire_gsa_subset_noDuplicates.fam
# 464 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/bernstein_gsa_subset_noDuplicates.fam
# 1 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/ccfa_gsa_subset_noDuplicates.fam
# 453 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_370k_old_gwas_subset_noDuplicates.fam
# 467 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_610k_old_gwas_subset_noDuplicates.fam
# 17 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_gsa_subset_noDuplicates.fam
# 673 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_omni_old_gwas_subset_noDuplicates.fam
# 68 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/farkkila_gsa_subset_noDuplicates.fam
# 425 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/finland_illugwas_subset_noDuplicates.fam
# 2722 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/franchimont_gsa_subset_noDuplicates.fam
# 857 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/franke_gsa_subset_noDuplicates.fam
# 2631 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/german_affy6_old_gwas_subset_noDuplicates.fam
# 4676 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas1_subset_noDuplicates.fam
# 5132 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas2_subset_noDuplicates.fam
# 936 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/helmsley_xavier_prism_gsa_subset_noDuplicates.fam
# 418 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/hyams_protect_gsa_subset_noDuplicates.fam
# 949 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/italy_gsa_subset_noDuplicates.fam
# 13777 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/kiel_austria_sibdcs_gsa_subset_noDuplicates.fam
# 2780 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/lewis_sparc_gsa_subset_noDuplicates.fam
# 2210 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/lithuania_gsa_subset_noDuplicates.fam
# 5 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mccauley_gsa_subset_noDuplicates.fam
# 1602 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mccauley_new_gsa_subset_noDuplicates.fam
# 5898 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mcgovern_gsa_subset_noDuplicates.fam
# 1107 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/moayyedi_imagine_gsa_subset_noDuplicates.fam
# 4554 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/netherlands_gsa_subset_noDuplicates.fam
# 860 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/newberry_share_gsa_subset_noDuplicates.fam
# 185 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_broad_gsa_subset_noDuplicates.fam
# 1448 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_cho_gsa_subset_noDuplicates.fam
# 1320 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_duerr_gsa_subset_noDuplicates.fam
# 5911 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_feinstein_gsa_subset_noDuplicates.fam
# 1346 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_old_gwas_subset_noDuplicates.fam
# 640 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_rioux_gsa_subset_noDuplicates.fam
# 1952 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_silverberg_gsa_subset_noDuplicates.fam
# 550 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/norway_affy6_old_gwas_subset_noDuplicates.fam
# 877 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/palotie_hus_gsa_subset_noDuplicates.fam
# 600 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/pekow_share_gsa_subset_noDuplicates.fam
# 2030 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/pittsburgh_gsa_subset_noDuplicates.fam
# 42 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/prism_nfe_gsa_subset_noDuplicates.fam
# 569 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/prism_nfe_gwas_subset_noDuplicates.fam
# 167 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/rioux_igenomed_gsa_subset_noDuplicates.fam
# 1080 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/sands_msccr_gsa_subset_noDuplicates.fam
# 262 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/slovenia_gsa_subset_noDuplicates.fam
# 3443 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/spain_gsa_subset_noDuplicates.fam
# 1467 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/stampfer_gsa_subset_noDuplicates.fam
# 1379 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/sweden_gsa_subset_noDuplicates.fam
# 1248 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/swedish_uc_old_gwas_subset_noDuplicates.fam
# 4678 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/vermeire_gsa_subset_noDuplicates.fam
# 13 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/weersma_gsa_subset_noDuplicates.fam
# 663 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/xavier_prism_gsa_subset_noDuplicates.fam
# 651 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/xavier_share_gsa_subset_noDuplicates.fam
# 114269 total
##########
# 114269 total


#######################
# 13.1 MERGE ALL FILES

#### /software/R-4.3.1/bin/R

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



#####

cohorts<-cohorts[-1]

dat<-matrix(ncol=3,nrow=length(cohorts))
dat<-as.data.frame(dat)

for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_withDuplicates.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_withDuplicates.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_withDuplicates.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/relatedness/list_cohorts_tomerge_withDuplicates.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#####

# helmsley_prism_gsa has no sample after excluding dups:
cohorts<-cohorts[which(!cohorts %in% "helmsley_prism_gsa")]

dat<-matrix(ncol=3,nrow=length(cohorts))
dat<-as.data.frame(dat)

for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_noDuplicates.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_noDuplicates.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_noDuplicates.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/relatedness/list_cohorts_tomerge_noDuplicates.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


#############################


path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/relatedness/logs/stderr_plink_39 \
-o ${path_gwas}pre_imputation/QC/relatedness/logs/stdout_plink_39 \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/relatedness/australia_omniexome_subset_withDuplicates \
--allow-no-sex \
--merge-list ${path_gwas}pre_imputation/QC/relatedness/list_cohorts_tomerge_withDuplicates.txt \
--make-bed --out ${path_gwas}pre_imputation/QC/relatedness/iibdgc_merged_withDuplicates"

bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/relatedness/logs/stderr_plink_40 \
-o ${path_gwas}pre_imputation/QC/relatedness/logs/stdout_plink_40 \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/relatedness/australia_omniexome_subset_noDuplicates \
--allow-no-sex \
--merge-list ${path_gwas}pre_imputation/QC/relatedness/list_cohorts_tomerge_noDuplicates.txt \
--make-bed --out ${path_gwas}pre_imputation/QC/relatedness/iibdgc_merged_noDuplicates"

wc -l ${path_gwas}pre_imputation/QC/relatedness/iibdgc_merged_withDuplicates.fam
# 146953 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/iibdgc_merged_withDuplicates.fam
wc -l ${path_gwas}pre_imputation/QC/relatedness/iibdgc_merged_noDuplicates.fam
# 114269 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/iibdgc_merged_noDuplicates.fam


# # remove intermediate files - at the very end
# for i in ${studies[@]}
# do rm ${path_gwas}pre_imputation/QC/relatedness/{i}_subset_noDuplicates.* \
# && rm ${path_gwas}pre_imputation/QC/relatedness/{i}_subset_withDuplicates.*
# done
#   
  
###########################################
# 13.1- REDEFINE LIST OF COMMON VARIANTS: #
###########################################

path_gwas=/path/to/ibdgwas/IIBDGC/
  
# 1.1 EXCLUDE VARIANTS IN HIGH LD REGIONS:
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/relatedness/iibdgc_merged_noDuplicates \
--allow-no-sex \
--exclude /nfs/team152/carl/protocols/QC-paper/high-LD-regions.txt --range \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD
# 3150 variants and 114269 people pass filters and QC.
# Among remaining phenotypes, 71892 are cases and 42377 are controls.


# 1.3 EXCLUDE ANY ASSOCIATED VARIANT - exclude from analysis samples with sex==0 - not confirmed that the sample is the right one

########## /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(viridis)

path<-"/path/to/ibdgwas/IIBDGC/"

a<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates.fam",sep=""),head=F)
table(a$V5)
#    0     1     2 
# 1297 54775 58197 

write.table(a[which(a$V5==0),c(1,2)],paste(path,"pre_imputation/QC/pca_1000gp/list_iibdgc_merged_noDuplicates_sex_0",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

#######


/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD \
--remove ${path_gwas}pre_imputation/QC/pca_1000gp/list_iibdgc_merged_noDuplicates_sex_0 \
--assoc --out ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD
# 114269 phenotype values loaded from .fam.
# --remove: 112972 people remaining.

#### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(viridis)

path<-"/path/to/ibdgwas/IIBDGC/"

as<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD.assoc",sep=""),head=T)

summary(as$P)
# 
nrow(as[which(as$P<=1E-4),])
# [1] 927

write.table(as[which(as$P<=1E-4),"SNP",drop=F],paste(path,"pre_imputation/QC/pca_1000gp/list_associated_variants_toexclude",sep=""),col.names=F
            ,row.names=F,quote=F,sep="\t")

nrow(as[which(as$P>1E-4),])
# [1] 2223
write.table(as[which(as$P>1E-4),"SNP",drop=F],paste(path,"pre_imputation/QC/pca_1000gp/list_noassociated_variants_tokeep",sep=""),col.names=F
            ,row.names=F,quote=F,sep="\t")

####

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD \
--exclude ${path_gwas}pre_imputation/QC/pca_1000gp/list_associated_variants_toexclude \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD_noassoclogP4
# 2223 variants and 114269 people pass filters and QC.
# Among remaining phenotypes, 71892 are cases and 42377 are controls.



# 1.2 PRUNE VARIANTS:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD_noassoclogP4 \
--allow-no-sex \
--indep-pairwise 50 5 0.2 \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD_noassoclogP4

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD_noassoclogP4 \
--extract ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD_noassoclogP4.prune.in \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD_noassoclogP4_pruned
# 2173 variants and 114269 people pass filters and QC.
# Among remaining phenotypes, 71892 are cases and 42377 are controls.


## remove intermed files:
rm ${path_gwas}pre_imputation/QC/relatedness/iibdgc_merged_noDuplicates.*
rm ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD.*
rm ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD_pruned.*
  
  
###########################################################################################################################################################################

######################################
# 13.2.- SUBSET LIST OF COMMON VARIANTS #
######################################

#################################################
# 13.2.1- FROM IIBDGC FILES WITH DUP SAMPLES:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/relatedness/iibdgc_merged_withDuplicates \
--extract ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD_noassoclogP4.prune.in \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned
# 2173 variants and 146953 people pass filters and QC.
# Among remaining phenotypes, 95320 are cases and 51633 are controls.

# remove intermed files:
rm ${path_gwas}pre_imputation/QC/relatedness/iibdgc_merged_withDuplicates.*
  
#################################################
# 13.2.1- FROM 1000GP FILES AND MERGE 

for chr in {1..22}; do /path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}pre_imputation/QC/1000gp/1000GP_chr${chr}_b37.bed \
--bim ${path_gwas}pre_imputation/QC/1000gp/1000GP_chr${chr}_b37_edited.bim \
--fam ${path_gwas}pre_imputation/QC/1000gp/1000GP_chr${chr}_b37.fam \
--extract ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD_noassoclogP4.prune.in \
--make-bed --out ${path_gwas}pre_imputation/QC/1000gp/1000GP_chr${chr}_b37_subset_iibdgc_noHighLD_noassoclogP4_pruned;done

#### R

# no chrx or chry
dat<-matrix(ncol=3,nrow=21)
dat<-as.data.frame(dat)
for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/1000gp/1000GP_chr",i+1,"_b37_subset_iibdgc_noHighLD_noassoclogP4_pruned.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/1000gp/1000GP_chr",i+1,"_b37_subset_iibdgc_noHighLD_noassoclogP4_pruned.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/1000gp/1000GP_chr",i+1,"_b37_subset_iibdgc_noHighLD_noassoclogP4_pruned.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/1000gp/list_1000GP_files_subset_iibdgc_noHighLD_noassoclogP4_pruned_tomerge.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/1000gp/1000GP_chr1_b37_subset_iibdgc_noHighLD_noassoclogP4_pruned \
--merge-list ${path_gwas}pre_imputation/QC/1000gp/list_1000GP_files_subset_iibdgc_noHighLD_noassoclogP4_pruned_tomerge.txt \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc_noHighLD_noassoclogP4_pruned
# 2173 variants and 2504 people pass filters and QC.
# Note: No phenotypes present.

# remove intermediate files:
rm ${path_gwas}pre_imputation/QC/1000gp/1000GP_chr*_b37_subset*
  
###################################################################################################################################################
# # EXPLORE PROJECTIONS - pre plink2 instructions
# 
# # 
# # Yes, it's likely that something like PLINK 1.9's --pca-clusters/--pca-cluster-names will eventually make it into 2.0.  
# # With that said, PCA projection is actually already supported, the workflow is just a bit more convoluted for now.
# # 
# # Step 1: Export allele frequencies and PCA variant weights from your reference dataset.  E.g.
# # plink2 --bfile hapmap --freq --pca var-wts --out pca_hapmap
# # 
# # Step 2: Use --score to compute the necessary dot products with the variant weights.  E.g.
# # plink2 --bfile mydata --read-freq pca_hapmap.afreq --score pca_hapmap.eigenvec.var 2 3 header-read no-mean-imputation variance-normalize 
# # --score-col-nums 5-14 --out pca_proj_mydata
# # 
# # These PCs will be scaled a bit differently from pca_hapmap.eigenvec (you need to multiply or divide the PCs by sqrt(eigenvalue) to put 
# them on the same scale).  One way to make the PCs directly comparable is to also run step 2 on the original dataset.
# # 

# # EXPLORE PROJECTIONS - plink2 instructions
# 
# PCA projection with --score
# Since --score's new 'variance-standardize' modifier applies the same transformation to G as --pca does, --score can now 
# execute the vector-matrix multiply corresponding to PCA projection.
# 
# The following command exports PCs to project onto, along with the allele frequencies needed to calibrate the 'variance-standardize' operation:
# 
# plink2 --pfile ref_data \
#        --freq counts \
#        --pca allele-wts \
#        --out ref_pcs
# 
# You can then project onto those PCs with
# 
# plink2 --pfile new_data \
#        --read-freq ref_pcs.acount \
#        --score ref_pcs.eigenvec.allele 2 5 header-read no-mean-imputation \
#                variance-standardize \
#        --score-col-nums 6-15 \
#        --out new_projection
# 
# Note that these PCs will be scaled a bit differently from ref_data.eigenvec; you need to multiply or divide the PCs by a multiple of sqrt(eigenvalue)
# to put them on the same scale.
###################################################################################################################################################


################################
# 13.3 - RUN PCA AND PROJECT DATA #
################################

#########################################################
## 13.3.1.- ESTIMATE PC, EXPORT PCS, AND ALLELE FREQUENCIES:

# The following command exports PCs to project onto, along with the allele frequencies needed to calibrate the 'variance-standardize' operation:

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=1000 # next time only this

bsub -J"bcftools" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 10 \
-e ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stderr_pca_1000gp \
-o ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_pca_1000gp \
"/path/to/software/./plink2  \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc_noHighLD_noassoclogP4_pruned \
--freq counts \
--pca allele-wts \
--out ${path_gwas}pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc_noHighLD_noassoclogP4_pruned_ref_pcs \
--threads 8 --allow-no-sex --memory $MEM"
# Job <2549848> is submitted to queue <normal>.

less ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_pca_1000gp
# 2504 samples (0 females, 0 males, 2504 ambiguous; 2504 founders) loaded from
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc_noHigh
# LD_noassoclogP4_pruned.fam.
# 2173 variants loaded from
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc_noHigh
# LD_noassoclogP4_pruned.bim.
# Note: No phenotype data present.
# Extracting eigenvalues and eigenvectors... done.
# --pca: Allele weights written to
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc_noHigh
# LD_noassoclogP4_pruned_ref_pcs.eigenvec.allele
# .
# --pca: Eigenvectors written to
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc_noHigh
# LD_noassoclogP4_pruned_ref_pcs.eigenvec
# , and eigenvalues written to
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc_noHigh
# LD_noassoclogP4_pruned_ref_pcs.eigenval

#########################################################
## 13.3.2 MERGE 1000GP AND IIBDGC (WITH DUPLICATED SAMPLES)

# plink1.9, 2.0 has not bmerge available yet

path=/path/to/ibdgwas/IIBDGC/
MEM=800 

bsub -J"bcftools" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stderr_pca_1000gp_merge_iibdgc \
-o ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_pca_1000gp_merge_iibdgc \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc_noHighLD_noassoclogP4_pruned \
--bmerge ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned.bed \
${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned.bim \
${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned.fam \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/1000GP_all_b37_plus_iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned"
# Job <2550329> is submitted to queue <normal>.

less ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_pca_1000gp_merge_iibdgc
# 2173 variants and 149457 people pass filters and QC.
# Among remaining phenotypes, 95319 are cases and 51633 are controls.  (2504
#                                                                       phenotypes are missing.)

146953+2504
# [1] 149457 OK


#########################################################
## 13.3.3.- PROJECT ONTO THOSE PCS 1000GP AND IIBDGC

path=/path/to/ibdgwas/IIBDGC/
MEM=1000 # next time only this

bsub -J"bcftools" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 10 \
-e ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stderr_pca_1000gp_iibdgc_pca \
-o ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_pca_1000gp_iibdgc_pca \
"/path/to/software/./plink2 \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/1000GP_all_b37_plus_iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned \
--read-freq ${path_gwas}pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc_noHighLD_noassoclogP4_pruned_ref_pcs.acount \
--score ${path_gwas}pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc_noHighLD_noassoclogP4_pruned_ref_pcs.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
--score-col-nums 6-15 \
--out ${path_gwas}pre_imputation/QC/pca_1000gp/1000GP_all_b37_plus_iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned_new_projection"
# Job <993401> is submitted to queue <normal>.

less ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_pca_1000gp_iibdgc_pca
# 149457 samples (74713 females, 70609 males, 4135 ambiguous; 149457 founders)
# 2173 variants loaded from
# --score: Results written to
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc_noHighL
# D_noassoclogP4_pruned_new_projection.sscore


########## /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(viridis)

path<-"/path/to/ibdgwas/IIBDGC/"

# no german_illu, or chop
cohorts<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
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

length(cohorts)
# [1] 58


for (i in 1:length(cohorts)){
  a<-read.table(paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_withDuplicates.fam",sep=""),head=F)
  a$cohort<-cohorts[i]
  if(i==1){
    fam<-a
  }else{
    fam<-rbind(a,fam)
  }
}
dim(fam)
# [1] 146953      7

table(fam$cohort)

samples1000<-read.table("/path/to/project",head=T)
table(samples1000$pop,samples1000$super_pop)
#      AFR AMR EAS EUR SAS
# ACB  96   0   0   0   0
# ASW  61   0   0   0   0
# BEB   0   0   0   0  86
# CDX   0   0  93   0   0
# CEU   0   0   0  99   0
# CHB   0   0 103   0   0
# CHS   0   0 105   0   0
# CLM   0  94   0   0   0
# ESN  99   0   0   0   0
# FIN   0   0   0  99   0
# GBR   0   0   0  91   0
# GIH   0   0   0   0 103
# GWD 113   0   0   0   0
# IBS   0   0   0 107   0
# ITU   0   0   0   0 102
# JPT   0   0 104   0   0
# KHV   0   0  99   0   0
# LWK  99   0   0   0   0
# MSL  85   0   0   0   0
# MXL   0  64   0   0   0
# PEL   0  85   0   0   0
# PJL   0   0   0   0  96
# PUR   0 104   0   0   0
# STU   0   0   0   0 102
# TSI   0   0   0 107   0
# YRI 108   0   0   0   0

fam<-fam[,c(1,7)]
colnames(fam)[1]<-c("sample")
fam$super_pop<-"IIBDGC"
fam$pop<-NA

samples1000$cohort<-"1000GP"
samples1000<-samples1000[,c("sample","cohort","super_pop","pop")]
fam<-rbind(fam,samples1000)
dim(fam)
# [1] 149457      4


pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/1000GP_all_b37_plus_iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned_new_projection.sscore",sep=""),head=F)
colnames(pca)<-c("FID","IID","PHENO1","ALLELE_CT","NAMED_ALLELE_DOSAGE_SUM","PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG"
                 ,"PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG")
dim(pca)
# [1] 149457     12


pca<-merge(pca,fam[,c("sample","cohort","super_pop","pop")],by.x="IID",by.y="sample",all.x=T)
dim(pca)
# [1] 149457     18

table(pca$super_pop)
# AFR    AMR    EAS    EUR IIBDGC    SAS 
# 661    347    504    503 146953    489 

# TOTAL VARIANCE EXPLAINED BY EACH PC:

eigenval<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc_noHighLD_noassoclogP4_pruned_ref_pcs.eigenval",sep=""),head=F)

eigenval$var_exp<-NA
for (i in 1:nrow(eigenval)){
  eigenval$var_exp[i]<-(eigenval$V1[i] / sum(eigenval$V1))*100
}

eigenval
#           V1   var_exp
# 1  213.76500 55.708561
# 2   88.39460 23.036212
# 3   29.43400  7.670693
# 4   23.35740  6.087092
# 5    6.00457  1.564830
# 6    5.35610  1.395835
# 7    4.73664  1.234399
# 8    4.44954  1.159579
# 9    4.14242  1.079542
# 10   4.07993  1.063257

# PC1 EUR - non EUR
# PC2 AFR  - non AFR
# 



pca$study<-NA

pca$study[which(pca$cohort=="australia_omniexome")]<-"Australia"
pca$study[which(pca$cohort %in% c("gwas1","gwas2","all_hce"))]<-"UK"
pca$study[which(pca$cohort %in% c("spain_gsa","basque_gsa"))]<-"Spain"
pca$study[which(pca$cohort %in% c("kiel_austria_sibdcs_gsa","german_affy6_old_gwas","franke_gsa"))]<-"Germany"                 
pca$study[which(pca$cohort %in% c("italy_gsa"))]<-"Italy"   
pca$study[which(pca$cohort %in% c("netherlands_gsa"))]<-"Netherlands"                     
pca$study[which(pca$cohort %in% c("slovenia_gsa"))]<-"Slovenia"                     
pca$study[which(pca$cohort %in% c("sweden_gsa","swedish_uc_old_gwas"))]<-"Sweden"                    
pca$study[which(pca$cohort %in% c("finland_illugwas","farkkila_gsa"))]<-"Finland"
pca$study[which(pca$cohort %in% c("chop_old_gwas"))]<-"CHOP"
pca$study[which(pca$cohort %in% c("norway_affy6_old_gwas"))]<-"Norway"
pca$study[which(pca$cohort %in% c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","belgium_louis_gsa","belgium_franchimont_gsa"
                                  ,"belgium_vermeire_gsa","franchimont_gsa"))]<-"Belgium"
pca$study[which(pca$cohort %in% c("lithuania_gsa"))]<-"Lithuania"
pca$study[which(pca$cohort %in% c("prism_nfe_gwas","prism_nfe_gsa","cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas",
                                  "cedars_gsa",
                                  "niddk_broad_gsa","niddk_feinstein_gsa","niddk_old_gwas","pittsburgh_gsa","mccauley_gsa","ccfa_gsa",
                                  "bernstein_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
                                  "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
                                  "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
                                  "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
                                  "xavier_share_gsa"))]<-"USA"


table(pca$study,useNA="ifany")
# Australia     Belgium     Finland     Germany       Italy   Lithuania 
# 1306       11489         514       18093         953        2231 
# Netherlands      Norway    Slovenia       Spain      Sweden          UK 
# 4581         550         264        4951        2665       36348 
# USA        <NA> 
# 63008        2504 


####################
# 4.- PLOT RESULTS #
####################

######################################################
# 4.1 PLOT 1000GP ALONE, COLOUR BY SUPER POPULATION:

# creat plots PC1 to PC6, 1000GP only
pca$pop<-as.factor(pca$pop)
pca$pop<-factor(pca$pop, levels=c("CEU","FIN","GBR","IBS","TSI",
                                  "CDX","CHB","CHS","JPT","KHV",
                                  "BEB","GIH","ITU","PJL","STU",
                                  "ACB","ASW","ESN","GWD","LWK","MSL","YRI",
                                  "CLM","MXL","PEL","PUR"))

# EUR = "CEU","FIN","GBR","IBS","TSI"
# EAS = "CDX","CHB","CHS","JPT","KHV"
# SAS =  "BEB","GIH","ITU","PJL","STU"
# AFR = "ACB","ASW","ESN","GWD","LWK","MSL","YRI"
# AMF = "CLM","MXL","PEL","PUR"

pca$super_pop<-as.factor(pca$super_pop)
pca$super_pop<-factor(pca$super_pop, levels=c("EUR","AFR","EAS","SAS","AMR","IIBDGC"))

colnames(pca)<-gsub("_AVG","",colnames(pca))

pna<-qplot()+theme(
  panel.background = element_rect(fill = "transparent") # bg of the panel
  , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
  , panel.grid.major = element_blank() # get rid of major grid
  , panel.grid.minor = element_blank() # get rid of minor grid
  , legend.background = element_rect(fill = "transparent") # get rid of legend bg
  , axis.title=element_blank()
  , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
)

p11<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p12<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p13<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p14<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p15<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p16<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")

p21<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p22<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p23<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p24<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p25<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p26<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")

p31<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p32<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p33<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p34<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p35<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p36<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")

p41<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p42<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p43<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p44<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p45<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p46<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")

p51<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p52<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p53<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p54<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p55<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p56<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")

p61<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p62<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p63<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p64<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p65<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p66<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")


r1<-ggarrange(pna,p12,p13,p14,p15,p16,ncol=6,legend=c("none"))
r2<-ggarrange(p21,pna,p23,p24,p25,p26,ncol=6,legend=c("none"))
r3<-ggarrange(p31,p32,pna,p34,p35,p36,ncol=6,legend=c("none"))
r4<-ggarrange(p41,p42,p43,pna,p45,p36,ncol=6,legend=c("none"))
r5<-ggarrange(p51,p52,p53,p54,pna,p56,ncol=6,legend=c("none"))
r6<-ggarrange(p61,p62,p63,p64,p65,pna,ncol=6,common.legend = TRUE,legend=c("bottom"))
dev.off()


pdf(paste(path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_3kvariants.pdf",sep=""),width=30,height=30)
ggarrange(r1,r2,r3,r4,r5,r6,nrow=6,common.legend = TRUE,legend=c("bottom"),heights = c(1,1,1,1,1,1.2))
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_3kvariants.pdf ~/tmp_plots/",sep=""))



#######################################################################
# 13.4.1 PLOT 1000GP ALONE, AND SEPARATELY , COLOUR BY SUPER POPULATION:


# PC1 EUR - non EUR
# PC2 AFR  - non AFR
# PC3 SAS - non SAS


pca$study<-as.factor(pca$study)
pca$study<-factor(pca$study, levels=c("Spain","Italy","Slovenia","Belgium","Germany","UK","Netherlands","Lithuania","Norway","Sweden","Finland",
                                      "CHOP","USA","Australia"))

# PC1 - PC2

ylims<-c(min(pca$PC1),max(pca$PC1))
xlims<-c(min(pca$PC2),max(pca$PC2))

p0<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC2,PC1)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims)

p1<-ggplot(pca[which(pca$cohort!="1000GP"),],aes(PC2,PC1)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims)


pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2.pdf",sep=""),width=14,height=6)
ggarrange(p0,p1,legend=c("right"),widths = c(1,1.1))
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2.pdf ~/tmp_plots/",sep=""))



# PC3 - PC4

ylims<-c(min(pca$PC3),max(pca$PC3))
xlims<-c(min(pca$PC4),max(pca$PC4))

p0<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC4,PC3)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims)

p1<-ggplot(pca[which(pca$cohort!="1000GP"),],aes(PC4,PC3)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims)


pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC3_PC4.pdf",sep=""),width=14,height=6)
ggarrange(p0,p1,legend=c("right"),widths = c(1,1))
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC3_PC4.pdf ~/tmp_plots/",sep=""))


############################################################################################################################################
# INFER POPULATIONS:

pop<-c("AMR","AFR","EAS","SAS","EUR")

pca$inferred_population<-NA

for (i in 1:length(pop)) {
  
  pca$inferred_population[which(pca$PC1>=min(pca$PC1[which(pca$super_pop==pop[i])]) & pca$PC1<=max(pca$PC1[which(pca$super_pop==pop[i])]) &
                                  pca$PC2>=min(pca$PC2[which(pca$super_pop==pop[i])]) & pca$PC2<=max(pca$PC2[which(pca$super_pop==pop[i])]) &
                                  pca$PC3>=min(pca$PC3[which(pca$super_pop==pop[i])]) & pca$PC3<=max(pca$PC3[which(pca$super_pop==pop[i])]) &
                                  pca$PC4>=min(pca$PC4[which(pca$super_pop==pop[i])]) & pca$PC4<=max(pca$PC4[which(pca$super_pop==pop[i])]) &
                                  pca$PC5>=min(pca$PC5[which(pca$super_pop==pop[i])]) & pca$PC5<=max(pca$PC5[which(pca$super_pop==pop[i])]) &
                                  pca$PC6>=min(pca$PC6[which(pca$super_pop==pop[i])]) & pca$PC6<=max(pca$PC6[which(pca$super_pop==pop[i])]) &
                                  pca$PC7>=min(pca$PC7[which(pca$super_pop==pop[i])]) & pca$PC7<=max(pca$PC7[which(pca$super_pop==pop[i])]) &
                                  pca$PC8>=min(pca$PC8[which(pca$super_pop==pop[i])]) & pca$PC8<=max(pca$PC8[which(pca$super_pop==pop[i])]) &
                                  pca$PC9>=min(pca$PC9[which(pca$super_pop==pop[i])]) & pca$PC9<=max(pca$PC9[which(pca$super_pop==pop[i])]) &
                                  pca$PC10>=min(pca$PC10[which(pca$super_pop==pop[i])]) & pca$PC10<=max(pca$PC10[which(pca$super_pop==pop[i])]) )]<-pop[i]
  
}


table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           AFR    AMR    EAS    EUR    SAS   <NA>
# EUR         0      0      0    503      0      0
# AFR       661      0      0      0      0      0
# EAS         0      0    504      0      0      0
# SAS         0      0      0      0    489      0
# AMR        66    250      0     31      0      0
# IIBDGC   2545   2753    989 137523   1312   1831


# PC1 - PC2

ylims<-c(min(pca$PC1),max(pca$PC1))
xlims<-c(min(pca$PC2),max(pca$PC2))

p0<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC2,PC1)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims)

p1<-ggplot(pca[which(is.na(pca$inferred_population)),],aes(PC2,PC1)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims)

# PC1 - PC3

ylims<-c(min(pca$PC1),max(pca$PC1))
xlims<-c(min(pca$PC3),max(pca$PC3))

p2<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC3,PC1)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims)

p3<-ggplot(pca[which(is.na(pca$inferred_population)),],aes(PC3,PC1)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims)


c1<-ggarrange(p0,p2,nrow=2,common.legend = TRUE,legend=c("right"))
c2<-ggarrange(p1,p3,nrow=2,common.legend = TRUE,legend=c("right"))
dev.off()

pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_PC3_noInferredEth.pdf",sep=""),width=14,height=11)
ggarrange(c1,c2,widths = c(1,1),ncol=2)
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_PC3_noInferredEth.pdf ~/tmp_plots/",sep=""))


####################

# EXPAND LIMITS FOR EUR SUBSET OF SAMPLES

pc1max<-max(pca$PC1[which(pca$super_pop=="EUR")])+((max(pca$PC1[which(pca$super_pop=="EUR")])-min(pca$PC1[which(pca$super_pop=="EUR")]))*0.3)
pc1min<-min(pca$PC1[which(pca$super_pop=="EUR")])-((max(pca$PC1[which(pca$super_pop=="EUR")])-min(pca$PC1[which(pca$super_pop=="EUR")]))*0.3)

pc2max<-max(pca$PC2[which(pca$super_pop=="EUR")])+((max(pca$PC2[which(pca$super_pop=="EUR")])-min(pca$PC2[which(pca$super_pop=="EUR")]))*0.3)
pc2min<-min(pca$PC2[which(pca$super_pop=="EUR")])-((max(pca$PC2[which(pca$super_pop=="EUR")])-min(pca$PC2[which(pca$super_pop=="EUR")]))*0.3)

pc3max<-max(pca$PC3[which(pca$super_pop=="EUR")])+((max(pca$PC3[which(pca$super_pop=="EUR")])-min(pca$PC3[which(pca$super_pop=="EUR")]))*0.3)
pc3min<-min(pca$PC3[which(pca$super_pop=="EUR")])-((max(pca$PC3[which(pca$super_pop=="EUR")])-min(pca$PC3[which(pca$super_pop=="EUR")]))*0.3)


print(paste("PC1 expanded limits",pc1max,pc1min))
# [1] "PC1 expanded limits 0.12924582 0.01132678"
print(paste("PC2 expanded limits",pc2max,pc2min))
# "PC2 expanded limits -0.05502758 -0.17508102"
print(paste("PC3 expanded limits",pc3max,pc3min))
# [1] "PC3 expanded limits 0.03398294 -0.06531114"


# PC1 - PC2 - with new limits


ylims<-c(min(pca$PC1),max(pca$PC1))
xlims<-c(min(pca$PC2),max(pca$PC2))

p0<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC2,PC1)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims) + 
  geom_hline(yintercept=pc1max, linetype="dashed", color = "grey") + 
  geom_hline(yintercept=pc1min, linetype="dashed", color = "grey") + 
  geom_vline(xintercept=pc2max, linetype="dashed", color = "grey") + 
  geom_vline(xintercept=pc2min, linetype="dashed", color = "grey")


p1<-ggplot(pca[which(is.na(pca$inferred_population)),],aes(PC2,PC1)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims) + 
  geom_hline(yintercept=pc1max, linetype="dashed", color = "grey") + 
  geom_hline(yintercept=pc1min, linetype="dashed", color = "grey") + 
  geom_vline(xintercept=pc2max, linetype="dashed", color = "grey") + 
  geom_vline(xintercept=pc2min, linetype="dashed", color = "grey")


# PC1 - PC3

ylims<-c(min(pca$PC1),max(pca$PC1))
xlims<-c(min(pca$PC3),max(pca$PC3))

p2<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC3,PC1)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims) + 
  geom_hline(yintercept=pc1max, linetype="dashed", color = "grey") + 
  geom_hline(yintercept=pc1min, linetype="dashed", color = "grey") + 
  geom_vline(xintercept=pc3max, linetype="dashed", color = "grey") + 
  geom_vline(xintercept=pc3min, linetype="dashed", color = "grey")

p3<-ggplot(pca[which(is.na(pca$inferred_population)),],aes(PC3,PC1)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims) + 
  geom_hline(yintercept=pc1max, linetype="dashed", color = "grey") + 
  geom_hline(yintercept=pc1min, linetype="dashed", color = "grey") + 
  geom_vline(xintercept=pc3max, linetype="dashed", color = "grey") + 
  geom_vline(xintercept=pc3min, linetype="dashed", color = "grey")


c1<-ggarrange(p0,p2,nrow=2,common.legend = TRUE,legend=c("right"))
c2<-ggarrange(p1,p3,nrow=2,common.legend = TRUE,legend=c("right"))
dev.off()

pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_PC3_noInferredEth_newLim.pdf",sep=""),width=14,height=11)
ggarrange(c1,c2,widths = c(1,1),ncol=2)
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_PC3_noInferredEth_newLim.pdf ~/tmp_plots/",sep=""))


#### reclasify

pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max & pca$PC1>=pc1min & pca$PC2<=pc2max & pca$PC2>=pc2min
                              & pca$PC3<=pc3max & pca$PC3>=pc3min)]<-"EUR"



table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           AFR    AMR    EAS    EUR    SAS   <NA>
# EUR         0      0      0    503      0      0
# AFR       661      0      0      0      0      0
# EAS         0      0    504      0      0      0
# SAS         0      0      0      0    489      0
# AMR        66    250      0     31      0      0
# IIBDGC   2545   2753    989 138901   1312    453


# 
# # PC1 - PC2 - with new limits
# 
# 
# ylims<-c(min(pca$PC1),max(pca$PC1))
# xlims<-c(min(pca$PC2),max(pca$PC2))
# 
# p0<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC2,PC1)) +
#   geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") +
#   xlim(xlims) + ylim(ylims) +
#   geom_hline(yintercept=pc1max, linetype="dashed", color = "grey") +
#   geom_hline(yintercept=pc1min, linetype="dashed", color = "grey") +
#   geom_vline(xintercept=pc2max, linetype="dashed", color = "grey") +
#   geom_vline(xintercept=pc2min, linetype="dashed", color = "grey")
# 
# 
# p1<-ggplot(pca[which(is.na(pca$inferred_population)),],aes(PC2,PC1)) +
#   geom_point(aes(color = study)) +
#   xlim(xlims) + ylim(ylims) +
#   geom_hline(yintercept=pc1max, linetype="dashed", color = "grey") +
#   geom_hline(yintercept=pc1min, linetype="dashed", color = "grey") +
#   geom_vline(xintercept=pc2max, linetype="dashed", color = "grey") +
#   geom_vline(xintercept=pc2min, linetype="dashed", color = "grey")
# 
# 
# # PC1 - PC3
# 
# ylims<-c(min(pca$PC1),max(pca$PC1))
# xlims<-c(min(pca$PC3),max(pca$PC3))
# 
# p2<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC3,PC1)) +
#   geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") +
#   xlim(xlims) + ylim(ylims) +
#   geom_hline(yintercept=pc1max, linetype="dashed", color = "grey") +
#   geom_hline(yintercept=pc1min, linetype="dashed", color = "grey") +
#   geom_vline(xintercept=pc3max, linetype="dashed", color = "grey") +
#   geom_vline(xintercept=pc3min, linetype="dashed", color = "grey")
# 
# p3<-ggplot(pca[which(is.na(pca$inferred_population)),],aes(PC3,PC1)) +
#   geom_point(aes(color = study)) +
#   xlim(xlims) + ylim(ylims) +
#   geom_hline(yintercept=pc1max, linetype="dashed", color = "grey") +
#   geom_hline(yintercept=pc1min, linetype="dashed", color = "grey") +
#   geom_vline(xintercept=pc3max, linetype="dashed", color = "grey") +
#   geom_vline(xintercept=pc3min, linetype="dashed", color = "grey")
# 
# 
# c1<-ggarrange(p0,p2,nrow=2,common.legend = TRUE,legend=c("right"))
# c2<-ggarrange(p1,p3,nrow=2,common.legend = TRUE,legend=c("right"))
# dev.off()
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_PC3_noInferredEth_newLim_2.pdf",sep=""),width=14,height=11)
# ggarrange(c1,c2,widths = c(1,1),ncol=2)
# dev.off()
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_PC3_noInferredEth_newLim_2.pdf ~/tmp_plots/",sep=""))
# 


####################

# EXPAND LIMITS FOR SAS SUBSET OF SAMPLES

pc1max_sas<-max(pca$PC1[which(pca$super_pop=="SAS")])+((max(pca$PC1[which(pca$super_pop=="SAS")])-min(pca$PC1[which(pca$super_pop=="SAS")]))*0.3)
pc1min_sas<-min(pca$PC1[which(pca$super_pop=="SAS")])-((max(pca$PC1[which(pca$super_pop=="SAS")])-min(pca$PC1[which(pca$super_pop=="SAS")]))*0.3)

pc2max_sas<-max(pca$PC2[which(pca$super_pop=="SAS")])+((max(pca$PC2[which(pca$super_pop=="SAS")])-min(pca$PC2[which(pca$super_pop=="SAS")]))*0.3)
pc2min_sas<-min(pca$PC2[which(pca$super_pop=="SAS")])-((max(pca$PC2[which(pca$super_pop=="SAS")])-min(pca$PC2[which(pca$super_pop=="SAS")]))*0.3)

pc3max_sas<-max(pca$PC3[which(pca$super_pop=="SAS")])+((max(pca$PC3[which(pca$super_pop=="SAS")])-min(pca$PC3[which(pca$super_pop=="SAS")]))*0.3)
pc3min_sas<-min(pca$PC3[which(pca$super_pop=="SAS")])-((max(pca$PC3[which(pca$super_pop=="SAS")])-min(pca$PC3[which(pca$super_pop=="SAS")]))*0.3)


print(paste("PC1 expanded limits",pc1max_sas,pc1min_sas))
# [1] "PC1 expanded limits 0.13417616 0.02434064"
print(paste("PC2 expanded limits",pc2max_sas,pc2min_sas))
# [1] "PC2 expanded limits 0.05007192 -0.12805192"
print(paste("PC3 expanded limits",pc3max_sas,pc3min_sas))
# [1] "PC3 expanded limits 0.14906881 0.01080849"

pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max_sas & pca$PC1>=pc1min_sas & pca$PC2<=pc2max_sas & pca$PC2>=pc2min_sas
                              & pca$PC3<=pc3max_sas & pca$PC3>=pc3min_sas)]<-"SAS"

table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           AFR    AMR    EAS    EUR    SAS   <NA>
# EUR         0      0      0    503      0      0
# AFR       661      0      0      0      0      0
# EAS         0      0    504      0      0      0
# SAS         0      0      0      0    489      0
# AMR        66    250      0     31      0      0
# IIBDGC   2545   2753    989 138901   1482    283

####################

# EXPAND LIMITS FOR EAS SUBSET OF SAMPLES

pc1max_eas<-max(pca$PC1[which(pca$super_pop=="EAS")])+((max(pca$PC1[which(pca$super_pop=="EAS")])-min(pca$PC1[which(pca$super_pop=="EAS")]))*0.3)
pc1min_eas<-min(pca$PC1[which(pca$super_pop=="EAS")])-((max(pca$PC1[which(pca$super_pop=="EAS")])-min(pca$PC1[which(pca$super_pop=="EAS")]))*0.3)

pc2max_eas<-max(pca$PC2[which(pca$super_pop=="EAS")])+((max(pca$PC2[which(pca$super_pop=="EAS")])-min(pca$PC2[which(pca$super_pop=="EAS")]))*0.3)
pc2min_eas<-min(pca$PC2[which(pca$super_pop=="EAS")])-((max(pca$PC2[which(pca$super_pop=="EAS")])-min(pca$PC2[which(pca$super_pop=="EAS")]))*0.3)

pc3max_eas<-max(pca$PC3[which(pca$super_pop=="EAS")])+((max(pca$PC3[which(pca$super_pop=="EAS")])-min(pca$PC3[which(pca$super_pop=="EAS")]))*0.3)
pc3min_eas<-min(pca$PC3[which(pca$super_pop=="EAS")])-((max(pca$PC3[which(pca$super_pop=="EAS")])-min(pca$PC3[which(pca$super_pop=="EAS")]))*0.3)


print(paste("PC1 expanded limits",pc1max_eas,pc1min_eas))
# [1] "PC1 expanded limits 0.16424364 0.08424556"
print(paste("PC2 expanded limits",pc2max_eas,pc2min_eas))
#[1] "PC2 expanded limits 0.206812 0.104236"
print(paste("PC3 expanded limits",pc3max_eas,pc3min_eas))
# [1] "PC3 expanded limits 0.06467009 -0.05431119"

pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max_eas & pca$PC1>=pc1min_eas & pca$PC2<=pc2max_eas & pca$PC2>=pc2min_eas
                              & pca$PC3<=pc3max_eas & pca$PC3>=pc3min_eas)]<-"EAS"

table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           AFR    AMR    EAS    EUR    SAS   <NA>
# EUR         0      0      0    503      0      0
# AFR       661      0      0      0      0      0
# EAS         0      0    504      0      0      0
# SAS         0      0      0      0    489      0
# AMR        66    250      0     31      0      0
# IIBDGC   2545   2753   1030 138901   1482    242

####################

####################

# EXPAND LIMITS FOR AFR SUBSET OF SAMPLES

pc1max_afr<-max(pca$PC1[which(pca$super_pop=="AFR")])+((max(pca$PC1[which(pca$super_pop=="AFR")])-min(pca$PC1[which(pca$super_pop=="AFR")]))*0.3)
pc1min_afr<-min(pca$PC1[which(pca$super_pop=="AFR")])-((max(pca$PC1[which(pca$super_pop=="AFR")])-min(pca$PC1[which(pca$super_pop=="AFR")]))*0.3)

pc2max_afr<-max(pca$PC2[which(pca$super_pop=="AFR")])+((max(pca$PC2[which(pca$super_pop=="AFR")])-min(pca$PC2[which(pca$super_pop=="AFR")]))*0.3)
pc2min_afr<-min(pca$PC2[which(pca$super_pop=="AFR")])-((max(pca$PC2[which(pca$super_pop=="AFR")])-min(pca$PC2[which(pca$super_pop=="AFR")]))*0.3)

pc3max_afr<-max(pca$PC3[which(pca$super_pop=="AFR")])+((max(pca$PC3[which(pca$super_pop=="AFR")])-min(pca$PC3[which(pca$super_pop=="AFR")]))*0.3)
pc3min_afr<-min(pca$PC3[which(pca$super_pop=="AFR")])-((max(pca$PC3[which(pca$super_pop=="AFR")])-min(pca$PC3[which(pca$super_pop=="AFR")]))*0.3)


print(paste("PC1 expanded limits",pc1max_afr,pc1min_afr))
# [1] "PC1 expanded limits 0.2187623 -0.4029193"
print(paste("PC2 expanded limits",pc2max_afr,pc2min_afr))
# [1] "PC2 expanded limits 0.08841407 -0.07793137"
print(paste("PC3 expanded limits",pc3max_afr,pc3min_afr))
# [1] "PC3 expanded limits 0.09617305 -0.18353655"

pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max_afr & pca$PC1>=pc1min_afr & pca$PC2<=pc2max_afr & pca$PC2>=pc2min_afr
                              & pca$PC3<=pc3max_afr & pca$PC3>=pc3min_afr)]<-"AFR"

table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           AFR    AMR    EAS    EUR    SAS   <NA>
# EUR         0      0      0    503      0      0
# AFR       661      0      0      0      0      0
# EAS         0      0    504      0      0      0
# SAS         0      0      0      0    489      0
# AMR        66    250      0     31      0      0
# IIBDGC   2699   2753   1030 138901   1482     88

####################

# EXPAND LIMITS FOR AMR SUBSET OF SAMPLES

pc1max_amr<-max(pca$PC1[which(pca$super_pop=="AMR")])+((max(pca$PC1[which(pca$super_pop=="AMR")])-min(pca$PC1[which(pca$super_pop=="AMR")]))*0.3)
pc1min_amr<-min(pca$PC1[which(pca$super_pop=="AMR")])-((max(pca$PC1[which(pca$super_pop=="AMR")])-min(pca$PC1[which(pca$super_pop=="AMR")]))*0.3)

pc2max_amr<-max(pca$PC2[which(pca$super_pop=="AMR")])+((max(pca$PC2[which(pca$super_pop=="AMR")])-min(pca$PC2[which(pca$super_pop=="AMR")]))*0.3)
pc2min_amr<-min(pca$PC2[which(pca$super_pop=="AMR")])-((max(pca$PC2[which(pca$super_pop=="AMR")])-min(pca$PC2[which(pca$super_pop=="AMR")]))*0.3)

pc3max_amr<-max(pca$PC3[which(pca$super_pop=="AMR")])+((max(pca$PC3[which(pca$super_pop=="AMR")])-min(pca$PC3[which(pca$super_pop=="AMR")]))*0.3)
pc3min_amr<-min(pca$PC3[which(pca$super_pop=="AMR")])-((max(pca$PC3[which(pca$super_pop=="AMR")])-min(pca$PC3[which(pca$super_pop=="AMR")]))*0.3)


print(paste("PC1 expanded limits",pc1max_amr,pc1min_amr))
# [1] "PC1 expanded limits 0.2300459 -0.2356709"
print(paste("PC2 expanded limits",pc2max_amr,pc2min_amr))
# [1] "PC2 expanded limits 0.14405716 -0.19010796"
print(paste("PC3 expanded limits",pc3max_amr,pc3min_amr))
# [1] "PC3 expanded limits 0.064479033 -0.295113623"

pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max_amr & pca$PC1>=pc1min_amr & pca$PC2<=pc2max_amr & pca$PC2>=pc2min_amr
                              & pca$PC3<=pc3max_amr & pca$PC3>=pc3min_amr)]<-"AMR"

table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           AFR    AMR    EAS    EUR    SAS
# EUR         0      0      0    503      0
# AFR       661      0      0      0      0
# EAS         0      0    504      0      0
# SAS         0      0      0      0    489
# AMR        66    250      0     31      0
# IIBDGC   2699   2841   1030 138901   1482


pca$inferred_population<-as.factor(pca$inferred_population)
pca$inferred_population<-factor(pca$inferred_population, levels=c("EUR","AFR","EAS","SAS","AMR"))

###########################################
# some inferred AFR in AMR, reclassify:

table(pca[which(pca$PC3< -0.05),"super_pop"])
# EUR    AFR    EAS    SAS    AMR IIBDGC 
#   0      3      0      0    251   1232  
table(pca[which(pca$PC3< -0.05),"inferred_population"])
# EUR AFR EAS SAS AMR 
#  17 602   0   0 867 

pca$inferred_population[which((pca$inferred_population=="AFR") & (pca$PC3< -0.05) )]<-"AMR"
table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           EUR    AFR    EAS    SAS    AMR
# EUR       503      0      0      0      0
# AFR         0    658      0      0      3
# EAS         0      0    504      0      0
# SAS         0      0      0    489      0
# AMR        31      7      0      0    309
# IIBDGC 138901   2159   1030   1482   3381


table(pca[which(pca$PC1< -0.05),"super_pop"])
# EUR    AFR    EAS    SAS    AMR IIBDGC 
#        0    658      0      0      4   2013 
table(pca[which(pca$PC1> -0.05),"super_pop"])
# EUR    AFR    EAS    SAS    AMR IIBDGC 
# 503      3    504    489    343 144940 

table(pca[which(pca$PC1< -0.05),"inferred_population"])
# EUR  AFR  EAS  SAS  AMR 
#  0 2459    0    0  216 

pca$inferred_population[which((pca$inferred_population=="AMR") & (pca$PC1< -0.05) )]<-"AFR"
table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           EUR    AFR    EAS    SAS    AMR
# EUR       503      0      0      0      0
# AFR         0    658      0      0      3
# EAS         0      0    504      0      0
# SAS         0      0      0    489      0
# AMR        31      8      0      0    308
# IIBDGC 138901   2374   1030   1482   3166


pca$inferred_population[which((pca$inferred_population=="AFR") & (pca$PC1> -0.05) )]<-"AMR"
table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           EUR    AFR    EAS    SAS    AMR
# EUR       503      0      0      0      0
# AFR         0    658      0      0      3
# EAS         0      0    504      0      0
# SAS         0      0      0    489      0
# AMR        31      4      0      0    312
# IIBDGC 138901   2013   1030   1482   3527

# PLOT FINAL RESULTS

# PC1 - PC2 - with new limits

ylims<-c(min(pca$PC1),max(pca$PC1))
xlims<-c(min(pca$PC2),max(pca$PC2))

p0<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC2,PC1)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims) 


p1<-ggplot(pca[which(pca$cohort!="1000GP"),],aes(PC2,PC1)) +
  geom_point(aes(color = inferred_population)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims) 


# PC1 - PC3

ylims<-c(min(pca$PC1),max(pca$PC1))
xlims<-c(min(pca$PC3),max(pca$PC3))

p2<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC3,PC1)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims) 

p3<-ggplot(pca[which(pca$cohort!="1000GP"),],aes(PC3,PC1)) +
  geom_point(aes(color = inferred_population)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims) 


c1<-ggarrange(p0,p2,nrow=2,common.legend = TRUE,legend=c("right"))
c1<-annotate_figure(c1,top = text_grob("1000GP"))
c2<-ggarrange(p1,p3,nrow=2,common.legend = TRUE,legend=c("right"))
c2<-annotate_figure(c2,top = text_grob("IIBDGC"))
dev.off()

pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_PC3_inferred_final.pdf",sep=""),width=14,height=11)
ggarrange(c1,c2,widths = c(1,1.1),ncol=2)
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_PC3_inferred_final.pdf ~/tmp_plots/",sep=""))


#########

x1<-as.data.frame.matrix(table(pca$cohort[which(pca$PHENO1=="2")],pca$inferred_population[which(pca$PHENO1=="2")],useNA="ifany"))
colnames(x1)<-paste(colnames(x1),"n_cases",sep="_")
x1$study<-rownames(x1)
x2<-as.data.frame.matrix(table(pca$cohort[which(pca$PHENO1=="1")],pca$inferred_population[which(pca$PHENO1=="1")],useNA="ifany"))
colnames(x2)<-paste(colnames(x2),"n_ctr",sep="_")
x2$study<-rownames(x2)



x<-merge(x1,x2,by="study",all=T)
write.table(x,paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_ancestries_per_study_case_control_withDuplicates.tsv",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")

write.table(pca,paste(path,"pre_imputation/QC/pca_1000gp/1000GP_all_b37_plus_iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned_new_projection_edited_withDuplicates.sscore",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")
               

#########################################################
# 5.- EXTRACT EUR SAMPLES AND IDENTIFY ASHKENAZI SUBSET #
#########################################################

pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/1000GP_all_b37_plus_iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned_new_projection_edited_withDuplicates.sscore",sep=""),head=T)

table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           AFR    AMR    EAS    EUR    SAS
# AFR       658      3      0      0      0
# AMR         4    312      0     31      0
# EAS         0      0    504      0      0
# EUR         0      0      0    503      0
# IIBDGC   2013   3527   1030 138901   1482
# SAS         0      0      0      0    489



write.table(pca[which(pca$inferred_population=="EUR"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eur_ancestry_samples_withDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca[which(pca$inferred_population=="AFR"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_afr_ancestry_samples_withDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca[which(pca$inferred_population=="EAS"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eas_ancestry_samples_withDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca[which(pca$inferred_population=="SAS"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_sas_ancestry_samples_withDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca[which(pca$inferred_population=="AMR"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_amr_ancestry_samples_withDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

#### create lists with no duplicates:

fam<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_noDuplicates.fam",sep=""),head=F)

pca_nodup<-pca[which(pca$IID %in% fam$V2),]
table(pca_nodup$super_pop,pca_nodup$inferred_population,useNA="ifany")
#           AFR    AMR    EAS    EUR    SAS
# AFR         0      0      0      0      0
# AMR         0      0      0      0      0
# EAS         0      0      0      0      0
# EUR         0      0      0      0      0
# IIBDGC   1573   2752    773 108003   1168
# SAS         0      0      0      0      0

write.table(pca_nodup[which(pca_nodup$inferred_population=="EUR"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eur_ancestry_samples_noDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca_nodup[which(pca_nodup$inferred_population=="AFR"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_afr_ancestry_samples_noDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca_nodup[which(pca_nodup$inferred_population=="EAS"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eas_ancestry_samples_noDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca_nodup[which(pca_nodup$inferred_population=="SAS"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_sas_ancestry_samples_noDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca_nodup[which(pca_nodup$inferred_population=="AMR"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_amr_ancestry_samples_noDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")


x1<-as.data.frame.matrix(table(pca_nodup$cohort[which(pca_nodup$PHENO1=="2")],pca_nodup$inferred_population[which(pca_nodup$PHENO1=="2")],useNA="ifany"))
colnames(x1)<-paste(colnames(x1),"n_cases",sep="_")
x1$study<-rownames(x1)
x2<-as.data.frame.matrix(table(pca_nodup$cohort[which(pca_nodup$PHENO1=="1")],pca_nodup$inferred_population[which(pca_nodup$PHENO1=="1")],useNA="ifany"))
colnames(x2)<-paste(colnames(x2),"n_ctr",sep="_")
x2$study<-rownames(x2)

x<-merge(x1,x2,by="study",all=T)
write.table(x,paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_ancestries_per_study_case_control_noDuplicates.tsv",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")

  
###########################################################################################################################################

# RUN PCA IN EUROEPAN ANCESTRY SAMPLES, ESTIMATE PCS IN NON DUPLICATED EUR SUBSET, AND PROJECT ALL:


#########################################################
## 13.5.1 - ESTIMATE PC, EXPORT PCS, AND ALLELE FREQUENCIES:

# The following command exports PCs to project onto, along with the allele frequencies needed to calibrate the 'variance-standardize' operation:

wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD_noassoclogP4_pruned.fam
# 114269 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD_noassoclogP4_pruned.fam
wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned.fam
# 146953 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned.fam

wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned.bim
# 2173
wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned.bim
# 2173

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD_noassoclogP4_pruned \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_ancestry_samples_noDuplicates \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD_noassoclogP4_pruned_EURonly_noDuplicates
# 114269 phenotype values loaded from .fam.
# 2173 variants and 108003 people pass filters and QC.
# Among remaining phenotypes, 67230 are cases and 40773 are controls.

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_ancestry_samples_withDuplicates \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned_EURonly_withDuplicates
# 146953 phenotype values loaded from .fam.
# 2173 variants and 138901 people pass filters and QC.
# Among remaining phenotypes, 89204 are cases and 49697 are controls.


path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=5000 # next time only this

bsub -J"bcftools" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 10 \
-e ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stderr_pca_eur_nodup \
-o ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_pca_eur_nodup \
"/path/to/software/./plink2  \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_noDuplicates_noHighLD_noassoclogP4_pruned_EURonly_noDuplicates \
--freq counts \
--pca approx allele-wts \
--out ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_EUR_nodup_ref_pcs \
--threads 8 --allow-no-sex --memory $MEM"
# Job <2583325> is submitted to queue <normal>.


#########################################################
## 1.2 - PROJECT ONTO THOSE PCS ALL IIBDGC

path=/path/to/ibdgwas/IIBDGC/
MEM=1000 # next time only this

bsub -J"bcftools" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 10 \
-e ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_pca_eur_nodup_project_withDuplicates \
-o ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_pca_eur_nodup_project_withDuplicates \
"/path/to/software/./plink2 \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned_EURonly_withDuplicates \
--read-freq ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_EUR_nodup_ref_pcs.acount \
--score ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_EUR_nodup_ref_pcs.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
--score-col-nums 6-15 \
--out ${path_gwas}pre_imputation/QC/pca_1000gp/iibdgc_EUR_withDuplicates_ref_pcs_new_projection"
# Job <995336> is submitted to queue <normal>.



###############################

### /software/R-4.3.1/bin/R


library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(viridis)

path<-"/path/to/ibdgwas/IIBDGC/"

# no german_illu or chop
cohorts<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
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

length(cohorts)
# [1] 58


for (i in 1:length(cohorts)){
  a<-read.table(paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_withDuplicates.fam",sep=""),head=F)
  a$cohort<-cohorts[i]
  if(i==1){
    fam<-a
  }else{
    fam<-rbind(a,fam)
  }
}

dim(fam)
# [1] 138901      7 # OK


pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_EUR_withDuplicates_ref_pcs_new_projection.sscore",sep=""),head=F)
colnames(pca)<-c("FID","IID","PHENO1","ALLELE_CT","NAMED_ALLELE_DOSAGE_SUM","PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG"
                 ,"PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG")
dim(pca)
# [1]  138901     15 # OK
colnames(pca)<-gsub("_AVG","",colnames(pca))

pca<-merge(pca,fam[,c("V1","cohort")],by.x="FID",by.y="V1",all.x=T)
dim(pca)
# [1] 138901 16

table(pca$cohort)

############
# all_hce       australia_omniexome                basque_gsa 
# 22705                      1252                      1497 
# belgium_franchimont_gsa     belgium_inf1_old_gwas     belgium_inf2_old_gwas 
# 1468                      1407                       271 
# belgium_louis_gsa      belgium_vermeire_gsa             bernstein_gsa 
# 1506                      3912                       484 
# ccfa_gsa      cedars_370k_old_gwas      cedars_610k_old_gwas 
# 1940                       513                       830 
# cedars_gsa      cedars_omni_old_gwas              farkkila_gsa 
# 2714                      1114                        68 
# finland_illugwas           franchimont_gsa                franke_gsa 
# 442                      2376                       856 
# german_affy6_old_gwas                     gwas1                     gwas2 
# 2813                      4673                      7763 
# helmsley_prism_gsa helmsley_xavier_prism_gsa         hyams_protect_gsa 
# 676                      1183                       325 
# italy_gsa   kiel_austria_sibdcs_gsa           lewis_sparc_gsa 
# 948                     14160                      2542 
# lithuania_gsa              mccauley_gsa          mccauley_new_gsa 
# 2219                       764                      1329 
# mcgovern_gsa      moayyedi_imagine_gsa           netherlands_gsa 
# 5190                      1018                      4358 
# newberry_share_gsa           niddk_broad_gsa             niddk_cho_gsa 
# 752                      5158                      1620 
# niddk_duerr_gsa       niddk_feinstein_gsa            niddk_old_gwas 
# 1875                      7013                      2740 
# niddk_rioux_gsa      niddk_silverberg_gsa     norway_affy6_old_gwas 
# 877                      2134                       548 
# palotie_hus_gsa           pekow_share_gsa            pittsburgh_gsa 
# 862                       540                      2717 
# prism_nfe_gsa            prism_nfe_gwas        rioux_igenomed_gsa 
# 430                       724                       172 
# sands_msccr_gsa              slovenia_gsa                 spain_gsa 
# 1129                       264                      3431 
# stampfer_gsa                sweden_gsa       swedish_uc_old_gwas 
# 1443                      1369                      1261 
# vermeire_gsa               weersma_gsa          xavier_prism_gsa 
# 4577                       686                       622 
# xavier_share_gsa 
# 641 
############

# TOTAL VARIANCE EXPLAINED BY EACH PC:


eigenval<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_EUR_nodup_ref_pcs.eigenval",sep=""),head=F)

eigenval$var_exp<-NA
for (i in 1:nrow(eigenval)){
  eigenval$var_exp[i]<-(eigenval$V1[i] / sum(eigenval$V1))*100
}

eigenval
#          V1   var_exp
# 1  219.7500 21.121475
# 2  128.9840 12.397417
# 3  101.1800  9.725009
# 4   88.4031  8.496946
# 5   86.4672  8.310875
# 6   85.8062  8.247342
# 7   84.1642  8.089520
# 8   83.1815  7.995067
# 9   81.7221  7.854795
# 10  80.7520  7.761553


### check loadings for PC1 - later used for infer ashkenazi ancestry:


score<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_EUR_nodup_ref_pcs.eigenvec.allele",sep=""),head=F)
colnames(score)<-c("chr","snp","alt","ref","allele",paste("PC",seq(1:10),sep=""))
score$index<-rownames(score)

ggplot(score, aes(chr, PC1)) + geom_point(aes(colour = chr))+ scale_color_viridis(discrete=F)

table(score$chr)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
# 360 398 258 268 312 310 196 286 216 228 230 140 164 134  94 146  94 154  82 108 
# 21  22 
# 70  98 

score$chr_tmp<-gsub(":.*","",score$snp)
table(score$chr,score$chr_tmp)
# all OK

score$pos<-gsub("[0-9]{1,2}:","",score$snp)
score$pos<-gsub("_.*","",score$pos)
score$pos<-as.numeric(score$pos)
score<-score[order(score$chr,score$pos),]
score$index<-seq(1:nrow(score))


for (i in 1){
  print(i)
  tmp<-score[,c("index","chr",paste("PC",i,sep=""))]
  colnames(tmp)<-c("index","chr","PC")
  pdf(paste("/path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/PC",i,"_loadings_iibdgc_EUR_nodup.pdf",sep=""),width=10,height=7)
  # print(ggplot(tmp, aes(index, PC))  + stat_binhex()+ scale_color_viridis(discrete=F))
  print(ggplot(tmp, aes(index, PC)) + geom_point(aes(colour = chr)) + scale_color_viridis(discrete=F))
  dev.off()
  rm(tmp)
}

# not driven by any specific locus


pca$study<-NA

pca$study[which(pca$cohort=="australia_omniexome")]<-"Australia"
pca$study[which(pca$cohort %in% c("gwas1","gwas2","all_hce"))]<-"UK"
pca$study[which(pca$cohort %in% c("spain_gsa","basque_gsa"))]<-"Spain"
pca$study[which(pca$cohort %in% c("kiel_austria_sibdcs_gsa","german_affy6_old_gwas","franke_gsa"))]<-"Germany"                 
pca$study[which(pca$cohort %in% c("italy_gsa"))]<-"Italy"   
pca$study[which(pca$cohort %in% c("netherlands_gsa"))]<-"Netherlands"                     
pca$study[which(pca$cohort %in% c("slovenia_gsa"))]<-"Slovenia"                     
pca$study[which(pca$cohort %in% c("sweden_gsa","swedish_uc_old_gwas"))]<-"Sweden"                    
pca$study[which(pca$cohort %in% c("finland_illugwas","farkkila_gsa"))]<-"Finland"
pca$study[which(pca$cohort %in% c("norway_affy6_old_gwas"))]<-"Norway"
pca$study[which(pca$cohort %in% c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","belgium_louis_gsa","belgium_franchimont_gsa"
                                  ,"belgium_vermeire_gsa","franchimont_gsa"))]<-"Belgium"
pca$study[which(pca$cohort %in% c("lithuania_gsa"))]<-"Lithuania"
pca$study[which(pca$cohort %in% c("prism_nfe_gwas","prism_nfe_gsa","cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas",
                                  "cedars_gsa",
                                  "niddk_broad_gsa","niddk_feinstein_gsa","niddk_old_gwas","pittsburgh_gsa","mccauley_gsa","ccfa_gsa",
                                  "bernstein_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
                                  "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
                                  "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
                                  "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
                                  "xavier_share_gsa"))]<-"USA"

table(pca$study,useNA="ifany")
# Australia     Belgium     Finland     Germany       Italy   Lithuania 
# 1252       10940         510       17829         948        2219 
# Netherlands      Norway    Slovenia       Spain      Sweden          UK 
# 4358         548         264        4928        2630       35141 
# USA 
# 57334 


pca$study<-as.factor(pca$study)
pca$study<-factor(pca$study, levels=c("Spain","Italy","Slovenia","Belgium","Germany","Netherlands","UK","Lithuania",
                                      "Sweden","Norway","Finland","Australia","USA"))


table(pca$study,useNA="ifany")
# Spain       Italy    Slovenia     Belgium     Germany Netherlands 
# 4928         948         264       10940       17829        4358 
# UK   Lithuania      Sweden      Norway     Finland   Australia 
# 35141        2219        2630         548         510        1252 
# USA 
# 57334 


########################
# 6.- PLOT EUR RESULTS #
########################

######################################################
# 4.1 PLOT 1000GP ALONE, COLOUR BY COUNTRY:

# create plots PC1 to PC10

# pna<-qplot()+theme(
#   panel.background = element_rect(fill = "transparent") # bg of the panel
#   , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
#   , panel.grid.major = element_blank() # get rid of major grid
#   , panel.grid.minor = element_blank() # get rid of minor grid
#   , legend.background = element_rect(fill = "transparent") # get rid of legend bg
#   , axis.title=element_blank()
#   , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
# )
# 
# 
# p11<-qplot(PC1,PC1, data = pca, colour = study)
# p12<-qplot(PC2,PC1, data = pca, colour = study)
# p13<-qplot(PC3,PC1, data = pca, colour = study)
# p14<-qplot(PC4,PC1, data = pca, colour = study)
# p15<-qplot(PC5,PC1, data = pca, colour = study)
# p16<-qplot(PC6,PC1, data = pca, colour = study)
# p17<-qplot(PC7,PC1, data = pca, colour = study)
# p18<-qplot(PC8,PC1, data = pca, colour = study)
# p19<-qplot(PC9,PC1, data = pca, colour = study)
# p110<-qplot(PC10,PC1, data = pca, colour = study)
# 
# p21<-qplot(PC1,PC2, data = pca, colour = study)
# p22<-qplot(PC2,PC2, data = pca, colour = study)
# p23<-qplot(PC3,PC2, data = pca, colour = study)
# p24<-qplot(PC4,PC2, data = pca, colour = study)
# p25<-qplot(PC5,PC2, data = pca, colour = study)
# p26<-qplot(PC6,PC2, data = pca, colour = study)
# p27<-qplot(PC7,PC2, data = pca, colour = study)
# p28<-qplot(PC8,PC2, data = pca, colour = study)
# p29<-qplot(PC9,PC2, data = pca, colour = study)
# p210<-qplot(PC10,PC2, data = pca, colour = study)
# 
# p31<-qplot(PC1,PC3, data = pca, colour = study)
# p32<-qplot(PC2,PC3, data = pca, colour = study)
# p33<-qplot(PC3,PC3, data = pca, colour = study)
# p34<-qplot(PC4,PC3, data = pca, colour = study)
# p35<-qplot(PC5,PC3, data = pca, colour = study)
# p36<-qplot(PC6,PC3, data = pca, colour = study)
# p37<-qplot(PC7,PC3, data = pca, colour = study)
# p38<-qplot(PC8,PC3, data = pca, colour = study)
# p39<-qplot(PC9,PC3, data = pca, colour = study)
# p310<-qplot(PC10,PC3, data = pca, colour = study)
# 
# p41<-qplot(PC1,PC4, data = pca, colour = study)
# p42<-qplot(PC2,PC4, data = pca, colour = study)
# p43<-qplot(PC3,PC4, data = pca, colour = study)
# p44<-qplot(PC4,PC4, data = pca, colour = study)
# p45<-qplot(PC5,PC4, data = pca, colour = study)
# p46<-qplot(PC6,PC4, data = pca, colour = study)
# p47<-qplot(PC7,PC4, data = pca, colour = study)
# p48<-qplot(PC8,PC4, data = pca, colour = study)
# p49<-qplot(PC9,PC4, data = pca, colour = study)
# p410<-qplot(PC10,PC4, data = pca, colour = study)
# 
# p51<-qplot(PC1,PC5, data = pca, colour = study)
# p52<-qplot(PC2,PC5, data = pca, colour = study)
# p53<-qplot(PC3,PC5, data = pca, colour = study)
# p54<-qplot(PC4,PC5, data = pca, colour = study)
# p55<-qplot(PC5,PC5, data = pca, colour = study)
# p56<-qplot(PC6,PC5, data = pca, colour = study)
# p57<-qplot(PC7,PC5, data = pca, colour = study)
# p58<-qplot(PC8,PC5, data = pca, colour = study)
# p59<-qplot(PC9,PC5, data = pca, colour = study)
# p510<-qplot(PC10,PC5, data = pca, colour = study)
# 
# p61<-qplot(PC1,PC6, data = pca, colour = study)
# p62<-qplot(PC2,PC6, data = pca, colour = study)
# p63<-qplot(PC3,PC6, data = pca, colour = study)
# p64<-qplot(PC4,PC6, data = pca, colour = study)
# p65<-qplot(PC5,PC6, data = pca, colour = study)
# p66<-qplot(PC6,PC6, data = pca, colour = study)
# p67<-qplot(PC7,PC6, data = pca, colour = study)
# p68<-qplot(PC8,PC6, data = pca, colour = study)
# p69<-qplot(PC9,PC6, data = pca, colour = study)
# p610<-qplot(PC10,PC6, data = pca, colour = study)
# 
# p61<-qplot(PC1,PC6, data = pca, colour = study)
# p62<-qplot(PC2,PC6, data = pca, colour = study)
# p63<-qplot(PC3,PC6, data = pca, colour = study)
# p64<-qplot(PC4,PC6, data = pca, colour = study)
# p65<-qplot(PC5,PC6, data = pca, colour = study)
# p66<-qplot(PC6,PC6, data = pca, colour = study)
# p67<-qplot(PC7,PC6, data = pca, colour = study)
# p68<-qplot(PC8,PC6, data = pca, colour = study)
# p69<-qplot(PC9,PC6, data = pca, colour = study)
# p610<-qplot(PC10,PC6, data = pca, colour = study)
# 
# p71<-qplot(PC1,PC7, data = pca, colour = study)
# p72<-qplot(PC2,PC7, data = pca, colour = study)
# p73<-qplot(PC3,PC7, data = pca, colour = study)
# p74<-qplot(PC4,PC7, data = pca, colour = study)
# p75<-qplot(PC5,PC7, data = pca, colour = study)
# p76<-qplot(PC6,PC7, data = pca, colour = study)
# p77<-qplot(PC7,PC7, data = pca, colour = study)
# p78<-qplot(PC8,PC7, data = pca, colour = study)
# p79<-qplot(PC9,PC7, data = pca, colour = study)
# p710<-qplot(PC10,PC7, data = pca, colour = study)
# 
# p81<-qplot(PC1,PC8, data = pca, colour = study)
# p82<-qplot(PC2,PC8, data = pca, colour = study)
# p83<-qplot(PC3,PC8, data = pca, colour = study)
# p84<-qplot(PC4,PC8, data = pca, colour = study)
# p85<-qplot(PC5,PC8, data = pca, colour = study)
# p86<-qplot(PC6,PC8, data = pca, colour = study)
# p87<-qplot(PC7,PC8, data = pca, colour = study)
# p88<-qplot(PC8,PC8, data = pca, colour = study)
# p89<-qplot(PC9,PC8, data = pca, colour = study)
# p810<-qplot(PC10,PC8, data = pca, colour = study)
# 
# p91<-qplot(PC1,PC9, data = pca, colour = study)
# p92<-qplot(PC2,PC9, data = pca, colour = study)
# p93<-qplot(PC3,PC9, data = pca, colour = study)
# p94<-qplot(PC4,PC9, data = pca, colour = study)
# p95<-qplot(PC5,PC9, data = pca, colour = study)
# p96<-qplot(PC6,PC9, data = pca, colour = study)
# p97<-qplot(PC7,PC9, data = pca, colour = study)
# p98<-qplot(PC8,PC9, data = pca, colour = study)
# p99<-qplot(PC9,PC9, data = pca, colour = study)
# p910<-qplot(PC10,PC9, data = pca, colour = study)
# 
# p101<-qplot(PC1,PC10, data = pca, colour = study)
# p102<-qplot(PC2,PC10, data = pca, colour = study)
# p103<-qplot(PC3,PC10, data = pca, colour = study)
# p104<-qplot(PC4,PC10, data = pca, colour = study)
# p105<-qplot(PC5,PC10, data = pca, colour = study)
# p106<-qplot(PC6,PC10, data = pca, colour = study)
# p107<-qplot(PC7,PC10, data = pca, colour = study)
# p108<-qplot(PC8,PC10, data = pca, colour = study)
# p109<-qplot(PC9,PC10, data = pca, colour = study)
# p1010<-qplot(PC10,PC10, data = pca, colour = study)
# 
# 
# r1<-ggarrange(pna,p12,p13,p14,p15,p16,p17,p18,p19,p110,ncol=10,legend=c("none"))
# r2<-ggarrange(p21,pna,p23,p24,p25,p26,p27,p28,p29,p210,ncol=10,legend=c("none"))
# r3<-ggarrange(p31,p32,pna,p34,p35,p36,p37,p38,p39,p310,ncol=10,legend=c("none"))
# r4<-ggarrange(p41,p42,p43,pna,p45,p46,p47,p48,p49,p410,ncol=10,legend=c("none"))
# r5<-ggarrange(p51,p52,p53,p54,pna,p56,p57,p58,p59,p510,ncol=10,legend=c("none"))
# r6<-ggarrange(p61,p62,p63,p64,p65,pna,p67,p68,p69,p610,ncol=10,legend=c("none"))
# r7<-ggarrange(p71,p72,p73,p74,p75,p76,pna,p78,p79,p710,ncol=10,legend=c("none"))
# r8<-ggarrange(p81,p82,p83,p84,p85,p86,p87,pna,p89,p810,ncol=10,legend=c("none"))
# r9<-ggarrange(p91,p92,p93,p94,p95,p96,p97,p98,pna,p910,ncol=10,legend=c("none"))
# r10<-ggarrange(p101,p102,p103,p104,p105,p106,p107,p108,p109,pna,ncol=10,common.legend = TRUE,legend=c("bottom"))
# 
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_PC7_PC8_PC9_PC10_3kvariants_iibdgc_EURonly.pdf",sep=""),width=45,height=45)
# ggarrange(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,nrow=10,common.legend = TRUE,legend=c("bottom"),heights = c(1,1,1,1,1,1,1,1,1,1.2))
# dev.off()
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_PC7_PC8_PC9_PC10_3kvariants_iibdgc_EURonly.pdf ~/tmp_plots/",sep=""))


### SELF-REPORTED ASKENAZI:

# only new cohorts and cedars will have data for this:

pheno<-fread("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes-4313226.csv.gz")
pheno<-as.data.frame(pheno)
colnames(pheno)[14]<-"self_jewish"
pheno_cd<-fread(paste(path,"pheno/list_cedars_old_gwas_ids_TH.txt",sep=""),head=T,sep="\t")

pheno<-rbind(pheno[,c("sample_id","self_jewish")],pheno_cd[,c("sample_id","self_jewish")])
table(pheno$self_jewish)
#           Jewish         No Non-Jewish    Unknown        Yes
# 65002       1217      33975       1513      16224       4546

pheno$self_jewish<-as.character(pheno$self_jewish)
pheno$self_jewish[which(pheno$self_jewish %in% c("","Unknown"))]<-"Unknown"
pheno$self_jewish[which(pheno$self_jewish %in% c("No","Non-Jewish"))]<-"Non-Jewish"
pheno$self_jewish[which(pheno$self_jewish %in% c("Yes","Jewish"))]<-"Jewish"
table(pheno$self_jewish)
# Jewish Non-Jewish    Unknown 
#   5763      35488      81226



pca$self_jewish<-"Unknown"
pca$self_jewish[which(pca$FID %in% pheno$sample_id[which(pheno$self_jewish=="Non-Jewish")])]<-"Non-Jewish"
pca$self_jewish[which(pca$FID %in% pheno$sample_id[which(pheno$self_jewish=="Jewish")])]<-"Jewish"
table(pca$self_jewish)
# Jewish Non-Jewish    Unknown 
#   5493      30184     103224 


# p11<-qplot(pca$PC1,pca$PC1, data = pca, colour = self_jewish)
# p12<-qplot(pca$PC2,pca$PC1, data = pca, colour = self_jewish)
# p13<-qplot(pca$PC3,pca$PC1, data = pca, colour = self_jewish)
# p14<-qplot(pca$PC4,pca$PC1, data = pca, colour = self_jewish)
# p15<-qplot(pca$PC5,pca$PC1, data = pca, colour = self_jewish)
# p16<-qplot(pca$PC6,pca$PC1, data = pca, colour = self_jewish)
# 
# p21<-qplot(pca$PC1,pca$PC2, data = pca, colour = self_jewish)
# p22<-qplot(pca$PC2,pca$PC2, data = pca, colour = self_jewish)
# p23<-qplot(pca$PC3,pca$PC2, data = pca, colour = self_jewish)
# p24<-qplot(pca$PC4,pca$PC2, data = pca, colour = self_jewish)
# p25<-qplot(pca$PC5,pca$PC2, data = pca, colour = self_jewish)
# p26<-qplot(pca$PC6,pca$PC2, data = pca, colour = self_jewish)
# 
# p31<-qplot(pca$PC1,pca$PC3, data = pca, colour = self_jewish)
# p32<-qplot(pca$PC2,pca$PC3, data = pca, colour = self_jewish)
# p33<-qplot(pca$PC3,pca$PC3, data = pca, colour = self_jewish)
# p34<-qplot(pca$PC4,pca$PC3, data = pca, colour = self_jewish)
# p35<-qplot(pca$PC5,pca$PC3, data = pca, colour = self_jewish)
# p36<-qplot(pca$PC6,pca$PC3, data = pca, colour = self_jewish)
# 
# p41<-qplot(pca$PC1,pca$PC4, data = pca, colour = self_jewish)
# p42<-qplot(pca$PC2,pca$PC4, data = pca, colour = self_jewish)
# p43<-qplot(pca$PC3,pca$PC4, data = pca, colour = self_jewish)
# p44<-qplot(pca$PC4,pca$PC4, data = pca, colour = self_jewish)
# p45<-qplot(pca$PC5,pca$PC4, data = pca, colour = self_jewish)
# p46<-qplot(pca$PC6,pca$PC4, data = pca, colour = self_jewish)
# 
# p51<-qplot(pca$PC1,pca$PC5, data = pca, colour = self_jewish)
# p52<-qplot(pca$PC2,pca$PC5, data = pca, colour = self_jewish)
# p53<-qplot(pca$PC3,pca$PC5, data = pca, colour = self_jewish)
# p54<-qplot(pca$PC4,pca$PC5, data = pca, colour = self_jewish)
# p55<-qplot(pca$PC5,pca$PC5, data = pca, colour = self_jewish)
# p56<-qplot(pca$PC6,pca$PC5, data = pca, colour = self_jewish)
# 
# p61<-qplot(pca$PC1,pca$PC6, data = pca, colour = self_jewish)
# p62<-qplot(pca$PC2,pca$PC6, data = pca, colour = self_jewish)
# p63<-qplot(pca$PC3,pca$PC6, data = pca, colour = self_jewish)
# p64<-qplot(pca$PC4,pca$PC6, data = pca, colour = self_jewish)
# p65<-qplot(pca$PC5,pca$PC6, data = pca, colour = self_jewish)
# p66<-qplot(pca$PC6,pca$PC6, data = pca, colour = self_jewish)
# 
# 
# r1<-ggarrange(pna,p12,p13,p14,p15,p16,ncol=6,legend=c("none"))
# r2<-ggarrange(p21,pna,p23,p24,p25,p26,ncol=6,legend=c("none"))
# r3<-ggarrange(p31,p32,pna,p34,p35,p36,ncol=6,legend=c("none"))
# r4<-ggarrange(p41,p42,p43,pna,p45,p36,ncol=6,legend=c("none"))
# r5<-ggarrange(p51,p52,p53,p54,pna,p56,ncol=6,legend=c("none"))
# r6<-ggarrange(p61,p62,p63,p64,p65,pna,ncol=6,common.legend = TRUE,legend=c("bottom"))
# dev.off()
# 
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_3kvariants_iibdgc_EURonly_sr_jewish.pdf",sep=""),width=30,height=30)
# ggarrange(r1,r2,r3,r4,r5,r6,nrow=6,common.legend = TRUE,legend=c("bottom"),heights = c(1,1,1,1,1,1.2))
# dev.off()
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_3kvariants_iibdgc_EURonly_sr_jewish.pdf ~/tmp_plots/",sep=""))

# ############################################
# # PC1 - JEWISH ANCESTRY:
# 
# ylims<-c(min(pca$PC1),max(pca$PC1))
# xlims<-c(min(pca$PC2),max(pca$PC2))
# 
# p1<-ggplot(pca,aes(PC2,PC1)) +
#   geom_point(aes(color = study)) + 
#   xlim(xlims) + ylim(ylims)
# 
# p2<-ggplot(pca[which(pca$self_jewish=="Jewish"),],aes(PC2,PC1,color = self_jewish)) +
#   geom_point() + scale_color_manual(values=c("#E69F00")) +
#   xlim(xlims) + ylim(ylims)
# 
# p3<-ggplot(pca[which(pca$self_jewish=="Non-Jewish"),],aes(PC2,PC1,color = self_jewish)) +
#   geom_point() + scale_color_manual(values=c("#56B4E9")) +
#   xlim(xlims) + ylim(ylims)
# 
# p4<-ggplot(pca[which(pca$self_jewish=="Unknown"),],aes(PC2,PC1,color = self_jewish)) +
#   geom_point() + scale_color_manual(values=c("#999999")) +
#   xlim(xlims) + ylim(ylims)
# 
# p100<-ggarrange(p1,pna,pna,ncol=3)
# p234<-ggarrange(p2,p3,p4,ncol=3)
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_PC2_EURonly_sr_jewish.pdf",sep=""),width=22,height=12)
# ggarrange(p100,p234,nrow=2)
# dev.off()
# 
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_PC2_EURonly_sr_jewish.pdf ~/tmp_plots/",sep=""))
# 
# 
# ########################
# DENSITY PLOT FOR PC1:

xlims<-c(min(pca$PC1),max(pca$PC1))

p1<-ggplot(pca[which(pca$self_jewish %in% c("Jewish","Non-Jewish")),], aes(x=PC1, fill=self_jewish)) +
  geom_density(alpha=0.4) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  xlim(xlims)

p2<-ggplot(pca[which(pca$self_jewish %in% c("Unknown")),], aes(x=PC1, fill=self_jewish)) +
  geom_density(alpha=0.4) + scale_fill_manual(values=c("#999999")) +
  xlim(xlims)

pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_density_EURonly_sr_jewish.pdf",sep=""),width=14,height=10)
ggarrange(p1,p2,nrow=2)
dev.off()

system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_density_EURonly_sr_jewish.pdf ~/tmp_plots/",sep=""))

/Users/username/Desktop/IIBDGC_GWAS_call_20230123.pptx

# #### histogram
# 
# fd=function(x) {
#   n=length(x)
#   r=IQR(x)
#   2*r/n^(1/3)
# }
# 
# fd_bin<-fd(pca$PC1)
# 
# p1<-ggplot(pca[which(pca$self_jewish %in% c("Jewish","Non-Jewish")),], aes(x=PC1,color=self_jewish, fill=self_jewish)) +
#   geom_histogram(position="identity", alpha=0.5, binwidth=fd_bin) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) + 
#   scale_color_manual(values=c("#E69F00", "#56B4E9")) +
#   xlim(xlims)
# 
# p2<-ggplot(pca[which(pca$self_jewish %in% c("Unknown")),], aes(x=PC1,color=self_jewish, fill=self_jewish)) +
#   geom_histogram(position="identity", alpha=0.5, binwidth=fd_bin) + scale_fill_manual(values=c("#999999")) + 
#   scale_color_manual(values=c("#999999")) +
#   xlim(xlims)
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_histogram_EURonly_sr_jewish.pdf",sep=""),width=14,height=10)
# ggarrange(p1,p2,nrow=2)
# dev.off()
# 
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_histogram_EURonly_sr_jewish.pdf ~/tmp_plots/",sep=""))
# 
# 
### ESTIMATE CUTOFF POINT:

cutpoint<- 0.045

### density:

p1<-ggplot(pca[which(pca$self_jewish %in% c("Jewish","Non-Jewish")),], aes(x=PC1, fill=self_jewish)) +
  geom_density(alpha=0.4) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) + geom_vline(xintercept = cutpoint, linetype="dotted", color = "darkgrey") +
  xlim(xlims)
p2<-ggplot(pca[which(pca$self_jewish %in% c("Unknown")),], aes(x=PC1, fill=self_jewish)) +
  geom_density(alpha=0.4) + scale_fill_manual(values=c("#999999")) + geom_vline(xintercept = cutpoint, linetype="dotted", color = "black") +
  xlim(xlims)

pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_density_EURonly_sr_jewish_cutoff.pdf",sep=""),width=14,height=10)
ggarrange(p1,p2,nrow=2)
dev.off()

system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_density_EURonly_sr_jewish_cutoff.pdf ~/tmp_plots/",sep=""))
# 
# 
# ### histogram:
# 
# p1<-ggplot(pca[which(pca$self_jewish %in% c("Jewish","Non-Jewish")),], aes(x=PC1,color=self_jewish, fill=self_jewish)) +
#   geom_histogram(position="identity", alpha=0.5, binwidth=fd_bin) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) + 
#   scale_color_manual(values=c("#E69F00", "#56B4E9")) + geom_vline(xintercept = cutpoint, linetype="dotted", color = "black") +
#   xlim(xlims)
# 
# p2<-ggplot(pca[which(pca$self_jewish %in% c("Unknown")),], aes(x=PC1,color=self_jewish, fill=self_jewish)) +
#   geom_histogram(position="identity", alpha=0.5, binwidth=fd_bin) + scale_fill_manual(values=c("#999999")) + 
#   scale_color_manual(values=c("#999999")) + geom_vline(xintercept = cutpoint, linetype="dotted", color = "black") +
#   xlim(xlims)
# 
# pca$tmp<-"All IIBDG EUR"
# 
# p3<-ggplot(pca, aes(x=PC1,color=tmp, fill=tmp)) +
#   geom_histogram(position="identity", alpha=0.5, binwidth=fd_bin) + scale_fill_manual(values=c("lightblue4")) + 
#   scale_color_manual(values=c("lightblue4")) + geom_vline(xintercept = cutpoint, linetype="dotted", color = "black") +
#   xlim(xlims)
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_histogram_EURonly_sr_jewish_cutoff.pdf",sep=""),width=14,height=15)
# ggarrange(p1,p2,p3,nrow=3)
# dev.off()
# 
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_histogram_EURonly_sr_jewish_cutoff.pdf ~/tmp_plots/",sep=""))
# 


# reclassify:

pca$pca_jewish<-NA
pca$pca_jewish[which(pca$PC1<cutpoint)]<-"Non-Jewish"
pca$pca_jewish[which(pca$PC1>=cutpoint)]<-"Jewish"
table(pca$pca_jewish)
# Jewish Non-Jewish 
# 10274     128627    

table(pca$pca_jewish,pca$self_jewish)
#            Jewish Non-Jewish Unknown
# Jewish       3858        712    5704
# Non-Jewish   1635      29472   97520


3858/sum(3858,1635)
# [1] 0.7023484 sensitivity

29472/sum(29472,712)
# [1] 0.9764113 specificity



3858/sum(3858,712)
# [1] 0.8442013# PPV 

29472/sum(29472,1635)
# [1] 0.9474395 # NPV




# ###################################################################
# # PLOT CASES AND CONTROLS PER STUDY, FOR JEWISH SR, no duplicates:

pca$pheno<-NA
pca$pheno[which(pca$FID %in% fam$V1[which(fam$V6==1)])]<-"1"
pca$pheno[which(pca$FID %in% fam$V1[which(fam$V6==2)])]<-"2"

table(pca$cohort[which(pca$pheno=="2")],pca$self_jewish[which(pca$pheno=="2")],useNA="ifany")
table(pca$cohort[which(pca$pheno=="1")],pca$self_jewish[which(pca$pheno=="1")],useNA="ifany")


#######################################################
# PLOT CASES AND CONTROLS PER STUDY, FOR JEWISH PCA:


eur_nodup<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/list_eur_ancestry_samples_noDuplicates",sep=""),head=F)

x1<-as.data.frame.matrix(table(pca$cohort[which(pca$pheno=="2" & pca$FID %in% eur_nodup$V1)],
                               pca$pca_jewish[which(pca$pheno=="2" & pca$FID %in% eur_nodup$V1)],useNA="ifany"))
colnames(x1)<-paste(colnames(x1),"n_cases",sep="_")
x1$study<-rownames(x1)

x2<-as.data.frame.matrix(table(pca$cohort[which(pca$pheno=="1" & pca$FID %in% eur_nodup$V1)],
                               pca$pca_jewish[which(pca$pheno=="1" & pca$FID %in% eur_nodup$V1)],useNA="ifany"))
colnames(x2)<-paste(colnames(x2),"n_ctr",sep="_")
x2$study<-rownames(x2)

x<-merge(x1,x2,by="study",all=T)
write.table(x,paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_ancestries_per_study_case_control_noDuplicates_EUR_jewish.tsv",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")
  
# # PC2 - PC3
# 
# xlims<-c(min(pca$PC3),max(pca$PC3))
# ylims<-c(min(pca$PC2),max(pca$PC2))
# 
# p1<-ggplot(pca,aes(PC3,PC2)) +
#   geom_point(aes(color = study)) + 
#   xlim(xlims) + ylim(ylims)
# 
# p2<-ggplot(pca,aes(PC3,PC2)) +
#   geom_point(aes(color = study)) + 
#   xlim(xlims) + ylim(ylims)
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_PC2_PC3_EURonly.pdf",sep=""),width=7,height=7)
# p23
# dev.off()
# 
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_PC2_PC3_EURonly.pdf ~/tmp_plots/",sep=""))



# # PC2 - PC3 -  all studies
# 
# table(pca$study)
# # Spain       Italy    Slovenia     Belgium     Germany Netherlands 
# # 4935         947         264        8534       16937        4352 
# # UK   Lithuania      Sweden      Norway     Finland   Australia 
# # 35114        2198        2584         547         443        1251 
# # CHOP         USA 
# # 8475       26531 
# 
# 
# # to emulate ggplot2 palette:
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# 
# n = length(levels(pca$study))
# cols = gg_color_hue(n)
# 
# pca$PC2_inv<-pca$PC2*-1
# 
# xlims<-c(min(pca$PC2_inv),max(pca$PC2_inv))
# ylims<-c(min(pca$PC3),max(pca$PC3))
# 
# p0<-ggplot(pca,aes(PC2_inv,PC3)) +
#   geom_point(aes(color = study)) + 
#   xlim(xlims) + ylim(ylims)
# 
# for (i in 1:length(levels(pca$study))) {
#   print(i)
#   assign(paste("p",i,sep=""),ggplot(pca[which(pca$study==levels(pca$study)[i]),],aes(PC2_inv,PC3,color=study)) +
#            geom_point() + scale_color_manual(values=cols[i]) + 
#            xlim(xlims) + ylim(ylims))
# }
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_PC2_PC3_EURonly_per_study.pdf",sep=""),width=30,height=27)
# print(ggarrange(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,legend=c("right"),ncol=4,nrow=5))
# dev.off()
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/ibdgc_PC2_PC3_EURonly_per_study.pdf ~/tmp_plots/",sep=""))
# 

##### create two additional lists of samples:
# - EUR samples all
# - EUR samples non-jewish
# - EUR samples Jewish


write.table(pca[which(pca$pca_jewish=="Jewish"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eur_jewish_ancestry_samples_withDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")
write.table(pca[which(pca$pca_jewish=="Non-Jewish"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eur_nonjewish_ancestry_samples_withDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")
write.table(pca,paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_EUR_withDuplicates_ref_pcs_new_projection_edited.sscore",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")


eur_nodup<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/list_eur_ancestry_samples_noDuplicates",sep=""),head=F)
pca_nondup<-pca[which(pca$FID %in% eur_nodup$V1),]
write.table(pca_nondup[which(pca_nondup$pca_jewish=="Non-Jewish"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eur_nonjewish_ancestry_samples_noDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")



### some old studies keep only small N of samples - repeat with non duplicates per cohort


# no german_illu, or chop
cohorts<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
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

length(cohorts)
# [1] 58

for (i in 1:length(cohorts)){
  a<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.fam",sep=""),head=F)
  a$cohort<-cohorts[i]
  if(i==1){
    fam<-a
  }else{
    fam<-rbind(a,fam)
  }
}
dim(fam)
# [1] 144356      7

pca_nondup_2<-pca[which(pca$FID %in% fam$V1),]
dim(pca_nondup_2)
# [1] 136492     20
write.table(pca_nondup_2[which(pca_nondup_2$pca_jewish=="Non-Jewish"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eur_nonjewish_ancestry_samples_noDuplicates_perstudy",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")


write.table(pca_nondup_2[which(pca_nondup_2$pca_jewish=="Jewish"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eur_jewish_ancestry_samples_noDuplicates_perstudy",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

#######################


#############################################################################################################################################################
# CREATE AN UNIQUE PCA FILE:

pca_eur<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_EUR_withDuplicates_ref_pcs_new_projection_edited.sscore",sep=""),head=T,sep="\t")
dim(pca_eur)
# [1] 138901    19

pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/1000GP_all_b37_plus_iibdgc_merged_withDuplicates_noHighLD_noassoclogP4_pruned_new_projection_edited_withDuplicates.sscore",sep=""),head=T,sep="\t")
dim(pca)
# [1] 149457     21

colnames(pca_eur)[6:15]<-paste(colnames(pca_eur)[6:15],"_EUR",sep="")

pca<-merge(pca,pca_eur[,c(1,6:15,18,19)],by="FID",all.x=T)
pca<-pca[,c(1:3,6:32)]

table(pca$super_pop)
# AFR    AMR    EAS    EUR IIBDGC    SAS 
# 661    347    504    503 146953    489

write.table(pca,paste(path,"pre_imputation/QC/pca_1000gp/pca_1000gp_iibdgc_pca_pcaEur_withDuplicates.tsv",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")

