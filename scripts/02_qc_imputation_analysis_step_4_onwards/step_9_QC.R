# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
##########
# 9 - QC #
##########


####################################
# 9.1 - REMOVE SAMPLES CallPP <0.80%

path_gwas=/path/to/ibdgwas/IIBDGC/
studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
  
MEM=500
  
for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_22_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_22_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry \
--allow-no-sex \
--mind 0.20 \
--make-bed --out ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8"
done
  
for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_22_${i} | grep -E "completed"
done


for i in ${studies[@]}
do echo ${i} && tail -150 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_22_${i} | grep -E "removed"
done

############

# all_hce
# 1 person removed due to missing genotype data (--mind).
# niddk_old_gwas
# 13 people removed due to missing genotype data (--mind).
# australia_omniexome
# 0 people removed due to missing genotype data (--mind).
# gwas1
# 0 people removed due to missing genotype data (--mind).
# gwas2
# 0 people removed due to missing genotype data (--mind).
# pittsburgh_gsa
# 0 people removed due to missing genotype data (--mind).
# spain_gsa
# 0 people removed due to missing genotype data (--mind).
# italy_gsa
# 18 people removed due to missing genotype data (--mind).
# kiel_austria_sibdcs_gsa
# 4 people removed due to missing genotype data (--mind).
# netherlands_gsa
# 34 people removed due to missing genotype data (--mind).
# slovenia_gsa
# 1 person removed due to missing genotype data (--mind).
# sweden_gsa
# 0 people removed due to missing genotype data (--mind).
# niddk_broad_gsa
# 0 people removed due to missing genotype data (--mind).
# niddk_feinstein_gsa
# 0 people removed due to missing genotype data (--mind).
# basque_gsa
# 2 people removed due to missing genotype data (--mind).
# lithuania_gsa
# 0 people removed due to missing genotype data (--mind).
# belgium_louis_gsa
# 0 people removed due to missing genotype data (--mind).
# belgium_franchimont_gsa
# 0 people removed due to missing genotype data (--mind).
# belgium_vermeire_gsa
# 0 people removed due to missing genotype data (--mind).
# prism_nfe_gsa
# 0 people removed due to missing genotype data (--mind).
# prism_nfe_gwas
# 13 people removed due to missing genotype data (--mind).
# finland_illugwas
# 0 people removed due to missing genotype data (--mind).
# german_affy6_old_gwas
# 0 people removed due to missing genotype data (--mind).
# norway_affy6_old_gwas
# 0 people removed due to missing genotype data (--mind).
# belgium_inf1_old_gwas
# 0 people removed due to missing genotype data (--mind).
# belgium_inf2_old_gwas
# 0 people removed due to missing genotype data (--mind).
# cedars_370k_old_gwas
# 0 people removed due to missing genotype data (--mind).
# cedars_610k_old_gwas
# 0 people removed due to missing genotype data (--mind).
# cedars_omni_old_gwas
# 0 people removed due to missing genotype data (--mind).
# swedish_uc_old_gwas
# 0 people removed due to missing genotype data (--mind).
# mccauley_gsa
# 0 people removed due to missing genotype data (--mind).
# ccfa_gsa
# 0 people removed due to missing genotype data (--mind).
# cedars_gsa
# 1 person removed due to missing genotype data (--mind).
# bernstein_gsa
# 0 people removed due to missing genotype data (--mind).
# farkkila_gsa
# 0 people removed due to missing genotype data (--mind).
# franchimont_gsa
# 0 people removed due to missing genotype data (--mind).
# franke_gsa
# 0 people removed due to missing genotype data (--mind).
# helmsley_prism_gsa
# 0 people removed due to missing genotype data (--mind).
# helmsley_xavier_prism_gsa
# 0 people removed due to missing genotype data (--mind).
# hyams_protect_gsa
# 0 people removed due to missing genotype data (--mind).
# lewis_sparc_gsa
# 0 people removed due to missing genotype data (--mind).
# mccauley_new_gsa
# 0 people removed due to missing genotype data (--mind).
# mcgovern_gsa
# 0 people removed due to missing genotype data (--mind).
# moayyedi_imagine_gsa
# 0 people removed due to missing genotype data (--mind).
# newberry_share_gsa
# 0 people removed due to missing genotype data (--mind).
# niddk_cho_gsa
# 0 people removed due to missing genotype data (--mind).
# niddk_duerr_gsa
# 0 people removed due to missing genotype data (--mind).
# niddk_rioux_gsa
# 0 people removed due to missing genotype data (--mind).
# niddk_silverberg_gsa
# 0 people removed due to missing genotype data (--mind).
# palotie_hus_gsa
# 0 people removed due to missing genotype data (--mind).
# pekow_share_gsa
# 0 people removed due to missing genotype data (--mind).
# rioux_igenomed_gsa
# 0 people removed due to missing genotype data (--mind).
# sands_msccr_gsa
# 0 people removed due to missing genotype data (--mind).
# stampfer_gsa
# 0 people removed due to missing genotype data (--mind).
# vermeire_gsa
# 0 people removed due to missing genotype data (--mind).
# weersma_gsa
# 0 people removed due to missing genotype data (--mind).
# xavier_prism_gsa
# 0 people removed due to missing genotype data (--mind).
# xavier_share_gsa
# 0 people removed due to missing genotype data (--mind).

############

####################################
# 9.2- REMOVE VARIANTS CallPP <0.80%

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_23_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_23_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8 \
--geno 0.20 \
--make-bed --out ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8"
done


for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_23_${i} | grep -E "completed"
done


for i in ${studies[@]}
do echo ${i} && tail -150 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_23_${i} | grep -E "removed"
done

############

# all_hce
# 612 variants removed due to missing genotype data (--geno).
# niddk_old_gwas
# 126 variants removed due to missing genotype data (--geno).
# australia_omniexome
# 41688 variants removed due to missing genotype data (--geno).
# gwas1
# 0 variants removed due to missing genotype data (--geno).
# gwas2
# 275 variants removed due to missing genotype data (--geno).
# pittsburgh_gsa
# 799 variants removed due to missing genotype data (--geno).
# spain_gsa
# 1219 variants removed due to missing genotype data (--geno).
# italy_gsa
# 2395 variants removed due to missing genotype data (--geno).
# kiel_austria_sibdcs_gsa
# 43025 variants removed due to missing genotype data (--geno).
# netherlands_gsa
# 1699 variants removed due to missing genotype data (--geno).
# slovenia_gsa
# 2629 variants removed due to missing genotype data (--geno).
# sweden_gsa
# 1554 variants removed due to missing genotype data (--geno).
# niddk_broad_gsa
# 3000 variants removed due to missing genotype data (--geno).
# niddk_feinstein_gsa
# 975 variants removed due to missing genotype data (--geno).
# basque_gsa
# 486 variants removed due to missing genotype data (--geno).
# lithuania_gsa
# 6995 variants removed due to missing genotype data (--geno).
# belgium_louis_gsa
# 2073 variants removed due to missing genotype data (--geno).
# belgium_franchimont_gsa
# 4704 variants removed due to missing genotype data (--geno).
# belgium_vermeire_gsa
# 4142 variants removed due to missing genotype data (--geno).
# prism_nfe_gsa
# 3516 variants removed due to missing genotype data (--geno).
# prism_nfe_gwas
# 73 variants removed due to missing genotype data (--geno).
# finland_illugwas
# 374 variants removed due to missing genotype data (--geno).
# german_affy6_old_gwas
# 2842 variants removed due to missing genotype data (--geno).
# norway_affy6_old_gwas
# 21975 variants removed due to missing genotype data (--geno).
# belgium_inf1_old_gwas
# 0 variants removed due to missing genotype data (--geno).
# belgium_inf2_old_gwas
# 0 variants removed due to missing genotype data (--geno).
# cedars_370k_old_gwas
# 342 variants removed due to missing genotype data (--geno).
# cedars_610k_old_gwas
# 622 variants removed due to missing genotype data (--geno).
# cedars_omni_old_gwas
# 263 variants removed due to missing genotype data (--geno).
# swedish_uc_old_gwas
# 0 variants removed due to missing genotype data (--geno).
# mccauley_gsa
# 3542 variants removed due to missing genotype data (--geno).
# ccfa_gsa
# 1616 variants removed due to missing genotype data (--geno).
# cedars_gsa
# 4876 variants removed due to missing genotype data (--geno).
# bernstein_gsa
# 1715 variants removed due to missing genotype data (--geno).
# farkkila_gsa
# 4377 variants removed due to missing genotype data (--geno).
# franchimont_gsa
# 3232 variants removed due to missing genotype data (--geno).
# franke_gsa
# 2857 variants removed due to missing genotype data (--geno).
# helmsley_prism_gsa
# 70 variants removed due to missing genotype data (--geno).
# helmsley_xavier_prism_gsa
# 56 variants removed due to missing genotype data (--geno).
# hyams_protect_gsa
# 2983 variants removed due to missing genotype data (--geno).
# lewis_sparc_gsa
# 1603 variants removed due to missing genotype data (--geno).
# mccauley_new_gsa
# 1716 variants removed due to missing genotype data (--geno).
# mcgovern_gsa
# 3345 variants removed due to missing genotype data (--geno).
# moayyedi_imagine_gsa
# 851 variants removed due to missing genotype data (--geno).
# newberry_share_gsa
# 1156 variants removed due to missing genotype data (--geno).
# niddk_cho_gsa
# 1816 variants removed due to missing genotype data (--geno).
# niddk_duerr_gsa
# 1448 variants removed due to missing genotype data (--geno).
# niddk_rioux_gsa
# 1773 variants removed due to missing genotype data (--geno).
# niddk_silverberg_gsa
# 2080 variants removed due to missing genotype data (--geno).
# palotie_hus_gsa
# 2122 variants removed due to missing genotype data (--geno).
# pekow_share_gsa
# 885 variants removed due to missing genotype data (--geno).
# rioux_igenomed_gsa
# 1595 variants removed due to missing genotype data (--geno).
# sands_msccr_gsa
# 1010 variants removed due to missing genotype data (--geno).
# stampfer_gsa
# 3640 variants removed due to missing genotype data (--geno).
# vermeire_gsa
# 4206 variants removed due to missing genotype data (--geno).
# weersma_gsa
# 699 variants removed due to missing genotype data (--geno).
# xavier_prism_gsa
# 3084 variants removed due to missing genotype data (--geno).
# xavier_share_gsa
# 3994 variants removed due to missing genotype data (--geno).

############

####################################
# 9.3 - REMOVE SAMPLES CallPP <0.95%

MEM=500

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_24_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_24_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8 \
--allow-no-sex \
--mind 0.05 \
--make-bed --out ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_24_${i} | grep -E "completed"
done


for i in ${studies[@]}
do echo ${i} && tail -150 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_24_${i} | grep -E "removed"
done

############
# 
# all_hce
# 113 people removed due to missing genotype data (--mind).
# niddk_old_gwas
# 42 people removed due to missing genotype data (--mind).
# australia_omniexome
# 0 people removed due to missing genotype data (--mind).
# gwas1
# 0 people removed due to missing genotype data (--mind).
# gwas2
# 0 people removed due to missing genotype data (--mind).
# pittsburgh_gsa
# 26 people removed due to missing genotype data (--mind).
# spain_gsa
# 0 people removed due to missing genotype data (--mind).
# italy_gsa
# 9 people removed due to missing genotype data (--mind).
# kiel_austria_sibdcs_gsa
# 243 people removed due to missing genotype data (--mind).
# netherlands_gsa
# 63 people removed due to missing genotype data (--mind).
# slovenia_gsa
# 5 people removed due to missing genotype data (--mind).
# sweden_gsa
# 0 people removed due to missing genotype data (--mind).
# niddk_broad_gsa
# 14 people removed due to missing genotype data (--mind).
# niddk_feinstein_gsa
# 22 people removed due to missing genotype data (--mind).
# basque_gsa
# 6 people removed due to missing genotype data (--mind).
# lithuania_gsa
# 0 people removed due to missing genotype data (--mind).
# belgium_louis_gsa
# 1 person removed due to missing genotype data (--mind).
# belgium_franchimont_gsa
# 0 people removed due to missing genotype data (--mind).
# belgium_vermeire_gsa
# 9 people removed due to missing genotype data (--mind).
# prism_nfe_gsa
# 0 people removed due to missing genotype data (--mind).
# prism_nfe_gwas
# 16 people removed due to missing genotype data (--mind).
# finland_illugwas
# 1 person removed due to missing genotype data (--mind).
# german_affy6_old_gwas
# 4 people removed due to missing genotype data (--mind).
# norway_affy6_old_gwas
# 0 people removed due to missing genotype data (--mind).
# belgium_inf1_old_gwas
# 0 people removed due to missing genotype data (--mind).
# belgium_inf2_old_gwas
# 0 people removed due to missing genotype data (--mind).
# cedars_370k_old_gwas
# 0 people removed due to missing genotype data (--mind).
# cedars_610k_old_gwas
# 0 people removed due to missing genotype data (--mind).
# cedars_omni_old_gwas
# 0 people removed due to missing genotype data (--mind).
# swedish_uc_old_gwas
# 0 people removed due to missing genotype data (--mind).
# mccauley_gsa
# 0 people removed due to missing genotype data (--mind).
# ccfa_gsa
# 0 people removed due to missing genotype data (--mind).
# cedars_gsa
# 1 person removed due to missing genotype data (--mind).
# bernstein_gsa
# 0 people removed due to missing genotype data (--mind).
# farkkila_gsa
# 0 people removed due to missing genotype data (--mind).
# franchimont_gsa
# 1 person removed due to missing genotype data (--mind).
# franke_gsa
# 0 people removed due to missing genotype data (--mind).
# helmsley_prism_gsa
# 7 people removed due to missing genotype data (--mind).
# helmsley_xavier_prism_gsa
# 0 people removed due to missing genotype data (--mind).
# hyams_protect_gsa
# 0 people removed due to missing genotype data (--mind).
# lewis_sparc_gsa
# 0 people removed due to missing genotype data (--mind).
# mccauley_new_gsa
# 0 people removed due to missing genotype data (--mind).
# mcgovern_gsa
# 1 person removed due to missing genotype data (--mind).
# moayyedi_imagine_gsa
# 0 people removed due to missing genotype data (--mind).
# newberry_share_gsa
# 2 people removed due to missing genotype data (--mind).
# niddk_cho_gsa
# 0 people removed due to missing genotype data (--mind).
# niddk_duerr_gsa
# 0 people removed due to missing genotype data (--mind).
# niddk_rioux_gsa
# 0 people removed due to missing genotype data (--mind).
# niddk_silverberg_gsa
# 0 people removed due to missing genotype data (--mind).
# palotie_hus_gsa
# 0 people removed due to missing genotype data (--mind).
# pekow_share_gsa
# 0 people removed due to missing genotype data (--mind).
# rioux_igenomed_gsa
# 0 people removed due to missing genotype data (--mind).
# sands_msccr_gsa
# 0 people removed due to missing genotype data (--mind).
# stampfer_gsa
# 8 people removed due to missing genotype data (--mind).
# vermeire_gsa
# 9 people removed due to missing genotype data (--mind).
# weersma_gsa
# 0 people removed due to missing genotype data (--mind).
# xavier_prism_gsa
# 0 people removed due to missing genotype data (--mind).
# xavier_share_gsa
# 0 people removed due to missing genotype data (--mind).

############

#####################################
# 9.4 - REMOVE VARIANTS CallPP <0.95%

MEM=500

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_25_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_25_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95 \
--allow-no-sex \
--geno 0.05 \
--make-bed --out ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_25_${i} | grep -E "completed"
done


for i in ${studies[@]}
do echo ${i} && tail -150 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_25_${i} | grep -E "removed"
done

############

# all_hce
# 6662 variants removed due to missing genotype data (--geno).
# niddk_old_gwas
# 6007 variants removed due to missing genotype data (--geno).
# australia_omniexome
# 4090 variants removed due to missing genotype data (--geno).
# gwas1
# 370 variants removed due to missing genotype data (--geno).
# gwas2
# 1311 variants removed due to missing genotype data (--geno).
# pittsburgh_gsa
# 14961 variants removed due to missing genotype data (--geno).
# spain_gsa
# 2498 variants removed due to missing genotype data (--geno).
# italy_gsa
# 5211 variants removed due to missing genotype data (--geno).
# kiel_austria_sibdcs_gsa
# 36031 variants removed due to missing genotype data (--geno).
# netherlands_gsa
# 6451 variants removed due to missing genotype data (--geno).
# slovenia_gsa
# 20243 variants removed due to missing genotype data (--geno).
# sweden_gsa
# 7213 variants removed due to missing genotype data (--geno).
# niddk_broad_gsa
# 22583 variants removed due to missing genotype data (--geno).
# niddk_feinstein_gsa
# 4013 variants removed due to missing genotype data (--geno).
# basque_gsa
# 2850 variants removed due to missing genotype data (--geno).
# lithuania_gsa
# 13360 variants removed due to missing genotype data (--geno).
# belgium_louis_gsa
# 7733 variants removed due to missing genotype data (--geno).
# belgium_franchimont_gsa
# 13973 variants removed due to missing genotype data (--geno).
# belgium_vermeire_gsa
# 11166 variants removed due to missing genotype data (--geno).
# prism_nfe_gsa
# 8624 variants removed due to missing genotype data (--geno).
# prism_nfe_gwas
# 780 variants removed due to missing genotype data (--geno).
# finland_illugwas
# 1722 variants removed due to missing genotype data (--geno).
# german_affy6_old_gwas
# 46686 variants removed due to missing genotype data (--geno).
# norway_affy6_old_gwas
# 67855 variants removed due to missing genotype data (--geno).
# belgium_inf1_old_gwas
# 1 variant removed due to missing genotype data (--geno).
# belgium_inf2_old_gwas
# 0 variants removed due to missing genotype data (--geno).
# cedars_370k_old_gwas
# 776 variants removed due to missing genotype data (--geno).
# cedars_610k_old_gwas
# 1472 variants removed due to missing genotype data (--geno).
# cedars_omni_old_gwas
# 921 variants removed due to missing genotype data (--geno).
# swedish_uc_old_gwas
# 0 variants removed due to missing genotype data (--geno).
# mccauley_gsa
# 8147 variants removed due to missing genotype data (--geno).
# ccfa_gsa
# 5899 variants removed due to missing genotype data (--geno).
# cedars_gsa
# 12710 variants removed due to missing genotype data (--geno).
# bernstein_gsa
# 7520 variants removed due to missing genotype data (--geno).
# farkkila_gsa
# 10553 variants removed due to missing genotype data (--geno).
# franchimont_gsa
# 12946 variants removed due to missing genotype data (--geno).
# franke_gsa
# 7574 variants removed due to missing genotype data (--geno).
# helmsley_prism_gsa
# 639 variants removed due to missing genotype data (--geno).
# helmsley_xavier_prism_gsa
# 651 variants removed due to missing genotype data (--geno).
# hyams_protect_gsa
# 6268 variants removed due to missing genotype data (--geno).
# lewis_sparc_gsa
# 6258 variants removed due to missing genotype data (--geno).
# mccauley_new_gsa
# 7458 variants removed due to missing genotype data (--geno).
# mcgovern_gsa
# 11166 variants removed due to missing genotype data (--geno).
# moayyedi_imagine_gsa
# 3667 variants removed due to missing genotype data (--geno).
# newberry_share_gsa
# 4967 variants removed due to missing genotype data (--geno).
# niddk_cho_gsa
# 8530 variants removed due to missing genotype data (--geno).
# niddk_duerr_gsa
# 7353 variants removed due to missing genotype data (--geno).
# niddk_rioux_gsa
# 8455 variants removed due to missing genotype data (--geno).
# niddk_silverberg_gsa
# 8745 variants removed due to missing genotype data (--geno).
# palotie_hus_gsa
# 6634 variants removed due to missing genotype data (--geno).
# pekow_share_gsa
# 3817 variants removed due to missing genotype data (--geno).
# rioux_igenomed_gsa
# 4896 variants removed due to missing genotype data (--geno).
# sands_msccr_gsa
# 4876 variants removed due to missing genotype data (--geno).
# stampfer_gsa
# 17208 variants removed due to missing genotype data (--geno).
# vermeire_gsa
# 11135 variants removed due to missing genotype data (--geno).
# weersma_gsa
# 3913 variants removed due to missing genotype data (--geno).
# xavier_prism_gsa
# 8489 variants removed due to missing genotype data (--geno).
# xavier_share_gsa
# 9487 variants removed due to missing genotype data (--geno).

############


#######################################################
# 9.5 - REMOVE VARIANTS FREQ <0.01 AND CallPP <0.98%


MEM=500

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_26_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_26_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95 \
--allow-no-sex \
--missing --out ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95"
done


for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_27_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_27_${i} \
"//path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95 \
--allow-no-sex \
--freq --out ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_26_${i} | grep -E "completed"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_27_${i} | grep -E "completed"
done


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
  sample_miss<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95.imiss",sep=""),head=T)
  var_miss<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95.lmiss",sep=""),head=T)
  frq<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95.frq",sep=""),head=T)
  
  var<-merge(frq[,c(2:6)],var_miss,by="SNP")
  var.1<-var[which(var$MAF<0.01 & var$F_MISS>0.02),]
  
  print(nrow(var.1))
  
  write.table(var.1[,"SNP",drop=F],paste(path,"pre_imputation/QC/",studies[j],"/list_monomorphic_vcr098_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

  rm(var,var.1,sample_miss,var_miss,frq)
}

q("no")

############

# [1] "all_hce"
# [1] 783
# [1] "niddk_old_gwas"
# [1] 10
# [1] "australia_omniexome"
# [1] 755
# [1] "gwas1"
# [1] 14
# [1] "gwas2"
# [1] 113
# [1] "pittsburgh_gsa"
# [1] 3863
# [1] "spain_gsa"
# [1] 259
# [1] "italy_gsa"
# [1] 1234
# [1] "kiel_austria_sibdcs_gsa"
# [1] 5600
# [1] "netherlands_gsa"
# [1] 4369
# [1] "slovenia_gsa"
# [1] 3951
# [1] "sweden_gsa"
# [1] 1375
# [1] "niddk_broad_gsa"
# [1] 3660
# [1] "niddk_feinstein_gsa"
# [1] 1601
# [1] "basque_gsa"
# [1] 983
# [1] "lithuania_gsa"
# [1] 2070
# [1] "belgium_louis_gsa"
# [1] 1757
# [1] "belgium_franchimont_gsa"
# [1] 2147
# [1] "belgium_vermeire_gsa"
# [1] 2291
# [1] "prism_nfe_gsa"
# [1] 1286
# [1] "prism_nfe_gwas"
# [1] 34
# [1] "finland_illugwas"
# [1] 38
# [1] "german_affy6_old_gwas"
# [1] 13823
# [1] "norway_affy6_old_gwas"
# [1] 7572
# [1] "belgium_inf1_old_gwas"
# [1] 9
# [1] "belgium_inf2_old_gwas"
# [1] 7
# [1] "cedars_370k_old_gwas"
# [1] 26
# [1] "cedars_610k_old_gwas"
# [1] 117
# [1] "cedars_omni_old_gwas"
# [1] 90
# [1] "swedish_uc_old_gwas"
# [1] 0
# [1] "mccauley_gsa"
# [1] 1372
# [1] "ccfa_gsa"
# [1] 1373
# [1] "cedars_gsa"
# [1] 2440
# [1] "bernstein_gsa"
# [1] 1270
# [1] "farkkila_gsa"
# [1] 931
# [1] "franchimont_gsa"
# [1] 2641
# [1] "franke_gsa"
# [1] 1323
# [1] "helmsley_prism_gsa"
# [1] 18
# [1] "helmsley_xavier_prism_gsa"
# [1] 10
# [1] "hyams_protect_gsa"
# [1] 838
# [1] "lewis_sparc_gsa"
# [1] 1485
# [1] "mccauley_new_gsa"
# [1] 1419
# [1] "mcgovern_gsa"
# [1] 2281
# [1] "moayyedi_imagine_gsa"
# [1] 1126
# [1] "newberry_share_gsa"
# [1] 1167
# [1] "niddk_cho_gsa"
# [1] 2018
# [1] "niddk_duerr_gsa"
# [1] 2041
# [1] "niddk_rioux_gsa"
# [1] 1763
# [1] "niddk_silverberg_gsa"
# [1] 2067
# [1] "palotie_hus_gsa"
# [1] 1179
# [1] "pekow_share_gsa"
# [1] 960
# [1] "rioux_igenomed_gsa"
# [1] 878
# [1] "sands_msccr_gsa"
# [1] 1133
# [1] "stampfer_gsa"
# [1] 4689
# [1] "vermeire_gsa"
# [1] 2254
# [1] "weersma_gsa"
# [1] 1154
# [1] "xavier_prism_gsa"
# [1] 1582
# [1] "xavier_share_gsa"
# [1] 1439

############


for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_28_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_28_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95 \
--allow-no-sex \
--exclude ${path_gwas}/pre_imputation/QC/${i}/list_monomorphic_vcr098_var_exclude \
--make-bed \
--out ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01"
done

##############

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_28_${i} | grep -E "completed"
done

