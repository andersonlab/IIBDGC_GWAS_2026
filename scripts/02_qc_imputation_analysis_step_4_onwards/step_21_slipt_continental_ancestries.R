# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
########################################################
# 21 - SPLIT COHORTS INTO BROAD CONTINENTAL ANCESTRIES #
########################################################

#####################################
# 21.1 SPLIT (INCLUDING Duplicates)


studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

ancestry=(eur)

# EUR

for j in ${ancestry[@]}
do
for i in ${studies[@]}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_split_ancestry_eur_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_split_ancestry_eur_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged \
--keep-allele-order --allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_${j}_ancestry_samples_withDuplicates \
--freq \
--make-bed --out ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${j}"
done
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_split_ancestry_eur_${i} | grep -E "completed"
done

# non-EUR

for j in ${ancestry[@]}
do
for i in ${studies[@]}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_split_ancestry_noneur_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_split_ancestry_noneur_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged \
--keep-allele-order --allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/pca_1000gp/list_${j}_ancestry_samples_withDuplicates \
--freq \
--make-bed --out ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_non${j}"
done
done


for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_split_ancestry_noneur_${i} | grep -E "completed"
done



for i in ${studies[@]}
do
echo ${i} && wc -l ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_*.fam
done

#########
# all_hce
# 22588 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/all_hce/all_hce_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 1175 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/all_hce/all_hce_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 23763 total
# niddk_old_gwas
# 2712 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 20 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 2732 total
# australia_omniexome
# 1245 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/australia_omniexome/australia_omniexome_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 54 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/australia_omniexome/australia_omniexome_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 1299 total
# gwas1
# 4652 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/gwas1/gwas1_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 7 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/gwas1/gwas1_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 4659 total
# gwas2
# 7758 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/gwas2/gwas2_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 15 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/gwas2/gwas2_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 7773 total
# pittsburgh_gsa
# 2701 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 8 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 2709 total
# spain_gsa
# 3396 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/spain_gsa/spain_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 12 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/spain_gsa/spain_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 3408 total
# italy_gsa
# 938 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/italy_gsa/italy_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 5 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/italy_gsa/italy_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 943 total
# kiel_austria_sibdcs_gsa
# 14086 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 243 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 14329 total
# netherlands_gsa
# 4315 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/netherlands_gsa/netherlands_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 222 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/netherlands_gsa/netherlands_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 4537 total
# slovenia_gsa
# 263 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/slovenia_gsa/slovenia_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# sweden_gsa
# 1364 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/sweden_gsa/sweden_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 32 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/sweden_gsa/sweden_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 1396 total
# niddk_broad_gsa
# 5130 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 296 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 5426 total
# niddk_feinstein_gsa
# 6965 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 1159 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 8124 total
# basque_gsa
# 1487 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/basque_gsa/basque_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 11 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/basque_gsa/basque_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 1498 total
# lithuania_gsa
# 2213 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/lithuania_gsa/lithuania_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 12 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/lithuania_gsa/lithuania_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 2225 total
# belgium_louis_gsa
# 1499 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 12 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 1511 total
# belgium_franchimont_gsa
# 1459 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 32 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 1491 total
# belgium_vermeire_gsa
# 3894 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 80 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 3974 total
# prism_nfe_gsa
# 426 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 36 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 462 total
# prism_nfe_gwas
# 717 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 93 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 810 total
# finland_illugwas
# 440 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/finland_illugwas/finland_illugwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 4 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/finland_illugwas/finland_illugwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 444 total
# german_affy6_old_gwas
# 2787 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 10 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 2797 total
# norway_affy6_old_gwas
# 542 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 2 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 544 total
# belgium_inf1_old_gwas
# 1396 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 10 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 1406 total
# belgium_inf2_old_gwas
# 270 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 1 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 271 total
# cedars_370k_old_gwas
# 509 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 88 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 597 total
# cedars_610k_old_gwas
# 825 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 56 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 881 total
# cedars_omni_old_gwas
# 1104 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 103 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 1207 total
# swedish_uc_old_gwas
# 1255 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 3 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 1258 total
# mccauley_gsa
# 758 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mccauley_gsa/mccauley_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 18 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mccauley_gsa/mccauley_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 776 total
# ccfa_gsa
# 1932 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/ccfa_gsa/ccfa_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 235 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/ccfa_gsa/ccfa_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 2167 total
# cedars_gsa
# 2699 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_gsa/cedars_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 360 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_gsa/cedars_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 3059 total
# bernstein_gsa
# 482 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/bernstein_gsa/bernstein_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 27 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/bernstein_gsa/bernstein_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 509 total
# farkkila_gsa
# 68 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/farkkila_gsa/farkkila_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# franchimont_gsa
# 2360 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/franchimont_gsa/franchimont_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 409 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/franchimont_gsa/franchimont_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 2769 total
# franke_gsa
# 853 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/franke_gsa/franke_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 11 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/franke_gsa/franke_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 864 total
# helmsley_prism_gsa
# 669 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 84 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 753 total
# helmsley_xavier_prism_gsa
# 1171 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 110 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 1281 total
# hyams_protect_gsa
# 323 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 93 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 416 total
# lewis_sparc_gsa
# 2529 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 317 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 2846 total
# mccauley_new_gsa
# 1323 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 287 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 1610 total
# mcgovern_gsa
# 5149 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 829 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 5978 total
# moayyedi_imagine_gsa
# 1011 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 109 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 1120 total
# newberry_share_gsa
# 747 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 111 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 858 total
# niddk_cho_gsa
# 1611 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 143 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 1754 total
# niddk_duerr_gsa
# 1868 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 69 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 1937 total
# niddk_rioux_gsa
# 871 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 42 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 913 total
# niddk_silverberg_gsa
# 2112 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 235 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 2347 total
# palotie_hus_gsa
# 853 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 16 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 869 total
# pekow_share_gsa
# 538 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 94 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 632 total
# rioux_igenomed_gsa
# 171 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 10 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 181 total
# sands_msccr_gsa
# 1120 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 299 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 1419 total
# stampfer_gsa
# 1439 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/stampfer_gsa/stampfer_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 25 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/stampfer_gsa/stampfer_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 1464 total
# vermeire_gsa
# 4540 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/vermeire_gsa/vermeire_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 119 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/vermeire_gsa/vermeire_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 4659 total
# weersma_gsa
# 682 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/weersma_gsa/weersma_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 23 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/weersma_gsa/weersma_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 705 total
# xavier_prism_gsa
# 619 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 67 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 686 total
# xavier_share_gsa
# 637 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 55 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 692 total
#########


# #####################################################
# # 21.2 SPLIT (Non Duplicates) - just to get exact N
# 
# ancestry=(eur amr afr sas eas)
# 
# for j in ${ancestry[@]}
# do
# for i in ${studies[@]}
# do
# bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
# -e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_split_ancestry_noneur_${i} \
# -o ${path_gwas}pre_imputation/QC/${i}/logs/stderr_split_ancestry_noneur_${i} \
# "/path/to/software/plink_linux_x86_64_20181202/./plink \
# --bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged \
# --keep-allele-order --allow-no-sex \
# --keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_${j}_ancestry_samples_noDuplicates \
# --make-bed --out ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${j}_nodup"
# done
# done
# 
# for i in ${studies[@]}
# do
# echo ${i} && wc -l ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_*_nodup.fam
# done

# rm after this



#####################################################
# 21.3 REMOVE MONOMORPHIC VARIANTS PER BATCH

### /software/R-4.3.1/bin/R

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

ancestry<-c("eur","noneur")

for (j in 1:length(cohorts)) {
  print(cohorts[j])
  for (i in 1:length(ancestry)) {
    file_freq<-paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_",ancestry[i],".frq",sep="")
    if(file.exists(file_freq)){
      frq<-read.table(file_freq,
                      head=T)
      mono<-frq[which(frq$MAF==0),]
      print(paste("N Monomorphic ",ancestry[i],": ",nrow(mono),sep=""))
      write.table(mono[,"SNP",drop=F],paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_list_variants_toexclude_",ancestry[i],"_monomorphic",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    }
    rm(mono,frq)
  }

}

#############
# [1] "australia_omniexome"
# [1] "N Monomorphic eur: 12010"
# [1] "N Monomorphic noneur: 65226"
# [1] "gwas1"
# [1] "N Monomorphic eur: 75"
# [1] "N Monomorphic noneur: 75357"
# [1] "gwas2"
# [1] "N Monomorphic eur: 0"
# [1] "N Monomorphic noneur: 62127"
# [1] "all_hce"
# [1] "N Monomorphic eur: 43412"
# [1] "N Monomorphic noneur: 16919"
# [1] "pittsburgh_gsa"
# [1] "N Monomorphic eur: 2796"
# [1] "N Monomorphic noneur: 129382"
# [1] "spain_gsa"
# [1] "N Monomorphic eur: 3"
# [1] "N Monomorphic noneur: 42667"
# [1] "italy_gsa"
# [1] "N Monomorphic eur: 1265"
# [1] "N Monomorphic noneur: 212293"
# [1] "kiel_austria_sibdcs_gsa"
# [1] "N Monomorphic eur: 12178"
# [1] "N Monomorphic noneur: 26634"
# [1] "netherlands_gsa"
# [1] "N Monomorphic eur: 20616"
# [1] "N Monomorphic noneur: 16932"
# [1] "slovenia_gsa"
# [1] "N Monomorphic eur: 1"
# [1] "sweden_gsa"
# [1] "N Monomorphic eur: 7192"
# [1] "N Monomorphic noneur: 77029"
# [1] "niddk_broad_gsa"
# [1] "N Monomorphic eur: 23878"
# [1] "N Monomorphic noneur: 19083"
# [1] "niddk_feinstein_gsa"
# [1] "N Monomorphic eur: 32312"
# [1] "N Monomorphic noneur: 13476"
# [1] "basque_gsa"
# [1] "N Monomorphic eur: 1052"
# [1] "N Monomorphic noneur: 155001"
# [1] "prism_nfe_gsa"
# [1] "N Monomorphic eur: 11260"
# [1] "N Monomorphic noneur: 75455"
# [1] "lithuania_gsa"
# [1] "N Monomorphic eur: 315"
# [1] "N Monomorphic noneur: 162306"
# [1] "belgium_louis_gsa"
# [1] "N Monomorphic eur: 1628"
# [1] "N Monomorphic noneur: 153169"
# [1] "belgium_franchimont_gsa"
# [1] "N Monomorphic eur: 2656"
# [1] "N Monomorphic noneur: 77818"
# [1] "belgium_vermeire_gsa"
# [1] "N Monomorphic eur: 49590"
# [1] "N Monomorphic noneur: 11623"
# [1] "prism_nfe_gwas"
# [1] "N Monomorphic eur: 1246"
# [1] "N Monomorphic noneur: 0"
# [1] "finland_illugwas"
# [1] "N Monomorphic eur: 1155"
# [1] "N Monomorphic noneur: 25650"
# [1] "german_affy6_old_gwas"
# [1] "N Monomorphic eur: 4930"
# [1] "N Monomorphic noneur: 83043"
# [1] "norway_affy6_old_gwas"
# [1] "N Monomorphic eur: 210"
# [1] "N Monomorphic noneur: 258566"
# [1] "belgium_inf1_old_gwas"
# [1] "N Monomorphic eur: 1"
# [1] "N Monomorphic noneur: 9442"
# [1] "belgium_inf2_old_gwas"
# [1] "N Monomorphic eur: 2"
# [1] "N Monomorphic noneur: 121400"
# [1] "niddk_old_gwas"
# [1] "N Monomorphic eur: 4"
# [1] "N Monomorphic noneur: 2142"
# [1] "cedars_370k_old_gwas"
# [1] "N Monomorphic eur: 1128"
# [1] "N Monomorphic noneur: 421"
# [1] "cedars_610k_old_gwas"
# [1] "N Monomorphic eur: 5024"
# [1] "N Monomorphic noneur: 1392"
# [1] "cedars_omni_old_gwas"
# [1] "N Monomorphic eur: 13002"
# [1] "N Monomorphic noneur: 1244"
# [1] "swedish_uc_old_gwas"
# [1] "N Monomorphic eur: 0"
# [1] "N Monomorphic noneur: 53065"
# [1] "mccauley_gsa"
# [1] "N Monomorphic eur: 1640"
# [1] "N Monomorphic noneur: 116573"
# [1] "ccfa_gsa"
# [1] "N Monomorphic eur: 25457"
# [1] "N Monomorphic noneur: 23662"
# [1] "cedars_gsa"
# [1] "N Monomorphic eur: 32495"
# [1] "N Monomorphic noneur: 16375"
# [1] "bernstein_gsa"
# [1] "N Monomorphic eur: 4781"
# [1] "N Monomorphic noneur: 82662"
# [1] "farkkila_gsa"
# [1] "N Monomorphic eur: 3"
# [1] "franchimont_gsa"
# [1] "N Monomorphic eur: 9699"
# [1] "N Monomorphic noneur: 22082"
# [1] "franke_gsa"
# [1] "N Monomorphic eur: 6191"
# [1] "N Monomorphic noneur: 157263"
# [1] "helmsley_prism_gsa"
# [1] "N Monomorphic eur: 1356"
# [1] "N Monomorphic noneur: 7"
# [1] "helmsley_xavier_prism_gsa"
# [1] "N Monomorphic eur: 648"
# [1] "N Monomorphic noneur: 0"
# [1] "hyams_protect_gsa"
# [1] "N Monomorphic eur: 18402"
# [1] "N Monomorphic noneur: 35544"
# [1] "lewis_sparc_gsa"
# [1] "N Monomorphic eur: 26593"
# [1] "N Monomorphic noneur: 20322"
# [1] "mccauley_new_gsa"
# [1] "N Monomorphic eur: 8643"
# [1] "N Monomorphic noneur: 18391"
# [1] "mcgovern_gsa"
# [1] "N Monomorphic eur: 33098"
# [1] "N Monomorphic noneur: 12981"
# [1] "moayyedi_imagine_gsa"
# [1] "N Monomorphic eur: 25193"
# [1] "N Monomorphic noneur: 32920"
# [1] "newberry_share_gsa"
# [1] "N Monomorphic eur: 19421"
# [1] "N Monomorphic noneur: 41073"
# [1] "niddk_cho_gsa"
# [1] "N Monomorphic eur: 21999"
# [1] "N Monomorphic noneur: 34450"
# [1] "niddk_duerr_gsa"
# [1] "N Monomorphic eur: 10262"
# [1] "N Monomorphic noneur: 63712"
# [1] "niddk_rioux_gsa"
# [1] "N Monomorphic eur: 12229"
# [1] "N Monomorphic noneur: 81781"
# [1] "niddk_silverberg_gsa"
# [1] "N Monomorphic eur: 23337"
# [1] "N Monomorphic noneur: 27226"
# [1] "palotie_hus_gsa"
# [1] "N Monomorphic eur: 2758"
# [1] "N Monomorphic noneur: 120236"
# [1] "pekow_share_gsa"
# [1] "N Monomorphic eur: 16976"
# [1] "N Monomorphic noneur: 49902"
# [1] "rioux_igenomed_gsa"
# [1] "N Monomorphic eur: 5246"
# [1] "N Monomorphic noneur: 147679"
# [1] "sands_msccr_gsa"
# [1] "N Monomorphic eur: 29532"
# [1] "N Monomorphic noneur: 13768"
# [1] "stampfer_gsa"
# [1] "N Monomorphic eur: 9216"
# [1] "N Monomorphic noneur: 116762"
# [1] "vermeire_gsa"
# [1] "N Monomorphic eur: 10194"
# [1] "N Monomorphic noneur: 45656"
# [1] "weersma_gsa"
# [1] "N Monomorphic eur: 10393"
# [1] "N Monomorphic noneur: 108512"
# [1] "xavier_prism_gsa"
# [1] "N Monomorphic eur: 15663"
# [1] "N Monomorphic noneur: 50676"
# [1] "xavier_share_gsa"
# [1] "N Monomorphic eur: 15581"
# [1] "N Monomorphic noneur: 62788"
#############


#####################################
# 21.1 SPLIT (INCLUDING Duplicates)


ancestry=(eur noneur)

for j in ${ancestry[@]}
do
for i in ${studies[@]}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_split_ancestry_${j}_2_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_split_ancestry_${j}_2_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${j} \
--keep-allele-order --allow-no-sex \
--exclude ${path_gwas}/pre_imputation/QC/${i}/${i}_list_variants_toexclude_${j}_monomorphic \
--freq \
--make-bed --out ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${j}_nomonom"
done
done

for j in ${ancestry[@]}
do
for i in ${studies[@]}
do
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_split_ancestry_${j}_2_${i}| grep -E "completed"
done
done

for i in ${studies[@]}
do
echo ${i} && wc -l ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_*_nomonom.fam \
&& wc -l ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_*_nomonom.bim
done


for i in ${studies[@]}
do
echo ${i} && wc -l ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_*_nomonom.fam \
&& wc -l ${path_gwas}/pre_imputation/QC/${i}/${i}_list_variants_toexclude_eur_monomorphic \
&& wc -l ${path_gwas}/pre_imputation/QC/${i}/${i}_list_variants_toexclude_noneur_monomorphic \
&& wc -l ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_*_nomonom.bim \
&& echo "      "
done


############
# all_hce
# 22588 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/all_hce/all_hce_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 1175 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/all_hce/all_hce_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 359209 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/all_hce/all_hce_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 385702 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/all_hce/all_hce_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim

# niddk_old_gwas
# 2712 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 20 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 296703 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 294565 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim

# australia_omniexome
# 1245 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/australia_omniexome/australia_omniexome_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 54 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/australia_omniexome/australia_omniexome_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 714842 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/australia_omniexome/australia_omniexome_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 661626 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/australia_omniexome/australia_omniexome_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim

# gwas1
# 4652 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/gwas1/gwas1_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 7 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/gwas1/gwas1_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 436931 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/gwas1/gwas1_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 361649 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/gwas1/gwas1_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# gwas2
# 7758 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/gwas2/gwas2_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 15 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/gwas2/gwas2_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
#
# 757166 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/gwas2/gwas2_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 695039 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/gwas2/gwas2_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# pittsburgh_gsa
# 2701 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 8 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 870503 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 743917 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# spain_gsa
# 3396 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/spain_gsa/spain_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 12 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/spain_gsa/spain_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 567250 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/spain_gsa/spain_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 524586 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/spain_gsa/spain_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# italy_gsa
# 938 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/italy_gsa/italy_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 5 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/italy_gsa/italy_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
#
# 501246 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/italy_gsa/italy_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 290218 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/italy_gsa/italy_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim

# kiel_austria_sibdcs_gsa
# 14086 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 243 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 494956 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 480500 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# netherlands_gsa
# 4315 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/netherlands_gsa/netherlands_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 222 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/netherlands_gsa/netherlands_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 499774 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/netherlands_gsa/netherlands_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 503458 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/netherlands_gsa/netherlands_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# slovenia_gsa
# 263 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/slovenia_gsa/slovenia_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 472855 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/slovenia_gsa/slovenia_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim

# sweden_gsa
# 1364 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/sweden_gsa/sweden_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 32 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/sweden_gsa/sweden_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
#
# 488408 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/sweden_gsa/sweden_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 418571 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/sweden_gsa/sweden_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# niddk_broad_gsa
# 5130 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 296 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 511878 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 516673 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# niddk_feinstein_gsa
# 6965 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 1159 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
#
# 516117 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 534953 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# basque_gsa
# 1487 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/basque_gsa/basque_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 11 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/basque_gsa/basque_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
#
# 502233 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/basque_gsa/basque_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 348284 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/basque_gsa/basque_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# lithuania_gsa
# 2213 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/lithuania_gsa/lithuania_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 12 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/lithuania_gsa/lithuania_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 496671 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/lithuania_gsa/lithuania_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 334680 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/lithuania_gsa/lithuania_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# belgium_louis_gsa
# 1499 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 12 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 496595 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 345054 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# belgium_franchimont_gsa
# 1459 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 32 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 492892 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 417730 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# belgium_vermeire_gsa
# 3894 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 80 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
#
# 509393 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 547360 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# prism_nfe_gsa
# 426 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 36 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 477019 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 412824 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# prism_nfe_gwas
# 717 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 93 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 238491 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 239737 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# finland_illugwas
# 440 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/finland_illugwas/finland_illugwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 4 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/finland_illugwas/finland_illugwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 230707 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/finland_illugwas/finland_illugwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 206212 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/finland_illugwas/finland_illugwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# german_affy6_old_gwas
# 2787 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 10 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 761510 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 683397 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# norway_affy6_old_gwas
# 542 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 2 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 707859 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 449503 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# belgium_inf1_old_gwas
# 1396 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 10 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
#
# 298228 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 288787 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# belgium_inf2_old_gwas
# 270 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 1 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
#
# 286778 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 165380 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# cedars_370k_old_gwas
# 509 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 88 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
#
# 332472 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 333179 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# cedars_610k_old_gwas
# 825 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 56 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 569781 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 573413 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim

# cedars_omni_old_gwas
# 1104 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 103 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 590026 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 601784 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# swedish_uc_old_gwas
# 1255 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 3 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 293205 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 240140 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# mccauley_gsa
# 758 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mccauley_gsa/mccauley_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 18 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mccauley_gsa/mccauley_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 493023 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mccauley_gsa/mccauley_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 378090 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mccauley_gsa/mccauley_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim

# ccfa_gsa
# 1932 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/ccfa_gsa/ccfa_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 235 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/ccfa_gsa/ccfa_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 503789 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/ccfa_gsa/ccfa_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 505584 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/ccfa_gsa/ccfa_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# cedars_gsa
# 2699 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_gsa/cedars_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 360 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_gsa/cedars_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 510760 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_gsa/cedars_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 526880 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_gsa/cedars_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# bernstein_gsa
# 482 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/bernstein_gsa/bernstein_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 27 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/bernstein_gsa/bernstein_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 418994 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/bernstein_gsa/bernstein_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 341113 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/bernstein_gsa/bernstein_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# farkkila_gsa
# 68 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/farkkila_gsa/farkkila_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 435687 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/farkkila_gsa/farkkila_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim

# franchimont_gsa
# 2360 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/franchimont_gsa/franchimont_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 409 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/franchimont_gsa/franchimont_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 501071 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/franchimont_gsa/franchimont_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 488688 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/franchimont_gsa/franchimont_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# franke_gsa
# 853 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/franke_gsa/franke_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 11 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/franke_gsa/franke_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 488896 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/franke_gsa/franke_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 337824 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/franke_gsa/franke_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# helmsley_prism_gsa
# 669 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 84 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 238576 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 239925 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# helmsley_xavier_prism_gsa
# 1171 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 110 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 239320 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 239968 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# hyams_protect_gsa
# 323 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 93 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 485724 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 468582 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# lewis_sparc_gsa
# 2529 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 317 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 508535 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 514806 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# mccauley_new_gsa
# 1323 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 287 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 502933 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 493185 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# mcgovern_gsa
# 5149 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 829 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 522701 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 542818 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# moayyedi_imagine_gsa
# 1011 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 109 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 493420 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 485693 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# newberry_share_gsa
# 747 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 111 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 488663 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 467011 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# niddk_cho_gsa
# 1611 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 143 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 503239 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 490788 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim

# niddk_duerr_gsa
# 1868 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 69 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 504360 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 450910 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim

# niddk_rioux_gsa
# 871 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 42 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 493739 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 424187 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# niddk_silverberg_gsa
# 2112 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 235 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 508459 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 504570 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# palotie_hus_gsa
# 853 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 16 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 485182 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 367704 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# pekow_share_gsa
# 538 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 94 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 487596 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 454670 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# rioux_igenomed_gsa
# 171 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 10 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam

# 475257 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 332824 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# sands_msccr_gsa
# 1120 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 299 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 494862 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 510626 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# stampfer_gsa
# 1439 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/stampfer_gsa/stampfer_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 25 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/stampfer_gsa/stampfer_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 500071 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/stampfer_gsa/stampfer_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 392525 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/stampfer_gsa/stampfer_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# vermeire_gsa
# 4540 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/vermeire_gsa/vermeire_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 119 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/vermeire_gsa/vermeire_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 510372 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/vermeire_gsa/vermeire_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 474910 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/vermeire_gsa/vermeire_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# weersma_gsa
# 682 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/weersma_gsa/weersma_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 23 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/weersma_gsa/weersma_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 483773 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/weersma_gsa/weersma_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 385654 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/weersma_gsa/weersma_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# xavier_prism_gsa
# 619 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 67 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 488695 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 453682 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim


# xavier_share_gsa
# 637 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 55 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 
# 490342 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 443135 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim

############