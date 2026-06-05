# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
######################
# 22.- CREATE  FILES #
######################

###########################################
# 22.1 Convert ped/map files to VCF files

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

for i in ${studies[@]}
do 
mkdir ${path_gwas}imputation_ready/${i}/
done
  
ancestry=(eur noneur)

for chr in {1..23}
do
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_${j}_${i}_chr${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_${j}_${i}_chr${chr} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${j}_nomonom \
--allow-no-sex \
--keep-allele-order \
--chr ${chr} \
--output-chr chr26 --recode vcf-iid --out ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022"
done
done
done
  
for i in ${studies[@]}
do 
echo ${i} && for j in ${ancestry[@]}
do 
ls -la ${path_gwas}imputation_ready/${i}/${i}_${j}_chr22_2022.vcf
done
done

###########################################
# 22.2 Create a sorted vcf.gz

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

ancestry=(eur noneur)

for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
for chr in 1
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_2_${j}_${i}_chr${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_2_${j}_${i}_chr${chr} \
"/path/to/software/bcftools-1.16/./bcftools sort ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022.vcf \
-Oz -o ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022.vcf.gz"
done
done
done

# rm intermediate files:

for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/${i}_*_chr*.vcf
done
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/*.hh
done
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/*.log
done
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/*.nosex
done

  