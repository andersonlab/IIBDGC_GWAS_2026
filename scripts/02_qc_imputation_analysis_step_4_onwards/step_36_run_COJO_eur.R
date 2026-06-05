# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# double check we are not leaving out any region by using COJO to determine GWAS significant hits:

# give bfile as input

# --bfile test 
# Note: --bgen, --pfile, --bpfile, --mbfile, --mbgen, --mpfile and --mbpfile are currently only supported in the GRM and fastGWA analyses. 
# More functions will be available after rewriting some of the legacy codes. All the QC flags (e.g. --keep, --extract, --maf) in GCTA 1.92.4 are 
# currently applicable to these two formats.


##############################################################################################################################################
## CREATE .BGEN FILES FROM ALL ARRAYS USING SAME FILES (and N samples) USED FOR REGENIE STEP2, BUT CONVERTED AND SUBSET BY VARIANTS INCLUDED IN FINAL ANALYSIS:


# singularity exec iibdgc_postprocess_10_singularity.sif

path_gwas="/path/to/ibdgwas/IIBDGC/"

MEM=200
array=(illumina370 humanomniexpress affymetrix500 humanomni1 quad610 illuminaexome affymetrix6 humancoreexome gsa)
array=(illumina370 humanomniexpress affymetrix500 humanomni1 quad610 illuminaexome affymetrix6)
pheno=(cd uc ibd) 

for i in ${array[@]}
do
for ph in ${pheno[@]}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${i}_${chr}_ldsc_ldsc_plink_eur_tier2_stderr \
-o ${path_gwas}post_imputation/log/${i}_${chr}_ldsc_ldsc_plink_eur_tier2_stdout \
"awk -F '\t' '\$3 != "NA" { print }' ${path_gwas}post_imputation/2022/${i}/phenotype_data/${i}_all_studies_merged_eur_all_phenotype_${ph} \
| cut -f1-2 > ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_list_samples_included_in_${ph}_analysis"
done
done

for i in ${array[@]}
do
echo ${i} && for ph in ${pheno[@]}
do 
echo ${ph} && wc -l ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_list_samples_included_in_${ph}_analysis
done

# illumina370
# cd
# 4098 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/illumina370_list_samples_included_in_cd_analysis
# uc
# 4594 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/illumina370_list_samples_included_in_uc_analysis
# ibd
# 3782 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/illumina370_list_samples_included_in_ibd_analysis
# humanomniexpress
# cd
# 657 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/humanomniexpress_list_samples_included_in_cd_analysis
# uc
# 657 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/humanomniexpress_list_samples_included_in_uc_analysis
# ibd
# 657 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/humanomniexpress_list_samples_included_in_ibd_analysis
# affymetrix500
# cd
# 1767 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/affymetrix500_list_samples_included_in_cd_analysis
# uc
# 1767 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/affymetrix500_list_samples_included_in_uc_analysis
# ibd
# 1767 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/affymetrix500_list_samples_included_in_ibd_analysis
# humanomni1
# cd
# 1429 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/humanomni1_list_samples_included_in_cd_analysis
# uc
# 1429 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/humanomni1_list_samples_included_in_uc_analysis
# ibd
# 1429 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/humanomni1_list_samples_included_in_ibd_analysis
# quad610
# cd
# 1944 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/quad610_list_samples_included_in_cd_analysis
# uc
# 1944 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/quad610_list_samples_included_in_uc_analysis
# ibd
# 1944 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/quad610_list_samples_included_in_ibd_analysis
# illuminaexome
# cd
# 2677 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/illuminaexome_list_samples_included_in_cd_analysis
# uc
# 2677 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/illuminaexome_list_samples_included_in_uc_analysis
# ibd
# 2677 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/illuminaexome_list_samples_included_in_ibd_analysis
# affymetrix6
# cd
# 8326 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/affymetrix6_list_samples_included_in_cd_analysis
# uc
# 6445 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/affymetrix6_list_samples_included_in_uc_analysis
# ibd
# 6445 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/affymetrix6_list_samples_included_in_ibd_analysis
# humancoreexome
# cd
# 15803 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/humancoreexome_list_samples_included_in_cd_analysis
# uc
# 14867 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/humancoreexome_list_samples_included_in_uc_analysis
# ibd
# 20564 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/humancoreexome_list_samples_included_in_ibd_analysis
# gsa
# cd
# 42600 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/gsa_list_samples_included_in_cd_analysis
# uc
# 30578 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/gsa_list_samples_included_in_uc_analysis
# ibd
# 57999 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/gsa_list_samples_included_in_ibd_analysis


# 1.2- subset by variants with MAF 0.001 and info 0.4 as in regenie output

for i in ${array[@]}
do
for ph in ${pheno[@]}
do
for chr in {1..22} X
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${i}_${chr}_ldsc_subset_eur_tier2_stderr \
-o ${path_gwas}post_imputation/log/${i}_${chr}_ldsc_subset_eur_tier2_stdout \
"zcat ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/chr${chr}_${i}_eur_all_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz | awk -F ' ' '\$6 >= 0.001 && \$6 <= 0.999 && \$9 >= 0.4 { print \$3 }' \
> ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_chr${chr}_list_variants_included_in_${ph}_analysis"
done
done
done


for i in ${array[@]}
do
echo ${i} && for ph in ${pheno[@]}
do 
echo ${ph} && for chr in {1..22} X
do
echo ${chr} && wc -l  ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_chr${chr}_list_variants_included_in_${ph}_analysis
done
done

for i in ${array[@]}
do
echo ${i} && for ph in ${pheno[@]}
do 
ls -la ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_chr*_list_variants_included_in_${ph}_analysis | wc -l
done
done


MEM=20000
for ph in ${pheno[@]}
do 
for i in ${array[@]}
do
for chr in {1..22} X
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 4 \
-e ${path_gwas}post_imputation/log/${i}_${chr}_${ph}_ldsc_subset_eur_tier2_stderr \
-o ${path_gwas}post_imputation/log/${i}_${chr}_${ph}_ldsc_subset_eur_tier2_stdout \
"plink2 \
--bgen ${path_gwas}post_imputation/2022/${i}/imputed_data/${i}_chr${chr}.dose.subset.bgen 'ref-first' \
--sample ${path_gwas}post_imputation/2022/${i}/imputed_data/${i}_chr${chr}.dose.subset.sample \
--hard-call-threshold 0.2 \
--keep ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_list_samples_included_in_${ph}_analysis \
--extract ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_chr${chr}_list_variants_included_in_${ph}_analysis \
--threads 4 \
--make-bed --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_chr${chr}_subset_included_in_${ph}_analysis"
done
done
done

for i in ${array[@]}
do
echo ${i} && for chr in {1..22} X
do
echo ${chr} && tail -30 ${path_gwas}post_imputation/log/${i}_${chr}_${ph}_ldsc_subset_eur_tier2_stdout  | grep -E "Successfully|Exited"
done
done

rm ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/chr*_subset_included_in_*_analysis_mergelist.txt

array=(illumina370 humanomniexpress affymetrix500 humanomni1 quad610 illuminaexome affymetrix6 humancoreexome)

# for ph in ibd
for ph in uc
do
for chr in {1..22}
do
for i in ${array[@]}
do
echo ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_chr${chr}_subset_included_in_${ph}_analysis >> \
${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/chr${chr}_subset_included_in_${ph}_analysis_mergelist.txt
done
done
done


array=(illumina370 affymetrix500 humanomni1 quad610 illuminaexome affymetrix6 humancoreexome)

for ph in cd
do
for chr in {1..22}
do
for i in ${array[@]}
do
echo ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_chr${chr}_subset_included_in_${ph}_analysis >> \
${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/chr${chr}_subset_included_in_${ph}_analysis_mergelist.txt
done
done
done

array=(humanomniexpress quad610 affymetrix6)
array=(quad610 affymetrix6)

# for ph in ibd uc
for ph in cd
do
for chr in X
do
for i in ${array[@]}
do
echo ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_chr${chr}_subset_included_in_${ph}_analysis >> \
${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/chr${chr}_subset_included_in_${ph}_analysis_mergelist.txt
done
done
done


for chr in {1..22} X
do
ls -la ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/chr${chr}_subset_included_in_${ph}_analysis_mergelist.txt
done

MEM=80000
for ph in ibd cd uc
do
for chr in {1..22} X
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 4 \
-e ${path_gwas}post_imputation/log/bmerge_${chr}_${ph}_cojo_eur_tier2_stderr \
-o ${path_gwas}post_imputation/log/bmerge_${chr}_${ph}_cojo_eur_tier2_stdout \
"/path/to/software/username/plink_linux_x86_64_20231211/./plink \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/gsa_chr${chr}_subset_included_in_${ph}_analysis \
--merge-list ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/chr${chr}_subset_included_in_${ph}_analysis_mergelist.txt \
--keep-allele-order --threads 4 \
--make-bed --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_${ph}_analysis"
done
done


for ph in ${pheno[@]}
do
for chr in {1..22} X
do
echo ${chr} && tail -30 ${path_gwas}post_imputation/log/bmerge_${chr}_${ph}_cojo_eur_tier2_stdout  | grep -E 'Successfully|Exite'
# echo ${chr} && tail -30 ${path_gwas}post_imputation/log/bmerge_${chr}_${ph}_cojo_eur_tier2_stdout  | grep -E 'Max Memory'
done
done



##### RUN COJO NON-INFORMED - FIND ANY CONDITIONALLY SIGNIF VARIANT AT 5e-8 - include Neff:

# subset the data by Neff>=0.5 and hetPval >=1E-15 

echo "SNP A1 A2 freq b se p N" \
> ${path_gwas}post_imputation/2022/analysis/conditional_analysis/header.ma

MEM=800
path_gwas=/path/to/ibdgwas/IIBDGC/
pheno=(ibd cd uc)

for chr in {1..22} X
do
for ph in ${pheno[@]}
do 
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 4 \
-e ${path_gwas}post_imputation/log/${chr}_${ph}_cojo_sumstats_eur_tier2_stderr \
-o ${path_gwas}post_imputation/log/${chr}_${ph}_cojo_sumstats_eur_tier2_stdout \
"zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/${chr}_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz \
 | awk -F ' ' '\$12 >= 1E-15 { print }' | awk -F ' ' '\$23 >= 0.5 { print }' | awk -F '\t' '{ print \$1,\$4,\$3,\$16,\$5,\$6,\$7,\$15}' > \
${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${chr}_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_Neff_cojo_format.ma"
done
done

for chr in {1..22} X
do
for ph in ${pheno[@]}
do 
ls -la ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${chr}_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_Neff_cojo_format.ma
done
done


for ph in ${pheno[@]}
do cat ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/*_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_Neff_cojo_format.ma \
| grep -v "MarkerName" \
> ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_Neff_cojo_format_tmp.ma
done

for ph in ${pheno[@]}
do 
cat ${path_gwas}post_imputation/2022/analysis/conditional_analysis/header.ma  ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_Neff_cojo_format_tmp.ma \
> ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_Neff_cojo_format.ma
done

for ph in ${pheno[@]}
do 
rm ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_Neff_cojo_format_tmp.ma
done

# for ph in ${pheno[@]}
# do 
# sed -i 's/chrX:/chr23:/g' ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_Neff_cojo_format.ma
# done

for ph in ${pheno[@]}
do echo ${ph} && wc -l ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_Neff_cojo_format.ma
done

# ibd
# 13656989 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/allchr_ibd_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_Neff_cojo_format.ma
# cd
# 13909224 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/allchr_cd_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_Neff_cojo_format.ma
# uc
# 13694540 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/allchr_uc_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_Neff_cojo_format.ma


path_gwas=/path/to/ibdgwas/IIBDGC/
pheno=(ibd cd uc)
MEM=55000

for ph in ${pheno[@]}
do
for chr in {1..22} 
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q long \
-o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_no_conditioning_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_no_conditioning_stderr \
"/path/to/software/username/gcta-1.94.1-linux-kernel-3-x86_64/./gcta64 \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_${ph}_analysis \
--chr ${chr} \
--cojo-file ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_Neff_cojo_format.ma \
--cojo-slct \
--cojo-p 5e-8 \
--out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${chr}_${ph}_independent_index_variants_eur_tier_2_ld_allarrays_no_conditioning"
done
done




for ph in ${pheno[@]}
do
for chr in X
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q long \
-o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_no_conditioning_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_no_conditioning_stderr \
"/path/to/software/username/gcta-1.94.1-linux-kernel-3-x86_64/./gcta64 \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_${ph}_analysis \
--chr 23 \
--cojo-file ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_Neff_cojo_format.ma \
--cojo-slct \
--cojo-p 5e-8 \
--out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${chr}_${ph}_independent_index_variants_eur_tier_2_ld_allarrays_no_conditioning"
done
done




for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22} X
do 
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_no_conditioning_stdout | grep -E "Successfully|Exited"
done
done


for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22}
do 
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_no_conditioning_stdout | grep -E "Max Memory"
done
done

###############

# run conditional on known but now with all samples as LD ref:

path_gwas=/path/to/ibdgwas/IIBDGC/
pheno=(ibd cd uc)
# MEM=35000
MEM=55000

for ph in ${pheno[@]}
do
# for chr in {1..22}
for chr in {1..4}
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G humgen-priority \
-o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_conditioning_known_ibd_ld_allarrays_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_conditioning_known_ibd_ld_allarrays_stderr \
"/path/to/software/username/gcta-1.94.1-linux-kernel-3-x86_64/./gcta64 \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_${ph}_analysis \
--chr ${chr} \
--cojo-file ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format.ma \
--cojo-cond ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/cond_chr${chr}.snplist \
--out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${chr}_${ph}_independent_index_variants_eur_tier_2_ld_allarrays_conditioning_on_known_ibd_index"
done



for ph in ${pheno[@]}
do
for chr in X
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G humgen-priority \
-o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_conditioning_known_ibd_ld_allarrays_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_conditioning_known_ibd_ld_allarrays_stderr \
"/path/to/software/username/gcta-1.94.1-linux-kernel-3-x86_64/./gcta64 \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_${ph}_analysis \
--chr 23 \
--cojo-file ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format.ma \
--cojo-cond ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/cond_chr${chr}.snplist \
--out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${chr}_${ph}_independent_index_variants_eur_tier_2_ld_allarrays_conditioning_on_known_ibd_index"
done
done



for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22} X
do 
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_conditioning_known_ibd_ld_allarrays_stdout | grep -E "Successfully|Exited"
done
done
