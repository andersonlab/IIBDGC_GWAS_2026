# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# TIER 2 EUR

# singularity exec iibdgc_postprocess_10_singularity.sif

# run FIXED EFFECTS meta-analysis with all datasets:

path_gwas="/path/to/ibdgwas/IIBDGC/"
pheno=(cd)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=30000

for ph in ${pheno[@]}
do
bsub -J"noGC_meta" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_eur_eas_sas_tier_2_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_eur_eas_sas_tier_2_stdout \
"sh ~/git/IIBDGC_GWAS/scripts/other/run_metal_noGC_PCs_firthse_eur_eas_sas_tier_2_2022.sh ${ph}"
done


for ph in ${pheno[@]}
do
echo ${ph} && tail -50 ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_eur_eas_sas_tier_2_stdout | grep "Successfully"
done



######################################
# reformat summary statistic results

# split file by chr:

pheno=(cd)

for ph in ${pheno[@]}
do 
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/allchr_eur_all_${ph}_meta_noGC_PCs_firthse_eur_eas_sas_tier_2_1.txt | head -1 > \
${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/header
done

j=eur_eas_sas
MEM=2000
for ph in ${pheno[@]}
do 
echo ${ph} && for chr in {1..22} X
do
bsub -J"noGC_meta" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/chr${chr}_metal_noGC_PCs_firthse_${ph}_${j}_tier_2_split_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/chr${chr}_metal_noGC_PCs_firthse_${ph}_${j}_tier_2_split_stdout \
"echo ${chr} && grep chr${chr}: ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/allchr_eur_all_${ph}_meta_noGC_PCs_firthse_eur_eas_sas_tier_2_1.txt \
| gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/chr${chr}_eur_eas_sas_${ph}_meta_noGC_PCs_firthse_1.txt.gz"
done
done

for ph in ${pheno[@]}
do 
echo ${ph} && for chr in {1..22} X
do echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/chr${chr}_metal_noGC_PCs_firthse_${ph}_${j}_tier_2_split_stdout | grep "Successfully"
done
done

for ph in ${pheno[@]}
do 
echo ${ph} && for chr in {1..22} X
do echo ${chr} && ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/chr${chr}_eur_eas_sas_${ph}_meta_noGC_PCs_firthse_1.txt.gz
done
done


###############################
# reformat the file

pheno=(ibd cd uc)
MEM=55000

for ph in ${pheno[@]}
do
for chr in {1..22} X
do
bsub -J"format" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${chr}_format_metal_eur_eas_sas_tier_2_withNeff_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${chr}_format_metal_eur_eas_sas_tier_2_withNeff_stdout \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/format_metal_files_eur_eas_sas_tier_2_withNeff.R ${ph} ${chr} \
> ${path_gwas}scripts/logs/format_metal_files_eur_eas_sas_tier_2_withNeff_${ph}_${chr}.Rout"
done
done

for ph in ${pheno[@]}
do
for chr in {1..22} X
do echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${chr}_format_metal_eur_eas_sas_tier_2_withNeff_stdout | grep "Successfully"
done
done

for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22} X
do echo ${chr} && cat ${path_gwas}scripts/logs/format_metal_files_eur_eas_sas_tier_2_withNeff_${ph}_${chr}.Rout
done
done

for ph in ${pheno[@]}
do
ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/*_${ph}_meta_eur_eas_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz | wc -l
done



# ## plot manhattan with sas summary stats:

MEM=6000

release="eur_eas_sas_tier_2"

bsub -J"plot_imp1" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -m "modern_hardware" -G ibdgwas -n 2 \
-e ${path_gwas}scripts/logs/draw_manhattan_plots_stderr_${release} \
-o ${path_gwas}scripts/logs/draw_manhattan_plots_stdout_${release} \
"/software/R-4.3.1/bin/Rscript ~/git/IIBDGC_GWAS/scripts/other/draw_manhattan_plots.R ${release} > \
${path_gwas}scripts/logs/draw_manhattan_plots_${release}.Rout"





