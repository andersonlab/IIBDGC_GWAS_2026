# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# before meta-analysis make sure to keep only variants per array that have info >0.4 and maf>0.001 in the analysis - controlled by metal

###################################################################################################

# combine all files from ibdbioresource

path_gwas="/path/to/ibdgwas/IIBDGC/"
pheno=(ibd cd uc)

for ph in ${pheno[@]}
do 
zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/5_interval_ibdbioresource_noneur_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz | head -1  > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/header_regenie_sas
zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/*_interval_ibdbioresource_noneur_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz | grep -v "CHROM" > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/tmp_sas.regenie
cat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/header_regenie_sas ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/tmp_sas.regenie | gzip \
> ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/allchr_interval_ibdbioresource_noneur_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz
rm ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/header_regenie_sas 
rm ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/tmp_sas.regenie
done

for ph in ${pheno[@]}
do 
echo ${ph} && zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/allchr_interval_ibdbioresource_noneur_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz | wc -l
done

# ibd
# 39660073
# cd
# 32797941
# uc
# 35512511


for ph in ${pheno[@]}
do 
ls -la ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/*_interval_ibdbioresource_noneur_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz
done

for ph in ${pheno[@]}
do 
ls -la ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/allchr_interval_ibdbioresource_noneur_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz
done



###############################
# combine all files from ukbb

path_gwas="/path/to/ibdgwas/IIBDGC/"
pheno=(ibd cd uc)

for ph in ${pheno[@]}
do 
zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/chr1_ukbb_sas_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz | head -1  > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/header_regenie
zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/chr*_ukbb_sas_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz | grep -v "CHROM" > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/tmp.regenie
cat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/header_regenie ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/tmp.regenie | sed 's/chr23:/chrX:/g' | gzip \
> ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/allchr_ukbb_sas_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz
rm ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/header_regenie 
rm ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/tmp.regenie
done


for ph in ${pheno[@]}
do 
echo ${ph} && zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/allchr_ukbb_sas_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz | wc -l
done

# ibd
# 33018883
# cd
# 32774007
# uc
# 32912932


for ph in ${pheno[@]}
do 
ls -la ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/allchr_ukbb_sas_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz 
done


for ph in ${pheno[@]}
do 
ls -la ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/chr*_ukbb_sas_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz
done


# run meta analysis

path_gwas="/path/to/ibdgwas/IIBDGC/"
pheno=(ibd cd uc)
ancestry=(sas)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=15000

for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
bsub -J"noGC_meta" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_stdout \
"sh ~/git/IIBDGC_GWAS/scripts/other/run_metal_noGC_PCs_firthse_sas_2022.sh ${ph} ${j}"
done
done


for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
echo ${ph} && \
ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_stdout && \
tail -50 ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_stdout | grep "Successfully"
done
done



######################################
# reformat summary statistic results

# split file by chr:

pheno=(ibd cd uc)

for ph in ${pheno[@]}
do 
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/allchr_sas_${ph}_meta_noGC_PCs_firthse_1.txt | head -1 > \
${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/header
done


for ph in ${pheno[@]}
do 
echo ${ph} && for chr in {1..22} X
do
echo ${chr} && grep chr${chr}: ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/allchr_sas_${ph}_meta_noGC_PCs_firthse_1.txt \
| gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/chr${chr}_sas_${ph}_meta_noGC_PCs_firthse_1.txt.gz
done
done


zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/allchr_interval_ibdbioresource_noneur_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz | head -1 && \
zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/allchr_interval_ibdbioresource_noneur_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz | grep chr22:37658671:G:GC
CHROM   GENPOS  ID      ALLELE0 ALLELE1 A1FREQ  A1FREQ_CASES    A1FREQ_CONTROLS INFO    N       N_CASES N_CONTROLS      TEST    BETA    SE      CHISQ   LOG10P  EXTRA
22      37658671        chr22:37658671:G:GC     G       GC      0.0141485       0.0157116       0.0121671       0.4     1873    1047    826     ADD     0.9607130.475595 4.08049 1.3627  NA

zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/allchr_ukbb_sas_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz | head -1 \
&& zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/allchr_ukbb_sas_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz | grep chr22:37658671:G:GC
CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ A1FREQ_CASES A1FREQ_CONTROLS INFO N N_CASES N_CONTROLS TEST BETA SE CHISQ LOG10P EXTRA
22 37658671 chr22:37658671:G:GC G GC 0.0194737 0.0166789 0.0195464 0.409643 6326 160 6166 ADD -0.243851 0.698617 0.121835 0.138434 NA

zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/chr${chr}_sas_${ph}_meta_noGC_PCs_firthse_1.txt.gz | head -1 && \
${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/chr${chr}_sas_${ph}_meta_noGC_PCs_firthse_1.txt.gz | grep chr22:37658671:G:GC

###############################
# reformat the file

pheno=(ibd cd uc)
MEM=20000

for ph in ${pheno[@]}
do
for chr in {1..22} X
do
bsub -J"format" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${chr}_format_metal_sas_tier_2_withNeff_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${chr}_format_metal_sas_tier_2_withNeff_stdout \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/format_metal_files_sas_tier_2_withNeff.R ${ph} ${chr} \
> ${path_gwas}scripts/logs/format_metal_files_sas_tier_2_withNeff_${ph}_${chr}.Rout"
done
done


for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22} X
do
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${chr}_format_metal_sas_tier_2_withNeff_stdout | grep "Successfully"
done
done


for ph in ${pheno[@]}
do
for chr in {1..22} X
do echo ${chr} \
&& zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/${chr}_${ph}_meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz | wc -l \
&& zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/${chr}_${ph}_meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz | wc -l
done
done



for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22} X
do echo ${chr} && cat ${path_gwas}scripts/logs/format_metal_files_sas_tier_2_${ph}_${chr}.Rout
done
done


## plot manhattan with sas summary stats:

MEM=6000

release="sas_tier_2"

bsub -J"plot_imp1" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -m "modern_hardware" -G ibdgwas -n 2 \
-e ${path_gwas}scripts/logs/draw_manhattan_plots_stderr_${release} \
-o ${path_gwas}scripts/logs/draw_manhattan_plots_stdout_${release} \
"/software/R-4.3.1/bin/Rscript ~/git/IIBDGC_GWAS/scripts/other/draw_manhattan_plots.R ${release} > \
${path_gwas}scripts/logs/draw_manhattan_plots_${release}.Rout"