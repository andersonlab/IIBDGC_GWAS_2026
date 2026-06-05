# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# Author: https://orcid.org/0000-0002-7256-9752
#
# RUN META-ANALYSIS OF EUROPEAN SAS and EAS SAMPLES WITH GENOTYPED DATA AVAILABLE

ph=$1

path_gwas=/path/to/ibdgwas/IIBDGC/
path=/path/to/ibdgwas_bioresource/

## CREATE INPUT FILE PER CHROMOSOME
# ALLELELABELS NON-EFFECT-ALLELE EFFECT-ALLELE 

echo "
MARKERLABEL ID
ALLELELABELS ALLELE0 ALLELE1
EFFECT BETA
STDERR SE

# GENOMICCONTROL ON

SCHEME STDERR

ADDFILTER A1FREQ >= 0.001
ADDFILTER A1FREQ <= 0.999
ADDFILTER INFO >= 0.4

CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

PROCESS ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/eas_iibdgc/${ph}/ibd_EAS_SiKJ_meta_${ph}.regenie.gz
PROCESS ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/allchr_interval_ibdbioresource_noneur_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz
PROCESS ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/allchr_ukbb_sas_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz
PROCESS ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/allchr_ukbb_eur_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz
PROCESS ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/finngen/${ph}/allchr_finngen_r10_K11_${ph}.regenie.gz
PROCESS ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/decode/${ph}/allchr_decode_${ph}_30062021_edited_noMult.regenie.gz
PROCESS ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/danish/${ph}/allchr_danish_gsa_${ph}_eur_sex_10PCs_saige_spa.regenie.gz
PROCESS ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/allchr_interval_ibdbioresource_eur_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/illumina370/${ph}/allchr_illumina370_eur_all_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/affymetrix6/${ph}/allchr_affymetrix6_eur_all_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/humanomniexpress/${ph}/allchr_humanomniexpress_eur_all_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/affymetrix500/${ph}/allchr_affymetrix500_eur_all_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/humancoreexome/${ph}/allchr_humancoreexome_eur_all_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/humanomni1/${ph}/allchr_humanomni1_eur_all_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/quad610/${ph}/allchr_quad610_eur_all_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/gsa/${ph}/allchr_gsa_eur_all_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/illuminaexome/${ph}/allchr_illuminaexome_eur_all_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz

OUTFILE ${path_gwas}/post_imputation/2022/analysis/metaanalysis/${ph}/allchr_eur_all_${ph}_meta_noGC_PCs_firthse_eur_eas_sas_tier_2_ .txt

ANALYZE HETEROGENEITY

QUIT" > ${path_gwas}/post_imputation/2022/analysis/metaanalysis/log/allchr_eur_all_${ph}_metal_input_noGC_PCs_firthse_eur_eas_sas_tier_2.txt

# EXECUTE METAL
/path/to/software/metal/./metal ${path_gwas}/post_imputation/2022/analysis/metaanalysis/log/allchr_eur_all_${ph}_metal_input_noGC_PCs_firthse_eur_eas_sas_tier_2.txt



