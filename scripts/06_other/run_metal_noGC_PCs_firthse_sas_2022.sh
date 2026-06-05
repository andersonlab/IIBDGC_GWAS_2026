# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# Author: https://orcid.org/0000-0002-7256-9752
#
# RUN META-ANALYSIS OF SOUTH ASIAN SAMPLES WITH GENOTYPED DATA AVAILABLE

ph=$1
j=$2

path_gwas=/path/to/ibdgwas/IIBDGC/

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
ADDFILTER INFO >= 0.400

CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

PROCESS ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/allchr_interval_ibdbioresource_noneur_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz
PROCESS ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/allchr_ukbb_sas_step2_${ph}_sas_sex_PCs_firthse_${ph}.regenie.gz

OUTFILE ${path_gwas}/post_imputation/2022/analysis/metaanalysis/${ph}/allchr_${j}_${ph}_meta_noGC_PCs_firthse_ .txt

ANALYZE HETEROGENEITY

QUIT" > ${path_gwas}/post_imputation/2022/analysis/metaanalysis/log/allchr_${j}_${ph}_metal_input_noGC_PCs_firthse.txt


# EXECUTE METAL
/path/to/software/metal/./metal ${path_gwas}/post_imputation/2022/analysis/metaanalysis/log/allchr_${j}_${ph}_metal_input_noGC_PCs_firthse.txt
