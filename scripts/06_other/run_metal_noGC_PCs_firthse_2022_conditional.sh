# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# Author: https://orcid.org/0000-0002-7256-9752
#
# RUN META-ANALYSIS OF EUROPEAN SAMPLES WITH GENOTYPED DATA AVAILABLE

ph=$1
j=$2
signal=$3
chr=$4

path_gwas=/path/to/ibdgwas/IIBDGC/
  
if [ $ph == "cd" ] ; then

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
ADDFILTER INFO > 0.4

CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/illumina370/${ph}/chr${chr}_illumina370_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/affymetrix6/${ph}/chr${chr}_affymetrix6_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie
# PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/humanomniexpress/${ph}/chr${chr}_humanomniexpress_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/affymetrix500/${ph}/chr${chr}_affymetrix500_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/humancoreexome/${ph}/chr${chr}_humancoreexome_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/humanomni1/${ph}/chr${chr}_humanomni1_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/quad610/${ph}/chr${chr}_quad610_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/gsa/${ph}/chr${chr}_gsa_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/illuminaexome/${ph}/chr${chr}_illuminaexome_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie


OUTFILE ${path_gwas}/post_imputation/2022/analysis/metaanalysis/${ph}/chr${chr}_${j}_${ph}_meta_noGC_PCs_${signal}_firthse_noGC_ .txt

ANALYZE HETEROGENEITY

QUIT" > ${path_gwas}/post_imputation/2022/analysis/metaanalysis/log/chr${chr}_${j}_${ph}_metal_input_noGC_PCs_${signal}_firthse_noGC.txt

else

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
ADDFILTER INFO > 0.4

CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/illumina370/${ph}/chr${chr}_illumina370_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/affymetrix6/${ph}/chr${chr}_affymetrix6_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/humanomniexpress/${ph}/chr${chr}_humanomniexpress_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/affymetrix500/${ph}/chr${chr}_affymetrix500_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/humancoreexome/${ph}/chr${chr}_humancoreexome_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/humanomni1/${ph}/chr${chr}_humanomni1_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/quad610/${ph}/chr${chr}_quad610_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/gsa/${ph}/chr${chr}_gsa_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/illuminaexome/${ph}/chr${chr}_illuminaexome_${j}_step2_${ph}_eur_sex_PCs_${signal}_firthse_${ph}.regenie

OUTFILE ${path_gwas}post_imputation/2022/analysis/conditional_analysis_post_fm/regenie/metaanalysis/${ph}/chr${chr}_${j}_${ph}_meta_noGC_PCs_${signal}_firthse_noGC_ .txt

ANALYZE HETEROGENEITY

QUIT" > ${path_gwas}/post_imputation/2022/analysis/metaanalysis/log/chr${chr}_${j}_${ph}_metal_input_noGC_PCs_${signal}_firthse_noGC.txt

fi

# EXECUTE METAL
/path/to/software/username/metal/./metal ${path_gwas}/post_imputation/2022/analysis/metaanalysis/log/chr${chr}_${j}_${ph}_metal_input_noGC_PCs_${signal}_firthse_noGC.txt


