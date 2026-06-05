# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# Author: https://orcid.org/0000-0002-7256-9752
#
# RUN META-ANALYSIS OF EUROPEAN SAMPLES WITH GENOTYPED DATA AVAILABLE

ph=$1
j=$2

path_gwas=/path/to/ibdgwas/IIBDGC/
  
if [ $ph == "cd" ] ; then

## CREATE INPUT FILE PER CHROMOSOME
# ALLELELABELS NON-EFFECT-ALLELE EFFECT-ALLELE 
  
echo "
MARKERLABEL ID
ALLELELABELS ALLELE0 ALLELE1
EFFECT BETA
STDERR SE

GENOMICCONTROL ON

SCHEME STDERR

ADDFILTER A1FREQ >= 0.001
ADDFILTER A1FREQ <= 0.999

CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/illumina370/${ph}/allchr_illumina370_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/affymetrix6/${ph}/allchr_affymetrix6_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie
# PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/humanomniexpress/${ph}/allchr_humanomniexpress_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/affymetrix500/${ph}/allchr_affymetrix500_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/humancoreexome/${ph}/allchr_humancoreexome_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/humanomni1/${ph}/allchr_humanomni1_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/quad610/${ph}/allchr_quad610_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/gsa/${ph}/allchr_gsa_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/illuminaexome/${ph}/allchr_illuminaexome_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie


OUTFILE ${path_gwas}/post_imputation/2022/analysis/metaanalysis/${ph}/${chr}_${j}_${ph}_meta_withGC_PCs_firthse_ .txt

ANALYZE HETEROGENEITY

QUIT" > ${path_gwas}/post_imputation/2022/analysis/metaanalysis/log/allchr_${j}_${ph}_metal_input_withGC_PCs_firthse.txt

else

## CREATE INPUT FILE PER CHROMOSOME
# ALLELELABELS NON-EFFECT-ALLELE EFFECT-ALLELE 
  
echo "
MARKERLABEL ID
ALLELELABELS ALLELE0 ALLELE1
EFFECT BETA
STDERR SE

GENOMICCONTROL ON

SCHEME STDERR

ADDFILTER A1FREQ >= 0.001
ADDFILTER A1FREQ <= 0.999

CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/illumina370/${ph}/allchr_illumina370_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/affymetrix6/${ph}/allchr_affymetrix6_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/humanomniexpress/${ph}/allchr_humanomniexpress_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/affymetrix500/${ph}/allchr_affymetrix500_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/humancoreexome/${ph}/allchr_humancoreexome_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/humanomni1/${ph}/allchr_humanomni1_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/quad610/${ph}/allchr_quad610_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/gsa/${ph}/allchr_gsa_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie
PROCESS ${path_gwas}post_imputation/2022/analysis/regenie/illuminaexome/${ph}/allchr_illuminaexome_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}_subset_maf_0.001_info_0.4.regenie

OUTFILE ${path_gwas}/post_imputation/2022/analysis/metaanalysis/${ph}/${chr}_${j}_${ph}_meta_withGC_PCs_firthse_ .txt

ANALYZE HETEROGENEITY

QUIT" > ${path_gwas}/post_imputation/2022/analysis/metaanalysis/log/allchr_${j}_${ph}_metal_input_withGC_PCs_firthse.txt

fi

# EXECUTE METAL
/path/to/software/metal/./metal ${path_gwas}/post_imputation/2022/analysis/metaanalysis/log/allchr_${j}_${ph}_metal_input_withGC_PCs_firthse.txt


