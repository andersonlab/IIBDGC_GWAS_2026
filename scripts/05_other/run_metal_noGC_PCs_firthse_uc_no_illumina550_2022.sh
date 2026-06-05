# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# Author: https://orcid.org/0000-0002-7256-9752
#
# RUN META-ANALYSIS OF EUROPEAN SAMPLES WITH GENOTYPED DATA AVAILABLE

chr=$1

path_gwas=/path/to/ibdgwas/IIBDGC/
  
## CREATE INPUT FILE PER CHROMOSOME
# ALLELELABELS NON-EFFECT-ALLELE EFFECT-ALLELE 
  
echo "
MARKERLABEL ID
ALLELELABELS ALLELE0 ALLELE1
EFFECT BETA.Y1
STDERR SE.Y1

GENOMICCONTROL ON

SCHEME STDERR

PROCESS ${path_gwas}post_imputation/analysis/regenie/illumina370/uc/chr${chr}_illumina370_step2_uc_eur_sex_PCs_firthse.regenie
# PROCESS ${path_gwas}post_imputation/analysis/regenie/illumina550/uc/chr${chr}_illumina550_step2_uc_eur_sex_PCs_firthse.regenie
PROCESS ${path_gwas}post_imputation/analysis/regenie/humancoreexome/uc/chr${chr}_humancoreexome_step2_uc_eur_sex_PCs_firthse.regenie
PROCESS ${path_gwas}post_imputation/analysis/regenie/humanomni1/uc/chr${chr}_humanomni1_step2_uc_eur_sex_PCs_firthse.regenie
PROCESS ${path_gwas}post_imputation/analysis/regenie/quad610/uc/chr${chr}_quad610_step2_uc_eur_sex_PCs_firthse.regenie
PROCESS ${path_gwas}post_imputation/analysis/regenie/gsa/uc/chr${chr}_gsa_step2_uc_eur_sex_PCs_firthse.regenie
PROCESS ${path_gwas}post_imputation/analysis/regenie/illuminaexome/uc/chr${chr}_illuminaexome_step2_uc_eur_sex_PCs_firthse.regenie
PROCESS ${path_gwas}post_imputation/analysis/regenie/affymetrix6/uc/chr${chr}_affymetrix6_step2_uc_eur_sex_PCs_firthse.regenie
PROCESS ${path_gwas}post_imputation/analysis/regenie/humanomniexpress/uc/chr${chr}_humanomniexpress_step2_uc_eur_sex_PCs_firthse.regenie

OUTFILE ${path_gwas}/post_imputation/analysis/metaanalysis/uc/${chr}_meta_noGC_PCs_firthse_no_illumina550_ .txt

ANALYZE HETEROGENEITY

QUIT" > ${path_gwas}/post_imputation/analysis/metaanalysis/uc/chr${chr}_metal_input_noGC_PCs_firthse_no_illumina550.txt


# EXECUTE METAL
/path/to/software/metal/./metal ${path_gwas}/post_imputation/analysis/metaanalysis/uc/chr${chr}_metal_input_noGC_PCs_firthse_no_illumina550.txt

