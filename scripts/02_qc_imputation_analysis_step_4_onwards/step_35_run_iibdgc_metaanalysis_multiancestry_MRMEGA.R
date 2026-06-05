# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# RUN ANALYSIS USING MR-MEGA:

cd /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/mrmega/
  
# SEE https://genomics.ut.ee/en/mr-mega

# Each GWA study file has mandatory column headers:
# 1) MARKERNAME – snp name
# 2) EA – effect allele
# 3) NEA – non effect allele
# 4) OR - odds ratio
# 5) OR_95L - lower confidence interval of OR
# 6) OR_95U - upper confidence interval of OR
# 7) EAF – effect allele frequency
# 8) N - sample size
# 9) CHROMOSOME  - chromosome of marker
# 10) POSITION - position of marker


# fist attempt using META-ANALYSIS (IIBDGC + DANISH + UKBB) with metal + MRMEGA with decode and finngen, but error saying not enough cohorts

##############################
# 1.- REFORMAT SUMMARY STATS #
##############################

################################################
# 1.1 META-ANALYSIS IIBDGC + DANISH + UKBB 

path_gwas=/path/to/ibdgwas/IIBDGC/

MEM=65000
pheno=(ibd cd uc)
ancestry=(eur)

MEM=25000
pheno=(ibd cd uc)
ancestry=(sas eas)

for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
bsub -J"for_MRMEGA" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/format_summary_stats_to_MRMEGA_${ph}__${j}_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/format_summary_stats_to_MRMEGA_${ph}__${j}_stdout \
"/software/R-4.3.1/bin/Rscript ~/git/IIBDGC_GWAS/scripts/other/format_summary_stats_files_for_MRMEGA.R ${ph} ${j} \
> ${path_gwas}scripts/logs/format_summary_stats_IIBDGC_danish_ukbb_to_MRMEGA_${ph}__${j}.Rout"
done
done

pheno=(ibd cd uc)
ancestry=(eur sas eas)
for j in ${ancestry[@]}
do
echo ${j} && for ph in ${pheno[@]}
do
echo ${ph} && tail -50 ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/format_summary_stats_to_MRMEGA_${ph}__${j}_stdout | grep "Successfully"
done
done

for chr in 22
do 


###################
# 2.- RUN MR-MEGA #
###################

# INPUT FILES FOR MR-MEGA
# For running MR-MEGA you have to create an input file (default name “mr-mega.in”), which contains the list of all study files. 
# The should have each results' file on separate row. 
# 
# Sample “MR-MEGA.in” file:
# Pop1.txt.gz
# Pop2.txt.gz
# Pop3.txt.gz
# Pop4.txt.gz
# Pop5.txt.gz
# Pop6.txt.gz
# Pop7.txt.gz
# Pop8.txt.gz

cd /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/mrmega/
  

pheno=(ibd cd uc)

for ph in ${pheno[@]}
do
ls ${path_gwas}post_imputation/2022/analysis/metaanalysis/mrmega/input_files/${ph}/allchr_*_MRMEGA_format.txt.gz > /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/mrmega/scripts/mrmega_eur_eas_sas_tier_2_${ph}_allchr.in
done


# pcs:
# --pc <int>
# This specifies the number od PC to use in regression. Default = 4. Please note that the PC count must be < cohort count - 2. 
# Therefore, if five cohorts have been used in the analyse, then the maximum number of PC-s can be two!

# we expect 4 groups (previous PCA - eur-non jewish, eur-jewish, sas and eas)

path_gwas=/path/to/ibdgwas/IIBDGC/

ph=cd

MEM=55000
for n_pc in {1..4}
do
bsub -J"for_MRMEGA" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/MRMEGA_test_pc${n_pc}_${ph}_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/MRMEGA_test_pc${n_pc}_${ph}_stdout \
"/path/to/software/MRMEGA/./MR-MEGA \
-i ${path_gwas}post_imputation/2022/analysis/metaanalysis/mrmega/scripts/mrmega_eur_eas_sas_tier_2_${ph}_allchr.in \
--pc ${n_pc} -o ${path_gwas}post_imputation/2022/analysis/metaanalysis/mrmega/output_files/${ph}/mrmega_eur_eas_sas_tier_2_${ph}_allchr_pc${n_pc}"
done

# see test_mrmeag_results.R, 2 pcs enough to discriminate into eur (eur fin) sas and eas
# based on results above, two pcs are enough to separate studies by broad continental ancestry, test now per phenotype

n_pc=1
MEM=55000
pheno=(ibd cd uc)

for ph in ${pheno[@]}
do
bsub -J"for_MRMEGA" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/MRMEGA_test_pc${n_pc}_${ph}_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/MRMEGA_test_pc${n_pc}_${ph}_stdout \
"/path/to/software/MRMEGA/./MR-MEGA \
-i ${path_gwas}post_imputation/2022/analysis/metaanalysis/mrmega/scripts/mrmega_eur_eas_sas_tier_2_${ph}_allchr.in \
--pc ${n_pc} -o ${path_gwas}post_imputation/2022/analysis/metaanalysis/mrmega/output_files/${ph}/mrmega_eur_eas_sas_tier_2_${ph}_allchr_pc${n_pc}"
done


pheno=(ibd cd uc)

for ph in ${pheno[@]}
do
tail -50 ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/MRMEGA_test_pc2_${ph}_stdout | grep "Successfully"
done
done


pheno=(ibd cd uc)

for ph in ${pheno[@]}
do
ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/mrmega/output_files/${ph}/mrmega_eur_eas_sas_tier_2_${ph}_allchr_pc${n_pc}*
done
done



#### compare results to metal:

pheno=(ibd cd uc)
path_gwas=/path/to/ibdgwas/IIBDGC/


MEM=45000

for ph in ${pheno[@]}
do
bsub -J"format" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 2 \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_compare_metal_mrmega_multiancestr_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_compare_metal_mrmega_multiancestr_stdout \
"/software/R-4.3.1/bin/Rscript ~/git/IIBDGC_GWAS/scripts/other/compare_metal_mrmega_multiancestry.R ${ph} \
> ${path_gwas}scripts/logs/compare_metal_mrmega_multiancestry_${ph}.Rout"
done


pheno=(ibd cd uc)

for ph in ${pheno[@]}
do
tail -50 ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_compare_metal_mrmega_multiancestr_stdout | grep "Successfully"
done



