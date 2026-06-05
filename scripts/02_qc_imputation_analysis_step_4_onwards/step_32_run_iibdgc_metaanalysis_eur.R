# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# before meta-analysis make sure to keep only variants per array that have info >0.4 and maf>0.001 in the analysis - controlled by metal
# singularity exec iibdgc_postprocess_10_singularity.sif

path_gwas=/path/to/ibdgwas/IIBDGC/
  
MEM=7000
bsub -J"plot_imp1" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -m "modern_hardware" -G ibdgwas -n 2 -q long \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/plot_number_variants_post_qc_old_new_merging_procedure_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/plot_number_variants_post_qc_old_new_merging_procedure_stdout \
"/software/R-4.3.1/bin/Rscript ${path_gwas}scripts/plot_number_variants_post_qc_old_new_merging_procedure.R > \
${path_gwas}scripts/logs/plot_number_variants_post_qc_old_new_merging_procedure.Rout"


######################################################################################################################################################

# plot final INFO distribution SNPs included in meta analysis


MEM=6000
bsub -J"plot_imp2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -m "modern_hardware" -G ibdgwas -n 2 \
-e ${path_gwas}post_imputation/analysis/metaanalysis/log/plot_info_rsqr_distribution_variants_stderr \
-o ${path_gwas}post_imputation/analysis/metaanalysis/log/plot_info_rsqr_distribution_variants_stdout \
"/software/R-4.3.1/bin/Rscript ${path_gwas}scripts/plot_info_rsqr_distribution_variants.R > \
${path_gwas}scripts/logs/plot_info_rsqr_distribution_variants.Rout"


~/Desktop/


####################################################################################################################################################################################
####################################################################################################################################################################################

###############
### with GC ###
###############

# pheno=(ibd cd uc)
# array=(illumina370 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)
# ancestry=(eur_nonjewish eur_all)


# files need to be combined first:

path_gwas=/path/to/ibdgwas/IIBDGC/
  
pheno=(ibd cd uc)
ancestry=(eur_all)
array=(illumina370 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)
array=(humancoreexome)

for j in ${ancestry[@]}
do
for i in ${array[@]}
do
for ph in ${pheno[@]}
do 
zcat ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/chr1_${i}_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz | head -1 > ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/header_regenie
zcat ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/chr*_${i}_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz | grep -v "CHROM" > ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/tmp.regenie
cat ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/header_regenie ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/tmp.regenie | gzip \
> ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/allchr_${i}_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz
rm ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/header_regenie 
rm ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/tmp.regenie
done
done
done

j=eur_all

for i in ${array[@]}
do
echo ${i} && for ph in ${pheno[@]}
do 
echo ${ph} && zcat ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/allchr_${i}_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz | wc -l
done
done


##############
# illumina370
# ibd
# 34513504
# cd
# 29173382
# uc
# 28142683
# affymetrix6
# ibd
# 55343003
# cd
# 34311790
# uc
# 55145212
# humanomniexpress
# ibd
# 25460379
# cd
# 0
# uc
# 25460379
# affymetrix500
# ibd
# 41157933
# cd
# 40961945
# uc
# 33274217
# humancoreexome
# ibd
# 81772404
# cd
# 76214736
# uc
# 71200106
# humanomni1
# ibd
# 32896125
# cd
# 30123115
# uc
# 29812510
# quad610
# ibd
# 41752845
# cd
# 37630965
# uc
# 35333915
# gsa
# ibd
# 134478741
# cd
# 121354779
# uc
# 106635958
# illuminaexome
# ibd
# 25400370
# cd
# 21747087
# uc
# 20233675

for j in ${ancestry[@]}
do
for i in ${array[@]}
do
for ph in ${pheno[@]}
do 
zcat ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/allchr_${i}_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz | tail -2 && \
zcat ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/chrX_${i}_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz | tail -2
done
done
done

for j in ${ancestry[@]}
do
echo ${j} && for i in ${array[@]}
do
echo ${i} && for ph in ${pheno[@]}
do 
echo ${ph} && 
zcat ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/allchr_${i}_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz | head -1 | wc -w && \
zcat ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/allchr_${i}_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz | tail -1 | wc -w
done
done
done


for j in ${ancestry[@]}
do
for i in ${array[@]}
do
echo ${i} && for ph in ${pheno[@]}
do 
echo ${ph} && for chr in {1..22} X
do 
echo ${chr} && \
zcat ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/chr${chr}_${i}_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz | head -1 | wc -w && \
zcat ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/chr${chr}_${i}_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz | tail -1 | wc -w
done
done
done
done

  
  
path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=30000

pheno=(ibd cd uc)
ancestry=(eur_all)

for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
bsub -J"ibd" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_withGC_PCs_firthse_${ph}_${j}_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_withGC_PCs_firthse_${ph}_${j}_stdout \
"sh ${path_gwas}scripts/run_metal_withGC_PCs_firthse_2022.sh ${ph} ${j}"
done
done


for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do 
echo ${ph} && less ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_withGC_PCs_firthse_${ph}_${j}_stdout | \
grep -E 'Processing file|Genomic control parameter'
done
done
 
#############################

# checkout lambda per study:

# ibd
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/illumina370/ibd/allchr_illumina370_eur_all_step2_ibd_eur_sex_PCs_firthse_ibd.regenie'
# ## Genomic control parameter is 1.083, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/affymetrix6/ibd/allchr_affymetrix6_eur_all_step2_ibd_eur_sex_PCs_firthse_ibd.regenie'
# ## Genomic control parameter is 1.083, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/humanomniexpress/ibd/allchr_humanomniexpress_eur_all_step2_ibd_eur_sex_PCs_firthse_ibd.regenie'
# ## Genomic control parameter is 1.086, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/affymetrix500/ibd/allchr_affymetrix500_eur_all_step2_ibd_eur_sex_PCs_firthse_ibd.regenie'
# ## Genomic control parameter is 1.051, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/humancoreexome/ibd/allchr_humancoreexome_eur_all_step2_ibd_eur_sex_PCs_firthse_ibd.regenie'
# ## Genomic control parameter is 1.109, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/humanomni1/ibd/allchr_humanomni1_eur_all_step2_ibd_eur_sex_PCs_firthse_ibd.regenie'
# ## Genomic control parameter is 1.031, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/quad610/ibd/allchr_quad610_eur_all_step2_ibd_eur_sex_PCs_firthse_ibd.regenie'
# ## Genomic control parameter is 1.063, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/gsa/ibd/allchr_gsa_eur_all_step2_ibd_eur_sex_PCs_firthse_ibd.regenie'
# ## Genomic control parameter is 1.104, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/illuminaexome/ibd/allchr_illuminaexome_eur_all_step2_ibd_eur_sex_PCs_firthse_ibd.regenie'
# ## Genomic control parameter is 1.146, adjusting test statistics
# 
# cd
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/illumina370/cd/allchr_illumina370_eur_all_step2_cd_eur_sex_PCs_firthse_cd.regenie'
# ## Genomic control parameter is 1.051, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/affymetrix6/cd/allchr_affymetrix6_eur_all_step2_cd_eur_sex_PCs_firthse_cd.regenie'
# ## Genomic control parameter is 0.746, no adjustment made
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/humanomniexpress/cd/allchr_humanomniexpress_eur_all_step2_cd_eur_sex_PCs_firthse_cd.regenie'
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/affymetrix500/cd/allchr_affymetrix500_eur_all_step2_cd_eur_sex_PCs_firthse_cd.regenie'
# ## Genomic control parameter is 1.051, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/humancoreexome/cd/allchr_humancoreexome_eur_all_step2_cd_eur_sex_PCs_firthse_cd.regenie'
# ## Genomic control parameter is 1.097, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/humanomni1/cd/allchr_humanomni1_eur_all_step2_cd_eur_sex_PCs_firthse_cd.regenie'
# ## Genomic control parameter is 1.016, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/quad610/cd/allchr_quad610_eur_all_step2_cd_eur_sex_PCs_firthse_cd.regenie'
# ## Genomic control parameter is 1.059, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/gsa/cd/allchr_gsa_eur_all_step2_cd_eur_sex_PCs_firthse_cd.regenie'
# ## Genomic control parameter is 1.079, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/illuminaexome/cd/allchr_illuminaexome_eur_all_step2_cd_eur_sex_PCs_firthse_cd.regenie'
# ## Genomic control parameter is 1.047, adjusting test statistics
# 
# uc
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/illumina370/uc/allchr_illumina370_eur_all_step2_uc_eur_sex_PCs_firthse_uc.regenie'
# ## Genomic control parameter is 0.770, no adjustment made
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/affymetrix6/uc/allchr_affymetrix6_eur_all_step2_uc_eur_sex_PCs_firthse_uc.regenie'
# ## Genomic control parameter is 1.082, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/humanomniexpress/uc/allchr_humanomniexpress_eur_all_step2_uc_eur_sex_PCs_firthse_uc.regenie'
# ## Genomic control parameter is 1.086, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/affymetrix500/uc/allchr_affymetrix500_eur_all_step2_uc_eur_sex_PCs_firthse_uc.regenie'
# ## Genomic control parameter is 0.681, no adjustment made
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/humancoreexome/uc/allchr_humancoreexome_eur_all_step2_uc_eur_sex_PCs_firthse_uc.regenie'
# ## Genomic control parameter is 1.066, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/humanomni1/uc/allchr_humanomni1_eur_all_step2_uc_eur_sex_PCs_firthse_uc.regenie'
# ## Genomic control parameter is 0.990, no adjustment made
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/quad610/uc/allchr_quad610_eur_all_step2_uc_eur_sex_PCs_firthse_uc.regenie'
# ## Genomic control parameter is 1.040, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/gsa/uc/allchr_gsa_eur_all_step2_uc_eur_sex_PCs_firthse_uc.regenie'
# ## Genomic control parameter is 1.075, adjusting test statistics
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/illuminaexome/uc/allchr_illuminaexome_eur_all_step2_uc_eur_sex_PCs_firthse_uc.regenie'
# ## Genomic control parameter is 1.151, adjusting test statistics



#############################

# remove gc files


pheno=(ibd cd uc)
ancestry=(eur_all)

for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do 
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/allchr_eur_all_${ph}_meta_withGC_PCs_firthse_withGC_1.txt
done
done


###################################################################################################

###############
### no GC ###
###############

  
pheno=(ibd cd uc)
ancestry=(eur_all)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=800

for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
for chr in {1..22}
do
bsub -J"noGC_meta" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/chr${chr}_metal_noGC_PCs_firthse_${ph}_${j}_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/chr${chr}_metal_noGC_PCs_firthse_${ph}_${j}_stdout \
"sh ${path_gwas}scripts/run_metal_noGC_PCs_firthse_2022.sh ${ph} ${j} ${chr}"
done
done
done

pheno=(ibd cd uc)
ancestry=(eur_all)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=800

for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
for chr in X
do
bsub -J"noGC_meta" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/chr${chr}_metal_noGC_PCs_firthse_${ph}_${j}_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/chr${chr}_metal_noGC_PCs_firthse_${ph}_${j}_stdout \
"sh ${path_gwas}scripts/run_metal_noGC_PCs_firthse_chrX_2022.sh ${ph} ${j} ${chr}"
done
done
done



for j in ${ancestry[@]}
do
echo ${j} && for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22} X
do
echo ${chr} && tail -50   ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/chr${chr}_metal_noGC_PCs_firthse_${ph}_${j}_stdout | grep "Successfully"
done
done
done


# compress output results:

pheno=(ibd cd uc)
ancestry=(eur_all)

path_gwas=/path/to/ibdgwas/IIBDGC/
  
for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
for chr in X
do gzip -f ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/chr${chr}_${j}_${ph}_meta_noGC_PCs_firthse_1.txt
done
done
done

######################################
# reformat summary statistic results

pheno=(ibd cd uc)
ancestry=(eur_all)
MEM=30000

for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
for chr in {1..22} X
do
bsub -J"format" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 2 \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${j}_${chr}_format_metal_files_eur_tier_1_withNeff_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${j}_${chr}_format_metal_files_eur_tier_1_withNeff_stdout \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/format_metal_files_eur_tier_1_withNeff.R ${ph} ${chr} \
> ${path_gwas}scripts/logs/format_metal_files_eur_tier_1_withNeff_${ph}_${j}_${chr}.Rout"
done
done
done


for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22} X
do
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${j}_${chr}_format_metal_files_eur_tier_1_withNeff_stdout | grep "Successfully"
done
done
done

for ph in ${pheno[@]}
do
echo ${ph} && ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/*_${ph}_meta_eur_tier_1_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz | wc -l
done


for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22} X
do
echo ${chr} \
&& zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/${chr}_${ph}_meta_eur_tier_1_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz | wc -l
done
done


######################################
# add rsid db154

for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
for chr in {1..22} X
do
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/${chr}_${ph}_meta_eur_tier_1_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz \
| cut -f2-4 | awk -v chr=${chr} 'BEGIN {OFS = ""} {print chr,"\t",$1,"\t",$2,",",$3}' | sed '1d' > \
${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/iibdgc/chr${chr}_${ph}_${j}_list_chr_position_alleles
done
done
done

for chr in {1..22} X
do 
cat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/iibdgc/chr${chr}_ibd_${j}_list_chr_position_alleles ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/iibdgc/chr${chr}_cd_${j}_list_chr_position_alleles ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/iibdgc/chr${chr}_uc_${j}_list_chr_position_alleles | \
 sort -nk2,2 | uniq > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/iibdgc/chr${chr}_${j}_list_chr_position_alleles
done


# remove previous files:
for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
for chr in {1..22} X
do 
rm ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/iibdgc/chr${chr}_${ph}_${j}_list_chr_position_alleles
done
done
done


MEM=500

for chr in {1..22} X
do 
bsub -J"format" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 2 \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/chr${chr}_format_metal_2_IIBDGC_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/chr${chr}_format_metal_2_IIBDGC_stdout \
"/path/to/software/bcftools-1.16/./bcftools view -T ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/iibdgc/chr${chr}_${j}_list_chr_position_alleles \
${path_gwas}resources/dbsnp154/GRCh38.dbSNP154_multiallelic_split.vcf.gz | \
cut -f1-5 \
> ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/iibdgc/chr${chr}_${j}_list_chr_position_alleles_with_rsid_dbsnp154"
done


path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=3000

for chr in {1..22} X
do 
bsub -J"rsid" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G ibdgwas \
-e ${path_gwas}post_imputation/analysis/metaanalysis/log/add_rsid_${chr}_stderr \
-o ${path_gwas}post_imputation/analysis/metaanalysis/log/add_rsid_${chr}_stdout \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/add_rsid_to_formatted_metal_files_2022.R ${chr} > \
${path_gwas}scripts/logs/add_rsid_to_formatted_metal_files_${chr}.Rout"
done


for ph in ${pheno[@]}
do
for chr in {1..22} X
do 
ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/${chr}_${ph}_meta_eur_tier_1_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff_with_rsid_dbsnp154.txt.gz 
done
done


# remove intermediate files 

for ph in ${pheno[@]}
do ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/allchr_*
done



################################################
# plot results

# MEM=10000
# bsub -J"plot_imp3" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -m "modern_hardware" -G ibdgwas -n 2 \
# -e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/draw_manhattan_plots_stderr \
# -o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/draw_manhattan_plots_stdout \
# "/software/R-4.3.1/bin/Rscript ${path_gwas}scripts/draw_manhattan_plots.R > \
# ${path_gwas}scripts/logs/draw_manhattan_plots.Rout"

## plot manhattan with sas summary stats:

MEM=10000

release="eur_tier_1"

bsub -J"plot_imp1" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -m "modern_hardware" -G ibdgwas -n 2 \
-e ${path_gwas}scripts/logs/draw_manhattan_plots_stderr_${release} \
-o ${path_gwas}scripts/logs/draw_manhattan_plots_stdout_${release} \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/draw_manhattan_plots.R ${release} > \
${path_gwas}scripts/logs/draw_manhattan_plots_${release}.Rout"

  
################################################
# get number of variants per % samples included

MEM=15000

bsub -J"plot_imp4" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -m "modern_hardware" -G ibdgwas -n 2 \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/plot_summary_number_variants_per_sample_size_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/plot_summary_number_variants_per_sample_size_stdout \
"/software/R-4.3.1/bin/Rscript ${path_gwas}scripts/plot_summary_number_variants_per_sample_size.R \
> ${path_gwas}scripts/logs/plot_summary_number_variants_per_sample_size.Rout"


~/Desktop/

  
#####################################
### no GC - sensitivity analysis ###
#####################################

# 
# for i in ${array[@]}
# do
# echo ${i} && less ${path_gwas}scripts/run_metal_noGC_PCs_firthse_no_${i}_2022.sh | grep '# PROCESS'
# done

array=(illumina370 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)
pheno=(ibd cd uc)
ancestry=(eur_all)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=30000

for i in ${array[@]}
do
for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
bsub -J"noGC_meta" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_no_${i}_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_no_${i}_stdout \
"sh ${path_gwas}scripts/run_metal_noGC_PCs_firthse_no_${i}_2022.sh ${ph} ${j} ${i}"
done
done
done


for i in ${array[@]}
do
echo ${i} && for j in ${ancestry[@]}
do
echo ${j} && for ph in ${pheno[@]}
do
echo ${ph} && tail -50 ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_no_${i}_stdout | grep "Successfully"
done
done
done
done

# file with high het driven by illumianexome saved as 'cd_list_variants_high_HetPval_sensitivity_discrepancy_illuminaexome.xlsx'

###########################
# # illumina370
# # eur_all
# # ibd
# # Successfully completed.
# # cd
# # uc
# # Successfully completed.
# 
# # affymetrix6
# # eur_all
# # ibd
# # Successfully completed.
# # cd
# # Successfully completed.
# # uc
# # Successfully completed.
# 
# # humanomniexpress
# # eur_all
# # ibd
# # Successfully completed.
# # cd
# # Successfully completed.
# # uc
# # Successfully completed.
# 
# # affymetrix500
# # eur_all
# # ibd
# # Successfully completed.
# # cd
# # Successfully completed.
# # uc
# # Successfully completed.
# 
# # humancoreexome
# # eur_all
# # ibd
# # Successfully completed.
# # cd
# # Successfully completed.
# # uc
# # Successfully completed.
# 
# # humanomni1
# # eur_all
# # ibd
# # Successfully completed.
# # cd
# # Successfully completed.
# # uc
# # Successfully completed.
# 
# # quad610
# # eur_all
# # ibd
# # Successfully completed.
# # cd
# # Successfully completed.
# # uc
# # Successfully completed.
# 
# # gsa
# # eur_all
# # ibd
# # Successfully completed.
# # cd
# # Successfully completed.
# # uc
# # Successfully completed.
# 
# # illuminaexome
# # eur_all
# # ibd
# # Successfully completed.
# # cd
# # Successfully completed.
# # uc
# # Successfully completed.
###########################

###################
# plot results

  
path_gwas=/path/to/ibdgwas/IIBDGC/
  
MEM=8000

array=(illumina370 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)
pheno=(ibd cd uc)

ancestry=(eur_all)

for i in ${array[@]}
do
for ph in ${pheno[@]}
do
bsub -J"senst" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${i}_${ph}_sensitivity_test_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${i}_${ph}_sensitivity_test_stdout \
"/software/R-4.3.1/bin/Rscript ${path_gwas}scripts/plot_sensitivity_test_meta_analysis_results.R ${ph} ${i} \
> ${path_gwas}scripts/logs/sensitivity_test_meta_analysis_results_${ph}_${i}.Rout"
done
done

for i in ${array[@]}
do
echo ${i} && for j in ${ancestry[@]}
do
echo ${j} && for ph in ${pheno[@]}
do
echo ${ph} && tail -50 ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${i}_${ph}_sensitivity_test_stdout | grep "Successfully"
done
done
done
done


  

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################




########################################################################################################################
# compress all regenie file

array=(illumina370 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)
pheno=(ibd cd uc)
ancestry=(eur_all)

path_gwas=/path/to/ibdgwas/IIBDGC/

for i in ${array[@]}
do
for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
gzip ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/allchr_${i}_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie
done
done
done

#########

array=(illumina370 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)
pheno=(ibd cd uc)
ancestry=(eur_all)

path_gwas=/path/to/ibdgwas/IIBDGC/
  
for i in ${array[@]}
do
for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
for chr in {1..22} X
do
gzip ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/chr${chr}_${i}_${j}_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie
done
done
done
done



########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

# TIER 2 EUR

# run sensitibity analysis with TIER 2 EUR datasets:

#####################################
### no GC - sensitivity analysis ###
#####################################

# 
# for i in ${array[@]}
# do
# echo ${i} && less ${path_gwas}scripts/run_metal_noGC_PCs_firthse_no_${i}_2022.sh | grep '# PROCESS'
# done

###############################
# combine all files from ukbb


path_gwas="/path/to/ibdgwas/IIBDGC/"
pheno=(ibd cd uc)

for ph in ${pheno[@]}
do 
zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/chr1_ukbb_eur_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz | head -1  > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/tailer_regenie
zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/chr*_ukbb_eur_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz | grep -v "CHROM" > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/tmp.regenie
cat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/tailer_regenie ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/tmp.regenie \
| gzip > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/allchr_ukbb_eur_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz
rm ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/tailer_regenie 
rm ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/tmp.regenie
done


for ph in ${pheno[@]}
do 
echo ${ph} && zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/${ph}/allchr_ukbb_eur_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz | wc -l
done

# 151487379
# 151104448
# 151298423

# run meta analysis

study=(ukbb)
pheno=(ibd cd uc)
ancestry=(eur_all)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=30000

for i in ${study[@]}
do
for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
bsub -J"noGC_meta" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_no_${i}_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_no_${i}_stdout \
"sh ${path_gwas}scripts/run_metal_noGC_PCs_firthse_tier_2_with_${i}_2022.sh ${ph} ${j} ${i}"
done
done
done



###############################
# combine all files from finngen

path_gwas="/path/to/ibdgwas/IIBDGC/"
pheno=(ibd cd uc)

for ph in ${pheno[@]}
do 
cat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/finngen/${ph}/chr1_finngen_r10_K11_${ph}.regenie | tail -1  > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/finngen/${ph}/tailer_regenie
cat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/finngen/${ph}/chr*_finngen_r10_K11_${ph}.regenie | grep -v "CHROM" > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/finngen/${ph}/tmp.regenie
cat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/finngen/${ph}/tailer_regenie ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/finngen/${ph}/tmp.regenie \
> ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/finngen/${ph}/allchr_finngen_r10_K11_${ph}.regenie
rm ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/finngen/${ph}/tailer_regenie 
rm ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/finngen/${ph}/tmp.regenie
done

for ph in ${pheno[@]}
do 
gzip ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/finngen/${ph}/allchr_finngen_r10_K11_${ph}.regenie
done


# run meta analysis

study=(finngen)
pheno=(ibd cd uc)
ancestry=(eur_all)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=30000

for i in ${study[@]}
do
for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
bsub -J"noGC_meta" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_no_${i}_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_no_${i}_stdout \
"sh ${path_gwas}scripts/run_metal_noGC_PCs_firthse_tier_2_with_${i}_2022.sh ${ph} ${j} ${i}"
done
done
done



###############################
# combine all files from the danish study


path_gwas="/path/to/ibdgwas/IIBDGC/"
pheno=(ibd cd uc)

for ph in ${pheno[@]}
do 
cat ${path_gwas}post_imputation/analysis/stage_2_summary_statistics/danish/${ph}/chr1_danish_gsa_${ph}_eur_sex_10PCs_saige_spa.regenie | tail -1  > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/danish/${ph}/tailer_regenie
cat ${path_gwas}post_imputation/analysis/stage_2_summary_statistics/danish/${ph}/chr*_danish_gsa_${ph}_eur_sex_10PCs_saige_spa.regenie | grep -v "CHROM" > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/danish/${ph}/tmp.regenie
cat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/danish/${ph}/tailer_regenie ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/danish/${ph}/tmp.regenie > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/danish/${ph}/allchr_danish_gsa_${ph}_eur_sex_10PCs_saige_spa.regenie
rm ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/danish/${ph}/tailer_regenie 
rm ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/danish/${ph}/tmp.regenie
done

for ph in ${pheno[@]}
do 
gzip ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/danish/${ph}/allchr_danish_gsa_${ph}_eur_sex_10PCs_saige_spa.regenie
done


# run meta analysis

study=(danish)
pheno=(ibd cd uc)
ancestry=(eur_all)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=30000

for i in ${study[@]}
do
for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
bsub -J"noGC_meta" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_no_${i}_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_no_${i}_stdout \
"sh ${path_gwas}scripts/run_metal_noGC_PCs_firthse_tier_2_with_${i}_2022.sh ${ph} ${j} ${i}"
done
done
done


###############################
# combine all files from decode


path_gwas="/path/to/ibdgwas/IIBDGC/"
pheno=(ibd cd uc)

for ph in ${pheno[@]}
do 
cat ${path_gwas}post_imputation/analysis/stage_2_summary_statistics/decode/${ph}/chr1_DECODE_${ph}_30062021_edited_noMult.regenie | tail -1  > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/decode/${ph}/tailer_regenie
cat ${path_gwas}post_imputation/analysis/stage_2_summary_statistics/decode/${ph}/chr*_DECODE_${ph}_30062021_edited_noMult.regenie | grep -v "CHROM" > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/decode/${ph}/tmp.regenie
cat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/decode/${ph}/tailer_regenie ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/decode/${ph}/tmp.regenie |
gzip > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/decode/${ph}/allchr_decode_${ph}_30062021_edited_noMult.regenie.gz
rm ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/decode/${ph}/tailer_regenie 
rm ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/decode/${ph}/tmp.regenie
done




# run meta analysis

study=(decode)
pheno=(ibd cd uc)
ancestry=(eur_all)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=10000

for i in ${study[@]}
do
for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
bsub -J"noGC_meta" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_no_${i}_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_no_${i}_stdout \
"sh ${path_gwas}scripts/run_metal_noGC_PCs_firthse_tier_2_with_${i}_2022.sh ${ph} ${j} ${i}"
done
done
done


###############################
# combine all files from ibdbioresource

# zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/${chr}_interval_ibdbioresource_eur_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz | \
# grep chr16:50729867:G:GC

path_gwas="/path/to/ibdgwas/IIBDGC/"
pheno=(ibd cd uc)

for ph in ${pheno[@]}
do 
zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/1_interval_ibdbioresource_eur_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz | head -1  > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/header_regenie
zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/*_interval_ibdbioresource_eur_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz | grep -v "CHROM" > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/tmp.regenie
cat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/header_regenie ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/tmp.regenie | gzip \
> ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/allchr_interval_ibdbioresource_eur_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz
rm ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/header_regenie 
rm ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/tmp.regenie
done

for ph in ${pheno[@]}
do 
echo ${ph} && zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/${ph}/allchr_interval_ibdbioresource_eur_step2_${ph}_eur_sex_PCs_firthse_${ph}.regenie.gz | wc -l
done

# ibd
# 77628758
# cd
# 75211407
# uc
# 75242819


# run meta analysis

study=(interval_ibdbioresource)
pheno=(ibd cd uc)
ancestry=(eur_all)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=18000

for i in ${study[@]}
do
for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
bsub -J"noGC_meta" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_no_${i}_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_no_${i}_stdout \
"sh ${path_gwas}scripts/run_metal_noGC_PCs_firthse_tier_2_with_${i}_2022.sh ${ph} ${j} ${i}"
done
done
done


#############################
# run meta analysis with eas

study=(eas)
pheno=(ibd cd uc)
ancestry=(eur_all)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=10000

for i in ${study[@]}
do
for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
bsub -J"noGC_meta" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_no_${i}_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_no_${i}_stdout \
"sh ${path_gwas}scripts/run_metal_noGC_PCs_firthse_tier_2_with_${i}_2022.sh ${ph} ${j} ${i}"
done
done
done


###############################

path_gwas=/path/to/ibdgwas/IIBDGC/
study=(ukbb finngen danish decode eas interval_ibdbioresource)
j=eur_all

for i in ${study[@]}
do
echo ${i} && for ph in ${pheno[@]}
do ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_no_${i}_stdout | tail -2 && \
tail -50  ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_no_${i}_stdout | grep "Successfully"
done
done


###################################
# # ukbb
# # Successfully completed.
# # Successfully completed.
# # Successfully completed.
# # finngen
# # Successfully completed.
# # Successfully completed.
# # Successfully completed.
# # danish
# # Successfully completed.
# # Successfully completed.
# # Successfully completed.
# # decode
# # Successfully completed.
# # Successfully completed.
# # Successfully completed.
# # eas
# # Successfully completed.
# # Successfully completed.
# # Successfully completed.
# interval_ibdbioresource
# Successfully completed.
# Successfully completed.
# Successfully completed.
# 
# ############# ukbb
# 
# ### ibd
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/illumina370/ibd/allchr_illumina370_eur_all_step2_ibd_eur_sex_PCs_firthse_ibd.regenie'
# # 15964870 lines filtered, for the following reasons:
# # 15654764 lines filtered due to constraint A1FREQ >= 0.001
# #       10 lines filtered due to constraint A1FREQ <= 0.999
# #   310096 lines filtered due to constraint INFO > 0.4
# ## Processed 18548633 markers ...
# ## Smallest p-value is 4.90e-267 at marker 'chr6:32333650:C:T'
# 
# ### cd
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/cd/allchr_ukbb_eur_step2_cd_eur_sex_PCs_firthse_cd.regenie'
# # 125316279 lines filtered, for the following reasons:
# # 125122120 lines filtered due to constraint A1FREQ >= 0.001
# #     1392 lines filtered due to constraint A1FREQ <= 0.999
# #   192767 lines filtered due to constraint INFO > 0.4
# ## Processed 13474871 markers ...
# ## Smallest p-value is 4.27e-311 at marker 'chr16:50729867:G:GC'
# 
# 
# ### uc
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/uc/allchr_ukbb_eur_step2_uc_eur_sex_PCs_firthse_uc.regenie'
# # 125493615 lines filtered, for the following reasons:
# # 125299476 lines filtered due to constraint A1FREQ >= 0.001
# #     1394 lines filtered due to constraint A1FREQ <= 0.999
# #   192745 lines filtered due to constraint INFO > 0.4
# ## Processed 13474946 markers ...
# ## Smallest p-value is 2.52e-220 at marker 'chr6:32333650:C:T'
# 
# ############# finngen
# 
# ### ibd
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/stage_2_summary_statistics/finngen/cd/allchr_finngen_r10_K11_ibd.regenie'
# #  6691700 lines filtered, for the following reasons:
# #  6669200 lines filtered due to constraint A1FREQ >= 0.001
# #    22500 lines filtered due to constraint A1FREQ <= 0.999
# ## Processed 14614427 markers ...
# ## Smallest p-value is 6.56e-234 at marker 'chr17:53882759:C:A'
# 
# 
# ### cd
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/stage_2_summary_statistics/finngen/cd/allchr_finngen_r10_K11_cd.regenie'
# #  6692091 lines filtered, for the following reasons:
# #  6669590 lines filtered due to constraint A1FREQ >= 0.001
# #    22501 lines filtered due to constraint A1FREQ <= 0.999
# ## Processed 14614252 markers ...
# ## Smallest p-value is 1.37e-305 at marker 'chr16:50729867:G:GC'
# 
# ### uc
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/stage_2_summary_statistics/finngen/uc/allchr_finngen_r10_K11_uc.regenie'
# #  6692070 lines filtered, for the following reasons:
# #  6669569 lines filtered due to constraint A1FREQ >= 0.001
# #    22501 lines filtered due to constraint A1FREQ <= 0.999
# ## Processed 14614268 markers ...
# ## Smallest p-value is 4.71e-195 at marker 'chr6:32625390:C:T'
# 
# 
# ############# danish
# 
# ### ibd
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/stage_2_summary_statistics/danish/ibd/allchr_danish_gsa_ibd_eur_sex_10PCs_saige_spa.regenie'
# ## WARNING: Marker 'chr10:34058838:C:CCG' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr11:48651855:G:A' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr15:101074678:GA:G' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr19:17189326:TAAAACAAAAC:T' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr1:154869854:TTGCTGCTGC:T' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr20:31482970:C:CA' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr2:52579459:G:A' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr2:55439239:G:A' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr3:3359130:C:CG' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr3:3359130:C:T' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr7:34643641:C:T' duplicated in input, first instance used, others skipped
# ## Although a filter is defined, no lines were excluded because of it
# ## Processed 13950789 markers ...
# ## Smallest p-value is 6.56e-234 at marker 'chr17:53882759:C:A'
# 
# ### cd
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/stage_2_summary_statistics/danish/cd/allchr_danish_gsa_cd_eur_sex_10PCs_saige_spa.regenie'
# ## WARNING: Marker 'chr10:34058838:C:CCG' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr11:48651855:G:A' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr15:101074678:GA:G' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr19:17189326:TAAAACAAAAC:T' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr1:154869854:TTGCTGCTGC:T' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr20:31482970:C:CA' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr2:52579459:G:A' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr2:55439239:G:A' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr3:3359130:C:CG' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr3:3359130:C:T' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr7:34643641:C:T' duplicated in input, first instance used, others skipped
# ## Although a filter is defined, no lines were excluded because of it
# ## Processed 13975773 markers ...
# ## Smallest p-value is 1.06e-304 at marker 'chr16:50729867:G:GC'
# 
# ### uc
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/stage_2_summary_statistics/danish/uc/allchr_danish_gsa_uc_eur_sex_10PCs_saige_spa.regenie'
# ## WARNING: Marker 'chr10:34058838:C:CCG' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr11:48651855:G:A' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr15:101074678:GA:G' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr19:17189326:TAAAACAAAAC:T' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr1:154869854:TTGCTGCTGC:T' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr20:31482970:C:CA' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr2:55439239:G:A' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr3:3359130:C:CG' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr3:3359130:C:T' duplicated in input, first instance used, others skipped
# ## WARNING: Marker 'chr7:34643641:C:T' duplicated in input, first instance used, others skipped
# ## Although a filter is defined, no lines were excluded because of it
# ## Processed 13967802 markers ...
# ## Smallest p-value is 1.31e-191 at marker 'chr6:32653132:A:G'
# 
# ############# decode
# 
# ### ibd
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/stage_2_summary_statistics/decode/ibd/allchr_decode_ibd_30062021_edited_noMult.regenie'
# ## Although a filter is defined, no lines were excluded because of it
# ## Processed 16146964 markers ...
# ## Smallest p-value is 6.56e-234 at marker 'chr17:53882759:C:A'
# 
# ### cd
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/stage_2_summary_statistics/decode/cd/allchr_decode_cd_30062021_edited_noMult.regenie'
# ## Although a filter is defined, no lines were excluded because of it
# ## Processed 16146964 markers ...
# ## Smallest p-value is 1.44e-292 at marker 'chr16:50729867:G:GC'
# 
# ### uc
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/stage_2_summary_statistics/decode/uc/allchr_decode_uc_30062021_edited_noMult.regenie'
# ## Although a filter is defined, no lines were excluded because of it
# ## Processed 16146964 markers ...
# ## Smallest p-value is 1.87e-190 at marker 'chr6:32653132:A:G'
# 
# ############# eas
# 
# ### ibd
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/illumina370/ibd/allchr_illumina370_eur_all_step2_ibd_eur_sex_PCs_firthse_ibd.regenie'
# # 15964870 lines filtered, for the following reasons:
# # 15654764 lines filtered due to constraint A1FREQ >= 0.001
# #       10 lines filtered due to constraint A1FREQ <= 0.999
# #   310096 lines filtered due to constraint INFO > 0.4
# ## Processed 18548633 markers ...
# 
# ### cd
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/stage_2_summary_statistics/eas_iibdgc/cd/ibd_EAS_SiKJ_meta_cd.regenie'
# #     4760 lines filtered, for the following reasons:
# #     4758 lines filtered due to constraint A1FREQ >= 0.001
# #        2 lines filtered due to constraint A1FREQ <= 0.999
# ## Processed 13493322 markers ...
# 
# ### uc
# ## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/stage_2_summary_statistics/eas_iibdgc/uc/ibd_EAS_SiKJ_meta_uc.regenie'
# #     5185 lines filtered, for the following reasons:
# #     5183 lines filtered due to constraint A1FREQ >= 0.001
# #        2 lines filtered due to constraint A1FREQ <= 0.999
# ## Processed 13489536 markers ...
# 
# 
# ############# interval_ibdbioresource
# ### ibd
# # Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/ibd/allchr_interval_ibdbioresource_eur_step2_ibd_eur_sex_PCs_firthse_ibd.regenie.gz'
# # WARNING: Marker '10_50134753_G_A' duplicated in input, first instance used, others skipped
# 61819120 lines filtered, for the following reasons:
# 61105451 lines filtered due to constraint A1FREQ >= 0.001
#   213856 lines filtered due to constraint A1FREQ <= 0.999
#   499813 lines filtered due to constraint INFO > 0.4
## Processed 15809636 markers ...
# 
# ### cd
# # Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/cd/allchr_interval_ibdbioresource_eur_step2_cd_eur_sex_PCs_firthse_cd.regenie.gz'
## WARNING: Marker 'chr10:50134753:G:A' duplicated in input, first instance used, others skipped
# 56813594 lines filtered, for the following reasons:
# 56243225 lines filtered due to constraint A1FREQ >= 0.001
#   208731 lines filtered due to constraint A1FREQ <= 0.999
#   361638 lines filtered due to constraint INFO > 0.4
## Processed 15339123 markers ...
#
# ### uc
## Processing file '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/uc/allchr_interval_ibdbioresource_eur_step2_uc_eur_sex_PCs_firthse_uc.regenie.gz'
## WARNING: Marker '10_50134753_G_A' duplicated in input, first instance used, others skipped
# 56853563 lines filtered, for the following reasons:
# 56281688 lines filtered due to constraint A1FREQ >= 0.001
#   208494 lines filtered due to constraint A1FREQ <= 0.999
#   363381 lines filtered due to constraint INFO > 0.4
## Processed 15333550 markers ...
# ###############
# 
# # inspect regions driven highest heterogeneity:
###################################

###################
# plot results


path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=35000

study=(ukbb finngen danish decode eas)
study=(interval_ibdbioresource)
pheno=(cd uc)
ancestry=(eur_all)

for i in ${study[@]}
do
for ph in ${pheno[@]}
do
bsub -J"senst" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${i}_${ph}_sensitivity_test_tier_2_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${i}_${ph}_sensitivity_test_tier_2_stdout \
"/software/R-4.3.1/bin/Rscript ${path_gwas}scripts/plot_sensitivity_test_meta_analysis_results_tier_2.R ${ph} ${i} \
> ${path_gwas}scripts/logs/sensitivity_test_meta_analysis_results_tier_2_${ph}_${i}.Rout"
done
done

# continue here


for i in ${study[@]}
do
echo ${i} && for ph in ${pheno[@]}
do
echo ${ph} && tail -50 ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${i}_${ph}_sensitivity_test_tier_2_stdout | grep "Successfully"
done
done




########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

# TIER 2 EUR

# run meta-analysis with all TIER 2 EUR datasets (minus IBD BioResource+Interval):
path_gwas="/path/to/ibdgwas/IIBDGC/"
pheno=(ibd cd uc)
ancestry=(eur_all)

MEM=15000

for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
bsub -J"noGC_meta" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_tier_2_no_ibdbioresource_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_tier_2_no_ibdbioresource_stdout \
"sh ${path_gwas}scripts/run_metal_noGC_PCs_firthse_eur_tier_2_2022_no_ibdbioresource.sh ${ph} ${j}"
done
done

for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
echo ${ph} && tail -50 ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_tier_2_no_ibdbioresource_stdout | grep "Successfully"
done
done



######################################################
# run meta-analysis with ALL TIER 2 EUR datasets :

# run meta-analysis with ALL TIER 2 EUR datasets - with GC:

path_gwas="/path/to/ibdgwas/IIBDGC/"
pheno=(ibd cd uc)
ancestry=(eur_all)

MEM=35000

for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
bsub -J"noGC_meta" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_withGC_PCs_firthse_${ph}_${j}_tier_2_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_withGC_PCs_firthse_${ph}_${j}_tier_2_stdout \
"sh ~/git/IIBDGC_GWAS/scripts/other/run_metal_withGC_PCs_firthse_eur_tier_2_2022.sh ${ph} ${j}"
done
done



path_gwas="/path/to/ibdgwas/IIBDGC/"
pheno=(ibd cd uc)
ancestry=(eur_all)

MEM=35000

for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
bsub -J"noGC_meta" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_tier_2_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_tier_2_stdout \
"sh ~/git/IIBDGC_GWAS/scripts/other/run_metal_noGC_PCs_firthse_eur_tier_2_2022.sh ${ph} ${j}"
done
done


for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do echo ${ph} && ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_tier_2_stdout | tail -2 && \
tail -50 ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/allchr_metal_noGC_PCs_firthse_${ph}_${j}_tier_2_stdout | grep "Successfully"
done
done



######################################
# reformat summary statistic results

# split file by chr:

pheno=(cd)

for ph in ${pheno[@]}
do 
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/allchr_eur_all_${ph}_meta_noGC_PCs_firthse_eur_tier_2_1.txt | head -1 > \
${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/header
done

MEM=2000
j=eur_all
for ph in ${pheno[@]}
do 
echo ${ph} && for chr in {1..22} X
do
bsub -J"noGC_meta" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/chr${chr}_metal_noGC_PCs_firthse_${ph}_${j}_tier_2_split_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/chr${chr}_metal_noGC_PCs_firthse_${ph}_${j}_tier_2_split_stdout \
"echo ${chr} && grep chr${chr}: ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/allchr_eur_all_${ph}_meta_noGC_PCs_firthse_eur_tier_2_1.txt \
| gzip >  ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/chr${chr}_eur_all_${ph}_meta_noGC_PCs_firthse_eur_tier_2_1_tmp.txt.gz"
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
do echo ${chr} && ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/chr${chr}_eur_all_${ph}_meta_noGC_PCs_firthse_eur_tier_2_1_tmp.txt.gz
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
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${chr}_format_metal_eur_tier_2_withNeff_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${chr}_format_metal_eur_tier_2_withNeff_stdout \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/format_metal_files_eur_tier_2_withNeff.R ${ph} ${chr} \
> ${path_gwas}scripts/logs/format_metal_files_eur_tier_2_withNeff_${ph}_${chr}.Rout"
done
done



for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22} X
do echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${chr}_format_metal_eur_tier_2_withNeff_stdout | grep -E "Successfully|Exited"
done
done

for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22} X
do echo ${chr} && cat ${path_gwas}scripts/logs/format_metal_files_eur_tier_2_withNeff_${ph}_${chr}.Rout
done
done


for ph in ${pheno[@]}
do
for chr in {1..22} X
do ls -la  ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/${chr}_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz
done
done

for ph in ${pheno[@]}
do
for chr in {1..22} X
do 
echo ${chr} && \
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/${chr}_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz | head -1 | wc -w \
&& zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/${chr}_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz | tail -1 | wc -w
done
done



# look for potential duplicates:

# for ph in ${pheno[@]}
# do
# echo ${ph} && for chr in {1..22} X
# do echo ${chr} && zcat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/${ph}/${chr}_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz | cut -f1 | sort | uniq --count | tail -3
# done
# done

## plot manhattan with sas summary stats:

MEM=6000

release="eur_tier_2"

bsub -J"plot_imp1" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -m "modern_hardware" -G ibdgwas -n 2 \
-e ${path_gwas}scripts/logs/draw_manhattan_plots_stderr_${release} \
-o ${path_gwas}scripts/logs/draw_manhattan_plots_stdout_${release} \
"/software/R-4.3.1/bin/Rscript ~/git/IIBDGC_GWAS/scripts/other/draw_manhattan_plots.R ${release} > \
${path_gwas}scripts/logs/draw_manhattan_plots_${release}.Rout"



###############################
# reformat file for coloc:

path_gwas="/path/to/ibdgwas/IIBDGC/"
pheno=(ibd cd uc)

MEM=1600

for ph in ${pheno[@]}
do
for chr in {1..22} X
do
bsub -J"formatcol" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${chr}_format_files_eur_tier_2_for_coloc_stderr \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${chr}_format_files_eur_tier_2_for_coloc_stdout \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/format_files_eur_tier_2_for_coloc_withNeff.R ${ph} ${chr} \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/format_summary_stats_files_eur_tier_2_for_coloc_withNeff_${ph}_${chr}.Rout"
done
done


for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22} X
do echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${chr}_format_files_eur_tier_2_for_coloc_stdout | grep -E "Successfully|Exited"
done
done


for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22} X
do echo ${chr} && cat  ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/format_summary_stats_files_eur_tier_2_for_coloc_withNeff_${ph}_${chr}.Rout
done
done



##### reformat final files for coloc:

cd /path/to/project

head -1 uc_chr22_summary_stats_metal_nogc_July2023_Neff_eur_tier2.tsv > header.tsv

pheno=(ibd cd uc)

for ph in ${pheno[@]}
do
cat ${ph}_chr*_summary_stats_metal_nogc_July2023_Neff_eur_tier2.tsv | grep -v "RSid" > tmp
sort -k2,2n -k3,3n tmp > tmp2
cat header.tsv tmp2 > ${ph}_allchr_summary_stats_nogc_July2023_Neff_eur_tier2_sorted.txt
bgzip ${ph}_allchr_summary_stats_nogc_July2023_Neff_eur_tier2_sorted.txt
tabix -s2 -b3 -e3 -S1 ${ph}_allchr_summary_stats_nogc_July2023_Neff_eur_tier2_sorted.txt.gz
done

for ph in ${pheno[@]}
do
echo ${ph} && zcat ${ph}_allchr_summary_stats_nogc_July2023_Neff_eur_tier2_sorted.txt.gz | wc -l
done

# ibd
# 13645384
# cd
# 13885188
# uc
# 13682870

for ph in ${pheno[@]}
do
mv ${ph}_allchr_summary_stats_nogc_July2023_Neff_eur_tier2_sorted.txt.gz.tbi ${ph}_allchr_summary_stats_nogc_Oct2023_Neff_eur_tier2_sorted.txt.gz.tbi
done

# remove intermediate files:
for chr in {1..22} X
do
rm *_chr*_summary_stats_metal_nogc_July2023_Neff_eur_tier2.tsv
done

rm tmp
rm tmp2

# create top hits 
for i in  *_Neff_eur_tier2_sorted.txt.gz; do g1=$(echo $i|cut -f1 -d "."); echo "awk '\$6<=1E-5' <(zcat ${i}) | bgzip -c > ${g1}.top_hits.txt.gz"; done
awk '$6<=1E-5' <(zcat cd_allchr_summary_stats_nogc_Oct2023_Neff_eur_tier2_sorted.txt.gz) | bgzip -c > cd_allchr_summary_stats_nogc_Oct2023_Neff_eur_tier2_sorted.top_hits.txt.gz
awk '$6<=1E-5' <(zcat ibd_allchr_summary_stats_nogc_Oct2023_Neff_eur_tier2_sorted.txt.gz) | bgzip -c > ibd_allchr_summary_stats_nogc_Oct2023_Neff_eur_tier2_sorted.top_hits.txt.gz
awk '$6<=1E-5' <(zcat uc_allchr_summary_stats_nogc_Oct2023_Neff_eur_tier2_sorted.txt.gz) | bgzip -c > uc_allchr_summary_stats_nogc_Oct2023_Neff_eur_tier2_sorted.top_hits.txt.gz

# update file names in ref file
/path/to/project


# cd /path/to/project


# save a copy of the full files for Jack to investigate:

pheno=(ibd cd uc)
path_gwas="/path/to/ibdgwas/IIBDGC/"

for ph in ${pheno[@]}
do
cp ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/*_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz /path/to/project
done