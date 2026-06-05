# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
#########################################################################################################################################################################

# singularity exec iibdgc_postprocess_10_singularity.sif
# singularity exec polyfun_2_singularity.sif

# path_gwas="/path/to/ibdgwas/IIBDGC/"
path_gwas="/path/to/ibdgwas/IIBDGC/"

pheno=(ibd cd uc)


# https://github.com/omerwe/polyfun/wiki/1.-Computing-prior-causal-probabilities-with-PolyFun#1-create-a-munged-summary-statistics-file-in-a-polyfun-friendly-parquet-format

# we will use step 3:
# Computing prior causal probabilities non-parametrically. This is the most robust approach, but it is computationally intensive and requires 
# access to individual-level genotypic data from a large reference panel (optimally >10,000 population-matched individuals), unless you are 
# analyzing summary statistics based only on UK Biobank British-ancestry individuals.

# USE EUR TIER 1 summary stat files for this

#########################################################
# 1.- PREPARE SUMMARY STATISTICS IN THE FORMAT REQUIRED #
#########################################################

# 1. Create a munged summary statistics file in a PolyFun-friendly parquet format.
# To do this, use the script munge_polyfun_sumstats.py, which takes an input summary statistics file and creates a munged output file. The script tries to be flexible and accommodate multiple file formats and column names. It generally requires only a sample size parameter (n) and a whitespace-delimited input file with SNP rsids, chromosome and base pair info, and either a p-value, an effect size estimate and its standard error, a Z-score or a p-value.

# Here is a usage example:

# python munge_polyfun_sumstats.py \
#   --sumstats example_data/boltlmm_sumstats.gz \
#   --n 327209 \
#   --out example_data/sumstats_munged.parquet \
#   --min-info 0.6 \
#   --min-maf 0.001
# This takes the input BOLT-LMM file example_data/boltlmm_sumstats.gz and converts it to the parquet file example_data/sumstats_munged.parquet, excluding SNPs with INFO score<0.6, with MAF<0.001 or in the MHC region. It will additionally compute the BOLT-LMM effective sample size. You can see other possible arguments with the command python munge_polyfun_sumstats.py --help. You can see the output file by opening the parquet file through python with the command df = pd.read_parquet('example_data/sumstats_munged.parquet')

# usage: munge_polyfun_sumstats.py [-h] --sumstats SUMSTATS --out OUT [--n N] [--min-info MIN_INFO] [--min-maf MIN_MAF] [--remove-strand-ambig]
#                                  [--chi2-cutoff CHI2_CUTOFF] [--keep-hla] [--no-neff]

# The sumstats_file should be either a parquet or a whitespace-delimited file (that can be gzipped) with a header line. It can accept 
# files with any combination of columns, as long as they include the following columns:
# CHR - chromosome
# BP - base pair position (in hg19 coordinates)
# A1 - The effect allele (i.e., the sign of the effect size is with respect to A1)
# A2 - the second allele

# alleles naming needs to be changed, effect allele here refers to A1

# BOL
# SNP	CHR	BP	GENPOS	ALLELE1	ALLELE0	A1FREQ	INFO	CHISQ_LINREG	P_LINREG	BETA	SE	CHISQ_BOLT_LMM_INF	P_BOLT_LMM_INF	CHISQ_BOLT_LMM	P_BOLT_LMM
# rs569799965	1	54945	0.0	C	A	0.9940129999999999	0.43888299999999997	1.69708	0.19	0.0235362	0.0216915	1.1773200000000001	0.28	1.2354399999999999	0.27
# rs138802168	1	362907	0.0	G	C	0.998551	0.436337	0.392685	0.53	0.0380167	0.049685599999999996	0.585449	0.44	0.445565	0.5


##################################################
# CREATE A REFERENCE UNION FILE FOR ALL VARIANTS

path_gwas="/path/to/ibdgwas/IIBDGC/"
pheno=(ibd cd uc)

for ph in ${pheno[@]}
do
for chr in {1..22} X
do
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/${chr}_${ph}_meta_eur_tier2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz | cut -f 1,3-5 | 
awk -v OFS="\t" -v myvar=${chr} '{print $0,myvar}' | sed '1d' \
>> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp_${ph}
done
done
  

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp_*
  # 34460879 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/tmp_cd
  # 34745784 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/tmp_ibd
  # 35873018 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/tmp_uc


# combine, sort and remove duplicates:
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp_cd \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp_ibd \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp_uc \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp_all

sort -k5,5n -k2,2n ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp_all \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp_all_sorted

cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp_all_sorted | \
cut -f 1 | uniq > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis
# 38680463

# remove intermediate files:
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp_*

# make a bed file:

awk 'BEGIN {FS = ":";OFS = "\t"} {print $1,$2,$2+1}' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp

paste -d'\t' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis.bed
  
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp

# 'V' sorts alphanumerically
sort -k1,1V -k2,2n ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis.bed \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_sorted.bed

wc -l  ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_sorted.bed
# 38680463 

rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis.bed


##############################################
# 1.1 - concatenate all summary stats files

ph=ibd
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/1_${ph}_meta_eur_tier2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz | head -1 \
| awk -F'\t' 'BEGIN {OFS=FS} {print $0,"chr"}' | sed 's/MarkerName/SNP/g' | sed 's/chr/CHR/g' | sed 's/Position_b38/BP/g' | sed 's/A1/ALLELE2/g' | sed 's/A2/ALLELE1/g' | sed 's/avgALLELE1FREQ/A1FREQ/g' > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/header_tmp

pheno=(ibd cd uc)

for ph in ${pheno[@]}
do 
echo ${ph} && rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_tmp
for chr in {1..22}
do
echo ${chr} && \
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/${chr}_${ph}_meta_eur_tier2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz | \
sed '1d' | awk -F'\t' 'BEGIN {OFS=FS} {print $0,'${chr}'}' \
>> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_tmp
done
done

# for ph in ${pheno[@]}
# do 
# cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_tmp | cut -f23 | sort | uniq -dc
# done

chr="X"
for ph in ${pheno[@]}
do 
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/${chr}_${ph}_meta_eur_tier2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz | \
sed '1d' | awk -F'\t' 'BEGIN {OFS=FS} {print $0,"23"}' \
>> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_tmp
done

for ph in ${pheno[@]}
do 
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_tmp | cut -f24 | sort | uniq -dc
done


for ph in ${pheno[@]}
do 
echo ${ph} && wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_tmp
done
# ibd
# 34745780 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_ibd_tmp
# cd
# 34460875 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_cd_tmp
# uc
# 35873014 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_uc_tmp

# add header
for ph in ${pheno[@]}
do 
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/header_tmp ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_tmp > \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_tmp2
done

# check out number of columns
for ph in ${pheno[@]}
do 
echo ${ph} && head -2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_tmp2 && \
tail -1 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_tmp2
done

for ph in ${pheno[@]}
do 
echo ${ph} && cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_tmp2 | cut -f24 | sort | uniq -dc
done

for ph in ${pheno[@]}
do
echo ${ph} && wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_tmp2
done
# ibd
# 34745781 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_ibd_tmp2
# cd
# 34460876 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_cd_tmp2
# uc
# 35873015 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_uc_tmp2


#################################################################
# # 1.2 - retain only those variants with rate_Neff >=0.5


for ph in ${pheno[@]}
do
awk -F'\t' '$23 >= 0.5 && $14 >= 0.9' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_tmp2 | cut -f1-7,14,16,21,23-25 | \
gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz
done


# for ph in ${pheno[@]}
# do
# echo ${ph} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz | cut -f11 | uniq -dc
# done

for ph in ${pheno[@]}
do 
echo ${ph} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz | wc -l
done
# # ibd
# 12674691
# cd
# 12840263
# uc
# 12833584

# summary Neff:
for ph in ${pheno[@]}
do 
echo ${ph} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz | sed '1d' | cut -f10 | sort -n | sed -n '1s/^/min=/p; $s/^/max=/p'
done
# ibd
# min=2e+05
# max=240386
# cd
# min=1e+05
# max=118922
# uc
# min=1e+05
# max=145099

# summary imputation info:
for ph in ${pheno[@]}
do 
echo ${ph} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz | sed '1d' | cut -f8 | sort -n | sed -n '1s/^/min=/p; $s/^/max=/p'
done
# ibd
# min=0.900000022395203
# max=1
# cd
# min=0.90000004003946
# max=1
# uc
# min=0.900000232774491
# max=1

# summary frequency:
for ph in ${pheno[@]}
do 
echo ${ph} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz | sed '1d' | cut -f9 | sort -n | sed -n '1s/^/min=/p; $s/^/max=/p'
done
# ibd
# min=0.00100031498623751
# max=0.998720295358052
# cd
# min=0.00100257038018178
# max=0.998801741664528
# uc
# min=0.00100058210613968
# max=0.9987490588975

# remove tmp files:
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_*_tmp2


## CREATE MUNGED.PARKET FILES

# files already thresholded at the levels we want, set up dummy thresholds so no more variants excluded - otherwise it will use defautls

pheno=(ibd cd uc)

# rename Neff to N so that munge can recognise it:
for ph in ${pheno[@]}
do
gunzip ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz && \
sed -i -e 's/Neff/N/g' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_eur_tier2_list_variants_rate_0.5_info_0.9.tsv && \
gzip ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_eur_tier2_list_variants_rate_0.5_info_0.9.tsv
done

for ph in ${pheno[@]}
do
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz | head 
done

# excluding HLA:
MEM=18000

for ph in ${pheno[@]}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${ph}_allchr_ldsc_sumstats_eur_tier2_rate_0.5_info_0.9_nohla_mungen_format_stderr \
-o ${path_gwas}post_imputation/log/${ph}_allchr_ldsc_sumstats_eur_tier2_rate_0.5_info_0.9_nohla_mungen_format_stdout \
"munge_polyfun_sumstats.py \
--sumstats ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz \
--min-info 0 \
--min-maf 0 \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/${ph}_eur_tier2_list_variants_rate_0.5_info_0.9_nohla_sumstats_munged.parquet"
done

for ph in ${pheno[@]}
do
tail -30 ${path_gwas}post_imputation/log/${ph}_allchr_ldsc_sumstats_eur_tier2_rate_0.5_info_0.9_nohla_mungen_format_stdout | grep -E "Successfully|Exited"
done

for ph in ${pheno[@]}
do
ls -la  ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/${ph}_eur_tier2_list_variants_rate_0.5_info_0.9_nohla_sumstats_munged.parquet
done


####################################################################################################
# CREATE A REFERENCE UNION FILE FOR EUR TIER 2 - 50% Neff samples contributing

path_gwas="/path/to/ibdgwas/IIBDGC/"
pheno=(ibd cd uc)

# combine, sort and remove duplicates:
cat <(zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_cd_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz | sed '1d') \
<(zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_ibd_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz | sed '1d') \
<(zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_uc_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz | sed '1d') \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp_all

sort -k12,12n -k2,2n ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp_all \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp_all_sorted

cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp_all_sorted | \
cut -f 1 | uniq > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9
# 12979489

# remove intermediate files:
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp_*


##############################################################################################################################################
## CREATE .BGEN FILES FROM GSA USING SAME FILES (and N samples) USED FOR REGENIE STEP2, BUT CONVERTED AND SUBSET BY VARIANTS INCLUDED IN FINAL ANALYSIS:

MEM=20000
array=(illumina370 humanomniexpress affymetrix500 humanomni1 quad610 illuminaexome affymetrix6 humancoreexome gsa)

pheno=(cd uc ibd)

for i in ${array[@]}
do
for ph in ${pheno[@]}
do
awk -F '\t' '$3 != "NA" { print }' ${path_gwas}post_imputation/2022/${i}/phenotype_data/${i}_all_studies_merged_eur_all_phenotype_${ph} \
| cut -f1-2 > ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_list_samples_included_in_${ph}_analysis
done
done

for i in ${array[@]}
do
echo ${i} && for ph in ${pheno[@]}
do 
echo ${ph} && wc -l ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_list_samples_included_in_${ph}_analysis
done
done


for i in ${array[@]}
do
for chr in {1..22} X
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 4 \
-e ${path_gwas}post_imputation/log/${i}_${chr}_ldsc_ldsc_plink_eur_tier2_stderr \
-o ${path_gwas}post_imputation/log/${i}_${chr}_ldsc_ldsc_plink_eur_tier2_stdout \
"plink2 \
--bgen ${path_gwas}post_imputation/2022/${i}/imputed_data/${i}_chr${chr}.dose.subset.bgen 'ref-first' \
--sample ${path_gwas}post_imputation/2022/${i}/imputed_data/${i}_chr${chr}.dose.subset.sample \
--keep ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_list_samples_included_in_ibd_analysis \
--extract ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9 \
--threads 4 \
--make-bed --out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/${i}_chr${chr}_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9"
done
done


for i in ${array[@]}
do
echo ${i} && for chr in {1..22} X
do
echo ${chr} && tail -30 ${path_gwas}post_imputation/log/${i}_${chr}_ldsc_ldsc_plink_eur_tier2_stdout | grep -E "Successfully|Exited"
done
done

for i in ${array[@]}
do
echo ${i} && ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/${i}_chr*_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9.bim | wc -l
done
done


array=(illumina370 humanomniexpress affymetrix500 humanomni1 quad610 illuminaexome affymetrix6 humancoreexome)

for chr in {1..22}
do
for i in ${array[@]}
do
echo ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/${i}_chr${chr}_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9 >> \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/${chr}_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_mergelist.txt
done
done

array=(humanomniexpress quad610 affymetrix6)

for chr in X
do
for i in ${array[@]}
do
echo ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/${i}_chr${chr}_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9 >> \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/${chr}_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_mergelist.txt
done
done

for chr in {1..22} X
do
ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/${chr}_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_mergelist.txt
done

MEM=45000
for chr in {1..22} X
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 4 \
-e ${path_gwas}post_imputation/log/bmerge_${chr}_ldsc_ldsc_plink_eur_tier2_stderr \
-o ${path_gwas}post_imputation/log/bmerge_${chr}_ldsc_ldsc_plink_eur_tier2_stdout \
"/path/to/software/username/plink_linux_x86_64_20231211/./plink \
--bfile ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/gsa_chr${chr}_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9 \
--merge-list ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/${chr}_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_mergelist.txt \
--keep-allele-order --threads 4 \
--make-bed --out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/allarrays_chr${chr}_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_rate_0.5_info_0.9"
done

for chr in {1..22} X
do
echo ${chr} && tail -30 ${path_gwas}post_imputation/log/bmerge_${chr}_ldsc_ldsc_plink_eur_tier2_stdout | grep -E 'Successfully|Exited'
done


## edit the list of variants to map the data we have genotyping array info for:

# MEM=8000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

path_gwas="/path/to/ibdgwas/IIBDGC/"
library(data.table)

pheno<-c("ibd","cd","uc")

for (chr in c(1:22)) {
  print(chr)
  
  tmp1<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/allarrays_chr",chr,"_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_rate_0.5_info_0.9.bim"))
   
  if(chr==1) {
    var<-tmp1
  }else{
    var<-rbind(var,tmp1)
  }
}

dim(var)

list<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9"),head=F)
dim(list)
# [1] 12979489        1

list<-list[which(list$V1 %in% var$V2),]
dim(list)
# [1] 12966826

write.table(list,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9"),
col.names=F,row.names=F,quote=F)

q("no")

######

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9
# 12966826

# remove intermediate plink files:

array=(illumina370 humanomniexpress affymetrix500 humanomni1 quad610 illuminaexome affymetrix6 humancoreexome)
for i in ${array[@]}
do
for chr in {1..22}
do
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/${i}_chr${chr}_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9.*
done
done


# make a bed file:

awk 'BEGIN {FS = ":";OFS = "\t"} {print $1,$2,$2+1}' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9 \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp

paste -d'\t' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9 \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9.bed
  
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp

# 'V' sorts alphanumerically
sort -k1,1V -k2,2n ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9.bed \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed

wc -l  ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed
#  12966826

rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9.bed



### estimate regression weights, no partitioned:

MEM=60000
for chr in 2
# MEM=55000
# for chr in 1 {3..6} 
# MEM=45000
# for chr in {7..8}
# MEM=40000
# for chr in {9..12}
# MEM=30000
# for chr in {13..16}
# MEM=25000
# for chr in {17..19}
# MEM=20000
# for chr in {20..22}
do
bsub -J"ldsc_w" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q week \
-e ${path_gwas}post_imputation/log/${chr}_ldsc_eur_tier2_rate_0.5_info_0.9_w_stderr \
-o ${path_gwas}post_imputation/log/${chr}_ldsc_eur_tier2_rate_0.5_info_0.9_w_stdout \
"ldsc.py \
--bfile ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/allarrays_chr${chr}_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9 \
--l2 \
--ld-wind-kb 1000 \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/eur_tier2_rate_0.5_info_0.9_allarrays/weights.april2025.${chr}"
done


for chr in {1..22}
do
echo ${chr} && tail -30 ${path_gwas}post_imputation/log/${chr}_ldsc_eur_tier2_rate_0.5_info_0.9_w_stdout | grep -E "Successfully|Exit"
done

for chr in {1..22}
do
echo ${chr} && tail -30 ${path_gwas}post_imputation/log/${chr}_ldsc_eur_tier2_rate_0.5_info_0.9_w_stdout | grep -E "Max Memory"
done



#########################################################################################################################################################################

#################################################################################
# 2.- DEFINE BASELINE ANNOTATION AND CELL TYPE SPEFICIC ANNOTATION USING S-LDSC #
#################################################################################

# https://github.com/omerwe/polyfun/wiki/5.-Estimating-functional-enrichment-using-S-LDSC

# ldsc.py -h

ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/ \
| grep -v list_all | grep -v ldsc.log | grep -v list_to_exclude | wc -l
# 205


# 2.1 CREATE BASELINE ANNOTATION FILE, FOLLOWING BASELINE USED IN WEISSBROD, NON CORRELATED FIELDS
# subset the original files, split by chr, and keep only the annotation:

# https://github.com/RHReynolds/LDSCforRyten
# We recommend that for estimating heritability enrichment (i.e., %h2/%SNPs) of any annotation, including tissue-specific annotations, 
# it is best to use baselineLD v2.2.

cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/

# exclude encode, remaining encode, Soskic et al, Zhang et al, and immune cell atlas, exclude as well vahedi, not part of baseline
files=($(ls | grep -v cre | grep -v ChM | grep -v cCRE | grep -v ATAC | grep -v atac | grep -v viestra | \
grep -v jarvis | grep -v Vahedi | grep -v 5utr | grep -v hickey |\
grep -v epithelial_peaks | grep -v immune_peaks | grep -v stromal_peaks | grep -v epithelial | grep -v ibd_regions | \
grep -v outside_ibd_regions | grep -v list_all | grep -v list_to_exclude | grep -v baseline | grep -v ENCFF | grep -v ldsc.log)) 

for i in ${files[@]}
do echo ${i}
done

files=($(echo ${files[*]} | sed 's/_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed.gz//g' ))
files=($(echo ${files[*]} | sed 's/allchr_//g'))

# list of annotations
for i in ${files[*]}
do
echo ${i}
done

# length array
echo ${#files[@]}
# 103

for i in ${files[*]}
do
echo ${i} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/allchr_${i}_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed.gz \
| cut -f2  | awk '{s+=$1}END{print s}'
done

for i in ${files[*]}
do
echo ${i} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/allchr_${i}_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed.gz \
| wc -l
done

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed
# 12966826

# concatennate all files to a template - A1 is effect allele for polyfun:
cut -f4 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed \
| awk 'BEGIN {FS = ":";OFS = "\t"} {print $0,$1,$2,$4,$3}' \
| awk 'BEGIN {FS = "\t";OFS = "\t"} {$6="...";print $0}' | awk -F'\t' -v OFS='\t' 'BEGIN {print "SNP","CHR","BP","A1","A2","..."}1' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/merged_files/template_variants_eur_tier2_rate_0.5_info_0.9


wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/merged_files/template_variants_eur_tier2_rate_0.5_info_0.9
# 12966827

wc -l  ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed
# 12966826 


# # subset the annotation files so they only include the 14K variants:

MEM=6000
for i in ${files[*]}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/allchr_${i}_subset_eur_tier2_rate_0.5_info_0.9_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_subset_eur_tier2_rate_0.5_info_0.9_stdout \
"awk 'BEGIN{FS=OFS=\"\\t\"} FNR==NR{arr[\$1];next} ((\$1) in arr)' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/merged_files/template_variants_eur_tier2_rate_0.5_info_0.9 \
<(zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/allchr_${i}_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed.gz) \
| awk -F'\t' -v OFS='\t' 'BEGIN {print \"variant\",\"${i}\"}1' | gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz"
done

for i in ${files[*]}
do echo ${i} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | wc -l
done
# 12966827 OK


# compare order of variants in file:
for i in ${files[*]}
do
echo ${i} && diff <(awk '{print $1}' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/merged_files/template_variants_eur_tier2_rate_0.5_info_0.9) \
<(zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | awk '{print $1}')
done

# 1c1
# < SNP
# ---

# match those that are not aligned
# MEM=8000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
path<-"/path/to/ibdgwas/IIBDGC/"

ref<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/annotation/matrix/merged_files/template_variants_eur_tier2_rate_0.5_info_0.9",sep=""),head=T)

sort_again<-c("missense","frameshift","splice","synonymous")

for (j in c(1:length(sort_again))) {
    
    print(sort_again[j])
    tmp<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_",sort_again[j],"_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz",sep=""),head=T)
    dim(tmp)
    tmp<-tmp[match(ref$SNP,tmp$variant),]
    print(table(tmp$variant==ref$SNP))

    fwrite(tmp,paste(path,"post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_",sort_again[j],"_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed",sep=""),
    col.names=T,row.names=F,quote=F,sep="\t")

    rm(tmp)

}

q("no")

###############################

files_2=("missense" "frameshift" "splice" "synonymous")
for i in ${files_2[*]}
do
echo ${i} && diff <(awk '{print $1}' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/merged_files/template_variants_eur_tier2_rate_0.5_info_0.9) \
<(cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed | awk '{print $1}')
done

for i in ${files_2[*]}
do
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz && \
gzip ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed
done

# compare order of variants in file:
for i in ${files[*]}
do
echo ${i} && diff <(awk '{print $1}' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/merged_files/template_variants_eur_tier2_rate_0.5_info_0.9) \
<(zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | awk '{print $1}')
done


###############################


# paste together
command='paste ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/merged_files/template_variants_eur_tier2_rate_0.5_info_0.9'
for i in ${files[*]}
do
command="$command <(gzip -cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | cut -f2)"
done
eval $command | gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.april2025.annot.gz \
&& unset command

zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.april2025.annot.gz \
| head

zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.april2025.annot.gz \
| wc -l
# 12966827


zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.april2025.annot.gz \
| head -1 | wc -w
# 109
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.april2025.annot.gz \
| tail -1 | wc -w
# 109

# subset by chr
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.april2025.annot.gz \
| awk '{print $0 >> "/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.april2025."$2"_tmp"}' 

# add header
for chr in {1..22} X
do
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.april2025.CHR_tmp \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.april2025.chr${chr}_tmp \
| gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.april2025.chr${chr}.gz
done

# rm intermed file
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.april2025.*_tmp


# #####################

# The .annot file must contain the same SNPs in the same order as the .bim file; 
# sort out file using as template the bim file (some multiallelic varaints causing trouble)

for chr in {1..22}
do
awk -F'\t' '{OFS=FS} NR==FNR {h[$1] = $0; next} {print h[$2]}'  \
<(zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.april2025.chr${chr}.gz) \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/allarrays_chr${chr}_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_rate_0.5_info_0.9.bim \
| cat <(zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.april2025.chr${chr}.gz | head -1) - \
| gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.sorted.april2025.chr${chr}.gz
done

for chr in {1..22}
do
echo ${chr} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.april2025.chr${chr}.gz | wc -l && \
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.sorted.april2025.chr${chr}.gz | wc -l
done

for chr in {1..22}
do
echo ${chr} && diff <(zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.sorted.april2025.chr${chr}.gz | cut -f1) \
<(cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/allarrays_chr${chr}_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_rate_0.5_info_0.9.bim | awk '{print $2}')
done


# create thin annotation files to run ldsc, script uses different annotation format that polyfun:
for chr in {1..22} 
do
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.sorted.april2025.chr${chr}.gz \
| cut -f 7-109 | gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.sorted.april2025.thin.${chr}.annot.gz
done

### LDSC to identify relevant annotation tracks:

# this command generates ld scores (not the ones for polyfun, but the ones for partitioned heritability, but also l2. files), 
# no cm in files, set by kb 
# from: https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial#univariate-ld-scores
# The --ld-wind-cm flag tells lsdc to use a 1 cM window to estimate LD Scores. The other options are --ld-wind-kb, which defines the window size in kilobases, 
# and --ld-wind-snp, which defines the window size in terms of a number of SNPs. We recommend using --ld-wind-cm, because this allows the window size to vary with 
# the range of LD. It is sensible to use a larger window (as measured in kb) in regions like the MHC where LD spans over tens of megabases than in regions with 
# high recombination rate, where LD doesn't extend beyond ~100kb


# used before 1000 KB = 1Mb but computationally took too long. After reviewing LD in EUR (see https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-10-338)
# Figure 1, average r2 +/-Kb is <0.1, and %SNP pairs r2>0.8 is 0

# path_gwas=${path_gwas}

# MEM=58000
# for chr in 2
# MEM=55000
# for chr in 1 {3..6} 
# MEM=45000
# for chr in {7..8}
# MEM=40000
# for chr in {9..12}
# MEM=30000
# for chr in {13..16}
# MEM=25000
# for chr in {17..19}
MEM=20000
for chr in {20..22}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q week \
-e ${path_gwas}post_imputation/log/${chr}_ldsc_eur_tier2_rate_0.5_info_0.9_stderr \
-o ${path_gwas}post_imputation/log/${chr}_ldsc_eur_tier2_rate_0.5_info_0.9_stdout \
"ldsc.py \
--bfile ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/allarrays_chr${chr}_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_rate_0.5_info_0.9 \
--l2 \
--ld-wind-kb 1000 \
--annot ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.sorted.april2025.thin.${chr}.annot.gz \
--thin-annot \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.sorted.april2025.thin.${chr}"
done

# COMPLETED

for chr in {1..22}
do
echo ${chr} && tail -30 ${path_gwas}post_imputation/log/${chr}_ldsc_eur_tier2_rate_0.5_info_0.9_stdout | grep -E "Successfully|Exited|Max Memory"
done

### run the enrichment analyses with the baseline model:

MEM=50000

pheno=(ibd cd uc)
for ph in ${pheno[@]}
do
bsub -J"ldsc_pf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G humgen-priority \
-e ${path_gwas}post_imputation/log/${ph}_ldsc_eur_tier2_rate_0.5_info_0.9_pf_stderr \
-o ${path_gwas}post_imputation/log/${ph}_ldsc_eur_tier2_rate_0.5_info_0.9_pf_stdout \
"ldsc.py \
--h2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/${ph}_eur_tier2_list_variants_rate_0.5_info_0.9_nohla_sumstats_munged.parquet \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/weights.april2025. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.sorted.april2025.thin. \
--overlap-annot \
--thin-annot \
--print-coefficients \
--print-delete-vals \
--not-M-5-50 \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/${ph}.baseline.files_version_april2025"
done

for ph in ${pheno[@]}
do
tail -50 ${path_gwas}post_imputation/log/${ph}_ldsc_eur_tier2_rate_0.5_info_0.9_pf_stdout | grep "Successfully"
done



##################################################################################################################################
##################################################################################################################################
### add the other annotations, one at a time

cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/

# exclude encode, remaining encode, Soskic et al, Zhang et al, and immune cell atlas, exclude as well vahedi, not part of baseline
files_baseline=($(ls | grep -v cre | grep -v ChM | grep -v cCRE | grep -v ATAC | grep -v atac | grep -v viestra | \
grep -v jarvis | grep -v Vahedi | grep -v 5utr | grep -v hickey |\
grep -v epithelial_peaks | grep -v immune_peaks | grep -v stromal_peaks | grep -v epithelial | grep -v ibd_regions | \
grep -v outside_ibd_regions | grep -v list_all | grep -v list_to_exclude | grep -v baseline | grep -v ENCFF)) 


# length array
echo ${#files_baseline[@]}
#103

##########################################
######## rest of annotations:

# viestra_tf - Viestra TF - 1 file
# cCRE_Accessibility - Zhang - 54 files
# ATAC - soskic - 55
# ATAC_counts - Immune Cell atlas - 45
# sc_atac_perianal - Levantovsky - 8 - not run yet
# peaks_b38_list_union - hickey - 180


cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/
files_viestra=($(ls | grep rate_0.5_info_0.9 | grep viestra)) 
echo ${#files_viestra[@]}
# 1

files_zhang=($(ls | grep rate_0.5_info_0.9 | grep cCRE)) 
echo ${#files_zhang[@]}
# 54

files_immcellatlas=($(ls | grep rate_0.5_info_0.9 | grep ATAC_counts)) 
echo ${#files_immcellatlas[@]}
# 45

files_soskic=($(ls | grep rate_0.5_info_0.9 | grep ChM)) 
echo ${#files_soskic[@]}
# 55

files_hickey=($(ls | grep rate_0.5_info_0.9 | grep hickey)) 
echo ${#files_hickey[@]}
# 180

files_encode=($(ls | grep rate_0.5_info_0.9 | grep ENCFF)) 
echo ${#files_encode[@]}
# 251

### all

files=($(ls | grep rate_0.5_info_0.9 | grep -E 'viestra|cCRE|ATAC_counts|ChM|hickey|ENCFF')) 
echo ${#files[@]}
# 586

files=($(echo ${files[*]} | sed 's/_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz//g' ))
files=($(echo ${files[*]} | sed 's/allchr_//g'))
echo ${#files[@]}
# 586

for i in ${files[*]}
do
echo ${i} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz \
| wc -l
done
# 25448334

# compare order of variants in file:
MEM=900

for i in ${files[*]}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/allchr_${i}_subset_check_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_subset_check_stdout \
"diff <(awk '{print \$1}' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/merged_files/template_variants_eur_tier2_rate_0.5_info_0.9) \
<(zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | awk '{print \$1}') \
>> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/var_check_${i}"
done

# all end as exited
for i in ${files[*]}
do
echo ${i} && tail -30 ${path_gwas}post_imputation/log/allchr_${i}_subset_check_stdout | grep Exited
done

for i in ${files[*]}
do
echo ${i} >> ~/tmp.txt && tail ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/var_check_${i} >> ~/tmp.txt 
done
# all OK
rm ~/tmp.txt 

for i in ${files[*]}
do
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/var_check_${i}
done



# paste together 
MEM=300
for i in ${files[*]}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q small \
-e ${path_gwas}post_imputation/log/allchr_${i}_paste_template_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_paste_template_stdout \
"paste ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/merged_files/template_variants_eur_tier2_rate_0.5_info_0.9 \
<(gzip -cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | cut -f2) \
| gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}.april2025.annot.gz"
done

for i in ${files[*]}
do
echo ${i} && tail -30  ${path_gwas}post_imputation/log/allchr_${i}_paste_template_stdout | grep "Successfully"
done

for i in ${files[*]}
do
echo ${i} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}.april2025.annot.gz \
| wc -l
done
# 12966827


# make a folder for each files (to avoid >1000 files per directory):
for i in ${files[*]}
do
mkdir ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/
done

# subset by chr - running for non-encode
MEM=200
for i in ${files[*]}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q small \
-e ${path_gwas}post_imputation/log/allchr_${i}_subset_per_chr_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_subset_per_chr_stdout \
"zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}.april2025.annot.gz \
| awk '{print \$0 >> \"/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.april2025.\"\$2\"_tmp\"}'"
done

for i in ${files[*]}
do
echo ${i} && tail -30 ${path_gwas}post_imputation/log/allchr_${i}_subset_per_chr_stdout | grep "Successfully" 
done

for i in ${files[*]}
do
echo ${i} && ls -la /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.april2025.* | wc -l
done


# # rm intermed file
# for i in ${files[*]}
# do
# rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_subset.bed.gz
# done


# The .annot file must contain the same SNPs in the same order as the .bim file; 
# sort out file using as template the bim file (some multiallelic varaints causing trouble)
# compare with bim file and add header
# create thin annotation files


MEM=1000
for i in ${files[*]}
do
for chr in {1..22}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q small \
-e ${path_gwas}post_imputation/log/allchr_${i}_${chr}_add_header_per_chr_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_${chr}_add_header_per_chr_stdout \
"awk -F'\t' '{OFS=FS} NR==FNR {h[\$1] = \$0; next} {print h[\$2]}' \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.april2025.chr${chr}_tmp \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/allarrays_chr${chr}_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_rate_0.5_info_0.9.bim \
| cat <(cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.april2025.CHR_tmp) - \
| cut -f 7 | gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.${chr}.annot.gz"
done
done


for i in ${files[*]}
do
echo ${i} >> ~/tmp.txt && for chr in {1..22}
do 
echo ${chr} >> ~/tmp.txt && tail -30 ${path_gwas}post_imputation/log/allchr_${i}_${chr}_add_header_per_chr_stdout | grep -E "Successfully|Exit" >> ~/tmp.txt
done
done

cat ~/tmp.txt | grep Exited | wc -l
# 0
cat ~/tmp.txt | grep Successfully | wc -l
# 12892

for i in ${files[*]}
do
ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.1.annot.gz
done

for i in ${files[*]}
do
echo ${i} && ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.*.annot.gz | wc -l
done


for chr in {1..22}
do
for i in ${files[*]}
do
echo ${i} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.${chr}.annot.gz | wc -l
done
done

# # rm intermed file
# for i in ${files[*]}
# do
# rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.april2025.*_tmp
# done

rm ~/tmp.txt

# start submitting zhang the rest (other than ENCODE)
# files=($(ls | grep rate_0.5_info_0.9 | grep -E 'cCRE|viestra|ATAC_counts|ChM|hickey')) 
# echo ${#files[@]}
# # 335

# files=($(ls | grep rate_0.5_info_0.9 | grep -E 'ENCFF')) 
# echo ${#files[@]}
# # 251


cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/
files=($(ls | grep rate_0.5_info_0.9 | grep -E 'cCRE|viestra|ATAC_counts|ChM|hickey|ENCFF')) 
echo ${#files[@]}

files=($(echo ${files[*]} | sed 's/_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz//g' ))
files=($(echo ${files[*]} | sed 's/allchr_//g'))
echo ${#files[@]}
#586

for i in ${files[*]}
do
# MEM=60000 && for chr in 2
# MEM=55000 && for chr in 1 {3..6} 
# MEM=45000 && for chr in {7..8}
# MEM=40000 && for chr in {9..12}
# MEM=30000 && for chr in {13..16}
MEM=25000 && for chr in {17..19}
# MEM=20000 && for chr in {20..22}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G humgen-priority -q week \
-e ${path_gwas}post_imputation/log/${chr}_${i}_ldsc_stderr \
-o ${path_gwas}post_imputation/log/${chr}_${i}_ldsc_stdout \
"ldsc.py \
--bfile ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/allarrays_chr${chr}_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_rate_0.5_info_0.9 \
--l2 \
--ld-wind-kb 1000 \
--annot ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.${chr}.annot.gz \
--thin-annot \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.${chr}"
done
done



rm ~/tmp_list.txt
for i in ${files[*]}
do
echo ${i} >> ~/tmp_list.txt && for chr in {1..22}
do 
echo ${chr} >> ~/tmp_list.txt && tail -30 ${path_gwas}post_imputation/log/${chr}_${i}_ldsc_stdout | grep -E "Successfully|Exited" >> ~/tmp_list.txt
done
done



# for i in ${files[*]}
# do
# echo ${i} && for chr in {1..22} X
# do
# echo ${chr} && tail -50 ${path_gwas}post_imputation/log/${chr}_${i}_ldsc_stdout | grep -E "Successfully|Exit"
# done
# done


for i in ${files[*]}
do
echo ${i} && ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.*.l2.M | wc -l
done


rm  ~/tmp
for i in ${files[*]}
do
echo ${i} >> ~/tmp 
for chr in {1..22}
do
echo ${chr} >> ~/tmp && \
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.${chr}.log | grep 'Analysis finished' >> ~/tmp
done
done




# resubmit manually files not completed:
# i=non_multiome_stromal_Crypt_Fibroblasts_2_hickey

# for chr in 6
# do

# if [[ "${chr}" == "2" ]] ; then
# MEM=60000
# elif [[ "${chr}" == "1" ]] ; then
# MEM=55000
# elif [[ "${chr}" -ge 3 && "${chr}" -le 6 ]] ; then
# MEM=55000
# elif [[ "${chr}" -ge 7 && "${chr}" -le 8 ]] ; then
# MEM=45000
# elif [[ "${chr}" -ge 9 && "${chr}" -le 12 ]] ; then
# MEM=40000
# elif [[ "${chr}" -ge 13 && "${chr}" -le 16 ]] ; then
# MEM=30000
# elif [[ "${chr}" -ge 17 && "${chr}" -le 19 ]] ; then
# MEM=25000
# elif [[ "${chr}" -ge 20 && "${chr}" -le 22 ]] ; then
# MEM=20000
# fi


# bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G humgen-priority -q week \
# -e ${path_gwas}post_imputation/log/${chr}_${i}_ldsc_stderr \
# -o ${path_gwas}post_imputation/log/${chr}_${i}_ldsc_stdout \
# "ldsc.py \
# --bfile ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/allarrays_chr${chr}_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_rate_0.5_info_0.9 \
# --l2 \
# --ld-wind-kb 1000 \
# --annot ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.${chr}.annot.gz \
# --thin-annot \
# --out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.${chr}"

# done

# duodenum_epithelial_Best4_Enterocytes_hickey - chr6 resubmitted on 16 J





# # submit the ldsc jobs

MEM=73000
pheno=(ibd cd uc)

for i in ${files[*]}
do
for ph in ${pheno[@]}
do
bsub -J"ldsc_pf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G humgen-priority -q long \
-e ${path_gwas}post_imputation/log/${ph}_${i}_ldsc_pf_stderr \
-o ${path_gwas}post_imputation/log/${ph}_${i}_ldsc_pf_stdout \
"ldsc.py \
--h2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/${ph}_eur_tier2_list_variants_rate_0.5_info_0.9_nohla_sumstats_munged.parquet \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/weights.april2025. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.sorted.april2025.thin.,${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin. \
--overlap-annot \
--thin-annot \
--print-coefficients \
--print-delete-vals \
--not-M-5-50 \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/${ph}.baseline.${i}.files_version_april2025"
done
done

/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/uc.baseline.ATAC_counts_Memory_Tregs-S_merged_samples.files_version_april2025


rm ~/tmp
for ph in ${pheno[@]}
do
for i in ${files[*]}
do
echo ${i} >> ~/tmp && tail -50 ${path_gwas}post_imputation/log/${ph}_${i}_ldsc_pf_stdout | grep -E "Successfully|Exit" >> ~/tmp && \
echo ${ph} >> ~/tmp
done
done



# all completed

# lfs quota -hg ibdgwas /path/to/project

# ###############################################

# # plot enrichment analysis results

# ~/git/IIBDGC_GWAS/scripts/other/plot_enrichment_analyses.R


# # cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/
# cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/

# # all but encode
# files=($(ls | grep -E 'sc_atac_perianal' | grep -v non_multiome  | grep -v jarvis | grep -v allchr_ATAC_counts_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis.bed.gz | grep -v allchr_ATAC_counts_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis.bed.gz | grep -v .april2025.annot.gz | uniq)) 

# # submitted for Kyle to check out
# for i in ${files[*]}
# do
# cp ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.chr*.gz ${path_gwas}to_upload_to_iibdgc_globus/
# done

# ##################################################################################################################################################################
# ##################################################################################################################################################################

### run the same analysis with a negative control - height
# see reformat_other_gwas_summary_statistics_to_munge.R

echo ${#files[@]}
# 586


MEM=73000

for i in ${files[*]}
do
bsub -J"ldsc_pf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G humgen-priority \
-e ${path_gwas}post_imputation/log/height_${i}_ldsc_pf_height_stderr \
-o ${path_gwas}post_imputation/log/height_${i}_ldsc_pf_height_stdout \
"ldsc.py \
--h2 ${path_gwas}resources/gwas_summary_statistics/sakaue_2021_height_list_variants_maf_0.001_munged.parquet \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/weights.april2025. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.sorted.april2025.thin.,${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin. \
--overlap-annot \
--thin-annot \
--print-coefficients \
--print-delete-vals \
--not-M-5-50 \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/other_traits/height.baseline.${i}.files_version_april2025"
done


# submitted
for i in ${files[*]}
do
echo ${i} && tail -50 ${path_gwas}post_imputation/log/height_${i}_ldsc_pf_height_stdout  | grep -E "Successfully|Exit"
done


# COMPLETED
