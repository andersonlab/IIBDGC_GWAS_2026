# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# ESTIMATE ENRICHMENT IN ANNOTATIONS:

# start getting ld tags from our in sample data - use r=0.5 as in Susie purity -> r2 0.25

# singularity exec iibdgc_postprocess_10_singularity.sif


path_gwas="/path/to/ibdgwas/IIBDGC/"

eqtl_study=pQTL_sparc

mkdir ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/
mkdir -p ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/

path="/path/to/project"

dataset=plasma_ibd_patients

MEM=4500

for chr in {1..22}
do
bsub -J"ld" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${dataset}_${chr}_ld_stdout \
-e ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${dataset}_${chr}_ld_stderr \
"/path/to/software/username/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_ibd_analysis \
--show-tags <(zcat /path/to/project | cut -f8 | sort | uniq | sed 's/_/:/g' | sed 's/^/chr/g') \
--tag-r2 0.25 \
--tag-kb 2000 \
--out ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${dataset}_chr${chr}_list_tier_2_ld_r2_0.25"
done

for chr in {1..22}
do
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${dataset}_${chr}_ld_stdout | grep Successfully
done


echo ${dataset} && ls -la ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${dataset}_chr*_list_tier_2_ld_r2_0.25.tags | wc -l


# concatenate all files

MEM=200

bsub -J"ld" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${dataset}_allchr_ld_stdout \
-e ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${dataset}_allchr_ld_stderr \
"cat ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${dataset}_chr*_list_tier_2_ld_r2_0.25.tags  > \
${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${dataset}_allchr_list_tier_2_ld_r2_0.25.tags"
 
echo ${dataset} && tail -50 ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${dataset}_allchr_ld_stdout | grep Successfully

echo ${dataset} && less ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${dataset}_allchr_list_tier_2_ld_r2_0.25.tags | wc -l



#  make a bed file per condition

MEM=200
bsub -J"ld" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${dataset}_allchr_bed1_ld_stdout \
-e ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${dataset}_allchr_bed1_ld_stderr \
"awk 'BEGIN {FS = \":\";OFS = \"\t\"} {print \$1,\$2,\$2+1}' ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${dataset}_allchr_list_tier_2_ld_r2_0.25.tags > ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${dataset}_allchr_list_tier_2_ld_r2_0.25_tmp"


# double check files:
echo ${dataset} && \
wc -l ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${dataset}_allchr_list_tier_2_ld_r2_0.25_tmp && \
wc -l ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${dataset}_allchr_list_tier_2_ld_r2_0.25.tags


MEM=200
bsub -J"ld" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${dataset}_allchr_bed2_ld_stdout \
-e ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${dataset}_allchr_bed2_ld_stderr \
"paste -d'\t' ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${dataset}_allchr_list_tier_2_ld_r2_0.25_tmp \
${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${dataset}_allchr_list_tier_2_ld_r2_0.25.tags \
> ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${dataset}_allchr_list_tier_2_ld_r2_0.25.bed"

echo ${dataset} && \
head ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${dataset}_allchr_list_tier_2_ld_r2_0.25.bed && \
tail ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${dataset}_allchr_list_tier_2_ld_r2_0.25.bed


# remove intermediate fiels
rm ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_*_allchr_list_tier_2_ld_r2_0.25_tmp
rm ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_*_list_tier_2_ld_r2_0.25.tags

# 'V' sorts alphanumerically
MEM=200

bsub -J"ld" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${dataset}_allchr_bed3_ld_stdout \
-e ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${dataset}_allchr_bed3_ld_stderr \
"sort -k1,1V -k2,2n ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${dataset}_allchr_list_tier_2_ld_r2_0.25.bed \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/${eqtl_study}_${dataset}_allchr_list_tier_2_ld_r2_0.25_sorted.bed"


echo ${dataset} && tail -50 ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${dataset}_allchr_bed3_ld_stdout | grep Successfully


# remove intermed files
rm ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_*_allchr_list_tier_2_ld_r2_0.25.bed

cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/
ls -la *_sorted.bed | wc -l
# 1

files_eqtl=($(ls | grep "_sorted.bed" | sed 's/_sorted\.bed//g' ))


# intersect each file with list of analysed variants, and add header

MEM=1000
for i in ${files_eqtl[@]}
do
bsub -J"eqlt" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${i}_eqtl_3_stderr \
-o ${path_gwas}post_imputation/log/${i}_eqtl_3_stdout \
"bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed \
-b ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/${i}_sorted.bed -c \
| awk -F'\t' 'BEGIN {OFS=FS} {if (\$5==\"2\") \$5=\"1\"; print \$0}' | awk -F'\t' 'BEGIN {OFS=FS} {print \$4,\$5}' | awk -F'\t' -v OFS='\t' 'BEGIN {print \"variant\",\"'\"${i}\"'\"}1' \
| sed 's/condition/\${i}/g' | gzip > \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz"
done



for i in ${files_eqtl[@]}
do 
echo ${i} && tail -50 ${path_gwas}post_imputation/log/${i}_eqtl_3_stdout | grep "Successfully"
done


# check out the length of the files:
for i in ${files_eqtl[@]}
do 
echo ${i} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | wc -l
done
# 12966827 - OK

# double check values are 0|1, and count:
for i in ${files_eqtl[@]}
do
echo ${i} && \
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | cut -f2 | sort | uniq -c
done

# pQTL_sparc_plasma_ibd_patients_allchr_list_tier_2_ld_r2_0.25
# 12910789 0
#   56037 1
#       1 pQTL_sparc_plasma_ibd_patients_allchr_list_tier_2_ld_r2_0.25


# remove tmp and other files
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/*_tmp.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/*_sorted.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/*_sorted.bed.gz


# mv data to final folder
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/allchr_*_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/




