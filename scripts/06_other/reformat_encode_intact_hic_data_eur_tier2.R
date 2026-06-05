# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# recode_encode_ccres_data

# singularity exec iibdgc_postprocess_10_singularity.sif

# # ENCODE ANNOTATION - sCRE - recode by type:
path_gwas="/path/to/ibdgwas/IIBDGC/"

files_path=/path/to/project
cd ${files_path}


files_intact_hic=($(ls | grep bedpe.gz | sed 's/\.bedpe\.gz//g' ))

for i in ${files_intact_hic[@]}
do echo ${i}
done

echo ${#files_intact_hic[@]}
# 44 - loop and interaction files, keep only the loop files

MEM=1200

for i in ${files_intact_hic[@]}
do
bsub -J"hic" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${i}_encode_intact_hic_stderr \
-o ${path_gwas}post_imputation/log/${i}_encode_intact_hic_stdout \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/split_bedpe_intact_hic_annotation_per_type.R ${i} > \
${path_gwas}post_imputation/2022/log/split_bedpe_intact_hic_annotation_per_type${i}.Rout"
done

for i in ${files_intact_hic[@]}
do 
echo ${i} && tail -50 ${path_gwas}post_imputation/log/${i}_encode_intact_hic_stdout | grep "Successfully"
done

for i in ${files_intact_hic[@]}
do 
echo ${i} && ls -la ${path_gwas}/post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/${i}_intact_hic.bed.gz
done

ls -la ${path_gwas}/post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/*_intact_hic.bed.gz | wc -l
#  22

cd ${path_gwas}/post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/
files_intact_hic=($(ls | grep "_intact_hic.bed.gz" | sed 's/\.bed\.gz//g' ))

for i in ${files_intact_hic[@]}
do echo ${i}
done

# this list now only includes the loops
echo ${#files_intact_hic[@]}
# 22

MEM=200

# sort files
for i in ${files_intact_hic[@]}
do 
bsub -J"hic" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${i}_encode_intact_hic_1_stderr \
-o ${path_gwas}post_imputation/log/${i}_encode_intact_hic_1_stdout \
"zcat ${path_gwas}/post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/${i}.bed.gz \
| sort -u \
| sort -k1,1V -k2,2n \
| gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/${i}_sorted.bed.gz"
done


for i in ${files_intact_hic[@]}
do 
echo ${i} && tail -30 ${path_gwas}post_imputation/log/${i}_encode_intact_hic_1_stdout | grep "Successfully"
done


for i in ${files_intact_hic[@]}
do 
echo ${i} && ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/${i}_sorted.bed.gz
done


# in case of overlapping regions, merge regions that overlap, estimate mean value and run:

for i in ${files_intact_hic[@]}
do
bsub -J"dnase_atac" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 2 \
-e ${path_gwas}post_imputation/log/${i}_encode_intact_hic_2_stderr \
-o ${path_gwas}post_imputation/log/${i}_encode_intact_hic_2_stdout \
"bedtools merge -i <(zcat ${path_gwas}/post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/${i}_sorted.bed.gz) \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/${i}_sorted_nooverlaps.bed"
done

for i in ${files_intact_hic[@]}
do 
echo ${i} && tail -50 ${path_gwas}post_imputation/log/${i}_encode_intact_hic_2_stdout | grep "Successfully"
done

ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/*_sorted_nooverlaps.bed | wc -l
# 22

# intersect each file with list of analysed variants, replace . by 0, and add header
MEM=850

for i in ${files_intact_hic[@]}
do
bsub -J"hic" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 2 \
-e ${path_gwas}post_imputation/log/${i}_encode_intact_hic_3_stderr \
-o ${path_gwas}post_imputation/log/${i}_encode_intact_hic_3_stdout \
"bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed \
-b ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/${i}_sorted_nooverlaps.bed -c \
| awk -F'\t' 'BEGIN {OFS=FS} {if (\$5==\"2\") \$5=\"1\"; print \$0}' | awk -F'\t' 'BEGIN {OFS=FS} {print \$4,\$5}' | awk -F'\t' -v OFS='\t' 'BEGIN {print \"variant\",\"'\"${i}\"'\"}1' \
| sed 's/\\t\./\\t0/gm' | gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz"
done

for i in ${files_intact_hic[@]}
do 
echo ${i} && tail -50 ${path_gwas}post_imputation/log/${i}_encode_intact_hic_3_stdout | grep "Successfully"
done

# check out the length of the files:
for i in ${files_intact_hic[@]}
do
echo ${i} && \
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | wc -l
done


# double check min value is always positive:
for i in ${files_intact_hic[@]}
do
echo ${i} && \
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | cut -f2  | sort -n |  uniq -c
done


# remove tmp and other files
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/*_tmp.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/*_sorted.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/*_sorted_nooverlaps.bed


# mv data to final folder
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/allchr_*_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/


