# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# recode_encode_atac_dnase

# singularity exec iibdgc_postprocess_10_singularity.sif

# ENCODE ANNOTATION - sCRE - recode by type:
path_gwas="/path/to/ibdgwas/IIBDGC/"


####################
# DNASE ASSAYS:

files_path=/path/to/project
cd ${files_path}

files_dnase=($(ls | grep bed | sed 's/\.bed\.gz//g' ))

for i in ${files_dnase[@]}
do echo ${i}
done

echo ${#files_dnase[@]}
# 183


MEM=500

# sort files
for i in ${files_dnase[@]}
do 
bsub -J"dnase_atac" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${i}_encode_dnase_atac_1_stderr \
-o ${path_gwas}post_imputation/log/${i}_encode_dnase_atac_1_stdout \
"zcat ${files_path}/${i}.bed.gz \
| sort -u \
| sort -k1,1V -k2,2n \
| gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/dnase_atacseq/bed_b38/${i}_sorted.bed.gz"
done

for i in ${files_dnase[@]}
do 
echo ${i} && tail -50 ${path_gwas}post_imputation/log/${i}_encode_dnase_atac_1_stdout | grep "Successfully"
done


####################
# ATACSEQ ASSAYS:

files_path=/path/to/project
cd ${files_path}

files_atac=($(ls | grep bed | sed 's/\.bed\.gz//g' ))

for i in ${files_atac[@]}
do echo ${i}
done

echo ${#files_atac[@]}
# 46

MEM=500

# sort files
for i in ${files_atac[@]}
do 
bsub -J"dnase_atac" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${i}_encode_dnase_atac_1_stderr \
-o ${path_gwas}post_imputation/log/${i}_encode_dnase_atac_1_stdout \
"zcat ${files_path}/${i}.bed.gz \
| sort -u \
| sort -k1,1V -k2,2n \
| gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/dnase_atacseq/bed_b38/${i}_sorted.bed.gz"
done


for i in ${files_atac[@]}
do 
echo ${i} && tail -50 ${path_gwas}post_imputation/log/${i}_encode_dnase_atac_1_stdout | grep "Successfully"
done


#############################################################################


cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/dnase_atacseq/bed_b38/
ls -la *_sorted.bed.gz | wc -l
# 229


files_dnase_atac=($(ls | grep "_sorted.bed" | sed 's/_sorted\.bed\.gz//g' ))

for i in ${files_dnase_atac[@]}
do echo ${i}
done

echo ${#files_dnase_atac[@]}
# 229

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed
# 12966826


# duplicated IDs, merge regions that overlap, estimate mean value and run:
MEM=200
for i in ${files_dnase_atac[@]}
do
bsub -J"dnase_atac" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 2 \
-e ${path_gwas}post_imputation/log/${i}_encode_dnase_atac_2_stderr \
-o ${path_gwas}post_imputation/log/${i}_encode_dnase_atac_2_stdout \
"bedtools merge -i <(zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/dnase_atacseq/bed_b38/${i}_sorted.bed.gz) -c 7 \
-o mean > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/dnase_atacseq/bed_b38/${i}_sorted_nooverlaps.bed"
done


for i in ${files_atac[@]}
do 
echo ${i} && tail -50 ${path_gwas}post_imputation/log/${i}_encode_dnase_atac_2_stdout | grep "Successfully"
done

for i in ${files_atac[@]}
do 
echo ${i} && ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/dnase_atacseq/bed_b38/${i}_sorted_nooverlaps.bed
done


# intersect each file with list of analysed variants, and add header
MEM=250

for i in ${files_dnase_atac[@]}
do
bsub -J"dnase_atac" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${i}_encode_dnase_atac_3_stderr \
-o ${path_gwas}post_imputation/log/${i}_encode_dnase_atac_3_stdout \
"bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed \
-b ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/dnase_atacseq/bed_b38/${i}_sorted_nooverlaps.bed \
-loj | awk -F'\t' 'BEGIN {OFS=FS} {print \$4,\$8}' | awk -F'\t' -v OFS='\t' 'BEGIN {print \"variant\",\"'\"${i}\"'\"}1' | \
sed 's/\\t\./\\t0/gm' | gzip > \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/dnase_atacseq/bed_b38/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz"
done

for i in ${files_dnase_atac[@]}
do 
echo ${i} && tail -50 ${path_gwas}post_imputation/log/${i}_encode_dnase_atac_3_stdout | grep "Successfully"
done


# check out the length of the files:
for i in ${files_dnase_atac[@]}
do 
echo ${i} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/dnase_atacseq/bed_b38/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | wc -l
done

# double check min value is always positive:
for i in ${files_dnase_atac[@]}
do
echo ${i} && \
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/dnase_atacseq/bed_b38/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | cut -f2 | sort -n | sed -n '1s/^/min=/p; $s/^/max=/p'
done

# remove tmp and other files
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/dnase_atacseq/bed_b38/*_tmp.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/dnase_atacseq/bed_b38/*_sorted.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/dnase_atacseq/bed_b38/*_sorted.bed.gz
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/dnase_atacseq/bed_b38/*_sorted_nooverlaps.bed


# mv data to final folder
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/dnase_atacseq/bed_b38/allchr_*_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/


