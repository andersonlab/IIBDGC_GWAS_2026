# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# ESTIMATE ENRICHMENT IN ANNOTATIONS:

# start getting ld tags from our in sample data - use r=0.5 as in Susie purity -> r2 0.25

# singularity exec iibdgc_postprocess_10_singularity.sif


path_gwas="/path/to/ibdgwas/IIBDGC/"

eqtl_study=macromap
condition=CIL_24

mkdir ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/
mkdir ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/

conditions=(IFNB_6 IFNB_24 P3C_6 P3C_24 CIL_6 CIL_24 IFNG_6 IFNG_24 PIC_6 PIC_24 IL4_6 IL4_24 Prec_D0 Prec_D2 Ctrl_6 Ctrl_24 LIL10_6 LIL10_24 R848_6 R848_24 MBP_6 MBP_24 sLPS_6 sLPS_24)

MEM=4500

for condition in ${conditions[@]}
do 
for chr in {1..22}
do
bsub -J"ld" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_${chr}_ld_stdout \
-e ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_${chr}_ld_stderr \
"/path/to/software/username/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_ibd_analysis \
--show-tags <(zcat /path/to/project | cut -f8 | sort | uniq | sed 's/_/:/g') \
--tag-r2 0.25 \
--tag-kb 2000 \
--out ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_chr${chr}_list_tier_2_ld_r2_0.25"
done
done


for condition in ${conditions[@]}
do 
echo ${condition} && for chr in {1..22}
do
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_${chr}_ld_stdout | grep Successfully
done
done

for condition in ${conditions[@]}
do 
echo ${condition} && ls -la ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_chr*_list_tier_2_ld_r2_0.25.tags | wc -l
done

# concatenate all files

MEM=200

for condition in ${conditions[@]}
do 
bsub -J"ld" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_ld_stdout \
-e ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_ld_stderr \
"cat ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_chr*_list_tier_2_ld_r2_0.25.tags  > \
${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_0.25.tags"
done

for condition in ${conditions[@]}
do 
echo ${condition} && tail -50 ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_ld_stdout | grep Successfully
done
done

for condition in ${conditions[@]}
do 
echo ${condition} && less ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_0.25.tags | wc -l
done


#  make a bed file per condition

MEM=200
for condition in ${conditions[@]}
do 
bsub -J"ld" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_bed1_ld_stdout \
-e ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_bed1_ld_stderr \
"awk 'BEGIN {FS = \":\";OFS = \"\t\"} {print \$1,\$2,\$2+1}' ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_0.25.tags > ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_0.25_tmp"
done

# double check files:
for condition in ${conditions[@]}
do echo ${condition} && \
wc -l ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_0.25_tmp && \
wc -l ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_0.25.tags
done


MEM=200
for condition in ${conditions[@]}
do 
bsub -J"ld" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_bed2_ld_stdout \
-e ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_bed2_ld_stderr \
"paste -d'\t' ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_0.25_tmp \
${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_0.25.tags \
> ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_0.25.bed"
done

for condition in ${conditions[@]}
do echo ${condition} && \
head ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_0.25.bed && \
tail ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_0.25.bed
done

# remove intermediate fiels
rm ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_*_allchr_list_tier_2_ld_r2_0.25_tmp
rm ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_*_list_tier_2_ld_r2_0.25.tags

# 'V' sorts alphanumerically
MEM=200
for condition in ${conditions[@]}
do 
bsub -J"ld" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_bed3_ld_stdout \
-e ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_bed3_ld_stderr \
"sort -k1,1V -k2,2n ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_0.25.bed \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38//${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_0.25_sorted.bed"
done

for condition in ${conditions[@]}
do 
echo ${condition} && tail -50 ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_bed3_ld_stdout | grep Successfully
done

# remove intermed files
rm ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_*_allchr_list_tier_2_ld_r2_0.25.bed

cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/
ls -la *_sorted.bed | wc -l
# 24

files_eqtl=($(ls | grep "_sorted.bed" | sed 's/_sorted\.bed//g' ))


# intersect each file with list of analysed variants, and add header

MEM=1000
for i in ${files_eqtl[@]}
do
bsub -J"macromap" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
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

# macromap_CIL_24_allchr_list_tier_2_ld_r2_0.25
# 12144201 0
#  822625 1
#       1 macromap_CIL_24_allchr_list_tier_2_ld_r2_0.25
# macromap_CIL_6_allchr_list_tier_2_ld_r2_0.25
# 12143608 0
#  823218 1
#       1 macromap_CIL_6_allchr_list_tier_2_ld_r2_0.25
# macromap_Ctrl_24_allchr_list_tier_2_ld_r2_0.25
# 12125230 0
#  841596 1
#       1 macromap_Ctrl_24_allchr_list_tier_2_ld_r2_0.25
# macromap_Ctrl_6_allchr_list_tier_2_ld_r2_0.25
# 12129871 0
#  836955 1
#       1 macromap_Ctrl_6_allchr_list_tier_2_ld_r2_0.25
# macromap_IFNB_24_allchr_list_tier_2_ld_r2_0.25
# 12130145 0
#  836681 1
#       1 macromap_IFNB_24_allchr_list_tier_2_ld_r2_0.25
# macromap_IFNB_6_allchr_list_tier_2_ld_r2_0.25
# 12131783 0
#  835043 1
#       1 macromap_IFNB_6_allchr_list_tier_2_ld_r2_0.25
# macromap_IFNG_24_allchr_list_tier_2_ld_r2_0.25
# 12135292 0
#  831534 1
#       1 macromap_IFNG_24_allchr_list_tier_2_ld_r2_0.25
# macromap_IFNG_6_allchr_list_tier_2_ld_r2_0.25
# 12128637 0
#  838189 1
#       1 macromap_IFNG_6_allchr_list_tier_2_ld_r2_0.25
# macromap_IL4_24_allchr_list_tier_2_ld_r2_0.25
# 12115190 0
#  851636 1
#       1 macromap_IL4_24_allchr_list_tier_2_ld_r2_0.25
# macromap_IL4_6_allchr_list_tier_2_ld_r2_0.25
# 12126728 0
#  840098 1
#       1 macromap_IL4_6_allchr_list_tier_2_ld_r2_0.25
# macromap_LIL10_24_allchr_list_tier_2_ld_r2_0.25
# 12140914 0
#  825912 1
#       1 macromap_LIL10_24_allchr_list_tier_2_ld_r2_0.25
# macromap_LIL10_6_allchr_list_tier_2_ld_r2_0.25
# 12138573 0
#  828253 1
#       1 macromap_LIL10_6_allchr_list_tier_2_ld_r2_0.25
# macromap_MBP_24_allchr_list_tier_2_ld_r2_0.25
# 12152916 0
#  813910 1
#       1 macromap_MBP_24_allchr_list_tier_2_ld_r2_0.25
# macromap_MBP_6_allchr_list_tier_2_ld_r2_0.25
# 12140540 0
#  826286 1
#       1 macromap_MBP_6_allchr_list_tier_2_ld_r2_0.25
# macromap_P3C_24_allchr_list_tier_2_ld_r2_0.25
# 12127863 0
#  838963 1
#       1 macromap_P3C_24_allchr_list_tier_2_ld_r2_0.25
# macromap_P3C_6_allchr_list_tier_2_ld_r2_0.25
# 12151095 0
#  815731 1
#       1 macromap_P3C_6_allchr_list_tier_2_ld_r2_0.25
# macromap_PIC_24_allchr_list_tier_2_ld_r2_0.25
# 12147801 0
#  819025 1
#       1 macromap_PIC_24_allchr_list_tier_2_ld_r2_0.25
# macromap_PIC_6_allchr_list_tier_2_ld_r2_0.25
# 12156807 0
#  810019 1
#       1 macromap_PIC_6_allchr_list_tier_2_ld_r2_0.25
# macromap_Prec_D0_allchr_list_tier_2_ld_r2_0.25
# 12178409 0
#  788417 1
#       1 macromap_Prec_D0_allchr_list_tier_2_ld_r2_0.25
# macromap_Prec_D2_allchr_list_tier_2_ld_r2_0.25
# 12115537 0
#  851289 1
#       1 macromap_Prec_D2_allchr_list_tier_2_ld_r2_0.25
# macromap_R848_24_allchr_list_tier_2_ld_r2_0.25
# 12126929 0
#  839897 1
#       1 macromap_R848_24_allchr_list_tier_2_ld_r2_0.25
# macromap_R848_6_allchr_list_tier_2_ld_r2_0.25
# 12125144 0
#  841682 1
#       1 macromap_R848_6_allchr_list_tier_2_ld_r2_0.25
# macromap_sLPS_24_allchr_list_tier_2_ld_r2_0.25
# 12123807 0
#  843019 1
#       1 macromap_sLPS_24_allchr_list_tier_2_ld_r2_0.25
# macromap_sLPS_6_allchr_list_tier_2_ld_r2_0.25
# 12155866 0
#  810960 1
#       1 macromap_sLPS_6_allchr_list_tier_2_ld_r2_0.25


# remove tmp and other files
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/*_tmp.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/*_sorted.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/*_sorted.bed.gz


# mv data to final folder
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/allchr_*_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/




