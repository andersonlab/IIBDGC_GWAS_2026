# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# ESTIMATE ENRICHMENT IN ANNOTATIONS:

# start getting ld tags from our in sample data - use r=0.5 as in Susie purity -> r2 0.25
# test here additional thresholds at 0.5 and ${ld}

# singularity exec iibdgc_postprocess_10_singularity.sif

path_gwas="/path/to/ibdgwas/IIBDGC/"

eqtl_study=macromap

mkdir ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/
mkdir ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/

conditions=(IFNB_6 IFNB_24 P3C_6 P3C_24 CIL_6 CIL_24 IFNG_6 IFNG_24 PIC_6 PIC_24 IL4_6 IL4_24 Prec_D0 Prec_D2 Ctrl_6 Ctrl_24 LIL10_6 LIL10_24 R848_6 R848_24 MBP_6 MBP_24 sLPS_6 sLPS_24)

MEM=4500

for condition in ${conditions[@]}
do 
for chr in {1..22}
do
for ld in 0.5 0.75
do
bsub -J"ld" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_${chr}_${ld}_ld_stdout \
-e ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_${chr}_${ld}_ld_stderr \
"/path/to/software/username/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_ibd_analysis \
--show-tags <(zcat /path/to/project | cut -f8 | sort | uniq | sed 's/_/:/g') \
--tag-r2 ${ld} \
--tag-kb 2000 \
--out ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_chr${chr}_list_tier_2_ld_r2_${ld}"
done
done
done

# SUBMITTED - both

for condition in ${conditions[@]}
do 
echo ${condition} && for chr in {1..22}
do
for ld in 0.5 0.75
do
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_${chr}_${ld}_ld_stdout | grep Successfully
done
done
done

for condition in ${conditions[@]}
do 
for ld in 0.5 0.75
do
echo ${condition} && ls -la ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_chr*_list_tier_2_ld_r2_${ld}.tags | wc -l
done
done

# concatenate all files

MEM=200

for condition in ${conditions[@]}
do 
for ld in 0.5 0.75
do
bsub -J"ld" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_${ld}_ld_stdout \
-e ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_${ld}_ld_stderr \
"cat ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_chr*_list_tier_2_ld_r2_${ld}.tags  > \
${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_${ld}.tags"
done
done


for condition in ${conditions[@]}
do 
for ld in 0.5 0.75
do
echo ${condition} && tail -50 ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_${ld}_ld_stdout | grep Successfully
done
done

for condition in ${conditions[@]}
do 
for ld in 0.5 0.75
do
echo ${condition} && less ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_${ld}.tags | wc -l
done
done

#  make a bed file per condition

MEM=200
for condition in ${conditions[@]}
do 
for ld in 0.5 0.75
do
bsub -J"ld" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_bed1_${ld}_ld_stdout \
-e ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_bed1_${ld}_ld_stderr \
"awk 'BEGIN {FS = \":\";OFS = \"\t\"} {print \$1,\$2,\$2+1}' ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_${ld}.tags > ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_${ld}_tmp"
done
done

# double check files:
for condition in ${conditions[@]}
do echo ${condition} && \
for ld in 0.5 0.75
do
wc -l ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_${ld}_tmp && \
wc -l ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_${ld}.tags
done
done

MEM=200
for condition in ${conditions[@]}
do 
for ld in 0.5 0.75
do
bsub -J"ld" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_bed2_${ld}_ld_stdout \
-e ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_bed2_${ld}_ld_stderr \
"paste -d'\t' ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_${ld}_tmp \
${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_${ld}.tags \
> ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_${ld}.bed"
done
done

for condition in ${conditions[@]}
do echo ${condition} && \
for ld in 0.5 0.75
do
head ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_${ld}.bed && \
tail ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_${ld}.bed
done
done

# remove intermediate fields
rm ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_*_allchr_list_tier_2_ld_r2_${ld}_tmp
rm ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_*_list_tier_2_ld_r2_${ld}.tags

# 'V' sorts alphanumerically
MEM=200
for condition in ${conditions[@]}
do 
for ld in 0.5 0.75
do
bsub -J"ld" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_bed3_${ld}_ld_stdout \
-e ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_bed3_${ld}_ld_stderr \
"sort -k1,1V -k2,2n ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_${ld}.bed \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38//${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_${ld}_sorted.bed"
done
done


for condition in ${conditions[@]}
do 
for ld in 0.5 0.75
do
echo ${condition} && tail -50 ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_bed3_${ld}_ld_stdout | grep Successfully
done
done

# remove intermed files
for ld in 0.5 0.75
do
rm ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_*_allchr_list_tier_2_ld_r2_${ld}.bed
done


cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/
ls -la *_sorted.bed | wc -l
# 48

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

# macromap_CIL_24_allchr_list_tier_2_ld_r2_0.5
# 12570117 0
#  396709 1
#       1 macromap_CIL_24_allchr_list_tier_2_ld_r2_0.5
# macromap_CIL_24_allchr_list_tier_2_ld_r2_0.75
# 12768338 0
#  198488 1
#       1 macromap_CIL_24_allchr_list_tier_2_ld_r2_0.75
# macromap_CIL_6_allchr_list_tier_2_ld_r2_0.5
# 12577361 0
#  389465 1
#       1 macromap_CIL_6_allchr_list_tier_2_ld_r2_0.5
# macromap_CIL_6_allchr_list_tier_2_ld_r2_0.75
# 12775820 0
#  191006 1
#       1 macromap_CIL_6_allchr_list_tier_2_ld_r2_0.75
# macromap_Ctrl_24_allchr_list_tier_2_ld_r2_0.5
# 12562891 0
#  403935 1
#       1 macromap_Ctrl_24_allchr_list_tier_2_ld_r2_0.5
# macromap_Ctrl_24_allchr_list_tier_2_ld_r2_0.75
# 12761563 0
#  205263 1
#       1 macromap_Ctrl_24_allchr_list_tier_2_ld_r2_0.75
# macromap_Ctrl_6_allchr_list_tier_2_ld_r2_0.5
# 12569673 0
#  397153 1
#       1 macromap_Ctrl_6_allchr_list_tier_2_ld_r2_0.5
# macromap_Ctrl_6_allchr_list_tier_2_ld_r2_0.75
# 12772614 0
#  194212 1
#       1 macromap_Ctrl_6_allchr_list_tier_2_ld_r2_0.75
# macromap_IFNB_24_allchr_list_tier_2_ld_r2_0.5
# 12559754 0
#  407072 1
#       1 macromap_IFNB_24_allchr_list_tier_2_ld_r2_0.5
# macromap_IFNB_24_allchr_list_tier_2_ld_r2_0.75
# 12759317 0
#  207509 1
#       1 macromap_IFNB_24_allchr_list_tier_2_ld_r2_0.75
# macromap_IFNB_6_allchr_list_tier_2_ld_r2_0.5
# 12560954 0
#  405872 1
#       1 macromap_IFNB_6_allchr_list_tier_2_ld_r2_0.5
# macromap_IFNB_6_allchr_list_tier_2_ld_r2_0.75
# 12759324 0
#  207502 1
#       1 macromap_IFNB_6_allchr_list_tier_2_ld_r2_0.75
# macromap_IFNG_24_allchr_list_tier_2_ld_r2_0.5
# 12560559 0
#  406267 1
#       1 macromap_IFNG_24_allchr_list_tier_2_ld_r2_0.5
# macromap_IFNG_24_allchr_list_tier_2_ld_r2_0.75
# 12758163 0
#  208663 1
#       1 macromap_IFNG_24_allchr_list_tier_2_ld_r2_0.75
# macromap_IFNG_6_allchr_list_tier_2_ld_r2_0.5
# 12564608 0
#  402218 1
#       1 macromap_IFNG_6_allchr_list_tier_2_ld_r2_0.5
# macromap_IFNG_6_allchr_list_tier_2_ld_r2_0.75
# 12768577 0
#  198249 1
#       1 macromap_IFNG_6_allchr_list_tier_2_ld_r2_0.75
# macromap_IL4_24_allchr_list_tier_2_ld_r2_0.5
# 12549697 0
#  417129 1
#       1 macromap_IL4_24_allchr_list_tier_2_ld_r2_0.5
# macromap_IL4_24_allchr_list_tier_2_ld_r2_0.75
# 12757862 0
#  208964 1
#       1 macromap_IL4_24_allchr_list_tier_2_ld_r2_0.75
# macromap_IL4_6_allchr_list_tier_2_ld_r2_0.5
# 12566567 0
#  400259 1
#       1 macromap_IL4_6_allchr_list_tier_2_ld_r2_0.5
# macromap_IL4_6_allchr_list_tier_2_ld_r2_0.75
# 12769819 0
#  197007 1
#       1 macromap_IL4_6_allchr_list_tier_2_ld_r2_0.75
# macromap_LIL10_24_allchr_list_tier_2_ld_r2_0.5
# 12573361 0
#  393465 1
#       1 macromap_LIL10_24_allchr_list_tier_2_ld_r2_0.5
# macromap_LIL10_24_allchr_list_tier_2_ld_r2_0.75
# 12769750 0
#  197076 1
#       1 macromap_LIL10_24_allchr_list_tier_2_ld_r2_0.75
# macromap_LIL10_6_allchr_list_tier_2_ld_r2_0.5
# 12567649 0
#  399177 1
#       1 macromap_LIL10_6_allchr_list_tier_2_ld_r2_0.5
# macromap_LIL10_6_allchr_list_tier_2_ld_r2_0.75
# 12766474 0
#  200352 1
#       1 macromap_LIL10_6_allchr_list_tier_2_ld_r2_0.75
# macromap_MBP_24_allchr_list_tier_2_ld_r2_0.5
# 12577430 0
#  389396 1
#       1 macromap_MBP_24_allchr_list_tier_2_ld_r2_0.5
# macromap_MBP_24_allchr_list_tier_2_ld_r2_0.75
# 12774128 0
#  192698 1
#       1 macromap_MBP_24_allchr_list_tier_2_ld_r2_0.75
# macromap_MBP_6_allchr_list_tier_2_ld_r2_0.5
# 12574816 0
#  392010 1
#       1 macromap_MBP_6_allchr_list_tier_2_ld_r2_0.5
# macromap_MBP_6_allchr_list_tier_2_ld_r2_0.75
# 12771968 0
#  194858 1
#       1 macromap_MBP_6_allchr_list_tier_2_ld_r2_0.75
# macromap_P3C_24_allchr_list_tier_2_ld_r2_0.5
# 12563145 0
#  403681 1
#       1 macromap_P3C_24_allchr_list_tier_2_ld_r2_0.5
# macromap_P3C_24_allchr_list_tier_2_ld_r2_0.75
# 12761947 0
#  204879 1
#       1 macromap_P3C_24_allchr_list_tier_2_ld_r2_0.75
# macromap_P3C_6_allchr_list_tier_2_ld_r2_0.5
# 12574325 0
#  392501 1
#       1 macromap_P3C_6_allchr_list_tier_2_ld_r2_0.5
# macromap_P3C_6_allchr_list_tier_2_ld_r2_0.75
# 12772584 0
#  194242 1
#       1 macromap_P3C_6_allchr_list_tier_2_ld_r2_0.75
# macromap_PIC_24_allchr_list_tier_2_ld_r2_0.5
# 12574431 0
#  392395 1
#       1 macromap_PIC_24_allchr_list_tier_2_ld_r2_0.5
# macromap_PIC_24_allchr_list_tier_2_ld_r2_0.75
# 12767909 0
#  198917 1
#       1 macromap_PIC_24_allchr_list_tier_2_ld_r2_0.75
# macromap_PIC_6_allchr_list_tier_2_ld_r2_0.5
# 12588488 0
#  378338 1
#       1 macromap_PIC_6_allchr_list_tier_2_ld_r2_0.5
# macromap_PIC_6_allchr_list_tier_2_ld_r2_0.75
# 12782878 0
#  183948 1
#       1 macromap_PIC_6_allchr_list_tier_2_ld_r2_0.75
# macromap_Prec_D0_allchr_list_tier_2_ld_r2_0.5
# 12598807 0
#  368019 1
#       1 macromap_Prec_D0_allchr_list_tier_2_ld_r2_0.5
# macromap_Prec_D0_allchr_list_tier_2_ld_r2_0.75
# 12784468 0
#  182358 1
#       1 macromap_Prec_D0_allchr_list_tier_2_ld_r2_0.75
# macromap_Prec_D2_allchr_list_tier_2_ld_r2_0.5
# 12557034 0
#  409792 1
#       1 macromap_Prec_D2_allchr_list_tier_2_ld_r2_0.5
# macromap_Prec_D2_allchr_list_tier_2_ld_r2_0.75
# 12761544 0
#  205282 1
#       1 macromap_Prec_D2_allchr_list_tier_2_ld_r2_0.75
# macromap_R848_24_allchr_list_tier_2_ld_r2_0.5
# 12559250 0
#  407576 1
#       1 macromap_R848_24_allchr_list_tier_2_ld_r2_0.5
# macromap_R848_24_allchr_list_tier_2_ld_r2_0.75
# 12761681 0
#  205145 1
#       1 macromap_R848_24_allchr_list_tier_2_ld_r2_0.75
# macromap_R848_6_allchr_list_tier_2_ld_r2_0.5
# 12566801 0
#  400025 1
#       1 macromap_R848_6_allchr_list_tier_2_ld_r2_0.5
# macromap_R848_6_allchr_list_tier_2_ld_r2_0.75
# 12767059 0
#  199767 1
#       1 macromap_R848_6_allchr_list_tier_2_ld_r2_0.75
# macromap_sLPS_24_allchr_list_tier_2_ld_r2_0.5
# 12562367 0
#  404459 1
#       1 macromap_sLPS_24_allchr_list_tier_2_ld_r2_0.5
# macromap_sLPS_24_allchr_list_tier_2_ld_r2_0.75
# 12768298 0
#  198528 1
#       1 macromap_sLPS_24_allchr_list_tier_2_ld_r2_0.75
# macromap_sLPS_6_allchr_list_tier_2_ld_r2_0.5
# 12583228 0
#  383598 1
#       1 macromap_sLPS_6_allchr_list_tier_2_ld_r2_0.5
# macromap_sLPS_6_allchr_list_tier_2_ld_r2_0.75
# 12773457 0
#  193369 1
#       1 macromap_sLPS_6_allchr_list_tier_2_ld_r2_0.75



# remove tmp and other files
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/*_tmp.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/*_sorted.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/*_sorted.bed.gz


# mv data to final folder
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/allchr_*_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/




