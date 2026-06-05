# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# reformat_zhang_data

# ### SC-ATAC SEQ FROM ZHANG ET AL

# # download scATAC-seq data from Zhang et al (https://www.biorxiv.org/content/biorxiv/early/2021/02/17/2021.02.17.431699.full.pdf; http://catlas.org/humantissue):
# # This is an updated link for the correct file: http://renlab.sdsc.edu/kai/data/other/cCRE_Accessibility.tsv.gz


# singularity exec iibdgc_postprocess_10_singularity.sif

path_gwas=/path/to/ibdgwas/IIBDGC/
cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/
# emacs README


# path_gwas=/path/to/ibdgwas/IIBDGC/
wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/cCRE_Accessibility.tsv 
# 756415 cCRE_Accessibility.tsv


# transfrom into bed files:
cut -f 2-55 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/cCRE_Accessibility.tsv \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/tmp

cut -f 1 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/cCRE_Accessibility.tsv | awk -F':' '{OFS="\t"; print $1,$2}' | \
awk -F'-' '{OFS="\t"; print $1,$2}' | paste - ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/tmp \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/cCRE_Accessibility_b38_tmp.bed


# KEEP ONLY CHR1 TO CHR22, EXCLUDE DUP LINES, and SORT FILES:
awk '$1 ~ /^chr(1?[0-9]|2[0-2])$/' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/cCRE_Accessibility_b38_tmp.bed \
| sort -u \
| sort -k1,1V -k2,2n \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/cCRE_Accessibility_b38_sorted.bed

head ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/cCRE_Accessibility_b38_sorted.bed | cut -f1-5
tail ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/cCRE_Accessibility_b38_sorted.bed | cut -f1-5
wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/cCRE_Accessibility_b38_sorted.bed
# 732908
head -1 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/cCRE_Accessibility_b38_sorted.bed | wc -w 
# 57


# intersect, header and replace . by 0
MEM=5000

bsub -J"ldsc_intersect" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/allchr_${i}_subset_eur_tier2_rate_0.5_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_subset_eur_tier2_rate_0.5_stdout \
"cat <(cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/cCRE_Accessibility.tsv | head -1 | cut -f2-55 \
| awk -F'\\t' '{OFS=FS} {print \"variant\",\$0}') <(bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed \
-b ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/cCRE_Accessibility_b38_sorted.bed \
-loj | cut -f 4,8-61) | sed 's/\\t\./\\t0/gm'| gzip > \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/allchr_cCRE_Accessibility_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz"

# CONTINUE HERE


less  ${path_gwas}post_imputation/log/allchr_${i}_subset_eur_tier2_rate_0.5_stdout

zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/allchr_cCRE_Accessibility_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | head -1 | wc -w
# 55
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/allchr_cCRE_Accessibility_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | tail -1 | wc -w
# 55

zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/allchr_cCRE_Accessibility_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | wc -l
# 12966827


#########################################
# split the file into individual files:


# split into different files and save in final folder:
for n in {2..55}
do
export TMP_TEST=$(zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/allchr_cCRE_Accessibility_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz \
                  | cut -f${n} | head -1) && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/allchr_cCRE_Accessibility_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz \
| cut -f1,${n} | gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/allchr_cCRE_Accessibility_${TMP_TEST}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz
done



cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/
files_ccre=($(ls | grep _b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | grep -v allchr_cCRE_Accessibility_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz))

for i in ${files_ccre[@]}
do
echo ${i} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/${i} | wc -l
done
# 12966827

echo ${#files_ccre[@]}
# 54

for i in ${files_ccre[@]}
do
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/zhang_sc_atacseq/${i} \
/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/
done

