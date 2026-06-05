# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#


# singularity exec iibdgc_postprocess_10_singularity.sif

path_gwas=/path/to/ibdgwas/IIBDGC/
cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/nucleotide_diversity/

# get CEU file from:
wget https://pophuman.uab.cat/files/wig/nuc_diversity_within_CEU_10kb.bw


# convert to bedgraph:
/path/to/software/username/./bigWigToBedGraph nuc_diversity_within_CEU_10kb.bw nuc_diversity_within_CEU_10kb.bg


# intersect:
/path/to/software/username/./bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed \
-b <(cut -f1-4 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/nucleotide_diversity/nuc_diversity_within_CEU_10kb.bg) \
-loj | awk -F'\t' 'BEGIN {OFS=FS} {print $4,$8}' | awk -F'\t' -v OFS='\t' 'BEGIN {print "variant","nuc_diversity_within_CEU_10kb"}1' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/nucleotide_diversity/allchr_nuc_diversity_within_CEU_10kb_b38_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed

head ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/nucleotide_diversity/allchr_nuc_diversity_within_CEU_10kb_b38_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed

# double check min value is always positive:
cut -f2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/nucleotide_diversity/allchr_nuc_diversity_within_CEU_10kb_b38_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed | sort -n | sed -n '1s/^/min=/p; $s/^/max=/p'

# # Replace '.' with '0' - this is not accurate for allele age, do not use allele age in analysis
sed 's/\t\./\t0/gm' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/nucleotide_diversity/allchr_nuc_diversity_within_CEU_10kb_b38_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/nucleotide_diversity/allchr_nuc_diversity_within_CEU_10kb_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed

cut -f2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/nucleotide_diversity/allchr_nuc_diversity_within_CEU_10kb_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed | sort -n | sed -n '1s/^/min=/p; $s/^/max=/p'

# mv data to final folder
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/nucleotide_diversity/allchr_nuc_diversity_within_CEU_10kb_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/

# remove intermed files
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/nucleotide_diversity/*tmp*