# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# ###########################################################################################################################################################################
# ######################################################################################################################################################################

# Viestra TF fingerprints:
# singularity exec iibdgc_postprocess_10_singularity.sif

path_gwas="/path/to/ibdgwas/IIBDGC/"

https://resources.altius.org/~jvierstra/projects/footprinting.2020/consensus.index/consensus_footprints_and_motifs_hg38.bed.gz

cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/viestra_tf/

wget https://resources.altius.org/~jvierstra/projects/footprinting.2020/consensus.index/consensus_footprints_and_collapsed_motifs_hg38.bed.gz
wget https://resources.altius.org/~jvierstra/projects/footprinting.2020/consensus.index/consensus_footprints_and_collapsed_motifs_hg38.bed.gz.tbi


# Column	Example	Description
# 1	contig	chr10	Chromosome
# 2	start	97320044	Start position (0-based)
# 3	stop	97320056	End position (start+1)
# 4	identifier	10.754379.4	Unique identifier (DHS_chom#.DHS_position%.fp_position%; DHS_chrom#.DHS_position% = DHS index identifier)
# 5	mean_signal	55.865317	Mean footprint -log(1-posterior) across biosamples (“confidence score”)
# 6	num_samples	9	Number of unique biosamples contributing to this index footprint (posterior >= 0.99)
# 7	num_fps	9	Number of unique footprints across all samples that contributed to this consensus footprint
# 8	width	12	Width of consensus footprint (column 3-column 2)
# 9	summit_pos	97320049	Estimated footprint summit position
# 10	core_start	97320042	Start position of core-region containing 95% of per-biosample summits
# 11	core_end	97320053	End position of core-region containing 95% of per-biosample summits
# 12	motif_clusters	CTCF;KLF/SP/2;ZNF563	Non-redundant motif archetype matches w/ 90% overlap, ; delimited


# # # some duplicated IDs, merge regions that overlap, keep max value and run:
MEM=500
bsub -J"ldsc_intersect" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/allchr_viestra_merge_eur_tier2_rate_0.5_stderr \
-o ${path_gwas}post_imputation/log/allchr_viestra_merge_eur_tier2_rate_0.5_stdout \
"bedtools merge -i <(zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/viestra_tf/consensus_footprints_and_collapsed_motifs_hg38.bed.gz | cut -f1-3,5) -c 4 \
-o max > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/viestra_tf/consensus_footprints_and_collapsed_motifs_hg38_nooverlaps.bed"


tail -30 ${path_gwas}post_imputation/log/allchr_viestra_merge_eur_tier2_rate_0.5_stdout | grep Successfully

# intersect, header and replace . by 0
MEM=4000
bsub -J"ldsc_intersect" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/allchr_viestra_intersect_eur_tier2_rate_0.5_stderr \
-o ${path_gwas}post_imputation/log/allchr_viestra_intersect_eur_tier2_rate_0.5_stdout \
"bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed \
-b ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/viestra_tf/consensus_footprints_and_collapsed_motifs_hg38_nooverlaps.bed \
-loj | cut -f 4,8 | awk -F'\\t' -v OFS='\\t' 'BEGIN {print \"variant\",\"viestra_tf\"}1' | sed 's/\\t\./\\t0/gm' | gzip \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/viestra_tf/allchr_viestra_tf_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz"

zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/viestra_tf/allchr_viestra_tf_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | wc -l
# 12966827

# mv data to final folder
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/viestra_tf/allchr_viestra_tf_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/


# remove intermediate files:
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/viestra_tf/*_tmp*.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/viestra_tf/*_nooverlaps*.bed

