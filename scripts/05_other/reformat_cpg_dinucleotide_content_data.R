# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# reformat CpG dinucleotide content ±50 kb for CpG content

cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/

# keep 0 any bp <0
awk -v OFS="\t" '{print $1, $2-50000, $3+50000, $4}' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed \
| awk -v OFS="\t" '$2<0{$2=0}1' > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_plus_minus_50kb.bed

awk -v OFS="\t" '{print $1, $2-100000, $3+100000, $4}' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed \
| awk -v OFS="\t" '$2<0{$2=0}1' > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_plus_minus_1Mb.bed

# ###################
# # CpG content

# # https://genome-euro.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=regulation&hgta_track=cpgIslandExt&hgta_table=cpgIslandExt&hgta_doSchema=describe+table+schema
# # cpgIslandExt

# # Database: hg38    Primary Table: cpgIslandExt    Row Count: 31,144   Data last updated: 2018-08-10
# # Format description: Describes the CpG Islands (includes observed/expected ratio)
# # field	example	SQL type	info	description
# # bin	585	smallint(6)	range	Indexing field to speed chromosome range queries.
# # chrom	chr1	varchar(255)	values	Reference sequence chromosome or scaffold
# # chromStart	28735	int(10) unsigned	range	Start position in chromosome
# # chromEnd	29737	int(10) unsigned	range	End position in chromosome
# # name	CpG: 111	varchar(255)	values	CpG Island
# # length	1002	int(10) unsigned	range	Island Length
# # cpgNum	111	int(10) unsigned	range	Number of CpGs in island      - USE THIS
# # gcNum	731	int(10) unsigned	range	Number of C and G in island
# # perCpg	22.2	float	range	Percentage of island that is CpG
# # perGc	73	float	range	Percentage of island that is C or G
# # obsExp	0.85	float	range	Ratio of observed(cpgNum) to expected(numC*numG/length) CpG in island

# # Obs/Exp CpG = Number of CpG * N / (Number of C * Number of G)

# wc -l cpgIslandExt
# # 31145 cpgIslandExt

# # split file by chr
for chr in X
do
grep -w chr${chr} ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_plus_minus_50kb.bed > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/cpgislands/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_plus_minus_50kb.bed
done

cut -f 2-11 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/cpgislands/cpgIslandExt \
| awk '{if (NR!=1) {print}}' > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/cpgislands/cpgIslandExt.bed


# # write the N overlap, to later estimate approx N CpG
for chr in {1..22}
do
bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/cpgislands/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_plus_minus_50kb.bed \
-b ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/cpgislands/cpgIslandExt.bed -wo \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/cpgislands/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_plus_minus_50kb_cpgIslandExt_merged.bed
done


# time consuming, restrict this step to list of variants in:

path_gwas=/path/to/ibdgwas/IIBDGC/
for chr in {1..21}
do 
bash ~/git/IIBDGC_GWAS/scripts/other/estimate_CpG.sh ${chr}
done

# Job <351927..351991> is submitted to queue <week>.

cat <(zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/allchr_plus_minus_50kb_withCpG_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed.gz) \
<(awk -F'\t' 'BEGIN {OFS=FS} {print $4,$5}' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/cpgislands/chrX_list_union_variants_cd_uc_ibd_metaanalysis_plus_minus_50kb_withCpG.bed | sed "/^V4/d") \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/cpgislands/allchr_plus_minus_50kb_withCpG_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed


# combine into one file, remove header, sort by position and chr, and keep relevant columns

path_gwas=/path/to/ibdgwas/IIBDGC/
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/cpgislands/chr*_list_union_variants_cd_uc_ibd_metaanalysis_plus_minus_50kb_withCpG.bed \
| sed "/^V1/d" | sort -k1,1V -k2,2n | awk -F'\t' 'BEGIN {OFS=FS} {print $4,$5}' > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/cpgislands/allchr_plus_minus_50kb_withCpG_list_union_variants_cd_uc_ibd_metaanalysis.bed

wc -l  ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/cpgislands/allchr_plus_minus_50kb_withCpG_list_union_variants_cd_uc_ibd_metaanalysis.bed
# 49246387 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed

# all OK

# rm intemediate files:
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/cpgislands/chr*
  
# move file to final folder
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/cpgislands/allchr_plus_minus_50kb_withCpG_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/

# add header
awk -F'\t' -v OFS='\t' 'BEGIN {print "variant","50kb_withCpG"}1' \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/allchr_plus_minus_50kb_withCpG_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/tmp && \
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/tmp ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/allchr_plus_minus_50kb_withCpG_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/allchr_plus_minus_50kb_withCpG_list_union_variants_cd_uc_ibd_metaanalysis.bed
# 49246388
head ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/allchr_plus_minus_50kb_withCpG_list_union_variants_cd_uc_ibd_metaanalysis.bed


# ###########################################################################################################################################################################
