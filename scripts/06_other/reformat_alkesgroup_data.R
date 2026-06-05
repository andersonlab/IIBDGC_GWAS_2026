# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# reformat alkesgroup_data
# ########################################################################################################################################################
# ########################################################################################################################################################

# # 2.- BED FILES PROVIDED BY ALKES GROUP

path_gwas=/path/to/ibdgwas/IIBDGC/
cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/baselineLD_v2.2_bedfiles.tgz
tar -xvf baselineLD_v2.2_bedfiles.tgz

cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b37/
files=($(ls | sed 's/\.bed//g' ))

for i in ${files[@]}
do 

${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b37/${i}.bed \
${path_gwas}previous_qced_b38/liftover/hg19ToHg38.over.chain.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/${i}_b38.bed \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/nolifted_${i}_b38

done

for i in ${files[@]}
do 
echo ${i} && cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/nolifted_${i}_b38 | sed "/^#/d" | wc -l
done

for i in ${files[@]}
do 
wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/${i}_b38.bed
done
# 

# # 0 /path/to/ibdgwas/IIBDGC/resources/annotation/alkesgroup_baseline_bedfiles/bed_b38/ASMC_b38.bed

# add chr to first column:
# cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b37/
# head ASMC.bed
# # # # 1	752795	752796	1.2280100000e+04	0.0301205
# # # # 1	752874	752875	1.2302500000e+04	0.998996
# # # # 1	752893	752894	1.2312900000e+04	0.172691

# sed -e 's/^/chr/' ASMC.bed > ASMC_edited.bed

# ${path_gwas}previous_qced_b38/liftover/liftOver \
# ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b37/ASMC_edited.bed \
# ${path_gwas}previous_qced_b38/liftover/hg19ToHg38.over.chain.gz \
# ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/ASMC_b38.bed \
# ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/nolifted_ASMC_b38

# wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/ASMC_b38.bed
# # 19536981 /path/to/ibdgwas/IIBDGC/resources/annotation/alkesgroup_baseline_bedfiles/bed_b38/ASMC_b38.bed

# rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/ASMC_b38.bed
# rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/ASMC_edited_b38.bed

ls ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b37/*.bed | wc -l
# 83
ls ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/*_b38.bed | wc -l
# 81 OK 

head -2000 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/Coding_UCSC_b38.bed | uniq -d
# this introduces some duplicated rows, remove duplicates after sorting out the files:

awk '!seen[$0]++' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/Coding_UCSC_b38.bed \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/Coding_UCSC_b38_sorted_nodup.bed
# 195710 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/Coding_UCSC_b38_sorted.bed
# 195547 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/Coding_UCSC_b38_sorted_nodup.bed

sort -u ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/Coding_UCSC_b38.bed \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/Coding_UCSC_b38_sorted_nodup.bed
# 195710 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/Coding_UCSC_b38.bed
# 195547 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/Coding_UCSC_b38_sorted_nodup.bed


rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/nolifted_*
cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/
files=($(ls *_b38.bed | sed 's/\.bed//g' ))


# EXCLUDE DUP LINES, and SORT FILES:
for i in ${files[@]}
do
# KEEP ONLY CHR1 TO CHR22,  awk '$1 ~ /^chr(1?[0-9]|2[0-2])$/' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/${i}.bed \
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/${i}.bed \
| sort -u \
| sort -k1,1V -k2,2n \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/${i}_sorted.bed
done


files=($(ls *_b38_sorted.bed | sed 's/_sorted\.bed//g' ))

# # remove intermed files:
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/*b38.bed
  
ls ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/*_b38_sorted.bed | wc -l
# # 81 OK

# create intersection files with 
# ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis.bed

# # # 
# bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed \
# -b ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/Coding_UCSC_b38_sorted.bed -c \
# > ~/test.bed

# cat ~/test.bed | cut -f5 | sort | uniq -c
# # 10140616 0
# # 163273 1
# # 15 2

# awk '{if ($5 == 2) print $0;}' ~/test.bed
# # chr1	13391965	13391966	chr1:13391965:A:G	2
# # chr1	13392079	13392080	chr1:13392079:T:C	2
# # chr1	13392130	13392131	chr1:13392130:C:T	2
# # chr1	13392142	13392143	chr1:13392142:G:A	2
# # chr1	13392184	13392185	chr1:13392184:C:T	2
# # chr1	13392592	13392593	chr1:13392592:T:G	2
# # chr10	50008871	50008872	chr10:50008871:G:A	2
# # chr10	50008914	50008915	chr10:50008914:CAA:C	2
# # chr10	50009016	50009017	chr10:50009016:G:A	2
# # chr10	50009037	50009038	chr10:50009037:C:T	2
# # chr10	50009298	50009299	chr10:50009298:G:A	2
# # chr10	50009802	50009803	chr10:50009802:G:A	2
# # chr10	50009819	50009820	chr10:50009819:G:A	2
# # chr10	50009910	50009911	chr10:50009910:C:T	2

# # this is caused by overlapping regions (after liftover)
# # chr1	13391943	13392604
# # chr1	13391943	13392631

# # for binomial annotation, not an issue, edit 2 to 1
# # explore how this affects quantitative ones

# # # ~/test.bed
# # # 
# awk -F'\t' 'BEGIN {OFS=FS} {if ($5=="2") $5="1"; print $0}' ~/test.bed > ~/test2.bed 
# awk '{if ($5 == 2) print $0;}' ~/test2.bed
# # # 
# cat ~/test.bed | cut -f5 | sort | uniq -c
# # 10140616 0
# # 163273 1
# # 15 2
# # # 
# cat ~/test2.bed | cut -f5 | sort | uniq -c
# # 10140616 0
# # 163288 1 
# # IT WORKS WELL

# bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed \
# -b ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/Coding_UCSC_b38_sorted.bed -c \
# | awk -F'\t' 'BEGIN {OFS=FS} {if ($5=="2") $5="1"; print $0}' | awk -F'\t' 'BEGIN {OFS=FS} {print $4,$5}' > ~/test.bed


# intersect each file with list of analysed variants, and add header
for i in ${files[@]}
do
/path/to/software/username/./bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed \
-b ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/${i}_sorted.bed -c \
| awk -F'\t' 'BEGIN {OFS=FS} {if ($5=="2") $5="1"; print $0}' | awk -F'\t' 'BEGIN {OFS=FS} {print $4,$5}' | awk -F'\t' -v OFS='\t' 'BEGIN {print "variant","'"${i}"'"}1' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/allchr_${i}_list_union_variants_cd_uc_ibd_metaanalysis.bed
done


for i in ${files[@]}
do
wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/allchr_${i}_list_union_variants_cd_uc_ibd_metaanalysis.bed
done
# 49283611

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed
# # 49283610


# ## FOR THE QUANTITATIVE ONES, re-do keeping the cuantitative value:

# # alleleage_b38_sorted.bed
# # Backgrd_Selection_Stat_b38_sorted.bed
# # Human_Enhancer_Villar_Species_Enhancer_Count_b38_sorted.bed
# # BLUEPRINT_FE_META_TISSUE_DNAMETH_MaxCPP_b38_sorted.bed
# # BLUEPRINT_FE_META_TISSUE_H3K27ac_MaxCPP_b38_sorted.bed
# # BLUEPRINT_FE_META_TISSUE_H3K4me1_MaxCPP_b38_sorted.bed

files=(alleleage_b38 Backgrd_Selection_Stat_b38 Human_Enhancer_Villar_Species_Enhancer_Count_b38 
       BLUEPRINT_FE_META_TISSUE_DNAMETH_MaxCPP_b38 BLUEPRINT_FE_META_TISSUE_H3K27ac_MaxCPP_b38 
       BLUEPRINT_FE_META_TISSUE_H3K4me1_MaxCPP_b38 GTEx_FE_META_TISSUE_GE_MaxCPP_b38)

for i in ${files[@]}
do
/path/to/software/username/./bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed \
-b ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/${i}_sorted.bed \
-loj | awk -F'\t' 'BEGIN {OFS=FS} {print $4,$8}' | awk -F'\t' -v OFS='\t' 'BEGIN {print "variant","'"${i}"'"}1' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/allchr_${i}_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed
done

/path/to/software/username/./bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed \
-b ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/${i}_sorted.bed \
-loj | awk -F'\t' 'BEGIN {OFS=FS} {print $4,$8}' | awk -F'\t' -v OFS='\t' 'BEGIN {print "variant","'"${i}"'"}1' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/allchr_${i}_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed


# # some duplicated IDs in Backgrd_Selection_Stat_b38, merge regions that overlap, estimate mean value and run:
i=Backgrd_Selection_Stat_b38

/path/to/software/username/./bedtools merge -i ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/${i}_sorted.bed -c 4 \
-o mean > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/${i}_sorted_nooverlaps.bed

/path/to/software/username/./bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed \
-b ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/${i}_sorted_nooverlaps.bed \
-loj | awk -F'\t' 'BEGIN {OFS=FS} {print $4,$8}' | awk -F'\t' -v OFS='\t' 'BEGIN {print "variant","'"${i}"'"}1' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/allchr_${i}_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/allchr_${i}_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed
# 10303905 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/allchr_Backgrd_Selection_Stat_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed


# double check min value is always positive:
for i in ${files[@]}
do
echo ${i} && \
cut -f2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/allchr_${i}_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed | sort -n | sed -n '1s/^/min=/p; $s/^/max=/p'
done
# min=.
# max=497505
# Backgrd_Selection_Stat_b38
# min=.
# max=0.8812939347
# Human_Enhancer_Villar_Species_Enhancer_Count_b38
# min=.
# max=9
# BLUEPRINT_FE_META_TISSUE_DNAMETH_MaxCPP_b38
# min=.
# max=1.000000
# BLUEPRINT_FE_META_TISSUE_H3K27ac_MaxCPP_b38
# min=.
# max=1.000000
# BLUEPRINT_FE_META_TISSUE_H3K4me1_MaxCPP_b38
# min=.
# max=1.000000
# GTEx_FE_META_TISSUE_GE_MaxCPP_b38
# min=.
# max=1.000000

# for i in ${files[@]}
# do
# wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/allchr_${i}_list_union_variants_cd_uc_ibd_metaanalysis.bed
# done


# # Replace '.' with '0' - this is not accurate for allele age, do not use allele age in analysis

for i in ${files[@]}
do
sed 's/\t\./\t0/gm' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/allchr_${i}_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/allchr_${i}_list_union_variants_cd_uc_ibd_metaanalysis.bed
done

# double check min value is the 0:
for i in ${files[@]}
do
echo ${i} && \
cut -f2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/allchr_${i}_list_union_variants_cd_uc_ibd_metaanalysis.bed | sort -n | sed -n '1s/^/min=/p; $s/^/max=/p'
done
# alleleage_b38
# min=0
# max=497505
# Backgrd_Selection_Stat_b38
# min=0
# max=1
# Human_Enhancer_Villar_Species_Enhancer_Count_b38
# min=0
# max=9
# BLUEPRINT_FE_META_TISSUE_DNAMETH_MaxCPP_b38
# min=0
# max=1.000000
# BLUEPRINT_FE_META_TISSUE_H3K27ac_MaxCPP_b38
# min=0
# max=1.000000
# BLUEPRINT_FE_META_TISSUE_H3K4me1_MaxCPP_b38
# min=0
# max=1.000000
# GTEx_FE_META_TISSUE_GE_MaxCPP_b38
# min=0
# max=1.000000

# remove tmp and other files
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/*_tmp.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/*_sorted.bed

# # remove allele age, only a subset of individual position will have information about it:
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/allchr_alleleage_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed

# compress file:
gzip ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/*.bed

# mv data to final folder
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/alkesgroup_baseline_bedfiles/bed_b38/allchr_*_list_union_variants_cd_uc_ibd_metaanalysis.bed.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/


