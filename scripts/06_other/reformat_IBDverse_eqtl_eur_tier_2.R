# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# ESTIMATE ENRICHMENT IN ANNOTATIONS:

# start getting ld tags from our in sample data - use r=0.5 as in Susie purity -> r2 0.25

# singularity exec iibdgc_postprocess_10_singularity.sif
path_gwas="/path/to/ibdgwas/IIBDGC/"

eqtl_study=IBDverse

mkdir ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/
mkdir -p ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/

path="/path/to/project"
conditions=($(ls /path/to/project

for condition in ${conditions[@]}
do
echo ${condition} 
done

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
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/${eqtl_study}_${condition}_allchr_list_tier_2_ld_r2_0.25_sorted.bed"
done


for condition in ${conditions[@]}
do 
echo ${condition} && tail -50 ${path_gwas}post_imputation/2022/log/allarray_${eqtl_study}_${condition}_allchr_bed3_ld_stdout | grep Successfully
done

# remove intermed files
rm ${path_gwas}post_imputation/2022/analysis/ld_data/eur_tier2/${eqtl_study}/${eqtl_study}_*_allchr_list_tier_2_ld_r2_0.25.bed

cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/
ls -la *_sorted.bed | wc -l
# 252

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


# IBDverse_dMean__B_0_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12054440 0
#  912386 1
#       1 IBDverse_dMean__B_0_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_0_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12008007 0
#  958819 1
#       1 IBDverse_dMean__B_0_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_0_r_all_allchr_list_tier_2_ld_r2_0.25
# 12225660 0
#  741163 1
#       3 3
#       1 IBDverse_dMean__B_0_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_0_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12122795 0
#  844031 1
#       1 IBDverse_dMean__B_0_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_10_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12236769 0
#  730057 1
#       1 IBDverse_dMean__B_10_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_10_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12273608 0
#  693218 1
#       1 IBDverse_dMean__B_10_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_10_r_all_allchr_list_tier_2_ld_r2_0.25
# 12403698 0
#  563128 1
#       1 IBDverse_dMean__B_10_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_10_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12375023 0
#  591803 1
#       1 IBDverse_dMean__B_10_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_12_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12357648 0
#  609178 1
#       1 IBDverse_dMean__B_12_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_12_r_all_allchr_list_tier_2_ld_r2_0.25
# 12464085 0
#  502741 1
#       1 IBDverse_dMean__B_12_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_12_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12376416 0
#  590410 1
#       1 IBDverse_dMean__B_12_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_13_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12424652 0
#  542174 1
#       1 IBDverse_dMean__B_13_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_13_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12419049 0
#  547777 1
#       1 IBDverse_dMean__B_13_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_1_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11987859 0
#  978967 1
#       1 IBDverse_dMean__B_1_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_1_r_all_allchr_list_tier_2_ld_r2_0.25
# 12144532 0
#  822294 1
#       1 IBDverse_dMean__B_1_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_1_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12039346 0
#  927480 1
#       1 IBDverse_dMean__B_1_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_2_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12305788 0
#  661038 1
#       1 IBDverse_dMean__B_2_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_2_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12021511 0
#  945315 1
#       1 IBDverse_dMean__B_2_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_2_r_all_allchr_list_tier_2_ld_r2_0.25
# 12222387 0
#  744439 1
#       1 IBDverse_dMean__B_2_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_2_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12063446 0
#  903380 1
#       1 IBDverse_dMean__B_2_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_3_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12254738 0
#  712088 1
#       1 IBDverse_dMean__B_3_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_3_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12167991 0
#  798835 1
#       1 IBDverse_dMean__B_3_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_3_r_all_allchr_list_tier_2_ld_r2_0.25
# 12296059 0
#  670767 1
#       1 IBDverse_dMean__B_3_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_3_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12268922 0
#  697904 1
#       1 IBDverse_dMean__B_3_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_4_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12304006 0
#  662820 1
#       1 IBDverse_dMean__B_4_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_4_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12159370 0
#  807456 1
#       1 IBDverse_dMean__B_4_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_4_r_all_allchr_list_tier_2_ld_r2_0.25
# 12352527 0
#  614299 1
#       1 IBDverse_dMean__B_4_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_4_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12188920 0
#  777906 1
#       1 IBDverse_dMean__B_4_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_5_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12281617 0
#  685209 1
#       1 IBDverse_dMean__B_5_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_5_r_all_allchr_list_tier_2_ld_r2_0.25
# 12450598 0
#  516228 1
#       1 IBDverse_dMean__B_5_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_5_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12284333 0
#  682493 1
#       1 IBDverse_dMean__B_5_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_6_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12161471 0
#  805355 1
#       1 IBDverse_dMean__B_6_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_6_r_all_allchr_list_tier_2_ld_r2_0.25
# 12274605 0
#  692221 1
#       1 IBDverse_dMean__B_6_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_6_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12187748 0
#  779078 1
#       1 IBDverse_dMean__B_6_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_8_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12122203 0
#  844623 1
#       1 IBDverse_dMean__B_8_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_8_r_all_allchr_list_tier_2_ld_r2_0.25
# 12242483 0
#  724343 1
#       1 IBDverse_dMean__B_8_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_8_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12153036 0
#  813790 1
#       1 IBDverse_dMean__B_8_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_blood_all_allchr_list_tier_2_ld_r2_0.25
# 11971047 0
#  995779 1
#       1 IBDverse_dMean__B_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11810407 0
# 1156419 1
#       1 IBDverse_dMean__B_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_r_all_allchr_list_tier_2_ld_r2_0.25
# 12026280 0
#  940546 1
#       1 IBDverse_dMean__B_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__B_ti_all_allchr_list_tier_2_ld_r2_0.25
# 11916555 0
# 1050271 1
#       1 IBDverse_dMean__B_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Colonocyte_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11659506 0
# 1307320 1
#       1 IBDverse_dMean__Colonocyte_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Colonocyte_r_all_allchr_list_tier_2_ld_r2_0.25
# 11636599 0
# 1330227 1
#       1 IBDverse_dMean__Colonocyte_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Colonocyte_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12151289 0
#  815534 1
#       3 3
#       1 IBDverse_dMean__Colonocyte_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Enterocyte_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11692049 0
# 1274777 1
#       1 IBDverse_dMean__Enterocyte_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Enterocyte_ti_all_allchr_list_tier_2_ld_r2_0.25
# 11695250 0
# 1271576 1
#       1 IBDverse_dMean__Enterocyte_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_0_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11746315 0
# 1220511 1
#       1 IBDverse_dMean__Epithelial_0_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_0_r_all_allchr_list_tier_2_ld_r2_0.25
# 11747635 0
# 1219191 1
#       1 IBDverse_dMean__Epithelial_0_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_10_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11876381 0
# 1090445 1
#       1 IBDverse_dMean__Epithelial_10_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_10_r_all_allchr_list_tier_2_ld_r2_0.25
# 11868669 0
# 1098157 1
#       1 IBDverse_dMean__Epithelial_10_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_12_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11943373 0
# 1023453 1
#       1 IBDverse_dMean__Epithelial_12_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_12_ti_all_allchr_list_tier_2_ld_r2_0.25
# 11954854 0
# 1011972 1
#       1 IBDverse_dMean__Epithelial_12_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_13_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11932714 0
# 1034112 1
#       1 IBDverse_dMean__Epithelial_13_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_13_r_all_allchr_list_tier_2_ld_r2_0.25
# 11938339 0
# 1028487 1
#       1 IBDverse_dMean__Epithelial_13_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_14_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11928693 0
# 1038133 1
#       1 IBDverse_dMean__Epithelial_14_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_14_ti_all_allchr_list_tier_2_ld_r2_0.25
# 11917924 0
# 1048902 1
#       1 IBDverse_dMean__Epithelial_14_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_15_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12073318 0
#  893508 1
#       1 IBDverse_dMean__Epithelial_15_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_15_r_all_allchr_list_tier_2_ld_r2_0.25
# 12151453 0
#  815370 1
#       3 3
#       1 IBDverse_dMean__Epithelial_15_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_15_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12156445 0
#  810381 1
#       1 IBDverse_dMean__Epithelial_15_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_16_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11932643 0
# 1034183 1
#       1 IBDverse_dMean__Epithelial_16_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_16_ti_all_allchr_list_tier_2_ld_r2_0.25
# 11938823 0
# 1028003 1
#       1 IBDverse_dMean__Epithelial_16_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_17_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11971860 0
#  994963 1
#       3 3
#       1 IBDverse_dMean__Epithelial_17_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_17_r_all_allchr_list_tier_2_ld_r2_0.25
# 11972065 0
#  994761 1
#       1 IBDverse_dMean__Epithelial_17_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_18_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12062465 0
#  904361 1
#       1 IBDverse_dMean__Epithelial_18_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_18_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12064947 0
#  901879 1
#       1 IBDverse_dMean__Epithelial_18_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_19_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12008117 0
#  958709 1
#       1 IBDverse_dMean__Epithelial_19_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_19_r_all_allchr_list_tier_2_ld_r2_0.25
# 12011760 0
#  955066 1
#       1 IBDverse_dMean__Epithelial_19_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_1_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11852102 0
# 1114721 1
#       3 3
#       1 IBDverse_dMean__Epithelial_1_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_1_ti_all_allchr_list_tier_2_ld_r2_0.25
# 11855591 0
# 1111232 1
#       3 3
#       1 IBDverse_dMean__Epithelial_1_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_20_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12161247 0
#  805579 1
#       1 IBDverse_dMean__Epithelial_20_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_20_r_all_allchr_list_tier_2_ld_r2_0.25
# 12159659 0
#  807167 1
#       1 IBDverse_dMean__Epithelial_20_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_22_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12120336 0
#  846490 1
#       1 IBDverse_dMean__Epithelial_22_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_22_r_all_allchr_list_tier_2_ld_r2_0.25
# 12180470 0
#  786356 1
#       1 IBDverse_dMean__Epithelial_22_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_22_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12116040 0
#  850786 1
#       1 IBDverse_dMean__Epithelial_22_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_23_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12033567 0
#  933259 1
#       1 IBDverse_dMean__Epithelial_23_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_23_r_all_allchr_list_tier_2_ld_r2_0.25
# 12095268 0
#  871558 1
#       1 IBDverse_dMean__Epithelial_23_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_23_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12147965 0
#  818861 1
#       1 IBDverse_dMean__Epithelial_23_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_24_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12135449 0
#  831377 1
#       1 IBDverse_dMean__Epithelial_24_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_24_r_all_allchr_list_tier_2_ld_r2_0.25
# 12145875 0
#  820951 1
#       1 IBDverse_dMean__Epithelial_24_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_25_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12125848 0
#  840978 1
#       1 IBDverse_dMean__Epithelial_25_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_25_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12134961 0
#  831865 1
#       1 IBDverse_dMean__Epithelial_25_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_26_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12149325 0
#  817501 1
#       1 IBDverse_dMean__Epithelial_26_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_26_r_all_allchr_list_tier_2_ld_r2_0.25
# 12171034 0
#  795792 1
#       1 IBDverse_dMean__Epithelial_26_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_26_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12283379 0
#  683447 1
#       1 IBDverse_dMean__Epithelial_26_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_27_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12168185 0
#  798641 1
#       1 IBDverse_dMean__Epithelial_27_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_27_r_all_allchr_list_tier_2_ld_r2_0.25
# 12213288 0
#  753538 1
#       1 IBDverse_dMean__Epithelial_27_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_27_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12279268 0
#  687558 1
#       1 IBDverse_dMean__Epithelial_27_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_28_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12160212 0
#  806614 1
#       1 IBDverse_dMean__Epithelial_28_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_28_r_all_allchr_list_tier_2_ld_r2_0.25
# 12286629 0
#  680197 1
#       1 IBDverse_dMean__Epithelial_28_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_28_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12205913 0
#  760913 1
#       1 IBDverse_dMean__Epithelial_28_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_29_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12195564 0
#  771262 1
#       1 IBDverse_dMean__Epithelial_29_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_29_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12194596 0
#  772230 1
#       1 IBDverse_dMean__Epithelial_29_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_2_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11802921 0
# 1163905 1
#       1 IBDverse_dMean__Epithelial_2_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_2_ti_all_allchr_list_tier_2_ld_r2_0.25
# 11802327 0
# 1164499 1
#       1 IBDverse_dMean__Epithelial_2_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_31_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12221158 0
#  745668 1
#       1 IBDverse_dMean__Epithelial_31_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_31_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12219709 0
#  747117 1
#       1 IBDverse_dMean__Epithelial_31_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_32_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12374479 0
#  592347 1
#       1 IBDverse_dMean__Epithelial_32_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_32_r_all_allchr_list_tier_2_ld_r2_0.25
# 12373464 0
#  593362 1
#       1 IBDverse_dMean__Epithelial_32_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_3_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11780637 0
# 1186189 1
#       1 IBDverse_dMean__Epithelial_3_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_3_r_all_allchr_list_tier_2_ld_r2_0.25
# 11767410 0
# 1199416 1
#       1 IBDverse_dMean__Epithelial_3_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_4_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11965664 0
# 1001162 1
#       1 IBDverse_dMean__Epithelial_4_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_4_ti_all_allchr_list_tier_2_ld_r2_0.25
# 11957194 0
# 1009632 1
#       1 IBDverse_dMean__Epithelial_4_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_5_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11839386 0
# 1127440 1
#       1 IBDverse_dMean__Epithelial_5_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_5_r_all_allchr_list_tier_2_ld_r2_0.25
# 11827892 0
# 1138934 1
#       1 IBDverse_dMean__Epithelial_5_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_6_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11995285 0
#  971541 1
#       1 IBDverse_dMean__Epithelial_6_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_6_r_all_allchr_list_tier_2_ld_r2_0.25
# 12200480 0
#  766346 1
#       1 IBDverse_dMean__Epithelial_6_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_6_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12023760 0
#  943066 1
#       1 IBDverse_dMean__Epithelial_6_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_7_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11860914 0
# 1105912 1
#       1 IBDverse_dMean__Epithelial_7_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_7_ti_all_allchr_list_tier_2_ld_r2_0.25
# 11859553 0
# 1107273 1
#       1 IBDverse_dMean__Epithelial_7_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_9_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11936016 0
# 1030810 1
#       1 IBDverse_dMean__Epithelial_9_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Epithelial_9_ti_all_allchr_list_tier_2_ld_r2_0.25
# 11930450 0
# 1036376 1
#       1 IBDverse_dMean__Epithelial_9_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mast_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12322042 0
#  644784 1
#       1 IBDverse_dMean__Mast_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mast_r_all_allchr_list_tier_2_ld_r2_0.25
# 12389566 0
#  577260 1
#       1 IBDverse_dMean__Mast_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mast_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12337087 0
#  629739 1
#       1 IBDverse_dMean__Mast_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_0_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12313679 0
#  653147 1
#       1 IBDverse_dMean__Mesenchymal_0_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_0_r_all_allchr_list_tier_2_ld_r2_0.25
# 12361354 0
#  605472 1
#       1 IBDverse_dMean__Mesenchymal_0_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_0_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12347943 0
#  618883 1
#       1 IBDverse_dMean__Mesenchymal_0_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_1_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12271376 0
#  695450 1
#       1 IBDverse_dMean__Mesenchymal_1_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_1_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12277286 0
#  689540 1
#       1 IBDverse_dMean__Mesenchymal_1_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_2_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12352465 0
#  614361 1
#       1 IBDverse_dMean__Mesenchymal_2_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_2_r_all_allchr_list_tier_2_ld_r2_0.25
# 12427694 0
#  539132 1
#       1 IBDverse_dMean__Mesenchymal_2_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_2_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12345440 0
#  621386 1
#       1 IBDverse_dMean__Mesenchymal_2_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_3_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12383755 0
#  583071 1
#       1 IBDverse_dMean__Mesenchymal_3_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_3_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12394156 0
#  572670 1
#       1 IBDverse_dMean__Mesenchymal_3_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_4_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12355625 0
#  611201 1
#       1 IBDverse_dMean__Mesenchymal_4_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_4_r_all_allchr_list_tier_2_ld_r2_0.25
# 12357416 0
#  609410 1
#       1 IBDverse_dMean__Mesenchymal_4_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_4_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12354445 0
#  612381 1
#       1 IBDverse_dMean__Mesenchymal_4_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_5_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12338497 0
#  628329 1
#       1 IBDverse_dMean__Mesenchymal_5_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_5_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12266025 0
#  700801 1
#       1 IBDverse_dMean__Mesenchymal_5_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_6_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12385882 0
#  580944 1
#       1 IBDverse_dMean__Mesenchymal_6_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_6_r_all_allchr_list_tier_2_ld_r2_0.25
# 12411600 0
#  555226 1
#       1 IBDverse_dMean__Mesenchymal_6_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_7_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12373909 0
#  592917 1
#       1 IBDverse_dMean__Mesenchymal_7_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_7_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12385139 0
#  581687 1
#       1 IBDverse_dMean__Mesenchymal_7_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12121885 0
#  844941 1
#       1 IBDverse_dMean__Mesenchymal_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_r_all_allchr_list_tier_2_ld_r2_0.25
# 12247969 0
#  718857 1
#       1 IBDverse_dMean__Mesenchymal_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Mesenchymal_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12161596 0
#  805230 1
#       1 IBDverse_dMean__Mesenchymal_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_0_blood_all_allchr_list_tier_2_ld_r2_0.25
# 11887382 0
# 1079444 1
#       1 IBDverse_dMean__Myeloid_0_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_0_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11908504 0
# 1058322 1
#       1 IBDverse_dMean__Myeloid_0_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_0_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12411797 0
#  555029 1
#       1 IBDverse_dMean__Myeloid_0_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_10_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12342822 0
#  624004 1
#       1 IBDverse_dMean__Myeloid_10_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_10_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12330070 0
#  636756 1
#       1 IBDverse_dMean__Myeloid_10_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_11_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12190890 0
#  775936 1
#       1 IBDverse_dMean__Myeloid_11_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_11_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12191009 0
#  775817 1
#       1 IBDverse_dMean__Myeloid_11_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_11_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12348774 0
#  618052 1
#       1 IBDverse_dMean__Myeloid_11_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_14_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12567929 0
#  398897 1
#       1 IBDverse_dMean__Myeloid_14_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_14_r_all_allchr_list_tier_2_ld_r2_0.25
# 12626227 0
#  340599 1
#       1 IBDverse_dMean__Myeloid_14_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_14_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12568512 0
#  398314 1
#       1 IBDverse_dMean__Myeloid_14_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_16_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12360877 0
#  605949 1
#       1 IBDverse_dMean__Myeloid_16_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_16_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12363852 0
#  602971 1
#       3 3
#       1 IBDverse_dMean__Myeloid_16_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_1_blood_all_allchr_list_tier_2_ld_r2_0.25
# 11905910 0
# 1060916 1
#       1 IBDverse_dMean__Myeloid_1_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_1_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11910970 0
# 1055856 1
#       1 IBDverse_dMean__Myeloid_1_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_2_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12380002 0
#  586824 1
#       1 IBDverse_dMean__Myeloid_2_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_2_r_all_allchr_list_tier_2_ld_r2_0.25
# 12577689 0
#  389137 1
#       1 IBDverse_dMean__Myeloid_2_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_2_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12387240 0
#  579586 1
#       1 IBDverse_dMean__Myeloid_2_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_3_blood_all_allchr_list_tier_2_ld_r2_0.25
# 11995366 0
#  971460 1
#       1 IBDverse_dMean__Myeloid_3_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_3_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11996498 0
#  970328 1
#       1 IBDverse_dMean__Myeloid_3_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_4_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12100376 0
#  866450 1
#       1 IBDverse_dMean__Myeloid_4_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_4_r_all_allchr_list_tier_2_ld_r2_0.25
# 12332011 0
#  634815 1
#       1 IBDverse_dMean__Myeloid_4_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_4_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12116516 0
#  850310 1
#       1 IBDverse_dMean__Myeloid_4_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_7_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12524386 0
#  442440 1
#       1 IBDverse_dMean__Myeloid_7_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_7_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12526651 0
#  440175 1
#       1 IBDverse_dMean__Myeloid_7_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_7_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12529492 0
#  437334 1
#       1 IBDverse_dMean__Myeloid_7_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_8_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12096391 0
#  870435 1
#       1 IBDverse_dMean__Myeloid_8_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_8_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12127166 0
#  839657 1
#       3 3
#       1 IBDverse_dMean__Myeloid_8_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_8_r_all_allchr_list_tier_2_ld_r2_0.25
# 12373196 0
#  593630 1
#       1 IBDverse_dMean__Myeloid_8_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_8_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12230026 0
#  736800 1
#       1 IBDverse_dMean__Myeloid_8_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_blood_all_allchr_list_tier_2_ld_r2_0.25
# 11756373 0
# 1210453 1
#       1 IBDverse_dMean__Myeloid_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11874660 0
# 1092166 1
#       1 IBDverse_dMean__Myeloid_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_r_all_allchr_list_tier_2_ld_r2_0.25
# 12248675 0
#  718151 1
#       1 IBDverse_dMean__Myeloid_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Myeloid_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12019684 0
#  947142 1
#       1 IBDverse_dMean__Myeloid_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Plasma_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12276203 0
#  690623 1
#       1 IBDverse_dMean__Plasma_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Plasma_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11942003 0
# 1024823 1
#       1 IBDverse_dMean__Plasma_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Plasma_r_all_allchr_list_tier_2_ld_r2_0.25
# 12137797 0
#  829029 1
#       1 IBDverse_dMean__Plasma_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Plasma_ti_all_allchr_list_tier_2_ld_r2_0.25
# 11964908 0
# 1001918 1
#       1 IBDverse_dMean__Plasma_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Platelet_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12346351 0
#  620475 1
#       1 IBDverse_dMean__Platelet_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Platelet_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12343851 0
#  622975 1
#       1 IBDverse_dMean__Platelet_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Secretory_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11642085 0
# 1324741 1
#       1 IBDverse_dMean__Secretory_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Secretory_r_all_allchr_list_tier_2_ld_r2_0.25
# 11728112 0
# 1238714 1
#       1 IBDverse_dMean__Secretory_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Secretory_ti_all_allchr_list_tier_2_ld_r2_0.25
# 11710161 0
# 1256665 1
#       1 IBDverse_dMean__Secretory_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Stem_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11723323 0
# 1243503 1
#       1 IBDverse_dMean__Stem_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Stem_r_all_allchr_list_tier_2_ld_r2_0.25
# 11940962 0
# 1025864 1
#       1 IBDverse_dMean__Stem_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__Stem_ti_all_allchr_list_tier_2_ld_r2_0.25
# 11759836 0
# 1206990 1
#       1 IBDverse_dMean__Stem_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_0_blood_all_allchr_list_tier_2_ld_r2_0.25
# 11902865 0
# 1063961 1
#       1 IBDverse_dMean__T_0_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_0_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12054528 0
#  912298 1
#       1 IBDverse_dMean__T_0_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_0_r_all_allchr_list_tier_2_ld_r2_0.25
# 12382404 0
#  584422 1
#       1 IBDverse_dMean__T_0_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_0_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12234920 0
#  731906 1
#       1 IBDverse_dMean__T_0_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_10_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12113128 0
#  853698 1
#       1 IBDverse_dMean__T_10_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_10_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12200175 0
#  766651 1
#       1 IBDverse_dMean__T_10_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_10_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12407892 0
#  558934 1
#       1 IBDverse_dMean__T_10_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_12_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12187964 0
#  778862 1
#       1 IBDverse_dMean__T_12_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_12_r_all_allchr_list_tier_2_ld_r2_0.25
# 12327164 0
#  639662 1
#       1 IBDverse_dMean__T_12_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_12_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12236862 0
#  729964 1
#       1 IBDverse_dMean__T_12_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_13_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12110879 0
#  855947 1
#       1 IBDverse_dMean__T_13_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_13_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12232759 0
#  734067 1
#       1 IBDverse_dMean__T_13_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_13_r_all_allchr_list_tier_2_ld_r2_0.25
# 12502487 0
#  464339 1
#       1 IBDverse_dMean__T_13_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_13_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12383664 0
#  583162 1
#       1 IBDverse_dMean__T_13_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_15_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12184228 0
#  782598 1
#       1 IBDverse_dMean__T_15_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_15_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12178763 0
#  788063 1
#       1 IBDverse_dMean__T_15_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_16_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12309887 0
#  656939 1
#       1 IBDverse_dMean__T_16_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_16_r_all_allchr_list_tier_2_ld_r2_0.25
# 12462077 0
#  504749 1
#       1 IBDverse_dMean__T_16_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_16_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12300357 0
#  666469 1
#       1 IBDverse_dMean__T_16_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_19_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12287360 0
#  679466 1
#       1 IBDverse_dMean__T_19_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_19_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12316848 0
#  649978 1
#       1 IBDverse_dMean__T_19_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_19_r_all_allchr_list_tier_2_ld_r2_0.25
# 12484694 0
#  482132 1
#       1 IBDverse_dMean__T_19_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_19_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12389022 0
#  577804 1
#       1 IBDverse_dMean__T_19_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_1_blood_all_allchr_list_tier_2_ld_r2_0.25
# 11854376 0
# 1112450 1
#       1 IBDverse_dMean__T_1_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_1_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12056711 0
#  910115 1
#       1 IBDverse_dMean__T_1_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_1_r_all_allchr_list_tier_2_ld_r2_0.25
# 12445052 0
#  521774 1
#       1 IBDverse_dMean__T_1_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_1_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12399204 0
#  567622 1
#       1 IBDverse_dMean__T_1_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_20_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12275571 0
#  691255 1
#       1 IBDverse_dMean__T_20_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_20_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12305124 0
#  661702 1
#       1 IBDverse_dMean__T_20_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_20_r_all_allchr_list_tier_2_ld_r2_0.25
# 12428178 0
#  538648 1
#       1 IBDverse_dMean__T_20_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_20_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12364033 0
#  602793 1
#       1 IBDverse_dMean__T_20_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_21_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12332072 0
#  634754 1
#       1 IBDverse_dMean__T_21_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_21_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12339600 0
#  627226 1
#       1 IBDverse_dMean__T_21_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_21_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12389927 0
#  576899 1
#       1 IBDverse_dMean__T_21_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_24_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12392999 0
#  573827 1
#       1 IBDverse_dMean__T_24_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_24_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12386553 0
#  580273 1
#       1 IBDverse_dMean__T_24_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_2_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12397043 0
#  569783 1
#       1 IBDverse_dMean__T_2_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_2_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12054317 0
#  912509 1
#       1 IBDverse_dMean__T_2_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_2_r_all_allchr_list_tier_2_ld_r2_0.25
# 12356469 0
#  610357 1
#       1 IBDverse_dMean__T_2_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_2_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12069626 0
#  897200 1
#       1 IBDverse_dMean__T_2_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_3_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12179382 0
#  787444 1
#       1 IBDverse_dMean__T_3_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_3_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12134803 0
#  832023 1
#       1 IBDverse_dMean__T_3_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_3_r_all_allchr_list_tier_2_ld_r2_0.25
# 12388850 0
#  577976 1
#       1 IBDverse_dMean__T_3_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_3_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12205631 0
#  761195 1
#       1 IBDverse_dMean__T_3_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_4_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12052705 0
#  914121 1
#       1 IBDverse_dMean__T_4_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_4_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12176753 0
#  790073 1
#       1 IBDverse_dMean__T_4_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_4_r_all_allchr_list_tier_2_ld_r2_0.25
# 12411326 0
#  555500 1
#       1 IBDverse_dMean__T_4_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_4_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12330299 0
#  636527 1
#       1 IBDverse_dMean__T_4_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_5_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12126234 0
#  840592 1
#       1 IBDverse_dMean__T_5_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_5_r_all_allchr_list_tier_2_ld_r2_0.25
# 12255105 0
#  711721 1
#       1 IBDverse_dMean__T_5_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_5_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12160035 0
#  806791 1
#       1 IBDverse_dMean__T_5_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_6_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12161861 0
#  804965 1
#       1 IBDverse_dMean__T_6_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_6_r_all_allchr_list_tier_2_ld_r2_0.25
# 12447599 0
#  519227 1
#       1 IBDverse_dMean__T_6_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_6_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12169324 0
#  797502 1
#       1 IBDverse_dMean__T_6_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_7_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12015258 0
#  951568 1
#       1 IBDverse_dMean__T_7_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_7_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12028397 0
#  938429 1
#       1 IBDverse_dMean__T_7_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_8_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12065104 0
#  901722 1
#       1 IBDverse_dMean__T_8_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_8_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12085612 0
#  881214 1
#       1 IBDverse_dMean__T_8_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_9_blood_all_allchr_list_tier_2_ld_r2_0.25
# 12034938 0
#  931888 1
#       1 IBDverse_dMean__T_9_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_9_ct_all_allchr_list_tier_2_ld_r2_0.25
# 12148713 0
#  818113 1
#       1 IBDverse_dMean__T_9_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_9_r_all_allchr_list_tier_2_ld_r2_0.25
# 12434393 0
#  532433 1
#       1 IBDverse_dMean__T_9_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_9_ti_all_allchr_list_tier_2_ld_r2_0.25
# 12352073 0
#  614753 1
#       1 IBDverse_dMean__T_9_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_blood_all_allchr_list_tier_2_ld_r2_0.25
# 11717297 0
# 1249529 1
#       1 IBDverse_dMean__T_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11731789 0
# 1235037 1
#       1 IBDverse_dMean__T_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_r_all_allchr_list_tier_2_ld_r2_0.25
# 12047490 0
#  919336 1
#       1 IBDverse_dMean__T_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_ti_all_allchr_list_tier_2_ld_r2_0.25
# 11835050 0
# 1131776 1
#       1 IBDverse_dMean__T_ti_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__unannotated_blood_all_allchr_list_tier_2_ld_r2_0.25
# 11653498 0
# 1313328 1
#       1 IBDverse_dMean__unannotated_blood_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__unannotated_ct_all_allchr_list_tier_2_ld_r2_0.25
# 11478505 0
# 1488321 1
#       1 IBDverse_dMean__unannotated_ct_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__unannotated_r_all_allchr_list_tier_2_ld_r2_0.25
# 11562973 0
# 1403853 1
#       1 IBDverse_dMean__unannotated_r_all_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__unannotated_ti_all_allchr_list_tier_2_ld_r2_0.25
# 11531968 0
# 1434858 1
#       1 IBDverse_dMean__unannotated_ti_all_allchr_list_tier_2_ld_r2_0.25

# remove tmp and other files
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/*_tmp.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/*_sorted.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/*_sorted.bed.gz


# mv data to final folder
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/${eqtl_study}/bed_b38/allchr_*_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/




