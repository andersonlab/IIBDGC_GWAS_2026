# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# ESTIMATE ENRICHMENT IN ANNOTATIONS:

# start getting ld tags from our in sample data - use r=0.5 as in Susie purity -> r2 0.25, and create bed files with those annotations
# see /path/to/user/git/IIBDGC_GWAS/scripts/other/reformat_${eqtl_study}_eqtl_eur_tier_2.R


# singularity exec iibdgc_postprocess_10_singularity.sif
# singularity exec polyfun_2_singularity.sif


path_gwas="/path/to/ibdgwas/IIBDGC/"



##########################################
######## eQTL annotations: add each eqtl annotation, one at a time

cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/
# files=($(ls | grep rate_0.5_info_0.9 | grep -E 'macromap')) 
# echo ${#files[@]}
# # 24 - completed


files=($(ls | grep rate_0.5_info_0.9 | grep -E 'hu_2021')) 
echo ${#files[@]}
# 1

files=($(ls | grep rate_0.5_info_0.9 | grep -E 'edQTL_GTEX')) 
echo ${#files[@]}
# 10

files=($(ls | grep rate_0.5_info_0.9 | grep -E 'eQTL_catalogue')) 
echo ${#files[@]}
# 224

files=($(ls | grep rate_0.5_info_0.9 | grep -E 'hu_2021|edQTL_GTEX|eQTL_catalogue')) 
echo ${#files[@]}
# 235

files=($(ls | grep rate_0.5_info_0.9 | grep -E 'IBDverse')) 
echo ${#files[@]}
# 252

files=($(ls | grep rate_0.5_info_0.9 | grep -E 'sparc')) 
echo ${#files[@]}

files=($(ls | grep rate_0.5_info_0.9 | grep -E 'decode')) 
echo ${#files[@]}

files=($(echo ${files[*]} | sed 's/_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz//g' ))
files=($(echo ${files[*]} | sed 's/\ballchr_//g'))
echo ${#files[@]}
# 235
# 252
# 2


for i in ${files[*]}
do
echo ${i}
done

for i in ${files[*]}
do
echo ${i} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz \
| wc -l
done
# 12966827

# compare order of variants in file:
MEM=900

for i in ${files[*]}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/allchr_${i}_subset_check_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_subset_check_stdout \
"diff <(awk '{print \$1}' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/merged_files/template_variants_eur_tier2_rate_0.5_info_0.9) \
<(zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | awk '{print \$1}') \
>> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/var_check_${i}"
done

# all end as exited
for i in ${files[*]}
do
echo ${i} && tail -30 ${path_gwas}post_imputation/log/allchr_${i}_subset_check_stdout | grep Exited
done

rm ~/tmp.txt
for i in ${files[*]}
do
echo ${i} >> ~/tmp.txt && tail ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/var_check_${i} >> ~/tmp.txt 
done
# all OK
rm ~/tmp.txt 

for i in ${files[*]}
do
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/var_check_${i}
done



# paste together 
MEM=300
for i in ${files[*]}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q small \
-e ${path_gwas}post_imputation/log/allchr_${i}_paste_template_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_paste_template_stdout \
"paste ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/merged_files/template_variants_eur_tier2_rate_0.5_info_0.9 \
<(gzip -cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | cut -f2) \
| gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}.april2025.annot.gz"
done

for i in ${files[*]}
do
echo ${i} && tail -30  ${path_gwas}post_imputation/log/allchr_${i}_paste_template_stdout | grep "Successfully"
done

for i in ${files[*]}
do
echo ${i} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}.april2025.annot.gz \
| wc -l
done
# 12966827


# make a folder for each files (to avoid >1000 files per directory):
for i in ${files[*]}
do
mkdir ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/
done

# subset by chr
MEM=200
for i in ${files[*]}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q small \
-e ${path_gwas}post_imputation/log/allchr_${i}_subset_per_chr_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_subset_per_chr_stdout \
"zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}.april2025.annot.gz \
| awk '{print \$0 >> \"/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.april2025.\"\$2\"_tmp\"}'"
done

for i in ${files[*]}
do
echo ${i} && tail -30 ${path_gwas}post_imputation/log/allchr_${i}_subset_per_chr_stdout | grep "Successfully" 
done

for i in ${files[*]}
do
echo ${i} && ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.april2025.* | wc -l
done


# rm intermed file
# for i in ${files[*]}
# do
# rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz
# done


# The .annot file must contain the same SNPs in the same order as the .bim file; 
# sort out file using as template the bim file (some multiallelic varaints causing trouble)
# compare with bim file and add header
# create thin annotation files


MEM=1000
for i in ${files[*]}
do
for chr in {1..22}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q small \
-e ${path_gwas}post_imputation/log/allchr_${i}_${chr}_add_header_per_chr_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_${chr}_add_header_per_chr_stdout \
"awk -F'\t' '{OFS=FS} NR==FNR {h[\$1] = \$0; next} {print h[\$2]}' \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.april2025.chr${chr}_tmp \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/allarrays_chr${chr}_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_rate_0.5_info_0.9.bim \
| cat <(cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.april2025.CHR_tmp) - \
| cut -f 7 | gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.${chr}.annot.gz"
done
done

rm ~/tmp.txt
for i in ${files[*]}
do
echo ${i} >> ~/tmp.txt && for chr in {1..22}
do 
echo ${chr} >> ~/tmp.txt && tail -30 ${path_gwas}post_imputation/log/allchr_${i}_${chr}_add_header_per_chr_stdout | grep -E "Successfully|Exit" >> ~/tmp.txt
done
done

cat ~/tmp.txt | grep Exited | wc -l
# 0
cat ~/tmp.txt | grep Successfully | wc -l
# 44

for i in ${files[*]}
do
ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.1.annot.gz
done

for i in ${files[*]}
do
echo ${i} && ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.*.annot.gz | wc -l
done


for chr in {1..22}
do
for i in ${files[*]}
do
echo ${i} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.${chr}.annot.gz | wc -l
done
done

# # rm intermed file
# for i in ${files[*]}
# do
# rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.april2025.*_tmp
# done

rm ~/tmp.txt


cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/

# files=($(ls | grep -E 'hu_2021|edQTL_GTEX|eQTL_catalogue|IBDverse' | grep '.april2025.annot.gz')) 
files=($(ls | grep -E 'pQTL' | grep '.april2025.annot.gz')) 
files=($(echo ${files[*]} | sed 's/.april2025.annot.gz//g'))
echo ${#files[@]}
# 3

for i in ${files[*]}
do
for chr in {1..22}
do
if [[ "${chr}" == "2" ]] ; then
MEM=60000 && queue=basement
elif [[ "${chr}" == "1" ]] ; then
MEM=55000 && queue=basement
elif [[ "${chr}" -ge 3 && "${chr}" -le 6 ]] ; then
MEM=55000 && queue=basement
elif [[ "${chr}" -ge 7 && "${chr}" -le 8 ]] ; then
MEM=45000 && queue=basement
elif [[ "${chr}" -ge 9 && "${chr}" -le 12 ]] ; then
MEM=40000 && queue=week
elif [[ "${chr}" -ge 13 && "${chr}" -le 16 ]] ; then
MEM=30000 && queue=week
elif [[ "${chr}" -ge 17 && "${chr}" -le 19 ]] ; then
MEM=25000 && queue=week
elif [[ "${chr}" -ge 20 && "${chr}" -le 22 ]] ; then
MEM=20000 && queue=week
fi
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G humgen-priority -q ${queue} \
-e ${path_gwas}post_imputation/log/${chr}_${i}_ldsc_stderr \
-o ${path_gwas}post_imputation/log/${chr}_${i}_ldsc_stdout \
"ldsc.py \
--bfile ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/allarrays_chr${chr}_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_rate_0.5_info_0.9 \
--l2 \
--ld-wind-kb 1000 \
--annot ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.${chr}.annot.gz \
--thin-annot \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.${chr}"
done
done

# CONTINUE HERE


# resubmit jobs:
# eQTL_catalogue_QTS000001_QTD000001_allchr_list_tier_2_ld_r2_0.25
# eQTL_catalogue_QTS000024_QTD000409_allchr_list_tier_2_ld_r2_0.25
# IBDverse_dMean__T_1_r_all_allchr_list_tier_2_ld_r2_0.25



rm ~/tmp_list.R
for i in ${files[*]}
do
echo ${i} >> ~/tmp_list.R && for chr in {1..22}
do 
echo ${chr} >> ~/tmp_list.R && tail -30 ${path_gwas}post_imputation/log/${chr}_${i}_ldsc_stdout | grep -E "Exited|Successfully" >> ~/tmp_list.R
done
done

cat  ~/tmp_list.R | grep Exited | wc -l
cat  ~/tmp_list.R | grep Successfully | wc -l

rm ~/tmp

cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/

files=($(ls | grep -E 'hu_2021|edQTL_GTEX|eQTL_catalogue|IBDverse|macromap' | grep 'r2_0.25' | grep '.april2025.annot.gz')) 
files=($(ls | grep -E 'pQTL' | grep 'r2_0.25' | grep '.april2025.annot.gz')) 
files=($(echo ${files[*]} | sed 's/.april2025.annot.gz//g'))
echo ${#files[@]}
# 511
# 3


# # submit the ldsc jobs

MEM=58000

pheno=(ibd cd uc)

for i in ${files[*]}
do
for ph in ${pheno[@]}
do
bsub -J"ldsc_pf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G humgen-priority -q long \
-e ${path_gwas}post_imputation/log/${ph}_${i}_ldsc_pf_stderr \
-o ${path_gwas}post_imputation/log/${ph}_${i}_ldsc_pf_stdout \
"ldsc.py \
--h2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/${ph}_eur_tier2_list_variants_rate_0.5_info_0.9_nohla_sumstats_munged.parquet \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/weights.april2025. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.sorted.april2025.thin.,${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin. \
--overlap-annot \
--thin-annot \
--print-coefficients \
--print-delete-vals \
--not-M-5-50 \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/${ph}.baseline.${i}.files_version_april2025"
done
done


rm ~/tmp_list.txt
for ph in ${pheno[@]}
do
for i in ${files[*]}
do
echo ${i} >> ~/tmp_list.txt && tail -50 ${path_gwas}post_imputation/log/${ph}_${i}_ldsc_pf_stdout | grep -E "Successfully|Exit" >> ~/tmp_list.txt && \
echo ${ph} >> ~/tmp_list.txt
done
done

cat  ~/tmp_list.txt | grep Exit | wc -l
cat  ~/tmp_list.txt | grep Successfully | wc -l


# # all completed


##################################################################
# ### run the same analysis with a negative control - height
# # see reformat_other_gwas_summary_statistics_to_munge.R


cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/

files=($(ls | grep -E 'macromap' | grep -E 'r2_0.25|r2_0.75' | grep '.april2025.annot.gz')) 
files=($(echo ${files[*]} | sed 's/.april2025.annot.gz//g'))
echo ${#files[@]}
# 48


MEM=58000
for i in ${files[*]}
do
bsub -J"ldsc_pf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q long \
-e ${path_gwas}post_imputation/log/height_${i}_ldsc_pf_height_stderr \
-o ${path_gwas}post_imputation/log/height_${i}_ldsc_pf_height_stdout \
"ldsc.py \
--h2 ${path_gwas}resources/gwas_summary_statistics/sakaue_2021_height_list_variants_maf_0.001_munged.parquet \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/weights.april2025. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/baseline_files/baseline.sorted.april2025.thin.,${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin. \
--overlap-annot \
--thin-annot \
--print-coefficients \
--print-delete-vals \
--not-M-5-50 \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/other_traits/height.baseline.${i}.files_version_april2025"
done



rm ~/tmp_list.txt
for i in ${files[*]}
do
echo ${i} >> ~/tmp_list.txt && tail -50 ${path_gwas}post_imputation/log/height_${i}_ldsc_pf_height_stdout | grep -E "Successfully|Exit" >> ~/tmp_list.txt && \
echo ${ph} >> ~/tmp_list.txt
done



################################################

# plot enrichment analysis results, and create report files 

# # ~/git/IIBDGC_GWAS/scripts/other/plot_enrichment_analyses.R


# # # cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/
# # cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/

# # # all but encode
# # files=($(ls | grep -E 'sc_atac_perianal' | grep -v non_multiome  | grep -v jarvis | grep -v allchr_ATAC_counts_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis.bed.gz | grep -v allchr_ATAC_counts_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis.bed.gz | grep -v .april2025.annot.gz | uniq)) 

# # # submitted for Kyle to check out
# # for i in ${files[*]}
# # do
# # cp ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.chr*.gz ${path_gwas}to_upload_to_iibdgc_globus/
# # done

### COMPLETED


cd_pQTL_decode_final_olink_ukb_bi_allchr_list_tier_2_ld_r2_0.25_ldsc_pf_stderr