# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# reformat_hickey_data

# download atac-seq peaks from:
# https://www.nature.com/articles/s41586-023-05915-x#Sec55
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.0zpc8672f



# ## extracted from README:

# # the ***_peak_matrix.mtx file contains the matrix data, 
# # the ***_peak_matrix_cells.tsv file contains the cell barcodes represented in the matrix, 
# # and the ***_peaks.bed contains the peaks represented in the matrix. 


# # https://kb.10xgenomics.com/hc/en-us/articles/360023561492-How-can-I-convert-the-peak-barcode-matrix-from-Cell-Ranger-ATAC-1-x-to-a-CSV-file

# # Print line number along with contents of barcodes.tsv and peaks.bed
# awk -F "\t" 'BEGIN { OFS = "," }; {print NR,$1}' barcodes.tsv | sort -t, -k 1b,1 > numbered_barcodes.csv
# awk -F "\t" 'BEGIN { OFS = "," }; {print NR,$1,$2,$3}' peaks.bed | sort -t, -k 1b,1 > numbered_peaks.csv

# # Skip the header lines and sort matrix.mtx
# tail -n +4 matrix.mtx | awk -F " " 'BEGIN { OFS = "," }; {print $1,$2,$3}' | sort -t, -k 1b,1 > peak_sorted_matrix.csv
# tail -n +4 matrix.mtx | awk -F " " 'BEGIN { OFS = "," }; {print $1,$2,$3}' | sort -t, -k 2b,2 > barcode_sorted_matrix.csv

# # Use join to replace line number with barcodes and peaks
# join -t, -1 1 -2 1 numbered_peaks.csv peak_sorted_matrix.csv | cut -d, -f 2,3,4,5,6 | sort -t, -k 4b,4 | join -t, -1 1 -2 4 numbered_barcodes.csv - | cut -d, -f 2,3,4,5,6 > final_matrix.csv

# # Remove temp files
# rm -f barcode_sorted_matrix.csv peak_sorted_matrix.csv numbered_barcodes.csv numbered_peaks.csv



path_gwas="/path/to/ibdgwas/IIBDGC/"
cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/hickey_2023/

files_hickey=($(ls | grep mtx | sed 's/_peak_matrix\.mtx//g' ))

ls -la *_matrix_cells.tsv | wc -l
ls -la *_peaks.bed | wc -l
ls -la *_matrix.mtx | wc -l
ls -la *_cell_types_*.tsv | wc -l
#12 OK

for i in ${files_hickey[@]}
do echo ${i}
done

# Print line number along with contents of barcodes.tsv and peaks.bed
MEM=250
for i in ${files_hickey[@]}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/allchr_${i}_numbered_peaks_eur_tier2_rate_0.5_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_numbered_peaks_eur_tier2_rate_0.5_stdout \
"awk -F \"\\t\" 'BEGIN { OFS = \",\" }; {print NR,\$1}' ${i}_peak_matrix_cells.tsv | sort -t, -k 1b,1 > ${i}_peak_matrix_cells_numbered_barcodes.csv && \
awk -F \" \" 'BEGIN { OFS = \",\" }; {print NR,\$1,\$2,\$3}' ${i}_peaks.bed | sort -t, -k 1b,1 > ${i}_peak_numbered_peaks.csv"
done


for i in ${files_hickey[@]}
do
echo ${i} && tail -50 ${path_gwas}post_imputation/log/allchr_${i}_numbered_peaks_eur_tier2_rate_0.5_stdout | grep "Successfully"
done

for i in ${files_hickey[@]}
do
ls -la ${i}_peak_numbered_peaks.csv
done


# Skip the header lines and sort matrix.mtx

MEM=400
for i in ${files_hickey[@]}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/allchr_${i}_sort_matrix_eur_tier2_rate_0.5_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_sort_matrix_eur_tier2_rate_0.5_stdout \
"tail -n +4 ${i}_peak_matrix.mtx | awk -F ' ' 'BEGIN { OFS = \",\" }; {print \$1,\$2,\$3}' | sort -t, -k 1b,1 > ${i}_peak_matrix_peak_sorted_matrix.csv && \
join -t, -1 1 -2 1 ${i}_peak_numbered_peaks.csv ${i}_peak_matrix_peak_sorted_matrix.csv | \
cut -d, -f 2,3,4,5,6 | sort -t, -k 4b,4 | join -t, -1 1 -2 4 ${i}_peak_matrix_cells_numbered_barcodes.csv - | \
cut -d, -f 2,3,4,5,6 > ${i}_hickney_peak_matrix_final_matrix.csv"
done


for i in ${files_hickey[@]}
do
echo ${i} && tail -50 ${path_gwas}post_imputation/log/allchr_${i}_sort_matrix_eur_tier2_rate_0.5_stderr | grep "Successfully"
done

for i in ${files_hickey[@]}
do
ls -la ${i}_hickney_peak_matrix_final_matrix.csv
done

# Use join to replace line number with barcodes and peaks

for i in ${files_hickey[@]}
do
ls -la ${i}_hickney_peak_matrix_final_matrix.csv
done

# Remove temp files
rm -f *_numbered_barcodes.csv *_numbered_peaks.csv *_sorted_matrix.csv

# create one file per cell type - recode cell type using 

MEM=40000

for i in ${files_hickey[@]}
do
bsub -J"ldsc_split" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/allchr_${i}_split_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_split_stdout \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/split_mtx_cell_ranger_per_cell_type.R ${i} > \
${path_gwas}post_imputation/log/split_mtx_cell_ranger_per_cell_type_${i}.Rout"
done

for i in ${files_hickey[@]}
do
echo ${i} && tail -50 ${path_gwas}post_imputation/log/allchr_${i}_split_stdout | grep -E 'Successfully|Exit'
done


cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/hickey_2023/
files_hickey=($(ls | grep -v peaks | grep bed | sed 's/\.bed//g' ))
echo ${#files_hickey[@]}
# files_hickey=($(echo ${files_hickey[*]} | sed 's/_hickey//g'))
# echo ${#files_hickey[@]}
# 180

for i in ${files_hickey[@]}
do echo ${i}
done

for i in ${files_hickey[@]}
do cut -f1 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/hickey_2023/${i}.bed | sort | uniq -c
done

for i in ${files_hickey[@]}
do ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/hickey_2023/${i}.bed
done

# intersect each file with list of analysed variants, and add header
MEM=1000
for i in ${files_hickey[@]}
do 
bsub -J"ldsc_intersect" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/allchr_${i}_subset_eur_tier2_rate_0.5_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_subset_eur_tier2_rate_0.5_stdout \
"bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed \
-b <(sed  '1d' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/hickey_2023/${i}.bed) -c \
| awk -F'\\t' 'BEGIN {OFS=FS} {if (\$5==\"2\") \$5=\"1\"; print \$0}' | awk -F'\\t' 'BEGIN {OFS=FS} {print \$4,\$5}' | awk -F'\\t' -v OFS='\\t' 'BEGIN {print \"variant\",\"'\"${i}\"'\"}1' \
| gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/hickey_2023/allchr_${i}_hickey_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz"
done

# CONTINUE HERE


for i in ${files_hickey[@]}
do
echo ${i} && tail -30 ${path_gwas}post_imputation/log/allchr_${i}_subset_eur_tier2_rate_0.5_stdout | grep "Successfully"
done

for i in ${files_hickey[@]}
do 
echo ${i} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/hickey_2023/allchr_${i}_hickey_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | wc -l
done


# mv data to final folder
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/hickey_2023/allchr_*_hickey_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/

