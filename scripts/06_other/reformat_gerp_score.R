# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# reformat_gerp

# data in b37 downloaded from:
https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=2157287288_bEA1MhRtv3WEQkDOtghJmNpbnUH8&clade=mammal&org=Human&db=hg19&hgta_group=compGeno&hgta_track=allHg19RS_BW&hgta_table=0&hgta_regionType=genome&position=chr7%3A155%2C592%2C223-155%2C605%2C565&hgta_outputType=wigData&hgta_outFileName=

path_gwas=/path/to/ibdgwas/IIBDGC/
cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/


# convert to bedgraph:
/path/to/software/username/./bigWigToBedGraph All_hg19_RS.bw All_hg19_RS.bg

MEM=8000
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/gerp_bedgraph_stderr \
-o ${path_gwas}post_imputation/log/gerp_bedgraph_stdout \
"/path/to/software/username/./bigWigToBedGraph \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/All_hg19_RS.bw \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/All_hg19_RS.bg"


# liftover:
${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/All_hg19_RS.bg \
${path_gwas}previous_qced_b38/liftover/hg19ToHg38.over.chain.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/All_hg19_RS_b38.bed \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/All_hg19_RS_nolifted_b38


# sort


# MEM=6500

# bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
# -e ${path_gwas}post_imputation/log/gerp_score_sort_stderr \
# -o ${path_gwas}post_imputation/log/gerp_score_sort_stdout \
# "cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/All_hg19_RS_b38.bed \
# | sort -u \
# | sort -k1,1V -k2,2n \
# > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/gerp_score.bed"



# # split file by chr
MEM=600
for chr in {1..21} X
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/gerp_score_chr${chr}_subset_sort_stderr \
-o ${path_gwas}post_imputation/log/gerp_score_chr${chr}_subset_sort_stdout \
"grep -w chr${chr} ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/All_hg19_RS_b38.bed \
| sort -k2,2n > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/chr${chr}_hg19_RS_b38_sorted.bed"
done

MEM=500
for chr in {1..22} X
do 
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/gerp_score_chr${chr}_merge_stderr \
-o ${path_gwas}post_imputation/log/gerp_score_chr${chr}_merge_stdout \
"/path/to/software/username/./bedtools merge -i ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/chr${chr}_hg19_RS_b38_sorted.bed -c 4 -d -1 \
-o mean -sorted > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/chr${chr}_hg19_RS_b38_sorted_nooverlaps.bed"
done

# intersect:

MEM=250
for chr in {1..22} X
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/gerp_score_chr${chr}_intersect_stderr \
-o ${path_gwas}post_imputation/log/gerp_score_chr${chr}_intersect_stdout \
"/path/to/software/username/./bedtools intersect -a <(grep -w chr${chr} ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed)  \
-b ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/chr${chr}_hg19_RS_b38_sorted_nooverlaps.bed \
-loj -sorted | awk -F'\t' 'BEGIN {OFS=FS} {print \$4,\$8}' | \
awk -F'\t' -v OFS='\t' 'BEGIN {print \"variant\",\"gerp_score\"}1' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/chr${chr}_gerp_score_list_union_variants_cd_uc_ibd_metaanalysis.bed"
done

rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/allchr_gerp_score_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed

# CONCATENATE
for chr in {1..22} X
do 
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/chr${chr}_gerp_score_list_union_variants_cd_uc_ibd_metaanalysis.bed | sed '1d' | uniq >> \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/allchr_gerp_score_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed
done

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/allchr_gerp_score_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed

cut -f2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/allchr_gerp_score_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed | sort -n | sed -n '1s/^/min=/p; $s/^/max=/p'
# min=-12.3
# max=9.95e-05
cut -f2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/allchr_gerp_score_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed | grep -w '\.' | wc -l
# 1512292

# distributed through the chrs by size (more likely in boundaries?)
for chr in {1..22} X 
do
echo ${chr} && cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/allchr_gerp_score_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed \
| grep -w '\.' | grep chr${chr}: | wc -l
done

# ADD HEADER AND Replace '.' with '0'
awk -F'\t' -v OFS='\t' 'BEGIN {print "variant","gerp_score"}1' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/allchr_gerp_score_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed \
| sed 's/\t\./\t0/gm' > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/allchr_gerp_score_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed

# double check min value is the 0:
cut -f2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/allchr_gerp_score_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed | sort -n | sed -n '1s/^/min=/p; $s/^/max=/p'
# min=-12.3
# max=9.95e-05

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/allchr_gerp_score_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed
49283611

# remove tmp and other files
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/All_hg19_RS.bg
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/*_sorted.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/All_hg19_RS_nolifted_b38
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/chr*_gerp_score_list_union_variants_cd_uc_ibd_metaanalysis.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/chr*_hg19_RS_b38_sorted_nooverlaps.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/*tmp.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/All_hg19_RS_b38.bed

gzip ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/allchr_gerp_score_list_union_variants_cd_uc_ibd_metaanalysis.bed

# mv data to final folder
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/gerp_score/allchr_gerp_score_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/


