# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# reformat_immune_cell_atlas
# ####################################
# ### IMMUNE CELL ATLAS

# # Immune cell atlas http://web.stanford.edu/group/pritchardlab/dataArchive/immune_atlas_web/index.html

# # atac-seq counts from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118189

# # Tn5 transposase adapters were trimmed with cutadapt version 1.13 with a minimum length of 20 and an overlap of 5 in paired-end mode.
# # Trimmed reads were aligned using bowtie2 version 2.2.9
# # Following removal of chrM reads we filtered with samtools version 1.4 and-f 2 and -F 1804 and MAPQ > 30, 
# # and removed duplicate reads with picard version 1.134
# # Used get_count function from nucleoATAC
# # Genome_build: hg19
# # Supplementary_files_format_and_content: tab delimited file including read count for each sample in each peak

# # The data used for downstream analysis thus typically consist of a count matrix, where each row corresponds to a genomic region (or “peak”)
# # and each column corresponds to a sample. The count in each cell of the matrix represents the number of read ends mapped to a particular
# # peak for a given sample, and is a proxy for the accessibility of the genomic region.

# wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118189/suppl/GSE118189%5FATAC%5Fcounts%2Etxt%2Egz

# singularity exec iibdgc_postprocess_10_singularity.sif
path_gwas=/path/to/ibdgwas/IIBDGC/

cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/
zcat GSE118189_ATAC_counts.txt.gz | wc -l
# 829943

# update b37 to b38:
zcat GSE118189_ATAC_counts.txt.gz | cut -f 1 | awk -F'_' '{print $1,$2,$3,$0}' | awk 'NR>1 {print}' > ATAC_counts_b37.bed
wc -l ATAC_counts_b37.bed
# 829942 ATAC_counts_b37.bed

${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b37.bed \
${path_gwas}previous_qced_b38/liftover/hg19ToHg38.over.chain.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b38.bed \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_nolifted

cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_nolifted | sed "/^#/d" | wc -l
# 180


# KEEP ONLY CHR1 TO CHR22, EXCLUDE DUP LINES, and SORT FILES:
awk '$1 ~ /^chr(1?[0-9]|2[0-2])$/' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b38.bed \
| sort -u \
| sort -k1,1V -k2,2n \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b38_sorted.bed


wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b37.bed
# 829942

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b38_sorted.bed
# 829577
# some regions lost after liftover

bedtools merge -i ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b38_sorted.bed | wc -l
# 829301
# some of the regions after liftover overlap


# create new file with chr start end B38
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/GSE118189_ATAC_counts.txt.gz | \
awk -F'\t' '{OFS=FS} NR==FNR {h[$1] = $0; next} {print h[$4],$0}' - ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b38_sorted.bed \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/tmp

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/tmp
# 829577 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/tmp

# get chr start end
cut -f 177-179 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/tmp > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/tmp1
cut -f 2-176 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/tmp > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/tmp2

paste ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/tmp1 \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/tmp2 \
| sort -k1,1V -k2,2n > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b38_edited.bed

head ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b38_edited.bed
tail ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b38_edited.bed
wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b38_edited.bed
# 829577

head -1 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b38_edited.bed | wc -w
# 178

seq -s, 178

# merge overlapping regions, and estimate mean counts for regions that are merged:
bedtools merge -i ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b38_edited.bed -c 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178 \
-o mean > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b38_edited_nooverlaps.bed

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b38_edited_nooverlaps.bed
# 829301 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b38_edited_nooverlaps.bed

zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/GSE118189_ATAC_counts.txt.gz | head -1 \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/header1


# # split by sample and cell type, combine by cell type

MEM=4000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

library(data.table)

head<-read.table("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/header1",head=T,check.names=F)
bed<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b38_edited_nooverlaps.bed",head=F)

head<-c("chr","start","end",as.character(colnames(head)))
colnames(bed)<-head

x<-colnames(bed)
x<-gsub("^[0-9]{4}-","",x)
x<-x[!duplicated(x)]
x<-x[!x %in% c("chr","start","end")]

bed<-as.data.frame(bed)
for (j in 1:length(x)) {
  tmp<-bed[,colnames(bed)[grep(x[j],colnames(bed))]]
  tmp[,x[j]]<-rowSums(tmp)
  tmp<-cbind(bed[,1:3],tmp[,ncol(tmp)])
  write.table(tmp,paste("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/split_files/ATAC_counts_",x[j],"_merged_samples_b38_edited_nooverlaps.bed",sep=""),
              col.names=F,row.names=F,quote=F,sep="\t")
}

q("no")

##############################

# rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/tmp*

# wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed
# # 25448333

# # interset original files:
# MEM=6000
# bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
# -e ${path_gwas}post_imputation/log/allchr_ATAC_subset_eur_tier2_rate_0.5_stderr \
# -o ${path_gwas}post_imputation/log/allchr_ATAC_subset_eur_tier2_rate_0.5_stdout \
# "bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed \
# -b ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/ATAC_counts_b38_edited_nooverlaps.bed \
# -loj | cut -f 4,8-182 > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/allchr_ATAC_counts_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_tmp.bed"


# wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/allchr_ATAC_counts_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_tmp.bed
# # 25448333 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/allchr_ATAC_counts_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_tmp.bed
# wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed
# # 25448333

# # add header:
# zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/GSE118189_ATAC_counts.txt.gz | head -1 \
# | awk -F'\t' -v OFS=FS '{print "variant",$0}' > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/header


# cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/header \
# ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/allchr_ATAC_counts_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_tmp.bed \
# > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/allchr_ATAC_counts_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_tmp1.bed

# wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/allchr_ATAC_counts_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_tmp1.bed
# # 25448334

# # Replace '.' with '0' - this is not accurate for allele age, do not use allele age in analysis
# sed 's/\t\./\t0/gm' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/allchr_ATAC_counts_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_tmp1.bed \
# > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/allchr_ATAC_counts_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5.bed


# # # mv data to final folder
# # mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/allchr_ATAC_counts_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5.bed \
# # /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/

# # wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_ATAC_counts_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5.bed
# # # 25448334


# ####################################################

# interset merged files:
cd /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/split_files/
files_atac=($(ls | grep ATAC_counts))

echo ${#files_atac[@]}
# 45
  
files_atac=($(echo ${files_atac[*]} | sed 's/\_b38\_edited\_nooverlaps\.bed//g' ))

for i in ${files_atac[@]}
do
echo ${i}
done

# intersect with list 25M variants, add header, Replace '.' with '0'
MEM=1000
for i in ${files_atac[@]}
do
bsub -J"ldsc_intersect" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/allchr_${i}_subset_eur_tier2_rate_0.5_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_subset_eur_tier2_rate_0.5_stdout \
"bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed \
-b ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/split_files/${i}_b38_edited_nooverlaps.bed \
-loj | cut -f 4,8 | awk -F'\t' -v OFS='\t' 'BEGIN {print \"variant\",\"'\"${i}\"'\"}1' \
| sed 's/\\t\./\\t0/gm' | gzip > \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/split_files/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz"
done


for i in ${files_atac[@]}
do
echo ${i} && tail -30 ${path_gwas}post_imputation/log/allchr_${i}_subset_eur_tier2_rate_0.5_stdout | grep "Successfully"
done

for i in ${files_atac[@]}
do
echo ${i} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/split_files/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | wc -l
done

# mv data to final folder
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/split_files/allchr_*_b38_list_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz \
/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/

# rm intermed files:
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/immune_cell_atlas/split_files/ATAC_counts_*


  
