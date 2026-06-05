# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# singularity exec iibdgc_postprocess_10_singularity.sif
# singularity exec polyfun_2_singularity.sif

# path_gwas="/path/to/ibdgwas/IIBDGC/"
path_gwas="/path/to/ibdgwas/IIBDGC/"
pheno=(ibd cd uc)

# https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
# https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#conversion-to-liability-scale

# see also polyfun.R
# read https://www.med.unc.edu/pgc/wp-content/uploads/sites/959/2019/01/pgc_stat_bulik_2015.pdf

# For heritability and genetic covariance, it is customary to report heritability on the liability scale,
#  because liability scale heritability is comparable across
# studies with different prevalences. By default, the --h2 and --rg flags in ldsc output observed scale heritability.
# To convert to the liability scale, 
# we need to tell ldsc the sample and population prevalence for each trait using the --samp-prev and --pop-prev flags, 
# respectively. The population prevalence of scz and bip are both around 1%, 
# and the sample prevalence in each of these studies was about 50%, so
# 
# ldsc.py \
# --rg scz.sumstats.gz,bip.sumstats.gz \
# --ref-ld-chr eur_w_ld_chr/ \
# --w-ld-chr eur_w_ld_chr/ \
# --out scz_bip \
# --samp-prev 0.5,0.5 \
# --pop-prev 0.01,0.01
# 
# Note that the --samp-prev and --pop-prev flags also work with the --h2 flag. Conversion to liability scale affects only the SNP-heritability estimate; 
# it does not affect the LD Score regression intercept. The output is the same as before, except 'Observed' is replaced with 'Liability', 
# and the numbers are reported on the liability scale. 

# ldsc.py \
# --h2 scz.sumstats.gz \
# --ref-ld-chr eur_w_ld_chr/ \
# --w-ld-chr eur_w_ld_chr/ \
# --out scz_h2
# less scz_h2.log

# https://crohnsandcolitis.org.uk/media/4e5ccomz/epidemiology-summary-final.pdf

# From a 2020 baseline population of 16.1 million individuals, over 131,000 people were identified
# with IBD which represented a prevalence of 0.8%. This equates to one in every 123 people having a
# diagnosis of IBD. Crohn’s & Colitis UK have calculated that this equates to over 500,000 people in the
# UK currently living with this condition. Ulcerative colitis continues to have a higher prevalence (0.4%)
# compared with Crohn’s disease (0.3%) and unclassified IBD (0.07%) in the UK population.

# https://gut.bmj.com/content/68/11/1953
# Point prevalence estimates on 31 December 2018 for IBD overall, CD and UC were 725, 276 and 397 per 100000 people, respectively.

round(725/100000,3)
# [1] 0.007
round(276/100000,3)
# [1] 0.003
round(397/100000,3)
# [1] 0.004

# 1.- Estimate h2 on the liability scale:

MEM=10000

pheno=(ibd cd uc)
ppreval=(0.007 0.003 0.004)
spreval=(0.5 0.5 0.5)
# using now Neff, so the sample prevalence = 0.5

for i in {0..2}  
do echo ${pheno[${i}]} && echo ${ppreval[${i}]} && echo ${spreval[${i}]}
done

MEM=13000
for i in {0..2}  
do
bsub -J"h2_py2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/log/${pheno[${i}]}_eur_tier2_rate_0.5_info_0.9_allarrays_h2_ldsc_py2_stderr \
-o ${path_gwas}post_imputation/log/${pheno[${i}]}_eur_tier2_rate_0.5_info_0.9_allarrays_h2_ldsc_py2_stdout \
"ldsc.py \
--h2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9_nohla_sumstats_munged.parquet \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/eur_tier2_rate_0.5_info_0.9_allarrays/weights.april2025. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/eur_tier2_rate_0.5_info_0.9_allarrays/weights.april2025. \
--not-M-5-50 \
--samp-prev ${spreval[${i}]} \
--pop-prev ${ppreval[${i}]} \
--no-check-alleles \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/heritability/${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9_allarrays_heritability_estimation_liability_scale"
done




for i in {0..2}  
do echo ${pheno[${i}]} && cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/heritability/${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9_allarrays_heritability_estimation_liability_scale.log | \
grep -E "Total Liability scale h2|Lambda GC|Intercept|Ratio|SNPs remain"
done
# ibd
# After merging with reference panel LD, 12605850 SNPs remain.
# After merging with regression SNP LD, 12605850 SNPs remain.
# Total Liability scale h2: 0.2362 (0.015)
# Lambda GC: 1.3437
# Intercept: 1.0683 (0.004)
# Ratio: 0.0953 (0.0056)
# cd
# After merging with reference panel LD, 12768994 SNPs remain.
# After merging with regression SNP LD, 12768994 SNPs remain.
# Total Liability scale h2: 0.3054 (0.0242)
# Lambda GC: 1.2521
# Intercept: 1.0413 (0.0034)
# Ratio: 0.0743 (0.0061)
# uc
# After merging with reference panel LD, 12764456 SNPs remain.
# After merging with regression SNP LD, 12764456 SNPs remain.
# Total Liability scale h2: 0.2442 (0.0151)
# Lambda GC: 1.2627
# Intercept: 1.0459 (0.0034)
# Ratio: 0.0928 (0.0068)


#########################################################################################################################################

# rationale to keep lower freq variants, ie NealeLab using 1% in:
# https://nealelab.github.io/UKBB_ldsc/methods.html

# however to test this:

# keeping only variants with at different MAF thresholdss, still excluding HLA:


### CREATE MUNGE.PARQUET FILES

MEM=8000

for i in {0..2}  
do
echo ${pheno[${i}]} 
for maf in 0.005 0.01 0.05
do 
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${pheno[${i}]}_allchr_ldsc_sumstats_eur_tier2_rate_0.5_info_0.9_${maf}_nohla_mungen_format_stderr \
-o ${path_gwas}post_imputation/log/${pheno[${i}]}_allchr_ldsc_sumstats_eur_tier2_rate_0.5_info_0.9_${maf}_nohla_mungen_format_stdout \
"munge_polyfun_sumstats.py \
--sumstats ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz \
--min-info 0 \
--min-maf ${maf} \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9_maf_${maf}_nohla_sumstats_munged.parquet"
done
done


for i in {0..2}  
do
echo ${pheno[${i}]} 
for maf in 0.005 0.01 0.05
do echo ${maf} && tail -50 ${path_gwas}post_imputation/log/${pheno[${i}]}_allchr_ldsc_sumstats_eur_tier2_rate_0.5_info_0.9_${maf}_nohla_mungen_format_stdout | grep "Successfully"
done
done

for i in {0..2}  
do
echo ${pheno[${i}]} 
for maf in 0.005 0.01 0.05
do echo ${maf} && ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9_maf_${maf}_nohla_sumstats_munged.parquet
done
done

### ESTIMATE HERITABILITY

MEM=13000
for i in {0..2}  
do
echo ${pheno[${i}]} 
for maf in 0.05
do 
bsub -J"h2_py2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/log/${pheno[${i}]}_h2_ldsc_py2_common_0.005_${maf}_stderr \
-o ${path_gwas}post_imputation/log/${pheno[${i}]}_h2_ldsc_py2_common_0.005_${maf}_stdout \
"ldsc.py \
--h2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9_maf_${maf}_nohla_sumstats_munged.parquet \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/eur_tier2_rate_0.5_info_0.9_allarrays/weights.april2025. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/eur_tier2_rate_0.5_info_0.9_allarrays/weights.april2025. \
--not-M-5-50 \
--samp-prev ${spreval[${i}]} \
--pop-prev ${ppreval[${i}]} \
--no-check-alleles \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/heritability/${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9_allarrays_maf_${maf}_heritability_estimation_liability_scale"
done
done


for i in {0..2}  
do 
echo ${pheno[${i}]} &&
for maf in 0.005 0.01 0.05
do echo ${maf} && cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/heritability/${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9_allarrays_maf_${maf}_heritability_estimation_liability_scale.log | \
grep -E "Total Liability scale h2|Lambda GC|Intercept|Ratio|SNPs remain"
done
done





########################################################################################################################################
########################################################################################################################################
# follow the same approach as in 
# https://www.nature.com/articles/s41586-022-05275-y#Sec4

# To estimate the fraction of heritability that is explained by common variants within the 21% of the genome overlapping GWS loci, we calculated two genomic relationship matrices (GRMs)—one for SNPs within these loci and one for SNPs outside these loci—and then used both matrices to estimate a stratified SNP-based heritability (
# ) of height in eight independent samples of all five population groups represented in our METAFE (Fig. 3 and Methods). Altogether, our stratified estimation of SNP-based heritability shows that SNPs within these 7,209 GWS loci explain around 100% of 
#  in EUR and more than 90% of 
#  across all non-EUR groups, despite being drawn from less than 21% of the genome (Fig. 3). We also varied the window size used to define GWS loci and found that 35 kb was the smallest window size for which this level of saturation of SNP-based heritability could be achieved (Supplementary Fig. 18).


MEM=2000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(rcartocolor)
library(qvalue)

# path_gwas="/path/to/ibdgwas/IIBDGC/"
path_gwas="/path/to/ibdgwas/IIBDGC/"

all<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_conditional_on_known_index_gwas_fm_2023_with_iibdgc_tier1_tier2_eas_sas_summary_stats_plus_unsupervised_conditional.tsv"),head=T)

all<-all[,c("updated_region","chr","min","max")]
all<-all[!duplicated(all),]
head(all)

# all$start<-all$pos-35000
# all$end<-all$pos+35000
all$chr<-as.numeric(all$chr)

all<-all[order(all$chr,all$min,decreasing=F),]
all$chr<-paste("chr",all$chr,sep="")
# all<-all[,c("chr","start","end")]
all<-all[which(!is.na(all$min)),]

all$interval<-all$max-all$min
sum(all$interval)
# 582867541

# source genome size: https://www.ncbi.nlm.nih.gov/grc/human/data
 sum(all$interval)/(3099734149) 
# [1] 0.1880379

write.table(all[,2:4],,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/ibd_regions.bed"),
col.names=F,row.names=F,quote=F,sep="\t")


# create a similar one with just the old regions:

old<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_merged_regions_list_known_index_gwas_fm_2023.bed"),head=F)
old$start<-old$V2-500000
old$end<-old$V3+500000

old[which(old$start<0)]<-0

summary(old$start)
summary(old$end)

old$chr<-paste0("chr",old$V1)

old<-old[which(old$chr!="chr0"),]

write.table(old[,c("chr","start","end")],paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/ibd_old_regions.bed"),
col.names=F,row.names=F,quote=F,sep="\t")


##############################################################################
# # intersect each file with list of analysed variants, and add header


# 1.- define the list of variants, note that cd/uc/ibd analysis will include different cohorts, thus the N of variants is not the same as the 
# intersection of these

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed
# 12966826

# path_gwas="/path/to/ibdgwas/IIBDGC/"
path_gwas="/path/to/ibdgwas/IIBDGC/"


####################
# ALL REGIONS:

bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed \
-b ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/ibd_regions.bed -c \
| awk -F'\t' 'BEGIN {OFS=FS} {if ($5=="2") $5="1"; print $0}' | awk -F'\t' 'BEGIN {OFS=FS} {print $4,$5}' | awk -F'\t' -v OFS='\t' 'BEGIN {print "variant","ibd_regions"}1' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_ibd_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed

sed 's/\t[1-9]/\t1/gm' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_ibd_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_ibd_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset_tmp.bed

# ## reverse

sed 's/\t1/\ta/gm' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_ibd_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset_tmp.bed \
| sed 's/\t0/\tb/gm' | sed 's/\ta/\t0/gm' | sed 's/\tb/\t1/gm' | sed 's/\tibd_regions/\toutside_ibd_regions/gm' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_outside_ibd_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed

mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_ibd_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset_tmp.bed \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_ibd_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed

cut -f2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_ibd_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed | sort | uniq -c
# 10204376 0
# 2762450 1
#       1 ibd_regions

cut -f2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_outside_ibd_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed | sort | uniq -c
# 2762450 0
# 10204376 1
#       1 outside_ibd_regions

gzip ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_*_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed

# mv data to final folder
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_*_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/

####################
# ONLY OLD REGIONS:


bedtools intersect -a ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9_sorted.bed \
-b ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/ibd_old_regions.bed -c \
| awk -F'\t' 'BEGIN {OFS=FS} {if ($5=="2") $5="1"; print $0}' | awk -F'\t' 'BEGIN {OFS=FS} {print $4,$5}' | awk -F'\t' -v OFS='\t' 'BEGIN {print "variant","ibd_regions"}1' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_ibd_old_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed

sed 's/\t[1-9]/\t1/gm' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_ibd_old_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_ibd_old_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset_tmp.bed

# ## reverse

sed 's/\t1/\ta/gm' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_ibd_old_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset_tmp.bed \
| sed 's/\t0/\tb/gm' | sed 's/\ta/\t0/gm' | sed 's/\tb/\t1/gm' | sed 's/\tibd_regions/\toutside_ibd_regions/gm' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_outside_ibd_old_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed

mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_ibd_old_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset_tmp.bed \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_ibd_old_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed

cut -f2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_ibd_old_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed | sort | uniq -c
# 11517690 0
# 1449136 1
#       1 ibd_regions

cut -f2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_outside_ibd_old_regions_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed | sort | uniq -c
# 1449136 0
# 11517690 1
#       1 outside_ibd_regions

gzip ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_*_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed

# mv data to final folder
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ibd_regions/allchr_*_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/




##########################

cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/

files_ibd=($(ls | grep -E 'allchr_ibd|allchr_outside_ibd' | grep 'rate_0.5_info_0.9')) 
echo ${#files_ibd[@]}
# 4

files_ibd=($(ls | grep -E 'allchr_ibd|allchr_outside_ibd' | grep 'rate_0.5_info_0.9')) 
echo ${#files_ibd[@]}


files_ibd=($(echo ${files_ibd[*]} | sed 's/_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz//g' ))
files_ibd=($(echo ${files_ibd[*]} | sed 's/allchr_//g'))
echo ${#files_ibd[@]}

for i in ${files_ibd[*]}
do
echo ${i} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz \
| wc -l
done
# 12966827

# compare order of variants in file:
MEM=900

for i in ${files_ibd[*]}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/allchr_${i}_subset_check_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_subset_check_stdout \
"diff <(awk '{print \$1}' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/merged_files/template_variants_eur_tier2_rate_0.5_info_0.9) \
<(zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | awk '{print \$1}') \
>> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/var_check_${i}"
done

for i in ${files_ibd[*]}
do
echo ${i} && tail -30 ${path_gwas}post_imputation/log/allchr_${i}_subset_check_stdout | grep Exited
done

rm ~/tmp.txt 
for i in ${files_ibd[*]}
do
echo ${i} >> ~/tmp.txt && tail ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/var_check_${i} >> ~/tmp.txt 
done
# all OK
rm ~/tmp.txt 

for i in ${files_ibd[*]}
do
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/var_check_${i}
done



# paste together 
MEM=300
for i in ${files_ibd[*]}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q small \
-e ${path_gwas}post_imputation/log/allchr_${i}_paste_template_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_paste_template_stdout \
"paste ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/merged_files/template_variants_eur_tier2_rate_0.5_info_0.9 \
<(gzip -cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_rate_0.5_info_0.9_subset.bed.gz | cut -f2) \
| gzip > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}.april2025.annot.gz"
done

for i in ${files_ibd[*]}
do
echo ${i} && tail -30  ${path_gwas}post_imputation/log/allchr_${i}_paste_template_stdout | grep "Successfully"
done

for i in ${files_ibd[*]}
do
echo ${i} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}.april2025.annot.gz \
| wc -l
done
# 12966827


# make a folder for each files (to avoid >1000 files per directory):
for i in ${files_ibd[*]}
do
mkdir ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/
done

# subset by chr
MEM=200
for i in ${files_ibd[*]}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q small \
-e ${path_gwas}post_imputation/log/allchr_${i}_subset_per_chr_stderr \
-o ${path_gwas}post_imputation/log/allchr_${i}_subset_per_chr_stdout \
"zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}.april2025.annot.gz \
| awk '{print \$0 >> \"/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.april2025.\"\$2\"_tmp\"}'"
done

for i in ${files_ibd[*]}
do
echo ${i} && tail -30 ${path_gwas}post_imputation/log/allchr_${i}_subset_per_chr_stdout | grep "Successfully" 
done

for i in ${files_ibd[*]}
do
echo ${i} && ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.april2025.* | wc -l
done

for i in ${files_ibd[*]}
do
echo ${i} && wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.april2025.*
done

# # rm intermed file
for i in ${files_ibd[*]}
do
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_april2025/allchr_${i}_b38_list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_eur_tier2_subset.bed.gz
done


# The .annot file must contain the same SNPs in the same order as the .bim file; 
# sort out file using as template the bim file (some multiallelic varaints causing trouble)
# compare with bim file and add header
# create thin annotation files

MEM=1000
for i in ${files_ibd[*]}
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
for i in ${files_ibd[*]}
do
echo ${i} >> ~/tmp.txt && for chr in {1..22}
do 
echo ${chr} >> ~/tmp.txt && tail -30 ${path_gwas}post_imputation/log/allchr_${i}_${chr}_add_header_per_chr_stdout | grep -E "Successfully|Exit" >> ~/tmp.txt
done
done

cat ~/tmp.txt | grep Exited | wc -l
# 0
cat ~/tmp.txt | grep Successfully | wc -l
# 88

rm ~/tmp.txt

for i in ${files_ibd[*]}
do
ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.1.annot.gz
done

for i in ${files_ibd[*]}
do
echo ${i} && ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.*.annot.gz | wc -l
done


for chr in {1..22}
do
for i in ${files_ibd[*]}
do
echo ${i} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.${chr}.annot.gz | wc -l
done
done

for i in ${files_ibd[*]}
do
# MEM=60000 && for chr in 2
# MEM=55000 && for chr in 1 {3..6}
# MEM=45000 && for chr in {7..9}
# MEM=40000 && for chr in {10..12}
# MEM=30000 && for chr in {13..16}
# MEM=25000 && for chr in {17..19}
MEM=20000 && for chr in {20..22}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q week \
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

for i in ${files_ibd[*]}
do
for chr in {1..22}
do
echo ${chr} && tail -30 ${path_gwas}post_imputation/log/${chr}_${i}_ldsc_stdout | grep -E "Successfully|Exited|Max Memory"
done
done


######################################################################


cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/

# all but encode
files=($(ls | grep -E 'outside_ibd_regions|outside_ibd_old_regions' | grep -v jarvis | grep -v 35kb | grep -v allchr_ATAC_counts_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed.gz | grep -v allchr_ATAC_counts_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed.gz | grep -v .april2025.annot.gz | uniq | sed 's/outside_//')) 

echo ${#files[@]}
# 2

# submit the ldsc jobs

# issue with running this with only one annotation at a time stated here:
# https://github.com/bulik/ldsc/issues/138

# This is why both annotations (outside + included) were required

MEM=75000
pheno=(ibd cd uc)

for i in ${files[*]}
do
for ph in ${pheno[@]}
do
bsub -J"ldsc_pf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${ph}_${i}_ldsc_pf_no_baseline_stderr \
-o ${path_gwas}post_imputation/log/${ph}_${i}_ldsc_pf_no_baseline_stdout \
"ldsc.py \
--h2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/${ph}_eur_tier2_list_variants_rate_0.5_info_0.9_nohla_sumstats_munged.parquet \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/eur_tier2_rate_0.5_info_0.9_allarrays/weights.april2025. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/${i}/${i}.sorted.april2025.thin.,${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/no_baseline_files/outside_${i}/outside_${i}.sorted.april2025.thin. \
--overlap-annot \
--thin-annot \
--print-coefficients \
--print-delete-vals \
--not-M-5-50 \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/${ph}.no_baseline_both_annotations.${i}.files_version_april2025"
done
done

for i in ${files[*]}
do
for ph in ${pheno[@]}
do
echo ${i} && tail -50 ${path_gwas}post_imputation/log/${ph}_ibd_regions_ldsc_pf_no_baseline_stdout | grep -E "Successfully|Exit"
done
done

#### explore the results

MEM=2000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(rcartocolor)
library(qvalue)

# path_gwas="/path/to/ibdgwas/IIBDGC/"
path_gwas="/path/to/ibdgwas/IIBDGC/"


files<-list.files(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",sep=""))
files<-files[grep(".results",files)]

files<-files[grep("no_baseline_both_annotations",files)]
files<-files[!grepl("annotations.files",files)]

for (i in 1:length(files)) {

    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files[i],sep=""),head=T)
    tmp$pheno<-toupper(gsub(".no_baseline_both_annotations.*","",files[i]))
    tmp$class<-gsub(".*baseline_both_annotations.","",files[i])
    tmp$class<-gsub(".files_version_april2025.results","",tmp$class)
    tmp$class<-paste0(c("inside_","outside_"),tmp$class)
    # tmp<-tmp[which(tmp$Category %in% c("L2_1","L2_0")),]
    # tmp$Category<-gsub("L2_0","",tmp$Category)
    # tmp$Category<-gsub("L2_0","",tmp$Category)

    if(i==1) {
        enrich<-tmp
    } else {
        enrich<-rbind(enrich,tmp)
    }

}

enrich
#     Category Prop._SNPs Prop._h2 Prop._h2_std_error Enrichment
#       <char>      <num>    <num>              <num>      <num>
#  1:     L2_0    0.11176  0.58336           0.021811    5.21990
#  2:     L2_1    0.88811  0.41670           0.021814    0.46920
#  3:     L2_0    0.21304  0.76805           0.016558    3.60520
#  4:     L2_1    0.78683  0.23199           0.016561    0.29484
#  5:     L2_0    0.11176  0.60842           0.019017    5.44420
#  6:     L2_1    0.88811  0.39164           0.019020    0.44098
#  7:     L2_0    0.21304  0.79452           0.013502    3.72950
#  8:     L2_1    0.78683  0.20552           0.013504    0.26120
#  9:     L2_0    0.11176  0.56706           0.020203    5.07400
# 10:     L2_1    0.88811  0.43301           0.020206    0.48756
# 11:     L2_0    0.21304  0.77177           0.015188    3.62270
# 12:     L2_1    0.78683  0.22827           0.015191    0.29012
#     Enrichment_std_error Enrichment_p Coefficient Coefficient_std_error
#                    <num>        <num>       <num>                 <num>
#  1:             0.195160   1.8066e-29  2.6855e-07            1.8146e-08
#  2:             0.024562   1.8580e-29  2.4135e-08            1.2509e-09
#  3:             0.077725   5.3065e-39  1.8806e-07            1.0231e-08
#  4:             0.021048   5.4234e-39  1.5377e-08            9.6919e-10
#  5:             0.170170   1.0263e-32  1.9961e-07            1.2744e-08
#  6:             0.021417   1.0520e-32  1.6166e-08            8.0478e-10
#  7:             0.063377   2.1875e-43  1.3861e-07            7.1049e-09
#  8:             0.017162   2.2274e-43  9.7063e-09            5.6240e-10
#  9:             0.180780   4.0709e-30  2.1082e-07            1.4078e-08
# 10:             0.022752   4.1859e-30  2.0255e-08            1.0420e-09
# 11:             0.071293   1.0774e-41  1.5311e-07            7.9656e-09
# 12:             0.019306   1.1005e-41  1.2260e-08            7.3403e-10
#     Coefficient_z-score  pheno           class
#                   <num> <char>          <char>
#  1:              14.799     CD ibd_old_regions
#  2:              19.295     CD ibd_old_regions
#  3:              18.380     CD     ibd_regions
#  4:              15.866     CD     ibd_regions
#  5:              15.663    IBD ibd_old_regions
#  6:              20.087    IBD ibd_old_regions
#  7:              19.509    IBD     ibd_regions
#  8:              17.259    IBD     ibd_regions
#  9:              14.976     UC ibd_old_regions
# 10:              19.439     UC ibd_old_regions
# 11:              19.222     UC     ibd_regions
# 12:              16.702     UC     ibd_regions



# below old
enrich
#     Category Prop._SNPs Prop._h2 Prop._h2_std_error Enrichment
#       <char>      <num>    <num>              <num>      <num>
#  1:     L2_0   0.214390  0.75214           0.019864    3.50830
#  2:     L2_1   0.785080  0.24805           0.019877    0.31596
#  3:     L2_0   0.017359  0.50717           0.023504   29.21600
#  4:     L2_1   0.982110  0.49312           0.023516    0.50210
#  5:     L2_0   0.214390  0.75704           0.016965    3.53120
#  6:     L2_1   0.785080  0.24316           0.016975    0.30972
#  7:     L2_0   0.017359  0.53687           0.023345   30.92700
#  8:     L2_1   0.982110  0.46341           0.023357    0.47185
#  9:     L2_0   0.214390  0.72221           0.019359    3.36870
# 10:     L2_1   0.785080  0.27800           0.019371    0.35410
# 11:     L2_0   0.017359  0.48760           0.026355   28.08900
# 12:     L2_1   0.982110  0.51270           0.026368    0.52203
#     Enrichment_std_error Enrichment_p Coefficient Coefficient_std_error
#                    <num>        <num>       <num>                 <num>
#  1:             0.092656   1.3039e-33  1.8086e-07            1.0869e-08
#  2:             0.025318   1.4244e-33  1.6279e-08            1.2104e-09
#  3:             1.354000   2.6318e-33  1.5119e-06            1.0180e-07
#  4:             0.023944   2.9592e-33  2.5971e-08            1.6453e-09
#  5:             0.079130   4.3009e-38  1.4412e-07            7.9335e-09
#  6:             0.021622   4.7081e-38  1.2633e-08            7.8760e-10
#  7:             1.344800   3.6964e-35  1.2692e-06            8.2157e-08
#  8:             0.023782   4.1465e-35  1.9355e-08            1.1890e-09
#  9:             0.090300   1.0709e-36  1.4866e-07            8.2464e-09
# 10:             0.024674   1.1918e-36  1.5616e-08            1.1080e-09
# 11:             1.518200   3.1508e-31  1.2408e-06            8.7247e-08
# 12:             0.026848   3.6123e-31  2.3050e-08            1.5515e-09
#     Coefficient_z-score
#                   <num>
#  1:              16.640
#  2:              13.449
#  3:              14.851
#  4:              15.785
#  5:              18.166
#  6:              16.040
#  7:              15.449
#  8:              16.279
#  9:              18.027
# 10:              14.094
# 11:              14.222
# 12:              14.856
#                                                                                        pheno
#                                                                                       <char>
#  1:                           CD.NO_BASELINE_BOTH_ANNOTATIONS.FILES_VERSION_april2025.RESULTS
#  2:                           CD.NO_BASELINE_BOTH_ANNOTATIONS.FILES_VERSION_april2025.RESULTS
#  3:  CD.NO_BASELINE_BOTH_ANNOTATIONS.OUTSIDE_IBD_REGIONS_35KB.FILES_VERSION_april2025.RESULTS
#  4:  CD.NO_BASELINE_BOTH_ANNOTATIONS.OUTSIDE_IBD_REGIONS_35KB.FILES_VERSION_april2025.RESULTS
#  5:                          IBD.NO_BASELINE_BOTH_ANNOTATIONS.FILES_VERSION_april2025.RESULTS
#  6:                          IBD.NO_BASELINE_BOTH_ANNOTATIONS.FILES_VERSION_april2025.RESULTS
#  7: IBD.NO_BASELINE_BOTH_ANNOTATIONS.OUTSIDE_IBD_REGIONS_35KB.FILES_VERSION_april2025.RESULTS
#  8: IBD.NO_BASELINE_BOTH_ANNOTATIONS.OUTSIDE_IBD_REGIONS_35KB.FILES_VERSION_april2025.RESULTS
#  9:                           UC.NO_BASELINE_BOTH_ANNOTATIONS.FILES_VERSION_april2025.RESULTS
# 10:                           UC.NO_BASELINE_BOTH_ANNOTATIONS.FILES_VERSION_april2025.RESULTS
# 11:  UC.NO_BASELINE_BOTH_ANNOTATIONS.OUTSIDE_IBD_REGIONS_35KB.FILES_VERSION_april2025.RESULTS
# 12:  UC.NO_BASELINE_BOTH_ANNOTATIONS.OUTSIDE_IBD_REGIONS_35KB.FILES_VERSION_april2025.RESULTS





######################################################################################################################################################################
######################################################################################################################################################################

## plot those results, plus output from:
IIBDGC_GWAS/scripts/other/investigate_effect_of_using_different_ld_hm3_in_genetic_correlation.R
IIBDGC_GWAS/scripts/other/plot_heritability_liability_scale_different_thresholds.R


# COMPLETED JUNE 25