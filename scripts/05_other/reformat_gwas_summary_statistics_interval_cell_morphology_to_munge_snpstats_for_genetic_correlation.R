# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
######################################################################################################################################################################
######################################################################################################################################################################
### MANUALLY ADD - Interval cell morphology GWAS - 40 traits

cd /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/

# If you use the data provided in this directory, please cite:
# Akbari, P., Vuckovic, D., Stefanucci, L. et al. A genome-wide association study of blood cell morphology identifies cellular proteins implicated in disease aetiology. Nat Commun 14, 5023 (2023). https://doi.org/10.1038/s41467-023-40679-y
# Published
# 18 August 2023
# DOI
# https://doi.org/10.1038/s41467-023-40679-y


# saved a map for this dataset based on supplementary data 2, and data shared by Klaudia 
# /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/Akbari_2023_supplementary_data_2.txt
# /path/to/project


# harmonise nomenclature and liftover to b38:
# singularity exec iibdgc_postprocess_10_singularity.sif
path_gwas="/path/to/ibdgwas/IIBDGC/"


gwas_id=($(ls /path/to/project | grep ".tsv.gz$" | \
sed 's/.tsv.gz//g'))
echo ${#gwas_id[@]}
# 63 - same as listed in summary file


MEM=28000

for i in {1..62}
do
bsub -J"subset_gwas" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/subset_gwas_variants_may2025_${gwas_id[i]}_stdout \
-e ${path_gwas}post_imputation/2022/log/subset_gwas_variants_may2025_${gwas_id[i]}_stderr \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/reformat_gwas_summary_statistics_interval_cell_morphology_gwas_summary_statistics_liftover_step1.R ${gwas_id[i]} > \
${path_gwas}post_imputation/2022/log/reformat_gwas_summary_statistics_interval_cell_morphology_gwas_summary_statistics_liftover_step1_${gwas_id[i]}.Rout"
done


for i in {1..62}
do
echo ${i} && tail -50 ${path_gwas}post_imputation/2022/log/subset_gwas_variants_may2025_${gwas_id[i]}_stdout  | grep "Successfully"
done

# (i=56 "RE-LYMP(L)%"" - manually ran)

# replace some special characters in the name:

gwas_id=($(ls /path/to/project | \
sed 's/-/_/g' | sed 's/(/_/g' | sed 's/)/_/g' | grep ".tsv.gz$" | sed 's/.tsv.gz//g'))
echo ${#gwas_id[@]}

for i in {0..62}; do echo ${i} && echo ${gwas_id[i]}; done

for i in {0..62}
do
ls -la ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_edited.tsv.gz
done

#########################################

# gunzip /path/to/ibdgwas/IIBDGC/resources/gnomad/gnomad_freq_edited.gz 

MEM=8200

for i in {0..62}
do
bsub -J"subset_gwas" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/subset_gwas_variants_2_${gwas_id[i]}_stdout \
-e ${path_gwas}post_imputation/2022/log/subset_gwas_variants_2_${gwas_id[i]}_stderr \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/subset_gwas_variants_for_genetic_correlation.R ${gwas_id[i]} > \
${path_gwas}post_imputation/2022/log/subset_gwas_variants_for_genetic_correlation_${gwas_id[i]}.Rout"
done

# continue here!

for i in  {0..62}
do
echo ${i} && tail -50 ${path_gwas}post_imputation/2022/log/subset_gwas_variants_2_${gwas_id[i]}_stdout | grep -E "Successfully|Exit"
done

for i in {0..62}
do
ls -la ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_edited_2.tsv.gz
done

for i in {0..62}
do
cat ${path_gwas}post_imputation/2022/log/subset_gwas_variants_for_genetic_correlation_${gwas_id[i]}.Rout
done

gzip /path/to/ibdgwas/IIBDGC/resources/gnomad/gnomad_freq_edited

#######################################################################################################
# reformat to munged file

path_gwas="/path/to/ibdgwas/IIBDGC/"
module unload HGI/softpack/groups/team152/iibdgc_postprocess/10

# version for polyfun/munge
# singularity exec polyfun_2_singularity.sif

MEM=7000

for i in {0..62}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${gwas_id[i]}_ldsc_sumstats_only_SNPs_polyfun_mungen_format_stderr \
-o ${path_gwas}post_imputation/log/${gwas_id[i]}_ldsc_sumstats_only_SNPs_polyfun_mungen_format_stdout \
"munge_polyfun_sumstats.py \
--sumstats ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_edited_2.tsv.gz \
--min-info 0 \
--min-maf 0 \
--remove-strand-ambig \
--out ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_sumstats_munged.parquet"
done

for i in {0..62}
do
echo ${i} && tail -50 ${path_gwas}post_imputation/log/${gwas_id[i]}_ldsc_sumstats_only_SNPs_polyfun_mungen_format_stdout  | grep -E "Successfully|Exit"
done

for i in {0..62}
do
echo ${i} && ls -la ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_sumstats_munged.parquet
done


# SAVE PARQUET FILE AS FLAT TEXT FILE AS INPUT FOR LDSC/LDSC.PY

module unload HGI/softpack/groups/team152/iibdgc_postprocess/10
module unload HGI/softpack/groups/team152/polyfun-2/1

MEM=3000
for i in {0..62}
do
bsub -J"subset_gwas" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/munge_parquet_to_sumstats_${gwas_id[i]}_stdout \
-e ${path_gwas}post_imputation/2022/log/munge_parquet_to_sumstats_${gwas_id[i]}_stderr \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/reformat_munge_parquet_to_sumtats_for_genetic_correlation.R ${gwas_id[i]} > \
${path_gwas}post_imputation/2022/log/reformat_munge_parquet_to_sumtats_for_genetic_correlation_${gwas_id[i]}.Rout"
done

for i in {0..62}
do
echo ${i} && tail -50 ${path_gwas}post_imputation/2022/log/munge_parquet_to_sumstats_${gwas_id[i]}_stdout | grep -E "Successfully|Exit"
done

for i in {0..62}
do
echo ${i} && ls -la ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_sumstats_munged.sumstats
done




# check out direction of effect sizes 

 zcat /path/to/project | head -1 && \
 zcat /path/to/project | \
 grep -E 'rs6993770|rs11774659'

# #VARIANT        ID_dbSNP        CHR     BP      GENPOS  REF     ALT     ALT_MINOR       DIRECTION       EFFECT  SE      P       MLOG10P GWSIG   ALT_FREQ        MA_FREQ R2    INFO
# 8:106581528_A_T rs6993770       8       106581528       1.26063 A       T       TRUE    +       8.44039e-02     7.93877e-03     2.11842e-26     2.56740e+01     TRUE    2.85702e-01    2.85702e-01     2.90769e-03     9.93255e-01
# 8:131190336_C_T rs11774659      8       131190336       1.53525 C       T       TRUE    -       -4.61969e-02    7.17871e-03     1.23258e-10     9.90918e+00     TRUE    4.91047e-01    4.91047e-01     1.06673e-03     9.90950e-01

# according to supplementary table 4 - effect size is from Allele0:
# rs6993770-T 0.08; AF T = 0.29
# rs11774659-T -0.42; AF T = 0.49

# first conversion 
path_gwas="/path/to/ibdgwas/IIBDGC/"

zcat ${path_gwas}resources/gwas_summary_statistics/H_IPF_edited.tsv.gz | head -1 && \
zcat ${path_gwas}resources/gwas_summary_statistics/H_IPF_edited.tsv.gz | grep -E '8:105569300|8:130178090'


# SNP     CHR     BP      A0      A1      BETA    SE      PVALUE  INFO    N       A1FREQ
# chr8:105569300:A:T      8       105569300       A       T       0.0844039       0.00793877      2.11842e-26     0.993255        41515   0.285702
# chr8:130178090:C:T      8       130178090       C       T       -0.0461969      0.00717871      1.23258e-10     0.99095 41515   0.491047

# COMPLETED JULY 25