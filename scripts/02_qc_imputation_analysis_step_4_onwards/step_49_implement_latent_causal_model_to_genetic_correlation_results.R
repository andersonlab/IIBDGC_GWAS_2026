# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# this stage follows the approached implemented by Luke Connor in:
# https://github.com/lukejoconnor/LCV/blob/master/R/RunLCV.R


# singularity exec iibdgc_postprocess_10_singularity.sif
# singularity exec iibdgc_postprocess_10_singularity.sif

path_gwas="/path/to/ibdgwas/IIBDGC/"

######################################################################################################
## take the list of significant rg results


gwas_id=($(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/genetic_correlation/list_traits_significantly_correlated_with_ibd_cd_uc)) 
echo ${#gwas_id[@]}

for i in {0..86}
do
echo ${gwas_id[i]}
done

# /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/genetic_correlation/list_traits_significantly_correlated_with_ibd_cd_uc


pheno=(ibd cd uc)

MEM=15000

for ph in ${pheno[@]}
do 
for i in {0..86}
do
bsub -J"lcv" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/lcv_${ph}_${gwas_id[i]}_stdout \
-e ${path_gwas}post_imputation/2022/log/lcv_${ph}_${gwas_id[i]}_stderr \
"Rscript /path/to/user/git/IIBDGC_GWAS/scripts/other/run_lcv.R ${ph} ${gwas_id[i]} > \
${path_gwas}post_imputation/2022/log/run_lcv_${ph}_${gwas_id[i]}.Rout"
done
done


for ph in ${pheno[@]}
do 
echo ${ph} && for i in {0..36}
do
echo ${gwas_id[i]} && tail -50  ${path_gwas}post_imputation/2022/log/lcv_${ph}_${gwas_id[i]}_stdout | grep -E "Successfully|TERM_MEMLIMIT"
done
done

# OUTPUT VARIABLES: 
#	lcv.output, a list with named entries:
#   "zscore", Z score for partial genetic causality. zscore>>0 implies gcp>0.
#	pval.gcpzero.2tailed, 2-tailed p-value for null that gcp=0.
#   "gcp.pm", posterior mean gcp (gcp=1: trait 1 -> trait 2; gcp=-1: trait 2-> trait 1); 
#   "gcp.pse", posterior standard error for gcp; 
#   "rho.est", estimated genetic correlation; 
#   "rho.err", standard error of rho estimate;
#   "pval.fullycausal" [2 entries], p-values for null that gcp=1 or that gcp=-1, respectively; 
#   "h2.zscore" [2 entries], z scores for trait 1 and trait 2 being heritable, respectively;
#   we recommend reporting results for h2.zscore > 7 (a very stringent threshold).

# plot output
~/git/IIBDGC_GWAS/scripts/other/plot_lcv.R
