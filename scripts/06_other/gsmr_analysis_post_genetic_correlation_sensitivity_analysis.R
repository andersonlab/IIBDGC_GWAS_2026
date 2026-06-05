# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# # singularity exec iibdgc_postprocess_10_singularity.sif
# path_gwas="/path/to/ibdgwas/IIBDGC/"

# to test
# MEM=10000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q normal R 

# gsmr requires the R-package(s)
# install.packages(c('survey'))
# install gsmr
# install.packages("https://yanglab.westlake.edu.cn/software/gsmr/static/gsmr_1.1.1.tar.gz",repos=NULL,type="source")

library(data.table)
library(R.utils)
library(dplyr)
library(ggplot2)
library("gsmr")
library(dplyr)

rm(list=ls())

# only multithreaded when writing
setDTthreads(2)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

# provide trait (in this case combination of database + trait) as input
args = commandArgs(trailingOnly=TRUE)

dataset<-args[2]
pheno<-args[1]

p_threshold<-5E-8

# # to test model
# pheno<-"ibd"
# pheno<-"uc"
# pheno<-"cd"
# dataset<-"GCST004988"
# dataset<-"GCST005531"
# dataset<-"GCST90029070"
# dataset<-"GCST005529"

print(dataset)
print(pheno)
library(stringr)

# load results and list of variants:
summary_data<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/plots/gsmr_eur_tier2_",dataset,"_",pheno,"_",p_threshold,"_data_source.tsv.gz",sep=""),head=T)

# load results with excluded (labelled as pleiotropic) variants:
data<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/gsmr_eur_tier2_",dataset,"_",pheno,"_",p_threshold,".tsv.gz",sep=""),head=T)
pleiotropic_snps<-unlist(strsplit(data$pleio_snps,split="\\|"))

# risk factor
gsmr_data<-fread(paste(path_gwas,"resources/gwas_summary_statistics/",dataset,"_edited_2.tsv.gz",sep=""),head=T)
gsmr_data<-as.data.frame(gsmr_data)

gsmr_data<-gsmr_data[which(gsmr_data$SNP %in% summary_data$SNP),]


# load the IBD|CD|UC summary stats - and retain
gsmr_data_ibd_traits<-fread(paste(path_gwas,"resources/gwas_summary_statistics/allchr_",pheno,"_eur_tier1_list_variants_only_SNPs_rate_0.5_with_pval.tsv.gz",sep=""),head=T)
gsmr_data_ibd_traits<-gsmr_data_ibd_traits[which(gsmr_data_ibd_traits$SNP %in% gsmr_data$SNP),]

gsmr_data<-gsmr_data[which(gsmr_data$SNP %in% gsmr_data_ibd_traits$SNP),]
gsmr_data<-gsmr_data[match(gsmr_data_ibd_traits$SNP,gsmr_data$SNP),]

table(gsmr_data$SNP==gsmr_data_ibd_traits$SNP)
table(gsmr_data$ALLELE1==gsmr_data_ibd_traits$ALLELE1)

# load LD matrix - created in
snp_coeff_id<-scan(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/ld_files/list_union_variants_gwas_significant_rg_significant_gsa_genotypes.xmat.gz",sep=""), what="", nlines=1)
snp_coeff_id<-snp_coeff_id[which(snp_coeff_id %in% gsmr_data$SNP)]

# we need >10 independent variants - >10 GWAS significant variants

if (length(snp_coeff_id)>10) {

    snp_coeff<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/ld_files/list_union_variants_gwas_significant_rg_significant_gsa_genotypes.xmat.gz",sep=""), header=T,select=snp_coeff_id)
    snp_coeff<-snp_coeff[-1,]

    snp_coeff_id<-colnames(snp_coeff)
    snp_coeff<-as.data.frame(snp_coeff)

    snp_coeff <- snp_coeff %>% mutate_if(is.character, as.numeric)

    # Match the SNP genotype data with the summary data
    snp_id<-Reduce(intersect, list(gsmr_data$SNP, snp_coeff_id))
    gsmr_data<-gsmr_data[match(snp_id, gsmr_data$SNP),]
    snp_order<-match(snp_id, snp_coeff_id)
    snp_coeff_id<-snp_coeff_id[snp_order]
    snp_coeff<-snp_coeff[, snp_coeff_id]

    # Calculate the LD correlation matrix
    ldrho<-cor(snp_coeff,use="pairwise.complete.obs")
    # the correlation or covariance between each pair of variables is computed using all complete pairs of observations on those variables

    # Check the size of the correlation matrix and double-check if the order of the SNPs in the LD correlation matrix is consistent with that in the GWAS summary data
    colnames(ldrho)<-snp_coeff_id
    rownames(ldrho)<-snp_coeff_id
    print(dim(ldrho))
    print(dim(gsmr_data))

    table(gsmr_data$SNP==gsmr_data_ibd_traits$SNP)
    table(gsmr_data$SNP==snp_coeff_id)


    # Define the paramaters to run GSMR:
    gwas_thresh = 5e-8    # GWAS threshold to select SNPs as the instruments for the GSMR analysis
    single_snp_heidi_thresh<-0.01    # p-value threshold for single-SNP-based HEIDI-outlier analysis
    multi_snps_heidi_thresh<-0.01    # p-value threshold for multi-SNP-based HEIDI-outlier analysis
    global_het_pval<-0.01  # p-value threshold when heidi is now on
    nsnps_thresh = 10   # the minimum number of instruments required for the GSMR analysis 
    heidi_outlier_flag = T    # flag for HEIDI-outlier analysis - keep same analysis without additioanl exclussions, just one out each time
    ld_r2_thresh = 0.05    # LD r2 threshold to remove SNPs in high LD
    ld_fdr_thresh = 0.05   # FDR threshold to remove the chance correlations between the SNP instruments
    gsmr2_beta = 1     # 0 - the original HEIDI-outlier method; 1 - the new global HEIDI-outlier method  

    # Sample size of the reference sample for LD
    if (pheno=="ibd") {
        n_ref<-57998
    } else if (pheno=="cd") {
        n_ref<-42599 
    } else if (pheno=="uc") {
        n_ref<-30577 
    }  

    # run the original analyses but without HEIDI (that is retaining the pleiotropic outliers)

    snpfreq<-gsmr_data$A1FREQ            # allele frequencies of the SNPs
    bzx<-gsmr_data$BETA     # effects of the instruments on risk factor
    bzx_se<-gsmr_data$SE       # standard errors of bzx
    bzx_n<-gsmr_data$N         # GWAS sample size for the risk factor
    bzx_pval<-gsmr_data$'P-value'  # p-values for bzx

    bzy<-gsmr_data_ibd_traits$BETA    # SNP effects on the disease
    bzy_se<-gsmr_data_ibd_traits$SE   # standard errors of bzy
    bzy_pval<-gsmr_data_ibd_traits$'P-value'    # p-values for bzy

    snp_coeff_id_tmp<-snp_coeff_id

    # extract results - with heidi:
    heidi_outlier_flag = F    # flag for HEIDI-outlier analysis - keep same analysis without additioanl exclussions, just one out each time

    gsmr_results_noheidi<-gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, single_snp_heidi_thresh, multi_snps_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)    # GSMR analysis 
    data$gsmr_no_heidi_bxy_pvalue<-gsmr_results_noheidi$bxy_pval

    # save main results:
    fwrite(data,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/gsmr_eur_tier2_",dataset,"_",pheno,"_",p_threshold,"_with_without_heidi.tsv.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")


    # leave one out approach - make some room for allowing this comparison
    nsnps_thresh = 10-3

    gsmr_data<-gsmr_data[which(!gsmr_data$SNP %in% pleiotropic_snps),]

    data<-as.data.frame(matrix(ncol=5,nrow=nrow(gsmr_data)))
    colnames(data)<-c("SNP_out","gsmr_with_heidi_bxy_pvalue","gsmr_with_heidi_N","gsmr_no_heidi_bxy_pvalue","gsmr_no_heidi_N")

    ldrho<-ldrho[gsmr_data$SNP,gsmr_data$SNP]
    snp_coeff_id<-snp_coeff_id[which(snp_coeff_id %in% gsmr_data$SNP)]

    gsmr_data_ibd_traits<-gsmr_data_ibd_traits[which(gsmr_data_ibd_traits$SNP %in% gsmr_data$SNP)]

    table(gsmr_data$SNP==gsmr_data_ibd_traits$SNP)
    table(gsmr_data$SNP==snp_coeff_id)

    for (n in 1:nrow(gsmr_data)) {

        data$SNP_out[n]<-gsmr_data$SNP[n]

        snpfreq<-gsmr_data$A1FREQ[-n]             # allele frequencies of the SNPs
        bzx<-gsmr_data$BETA[-n]     # effects of the instruments on risk factor
        bzx_se<-gsmr_data$SE[-n]       # standard errors of bzx
        bzx_n<-gsmr_data$N[-n]          # GWAS sample size for the risk factor
        bzx_pval<-gsmr_data$'P-value'[-n]   # p-values for bzx

        bzy<-gsmr_data_ibd_traits$BETA[-n]    # SNP effects on the disease
        bzy_se<-gsmr_data_ibd_traits$SE[-n]    # standard errors of bzy
        bzy_pval<-gsmr_data_ibd_traits$'P-value'[-n]    # p-values for bzy

        ldrho_tmp<-ldrho[-n,-n]
        snp_coeff_id_tmp<-snp_coeff_id[-n]

        # extract results - with heidi:
        heidi_outlier_flag = T    # flag for HEIDI-outlier analysis - keep same analysis without additioanl exclussions, just one out each time

        gsmr_results<-gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho_tmp, snp_coeff_id_tmp, n_ref, heidi_outlier_flag, gwas_thresh, single_snp_heidi_thresh, multi_snps_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)    # GSMR analysis 
        data$gsmr_with_heidi_bxy_pvalue[n]<-gsmr_results$bxy_pval
        data$gsmr_with_heidi_N[n]<-length(gsmr_results$used_index)

        # extract results - no heidi:
        heidi_outlier_flag = F    # flag for HEIDI-outlier analysis - keep same analysis without additioanl exclussions, just one out each time
        global_het_pval<-0.01
        gsmr_results<-gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho_tmp, snp_coeff_id_tmp, n_ref, heidi_outlier_flag, gwas_thresh, single_snp_heidi_thresh, multi_snps_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)    # GSMR analysis 
        data$gsmr_no_heidi_bxy_pvalue[n]<-gsmr_results$bxy_pval
        data$gsmr_no_heidi_N[n]<-length(gsmr_results$used_index)

    }  


    # save main results:
    fwrite(data,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/gsmr_eur_tier2_",dataset,"_",pheno,"_",p_threshold,"_sensitivity_test.tsv.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")


} 


q("no")


