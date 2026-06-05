# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# # singularity exec iibdgc_postprocess_10_singularity.sif
# path_gwas="/path/to/ibdgwas/IIBDGC/"

# to test
# MEM=25000
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
# dataset<-"GCST004988"
# dataset<-"GCST005531"
# dataset<-"GCST90029070"

print(dataset)
print(pheno)


gsmr_data<-fread(paste(path_gwas,"resources/gwas_summary_statistics/",dataset,"_edited_2.tsv.gz",sep=""),head=T)
gsmr_data<-as.data.frame(gsmr_data)

gsmr_data_to_retain<-fread(paste(path_gwas,"resources/gwas_summary_statistics/",dataset,"_only_SNPs_sumstats_munged.sumstats",sep=""),head=T)
gsmr_data_to_retain<-as.data.frame(gsmr_data_to_retain)

gsmr_data<-gsmr_data[which(gsmr_data$SNP %in% gsmr_data_to_retain$SNP),]

gsmr_data<-gsmr_data[match(gsmr_data_to_retain$SNP,gsmr_data$SNP),]
print(table(gsmr_data$ALLELE1==gsmr_data_to_retain$A1))
rm(gsmr_data_to_retain)


# keep only relatively freq samples
gsmr_data$A1FREQ<-as.numeric(gsmr_data$A1FREQ)
gsmr_data<-gsmr_data[which(gsmr_data$A1FREQ>=0.001 & gsmr_data$A1FREQ<=1-0.001),]
dim(gsmr_data)
#[1] 104809     11- GCST005531 - uc


# load the IBD|CD|UC summary stats - and retain
gsmr_data_ibd_traits<-fread(paste(path_gwas,"resources/gwas_summary_statistics/allchr_",pheno,"_eur_tier1_list_variants_only_SNPs_rate_0.5_with_pval.tsv.gz",sep=""),head=T)
gsmr_data_ibd_traits<-gsmr_data_ibd_traits[which(gsmr_data_ibd_traits$SNP %in% gsmr_data$SNP),]

# retain only strong instruments:
gsmr_data_ibd_traits$'P-value'<-as.numeric(gsmr_data_ibd_traits$'P-value')
gsmr_data_ibd_traits<-gsmr_data_ibd_traits[which(gsmr_data_ibd_traits$'P-value'<=p_threshold),]
dim(gsmr_data_ibd_traits)
# [1] 2831   11

gsmr_data<-gsmr_data[which(gsmr_data$SNP %in% gsmr_data_ibd_traits$SNP),]
gsmr_data_ibd_traits<-gsmr_data_ibd_traits[which(gsmr_data_ibd_traits$SNP %in% gsmr_data$SNP),]

gsmr_data<-gsmr_data[match(gsmr_data_ibd_traits$SNP,gsmr_data$SNP),]

table(gsmr_data$SNP==gsmr_data_ibd_traits$SNP)
table(gsmr_data$ALLELE1==gsmr_data_ibd_traits$ALLELE1)


# load LD matrix - created in
snp_coeff_id<-scan(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/ld_files/list_union_variants_ibd_cd_uc_gwas_significant_gsa_genotypes.xmat.gz",sep=""), what="", nlines=1)
snp_coeff_id<-snp_coeff_id[which(snp_coeff_id %in% gsmr_data$SNP)]

# we need >10 independent variants - >10 GWAS significant variants

if (length(snp_coeff_id)>10) {

    snp_coeff<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/ld_files/list_union_variants_ibd_cd_uc_gwas_significant_gsa_genotypes.xmat.gz",sep=""), header=T,select=snp_coeff_id)
    snp_coeff<-snp_coeff[-1,]

    snp_coeff_id<-colnames(snp_coeff)
    snp_coeff<-as.data.frame(snp_coeff)

    snp_coeff <- snp_coeff %>% mutate_if(is.character, as.numeric)

    # Match the SNP genotype data with the summary data
    snp_id<-Reduce(intersect, list(gsmr_data_ibd_traits$SNP, snp_coeff_id))
    gsmr_data_ibd_traits<-gsmr_data_ibd_traits[match(snp_id, gsmr_data_ibd_traits$SNP),]
    gsmr_data<-gsmr_data[match(snp_id, gsmr_data$SNP),]
    snp_order<-match(snp_id, snp_coeff_id)
    snp_coeff_id<-snp_coeff_id[snp_order]
    snp_coeff<-snp_coeff[, snp_order]

    # Calculate the LD correlation matrix
    ldrho<-cor(snp_coeff,use="pairwise.complete.obs")
    # the correlation or covariance between each pair of variables is computed using all complete pairs of observations on those variables

    # Check the size of the correlation matrix and double-check if the order of the SNPs in the LD correlation matrix is consistent with that in the GWAS summary data
    colnames(ldrho)<-snp_coeff_id
    rownames(ldrho)<-snp_coeff_id
    print(dim(ldrho))
    print(dim(gsmr_data))
    print(dim(gsmr_data_ibd_traits))


    # Standardise the effect size
    # This is an optional process. If the risk factor was not standardised in GWAS, the effect sizes can be scaled using the method below. Note that this process requires allele 
    # frequencies, z-statistics, and sample size. After scaling, bzx is interpreted as the per-allele effect of a SNP on the exposure in standard deviation units.

    snpfreq<-gsmr_data_ibd_traits$A1FREQ             # allele frequencies of the SNPs
    bzx<-gsmr_data_ibd_traits$BETA     # effects of the instruments on risk factor
    bzx_se<-gsmr_data_ibd_traits$SE       # standard errors of bzx
    bzx_n<-gsmr_data_ibd_traits$N          # GWAS sample size for the risk factor
    # std_zx<-std_effect(snpfreq, bzx, bzx_se, bzx_n)    # perform standardisation
    # gsmr_data_ibd_traits$std_bzx<-std_zx$b    # standardized bzx
    # gsmr_data_ibd_traits$std_bzx_se<-std_zx$se    # standardized bzx_se
    # head(gsmr_data)
    # bzx<-gsmr_data_ibd_traits$std_bzx    # SNP effects on the risk factor
    # bzx_se<-gsmr_data_ibd_traits$std_bzx_se    # standard errors of bzx
    bzx_pval<-gsmr_data_ibd_traits$'P-value'   # p-values for bzx

    bzy<-gsmr_data$BETA    # SNP effects on the disease
    bzy_se<-gsmr_data$SE    # standard errors of bzy
    bzy_pval<-gsmr_data$'P-value'    # p-values for bzy

    gsmr_data_ibd_traits$bzx<-bzx
    gsmr_data_ibd_traits$bzx_se<-bzx_se
    gsmr_data_ibd_traits$bzx_pval<-bzx_pval
    gsmr_data_ibd_traits$bzy_pval<-bzy_pval

    # Sample size of the reference sample for LD
    if (pheno=="ibd") {
        n_ref<-57998
    } else if (pheno=="cd") {
        n_ref<-42599 
    } else if (pheno=="uc") {
        n_ref<-30577 
    }


    # Define the paramaters to run GSMR:
    gwas_thresh = 5e-8    # GWAS threshold to select SNPs as the instruments for the GSMR analysis
    single_snp_heidi_thresh<-0.01    # p-value threshold for single-SNP-based HEIDI-outlier analysis
    multi_snps_heidi_thresh<-0.01    # p-value threshold for multi-SNP-based HEIDI-outlier analysis
    nsnps_thresh = 10   # the minimum number of instruments required for the GSMR analysis
    heidi_outlier_flag = T    # flag for HEIDI-outlier analysis
    ld_r2_thresh = 0.05    # LD r2 threshold to remove SNPs in high LD
    ld_fdr_thresh = 0.05   # FDR threshold to remove the chance correlations between the SNP instruments
    gsmr2_beta = 1     # 0 - the original HEIDI-outlier method; 1 - the new global HEIDI-outlier method 

    if (nrow(gsmr_data_ibd_traits[which(gsmr_data_ibd_traits$'P-value'<=p_threshold),])<10) {

        print("Not enough variants after pvalue thresholds applied")

    } else {

        gsmr_results<-gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, single_snp_heidi_thresh, multi_snps_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)    # GSMR analysis 
        filtered_index<-gsmr_results$used_index

        # extract the results:
        data<-as.data.frame(matrix(ncol=8,nrow=1))
        colnames(data)<-c("dataset","bxy","bxy_se","bxy_pval","index_snps","pleio_snps","pvalue_on_risk_factor","pvalue_on_disease")
        data$dataset<-dataset
        data$bxy<-gsmr_results$bxy
        data$bxy_se<-gsmr_results$bxy_se
        data$bxy_pval<-gsmr_results$bxy_pval
        data$index_snps<-paste(gsmr_data$SNP[gsmr_results$used_index],collapse="|")
        data$pleio_snps<-paste(gsmr_results$pleio_snps,collapse="|")

        data$pvalue_on_risk_factor<-paste(gsmr_data_ibd_traits$'P-value'[which(gsmr_data$SNP %in% gsmr_data$SNP[gsmr_results$used_index])],collapse="|")

        # extract pvalue of the SNPs on the disease:
        data$pvalue_on_disease<-paste(gsmr_data$'P-value'[which(gsmr_data$SNP %in% gsmr_data$SNP[gsmr_results$used_index])],collapse="|")

        # # retrieve the list of SNPs used in the analysis - after excluding the outliers:
        list_snps<-gsmr_data$SNP[gsmr_results$used_index]


        effect_col = colors()[75]
        # vals_all<-max(abs(c(bzx[filtered_index]-bzx_se[filtered_index], bzx[filtered_index]+bzx_se[filtered_index],bzy[filtered_index]-bzy_se[filtered_index], bzy[filtered_index]+bzy_se[filtered_index])))
        # vals = c(bzx[filtered_index]-bzx_se[filtered_index], bzx[filtered_index]+bzx_se[filtered_index])
        # xmin = min(vals); xmax = max(vals)
        # vals = c(bzy[filtered_index]-bzy_se[filtered_index], bzy[filtered_index]+bzy_se[filtered_index])
        # ymin = min(vals); ymax = max(vals)

        # xmin = vals_all*-1; xmax = vals_all
        # ymin = vals_all*-1; ymax = vals_all
        # par(mar=c(5,5,4,2))

        summary_data<-cbind(bzx[filtered_index],bzy[filtered_index],bzx_se[filtered_index],bzy_se[filtered_index])
        colnames(summary_data)<-c("bzx","bzy","bzx_se","bzy_se")
        summary_data<-as.data.frame(summary_data)


        p<-ggplot(summary_data,aes(x=bzx,y=bzy)) + 
        geom_point(size = 3) + xlab(paste("Beta ",toupper(pheno)," susceptibility (bzx)",sep="")) + 
        geom_abline(slope=data$bxy,colour="blue",lty=2) + 
        ylab(paste("Beta ",dataset," (bzy)",sep=""))  + 
        scale_color_brewer(palette = "Paired") + 
        # ylim(-1,1) + xlim(-0.1,0.1)+ 
        geom_errorbar(aes(ymin=bzy-bzy_se, ymax=bzy+bzy_se)) +  
        geom_errorbar(aes(xmin=bzx-bzx_se, xmax=bzx+bzx_se)) +
        # facet_grid(. ~ source,scales = "free") + theme_bw() +
        theme(legend.title=element_blank(),
        axis.text=element_text(size=12),axis.title=element_text(size=12),strip.text.x = element_text(size = 12)) +
        ggtitle(paste(toupper(pheno),"\nbxy = ",format(gsmr_results$bxy, scientific=F,digits=2),
                "; P-value =",format(gsmr_results$bxy_pval, scientific=TRUE,digits=3),"\nN instruments =",format(length(gsmr_results$used_index), scientific=F,digits=2)))


        ggsave(
            paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/plots/gsmr_eur_tier2_",dataset,"_",pheno,"_",p_threshold,"_reverse_unsupervised.png",sep=""),
            p,
            width = 45,
            height = 45,
            dpi = 400,
            units = c("mm"),
            limitsize = T,scale=4
        )

        # nsnps = length(bzx[filtered_index])
        # for( i in 1:nsnps ) {
        #     # x axis
        #     xstart = bzx[filtered_index [i]] - bzx_se[filtered_index[i]]; xend = bzx[filtered_index[i]] + bzx_se[filtered_index[i]]
        #     ystart = bzy[filtered_index[i]]; yend = bzy[filtered_index[i]]
        #     segments(xstart, ystart, xend, yend, lwd=1.5, col=effect_col)
        #     # y axis
        #     xstart = bzx[filtered_index[i]]; xend = bzx[filtered_index[i]] 
        #     ystart = bzy[filtered_index[i]] - bzy_se[filtered_index[i]]; yend = bzy[filtered_index[i]] + bzy_se[filtered_index[i]]
        #     segments(xstart, ystart, xend, yend, lwd=1.5, col=effect_col)
        # }

        # save main results:
        fwrite(data,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/gsmr_eur_tier2_",dataset,"_",pheno,"_",p_threshold,"_reverse_unsupervised.tsv.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

        summary_data$SNP<-list_snps
        
        # save results for the plot
        fwrite(summary_data,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/plots/gsmr_eur_tier2_",dataset,"_",pheno,"_",p_threshold,"_data_source_reverse_unsupervised.tsv.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

    }

} else {
    print(paste("Not enough instruments at pvalue",p_threshold))
}


q("no")
