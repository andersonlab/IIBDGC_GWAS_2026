# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# singularity exec iibdgc_postprocess_10_singularity.sif
path_gwas="/path/to/ibdgwas/IIBDGC/"

MEM=1000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(ggplot2)
library(colorspace)
library(dplyr)
library(tibble)
library(ggrepel)
library(ggpubr)
library(qvalue)

path_gwas="/path/to/ibdgwas/IIBDGC/"

files<-list.files(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/",sep=""))

# pvalue_threshold<-0.05/30

# relabel to phenotypes the study accessions
dat<-fread(paste(path_gwas,"resources/gwas_summary_statistics/gwas-catalog-v1.0.3.1-studies-r2024-06-17.tsv.gz",sep=""),quote="")
dat$study_accession<-dat$'STUDY ACCESSION'
dat$disease_trait<-dat$'DISEASE/TRAIT'

## reverse supervised results:
files_superv<-files[grep("_reverse_supervised_with_without_heidi.tsv.gz",files)]

for (i in 1:length(files_superv)){
    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/",files_superv[i],sep=""),head=T)
    tmp$pheno<-gsub("gsmr_eur_tier2_GCST[0-9]*_","",files_superv[i])
    tmp$pheno<-toupper(gsub("_5e-08_.*_with_without_heidi.tsv.gz","",tmp$pheno))

    if(i==1){
        superv<-tmp
    } else {
        superv<-rbind(superv,tmp) 
    }
}

superv<-merge(superv,dat[,c("study_accession","disease_trait")],by.x="dataset",by.y="study_accession")
superv$method<-"reverse_supervised"

## reverse unsupervised results:
files_unsuperv<-files[grep("_reverse_unsupervised_with_without_heidi.tsv.gz",files)]

for (i in 1:length(files_unsuperv)){
    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/",files_unsuperv[i],sep=""),head=T)
    tmp$pheno<-gsub("gsmr_eur_tier2_GCST[0-9]*_","",files_superv[i])
    tmp$pheno<-toupper(gsub("_5e-08_.*_with_without_heidi.tsv.gz","",tmp$pheno))

    if(i==1){
        unsuperv<-tmp
    } else {
        unsuperv<-rbind(unsuperv,tmp) 
    }
}

unsuperv<-merge(unsuperv,dat[,c("study_accession","disease_trait")],by.x="dataset",by.y="study_accession")
unsuperv$method<-"reverse_unsupervised"

## forward results:
files_forw<-files[grep("_5e-08_with_without_heidi.tsv.gz",files)]

for (i in 1:length(files_forw)){

    file_tmp<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/",files_forw[i],sep="")
    
    if(file.exists(file_tmp)) {
        tmp<-fread(file_tmp,head=T)
        tmp$pheno<-gsub("gsmr_eur_tier2_GCST[0-9]*_","",files_forw[i])
        tmp$pheno<-toupper(gsub("_5e-08_with_without_heidi.tsv.gz","",tmp$pheno))

        if(i==1){
            forw<-tmp
        } else {
            forw<-rbind(forw,tmp) 
        }
        rm(tmp)

    }

}

forw<-merge(forw,dat[,c("study_accession","disease_trait")],by.x="dataset",by.y="study_accession")
forw$method<-"forward"

# get the list  of significant GSMR in at least one approach:
all<-rbind(forw,unsuperv,superv)
dim(all)
# [1] 270  11


# select only the pairs that are significantly correlated in genetic correlation step - see 
# /path/to/user/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_45_compute_genetic_correlation.R

list_signif<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/list_pairs_traits_significantly_correlated_with_ibd_cd_uc",sep=""),head=F)
colnames(list_signif)<-c("dataset","disease_trait","pheno")
list_signif
list_signif$pair<-paste(list_signif$dataset,toupper(list_signif$pheno),sep="_")
dim(list_signif)
# [1] 42  4

all$pair<-paste(all$dataset,toupper(all$pheno),sep="_")
all<-all[which(all$pair %in% list_signif$pair),]
dim(all)
# [1] 126  13

42*3
# [1] 126 - still one missing, still running

# define which MR are significant after controlling for false positives:
all$bxy_qval<-qvalue(all$bxy_pval)$qvalues

tmp<-all[,c("dataset","disease_trait","pheno")]
tmp<-tmp[!duplicated(tmp),]
dim(tmp)
# [1] 42  2

tmp$pheno<-tolower(tmp$pheno)

# define the set of limits:

min_lim1<-min(all$bxy[!is.na(all$bxy) & !is.na(all$bxy)]-all$bxy_se[!is.na(all$bxy) & !is.na(all$bxy)],na.rm=T)-0.1
max_lim1<-max(all$bxy[!is.na(all$bxy) & !is.na(all$bxy)]+all$bxy_se[!is.na(all$bxy) & !is.na(all$bxy)],na.rm=T)+0.1


tmp$forward_supervised_pvalue<-NA
tmp$forward_supervised_bxy<-NA
tmp$forward_supervised_sensitivity_pass<-NA
tmp$forward_unsupervised_pvalue<-NA
tmp$forward_unsupervised_bxy<-NA
tmp$forward_unsupervised_sensitivity_pass<-NA
tmp$reverse_unsupervised_pvalue<-NA
tmp$reverse_unsupervised_bxy<-NA
tmp$reverse_unsupervised_sensitivity_pass<-NA


# how to plot all together:

for (i in 1:nrow(tmp)) {
        
    ph<-tolower(tmp$pheno[i])

    # per SNP summary files
    file_sum_superv<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/plots/gsmr_eur_tier2_",tmp$dataset[i],"_",ph,"_5e-08_data_source_reverse_supervised.tsv.gz",sep="")
    file_sum_unsuperv<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/plots/gsmr_eur_tier2_",tmp$dataset[i],"_",ph,"_5e-08_data_source_reverse_unsupervised.tsv.gz",sep="")
    file_sum_forw<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/plots/gsmr_eur_tier2_",tmp$dataset[i],"_",ph,"_5e-08_data_source.tsv.gz",sep="")

    # GSMR results
    file_sum_superv_2<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/gsmr_eur_tier2_",tmp$dataset[i],"_",ph,"_5e-08_reverse_supervised_with_without_heidi.tsv.gz",sep="")
    file_sum_unsuperv_2<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/gsmr_eur_tier2_",tmp$dataset[i],"_",ph,"_5e-08_reverse_unsupervised_with_without_heidi.tsv.gz",sep="")
    file_sum_forw_2<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/gsmr_eur_tier2_",tmp$dataset[i],"_",ph,"_5e-08_with_without_heidi.tsv.gz",sep="")

    # sensitivity analyses
    file_sum_superv_3<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/gsmr_eur_tier2_",tmp$dataset[i],"_",ph,"_5e-08_reverse_supervised_sensitivity_test.tsv.gz",sep="")
    file_sum_unsuperv_3<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/gsmr_eur_tier2_",tmp$dataset[i],"_",ph,"_5e-08_reverse_unsupervised_sensitivity_test.tsv.gz",sep="")
    file_sum_forw_3<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/gsmr_eur_tier2_",tmp$dataset[i],"_",ph,"_5e-08_sensitivity_test.tsv.gz",sep="")


    if ( file.exists(file_sum_superv) & file.exists(file_sum_unsuperv) & file.exists(file_sum_forw) & 
    file.exists(file_sum_superv_2) & file.exists(file_sum_unsuperv_2) & file.exists(file_sum_forw_2) & 
    file.exists(file_sum_superv_3) & file.exists(file_sum_unsuperv_3) & file.exists(file_sum_forw_3)) {

        # print("All data ready")
        # print(tmp[i,1:3],ph)
    
        sum_superv<-fread(file_sum_superv)
        sum_superv$method<-"reverse_supervised"
        sum_superv$dataset<-tmp$dataset[i]
        sum_superv$disease_trait<-tmp$disease_trait[i]
        sum_superv$bxy_pval<-all$bxy_pval[which(all$method=="reverse_supervised" & all$dataset==tmp$dataset[i] & all$pheno==toupper(ph))]
        sum_superv$bxy_qval<-all$bxy_qval[which(all$method=="reverse_supervised" & all$dataset==tmp$dataset[i] & all$pheno==toupper(ph))]
        sum_superv$bxy<-all$bxy[which(all$method=="reverse_supervised" & all$dataset==tmp$dataset[i] & all$pheno==toupper(ph))]

        model<-fread(file_sum_superv_2)
        pleiotropic_snps<-unlist(strsplit(model$pleio_snps,split="\\|"))

        # highliht pleiotropic SNPs:
        sum_superv$highlight<-"normal"
        sum_superv$highlight[which(sum_superv$SNP %in% pleiotropic_snps)]<-"highlight_grey"

        tmp$forward_supervised_pvalue[i]<-sum_superv$bxy_pval[1]
        tmp$forward_supervised_bxy[i]<-sum_superv$bxy[1]

        if (sum_superv$bxy_pval[1]<0.05) {

            colour_plot<-"black"

            # incorporate sensitivity analyses - leave one out from model with Heidi=T
            sa_superv<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/gsmr_eur_tier2_",tmp$dataset[i],"_",ph,"_5e-08_reverse_supervised_sensitivity_test.tsv.gz",sep=""))
            sa_superv<-sa_superv[which(sa_superv$gsmr_no_heidi_bxy_pvalue>0.05),]

            # highliht SNPs that do not pass sensitivity analysis:
            sum_superv$highlight[which(sum_superv$SNP %in% sa_superv$SNP_out)]<-"highlight_blue"

            if (nrow(sa_superv)>0) {
                tmp$forward_supervised_sensitivity_pass[i]<-"N"
            }
        

        } else {
            colour_plot<-"darkgrey"


            # incorporate sensitivity analyses - leave one out from model with Heidi=T
            sa_superv<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/gsmr_eur_tier2_",tmp$dataset[i],"_",ph,"_5e-08_reverse_supervised_sensitivity_test.tsv.gz",sep=""))
            sa_superv<-sa_superv[which(sa_superv$gsmr_no_heidi_bxy_pvalue<0.05),]

            # highliht SNPs that do not pass sensitivity analysis:
            sum_superv$highlight[which(sum_superv$SNP %in% sa_superv$SNP_out)]<-"highlight_red"

            if (nrow(sa_superv)>0) {
                tmp$forward_supervised_sensitivity_pass[i]<-"N"
            }

        }


        mycolours <- c("highlight_red" = "#ef8a62", "highlight_blue" = "#67a9cf","highlight_grey" = "lightgrey", "normal" = colour_plot)

        p1<-ggplot(sum_superv,aes(x=bzx,y=bzy)) + 
        geom_point(size = 3,aes(colour = highlight)) + xlab(paste("Beta ",toupper(ph)," susceptibility (bzx)",sep="")) + 
        geom_abline(slope=sum_superv$bxy[1],colour="blue",lty=2) + 
        ylab(paste("Beta ",tmp$disease_trait[i]," (bzy)",sep=""))  + 
        # scale_color_brewer(palette = "Paired") + 
        # ylim(-1,1) + xlim(-0.1,0.1)+ 
        geom_errorbar(aes(ymin=bzy-bzy_se, ymax=bzy+bzy_se, colour = highlight)) +  
        geom_errorbar(aes(xmin=bzx-bzx_se, xmax=bzx+bzx_se, colour = highlight)) +
        # facet_grid(. ~ source,scales = "free") + theme_bw() +
        scale_color_manual(values = mycolours) +
        theme(legend.title=element_blank(),legend.position="none",
        axis.text=element_text(size=12),axis.title=element_text(size=12),strip.text.x = element_text(size = 12)) +
        ggtitle(paste("Risk factor ",toupper(ph), ", on ",tmp$disease_trait[i],", supervised",
        "\nbxy = ",format(sum_superv$bxy[2], scientific=F,digits=2),
                "; P-value =",format(sum_superv$bxy_pval[1], scientific=TRUE,digits=3),
                # "; q-value =",format(sum_superv$bxy_qval[1], scientific=TRUE,digits=3),
                "\nN instruments = ",format(length(sum_superv$SNP)-length(pleiotropic_snps), scientific=F,digits=2),
                "; N pleitropic = ",format(length(pleiotropic_snps), scientific=F,digits=2),
                "\nP-value (no HEIDI) =",format(model$gsmr_no_heidi_bxy_pvalue[1], scientific=TRUE,digits=3),
                sep=""))


        ##############################

        sum_unsuperv<-fread(file_sum_unsuperv)
        sum_unsuperv$method<-"reverse_unsupervised"
        sum_unsuperv$dataset<-tmp$dataset[i]
        sum_unsuperv$disease_trait<-tmp$disease_trait[i]
        sum_unsuperv$bxy_pval<-all$bxy_pval[which(all$method=="reverse_unsupervised" & all$dataset==tmp$dataset[i] & all$pheno==toupper(ph))]
        sum_unsuperv$bxy_qval<-all$bxy_qval[which(all$method=="reverse_unsupervised" & all$dataset==tmp$dataset[i] & all$pheno==toupper(ph))]
        sum_unsuperv$bxy<-all$bxy[which(all$method=="reverse_unsupervised" & all$dataset==tmp$dataset[i] & all$pheno==toupper(ph))]

        model<-fread(file_sum_unsuperv_2)
        pleiotropic_snps<-unlist(strsplit(model$pleio_snps,split="\\|"))

        # highliht pleiotropic SNPs:
        sum_unsuperv$highlight<-"normal"
        sum_unsuperv$highlight[which(sum_unsuperv$SNP %in% pleiotropic_snps)]<-"highlight_grey"

        tmp$forward_unsupervised_pvalue[i]<-sum_unsuperv$bxy_pval[1]
        tmp$forward_unsupervised_bxy[i]<-sum_unsuperv$bxy[1]


        if (sum_unsuperv$bxy_pval[1]<0.05) {

            colour_plot<-"black"
                
            # incorporate sensitivity analyses - leave one out
            sa_unsuperv<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/gsmr_eur_tier2_",tmp$dataset[i],"_",ph,"_5e-08_reverse_unsupervised_sensitivity_test.tsv.gz",sep=""))
            sa_unsuperv<-sa_unsuperv[which(sa_unsuperv$gsmr_no_heidi_bxy_pvalue>0.05),]
                
            # highliht SNPs that do not pass sensitivity analysis:
            sum_unsuperv$highlight[which(sum_unsuperv$SNP %in% sa_unsuperv$SNP_out)]<-"highlight_blue"

            if (nrow(sa_unsuperv)>0) {
                tmp$forward_unsupervised_sensitivity_pass[i]<-"N"
            }

        } else {

            colour_plot<-"darkgrey"
                
            # incorporate sensitivity analyses - leave one out
            sa_unsuperv<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/gsmr_eur_tier2_",tmp$dataset[i],"_",ph,"_5e-08_reverse_unsupervised_sensitivity_test.tsv.gz",sep=""))
            sa_unsuperv<-sa_unsuperv[which(sa_unsuperv$gsmr_no_heidi_bxy_pvalue<0.05),]
                
            # highliht SNPs that do not pass sensitivity analysis:
            sum_unsuperv$highlight[which(sum_unsuperv$SNP %in% sa_unsuperv$SNP_out)]<-"highlight_red"

            if (nrow(sa_unsuperv)>0) {
                tmp$forward_unsupervised_sensitivity_pass[i]<-"N"
            }

        }

        mycolours <- c("highlight_red" = "#ef8a62", "highlight_blue" = "#67a9cf","highlight_grey" = "lightgrey", "normal" = colour_plot)

        p2<-ggplot(sum_unsuperv,aes(x=bzx,y=bzy)) + 
        geom_point(size = 3,aes(colour = highlight)) + xlab(paste("Beta ",toupper(ph)," susceptibility (bzx)",sep="")) + 
        geom_abline(slope=sum_unsuperv$bxy[1],colour="blue",lty=2) + 
        ylab(paste("Beta ",tmp$disease_trait[i]," (bzy)",sep=""))  + 
        # scale_color_brewer(palette = "Paired") + 
        # ylim(-1,1) + xlim(-0.1,0.1)+ 
        geom_errorbar(aes(ymin=bzy-bzy_se, ymax=bzy+bzy_se, colour = highlight)) +  
        geom_errorbar(aes(xmin=bzx-bzx_se, xmax=bzx+bzx_se, colour = highlight)) +
        # facet_grid(. ~ source,scales = "free") + theme_bw() +
        scale_color_manual(values = mycolours) +
        theme(legend.title=element_blank(),legend.position="none",
        axis.text=element_text(size=12),axis.title=element_text(size=12),strip.text.x = element_text(size = 12)) +
        ggtitle(paste("Risk factor ",toupper(ph), ", on ",tmp$disease_trait[i],", unsupervised",
                "\nbxy = ",format(sum_unsuperv$bxy[1], scientific=F,digits=2),
                "; P-value =",format(sum_unsuperv$bxy_pval[1], scientific=TRUE,digits=3),
                # "; q-value =",format(sum_unsuperv$bxy_qval[1], scientific=TRUE,digits=3),
                "\nN instruments = ",format(length(sum_unsuperv$SNP)-length(pleiotropic_snps), scientific=F,digits=2),
                "; N pleitropic = ",format(length(pleiotropic_snps), scientific=F,digits=2),
                "\nP-value (no HEIDI) =",format(model$gsmr_no_heidi_bxy_pvalue[1], scientific=TRUE,digits=3),
                sep=""))


        ##############################

        sum_forw<-fread(file_sum_forw)
        sum_forw$method<-"forward"
        sum_forw$dataset<-tmp$dataset[i]
        sum_forw$disease_trait<-tmp$disease_trait[i]
        sum_forw$bxy_pval<-all$bxy_pval[which(all$method=="forward" & all$dataset==tmp$dataset[i] & all$pheno==toupper(ph))]
        sum_forw$bxy_qval<-all$bxy_qval[which(all$method=="forward" & all$dataset==tmp$dataset[i] & all$pheno==toupper(ph))]
        sum_forw$bxy<-all$bxy[which(all$method=="forward" & all$dataset==tmp$dataset[i] & all$pheno==toupper(ph))]

        model<-fread(file_sum_forw_2)
        pleiotropic_snps<-unlist(strsplit(model$pleio_snps,split="\\|"))

        # highliht pleiotropic SNPs:
        sum_forw$highlight<-"normal"
        sum_forw$highlight[which(sum_forw$SNP %in% pleiotropic_snps)]<-"highlight_grey"

        tmp$reverse_unsupervised_pvalue[i]<-sum_forw$bxy_pval[1]
        tmp$reverse_unsupervised_bxy[i]<-sum_forw$bxy[1]

        if (sum_forw$bxy_pval[1]<0.05) {

            colour_plot<-"black"
                
            # incorporate sensitivity analyses - leave one out
            sa_forw<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/gsmr_eur_tier2_",tmp$dataset[i],"_",ph,"_5e-08_sensitivity_test.tsv.gz",sep=""))
            sa_forw<-sa_forw[which(sa_forw$gsmr_no_heidi_bxy_pvalue>0.05),]
                
            # highliht SNPs that do not pass sensitivity analysis:
            sum_forw$highlight[which(sum_forw$SNP %in% sa_forw$SNP_out)]<-"highlight_blue"

            if (nrow(sa_forw)>0) {
                tmp$forward_unsupervised_sensitivity_pass[i]<-"N"
            }

        } else {

            colour_plot<-"darkgrey"
                
            # incorporate sensitivity analyses - leave one out
            sa_forw<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/mendelian_randomization/results/gsmr_eur_tier2_",tmp$dataset[i],"_",ph,"_5e-08_sensitivity_test.tsv.gz",sep=""))
            sa_forw<-sa_forw[which(sa_forw$gsmr_no_heidi_bxy_pvalue<0.05),]
                
            # highliht SNPs that do not pass sensitivity analysis:
            sum_forw$highlight[which(sum_forw$SNP %in% sa_forw$SNP_out)]<-"highlight_red"

            if (nrow(sa_forw)>0) {
                tmp$forward_unsupervised_sensitivity_pass[i]<-"N"
            }

        }

        mycolours <- c("highlight_red" = "#ef8a62", "highlight_blue" = "#67a9cf","highlight_grey" = "lightgrey", "normal" = colour_plot)

        p3<-ggplot(sum_forw,aes(x=bzx,y=bzy)) + 
        geom_point(size = 3,aes(colour = highlight)) + ylab(paste("Beta ",toupper(ph)," susceptibility (bzy)",sep="")) + 
        geom_abline(slope=sum_forw$bxy[1],colour="blue",lty=2) + 
        xlab(paste("Beta ",tmp$disease_trait[i]," (bzx)",sep=""))  + 
        # scale_color_brewer(palette = "Paired") + 
        # ylim(-1,1) + xlim(-0.1,0.1)+ 
        geom_errorbar(aes(ymin=bzy-bzy_se, ymax=bzy+bzy_se, colour = highlight)) +  
        geom_errorbar(aes(xmin=bzx-bzx_se, xmax=bzx+bzx_se, colour = highlight)) +
        # facet_grid(. ~ source,scales = "free") + theme_bw() +
        scale_color_manual(values = mycolours) +
        theme(legend.title=element_blank(),legend.position="none",
        axis.text=element_text(size=12),axis.title=element_text(size=12),strip.text.x = element_text(size = 12)) +
        ggtitle(paste("Risk factor ",tmp$disease_trait[i], ", on ",toupper(ph),", unsupervised",
                "\nbxy = ",format(sum_forw$bxy[1], scientific=F,digits=2),
                "; P-value =",format(sum_forw$bxy_pval[1], scientific=TRUE,digits=3),
                # "; q-value =",format(sum_forw$bxy_qval[1], scientific=TRUE,digits=3),
                "\nN instruments = ",format(length(sum_forw$SNP)-length(pleiotropic_snps), scientific=F,digits=2),
                "; N pleitropic = ",format(length(pleiotropic_snps), scientific=F,digits=2),
                "\nP-value (no HEIDI) =",format(model$gsmr_no_heidi_bxy_pvalue[1], scientific=TRUE,digits=3),
                sep=""))

        p<-ggarrange(p2,p1,p3,nrow=2,ncol=2)

        trait<-tolower(gsub(" ","_",tmp$disease_trait[i]))

        ggsave(
            paste("~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_",tmp$dataset[i],"_",trait,"_",ph,"_gsmr.pdf",sep=""),
            p,
            width = 180,
            height = 180,
            dpi = 300,
            units = c("mm"),
            limitsize = T,scale=2
        )

        rm(p,p1,p2,p3,sum_superv,sum_unsuperv,sum_forw,sa_superv,sa_unsuperv,sa_forw)
            
    } else {
        print("Not all data ready")
        print(tmp[i,1:3],ph)
    }

    rm(file_sum_superv,file_sum_unsuperv,file_sum_forw,file_sum_superv_2,file_sum_unsuperv_2,
    file_sum_forw_2,file_sum_superv_3,file_sum_unsuperv_3,file_sum_forw_3)

}


# [1] "Not all data ready"
#         dataset               disease_trait  pheno
#          <char>                      <char> <char>
# 1: GCST90002322 Mean corpuscular hemoglobin     cd
# [1] "Not all data ready"
#         dataset               disease_trait  pheno
#          <char>                      <char> <char>
# 1: GCST90002322 Mean corpuscular hemoglobin    ibd
# [1] "Not all data ready"
#         dataset           disease_trait  pheno
#          <char>                  <char> <char>
# 1: GCST90002334 Mean corpuscular volume     cd
# [1] "Not all data ready"
#         dataset           disease_trait  pheno
#          <char>                  <char> <char>
# 1: GCST90002334 Mean corpuscular volume    ibd
# [1] "Not all data ready"
#         dataset        disease_trait  pheno
#          <char>               <char> <char>
# 1: GCST90002363 Red blood cell count    ibd
# [1] "Not all data ready"
#         dataset        disease_trait  pheno
#          <char>               <char> <char>
# 1: GCST90002363 Red blood cell count     uc
# [1] "Not all data ready"
#         dataset     disease_trait  pheno
#          <char>            <char> <char>
# 1: GCST90129505 Colorectal cancer     cd


data1<-tmp[,c("dataset","disease_trait","pheno","forward_supervised_pvalue","forward_supervised_bxy","forward_supervised_sensitivity_pass")]
data1$mr<-"forward_supervised"
data2<-tmp[,c("dataset","disease_trait","pheno","forward_unsupervised_pvalue","forward_unsupervised_bxy","forward_unsupervised_sensitivity_pass")]
data2$mr<-"forward_unsupervised"
data3<-tmp[,c("dataset","disease_trait","pheno","reverse_unsupervised_pvalue","reverse_unsupervised_bxy","reverse_unsupervised_sensitivity_pass")]
data3$mr<-"reverse_unsupervised"

colnames(data1)<-c("dataset","disease_trait","pheno","pvalue","bxy","sensitivity_pass","mr")
colnames(data2)<-c("dataset","disease_trait","pheno","pvalue","bxy","sensitivity_pass","mr")
colnames(data3)<-c("dataset","disease_trait","pheno","pvalue","bxy","sensitivity_pass","mr")

data<-rbind(data1,data2,data3)

data$sensitivity_pass[which(!is.na(data$pvalue) & is.na(data$sensitivity_pass))]<-"Y"
data<-data[which(!is.na(data$pvalue)),]

ggplot(data, aes(x=mr, y=disease_trait, fill=bxy)) + 
  geom_tile() + facet_grid( ~ pheno, scales = "free",space = "free") + scale_fill_distiller(palette = "Spectral", direction = 1) + 
  scale_colour_gradient2(limits = c(-4,4),midpoint = 0,low ="yellow",mid = "white",high = "blue")




