# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# To launch interactively:
# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM=8000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q cpu-interactive R

library(data.table)
library(rtracklayer)
library(stringr)
library(qvalue)
library(ggplot2)
library(ggpubr)

rm(list=ls())

pheno<-c("ibd","cd","uc")

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
path_eqtl<-"/path/to/project"


# LOAD THE MAP FOR COLOCALIZATION RESULTS - SEE ~/git/snakemake_colocalisation/scripts/create_maps_for_IIBDGC.R
map<-fread("/path/to/project",head=T)
map<-as.data.frame(map)
map<-map[which(map$cohort!="edQTL_GTEX"),]
map<-map[which(map$quant_method %in% c("aptamer","ge","microarray")),]

# LOAD THE INDEPENDENT SIGNALS - see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_38_2_define_independent_signals_using_as_ref_cojo_agnostic_results.R
all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 619 192

# LOAD THE COLOCALIZATION RESULTS - SEE ~/git/IIBDGC_GWAS/scripts/other/evaluate_direction_effect_QTL_risk_variant.R

for (chr in c(1:22)) {
    for (ph in pheno) {
        tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/tmp_coloc_filtered_by_gwas_r2_0.6_allpheno_with_beta_with_wald_ratio_",chr,"_",ph,".tsv"),head=T,sep="\t")
        if (ph=="ibd" & chr==1) {
            dat<-tmp
        }else{
            dat<-rbind(dat,tmp,fill=T)
        }
    }
}

dat<-as.data.frame(dat)
dat$se_wald<-dat$sd_wald_ratio/dat$n_variants_in_ld_for_wr
dat$ci_low_wald_ratio<-NA
dat$ci_up_wald_ratio<-NA

betas_list <- strsplit(dat$betas_wald_ratio, "\\|")
dat$ci_low_wald_ratio <- vapply(betas_list, function(x) quantile(as.numeric(x), 0.025, na.rm=T), numeric(1))
dat$ci_up_wald_ratio  <- vapply(betas_list, function(x) quantile(as.numeric(x), 0.975, na.rm=T), numeric(1))

table(dat$cohort[!grepl("ENSG",dat$'gene')])
# pQTL_decode 
#          27

dat$gene[!grepl("ENSG",dat$gene)]
dat$gene_name[which(dat$gene_name=="")]<-dat$gene[which(dat$gene_name=="")]

table(dat$pheno_coloc)
#   cd  ibd   uc 
# 6130 6969 5377

coloc_results<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_filtered_by_gwas_r2_0.6_allpheno_noFM_with_beta.tsv",sep=""))
coloc_results<-as.data.frame(coloc_results)

dat$X<-paste(dat$condition_name,dat$gene,dat$MarkerName_eqtl_lead,sep="_")
coloc_results$X<-paste(coloc_results$condition_name,coloc_results$gene,coloc_results$MarkerName_eqtl_lead,sep="_")

dim(dat[which(!dat$X %in% coloc_results$X),])
# 0

dim(coloc_results[which(!coloc_results$X %in% dat$X),])
# 0

length(names(table(dat$gene_name)))
# [1] 625


# for each gene, plot the wald_test results:
summary(dat$n_variants_in_ld_for_wr)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1.00   15.00   42.00   70.85   99.00  516.00

list_genes<-names(table(dat$gene_name))
length(list_genes)
# [1] 625

dat<-merge(dat[,which(!colnames(dat) %in% colnames(map))],map,by.x="condition_name",by.y="id_map")
dim(dat)
# [1] 18476    30

fwrite(dat,"~/git/IIBDGC_GWAS/plots/paper_tables/Source_data_coloc_results_beta_wald.tsv.gz",col.names=T,row.names=F,quote=F,sep="\t")

dat$direction_wald<-NA
dat$direction_wald[which(dat$mean_wald_ratio<0)]<-"negative"
dat$direction_wald[which(dat$mean_wald_ratio>0)]<-"positive"


xx<-as.data.frame.matrix(table(dat$gene_name,dat$direction_wald))
dim(xx)
# [1] 625   3


for (gene in list_genes) {

    gene_dat<-dat[which(dat$gene_name==gene),]
    gene_dat$direction_wald<-factor(gene_dat$direction_wald,levels=c("negative","positive"))

    variants<-names(table(gene_dat$IIBDGC_GWAS_index_variant))
    plot_list<-vector("list", length(variants))

    for (i in seq_along(variants)) {

        var<-variants[i]

        tmp<-gene_dat[which(gene_dat$IIBDGC_GWAS_index_variant==var),]
        tmp$direction_wald<-factor(tmp$direction_wald,levels=c("negative","positive"))
        tmp$tissue_cell_condition_label<-paste(tmp$study_label.x,tmp$tissue_cell_condition_label)

        xlim_abs<-max(abs(tmp$ci_low_wald_ratio), abs(tmp$ci_up_wald_ratio), na.rm=T)
        xlim<-c(round(xlim_abs)*-1-10, round(xlim_abs)+10)

        plot_list[[i]]<-ggplot(tmp, aes(x=mean_wald_ratio,y=tissue_cell_condition_label,color=direction_wald)) +
        geom_point() + theme(axis.title.y=element_blank(),) +
        xlim(xlim) +
        geom_errorbar(aes(xmin=ci_low_wald_ratio, xmax=ci_up_wald_ratio), width=0.2) +
        scale_size_continuous(limits=c(0,520), breaks=seq(0,520, by=40),name="Number variants (r2≥0.6)") +
        theme(plot.margin = unit(c(1,1,1,0.5), "cm"))+ geom_vline(xintercept=c(0), linetype="dotted")+
        # scale_color_manual(values=c("red","blue"),breaks = c("negative", "positive"),labels = c("Opposite direction", "Same direction"),drop = F)
        scale_color_manual(values=c("blue","red"),drop=F) +
        guides(color="none", size=guide_legend()) + ggtitle(paste(gene, var)) + facet_grid(cell_label~ pheno_coloc,scales = "free",, space = "free") +
        labs(x="Mean Beta-wald (95% quantile)") + theme(strip.text.y.right = element_text(angle = 0))

    }

    n_labels<-length(levels(as.factor(gene_dat$tissue_cell_condition_label)))
    hei<-ceiling(n_labels*4)+20

    if (n_labels<200) {

        if (length(variants)>1) {
            hei<-hei*length(variants)
            p<-ggarrange(plotlist=plot_list, nrow=length(variants))
        } else {
            p<-plot_list[[1]]
        }

        ggsave(
        paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/plots/",gene,"_coloc_direction_effect.png"),
        p,
        width = 120,
        height = hei,
        dpi = 400,
        units = c("mm"),
        limitsize = T,scale=2
        )

    } else {
        print(paste("Plot too large:",gene))
    }

    rm(p, gene_dat, plot_list)
}
# [1] "Plot too large: ERAP2"
# [1] "Plot too large: RHD"

xx$list_genes<-rownames(xx)


# number of genes with same / opposite direction effects:
dim(xx)
# [1] 625   3

# evidence from only one QTL
dim(xx[which( (xx$negative==1 & xx$positive==0) | (xx$negative==0 & xx$positive==1)),])
# [1] 21  2

# nominated by >1 QTL
625-21
# [1] 604

# opposite directions of effect
dim(xx[which(xx$negative!=0 & xx$positive!=0),]) 
# [1] 78  2

604-78
# [1] 526 # number of coloc wiht > 1 line of evidence and in the same direction

## annotate the list of genes with this information:

list_genes<-fread(paste0(path_gwas,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence.tsv.gz"))
list_genes<-merge(list_genes,xx,by="list_genes",all.x=T)

# update the list
write.table(list_genes,paste0(path_gwas,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence_with_coloc_direction.tsv.gz"),col.names=T,row.names=F,quote=F,sep="\t")

xx

q("no")

# cd /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/coloc_results/plots/
# tar -czvf plot_coloc_results_624.tar.gz *png

