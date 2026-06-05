# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
#############################################################################################

# takes output from step_40_add_results_from_coloc_and_MR.R

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=65000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(stringr)
library(rtracklayer)
library(ggh4x)
library(qvalue)
library(colorspace)
library(R.utils)

rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
path<-"/path/to/project"

# only multithreaded when writing
setDTthreads(2)
path_gwas<-"/path/to/ibdgwas/IIBDGC/"

args = commandArgs(trailingOnly=TRUE)
id<-args[1]

gene_symbol<-read.table(paste(path_gwas,"post_imputation/2022/analysis/colocalization/list_coloc_genes_pph4_80_feb24.tsv.gz",sep=""),head=T)
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="TSPO")]
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="SELL")]
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="WNK1")]
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="NCF4")]
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol==id)]
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="BCL2")]
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="VSIR")]
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="BCL2L11")]
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="LILRA6")]
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="CD200R1")]

gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="VSIR")]
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="ITGAV")]
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="ITGA4")]
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="ITGAL")]
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="ICAM5")]
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="PLAUR")]
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="CD28")]
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="MMP10")]
# gene_symbol<-gene_symbol$gene_symbol[which(gene_symbol$gene_symbol=="LY75")]

coloc<-fread(paste(path_gwas,"post_imputation/2022/analysis/colocalization/list_coloc_pph4_80_feb24.tsv.gz",sep=""),head=T)
coloc<-as.data.frame(coloc)
coloc<-coloc[which(coloc$gene_symbol %in% gene_symbol),]

# coloc<-coloc[which(coloc$study_label %in% c("CAP","GEUVADIS")),]


for (i in 1:nrow(coloc)) {

    if (coloc$study_label[i]=="Yazar_2022") {
        tmp<-fread(paste(path,"QTL_managed_access/yazar_2022/OneK1K_matrix_eQTL_results/nominal/",coloc$condition_name[i],"/",coloc$condition_name[i],".tsv.gz",sep=""),head=F)
        colnames(tmp)<-c("variant","r2","pvalue","molecular_trait_object_id","molecular_trait_id","maf",
                "gene_id","median_tpm","beta","se","an","ac","chromosome","position","ref","alt","type","rsid")


    } else if (coloc$study_label[i]=="Hu_2021") {
        tmp<-fread(paste(path,"QTL_managed_access/hu_2021/nominal/",coloc$condition_name[i],"/",coloc$condition_name[i],".tsv.gz",sep=""),head=F)
        colnames(tmp)<-c("variant","r2","pvalue","molecular_trait_object_id","molecular_trait_id","maf",
            "gene_id","median_tpm","beta","se","an","ac","chromosome","position","ref","alt","type","rsid")

    } else if (coloc$study_label[i]=="Panousis_2024") {
        tmp<-fread(paste(path,"macromap/nominal/",coloc$condition_name[i],"/",coloc$condition_name[i],".tsv.gz",sep=""),head=F)
        colnames(tmp)<-c("variant","r2","pvalue","molecular_trait_object_id","molecular_trait_id","maf",
            "gene_id","median_tpm","beta","se","an","ac","chromosome","position","ref","alt","type","rsid")

    } else if (coloc$study_label[i]=="eQTLGen_2021") {
        tmp<-fread(paste(path,"eQTLGen/nominal/",coloc$condition_name[i],"/",coloc$condition_name[i],".tsv.gz",sep=""),head=F)
        colnames(tmp)<-c("variant","r2","pvalue","molecular_trait_object_id","molecular_trait_id","maf",
            "gene_id","median_tpm","beta","se","an","ac","chromosome","position","ref","alt","type","rsid")

    } else if (coloc$study_label[i]=="blueprint") {
        tmp<-fread(paste(path,"blueprint//nominal/",coloc$condition_name[i],"/",coloc$condition_name[i],".tsv.gz",sep=""),head=F)
        colnames(tmp)<-c("variant","r2","pvalue","molecular_trait_object_id","molecular_trait_id","maf",
            "gene_id","median_tpm","beta","se","an","ac","chromosome","position","ref","alt","type","rsid")

    } else {
        tmp<-fread(paste(path,"eQTL_catalog_V6//nominal/",coloc$condition_name[i],"/",coloc$condition_name[i],".tsv.gz",sep=""),head=F)
        colnames(tmp)<-c("variant","r2","pvalue","molecular_trait_object_id","molecular_trait_id","maf",
            "gene_id","median_tpm","beta","se","an","ac","chromosome","position","ref","alt","type","rsid")
    }


    tmp<-tmp[which(tmp$molecular_trait_object_id==coloc$phenotype_id[i]),]
    tmp$MarkerName<-gsub("_",":",tmp$variant)

    tmp<-tmp[,c("MarkerName","molecular_trait_object_id","beta","se","pvalue")]
    tmp$pvalue<-as.numeric(tmp$pvalue)

    tmp$qvalue<-qvalue(tmp$pvalue)$qvalue
    tmp<-tmp[order(tmp$qvalue,decreasing=F),]
    tmp<-tmp[which(tmp$qvalue<0.05),]

    tmp$gene_symbol<-coloc$gene_symbol[i]
    tmp$study_label<-coloc$study_label[i]
    tmp$sample_group<-coloc$sample_group[i]

    if(i==1) {
        eqtl<-tmp
    } else {
        eqtl<-rbind(eqtl,tmp)
    }
    rm(tmp)

}

risk<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",tolower(coloc$gwas_trait[1]),"/",coloc$chr[1],"_",tolower(coloc$gwas_trait[1]),"_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),head=T)
risk<-risk[,c("MarkerName","BETA","SE","P-value")]
colnames(risk)[4]<-"pvalue_disease"

eqtl<-as.data.frame(eqtl)

eqtl<-merge(eqtl,risk,by="MarkerName")
eqtl$beta<-as.numeric(eqtl$beta)

# # estimate ld between variants to define independent signals
# write.table(eqtl$MarkerName[!duplicated(eqtl$MarkerName)],paste(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/",gene_symbol,"_eqtl_risk_variants_coloc_pph4_80_feb24_for_ld",sep=""),
# col.names=F,row.names=F,quote=F,sep="\t")

# # estimate LD between those and index variants above:
# system(paste("plink2 --bfile ",path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/gsa_allchr_subset_included_in_",tolower(coloc$gwas_trait[1]),"_analysis --extract ",path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/",gene_symbol,"_eqtl_risk_variants_coloc_pph4_80_feb24_for_ld --make-bed --out ",path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/gsa_allchr_list_",gene_symbol,"_eqtl_risk_variants_coloc_pph4_80_feb24_for_ld",sep=""))
# system(paste("plink --noweb --bfile ",path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/gsa_allchr_list_",gene_symbol,"_eqtl_risk_variants_coloc_pph4_80_feb24_for_ld --r2 --ld-window-kb 10000000 --ld-window-r2 0 --out ",path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/gsa_allchr_list_",gene_symbol,"_eqtl_risk_variants_coloc_pph4_80_feb24_for_ld",sep=""))

# # keep the most significant variant per signal/study
# ld<-fread(paste(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/gsa_allchr_list_",gene_symbol,"_eqtl_risk_variants_coloc_pph4_80_feb24_for_ld.ld",sep=""),head=T)

# eqtl$index<-NA
# eqtl$study_cell<-paste(eqtl$study_label,eqtl$sample_group,sep="_")

# study_cell<-names(table(eqtl$study_cell))
# for (i in 1:length(study_cell)) {
#     eqtl$index[which(eqtl$study_cell==study_cell[i] & eqtl$pvalue==min(eqtl$pvalue[which(eqtl$study_cell==study_cell[i])]))]<-1
# }

# for (i in 1:length(study_cell)) {

#     tmp<-eqtl[which(eqtl$study_cell==study_cell[i]),]

#     # get variants not in high ld:
#     tmp_ld<-ld[which(ld$SNP_A %in% tmp$MarkerName & ld$SNP_B %in% tmp$MarkerName),]

#     tmp_ld<-tmp_ld[which(tmp_ld$SNP_A==tmp$MarkerName[which(tmp$index==1)] | tmp_ld$SNP_B==tmp$MarkerName[which(tmp$index==1)] ),]
#     tmp_ld<-tmp_ld[which(tmp_ld$R2<0.2),]

#     if(nrow(tmp_ld)>0) {

#         tmp_ld<-tmp_ld[which(tmp_ld$R2==min(tmp_ld$R2)),]
#         tmp<-tmp[which(tmp$MarkerName %in% c(as.character(tmp_ld$SNP_A),as.character(tmp_ld$SNP_B))),]

#     } else {
#         tmp<-tmp[which(tmp$index==1),]
#     }
    
#     if(i==1) {
#         eqtl_tmp<-tmp
#     } else {
#         eqtl_tmp<-rbind(eqtl_tmp,tmp)
#     }

# }

# # keep only those in high LD wiht index variant??

eqtl_tmp<-eqtl[which(eqtl$MarkerName %in% coloc$MarkerName),]

eqtl_tmp$BETA<-as.numeric(eqtl_tmp$BETA)
eqtl_tmp$SE<-as.numeric(eqtl_tmp$SE)
eqtl_tmp$beta<-as.numeric(eqtl_tmp$beta)
eqtl_tmp$se<-as.numeric(eqtl_tmp$se)

eqtl_tmp$pvalue_disease<-as.numeric(eqtl_tmp$pvalue_disease)
eqtl_tmp<-eqtl_tmp[which(eqtl_tmp$pvalue_disease<1E-3),]

# eqtl_tmp$study_cell<-paste(e)

# keep a copy of the table:



# plot them by direction of effect
same_dir<-eqtl_tmp[which( (eqtl_tmp$BETA>0 & eqtl_tmp$beta>0) | (eqtl_tmp$BETA<0 & eqtl_tmp$beta<0)),]
same_dir$slope<-1
opp_dir<-eqtl_tmp[which( (eqtl_tmp$BETA<0 & eqtl_tmp$beta>0) | (eqtl_tmp$BETA>0 & eqtl_tmp$beta<0)),]
opp_dir$slope<-(-1)

dir<-rbind(same_dir,opp_dir)
order_study<-names(table(dir$study_label))
order_study<-order_study[!duplicated(order_study)]

dir$study_label<-factor(dir$study_label,levels=order_study)
study_label<-levels(dir$study_label)
same_dir$study_label<-factor(same_dir$study_label,levels=order_study)
opp_dir$study_label<-factor(opp_dir$study_label,levels=order_study)

numbercols_samedir<-length(names(table(as.character(same_dir$study_label))))
study_label_samedir<-names(table(as.character(same_dir$study_label)))

rm(scales_samedir)
for (i in 1:numbercols_samedir) {
    if(i==1) {
        scales_samedir<-list(scale_x_continuous(limits = c(-1*max(abs(same_dir$beta[which(same_dir$study_label==study_label_samedir[i])])+same_dir$se[which(same_dir$study_label==study_label_samedir[i])]),
        max(abs(same_dir$beta[which(same_dir$study_label==study_label_samedir[i])])+same_dir$se[which(same_dir$study_label==study_label_samedir[i])]))))
    } else{
        scales_samedir<-append(scales_samedir,list(scale_x_continuous(limits = c(-1*max(abs(same_dir$beta[which(same_dir$study_label==study_label_samedir[i])])+same_dir$se[which(same_dir$study_label==study_label_samedir[i])]),
        max(abs(same_dir$beta[which(same_dir$study_label==study_label_samedir[i])])+same_dir$se[which(same_dir$study_label==study_label_samedir[i])])))))
    }
}

numbercols_oppdir<-length(names(table(as.character(opp_dir$study_label))))
study_label_oppdir<-names(table(as.character(opp_dir$study_label)))

rm(scales_oppdir)
for (i in 1:numbercols_oppdir) {
    if(i==1) {
        scales_oppdir<-list(scale_x_continuous(limits = c(-1*max(abs(opp_dir$beta[which(opp_dir$study_label==study_label_oppdir[i])])+opp_dir$se[which(opp_dir$study_label==study_label_oppdir[i])]),
        max(abs(opp_dir$beta[which(opp_dir$study_label==study_label_oppdir[i])])+opp_dir$se[which(opp_dir$study_label==study_label_oppdir[i])]))))
    } else{
        scales_oppdir<-append(scales_oppdir,list(scale_x_continuous(limits = c(-1*max(abs(opp_dir$beta[which(opp_dir$study_label==study_label_oppdir[i])])+opp_dir$se[which(opp_dir$study_label==study_label_oppdir[i])]),
        max(abs(opp_dir$beta[which(opp_dir$study_label==study_label_oppdir[i])])+opp_dir$se[which(opp_dir$study_label==study_label_oppdir[i])])))))
    }

}

# for (i in 1:numbercols) {
#     dir$intercept[which(dir$study_label==study_label[i])]<-(-1*max(abs(dir$beta[which(dir$study_label==study_label[i])])))

#     dir$x1[which(dir$study_label==study_label[i])]<-(max(abs(dir$beta[which(dir$study_label==study_label[i])])))
#     dir$y1[which(dir$study_label==study_label[i])]<-(max(abs(dir$BETA[which(dir$study_label==study_label[i])])))

# }
# dir$x<-dir$slope*dir$x
# dir$y<-dir$slope*dir$y

ylimits<-c(-1*max(abs(eqtl_tmp$BETA)+eqtl_tmp$SE),max(abs(eqtl_tmp$BETA)+eqtl_tmp$SE))

# my_colors <- RColorBrewer::brewer.pal(12, "Paired")[1:12]
# col_ligth<-darken(my_colors, 0.80, space = "HCL")
# my_colors<-c(my_colors,col_ligth)
# scales::show_col(my_colors)

if (nrow(opp_dir)>0 & nrow(same_dir)>0) {


    p1<-ggplot(same_dir, aes(y=BETA, x=beta,color=sample_group)) +
        # annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, fill= alpha("lightcyan",.4))  + 
        # annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0 , fill= alpha("lightcyan",.4)) +
        geom_point(size = 3) + ylab(paste("Beta ",coloc$gwas_trait[1]," susceptibility",sep="")) + 
        ylim(ylimits) + 
        # xlim(xlimits) +
        xlab(paste("Beta ",gene_symbol," eQTL",sep=""))  + 
        scale_color_brewer(palette = "Paired") + 
        # scale_colour_manual(values = my_colors) + 
        geom_hline(aes(yintercept=0), color = "darkgrey", linetype = "dashed") +
        geom_vline(aes(xintercept=0), color = "darkgrey", linetype = "dashed") + theme_bw() +
        geom_errorbar(aes(xmin=beta-se, xmax=beta+se), width = 0) + geom_errorbar(aes(ymin=BETA-SE, ymax=BETA+SE), width = 0) +
        facet_grid(. ~ study_label,scales = "free")  +
        theme(plot.title=element_text(face = "italic",size=12),legend.title=element_blank(),
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.text=element_text(size=12),legend.position="top",strip.text.x = element_text(size = 12)) +
        facetted_pos_scales(x = scales_samedir) 
            
    p2<-ggplot(opp_dir, aes(y=BETA, x=beta,color=sample_group)) +
        # annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = -Inf, fill= alpha("lightcoral",.2))  + 
        # annotate("rect", xmin = 0, xmax = -Inf, ymin = Inf, ymax = 0, fill= alpha("lightcoral",.2)) +
        geom_point(size = 3) + ylab(paste("Beta ",coloc$gwas_trait[1]," susceptibility",sep="")) + 
        ylim(ylimits) + 
        # xlim(xlimits) +
        xlab(paste("Beta ",gene_symbol," eQTL",sep=""))  + 
        scale_color_brewer(palette = "Paired") + 
        # scale_colour_manual(values = my_colors) + 
        geom_hline(aes(yintercept=0), color = "darkgrey", linetype = "dashed") +
        geom_vline(aes(xintercept=0), color = "darkgrey", linetype = "dashed") + theme_bw() +
        geom_errorbar(aes(xmin=beta-se, xmax=beta+se), width = 0) + geom_errorbar(aes(ymin=BETA-SE, ymax=BETA+SE), width = 0) +
        facet_grid(. ~ study_label,scales = "free")  +
        theme(plot.title=element_text(face = "italic",size=12),legend.title=element_blank(),
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.text=element_text(size=12),legend.position="top",strip.text.x = element_text(size = 12)) +
        facetted_pos_scales(x = scales_oppdir) 
    p<-ggarrange(p1,p2,widths=c(numbercols_samedir,numbercols_oppdir),ncol=2,align = "h")
    # p

    pdf(paste("~/git/IIBDGC_GWAS/plots/",gene_symbol,"_correlation_beta_susceptibility.pdf",sep=""),height=4, width=sum(numbercols_oppdir,numbercols_samedir)*3)
    print(p)
    dev.off()



    p1<-ggplot(same_dir, aes(y=BETA, x=beta,color=sample_group)) +
        annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, fill= alpha("lightcyan",.6))  + 
        annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0 , fill= alpha("lightcyan",.6)) +
        geom_point(size = 3) + ylab(paste("Beta ",coloc$gwas_trait[1]," susceptibility",sep="")) + 
        ylim(ylimits) + 
        # xlim(xlimits) +
        xlab(paste("Beta ",gene_symbol," eQTL",sep=""))  + 
        scale_color_brewer(palette = "Paired") + 
        # scale_colour_manual(values = my_colors) + 
        geom_hline(aes(yintercept=0), color = "darkgrey", linetype = "dashed") +
        geom_vline(aes(xintercept=0), color = "darkgrey", linetype = "dashed") + theme_bw() +
        geom_errorbar(aes(xmin=beta-se, xmax=beta+se), width = 0) + geom_errorbar(aes(ymin=BETA-SE, ymax=BETA+SE), width = 0) +
        facet_grid(. ~ study_label,scales = "free")  +
        theme(plot.title=element_text(face = "italic",size=12),legend.title=element_blank(),
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.text=element_text(size=12),legend.position="top",strip.text.x = element_text(size = 12)) +
        facetted_pos_scales(x = scales_samedir) 
            
    p2<-ggplot(opp_dir, aes(y=BETA, x=beta,color=sample_group)) +
        annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = -Inf, fill= alpha("lightcoral",.2))  + 
        annotate("rect", xmin = 0, xmax = -Inf, ymin = Inf, ymax = 0, fill= alpha("lightcoral",.2)) +
        geom_point(size = 3) + ylab(paste("Beta ",coloc$gwas_trait[1]," susceptibility",sep="")) + 
        ylim(ylimits) + 
        # xlim(xlimits) +
        xlab(paste("Beta ",gene_symbol," eQTL",sep=""))  + 
        scale_color_brewer(palette = "Paired") + 
        # scale_colour_manual(values = my_colors) + 
        geom_hline(aes(yintercept=0), color = "darkgrey", linetype = "dashed") +
        geom_vline(aes(xintercept=0), color = "darkgrey", linetype = "dashed") + theme_bw() +
        geom_errorbar(aes(xmin=beta-se, xmax=beta+se), width = 0) + geom_errorbar(aes(ymin=BETA-SE, ymax=BETA+SE), width = 0) +
        facet_grid(. ~ study_label,scales = "free")  +
        theme(plot.title=element_text(face = "italic",size=12),legend.title=element_blank(),
        axis.text=element_text(size=12),axis.title=element_text(size=12),
        legend.text=element_text(size=12),legend.position="top",strip.text.x = element_text(size = 12)) +
        facetted_pos_scales(x = scales_oppdir) 
    p<-ggarrange(p1,p2,widths=c(numbercols_samedir,numbercols_oppdir),ncol=2,align = "h")
    # p



    pdf(paste("~/git/IIBDGC_GWAS/plots/",gene_symbol,"_correlation_beta_susceptibility_with_direction_lines.pdf",sep=""),height=4, width=sum(numbercols_oppdir,numbercols_samedir)*4)
    print(p)
    dev.off()



} else if (nrow(opp_dir)==0 & nrow(same_dir)>0) {

        p<-ggplot(same_dir, aes(y=BETA, x=beta,color=sample_group)) +
                # annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, fill= alpha("lightcyan",.4))  + 
                # annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0 , fill= alpha("lightcyan",.4)) +
                geom_point(size = 3) + ylab(paste("Beta ",coloc$gwas_trait[1]," susceptibility",sep="")) + 
                ylim(ylimits) + 
                # xlim(xlimits) +
                xlab(paste("Beta ",gene_symbol," eQTL",sep=""))  + 
                scale_color_brewer(palette = "Paired") + 
                # scale_colour_manual(values = my_colors) + 
                geom_hline(aes(yintercept=0), color = "darkgrey", linetype = "dashed") +
                geom_vline(aes(xintercept=0), color = "darkgrey", linetype = "dashed") + theme_bw() +
                geom_errorbar(aes(xmin=beta-se, xmax=beta+se), width = 0) + geom_errorbar(aes(ymin=BETA-SE, ymax=BETA+SE), width = 0) +
                facet_grid(. ~ study_label,scales = "free")  +
                theme(plot.title=element_text(face = "italic",size=12),legend.title=element_blank(),
                axis.text=element_text(size=12),axis.title=element_text(size=12),
                legend.text=element_text(size=12),legend.position="top",strip.text.x = element_text(size = 12)) +
                facetted_pos_scales(x = scales_samedir) 
            
        pdf(paste("~/git/IIBDGC_GWAS/plots/",gene_symbol,"_correlation_beta_susceptibility.pdf",sep=""),height=4, width=sum(numbercols_oppdir,numbercols_samedir)*4)
            print(p)
        dev.off()

        p<-ggplot(same_dir, aes(y=BETA, x=beta,color=sample_group)) +
                annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, fill= alpha("lightcyan",.4))  + 
                annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0 , fill= alpha("lightcyan",.4)) +
                geom_point(size = 3) + ylab(paste("Beta ",coloc$gwas_trait[1]," susceptibility",sep="")) + 
                ylim(ylimits) + 
                # xlim(xlimits) +
                xlab(paste("Beta ",gene_symbol," eQTL",sep=""))  + 
                scale_color_brewer(palette = "Paired") + 
                # scale_colour_manual(values = my_colors) + 
                geom_hline(aes(yintercept=0), color = "darkgrey", linetype = "dashed") +
                geom_vline(aes(xintercept=0), color = "darkgrey", linetype = "dashed") + theme_bw() +
                geom_errorbar(aes(xmin=beta-se, xmax=beta+se), width = 0) + geom_errorbar(aes(ymin=BETA-SE, ymax=BETA+SE), width = 0) +
                facet_grid(. ~ study_label,scales = "free")  +
                theme(plot.title=element_text(face = "italic",size=12),legend.title=element_blank(),
                axis.text=element_text(size=12),axis.title=element_text(size=12),
                legend.text=element_text(size=12),legend.position="top",strip.text.x = element_text(size = 12)) +
                facetted_pos_scales(x = scales_samedir) 

        pdf(paste("~/git/IIBDGC_GWAS/plots/",gene_symbol,"_correlation_beta_susceptibility_with_direction_lines.pdf",sep=""),height=4, width=sum(numbercols_oppdir,numbercols_samedir)*3)
            print(p)
        dev.off()

} else if (nrow(opp_dir)>0 & nrow(same_dir)==0) {

        p<-ggplot(opp_dir, aes(y=BETA, x=beta,color=sample_group)) +
                # annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = -Inf, fill= alpha("lightcoral",.2))  + 
                # annotate("rect", xmin = 0, xmax = -Inf, ymin = Inf, ymax = 0, fill= alpha("lightcoral",.2)) +
                geom_point(size = 3) + ylab(paste("Beta ",coloc$gwas_trait[1]," susceptibility",sep="")) + 
                ylim(ylimits) + 
                # xlim(xlimits) +
                xlab(paste("Beta ",gene_symbol," eQTL",sep=""))  + 
                scale_color_brewer(palette = "Paired") + 
                # scale_colour_manual(values = my_colors) + 
                geom_hline(aes(yintercept=0), color = "darkgrey", linetype = "dashed") +
                geom_vline(aes(xintercept=0), color = "darkgrey", linetype = "dashed") + theme_bw() +
                geom_errorbar(aes(xmin=beta-se, xmax=beta+se), width = 0) + geom_errorbar(aes(ymin=BETA-SE, ymax=BETA+SE), width = 0) +
                facet_grid(. ~ study_label,scales = "free")  +
                theme(plot.title=element_text(face = "italic",size=12),legend.title=element_blank(),
                axis.text=element_text(size=12),axis.title=element_text(size=12),
                legend.text=element_text(size=12),legend.position="top",strip.text.x = element_text(size = 12)) +
                facetted_pos_scales(x = scales_oppdir)

         pdf(paste("~/git/IIBDGC_GWAS/plots/",gene_symbol,"_correlation_beta_susceptibility.pdf",sep=""),height=4, width=sum(numbercols_oppdir,numbercols_samedir)*3)
            print(p)
        dev.off()
    
        p<-ggplot(opp_dir, aes(y=BETA, x=beta,color=sample_group)) +
            annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = -Inf, fill= alpha("lightcoral",.2))  + 
            annotate("rect", xmin = 0, xmax = -Inf, ymin = Inf, ymax = 0, fill= alpha("lightcoral",.2)) +
            geom_point(size = 3) + ylab(paste("Beta ",coloc$gwas_trait[1]," susceptibility",sep="")) + 
            ylim(ylimits) + 
            # xlim(xlimits) +
            xlab(paste("Beta ",gene_symbol," eQTL",sep=""))  + 
            scale_color_brewer(palette = "Paired") + 
            # scale_colour_manual(values = my_colors) + 
            geom_hline(aes(yintercept=0), color = "darkgrey", linetype = "dashed") +
            geom_vline(aes(xintercept=0), color = "darkgrey", linetype = "dashed") + theme_bw() +
            geom_errorbar(aes(xmin=beta-se, xmax=beta+se), width = 0) + geom_errorbar(aes(ymin=BETA-SE, ymax=BETA+SE), width = 0) +
            facet_grid(. ~ study_label,scales = "free")  +
            theme(plot.title=element_text(face = "italic",size=12),legend.title=element_blank(),
            axis.text=element_text(size=12),axis.title=element_text(size=12),
            legend.text=element_text(size=12),legend.position="top",strip.text.x = element_text(size = 12)) +
            facetted_pos_scales(x = scales_oppdir)

        pdf(paste("~/git/IIBDGC_GWAS/plots/",gene_symbol,"_correlation_beta_susceptibility_with_direction_lines.pdf",sep=""),height=4, width=sum(numbercols_oppdir,numbercols_samedir)*3)
            print(p)
        dev.off()

}


q("no")


# same_dir
#           MarkerName molecular_trait_object_id     beta        se      pvalue
# 6 chr18:63257737:C:T           ENSG00000171791 0.226567 0.0480167 2.99382e-06
#        qvalue gene_symbol study_label  sample_group   BETA     SE
# 6 0.005011946        BCL2        GTEx artery_tibial 0.0329 0.0092
#   pvalue_disease index         study_cell slope
# 6      0.0003685     1 GTEx_artery_tibial     1


# same_dir
#             MarkerName molecular_trait_object_id     beta        se      pvalue
# 227 chr2:111150990:C:T           ENSG00000153094 0.386557 0.0972897 9.13600e-05
# 228 chr2:111150990:C:T           ENSG00000153094 0.386557 0.0972897 9.13600e-05
# 120 chr2:111119036:G:C           ENSG00000153094 0.460179 0.0691987 8.39828e-11
#           qvalue gene_symbol  study_label sample_group   BETA     SE
# 227 2.009654e-02     BCL2L11 Fairfax_2012  B-cell_CD19 0.0978 0.0130
# 228 2.009654e-02     BCL2L11 Fairfax_2012  B-cell_CD19 0.0978 0.0130
# 120 1.982011e-07     BCL2L11   Lepik_2017        blood 0.0963 0.0146
#     pvalue_disease index               study_cell slope
# 227      4.620e-14    NA Fairfax_2012_B-cell_CD19     1
# 228      4.620e-14    NA Fairfax_2012_B-cell_CD19     1
# 120      4.005e-11     1         Lepik_2017_blood     1