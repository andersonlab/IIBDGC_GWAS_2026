# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=4000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q cpu-interactive R


library(data.table)
library(rtracklayer)
library(stringr)
library(qvalue)
library(ggplot2)
library(ggpubr)

rm(list=ls())

# args<-commandArgs(trailingOnly=TRUE)
pheno<-c("ibd","cd","uc")

# for testing purposes

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
# [1] 653 182

# LOAD THE COLOCALIZATION RESULTS - SEE ~/git/IIBDGC_GWAS/scripts/other/evaluate_direction_effect_QTL_risk_variant.R

for (chr in c(1:22)) {
    for (ph in pheno) {
        tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/tmp_coloc_filtered_by_gwas_r2_0.6_allpheno_with_beta_with_wald_ratio_",chr,"_",ph,".tsv"),head=T)
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

for (i in 1:nrow(dat)) {
    dat$ci_low_wald_ratio[i]<-quantile(as.numeric(str_split(dat$betas_wald_ratio[i],"\\|",simplify=T)), c(0.025),na.rm=T)
    dat$ci_up_wald_ratio[i]<-quantile(as.numeric(str_split(dat$betas_wald_ratio[i],"\\|",simplify=T)), c(0.975),na.rm=T)
}

table(dat$cohort[!grepl("ENSG",dat$gene)])
# pQTL_decode  pQTL_sparc 
#          27           6

dat$gene[!grepl("ENSG",dat$gene)]
dat$gene_name[which(dat$gene_name=="")]<-dat$gene[which(dat$gene_name=="")]

table(dat$pheno_coloc)
#   cd  ibd   uc 
# 6233 7099 5490 

coloc_results<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_filtered_by_gwas_r2_0.6_allpheno_noFM_with_beta.tsv",sep=""))
coloc_results<-as.data.frame(coloc_results)

dat$X<-paste(dat$condition_name,dat$gene,dat$MarkerName_eqtl_lead,sep="_")
coloc_results$X<-paste(coloc_results$condition_name,coloc_results$gene,coloc_results$MarkerName_eqtl_lead,sep="_")

dim(dat[which(!dat$X %in% coloc_results$X),])
# 0

dim(coloc_results[which(!coloc_results$X %in% dat$X),])
# 0

length(names(table(dat$gene_name)))
# [1] 680


# for each gene, plot the wald_test results:
summary(dat$n_variants_in_ld_for_wr)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1.00   15.00   42.00   70.85   99.00  516.00

list_genes<-names(table(dat$gene_name))
length(list_genes)
# [1] 680

dat<-merge(dat,map,by.x="condition_name",by.y="id_map")
dim(dat)
# [1] 18822    27

dat$direction_wald<-NA
dat$direction_wald[which(dat$mean_wald_ratio<0)]<-"negative"
dat$direction_wald[which(dat$mean_wald_ratio>0)]<-"positive"

xx<-as.data.frame.matrix(table(dat$gene_name,dat$direction_wald))


# evidence from only one QTL
dim(xx[which( (xx$negative==1 & xx$positive==0) | (xx$negative==0 & xx$positive==1)),])
# [1] 59  2

# nominated by >1 QTL
> 680-59
# [1] 621

dim(xx[which(xx$negative!=0 & xx$positive!=0),]) 
# [1] 86  2

dat<-dat[which(!dat$cell_label %in% c("unannotated","","artery (aorta)","artery (coronary)","artery (tibial)","thyroid","plasma","stomach","esophagus (gej)","esophagus (mucosa)","esophagus (muscularis)",
"hematopoietic precursor cell","sigmoid colon","small intestine","transverse colon","LCL")),]


dat$cell_label[which(dat$cell_label=="B")]<-"B cell"
dat$cell_label[which(dat$cell_label=="memory B cell")]<-"B cell memory"
dat$cell_label[which(dat$cell_label=="gdT cell")]<-"T cell gd"
dat$cell_label[which(dat$cell_label=="blood")]<-"Blood"
dat$cell_label[which(dat$cell_label %in% c("monocyte","CD16+ monocyte"))]<-"Monocyte"
dat$cell_label[which(dat$cell_label %in% c("NK cell","CD56+ NK cell" ))]<-"NK cell"
dat$cell_label[which(dat$cell_label=="neutrophil")]<-"Neutrophil"
dat$cell_label[which(dat$cell_label=="macrophage")]<-"Macrophage"
dat$cell_label[which(dat$cell_label %in% c("T","T cell gd","Tfh cell","Th1 cell","Th17 cell","Th2 cell","MAIT cell",
"CD4+ CTL cell","CD4+ memory T cell","CD4+ T cell","Treg memory","Treg naive",
"CD4+ TCM cell","CD4+ TEM cell","CD8+ T cell","CD8+ TCM cell","CD8+ TEM cell"))]<-"T cell"
dat$cell_label[which(dat$cell_label=="transverse colon")]<-"Transverse colon"
dat$cell_label[which(dat$cell_label=="fibroblast")]<-"Fibroblast"
dat$cell_label[which(dat$cell_label=="platelet")]<-"Platelet"
dat$cell_label[which(dat$cell_label=="ileum")]<-"Ileum"
dat$cell_label[which(dat$cell_label=="rectum")]<-"Rectum"
dat$cell_label[which(dat$cell_label=="plasmablast")]<-"Plasmablast"

names(table(dat$cell_label))

dat$cell_label<-factor(dat$cell_label,levels=c("Blood","LCL","B cell","B cell memory",
"T cell","Macrophage","Monocyte","Myeloid","Neutrophil","NK cell","Platelet","Plasmablast","Plasma",
"Fibroblast","Transverse colon" ,"Ileum","Rectum","Colonocyte","Enterocyte","Mesenchymal","Secretory","Stem"))

names(table(dat$cell_label))

# Load pathways:
pathways<-fread(paste0("~/git/IIBDGC_GWAS/plots/pathways/tables/gene_enrichment_analysis_pathways_go_all_genes.txt"))


# Not selected:
# ,"INTERLEUKIN-23 RECEPTOR COMPLEX" - not enought coloc
# "ENDOTHELINS" - not informatice
# "INTEGRIN BINDING" - too large
# ,"INTEGRIN-MEDIATED SIGNALING PATHWAY"
# "LEUKOCYTE MIGRATION","RESPONSE TO TYPE II INTERFERON",
# "POSITIVE REGULATION OF LEUKOCYTE CELL-CELL ADHESION"
# "NEGATIVE REGULATION OF LYMPHOCYTE ACTIVATION"
# "AUTOPHAGOSOME ORGANIZATION"
# "RAF MAP KINASE CASCADE")
# "B CELL MEDIATED IMMUNITY"
# "NEGATIVE REGULATION OF CELL MIGRATION"
# "B CELL ACTIVATION"
# T CELL ACTIVATION
# CELLULAR EXTRAVASATION
# "IFN-GAMMA"
# "RESPONSE TO GROWTH FACTOR"
# "CELL ADHESION MEDIATED BY INTEGRIN"
# "MYELOID LEUKOCYTE MIGRATION",

# "INTERFERON SIGNALING",
# # 
# "NEGATIVE REGULATION OF PHOSPHATIDYLINOSITOL BIOSYNTHETIC PROCESS",
# ",,,,
# # ,
# "REACTIVE OXYGEN SPECIES METABOLIC PROCESS","REGULATION OF STRESS-ACTIVATED MAPK CASCADE","PHOSPHOLIPID METABOLIC PROCESS")



pathways_to_plot<-c("GRANULOCYTE MIGRATION","T CELL MIGRATION","INTEGRIN CELL SURFACE INTERACTIONS",
"UROKINASE PLASMINOGEN ACTIVATOR SIGNALING PATHWAY",
"PLATELET-DERIVED GROWTH FACTOR RECEPTOR SIGNALING PATHWAY","ENDOTHELINS","PID_ENDOTHELIN_PATHWAY","UNTITLED",
"INTERFERON SIGNALING","PHOSPHOLIPID METABOLIC PROCESS","PHOSPHATIDYLSERINE BINDING")

# pathways_to_plot<-c("REACTIVE OXYGEN SPECIES METABOLIC PROCESS")
xlim<-c(-40,40)

i=10

for (i in 1:length(pathways_to_plot)) {

    gene_plot<-unlist(strsplit(pathways$ibd[which(pathways$path==pathways_to_plot[i])],";"))

    tmp<-dat[which(dat$gene_name %in% gene_plot),]
        tmp$direction_wald<-factor(tmp$direction_wald,levels=c("negative","positive"))

    if (i %in% c(6,7)) {
        xlim<-c(round(max(abs(tmp$ci_low_wald_ratio),abs(tmp$ci_low_wald_ratio))*-1)-10,round(max(abs(tmp$ci_low_wald_ratio),abs(tmp$ci_low_wald_ratio)))+10)
    } else {
        xlim<-c(-40,40)
    }


    # sort order for some figures

    # if (pathways_to_plot[i]==c("INTERFERON SIGNALING")) {
    #     tmp$IIBDGC_GWAS_index_variant<-factor(tmp$IIBDGC_GWAS_index_variant,levels=c("chr5:132435113:C:T","chr5:132461855:A:G","chr7:128941781:A:G","chr16:85980635:T:C","chr16:85985910:G:C","chr21:33397315:TAA:T","chr17:42375526:A:G","chr9:4984530:G:C","chr18:12853459:T:G","chr10:102465559:C:G","chr10:73853123:T:C","chr10:73934322:G:GCT","chr11:308059:T:C","chr17:56822584:G:A","chr17:56830770:A:T","chr2:230249416:G:A", "chr3:122386581:A:G"))
    # } else if (pathways_to_plot[i]==c("PHOSPHOLIPID METABOLIC PROCESS")) {
    #     tmp$IIBDGC_GWAS_index_variant<-factor(tmp$IIBDGC_GWAS_index_variant,levels=c("chr18:49639885:T:C","chr19:40712706:G:A",
    #     "chr3:16883833:A:G","chr6:122400365:A:G","chr11:63929911:G:A","chr11:63796342:TC:T","chr16:28517460:G:A","chr16:28530851:T:C","chr11:61820833:A:G","chr17:56822584:G:A","chr17:56830770:A:T",
    #     "chr1:212839129:T:TTA","chr1:56456450:C:A","chr1:56463853:C:A","chr10:50592671:C:T","chr11:64248783:C:A",
    #     "chr2:233261680:A:T","chr20:6084192:T:C"))
    # } else if (pathways_to_plot[i]==c("REACTIVE OXYGEN SPECIES METABOLIC PROCESS")) {
    #     tmp$IIBDGC_GWAS_index_variant<-factor(tmp$IIBDGC_GWAS_index_variant,levels=c("chr1:200906769:T:C","chr16:69627624:TCAGGCTGCATAGTG:T","chr11:64248783:C:A","chr1:8020747:A:G","chr22:37470025:G:A","chr22:39264824:T:C","chr3:49684099:G:A","chr6:109055444:C:G","chr7:75073650:C:T"))
    # }

     


    # p<-
    ggplot(tmp, aes(x=mean_wald_ratio,y=gene_name,color=direction_wald)) +
            geom_point() + theme(axis.title.y=element_blank(),) + 
            xlim(xlim) + theme_bw() +
            theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
            legend.position="none",
            strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0),hjust = 0),
            strip.text.x = element_text(margin = margin(0,0,0,0),vjust = 0)) + 
            geom_vline(xintercept=c(0), linetype="dotted")+
            scale_color_manual(values=c("blue","red"),drop=F) +
            guides(color="none", size=guide_legend()) + facet_grid(IIBDGC_GWAS_index_variant ~ cell_label,scales = "free",space = "free") +
            labs(x="Mean Beta-wald (95% quantile)") + 
            ggtitle(pathways_to_plot[i])


      assign(paste0("p",i),p)

        # # save tmp plots for inicial visualization 
        # ggsave(
        # paste0("~/git/IIBDGC_GWAS/plots/tmp_figures/p",i,"_",pathways_to_plot[i],".png"),
        # p,
        # width = 10+(length(table(as.character(tmp$cell_label)))*15),
        # height = 5+(length(table(tmp$gene_name))*5),
        # dpi = 300,
        # units = c("mm"),
        # limitsize = T,scale=2
        # )

}





## plot an extended version for CARD9
list_genes<-"CARD9"

for (gene in list_genes) {

    # print(gene)

    tmp<-dat[which(dat$gene_name==gene),]
    tmp$direction_wald<-factor(tmp$direction_wald,levels=c("negative","positive"))

    variants<-names(table(tmp$IIBDGC_GWAS_index_variant))

    for (i in c(1:length(variants))) {

        var<-variants[i]

        tmp<-dat[which(dat$gene_name==gene & dat$IIBDGC_GWAS_index_variant==var),]
        tmp$direction_wald<-factor(tmp$direction_wald,levels=c("negative","positive"))

        tmp$tissue_cell_condition_label<-paste(tmp$study_label.x,tmp$tissue_cell_condition_label)

        # xlim<-c(round(max(abs(tmp$ci_low_wald_ratio),abs(tmp$ci_low_wald_ratio))*-1)-10,round(max(abs(tmp$ci_low_wald_ratio),abs(tmp$ci_low_wald_ratio)))+10)

        # keep only macrophages data:

        tmp<-tmp[which(tmp$cell_label=="Macrophage"),]

        p9<-ggplot(tmp, aes(x=mean_wald_ratio,y=tissue_cell_condition_label,color=direction_wald)) +
        geom_point() +  xlim(xlim) +
        geom_errorbar(aes(xmin=ci_low_wald_ratio, xmax=ci_up_wald_ratio), width=0.2) + theme_bw() +
        scale_size_continuous(limits=c(0,520), breaks=seq(0,520, by=40),name="Number variants (r2≥0.6)") + 
        scale_color_manual(values=c("blue","red"),drop=F) +
        guides(color="none", size=guide_legend()) + ggtitle(paste(gene, var)) + facet_grid(cell_label~ pheno_coloc,scales = "free",, space = "free") +
        labs(x="Mean Beta-wald (95% quantile)") + 
        theme(axis.title.y=element_blank(),strip.text.y.right = element_text(angle = 0),plot.margin = unit(c(0,0.5,0,1), "cm"))+ geom_vline(xintercept=c(0), linetype="dotted")

    }
}


p4r<-ggarrange(p4,p5,labels=c("","e"),nrow=2,heights=c(2.5,3))

p1r<-ggarrange(p1,p2,p3,p4r,labels=c("a","b","c","d"),nrow=2,ncol=2)
p2r<-ggarrange(p6,p7,labels=c("f","g"),nrow=2,heights=c(2,2.5))

p8r<-ggarrange(p8,p9,ncol=2,labels=c("h","i"))

p<-ggarrange(p1r,p8r,p2r,nrow=3,heights=c(2,0.9,4.2))


ggsave(
        paste0("~/git/IIBDGC_GWAS/plots/tmp_figures/Figure4.png"),
        p,
        width = 280,
        height = 320,
        dpi = 300,
        units = c("mm"),
        limitsize = T,scale=2
)


# Core Signaling Pathway: NOD2-INAVA-CARD9-MAPK 
# # NOD2 Activation: NOD2 recognizes muramyl dipeptide (MDP) and, upon activation, forms a complex with the adaptor protein CARD9. While traditionally associated with NF-κB, this interaction drives crucial p38 and JNK MAPK activation.
# # INAVA Amplifier: INAVA (Inflammatory Bowel Disease-Associated Protein) acts as a scaffold that recruits 14-3-3τ to the NOD2-RIP2-CARD9 complex. This recruitment is essential for optimal activation of p-ERK, p-p38, and p-IκBα.
# # CARD9 Dependency: CARD9 is specifically required for NOD2-mediated activation of the MAPK pathway (p38/JNK), particularly in the context of bacterial clearance (e.g., AIEC, S. aureus).
# # Inflammatory Output: This complex drives pro-inflammatory cytokine secretion (e.g., IL-6) and enhances ROS, RNS, and autophagy pathways to clear bacteria. 
# # Regulatory Interplay with FOXO1 and MYC 
# # The MAPK cascade initiated by NOD2/INAVA/CARD9 intersects with metabolic and proliferation checkpoints:
# # FOXO1-MYC Antagonism: FOXO1 typically promotes quiescence by antagonizing MYC, reducing metabolic activity and proliferation in cells.
# # Signaling Context: In scenarios of high stress or infection, where NOD2-INAVA is highly active, this likely occurs in a environment where FOXO1's role in inhibiting MYC is diminished (i.e., high Akt phosphorylation, which degrades FOXO1, or low nutrient stress).
# # IBD Connection: Lower INAVA expression (observed in IBD rs7554511 risk carriers) leads to decreased MAPK/NF-κB signaling and impaired bacterial clearance. This, combined with potential dysregulation of the FOXO1/MYC metabolic axis, likely contributes to the chronic inflammation seen in Crohn's disease. 


# # NOD2 - activates
# # INAVA - amplifies

# # FOXO1-MYC