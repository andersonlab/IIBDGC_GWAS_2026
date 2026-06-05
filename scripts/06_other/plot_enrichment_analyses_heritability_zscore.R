# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=2000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R 

library(data.table)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(rcartocolor)
library(qvalue)

path_gwas="/path/to/ibdgwas/IIBDGC/"

586+103
# [1] 689 - OK, 103 from baseline

threshold<-0.05/((1+54+45+55+180+251+103)*3)
threshold
# [1] 2.418965e-05

col<-carto_pal(12, "Safe")
col_ligth_1<-lighten(col, 0.90, space = "HCL")
col_ligth_2<-lighten(col, 0.70, space = "HCL")
col_ligth_3<-lighten(col, 0.50, space = "HCL")
col_ligth_4<-lighten(col, 0.30, space = "HCL")

colo<-c(col[0:11],col_ligth_1[0:11],col_ligth_2[0:11],col_ligth_3[0:11],col_ligth_4[0:11])


files<-list.files(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",sep=""))
files<-files[grep(".results",files)]
length(files)
# [1] 1770

files_2<-list.files(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/other_traits/",sep=""))
files_2<-files_2[grep(".results",files_2)]
length(files_2)
# [1] 585

files<-c(files,files_2)


###################
# baseline model:

files_baseline<-files[grep(".baseline.files_version_april2025.results",files)]

for (i in 1:length(files_baseline)) {

    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_baseline[i],sep=""),head=T)
    tmp$pheno<-toupper(gsub(".baseline.files_version_april2025.results","",files_baseline[i]))
    tmp$Category<-gsub("L2_0","",tmp$Category)

    if(i==1) {
        enrich<-tmp
    } else {
        enrich<-rbind(enrich,tmp)
    }

}

dim(enrich)
# [1] 309  11
309/3
# [1] 103

# only colour the significant results:
enrich$signif<-NA
enrich$signif[which(enrich$Enrichment_p<threshold)]<-"Signif"
enrich$zscore<-enrich$'Coefficient_z-score'

p1<-ggplot(enrich[!grepl("maf_bin",enrich$Category)], aes(y=Category, x=zscore,fill=signif)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) +
  # geom_text(vjust = 0,hjust=0.9,size=3) +
   scale_fill_manual(limits = levels(enrich$signif),values=c("lightblue")) +
#   geom_errorbar(aes(xmin=zscore-Coefficient_std_error, xmax=zscore+Coefficient_std_error), width=.2,
#                  position=position_dodge(.9)) + 
                 theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
                 plot.margin = unit(c(1,1,1,1), "cm")) + facet_grid(~ pheno, scales = "free",space = "free")+ theme(legend.position="none")


ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_enrichment_analyses_baseline_model_zscore.png",
  p1,
  width = 120,
  height = ceiling(length(levels(as.factor(enrich$Category)))*1.8),
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)


enrich$source<-"baseline_model"
enrich_baseline<-enrich
rm(enrich)

# ##########
# immune cell atlas:

files_atac<-files[grep("ATAC_counts",files)]
length(files_atac)
# 180

length(files_atac)/4
# [1] 45

data<-gsub(".files_version_april2025.results","",files_atac)

for (i in 1:length(files_atac)) {

    if(file.exists(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_atac[i],sep=""))) {
      tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_atac[i],sep=""),head=T)
    } else {
      tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/other_traits/",files_atac[i],sep=""),head=T)
    }

    tmp<-tmp[which(tmp$Category=="L2_1",)]
    tmp$Category<-data[i]
    tmp$pheno<-toupper(gsub(".baseline.*","",tmp$Category))
    tmp$Category<-(gsub(".*baseline.ATAC_counts_","",tmp$Category))
    tmp$Category<-(gsub("_merged_samples","",tmp$Category))

    if(i==1) {
        enrich<-tmp
    } else {
        enrich<-rbind(enrich,tmp)
    }
}

enrich$pheno<-factor(enrich$pheno,levels=c("CD","IBD","UC","HEIGHT"))
enrich$zscore<-enrich$'Coefficient_z-score'

enrich$cell_type<-gsub("_merged_samples","",enrich$Category)
enrich$cell_type<-gsub("-U$","",enrich$cell_type)
enrich$cell_type<-gsub("-S$","",enrich$cell_type)

enrich$cell_type<-factor(enrich$cell_type,levels = (c("pDCs","Myeloid_DCs","Monocytes",
                                                  "Immature_NK","Mature_NK","Memory_NK",
                                                  "Gamma_delta_T",
                                                  "Effector_CD4pos_T",
                                                  "Naive_Teffs","Memory_Teffs","Th1_precursors","Th2_precursors","Th17_precursors","Follicular_T_Helper",
                                                  "Naive_Tregs","Memory_Tregs",
                                                  "CD8pos_T","Naive_CD8_T","Central_memory_CD8pos_T","Effector_memory_CD8pos_T",
                                                  "Regulatory_T",
                                                  "Bulk_B","Naive_B","Mem_B","Plasmablasts")))



enrich$cell<-"T-cell"
enrich$cell[enrich$cell_type %in% c("Mature_NK","Mature_NK","Memory_NK","Immature_NK")]<-"NK"
enrich$cell[enrich$cell_type %in% c("Mem_B","Plasmablasts","Naive_B","Bulk_B")]<-"B-cell"
enrich$cell[enrich$cell_type %in% c("Monocytes")]<-"Monocyte"
enrich$cell[enrich$cell_type %in% c("Myeloid_DCs","pDCs")]<-"Dendritic cells"

enrich$cell<-factor(enrich$cell,levels=c("B-cell","T-cell","Monocyte","NK","Dendritic cells"))

enrich$class<-"ATAC-seq"
enrich_atac<-enrich
dim(enrich_atac)
# [1] 180 14


# order by cell
enrich<-enrich[order(enrich$cell,enrich$cell_type),]
levels_sort<-enrich$Category[!duplicated(enrich$Category)]
enrich$Category_ed<-enrich$Category
enrich$Category_ed<-factor(enrich$Category_ed,levels=levels_sort)

# only colour the significant results:
enrich$signif<-enrich$cell
enrich$signif[which(enrich$Enrichment_p>threshold)]<-NA

# p3<-ggplot(enrich, aes(x=Category, y=zscore,fill=cell,label = ifelse(Enrichment_p < threshold, "*", ""))) + 
p3<-ggplot(enrich, aes(y=Category, x=zscore,fill=signif)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) + 
  # geom_text(vjust = 0,hjust=0.9,size=3) +
  scale_fill_manual(limits = levels(enrich$signif),values=c(colo)) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
                 plot.margin = unit(c(1,1,1,1), "cm")) + facet_grid(~ pheno, space = "free") + theme(legend.position="right",legend.title=element_blank())

ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_enrichment_analyses_immune_cell_atlas_zscore.png",
  p3,
  width = 120,
  height = ceiling(length(levels(as.factor(enrich$Category)))*2),
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)


enrich$source<-"immune_cell_atlas"
enrich_immune<-enrich
rm(enrich)

#######################
# zhang

files_cre<-files[grep("cCRE",files)]
length(files_cre)
# [1] 216

length(files_cre)/4
# 54


data<-gsub(".files_version_april2025.results","",files_cre)


for (i in 1:length(files_cre)) {

    if(file.exists(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_cre[i],sep=""))) {
      tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_cre[i],sep=""),head=T)
    } else {
      tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/other_traits/",files_cre[i],sep=""),head=T)
    }

    tmp<-tmp[which(tmp$Category=="L2_1",)]
    tmp$Category<-data[i]
    tmp$pheno<-toupper(gsub(".baseline.*","",tmp$Category))
    tmp$Category<-(gsub(".*baseline.cCRE_Accessibility_","",tmp$Category))

    if(i==1) {
        enrich<-tmp
    } else {
        enrich<-rbind(enrich,tmp)
    }
}

enrich$pheno<-factor(enrich$pheno,levels=c("CD","IBD","UC","HEIGHT"))
enrich$zscore<-enrich$'Coefficient_z-score'

enrich$cell_type<-enrich$Category

enrich$cell<-enrich$cell_type
enrich$cell<-gsub("\\.[0-9]{1}$","",enrich$cell)
enrich$cell<-gsub("[0-9]","",enrich$cell)
enrich$cell[which(enrich$cell=="Adi")]<-"Adipocyte"
enrich$cell[which(enrich$cell=="Mes")]<-"Mesothelial cell"
enrich$cell[which(enrich$cell=="Grn")]<-"Granulosa cell"
enrich$cell[which(enrich$cell=="Smm")]<-"Smooth muscle"
enrich$cell[which(enrich$cell=="Vsm")]<-"Vascular smooth muscle"
enrich$cell[which(enrich$cell=="Esm")]<-"Esophageal smooth muscle"
enrich$cell[which(enrich$cell=="Mfb")]<-"Myofibroblast"
enrich$cell[which(enrich$cell=="Swn")]<-"Schwann cell"
enrich$cell[which(enrich$cell=="Fib")]<-"Fibroblast"
enrich$cell[which(enrich$cell=="Stl")]<-"Myosatellite cell"
enrich$cell[which(enrich$cell=="Skm")]<-"Skeletal myocyte"
enrich$cell[which(enrich$cell=="Cam")]<-"Cardiomyocyte"
enrich$cell[which(enrich$cell=="Mac")]<-"Macrophage"
enrich$cell[which(enrich$cell=="Tly")]<-"T-cell"
enrich$cell[which(enrich$cell=="Bly")]<-"B-cell"
enrich$cell[which(enrich$cell=="Mst")]<-"Mast cell"
enrich$cell[which(enrich$cell=="End")]<-"Endothelial cell"
enrich$cell[which(enrich$cell=="Msc")]<-"Miscellaneous stromal cell"
enrich$cell[which(enrich$cell=="Pal")]<-"Pheumocyte"
enrich$cell[which(enrich$cell=="Krt")]<-"Keratinocyte"
enrich$cell[which(enrich$cell=="Eep")]<-"Esophageal epithelial cell"
enrich$cell[which(enrich$cell=="Enc")]<-"Enterocyte"
enrich$cell[which(enrich$cell=="Gbl")]<-"Goblet cell"
enrich$cell[which(enrich$cell=="Fol")]<-"Follicular cell of thyroid"
enrich$cell[which(enrich$cell=="Lue")]<-"Luminal epithelial cell"
enrich$cell[which(enrich$cell=="Epc")]<-"Epithelial cells"
enrich$cell[which(enrich$cell=="Bas")]<-"Basal cell"
enrich$cell[which(enrich$cell=="Agb")]<-"Airway goblet cell"
enrich$cell[which(enrich$cell=="Dut")]<-"Ductal cell"
enrich$cell[which(enrich$cell=="Acn")]<-"Acinar cell of pancreas"
enrich$cell[which(enrich$cell=="Prt")]<-"Gastric parietal cell"
enrich$cell[which(enrich$cell=="Gcf")]<-"Gastric chief cell"
enrich$cell[which(enrich$cell=="Adc")]<-"Cortical cell of adrenal gland"
enrich$cell[which(enrich$cell=="Nec")]<-"Neuroendocrine cell"
enrich$cell[which(enrich$cell=="Hpc")]<-"Hepatocyte"

enrich$cell<-as.factor(enrich$cell)
# levels(enrich$cell)


enrich$cell<-factor(enrich$cell,levels=c("B-cell","T-cell","Macrophage",
"Endothelial cell","Enterocyte","Epithelial cells","Goblet cell",
"Acinar cell of pancreas","Adipocyte","Airway goblet cell","Basal cell","Cardiomyocyte",
"Cortical cell of adrenal gland","Ductal cell",
"Esophageal epithelial cell","Esophageal smooth muscle","Fibroblast","Follicular cell of thyroid",
"Gastric chief cell","Gastric parietal cell","Granulosa cell","Hepatocyte","Keratinocyte",
"Luminal epithelial cell","Mast cell","Mesothelial cell","Miscellaneous stromal cell",
"Myofibroblast","Myosatellite cell","Neuroendocrine cell","Pheumocyte","Schwann cell","Skeletal myocyte",
"Smooth muscle","Vascular smooth muscle"))


enrich$class<-paste(enrich$cell," (",enrich$Category,")",sep="")

# order by cell
enrich<-enrich[order(enrich$cell),]
levels_sort<-enrich$class[!duplicated(enrich$class)]
enrich$class_ed<-enrich$class
enrich$class_ed<-factor(enrich$class_ed,levels=levels_sort)

# only colour the significant results:
enrich$signif<-enrich$cell
enrich$signif[which(enrich$Enrichment_p>threshold)]<-NA


p4<-ggplot(enrich, aes(y=class_ed, x=zscore,fill=signif)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) + 
  # geom_text(vjust = 0,hjust=0.9,size=3) +
  scale_fill_manual(limits = levels(enrich$signif),values=c(colo)) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
                 plot.margin = unit(c(1,1,1,1), "cm")) + facet_grid(~ pheno, space = "free") + 
                 theme(legend.position="bottom",legend.title=element_blank(),legend.justification = c("right", "bottom"),)


ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_enrichment_analyses_ciselement_atlas_zscore.png",
  p4,
  width = 140,
  height = ceiling(length(levels(as.factor(enrich$Category)))*3),
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)




enrich$source<-"cis_element_atlas"
enrich_cis_element_atlas<-enrich
rm(enrich)


#######################
# hickey

files_hic<-files[grep("hickey|colon_epithelial",files)]
files_hic<-files_hic[!grepl("peak|both_annotations",files_hic)]
length(files_hic)
# [1] 719

data<-gsub(".files_version_april2025.results","",files_hic)


for (i in 1:length(files_hic)) {

    if(file.exists(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_hic[i],sep=""))) {
      tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_hic[i],sep=""),head=T)
    } else {
      tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/other_traits/",files_hic[i],sep=""),head=T)
    }
    tmp<-tmp[which(tmp$Category=="L2_1",)]
    tmp$Category<-data[i]
    tmp$pheno<-toupper(gsub(".baseline.*","",tmp$Category))
    tmp$Category<-(gsub(".*baseline.","",tmp$Category))
    tmp$Category<-(gsub("_peaks","",tmp$Category))

    if(i==1) {
        enrich<-tmp
    } else {
        enrich<-rbind(enrich,tmp)
    }
}

enrich$pheno<-factor(enrich$pheno,levels=c("CD","IBD","UC","HEIGHT"))
enrich<-enrich[!grepl("Unknown",enrich$Category),]
enrich$zscore<-enrich$'Coefficient_z-score'

enrich$cell_type<-enrich$Category

enrich$tissue<-NA
enrich$tissue[grep("colon_epithelial",enrich$cell_type)]<-"colon"
enrich$tissue[grep("ileum_epithelial",enrich$cell_type)]<-"ileum_epithelial"
enrich$tissue[grep("duodenum_epithelial",enrich$cell_type)]<-"duodenum_epithelial"
enrich$tissue[grep("immune",enrich$cell_type)]<-"immune"
enrich$tissue[grep("jejunum_epithelial",enrich$cell_type)]<-"jejunum_epithelial"
enrich$tissue[grep("stromal",enrich$cell_type)]<-"stromal"
enrich$tissue[grep("colon",enrich$cell_type)]<-"colon"
enrich$tissue<-factor(enrich$tissue,levels=c("immune","stromal","duodenum_epithelial","jejunum_epithelial","ileum_epithelial","colon"))

table(enrich$tissue,useNA="ifany")
  #            immune             stromal duodenum_epithelial  jejunum_epithelial 
  #                72                 116                 148                 132 
  #  ileum_epithelial               colon 
  #               132                 116 


enrich$cell_type<-gsub("non_multiome_","",enrich$cell_type)
enrich$cell_type<-gsub("colon_epithelial_","",enrich$cell_type)
enrich$cell_type<-gsub("ileum_epithelial_","",enrich$cell_type)
enrich$cell_type<-gsub("duodenum_epithelial_","",enrich$cell_type)
enrich$cell_type<-gsub("immune_","",enrich$cell_type)
enrich$cell_type<-gsub("jejunum_epithelial_","",enrich$cell_type)
enrich$cell_type<-gsub("stromal_","",enrich$cell_type)
enrich$cell_type<-gsub("_hickey","",enrich$cell_type)
table(enrich$cell_type)
            #         Adipocytes                        B_Cells 
            #                  4                              8 
            #  Best4_Enterocytes                            CD4 
            #                 32                              8 
            #                CD8      Crypt_Fibroblasts_1_WNT2B 
            #                  8                              8


enrich$assay<-"multiome"
enrich$assay[grep("non_multiome",enrich$Category)]<-"non_multiome"

table(enrich$assay)

enrich$cell<-enrich$cell_type
enrich$cell[which(enrich$cell_type %in% c("Best4_Enterocytes","Enterocytes","Immature_Enterocytes"))]<-"Enterocyte"
enrich$cell[which(enrich$cell_type %in% c("Adipocytes"))]<-"Adipocytes"
enrich$cell[which(enrich$cell_type %in% c("B_Cells"))]<-"B-cell"
enrich$cell[which(enrich$cell_type %in% c("NK","NK_Cells"))]<-"NK"
enrich$cell[which(enrich$cell_type %in% c("T_Cells","CD8","CD4"))]<-"T-cell"
enrich$cell[which(enrich$cell_type %in% c("Myofibroblasts_SM_1","Myofibroblasts_SM_2","Myofibroblasts_SM_3","Myofibroblasts_SM_DES_High"))]<-"Myofibroblast"
enrich$cell[which(enrich$cell_type %in% c("Secretory_Specialized_MUC5B","Secretory_Specialized_MUC6"))]<-"Secretory cell"
enrich$cell[which(enrich$cell_type %in% c("Crypt_Fibroblasts_1_WNT2B","Crypt_Fibroblasts_2","Crypt_Fibroblasts_3_RSPO3","Villus_Fibroblasts_WNT5B"))]<-"Fibroblast"
enrich$cell[which(enrich$cell_type %in% c("Epithelial"))]<-"Epithelial cell"
enrich$cell[which(enrich$cell_type %in% c("Enteroendocrine","EnteroendocrineUn","EnteroendocrineUn_1"))]<-"Enteroendocrine" 
enrich$cell[which(enrich$cell_type %in% c("Endothelial","Endothelial_CD36_Microvascular","Endothelial_Venules","Lymphatic_endothelial_cells",
"Lymphatic_Endothelial_Cells"))]<-"Endothelial cell"  
enrich$cell[which(enrich$cell_type %in% c("Immature_Goblet","Goblet"))]<-"Goblet cell" 
enrich$cell[which(enrich$cell_type %in% c("NEUROG3high","Neurons"))]<-"Neurons"
enrich$cell[which(enrich$cell_type %in% c("Mono_Macrophages"))]<-"Monocyte/Macrophage" 
enrich$cell[which(enrich$cell_type %in% c("D_Cells"))]<-"D-cell"
enrich$cell[which(enrich$cell_type %in% c("DC"))]<-"Dendritic cells" 
enrich$cell[which(enrich$cell_type %in% c("S_Cells"))]<-"S-cell" 

table(enrich$cell,useNA="ifany")
  #        Adipocytes              B-cell       CyclingImmune         CyclingTA_1 
  #                 4                   8                   4                  16 
  #       CyclingTA_2              D-cell     Dendritic cells    Endothelial cell 
  #                 4                  16                   4                  20 
  #  Enterochromaffin          Enterocyte     Enteroendocrine     Epithelial cell 
  #                16                  96                  48                  12 
  #          Exocrine          Fibroblast                Glia         Goblet cell 
  #                 4                  32                   8                  64 
  #           I_Cells                 ICC                 ILC             K_Cells 
  #                16                   8                   4                  12 
  #           L_Cells                Mast            Mo_Cells Monocyte/Macrophage 
  #                16                   8                  16                   8 
  #     Myofibroblast             Neurons                  NK              Paneth 
  #                28                  24                   8                  24 
  #         Pericytes              Plasma              S-cell      Secretory cell 
  #                 8                   8                  16                  12 
  #              Stem              T-cell                 TA1                 TA2 
  #                32                  20                  32                  32 
  #              Tuft 
  #                28

enrich$cell<-factor(enrich$cell,levels=c("B-cell","T-cell","Monocyte/Macrophage","Dendritic cells","NK",
"Endothelial cell","Enterocyte","Epithelial cell","Goblet cell","Adipocytes","CyclingImmune","CyclingTA_1","CyclingTA_2","D-cell",
"Enterochromaffin","Enteroendocrine","Exocrine","Fibroblast",
"Glia","I_Cells","ICC","ILC","K_Cells","L_Cells",
"Mast","Mo_Cells","Myofibroblast","Neurons","Paneth","Pericytes","Plasma","S-cell","Secretory cell","Stem",
"TA1","TA2","Tuft"))
table(enrich$cell,useNA="ifany")


# order by cell
enrich<-enrich[order(enrich$cell),]
levels_sort<-enrich$cell[!duplicated(enrich$cell)]

enrich$cell<-factor(enrich$cell,levels=levels_sort)

# only colour the significant results:
enrich$signif<-enrich$cell
enrich$signif[which(enrich$Enrichment_p>threshold)]<-NA


multiome<-enrich[which(enrich$assay=="multiome"),]
nomultiome<-enrich[which(enrich$assay=="non_multiome"),]

p5.1<-ggplot(multiome, aes(y=cell_type, x=zscore,fill=signif)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) + 
  scale_fill_manual(limits = levels(enrich$signif),values=c(colo)) +
                 theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
                 plot.margin = unit(c(1,1,1,1), "cm")) + facet_grid(tissue ~ pheno, scales = "free_y", space = "free") + 
                 theme(legend.position="bottom",legend.title=element_blank()) + xlim(0,max(enrich$zscore))

p5.2<-ggplot(nomultiome, aes(y=cell_type, x=zscore,fill=signif)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) + 
  scale_fill_manual(limits = levels(enrich$signif),values=c(colo)) +
theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
                 plot.margin = unit(c(1,1,1,1), "cm")) + facet_grid(tissue ~ pheno, scales = "free_y", space = "free") + 
                 theme(legend.position="bottom",legend.title=element_blank()) + xlim(0,max(enrich$zscore))

p5<-ggarrange(p5.1,p5.2,ncol=2,common.legend=T,legend="bottom",labels = c("a", "b"),align = "v",widths=c(1.1,1))

ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_enrichment_analyses_hickey_zscore.png",
  p5,
  width = 180,
  height = ceiling(length(levels(as.factor(enrich$cell_type)))*3.5),
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2.2
)



enrich$source<-"hickey"
enrich_hickey<-enrich
rm(enrich)


#######################
# soskic

files_sos<-files[grep("ChM",files)]
length(files_sos)
# [1] 220

data<-gsub(".files_version_april2025.results","",files_sos)

for (i in 1:length(files_sos)) {

     if(file.exists(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_sos[i],sep=""))) {
      tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_sos[i],sep=""),head=T)
    } else {
      tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/other_traits/",files_sos[i],sep=""),head=T)
    }
    tmp<-tmp[which(tmp$Category=="L2_1",)]
    tmp$Category<-data[i]
    tmp$pheno<-toupper(gsub(".baseline.*","",tmp$Category))
    tmp$Category<-(gsub(".*baseline.","",tmp$Category))
    tmp$Category<-(gsub("_peaks","",tmp$Category))

    if(i==1) {
        enrich<-tmp
    } else {
        enrich<-rbind(enrich,tmp)
    }
}


enrich$pheno<-factor(enrich$pheno,levels=c("CD","IBD","UC","HEIGHT"))
enrich$zscore<-enrich$'Coefficient_z-score'

enrich$cell_type<-gsub("ChM_","",enrich$Category)
enrich$cell_type<-gsub("_.*","",enrich$cell_type)
table(enrich$cell_type)
# macrophage     memory      naive 
#         44         88         88 


enrich$cell<-"T-cell"
enrich$cell[which(enrich$cell_type=="macrophage")]<-"Macrophage"
enrich$cell<-factor(enrich$cell,levels=c("T-cell","Macrophage"))

enrich$Category<-gsub("ChM_","",enrich$Category)

vec<-c("macrophage_6H_UNS","macrophage_24H_UNS",
"macrophage_6H_IFNG","macrophage_24H_IFNG",
"macrophage_24H_IL23",
"macrophage_6H_IL26","macrophage_24H_IL26",
"macrophage_6H_IL4","macrophage_24H_IL4",
"macrophage_6H_TNFA","macrophage_24H_TNFA",
"naive_16H_UNS","naive_D5_UNS",
"naive_16H_IFNB","naive_D5_IFNB",     
"naive_16H_IL10","naive_D5_IL10",
"naive_16H_IL21","naive_D5_IL21",
"naive_16H_IL27","naive_D5_IL27",
"naive_16H_ITREG","naive_D5_ITREG",
"naive_16H_TH0","naive_D5_TH0",
"naive_16H_TH1","naive_D5_TH1",
"naive_16H_TH17","naive_D5_TH17",
"naive_16H_TH2","naive_D5_TH2",
"naive_16H_TNFA","naive_D5_TNFA",
"memory_16H_UNS","memory_D5_UNS",
"memory_16H_IFNB","memory_D5_IFNB",
"memory_16H_IL10","memory_D5_IL10",
"memory_16H_IL21","memory_D5_IL21",    
"memory_16H_IL27","memory_D5_IL27",
"memory_16H_ITREG","memory_D5_ITREG",   
"memory_16H_TH0","memory_D5_TH0",
"memory_16H_TH1","memory_D5_TH1",     
"memory_16H_TH17","memory_D5_TH17",
"memory_16H_TH2","memory_D5_TH2",     
"memory_16H_TNFA","memory_D5_TNFA")

enrich$Category<-factor(enrich$Category,levels=vec)


# only colour the significant results:
enrich$signif<-enrich$cell
enrich$signif[which(enrich$Enrichment_p>threshold)]<-NA

p6<-ggplot(enrich, aes(y=Category, x=zscore,fill=signif)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) + 
  scale_fill_manual(limits = levels(enrich$signif),values=c(colo)) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
                 plot.margin = unit(c(1,1,1,1), "cm")) + facet_grid(cell ~ pheno, scales = "free_y", space = "free") + 
                 theme(legend.position="bottom",legend.title=element_blank())


ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_enrichment_analyses_soskic_zscore.png",
  p6,
  width = 140,
  height = ceiling(length(levels(as.factor(enrich$Category)))*3),
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)



enrich$source<-"soskic"
enrich_sos<-enrich
rm(enrich)

#######################
# encode:


files_enc<-files[grep("ENCFF",files)]
length(files_enc)
# [1] 1004

# maps:

dnase<-fread(paste0(path_gwas,"resources/ENCODE/encode_dnase_220250331_experiment_report.tsv"),head=T)
ihic<-fread(paste0(path_gwas,"resources/ENCODE/encode_intact_hic_20250331_experiment_report.tsv"),head=T)
atacs<-fread(paste0(path_gwas,"resources/ENCODE/encode_atacseq_20250331_experiment_report.tsv"),head=T)

mapecn<-rbind(dnase,ihic,atacs)
rm(dnase,ihic,atacs)

data<-gsub(".files_version_april2025.results","",files_enc)

for (i in 1:length(files_enc)) {

     if(file.exists(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_enc[i],sep=""))) {
      tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_enc[i],sep=""),head=T)
    } else {
      tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/other_traits/",files_enc[i],sep=""),head=T)
    }
    tmp<-tmp[which(tmp$Category=="L2_1",)]
    tmp$Category<-data[i]
    tmp$pheno<-toupper(gsub(".baseline.*","",tmp$Category))
    tmp$Category<-(gsub(".*baseline.","",tmp$Category))
    tmp$Category<-(gsub("_peaks","",tmp$Category))
    tmp$Category<-(gsub("_intact_hic","",tmp$Category))

    
    # get info from map
    tmp1<-mapecn[grep(tmp$Category,mapecn$Files),]

    tmp$assay<-tmp1$'Assay title'
    tmp$cell_type<-tmp1$'Biosample term name'
    tmp$cell_treatment<-tmp1$'Biosample treatment'

    if(i==1) {
        enrich<-tmp
    } else {
        enrich<-rbind(enrich,tmp)
    }
}

enrich$pheno<-factor(enrich$pheno,levels=c("CD","IBD","UC","HEIGHT"))
enrich$zscore<-enrich$'Coefficient_z-score'

table(enrich$assay,useNA="ifany")
  #  ATAC-seq   DNase-seq intact Hi-C 
  #       184         732          88 


table(enrich$assay,useNA="ifany")/4
  #  ATAC-seq   DNase-seq intact Hi-C 
  #        46         183          22 


table(enrich$cell_type)
#                                                     activated B cell 
#                                                                   12 
#                     activated CD4-positive, alpha-beta memory T cell 
#                                                                   12 
#                            activated CD4-positive, alpha-beta T cell 
#                                                                   76 
#  activated CD4-positive, CD25-positive, alpha-beta regulatory T cell 
#                                                                   12 
#                     activated CD8-positive, alpha-beta memory T cell 
#                                                                   20 
#                            activated CD8-positive, alpha-beta T cell 
#                                                                   20 
#            activated effector memory CD4-positive, alpha-beta T cell 
#                                                                    8 
#                                         activated gamma-delta T cell 
#                                                                    4 
#                      activated naive CD4-positive, alpha-beta T cell 
#                                                                   28 
#                      activated naive CD8-positive, alpha-beta T cell 
#                                                                   20 
#                                                     activated T-cell 
#                                                                   52 
#                                            activated T-helper 1 cell 
#                                                                    4 
#                                           activated T-helper 17 cell 
#                                                                    4 
#                                            activated T-helper 2 cell 
#                                                                    4 
#                                                               B cell 
#                                                                   12 
#                                               CD14-positive monocyte 
#                                                                   12 
#                               CD4-positive, alpha-beta memory T cell 
#                                                                   12 
#                                      CD4-positive, alpha-beta T cell 
#                                                                  100 
#            CD4-positive, CD25-positive, alpha-beta regulatory T cell 
#                                                                   16 
#                               CD8-positive, alpha-beta memory T cell 
#                                                                   12 
#                                      CD8-positive, alpha-beta T cell 
#                                                                   40 
#                       central memory CD4-positive, alpha-beta T cell 
#                                                                    4 
#                                                       dendritic cell 
#                                                                    4 
#                             effector CD4-positive, alpha-beta T cell 
#                                                                    4 
#                      effector memory CD4-positive, alpha-beta T cell 
#                                                                   12 
#                      effector memory CD8-positive, alpha-beta T cell 
#                                                                   16 
#                                         immature natural killer cell 
#                                                                    8 
#                                              inflammatory macrophage 
#                                                                   20 
#                                                        memory B cell 
#                                                                    4 
#                                                         naive B cell 
#                                                                   12 
#                 naive thymus-derived CD4-positive, alpha-beta T cell 
#                                                                   32 
#                 naive thymus-derived CD8-positive, alpha-beta T cell 
#                                                                   44 
#                                                  natural killer cell 
#                                                                   28 
#                 stimulated activated CD4-positive, alpha-beta T cell 
#                                                                   44 
#          stimulated activated CD8-positive, alpha-beta memory T cell 
#                                                                    4 
# stimulated activated effector memory CD8-positive, alpha-beta T cell 
#                                                                    4 
#                                   stimulated activated memory B cell 
#                                                                    4 
#                                    stimulated activated naive B cell 
#                                                                   12 
#           stimulated activated naive CD8-positive, alpha-beta T cell 
#                                                                   20 
#                                                suppressor macrophage 
#                                                                   32 
#                                                               T-cell 
#                                                                  156 
#                                                      T-helper 1 cell 
#                                                                   20 
#                                                     T-helper 17 cell 
#                                                                   16 
#                                                      T-helper 2 cell 
#                                                                   20 
#                                                      T-helper 9 cell 
#                                                                    4 



enrich$cell<-"T-cell"
enrich$cell[enrich$cell_type %in% c("B cell","activated B cell","memory B cell","naive B cell","stimulated activated memory B cell","stimulated activated naive B cell")]<-"B-cell"
enrich$cell[enrich$cell_type %in% c("CD14-positive monocyte")]<-"Monocyte"
enrich$cell[enrich$cell_type %in% c("dendritic cell")]<-"Dendritic cells"
enrich$cell[enrich$cell_type %in% c("immature natural killer cell","natural killer cell")]<-"NK"
enrich$cell[enrich$cell_type %in% c("inflammatory macrophage","suppressor macrophage")]<-"Macrophage"

table(enrich$cell_type,enrich$cell)

enrich$cell<-factor(enrich$cell,levels=c("B-cell","T-cell","Macrophage","Monocyte","NK","Dendritic cells"))


# only colour the significant results:
enrich$signif<-enrich$cell
enrich$signif[which(enrich$Enrichment_p>threshold)]<-NA

# duplicated cell types, safe the replicated ID:
enrich$cell_type_edited<-paste0(enrich$cell_type," (",enrich$Category,")")


atacs<-enrich[which(enrich$assay=="ATAC-seq"),]
dnase<-enrich[which(enrich$assay=="DNase-seq"),]
ihic<-enrich[which(enrich$assay=="intact Hi-C"),]


p7.1<-ggplot(atacs, aes(y=cell_type_edited, x=zscore,fill=signif)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) + 
  scale_fill_manual(limits = levels(enrich$signif),values=c(colo)) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
                 plot.margin = unit(c(1,1,1,1), "cm")) + facet_grid(cell ~ pheno, scales = "free_y", space = "free") + 
                 theme(legend.position="bottom",legend.title=element_blank()) + xlim(0,max(enrich$zscore))


p7.2<-ggplot(dnase, aes(y=cell_type_edited, x=zscore,fill=signif)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) + 
  scale_fill_manual(limits = levels(enrich$signif),values=c(colo)) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
                 plot.margin = unit(c(1,1,1,1), "cm")) + facet_grid(cell ~ pheno, scales = "free_y", space = "free") + 
                 theme(legend.position="none") + xlim(0,max(enrich$zscore))

p7.3<-ggplot(ihic, aes(y=cell_type_edited, x=zscore,fill=signif)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) + 
  scale_fill_manual(limits = levels(enrich$signif),values=c(colo)) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
                 plot.margin = unit(c(1,1,1,1), "cm")) + facet_grid(cell ~ pheno, scales = "free_y", space = "free") + 
                 theme(legend.position="bottom",legend.title=element_blank()) + xlim(0,max(enrich$zscore))



p7<-ggarrange(p7.2,
          ggarrange(p7.3,p7.1,nrow=2, labels = c("b","c"),heights=c(1,1.2),common.legend=T,legend="bottom",align="v"),
          ncol=2,labels = c("a"))

ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_enrichment_analyses_encode_zscore.png",
  p7,
  width = 240,
  height = ceiling(length(levels(as.factor(enrich$cell_type_edited)))*1.5),
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2.2
)



enrich$source<-"encode"
enrich_enc<-enrich
rm(enrich)



#### combine and save table

enrich_baseline$cell<-NA
enrich_baseline$cell_type<-NA

vec<-colnames(enrich_baseline)

enrich<-rbind(enrich_baseline,enrich_immune[,..vec])
enrich<-rbind(enrich,enrich_cis_element_atlas[,..vec])
# enrich<-rbind(enrich,enrich_levans[,..vec])
enrich<-rbind(enrich,enrich_hickey[,..vec])
enrich<-rbind(enrich,enrich_enc[,..vec])
enrich<-rbind(enrich,enrich_sos[,..vec])

nrow(enrich_baseline)/3
nrow(enrich_immune)/4
nrow(enrich_cis_element_atlas)/4
nrow(enrich_hickey)/4
nrow(enrich_enc)/4
nrow(enrich_sos)/4


write.table(enrich,"~/git/IIBDGC_GWAS/plots/paper_tables/Source_data_enrichment_analyses.tsv",col.names=T,row.names=F,quote=F,sep="\t")

q("no")

