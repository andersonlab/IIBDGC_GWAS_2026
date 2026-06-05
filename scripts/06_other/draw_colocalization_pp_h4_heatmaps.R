# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# DRAW COLOC PP:

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=8000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(ggplot2)
library(ggpubr)
# library(gtools)
library(stringr)
library(rtracklayer)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_conditional_on_known_index_gwas_fm_2023_with_iibdgc_tier1_tier2_eas_sas_summary_stats_plus_unsupervised_conditional.tsv",sep=""))
all<-as.data.frame(all)

coloc<-fread(paste(path_gwas,"post_imputation/2022/analysis/colocalization/list_coloc_pph4_80_feb24.tsv.gz",sep=""),head=T)

# subset for conferences - ASHG, initial ECCO
# regions_prioritize<-c("4_155142873_156183022","6_108456705_109568306","8_51925612_52994334","10_13504179_14574118","10_71205481_72205481","11_1_808060","12_288145_1464614","17_82571123_83636754","22_42218922_43668958")

# subset for single cell:
# regions_prioritize<-c("9_3580755_6322802","5_131190525_132974928")

# subset for ecco:
# genes_prioritize<-c("SELL","ITGAL","IL23R","VSIR","WNK1","LACC1","STAT3","METRNL","PTPN2","ICAM5","LILRB3","CD40","NCF4","MFNG","TSPO","KPNA1","NFKB1","TNFSF15")
genes_prioritize<-c("BCL2","BCL2L11")

coloc<-coloc[which(coloc$gene  %in% genes_prioritize),]
genes_prioritize[which(!genes_prioritize %in% coloc$gene)]
#0

table(coloc$updated_region)
#  10_71205481_72205481 1_169049811_170662153     12_288145_1464614 
#                     4                     3                     4 
#  13_43410570_44410571  16_29971173_31965789   1_65786524_68054183 
#                     7                    28                     4 
#  17_38829907_43388282  17_82571123_83636754  18_12298655_13403349 
#                     3                     2                     2 
#  19_53703491_55371596   19_9626168_11103383  20_45605657_46620613 
#                     1                     2                     7 
#  22_36362461_37975679  22_42218922_43668958 3_121757685_123015485 
#                    27                    12                     7 
# 4_101580492_103494462 9_114022527_115448648 
#                     9                    20 

table(coloc$sample_group)
names(table(coloc$sample_group))

# remove some - to make plot more readable
# list_tissue_toexclude<-c("artery_aorta","artery_coronary","artery_tibial","esophagus_gej","esophagus_muscularis")
# coloc<-coloc[which(!coloc$sample_group %in% list_tissue_toexclude)]
list_tissue_toexclude<-""

# map to rename conditions:
# extracted from eQTL catalog V6
map<-fread("/path/to/project")
map$condition_name<-paste(map$study_id,map$dataset,sep="_")
map<-map[,c("condition_name","sample_group","study_label")]

# add map for macromap data:
map2<-as.data.frame(matrix(c("IFNB_6","IFNB_24","P3C_6","P3C_24","CIL_6","CIL_24","IFNG_6","IFNG_24","PIC_6","PIC_24","IL4_6","IL4_24","Prec_D0","Prec_D2","Ctrl_6","Ctrl_24","LIL10_6","LIL10_24","R848_6","R848_24","MBP_6","MBP_24","sLPS_6","sLPS_24"),ncol=1))
colnames(map2)<-"condition_name"
map2$sample_group<-paste("macrophage_iPSC_derived_",map2$condition_name,sep="")
map2$study_label<-"Panousis_2024"

map<-rbind(map,map2)
head(map)
rm(map2)

# add map for yazar:
map2<-as.data.frame(matrix(c("BIN","BMem","CD4ET","CD4NC","CD4SOX4","CD8ET","CD8NC","CD8S100B","DC","MonoC","MonoNC","NK","NKR","Plasma"),ncol=1))
colnames(map2)<-"condition_name"
map2$sample_group<-map2$condition_name
map2$study_label<-"Yazar_2022"

map<-rbind(map,map2)
rm(map2)


# add map for Hu:
map2<-as.data.frame(matrix(c("intestinal_mucosa_colon_ileum"),ncol=1))
colnames(map2)<-"condition_name"
map2$sample_group<-map2$condition_name
map2$study_label<-"Hu_2021"

map<-rbind(map,map2)
head(map)
rm(map2)

# add map for eQTL gen:
map2<-as.data.frame(matrix(c("Blood"),ncol=1))
colnames(map2)<-"condition_name"
map2$sample_group<-map2$condition_name
map2$study_label<-"eQTLGen_2021"

map<-rbind(map,map2)
head(map)
rm(map2)

# create a map of gene IDs
gtf<-rtracklayer::import(paste(path_gwas,'post_imputation/2022/analysis/metaanalysis/annotation/gencode.v39.annotation.gtf.gz',sep=""))
gene<-as.data.frame(gtf)
gene<-gene[,c("gene_id","gene_name")]
gene$gene<-gsub("\\.[0-9]*","",gene$gene_id)
gene<-gene[!duplicated(gene$gene),]
colnames(gene)<-c("gene_id","gene_symbol","gene_id_2")
rm(gtf)

# re-extract the data, just for those regions and genes, but withough PPH4 restriction
pheno<-c("ibd","cd","uc")

rm(coloc_all)
for (i in 1:length(pheno)) {

    tmp2<-fread(paste("/path/to/project",toupper(pheno[i]),".gz",sep=""),head=T)
    tmp2$gene_symbol<-tmp2$phenotype_id
    vec<-colnames(tmp2)

    tmp1<-fread(paste("/path/to/project",toupper(pheno[i]),".gz",sep=""),head=T)
    tmp1$gene_id<-gsub("\\.[0-9]{1,2}","",tmp1$phenotype_id)
    tmp1<-merge(tmp1,gene[,c("gene_id_2","gene_symbol")],by.x="gene_id",by.y="gene_id_2",all.x=T)
    tmp1<-tmp1[,..vec]

    tmp3<-fread(paste("/path/to/project",toupper(pheno[i]),".gz",sep=""),head=T)
    tmp3$gene_id<-tmp3$phenotype_id
    tmp3<-merge(tmp3,gene[,c("gene_id_2","gene_symbol")],by.x="gene_id",by.y="gene_id_2",all.x=T)
    tmp3<-tmp3[,..vec]

    tmp<-rbind(tmp1,tmp2,tmp3)
    rm(tmp1,tmp2,tmp3)

    tmp$MarkerName<-gsub("_",":",tmp$gwas_lead)
    
    if (pheno[i]=="cd") {
        tmp<-tmp[which(tmp$MarkerName %in% all$MarkerName[which(all$phenotype=="CD")] | tmp$MarkerName %in% all$MarkerName_ld[which(all$phenotype=="CD")] | tmp$MarkerName %in% all$MarkerName[which(all$phenotype=="IBD_saturated")] | tmp$MarkerName %in% all$MarkerName_ld[which(all$phenotype=="IBD_saturated")] | 
        tmp$MarkerName %in% all$MarkerName[which(all$phenotype=="IBD_unsaturated")] | tmp$MarkerName %in% all$MarkerName_ld[which(all$phenotype=="IBD_unsaturated")] ),]
    } else if (pheno[i]=="uc") {
        tmp<-tmp[which(tmp$MarkerName %in% all$MarkerName[which(all$phenotype=="UC")] | tmp$MarkerName %in% all$MarkerName_ld[which(all$phenotype=="UC")]  | tmp$MarkerName %in% all$MarkerName[which(all$phenotype=="IBD_saturated")] | tmp$MarkerName %in% all$MarkerName_ld[which(all$phenotype=="IBD_saturated")] | 
        tmp$MarkerName %in% all$MarkerName[which(all$phenotype=="IBD_unsaturated")] | tmp$MarkerName %in% all$MarkerName_ld[which(all$phenotype=="IBD_unsaturated")] ),]
    } else {
        tmp<-tmp[which(tmp$MarkerName %in% all$MarkerName | tmp$MarkerName %in% all$MarkerName_ld),]
    }

    tmp<-merge(tmp,map,by="condition_name",all.x=T)

    # only those subset above
    tmp<-tmp[which(tmp$gene_symbol %in% coloc$gene_symbol),]

    # only those subset above
    tmp<-tmp[which(tmp$MarkerName %in% coloc$MarkerName),]

    # only tissues subset above:
    tmp<-tmp[which(!tmp$sample_group %in% list_tissue_toexclude),]

    
    if(i==1) {
        coloc_all<-tmp
    } else {
        coloc_all<-rbind(coloc_all,tmp)
    }
}

coloc<-coloc_all
rm(coloc_all)

# reorganise names:
names(table(as.character(coloc$sample_group)))

coloc$condition_name_read<-NA
coloc$cell_tissue<-NA

coloc$condition_name_read[which(coloc$sample_group %in% c("Blood"))]<-"Blood (eQTLGen)"
coloc$condition_name_read[which(coloc$sample_group %in% c("blood"))]<-"Blood (GTEx)"
coloc$condition_name_read[which(coloc$sample_group %in% c("esophagus_mucosa"))]<-"Esophagus (GTEx)"
coloc$condition_name_read[which(coloc$sample_group %in% c("stomach"))]<-"Stomach (GTEx)"
coloc$condition_name_read[which(coloc$sample_group %in% c("small_intestine"))]<-"Small intestine (GTEx)"
coloc$condition_name_read[which(coloc$sample_group %in% c("colon_transverse","transverse_colon"))]<-"Colon transverse (GTEx)"
coloc$condition_name_read[which(coloc$sample_group %in% c("colon_sigmoid"))]<-"Colon sigmoid (GTEx)"
coloc$condition_name_read[which(coloc$sample_group %in% c("rectum"))]<-"Rectum (GTEx)"
coloc$condition_name_read[which(coloc$sample_group %in% c("artery_tibial"))]<-"Artery (GTEx)"
coloc$cell_tissue[which(coloc$sample_group %in% c("artery_tibial","blood","Blood","esophagus_mucosa","stomach","small_intestine","colon_transverse","transverse_colon","rectum","colon_sigmoid"))]<-"Tissue"


coloc$condition_name_read[which(coloc$sample_group %in% c("BMem"))]<-"Memory B cell (Yazar 2020)"
coloc$condition_name_read[which(coloc$sample_group %in% c("B-cell_CD19"))]<-"CD19 B cell"
coloc$condition_name_read[which(coloc$sample_group %in% c("B-cell_naive"))]<-"Naive B cell"
coloc$condition_name_read[which(coloc$sample_group %in% c("BIN"))]<-"Naive B cell (Yazar 2020)"
coloc$cell_tissue[which(coloc$sample_group %in% c("BMem","B-cell_CD19","B-cell_naive","BIN"))]<-"B cell"

coloc$condition_name_read[which(coloc$sample_group %in% c("fibroblast"))]<-"Fibroblast"
coloc$cell_tissue[which(coloc$sample_group %in% c("fibroblast"))]<-"Fb"

coloc$condition_name_read[which(coloc$sample_group %in% c("NK"))]<-"NK cell"
coloc$condition_name_read[which(coloc$sample_group %in% c("NK-cell_naive"))]<-"Naive NK cell"
coloc$cell_tissue[which(coloc$sample_group %in% c("NK","NK-cell_naive"))]<-"NK cell"

coloc$condition_name_read[which(coloc$sample_group %in% c("plasma"))]<-"plasma"
coloc$cell_tissue[which(coloc$sample_group %in% c("plasma"))]<-"Blood plasma"

coloc$condition_name_read[which(coloc$sample_group %in% c("iPSC"))]<-"iPSC naive"
coloc$cell_tissue[which(coloc$sample_group %in% c("iPSC"))]<-"iPSC"

coloc$condition_name_read[which(coloc$sample_group %in% c("platelet"))]<-"Platelet"
coloc$cell_tissue[which(coloc$sample_group %in% c("platelet"))]<-"Platelet"


coloc$condition_name_read[which(coloc$sample_group %in% c("CD4NC","CD4ET"))]<-"CD4 T cell (Yazar 2022)"
coloc$condition_name_read[which(coloc$sample_group %in% c("T-cell"))]<-"T cell"
coloc$condition_name_read[which(coloc$sample_group %in% c("T-cell_CD8","CD8ET","CD8NC"))]<-"CD8 T cell" 
coloc$condition_name_read[which(coloc$sample_group %in% c("CD8_T-cell_naive"))]<-"Naive CD8 T cell"
coloc$condition_name_read[which(coloc$sample_group %in% c("CD8_T-cell_anti-CD3-CD28"))]<-"Activated CD8 T cell" 
coloc$condition_name_read[which(coloc$sample_group %in% c("Tfh_memory"))]<-"Memory folicular T helper"
coloc$condition_name_read[which(coloc$sample_group %in% c("Th1_memory"))]<-"Memory T helper 1"
coloc$condition_name_read[which(coloc$sample_group %in% c("Th2_memory"))]<-"Memory T helper 2"
coloc$condition_name_read[which(coloc$sample_group %in% c("Th17_memory"))]<-"Memory T helper 17" 
coloc$condition_name_read[which(coloc$sample_group %in% c("Th1-17_memory"))]<-"Memory T helper 1/17" 
coloc$condition_name_read[which(coloc$sample_group %in% c("Treg_naive"))]<-"Naive Treg" 
coloc$condition_name_read[which(coloc$sample_group %in% c("Treg_memory"))]<-"Memory Treg"  
coloc$cell_tissue[which(coloc$sample_group %in% c("CD4NC","CD4ET","CD8ET","CD8NC","T-cell","CD8_T-cell_naive","CD8_T-cell_anti-CD3-CD28","T-cell_CD8",
"Tfh_memory","Th1-17_memory","Th1_memory","Th17_memory","Th2_memory","Treg_naive","Treg_memory"))]<-"T cell"

coloc$condition_name_read[which(coloc$sample_group %in% c("macrophage_IFNg+Salmonella"))]<-"Macrophage + IFNg Salmonella" 
coloc$condition_name_read[which(coloc$sample_group %in% c("macrophage_IFNg","macrophage_iPSC_derived_IFNG_6","macrophage_iPSC_derived_IFNG_24"))]<-"Macrophage + IFNg" 
coloc$condition_name_read[which(coloc$sample_group %in% c("macrophage_Salmonella"))]<-"Macrophage + Salmonella"
coloc$condition_name_read[which(coloc$sample_group %in% c("macrophage_Listeria"))]<-"Macrophage + Listeria"

coloc$condition_name_read[which(coloc$sample_group %in% c("macrophage_iPSC_derived_Ctrl_6","macrophage_iPSC_derived_Ctrl_24"))]<-"Macrophage Ctr" 
coloc$condition_name_read[which(coloc$sample_group %in% c("macrophage_iPSC_derived_MBP_24","macrophage_iPSC_derived_MBP_6"))]<-"Macrophage + MBP" 
coloc$condition_name_read[which(coloc$sample_group %in% c("macrophage_iPSC_derived_Prec_D0","macrophage_iPSC_derived_Prec_D2"))]<-"Macrophage Prec" 


coloc$condition_name_read[which(coloc$sample_group %in% c("macrophage_iPSC_derived_IFNB_6","macrophage_iPSC_derived_IFNB_24"))]<-"Macrophage + IFNb" 
coloc$condition_name_read[which(coloc$sample_group %in% c("macrophage_iPSC_derived_CIL_24","macrophage_iPSC_derived_CIL_6"))]<-"Macrophage + CIL" 
coloc$condition_name_read[which(coloc$sample_group %in% c("macrophage_iPSC_derived_IL4_24","macrophage_iPSC_derived_IL4_6"))]<-"Macrophage + IL4" 
coloc$condition_name_read[which(coloc$sample_group %in% c("macrophage_iPSC_derived_LIL10_24","macrophage_iPSC_derived_LIL10_6"))]<-"Macrophage + LIL10"                     
coloc$condition_name_read[which(coloc$sample_group %in% c("macrophage_iPSC_derived_P3C_24","macrophage_iPSC_derived_P3C_6"))]<-"Macrophage + P3C"
coloc$condition_name_read[which(coloc$sample_group %in% c("macrophage_iPSC_derived_PIC_24","macrophage_iPSC_derived_PIC_6"))]<-"Macrophage + PIC"
coloc$condition_name_read[which(coloc$sample_group %in% c("macrophage_iPSC_derived_R848_24","macrophage_iPSC_derived_R848_6"))]<-"Macrophage + R848"
coloc$condition_name_read[which(coloc$sample_group %in% c("macrophage_iPSC_derived_sLPS_24","macrophage_iPSC_derived_sLPS_6"))]<-"Macrophage + sLPS"
coloc$condition_name_read[which(coloc$sample_group %in% c("macrophage_naive"))]<-"Naive Macrophage"
coloc$cell_tissue[which(coloc$sample_group %in% c("macrophage_iPSC_derived_Ctrl_6","macrophage_iPSC_derived_Ctrl_24","macrophage_IFNg+Salmonella","macrophage_IFNg","macrophage_iPSC_derived_IFNG_6","macrophage_iPSC_derived_MBP_24","macrophage_iPSC_derived_MBP_6","macrophage_iPSC_derived_Prec_D0","macrophage_iPSC_derived_Prec_D2",
"macrophage_Salmonella","macrophage_Listeria","macrophage_iPSC_derived_IFNB_6","macrophage_iPSC_derived_IFNB_24","macrophage_iPSC_derived_IFNG_24",
"macrophage_iPSC_derived_CIL_24","macrophage_iPSC_derived_CIL_6","macrophage_iPSC_derived_IL4_24","macrophage_iPSC_derived_IL4_6",
"macrophage_iPSC_derived_LIL10_24","macrophage_iPSC_derived_LIL10_6","macrophage_iPSC_derived_MBP_24","macrophage_iPSC_derived_MBP_6",
"macrophage_iPSC_derived_P3C_24","macrophage_iPSC_derived_P3C_6","macrophage_iPSC_derived_PIC_24","macrophage_iPSC_derived_PIC_6",
"macrophage_iPSC_derived_R848_24","macrophage_iPSC_derived_R848_6","macrophage_iPSC_derived_sLPS_24","macrophage_iPSC_derived_sLPS_6",
"macrophage_naive"))]<-"Macrophage"

coloc$condition_name_read[which(coloc$sample_group %in% c("monocyte","monocyte_naive","MonoC"))]<-"Naive Monocyte"
coloc$condition_name_read[which(coloc$sample_group %in% c("monocyte_CD14"))]<-"CD14 Monocyte"
coloc$condition_name_read[which(coloc$sample_group %in% c("monocyte_IAV"))]<-"Monocyte + IAV"
coloc$condition_name_read[which(coloc$sample_group %in% c("monocyte_IFN24"))]<-"Monocyte + IFN"
coloc$condition_name_read[which(coloc$sample_group %in% c("monocyte_LPS","monocyte_LPS2","monocyte_LPS24"))]<-"Monocyte + LPS"
coloc$condition_name_read[which(coloc$sample_group %in% c("monocyte_Pam3CSK4"))]<-"Monocyte + Pam3CSK4"
coloc$condition_name_read[which(coloc$sample_group %in% c("monocyte_R848"))]<-"Monocyte + R848"
coloc$cell_tissue[which(coloc$sample_group %in% c("monocyte","monocyte_naive","monocyte_CD14","monocyte_IAV","monocyte_IFN24","monocyte_LPS","monocyte_LPS2",
"monocyte_LPS24","monocyte_Pam3CSK4","monocyte_R848","MonoC"))]<-"Monocyte"

coloc$condition_name_read[which(coloc$sample_group %in% c("neutrophil"))]<-"Naive Neutrophil"
coloc$condition_name_read[which(coloc$sample_group %in% c("neutrophil_CD16"))]<-"CD16 Neutrophil" 
coloc$condition_name_read[which(coloc$sample_group %in% c("neutrophil_CD15"))]<-"CD15 Neutrophil"
coloc$cell_tissue[which(coloc$sample_group %in% c("neutrophil","neutrophil_CD15","neutrophil_CD16"))]<-"Neutroph"

coloc$condition_name_read[which(coloc$sample_group %in% c("LCL_statin"))]<-"LCL + statin"
coloc$condition_name_read[which(coloc$sample_group %in% c("LCL_naive","LCL"))]<-"LCL naive"
coloc$cell_tissue[which(coloc$sample_group %in% c("LCL_statin","LCL_naive","LCL"))]<-"LCL"

table(coloc$sample_group,coloc$cell_tissue,useNA="ifany")
table(coloc$sample_group,coloc$condition_name_read,useNA="ifany")

coloc[which(is.na(coloc$sample_group)),]
coloc[which(is.na(coloc$condition_name_read)),]
coloc[which(is.na(coloc$cell_tissue)),]

coloc$condition_name_read<-paste(coloc$condition_name_read," (",coloc$study_label,")",sep="")

# sort out data:

coloc1<-merge(coloc,all[,c("MarkerName","updated_region","class","class_signal")],by="MarkerName")
coloc2<-merge(coloc,all[,c("MarkerName_ld","updated_region","class","class_signal")],by.x="MarkerName",by.y="MarkerName_ld")
coloc<-rbind(coloc1,coloc2)
rm(coloc1,coloc2)


coloc$updated_region<-as.factor(coloc$updated_region)


coloc$condition_name_read<-as.factor(coloc$condition_name_read)
coloc$cell_tissue<-as.factor(coloc$cell_tissue)
coloc$gene_name<-as.factor(coloc$gene_name)


names(table(coloc$condition_name_read))
# coloc$condition_name_read<-factor(coloc$condition_name_read, levels = c("Blood","Esophagus","Stomach","Colon transverse","Rectum",
# "Fibroblast","Naive B cell","CD19 B cell",
# "T cell", "CD8 T cell","Naive CD8 T cell","Activated CD8 T cell",
# "Naive Treg","Memory Treg",
# "Memory folicular T helper", "Memory T helper 1/17","Memory T helper 17","Memory T helper 2",
# "Naive Monocyte","CD14 Monocyte","Monocyte + IAV","Monocyte + IFN","Monocyte + LPS","Monocyte + Pam3CSK4","Monocyte + R848", 
# "Naive Neutrophil","CD15 Neutrophil","CD16 Neutrophil",
# "Naive Macrophage",
# "Macrophage + IFNb","Macrophage + IFNg","Macrophage + IL4",
# "Macrophage + IFNg Salmonella","Macrophage + Salmonella","Macrophage + Listeria",
# "Macrophage + LIL10","Macrophage + sLPS","Macrophage + CIL","Macrophage + P3C",
# "Macrophage + R848", "Macrophage + PIC"))


names(table(coloc$cell_tissue))
coloc$cell_tissue<-factor(coloc$cell_tissue, levels = c("Tissue","Fb","B cell","T cell","Monocyte","Neutroph","NK cell","Macrophage","Platelet","iPSC","LCL","Blood plasma"))
names(table(coloc$cell_tissue))



pdf("~/git/IIBDGC_GWAS/plots/metaanalysis/colocalization_pp_h4_heatmap.pdf",width=20, height=10)
ggplot(coloc, mapping = aes(x = condition_name_read,y = gene_symbol,fill = PP.H4.abf)) + 
  geom_tile() + xlab(label = "Sample") + scale_fill_continuous(type = "viridis",name = "H4 Posterior\nProbalility") +
  #geom_text(aes(label = round(PP.H4.abf, 2))) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=13),plot.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y = element_text(face = "italic",size=12),
        strip.text.x = element_text(size = 12),strip.text.y = element_text(size = 12)) + 
  facet_grid(updated_region ~ cell_tissue,scales = "free", space = "free") +
  theme(strip.text.y = element_text(angle = 0, hjust = 0),strip.background = element_blank(),
  legend.position="bottom",axis.title.x =element_blank(),axis.title.y=element_blank(), plot.margin = margin(1, 0.5, 1, 3, "cm"))
dev.off()



q("no")