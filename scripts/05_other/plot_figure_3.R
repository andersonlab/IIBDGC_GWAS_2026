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

rm(list=ls())

path_gwas="/path/to/ibdgwas/IIBDGC/"

#################################################################################################
# compare significant results:

enrich<-fread("~/git/IIBDGC_GWAS/plots/paper_tables/Source_data_enrichment_analyses.tsv",head=T)
enrich<-as.data.frame(enrich)

threshold<-0.05/((1+54+45+55+180+251+103)*3)
threshold

enrich$signif<-as.character(enrich$signif)
enrich$signif<-"signif"
enrich$signif[which(enrich$Enrichment_p>threshold)]<-"ns"


cd<-enrich[which(enrich$pheno=="CD" & enrich$source!="baseline_model"),]
uc<-enrich[which(enrich$pheno=="UC" & enrich$source!="baseline_model"),]
height<-enrich[which(enrich$pheno=="HEIGHT" & enrich$source!="baseline_model"),]

table(cd$Category==uc$Category)
# TRUE 
#  584

table(cd$Category==height$Category)
# TRUE 
#  584

uc$z<-uc$"Coefficient_z-score"
colnames(uc)[ncol(uc)]<-"z_uc"
cd$z<-cd$"Coefficient_z-score"
colnames(cd)[ncol(cd)]<-"z_cd"
height$z<-height$"Coefficient_z-score"
colnames(height)[ncol(height)]<-"z_height"


table(cd$Category==uc$Category)
# TRUE 
#  584

table(cd$Category==height$Category)
# TRUE 
#  584

comp<-cbind(cd[,c("Category","source","z_cd","cell","cell_type")],uc[,c("z_uc"),drop=F])
comp<-cbind(comp,height[,c("z_height"),drop=F])

comp<-comp[which(comp$source!="baseline_model"),]

dim(comp)
# 584

table(cd$signif)
  #   ns signif 
  #  228    356 
table(uc$signif)
  #   ns signif 
  #  155    429 

comp$signif<-"signif"
comp$signif[which(comp$Category %in% cd$Category[which(cd$signif=="ns")])]<-"ns"
comp$signif[which(comp$Category %in% uc$Category[which(uc$signif=="ns")])]<-"ns"

table(comp$signif)
  #   ns signif 
  #  229    355 
comp$signif[which((comp$Category %in% cd$Category[which(cd$signif=="ns")]) & (comp$Category %in% uc$Category[which(uc$signif!="ns")]))]<-"ns_cd"
comp$signif[which((comp$Category %in% uc$Category[which(uc$signif=="ns")]) & (comp$Category %in% cd$Category[which(cd$signif!="ns")]))]<-"ns_uc"

comp<-comp[which(comp$signif!="ns"),]


table(comp$signif)
#  ns_cd  ns_uc signif 
#     74      1    355 

comp$source[grep("ENC",comp$Category)]<-"Encode"
comp$source[which(comp$source=="soskic")]<-"Soskic 2019"
comp$source[which(comp$source=="hickey")]<-"Hickey 2023"
comp$source[which(comp$source=="immune_cell_atlas")]<-"Calderon 2019"
comp$source[which(comp$source=="cis_element_atlas")]<-"Zhang 2021"

comp$source<-factor(comp$source,levels=c("Hickey 2023","Zhang 2021","Calderon 2019","Soskic 2019","Encode"))

comp$cell<-as.factor(as.character(comp$cell))

names(table(comp$cell))
#  [1] "Adipocytes"          "B-cell"              "CyclingImmune"      
#  [4] "CyclingTA_1"         "CyclingTA_2"         "D-cell"             
#  [7] "Dendritic cells"     "Endothelial cell"    "Enterochromaffin"   
# [10] "Enterocyte"          "Enteroendocrine"     "Epithelial cell"    
# [13] "Exocrine"            "Fibroblast"          "Glia"               
# [16] "Goblet cell"         "I_Cells"             "ICC"                
# [19] "ILC"                 "K_Cells"             "L_Cells"            
# [22] "Macrophage"          "Mast"                "Mast cell"          
# [25] "Mo_Cells"            "Monocyte"            "Monocyte/Macrophage"
# [28] "Myofibroblast"       "Neurons"             "NK"                 
# [31] "Paneth"              "Pericytes"           "Plasma"             
# [34] "S-cell"              "Secretory cell"      "Stem"               
# [37] "T-cell"              "TA1"                 "TA2"                
# [40] "Tuft" 

## create a plot with fold changes 
comp$cell2<-NA

comp$cell2[grepl("D-cell|I_Cells|K_Cells|L_Cells|Mo_Cells|S-cell|Enteroendocrine",comp$cell)]<-"Enteroendocrine cells"
comp$cell2[grepl("NK",comp$cell)]<-"NK cells"
comp$cell2[grepl("Enterocyte",comp$cell)]<-"Enterocytes"
comp$cell2[grepl("Endothelial cell",comp$cell)]<-"Endothelial cell"
comp$cell2[grepl("Goblet cell",comp$cell)]<-"Goblet cells"
comp$cell2[grepl("T-cell",comp$cell)]<-"T cells"
comp$cell2[grepl("Monocyte|Monocyte/Macrophage|Macrophage",comp$cell)]<-"Myeloid cells"
comp$cell2[grepl("Tuft",comp$cell)]<-"Tuft cells"
comp$cell2[grepl("Paneth",comp$cell)]<-"Paneth cells"
comp$cell2[grepl("B-cell",comp$cell)]<-"B cells"
comp$cell2[grepl("Fibroblast|Myofibroblast",comp$cell)]<-"Fibroblasts"
comp$cell2[grepl("Dendritic cells",comp$cell)]<-"Dendritic cells"
comp$cell2[grepl("Secretory cell",comp$cell)]<-"Secretory cells, MUC"
comp$cell2[grepl("Stem",comp$cell)]<-"Stem cells"
comp$cell2[grepl("Enterochromaffin",comp$cell)]<-"Enterochromaffin cells"
comp$cell2[grepl("Epithelial",comp$cell)]<-"Epithelial cells"
comp$cell2[grepl("TA2|TA1|CyclingTA_1|CyclingTA_2",comp$cell)]<-"TA1 TA2 Enterocytes"
# We divided this differentiation trajectory into five cell types using similar annotations to other studies (stem > TA2 > TA1 > immature enterocyte > enterocyte)

table(comp$cell2,useNA="ifany")
  #                 B cells           Dendritic cells          Endothelial cell 
  #                      24                         4                         4 
  #  Enterochromaffin cells               Enterocytes     Enteroendocrine cells 
  #                       4                        26                        24 
  #        Epithelial cells Fibroblasts|Myofibroblast              Goblet cells 
  #                       3                        15                        16 
  #           Myeloid cells                  NK cells              Paneth cells 
  #                      33                        15                         6 
  #    Secretory cells, MUC                Stem cells                   T cells 
  #                       3                         8                       201 
  #     TA1 TA2 Enterocytes                Tuft cells                      <NA> 
  #                      16                         7                        21 

comp[which(is.na(comp$cell2)),]

# small N<5 exclude
comp<-comp[which(comp$cell2 %in% c("Enterocytes","Goblet cells","Paneth cells","Tuft cells","Enteroendocrine cells","B cells","T cells","Myeloid cells","NK cells",
"Stem cells","TA1 TA2 Enterocytes","Fibroblasts")),]

table(comp$cell2,useNA="ifany")
        #       B cells           Enterocytes Enteroendocrine cells 
        #            24                    26                    24 
        #   Fibroblasts          Goblet cells         Myeloid cells 
        #            15                    16                    33 
        #      NK cells          Paneth cells            Stem cells 
        #            15                     6                     8 
        #       T cells   TA1 TA2 Enterocytes            Tuft cells 
        #           201                    16                     7 
  
comp$cell2<-factor(comp$cell2,levels=c("Stem cells","TA1 TA2 Enterocytes","Enterocytes","Goblet cells","Paneth cells","Tuft cells","Enteroendocrine cells",
"B cells","T cells","Myeloid cells","NK cells","Fibroblasts"))
table(comp$cell2,useNA="ifany")
#             Stem cells   TA1 TA2 Enterocytes           Enterocytes 
#                     8                    16                    26 
#          Goblet cells          Paneth cells            Tuft cells 
#                    16                     6                     7 
# Enteroendocrine cells               B cells               T cells 
#                    24                    24                   201 
#         Myeloid cells              NK cells           Fibroblasts 
#                    33                    15                    15

tmp<-comp[which(comp$cell2=="Myeloid cells"),]
# binomial distribution driven by Hickey/Soskic vs Encode


comp<-comp[which(!is.na(comp$cell2)),]
table(comp$cell2)

table(comp$signif)
    

table(comp$cell2,comp$source,useNA="ifany")
#                         Hickey 2023 Zhang 2021 Calderon 2019 Soskic 2019 Encode
#   Stem cells                      8          0             0           0      0
#   TA1 TA2 Enterocytes            16          0             0           0      0
#   Enterocytes                    24          2             0           0      0
#   Goblet cells                   16          0             0           0      0
#   Paneth cells                    6          0             0           0      0
#   Tuft cells                      7          0             0           0      0
#   Enteroendocrine cells          24          0             0           0      0
#   B cells                         2          1             7           0     14
#   T cells                         5          2            30          44    120
#   Myeloid cells                   2          2             2          11     16
#   NK cells                        2          0             4           0      9
#   Fibroblasts                    15          0             0           0      0


minlim<-min(comp$z_cd,comp$z_uc)
maxlim<-max(comp$z_cd,comp$z_uc)

names(table(comp$cell2,useNA="ifany"))
#  [1] "Stem cells"            "TA1 TA2 Enterocytes"   "Enterocytes"          
#  [4] "Goblet cells"          "Paneth cells"          "Tuft cells"           
#  [7] "Enteroendocrine cells" "B cells"               "T cells"              
# [10] "Myeloid cells"         "NK cells"              "Fibroblasts"

   
comp$class<-NA
comp$class[which(comp$cell2 %in% c("B cells","T cells","Myeloid cells","NK cells","Dendritic cells"))]<-"Immune cells"
comp$class[which(comp$cell2 %in% c("Fibroblasts"))]<-"Fibroblasts"
comp$class[which(comp$cell2 %in% c("Stem cells","TA1 TA2 Enterocytes","Enterocytes","Secretory cells, MUC","Goblet cells","Paneth cells","Tuft cells","Enteroendocrine cells"))]<-"Epithelial cells"

# comp$class[which(comp$cell2 %in% c("Stem cells","TA1 TA2 Enterocytes","Enterocytes","Epithelial cells","Secretory cells, MUC","Goblet cells","Paneth cells","Tuft cells","Enteroendocrine cells","Enterochromaffin cells"))]<-"Epithelial cells"

# comp$class2<-comp$class
# comp$class2[grep("ileum",comp$Category)]<-"Epithelial cells (Ileum Enteroendocrine cells)"
# comp$class2<-as.factor(comp$class2)

table(comp$cell2,comp$class,useNA="ifany")
#                         Epithelial cells Fibroblasts Immune cells
#   Stem cells                           8           0            0
#   TA1 TA2 Enterocytes                 16           0            0
#   Enterocytes                         26           0            0
#   Goblet cells                        16           0            0
#   Paneth cells                         6           0            0
#   Tuft cells                           7           0            0
#   Enteroendocrine cells               24           0            0
#   B cells                              0           0           24
#   T cells                              0           0          201
#   Myeloid cells                        0           0           33
#   NK cells                             0           0           15
#   Fibroblasts                          0          15            0


immune<-comp[which(comp$class %in% c("Immune cells")),]
epi<-comp[which(comp$class %in% c("Epithelial cells")),]
fb<-comp[which(comp$class %in% c("Fibroblasts")),]

summary(immune$z_uc)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.233   3.161   4.382   4.718   6.486   8.858 
summary(immune$z_cd)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   -0.07166  2.72390  4.55390  4.53932  6.44440  9.90920

summary(epi$z_uc)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.990   3.269   3.413   3.313   3.520   4.253 
summary(epi$z_cd)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.02974 1.60640 2.13840 1.84298 2.22440 2.53720


summary(fb$z_uc)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  -0.1452  0.2502  1.3838  0.9632  1.4548  1.6557
summary(fb$z_cd)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1888  0.6103  1.7478  1.2847  1.7601  2.0495


non_ep<-comp[which(comp$class %in% c("Immune cells","Fibroblasts")),]

p8.1<-ggplot(non_ep, aes(y=z_uc, x=z_cd,color=class)) + 
  geom_point() +  
  geom_smooth(method='lm')   +
  theme_bw() + theme(axis.title.y=element_text(vjust=1.1),axis.title.x=element_blank(),text = element_text(size = 15),panel.background = element_blank(),legend = element_blank(),
                 axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(0.2,0.1,0.2,1), "cm")) + 
    labs(y="UC Coefficient Z-score") + 
    facet_grid(~cell2, scales = "free", space = "free") + 
    # xlim(-1.1,maxlim+1) + ylim(-1.1,maxlim+1) + 
    scale_x_continuous(limits = c(0,10), breaks = seq(0,10,by=2)) +
    scale_y_continuous(limits = c(0,10), breaks = seq(0,10,by=2)) +
    geom_abline(color="grey",linetype = "dashed") + 
    scale_colour_manual(values = c(
      "Epithelial cells" = "#BB5566",
      "Immune cells" = "#1b9e77",
      "Fibroblasts" = "#db7107"
    ),drop = FALSE)



p8.2<-ggplot(epi, aes(y=z_uc, x=z_cd,color=class)) + 
  geom_point() +  
  geom_smooth(method='lm',aes(y=z_uc, x=z_cd,color=class))  + 
  scale_colour_manual(values = c(
      "Epithelial cells" = "#DDAA33",
      "Immune cells" = "#1b9e77",
      "Fibroblasts" = "#db7107"
    )) +
  theme_bw() + theme(axis.title.y=element_text(vjust=3),text = element_text(size = 15),panel.background = element_blank(),legend = element_blank(),
                 axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(0.2,0.1,0.2,1.2), "cm")) + 
    labs(x="CD Coefficient Z-score",y="UC Coefficient Z-score") + 
    facet_grid(~cell2, scales = "free", space = "free") + 
    # xlim(-1.1,4.5) + ylim(-1.1,4.5) + 
    scale_x_continuous(limits = c(0,4), breaks = seq(0,4,by=1)) +
    scale_y_continuous(limits = c(0,4), breaks = seq(0,4,by=1)) +
    geom_abline(color="grey",linetype = "dashed")


top_row<-ggarrange(p8.1, NULL, ncol = 2, widths=c(7.3,2.7))
p<-ggarrange(top_row,p8.2,nrow=2,labels=c("a","b"),heights=c(0.92,1))

ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Figure_3_version_4.pdf",
  p,
  width = 170,
  height = 62,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)




# q("no")

