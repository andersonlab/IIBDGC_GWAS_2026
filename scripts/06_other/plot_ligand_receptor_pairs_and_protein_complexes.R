# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#


# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=2000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R -q "cpu-interactive"

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(stringr)

rm(list=ls())
path_gwas<-"/path/to/ibdgwas/IIBDGC/"


# add direction of effect
eff<-fread(paste0(path_gwas,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_with_direction_effect.tsv.gz"))
eff<-eff[which(eff$list_genes!=""),]


eff$eff_direction<-NA
eff$eff_direction[which(eff$direction_effect_gsmr=="same" | eff$direction_effect_colocalization=="same")]<-"same"
eff$eff_direction[which(eff$direction_effect_gsmr=="opposite" | eff$direction_effect_colocalization=="opposite")]<-"opposite"
eff$eff_direction[which( (eff$direction_effect_gsmr=="same" & eff$direction_effect_colocalization=="opposite") | 
(eff$direction_effect_gsmr=="opposite" & eff$direction_effect_colocalization=="same") | 
(eff$direction_effect_gsmr=="opposite;same" | eff$direction_effect_colocalization=="opposite;same"))]<-"opposite;same"


table(eff$direction_effect_gsmr)
table(eff$direction_effect_gsmr,eff$eff_direction)


table(eff$direction_effect_colocalization)
table(eff$direction_effect_colocalization,eff$eff_direction)


eff$eff<-NA
eff$eff[which(eff$eff_direction=="same")]<-1
eff$eff[which(eff$eff_direction=="opposite")]<-(-1)
eff$eff[which(eff$eff_direction=="opposite;same")]<-0



#########################
## pairs ligands receptors - create a file with uniq combinations gene-ld_pair, using the file with additional drug information added by Kyle Gettler

reli<-fread("~/git/IIBDGC_GWAS/plots/paper_tables/Source_data_human_lr_pair_edited_with_OTAR_drug_evidence.txt")
reli<-as.data.frame(reli)
reli<-reli[which(reli$ligand_ibd_mr_coloc_exonic_coding==1 & reli$receptor_ibd_mr_coloc_exonic_coloc_coding==1),]

tmp1<-reli[,c("ligand_gene_symbol","lr_pair")]
colnames(tmp1)<-c("gene_symbol","lr_pair")

tmp2<-reli[,c("receptor_gene_symbol","lr_pair")]
colnames(tmp2)<-c("gene_symbol","lr_pair")

reli<-rbind(tmp1,tmp2)


reli<-merge(reli,eff[,c("list_genes","eff_direction","eff")],by.x='gene_symbol',by.y="list_genes")
reli<-reli[order(reli$lr_pair,decreasing=F),]

pairs_to_keep<-intersect(tmp1$lr_pair,tmp2$lr_pair)

reli<-reli[which(reli$lr_pair %in% pairs_to_keep),]
table(reli$lr_pair)
    # ADAM15_ITGAV      ADAM17_MUC1       CCL13_CCR2       CCL13_CCR5 
    #            2                2                2                2 
    #    CCL2_CCR2        CCL2_CCR5      CXCL5_CXCR1      CXCL5_CXCR2 
    #            2                2                2                2 
    #   GNAS_ADCY7      ICAM4_ITGA4      ICAM4_ITGAL      ICAM4_ITGAV 
    #            2                2                2                2 
    #  ICAM5_ITGAL      IL12B_IL23R      PDGFB_ITGAV       PLAU_ITGAV 
    #            2                2                2                2 
    #   PLAU_PLAUR TNFSF15_TNFRSF25      VCAM1_ITGA4       VEGFB_NRP1 
    #            2                2                2                2 
 


table(reli$eff_direction,reli$lr_pair,useNA="ifany")     
  #               ADAM15_ITGAV ADAM17_MUC1 CCL13_CCR2 CCL13_CCR5 CCL2_CCR2
  # opposite                 2           0          1          2         1
  # opposite;same            0           1          0          0         0
  # same                     0           1          1          0         1
  # <NA>                     0           0          0          0         0
               
  #               CCL2_CCR5 CXCL5_CXCR1 CXCL5_CXCR2 GNAS_ADCY7 ICAM4_ITGA4
  # opposite              2           2           2          1           0
  # opposite;same         0           0           0          0           0
  # same                  0           0           0          0           2
  # <NA>                  0           0           0          1           0
               
  #               ICAM4_ITGAL ICAM4_ITGAV ICAM5_ITGAL IL12B_IL23R PDGFB_ITGAV
  # opposite                0           1           0           0           2
  # opposite;same           0           0           0           0           0
  # same                    2           1           2           2           0
  # <NA>                    0           0           0           0           0
               
  #               PLAU_ITGAV PLAU_PLAUR TNFSF15_TNFRSF25 VCAM1_ITGA4 VEGFB_NRP1
  # opposite               1          0                0           0          1
  # opposite;same          0          0                1           0          0
  # same                   1          2                1           2          1
  # <NA>                   0          0                0           0          0



#########################
## protein complexes, exclude those that are already been addressed in th 

prc<-fread("~/git/IIBDGC_GWAS/plots/paper_tables/Source_data_human_protein_complexes_corum.txt.gz")
table(prc$n_ibd_genes_complex)
#    0    1    2    3 
# 4733  692   79    3 

prc<-prc[which(prc$n_ibd_genes_complex>1),]
table(prc$n_ibd_genes_complex)
#  2  3 
# 79  3

# extract the list of genes that contribute to those protein complexes
lgprc<-fread(paste0(path_gwas,"post_imputation/2022/analysis/final_tables/protein_complexes_ibd_genes_corum.txt.gz"),head=T)
lgprc<-lgprc[which(lgprc$protein_complex %in% prc$protein_complex),]
lgprc<-lgprc[which(lgprc$ibd_gene==1),]

dim(lgprc)
lgprc<-lgprc[!duplicated(lgprc),]
dim(lgprc)

dim(lgprc[!duplicated(lgprc$protein_complex),])
# [1] 82  7
dim(lgprc[!duplicated(lgprc$list_genes),])
# [1] 80  7


table(lgprc$n_ibd_genes_complex)
#   2   3 
# 158   9

table(lgprc$n_genes_complex)

# the following step excludes 3 complexes with 3 IBD genes each one:
#    list_genes ibd_gene n_ibd_genes_complex n_genes_complex
#        <char>    <int>               <int>           <int>
# 1:       RPS9        1                   3             104
# 2:      SNU13        1                   3             104
# 3:       TUFM        1                   3             104
# 4:      PHF5A        1                   3             112
# 5:        PNN        1                   3             112
# 6:     ZNF207        1                   3             112
# 7:     GEMIN2        1                   3             129
# 8:      PHF5A        1                   3             129
# 9:       RBM7        1                   3             129

lgprc<-lgprc[which(lgprc$n_genes_complex<5),]



dim(lgprc)
lgprc<-merge(lgprc,eff[,c("list_genes","eff_direction","eff")],by.x='list_genes',by.y="list_genes")
dim(lgprc)
# [1] 104 13

lgprc<-lgprc[!duplicated(lgprc),]
table(lgprc$n_ibd_genes_complex)
#   2 
# 104

list_uniq_genes<-names(table(lgprc$protein_complex)[which(table(lgprc$protein_complex)==1)])
list_uniq_genes
# 0
lgprc<-lgprc[which(!lgprc$protein_complex %in% c(list_uniq_genes)),]


lgprc<-lgprc[order(lgprc$protein_complex,decreasing=F),]

length(names(table(lgprc$protein_complex)))
# [1] 52

dim(lgprc)
# [1] 104  15

# retain only those that are not already represented in the ligand receptor panel:
lgprc$protein_complex<-gsub(";","_",lgprc$protein_complex)
lgprc<-lgprc[which(!lgprc$protein_complex %in% reli$lr_pair),]
dim(lgprc)
# [1] 102  15

# only show once pairs that are already represented in another figure:
lgprc<-lgprc[order(lgprc$protein_complex,lgprc$list_genes,decreasing=F),]

ids<-names(table(lgprc$protein_complex))

lgprc$id2<-NA
for (i in 1:length(ids)) {

  lgprc$id2[which(lgprc$protein_complex==ids[i])]<-paste(lgprc$list_genes[which(lgprc$protein_complex==ids[i])],collapse="_")

}

# exclude again those present in the lr_pairs
lgprc<-lgprc[which(!lgprc$id2 %in% reli$lr_pair),]
dim(lgprc)
# [1] 90 17

table(lgprc$id2)
    #   CCR2_CCR5      CD81_ITGA4     CEP19_CEP43       CRK_ELMO1       CUL2_ELOC 
    #           2               2               2               2               2 
    # CXCR1_CXCR2     CXCR2_PLCB3       ETS1_ETS2     FOXO1_SMAD3     FOXO3_SMAD3 
    #           2               2               2               2               4 
    # GNA12_ITGAV       GNAS_GNB2 GUCY1A1_GUCY1B1      IL6ST_JAK2     IL6ST_STAT3 
    #           2              22               4               2               2 
    # ITGAV_ITGB8     ITGAV_PLAUR     ITGAV_PLPP3     ITGAV_PTK2B         MAX_MYC 
    #           2               2               2               2               4 
    #   NCF1_NCF4   NDFIP1_UBE2L3       ORC3_ORC4       RBPJ_SPEN 
    #           2               2               4               6 


# exclude as well internal duplicates:
ids2<-names(table(lgprc$id2))

for (i in 1:length(ids2)) {


  if(nrow(lgprc[which(lgprc$id2==ids2[i])])>2) {
    
    tmp<-as.data.frame(names(table(lgprc$protein_complex[which(lgprc$id2==ids2[i])])))
    colnames(tmp)<-"protein_complex"
    tmp$nchar<-nchar(tmp$protein_complex)
    print(i)
    print(tmp)

    if (nrow(tmp[which(tmp$nchar==min(tmp$nchar)),])>1) {

      tmp<-tmp[order(tmp$nchar,decreasing=F),]
      lgprc$protein_complex[which(lgprc$protein_complex %in% tmp$protein_complex[2:nrow(tmp)])]<-"remove"

    }

    lgprc$protein_complex[which(lgprc$protein_complex %in% tmp$protein_complex[which(tmp$nchar!=min(tmp$nchar))])]<-"remove"
    
    rm(tmp)

  }

}
lgprc<-lgprc[which(lgprc$protein_complex!="remove"),]
dim(lgprc)
# [1] 56 10

table(lgprc$id2)
    #   CCR2_CCR5      CD81_ITGA4     CEP19_CEP43       CRK_ELMO1       CUL2_ELOC 
    #           2               2               2               2               2 
    # CXCR1_CXCR2     CXCR2_PLCB3       ETS1_ETS2     FOXO1_SMAD3     FOXO3_SMAD3 
    #           2               2               2               2               2 
    # GNA12_ITGAV       GNAS_GNB2 GUCY1A1_GUCY1B1     IL2RA_IL2RB      IL6ST_JAK2 
    #           2               2               2               2               2 
    # IL6ST_STAT3     ITGA4_VCAM1     ITGAV_ITGB8     ITGAV_PLAUR     ITGAV_PLPP3 
    #           2               2               2               2               2 
    # ITGAV_PTK2B       JAK2_TYK2         MAX_MYC       NCF1_NCF4   NDFIP1_UBE2L3 
    #           2               2               2               2               2 
    #  NOD2_RIPK2       ORC3_ORC4       RBPJ_SPEN 
    #           2               2               2 


####################################################
# create an unique file to share as source data:

tmp1<-reli[,c("lr_pair","gene_symbol","eff_direction","eff")]
colnames(tmp1)<-c("protein_complex_lr_pair","gene_symbol","eff_direction","eff")
tmp1$class<-"ligand_receptor_pair"
tmp2<-lgprc[,c("protein_complex","list_genes","eff_direction","eff")]
colnames(tmp2)<-c("protein_complex_lr_pair","gene_symbol","eff_direction","eff")
tmp2$class<-"protein_complex"

tmp<-rbind(tmp1,tmp2)

fwrite(tmp,"/path/to/user/git/IIBDGC_GWAS/plots/paper_tables/Source_data_ligand_receptors_protein_complexes.tsv",col.names=T,row.names=F,quote=F,sep="\t")



# visualise all the results:

p1<-ggplot(tmp[which(tmp$class=="ligand_receptor_pair"),], aes(y = gene_symbol, x=1, fill=eff)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + theme_bw()  + scale_fill_gradient2(low = 'darkblue', mid = 'white', high = 'darkred') +
  # geom_text(aes(label = eff_direction), size = 12 / .pt) +
    facet_grid(rows=vars(protein_complex_lr_pair), scales = "free",space = "free") + theme(legend.position="none",strip.text.y.right = element_text(angle = 0),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.title.y=element_blank())


p2<-ggplot(tmp[which(tmp$class=="protein_complex"),], aes(y = gene_symbol, x=1, fill=eff)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + theme_bw()  + scale_fill_gradient2(low = 'darkblue', mid = 'white', high = 'darkred') +
  # geom_text(aes(label = eff_direction), size = 12 / .pt) +
    facet_grid(rows=vars(protein_complex_lr_pair), scales = "free",space = "free") + theme(legend.position="none",strip.text.y.right = element_text(angle = 0),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.title.y=element_blank())





# exclude for the moment those pairs with coding, where direction of effect is less clear
to_remove<-tmp$protein_complex_lr_pair[which(is.na(tmp$eff))]
tmp<-tmp[which(!tmp$protein_complex_lr_pair %in% to_remove),]

# retain for the plot only those pairs with the same direction of effect:

ids<-names(table(tmp$protein_complex_lr_pair))

tmp$label<-NA

for (i in 1:length(ids)) {

  if (tmp$eff[which(tmp$protein_complex_lr_pair==ids[i])][1]!=tmp$eff[which(tmp$protein_complex_lr_pair==ids[i])][2]) {
    tmp$label[which(tmp$protein_complex_lr_pair==ids[i])]<-"remove"
  }

}

tmp<-tmp[which(is.na(tmp$label)),]
tmp<-tmp[which(tmp$eff_direction!="opposite;same"),]


p1<-ggplot(tmp[which(tmp$class=="ligand_receptor_pair"),], aes(y = gene_symbol, x=1, fill=eff)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + theme_bw()  + scale_fill_gradient2(low = 'darkblue', mid = 'white', high = 'darkred') +
  # geom_text(aes(label = eff_direction), size = 12 / .pt) +
    facet_grid(rows=vars(protein_complex_lr_pair), scales = "free",space = "free") + theme(legend.position="none",strip.text.y.right = element_text(angle = 0),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.title.y=element_blank())


p2<-ggplot(tmp[which(tmp$class=="protein_complex"),], aes(y = gene_symbol, x=1, fill=eff)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + theme_bw()  + scale_fill_gradient2(low = 'darkblue', mid = 'white', high = 'darkred') +
  # geom_text(aes(label = eff_direction), size = 12 / .pt) +
    facet_grid(rows=vars(protein_complex_lr_pair), scales = "free",space = "free") + theme(legend.position="none",strip.text.y.right = element_text(angle = 0),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.title.y=element_blank())


p2_bottom_margin <- p2 +
    theme(plot.margin = margin(5.5, 5.5, 40, 5.5))


# p<-ggarrange(p1,p2_bottom_margin,ncol=2,widths=c(0.85,1.1),labels = c("A", "B"))
p<-ggarrange(p1,p2_bottom_margin,ncol=2,widths=c(0.85,1.1))


ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Figure_4a_ligand_receptor_protein_complex_direction.png",
  p,
  width = 60,
  height = 60,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)


# LOAD THE FILE WITH THE PROTEIN COMPLEXES:
pcom<-fread("/path/to/ibdgwas/IIBDGC/resources/corum/corum_humanComplexes.txt",head=T)
pcom<-as.data.frame(pcom)
dim(pcom)

pcom$protein_complex<-gsub(";","_",pcom$subunits_gene_name)

pcom[which(pcom$protein_complex %in% tmp$protein_complex),]
#        protein_complex
#                 <char>
#  1:        CXCR1_CXCR2
#  2:        CXCR1_CXCR2
#  3: CXCR2_NHERF1_PLCB3
#  4: CXCR2_NHERF1_PLCB3
#  5:  FOXO3_SMAD3_SMAD4
#  6:  FOXO3_SMAD3_SMAD4
#  7:  GNA12_ITGAV_P2RY2
#  8:  GNA12_ITGAV_P2RY2
#  9:    GUCY1A1_GUCY1B1
# 10:    GUCY1A1_GUCY1B1
# 11:  ITGAV_ITGB3_PLPP3
# 12:  ITGAV_ITGB3_PLPP3
# 13:     NCF1_NCF2_NCF4
# 14:     NCF1_NCF2_NCF4
# 15:     ORC2_ORC3_ORC4
# 16:     ORC2_ORC3_ORC4
# 17:          RBPJ_SPEN
# 18:          RBPJ_SPEN
#                                                                                                                                                                      protein_complex_drug
#                                                                                                                                                                                    <char>
#  1:                                                                                                                                                                                      
#  2:                                                                                                                                                                                      
#  3:                                                                                                                                                                                      
#  4:                                                                                                                                                                                      
#  5:                                                                                                                                                                                      
#  6:                                                                                                                                                                                      
#  7:                                                                                                                                                                            diquafosol
#  8:                                                                                                                                                                            diquafosol
#  9: isosorbide dinitrate;isosorbide mononitrate;nitric oxide;glyceryl trinitrate;nitroprusside;isosorbide dinitrate;isosorbide mononitrate;nitric oxide;glyceryl trinitrate;nitroprusside
# 10: isosorbide dinitrate;isosorbide mononitrate;nitric oxide;glyceryl trinitrate;nitroprusside;isosorbide dinitrate;isosorbide mononitrate;nitric oxide;glyceryl trinitrate;nitroprusside
# 11:                                                                                                                                                      eptifibatide;tirofiban;abciximab
# 12:                                                                                                                                                      eptifibatide;tirofiban;abciximab
# 13:                                                                                                                                                                                      
# 14:                                                                                                                                                                                      
# 15:                                                                                                                                                                                      
# 16:                                                                                                                                                                                      
# 17:                                                                                                                                                                                      
# 18:                        

lgprc[which(lgprc$protein_complex %in% tmp$protein_complex & lgprc$protein_complex_drug!=""),c("list_genes","protein_complex","protein_complex_drug")]
lgprc[which(lgprc$protein_complex %in% tmp$protein_complex & lgprc$protein_complex_drug!=""),c("list_genes","protein_complex","protein_complex_drug")]
#    list_genes   protein_complex
#        <char>            <char>
# 1:      GNA12 GNA12_ITGAV_P2RY2
# 2:      ITGAV GNA12_ITGAV_P2RY2
# 3:    GUCY1A1   GUCY1A1_GUCY1B1
# 4:    GUCY1B1   GUCY1A1_GUCY1B1
# 5:      ITGAV ITGAV_ITGB3_PLPP3
# 6:      PLPP3 ITGAV_ITGB3_PLPP3
#                                                                                                                                                                     protein_complex_drug
#                                                                                                                                                                                   <char>
# 1:                                                                                                                                                                            diquafosol -> targets P2RY2
# 2:                                                                                                                                                                            diquafosol -> targets P2RY2
# 3: isosorbide dinitrate;isosorbide mononitrate;nitric oxide;glyceryl trinitrate;nitroprusside;isosorbide dinitrate;isosorbide mononitrate;nitric oxide;glyceryl trinitrate;nitroprusside
# 4: isosorbide dinitrate;isosorbide mononitrate;nitric oxide;glyceryl trinitrate;nitroprusside;isosorbide dinitrate;isosorbide mononitrate;nitric oxide;glyceryl trinitrate;nitroprusside
# 5:                                                                                                                                                      eptifibatide;tirofiban;abciximab
# 6:                                                                                                                                                      eptifibatide;tirofiban;abciximab


# encodes Lipid Phosphate Phosphatase 3 (LPP3), a cell-surface enzyme crucial for regulating bioactive lipid signals like


q("no")
