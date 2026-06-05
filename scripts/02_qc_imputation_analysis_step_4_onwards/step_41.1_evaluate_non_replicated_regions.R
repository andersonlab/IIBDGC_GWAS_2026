# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
############################################################################################################################
############################################################################################################################


# 1.- EVALUATE THE REGIONS WHERE WE CANNOT IDENTIFY ANY COJO UNSUPERVISED GENOME-WIDE SIGNIFICANT SIGNALS, OR GENOME-WIDE SIGNIFICANT SIGNALS IN MULTIANCESTRY TIER 2 DATASET

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=15000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

library(data.table)
library(ggplot2)
library(colorspace)
library(dplyr)
library(tibble)
library(ggrepel)
library(ggpubr)
rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

all_recap<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_no_subset_final_2.tsv.gz",sep=""))
dim(all_recap)
# [1] 631 163

reg<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_regions_420_cojo_supervised_plus_unsupervised.tsv.gz",sep=""))
reg<-reg[which(!reg$updated_region %in% all_recap$updated_region),]
dim(reg)
# [1] 49  2

all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_regions_no_conditional_significant_eur_tier2_no_fm_imput.tsv.gz",sep=""))
all<-all[which(all$updated_region %in% reg$updated_region)]
dim(all)
# [1]  47 168

# add the summary stats from de Lange:
pheno<-c("cd","uc","ibd")

for (ph in pheno) {

  # delange:

  tmp<-fread(paste0(path_gwas,"summary_files/summary_stats_",ph,"_delange_sorted_noheader_chr_pos_b37_pval_lifted_hg38.gz",sep=""),head=F)
  tmp<-tmp[,1:5]
  colnames(tmp)<-c("chr","position_b38","end","ID_b37","pvalue")
  tmp<-tmp[,c("chr","position_b38","ID_b37","pvalue")]

  tmp<-tmp[,c("ID_b37","chr","position_b38","pvalue")]

  colnames(tmp)[ncol(tmp)]<-paste(colnames(tmp)[ncol(tmp)],"_delange_",ph,sep="")

  if (!exists("old_delange")) {
    old_delange<-tmp
  } else {
    old_delange<-merge(old_delange,tmp[,c(1,4)],by="ID_b37",all=T)
  }
}


old_delange$ref<-gsub("[0-9]{1,2}:[0-9]*_","",old_delange$ID_b37)
old_delange$alt<-gsub(".*_","",old_delange$ref)
old_delange$ref<-gsub("_.*","",old_delange$ref)

old_delange$MarkerName_1<-paste(old_delange$chr,old_delange$position_b38,old_delange$ref,old_delange$alt,sep=":")
old_delange$MarkerName_2<-paste(old_delange$chr,old_delange$position_b38,old_delange$alt,old_delange$ref,sep=":")

old_delange1<-old_delange[which(old_delange$MarkerName_1 %in% all$MarkerName),]
old_delange2<-old_delange[which(old_delange$MarkerName_2 %in% all$MarkerName),]

dim(old_delange1)
# [1] 39 10
dim(old_delange2)
# [1]  0 10

old_delange1<-old_delange1[,c("ID_b37","chr","position_b38","pvalue_delange_cd","pvalue_delange_uc","pvalue_delange_ibd","ref","alt","MarkerName_1")]
old_delange2<-old_delange2[,c("ID_b37","chr","position_b38","pvalue_delange_cd","pvalue_delange_uc","pvalue_delange_ibd","ref","alt","MarkerName_2")]

colnames(old_delange1)<-c("ID_b37","chr","position_b38","pvalue_delange_cd","pvalue_delange_uc","pvalue_delange_ibd","ref","alt","MarkerName")
colnames(old_delange2)<-c("ID_b37","chr","position_b38","pvalue_delange_cd","pvalue_delange_uc","pvalue_delange_ibd","ref","alt","MarkerName")

old_delange<-rbind(old_delange1,old_delange2)
rm((old_delange1,old_delange2))

all<-merge(all,old_delange[,c("MarkerName","pvalue_delange_cd","pvalue_delange_uc","pvalue_delange_ibd")],by="MarkerName",all.x=T)

all$gwas_signif_delange<-"N"
all$gwas_signif_delange[which(all$pvalue_delange_cd<=5E-8 | all$pvalue_delange_ibd<=5E-8 | all$pvalue_delange_uc<=5E-8)]<-"Y"

table(all$'FM study',all$gwas_signif_delange)
#            N  Y
#            8  0
#   delange  0  5
#   Liu     17  5
#   noFM    12  0


# add the summary stats from de Liu Liu Gao et al.):

pheno<-c("cd","uc","ibd")

rm(old_liu)
for (ph in pheno) {

  # liu liu gao:

  tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/eas_iibdgc/liu-2022-east-asian-gwas/summary-stats/ibd_EAS_EUR_SiKJEF_meta_",toupper(ph),".TBL.txt.gz",sep=""),head=T)
  tmp<-tmp[,c("MarkerName","Allele1","Allele2","CHR","BP","P-value")]
  tmp<-tmp[which(tmp$MarkerName %in% all$MarkerName),]
  
  tmp<-tmp[,c("MarkerName","CHR","BP","P-value")]
  colnames(tmp)<-c("MarkerName","chr","position_b38","pvalue")
  
  colnames(tmp)[ncol(tmp)]<-paste(colnames(tmp)[ncol(tmp)],"_liu_",ph,sep="")

  if (!exists("old_liu")) {
    old_liu<-tmp
  } else {
    old_liu<-merge(old_liu,tmp[,c(1,4)],by="MarkerName",all=T)
  }
}

dim(old_liu)
# [1] 46  6

all<-merge(all,old_liu[,c("MarkerName","pvalue_liu_cd","pvalue_liu_uc","pvalue_liu_ibd")],by="MarkerName",all.x=T)


all$gwas_signif_liu<-"N"
all$gwas_signif_liu[which(all$pvalue_liu_cd<=5E-8 | all$pvalue_liu_ibd<=5E-8 | all$pvalue_liu_uc<=5E-8)]<-"Y"

table(all$'FM study',all$gwas_signif_liu)
  #          N  Y
  #          8  0
  # delange  2  3
  # Liu      7 15
  # noFM    11  1

all$old_region_signal_name

# keep all signals identified by Liu as Liu, all the rest as de Lange:

table(all$region_genomewide_significance_eur_eas_sas_tier_2)
#  N  Y 
# 42  5 

dim(all)
# [1]  47 225

dim(all[!duplicated(all$updated_region),])
# [1]  47 176

# region 20_45605657_46620613 with two signals


### DESCRIBE REGIONS THAT NOT REPLICATE IN EUR_SAS_EAS TIER 2:

ma<-all[which(all$region_genomewide_significance_eur_eas_sas_tier_2=="Y"),]
dim(ma)
# [1]  5 225

table(ma$'FM study')
      #   delange     Liu 
      # 3       1       1

table(ma$class_signal)
# new_cojo_supervised                 old 
#                   3                   2 

table(ma$'FM study',ma$class_signal)
#           old
#   delange   1
#   Liu       1


# is any of these variants in LD with any variant prioritise via FM or supervised conditional analyses??


### RETAIN ONLY SIGNALS/REGIONS THAT DO NOT REPLICATE IN EUR_SAS_EAS TIER 2:
all<-all[which(all$region_genomewide_significance_eur_eas_sas_tier_2=="N"),]
dim(all)
# [1]  42 225

table(all$'FM study')
      #   delange     Liu    noFM 
      # 5       4      21      12 

all$min_pval_iibdgc_multiancestry_tier_2<-pmin(all$"P-value_cd_eur_tier_2_eas_sas",all$"P-value_uc_eur_tier_2_eas_sas",all$"P-value_ibd_eur_tier_2_eas_sas")
all$min_pval_iibdgc_eur_tier_2<-pmin(all$"P-value_cd_eur_tier_2",all$"P-value_uc_eur_tier_2",all$"P-value_ibd_eur_tier_2")


all$min_pval_liu<-pmin(all$pvalue_liu_ibd,all$pvalue_liu_cd,all$pvalue_liu_uc)
all$min_pval_delange<-pmin(all$pvalue_delange_ibd,all$pvalue_delange_cd,all$pvalue_delange_uc)

all$region_genomewide_significance_eur_eas_sas_tier_2<-as.factor(all$region_genomewide_significance_eur_eas_sas_tier_2)
all$gwas_signif_delange<-as.factor(all$gwas_signif_delange)


col_all<-c("grey","blue")


################################################
## VARIANTS WE CAN ASSESS WITH DE LANGE:

dl<-all[which(all$'FM study'!="Liu"),]


dl$variant_id<-dl$MarkerName
dl$variant_id[which(dl$region_genomewide_significance_eur_eas_sas_tier_2=="N" & dl$gwas_signif_delange=="N")]<-""
dl$variant_id
#  [1] "chr10:124750812:G:A" ""                    "chr10:26890667:T:G" 
#  [4] ""                    ""                    ""                   
#  [7] ""                    ""                    ""                   
# [10] ""                    ""                    ""                   
# [13] ""                    ""                    "chr1:209797265:G:A" 
# [16] "chr20:6113242:T:G"   ""                    ""                   
# [19] ""                    ""                    ""     

dl[which(dl$variant_id!=""),c("variant_id","min_pval_iibdgc_multiancestry_tier_2")]
#             variant_id min_pval_iibdgc_multiancestry_tier_2
#                 <char>                                <num>
# 1: chr10:124750812:G:A                            4.410e-07
# 2:  chr10:26890667:T:G                            6.655e-02
# 3:  chr1:209797265:G:A                            1.873e-05
# 4:   chr20:6113242:T:G                            1.095e-04

## MULTIANCESTRY:

dl$label_ma<-""
dl$label_ma[which(dl$region_genomewide_significance_eur_eas_sas_tier_2=="N" & dl$gwas_signif_delange=="Y")]<-"Genome-wide significant in de Lange 2017"

dl$label_ma<-factor(dl$label_ma,levels=c("","Genome-wide significant in de Lange 2017"))
table(dl$label_ma)


p1<-ggplot(dl, aes(x=-log10(min_pval_delange),y=-log10(min_pval_iibdgc_multiancestry_tier_2),color=label_ma)) + 
  geom_point() + xlim(0,14) +
  ylim(0,8) +
  theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
 panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
  scale_color_manual(values=col_all) +
  geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=dl %>% filter(variant_id!=""),
  aes(label=variant_id),show.legend = FALSE)


## EUR:

dl$label_eur<-""
dl$label_eur[which(dl$region_genomewide_significance_eur_tier_2=="N" & dl$gwas_signif_delange=="Y")]<-"Genome-wide significant in de Lange 2017"

dl$label_eur<-factor(dl$label_eur,levels=c("","Genome-wide significant in de Lange 2017"))
table(dl$label_eur)
#                           
#                                                                                       17 
# No genome-wide significant in IIBDGC 2025 GWAS; Genome-wide significant in de Lange 2017 
#                                                                                        3 

p2<-ggplot(dl, aes(x=-log10(min_pval_delange),y=-log10(min_pval_iibdgc_eur_tier_2),color=label_eur)) + 
  geom_point() + xlim(0,14) +
  ylim(0,8) +
  theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
 panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
  scale_color_manual(values=col_all) +
  geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=dl %>% filter(variant_id!=""),
  aes(label=variant_id),show.legend = FALSE)
  
  
################################################
## VARIANTS WE CAN ASSESS WITH LIU:

liu<-all[which(all$'FM study'=="Liu"),]

liu$variant_id<-liu$MarkerName
liu$variant_id[which(liu$region_genomewide_significance_eur_eas_sas_tier_2=="N" & liu$gwas_signif_liu=="N")]<-""
liu$variant_id
#  [1] "chr10:101436625:AATAGATAGATAGATAG:AATAGATAGATAGATAGATAG"
#  [2] "chr11:17311651:A:G"                                     
#  [3] "chr11:36410895:C:T"                                     
#  [4] ""                                                       
#  [5] ""                                                       
#  [6] "chr12:69348988:AT:A"                                    
#  [7] "chr13:52062509:C:T"                                     
#  [8] "chr16:2117121:C:T"                                      
#  [9] ""                                                       
# [10] "chr18:2734945:C:G"                                      
# [11] "chr18:50977537:G:A"                                     
# [12] "chr19:36991870:G:A"                                     
# [13] "chr1:108112400:A:C"                                     
# [14] ""                                                       
# [15] "chr2:127092119:T:C"                                     
# [16] ""                                                       
# [17] "chr3:52013158:TGGA:T"                                   
# [18] "chr5:112513193:C:T"                                     
# [19] "chr5:80229545:A:G"                                      
# [20] ""                                                       
# [21] ""   

liu[which(liu$variant_id!=""),c("variant_id","min_pval_iibdgc_multiancestry_tier_2")]
#                                                  variant_id
#                                                      <char>
#  1: chr10:101436625:AATAGATAGATAGATAG:AATAGATAGATAGATAGATAG
#  2:                                      chr11:17311651:A:G
#  3:                                      chr11:36410895:C:T
#  4:                                     chr12:69348988:AT:A
#  5:                                      chr13:52062509:C:T
#  6:                                       chr16:2117121:C:T
#  7:                                       chr18:2734945:C:G
#  8:                                      chr18:50977537:G:A
#  9:                                      chr19:36991870:G:A
# 10:                                      chr1:108112400:A:C
# 11:                                      chr2:127092119:T:C
# 12:                                    chr3:52013158:TGGA:T
# 13:                                      chr5:112513193:C:T
# 14:                                       chr5:80229545:A:G
#     min_pval_iibdgc_multiancestry_tier_2
#                                    <num>
#  1:                                   NA
#  2:                            5.726e-04
#  3:                            5.198e-06
#  4:                            3.137e-05
#  5:                            4.815e-01
#  6:                            1.116e-01
#  7:                            1.157e-03
#  8:                            1.628e-07
#  9:                            1.074e-03
# 10:                            2.966e-04
# 11:                            3.838e-05
# 12:                            6.920e-07
# 13:                            1.220e-02
# 14:                            8.023e-06



## MULTIANCESTRY:

liu$label_ma<-""
liu$label_ma[which(liu$region_genomewide_significance_eur_eas_sas_tier_2=="N" & liu$gwas_signif_liu=="Y")]<-"Genome-wide significant in Liu 2023"

liu$label_ma<-factor(liu$label_ma,levels=c("","Genome-wide significant in Liu 2023"))
table(liu$label_ma)
#                          No genome-wide significant in IIBDGC 2025 GWAS or Liu 2023 
#                                                                                   7 
# No genome-wide significant in IIBDGC 2025 GWAS; Genome-wide significant in Liu 2023 
#                                                                                  14


p3<-ggplot(liu, aes(x=-log10(min_pval_liu),y=-log10(min_pval_iibdgc_multiancestry_tier_2),color=label_ma)) + 
  geom_point() + xlim(0,14) +
  ylim(0,8) +
  theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
 panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
  scale_color_manual(values=col_all) +
  geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=liu %>% filter(variant_id!=""),
  aes(label=variant_id),show.legend = FALSE)
  

## EUR:


liu$label_eur<-""
liu$label_eur[which(liu$region_genomewide_significance_eur_2=="N" & liu$gwas_signif_liu=="Y")]<-"Genome-wide significant in Liu 2023"

liu$label_eur<-factor(liu$label_ma,levels=c("","Genome-wide significant in Liu 2023"))
table(liu$label_eur)
#                           
#                                                                                   7 
# No genome-wide significant in IIBDGC 2025 GWAS; Genome-wide significant in Liu 2023 
#                                                                                  14


p4<-ggplot(liu, aes(x=-log10(min_pval_liu),y=-log10(min_pval_iibdgc_eur_tier_2),color=label_eur)) + 
  geom_point() + xlim(0,14) +
  ylim(0,8) +
  theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
 panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
  scale_color_manual(values=col_all) +
  geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=liu %>% filter(variant_id!=""),
  aes(label=variant_id),show.legend = FALSE)
  

p21<-ggarrange(p2,p1,nrow=2,common.legend=T,legend="bottom",labels = c("a", "c"))
p21<-annotate_figure(p21, top = text_grob("Association signals identified before Liu 2023"))
# p21


p43<-ggarrange(p4,p3,nrow=2,common.legend=T,legend="bottom",labels = c("b", "d"))
p43<-annotate_figure(p43, top = text_grob("Association signals identified in Liu 2023"))
# p43

p<-ggarrange(p21,p43,ncol=2)

ggsave(paste0("~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_non_replicated_signals.jpeg",sep=""),
    p,
    width = 15,
    height = 10,
    dpi = 300,
    units = c("in"),
    limitsize = T
)

q("no")