# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# singularity exec iibdgc_postprocess_10_singularity.sif

path_gwas="/path/to/ibdgwas/IIBDGC/"


MEM=2000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -G humgen-priority -q normal R 

library(data.table)
library(ggplot2)
library(qvalue)
library(dplyr)
library(tibble)

rm(list=ls())
path_gwas<-"/path/to/ibdgwas/IIBDGC/"


# load output files:

list_files<-list.files(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/genetic_correlation/lcv/"))
list_files<-list_files[!list_files %in% "plots"]
list_files[duplicated(list_files)]

rm(dat)
for (i in 1:length(list_files)) {

  tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/genetic_correlation/lcv/",list_files[i]))

  tmp$file<-NA
  tmp$file<-gsub("_lcv_output.tsv.gz","",list_files[i])

  if(i==1) {
    dat<-tmp
  }else{
    dat<-rbind(dat,tmp)
  }
  rm(tmp)
}

# 86*3
# [1] 261

length(list_files)
# [1] 261

dim(dat)
# [1] 261  13

# retain only pairs with significant rg:

list_signif_pairs<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/genetic_correlation/list_pairs_traits_significantly_correlated_with_ibd_cd_uc_with_qval",sep=""),head=T)
list_signif_pairs$pair_ibd_noibd<-paste(list_signif_pairs$p1,list_signif_pairs$p2,sep="_")

dim(list_signif_pairs)
# [1] 158  9

dat$pair_ibd_noibd<-paste(dat$d1,dat$d2,sep="_")
dim(dat[which(dat$pair_ibd_noibd %in% list_signif_pairs$pair_ibd_noibd),])
# [1] 158 14

dat<-dat[which(dat$pair_ibd_noibd %in% list_signif_pairs$pair_ibd_noibd),]

dat<-merge(dat,list_signif_pairs,by="pair_ibd_noibd")
colnames(dat)[22]<-"qval_rg"

# estimate significant resutls
dat$qval<-qvalue(dat$pval.gcpzero.2tailed)$qvalues

# set to NA those gcp with qval >0.05:
dat$gcp<-round(dat$gcp.pm,digits=1)

dim(dat[which(dat$qval>0.05),])
# 104
dat$gcp[which(dat$qval>0.05)]<-NA




# OUTPUT VARIABLES: 
#	lcv.output, a list with named entries:
#   "zscore", Z score for partial genetic causality. zscore>>0 implies gcp>0.
#	pval.gcpzero.2tailed, 2-tailed p-value for null that gcp=0.
#   "gcp.pm", posterior mean gcp (gcp=1: trait 1 -> trait 2; gcp=-1: trait 2-> trait 1); 
#   "gcp.pse", posterior standard error for gcp; 
#   "rho.est", estimated genetic correlation; 
#   "rho.err", standard error of rho estimate;
#   "pval.fullycausal" [2 entries], p-values for null that gcp=1 or that gcp=-1, respectively; 
#   "h2.zscore" [2 entries], z scores for trait 1 and trait 2 being heritable, respectively;
#   we recommend reporting results for h2.zscore > 7 (a very stringent threshold).



# set to NA traits with h2<7 - as recommended above
dim(dat[which(dat$h2.zscore1<7 | dat$h2.zscore2<7),])
# 9

names(table(dat$disease_trait[which(dat$h2.zscore1<7 | dat$h2.zscore2<7)]))
# [1] "Ankylosing spondylitis"                           
# [2] "Estimated degree of unsaturation"                 
# [3] "Primary sclerosing cholangitis"                   
# [4] "Ratio of omega-3 fatty acids to total fatty acids"

dat$gcp[which(dat$h2.zscore2<7)]<-NA

# report only pairs with GPC > 0.6:
dat[which(!is.na(dat$gcp) & abs(dat$gcp.pm)<0.6),c("p1","disease_trait")]
#         p1                                        disease_trait
#     <char>                                               <char>
#  1:     cd                                               Eczema
#  2:     cd                            C-reactive protein levels
#  3:     cd  Free cholesterol to total lipids ratio in large HDL
#  4:     cd                           Triglycerides in large HDL
#  5:     cd               Total cholesterol levels in medium HDL
#  6:     cd                     Cholesterol esters in medium HDL
#  7:     cd                       Free cholesterol in medium HDL
#  8:     cd                                 Phenylalanine levels
#  9:    ibd                         Systemic lupus erythematosus
# 10:    ibd                            Apolipoprotein A-I levels
# 11:    ibd                            Total cholesterol in HDL2
# 12:    ibd                      Total cholesterol levels in HDL
# 13:    ibd                       Mean diameter of HDL particles
# 14:    ibd                       Total cholesterol in large HDL
# 15:    ibd                      Cholesterol esters in large HDL
# 16:    ibd                        Free cholesterol in large HDL
# 17:    ibd                            Total lipids in large HDL
# 18:    ibd                 Concentration of large HDL particles
# 19:    ibd                           Phospholipids in large HDL
# 20:    ibd     Phospholipids to total lipids ratio in large HDL
# 21:    ibd Free cholesterol to total lipids ratio in medium HDL
# 22:    ibd                                 Phenylalanine levels
# 23:    ibd                      Phospholipids in very large HDL
#         p1                                        disease_trait


dat$gcp[which(abs(dat$gcp)<0.6)]<-NA

write.table(dat[,c("p1","p2","disease_trait","group","zscore","pval.gcpzero.2tailed","qval","gcp.pm","gcp.pse","rho.est","rho.err","pval.fullycausal1","pval.fullycausal2","h2.zscore1","h2.zscore2")],
"~/git/IIBDGC_GWAS/plots/paper_tables/Source_data_lcv.tsv",col.names=T,row.names=F,quote=F,sep="\t")



# keep only the significant ones to plot:
dat<-dat[which(dat$qval<0.05 & dat$h2.zscore2>7),]

dat$class<-NA
dat$class[grep("HDL",dat$disease_trait)]<-"HDL levels"
dat$class[grep("Apolipoprotein A-I levels",dat$disease_trait)]<-"Apolipoprotein A levels"
dat$class[grep("Phenylalanine levels",dat$disease_trait)]<-"Phenylalanine levels"
dat$class[grep("small VLDL|onounsaturated fatty acids",dat$disease_trait)]<-"Other Lipid levels"
dat$class[which(is.na(dat$gcp))]<-"|gcp| < 0.6"

table(dat$class)
# Apolipoprotein A levels              HDL levels      Other Lipid levels 
#                       1                      23                       3 
#    Phenylalanine levels             |gcp| < 0.6 
#                       1                      25


dat$class<-factor(dat$class,levels=c("Apolipoprotein A levels","HDL levels","Other Lipid levels","Phenylalanine levels","|gcp| < 0.6"))
cols<-c("#7eb0d5","#7c1158", "#4421af", "#1a53ff","grey")



# PLOT OPTION 1

# dat$label<-NA
# dat$label[which(dat$disease_trait %in% c("Apolipoprotein A-I levels","Phenylalanine levels"))]<-dat$disease_trait[which(dat$disease_trait %in% c("Apolipoprotein A-I levels","Phenylalanine levels"))]

# p<-ggplot(dat1, aes(x = rho.est, y = gcp.pm, color = class)) +
#     geom_point(shape = 16, size = 2) + xlim(-0.2,0.2) + ylim(-1,1) +
#     facet_grid(. ~ toupper(p1), scales = "free",space = "free") + 
#     ylab("Genetic Causal Effect\nnon-IBD trait on IBD <-> IBD trait on non-IBD\n(gcp)") + 
#     xlab("Genetic Correlation (rg)") +
#     theme_bw() + scale_colour_manual(values=cols) + theme(legend.title=element_blank()) 
#     # + 
#     # geom_label(aes(label=label))


# ggsave(
#   "~/git/IIBDGC_GWAS/plots/paper_figures/Figure_3_lcv_genetic_causal_effect.pdf",
#   p,
#   width = 150,
#   height = 50,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=2
# )

# # PLOT OPTION 1

# dat$rg_direction<-NA
# dat$rg_direction[which(dat$rg>0)]<-"+"
# dat$rg_direction[which(dat$rg<0)]<-"-"

# table(dat$rg_direction,useNA="ifany")
#  -  + 
# 45  8 

# palette_red_blue <- colorRampPalette(colors = c("#B2182B","grey", "#2166AC"))
# palette_red_blue(21)
# cols<-c("#B2182B","#B32839","#B43948","#B54957","#B65A65","grey","grey",
# "grey","grey","grey","grey","grey","grey","grey","grey","grey","#5F89B3","#5080B1","#4077AF","#306EAD","#2166AC")


# p<-ggplot(dat, aes(x = p1, y = disease_trait, fill = gcp)) +
#   geom_tile(color = "white",
#             lwd = 1.5,
#             linetype = 1) + theme_bw() +
#     facet_grid(group ~ toupper(p1) , scales = "free",space = "free") +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) +
#     scale_fill_stepsn(breaks=c(seq(-1,1,by=0.1)),colours=cols,limits=c(-1,1)) +
#   geom_text(aes(label = rg_direction), size = 12 / .pt)  +
#   theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())


# PLOT OPTION 2 - just export the data to later make a panel of it

panel_name<-"Figure_2b_lcv_genetic_causal_effect_heatmap"

dat$rg_direction<-NA
dat$rg_direction[which(dat$rg>0)]<-"+"
dat$rg_direction[which(dat$rg<0)]<-"-"
dat$rg_direction[which(dat$qval>0.05)]<-NA

dat$lcv_direction<-NA
dat$lcv_direction[which(dat$gcp>0)]<-"+"
dat$lcv_direction[which(dat$gcp<0)]<-"-"
dat$rg_direction[which(dat$qval>0.05)]<-NA


table(dat$lcv_direction,useNA="ifany")
  #  -    + <NA> 
  # 27    1   25 

list_traits_tokeep<-dat$disease_trait[which(dat$group %in% c("Metabolic biomarkers - Apolipoproteins","Metabolic biomarkers - Amino acids"))]
list_traits_tokeep<-c(list_traits_tokeep,"Total cholesterol in large HDL")

dat1<-dat[which(dat$disease_trait %in% list_traits_tokeep),]
dat1$disease_trait<-factor(dat1$disease_trait,levels=rev(c("Apolipoprotein A-I levels","Total cholesterol in large HDL","Phenylalanine levels")))
colnames(dat1)

dat1<-dat1[,c("disease_trait","p1","group","rg","qval","gcp","rg_direction","lcv_direction")]

# force keeping an uc category:
dat1_tmp<-dat1[which(dat1$p1=="ibd"),]
dat1_tmp[,4:ncol(dat1_tmp)]<-NA
dat1_tmp$p1<-"uc"

dat1<-rbind(dat1,dat1_tmp)

write.table(dat1,paste0("~/git/IIBDGC_GWAS/plots/paper_figures/",panel_name,".tsv"),col.names=T,row.names=F,sep="\t")

q("no")