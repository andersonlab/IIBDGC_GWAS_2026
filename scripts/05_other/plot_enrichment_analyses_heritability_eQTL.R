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
library(reshape2)


rm(list=ls())

path_gwas="/path/to/ibdgwas/IIBDGC/"

col<-carto_pal(12, "Safe")
col_ligth_1<-lighten(col, 0.90, space = "HCL")
col_ligth_2<-lighten(col, 0.70, space = "HCL")
col_ligth_3<-lighten(col, 0.50, space = "HCL")
col_ligth_4<-lighten(col, 0.30, space = "HCL")
col_ligth_5<-lighten(col, 0.10, space = "HCL")

colo<-c(col,col_ligth_1[0:11],col_ligth_2[0:11],col_ligth_3[0:11],col_ligth_4[0:11],col_ligth_5[0:11])
colo<-sample(colo)


#################################################################################################
# map to rename conditions - see /path/to/user/git/snakemake_colocalisation/scripts/create_maps_for_IIBDGC.R

map<-fread("/path/to/project",head=T)
map<-as.data.frame(map)
map<-map[which(map$quant_method %in% c("aptamer","ge","microarray")),]

########


files<-list.files(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",sep=""))
files<-files[grep(".results",files)]
# files<-files[grep("r2_0.25",files)]

files_2<-list.files(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/other_traits/",sep=""))
files_2<-files_2[grep(".results",files_2)]
# files_2<-files_2[grep("r2_0.25",files_2)]

files<-c(files,files_2)

length(files)
# [1] 4580

threshold<-0.05

###############################################################################################

# ###################
# # test r2 threshodl used in Macromap


# files_macromap<-files[grep("macromap",files)]

# for (i in 1:length(files_macromap)) {

#     file_tmp<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_macromap[i],sep="")

#     if(file.exists(file_tmp)) {
#       tmp<-fread(file_tmp,head=T)
#     } else {
#       tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/other_traits/",files_macromap[i],sep=""),head=T)
#     }    
    
#     tmp$pheno<-toupper(gsub(".baseline.files_version_april2025.results","",files_macromap[i]))

#     tmp<-tmp[which(tmp$Category=="L2_1",)]
#     tmp$Category<-files_macromap[i]
#     tmp$pheno<-toupper(gsub(".baseline.*","",tmp$Category))
#     tmp$Category<-(gsub(".*baseline.macromap_","",tmp$Category))
#     tmp$Category<-(gsub("_allchr_list_tier_2","",tmp$Category))
#     tmp$Category<-(gsub(".files_version_april2025.results","",tmp$Category))
#     tmp$ld<-(gsub(".*_ld_r2_","",tmp$Category))
#     tmp$Category<-(gsub("_ld_r2.*","",tmp$Category))

#     if(i==1) {
#         enrich<-tmp
#     } else {
#         enrich<-rbind(enrich,tmp)
#     }

# }

# dim(enrich)
# # [1] 264  11


# enrich<-enrich[which(enrich$pheno %in% c("IBD","CD","UC")),]
# enrich$Coefficient_zscore<-enrich$'Coefficient_z-score'

# enrich_resh_cd<-reshape(enrich[which(enrich$pheno=="CD"),c("Category","ld","Coefficient_zscore")], idvar = "Category", timevar = "ld", direction = "wide")
# enrich_resh_uc<-reshape(enrich[which(enrich$pheno=="UC"),c("Category","ld","Coefficient_zscore")], idvar = "Category", timevar = "ld", direction = "wide")
# enrich_resh_ibd<-reshape(enrich[which(enrich$pheno=="IBD"),c("Category","ld","Coefficient_zscore")], idvar = "Category", timevar = "ld", direction = "wide")


# ibd <- melt(round(cor(enrich_resh_ibd[,c(2:4)]),2))
# cd <- melt(round(cor(enrich_resh_cd[,c(2:4)]),2))
# uc <- melt(round(cor(enrich_resh_uc[,c(2:4)]),2))

# p1<-ggplot(data = cd, aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile() + geom_text(aes(label = value), color = "white", size = 3) + 
#   theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.title = element_blank()) +
#   scale_fill_gradient(high = "blue", low = "white", 
#    limit = c(0,1), space = "Lab", 
#     name="Pearson\nCorrelation") + ggtitle("CD")

# p2<-ggplot(data = ibd, aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile() + geom_text(aes(label = value), color = "white", size = 3) + 
#   theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.title = element_blank()) +
#   scale_fill_gradient(high = "blue", low = "white", 
#    limit = c(0,1), space = "Lab", 
#     name="Pearson\nCorrelation") + ggtitle("IBD")

# p3<-ggplot(data = uc, aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile() + geom_text(aes(label = value), color = "white", size = 3) + 
#   theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.title = element_blank()) +
#   scale_fill_gradient(high = "blue", low = "white", 
#    limit = c(0,1), space = "Lab", 
#     name="Pearson\nCorrelation") + ggtitle("UC")

# p_row2<-ggarrange(p1,p2,p3,ncol=3,common.legend=T,legend="right",labels=("c"))

# # compare pvalues:
# p1<-ggplot(enrich, aes(y=Category, x=-log10(Enrichment_p),color=Category)) +
#   geom_point() + scale_y_discrete(limits=rev) + 
#   scale_fill_manual(values=c(colo)) +
#   theme(axis.title.x = element_blank(),axis.title.y = element_blank()) + 
#   facet_grid(pheno ~ ld, space = "free") + theme(legend.position="right",legend.title=element_blank())


# p2<-ggplot(enrich, aes(y=Category, x=Coefficient_zscore,fill=Category)) +
#   geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) + 
#   scale_fill_manual(values=c(colo)) +
#   facet_grid(pheno ~ ld, scales = "free",space = "free")

# p_row1<-ggarrange(p2,p1,ncol=2,common.legend=T,legend="right",labels=c("a","b"))

# p<-ggarrange(p_row1,p_row2,nrow=2,heights=c(2,0.8))


# ggsave(
#   "~/git/IIBDGC_GWAS/plots/eqtl_enrichment/Supplementary_Figure_eqtl_enrichment_analyses_ld_thresholds_concordance.png",
#   p,
#   width = 240,
#   height = 240,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=2
# )


###############################################################################################

files<-files[grep("r2_0.25",files)]

###############################################################################################

###################
# Macromap

files_macromap<-files[grep("macromap",files)]


for (i in 1:length(files_macromap)) {

    file_tmp<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_macromap[i],sep="")

    if(file.exists(file_tmp)) {
      tmp<-fread(file_tmp,head=T)
    } else {
      tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/other_traits/",files_macromap[i],sep=""),head=T)
    }    
    
    tmp$pheno<-toupper(gsub(".baseline.files_version_april2025.results","",files_macromap[i]))

    tmp<-tmp[which(tmp$Category=="L2_1",)]
    tmp$Category<-files_macromap[i]
    tmp$pheno<-toupper(gsub(".baseline.*","",tmp$Category))
    tmp$Category<-(gsub(".*baseline.macromap_","",tmp$Category))
    tmp$Category<-(gsub("_allchr_list_tier_2","",tmp$Category))
    tmp$Category<-(gsub(".files_version_april2025.results","",tmp$Category))
    tmp$Category<-(gsub("_ld_r2.*","",tmp$Category))

    if(i==1) {
        enrich<-tmp
    } else {
        enrich<-rbind(enrich,tmp)
    }

}

dim(enrich)
# [1] 96  11

enrich$Category<-as.factor(enrich$Category)
enrich$Coefficient_zscore<-enrich$'Coefficient_z-score'
enrich$pheno<-factor(enrich$pheno,levels=c("CD","IBD","UC","HEIGHT"))

enrich<-merge(enrich,map,by.x="Category",by.y="id_map")

# # only colour the significant results:
enrich$signif<-enrich$tissue_cell_condition_label
enrich$signif[which(enrich$Enrichment_p>threshold)]<-NA


# To compare cell-type groups, we use the coefficient z-score in an analysis containing the full baseline model - source:https://github.com/bulik/ldsc/wiki/Partitioned-Heritability

xmin<-min(0,enrich$Coefficient_zscore-enrich$Coefficient_std_error)
xmax<-max(enrich$Coefficient_zscore+enrich$Coefficient_std_error)

p1<-
ggplot(enrich, aes(y=tissue_cell_condition_label, x=Coefficient_zscore,fill=signif)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) + 
  scale_fill_manual(values=c(colo)) +
  facet_grid( ~ pheno, scales = "free",space = "free") +
  xlim(xmin,xmax) + theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.title = element_blank(),legend.position="none",
   plot.margin = unit(c(1,1,1,1), "cm")) 


cd<-enrich[which(enrich$pheno=="CD"),c("tissue_cell_condition_label","Coefficient_zscore")]
colnames(cd)[2]<-paste(colnames(cd)[2],"cd",sep="_")
uc<-enrich[which(enrich$pheno=="UC"),c("tissue_cell_condition_label","Coefficient_zscore")]
colnames(uc)[2]<-paste(colnames(uc)[2],"uc",sep="_")
height<-enrich[which(enrich$pheno=="HEIGHT"),c("tissue_cell_condition_label","Coefficient_zscore")]
colnames(height)[2]<-paste(colnames(height)[2],"height",sep="_")

comp<-merge(cd,uc,by="tissue_cell_condition_label")
comp<-merge(comp,height,by="tissue_cell_condition_label")

rm(cd,uc,height)

minlim<-min(comp$Coefficient_zscore_cd,comp$Coefficient_zscore_uc)
maxlim<-max(comp$Coefficient_zscore_cd,comp$Coefficient_zscore_uc)

p2<-ggplot(comp, aes(y=Coefficient_zscore_cd, x=Coefficient_zscore_height,color=tissue_cell_condition_label)) + 
  geom_point(size=2) +  geom_smooth(method='lm') +
  scale_colour_manual(values=c(colo)) + 
  theme_bw() + theme(panel.background = element_blank(),legend.title = element_blank(),
                 axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(1,1,1,1), "cm")) + 
                 labs(x="Height Coefficient Z-score",y="CD Coefficient Z-score") + 
                 xlim(minlim-1,maxlim+1) + ylim(minlim-1,maxlim+1) + geom_abline(color="grey",linetype = "dashed")

p3<-ggplot(comp, aes(y=Coefficient_zscore_uc, x=Coefficient_zscore_height,color=tissue_cell_condition_label)) + 
  geom_point(size=2) +  geom_smooth(method='lm') +
  scale_colour_manual(values=c(colo)) + 
  theme_bw() + theme(panel.background = element_blank(),legend.title = element_blank(),
                 axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(1,1,1,1), "cm")) + 
                 labs(x="Height Coefficient Z-score",y="UC Coefficient Z-score") + 
                 xlim(minlim-1,maxlim+1) + ylim(minlim-1,maxlim+1) + geom_abline(color="grey",linetype = "dashed")


p23<-ggarrange(p2,p3,nrow=2,labels=c("b","c"))

p123<-ggarrange(p1,p23,ncol=2,widths=c(4,1.5),labels=c("a"))


ggsave(
  "~/git/IIBDGC_GWAS/plots/eqtl_enrichment/Supplementary_Figure_eqtl_enrichment_analyses_macromap.png",
  p123,
  width = 120,
  height = ceiling(length(levels(as.factor(enrich$Category)))*2.8),
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)

enrich_macromap<-enrich
dim(enrich_macromap)
# [1] 96 18


############################################################################
# edQTL GTEx

files_ed<-files[grep("edQTL_GTEX",files)]


for (i in 1:length(files_ed)) {

    file_tmp<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_ed[i],sep="")

    if(file.exists(file_tmp)) {
      tmp<-fread(file_tmp,head=T)
    } else {
      tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/other_traits/",files_ed[i],sep=""),head=T)
    }    
  

    tmp<-tmp[which(tmp$Category=="L2_1",)]
    tmp$Category<-files_ed[i]
    tmp$Category<-gsub("-","_",tmp$Category)
    tmp$pheno<-toupper(gsub(".baseline.*","",tmp$Category))
    tmp$Category<-(gsub(".*baseline.edQTL_GTEX_","",tmp$Category))
    tmp$Category<-(gsub("_allchr_list_tier_2","",tmp$Category))
    tmp$Category<-(gsub(".files_version_april2025.results","",tmp$Category))
    tmp$Category<-(gsub("_ld_r2.*","",tmp$Category))

    if(i==1) {
        enrich<-tmp
    } else {
        enrich<-rbind(enrich,tmp)
    }

}

dim(enrich)
# [1] 40  11

enrich$Category<-as.factor(enrich$Category)
enrich$Coefficient_zscore<-enrich$'Coefficient_z-score'
enrich$pheno<-factor(enrich$pheno,levels=c("CD","IBD","UC","HEIGHT"))

enrich<-merge(enrich,map,by.x="Category",by.y="id_map")

# # only colour the significant results:
enrich$signif<-enrich$tissue_cell_condition_label
enrich$signif[which(enrich$Enrichment_p>threshold)]<-NA


# To compare cell-type groups, we use the coefficient z-score in an analysis containing the full baseline model - source:https://github.com/bulik/ldsc/wiki/Partitioned-Heritability

xmin<-min(0,enrich$Coefficient_zscore-enrich$Coefficient_std_error)
xmax<-max(enrich$Coefficient_zscore+enrich$Coefficient_std_error)

p1<-
ggplot(enrich, aes(y=tissue_cell_condition_label, x=Coefficient_zscore,fill=signif)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) + 
  scale_fill_manual(values=c(colo)) +
  facet_grid( ~ pheno, scales = "free",space = "free") +
  xlim(xmin,xmax) + theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.title = element_blank(),legend.position="none",
   plot.margin = unit(c(1,1,1,1), "cm")) 


cd<-enrich[which(enrich$pheno=="CD"),c("tissue_cell_condition_label","Coefficient_zscore")]
colnames(cd)[2]<-paste(colnames(cd)[2],"cd",sep="_")
uc<-enrich[which(enrich$pheno=="UC"),c("tissue_cell_condition_label","Coefficient_zscore")]
colnames(uc)[2]<-paste(colnames(uc)[2],"uc",sep="_")
height<-enrich[which(enrich$pheno=="HEIGHT"),c("tissue_cell_condition_label","Coefficient_zscore")]
colnames(height)[2]<-paste(colnames(height)[2],"height",sep="_")

comp<-merge(cd,uc,by="tissue_cell_condition_label")
comp<-merge(comp,height,by="tissue_cell_condition_label")

rm(cd,uc,height)

minlim<-min(comp$Coefficient_zscore_cd,comp$Coefficient_zscore_uc)
maxlim<-max(comp$Coefficient_zscore_cd,comp$Coefficient_zscore_uc)

p2<-ggplot(comp, aes(y=Coefficient_zscore_cd, x=Coefficient_zscore_height,color=tissue_cell_condition_label)) + 
  geom_point(size=2) +  geom_smooth(method='lm') +
  scale_colour_manual(values=c(colo)) + 
  theme_bw() + theme(panel.background = element_blank(),legend.title = element_blank(),
                 axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(1,1,1,1), "cm")) + 
                 labs(x="Height Coefficient Z-score",y="CD Coefficient Z-score") + 
                 xlim(minlim-1,maxlim+1) + ylim(minlim-1,maxlim+1) + geom_abline(color="grey",linetype = "dashed")

p3<-ggplot(comp, aes(y=Coefficient_zscore_uc, x=Coefficient_zscore_height,color=tissue_cell_condition_label)) + 
  geom_point(size=2) +  geom_smooth(method='lm') +
  scale_colour_manual(values=c(colo)) + 
  theme_bw() + theme(panel.background = element_blank(),legend.title = element_blank(),
                 axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(1,1,1,1), "cm")) + 
                 labs(x="Height Coefficient Z-score",y="UC Coefficient Z-score") + 
                 xlim(minlim-1,maxlim+1) + ylim(minlim-1,maxlim+1) + geom_abline(color="grey",linetype = "dashed")



p23<-ggarrange(p2,p3,nrow=2,labels=c("b","c"))
p123<-ggarrange(p1,p23,ncol=2,widths=c(4,1.5),labels=c("a"))

hg<-ceiling(length(levels(as.factor(enrich$Category)))*2.8)

if (hg<50) {
  p23<-ggarrange(p2,p3,ncol=2,labels=c("b","c"))
  p123<-ggarrange(p1,p23,nrow=2,labels=c("a"))
}

ggsave(
  "~/git/IIBDGC_GWAS/plots/eqtl_enrichment/Supplementary_Figure_eqtl_enrichment_analyses_edQTL_GTEX.png",
  p123,
  width = 120,
  height = hg+40,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)


enrich_edQTL_GTEX<-enrich
dim(enrich_edQTL_GTEX)
# [[1]  8 18


############################################################################
# Hu 2021

files_hu<-files[grep("hu_2021",files)]


for (i in 1:length(files_ed)) {

    file_tmp<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_hu[i],sep="")

    if(file.exists(file_tmp)) {
      tmp<-fread(file_tmp,head=T)
    } else {
      tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/other_traits/",files_hu[i],sep=""),head=T)
    }    
  

    tmp<-tmp[which(tmp$Category=="L2_1",)]
    tmp$Category<-files_hu[i]
    tmp$Category<-gsub("-","_",tmp$Category)
    tmp$pheno<-toupper(gsub(".baseline.*","",tmp$Category))
    tmp$Category<-(gsub(".*baseline.hu_2021_","",tmp$Category))
    tmp$Category<-(gsub("_allchr_list_tier_2","",tmp$Category))
    tmp$Category<-(gsub(".files_version_april2025.results","",tmp$Category))
    tmp$Category<-(gsub("_ld_r2.*","",tmp$Category))

    if(i==1) {
        enrich<-tmp
    } else {
        enrich<-rbind(enrich,tmp)
    }

}

dim(enrich)
# [1] 96  11

enrich$Category<-as.factor(enrich$Category)
enrich$Coefficient_zscore<-enrich$'Coefficient_z-score'
enrich$pheno<-factor(enrich$pheno,levels=c("CD","IBD","UC","HEIGHT"))

enrich<-merge(enrich,map,by.x="Category",by.y="id_map")

# # only colour the significant results:
enrich$signif<-enrich$tissue_cell_condition_label
enrich$signif[which(enrich$Enrichment_p>threshold)]<-NA



# To compare cell-type groups, we use the coefficient z-score in an analysis containing the full baseline model - source:https://github.com/bulik/ldsc/wiki/Partitioned-Heritability

xmin<-min(0,enrich$Coefficient_zscore-enrich$Coefficient_std_error)
xmax<-max(enrich$Coefficient_zscore+enrich$Coefficient_std_error)

p1<-
ggplot(enrich, aes(y=tissue_cell_condition_label, x=Coefficient_zscore,fill=tissue_cell_condition_label)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) + 
  scale_fill_manual(values=c(colo)) +
  facet_grid( ~ pheno, scales = "free",space = "free") +
  xlim(xmin,xmax) + theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.title = element_blank(),legend.position="none",
   plot.margin = unit(c(1,1,1,1), "cm")) 


cd<-enrich[which(enrich$pheno=="CD"),c("tissue_cell_condition_label","Coefficient_zscore")]
colnames(cd)[2]<-paste(colnames(cd)[2],"cd",sep="_")
uc<-enrich[which(enrich$pheno=="UC"),c("tissue_cell_condition_label","Coefficient_zscore")]
colnames(uc)[2]<-paste(colnames(uc)[2],"uc",sep="_")
height<-enrich[which(enrich$pheno=="HEIGHT"),c("tissue_cell_condition_label","Coefficient_zscore")]
colnames(height)[2]<-paste(colnames(height)[2],"height",sep="_")

comp<-merge(cd,uc,by="tissue_cell_condition_label")
comp<-merge(comp,height,by="tissue_cell_condition_label")

rm(cd,uc,height)

minlim<-min(comp$Coefficient_zscore_cd,comp$Coefficient_zscore_uc)
maxlim<-max(comp$Coefficient_zscore_cd,comp$Coefficient_zscore_uc)

p2<-ggplot(comp, aes(y=Coefficient_zscore_cd, x=Coefficient_zscore_height,color=tissue_cell_condition_label)) + 
  geom_point(size=2) +  geom_smooth(method='lm') +
  scale_colour_manual(values=c(colo)) + 
  theme_bw() + theme(panel.background = element_blank(),legend.title = element_blank(),
                 axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(1,1,1,1), "cm")) + 
                 labs(x="Height Coefficient Z-score",y="CD Coefficient Z-score") + 
                 xlim(minlim-1,maxlim+1) + ylim(minlim-1,maxlim+1) + geom_abline(color="grey",linetype = "dashed")

p3<-ggplot(comp, aes(y=Coefficient_zscore_uc, x=Coefficient_zscore_height,color=tissue_cell_condition_label)) + 
  geom_point(size=2) +  geom_smooth(method='lm') +
  scale_colour_manual(values=c(colo)) + 
  theme_bw() + theme(panel.background = element_blank(),legend.title = element_blank(),
                 axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(1,1,1,1), "cm")) + 
                 labs(x="Height Coefficient Z-score",y="UC Coefficient Z-score") + 
                 xlim(minlim-1,maxlim+1) + ylim(minlim-1,maxlim+1) + geom_abline(color="grey",linetype = "dashed")


p23<-ggarrange(p2,p3,nrow=2,labels=c("b","c"))
p123<-ggarrange(p1,p23,ncol=2,widths=c(4,1.5),labels=c("a"))

hg<-ceiling(length(levels(as.factor(enrich$Category)))*2.8)

if (hg<50) {
  p23<-ggarrange(p2,p3,ncol=2,labels=c("b","c"))
  p123<-ggarrange(p1,p23,nrow=2,labels=c("a"))
}

ggsave(
  "~/git/IIBDGC_GWAS/plots/eqtl_enrichment/Supplementary_Figure_eqtl_enrichment_analyses_hu_2021.png",
  p123,
  width = 120,
  height = hg+40,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)

enrich_hu_2021<-enrich
dim(enrich_hu_2021)
# [1] 4 18


############################################################################
# eQTL_catalogue

files_cat<-files[grep("eQTL_catalogue",files)]
length(files_cat)
# 896

for (i in 1:length(files_cat)) {

    file_tmp<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_cat[i],sep="")
    if(file.exists(file_tmp)) {
      tmp<-fread(file_tmp,head=T)
    } else {
      tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/other_traits/",files_cat[i],sep=""),head=T)
    }

    tmp<-tmp[which(tmp$Category=="L2_1",)]
    tmp$Category<-files_cat[i]
    tmp$pheno<-toupper(gsub(".baseline.*","",tmp$Category))

    tmp$Category<-(gsub(".*baseline.","",tmp$Category))
    tmp$Category<-(gsub("_allchr_list_tier_2","",tmp$Category))
    tmp$Category<-(gsub(".files_version_april2025.results","",tmp$Category))
    tmp$Category<-(gsub("_ld_r2.*","",tmp$Category))
    tmp$Category<-(gsub("eQTL_catalogue_","",tmp$Category))

    if(i==1) {
        enrich<-tmp
    } else {
        enrich<-rbind(enrich,tmp)
    }

}

dim(enrich)
# [1] 896  11

enrich$Category<-as.factor(enrich$Category)
enrich$Coefficient_zscore<-enrich$'Coefficient_z-score'
enrich$pheno<-factor(enrich$pheno,levels=c("CD","IBD","UC","HEIGHT"))

enrich<-merge(enrich,map,by.x="Category",by.y="id_map")

# # only colour the significant results:
enrich$signif<-enrich$tissue_cell_condition_label
enrich$signif[which(enrich$Enrichment_p>threshold)]<-NA


# To compare cell-type groups, we use the coefficient z-score in an analysis containing the full baseline model - source:https://github.com/bulik/ldsc/wiki/Partitioned-Heritability

studies<-names(table(enrich$study_label))

for (i in 1:length(studies)) {

  tmp<-enrich[which(enrich$study_label==studies[i]),]

  xmin<-min(0,tmp$Coefficient_zscore-tmp$Coefficient_std_error)
  xmax<-max(tmp$Coefficient_zscore+tmp$Coefficient_std_error)

  p1<-
  ggplot(tmp, aes(y=tissue_cell_condition_label, x=Coefficient_zscore,fill=signif)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) + 
    scale_fill_manual(values=c(colo)) +
    facet_grid( ~ pheno, scales = "free",space = "free") +
    xlim(xmin,xmax) + theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.title = element_blank(),legend.position="none",
    plot.margin = unit(c(1,1,1,1), "cm")) 


  cd<-tmp[which(tmp$pheno=="CD"),c("tissue_cell_condition_label","Coefficient_zscore")]
  colnames(cd)[2]<-paste(colnames(cd)[2],"cd",sep="_")
  uc<-tmp[which(tmp$pheno=="UC"),c("tissue_cell_condition_label","Coefficient_zscore")]
  colnames(uc)[2]<-paste(colnames(uc)[2],"uc",sep="_")
  height<-tmp[which(tmp$pheno=="HEIGHT"),c("tissue_cell_condition_label","Coefficient_zscore")]
  colnames(height)[2]<-paste(colnames(height)[2],"height",sep="_")

  comp<-merge(cd,uc,by="tissue_cell_condition_label")
  comp<-merge(comp,height,by="tissue_cell_condition_label")

  rm(cd,uc,height)

  minlim<-floor(min(comp$Coefficient_zscore_cd,comp$Coefficient_zscore_uc,comp$Coefficient_zscore_height))
  maxlim<-ceiling(max(comp$Coefficient_zscore_cd,comp$Coefficient_zscore_uc,comp$Coefficient_zscore_height))

  p2<-ggplot(comp, aes(y=Coefficient_zscore_cd, x=Coefficient_zscore_height,color=tissue_cell_condition_label)) + 
    geom_point(size=2) +  geom_smooth(method='lm') +
    scale_colour_manual(values=c(colo)) + 
    theme_bw() + theme(panel.background = element_blank(),legend.title = element_blank(),
                  axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(1,1,1,1), "cm")) + 
                  labs(x="Height Coefficient Z-score",y="CD Coefficient Z-score") + 
                  xlim(minlim-1,maxlim+1) + ylim(minlim-1,maxlim+1) + geom_abline(color="grey",linetype = "dashed")

  p3<-ggplot(comp, aes(y=Coefficient_zscore_uc, x=Coefficient_zscore_height,color=tissue_cell_condition_label)) + 
    geom_point(size=2) +  geom_smooth(method='lm') +
    scale_colour_manual(values=c(colo)) + 
    theme_bw() + theme(panel.background = element_blank(),legend.title = element_blank(),
                  axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(1,1,1,1), "cm")) + 
                  labs(x="Height Coefficient Z-score",y="UC Coefficient Z-score") + 
                  xlim(minlim-1,maxlim+1) + ylim(minlim-1,maxlim+1) + geom_abline(color="grey",linetype = "dashed")


  p23<-ggarrange(p2,p3,nrow=2,labels=c("b","c"))
  p123<-ggarrange(p1,p23,ncol=2,widths=c(4,1.5),labels=c("a"))

  hg<-ceiling(length(levels(as.factor(tmp$tissue_cell_condition_label)))*2.8)

  if (hg<50) {
    p23<-ggarrange(p2,p3,ncol=2,labels=c("b","c"))
    p123<-ggarrange(p1,p23,nrow=2,labels=c("a"))
  }

  ggsave(
    paste0("~/git/IIBDGC_GWAS/plots/eqtl_enrichment/Supplementary_Figure_eqtl_enrichment_analyses_",studies[i],".png"),
    p123,
    width = 120,
    height = hg+40,
    dpi = 400,
    units = c("mm"),
    limitsize = T,scale=2
  )


}



enrich_gwas_catalogue<-enrich
dim(enrich_gwas_catalogue)
# [1] 860  18

############################################################################
# IBDverse

files_ibdverse<-files[grep("IBDverse",files)]
length(files_ibdverse)
# 1008

for (i in 1:length(files_ibdverse)) {

    file_tmp<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_ibdverse[i],sep="")
    if(file.exists(file_tmp)) {
      tmp<-fread(file_tmp,head=T)
    } else {
      tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/other_traits/",files_ibdverse[i],sep=""),head=T)
    }

    tmp<-tmp[which(tmp$Category=="L2_1",)]
    tmp$Category<-files_ibdverse[i]
    tmp$pheno<-toupper(gsub(".baseline.*","",tmp$Category))

    tmp$Category<-(gsub(".*baseline.","",tmp$Category))
    tmp$Category<-(gsub("_allchr_list_tier_2","",tmp$Category))
    tmp$Category<-(gsub(".files_version_april2025.results","",tmp$Category))
    tmp$Category<-(gsub("_ld_r2.*","",tmp$Category))
    tmp$Category<-(gsub("IBDverse_","",tmp$Category))

    if(i==1) {
        enrich<-tmp
    } else {
        enrich<-rbind(enrich,tmp)
    }

}

dim(enrich)
# [1] 1008  11


enrich$Category<-as.factor(enrich$Category)
enrich$Coefficient_zscore<-enrich$'Coefficient_z-score'
enrich$pheno<-factor(enrich$pheno,levels=c("CD","IBD","UC","HEIGHT"))

enrich<-merge(enrich,map,by.x="Category",by.y="id_map")

# # only colour the significant results:
enrich$signif<-enrich$tissue_cell_condition_label
enrich$signif[which(enrich$Enrichment_p>threshold)]<-NA


# To compare cell-type groups, we use the coefficient z-score in an analysis containing the full baseline model - source:https://github.com/bulik/ldsc/wiki/Partitioned-Heritability

studies<-names(table(enrich$cell_label))

for (i in 1:length(studies)) {

  tmp<-enrich[which(enrich$cell_label==studies[i]),]

  xmin<-min(0,tmp$Coefficient_zscore-tmp$Coefficient_std_error)
  xmax<-max(tmp$Coefficient_zscore+tmp$Coefficient_std_error)

  p1<-
  ggplot(tmp, aes(y=tissue_cell_condition_label, x=Coefficient_zscore,fill=signif)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) + 
    scale_fill_manual(values=c(colo)) +
    facet_grid(tissue_label ~ pheno, scales = "free",space = "free") +
    xlim(xmin,xmax) + theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.title = element_blank(),legend.position="none",
    plot.margin = unit(c(1,1,1,1), "cm")) 


  cd<-tmp[which(tmp$pheno=="CD"),c("Category","tissue_cell_condition_label","Coefficient_zscore")]
  colnames(cd)[3]<-paste(colnames(cd)[3],"cd",sep="_")
  uc<-tmp[which(tmp$pheno=="UC"),c("Category","Coefficient_zscore")]
  colnames(uc)[2]<-paste(colnames(uc)[2],"uc",sep="_")
  height<-tmp[which(tmp$pheno=="HEIGHT"),c("Category","Coefficient_zscore")]
  colnames(height)[2]<-paste(colnames(height)[2],"height",sep="_")

  comp<-merge(cd,uc,by="Category")
  comp<-merge(comp,height,by="Category")

  rm(cd,uc,height)

  minlim<-floor(min(comp$Coefficient_zscore_cd,comp$Coefficient_zscore_uc,comp$Coefficient_zscore_height))
  maxlim<-ceiling(max(comp$Coefficient_zscore_cd,comp$Coefficient_zscore_uc,comp$Coefficient_zscore_height))

  p2<-ggplot(comp, aes(y=Coefficient_zscore_cd, x=Coefficient_zscore_height,color=tissue_cell_condition_label)) + 
    geom_point(size=2) +  
    scale_colour_manual(values=c(colo)) + 
    theme_bw() + theme(panel.background = element_blank(),legend.title = element_blank(),
                  axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(1,1,1,1), "cm")) + 
                  labs(x="Height Coefficient Z-score",y="CD Coefficient Z-score") + 
                  xlim(minlim-1,maxlim+1) + ylim(minlim-1,maxlim+1) + geom_abline(color="grey",linetype = "dashed")

  p3<-ggplot(comp, aes(y=Coefficient_zscore_uc, x=Coefficient_zscore_height,color=tissue_cell_condition_label)) + 
    geom_point(size=2) +  
    scale_colour_manual(values=c(colo)) + 
    theme_bw() + theme(panel.background = element_blank(),legend.title = element_blank(),
                  axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(1,1,1,1), "cm")) + 
                  labs(x="Height Coefficient Z-score",y="UC Coefficient Z-score") + 
                  xlim(minlim-1,maxlim+1) + ylim(minlim-1,maxlim+1) + geom_abline(color="grey",linetype = "dashed")


  p23<-ggarrange(p2,p3,nrow=2,labels=c("b","c"))
  p123<-ggarrange(p1,p23,ncol=2,widths=c(4,1.5),labels=c("a"))

  hg<-ceiling(length(levels(as.factor(as.character(tmp$tissue_cell_condition_label))))*3)

  if (hg<50) {
    p23<-ggarrange(p2,p3,ncol=2,labels=c("b","c"))
    p123<-ggarrange(p1,p23,nrow=2,labels=c("a"),heights=c(1.2,0.8))
    hg<-hg+60
  } 

  ggsave(
    paste0("~/git/IIBDGC_GWAS/plots/eqtl_enrichment/Supplementary_Figure_eqtl_enrichment_analyses_IBDverse_",studies[i],".png"),
    p123,
    width = 120,
    height = hg,
    dpi = 400,
    units = c("mm"),
    limitsize = T,scale=2
  )


}


enrich_IBDverse<-enrich
dim(enrich_IBDverse)
# [1] 1008   19



############################################################################
# Sparc

files_sparc<-files[grep("sparc",files)]


for (i in 1:length(files_sparc)) {

    file_tmp<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_sparc[i],sep="")

    if(file.exists(file_tmp)) {
      tmp<-fread(file_tmp,head=T)
    } else {
      tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/other_traits/",files_sparc[i],sep=""),head=T)
    }    
  

    tmp<-tmp[which(tmp$Category=="L2_1",)]
    tmp$Category<-files_sparc[i]
    tmp$Category<-gsub("-","_",tmp$Category)
    tmp$pheno<-toupper(gsub(".baseline.*","",tmp$Category))
    tmp$Category<-(gsub(".*baseline.pQTL_sparc_","",tmp$Category))
    tmp$Category<-(gsub("_allchr_list_tier_2","",tmp$Category))
    tmp$Category<-(gsub(".files_version_april2025.results","",tmp$Category))
    tmp$Category<-(gsub("_ld_r2.*","",tmp$Category))

    if(i==1) {
        enrich<-tmp
    } else {
        enrich<-rbind(enrich,tmp)
    }

}

dim(enrich)
# [1] 4  11

enrich$Category<-as.factor(enrich$Category)
enrich$Coefficient_zscore<-enrich$'Coefficient_z-score'
enrich$pheno<-factor(enrich$pheno,levels=c("CD","IBD","UC","HEIGHT"))

enrich<-merge(enrich,map,by.x="Category",by.y="id_map")

# # only colour the significant results:
enrich$signif<-enrich$tissue_cell_condition_label
enrich$signif[which(enrich$Enrichment_p>threshold)]<-NA



# To compare cell-type groups, we use the coefficient z-score in an analysis containing the full baseline model - source:https://github.com/bulik/ldsc/wiki/Partitioned-Heritability

xmin<-min(0,enrich$Coefficient_zscore-enrich$Coefficient_std_error)
xmax<-max(enrich$Coefficient_zscore+enrich$Coefficient_std_error)

p1<-
ggplot(enrich, aes(y=tissue_cell_condition_label, x=Coefficient_zscore,fill=tissue_cell_condition_label)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) + 
  scale_fill_manual(values=c(colo)) +
  facet_grid( ~ pheno, scales = "free",space = "free") +
  xlim(xmin,xmax) + theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.title = element_blank(),legend.position="none",
   plot.margin = unit(c(1,1,1,1), "cm")) 


cd<-enrich[which(enrich$pheno=="CD"),c("tissue_cell_condition_label","Coefficient_zscore")]
colnames(cd)[2]<-paste(colnames(cd)[2],"cd",sep="_")
uc<-enrich[which(enrich$pheno=="UC"),c("tissue_cell_condition_label","Coefficient_zscore")]
colnames(uc)[2]<-paste(colnames(uc)[2],"uc",sep="_")
height<-enrich[which(enrich$pheno=="HEIGHT"),c("tissue_cell_condition_label","Coefficient_zscore")]
colnames(height)[2]<-paste(colnames(height)[2],"height",sep="_")

comp<-merge(cd,uc,by="tissue_cell_condition_label")
comp<-merge(comp,height,by="tissue_cell_condition_label")

rm(cd,uc,height)

minlim<-min(comp$Coefficient_zscore_cd,comp$Coefficient_zscore_uc,comp$Coefficient_zscore_height)
maxlim<-max(comp$Coefficient_zscore_cd,comp$Coefficient_zscore_uc,comp$Coefficient_zscore_height)

p2<-ggplot(comp, aes(y=Coefficient_zscore_cd, x=Coefficient_zscore_height,color=tissue_cell_condition_label)) + 
  geom_point(size=2) +  geom_smooth(method='lm') +
  scale_colour_manual(values=c(colo)) + 
  theme_bw() + theme(panel.background = element_blank(),legend.title = element_blank(),
                 axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(1,1,1,1), "cm")) + 
                 labs(x="Height Coefficient Z-score",y="CD Coefficient Z-score") + 
                 xlim(minlim-1,maxlim+1) + ylim(minlim-1,maxlim+1) + geom_abline(color="grey",linetype = "dashed")

p3<-ggplot(comp, aes(y=Coefficient_zscore_uc, x=Coefficient_zscore_height,color=tissue_cell_condition_label)) + 
  geom_point(size=2) +  geom_smooth(method='lm') +
  scale_colour_manual(values=c(colo)) + 
  theme_bw() + theme(panel.background = element_blank(),legend.title = element_blank(),
                 axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(1,1,1,1), "cm")) + 
                 labs(x="Height Coefficient Z-score",y="UC Coefficient Z-score") + 
                 xlim(minlim-1,maxlim+1) + ylim(minlim-1,maxlim+1) + geom_abline(color="grey",linetype = "dashed")


p23<-ggarrange(p2,p3,nrow=2,labels=c("b","c"))
p123<-ggarrange(p1,p23,ncol=2,widths=c(4,1.5),labels=c("a"))

hg<-ceiling(length(levels(as.factor(enrich$Category)))*2.8)

if (hg<50) {
  p23<-ggarrange(p2,p3,ncol=2,labels=c("b","c"))
  p123<-ggarrange(p1,p23,nrow=2,labels=c("a"))
}

ggsave(
  "~/git/IIBDGC_GWAS/plots/eqtl_enrichment/Supplementary_Figure_eqtl_enrichment_analyses_sparc.png",
  p123,
  width = 120,
  height = hg+40,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)

enrich_sparc<-enrich
dim(enrich_sparc)
# [1] 4 20


############################################################################
# DeCODE

files_decode<-files[grep("decode",files)]


for (i in 1:length(files_decode)) {

    file_tmp<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files_decode[i],sep="")

    if(file.exists(file_tmp)) {
      tmp<-fread(file_tmp,head=T)
    } else {
      tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/other_traits/",files_decode[i],sep=""),head=T)
    }    
  

    tmp<-tmp[which(tmp$Category=="L2_1",)]
    tmp$Category<-files_decode[i]
    tmp$Category<-gsub("-","_",tmp$Category)
    tmp$pheno<-toupper(gsub(".baseline.*","",tmp$Category))
    tmp$Category<-(gsub(".*baseline.pQTL_decode_","",tmp$Category))
    tmp$Category<-(gsub("_allchr_list_tier_2","",tmp$Category))
    tmp$Category<-(gsub(".files_version_april2025.results","",tmp$Category))
    tmp$Category<-(gsub("_ld_r2.*","",tmp$Category))

    if(i==1) {
        enrich<-tmp
    } else {
        enrich<-rbind(enrich,tmp)
    }

}

dim(enrich)
# [1] 8  11

enrich$Category<-as.factor(enrich$Category)
enrich$Coefficient_zscore<-enrich$'Coefficient_z-score'
enrich$pheno<-factor(enrich$pheno,levels=c("CD","IBD","UC","HEIGHT"))

enrich<-merge(enrich,map,by.x="Category",by.y="id_map")

# # only colour the significant results:
enrich$signif<-enrich$Category
enrich$signif[which(enrich$Enrichment_p>threshold)]<-NA



# To compare cell-type groups, we use the coefficient z-score in an analysis containing the full baseline model - source:https://github.com/bulik/ldsc/wiki/Partitioned-Heritability

xmin<-min(0,enrich$Coefficient_zscore-enrich$Coefficient_std_error)
xmax<-max(enrich$Coefficient_zscore+enrich$Coefficient_std_error)

p1<-
ggplot(enrich, aes(y=Category, x=Coefficient_zscore,fill=Category)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) + 
  scale_fill_manual(values=c(colo)) +
  facet_grid( ~ pheno, scales = "free",space = "free") +
  xlim(xmin,xmax) + theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.title = element_blank(),legend.position="none",
   plot.margin = unit(c(1,1,1,1), "cm")) 


cd<-enrich[which(enrich$pheno=="CD"),c("Category","Coefficient_zscore")]
colnames(cd)[2]<-paste(colnames(cd)[2],"cd",sep="_")
uc<-enrich[which(enrich$pheno=="UC"),c("Category","Coefficient_zscore")]
colnames(uc)[2]<-paste(colnames(uc)[2],"uc",sep="_")
height<-enrich[which(enrich$pheno=="HEIGHT"),c("Category","Coefficient_zscore")]
colnames(height)[2]<-paste(colnames(height)[2],"height",sep="_")

comp<-merge(cd,uc,by="Category")
comp<-merge(comp,height,by="Category")

rm(cd,uc,height)

minlim<-min(comp$Coefficient_zscore_cd,comp$Coefficient_zscore_uc,comp$Coefficient_zscore_height)
maxlim<-max(comp$Coefficient_zscore_cd,comp$Coefficient_zscore_uc,comp$Coefficient_zscore_height)

p2<-ggplot(comp, aes(y=Coefficient_zscore_cd, x=Coefficient_zscore_height,color=Category)) + 
  geom_point(size=2) +  geom_smooth(method='lm') +
  scale_colour_manual(values=c(colo)) + 
  theme_bw() + theme(panel.background = element_blank(),legend.title = element_blank(),
                 axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(1,1,1,1), "cm")) + 
                 labs(x="Height Coefficient Z-score",y="CD Coefficient Z-score") + 
                 xlim(minlim-1,maxlim+1) + ylim(minlim-1,maxlim+1) + geom_abline(color="grey",linetype = "dashed")

p3<-ggplot(comp, aes(y=Coefficient_zscore_uc, x=Coefficient_zscore_height,color=Category)) + 
  geom_point(size=2) +  geom_smooth(method='lm') +
  scale_colour_manual(values=c(colo)) + 
  theme_bw() + theme(panel.background = element_blank(),legend.title = element_blank(),
                 axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(1,1,1,1), "cm")) + 
                 labs(x="Height Coefficient Z-score",y="UC Coefficient Z-score") + 
                 xlim(minlim-1,maxlim+1) + ylim(minlim-1,maxlim+1) + geom_abline(color="grey",linetype = "dashed")


p23<-ggarrange(p2,p3,nrow=2,labels=c("b","c"))
p123<-ggarrange(p1,p23,ncol=2,widths=c(4,1.5),labels=c("a"))

hg<-ceiling(length(levels(as.factor(enrich$Category)))*2.8)

if (hg<50) {
  p23<-ggarrange(p2,p3,ncol=2,labels=c("b","c"))
  p123<-ggarrange(p1,p23,nrow=2,labels=c("a"))
}

ggsave(
  "~/git/IIBDGC_GWAS/plots/eqtl_enrichment/Supplementary_Figure_eqtl_enrichment_analyses_decode.png",
  p123,
  width = 120,
  height = hg+40,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)

enrich_decode<-enrich
dim(enrich_decode)
# [1] 8 20


##############################
# combine all results:

enrich<-rbind(enrich_edQTL_GTEX,enrich_gwas_catalogue,enrich_hu_2021,enrich_IBDverse,enrich_macromap,enrich_decode,enrich_sparc)

dim(enrich)

fwrite(enrich,"~/git/IIBDGC_GWAS/plots/paper_tables/Source_data_enrichment_analyses_eqtl_pqtl.tsv",col.names=T,row.names=F,quote=F,sep="\t")

q("no")



# # compare pvalues between ld at 0.25 and 0.75


# # compare pvalues:
# ggplot(enrich, aes(y=Category, x=-log10(Enrichment_p),color=signif)) +
#   geom_point() + scale_y_discrete(limits=rev) + 
#   # geom_text(vjust = 0,hjust=0.9,size=3) +
#   scale_fill_manual(values=c(colo)) +
#   theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
#                  plot.margin = unit(c(1,1,1,1), "cm")) + facet_grid(pheno ~ ld, space = "free") + theme(legend.position="right",legend.title=element_blank())


# ggplot(enrich, aes(y=Category, x=Coefficient_zscore,fill=signif)) +
#   geom_bar(stat="identity", color="black", position=position_dodge()) + scale_y_discrete(limits=rev) + 
#   scale_fill_manual(values=c(colo)) +
#   facet_grid(pheno ~ ld, scales = "free",space = "free")


# cd<-ggplot(enrich[which(enrich$pheno=="CD")], aes(y=Category, x=-log10(Enrichment_p),color=signif)) +
#   geom_point() + scale_y_discrete(limits=rev) + 
#   # geom_text(vjust = 0,hjust=0.9,size=3) +
#   scale_fill_manual(values=c(colo)) +
#   theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
#                  plot.margin = unit(c(1,1,1,1), "cm")) + facet_grid(~ ld, space = "free") + theme(legend.position="right",legend.title=element_blank())


# ibd<-ggplot(enrich[which(enrich$pheno=="IBD")], aes(y=Category, x=-log10(Enrichment_p),color=signif)) +
#   geom_point() + scale_y_discrete(limits=rev) + 
#   # geom_text(vjust = 0,hjust=0.9,size=3) +
#   scale_fill_manual(values=c(colo)) +
#   theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
#                  plot.margin = unit(c(1,1,1,1), "cm")) + facet_grid(~ ld, space = "free") + theme(legend.position="right",legend.title=element_blank())


# uc<-ggplot(enrich[which(enrich$pheno=="UC")], aes(y=Category, x=-log10(Enrichment_p),color=signif)) +
#   geom_point() + scale_y_discrete(limits=rev) + 
#   # geom_text(vjust = 0,hjust=0.9,size=3) +
#   scale_fill_manual(values=c(colo)) +
#   theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
#                  plot.margin = unit(c(1,1,1,1), "cm")) + facet_grid(~ ld, space = "free") + theme(legend.position="right",legend.title=element_blank())

# p1<-ggarrange(cd,ibd,uc,nrow=3,common.legend=T,legend="right",labels = c("cd", "ibd","uc"),align = "v")


# ggsave(
#   "~/git/IIBDGC_GWAS/plots/test_ld_thresholds_eqtl_enrichent.png",
#   p1,
#   width = 90,
#   height = 150,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=2
# )



