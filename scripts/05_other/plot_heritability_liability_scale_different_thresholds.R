# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# # singularity exec iibdgc_postprocess_10_singularity.sif
# updated with eur tier 2 results 

MEM=2000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R 

library(data.table)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(rcartocolor)
library(qvalue)
library(pwr)


path_gwas="/path/to/ibdgwas/IIBDGC/"


files<-list.files(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/heritability/",sep=""))
files<-files[grepl("_info_0.9",files)]
files<-files[grep("rate_0.5",files)]
files<-files[grep("eur_tier2",files)]
files<-files[grep("allarrays",files)]
files<-files[!grepl("no_liability_scale",files)]


rm(dat_final)
for (i in c(1:length(files))) {

    dat<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/heritability/",files[i],sep=""),head=F,fill=T)

    dat<-dat[grep("Total Liability scale",dat$V1),]
    dat$file<-files[i]

    dat$maf_threshold<-gsub(".*_eur_tier2_list_variants_rate_0.5_info_0.9.*_allarrays","",files[i])
    dat$maf_threshold<-gsub("_maf_","",dat$maf_threshold)
    dat$maf_threshold<-gsub("_heritability_estimation_liability_scale.log","",dat$maf_threshold)
    dat$maf_threshold
    if(dat$maf_threshold=="") {
        dat$maf_threshold<-"0.001"
    }

    dat$weights_reference<-gsub(".*_eur_tier2_list_variants_rate_0.5_info_0.9_","",dat$file)
    dat$weights_reference<-gsub("_.*","",dat$weights_reference)


    dat$V1<-gsub("Total Liability scale h2: ","",dat$V1)
    dat$h2<-unlist(strsplit(dat$V1," "))[1]
    dat$h2_se<-unlist(strsplit(dat$V1," "))[2]
    dat$h2_se<-gsub("\\(","",dat$h2_se)
    dat$h2_se<-gsub(")","",dat$h2_se)
    dat$pheno<-toupper(gsub("_eur_tier2.*","",files[i]))

    if (exists("dat_final")) {
        dat_final<-rbind(dat_final,dat)
    } else {
        dat_final<-dat
    }
    rm(dat)
}

dat_final$maf_threshold<-factor(dat_final$maf_threshold,levels=c("0.001","0.005","0.01","0.05"))
dat_final$h2<-as.numeric(dat_final$h2)
dat_final$h2_se<-as.numeric(dat_final$h2_se)


dat_final$pheno<-factor(dat_final$pheno,levels=c("IBD","CD","UC"))


rm(dat_final_2)
for (i in c(1:length(files))) {

    dat<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/heritability/",files[i],sep=""),head=F,fill=T)

    dat<-dat[grep("After merging with regression SNP LD",dat$V1),]

    dat$nsnps<-gsub("After merging with regression SNP LD, ","",dat$V1)
    dat$nsnps<-gsub(" SNPs remain.","",dat$nsnps)

    dat$file<-files[i]

    if (exists("dat_final_2")) {
        dat_final_2<-rbind(dat_final_2,dat)
    } else {
        dat_final_2<-dat
    }
    rm(dat)
}

dat_final<-merge(dat_final,dat_final_2[,c("file","nsnps")],by="file")

dat_final$pheno<-factor(dat_final$pheno,levels=c("IBD","CD","UC"))
dat_final$nsnps<-as.numeric(dat_final$nsnps)


rm(dat_final_3)
for (i in c(1:length(files))) {

    dat<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/heritability/",files[i],sep=""),head=F,fill=T)

    dat<-dat[grep("Intercept:",dat$V1),]
    dat$file<-files[i]

    dat$V1<-gsub("Intercept: ","",dat$V1)
    dat$intercept<-unlist(strsplit(dat$V1," "))[1]
    dat$intercept_se<-unlist(strsplit(dat$V1," "))[2]
    dat$intercept_se<-gsub("\\(","",dat$intercept_se)
    dat$intercept_se<-gsub(")","",dat$intercept_se)

    if (exists("dat_final_3")) {
        dat_final_3<-rbind(dat_final_3,dat)
    } else {
        dat_final_3<-dat
    }
    rm(dat)
}

dat_final<-merge(dat_final,dat_final_3[,c("file","intercept","intercept_se")],by="file")

dat_final$intercept<-as.numeric(dat_final$intercept)
dat_final$intercept_se<-as.numeric(dat_final$intercept_se)

dat_final$intercept_pvalue<-pnorm((dat_final$intercept-1)/dat_final$intercept_se, lower.tail=FALSE)



rm(dat_final_4)
for (i in c(1:length(files))) {

    dat<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/heritability/",files[i],sep=""),head=F,fill=T)

    dat<-dat[grep("Mean Chi",dat$V1),]
    dat$file<-files[i]

    dat$chi_sqr<-gsub("Mean Chi\\^2: ","",dat$V1)


    if (exists("dat_final_4")) {
        dat_final_4<-rbind(dat_final_4,dat)
    } else {
        dat_final_4<-dat
    }
    rm(dat)
}

dat_final<-merge(dat_final,dat_final_4[,c("file","chi_sqr")],by="file")

rm(dat_final_5)
for (i in c(1:length(files))) {

    dat<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/heritability/",files[i],sep=""),head=F,fill=T)

    dat<-dat[grep("Ratio",dat$V1),]
    dat$file<-files[i]

    dat$V1<-gsub("Ratio: ","",dat$V1)
    dat$ratio<-unlist(strsplit(dat$V1," "))[1]
    dat$ratio_se<-unlist(strsplit(dat$V1," "))[2]
    dat$ratio_se<-gsub("\\(","",dat$ratio_se)
    dat$ratio_se<-gsub(")","",dat$ratio_se)

    if(dat$ratio=="Ratio") {
        dat$ratio<-"<0"
        dat$ratio_se<-"NA"
    }

    if (exists("dat_final_5")) {
        dat_final_5<-rbind(dat_final_5,dat)
    } else {
        dat_final_5<-dat
    }
    rm(dat)
}

dat_final<-merge(dat_final,dat_final_5[,c("file","ratio","ratio_se")],by="file")


rm(dat_final_6)
for (i in c(1:length(files))) {

    dat<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/heritability/",files[i],sep=""),head=F,fill=T)

    dat<-dat[grep("Lambda GC",dat$V1),]
    dat$file<-files[i]

    dat$lambda_gc<-gsub("Lambda GC: ","",dat$V1)

    if (exists("dat_final_6")) {
        dat_final_6<-rbind(dat_final_6,dat)
    } else {
        dat_final_6<-dat
    }
    rm(dat)
}

dat_final<-merge(dat_final,dat_final_6[,c("file","lambda_gc")],by="file")
# save source data:
dat_final<-dat_final[,c("pheno","weights_reference","maf_threshold","h2","h2_se","nsnps","intercept","intercept_se","lambda_gc","ratio","ratio_se")]


# Consistent colours
phenotype_colors <- c(
  CD    = "#004488",
  UC    = "#BB5566",
  IBD = "#DDAA33"
)


p1<-ggplot(dat_final[which(dat_final$weights_reference=="allarrays"),], aes(x=maf_threshold, y=h2,fill=pheno)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  geom_errorbar(aes(ymin=h2-h2_se, ymax=h2+h2_se), width=.2,
                 position=position_dodge(.9)) + theme(axis.text.x=element_text(angle=45,hjust=1),plot.margin = unit(c(1,1,1,1), "cm")) +
                 facet_grid(,vars(pheno)) + theme(legend.position="none") + 
                 labs(y="Variant heritability\n(liability scale)",x="MAF threshold") + 
                 scale_fill_manual(values=phenotype_colors) +
                 scale_color_manual(values=phenotype_colors)

p2<-ggplot(dat_final[which(dat_final$weights_reference=="allarrays"),], aes(x=maf_threshold, y=nsnps,fill=pheno)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  geom_errorbar(aes(ymin=h2-h2_se, ymax=h2+h2_se), width=.2,
                 position=position_dodge(.9)) + theme(axis.text.x=element_text(angle=45,hjust=1),plot.margin = unit(c(1,1,1,1), "cm")) +
                 facet_grid(,vars(pheno)) + theme(legend.position="none") + 
                 labs(y="log10(Number variants)",x="MAF threshold") + 
                 scale_fill_manual(values=phenotype_colors)+
                 scale_color_manual(values=phenotype_colors)






##################################
### enrich per region:

path_gwas="/path/to/ibdgwas/IIBDGC/"

files<-list.files(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",sep=""))
files<-files[grep(".results",files)]

files<-files[grep("no_baseline_both_annotations",files)]
files<-files[!grepl("annotations.files",files)]

for (i in 1:length(files)) {

    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/output/",files[i],sep=""),head=T)
    tmp$pheno<-toupper(gsub(".no_baseline_both_annotations.*","",files[i]))
    tmp$class<-gsub(".*baseline_both_annotations.","",files[i])
    tmp$class<-gsub(".files_version_april2025.results","",tmp$class)
    tmp$class<-paste0(c("inside_","outside_"),tmp$class)
    # tmp<-tmp[which(tmp$Category %in% c("L2_1","L2_0")),]
    # tmp$Category<-gsub("L2_0","",tmp$Category)
    # tmp$Category<-gsub("L2_0","",tmp$Category)

    if(i==1) {
        enrich<-tmp
    } else {
        enrich<-rbind(enrich,tmp)
    }

}

enrich

enrich$pheno<-factor(enrich$pheno,levels=c("IBD","CD","UC"))

enrich$regions<-NA
enrich$regions<-gsub("inside_ibd_","",enrich$class)
enrich$regions<-gsub("outside_ibd_","",enrich$regions)
enrich$regions[which(enrich$regions=="regions")]<-"all_regions"

enrich$localization<-NA
enrich$localization<-gsub("_ibd.*","",enrich$class)
enrich$localization[which(enrich$localization=="inside")]<-paste(enrich$localization[which(enrich$localization=="inside")],enrich$pheno[which(enrich$localization=="inside")],sep="_")


# save source data:
write.table(enrich,"/path/to/user/git/IIBDGC_GWAS/plots/paper_tables/Source_data_enrichment_analyses_ibd_regions.tsv",col.names=T,row.names=F,quote=F,sep="\t")



# INCLUDE THE IBD REGIONS HERITABILITY in the plot

enrich$regions<-factor(enrich$regions,levels=c("old_regions","all_regions"))

# Consistent colours
loc_colors <- c(
  inside_CD    = "#004488",
  inside_UC    = "#BB5566",
  inside_IBD = "#DDAA33",
  outside = "grey"
)


p3<-ggplot(enrich, aes(x=regions, y=Prop._h2,fill=localization)) +
  geom_bar(stat="identity", color="black",position = position_fill(reverse = TRUE)) + 
                 facet_grid(,vars(pheno)) + theme(legend.position="none") + 
                 labs(y="Proportion of variant heritability\nin IBD regions")  + 
                 scale_fill_manual(values=loc_colors)


p<-ggarrange(p1,p2,p3,nrow=3,align="v",labels=c("a","b","c"))


ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_heritability_eur_tier2.png",
  p,
  width = 90,
  height = 120,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)

q("no")
