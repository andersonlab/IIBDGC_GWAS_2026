# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# compare ldscore from GSA, 1000GP and UKB

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=15000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G ibdgwas R

library(data.table)
library(ggplot2)
library(ggpubr)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("ibd","cd","uc")

# get hapmap variants:

list_hapmap_variants<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/w_hm3.snplist",sep=""))

for (chr in c(1:22)) {

    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/GRCh38/plink_files/1000G.EUR.hg38.",chr,".bim",sep=""))
    tmp<-tmp[which(tmp$V2 %in% list_hapmap_variants$SNP),]

    if(chr==1) {
        all<-tmp
    } else {
        all<-rbind(all,tmp)
    }
    rm(tmp)
}

# double check that Alele notation in bim does not reflect frequency
dim(all)
# [1] 1189841       6

all$chrposid<-paste("chr",all$V1,":",all$V4,":",all$V6,":",all$V5,sep="")
all$X<-paste("chr",all$V1,":",all$V4,sep="")

colnames(all)[2]<-"SNP"

dim(all)
# [1] 1189841       8

# UKB ldsc
ukbldsc<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.ldscore.gz")
ukbldsc<-ukbldsc[which(ukbldsc$SNP %in% all$SNP),]
dim(ukbldsc)
# [1] 1075486       8

ukbldsc<-as.data.frame(ukbldsc)
colnames(ukbldsc)[3:4]<-paste(colnames(ukbldsc)[3:4],"_ukb",sep="")

# GSA ldsc
for (chr in c(1:22)) {
    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/eur_tier2_rate_0.5_info_0.9_gsa/weights.april2025.",chr,".l2.ldscore.gz",sep=""))
    tmp<-tmp[which(tmp$SNP %in% all$chrposid),]

    if(chr==1) {
        gsaldsc<-tmp
    } else {
        gsaldsc<-rbind(gsaldsc,tmp)
    }

}
dim(gsaldsc)
# [1] 832715      4
gsaldsc<-as.data.frame(gsaldsc)
colnames(gsaldsc)[2]<-"chrposid"
gsaldsc<-merge(gsaldsc,all[,c("SNP","chrposid")],by="chrposid",all.x=T)
dim(gsaldsc)
[1] 832715      5

colnames(gsaldsc)[3:4]<-paste(colnames(gsaldsc)[3:4],"_gsa",sep="")


# allarrays ldsc
for (chr in c(1:22)) {
    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/eur_tier2_rate_0.5_info_0.9_allarrays/weights.april2025.",chr,".l2.ldscore.gz",sep=""))
    tmp<-tmp[which(tmp$SNP %in% all$chrposid),]

    if(chr==1) {
        allarraysldsc<-tmp
    } else {
        allarraysldsc<-rbind(allarraysldsc,tmp)
    }

}
dim(allarraysldsc)
# [1] 832734      4
allarraysldsc<-as.data.frame(allarraysldsc)
colnames(allarraysldsc)[2]<-"chrposid"
allarraysldsc<-merge(allarraysldsc,all[,c("SNP","chrposid")],by="chrposid",all.x=T)
dim(allarraysldsc)
# [1] 832734      5

colnames(allarraysldsc)[3:4]<-paste(colnames(allarraysldsc)[3:4],"_allarrays",sep="")



all_hm3<-merge(ukbldsc[,c("SNP","L2_ukb")],gsaldsc[,c("SNP","L2_gsa")],by="SNP")
all_hm3<-merge(all_hm3,allarraysldsc[,c("SNP","L2_allarrays")],by="SNP")

dim(all_hm3)
# [1] 750676      4

cor.test(all_hm3$L2_ukb, all_hm3$L2_gsa, method=c("pearson"))
# 0.9857921 

cor.test(all_hm3$L2_ukb, all_hm3$L2_allarrays, method=c("pearson"))
# 0.9921305

cor.test(all_hm3$L2_gsa, all_hm3$L2_allarrays, method=c("pearson"))
# 0.9982494


p1<-ggplot(all_hm3, aes(y=L2_gsa,x=L2_ukb)) + xlim(0,4000) + ylim(0,4000) + 
  geom_point() + stat_cor(method="pearson",label.x = 2500,label.y = 0) +
  geom_smooth(method=lm)


# compare GSA LD vs allarrays LD - with no variant constraint:


# GSA ldsc
for (chr in c(1:22)) {
    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/eur_tier2_rate_0.5_info_0.9_gsa/weights.april2025.",chr,".l2.ldscore.gz",sep=""))
    # tmp<-tmp[which(tmp$SNP %in% all$chrposid),]

    if(chr==1) {
        gsaldsc<-tmp
    } else {
        gsaldsc<-rbind(gsaldsc,tmp)
    }

}
dim(gsaldsc)
# [1] 12960265        4
gsaldsc<-as.data.frame(gsaldsc)
colnames(gsaldsc)[2]<-"chrposid"
colnames(gsaldsc)[3:4]<-paste(colnames(gsaldsc)[3:4],"_gsa",sep="")


# allarrays ldsc
for (chr in c(1:22)) {
    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/eur_tier2_rate_0.5_info_0.9_allarrays/weights.april2025.",chr,".l2.ldscore.gz",sep=""))
    # tmp<-tmp[which(tmp$SNP %in% all$chrposid),]

    if(chr==1) {
        allarraysldsc<-tmp
    } else {
        allarraysldsc<-rbind(allarraysldsc,tmp)
    }

}
dim(allarraysldsc)
# [1] 12965097      4
allarraysldsc<-as.data.frame(allarraysldsc)
colnames(allarraysldsc)[2]<-"chrposid"
colnames(allarraysldsc)[3:4]<-paste(colnames(allarraysldsc)[3:4],"_allarrays",sep="")

all<-merge(allarraysldsc[,c("chrposid","L2_allarrays")],gsaldsc[,c("chrposid","L2_gsa")],by="chrposid")
dim(all)
# [1] 12960265        3

cor.test(all$L2_gsa, all$L2_allarrays, method=c("pearson"))
#       cor 
# 0.9956192 

max(all$L2_gsa,all$L2_allarrays)


p2<-ggplot(all, aes(y=L2_gsa,x=L2_allarrays)) + xlim(0,4000) + ylim(0,4000) + 
  geom_point() + stat_cor(method="pearson",label.x = 2500,label.y = 0) +
  geom_smooth(method=lm)

p<-ggarrange(p1,p2,ncol=2)

ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_ldsc_reference_comparison.png",
  p,
  width = 150,
  height = 90,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=4
)

q("no")