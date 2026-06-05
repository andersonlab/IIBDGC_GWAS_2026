# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=25000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R


library(data.table)
library(ggplot2)
library(ggpubr)

path_gwas="/path/to/ibdgwas/IIBDGC/"


ph<-"ibd"
ph<-"cd"

for (chr in 1:22) {
    tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/",ph,"/",chr,"_",ph,"_meta_eur_tier_1_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz"))

    if(chr==1) {
        dat<-tmp
    } else {
        dat<-rbind(dat,tmp)
    }
    rm(tmp)

}

cor.test(dat$N, dat$Neff)
# 0.9995283
rm(dat)

ph<-"uc"

for (chr in 1:22) {
    tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/",ph,"/",chr,"_",ph,"_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz"))

    if(chr==1) {
        dat<-tmp
    } else {
        dat<-rbind(dat,tmp)
    }
    rm(tmp)

}

cor.test(dat$N, dat$Neff)
# 0.9529676 - IBD
# 0.8973351 - CD

# plot comparison N vs Neff
p<-ggplot(dat, aes(x=N, y=Neff)) + geom_point()
ggsave(
  paste0("~/git/IIBDGC_GWAS/plots/test_variants_ldsc/comparison_N_Neff_tier_2_",ph,".png"),
  p,
  width = 90,
  height = 90,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)

# plot comparison rate N vs Neff
p<-ggplot(dat, aes(x=rate_total_sample, y=rate_Neff)) + geom_point()

ggsave(
  "~/git/IIBDGC_GWAS/plots/test_variants_ldsc/comparison_rate_N_vs_Neff_eur_tier_2.png",
  p,
  width = 90,
  height = 90,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)

dat$MAF<-pmin(dat$avgA2FREQ_CONTROLS,1-dat$avgA2FREQ_CONTROLS)

# do we miss genomewide significant variants at differnet thresholds?
dat$pvalue<-as.numeric(dat$'P-value')

dim(dat[which(dat$pvalue<=5E-8),])
# [1] 63675    24
# [1] 62416    25 - CD

dim(dat[which(dat$pvalue<=5E-8 & dat$rate_Neff>=0.5),])
# [1] 52404    24
# [1] 40384    25 - CD


dim(dat[which(dat$pvalue<=5E-8 & dat$INFO>=0.9 & dat$rate_Neff>=0.5),])
# [1] 52328    24
# [1] 40319    25 - CD

tmp<-dat[which(dat$pvalue<=5E-8 & dat$INFO>=0.9 & dat$rate_Neff>=0.5),]

dim(dat[which(dat$pvalue<=5E-8 & dat$rate_Neff>=0.5 & !dat$MarkerName %in% tmp$MarkerName),])
# [1] 76 24
# [1] 65 25  - CD

range(dat$INFO[which(dat$pvalue<=5E-8 & dat$rate_Neff>=0.5 & !dat$MarkerName %in% tmp$MarkerName)])
# [1] 0.6003971 0.8996510
# [1] 0.6009698 0.8994308



# get number of SNPs at different Neff thresholds for ldsc
sum_dat<-as.data.frame(matrix(ncol=9,nrow=5))

colnames(sum_dat)<-c("Neff_rate_threshold","Nsnps","mean_info","sd_info","min_info","max_info","mean_MAF","min_MAF","max_MAF")
sum_dat$Neff_rate_threshold<-c(0.5,0.6,0.7,0.8,0.9)



for (i in 1:nrow(sum_dat)) {
    sum_dat$Nsnps[i]<-nrow(dat[which(dat$rate_Neff>=sum_dat$Neff_rate_threshold[i])])
    sum_dat$mean_info[i]<-mean(dat$INFO[which(dat$rate_Neff>=sum_dat$Neff_rate_threshold[i])])
    sum_dat$sd_info[i]<-sd(dat$INFO[which(dat$rate_Neff>=sum_dat$Neff_rate_threshold[i])])
    sum_dat$min_info[i]<-min(dat$INFO[which(dat$rate_Neff>=sum_dat$Neff_rate_threshold[i])])
    sum_dat$max_info[i]<-max(dat$INFO[which(dat$rate_Neff>=sum_dat$Neff_rate_threshold[i])])
    sum_dat$mean_MAF[i]<-mean(dat$MAF[which(dat$rate_Neff>=sum_dat$Neff_rate_threshold[i])])
    sum_dat$min_MAF[i]<-min(dat$MAF[which(dat$rate_Neff>=sum_dat$Neff_rate_threshold[i])])
    sum_dat$max_MAF[i]<-max(dat$MAF[which(dat$rate_Neff>=sum_dat$Neff_rate_threshold[i])])
}

sum_dat
#   Neff_rate_threshold    Nsnps mean_info    sd_info  min_info  max_info
# 1                 0.5 13649273 0.9679816 0.04783721 0.4078756 1.0000000
# 2                 0.6 13174214 0.9713733 0.03929956 0.4211392 1.0000000
# 3                 0.7 12690664 0.9728186 0.03721032 0.4306106 0.9999987
# 4                 0.8 10904574 0.9779091 0.02705553 0.6608554 0.9999987
# 5                 0.9  9504812 0.9829610 0.02036634 0.6998427 0.9999987
#    mean_MAF      min_MAF   max_MAF
# 1 0.1065753 0.0001082629 0.4999998
# 2 0.1095151 0.0001455672 0.4999998
# 3 0.1109274 0.0001853431 0.4999998
# 4 0.1098416 0.0005725001 0.4999998
# 5 0.1228781 0.0007170901 0.4999998

p1<-ggplot(sum_dat, aes(x=Neff_rate_threshold, y=Nsnps)) +
  geom_bar(stat="identity", color="black", position=position_dodge())

p2<-ggplot(sum_dat, aes(x=Neff_rate_threshold, y=mean_info)) +
geom_pointrange(aes(ymin=mean_info-sd_info, ymax=mean_info+sd_info)) + ylim(0.4,1.05) + ylab("Mean Info (sd)")

p3<-ggplot(sum_dat, aes(x=Neff_rate_threshold, y=mean_info)) +
geom_pointrange(aes(ymin=min_info, ymax=max_info)) + ylim(0.4,1.05) + ylab("Mean Info (range)")

p4<-ggplot(sum_dat, aes(x=Neff_rate_threshold, y=mean_MAF)) +
geom_pointrange(aes(ymin=min_MAF, ymax=max_MAF))+ ylab("Mean MAF (range)")

p_na<-ggarrange(p1,p2,p3,p4,ncol=1)
p_na

# ggsave(
#   "~/git/IIBDGC_GWAS/plots/test_variants_ldsc/comparison_Neff_threshold_on_NSNPS_and_INFO_eur_tier_2.png",
#   p,
#   width = 90,
#   height = 90,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=2
# )


# get number of SNPs at different Neff thresholds for ldsc
sum_dat<-as.data.frame(matrix(ncol=9,nrow=5))

colnames(sum_dat)<-c("Neff_rate_threshold","Nsnps","mean_info","sd_info","min_info","max_info","mean_MAF","min_MAF","max_MAF")
sum_dat$Neff_rate_threshold<-c(0.5,0.6,0.7,0.8,0.9)



for (i in 1:nrow(sum_dat)) {
    sum_dat$Nsnps[i]<-nrow(dat[which(dat$rate_Neff>=sum_dat$Neff_rate_threshold[i] & dat$INFO>0.9)])
    sum_dat$mean_info[i]<-mean(dat$INFO[which(dat$rate_Neff>=sum_dat$Neff_rate_threshold[i] & dat$INFO>0.9)])
    sum_dat$sd_info[i]<-sd(dat$INFO[which(dat$rate_Neff>=sum_dat$Neff_rate_threshold[i] & dat$INFO>0.9)])
    sum_dat$min_info[i]<-min(dat$INFO[which(dat$rate_Neff>=sum_dat$Neff_rate_threshold[i] & dat$INFO>0.9)])
    sum_dat$max_info[i]<-max(dat$INFO[which(dat$rate_Neff>=sum_dat$Neff_rate_threshold[i] & dat$INFO>0.9)])
    sum_dat$mean_MAF[i]<-mean(dat$MAF[which(dat$rate_Neff>=sum_dat$Neff_rate_threshold[i] & dat$INFO>0.9)])
    sum_dat$min_MAF[i]<-min(dat$MAF[which(dat$rate_Neff>=sum_dat$Neff_rate_threshold[i] & dat$INFO>0.9)])
    sum_dat$max_MAF[i]<-max(dat$MAF[which(dat$rate_Neff>=sum_dat$Neff_rate_threshold[i] & dat$INFO>0.9)])
}

sum_dat
#   Neff_rate_threshold    Nsnps mean_info    sd_info  min_info  max_info
# 1                 0.5 12666068 0.9783933 0.02307997 0.9000000 1.0000000
# 2                 0.6 12397221 0.9789356 0.02267392 0.9000000 1.0000000
# 3                 0.7 12020316 0.9794963 0.02226886 0.9000000 0.9999987
# 4                 0.8 10635315 0.9805992 0.02095939 0.9000000 0.9999987
# 5                 0.9  9421481 0.9838892 0.01778440 0.9000003 0.9999987
#    mean_MAF      min_MAF   max_MAF
# 1 0.1137476 0.0001082629 0.4999998
# 2 0.1153779 0.0001455672 0.4999998
# 3 0.1163122 0.0001853431 0.4999998
# 4 0.1119665 0.0005725001 0.4999998
# 5 0.1234498 0.0007170901 0.4999998

p1<-ggplot(sum_dat, aes(x=Neff_rate_threshold, y=Nsnps)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +  xlab("Neff_rate_threshold (INFO>0.9)")

p2<-ggplot(sum_dat, aes(x=Neff_rate_threshold, y=mean_info)) +
geom_pointrange(aes(ymin=mean_info-sd_info, ymax=mean_info+sd_info)) + ylim(0.4,1.05) + ylab("Mean Info (sd)") + xlab("Neff_rate_threshold (INFO>0.9)")

p3<-ggplot(sum_dat, aes(x=Neff_rate_threshold, y=mean_info)) +
geom_pointrange(aes(ymin=min_info, ymax=max_info)) + ylim(0.4,1.05) + ylab("Mean Info (range)") + xlab("Neff_rate_threshold (INFO>0.9)")

p4<-ggplot(sum_dat, aes(x=Neff_rate_threshold, y=mean_MAF)) +
geom_pointrange(aes(ymin=min_MAF, ymax=max_MAF))+ ylab("Mean MAF (range)") +  xlab("Neff_rate_threshold (INFO>0.9)")

p_09<-ggarrange(p1,p2,p3,p4,ncol=1)

p<-ggarrange(p_na,p_09,ncol=2)
p

ggsave(
  "~/git/IIBDGC_GWAS/plots/test_variants_ldsc/comparison_Neff_threshold_and_INFO_0.9_on_NSNPS_and_INFO_eur_tier_2.png",
  p,
  width = 90,
  height = 90,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)

