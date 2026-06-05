# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# singularity exec iibdgc_postprocess_10_singularity.sif
MEM=8000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)


path_gwas="/path/to/ibdgwas/IIBDGC/"

ph<-"ibd"
tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_",ph,"_eur_tier2_list_variants_rate_0.5.tsv.gz"),head=T)

# get summary of N SNPs at different thresolds:
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

p1<-ggplot(sum_dat, aes(x=Neff_rate_threshold, y=Nsnps)) +
  geom_bar(stat="identity", color="black", position=position_dodge())
p2<-ggplot(sum_dat, aes(x=Neff_rate_threshold, y=mean_info)) +
geom_pointrange(aes(ymin=mean_info-sd_info, ymax=mean_info+sd_info)) + ylim(0.4,1.05) + ylab("Mean Info (sd)")
p3<-ggplot(sum_dat, aes(x=Neff_rate_threshold, y=mean_info)) +
geom_pointrange(aes(ymin=min_info, ymax=max_info)) + ylim(0.4,1.05) + ylab("Mean Info (range)")
p4<-ggplot(sum_dat, aes(x=Neff_rate_threshold, y=mean_MAF)) +
geom_pointrange(aes(ymin=min_MAF, ymax=max_MAF))+ ylab("Mean MAF (range)")

p_na<-ggarrange(p1,p2,p3,p4,ncol=1)

# get summary of N SNPs at different thresolds, by subsetting those at INFO > 0.9 as recommended by LDSC:
# Imputation quality is correlated with LD Score, and low imputation quality yields lower test statistics, so imputation quality is a confounder for LD Score regression. To prevent bias from variable imputation quality, we usually remove poorly-imputed SNPs by filtering on INFO > 0.9. T
# source: https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation

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
  geom_bar(stat="identity", color="black", position=position_dodge())
p2<-ggplot(sum_dat, aes(x=Neff_rate_threshold, y=mean_info)) +
geom_pointrange(aes(ymin=mean_info-sd_info, ymax=mean_info+sd_info)) + ylim(0.4,1.05) + ylab("Mean Info (sd)")
p3<-ggplot(sum_dat, aes(x=Neff_rate_threshold, y=mean_info)) +
geom_pointrange(aes(ymin=min_info, ymax=max_info)) + ylim(0.4,1.05) + ylab("Mean Info (range)")
p4<-ggplot(sum_dat, aes(x=Neff_rate_threshold, y=mean_MAF)) +
geom_pointrange(aes(ymin=min_MAF, ymax=max_MAF))+ ylab("Mean MAF (range)")
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

q("no")