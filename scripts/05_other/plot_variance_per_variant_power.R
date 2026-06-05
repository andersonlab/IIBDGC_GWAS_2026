# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

###############################################################################################################
### ADD A PLOT WITH THE CUMMULATIVE VARIABILITY PER SNP

# singularity exec iibdgc_postprocess_10_singularity.sif
# updated with eur tier 2 results 

MEM=2000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q cpu-interactive R 

library(data.table)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(rcartocolor)
library(qvalue)
library(pwr)


path_gwas="/path/to/ibdgwas/IIBDGC/"


# LOAD THE INDEPENDENT SIGNALS - see /Z/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_60_consolidate_table_with_all_independent_signals.R
all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants_join_conditinal_model.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 653 182


##################################################################
### additive contribution of proportion of variance explained:

################################################################
## evaluate consistency of effect size between EUR and non-EUR

all$Effect_eur<-NA
all$StdErr_eur<-NA
all$pvalue_eur<-NA

all$Effect_eur[which(all$phenotype=="IBD_unsaturated")]<-all$BETA_ibd_eur_tier_2[which(all$phenotype=="IBD_unsaturated")]
all$StdErr_eur[which(all$phenotype=="IBD_unsaturated")]<-all$SE_ibd_eur_tier_2[which(all$phenotype=="IBD_unsaturated")]
all$pvalue_eur[which(all$phenotype=="IBD_unsaturated")]<-all$"P-value_ibd_eur_tier_2"[which(all$phenotype=="IBD_unsaturated")]

all$Effect_eur[which(all$phenotype=="CD")]<-all$BETA_cd_eur_tier_2[which(all$phenotype=="CD")]
all$StdErr_eur[which(all$phenotype=="CD")]<-all$SE_cd_eur_tier_2[which(all$phenotype=="CD")]
all$pvalue_eur[which(all$phenotype=="CD")]<-all$"P-value_cd_eur_tier_2"[which(all$phenotype=="CD")]

all$Effect_eur[which(all$phenotype=="UC")]<-all$BETA_uc_eur_tier_2[which(all$phenotype=="UC")]
all$StdErr_eur[which(all$phenotype=="UC")]<-all$SE_uc_eur_tier_2[which(all$phenotype=="UC")]
all$pvalue_eur[which(all$phenotype=="UC")]<-all$"P-value_uc_eur_tier_2"[which(all$phenotype=="UC")]


test<-all[which(all$phenotype=="IBD_saturated"),c("MarkerName","P-value_ibd_eur_tier_2","P-value_cd_eur_tier_2","P-value_uc_eur_tier_2")]
test$model<-NA
test$model[which(test$"P-value_cd_eur_tier_2"<test$"P-value_uc_eur_tier_2")]<-"CD"
test$model[which(test$"P-value_cd_eur_tier_2">test$"P-value_uc_eur_tier_2")]<-"UC"
table(test$model)
# CD UC 
# 70 44  
cd_model<-test$MarkerName[which(test$model=="CD")]
uc_model<-test$MarkerName[which(test$model=="UC")]


# list of IBD-sat variants more significant for CD
all$Effect_eur[which(all$MarkerName %in% cd_model)]<-all$BETA_cd_eur_tier_2[which(all$MarkerName %in% cd_model)]
all$StdErr_eur[which(all$MarkerName %in% cd_model)]<-all$SE_cd_eur_tier_2[which(all$MarkerName %in% cd_model)]
all$pvalue_eur[which(all$MarkerName %in% cd_model)]<-all$"P-value_cd_eur_tier_2"[which(all$MarkerName %in% cd_model)]

# list of IBD-sat variants more significant for UC
all$Effect_eur[which(all$MarkerName %in% uc_model)]<-all$BETA_uc_eur_tier_2[which(all$MarkerName %in% uc_model)]
all$StdErr_eur[which(all$MarkerName %in% uc_model)]<-all$SE_uc_eur_tier_2[which(all$MarkerName %in% uc_model)]
all$pvalue_eur[which(all$MarkerName %in% uc_model)]<-all$"P-value_uc_eur_tier_2"[which(all$MarkerName %in% uc_model)]


#### no liability:
# ref: https://cloufield.github.io/gwaslab/PerSNPh2/

# all$snp_variance<-NA
# all$snp_variance<-(2*(all$A2FREQ_CONTROLS_ibd_eur_tier_2))*(1-(all$A2FREQ_CONTROLS_ibd_eur_tier_2))*(all$Effect_eur)^2
# all$snp_proportion_variance<-all$snp_variance/(all$snp_variance+3.29)
# all<-all[order(all$snp_proportion_variance,decreasing=T),]

# all$cum_snp_variance<-cumsum(all$snp_proportion_variance)
# all$order<-seq(1:nrow(all))

#### liability

# b_logit: beta estimate from logistic regression
# maf: minor allele frequency of the variant
# k: Population prevalence of the disease
# p: Sample case proportion (case / total in your sample)

all$maf<-pmin(all$A2FREQ_CONTROLS_ibd_eur_tier_2,1-all$A2FREQ_CONTROLS_ibd_eur_tier_2)

b_logit_to_r2_lib<-function(b_logit,maf,k,p){
 R2obs <- b_logit^2*p*(1-p)*2*maf*(1-maf)
 t<-qnorm(1-k)
 z<-dnorm(t)
 C <- (k*(1-k))^2/(z^2 * p*(1-p))
 theta <- z*(p-k)/(k*(1-k))*(z*(p-k)/(k*(1-k))-t)
 R2lib <- R2obs*C/(1 + R2obs*C*theta)
 return(R2lib)
}

all$snp_variance<-b_logit_to_r2_lib(all$Effect_eur,all$maf,0.007,0.5)
all$snp_proportion_variance<-all$snp_variance/(all$snp_variance+3.29)
all<-all[order(all$snp_proportion_variance,decreasing=T),]

vec_new<-c("new_cojo_supervised_gw_significant_multiancestry_new_signal","new_cojo_supervised_gw_significant_multiancestry_new_signal_new_exonic_variant",
"new_cojo_unsupervised_new_signal","new_cojo_unsupervised_new_signal_new_exonic_variant")

all$class_signal_simplified<-"Known"
all$class_signal_simplified[which(all$class_signal_final_exome %in% vec_new)]<-"New"

table(all$class_signal_simplified)
# Known   New 
#   246   407 

signal_colors <- c(
  Known    = "darkgrey",
  New    = "#004488"
)

all$cum_snp_variance<-cumsum(all$snp_proportion_variance)
all<-all[order(all$snp_proportion_variance,decreasing=T),]
all$order<-seq(1:nrow(all))

tmp1<-all[,c("cum_snp_variance","order")]

tmp2<-tmp1[1:250,]

tmp2$cum_snp_variance<-NA
tmp2$order<-seq(1:250)+653


tmp<-rbind(tmp1,tmp2)

# Fit logarithmic curve using nls()
fit <- nls(order ~ a + b*(cum_snp_variance),start = list(a = 0, b = -20), data = tmp)



# Generate points for the fitted curve
curve_data <- data.frame(x = seq(min(tmp$order), max(tmp$order), length.out = 653))
colnames(curve_data)<-"order"
curve_data$cum_snp_variance <- predict(fit, newdata = curve_data)


na <- all[is.na(all$cum_snp_variance), "order", drop = FALSE]
predict(fit, newdata = na)


# Plot the data and the fitted curve
curve_data$y_normalized<-(max(tmp$cum_snp_variance,na.rm=T)-min(tmp$cum_snp_variance,na.rm=T))*((curve_data$cum_snp_variance)-min(curve_data$cum_snp_variance))/(max(curve_data$cum_snp_variance)-min(curve_data$cum_snp_variance))+min(tmp$cum_snp_variance,na.rm=T)


p1<-ggplot(all, aes(y=cum_snp_variance,x=order,color=class_signal_simplified)) +
    geom_point(size=1.5) + 
    # ylim(0,0.6) +
    xlim(0,753) +
  #  geom_line(data = curve_data, aes(x = order, y = y_normalized), color = "red", linetype = "dashed") +
    theme(legend.position="right",legend.title=element_blank(),panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank(),
        plot.title = element_text(face = "bold",size = 12),plot.margin = margin(t = 0.5, b = 0.01, r = 0.8, l = 1,unit = "cm")) + 
        xlab("Number of independent signals") + ylab("Cummulative SNP heritability") +
   scale_color_manual(values=signal_colors)





# all<-all[order(all$maf,decreasing=T),]
# all$cum_snp_variance_by_MAF<-cumsum(all$snp_proportion_variance)

# all1<-all[which(all$class_signal_simplified=="Known"),]
# all1<-all1[order(all1$maf,decreasing=T),]
# all1$cum_snp_variance_by_MAF<-cumsum(all1$snp_proportion_variance)


# all2<-all[which(all$class_signal_simplified!="Known"),]
# all2<-all2[order(all2$maf,decreasing=T),]


# all$cum_snp_variance_by_MAF<-cumsum(all$snp_proportion_variance)

# all<-rbind(all1,all2)

# p1<-ggplot(all, aes(x=-log10(maf),y=(cum_snp_variance_by_MAF),color=class_signal_simplified)) +
#     geom_point() + 
#     # ylim(0,0.6) +
#     theme(legend.position="none",legend.title=element_blank(),panel.grid.major = element_blank(), 
#     panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank(),
#         plot.title = element_text(face = "bold",size = 12)) + 
#         xlab("-log10(MAF)") + ylab("Cummulative proportion of \nheritability per variant (R2)") +
#         scale_color_manual(values=signal_colors)+facet_grid(~class_signal_simplified)

### plot SNP variance distribution:

# summary(all$snp_variance)
# #      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# # 1.254e-06 1.210e-04 1.861e-04 3.420e-04 3.432e-04 8.220e-03 

# all$maf_range<-NA
# all$maf_range[which(all$maf<=0.01)]<-'0.001-0.01'
# all$maf_range[which(all$maf>0.01 & all$maf<=0.05)]<-'0.01-0.05'
# all$maf_range[which(all$maf>0.05 & all$maf<=0.5)]<-'0.05-0.5'
# table(all$maf_range,useNA="ifany")

# all$maf_range<-as.factor(all$maf_range)


# p2<-ggplot(all, aes(x=maf_range, y=(snp_variance))) + 
#   ylim(0,0.002) +
#   geom_boxplot()+theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Heritability per variant")

p2<-ggplot(all, aes(y=snp_variance)) + 
  ylim(0,0.002) +
  geom_boxplot()+ theme_bw() + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) + ylab("SNP heritability")



# https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS3.html

# ibd
# max=240386
# cd
# max=118922
# uc
# max=145099

variance<-c(5e-6,1E-5,5e-5,1e-4,5e-4,1e-3)
variance<-seq(0.000025, 0.000175, 0.00005)

# f<-c(0.001,0.005,0.01,0.05)

f<-c(0.001)

dat<-as.data.frame(matrix(ncol=4,nrow=0))
colnames(dat)<-c("sample_size","power","variance","frequency")

for (i in 1:length(variance)) {

    for (j in 1:length(f)) {

        b.alt = sqrt(variance[i] / (2*f[j]*(1 - f[j])) ) # this is beta that explains 0.5%
        sigma = sqrt(1 - variance[i]) # error sd after SNP effect is accounted for
        ns = seq(50000, 1000000, 5000) # candidate n
        ses = sigma / sqrt( ns*2*f[j]*(1 - f[j]) ) # SE corresponding to each n

        q.thresh = qchisq(5e-8, df = 1, ncp = 0, lower = FALSE) # threshold corresp. alpha = 5e-8
        pwr = pchisq(q.thresh, df = 1, ncp = (b.alt/ses)^2, lower = FALSE) # power at alpha = 5e-8

        tmp<-as.data.frame(cbind(ns,pwr))
        colnames(tmp)<-c("sample_size","power")
        tmp$variance<-variance[i]
        tmp$frequency<-f[j]

        dat<-rbind(dat,tmp)
        rm(tmp)
    }

}



dat$variance<-formatC(dat$variance,format = "E",digits=2)
names(table(dat$variance))
dat$variance<-as.factor(dat$variance)
dat$variance<-factor(dat$variance,levels=c("2.50E-05","7.50E-05","1.25E-04","1.75E-04","2.25E-04","2.75E-04","3.25E-04","3.75E-04","4.25E-04","4.75E-04"))


# plot Neff:
  # CD    = "#004488",
  # UC    = "#BB5566",
  # IBD = "#DDAA33"

  # Neff IBD 240386, Neff CD 118922; Neff UC 145099


p3<-ggplot(dat, aes(x = sample_size, y = power)) +
  geom_line(color = "darkblue", linewidth = 1) +
  labs(
    x = "sample_size",
    y = "power",
  ) +
  scale_x_continuous(
    limits = c(50000, 1000000),
    breaks = seq(50000, 1000000, by = 100000)
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),strip.text.x=element_blank()) +  
  geom_hline(yintercept = 0.9,
           linetype = "dotted",
           color = "grey50") + 
    facet_grid(variance ~ 1) + 
    geom_vline(xintercept=240386,color="#DDAA33") + 
    geom_vline(xintercept=145099,color="#BB5566") +
    geom_vline(xintercept=118922,color="#004488")



p12<-ggarrange(p1,p2,ncol=2,widths=c(5,2),labels=c("a","b"))
pall<-ggarrange(p12,p3,ncol=2,widths=c(5,1.5),labels=c("","c"))

ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_variance_per_lead_cummulative_power.png",
  pall,
  width = 160,
  height = 50,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)


q("no")
