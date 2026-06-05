# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# ~/scripts/other/plot_effect_sizes_eur_multiancestry.R

# LD estimated in UKBB between rare varians and GWAS index:

# Setup (run in shell before starting R):
# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM=12000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -q cpu-interactive -G your_hpc_group R

library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr)

rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

# LOAD THE INDEPENDENT SIGNALS - see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_38_2_define_independent_signals_using_as_ref_cojo_agnostic_results.R
all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants_join_conditinal_model.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 619 182

pval_threshold<-0.05/nrow(all)
# [1] 7.92393e-05 - P-value threshold


#########################################################################################

## get numbers replicated in SAS:

class(all$"P-value_ibd_sas")
class(all$"P-value_cd_sas")
class(all$"P-value_uc_sas")

dim(all[which(all$"P-value_ibd_sas"<5E-8 | all$"P-value_cd_sas"<5E-8 | all$"P-value_uc_sas"<5E-8 ),])
# [1]  0 205

all[which(all$"P-value_ibd_sas"<5E-8 | all$"P-value_cd_sas"<5E-8 | all$"P-value_uc_sas"<5E-8),c("MarkerName","phenotype",
"BETA_cd_eur_tier_2","BETA_uc_eur_tier_2",
"BETA_cd_sas","P-value_cd_sas","BETA_uc_sas","P-value_uc_sas")]
#

#########################################################################################

## get numbers replicated in EAS:


class(all$"P-value_ibd_eas")
class(all$"P-value_cd_eas")
class(all$"P-value_uc_eas")

dim(all[which(all$"P-value_ibd_eas"<5E-8 | all$"P-value_cd_eas"<5E-8 | all$"P-value_uc_eas"<5E-8 ),])
# [1] 36 205



################################################################
## evaluate consistency of effect size between EUR and non-EUR, retain the MR mega results

all$"P-value_cd_eur_tier_2"<-as.numeric(all$"P-value_cd_eur_tier_2")
all$"P-value_uc_eur_tier_2"<-as.numeric(all$"P-value_uc_eur_tier_2")
all$"P-value_ibd_eur_tier_2"<-as.numeric(all$"P-value_ibd_eur_tier_2")

all$"P-value_cd_eur_tier_2_eas_sas"<-as.numeric(all$"P-value_cd_eur_tier_2_eas_sas")
all$"P-value_uc_eur_tier_2_eas_sas"<-as.numeric(all$"P-value_uc_eur_tier_2_eas_sas")
all$"P-value_ibd_eur_tier_2_eas_sas"<-as.numeric(all$"P-value_ibd_eur_tier_2_eas_sas")



# Helper: assign phenotype-specific source columns into destination columns
assign_by_pheno <- function(df, rows, dest_cols, src_cols) {
  for (i in seq_along(dest_cols)) {
    df[[dest_cols[i]]][rows] <- df[[src_cols[i]]][rows]
  }
  df
}

all$Effect_eur    <- NA
all$StdErr_eur    <- NA
all$pvalue_eur    <- NA
all$pvalue_eur_het <- NA

eur_dest <- c("Effect_eur", "StdErr_eur", "pvalue_eur", "pvalue_eur_het")
all <- assign_by_pheno(all, all$phenotype=="IBD_unsaturated", eur_dest,
  c("BETA_ibd_eur_tier_2", "SE_ibd_eur_tier_2", "P-value_ibd_eur_tier_2", "HetPVal_ibd_eur_tier_2"))
all <- assign_by_pheno(all, all$phenotype=="CD", eur_dest,
  c("BETA_cd_eur_tier_2",  "SE_cd_eur_tier_2",  "P-value_cd_eur_tier_2",  "HetPVal_cd_eur_tier_2"))
all <- assign_by_pheno(all, all$phenotype=="UC", eur_dest,
  c("BETA_uc_eur_tier_2",  "SE_uc_eur_tier_2",  "P-value_uc_eur_tier_2",  "HetPVal_uc_eur_tier_2"))


test<-all[which(all$phenotype=="IBD_saturated"),c("MarkerName","P-value_ibd_eur_tier_2","P-value_cd_eur_tier_2","P-value_uc_eur_tier_2")]
test$model<-NA
test$model[which(test$"P-value_cd_eur_tier_2"<test$"P-value_uc_eur_tier_2")]<-"CD"
test$model[which(test$"P-value_cd_eur_tier_2">test$"P-value_uc_eur_tier_2")]<-"UC"
table(test$model)
# CD UC 
# 67 42   
cd_model<-test$MarkerName[which(test$model=="CD")]
uc_model<-test$MarkerName[which(test$model=="UC")]


# list of IBD-sat variants more significant for CD
all <- assign_by_pheno(all, all$MarkerName %in% cd_model, eur_dest,
  c("BETA_cd_eur_tier_2", "SE_cd_eur_tier_2", "P-value_cd_eur_tier_2", "HetPVal_cd_eur_tier_2"))

# list of IBD-sat variants more significant for UC
all <- assign_by_pheno(all, all$MarkerName %in% uc_model, eur_dest,
  c("BETA_uc_eur_tier_2", "SE_uc_eur_tier_2", "P-value_uc_eur_tier_2", "HetPVal_uc_eur_tier_2"))

all$MarkerName[which(all$pvalue_eur_het<=pval_threshold)]
# [1] "chr1:25386842:AACACAC:A" "chr1:67240275:G:A"      
# [3] "chr10:125794381:G:A"     "chr16:50712015:C:T"     
# [5] "chr2:233274822:G:A"      "chr5:159402286:C:T"     
# [7] "chr7:5433979:G:C"        "chr9:136365140:C:G"   



#### multiancestry:

all$Effect_eur_eas_sas              <- NA
all$StdErr_eur_eas_sas              <- NA
all$pvalue_eur_eas_sas              <- NA
all$chisq_eur_eas_sas_ancestry_het  <- NA
all$pvalue_eur_eas_sas_ancestry_het <- NA
all$pvalue_eur_eas_sas_het          <- NA

ma_dest <- c("Effect_eur_eas_sas", "StdErr_eur_eas_sas", "pvalue_eur_eas_sas",
             "chisq_eur_eas_sas_ancestry_het", "pvalue_eur_eas_sas_het")
all <- assign_by_pheno(all, all$phenotype=="IBD_unsaturated", ma_dest,
  c("BETA_ibd_eur_tier_2_eas_sas", "SE_ibd_eur_tier_2_eas_sas", "P-value_ibd_eur_tier_2_eas_sas",
    "chisq_ancestry_het_ibd_eur_tier_2_eas_sas_mrmega", "P-value_residual_het_ibd_eur_tier_2_eas_sas_mrmega"))
all <- assign_by_pheno(all, all$phenotype=="CD", ma_dest,
  c("BETA_cd_eur_tier_2_eas_sas", "SE_cd_eur_tier_2_eas_sas", "P-value_cd_eur_tier_2_eas_sas",
    "chisq_ancestry_het_cd_eur_tier_2_eas_sas_mrmega", "P-value_residual_het_cd_eur_tier_2_eas_sas_mrmega"))
all <- assign_by_pheno(all, all$phenotype=="UC", ma_dest,
  c("BETA_uc_eur_tier_2_eas_sas", "SE_uc_eur_tier_2_eas_sas", "P-value_uc_eur_tier_2_eas_sas",
    "chisq_ancestry_het_uc_eur_tier_2_eas_sas_mrmega", "P-value_residual_het_uc_eur_tier_2_eas_sas_mrmega"))



test<-all[which(all$phenotype=="IBD_saturated"),c("MarkerName","P-value_ibd_eur_tier_2_eas_sas","P-value_cd_eur_tier_2_eas_sas","P-value_uc_eur_tier_2_eas_sas")]
test$model<-NA
test$model[which(test$"P-value_cd_eur_tier_2_eas_sas"<test$"P-value_uc_eur_tier_2_eas_sas")]<-"CD"
test$model[which(test$"P-value_cd_eur_tier_2_eas_sas">test$"P-value_uc_eur_tier_2_eas_sas")]<-"UC"
table(test$model)
# CD UC 
# 66 38  
cd_model<-test$MarkerName[which(test$model=="CD")]
uc_model<-test$MarkerName[which(test$model=="UC")]


# list of IBD-sat variants more significant for CD
all <- assign_by_pheno(all, all$MarkerName %in% cd_model, ma_dest,
  c("BETA_cd_eur_tier_2_eas_sas", "SE_cd_eur_tier_2_eas_sas", "P-value_cd_eur_tier_2_eas_sas",
    "chisq_ancestry_het_cd_eur_tier_2_eas_sas_mrmega", "P-value_residual_het_cd_eur_tier_2_eas_sas_mrmega"))

# list of IBD-sat variants more significant for UC
all <- assign_by_pheno(all, all$MarkerName %in% uc_model, ma_dest,
  c("BETA_uc_eur_tier_2_eas_sas", "SE_uc_eur_tier_2_eas_sas", "P-value_uc_eur_tier_2_eas_sas",  # BUG FIX: was pvalue_eur_eas_sas_ancestry_het
    "chisq_ancestry_het_uc_eur_tier_2_eas_sas_mrmega", "P-value_residual_het_uc_eur_tier_2_eas_sas_mrmega"))



all$pvalue_eur_eas_sas_ancestry_het<-pchisq(all$chisq_eur_eas_sas_ancestry_het,1,lower.tail = FALSE)
summary(all$pvalue_eur_eas_sas_ancestry_het)


all$zscore_eur_eas_sas<-all$Effect_eur_eas_sas/all$StdErr_eur_eas_sas
all$zscore_eur<-all$Effect_eur/all$StdErr_eur

all$zscore_diff<-(all$zscore_eur_eas_sas-all$zscore_eur)

low_lim<-mean(all$zscore_diff)-4*sd(all$zscore_diff)
upp_lim<-mean(all$zscore_diff)+4*sd(all$zscore_diff)

all$zscore_variants<-NA
all$zscore_variants[which(all$zscore_diff>upp_lim | all$zscore_diff<low_lim)]<-all$MarkerName[which(all$zscore_diff>upp_lim | all$zscore_diff<low_lim)]

all$zscore_variants[which(!is.na(all$zscore_variants))]
# [1] "chr1:67182913:G:A"  "chr10:62710915:C:T" "chr4:38333446:G:A"
# [4] "chr7:5433979:G:C"   "chr9:114801407:T:G"

all[which(!is.na(all$zscore_variants)),]


#### number of index variants with significat heterogeneity due to ancestry:

dim(all[which(all$pvalue_eur_eas_sas_ancestry_het<pval_threshold),])
# [1]  20 219

619-20
# [1] 599


all$het_variants_eur<-NA
all$het_variants_eur[which(all$pvalue_eur_het<pval_threshold)]<-all$MarkerName[which(all$pvalue_eur_het<pval_threshold)]
all$MarkerName[which(all$pvalue_eur_het<pval_threshold)]
# [1] "chr1:25386842:AACACAC:A" "chr1:67240275:G:A"      
# [3] "chr10:125794381:G:A"     "chr16:50712015:C:T"     
# [5] "chr2:233274822:G:A"      "chr5:159402286:C:T"     
# [7] "chr7:5433979:G:C"        "chr9:136365140:C:G"  

all[which(all$pvalue_eur_het<pval_threshold),c("MarkerName","pvalue_eur_het","phenotype")]
#                  MarkerName pvalue_eur_het       phenotype
# 51  chr1:25386842:AACACAC:A      8.468e-09 IBD_unsaturated
# 68        chr1:67240275:G:A      1.271e-08   IBD_saturated
# 80      chr10:125794381:G:A      2.717e-13 IBD_unsaturated
# 221      chr16:50712015:C:T      3.668e-05   IBD_saturated
# 344      chr2:233274822:G:A      1.572e-05              CD
# 473      chr5:159402286:C:T      6.539e-05 IBD_unsaturated
# 560        chr7:5433979:G:C      9.217e-18 IBD_unsaturated
# 601      chr9:136365140:C:G      1.324e-10 IBD_unsaturated


all$het_variants_multiancestry<-NA
all$het_variants_multiancestry[which(all$pvalue_eur_eas_sas_het<pval_threshold)]<-all$MarkerName[which(all$pvalue_eur_eas_sas_het<pval_threshold)]

all$MarkerName[which(all$pvalue_eur_eas_sas_het<pval_threshold)]
# [1] "chr1:25386842:AACACAC:A" "chr1:67240275:G:A"      
# [3] "chr10:125794381:G:A"     "chr16:50712015:C:T"     
# [5] "chr2:233274822:G:A"      "chr7:5433979:G:C"       
# [7] "chr9:136365140:C:G" 

all[which(all$pvalue_eur_eas_sas_ancestry_het<pval_threshold),"MarkerName"]
#  [1] "chr1:161547427:G:A"  "chr1:67214256:C:CT"  "chr10:110425838:G:A"
#  [4] "chr10:62710915:C:T"  "chr12:40337211:A:G"  "chr13:26957130:T:C" 
#  [7] "chr16:50500888:G:T"  "chr17:42345821:G:C"  "chr19:1127616:G:C"  
# [10] "chr4:38333446:G:A"   "chr5:150872803:G:A"  "chr5:159402286:C:T" 
# [13] "chr5:40498420:T:TG"  "chr5:59021173:A:G"   "chr7:27192143:C:T"  
# [16] "chr7:5433979:G:C"    "chr9:114801407:T:G"  "chr9:114845889:A:G" 
# [19] "chr9:114851353:GA:G" "chr9:136395559:G:A" 

all$het_variants_multiancestry[which(!is.na(all$het_variants_multiancestry))]
# [1] "chr1:25386842:AACACAC:A" "chr1:67240275:G:A"      
# [3] "chr10:125794381:G:A"     "chr16:50712015:C:T"     
# [5] "chr2:233274822:G:A"      "chr7:5433979:G:C"       
# [7] "chr9:136365140:C:G"  

all$het_ancestry_variants_multiancestry<-NA
all$het_ancestry_variants_multiancestry[which(all$pvalue_eur_eas_sas_ancestry_het<pval_threshold)]<-all$MarkerName[which(all$pvalue_eur_eas_sas_ancestry_het<pval_threshold)]

all[which(!is.na(all$het_ancestry_variants_multiancestry)),c("MarkerName","pvalue_eur_eas_sas_ancestry_het")]
#               MarkerName pvalue_eur_eas_sas_ancestry_het
# 19   chr1:161547427:G:A                    4.478709e-05
# 67   chr1:67214256:C:CT                    9.869727e-16
# 80  chr10:110425838:G:A                    1.179665e-05
# 97   chr10:62710915:C:T                    5.157739e-06
# 149  chr12:40337211:A:G                    1.348097e-05
# 163  chr13:26957130:T:C                    1.302431e-07
# 221  chr16:50500888:G:T                    7.702614e-06
# 259  chr17:42345821:G:C                    4.036603e-05
# 293   chr19:1127616:G:C                    3.407317e-05
# 459   chr4:38333446:G:A                    2.315243e-10
# 479  chr5:150872803:G:A                    4.577539e-05
# 482  chr5:159402286:C:T                    4.724379e-06
# 496  chr5:40498420:T:TG                    1.445898e-05
# 499   chr5:59021173:A:G                    1.055971e-05
# 559   chr7:27192143:C:T                    8.142125e-06
# 570    chr7:5433979:G:C                    9.306888e-18
# 591   chr8:15592628:A:C                    3.128484e-06
# 605  chr9:114801407:T:G                    1.490668e-36
# 606  chr9:114845889:A:G                    4.237494e-06
# 607 chr9:114851353:GA:G                    5.319543e-08
# 614  chr9:136395559:G:A                    3.476275e-07


dim(all[which(!is.na(all$het_variants_eur) | !is.na(all$het_variants_multiancestry) | !is.na(all$het_ancestry_variants_multiancestry)),c("MarkerName","pvalue_eur_het","pvalue_eur_eas_sas_het","pvalue_eur_eas_sas_ancestry_het","phenotype")])
# [1] 26  5


### PLOT 1, effect size and stderr estimates for EUR and multiancestry study

# Consistent colours
phenotype_colors <- c(
  CD    = "#004488",
  UC    = "#BB5566",
  'IBD_saturated' = "#DDAA33",
  'IBD_unsaturated'      = "#db7107"
)

theme_gwas <- theme(
  legend.position="right", legend.title=element_blank(), legend.key=element_blank(),
  panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
  panel.background=element_blank(), axis.line=element_line(colour="darkgrey"),
  plot.title=element_text(face="bold", size=12),
  plot.margin=margin(t=0.5, b=0.2, r=0.2, l=0.6, unit="cm")
)

## multiancestry heterogeneity
p1 <- ggplot(all, aes(x=zscore_eur_eas_sas, y=-log10(pvalue_eur_eas_sas_ancestry_het), color=phenotype)) +
  geom_point(shape=20) +
  xlim(-50, 50) +
  ylim(0, 40) +
  geom_hline(yintercept=-log10(pval_threshold), linetype="dashed", color="grey", linewidth=0.5) +
  geom_label_repel(aes(label=het_ancestry_variants_multiancestry), show.legend=FALSE, fontface="italic") +
  theme_gwas +
  scale_color_manual(values=phenotype_colors) +
  xlab("Zscore EUR EAS SAS") +
  ylab("-log10(Ancestry Het Pvalue)")

# ## EUR tier 2 heterogeneity
# p2 <- ggplot(all, aes(x=zscore_eur, y=-log10(pvalue_eur_het), color=phenotype)) +
#   geom_point(shape=20) +
#   xlim(-50, 50) +
#   ylim(0, 40) +
#   geom_hline(yintercept=-log10(pval_threshold), linetype="dashed", color="grey", linewidth=0.5) +
#   geom_label_repel(aes(label=het_variants_eur), show.legend=FALSE, fontface="italic") +
#   theme_gwas +
#   scale_color_manual(values=phenotype_colors) +
#   xlab("Zscore EUR EAS SAS") +
#   ylab("-log10(Ancestry Het Pvalue)")


# no new variants
p3 <- ggplot(all, aes(x=zscore_eur_eas_sas, y=zscore_eur, color=phenotype)) +
  geom_point(shape=20) +
  xlim(-50, 50) +
  ylim(-50, 50) +
  geom_abline(intercept=0,       linetype=2, color="lightgrey") +
  geom_abline(intercept=low_lim, linetype=2, color="lightgrey") +
  geom_abline(intercept=upp_lim, linetype=2, color="lightgrey") +
  geom_label_repel(aes(label=zscore_variants), show.legend=FALSE, seed=13, point.padding=unit(8, "lines")) +
  theme_gwas +
  scale_color_manual(values=phenotype_colors) +
  ylab("Zscore EUR") +
  xlab("Zscore EUR EAS SAS")


p<-ggarrange(p3,p1,ncol=2,labels=c("a","b"),common.legend=T,legend="bottom")


ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_effect_size_eur_multiancestry.png",
  p,
  width = 180,
  height = 80,
  dpi = 600,
  units = c("mm"),
  limitsize = T,scale=2
)


#####################################
### explore some of those results:


################################################################################################################################################################################
# TNRC18 variant 7-5397122-C-T - variant dropped from analyses, but tag by the following index as reported by Mark


finngen<-c("chr7:4468400:T:G","chr7:4827597:G:C","chr7:4881854:T:C","chr7:4981470:C:G","chr7:5139807:G:A","chr7:5147785:G:A","chr7:5266889:G:T","chr7:5285792:A:T","chr7:5319942:C:T","chr7:5333452:G:A","chr7:5433979:G:C","chr7:5453006:C:T","chr7:5611768:C:T","chr7:6161112:CCAAT:C","chr7:6352123:G:C")
all[which(all$MarkerName %in% finngen),c("MarkerName","phenotype","pvalue_eur_eas_sas_ancestry_het","P-value_residual_het_ibd_eur_tier_2_eas_sas_mrmega","P-value_residual_het_cd_eur_tier_2_eas_sas_mrmega","P-value_residual_het_uc_eur_tier_2_eas_sas_mrmega")]
#           MarkerName phenotype pvalue_eur_eas_sas_ancestry_het
# 570 chr7:5433979:G:C        CD                    9.306888e-18
#     P-value_residual_het_ibd_eur_tier_2_eas_sas_mrmega
# 570                                        1.84314e-16
#     P-value_residual_het_cd_eur_tier_2_eas_sas_mrmega
# 570                                       2.28257e-05
#     P-value_residual_het_uc_eur_tier_2_eas_sas_mrmega
# 570                                       5.26792e-14


# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM=2000
# bsub -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q normal \
#   -e ${path_gwas}post_imputation/2022/log/variant_lookup_stderr \
#   -o ${path_gwas}post_imputation/2022/log/variant_lookup_stdout \
#   "Rscript /path/to/user/git/IIBDGC_GWAS/requests/lookup_requests.R chr7:5397122:C:T"

# bsub -J"ph" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
#   -e ${path_gwas}post_imputation/2022/log/subset_files_by_ld_leads_allpheno_stderr \
#   -o ${path_gwas}post_imputation/2022/log/subset_files_by_ld_leads_allpheno_stdout \
#   "Rscript ~/git/IIBDGC_GWAS/scripts/other/coloc_subset_results_by_ld_leads_noFM.R ${ph} > \
#   ${path_gwas}post_imputation/2022/log/coloc_subset_results_by_ld_leads_noFM.Rout"

################################################################################################################################################################################

# define those variants with high heterogenetiy p-val and extract the sumamry statistics for each study:
high_het<-all[which(!is.na(all$het_variants_eur) | !is.na(all$het_variants_multiancestry) | !is.na(all$het_ancestry_variants_multiancestry)),c("MarkerName","chr","pvalue_eur_het","pvalue_eur_eas_sas_het","pvalue_eur_eas_sas_ancestry_het","phenotype")]
dim(high_het)

array<-c("illumina370","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome","ukbb","finngen","decode","danish","eas","ukbb_sas","ibdbioresource","ibdbioresource_sas")

pheno<-c("ibd","cd","uc")

chrs<-high_het$chr
chrs<-chrs[!duplicated(chrs)]

rm(list = intersect(ls(), c("tmp", "tmp1", "tmp11", "ibd", "cd", "uc")))

for (ph in pheno) {

  print(ph)
  
  tmp <- NULL

  for (chr in chrs) {

    print(chr)

    tmp1 <- NULL

    for (jj in 1:length(array)) {
      
      print(array[jj])

      file_tmp11<-paste(path_gwas,"post_imputation/2022/analysis/regenie/",array[jj],"/",ph,"/chr",chr,"_",array[jj],"_eur_all_step2_",ph,"_eur_sex_PCs_firthse_",ph,".regenie.gz",sep="")
      
      if (array[jj]=="ukbb") {
        file_tmp11<-paste(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/",ph,"/chr",chr,"_ukbb_eur_step2_",ph,"_eur_sex_PCs_firthse_",ph,".regenie.gz",sep="")
      } else if (array[jj]=="ukbb_sas") {
        file_tmp11<-paste(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/",ph,"/chr",chr,"_ukbb_sas_step2_",ph,"_sas_sex_PCs_firthse_",ph,".regenie.gz",sep="")
      } else if (array[jj]=="finngen") {
        file_tmp11<-paste(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/finngen/",ph,"/allchr_finngen_r10_K11_",ph,".regenie.gz",sep="")
      } else if (array[jj]=="ibdbioresource") {
        file_tmp11<-paste(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/",ph,"/",chr,"_interval_ibdbioresource_eur_step2_",ph,"_eur_sex_PCs_firthse_",ph,".regenie.gz",sep="")
      } else if (array[jj]=="ibdbioresource_sas") {
        file_tmp11<-paste(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/",ph,"/",chr,"_interval_ibdbioresource_noneur_step2_",ph,"_sas_sex_PCs_firthse_",ph,".regenie.gz",sep="")
      } else if (array[jj]=="decode") {
        file_tmp11<-paste(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/decode/",ph,"/allchr_decode_",ph,"_30062021_edited_noMult.regenie.gz",sep="")
      } else if (array[jj]=="danish") {
        file_tmp11<-paste(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/danish/",ph,"/allchr_danish_gsa_",ph,"_eur_sex_10PCs_saige_spa.regenie.gz",sep="")
      } else if (array[jj]=="eas") {
        file_tmp11<-paste(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/eas_iibdgc/",ph,"/ibd_EAS_SiKJ_meta_",ph,".regenie.gz",sep="")
      }

      if(!file.exists(file_tmp11)) next

      tmp11<-as.data.frame(fread(file_tmp11))
      tmp11<-tmp11[which(tmp11$ID %in% high_het$MarkerName),c("ID","BETA","SE","LOG10P")]

      if (array[jj]=="decode") {
        tmp11$pvalue<-10^(tmp11$LOG10P)
      } else {
        tmp11$pvalue<-10^(-tmp11$LOG10P)
      }
      tmp11<-tmp11[,c("ID","BETA","SE","pvalue")]
      colnames(tmp11)<-c("MarkerName","beta","se","pvalue")
      colnames(tmp11)[2:ncol(tmp11)]<-paste(colnames(tmp11)[2:ncol(tmp11)],ph,array[jj],sep="_")

      if(is.null(tmp1)) {
        tmp1<-tmp11
      } else {
        tmp1<-merge(tmp11,tmp1,by="MarkerName",all=T)
      }

      rm(tmp11)
    }


    if(is.null(tmp1)) next

    if(is.null(tmp)) {
      tmp<-tmp1
    } else {
      all_cols <- union(names(tmp), names(tmp1))
      tmp[setdiff(all_cols, names(tmp))] <- NA
      tmp1[setdiff(all_cols, names(tmp1))] <- NA
      tmp<-rbind(tmp,tmp1)
    }
    rm(tmp1)

  }
  
  assign(ph,tmp)

}

cd<-cd[!duplicated(cd),]
uc<-uc[!duplicated(uc),]
ibd<-ibd[!duplicated(ibd),]

dat<-merge(cd,uc,by="MarkerName")
dat<-merge(dat,ibd,by="MarkerName")

vec_cd<-colnames(dat)[grep("_cd_",colnames(dat))]
vec_cd<-c("MarkerName",vec_cd)

cd<-dat[,vec_cd]
cd<-cd[!duplicated(cd),]



dat<-dat[!duplicated(dat),]

fwrite(dat,"/path/to/user/git/IIBDGC_GWAS/plots/tmp_figures/index_variants_high_het_multiancestry_summary_stats.txt",col.names=T,row.names=F,quote=F,sep="\t")
fwrite(high_het,"/path/to/user/git/IIBDGC_GWAS/plots/tmp_figures/index_variants_high_het_multiancestry.txt",col.names=T,row.names=F,quote=F,sep="\t")

################################################################################################################################################################################

# inspect variants that have their significance dropped in the multiancestry analyses

dim(all[which(all$pvalue_eur_eas_sas>all$pvalue_eur),])
# [1]195 222


tmp<-all[which(all$pvalue_eur_eas_sas>all$pvalue_eur),
c("MarkerName","phenotype",
"bJ_ibd","bJ_se_ibd","pJ_ibd","bJ_cd","bJ_se_cd","pJ_cd","bJ_uc","bJ_se_uc","pJ_uc","Effect_eur","pvalue_eur","Effect_eur_eas_sas","pvalue_eur_eas_sas",
"BETA_ibd_eur_tier_2","P-value_ibd_eur_tier_2","BETA_cd_eur_tier_2","P-value_cd_eur_tier_2","BETA_uc_eur_tier_2","P-value_uc_eur_tier_2",
"BETA_ibd_eur_tier_2_eas_sas","P-value_ibd_eur_tier_2_eas_sas","BETA_cd_eur_tier_2_eas_sas","P-value_cd_eur_tier_2_eas_sas","BETA_uc_eur_tier_2_eas_sas","P-value_uc_eur_tier_2_eas_sas")]

fwrite(tmp,"/path/to/user/git/IIBDGC_GWAS/plots/tmp_figures/index_variants_lower_pvalue_multiancestry.txt",col.names=T,row.names=F,quote=F,sep="\t")

q("no")


################################################################################################################################################################################
# PART 2: LD between index variants and exonic variants
# (Run as a separate session, or remove the q("no") above)
#
# Setup (run in shell before starting R):
# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM=2000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -q cpu-interactive -G your_hpc_group R

library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr)

rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

# LOAD THE INDEPENDENT SIGNALS - see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_38_2_define_independent_signals_using_as_ref_cojo_agnostic_results.R
all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants_join_conditinal_model.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 631 182


ld<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/final_tables/ukbb_wgs_5mb.ld")
ld<-ld[which(ld$SNP_A %in% all$MarkerName | ld$SNP_B %in% all$MarkerName),]
ld<-ld[which(! (ld$SNP_A %in% all$MarkerName & ld$SNP_B %in% all$MarkerName)),]
dim(ld)
# [1] 1107    7



tmp<-all[,c("MarkerName","updated_region")]
tmp$best_exonic_variant_ld_ukb<-NA
tmp$best_exonic_variant_ld_r2_ukb<-NA


for (i in 1:nrow(tmp)) {

  ld_tmp<-ld[which(ld$SNP_A==tmp$MarkerName[i] | ld$SNP_B==tmp$MarkerName[i]),]

  
  if(nrow(ld_tmp)>0) {

    ld_tmp<-ld_tmp[which(ld_tmp$R2==max(ld_tmp$R2)),]

    vars<-c(ld_tmp$SNP_A,ld_tmp$SNP_B)
    vars<-vars[which(!vars %in% tmp$MarkerName)]

    tmp$best_exonic_variant_ld_ukb[i]<-paste(vars,collapse="|")
    tmp$best_exonic_variant_ld_r2_ukb[i]<-ld_tmp$R2[1]

  }

}


dim(tmp[which(!is.na(tmp$best_exonic_variant_ld_ukb)),])
# [1] 298   4

exo<-all[,c("MarkerName","exonic_variant_in_ld","exonic_variant_in_ld_gene_aac")]

tmp<-merge(tmp,exo,by="MarkerName",all=T)

fwrite(tmp,
paste0(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants_join_conditinal_model_UKB_LD_exome.tsv.gz"),
col.names=T,row.names=F,quote=F,sep="\t")