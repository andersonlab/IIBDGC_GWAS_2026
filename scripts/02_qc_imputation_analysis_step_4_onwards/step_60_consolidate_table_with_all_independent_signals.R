# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# To run interactively:
  # singularity exec iibdgc_postprocess_10_singularity.sif
  MEM=8000
  bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -q cpu-interactive -G your_hpc_group R

###############################################################

# 1.-  COMBINE THE RESULTS

library(data.table)

rm(list = ls())

path_gwas <- "/path/to/ibdgwas/IIBDGC/"

# LOAD THE INDEPENDENT SIGNALS
# see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_38_2_define_independent_signals_using_as_ref_cojo_agnostic_results.R
all <- fread(paste0(path_gwas, "post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz"))
all <- as.data.frame(all)
dim(all)
# [1] 619 192

cols_to_drop <- c("exclude", "var_id", "deviance_model_0", "deviance_model_1", "deviance_model_2",
  "deviance_model_3", "deviance_model_4", "pvalue_m1_m4", "pvalue_m2_m4", "pvalue_m3_m4",
  "pvalue_m1_m0", "pvalue_m2_m0", "pvalue_m3_m0", "model_1", "class", "class_signal",
  "class_signal_final", "variants_in_ld_r2_0.1")

all <- all[, !colnames(all) %in% cols_to_drop]


##############################################
# Helper: load allele frequencies across chromosomes and phenotypes

load_freq_by_ancestry <- function(all, path_gwas, ancestry, file_pattern, suffix) {
  pheno <- c("ibd", "cd", "uc")
  chrs <- sort(unique(all$chr[!is.na(all$chr)]))

  dat <- NULL
  for (i in seq_along(pheno)) {
    print(pheno[i])

    tmp <- rbindlist(lapply(chrs, function(j) {
      f <- sprintf("%spost_imputation/2022/analysis/metaanalysis/%s/%d_%s_%s",
        path_gwas, pheno[i], j, pheno[i], file_pattern)
      tmp1 <- fread(f, header = TRUE)
      tmp1 <- tmp1[tmp1$MarkerName %in% all$MarkerName,
                   .(MarkerName, A2FREQ = avgA2FREQ, A2FREQ_CASES = avgA2FREQ_CASES, A2FREQ_CONTROLS = avgA2FREQ_CONTROLS)]
      tmp1
    }))

    colnames(tmp)[2:4] <- paste(colnames(tmp)[2:4], pheno[i], suffix, sep = "_")

    if (i == 1) {
      dat <- tmp
    } else {
      dat <- merge(dat, tmp, by = "MarkerName", sort = FALSE)
    }
  }
  dat
}

# SAS tier2
dat_sas <- load_freq_by_ancestry(all, path_gwas,
  ancestry = "sas",
  file_pattern = "meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz",
  suffix = "sas_tier_2")
all <- merge(all, dat_sas, by = "MarkerName", all.x = TRUE)

# EUR tier2
dat_eur <- load_freq_by_ancestry(all, path_gwas,
  ancestry = "eur",
  file_pattern = "meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz",
  suffix = "eur_tier_2")
all <- merge(all, dat_eur, by = "MarkerName", all.x = TRUE)


##############################################
### Add EAS summary statistics

pheno <- c("ibd", "cd", "uc")

eas <- NULL
for (i in seq_along(pheno)) {
  print(pheno[i])

  tmp <- fread(paste0(path_gwas, "post_imputation/2022/analysis/stage_2_summary_statistics/eas_iibdgc/",
    pheno[i], "/ibd_EAS_SiKJ_meta_", pheno[i], ".regenie.gz"), header = TRUE)
  tmp <- tmp[tmp$ID %in% all$MarkerName, ]

  # LOG10P is -log10(p), so convert back to p-value correctly
  tmp$pvalue <- 10^(-tmp$LOG10P)

  if (i == 1) {
    tmp <- tmp[, .(MarkerName = ID, Position_b38 = GENPOS, A1 = ALLELE0, A2 = ALLELE1,
                   A2FREQ = A1FREQ, BETA, SE, `P-value` = pvalue)]
    colnames(tmp)[2:4] <- paste(colnames(tmp)[2:4], "eas", sep = "_")
    colnames(tmp)[5:ncol(tmp)] <- paste(colnames(tmp)[5:ncol(tmp)], pheno[i], "eas", sep = "_")
    eas <- tmp
  } else {
    tmp <- tmp[, .(MarkerName = ID, BETA, SE, `P-value` = pvalue)]
    colnames(tmp)[2:ncol(tmp)] <- paste(colnames(tmp)[2:ncol(tmp)], pheno[i], "eas", sep = "_")
    eas <- merge(eas, tmp, by = "MarkerName", sort = FALSE, all = TRUE)
  }
}

dim(eas)

all <- merge(all, eas, by = "MarkerName", all.x = TRUE)


fwrite(all,
  paste0(path_gwas, "post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants_join_conditinal_model.tsv.gz"),
  col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

q("no")

