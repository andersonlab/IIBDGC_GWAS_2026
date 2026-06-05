# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# Run before starting R (in shell, not here):
# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM=6000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q cpu-interactive R

library(rcartocolor)
library(qvalue)
library(data.table)
library(ggplot2)
library(colorspace)
library(dplyr)
library(tibble)
library(ggrepel)
library(ggpubr)
library(tidyverse)
library(scales)
library(png)
library(grid)

rm(list = ls())

path_gwas <- "/path/to/ibdgwas/IIBDGC/"

all <- fread(paste0(path_gwas, "post_imputation/2022/analysis/final_tables/",
  "list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_",
  "with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz"))
all <- as.data.frame(all)

chr  <- sort(unique(na.omit(all$chr)))
pheno <- c("ibd", "cd", "uc")

# ── Load per-chromosome frequency files and merge ──────────────────────────────
for (i in seq_along(pheno)) {

  tmp_list <- lapply(chr, function(j) {
    tmp1 <- fread(paste0(path_gwas,
      "post_imputation/2022/analysis/metaanalysis/", pheno[i], "/",
      j, "_", pheno[i],
      "_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_",
      "with_per_sample_rate_and_avgA2freq_withNeff.txt.gz"),
      header = TRUE)
    tmp1[tmp1$MarkerName %in% all$MarkerName,
         c("MarkerName", "avgA2FREQ", "avgA2FREQ_CASES", "avgA2FREQ_CONTROLS")]
  })
  tmp <- rbindlist(tmp_list)

  suffix <- paste0("_", pheno[i], "_eur_tier_2")
  setnames(tmp, 2:4, paste0(names(tmp)[2:4], suffix))

  if (i == 1) {
    dat <- tmp
  } else {
    dat <- merge(dat, tmp, by = "MarkerName", sort = FALSE)
  }
}

all <- merge(all, dat, by = "MarkerName", all.x = TRUE)

# MAF from IBD controls; pmin handles > 0.5 frequencies
all$MAF <- pmin(all$avgA2FREQ_CONTROLS_ibd_eur_tier_2,
                1 - all$avgA2FREQ_CONTROLS_ibd_eur_tier_2)

# ── Classify signals ───────────────────────────────────────────────────────────
known_classes <- c(
  "new_cojo_supervised_gw_significant_multiancestry_known_signal",
  "new_cojo_unsupervised_known_signal",
  "new_cojo_unsupervised_known_signal_known_exonic_variant",
  "new_cojo_unsupervised_new_signal_known_exonic_variant"
)
new_classes <- c(
  "new_cojo_supervised_gw_significant_multiancestry_new_signal",
  "new_cojo_supervised_gw_significant_multiancestry_new_signal_new_exonic_variant",
  "new_cojo_unsupervised_new_signal",
  "new_cojo_unsupervised_new_signal_new_exonic_variant"
)

all$class_region <- NA_character_
all$class_region[all$class_signal_final_exome %in% known_classes] <- "Known signals"
all$class_region[all$class_signal_final_exome %in% new_classes]   <- "New signals"

# ── Assign effect/SE/p from the most-relevant phenotype ───────────────────────
all$Effect <- NA_real_
all$StdErr <- NA_real_
all$pvalue <- NA_real_

assign_effect <- function(dat, rows, beta_col, se_col, p_col) {
  dat$Effect[rows] <- dat[[beta_col]][rows]
  dat$StdErr[rows] <- dat[[se_col]][rows]
  dat$pvalue[rows] <- dat[[p_col]][rows]
  dat
}

all <- assign_effect(all, all$phenotype == "IBD_unsaturated",
                     "BETA_ibd_eur_tier_2", "SE_ibd_eur_tier_2", "P-value_ibd_eur_tier_2")
all <- assign_effect(all, all$phenotype == "CD",
                     "BETA_cd_eur_tier_2",  "SE_cd_eur_tier_2",  "P-value_cd_eur_tier_2")
all <- assign_effect(all, all$phenotype == "UC",
                     "BETA_uc_eur_tier_2",  "SE_uc_eur_tier_2",  "P-value_uc_eur_tier_2")

# For IBD_saturated, use the sub-phenotype with the smaller p-value
ibd_sat <- all$phenotype == "IBD_saturated"
cd_better <- ibd_sat & (all$`P-value_cd_eur_tier_2` < all$`P-value_uc_eur_tier_2`)
uc_better <- ibd_sat & (all$`P-value_cd_eur_tier_2` > all$`P-value_uc_eur_tier_2`)

all <- assign_effect(all, cd_better,
                     "BETA_cd_eur_tier_2", "SE_cd_eur_tier_2", "P-value_cd_eur_tier_2")
all <- assign_effect(all, uc_better,
                     "BETA_uc_eur_tier_2", "SE_uc_eur_tier_2", "P-value_uc_eur_tier_2")

# ── Build pheno_class labels ───────────────────────────────────────────────────
make_pheno_class <- function(dat, short_labels = TRUE) {
  x <- table(dat$class_region, dat$phenotype, useNA = "ifany")

  fmt <- function(prefix, n) sprintf("%s (N=%d)", prefix, n)

  if (short_labels) {
    new_cd  <- fmt("CD",              x["New signals",   "CD"])
    kno_cd  <- fmt("CD",              x["Known signals", "CD"])
    new_uc  <- fmt("UC",              x["New signals",   "UC"])
    kno_uc  <- fmt("UC",              x["Known signals", "UC"])
    new_ibs <- fmt("IBD saturated",   x["New signals",   "IBD_saturated"])
    kno_ibs <- fmt("IBD saturated",   x["Known signals", "IBD_saturated"])
    new_ibu <- fmt("IBD unsaturated", x["New signals",   "IBD_unsaturated"])
    kno_ibu <- fmt("IBD unsaturated", x["Known signals", "IBD_unsaturated"])
  } else {
    new_cd  <- fmt("New CD signal",              x["New signals",   "CD"])
    kno_cd  <- fmt("Known CD signal",            x["Known signals", "CD"])
    new_uc  <- fmt("New UC signal",              x["New signals",   "UC"])
    kno_uc  <- fmt("Known UC signal",            x["Known signals", "UC"])
    new_ibs <- fmt("New IBD saturated signal",   x["New signals",   "IBD_saturated"])
    kno_ibs <- fmt("Known IBD saturated signal", x["Known signals", "IBD_saturated"])
    new_ibu <- fmt("New IBD unsaturated signal", x["New signals",   "IBD_unsaturated"])
    kno_ibu <- fmt("Known IBD unsaturated signal", x["Known signals", "IBD_unsaturated"])
  }

  lvls <- c(new_cd, kno_cd, new_uc, kno_uc, new_ibs, kno_ibs, new_ibu, kno_ibu)

  dat$pheno_class <- NA_character_
  dat$pheno_class[dat$class_region == "New signals"   & dat$phenotype == "CD"]             <- new_cd
  dat$pheno_class[dat$class_region == "Known signals" & dat$phenotype == "CD"]             <- kno_cd
  dat$pheno_class[dat$class_region == "New signals"   & dat$phenotype == "UC"]             <- new_uc
  dat$pheno_class[dat$class_region == "Known signals" & dat$phenotype == "UC"]             <- kno_uc
  dat$pheno_class[dat$class_region == "New signals"   & dat$phenotype == "IBD_saturated"]  <- new_ibs
  dat$pheno_class[dat$class_region == "Known signals" & dat$phenotype == "IBD_saturated"]  <- kno_ibs
  dat$pheno_class[dat$class_region == "New signals"   & dat$phenotype == "IBD_unsaturated"] <- new_ibu
  dat$pheno_class[dat$class_region == "Known signals" & dat$phenotype == "IBD_unsaturated"] <- kno_ibu

  dat$pheno_class <- factor(dat$pheno_class, levels = lvls)
  dat
}

all <- make_pheno_class(all, short_labels = FALSE)  # used only for the saved TSV

all$OR          <- exp(all$Effect)
all$confint_low <- exp(all$Effect - 1.96 * all$StdErr)
all$confint_up  <- exp(all$Effect + 1.96 * all$StdErr)

all <- all[order(all$MAF, decreasing = TRUE), ]
all$gene <- NA_character_

panel_name <- "Figure_1_independent_signals_effect_phenotype"
fwrite(all, paste0("~/git/IIBDGC_GWAS/plots/paper_figures/", panel_name, ".tsv.gz"),
       col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# ── Reload and rebuild plotting labels ────────────────────────────────────────
all <- fread(paste0("~/git/IIBDGC_GWAS/plots/paper_figures/", panel_name, ".tsv.gz"))
all <- as.data.frame(all)
all <- make_pheno_class(all, short_labels = TRUE)

# ── Colour palette ────────────────────────────────────────────────────────────
cols <- c('#004488', '#BB5566', '#DDAA33', "#db7107")

# ── y-axis limits (shared across both panels) ─────────────────────────────────
ymax <- ceiling(max(abs(log(all$confint_low)), abs(log(all$confint_up)), na.rm = TRUE) + 0.2)
ymin <- -ymax

maf_breaks <- as.numeric(c("0.001","0.002","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))

# ── Shared ggplot theme & axes ─────────────────────────────────────────────────
base_theme <- theme(
  text             = element_text(size = 14),
  legend.position  = "right",
  legend.title     = element_blank(),
  legend.key       = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line        = element_line(colour = "darkgrey"),
  axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1),
  plot.title       = element_text(face = "bold", size = 12),
  plot.margin      = margin(t = 0.5, b = 0.2, r = 0.2, l = 0.6, unit = "cm")
)

make_signal_plot <- function(dat, title) {
  ggplot(dat, aes(x = -log10(MAF), y = log(OR), color = pheno_class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
    geom_errorbar(aes(ymin = log(confint_low), ymax = log(confint_up))) +
    geom_point() +
    geom_label_repel(
      data = filter(dat, !is.na(gene) & gene != ""),
      aes(label = gene),
      force = 20, min.segment.length = 0, seed = 42,
      nudge_x = 0.07, nudge_y = 0.07,
      show.legend = FALSE, fontface = "italic"
    ) +
    scale_color_manual(values = cols) +
    scale_x_continuous(
      trans   = "reverse",
      limits  = c(-log10(0.001), -log10(0.6)),
      breaks  = -log10(maf_breaks),
      labels  = as.character(maf_breaks)
    ) +
    scale_y_continuous(limits = c(ymin, ymax)) +
    labs(y = "log(Odds Ratio)", x = "Minor Allele Frequency", title = title) +
    base_theme
}

# ── Known signals panel ────────────────────────────────────────────────────────
known <- all[all$class_region == "Known signals", ]
known$gene[known$exonic_variant_in_ld_gene_aac %in%
             c("IL10RA:P295L", "NOD2:V793M", "CARD9:.", "ADCY7:D439E", "IL23R:G149R")] <-
  known$exonic_variant_in_ld_gene_aac[known$exonic_variant_in_ld_gene_aac %in%
             c("IL10RA:P295L", "NOD2:V793M", "CARD9:.", "ADCY7:D439E", "IL23R:G149R")]

plot_signals_1 <- make_signal_plot(known,
  sprintf("Known Independent signals, N = %d", nrow(known)))

# ── New signals panel ──────────────────────────────────────────────────────────
new <- all[all$class_region == "New signals", ]
new$gene[new$MarkerName == "chr18:49639885:T:C"] <- "LIPG"
new$gene[new$MarkerName == "chr19:43648948:A:G"] <- "PLAUR"
new$gene[new$MarkerName == "chr7:592301:G:A"]    <- "PDGFA"
new$gene[new$MarkerName == "chr3:122386581:A:G"] <- "KPNA1"
new$gene[new$exonic_variant_in_ld_gene_aac %in%
           c("PTCD1:R113W", "PTGER2:C83G", "JAK2:R1063H")] <-
  new$exonic_variant_in_ld_gene_aac[new$exonic_variant_in_ld_gene_aac %in%
           c("PTCD1:R113W", "PTGER2:C83G", "JAK2:R1063H")]

plot_signals_2 <- make_signal_plot(new,
  sprintf("New Independent signals, N = %d", nrow(new)))

plot_signals <- ggarrange(plot_signals_1, plot_signals_2, nrow = 2, labels = c("b", "c"))

# ── Sample composition bar chart ───────────────────────────────────────────────
df_totals <- tribble(
  ~ancestry, ~IBD_total, ~Controls, ~CD,    ~UC,
  "EUR",      110392,     1239781,   51983,  50070,
  "SAS",        1207,        6992,     379,    633,
  "EAS",       14393,       15456,    7372,   6862
) %>%
  mutate(
    ancestry = factor(ancestry, levels = c("EUR", "EAS", "SAS")),
    IBDu     = IBD_total - (CD + UC)
  )

df_bar <- df_totals %>%
  select(ancestry, CD, UC, IBDu) %>%
  pivot_longer(cols = c(CD, UC, IBDu), names_to = "subtype", values_to = "n") %>%
  mutate(group = "Cases") %>%
  bind_rows(
    df_totals %>% transmute(ancestry, subtype = "Controls", n = Controls, group = "Controls")
  ) %>%
  group_by(group, ancestry) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(subtype = factor(subtype, levels = c("CD", "UC", "IBDu", "Controls")))

facet_totals <- df_bar %>%
  group_by(group) %>%
  summarise(N = sum(n), .groups = "drop") %>%
  mutate(group_label = paste0(group, " (N = ", comma(N), ")"))

df_bar <- df_bar %>%
  left_join(facet_totals, by = "group") %>%
  mutate(label_n = ifelse(prop >= 0.06 & n > 0, comma(n), NA_character_))

p_bar <- ggplot(df_bar, aes(x = ancestry, y = prop, fill = subtype)) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = label_n),
    color    = "white",
    position = position_stack(vjust = 0.5),
    size     = 3.3,
    na.rm    = TRUE
  ) +
  facet_grid(group_label ~ ., scales = "free_y") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(values = c(
    "CD"       = "#004488",
    "UC"       = "#BB5566",
    "IBDu"     = "#DDAA33",
    "Controls" = "#1b9e77"
  )) +
  labs(y = "Proportion of individuals per population ancestry", x = "Population ancestry", fill = NULL) +
  theme(
    text             = element_text(size = 14),
    legend.title     = element_blank(),
    legend.position  = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.y     = element_text(size = 14),
    axis.text        = element_text(color = "black"),
    axis.line        = element_line(colour = "darkgrey"),
    strip.text       = element_text(size = 14),
    plot.margin      = margin(5, 5, 5, 5)
  )

# ── Combine and save ───────────────────────────────────────────────────────────
p <- ggarrange(p_bar, plot_signals, widths = c(3, 7.5), labels = c("a", ""))

ggsave(
  paste0("~/git/IIBDGC_GWAS/plots/paper_figures/Figure_1_version_5.pdf"),
  p,
  width     = 170,
  height    = 90,
  dpi       = 600,
  units     = "mm",
  limitsize = TRUE,
  scale     = 2
)

q("no")
