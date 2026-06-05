# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
####################################
## DRAW A LOCUS ZOOM PLOT (locuszoomr + plink LD)

# # singularity exec iibdgc_postprocess_10_singularity.sif

# Usage:
# Rscript draw_locus_zoom_plot.R \
#   <pheno> \
#   <release> \
#   <region> \
#   [index_snp] \
#   [output_pdf] \
#   [study] \
#   [condition] \
#   [gene_id] \
#   [qtl_path]
#
# Arguments:
#   pheno      : phenotype, one of: ibd, cd, uc
#   release    : release label e.g. eur_eas_sas_tier_2
#   region     : region string CHR_START_END e.g. 1_52976954_54020750
#   index_snp  : (optional) index variant MarkerName e.g. chr1:12345:A:G
#                defaults to top SNP (lowest p-value)
#   output_pdf : (optional) output PDF path
#                defaults to locus_zoom_<pheno>_<region>.pdf in working dir
#   study      : (optional) study label, one of:
#                Yazar_2022, Hu_2021, Panousis_2024, eQTLGen_2021, blueprint, or eQTL_catalog_V6
#   condition  : (optional) condition_name to plot e.g. IFNB_6
#                required if study is provided
#   gene_id    : (optional) gene ID to subset nominal eQTL files (matched against gene_id column)
#                e.g. PDCD1; required if study is provided
#   qtl_path   : (optional) base path for nominal QTL files
#                required if study is provided

# for testing purposes:
# MEM=15000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R


library(data.table)
library(locuszoomr)
library(ggpubr)
library(qvalue)
library(GenomicRanges)
# library(Rsamtool)
library(dbplyr)
library(BiocFileCache)
library("ensembldb")
library(AnnotationHub)
library(EnsDb.Hsapiens.v86)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript draw_locus_zoom_plot.R <pheno> <release> <region> [index_snp] [output_pdf] [study] [condition] [gene_id] [qtl_path]")
}

pheno         <- args[1]
release       <- args[2]
region        <- args[3]
index_snp_arg <- if (length(args) >= 4 && nchar(args[4]) > 0) args[4] else NULL
output_pdf    <- if (length(args) >= 5 && nchar(args[5]) > 0) args[5] else NULL
study         <- if (length(args) >= 6 && nchar(args[6]) > 0) args[6] else NULL
condition     <- if (length(args) >= 7 && nchar(args[7]) > 0) args[7] else NULL
gene_id       <- if (length(args) >= 8 && nchar(args[8]) > 0) args[8] else NULL
path          <- if (length(args) >= 9 && nchar(args[9]) > 0) args[9] else NULL


##########
# # test:

# pheno         <- "ibd"
# release       <- "eur_tier_2"
# region        <- "2_201853808_204267954"
# index_snp_arg <- "chr2:203727298:G:A"
# study         <- "Schmiedel_2018"
# condition     <- "Treg_naive"
# gene_id       <- "CD28"
# path          <- "/path/to/project"
# output_pdf    <- paste0("~/git/IIBDGC_GWAS/plots/regional_manhattan_plots/locus_zoom_",
#                         pheno, "_", region, "_", gene_id, "_", study, "_", condition, ".pdf")


# path_gwas<-"/path/to/ibdgwas/IIBDGC/"
# index_snp_arg<-"chr11:85454082:A:G"
# region<-"11_84823557_86153989"
# pheno<-"ibd"
# release<-"eur_tier_2"

# path_gwas<-"/path/to/ibdgwas/IIBDGC/"
# index_snp_arg<-"chr7:5433979:G:C"
# region<-"7_3901126_7438156"
# pheno<-"ibd"
# release<-"eur_tier_2"

# path_gwas<-"/path/to/ibdgwas/IIBDGC/"
# index_snp_arg<-"chr1:53476954:G:A"
# region<-"1_52976954_54020750"
# pheno<-"cd"
# release<-"eur_tier_2"



# Build default output PDF name (includes gene/study/condition when provided)
if (is.null(output_pdf)) {
  qtl_suffix <- if (!is.null(gene_id) && !is.null(study) && !is.null(condition))
    paste0("_", gene_id, "_", study, "_", condition) else ""
  output_pdf <- paste0("~/git/IIBDGC_GWAS/plots/regional_manhattan_plots/locus_zoom_",
                       pheno, "_", region, qtl_suffix, ".pdf")
}

print(pheno)
print(release)
print(region)
print(index_snp_arg)

path_gwas <- "/path/to/ibdgwas/IIBDGC/"

# --- Parse region ---
chr     <- gsub("_.*", "", region)
tmp     <- gsub("^[0-9]{1,2}_", "", region)
max_pos <- as.numeric(gsub("^[0-9]*_", "", tmp))
min_pos <- as.numeric(gsub("_[0-9]*$", "", tmp))

cat("Region:", region, "\n")
cat("Chr:", chr, " Start:", min_pos, " End:", max_pos, "\n")

plink_prefix <- paste0(path_gwas, "post_imputation/2022/analysis/conditional_analysis/eur/",
                       "allarrays_chr", chr, "_subset_included_in_", pheno, "_analysis")
cat("Plink prefix:", plink_prefix, "\n")

# --- Read and filter summary stats ---
sumstats_file <- paste0(path_gwas, "post_imputation/2022/analysis/metaanalysis/", pheno, "/",
                        chr, "_", pheno, "_meta_", release,
                        "_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz")
cat("Reading summary stats:", sumstats_file, "\n")
risk <- fread(sumstats_file, header = TRUE)

# Extract position from MarkerName (format: chr1:12345:A:G)
risk$position <- gsub("chr[0-9]{1,2}:", "", risk$MarkerName)
risk$position <- as.numeric(gsub(":[A-Z]*:[A-Z]*$", "", risk$position))
risk$pval     <- as.numeric(risk[["P-value"]])
risk$chr      <- as.numeric(chr)
risk$variant  <- risk$MarkerName

rate_Neff_threshold <- if (release == "eur_tier_1") 0.8 else 0.5
risk <- risk[!is.na(risk$position) & !is.na(risk$pval) & risk$rate_Neff > rate_Neff_threshold, ]
risk <- risk[risk$position >= min_pos & risk$position <= max_pos, ]
# risk <- risk[which(risk$pval>1E-50),]

if (nrow(risk) == 0) stop("No variants found in region after filtering.")
cat("Variants in region:", nrow(risk), "\n")

# --- Identify index SNP ---
if (!is.null(index_snp_arg)) {
  if (!index_snp_arg %in% risk$variant) {
    warning("Specified index_snp not found in region; falling back to top SNP.")
    index_snp_arg <- NULL
  }
}
if (is.null(index_snp_arg)) {
  index_snp_arg <- risk$variant[which.min(risk$pval)]
}
cat("Index SNP:", index_snp_arg, "\n")

index_pos <- risk$position[risk$variant == index_snp_arg][1]

# Plot window: 200kb each side of index SNP (LD region stays as defined by region arg)
plot_min <- index_pos - 200000
plot_max <- index_pos + 200000

# --- Compute LD with plink ---
ld_dir    <- paste0(path_gwas, "post_imputation/2022/analysis/ld_data/eur_tier2/locus_zoom/")
ld_prefix <- paste0(ld_dir, pheno, "_", region)
ld_file   <- paste0(ld_prefix, ".ld")

if (file.exists(ld_file)) {
  cat("LD file already exists, skipping plink:", ld_file, "\n")
} else {
  bim <- fread(paste0(plink_prefix, ".bim"), header = FALSE,
               col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  bim <- bim[bim$CHR == as.numeric(chr), ]

  matched_bim <- bim[bim$BP == index_pos, ]
  if (nrow(matched_bim) == 0) {
    nearest_idx <- which.min(abs(bim$BP - index_pos))
    matched_bim <- bim[nearest_idx, ]
    warning(paste("Index SNP not found by position in bim; using nearest:", matched_bim$SNP))
  }
  plink_index_snp <- matched_bim$SNP[1]
  cat("Plink index SNP ID:", plink_index_snp, "\n")

  region_kb <- ceiling((max_pos - min_pos) / 1000) + 1

  plink_cmd <- paste(
    "plink",
    "--bfile", plink_prefix,
    "--chr", chr,
    "--from-bp", min_pos,
    "--to-bp", max_pos,
    "--r2",
    "--ld-snp", plink_index_snp,
    "--ld-window-kb", region_kb,
    "--ld-window 999999",
    "--ld-window-r2 0",
    "--out", ld_prefix,
    "--silent"
  )

  cat("Running plink for LD...\n")
  exit_code <- system(plink_cmd)
  if (exit_code != 0) stop("plink failed. Check plink_prefix and that plink is on PATH.")
  if (!file.exists(ld_file)) stop("plink LD output not found: ", ld_file)
}

ld       <- fread(ld_file, header = TRUE)
ld$BP_B  <- as.integer(trimws(ld$BP_B))
ld$R2    <- as.numeric(trimws(ld$R2))
ld$SNP_B <- trimws(ld$SNP_B)
cat("LD file:", nrow(ld), "pairs, R2 range:",
    round(min(ld$R2, na.rm=T), 3), "-", round(max(ld$R2, na.rm=T), 3), "\n")

# Build position-keyed lookup (most reliable since positions are unambiguous integers)
ld_by_pos <- setNames(ld$R2, as.character(ld$BP_B))

# --- Load Ensembl DB for gene track ---
ens_db <- tryCatch({
  library(EnsDb.Hsapiens.v86)
  EnsDb.Hsapiens.v86
}, error = function(e) {
  cat("EnsDb.Hsapiens.v86 not available; trying AnnotationHub...\n")
  library(AnnotationHub)
  ah <- AnnotationHub()
  q <- query(ah, c("EnsDb", "Homo sapiens", "GRCh38"))
  ah[[tail(names(q), 1)]]
})

# --- Create GWAS locus object ---
# Add r2 directly as a column so locuszoomr uses it without any name-matching
sumstats_df      <- as.data.frame(risk[, c("chr", "position", "pval", "variant")])
sumstats_df$r2   <- as.numeric(ld_by_pos[as.character(sumstats_df$position)])
sumstats_df$r2[sumstats_df$variant == index_snp_arg] <- 1
cat("LD matched for", sum(!is.na(sumstats_df$r2)), "/", nrow(sumstats_df),
    "variants (range:", round(min(sumstats_df$r2, na.rm=T), 2),
    "-", round(max(sumstats_df$r2, na.rm=T), 2), ")\n")

loc_gwas <- locus(
  data      = sumstats_df,
  seqname   = as.numeric(chr),
  xrange    = c(plot_min, plot_max),
  chrom     = "chr",
  pos       = "position",
  p         = "pval",
  labs      = "variant",
  LD        = "r2",
  index_snp = index_snp_arg,
  ens_db    = ens_db
)


# Get gene names
gtf<-rtracklayer::import(paste(path_gwas,'post_imputation/2022/analysis/metaanalysis/annotation/gencode.v39.annotation.gtf.gz',sep=""))
gene<-as.data.frame(gtf)
gene<-gene[,c("gene_id","gene_name")]
gene$gene<-gsub("\\.[0-9]*","",gene$gene_id)
gene<-gene[!duplicated(gene$gene),]
rm(gtf)


# --- Load eQTL data (if study, condition, gene_id, and path provided) ---
eqtl_loci <- list()

if (!is.null(study) && !is.null(condition) && !is.null(gene_id) && !is.null(path)) {

  cat("Loading eQTL data: study =", study, "| condition =", condition, "| gene =", gene_id, "\n")

  # Build nominal file path from path + study + condition
  qtl_file <- if (study == "Yazar_2022") {
    paste0(path, "QTL_managed_access/yazar_2022/OneK1K_matrix_eQTL_results/nominal/", condition, "/", condition, ".tsv.gz")
  } else if (study == "Hu_2021") {
    paste0(path, "QTL_managed_access/hu_2021/nominal/", condition, "/", condition, ".tsv.gz")
  } else if (study == "Panousis_2024") {
    paste0(path, "macromap/nominal/", condition, "/", condition, ".tsv.gz")
  } else if (study == "eQTLGen_2021") {
    paste0(path, "eQTLGen/nominal/", condition, "/", condition, ".tsv.gz")
  } else if (study == "blueprint") {
    paste0(path, "blueprint/nominal/", condition, "/", condition, ".tsv.gz")
  } else {
    paste0(path, "eQTL_catalog_V6/nominal/", condition, "/", condition, ".tsv.gz")
  }

  cat("Reading:", qtl_file, "\n")
  param <- GRanges(chr, IRanges(min_pos, max_pos))
  tbx   <- Rsamtools::TabixFile(qtl_file)
  res   <- Rsamtools::scanTabix(tbx, param = param)
  eqtl  <- as.data.frame(Map(function(elt) read.csv(textConnection(elt), sep = "\t", header = FALSE), res)[[1]])
  colnames(eqtl) <- c("variant", "r2", "pvalue", "molecular_trait_object_id", "molecular_trait_id", "maf",
                       "gene_id", "median_tpm", "beta", "se", "an", "ac", "chromosome", "position",
                       "ref", "alt", "type", "rsid")

  eqtl$gene_id<-gsub("\\.[0-9]*","",eqtl$gene_id)
  eqtl<-merge(eqtl,gene,by.x="gene_id",by.y="gene",all.x=T)

  eqtl <- eqtl[which(eqtl$gene_name == gene_id), ]
  eqtl$MarkerName <- gsub("_", ":", eqtl$variant)
  eqtl$pvalue     <- as.numeric(eqtl$pvalue)
  eqtl$position   <- as.numeric(eqtl$position)
  eqtl$chr        <- as.numeric(chr)

  if (nrow(eqtl) >= 2) {
    eqtl_df <- as.data.frame(eqtl[, c("chr", "position", "pvalue", "MarkerName")])
    colnames(eqtl_df)[3] <- "pval"

    eqtl_df$r2 <- as.numeric(ld_by_pos[as.character(eqtl_df$position)])
    eqtl_df$r2[eqtl_df$MarkerName == index_snp_arg] <- 1

    eqtl_index_snp <- eqtl_df$MarkerName[which.min(eqtl_df$pval)]
    cat("eQTL lead SNP:", eqtl_index_snp, "\n")

    loc_eqtl <- tryCatch(
      locus(
        data      = eqtl_df,
        seqname   = as.numeric(chr),
        xrange    = c(plot_min, plot_max),
        chrom     = "chr",
        pos       = "position",
        p         = "pval",
        labs      = "MarkerName",
        LD        = "r2",
        index_snp = eqtl_index_snp,
        ens_db    = ens_db
      ),
      error = function(e) { warning("Failed to create eQTL locus: ", e$message); NULL }
    )

    if (!is.null(loc_eqtl)) {
      attr(loc_eqtl, "panel_title") <- paste0(gene_id, " | ", study, " | ", condition)
      eqtl_loci[[1]] <- loc_eqtl
    }
  } else {
    warning("Fewer than 2 eQTL variants in region for gene ", gene_id, " — skipping eQTL panel.")
  }
}

# --- Plot ---
cat("Saving plot to:", output_pdf, "\n")

n_qtl_panels <- length(eqtl_loci)

gwas_title <- paste(toupper(pheno), "|", index_snp_arg)
main_title <- if (!is.null(study) && !is.null(condition) && !is.null(gene_id)) {
  paste0(toupper(pheno), " | ", index_snp_arg, " | ", gene_id, " | ", study, " | ", condition)
} else {
  paste0(toupper(pheno), " | ", index_snp_arg)
}

if (n_qtl_panels == 0) {
  # GWAS only
  pdf(output_pdf, height = 7, width = 9)
  print(locus_ggplot(loc_gwas, labels = "index",
               highlight = gene_id, highlight_col = "red") + ggtitle(gwas_title))
  dev.off()
} else {
  # Multi-panel: GWAS scatter + eQTL scatters + shared gene track
  p_gwas_leg <- gg_scatter(loc_gwas, labels = "index") +
    labs(title = gwas_title) +
    guides(colour = guide_legend(title = expression(r^2))) +
    theme(legend.position = "right")

  # Extract shared r2 legend from GWAS panel, then strip legend from all scatter panels
  shared_legend <- get_legend(p_gwas_leg)

  p_gwas <- p_gwas_leg + theme(legend.position = "none")

  p_eqtl <- lapply(eqtl_loci, function(loc) {
    gg_scatter(loc, labels = "index") +
      labs(title = attr(loc, "panel_title")) +
      theme(legend.position = "none")
  })

  p_genes <- gg_genetracks(loc_gwas, highlight = gene_id, highlight_col = "red")

  # Stack: GWAS on top, then eQTL panels, gene track at bottom (all without legend)
  all_panels <- c(list(p_gwas), p_eqtl, list(p_genes))
  panel_heights <- c(rep(3, 1 + n_qtl_panels), 1.5)

  panels_stacked <- ggarrange(
    plotlist = all_panels,
    ncol     = 1,
    heights  = panel_heights,
    align    = "hv"
  )

  # Attach single legend to the right of the stacked panels
  panels_with_legend <- ggarrange(
    panels_stacked,
    as_ggplot(shared_legend),
    ncol   = 2,
    widths = c(10, 1)
  )

  total_height <- 2 + sum(panel_heights)
  pdf(output_pdf, height = total_height, width = 9)
  print(
    annotate_figure(
      panels_with_legend,
      top = text_grob(main_title, face = "bold", size = 12)
    )
  )
  dev.off()
}

cat("Done.\n")
q("no")
