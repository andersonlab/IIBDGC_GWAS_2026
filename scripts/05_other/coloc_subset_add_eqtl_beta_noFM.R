# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM=12000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q cpu-interactive R

library(data.table)
library(rtracklayer)
library(stringr)

path_gwas <- "/path/to/ibdgwas/IIBDGC/"
path_eqtl <- "/path/to/project"

EQTL_COLS <- c("variant", "r2", "pvalue", "molecular_trait_object_id", "molecular_trait_id",
                "maf", "gene_id", "median_tpm", "beta", "se", "an", "ac",
                "chromosome", "position", "ref", "alt", "type", "rsid")

# Resolve the path to the nominal tabix file for a given cohort/condition
get_file_name <- function(cohort, condition_name) {
    if (cohort == "IBDverse") {
        paste0("/path/to/project",
               "IBDverse/snakemake_coloc/input/2025_06_11_IBDverse_coloc_all_gwas/nominal/",
               condition_name, "/", condition_name, ".tsv.gz")
    } else if (cohort == "hu_2021") {
        paste0(path_eqtl, "QTL_managed_access/", cohort,
               "/nominal/", condition_name, "/", condition_name, ".tsv.gz")
    } else {
        paste0(path_eqtl, cohort, "/nominal/", condition_name, "/", condition_name, ".tsv.gz")
    }
}

# Build the X key used to match rows back to coloc_results
make_x_key <- function(condition_name, id_col, variant, cohort) {
    if (cohort == "pQTL_sparc") {
        paste0(condition_name, "_", id_col, "_chr", variant)
    } else {
        paste(condition_name, id_col, variant, sep = "_")
    }
}

# LOAD THE COLOCALIZATION RESULTS
# See ~/git/IIBDGC_GWAS/scripts/other/coloc_subset_results_by_ld_leads_noFM.R
coloc_results <- as.data.frame(fread(paste0(path_gwas,
    "post_imputation/2022/analysis/metaanalysis/coloc_results/",
    "coloc_filtered_by_gwas_r2_0.6_allpheno_noFM.tsv")))

# Unique set of gene and eQTL variants
index_col <- unique(coloc_results[, c("condition_name", "MarkerName_eqtl_lead", "cohort")])

index_col$chr <- as.numeric(gsub(":.*", "",
    gsub("^chr", "", index_col$MarkerName_eqtl_lead)))

index_col$position <- as.numeric(gsub("^chr[0-9]{1,2}:", "",
    gsub(":[A-Z]*:[A-Z]*$", "", index_col$MarkerName_eqtl_lead)))

# Iterate over unique eQTL leads and extract beta/SE/p-value from tabix files
results_list <- vector("list", nrow(index_col))

for (i in seq_len(nrow(index_col))) {

    file_name <- get_file_name(index_col$cohort[i], index_col$condition_name[i])

    if (!file.exists(file_name)) {
        warning(sprintf("File not found (row %d): %s", i, file_name))
        next
    }

    param <- GRanges(as.character(index_col$chr[i]),
                     IRanges((index_col$position[i] - 5):(index_col$position[i] + 5)))
    tbx  <- Rsamtools::TabixFile(file_name)
    res  <- Rsamtools::scanTabix(tbx, param = param)

    if (length(res[[1]]) == 0) next

    tmp <- read.csv(textConnection(res[[1]]), sep = "\t", header = FALSE,
                    col.names = EQTL_COLS)

    # Match on gene_id key
    tmp$X <- make_x_key(index_col$condition_name[i], tmp$gene_id, tmp$variant, index_col$cohort[i])
    matched_gene <- tmp[tmp$X %in% coloc_results$X, ]

    # Match on molecular_trait_object_id key
    tmp$X <- make_x_key(index_col$condition_name[i], tmp$molecular_trait_object_id, tmp$variant, index_col$cohort[i])
    matched_trait <- tmp[tmp$X %in% coloc_results$X, ]

    combined <- unique(rbind(matched_gene, matched_trait))
    results_list[[i]] <- combined

    if (i %% 100 == 0) message(sprintf("Processed %d / %d", i, nrow(index_col)))
}

dat <- do.call(rbind, results_list)

unmatched <- sum(!coloc_results$X %in% dat$X)
if (unmatched > 0) warning(sprintf("%d coloc_results rows have no eQTL beta match", unmatched))

coloc_results <- merge(
    coloc_results,
    dat[, c("X", "beta", "se", "pvalue")],
    by = "X",
    all.x = TRUE
)
colnames(coloc_results)[colnames(coloc_results) == "beta"]   <- "Beta_eqtl"
colnames(coloc_results)[colnames(coloc_results) == "se"]     <- "SE_eqtl"
colnames(coloc_results)[colnames(coloc_results) == "pvalue"] <- "Pvalue_eqtl"

write.table(coloc_results,
    paste0(path_gwas, "post_imputation/2022/analysis/metaanalysis/coloc_results/",
           "coloc_filtered_by_gwas_r2_0.6_allpheno_noFM_with_beta.tsv"),
    col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
