# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# retrieve a curated list of receptors and ligands, and explore how many of these pairs are in our list of nearby genes/genes in locus - later explore likely target genes

# database used as ref: https://mips.helmholtz-muenchen.de/corum/
# paper describing database: https://academic.oup.com/nar/article/53/D1/D651/7889246

# version downloaded: News (03.09.2024): The latest CORUM 5.0 release is now available, featuring 7,193 mammalian protein complexes derived from 5,299 unique genes. In this update, we have integrated drug targets from approved drugs,
# as provided by DrugCentral. Our current analysis includes 725 drug targets, identifying 1975 instances where drugs influence the formation or function of protein complexes.


# To run interactively on the cluster:
# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM=6000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -q cpu-interactive -G your_hpc_group R


library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(stringr)

rm(list = ls())
path_gwas <- "/path/to/ibdgwas/IIBDGC/"

# LOAD THE INDEPENDENT SIGNALS - see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_38_2_define_independent_signals_using_as_ref_cojo_agnostic_results.R
all <- fread(paste0(path_gwas, "post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz"))
all <- as.data.frame(all)
dim(all)
# [1]  619 192

# LOAD THE FILE WITH THE PROTEIN COMPLEXES:
pcom <- fread("/path/to/ibdgwas/IIBDGC/resources/corum/corum_humanComplexes.txt", head = T)
pcom <- as.data.frame(pcom)
dim(pcom)
# [1] 5513   35

# LOAD THE LIST OF EFFECTOR GENES:
list_genes <- fread(paste0(path_gwas, "post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence_with_coloc_mr_direction_with_ligand_receptor.tsv.gz"))
list_genes <- as.data.frame(list_genes)

# LINK THE LIST OF GENES TO THE PROTEIN COMPLEXES AND DRUG TARGETS:

# number of genes in the union of all annotation sets identifiable in the symbol namespace (background)
p_genes <- unique(unlist(str_split(pcom$subunits_gene_name, "\\;")))
p_genes <- p_genes[p_genes != ""]
Np <- length(p_genes)

# number of effector genes that are part of at least one annotation set:
mp <- sum(list_genes$list_genes %in% p_genes)

# build per-complex gene membership table via lapply, then combine
dat_list <- lapply(seq_len(nrow(pcom)), function(i) {

    all_genes <- unique(unlist(str_split(pcom$subunits_gene_name[i], "\\;")))
    ibd_genes     <- all_genes[all_genes %in% list_genes$list_genes]
    non_ibd_genes <- all_genes[!all_genes %in% list_genes$list_genes]

    if (length(ibd_genes) == 0) return(NULL)

    n_ibd  <- length(ibd_genes)
    n_all  <- length(all_genes)
    pval   <- phyper(n_ibd - 1, n_all, Np - n_all, mp, lower.tail = FALSE, log.p = FALSE)

    rows <- data.frame(
        list_genes           = c(ibd_genes, non_ibd_genes),
        ibd_gene             = c(rep(1L, n_ibd), rep(0L, length(non_ibd_genes))),
        n_ibd_genes_complex  = n_ibd,
        n_genes_complex      = n_all,
        protein_complex      = pcom$subunits_gene_name[i],
        protein_complex_drug = pcom$subunits_drugs[i],
        p_value              = pval,
        stringsAsFactors     = FALSE
    )
    rows
})

dat <- rbindlist(dat_list)
dim(dat)
# [1] 3792    7

# total number of complexes
dim(pcom)
# [1] 5513   34

## how many different protein complexes with at least 1 IBD gene:
length(unique(dat$protein_complex[dat$n_ibd_genes_complex > 0]))
# [1] 765

## how many different IBD proteins in any complex:
length(unique(dat$list_genes[dat$n_ibd_genes_complex > 0 & dat$ibd_gene == 1]))
# [1] 206

# how many complexes with at least 2 IBD genes:
length(unique(dat$protein_complex[dat$n_ibd_genes_complex > 1]))
# [1] 101

fwrite(dat, paste0(path_gwas, "post_imputation/2022/analysis/final_tables/protein_complexes_ibd_genes_corum.txt.gz"),
       col.names = T, row.names = F, quote = F, sep = "\t")


### complex_centric analyses:
prc <- unique(dat[, c("protein_complex", "n_ibd_genes_complex", "n_genes_complex", "protein_complex_drug", "p_value")])
dim(prc)
# [1] 765   5

# add complexes with no IBD genes
prc0 <- pcom[!pcom$subunits_gene_name %in% prc$protein_complex, c("subunits_gene_name", "subunits_drugs")]
colnames(prc0) <- c("protein_complex", "protein_complex_drug")
prc0$n_ibd_genes_complex <- 0L

# vectorised gene count and p-value for complexes with no IBD genes
prc0$n_genes_complex <- lengths(str_split(prc0$protein_complex, "\\;"))
prc0$p_value <- phyper(
    prc0$n_ibd_genes_complex - 1,
    prc0$n_genes_complex,
    Np - prc0$n_genes_complex,
    mp,
    lower.tail = FALSE,
    log.p = FALSE
)
dim(prc0)
# [1] 4743    5

prc <- rbind(prc, prc0)
dim(prc)
# [1] 5508    6

prc$q_value <- p.adjust(prc$p_value, method = "BH", n = nrow(prc))

summary(prc$q_value)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.3537  1.0000  1.0000  0.9938  1.0000  1.0000 

# NO SIGNIFICANT RESULTS AFTER MULTIPLE TEST CORRECTION

# dim(prc[prc$p_value < 0.05 & prc$n_ibd_genes_complex > 1, ])
# [1] 83  6

# dim(prc[prc$p_value < 0.05 & prc$n_ibd_genes_complex > 1 & prc$n_genes_complex > 2, ])
# [1] 73  6

dim(prc[prc$n_ibd_genes_complex > 1 & prc$n_genes_complex > 2, ])
# [1] 91  6

table(prc$n_ibd_genes_complex)
#    0    1    2    3    4
# 4743  664   96    4    1

fwrite(prc, paste0(path_gwas, "post_imputation/2022/analysis/final_tables/protein_complexes_corum.txt.gz"),
       col.names = T, row.names = F, quote = F, sep = "\t")
fwrite(prc, "~/git/IIBDGC_GWAS/plots/paper_tables/Source_data_human_protein_complexes_corum.txt.gz",
       col.names = T, row.names = F, quote = F, sep = "\t")

# add protein complex membership flag to the effector genes table:
list_genes$protein_complex_gene <- as.integer(list_genes$list_genes %in% dat$list_genes[dat$ibd_gene == 1])

table(list_genes$protein_complex_gene)
#   0   1
# 458 206

fwrite(list_genes, paste0(path_gwas, "post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence_with_coloc_mr_direction_with_ligand_receptor_with_protein_complexes.tsv.gz"),
       col.names = T, row.names = F, quote = F, sep = "\t")

q("no")

##################

# plot these and /path/to/user/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_58_define_protein_complexes.R outputs:
# ~/git/IIBDGC_GWAS/scripts/other/plot_ligand_receptor_pairs_and_protein_complexes.R
