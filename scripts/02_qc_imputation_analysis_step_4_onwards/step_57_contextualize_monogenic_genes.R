# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# retrieve a curated list of receptors and ligands, and explore how many of these pairs are in our list of nearby genes/genes in locus - later explore likely target genes

# database used as ref: http://tcm.zju.edu.cn/celltalkdb/
# paper describing database: https://academic.oup.com/bib/article-abstract/22/4/bbaa269/5955941?redirectedFrom=fulltext&login=true#no-access-message

# version downloaded:
#
# File	Description	Statistics	Download times	Latest Download	License	Size
# human_lr_pair.txt Human ligand-receptor interaction pairs	3,399 lines, 10 columns.	2624 times	2022-10-04 18:24	GPL-3.0	355.15 Kb


# bash commands to run before opening R:
# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM=8000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -q cpu-interactive -G your_hpc_group R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(stringr)
library(rtracklayer)

set.seed(42)

path <- "/path/to/ibdgwas/IIBDGC/"


### probability of getting same N genes when randomly sampling from:
gtf  <- import(paste0(path, "post_imputation/2022/analysis/metaanalysis/annotation/gencode.v44.annotation.gtf.gz"))
gene <- as.data.frame(gtf)
gene <- gene[!gene$seqnames %in% c("chrX", "chrY", "chrM"), ]
gene <- gene[gene$type == "exon", c("gene_name", "gene_type", "gene_id")]
gene <- gene[!duplicated(gene), ]
gene <- gene[gene$gene_type == "protein_coding", ]

nrow(gene)
# [1] 19111

all_gene_names <- unique(gene$gene_name)


# LOAD THE INDEPENDENT SIGNALS - see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_38_2_define_independent_signals_using_as_ref_cojo_agnostic_results.R
list_genes <- as.data.frame(fread(paste0(path, "post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence_with_coloc_mr_direction_with_ligand_receptor_with_protein_complexes.tsv.gz")))
n_sample   <- nrow(list_genes)


#### list of monogenic IBD genes, extracted from Bolton, high and moderate penetrance:
mono <- fread(paste0(path, "resources/monogenic_pid_genes/list_monogenic_ibd_list_20260201.csv"), head = TRUE)
mono <- mono[mono$'list_2026 (based on Bolton high, moderate penetrance)' == "bolton_high_moderate_penetrance", ]
dim(mono)
# [1] 64 ...

mono_genes    <- mono$Gene
observed_mono <- sum(mono_genes %in% list_genes$list_genes)
cat("Monogenic IBD genes in effector list:", observed_mono, "\n")
# 7


#### list of IBD_PID genes by GEL:
pid <- fread(paste0(path, "resources/monogenic_pid_genes/IBD_PID_GenomicsEngland_v8.75_20260201.csv"), head = TRUE)
pid <- pid[, 1, drop = FALSE]
colnames(pid) <- "Gene"

pid_genes    <- pid$Gene
observed_pid <- sum(pid_genes %in% list_genes$list_genes)
cat("IBD-PID genes in effector list:", observed_pid, "\n")
# 42


#### Permutation test (500,000 iterations) ####
# Uses replicate() instead of a for loop for ~10x speed improvement
perm_counts <- replicate(500000, {
  sampled <- sample(all_gene_names, n_sample, replace = FALSE)
  c(mono = sum(mono_genes %in% sampled),
    pid  = sum(pid_genes  %in% sampled))
})

iter <- data.frame(
  N_monogenic = perm_counts["mono", ],
  N_PID       = perm_counts["pid",  ]
)

table(iter$N_monogenic)
table(iter$N_PID)


#### Monogenic enrichment
p_mono <- mean(iter$N_monogenic >= observed_mono)
or_mono <- observed_mono / mean(iter$N_monogenic)
cat("Monogenic p-value:", p_mono, "\n")   # ~0.002358
cat("Monogenic OR:",     or_mono, "\n")   # ~3.799536

#### IBD_PID enrichment
p_pid  <- mean(iter$N_PID >= observed_pid)
or_pid <- observed_pid / mean(iter$N_PID)
cat("IBD-PID p-value:", p_pid,  "\n")   # 2X-06
cat("IBD-PID OR:",     or_pid,  "\n")   # ~2.291823


#### Annotate list_genes with monogenic status ####
list_genes$monogenic_ibd <- as.integer(list_genes$list_genes %in% mono_genes)

# Reload full mono table to get effect direction columns
mono_full <- fread(paste0(path, "resources/monogenic_pid_genes/list_monogenic_ibd_list_20260201.csv"), head = TRUE)
mono_full <- mono_full[mono_full$'list_2026 (based on Bolton high, moderate penetrance)' == "bolton_high_moderate_penetrance", ]

table(mono_full$'Effect  (Uhlig JPGN 2020)', mono_full$'Effect Bolton (pre-publication)')
#               ? GOF LOF
#            1  0   0  11
# GOF        0  0   1   4
# LOF        0  2   3  41
# LOF & GOF  0  0   0   1

# Prefer Uhlig effect; fall back to Bolton where Uhlig is blank
mono_full$monogenic_effect <- mono_full$'Effect  (Uhlig JPGN 2020)'
blank <- mono_full$monogenic_effect == ""
mono_full$monogenic_effect[blank] <- mono_full$'Effect Bolton (pre-publication)'[blank]

table(mono_full$monogenic_effect)
#         GOF       LOF LOF & GOF
#   1       5        57         1

list_genes <- merge(list_genes, mono_full[, c("Gene", "monogenic_effect")],
                    by.x = "list_genes", by.y = "Gene", all.x = TRUE)

list_genes$monogenic_ibd_pid <- as.integer(list_genes$list_genes %in% pid_genes)


fwrite(list_genes,
       paste0(path, "post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence_with_coloc_mr_direction_with_ligand_receptor_with_protein_complexes_with_monogenic.tsv.gz"),
       col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

q("no")
