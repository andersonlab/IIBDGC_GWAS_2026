# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
##########################################################################
### ADD MR RESULTS
###
### Shell setup (run before launching R):
###
###   # singularity exec iibdgc_postprocess_10_singularity.sif
###   MEM=2000
###   bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R
###
##########################################################################

library(data.table)
library(stringr)
library(qvalue)

rm(list = ls())

path_gwas <- "/path/to/ibdgwas/IIBDGC/"

# -------------------------------------------------------------------------
# Helper: clean up pipe-separated strings produced by paste(..., NA, ...)
# -------------------------------------------------------------------------
clean_pipe_string <- function(x) {
    x <- gsub("\\|\\|",   "",  x)
    x <- gsub("\\|NA\\|", "|", x)
    x <- gsub("\\|NA$",   "",  x)
    x <- gsub("\\|$",     "",  x)
    x <- gsub("^\\|",     "",  x)
    x[x == "" | x == "NA"] <- NA
    x
}

# -------------------------------------------------------------------------
# Load pQTL protein lists
# -------------------------------------------------------------------------
pt1 <- fread(paste0(path_gwas, "post_imputation/2022/analysis/mendelian_randomization/raw_pqtl_decode/list_proteins_somascan_in_ibd_regions.tsv.gz"))
pt2 <- fread(paste0(path_gwas, "post_imputation/2022/analysis/mendelian_randomization/raw_pqtl_decode/list_proteins_olink_in_ibd_regions.tsv.gz"))

pt <- rbind(pt1, pt2)
rm(pt1, pt2)

pt$condition <- gsub(".txt.gz", "", pt$files)
pt <- pt[, c("condition", "gene_symbol", "source")]

datasets <- c("decode", "sparc")
pheno    <- c("ibd", "cd", "uc")

# -------------------------------------------------------------------------
# Load forward MR results (pQTL -> IBD)
# -------------------------------------------------------------------------
mr_list <- list()
for (ph in pheno) {
    for (ds in datasets) {
        tmp <- fread(paste0(path_gwas,
            "post_imputation/2022/analysis/mendelian_randomization/",
            ph, "/pqtl_", ds, "_results/output_files/", ph, "_ld_0.01.tsv.gz"))
        tmp$pheno <- ph
        tmp <- tmp[tmp$N_instruments >= 10, ]
        if (nrow(tmp) > 0) {
            tmp$qval <- qvalue(tmp$bxy_pval)$qval
            mr_list[[length(mr_list) + 1]] <- tmp
        }
    }
}
mr <- rbindlist(mr_list)
rm(mr_list)

dim(mr)
# [1] 1904   11

table(mr$pheno)
#  cd ibd  uc
# 639 632 633

# -------------------------------------------------------------------------
# Load reverse MR results (IBD -> pQTL)
# -------------------------------------------------------------------------
rmr_list <- list()
for (ph in pheno) {
    tmp <- fread(paste0(path_gwas,
        "post_imputation/2022/analysis/mendelian_randomization/reverse_mr/",
        ph, "/pqtl_results/output_files/", ph, "_reverse_mr_ld_0.01.tsv.gz"))
    tmp$pheno <- ph
    tmp$qval  <- qvalue(tmp$bxy_pval)$qval
    rmr_list[[length(rmr_list) + 1]] <- tmp
}
rmr <- rbindlist(rmr_list)
rmr <- unique(rmr)
rm(rmr_list)

dim(rmr)
# [1] 1904   10

# -------------------------------------------------------------------------
# Build condition_pheno key and filter rmr to conditions present in mr
# -------------------------------------------------------------------------
mr$condition_pheno  <- paste(mr$condition,  mr$pheno, sep = "_")
rmr$condition_pheno <- paste(rmr$condition, rmr$pheno, sep = "_")

# Keep only reverse-MR rows that match a forward-MR row
rmr <- rmr[rmr$condition_pheno %in% mr$condition_pheno, ]

rmr <- rmr[, c("condition_pheno", "N_instruments", "bxy", "bxy_se", "bxy_pval",
               "index_snps", "pleio_snps", "pvalue_on_risk_factor",
               "pvalue_on_pqtl", "qval")]

colnames(rmr)[2:ncol(rmr)] <- paste0(colnames(rmr)[2:ncol(rmr)], "_reverse_mr")

# -------------------------------------------------------------------------
# Merge forward + reverse MR; attach gene/source annotation
# -------------------------------------------------------------------------
mr <- merge(mr,  rmr, by = "condition_pheno")
mr <- merge(pt,  mr,  by = "condition", all.y = TRUE)

x <- as.data.frame.matrix(table(mr$gene_symbol, mr$source))
print(nrow(x[x$decode != 0 & x$ukb != 0, ]))
# [1] 56

# Explore overlap across q-value threshold combinations
qval_1 <- c(0.05, 0.01, 0.005, 0.05, 0.01)
qval_2 <- c(0.05, 0.01, 0.005, 0.01, 0.05)

for (i in seq_along(qval_1)) {
    print(paste("MR qval <", qval_1[i], "; reverse MR qval >=", qval_2[i]))
    sel <- mr$qval < qval_1[i] & mr$qval_reverse_mr >= qval_2[i]
    x   <- as.data.frame.matrix(table(mr$gene_symbol[sel], mr$source[sel]))
    print(nrow(x[x$decode != 0 & x$ukb != 0, ]))
}
# [1] "MR qval < 0.05 ; reverse MR qval >= 0.05"   -> 1 ***** selected
# [1] "MR qval < 0.01 ; reverse MR qval >= 0.01"   -> 2
# [1] "MR qval < 0.005 ; reverse MR qval >= 0.005" -> 2
# [1] "MR qval < 0.05 ; reverse MR qval >= 0.01"   -> 4 
# [1] "MR qval < 0.01 ; reverse MR qval >= 0.05"   -> 0

mr <- mr[mr$qval < 0.05 & mr$qval_reverse_mr >= 0.05, ]
dim(mr)
# 166

# -------------------------------------------------------------------------
# Link MR genes to GWAS regions (by gene name overlap)
# -------------------------------------------------------------------------
reg  <- fread(paste0(path_gwas, "post_imputation/2022/analysis/metaanalysis/forward_regression/list_regions_420_cojo_supervised_plus_unsupervised_with_gene_names.tsv.gz"))
pqtl <- names(table(mr$gene_symbol))

length(pqtl)
# [1] 117

#  [1] "ADAMTSL5" "AKR7L"    "AMBP"     "AOC3"     "APOBR"    "BCL2L11" 
#   [7] "BLVRB"    "BPIFA2"   "BRICD5"   "CAPS"     "CCL11"    "CCL13"   
#  [13] "CCL19"    "CCL25"    "CD200R1"  "CD28"     "CD300LG"  "CD302"   
#  [19] "CD3E"     "CD4"      "CD48"     "CD59"     "CD84"     "CDH1"    
#  [25] "CELA3A"   "CRTAM"    "CST6"     "CTRL"     "CTSO"     "CXCL16"  
#  [31] "CXCL5"    "CXCL8"    "DDC"      "DLL1"     "EDC4"     "ENTR1"   
#  [37] "ERBB2"    "ERLIN1"   "FCGR2B"   "FGF2"     "FGF8"     "FLT3LG"  
#  [43] "FURIN"    "GAL3ST2"  "GPN1"     "GPT"      "HGFAC"    "HMOX1"   
#  [49] "ICOSLG"   "IFNA2"    "IFNAR1"   "IFNGR1"   "IGSF9"    "IL10RA"  
#  [55] "IL12B"    "IL1R1"    "IL1R2"    "IL1RL1"   "IL2RB"    "INHBC"   
#  [61] "INSL5"    "ITGAV"    "ITGB6"    "ITIH3"    "KIR2DL1"  "KLK1"    
#  [67] "KLRB1"    "KPNA1"    "LDLR"     "LINGO1"   "LRPAP1"   "LRRC25"  
#  [73] "LTBR"     "LYPD3"    "MACROD1"  "MELTF"    "MET"      "MXRA8"   
#  [79] "MYL3"     "MYRF"     "NCR1"     "NOTCH1"   "NRBP1"    "NRP1"    
#  [85] "OMD"      "PAM"      "PDCD1"    "PLB1"     "PLTP"     "PLXNB2"  
#  [91] "PPY"      "PROCR"    "PSMB1"    "PVALB"    "RGS19"    "RNASET2" 
#  [97] "ROS1"     "SCARA3"   "SDF2L1"   "SELL"     "SELP"     "SEPTIN3" 
# [103] "SMAD3"    "SPINT1"   "STAT5B"   "TDRKH"    "TFF3"     "TFIP11"  
# [109] "TIMD4"    "TNFRSF1B" "TNFRSF4"  "TNFSF11"  "TNFSF12"  "TRADD"   
# [115] "VAV3"     "VCAM1"    "VWF" 


reg <- as.data.frame(reg)
reg$pqtl_mr_cd  <- NA_character_
reg$pqtl_mr_uc  <- NA_character_
reg$pqtl_mr_ibd <- NA_character_

for (gene in pqtl) {
    tmp      <- mr[mr$gene_symbol == gene, ]
    tmp$mr_results <- paste(tmp$gene_symbol, tmp$source,
                            formatC(tmp$bxy,              format = "f", digits = 2),
                            formatC(tmp$bxy_pval,         format = "E", digits = 2),
                            formatC(tmp$qval,             format = "E", digits = 2),
                            formatC(tmp$bxy_pval_reverse_mr, format = "E", digits = 2),
                            formatC(tmp$qval_reverse_mr,  format = "E", digits = 2),
                            sep = ",")

    tmp_reg <- reg[grep(gene, reg$genes_coding), c("updated_region", "genes_coding")]

    # Narrow to rows where gene is an exact element (not a substring)
    exact_match <- sapply(tmp_reg$genes_coding, function(gc) {
        gene %in% str_split(gc, "\\|", simplify = TRUE)
    })
    tmp_reg <- tmp_reg[exact_match, , drop = FALSE]

    if (nrow(tmp_reg) > 0) {
        rows <- reg$updated_region == tmp_reg$updated_region[1]
        reg$pqtl_mr_cd[rows]  <- paste(tmp$mr_results[tmp$pheno == "cd"],  reg$pqtl_mr_cd[rows],  sep = "|")
        reg$pqtl_mr_uc[rows]  <- paste(tmp$mr_results[tmp$pheno == "uc"],  reg$pqtl_mr_uc[rows],  sep = "|")
        reg$pqtl_mr_ibd[rows] <- paste(tmp$mr_results[tmp$pheno == "ibd"], reg$pqtl_mr_ibd[rows], sep = "|")
    }
}

reg <- reg[, c("updated_region", "genes_coding", "pqtl_mr_cd", "pqtl_mr_uc", "pqtl_mr_ibd")]
reg$pqtl_mr_cd  <- clean_pipe_string(reg$pqtl_mr_cd)
reg$pqtl_mr_uc  <- clean_pipe_string(reg$pqtl_mr_uc)
reg$pqtl_mr_ibd <- clean_pipe_string(reg$pqtl_mr_ibd)

## How many regions have at least one MR hit?
dim(reg[!is.na(reg$pqtl_mr_ibd) | !is.na(reg$pqtl_mr_cd) | !is.na(reg$pqtl_mr_uc), ])
# [1] 89  5

fwrite(reg,
    paste0(path_gwas, "post_imputation/2022/analysis/final_tables/list_regions_420_cojo_supervised_plus_unsupervised_with_gene_names_with_GSMR_results.tsv.gz"),
    col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# -------------------------------------------------------------------------
# Extract instrument SNPs for LD calculation
# -------------------------------------------------------------------------
all <- fread(paste0(path_gwas,
    "post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz"))
all <- as.data.frame(all)
dim(all)
# [1] 619 182

snps <- unique(c(
    str_split(mr$index_snps, "\\|", simplify = FALSE) |> unlist(),
    all$MarkerName
))
snps <- snps[snps != ""]

length(snps)
# [1] 3118

write.table(snps,
    paste0(path_gwas, "post_imputation/2022/analysis/metaanalysis/forward_regression/list_strong_instruments_for_ld.tsv"),
    col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# -------------------------------------------------------------------------
# Shell: extract subsets and compute LD (run in bash, not R)
# -------------------------------------------------------------------------
###
###   path_gwas=/path/to/ibdgwas/IIBDGC/
###   ph=ibd
###   MEM=2000
###
###   # Step 1: extract variants of interest per chromosome
###   for chr in {1..22}; do
###     bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" \
###       -G your_hpc_group -q normal \
###       -o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld1_stdout \
###       -e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld1_stderr \
###       "/path/to/software/username/./plink2 \
###         --bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_${ph}_analysis \
###         --extract ${path_gwas}post_imputation/2022/analysis/metaanalysis/forward_regression/list_strong_instruments_for_ld.tsv \
###         --make-bed --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_strong_instruments_for_ld"
###   done
###
###   # Step 2: compute pairwise LD
###   for chr in {1..22}; do
###     bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" \
###       -G your_hpc_group -q normal \
###       -o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld2_stdout \
###       -e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld2_stderr \
###       "/path/to/software/username/plink_linux_x86_64_20181202/./plink \
###         --bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_strong_instruments_for_ld \
###         --r2 --ld-window-kb '1000000000' --ld-window '1000000000' --ld-window-r2 '0' \
###         --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_strong_instruments_for_ld"
###   done
###

# -------------------------------------------------------------------------
# Read LD files back into R
# -------------------------------------------------------------------------
ld_list <- lapply(1:22, function(chr) {
    f <- paste0(path_gwas,
        "post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr",
        chr, "_list_strong_instruments_for_ld.ld")
    tmp <- read.table(f, header = TRUE)
    tmp[tmp$R2 >= 0.1, ]
})
ld <- do.call(rbind, ld_list)
rm(ld_list)

dim(ld)
# [1] 120  7

ld1 <- ld[ld$SNP_A %in% all$MarkerName, c("SNP_A", "SNP_B", "R2")]
ld2 <- ld[ld$SNP_B %in% all$MarkerName, c("SNP_B", "SNP_A", "R2")]
colnames(ld1) <- colnames(ld2) <- c("MarkerName", "SNP_mr", "R2")
ld <- rbind(ld1, ld2)

dim(ld)
# [1] 86   3

ld <- ld[ld$R2 > 0.6, ]
dim(ld)
# [1] 14  3

# -------------------------------------------------------------------------
# Tag MR instruments with matching GWAS index variants
# -------------------------------------------------------------------------
mr <- as.data.frame(mr)
mr$index_gwas_variant    <- NA_character_
mr$index_gwas_variant_ld <- NA_character_

for (i in seq_len(nrow(ld))) {
    tomatch <- c(ld$MarkerName[i], ld$SNP_mr[i])
    pattern  <- paste(tomatch, collapse = "|")
    hits     <- grep(pattern, mr$index_snps)
    mr$index_gwas_variant[hits]    <- ld$MarkerName[i]
    mr$index_gwas_variant_ld[hits] <- paste(ld$MarkerName[i], ld$SNP_mr[i], ld$R2[i], sep = "|")
}

dim(mr[!is.na(mr$index_gwas_variant), ])
# [1] 20 24

table(mr$gene_symbol[!is.na(mr$index_gwas_variant)])

# -------------------------------------------------------------------------
# Attach MR results to the main signals table
# -------------------------------------------------------------------------
all <- fread(paste0(path_gwas,
    "post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz"))
all <- as.data.frame(all)
dim(all)
# [1] 631 182

all <- all[, c("MarkerName", "updated_region", "phenotype", "class_signal_final")]

all$mr_ibd      <- NA_character_
all$mr_cd       <- NA_character_
all$mr_uc       <- NA_character_
all$mr_gene_ibd <- NA_character_
all$mr_gene_cd  <- NA_character_
all$mr_gene_uc  <- NA_character_

mr_hit <- mr[!is.na(mr$index_gwas_variant), ]

for (i in seq_len(nrow(all))) {
    tmp_region <- mr_hit[mr_hit$index_gwas_variant == all$MarkerName[i], ]
    if (nrow(tmp_region) == 0) next

    for (ph in pheno) {
        tmp <- tmp_region[tmp_region$pheno == ph, ]
        if (nrow(tmp) == 0) next

        tmp$mr_results <- paste(tmp$gene_symbol, tmp$source,
                                formatC(tmp$bxy,              format = "f", digits = 2),
                                formatC(tmp$bxy_pval,         format = "E", digits = 2),
                                formatC(tmp$qval,             format = "E", digits = 2),
                                formatC(tmp$bxy_pval_reverse_mr, format = "E", digits = 2),
                                formatC(tmp$qval_reverse_mr,  format = "E", digits = 2),
                                sep = ",")

        all[i, paste0("mr_", ph)]      <- paste(tmp$mr_results, collapse = "|")
        all[i, paste0("mr_gene_", ph)] <- paste(unique(tmp$gene_symbol), collapse = "|")
    }
}

dim(all[!is.na(all$mr_ibd) | !is.na(all$mr_cd) | !is.na(all$mr_uc), ])
# [1] 12 10

list_genes_linked_by_gsmr <- unique(na.omit(c(
    unlist(strsplit(all$mr_gene_cd,  "\\|")),
    unlist(strsplit(all$mr_gene_ibd, "\\|")),
    unlist(strsplit(all$mr_gene_uc,  "\\|"))
)))
list_genes_linked_by_gsmr

fwrite(all,
    paste0(path_gwas, "post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_mr_results.tsv.gz"),
    col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

all[!is.na(all$mr_ibd) | !is.na(all$mr_cd) | !is.na(all$mr_uc),
    c("MarkerName", "phenotype", "mr_gene_cd", "mr_gene_ibd", "mr_gene_uc", "class_signal_final")]

#              MarkerName       phenotype mr_gene_cd mr_gene_ibd mr_gene_uc
# 59    chr1:62795685:C:T IBD_unsaturated       LDLR        <NA>       NRP1
# 182  chr14:68794755:A:G IBD_unsaturated      MELTF        <NA>       <NA>
# 208  chr16:23242699:T:C              UC       <NA>       INSL5       <NA>
# 250  chr17:39932164:T:C IBD_unsaturated       SELL        NCR1       SELL
# 313  chr2:110858564:C:T              UC      IL2RB       CRTAM       <NA>
# 328 chr2:181459315:G:GT   IBD_saturated       <NA>       VCAM1       <NA>
# 359   chr2:62333197:A:G              CD CD28|PDCD1        <NA>       CD48
# 363    chr2:9521321:A:G IBD_unsaturated       <NA>    TNFRSF1B   TNFRSF1B
# 494  chr5:96970610:A:AT   IBD_saturated TIMD4|LTBR        <NA>       <NA>
# 517   chr6:19787512:G:A IBD_unsaturated    TNFSF11        <NA>       <NA>
# 577  chr8:129592317:G:A IBD_unsaturated       <NA>        <NA>     FLT3LG
# 586   chr8:58830080:G:C              UC      CCL13        <NA>       <NA>

pqtl[which(!pqtl %in% c(all$mr_gene_cd,all$mr_gene_uc,all$mr_gene_ibd,"CD28","PDCD1","TIMD4","LTBR"))]
#  [1] "ADAMTSL5" "AKR7L"    "AMBP"     "AOC3"     "APOBR"    "BCL2L11" 
#  [7] "BLVRB"    "BPIFA2"   "BRICD5"   "CAPS"     "CCL11"    "CCL19"   
# [13] "CCL25"    "CD200R1"  "CD300LG"  "CD302"    "CD3E"     "CD4"     
# [19] "CD59"     "CD84"     "CDH1"     "CELA3A"   "CST6"     "CTRL"    
# [25] "CTSO"     "CXCL16"   "CXCL5"    "CXCL8"    "DDC"      "DLL1"    
# [31] "EDC4"     "ENTR1"    "ERBB2"    "ERLIN1"   "FCGR2B"   "FGF2"    
# [37] "FGF8"     "FURIN"    "GAL3ST2"  "GPN1"     "GPT"      "HGFAC"   
# [43] "HMOX1"    "ICOSLG"   "IFNA2"    "IFNAR1"   "IFNGR1"   "IGSF9"   
# [49] "IL10RA"   "IL12B"    "IL1R1"    "IL1R2"    "IL1RL1"   "INHBC"   
# [55] "ITGAV"    "ITGB6"    "ITIH3"    "KIR2DL1"  "KLK1"     "KLRB1"   
# [61] "KPNA1"    "LINGO1"   "LRPAP1"   "LRRC25"   "LYPD3"    "MACROD1" 
# [67] "MET"      "MXRA8"    "MYL3"     "MYRF"     "NOTCH1"   "NRBP1"   
# [73] "OMD"      "PAM"      "PLB1"     "PLTP"     "PLXNB2"   "PPY"     
# [79] "PROCR"    "PSMB1"    "PVALB"    "RGS19"    "RNASET2"  "ROS1"    
# [85] "SCARA3"   "SDF2L1"   "SELP"     "SEPTIN3"  "SMAD3"    "SPINT1"  
# [91] "STAT5B"   "TDRKH"    "TFF3"     "TFIP11"   "TNFRSF4"  "TNFSF12" 
# [97] "TRADD"    "VAV3"     "VWF"   

pqtl[which(pqtl %in% c(all$mr_gene_cd,all$mr_gene_uc,all$mr_gene_ibd))]
#  [1] "CCL13"    "CD48"     "CRTAM"    "FLT3LG"   "IL2RB"    "INSL5"   
#  [7] "LDLR"     "MELTF"    "NCR1"     "NRP1"     "SELL"     "TNFRSF1B"
# [13] "TNFSF11"  "VCAM1"


q("no")
