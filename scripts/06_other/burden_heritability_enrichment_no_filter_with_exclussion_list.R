# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# Originally writen by Qian Zhang, modified by Laura Fachal

set.seed(314)

args <- commandArgs(T)
infile <- as.character(args[1])
outfile <- as.character(args[2])
n_perm <- as.double(args[3])

ph <- as.character(args[4])
exclussion_list<-args[5]

# test purposes:
# ph <- "ibd"
# infile <- paste0("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/coloc_results/test/allgenes_",ph,"_coloc_rank1_mr_exonic")
# outfile <-  paste0("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/coloc_results/test/allgenes_",ph,"_coloc_rank1_mr_exonic_nod2_exclussion_list_output")
# n_perm <- 500000
# exclussion_list <- paste0("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/coloc_results/test/nod2_exclussion_list")


#reference, IBD burden test summary statistics based on NSYN_AM mask
gene_reference <- readRDS(paste0("/path/to/project",ph,"_nsyn_am.rds"))

# load the exclussion list so it does not get included for the null distribution estimation
if (file.exists(exclussion_list)) {
  geneexclussionlist <- read.table(exclussion_list,header=F)
  gene_reference <- gene_reference[which(!gene_reference$geneid %in% geneexclussionlist$V1),]
}

ibd_split_by_chr <- split(gene_reference, gene_reference$CHROM)
gene_pos_lookup <- setNames(gene_reference$GENPOS, gene_reference$geneid)
gene_chr_lookup <- setNames(gene_reference$CHROM, gene_reference$geneid)
gene_r2_lookup <- setNames(gene_reference$r2, gene_reference$geneid)


#get target gene list
genelist <- read.table(infile,header=F)
geneid <- genelist[,1]


#keep genes available in the reference
geneid <- geneid[geneid %in% gene_reference$geneid]


#generate background gene pool, to avoid the bias caused by enrichment of GWAS signal in specific GWAS region, only select genes within 1Mbp region of target genes as background.

background_gene_pool <- lapply(geneid, function(gene) {
  chr <- gene_chr_lookup[gene]
  pos <- gene_pos_lookup[gene]
  chr_df <- ibd_split_by_chr[[as.character(chr)]]
  nearby_genes <- chr_df$geneid[chr_df$GENPOS > (pos - 1e6) & chr_df$GENPOS < (pos + 1e6)]
  nearby_genes
})


sampled_matrix <- sapply(background_gene_pool, function(glist) {
  if (length(glist) == 0) return(rep(NA, n_perm))
  sample(glist, size = n_perm, replace = TRUE)
})

sampled_r2_matrix <- matrix(gene_r2_lookup[as.vector(sampled_matrix)], nrow = n_perm)

mean_r2_perm <- rowMeans(sampled_r2_matrix, na.rm = TRUE)

mean_r2 <- mean(gene_r2_lookup[geneid])

print(paste0("Average R2 of target genes = ",mean_r2,"; Permutation P-value = ",sum(mean_r2<mean_r2_perm)/n_perm))

smy <- data.frame(permutate = c("Target_genes",1:n_perm),mean_r2 = c(mean_r2,mean_r2_perm))

write.csv(smy,file=outfile,quote=F,row.names=F)



