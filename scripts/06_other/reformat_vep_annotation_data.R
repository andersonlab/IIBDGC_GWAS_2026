# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# reformat_vep_annotation_data

# ## VEP

unset BCFTOOLS_PLUGINS
# singularity exec iibdgc_postprocess_10_singularity.sif

# get annotation data:

path_gwas=/path/to/ibdgwas/IIBDGC/

# # create new annotation files using sanger pipeline
  
sed 's/:\+/\t/g' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed | \
awk -v OFS="\t" '{print $6,$7}' > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp

awk -v OFS="\t" '{$10="0";print $1, $4, $10, $2}' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/list_union_variants_cd_uc_ibd_metaanalysis_sorted.map

paste -d'\t' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/list_union_variants_cd_uc_ibd_metaanalysis_sorted.map \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp> \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bim


# equals tped file:
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bim > \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/list_union_variants_cd_uc_ibd_metaanalysis_sorted.tped

# tfam:
echo '1  1  0  0  1  1' > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/list_union_variants_cd_uc_ibd_metaanalysis_sorted.tfam



MEM=50000
MEM1=100000

bsub -J"annota" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM:tmp=$MEM1]  span[hosts=1]" -G your_hpc_group -n 4 \
-e ${path_gwas}post_imputation/log/allchr_create_het_vcf_stderr \
-o ${path_gwas}post_imputation/log/allchr_create_het_vcf_stdout \
"/path/to/software/username/plink_linux_x86_64_20181202/./plink \
--tfile ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/list_union_variants_cd_uc_ibd_metaanalysis_sorted \
--memory $MEM1 \
--threads 4 \
--export vcf --out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/list_union_variants_cd_uc_ibd_metaanalysis_sorted.vcf"

cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/list_union_variants_cd_uc_ibd_metaanalysis_sorted.vcf | sed 's/^23/X/g' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/list_union_variants_cd_uc_ibd_metaanalysis_sorted_edited.vcf

MEM=300
bsub -J"annota" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]  span[hosts=1]" -G your_hpc_group -n 4 \
-e ${path_gwas}post_imputation/log/create_het_vcf_gz_stderr \
-o ${path_gwas}post_imputation/log/create_het_vcf_gz_stdout \
"bcftools view ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/list_union_variants_cd_uc_ibd_metaanalysis_sorted_edited.vcf -Oz \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/list_union_variants_cd_uc_ibd_metaanalysis_sorted.vcf.gz"

MEM=500
bsub -J"annota" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]  span[hosts=1]" -G your_hpc_group -n 4 \
-e ${path_gwas}post_imputation/log/create_het_vcf_gz_index_stderr \
-o ${path_gwas}post_imputation/log/create_het_vcf_gz_index_stdout \
"bcftools index -f ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/list_union_variants_cd_uc_ibd_metaanalysis_sorted.vcf.gz --threads 4"


# split into chr:

MEM=150

for chr in {1..23} X
do
bsub -J"annota" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]  span[hosts=1]" -G your_hpc_group -n 2 \
-e ${path_gwas}post_imputation/log/chr${chr}_create_het_vcf_stderr \
-o ${path_gwas}post_imputation/log/chr${chr}_create_het_vcf_stdout \
"bcftools view ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/list_union_variants_cd_uc_ibd_metaanalysis_sorted.vcf.gz \
--regions ${chr} -Oz -o ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_sorted.vcf.gz --threads 4"
done


for chr in {1..22} X
do
ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_sorted.vcf.gz
done


MEM=3000
for chr in {1..22} X
do
bsub -J"annota" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q long \
-e ${path_gwas}post_imputation/log/${chr}_annotation_vep_sanger_stderr \
-o ${path_gwas}post_imputation/log/${chr}_annotation_vep_sanger_stdout \
"/path/to/project ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_sorted.vcf.gz"
done

mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/*vep.* vep_110.1/


# index the files
path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=800

for chr in 22
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 2 \
-e ${path_gwas}post_imputation/log/${chr}_vep_index_stderr \
-o ${path_gwas}post_imputation/log/${chr}_vep_index_stdout \
"bcftools index -f ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/vep_110.1/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_sorted.vep.vcf.gz --threads 8"
done

for chr in {1..22} X
do 
echo ${chr} && tail -50 ${path_gwas}post_imputation/log/${chr}_vep_index_stdout | grep "Successfully"
done


# see headers - v101
bcftools +split-vep \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/vep_110.1/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_sorted.vep.vcf.gz \
-l
# 0       Allele
# 1       Consequence
# 2       IMPACT
# 3       SYMBOL
# 4       Gene
# 5       Feature_type
# 6       Feature
# 7       BIOTYPE
# 8       EXON
# 9       INTRON
# 10      HGVSc
# 11      HGVSp
# 12      cDNA_position
# 13      CDS_position
# 14      Protein_position
# 15      Amino_acids
# 16      Codons
# 17      Existing_variation
# 18      ALLELE_NUM
# 19      DISTANCE
# 20      STRAND
# 21      FLAGS
# 22      VARIANT_CLASS
# 23      SYMBOL_SOURCE
# 24      HGNC_ID
# 25      CANONICAL
# 26      MANE_SELECT
# 27      MANE_PLUS_CLINICAL
# 28      TSL
# 29      APPRIS
# 30      CCDS
# 31      ENSP
# 32      SWISSPROT
# 33      TREMBL
# 34      UNIPARC
# 35      UNIPROT_ISOFORM
# 36      SOURCE
# 37      GENE_PHENO
# 38      SIFT
# 39      PolyPhen
# 40      DOMAINS
# 41      miRNA
# 42      HGVS_OFFSET
# 43      AF
# 44      AFR_AF
# 45      AMR_AF
# 46      EAS_AF
# 47      EUR_AF
# 48      SAS_AF
# 49      gnomADe_AF
# 50      gnomADe_AFR_AF
# 51      gnomADe_AMR_AF
# 52      gnomADe_ASJ_AF
# 53      gnomADe_EAS_AF
# 54      gnomADe_FIN_AF
# 55      gnomADe_NFE_AF
# 56      gnomADe_OTH_AF
# 57      gnomADe_SAS_AF
# 58      gnomADg_AF
# 59      gnomADg_AFR_AF
# 60      gnomADg_AMI_AF
# 61      gnomADg_AMR_AF
# 62      gnomADg_ASJ_AF
# 63      gnomADg_EAS_AF
# 64      gnomADg_FIN_AF
# 65      gnomADg_MID_AF
# 66      gnomADg_NFE_AF
# 67      gnomADg_OTH_AF
# 68      gnomADg_SAS_AF
# 69      MAX_AF
# 70      MAX_AF_POPS
# 71      CLIN_SIG
# 72      SOMATIC
# 73      PHENO
# 74      PUBMED
# 75      MOTIF_NAME
# 76      MOTIF_POS
# 77      HIGH_INF_POS
# 78      MOTIF_SCORE_CHANGE
# 79      TRANSCRIPTION_FACTORS
# 80      SpliceRegion
# 81      GeneSplicer
# 82      existing_InFrame_oORFs
# 83      existing_OutOfFrame_oORFs
# 84      existing_uORFs
# 85      five_prime_UTR_variant_annotation
# 86      five_prime_UTR_variant_consequence
# 87      CADD_PHRED
# 88      CADD_RAW
# 89      AlphaMissense_pred
# 90      AlphaMissense_rankscore
# 91      AlphaMissense_score
# 92      Ensembl_transcriptid
# 93      LRT_pred
# 94      MutationTaster_pred
# 95      Polyphen2_HDIV_pred
# 96      Polyphen2_HVAR_pred
# 97      PrimateAI_pred
# 98      PrimateAI_rankscore
# 99      PrimateAI_score
# 100     SIFT_pred
# 101     Uniprot_acc
# 102     VEP_canonical
# 103     DisGeNET_PMID
# 104     DisGeNET_SCORE
# 105     PHENOTYPES
# 106     Conservation
# 107     LoF
# 108     LoF_filter
# 109     LoF_flags
# 110     LoF_info
# 111     REVEL
# 112     SpliceAI_pred_DP_AG
# 113     SpliceAI_pred_DP_AL
# 114     SpliceAI_pred_DP_DG
# 115     SpliceAI_pred_DP_DL
# 116     SpliceAI_pred_DS_AG
# 117     SpliceAI_pred_DS_AL
# 118     SpliceAI_pred_DS_DG
# 119     SpliceAI_pred_DS_DL
# 120     SpliceAI_pred_SYMBOL
# 121     EVE_CLASS
# 122     EVE_SCORE
# 123     ClinVar
# 124     ClinVar_CLNDN
# 125     ClinVar_CLNDNINCL
# 126     ClinVar_CLNDISDB
# 127     ClinVar_CLNDISDBINCL
# 128     ClinVar_AF_ESP
# 129     ClinVar_AF_EXAC
# 130     ClinVar_AF_TGP
# 131     ClinVar_ALLELEID
# 132     ClinVar_CLNDN
# 133     ClinVar_CLNDNINCL
# 134     ClinVar_CLNDISDB
# 135     ClinVar_CLNDISDBINCL
# 136     ClinVar_CLNHGVS
# 137     ClinVar_CLNREVSTAT
# 138     ClinVar_CLNSIG
# 139     ClinVar_CLNSIGCONF
# 140     ClinVar_CLNSIGINCL
# 141     ClinVar_CLNVC
# 142     ClinVar_CLNVCSO
# 143     ClinVar_CLNVI
# 144     ClinVar_DBVARID
# 145     ClinVar_GENEINFO
# 146     ClinVar_MC
# 147     ClinVar_ORIGIN
# 148     ClinVar_RS
# 149     ClinVar_SSR

# bcftools +split-vep -S -
  # Default consequence substrings ordered in ascending order by severity.
  # Consequences with the same severity can be put on the same line in arbitrary order.
  # See also https://m.ensembl.org/info/genome/variation/prediction/predicted_data.htm
# intergenic
# feature_truncation feature_elongation
# regulatory
# TF_binding_site TFBS
# downstream upstream
# non_coding_transcript non_coding
# intron NMD_transcript
# non_coding_transcript_exon
# 5_prime_utr 3_prime_utr
# coding_sequence mature_miRNA
# stop_retained start_retained synonymous
# incomplete_terminal_codon
# splice_region
# missense inframe protein_altering
# transcript_amplification
# exon_loss
# disruptive
# start_lost stop_lost stop_gained frameshift
# splice_acceptor splice_donor
# transcript_ablation


# ##################################################
# # # retrieve inframe del and ins plus missense:

bcftools \
+split-vep ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/vep_110.1/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_sorted.vep.vcf.gz \
-c Consequence -s :inframe \
| bcftools query -f'%ID\t%Consequence\n' | cut -f2 | sort | uniq -c

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=100

for chr in {1..22} X
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q small \
-e ${path_gwas}post_imputation/log/${chr}_vep_missense_stderr \
-o ${path_gwas}post_imputation/log/${chr}_vep_missense_stdout \
"bcftools \
+split-vep ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/vep_110.1/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_sorted.vep.vcf.gz \
-c Consequence -s :inframe \
| bcftools query -f'%CHROM\t%POS\t%ID\t%Consequence\n' > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_chr${chr}_missense"
done

# concatenate all files, sort, add code=1 when missense, and add header:
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_chr*_missense \
| sort -k1,1V -k2,2n | awk -F'\t' 'BEGIN {OFS=FS} {if ($4==".") print $3,"0";else print $3,"1"}' \
| awk -F'\t' -v OFS='\t' 'BEGIN {print "variant","missense"}1' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/allchr_missense_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/allchr_missense_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed
# 49283611
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/allchr_missense_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed | cut -f2 | sort | uniq -c
# 48983350 0
# 300260 1
# 1 missense

# remove tmp files
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_*
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/*_sorted.vcf.gz

############################################################
# # retrieves start_lost stop_lost stop_gained frameshift
# bcftools \
# +split-vep ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/vep_110.1/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_sorted.vep.vcf.gz \
# -c Consequence -s :frameshift \
# | bcftools query -f'%ID\t%Consequence\n' | cut -f2 | sort | uniq -c

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=100
for chr in {1..22} X
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q small \
-e ${path_gwas}post_imputation/log/${chr}_vep_frameshift_stderr \
-o ${path_gwas}post_imputation/log/${chr}_vep_frameshift_stdout \
"bcftools \
+split-vep ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/vep_110.1/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_sorted.vep.vcf.gz \
-c Consequence -s :frameshift \
| bcftools query -f'%CHROM\t%POS\t%ID\t%Consequence\n' > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_chr${chr}_frameshift"
done

# concatenate all files, sort, add code=1 when missense, and add header:
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_chr*_frameshift \
| sort -k1,1V -k2,2n | awk -F'\t' 'BEGIN {OFS=FS} {if ($4==".") print $3,"0";else print $3,"1"}' \
| awk -F'\t' -v OFS='\t' 'BEGIN {print "variant","frameshift"}1' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/allchr_frameshift_list_union_variants_cd_uc_ibd_metaanalysis.bed
wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/allchr_frameshift_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed
# 49283611

cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/allchr_frameshift_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed \
| cut -f2 | sort | uniq -c
# 49269073 0
# 14537 1
# 1 frameshift

# remove tmp files
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_*
  
# #############################################
# # # retrieves splice_acceptor splice_donor
# # bcftools \
# # +split-vep ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/vep_110.1/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_sorted.vep.vcf.gz \
# # -c Consequence -s :splice_donor \
# # | bcftools query -f'%ID\t%Consequence\n' | cut -f2 | sort | uniq -c

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=100
for chr in {1..22} X
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q small \
-e ${path_gwas}post_imputation/log/${chr}_vep_splice_stderr \
-o ${path_gwas}post_imputation/log/${chr}_vep_splice_stdout \
"bcftools \
+split-vep ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/vep_110.1/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_sorted.vep.vcf.gz \
-c Consequence -s :splice_donor \
| bcftools query -f'%CHROM\t%POS\t%ID\t%Consequence\n' > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_chr${chr}_splice"
done

# cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_chr1_splice \
# | cut -f4 | sort | uniq -c

# concatenate all files, sort, add code=1 when missense, and add header:
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_chr*_splice \
| sort -k1,1V -k2,2n | awk -F'\t' 'BEGIN {OFS=FS} {if ($4==".") print $3,"0";else print $3,"1"}' \
| awk -F'\t' -v OFS='\t' 'BEGIN {print "variant","splice"}1' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/allchr_splice_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed
wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/allchr_splice_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed
# 49283611
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/allchr_splice_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed \
| cut -f2 | sort | uniq -c
# 49257638 0
#   25972 1
#       1 splice

# remove tmp files
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_*
  
# #########################################################
# # # retrieves stop_retained start_retained synonymous
# # bcftools \
# # +split-vep ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/vep_110.1/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_sorted.vep.vcf.gz \
# # -c Consequence -s :synonymous \
# # | bcftools query -f'%ID\t%Consequence\n' | cut -f2 | sort | uniq -c

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=100
for chr in {1..22} X
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q small \
-e ${path_gwas}post_imputation/log/${chr}_vep_synonymous_stderr \
-o ${path_gwas}post_imputation/log/${chr}_vep_synonymous_stdout \
"bcftools \
+split-vep ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/vep_110.1/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_sorted.vep.vcf.gz \
-c Consequence -s :synonymous \
| bcftools query -f'%CHROM\t%POS\t%ID\t%Consequence\n' > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_chr${chr}_synonymous"
done

cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_chr1_synonymous \
| cut -f4 | sort | uniq -c

# concatenate all files, sort, add code=1 when missense, and add header:
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_chr*_synonymous \
| sort -k1,1V -k2,2n | awk -F'\t' 'BEGIN {OFS=FS} {if ($4==".") print $3,"0";else print $3,"1"}' \
| awk -F'\t' -v OFS='\t' 'BEGIN {print "variant","synonymous"}1' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/allchr_synonymous_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed
wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/allchr_synonymous_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed
# 49283611
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/allchr_synonymous_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed \
| cut -f2 | sort | uniq -c
# 49092050 0
#  191560 1
#       1 synonymous

# remove tmp files
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_*
  
# ##################################################################################################################

# variants that affect ORF:

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=100
for chr in {1..22} X
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q small \
-e ${path_gwas}post_imputation/log/${chr}_vep_5utr_stderr \
-o ${path_gwas}post_imputation/log/${chr}_vep_5utr_stdout \
"bcftools \
+split-vep ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/vep_110.1/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_sorted.vep.vcf.gz \
-f '%CHROM\t%POS\t%ID\t%five_prime_UTR_variant_consequence\n' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_chr${chr}_5utr"
done

cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_chr*_5utr \
| sort -k1,1V -k2,2n | awk -F'\t' 'BEGIN {OFS=FS} {print $3,"1"}' > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_allchr_5utr

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_allchr_5utr
# 20164 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_allchr_5utr

awk 'FNR == NR {h[$1] = $1; next;} {FS=OFS="\t"} FNR>0{if ($4 in h) {($5 = "1")} else {($5 = "0")} } 1' \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_allchr_5utr \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed \
| cut -f4-5 | awk -F'\t' -v OFS='\t' 'BEGIN {print "variant","5utr"}1' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/allchr_5utr_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed
cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/allchr_5utr_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed \
| cut -f2 | sort | uniq -c
# 49263446 0
#   20164 1
#       1 5utr

# remove tmp files
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/tmp_*
  
# #################################################

# mv data to final folder
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/allchr_*_list_union_variants_cd_uc_ibd_metaanalysis.bed \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/
  

# EXTRACT maf bins:

MEM=25000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(mltools)
library(Hmisc)

path<-"/path/to/ibdgwas/IIBDGC/"

dat<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed",sep=""),head=F)

pheno<-c("ibd","cd","uc")

chr10:50134753:G:A

for (chr in c(1:23)) {


    if (chr==23) {
        chr<-"X"
    } 

    tmp1<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/",pheno[1],"/",chr,"_",pheno[1],"_meta_eur_eas_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),head=T)
    tmp1<-tmp1[which(tmp1$MarkerName %in% dat$V4),]
    tmp1$MAF<-pmin(tmp1$avgA2FREQ_CONTROLS,1-tmp1$avgA2FREQ_CONTROLS)
    tmp1$MAF[which(is.na(tmp1$avgA2FREQ_CONTROLS))]<-pmin(tmp1$avgA2FREQ[which(is.na(tmp1$avgA2FREQ_CONTROLS))],1-tmp1$avgA2FREQ[which(is.na(tmp1$avgA2FREQ_CONTROLS))])
    tmp1<-tmp1[,c("MarkerName","MAF")]

    for (i in c(2:3)) {

        tmp2<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_eur_eas_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),head=T)
        tmp2<-tmp2[which((tmp2$MarkerName %in% dat$V4) & (!tmp2$MarkerName %in% tmp1$MarkerName)),]
        tmp2$MAF<-pmin(tmp2$avgA2FREQ_CONTROLS,1-tmp2$avgA2FREQ_CONTROLS)
        tmp2$MAF[which(is.na(tmp2$avgA2FREQ_CONTROLS))]<-pmin(tmp2$avgA2FREQ[which(is.na(tmp2$avgA2FREQ_CONTROLS))],1-tmp2$avgA2FREQ[which(is.na(tmp2$avgA2FREQ_CONTROLS))])
        tmp2<-tmp2[,c("MarkerName","MAF")]

        if (nrow(tmp2)>0) {
            tmp1<-rbind(tmp1,tmp2)
        }
        rm(tmp2)

    }
    if (chr==1) {
        tmp<-tmp1
    } else {
        tmp<-rbind(tmp,tmp1)
    }
    rm(tmp1)
}






dim(dat)
# [1] 49283610        4
dim(tmp)
# [1] 49283627        2

tmp<-tmp[!duplicated(tmp),]
tmp<-tmp[which(!duplicated(tmp$MarkerName)),]

dim(tmp[which(!tmp$MarkerName %in% dat$V4),])
# [1] 0   4 # 

summary(tmp$MAF,useNA="ifany")
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   0.000   0.002   0.035   0.006   0.500    3860

tmp$MAF[which(is.na(tmp$MAF))]<-0
# set of variants in chrX with only one study with samples contributing to the analysis


## Define the groups based on the list of variants we will use in fine mapping instead, and reassing all variants to these groups:
dat1<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier_1_metaanalysis_sorted.bed",sep=""),head=F)
tmp1<-tmp[which(tmp$MarkerName %in% dat1$V4),]

tmp1$maf_bin2 <- as.numeric(cut2(tmp1$MAF, g=20))
table(tmp1$maf_bin2)
#       1       2       3       4       5       6       7       8       9      10 
# 1556046 1428105 1496644 1494681 1484859 1492069 1492065 1492067 1492077 1492056 
#      11      12      13      14      15      16      17      18      19      20 
# 1492067 1492478 1491656 1492067 1492067 1492067 1492067 1492067 1492067 1492066 

tapply(tmp1$MAF, tmp1$maf_bin2, summary)
# # $`1`
# # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# # 0.001316 0.003219 0.003868 0.003783 0.004427 0.004958 
# # 
# # $`2`
# # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# # 0.004958 0.005486 0.006018 0.006035 0.006576 0.007172 
# # 
# # $`3`
# # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# # 0.007172 0.007801 0.008482 0.008512 0.009208 0.009975 
# # 
# # $`4`
# # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# # 0.009975 0.010805 0.011709 0.011765 0.012697 0.013779 
# # 
# # $`5`
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.01378 0.01494 0.01620 0.01627 0.01756 0.01901 
# # 
# # $`6`
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.01901 0.02063 0.02241 0.02252 0.02435 0.02645 
# # 
# # $`7`
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.02645 0.02865 0.03112 0.03123 0.03374 0.03652 
# # 
# # $`8`
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.03652 0.03951 0.04273 0.04293 0.04627 0.05000 
# # 
# # $`9`
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.05000 0.05398 0.05826 0.05834 0.06263 0.06724 
# # 
# # $`10`
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.06724 0.07217 0.07750 0.07764 0.08302 0.08876 
# # 
# # $`11`
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.08876 0.09484 0.10118 0.10133 0.10776 0.11447 
# # 
# # $`12`
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.1145  0.1215  0.1288  0.1289  0.1363  0.1441 
# # 
# # $`13`
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.1441  0.1522  0.1604  0.1607  0.1691  0.1780 
# # 
# # $`14`
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.1780  0.1870  0.1962  0.1964  0.2057  0.2153 
# # 
# # $`15`
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.2153  0.2254  0.2353  0.2356  0.2459  0.2566 
# # 
# # $`16`
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.2566  0.2675  0.2785  0.2786  0.2896  0.3012 
# # 
# # $`17`
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.3012  0.3130  0.3245  0.3245  0.3360  0.3480 
# # 
# # $`18`
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.3480  0.3602  0.3723  0.3725  0.3848  0.3972 
# # 
# # $`19`
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.3972  0.4097  0.4225  0.4226  0.4354  0.4483 
# # 
# # $`20`
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.4483  0.4612  0.4741  0.4740  0.4869  0.5000


tmp$maf_bin2<-NA

for (i in 1:20) {

    val_min<-min(tmp1$MAF[which(tmp1$maf_bin2==i)])
    val_max<-max(tmp1$MAF[which(tmp1$maf_bin2==i)])

    if(i==1) {
        tmp$maf_bin2[which(tmp$MAF<=val_max)]<-i
    } else {
        tmp$maf_bin2[which(tmp$MAF>=val_min & tmp$MAF<=val_max)]<-i
    }

}

summary(tmp$MAF[which(is.na(tmp$maf_bin2))])
tapply(tmp$MAF, tmp$maf_bin2, summary)

tmp[which(is.na(tmp$maf_bin2)),]
#                Variant     V1        V2        V3 maf_bin2
#                 <char> <char>     <int>     <int>    <num>
# 1:  chr11:74295061:C:T  chr11  74295061  74295062       NA
# 2: chr12:125574419:C:T  chr12 125574419 125574420       NA
# 3:    chr19:613631:A:C  chr19    613631    613632       NA
# 4:    chr19:613632:T:C  chr19    613632    613633       NA

tmp$maf_bin2[which(tmp$MarkerName %in% c("chr1:190483560:G:GA","chr11:74295061:C:T","chr12:125574419:C:T","chr12:6429051:T:G"))]<-17
tmp$maf_bin2[which(tmp$MarkerName %in% c("chr19:613631:A:C","chr19:613632:T:C"))]<-7

table(tmp$maf_bin2,useNA="ifany")
#        1        2        3        4        5        6        7        8 
# 10101540  1488659  2167197  2116201  2399054  2365528  2279943  2119603 
#        9       10       11       12       13       14       15       16 
#  2321466  2334237  2263001  2316297  2120551  2213719  1975021  1807211 
#       17       18       19       20 
#  1724031  1727810  1728808  1713733 


dat<-merge(dat,tmp[,c("MarkerName","maf_bin2")],by.x="V4",by.y="MarkerName",sort=F)
rm(dat1,dat2,tmp)
colnames(dat)[1]<-"Variant"

dat<-as.data.frame(dat)

for (i in 1:20) {
  dat[,paste("maf_bin_",i,sep="")]<-0
  dat[which(dat$maf_bin2==i),ncol(dat)]<-1
  table(dat[,ncol(dat)])
  
  fwrite(dat[,c(1,ncol(dat))],
              paste(path,"post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/allchr_",paste("maf_bin_",i,sep=""),"_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed",sep=""),
              col.names=T,row.names=F,quote=F,sep="\t")
}

q("no")



##########################################################################################################################################
##########################################################################################################################################

# add new consequences from vep_111


rm *vep*

MEM=2000
for chr in {1..22} X
do
bsub -J"annota" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q long \
-e ${path_gwas}post_imputation/log/${chr}_annotation_vs2_vep_sanger_stderr \
-o ${path_gwas}post_imputation/log/${chr}_annotation_vs2_vep_sanger_stdout \
"/path/to/project \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_sorted.vcf.gz"
done

# CONTINUE HERE!!!!


for chr in {1..22} X
do
echo ${chr} &&  tail -100 ${path_gwas}post_imputation/log/${chr}_annotation_vep_sanger_stdout | grep "Successfully completed."
done

for chr in {1..22} X
do
echo ${chr} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_sorted.vcf.gz | sed "/^#/d" | wc -l
done
# # ok sum = 49283610

# use bcftools +split-vep plugin to extract information"
#Version: 1.12-57-g0c2765b - updated bftools version




mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/*vep.* vep_111.1/

# index the files
path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=800

for chr in 22
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 2 \
-e ${path_gwas}post_imputation/log/${chr}_vep_index_stderr \
-o ${path_gwas}post_imputation/log/${chr}_vep_index_stdout \
"bcftools index -f ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/vep/vep_111.1/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_sorted.vep.vcf.gz --threads 8"
done

for chr in {1..22} X
do 
echo ${chr} && tail -50 ${path_gwas}post_imputation/log/${chr}_vep_index_stdout | grep "Successfully"
done

