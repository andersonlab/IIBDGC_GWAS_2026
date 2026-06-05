# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
######################################################################################################################################################################
######################################################################################################################################################################
### MANUALLY ADD - partially harmonised by GWAS catalog:

# 5. Successfully harmonised variants

# 97.78% ( 7790531 of 7967643 ) sites successfully harmonised.
# hm_code	Number	Percentage	Explanation
# 10	1372090	17.22%	Forward strand; Correct orientation; Already harmonised
# 11	5285453	66.34%	Forward strand; Flipped orientation; Requires harmonisation
# 12	1488	0.02%	Reverse strand; Correct orientation; Already harmonised
# 13	4252	0.05%	Reverse strand; Flipped orientation; Requires harmonisation
# 5	233936	2.94%	Palindromic; Assume forward strand; Correct orientation; Already harmonised
# 6	893312	11.21%	Palindromic; Assume forward strand; Flipped orientation; Requires harmonisation

# path_gwas="/path/to/ibdgwas/IIBDGC/"
path_gwas="/path/to/ibdgwas/IIBDGC/"

cd /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/

# download clonal hematopoyesis:
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476007/GCST90476007.tsv.gz
# in b38

gwas_id="GCST90476007"

# singularity exec iibdgc_postprocess_10_singularity.sif
# singularity exec polyfun_2_singularity.sif

# harmonise nomenclature based on the instructructions provided by GWAS catalog

# MEM=45000
# bsub -Is -M"$MEM" -R"select[model==Intel_Platinum && mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group R 

library(data.table)
library(R.utils)
library(tidyr)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

gwas_id<-"GCST90476007"

df<-fread(paste0(path_gwas,"resources/gwas_summary_statistics/",gwas_id,".tsv.gz"),)
dim(df)
#[1] 19709213 6

df$variant<-paste("chr",df$chromosome,":",df$base_pair_location,":",df$other_allele,":",df$effect_allele,sep="")
df<-df[which(df$effect_allele!=df$other_allele),]

# keep only SNPs
df<-df[which(df$effect_allele %in% c("A","T","C","G") & df$other_allele %in% c("A","T","C","G"))]
dim(df)
# [1] 18690322       13

df<-df[which(!is.na(df$base_pair_location)),]
dim(df)
# 18690304

bed<-df[,c("chromosome","base_pair_location","variant","other_allele","effect_allele")]

bed<-bed[!duplicated(bed),]
dim(bed)
#[1]   18685916      5

bed$Pos<-format(bed$base_pair_location, scientific=F)
bed$chr<-paste("chr",bed$chromosome,sep="")
setorder(bed,cols="chr","base_pair_location")

dim(bed)
# [1]  18685916      7


# create tfam and tped like file to convert to vcf file
bed$morgan<-0
bed$genotype<-paste(bed$other_allele,bed$effect_allele)

tped<-bed[,c("chr","variant","morgan","Pos","genotype")]

file_out_tped<-paste(path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_b38.tped",sep="")
write.table(tped,file_out_tped,
            col.names=F,row.names=F,sep="\t",quote=F)

tfam<-c("ID_1","ID_1","0","0","1","1")
tfam<-t(as.data.frame(tfam))
file_out_tfam<-paste(path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_b38.tfam",sep="")
write.table(tfam,file_out_tfam,
            col.names=F,row.names=F,sep="\t",quote=F)
  
# create A1 allele
file_out_allele<-paste(path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_Ref",sep="")
write.table(bed[,c("variant","other_allele")],file_out_allele,
            col.names=F,row.names=F,sep="\t",quote=F)



# create vcf:
system(paste("/path/to/software/username/plink_linux_x86_64_20181202/./plink ",
             "--tfile ",path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_b38 ",
             "--allow-no-sex ",
             "--a2-allele ",path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_Ref ",
             "--keep-allele-order --output-chr M --recode vcf-iid ",
             "--out ",path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_b38",sep=""))

# # Note that most PLINK analyses treat the A1 (usually minor) allele as the reference allele, which makes sense when only biallelic variants are involved.
# # However, since it is conventional for VCF files to set the major allele as the reference allele instead

# # # double check alleles, variants:

system(paste0("bcftools +fixref ",path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_b38.vcf ",
             "-Oz -o ",path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_b38.vcf.gz -- -f /path/to/ibdgwas/IIBDGC/resources/hg38/hg38_edited.fa -m top"))
# # SC, guessed strand convention
# SC      TOP-compatible  0
# SC      BOT-compatible  0
# # ST, substitution types
# ST      A>C     832256  4.5%
# ST      A>G     3824793 20.5%
# ST      A>T     635128  3.4%
# ST      C>A     671134  3.6%
# ST      C>G     775876  4.2%
# ST      C>T     2602317 13.9%
# ST      G>A     2608887 14.0%
# ST      G>C     777406  4.2%
# ST      G>T     672601  3.6%
# ST      T>A     634263  3.4%
# ST      T>C     3818879 20.4%
# ST      T>G     832376  4.5%
# # NS, Number of sites:
# NS      total           18685916
# NS      ref match       0       0.0%
# NS      ref mismatch    18685916        100.0%
# NS      flipped         1411927 7.6%
# NS      swapped         17273885        92.4%
# NS      flip+swap       0       0.0%
# NS      unresolved      104     0.0%
# NS      fixed pos       0       0.0%
# NS      errors          0
# NS      skipped         0
# NS      non-ACGT        0
# NS      non-SNP         0
# NS      non-biallelic   0


# # Double check conversion:
system(paste0("bcftools +fixref ",path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_b38.vcf.gz -- -f /path/to/ibdgwas/IIBDGC/resources/hg38/hg38_edited.fa"))
# # SC, guessed strand convention
# SC      TOP-compatible  0
# SC      BOT-compatible  0
# # ST, substitution types
# ST      A>C     671134  3.6%
# ST      A>G     2608887 14.0%
# ST      A>T     634283  3.4%
# ST      C>A     832256  4.5%
# ST      C>G     777404  4.2%
# ST      C>T     3818879 20.4%
# ST      G>A     3824793 20.5%
# ST      G>C     775878  4.2%
# ST      G>T     832376  4.5%
# ST      T>A     635108  3.4%
# ST      T>C     2602317 13.9%
# ST      T>G     672601  3.6%
# # NS, Number of sites:
# NS      total           18685916
# NS      ref match       18685812        100.0%
# NS      ref mismatch    104     0.0%
# NS      errors          0
# NS      skipped         0
# NS      non-ACGT        0
# NS      non-SNP         0
# NS      non-biallelic   0


 #### VCF to BED
system(paste0("/path/to/software/username/plink_linux_x86_64_20181202/./plink --vcf ",path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_b38.vcf.gz --keep-allele-order --allow-no-sex --double-id --make-bed --out ",path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_b38_2"))


# # ##############################################################
# # # 16.2 UPDATE NAME VARIANTS TO CHR:POSITION_REF_ALT in b38

system(paste0("zcat ",path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_b38.vcf.gz | cut -f '1-5' | awk '{print $3,$1\":\"$2\"_\"$4\"_\"$5}' > ",path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_list_variants_b38"))

########################


bim<-fread(paste0(path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_b38_2.bim"),head=F)
colnames(bim)[2]<-"variant"

ids<-fread(paste0(path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_list_variants_b38"),head=F)
ids<-ids[-(1:2),]
colnames(ids)<-c("variant","variant_b38")

bim<-merge(bim,ids,by="variant",all.y=T)
bim<-bim[,c("variant","variant_b38")]

# update alleles
bim$A1<-gsub(".*_","",bim$variant_b38)
bim$A0<-gsub("[0-9]*:[0-9]*_","",bim$variant_b38)
bim$A0<-gsub("_[A-Z]*","",bim$A0)
head(bim)
head(bed)

bed<-merge(bed,bim,by="variant",all.y=T)
dim(bed)

rm(bim,tped,tfam)

### add freq from gnomad:

# # GNOMAD:
# system(paste0("awk 'NR==FNR{vals[$2];next} ($1) in vals' ",path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_list_variants_b38 /path/to/ibdgwas/IIBDGC/resources/gnomad/gnomad_freq_edited > ",path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_list_variants_b38_gnomad"))
# system(paste0("gzip ",path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_list_variants_b38_gnomad"))


##############

# No need to add frequency"

# gnomad<-fread(paste0(path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_list_variants_b38_gnomad.gz"),head=F,sep=" ")
# gnomad<-separate_wider_delim(gnomad, cols = "V1", delim = "\t", names = c("SNP", "CHROM"))
# gnomad<-as.data.frame(gnomad)

# colnames(gnomad)<-c("SNP","CHROM","POS","REF","ALT","AF","AF_nfe","AF_afr","AF_amr","AF_eas","AF_sas","AF_asj")
# gnomad$AF_nfe<-as.numeric(as.character(gnomad$AF_nfe))

# bed<-merge(bed,gnomad[,c("SNP","AF_nfe")],by.x="variant_b38",by.y="SNP",all.x=T,sort=F)
# rm(gnomad)


dim(df[which(df$variant %in% bed$variant),])
# [1] 18690304       10

dim(bed)
# [1] 18685916       13

head(bed)
head(df)

df$beta<-log(df$odds_ratio)
df$SE<-(log(df$ci_upper)-log(df$ci_lower))/3.919928

# some beta set up to 0 - exclude
df<-df[which(df$beta!=0),]

df<-merge(df[,c("variant","beta","SE","p_value","effect_allele_frequency","r2")],bed,by="variant",all.y=T)
rm(bed)

df$N<-438632
df$INFO<-df$r2

# only retain those with info>0.4 and MAF >0.001:
df<-df[which(df$INFO>=0.4 & df$effect_allele_frequency>=0.001 & df$effect_allele_frequency<=0.999),]

length(df$beta[which(df$other_allele==df$A1 & df$effect_allele==df$A0)])
# 12825740
length(df$beta[which(df$effect_allele==df$A1 & df$other_allele==df$A0)])
# 86

# exclude those 86:
df<-df[which(df$other_allele==df$A1 & df$effect_allele==df$A0),]

dim(df)
# [1] 12713680       21

df<-df[,c("variant_b38","chromosome","base_pair_location","A0","A1","beta","SE","p_value","INFO","N","effect_allele_frequency")]

df<-df[which(!is.na(df$effect_allele_frequency)),]

colnames(df)<-c("SNP","CHR","BP","A0","A1","BETA","SE","PVALUE","INFO","N","A1FREQ")
dim(df)
# [1] 12713680       11

df$SNP<-gsub("_",":",df$SNP)
df$SNP<-paste("chr",df$SNP,sep="")

# remove duplicated ids:
dups<-df$SNP[duplicated(df$SNP)]
length(dups)
# [1] 3032
df<-df[!df$SNP %in% dups,]
dim(df)
# [1] 12707616       11

fwrite(df,paste(path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_edited.tsv.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

q("no")

###########################################################################################################################
###########################################################################################################################


MEM=8200

for i in 0
do
bsub -J"subset_gwas" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/subset_gwas_variants_${gwas_id[i]}_stdout \
-e ${path_gwas}post_imputation/2022/log/subset_gwas_variants_${gwas_id[i]}_stderr \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/subset_gwas_variants_for_genetic_correlation.R ${gwas_id[i]} > \
${path_gwas}post_imputation/2022/log/subset_gwas_variants_for_genetic_correlation_${gwas_id[i]}.Rout"
done


for i in  0
do
echo ${i} && tail -50 ${path_gwas}post_imputation/2022/log/subset_gwas_variants_${gwas_id[i]}_stdout  | grep -E "Successfully|Exit"
done

for i in 0
do
echo ${i} && ls -la ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_edited_2.tsv.gz
done

for i in 0
do
cat ${path_gwas}post_imputation/2022/log/subset_gwas_variants_for_genetic_correlation_${gwas_id[i]}.Rout
done

#######################################################################################################
# reformat to munged file

MEM=7000

for i in 0
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${gwas_id[i]}_ldsc_sumstats_only_SNPs_polyfun_mungen_format_stderr \
-o ${path_gwas}post_imputation/log/${gwas_id[i]}_ldsc_sumstats_only_SNPs_polyfun_mungen_format_stdout \
"munge_polyfun_sumstats.py \
--sumstats ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_edited_2.tsv.gz \
--min-info 0 \
--min-maf 0 \
--remove-strand-ambig \
--out ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_sumstats_munged.parquet"
done

# continue here


for i in 0
do
echo ${i} && tail -50 ${path_gwas}post_imputation/log/${gwas_id[i]}_ldsc_sumstats_only_SNPs_polyfun_mungen_format_stdout  | grep -E "Successfully|Exit"
done

for i in 0
do
echo ${i} && ls -la ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_sumstats_munged.parquet
done


# SAVE PARQUET FILE AS FLAT TEXT FILE AS INPUT FOR LDSC/LDSC.PY

module unload HGI/softpack/groups/team152/iibdgc_postprocess/10
module unload HGI/softpack/groups/team152/polyfun-2/1

MEM=3000
for i in 0
do
bsub -J"subset_gwas" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/munge_parquet_to_sumstats_${gwas_id[i]}_stdout \
-e ${path_gwas}post_imputation/2022/log/munge_parquet_to_sumstats_${gwas_id[i]}_stderr \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/reformat_munge_parquet_to_sumtats_for_genetic_correlation.R ${gwas_id[i]} > \
${path_gwas}post_imputation/2022/log/reformat_munge_parquet_to_sumtats_for_genetic_correlation_${gwas_id[i]}.Rout"
done


for i in 0
do
echo ${i} && tail -50 ${path_gwas}post_imputation/2022/log/munge_parquet_to_sumstats_${gwas_id[i]}_stdout | grep -E "Successfully|Exit"
done

for i in 0
do
echo ${i} && ls -la ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_sumstats_munged.sumstats
done

# COMPLETED - JUNE 25

