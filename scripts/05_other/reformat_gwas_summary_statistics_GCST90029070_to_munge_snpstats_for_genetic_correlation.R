# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
######################################################################################################################################################################
######################################################################################################################################################################
### MANUALLY ADD - CRP

# Saredo Said, Raha Pazoki, Ville Karhunen, Urmo Võsa, Symen Ligthart, Barbra Bodinier, Fotios Koskeridis, Paul Welsh,
#  Inflammation working group of CHARGE Consortium, Behrooz Z Alizadeh, Daniel I Chasman, Naveed Sattar, Marc Chadeau-Hyam, Evangelos Evangelou, Marjo-Riitta Jarvelin, 
# Paul Elliott, Ioanna Tzoulaki, Abbas Dehghan. 


# path_gwas="/path/to/ibdgwas/IIBDGC/"
path_gwas="/path/to/ibdgwas/IIBDGC/"
gwas_id="GCST90029070"

cd /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/

wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90029001-GCST90030000/GCST90029070/GCST90029070_buildGRCh37.tsv.gz

# singularity exec iibdgc_postprocess_10_singularity.sif
# singularity exec polyfun_2_singularity.sif


# harmonise nomenclature and liftover to b38:

# MEM=35000
# bsub -Is -M"$MEM" -R"select[model==Intel_Platinum && mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group R \

library(data.table)
library(R.utils)
library(tidyr)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

gwas_id<-"GCST90029070"

df<-fread(paste0(path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_buildGRCh37.tsv.gz"),)

dim(df)
#[1] 11106737 6

df$variant<-paste("chr",df$chromosome,":",df$base_pair_location,":",df$other_allele,":",df$effect_allele,sep="")
df<-df[which(df$effect_allele!=df$other_allele),]
dim(df)
#[1][1] 11106737       10 

# keep only SNPs
df<-df[which(df$effect_allele %in% c("A","T","C","G") & df$other_allele %in% c("A","T","C","G"))]
dim(df)
# [1] 9973021       13

bed<-df[,c("chromosome","base_pair_location","variant","other_allele","effect_allele")]

bed<-bed[!duplicated(bed),]
dim(bed)
#[1] 9973021        5

bed$end<-bed$base_pair_location+1

bed$Pos<-format(bed$base_pair_location, scientific=F)
bed$end<-format(bed$end, scientific=F)

bed$chr<-paste("chr",bed$chromosome,sep="")
setorder(bed,cols="chr","base_pair_location")

file_out<-paste(path_gwas,"/resources/gwas_summary_statistics/",gwas_id,"_buildGRCh37_b37.bed",sep="")
write.table(bed[,c("chr","base_pair_location","end","variant")],file_out,col.names=F,row.names=F,quote=F)


#### lift positions
system(paste("/path/to/software/username/liftOver ",
             file_out," /path/to/ibdgwas/IIBDGC/previous_qced_b38/liftover/hg19ToHg38.over.chain.gz ",
             file_out,"_lifted_hg38 ",file_out,"_no_lifted_hg38",sep=""))

system(paste("cut -f 4 ",file_out,"_no_lifted_hg38 | sed '/^#/d' > ",file_out,"_no_lifted_hg38_variants_to_exclude_tmp",sep=""))
system(paste("grep alt ",file_out,"_lifted_hg38 | cut -f 4 | cat - ",file_out,"_no_lifted_hg38_variants_to_exclude_tmp > ",
             file_out,"_no_lifted_hg38_variants_to_exclude",sep=""))

system(paste("wc -l ",file_out,"_no_lifted_hg38_variants_to_exclude",sep=""))
# 4627 /path/to/ibdgwas/IIBDGC//resources/gwas_summary_statistics/GCST90029070_buildGRCh37_b37.bed_no_lifted_hg38_variants_to_exclude


### exclude non lifted:

bed_up<-fread(paste(file_out,"_lifted_hg38",sep=""),head=F)
colnames(bed_up)[2]<-"pos_b38"
colnames(bed_up)[4]<-"variant"

dim(bed)
# [1] 9973021       7
bed<-merge(bed,bed_up[,c("variant","pos_b38")],by="variant",all.y=T,sort=F)
dim(bed)
# [1] 9970115       8


# create tfam and tped like file to convert to vcf file
bed$morgan<-0
bed$genotype<-paste(bed$other_allele,bed$effect_allele)

tped<-bed[,c("chr","variant","morgan","pos_b38","genotype")]

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
# ST      A>C     0       0.0%
# ST      A>G     0       0.0%
# ST      A>T     0       0.0%
# ST      C>A     823783  8.3%
# ST      C>G     0       0.0%
# ST      C>T     3385805 34.0%
# ST      G>A     3389701 34.0%
# ST      G>C     846711  8.5%
# ST      G>T     822129  8.2%
# ST      T>A     701484  7.0%
# ST      T>C     0       0.0%
# ST      T>G     0       0.0%
# NS, Number of sites:
# NS      total           9970115
# NS      ref match       5461767 54.8%
# NS      ref mismatch    4507846 45.2%
# NS      flipped         414715  4.2%
# NS      swapped         4085076 41.0%
# NS      flip+swap       368742  3.7%
# NS      unresolved      2082    0.0%
# NS      fixed pos       0       0.0%
# NS      errors          0
# NS      skipped         502
# NS      non-ACGT        502
# NS      non-SNP         0
# NS      non-biallelic   0


# # Double check:
system(paste0("bcftools +fixref ",path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_b38.vcf.gz -- -f /path/to/ibdgwas/IIBDGC/resources/hg38/hg38_edited.fa"))

# # SC, guessed strand convention
# SC      TOP-compatible  0
# SC      BOT-compatible  0
# # ST, substitution types
# ST      A>C     381585  3.8%
# ST      A>G     1482779 14.9%
# ST      A>T     350564  3.5%
# ST      C>A     442132  4.4%
# ST      C>G     422799  4.2%
# ST      C>T     1903695 19.1%
# ST      G>A     1906962 19.1%
# ST      G>C     423912  4.3%
# ST      G>T     440707  4.4%
# ST      T>A     350920  3.5%
# ST      T>C     1482070 14.9%
# ST      T>G     381488  3.8%
# # NS, Number of sites:
# NS      total           9970115
# NS      ref match       9967567 100.0%
# NS      ref mismatch    2046    0.0%
# NS      errors          0
# NS      skipped         502
# NS      non-ACGT        502
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

# GNOMAD:
system(paste0("awk 'NR==FNR{vals[$2];next} ($1) in vals' ",path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_list_variants_b38 /path/to/ibdgwas/IIBDGC/resources/gnomad/gnomad_freq_edited > ",path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_list_variants_b38_gnomad"))
system(paste0("gzip ",path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_list_variants_b38_gnomad"))


##############

gnomad<-fread(paste0(path_gwas,"resources/gwas_summary_statistics/list_variants_",gwas_id,"_b38_gnomad.gz"),head=F,sep=" ")
gnomad<-separate_wider_delim(gnomad, cols = "V1", delim = "\t", names = c("SNP", "CHROM"))
gnomad<-as.data.frame(gnomad)

colnames(gnomad)<-c("SNP","CHROM","POS","REF","ALT","AF","AF_nfe","AF_afr","AF_amr","AF_eas","AF_sas","AF_asj")
gnomad$AF_nfe<-as.numeric(as.character(gnomad$AF_nfe))

bed<-merge(bed,gnomad[,c("SNP","AF_nfe")],by.x="variant_b38",by.y="SNP",all.x=T,sort=F)
rm(gnomad)

dim(df[which(df$variant %in% bed$variant),])
# [1] 9970115       10

dim(bed)
# [1] 9970115       13

head(bed)
head(df)

df<-merge(df[,c("variant","beta","standard_error","p_value")],bed,by="variant",all.y=T)
rm(bed)

df$N<-575531
df$INFO<-1


length(df$beta[which(df$other_allele==df$A1 & df$effect_allele==df$A0)])
# 4495276
length(df$beta[which(df$effect_allele==df$A1 & df$other_allele==df$A0)])
# 5464315
4495276+5464315

dim(df)
# [1] 11102002       21

length(df$beta[which(df$effect_allele!=df$A1 & df$other_allele!=df$A0 & df$effect_allele!=df$A0 & df$other_allele!=df$A1)])
# [1] 10524
# exclude those 10K
df<-df[which(!(df$effect_allele!=df$A1 & df$other_allele!=df$A0 & df$effect_allele!=df$A0 & df$other_allele!=df$A1)),]

df$beta_ed<-df$beta
df$beta_ed[which(df$other_allele==df$A1 & df$effect_allele==df$A0)]<-df$beta[which(df$other_allele==df$A1 & df$effect_allele==df$A0)]*-1


# length(df$beta_ed[which(df$other_allele==df$A1 & df$effect_allele==df$A0)])

df<-df[,c("variant_b38","chromosome","pos_b38","A0","A1","beta_ed","standard_error","p_value","INFO","N","AF_nfe")]

# exclude those with no AF
dim(df[which(is.na(df$AF_nfe)),])
# [1] 83208      11

df<-df[which(!is.na(df$AF_nfe)),]

colnames(df)<-c("SNP","CHR","BP","A0","A1","BETA","SE","PVALUE","INFO","N","A1FREQ")
df<-df[which(df$A1FREQ>=0.001 & df$A1FREQ<=0.999),]
dim(df)
# [1] 9551888       11

df$SNP<-gsub("_",":",df$SNP)
df$SNP<-paste("chr",df$SNP,sep="")

# remove duplicated ids:
dups<-df$SNP[duplicated(df$SNP)]
length(dups)
# [1] 107
df<-df[!df$SNP %in% dups,]
dim(df)
# [1] 9551675       11

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

# CONTINUE HERE

for i in 0
do
echo ${i} && tail -50 ${path_gwas}post_imputation/2022/log/munge_parquet_to_sumstats_${gwas_id[i]}_stdout | grep -E "Successfully|Exit"
done

for i in 0
do
echo ${i} && ls -la ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_sumstats_munged.sumstats
done

# COMPLETED - JUNE 25

