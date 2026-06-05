# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# edit frequency file provided by TOPMed followin:

# https://www.well.ox.ac.uk/~wrayner/tools/

# The TOPMed reference panel is not available for direct download from this site, it needs to be created from the VCF of dbSNP submitted sites (currently ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz). 
# This can be downloaded from the Bravo Website https://bravo.sph.umich.edu/freeze5/hg38/ 

## SAVED IN:
#   /path/to/ibdgwas/IIBDGC/resources/TOPMed

# Once downloaded the VCF can be converted to an HRC formatted reference legend using the code here: CreateTOPMed.zip 
# Usage: ./CreateTOPMed.pl -i ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz 
# 
# By default this will create a file filtered for variants flagged as PASS only, if you wish to use all variants the -a flag overrides this.
# To override the default output file naming use -o filename. 


# ls -la /path/to/ibdgwas/IIBDGC/resources/TOPMed/
# -rw-rw-r-- 1 username ibdgwas 6904593910 Jun 10 13:49 bravo-dbsnp-all.vcf.gz
# -rwxrwxrwx 1 username ibdgwas       2654 May 15 11:14 CreateTOPMed.pl
# -rw-rw-r-- 1 username ibdgwas       1298 May 18 14:22 CreateTOPMed.zip
# -rw-rw-r-- 1 username ibdgwas 3983735769 Jun 10 18:36 PASS.Variantsbravo-dbsnp-all.tab.gz
  
# use file to identify which variants will be available in TOPMEd:

path_gwas=/path/to/ibdgwas/IIBDGC/
  

zcat PASS.Variantsbravo-dbsnp-all.tab.gz | head -10
# #CHROM	POS	ID	REF	ALT	AC	AN	AF
# 1	10498	TOPMed_freeze_5?chr1:10,498'01	G	A	1	125568	7.96381e-06
# 1	10498	TOPMed_freeze_5?chr1:10,498'02	G	T	1	125568	7.96381e-06
# 1	10509	TOPMed_freeze_5?chr1:10,509	G	A	1	125568	7.96381e-06
# 1	10527	TOPMed_freeze_5?chr1:10,527	C	T	1	125568	7.96381e-06
# 1	10531	TOPMed_freeze_5?chr1:10,531	C	G	2	125568	1.59276e-05
# 1	10534	TOPMed_freeze_5?chr1:10,534	A	G	1	125568	7.96381e-06
# 1	10535	TOPMed_freeze_5?chr1:10,535	G	C	2	125568	1.59276e-05
# 1	10548	TOPMed_freeze_5?chr1:10,548'01	C	T	3	125568	2.38914e-05
# 1	10550	TOPMed_freeze_5?chr1:10,550	G	A	1	125568	7.96381e-06


zcat ${path_gwas}resources/TOPMed/PASS.Variantsbravo-dbsnp-all.tab.gz | awk -v OFS='\t' '{print $0,$1":"$2"_"$4"_"$5}' \
> ${path_gwas}resources/TOPMed/PASS.Variantsbravo-dbsnp-all_edited.vcf

head ${path_gwas}resources/TOPMed/PASS.Variantsbravo-dbsnp-all_edited.vcf
# #CHROM	POS	ID	REF	ALT	AC	AN	AF	#CHROM:POS_REF_ALT
# 1	10498	TOPMed_freeze_5?chr1:10,498'01	G	A	1	125568	7.96381e-06	1:10498_G_A
# 1	10498	TOPMed_freeze_5?chr1:10,498'02	G	T	1	125568	7.96381e-06	1:10498_G_T
# 1	10509	TOPMed_freeze_5?chr1:10,509	G	A	1	125568	7.96381e-06	1:10509_G_A
# 1	10527	TOPMed_freeze_5?chr1:10,527	C	T	1	125568	7.96381e-06	1:10527_C_T
# 1	10531	TOPMed_freeze_5?chr1:10,531	C	G	2	125568	1.59276e-05	1:10531_C_G
# 1	10534	TOPMed_freeze_5?chr1:10,534	A	G	1	125568	7.96381e-06	1:10534_A_G
# 1	10535	TOPMed_freeze_5?chr1:10,535	G	C	2	125568	1.59276e-05	1:10535_G_C
# 1	10548	TOPMed_freeze_5?chr1:10,548'01	C	T	3	125568	2.38914e-05	1:10548_C_T
# 1	10550	TOPMed_freeze_5?chr1:10,550	G	A	1	125568	7.96381e-06	1:10550_G_A

wc -l ${path_gwas}resources/TOPMed/PASS.Variantsbravo-dbsnp-all_edited.vcf 
# 463071133 PASS.Variantsbravo-dbsnp-all_edited.vcf


# edit chr name:

sed 's/X:/23:/g' ${path_gwas}resources/TOPMed/PASS.Variantsbravo-dbsnp-all_edited.vcf \
> ${path_gwas}resources/TOPMed/PASS.Variantsbravo-dbsnp-all_edited_2.vcf

awk -v OFS='\t' '{print $0,$1":"$2}' ${path_gwas}resources/TOPMed/PASS.Variantsbravo-dbsnp-all_edited_2.vcf \
> ${path_gwas}resources/TOPMed/PASS.Variantsbravo-dbsnp-all_edited_3.vcf



########## double check if some of the variants that have low ER2 in Kyle's analysis are in original TOPmed file:

# reformat Kyle's file:

### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"


# VARIANTS LOW R2
r2<-read.csv("/path/to/ibdgwas/IIBDGC/from_kyle/LowER2snpsNIDDK.csv",head=T)
r2$chr<-gsub(":.*","",r2$SNP)
r2$chr<-paste("chr",r2$chr,sep="")

r2$pos<-gsub("^[0-9]{1,2}:","",r2$SNP)
r2$pos<-gsub(":[A-Z]*:[A-Z]*$","",r2$pos)

r2$chrpos<-paste(r2$chr,r2$pos,sep=":")

write.table(r2[,c("chr","pos")],"/path/to/ibdgwas/IIBDGC/from_kyle/LowER2snpsNIDDK_edited",col.names=F,row.names=F,sep="\t",quote=F)

######

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=1500

bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/log/stderr_topmed_kyle \
-o ${path_gwas}pre_imputation/QC/log/stdout_topmed_kyle \
"/path/to/software/bcftools-1.9/./bcftools filter \
-T ${path_gwas}from_kyle/LowER2snpsNIDDK_edited \
${path_gwas}resources/TOPMed/bravo-dbsnp-all.vcf.gz \
-o ${path_gwas}from_kyle/LowER2snpsNIDDK_edited_inTOPMed"


wc -l ${path_gwas}from_kyle/LowER2snpsNIDDK_edited_inTOPMed
# 953 /path/to/ibdgwas/IIBDGC/from_kyle/LowER2snpsNIDDK_edited_inTOPMed


##### R

a<-read.table("/path/to/ibdgwas/IIBDGC/from_kyle/LowER2snpsNIDDK_edited_inTOPMed",sep="\t")
a$SNP<-paste(a$V1,a$V2,a$V4,a$V5,sep=":")
a$SNP<-gsub("chr","",a$SNP)
  
a$chrpos<-paste(a$V1,a$V2,sep=":")
a$chrpos<-gsub("chr","",a$chrpos)


dim(a[which(a$SNP %in% r2$SNP),])
# [1] 65 11

dim(a[which(a$chrpos %in% r2$chrpos),])
# [1] 379  12


head(a[which(!a$SNP %in% r2$SNP),])
# V1      V2                             V3 V4 V5  V6   V7
# 1 chr1 1450947 TOPMed_freeze_5?chr1:1,450,947  C  T 255 PASS
# 2 chr1 1766856 TOPMed_freeze_5?chr1:1,766,856  C  A 255 PASS
# 4 chr1 2757675 TOPMed_freeze_5?chr1:2,757,675  A  C 157  SVM
# 5 chr1 3072956 TOPMed_freeze_5?chr1:3,072,956  T  G 119  SVM
# 6 chr1 3630630 TOPMed_freeze_5?chr1:3,630,630  C  A 255 PASS
# 7 chr1 3657308 TOPMed_freeze_5?chr1:3,657,308  G  T 160 PASS
# 
# V10           SNP       chrpos
# 1 125568:0.000804345 1:1450947:C:T chr1:1450947
# 2 125568:3.18552e-05 1:1766856:C:A chr1:1766856
# 4 125568:3.18552e-05 1:2757675:A:C chr1:2757675
# 5 125568:7.96381e-06 1:3072956:T:G chr1:3072956
# 6 125568:7.96381e-06 1:3630630:C:A chr1:3630630
# 7 125568:7.96381e-06 1:3657308:G:T chr1:3657308


r2[which(r2$chrpos %in% c("chr1:1450947","chr1:1766856","chr1:2757675","chr1:3072956")),]
# SNP     MAF      R2     ER2  chr     pos       chrpos
# 669 1:1450947:A:C 0.01115 0.98933 0.47504 chr1 1450947 chr1:1450947






