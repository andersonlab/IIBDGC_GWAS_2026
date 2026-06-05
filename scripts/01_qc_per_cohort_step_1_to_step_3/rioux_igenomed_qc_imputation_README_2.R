# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
anvil_ccdg_broad_ai_ibd_daly_rioux_igenomed_gsa.txt.gz

# new batch of text files downloaded to /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/

cd /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/
md5sum anvil_ccdg_broad_ai_ibd_daly_rioux_igenomed_gsa.txt.gz
# 17b76da6a20e89aa4fbfc9a085304953  anvil_ccdg_broad_ai_ibd_daly_rioux_igenomed_gsa.txt.gz

less /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/file-manifest.json 
#  "md5sum": "17b76da6a20e89aa4fbfc9a085304953", - OK


zcat /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/anvil_ccdg_broad_ai_ibd_daly_rioux_igenomed_gsa.txt.gz \
| head -10

# create folders for the study:
mkdir /path/to/ibdgwas/IIBDGC/pre_imputation/raw/rioux_igenomed/
  
# create a tmp copy with final file ID, # to keep consistency in naming - only copy while files are still downloading (hpc-server client checks out 
# downloaded files to evaluate whether those have been downloaded), remove afterwards
cp /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/anvil_ccdg_broad_ai_ibd_daly_rioux_igenomed_gsa.txt.gz \
/path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/rioux_igenomed_gsa.txt.gz

############################################################################## 
# make a temporary sorted fcr file:

# farm5:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/sort_fcr.sh rioux_igenomed normal gsa
# Job <798560> is submitted to queue <normal>.

study=rioux_igenomed
less ${path_gwas}scripts/logs/fcrsort_stderr_${study} 

zcat /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/rioux_igenomed_gsa.txt.gz | wc -l
# 127973953
zcat /path/to/ibdgwas/IIBDGC/pre_imputation/raw/rioux_igenomed/rioux_igenomed_gsa_sorted.txt.gz | wc -l
# 127973953

##############################################################################  
# generate map, fam and lgen files from this:

#############
# MAP FILES #
#############

# MAP file describes a single marker and must contain exactly 4 columns:
# chromosome (1-22, X, Y or 0 if unplaced)
# rs# or snp identifier
# Genetic distance (morgans)
# Base-pair position (bp units)

# farm5:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_map_from_fcr.sh rioux_igenomed long gsa
# Job <840263> is submitted to queue <long>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_map_from_fcr* /path/to/user/scripts/IIBDGC/


#############
# FAM FILES #
#############

# FAM (6 first columns of PED) PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:
# Family ID
# Individual ID
# Paternal ID
# Maternal ID
# Sex (1=male; 2=female; other=unknown)
# Phenotype

# farm5:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_fam_from_fcr.sh rioux_igenomed normal gsa
#Job <840298> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_fam_from_fcr* /path/to/user/scripts/IIBDGC/


##############
# LGEN FILES #
##############

# LGEN file, test.lgen
# family ID
# individual ID
# snp ID
# allele 1 of this genotype
# allele 2 of this genotype

# farm5:
path_gwas=/path/to/ibdgwas/IIBDGC/
bash ${path_gwas}scripts/create_lgen_from_fcr.sh rioux_igenomed normal gsa
# Job <840316> is submitted to queue <normal>.

# keep updated back-up copy:
# cp /path/to/ibdgwas/IIBDGC/scripts/create_lgen_from_fcr* /path/to/user/scripts/IIBDGC/


###########################

wc -l /path/to/ibdgwas/IIBDGC/pre_imputation/raw/rioux_igenomed/*.fam
# 186 /path/to/ibdgwas/IIBDGC/pre_imputation/raw/rioux_igenomed/rioux_igenomed_gsa.fam

wc -l /path/to/ibdgwas/IIBDGC/pre_imputation/raw/rioux_igenomed/*.map
# 688032 /path/to/ibdgwas/IIBDGC/pre_imputation/raw/rioux_igenomed/rioux_igenomed_gsa.map

zcat /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/rioux_igenomed_gsa.txt.gz \
| awk 'END { print NR }'

# 127973953/688032
# [1] 186 OK

# remove initial copy:
rm /path/to/ibdgwas/IIBDGC/raw_data_from_data_commons/downloads_March2023/rioux_igenomed_gsa.txt.gz

##########################################################################################################################################
##########################################################################################################################################


##########################################################################################################################################
##########################################################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
mkdir ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/
mkdir ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/logs/
  
MEM=800

bsub -J"pl" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/logs/stderr_plink_0_rioux_igenomed_gsa \
-o ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/logs/stdout_plink_0_rioux_igenomed_gsa \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--lfile /path/to/ibdgwas/IIBDGC/pre_imputation/raw/rioux_igenomed/rioux_igenomed_gsa \
--missing-genotype 0 \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19"
# Job <699751> is submitted to queue <normal>.

tail -200 ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/logs/stdout_plink_0_rioux_igenomed_gsa
# 688032 variants and 186 people pass filters and QC.
# Note: No phenotypes present.


#######################################################################################################################################

#####################
# 1.- REMOVE INDELS #
#####################

### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

bim<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19.bim",sep=""),head=F)
table(bim$V5,bim$V6)
#

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(indels)
#[1] 136166     6

table(indels$V1)
# 

# X    X chromosome                    -> 23
# Y    Y chromosome                    -> 24
# XY   Pseudo-autosomal region of X    -> 25
# MT   Mitochondrial                   -> 26

mt<-bim[which(bim$V1 %in% c(0,26)),]
table(mt$V1)
# 26 
# 133 

chry<-bim[which(bim$V1 %in% c(24)),]
table(chry$V1)
# 24 
# 1473 

all_remove<-rbind(indels,mt)
all_remove<-all_remove[!duplicated(all_remove),]

table(bim$V1[which(bim$V2 %in% all_remove$V2)])
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 8504 8272 6512 4773 5162 5269 5225 4116 3779 4036 5308 4594 2555 2957 3228 3302 
# 17   18   19   20   21   22   23   24   26 
# 4090 1957 3402 2036 1216 1536 2618  744  133 

dim(all_remove)
#[1] 136229      6

write.table(all_remove[,"V2",drop=F],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_indel_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



#######################################################################################################################################

###################################################
# 2.- EDIT FAM FILES AND REMOVE NON PHENO SAMPLES #
###################################################

pheno<-fread("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes-4313226.csv.gz")
pheno<-as.data.frame(pheno)
fam<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19.fam",sep=""),head=F)

table(as.character(pheno[which(pheno$sample_id %in% fam$V1),"batch"]))
# anvil_ccdg_broad_ai_ibd_daly_rioux_igenomed_gsa 
# 186 

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]))
# Crohn's Disease
#   anvil_ccdg_broad_ai_ibd_daly_rioux_igenomed_gsa   4             140
#                                                  
#                                                   Indeterminate
#   anvil_ccdg_broad_ai_ibd_daly_rioux_igenomed_gsa             4
#                                                  
#                                                   Ulcerative Colitis
#   anvil_ccdg_broad_ai_ibd_daly_rioux_igenomed_gsa                 38

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                                                     Affected
# anvil_ccdg_broad_ai_ibd_daly_rioux_igenomed_gsa   4      182

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$control[which(pheno$sample_id %in% fam$V1)]))
#                                                   0
# anvil_ccdg_broad_ai_ibd_daly_rioux_igenomed_gsa 182

table(as.character(pheno$diag[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$affection[which(pheno$sample_id %in% fam$V1)]))
#                        Affected
#                      4        0
# Crohn's Disease      0      140
# Indeterminate        0        4
# Ulcerative Colitis   0       38

table(as.character(pheno$batch[which(pheno$sample_id %in% fam$V1)]),as.character(pheno$sex[which(pheno$sample_id %in% fam$V1)]),useNA="ifany")
#                                                    Female Male
# anvil_ccdg_broad_ai_ibd_daly_rioux_igenomed_gsa  4     95   87


################################################
# Affection status, by default, should be coded:
# -9 missing 
# 0 missing
# 1 unaffected
# 2 affected

fam[,"V6"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$control=="1")]),"V6"]<-"1"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$affection=="Affected")]),"V6"]<-"2"

table(fam$V6)
# 2  -9 
# 182   4 

################################################
# Sex (1=male; 2=female; other=unknown)

fam[,"V5"]<-"-9"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Female")]),"V5"]<-"2"
fam[which(fam$V1 %in% pheno$sample_id[which(pheno$sex=="Male")]),"V5"]<-"1"
table(fam$V5)
# 1  2 -9 
# 87 95  4 

fam_test<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19.fam",sep=""),head=F)
table(fam$V1==fam_test$V1)
# TRUE 
# 186

write.table(fam,paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_edited.fam",sep=""),col.names = F,row.names = F,quote = F,sep="\t")

################################################

samples_remove<-fam[which(fam$V6=="-9"),c(1,2)]
dim(samples_remove)
# [1] 4 2

write.table(samples_remove,paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_samples_nopheno_remove",sep=""),col.names = T,row.names = F,quote = F,sep="\t")


######

path_gwas=/path/to/ibdgwas/IIBDGC/
  
  /path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19.bed \
--bim ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19.bim \
--fam ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_edited.fam \
--exclude ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/list_indel_var_exclude \
--remove ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/list_samples_nopheno_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind
# 551803 variants and 182 people pass filters and QC.
# Among remaining phenotypes, 182 are cases and 0 are controls.



#######################################################################################################################################

#####################
# 3.- UPDATE TO B37 #
#####################

# data already in b37



#######################################################################################################################################

#########################
# 4.- ALIGN TO + STRAND #
#########################

################################################################################################
# 4.0 create new edited .bim file to force chr25 (chrX pseudoautosom) to be X for this exercise:


### R

bim<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind.bim",sep=""),head=F)
table(bim$V1)
# 

# X    X chromosome                    -> 23
# Y    Y chromosome                    -> 24
# XY   Pseudo-autosomal region of X    -> 25
# MT   Mitochondrial                   -> 26

bim$V1[which(bim$V1==25)]<-"23"

table(bim$V1)
# 

write.table(bim,paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_edited.bim",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


##########################################################
# 4.1 GENERATE VCF (output-chr M - reformat chr23 to X)

# intermediate step to sort chr23:
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind.bed \
--bim ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_edited.bim \
--fam ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind.fam \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_edited2

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_edited2 \
--allow-no-sex \
--output-chr M --recode vcf-iid --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19


###########################################################
# 4.2 FLIP ALLELES TO + STRAND:

path_gwas=/path/to/ibdgwas/IIBDGC/
export BCFTOOLS_PLUGINS=/path/to/software/bcftools-1.9/plugins
MEM=2000

bsub -J"bcf" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/logs/stderr_bcftools_1_rioux_igenomed_gsa \
-o ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/logs/stdout_bcftools_1_rioux_igenomed_gsa \
"/path/to/software/bcftools-1.9/./bcftools +fixref ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19.vcf \
-Oz -o ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_posstrandaligned.vcf.gz --threads 4 \
-- -f /path/to/project -m top"
# Job <705296> is submitted to queue <normal>.


#### CONTINUE HERE

# # SC, guessed strand convention
# SC	TOP-compatible	0
# SC	BOT-compatible	0
# # ST, substitution types
# ST	A>C	24536	4.3%
# ST	A>G	97703	17.1%
# ST	A>T	815	0.1%
# ST	C>A	28998	5.1%
# ST	C>G	1288	0.2%
# ST	C>T	132446	23.2%
# ST	G>A	132238	23.2%
# ST	G>C	1270	0.2%
# ST	G>T	28678	5.0%
# ST	T>A	873	0.2%
# ST	T>C	97832	17.1%
# ST	T>G	24496	4.3%
# # NS, Number of sites:
# NS	total        	571173
# NS	ref match    	484679	84.9%
# NS	ref mismatch 	86494	15.1%
# NS	flipped      	358	0.1%
# NS	swapped      	86136	15.1%
# NS	flip+swap    	1844	0.3%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0


###########################################################
# 4.3 DOUBLE-CHECK THE CONVERSION

/path/to/software/bcftools-1.9/./bcftools +fixref ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_posstrandaligned.vcf.gz \
-- -f /path/to/project

# # SC, guessed strand convention
# SC	TOP-compatible	0
# SC	BOT-compatible	0
# # ST, substitution types
# ST	A>C	24589	4.3%
# ST	A>G	99283	17.4%
# ST	A>T	823	0.1%
# ST	C>A	28945	5.1%
# ST	C>G	1301	0.2%
# ST	C>T	130820	22.9%
# ST	G>A	130658	22.9%
# ST	G>C	1257	0.2%
# ST	G>T	28658	5.0%
# ST	T>A	865	0.2%
# ST	T>C	99458	17.4%
# ST	T>G	24516	4.3%
# # NS, Number of sites:
# NS	total        	571173
# NS	ref match    	571173	100.0%
# NS	ref mismatch 	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0



###########################################################
# 4.4 VCF to BED

/path/to/software/plink_linux_x86_64_20181202/./plink \
--vcf ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_posstrandaligned.vcf.gz \
--double-id \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr \
--missing \
--out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr



########################################################################################################################################

################################################################################
# 5.- UPDATE NAME VARIANTS TO CHR:POSITION_REF_ALT; REMOVE DUPLICATED VARIANTS #
################################################################################

zcat ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_posstrandaligned.vcf.gz | cut -f '1-5' | awk '{print $3,$1":"$2"_"$4"_"$5}' \
> ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_variants_rioux_igenomed_gsa_hg19_posstrandaligned

###### R

bim<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr.bim",sep=""),sep="\t",head=F)
ids<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_variants_rioux_igenomed_gsa_hg19_posstrandaligned",sep=""),sep=" ",head=F
                ,skip=33)

varmiss<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr.lmiss",sep=""),head=T)

colnames(ids)[2]<-"ids"

table(bim$V2==ids$V1)
# TRUE 
# 571173 


#the major allele is set to A2 by default by Plink, keep ids with real ref/alt as in vcf using the ids file
bim.1<-cbind(bim,ids[,"ids",drop=F])

bim.1<-merge(bim.1,varmiss[,c("SNP","F_MISS")],by.x="V2",by.y="SNP",sort=F)
table(bim.1$V2==bim$V2)
# TRUE 
# 571173
bim.1$ids<-as.character(bim.1$ids)

# identify duplicated variants (same chr position ref and alt)
dups<-bim.1[which(duplicated(bim.1$ids)),"ids"]

length(dups)
#[1] 294

bim.1[which(duplicated(bim.1$V2)),]
#<0 rows> (or 0-length row.names) OK

for (i in 1:length(dups)){
  
  tmp<-bim.1[which(bim.1$ids %in% dups[i]),]
  
  keep<-tmp[which(tmp$F_MISS==min(tmp$F_MISS)),]
  
  if (nrow(keep)>1){
    keep<-keep[1,]
  }
  
  exclude<-tmp[which(!tmp$V2 %in% keep$V2),]
  
  bim.1$ids[which(bim.1$V2 %in% exclude$V2)]<-paste(bim.1$ids[which(bim.1$V2 %in% exclude$V2)],"_rm",sep="")
  
}

dim(bim.1[which(duplicated(bim.1$ids)),])
#[1] 4 8

bim.1[which(duplicated(bim.1$ids)),]
# 0

duplicated_variants<-bim.1[grep("_rm",bim.1$ids),"ids"]
length(duplicated_variants)
# [1] 295

write.table(duplicated_variants,paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_duplicated_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

table(bim.1$V2==bim$V2)
# TRUE 
# 571173 

table(bim.1$V1)
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 44189 45940 38771 36179 33388 40708 30993 28847 24310 28036 27792 26761 19855 
# 14    15    16    17    18    19    20    21    22    23    24 
# 18230 17103 18562 16672 16415 12952 13647  7985  8214 14926   698 

# NOTE chrx pseudoautosomal now in chr23, recode this

bim_old<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind.bim",sep=""),head=F)
table(bim_old$V1)
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 44189 45940 38771 36179 33388 40708 30993 28847 24310 28036 27792 26761 19855 
# 14    15    16    17    18    19    20    21    22    23    24 
# 18230 17103 18562 16672 16415 12952 13647  7985  8214 14926   698 

# https://www.ncbi.nlm.nih.gov/grc/human
# Name Chr Start Stop
# PAR#1	X	60,001	2,699,520
# PAR#2	X	154,931,044	155,260,560
# PAR#1	Y	10,001	2,649,520
# PAR#2	Y	59,034,050	59,363,566

dim(bim.1[which( (bim.1$V1==23) & ((bim.1$V4>=60001 & bim$V4<=2699520) | (bim.1$V4>=154931044 & bim$V4<=155260560)) ),])
# [1] 472   8

dim(bim.1[which( (bim.1$V1==23) & ((bim.1$V4>=60001 & bim$V4<=2699520) | (bim.1$V4>=154931044 & bim$V4<=155260560)) & (!bim.1$ids %in% duplicated_variants) ),])
#[1] 471   8

dim(bim.1[which( (bim.1$V1==24) & ((bim.1$V4>=10001 & bim$V4<=2649520) | (bim.1$V4>=59034050 & bim$V4<=59363566)) ),])
# 0 8

write.table(bim.1[,c(2,7,3:6)],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_edited.bim",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


##########

# NOTE: use pre-vcf fam file which keeps the ca/ctr info:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr.bed \
--bim ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_edited.bim \
--fam ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind.fam \
--exclude ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_duplicated_var_exclude \
--split-x 'b37' \
--make-bed --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup
# 570878 variants and 788 people pass filters and QC.
# Among remaining phenotypes, 788 are cases and 0 are controls.
# --split-x: 471 chromosome codes changed.



##############################################################################################################################################

##############################################
# 6.- COMPARE ALLELE FREQUENCIES WITH 1000GP #
##############################################

#########################################################################
# 6.1 COMPARE ALLELE FREQUENCIES WITH 1000GP, GENERATE LIST OF VARIANTS:

cat ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup.bim | cut -f 2 > ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/list_variants_rioux_igenomed_gsa_hg19_noind_posstr_nodup

##############################################
# 6.2 EXTRACT VARIANTS FROM 1000GP

# see script 1000gp_data_to_compare.R about how 1000GP data was generated

for chr in {1..24}; do /path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/1000gp/1000GP_EUR_chr${chr}_b37 \
--extract ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/list_variants_rioux_igenomed_gsa_hg19_noind_posstr_nodup \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/1000GP_EUR_chr${chr}_b37_rioux_igenomed_gsa_variants;done


##############################################
# 6.3 COMBINE FILES

#### R

dat<-matrix(ncol=3,nrow=23)
dat<-as.data.frame(dat)
for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/rioux_igenomed_gsa/1000GP_EUR_chr",i+1,"_b37_rioux_igenomed_gsa_variants.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/rioux_igenomed_gsa/1000GP_EUR_chr",i+1,"_b37_rioux_igenomed_gsa_variants.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/rioux_igenomed_gsa/1000GP_EUR_chr",i+1,"_b37_rioux_igenomed_gsa_variants.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_1000GP_files_rioux_igenomed_gsa_variants_tomerge.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/1000GP_EUR_chr1_b37_rioux_igenomed_gsa_variants \
--merge-list ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_1000GP_files_rioux_igenomed_gsa_variants_tomerge.txt \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/1000GP_EUR_b37_rioux_igenomed_gsa_variants
# 567549 variants and 503 people pass filters and QC.

# remove intermediate files:
rm ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/1000GP_EUR_chr*_b37_rioux_igenomed_gsa_variants*
  
  
  ##############################################
# 6.4 ESTIMATE FREQUENCIES:

# 1000GP

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/1000GP_EUR_b37_rioux_igenomed_gsa_variants \
--freq \
--out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/1000GP_EUR_b37_rioux_igenomed_gsa_variants

# rioux_igenomed_gsa 

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup \
--freq \
--allow-no-sex \
--out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup


#### R

gp<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/1000GP_EUR_b37_rioux_igenomed_gsa_variants.frq",sep=""),head=T)
g1<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup.frq",sep=""),head=T)

dim(g1)
#[1] 570878      6
dim(gp)
#[1] 567549      6

dim(g1[which(!g1$SNP %in% gp$SNP),])
#[1] 3329     6 #

colnames(gp)[3:6]<-paste(colnames(gp)[3:6],"_gp",sep="")
colnames(g1)[3:6]<-paste(colnames(g1)[3:6],"_g1",sep="")

all<-merge(g1,gp[,2:6],by="SNP",all=T)

check<-all[which(all$A1_g1!=all$A1_gp),]
dim(check)
#[1] 10857   10

# keep only A/T C/G
check<-check[which( (check$A1_g1=="G" & check$A2_g1=="C") | (check$A1_g1=="C" & check$A2_g1=="G") | (check$A1_g1=="A" & check$A2_g1=="T") | (check$A1_g1=="T" & check$A2_g1=="A")),]
dim(check)
#[1] 2063    10

# LIST OF VARIANTS TO REMOVE, WE CANNOT REALLY BE SURE WHETHER THERE IS STRAND ISSUE OR NOT
remove<-check[which(check$MAF_g1>=0.45),]
dim(remove)
#[1] 107 10 - remove


# LIST OF VARIANTS TO FLIP:
flip<-check[which(check$MAF_g1<0.45),]
dim(flip)
#[1] 1956 10

flip<-flip[order(flip$MAF_g1,decreasing=T),]

pdf(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/plot_maf_AT_CG_rioux_igenomed_gsa_1000gp.pdf",sep=""))
plot(flip$MAF_g1,flip$MAF_gp)
dev.off()

remove_2<-flip[which(flip$MAF_g1>0.2 & flip$MAF_gp<0.1),]
remove_3<-flip[which(flip$MAF_g1<0.1 & flip$MAF_gp>0.2),]

remove<-rbind(remove,remove_2,remove_3)
dim(remove)
#[1] 123 10

write.table(remove[,"SNP"],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_variants_to_remove_AT_CG",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

flip<-flip[which(!flip$SNP %in% remove$SNP),]
dim(flip)
#[1] 1940    10

write.table(flip[,"SNP"],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_variants_to_flip_AT_CG",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


##########################################################################################
# 6.5 REMOVE VARIANTS WE CANNOT BE SURE ARE IN THE RIGHT STRAND, AND FLIP THE OTHERS

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup \
--exclude ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_variants_to_remove_AT_CG \
--flip ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_variants_to_flip_AT_CG \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip
# --flip: 1940 SNPs flipped.
# 570755 variants and 788 people pass filters and QC.

########################################################################################################################################

########################################################################################################################################

####################
# 7.- CHECK GENDER #
####################

#######################################
# 7.1 CREATE LIST OF MALES AND FEMALES

########################################################################################################################################

####################
# 7.- CHECK GENDER #
####################

#######################################
# 7.1 CREATE LIST OF MALES AND FEMALES

###### R 

fam<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip.fam",sep=""),head=F)
table(fam$V5,useNA="ifany")
# 1    2 
# 404 384 

fam.male<-fam[which(fam$V5==1),1:2] 
colnames(fam.male)<-c("FID","IID")

fam.female<-fam[which(fam$V5==2),1:2] 
colnames(fam.female)<-c("FID","IID")

write.table(fam.male,paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_male_samples",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
write.table(fam.female,paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_female_samples",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

###################################
# 7.2 KEEP ONLY FEMALES AND CHR23

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip \
--keep ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_female_samples --chr 23 --make-bed \
--out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_check_sex_females_only

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_check_sex_females_only \
--missing --hardy --freq --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_check_sex_females_only

###### R

bim<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_check_sex_females_only.bim",sep=""),sep="\t",head=F)
table(bim$V1)
# 23 
# 14417

hwe<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_check_sex_females_only.hwe",sep=""),head=T)
hwe<-hwe[which(hwe$TEST=="UNAFF"),]

frq<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_check_sex_females_only.frq",sep=""),head=T)
var_miss<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_check_sex_females_only.lmiss",sep=""),head=T)

all<-merge(hwe[,c("SNP","P")],frq[,c("SNP","MAF")],by="SNP",sort=F)
all<-merge(all,var_miss[,c("SNP","F_MISS")],by="SNP",sort=F)

all<-all[which(all$MAF>0.05 & all$F_MISS<0.01 & all$P>1E-4),]
dim(all)
# [1] 10492    4

write.table(all,paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_good_chrX_variants",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


##################################
# 7.3 KEEP ONLY MALES AND CHR24

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip \
--keep ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_male_samples --chr 24 --make-bed \
--out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_check_sex_males_only

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_check_sex_males_only \
--freq --missing \
--out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_check_sex_males_only


#  AND FIND VARIANTS WITH NO CALLS IN MOST OF FEMALES

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip \
--keep ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_female_samples --chr 24 --recode A \
--out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_check_sex_females_only_chr24


###### R

## males

frq<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_check_sex_males_only.frq",sep=""),head=T)
var_miss<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_check_sex_males_only.lmiss",sep=""),head=T)

all<-merge(frq[,c("SNP","MAF","NCHROBS")],var_miss[,c("SNP","F_MISS")],by="SNP",sort=F)

## females, find variatns with less % of calls:
ped<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_check_sex_females_only_chr24.raw",sep=""),head=T,check.names=F)

dat<-matrix(nrow=nrow(all),ncol=2)
dat<-as.data.frame(dat)
colnames(dat)<-c("variant","percentage_NA")

for (i in 1:nrow(dat)) {
  tmp<-ped[,6+i,drop=F]
  dat$variant[i]<-gsub("_[A-Z]{1}$","",colnames(tmp))
  dat$percentage_NA[i]<-nrow(tmp[which(is.na(tmp)),,drop=F])/nrow(tmp)
}

all<-all[which(all$SNP %in% dat$variant[which(dat$percentage_NA>0.95)]),]


# keep variants with large number of calls in males
all<-all[which(all$NCHROBS>=max(all$NCHROBS,na.rm=T)-(max(all$NCHROBS,na.rm=T)*0.005)),]

dim(all)
# [1] 128   4

write.table(all,paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_good_chrY_variants",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

##############

cat ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_good_chrX_variants ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_good_chrY_variants > \
${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_good_chrXY_variants

wc -l ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_good_chrXY_variants
# 10620 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/list_good_chrXY_variants

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip \
--extract ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_good_chrXY_variants --make-bed \
--out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_check_sex_good_chrXY

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_check_sex_good_chrXY \
--check-sex ycount --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_check_sex_good_chrXY

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_check_sex_good_chrXY \
--check-sex --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_check_sex_good_chrX

# In 'ycount' mode, gender is still imputed from the X chromosome, but female calls are downgraded to ambiguous whenever more than 0 
# nonmissing Y genotypes are present, and male calls are downgraded when fewer than 0 are present. (Note that these are counts, not rates.) 
# These thresholds are controllable with --check-sex ycount's optional 3rd and 4th numeric parameters.


###### R

sex<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_check_sex_good_chrXY.sexcheck",sep=""),head=T)

sex$PEDSEX<-as.factor(sex$PEDSEX)

# X thresholds:

fmin_male<-mean(sex[which(sex$PEDSEX==1),"F"])-(4*sd(sex[which(sex$PEDSEX==1),"F"]))
fmin_male
# [1] 0.6574553
fmax_female<-mean(sex[which(sex$PEDSEX==2),"F"])+(4*sd(sex[which(sex$PEDSEX==2),"F"]))
fmax_female
# [1] 0.4833173

# # Y thresholds: 

ymin_male<-mean(sex[which(sex$PEDSEX==1),"YCOUNT"])-(4*sd(sex[which(sex$PEDSEX==1),"YCOUNT"]))
ymin_male
#[1] 102.4786
ymax_female<-mean(sex[which(sex$PEDSEX==2),"YCOUNT"])+(4*sd(sex[which(sex$PEDSEX==2),"YCOUNT"]))
ymax_female
#[1] 60.9409


p1n<-ggplot(sex[which(sex$PEDSEX==1),],aes(x=F,y=YCOUNT)) + geom_point(colour = "#00AFBB") + 
  xlim(min(sex$F),max(sex$F)) + 
  ylim(0,max(sex$YCOUNT)) +  
  geom_hline(yintercept=ymin_male, linetype="dashed", color = "#00AFBB")+
  geom_vline(xintercept=fmin_male, linetype="dashed", color = "#00AFBB")

p2n<-ggplot(sex[which(sex$PEDSEX==2),],aes(x=F,y=YCOUNT)) + geom_point(colour = "#E7B800") + 
  xlim(min(sex$F),max(sex$F)) + 
  ylim(0,max(sex$YCOUNT)) +  
  geom_hline(yintercept=ymax_female, linetype="dashed", color = "#E7B800")+
  geom_vline(xintercept=fmax_female, linetype="dashed", color = "#E7B800")
# + geom_vline(xintercept=fmin_female, linetype="dashed", color = "#E7B800")

p3n<-ggplot(sex[which(sex$PEDSEX==0),],aes(x=F,y=YCOUNT)) + geom_point(colour = "#FC4E07") + 
  xlim(min(sex$F),max(sex$F)) + 
  ylim(0,max(sex$YCOUNT)) +  
  geom_hline(yintercept=ymin_male, linetype="dashed", color = "#00AFBB")+
  geom_vline(xintercept=fmin_male, linetype="dashed", color = "#00AFBB")+ 
  geom_hline(yintercept=ymax_female, linetype="dashed", color = "#E7B800")+
  geom_vline(xintercept=fmax_female, linetype="dashed", color = "#E7B800")
# + geom_vline(xintercept=fmin_female, linetype="dashed", color = "#E7B800")

pn<-ggarrange(p1n,p2n,p3n,ncol=3,labels=c("Male","Female","NA"))
dev.off()

pdf(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_histogram_homozigosity_chrXY_with_newlimits.pdf",sep=""),width =21)
pn
dev.off()


##############################

sex$SNPSEX2<-0
sex$SNPSEX2[which(sex$F<=fmax_female & sex$YCOUNT<=ymax_female)]<-2
sex$SNPSEX2[which(sex$F>=fmin_male & sex$YCOUNT>=ymin_male)]<-1

sex_chrx<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_check_sex_good_chrX.sexcheck",sep=""),head=T)
colnames(sex_chrx)[4]<-"SNPSEX_chrx"
sexm<-merge(sex,sex_chrx[,c("FID","SNPSEX_chrx")],by="FID",sort=F)


table(sexm$PEDSEX,sexm$SNPSEX2)
#     0   1   2
# 1   3 400   1
# 2   0   5 379


samples_exclude<-sex[which((sex$PEDSEX==2 & sex$SNPSEX2==1) | (sex$PEDSEX==1 & sex$SNPSEX2==2)),]
table(samples_exclude$PEDSEX,samples_exclude$SNPSEX2)
#   1 2
# 1 0 1
# 2 5 0

dim(samples_exclude)
#[1] 6 8

write.table(samples_exclude[,c("FID","IID")],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_samples_wrong_gender",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
write.table(samples_exclude[,c("FID","PEDSEX","SNPSEX2","F","YCOUNT")],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_list_samples_sex_discrepancy",sep=""),col.names=T,row.names=F,quote=F,sep="\t")


# recode the rest:

fam<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip.fam",sep=""),head=F)

fam_ed<-merge(fam,sex[,c("FID","SNPSEX2")],by.x="V1",by.y="FID",sort=F)

table(fam_ed$V2==fam$V2)
# TRUE 
# 788  

fam_ed<-fam_ed[,c("V1","V2","V3","V4","SNPSEX2","V6")]

write.table(fam_ed,paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_edited.fam",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

# N cases and ctr excluded:

table(fam$V6[which(fam$V1 %in% samples_exclude$FID)])
# 2 
# 6 

##############

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip.bed \
--bim ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip.bim \
--fam ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_edited.fam \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_samples_wrong_gender \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck
# 570755 variants and 782 people pass filters and QC.
# Among remaining phenotypes, 782 are cases and 0 are controls.

####################################################################################
# 7.4 SET TO MISSING MALE HET CHRX CALLS  - Finish setting hh calls to missing

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck \
--allow-no-sex \
--set-hh-missing \
--make-bed --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh



########################################################################################################################################

###################
# 8.- REMOVE ChrY #
###################

### R

bim<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh.bim",sep=""),head=F)

table(bim$V1)
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 44155 45905 38756 36141 33374 40627 30977 28835 24291 28023 27775 26752 19850 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18222 17091 18556 16664 16409 12941 13644  7983  8203 14417   694   470 

list_toremove<-bim[which(bim$V1==24),"V2",drop=F]

write.table(list_toremove,paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_chry_variants_toremove",sep=""),col.names=F,row.names=F,sep="\t",quote=F)

######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh \
--allow-no-sex \
--exclude ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_chry_variants_toremove \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry
# 570061 variants and 782 people pass filters and QC.



#########################################################################################################################################

##########
# 9 - QC #
##########

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry \
--allow-no-sex \
--missing --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry

### R

sample_miss<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.imiss",sep=""),head=T)
var_miss<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.lmiss",sep=""),head=T)

summary(sample_miss$F_MISS)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001435 0.003144 0.004202 0.004480 0.005496 0.015580 

summary(var_miss$F_MISS)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00000 0.00000 0.00000 0.00448 0.00000 0.99740 

####################################
# 9.1 - REMOVE SAMPLES CallPP <0.80%

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry \
--allow-no-sex \
--mind 0.20 \
--make-bed --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8
# 0 people removed due to missing genotype data (--mind).
# 570061 variants and 782 people pass filters and QC.


####################################
# 9.2- REMOVE VARIANTS CallPP <0.80%

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8 \
--geno 0.20 \
--make-bed --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8
# 3542 variants removed due to missing genotype data (--geno).
# 566519 variants and 782 people pass filters and QC.


####################################
# 9.3 - REMOVE SAMPLES CallPP <0.95%

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8 \
--allow-no-sex \
--mind 0.05 \
--make-bed --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95
# 0 people removed due to missing genotype data (--mind).
# 566519 variants and 782 people pass filters and QC.


#####################################
# 9.4 - REMOVE VARIANTS CallPP <0.95%

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95 \
--allow-no-sex \
--geno 0.05 \
--make-bed --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95
# 8147 variants removed due to missing genotype data (--geno).
# 558372 variants and 782 people pass filters and QC.



#######################################################
# 9.5 - REMOVE VARIANTS FREQ <0.01 AND CallPP <0.98%

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95 \
--allow-no-sex \
--missing --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95


/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95 \
--allow-no-sex \
--freq --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95

###### R

sample_miss<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95.imiss",sep=""),head=T)
var_miss<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95.lmiss",sep=""),head=T)
frq<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95.frq",sep=""),head=T)

var<-merge(frq[,c(2:6)],var_miss,by="SNP")
var.1<-var[which(var$MAF<0.01 & var$F_MISS>0.02),]
dim(var.1)
#[1] 1372     9

write.table(var.1[,"SNP",drop=F],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_monomorphic_vcr098_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


##############

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95 \
--allow-no-sex \
--exclude ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/list_monomorphic_vcr098_var_exclude \
--make-bed \
--out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01
# 557000 variants and 782 people pass filters and QC.



#########################################################################################################################################

##############################################################################
# 10.- REMOVE VARIANTS WITH SIGNIF DISCREPANCY IN CALL RATE BETWEN CA VS CTR #
##############################################################################

###### R

fam<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01.fam",sep=""),head=F)

table(fam$V6,useNA="ifany")
# 2 
# 782 

write.table(fam[which(fam$V6==1),1:2],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_controls",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
write.table(fam[which(fam$V6==2),1:2],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_cases",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

# ##############
# 
# /path/to/software/plink_linux_x86_64_20181202/./plink \
# --bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01 \
# --allow-no-sex \
# --keep ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_controls \
# --missing \
# --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_ctr
# 
# /path/to/software/plink_linux_x86_64_20181202/./plink \
# --bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01 \
# --allow-no-sex \
# --keep ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_cases \
# --missing \
# --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_cases
# 
# ###### R
# 
# ca_var_miss<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_cases.lmiss",sep=""),head=T)
# ctr_var_miss<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_ctr.lmiss",sep=""),head=T)
# 
# ca_var_miss$N_NO_MISS<-ca_var_miss$N_GENO-ca_var_miss$N_MISS
# ctr_var_miss$N_NO_MISS<-ctr_var_miss$N_GENO-ctr_var_miss$N_MISS
# 
# colnames(ca_var_miss)[c(3,6)]<-paste(colnames(ca_var_miss)[c(3,6)],"_cases",sep="")
# colnames(ctr_var_miss)[c(3,6)]<-paste(colnames(ctr_var_miss)[c(3,6)],"_ctr",sep="")
# 
# all<-merge(ca_var_miss[,c(2,3,6)],ctr_var_miss[,c(2,3,6)],by="SNP",sort=F)
# 
# 
# foo <- function(y){
#   # include here as.numeric to be sure that your values are numeric:
#   table <-  matrix(as.numeric(c(y[2], y[3], y[4], y[5])), ncol = 2, byrow = TRUE)
#   if(any(is.na(table))) p <- "error" else p <- fisher.test(table, alternative="two.sided")$p.value
#   p
# } 
# all$fisher_pvalue <- apply(all, 1, foo)
# 
# summary(all$fisher_pvalue)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #  0.00    1.00    1.00    0.87    1.00    1.00 
# 
# 
# dim(all[which(all$fisher_pvalue<1E-4),])
# #[1] 51  6
# 
# 
# pdf(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/hist_missingness_pvalue_ca_ctr.pdf",sep=""),height = 5,width = 10)
# ggplot(all, aes(x=(-log10(fisher_pvalue)))) + geom_histogram(binwidth=1) + geom_vline(xintercept=(-log10(1E-4)), linetype="dashed",color = "red", size=1)
# dev.off()
# 
# write.table(all,paste(path,"pre_imputation/QC/rioux_igenomed_gsa/comparison_missingness_cases_ctr.lmiss",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
# write.table(all[which(all$fisher_pvalue<1E-4),1,drop=F],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_variants_exclude_by_missingness_cases_ctr.lmiss",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

##############

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01 \
--allow-no-sex \
--make-bed \
--out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep
# 557000 variants and 782 people pass filters and QC.
# Among remaining phenotypes, 782 are cases and 0 are controls.

# --exclude ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_variants_exclude_by_missingness_cases_ctr.lmiss \

############################################################################################################################################

#############################################################################
# 11.- REMOVE VARIANTS WITH SIGNIF DISCREPANCY IN CALL RATE BETWEEN STUDIES #
#############################################################################

# This cohort is composed by just one study

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep \
--allow-no-sex \
--missing \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy
# 557000 variants and 782 people pass filters and QC.
# Among remaining phenotypes, 782 are cases and 0 are controls.



############################################################################################################################################

################################################
# 12.1- REMOVE INTRA-COHORT DUPLICATED SAMPLES #
################################################

# The amount of computer memory required by KING analysis is modest, at ~N ✕ M / 4 (where N is the number of samples and M is the number of SNPs)
# plus a small percentage of overhead cost. E.g., for a dataset consisting of 100,000 samples each genotyped at 1,000,000 SNPs, 
# the required memory size is ~25GB.


# (557000*782/4)/1E+6
# [1] 108.8935

path_gwas=/path/to/ibdgwas/IIBDGC/
  MEM=675

# recommended method for biobank-level estimations, faster than --kinship
bsub -J"kg" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/logs/stderr_king_rioux_igenomed_gsa \
-o ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/logs/stdout_king_rioux_igenomed_gsa \
"/path/to/software/./king \
-b ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy.bed \
--related --cpus 1 --prefix ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_king"
# Job <250858> is submitted to queue <normal>.


# The first algorithm estimates pair-wise kinship coefficients (through option --kinship), and the second algorithm infers pairwise IBD 
# (identical by descent) segments (through option --ibdseg). Both algorithms can also be integrated in a single inference procedure through 
# option --related.

### R

kin<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_king.kin0",sep=""),head=T)
table(kin$InfType)
# Dup/MZ     FS     PO 
# 1      5      4 

# Duplicated or Monozygotic twin  Dup/MZ
# Parent–offspring                PO
# Full sib                        FS
# 2nd Degree                      2nd

table(cut(kin$Kinship,breaks=c(1,0.354,0.177)))
# (0.177,0.354]     (0.354,1] 
#             9            1 

# however inferred reationship does not match kinship coefficient, two PO and one FS with Kinship <0.177, because it comes from ibdseg method,
# keep using kinship:

#kinship coefficient range
# >0.354 duplicate/MZ twin 
# [0.177, 0.354] 1st-degree 
# [0.0884, 0.177] 2nd-degree

# exclude then one that does not match with IBS, see below
# dup<-kin[which(kin$Kinship>0.354),c("FID1","FID2","Kinship")]

dup<-kin[which(kin$Kinship>0.4),c("FID1","FID2","Kinship")]

dup_ids<-c(as.character(dup$FID1),as.character(dup$FID2))
length(dup_ids)
# [1] 2
table(length(dup_ids)==(nrow(dup)*2))
# TRUE 
# 1 

dup_ids<-dup_ids[!duplicated(dup_ids)]
length(dup_ids)
#[1] 2

# CHECK ALSO THAT BOTH SAMPLES HAVE SAME PHENOTYPE AND GENDER, AND REMOVE SAMPLES LOWER CALL RATES

fam<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy.fam",sep=""),head=F)
colnames(fam)[5:6]<-c("sex","pheno")
sample_miss<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy.imiss",sep=""),head=T)

all<-merge(fam[,c(1,5,6)],sample_miss[,c("FID","F_MISS")],by.x="V1",by.y="FID",all.x=T,sort=F)

tmp<-all[which(all$V1 %in% dup_ids),]
dup_ids<-dup_ids[match(tmp$V1,dup_ids)]
rm(tmp)


for (i in 1:length(dup_ids)) {
  
  tmp1<-dup[which(dup$FID1==dup_ids[i]),]
  tmp2<-dup[which(dup$FID2==dup_ids[i]),]
  colnames(tmp2)[1:2]<-colnames(tmp2)[2:1]
  
  tmp<-rbind(tmp1,tmp2)
  
  ids_tmp<-c(as.character(tmp$FID1),as.character(tmp$FID2))
  ids_tmp<-ids_tmp[!duplicated(ids_tmp)]
  
  #keep same order as in dup_ids and in all
  ids_tmp<-ids_tmp[match(dup_ids[which(dup_ids %in% ids_tmp)],ids_tmp)]
  
  # get number of possible combinations
  n_possible_combinations<-nrow(permutations(length(ids_tmp), 2, letters[1:length(ids_tmp)]))/2
  
  
  if ( nrow(dup[which((dup$FID1 %in% ids_tmp) | (dup$FID2 %in% ids_tmp)),])==n_possible_combinations ) {
    
    # for duplicated samples that have same phenotype and sex, remove the ones with smaller call rate
    
    data<-as.data.frame(matrix(ncol=1,nrow=length(ids_tmp)))
    data$V1<-ids_tmp
    data<-merge(data,all,by="V1",all.x=T,sort=F)
    
    if ( (dim(table(data$sex))==1 & dim(table(data$pheno))==1) ) {
      
      if(!exists("data_remove")) {
        keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]
        data_remove<-data[which(!data$V1 %in% keep_sample),]
      } else {
        keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]
        data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
      }
      
    } else {
      
      data<-as.data.frame(matrix(ncol=1,nrow=length(ids_tmp)))
      data$V1<-ids_tmp
      data<-merge(data,all,by="V1",all.x=T,sort=F)
      
      # for duplicated samples that do not have same pheno and sex, remove all
      
      # print("Samples with different sex/pheno:")
      # print(data)
      
      if(!exists("data_remove")) {
        
        data_remove<-data
        
      } else {
        
        data_remove<-rbind(data_remove,data)
        
      }
      
      if(!exists("data_inconsist")){
        
        jj<-1
        data_inconsist<-data
        data_inconsist$group<-jj
        
      } else {
        
        jj<-jj+1
        data$group<-jj
        data_inconsist<-rbind(data_inconsist,data)
        
      }
    }
    
  } else {
    
    # not all combinations of duplicated pairs found, show list of IDs, to manually inspect issues
    
    print(paste("Number of expected combinations: ",n_possible_combinations,sep=""))
    print(paste("Number of observed combinations: ",nrow(dup[which((dup$FID1 %in% ids_tmp) | (dup$FID2 %in% ids_tmp)),]),sep=""))
    print(ids_tmp)
    
  }
  
}


dim(data_remove)
# [1] 2  4

data_remove<-data_remove[!duplicated(data_remove$V1),]
dim(data_remove)
# [1] 1   4

# FINAL DOUBLE CHECK
# all[which(all$V1 %in% c(as.character(dup$FID1),as.character(dup$FID2))),]
#

data_remove<-data_remove[,c(1,1)]
colnames(data_remove)<-c("FID","IID")

write.table(data_remove,paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_duplicated_samples",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

#### create a file to report inconsitend data:

# dim(data_inconsist)
# # [1] 
# data_inconsist<-data_inconsist[!duplicated(data_inconsist$V1),]
# dim(data_inconsist)
# # 
# # 
# data_inconsist$cohort<-"rioux_igenomed_gsa"
# dup[which((dup$FID1 %in% data_inconsist$V1) & (dup$FID2 %in% data_inconsist$V1)),]
# # FID1       FID2 Kinship
# # 30 sample_id sample_id  0.3705
# do not keep as dup!

# 
# # cohorts<-c("")
# # 
# # for (i in 1:length(cohorts)){
# #   fam<-read.table(paste(path,"pre_imputation/QC/XXX/",cohorts[i],"_hg19.fam",sep=""),head=F)
# #   data_inconsist$cohort[which(data_inconsist$V1 %in% fam$V1)]<-cohorts[i]
# # }
# 
# write.table(data_inconsist,paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_list_duplicated_samples_inconsistent_phenotype",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

######

path_gwas=/path/to/ibdgwas/IIBDGC/
  MEM=500

bsub -J"kg" -M"$MEM" -n4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/logs/stderr_plink_rioux_igenomed_gsa_12 \
-o ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/logs/stdout_plink_rioux_igenomed_gsa_12 \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_duplicated_samples \
--threads 4 \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample"
# Job <908954> is submitted to queue <normal>.

less ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/logs/stdout_plink_rioux_igenomed_gsa_12
# 557000 variants and 781 people pass filters and QC.
# Among remaining phenotypes, 781 are cases and 0 are controls.


############################################################################################################################################

################################################
# 12.2 - REMOVE INTRA-STUDY DUPLICATED SAMPLES #
################################################

# no other studies from same supplier

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy

# 557000 variants and 781 people pass filters and QC.
# Among remaining phenotypes, 781 are cases and 0 are controls.



############################################################################################################################################

#################################################
# 12.3 - REMOVE INTER-COHORT DUPLICATED SAMPLES #
#################################################

# see script evaluate_relatedness_IIBDGC.R

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_iibdgc_data_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy

# 557000 variants and 769 people pass filters and QC.
# Among remaining phenotypes, 769 are cases and 0 are controls.



############################################################################################################################################

########################################
# 13 - PCA TO DEFINE BROAD POPULATIONS #
########################################

# see script pca_iibdgc_1000GP_iibdgc_projection.R


############################################################################################################################################


################################
# 14 - HWE IN EUR ONLY SAMPLES #
################################

# samples genotyped in one batch
# evaluate HWE per batch and within EUR non-Jewish ancestry samples

### R

cohorts=("rioux_igenomed_gsa")

for (i in 1:length(cohorts)) {
  tmp<-read.table(paste(path,
                        "pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy.fam",sep=""),head=F)
  colnames(tmp)<-c("FID","IID")
  write.table(tmp,paste(path,"pre_imputation/QC/",cohorts[i],"/list_ids_",cohorts[i],sep=""),col.names=F,row.names=F,quote=F,sep="\t")
}



######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy \
--allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_nonjewish_ancestry_samples_noDuplicates \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/hwe_eur_nonjewish_autosomal
# 557000 variants and 464 people pass filters and QC.


cohorts=(rioux_igenomed_gsa)

for i in ${cohorts[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/hwe_eur_nonjewish_autosomal \
--allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_ids_$i \
--hardy \
--out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/hwe_eur_nonjewish_autosomal_$i
done


###

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy \
--allow-no-sex \
--filter-females --chr 23 25 \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/hwe_chr23_females

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/hwe_chr23_females \
--allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_nonjewish_ancestry_samples_noDuplicates --chr 23 25 \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/hwe_eur_nonjewish_chr23_females

cohorts=(rioux_igenomed_gsa)

for i in ${cohorts[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/hwe_eur_nonjewish_chr23_females \
--allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_ids_$i \
--hardy \
--out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/hwe_eur_nonjewish_chr23_females_$i
done


#### R

cohorts<-c("rioux_igenomed_gsa")

for (i in 1:length(cohorts)){
  if(i==1){
    aut_hwe<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/hwe_eur_nonjewish_autosomal_",cohorts[i],".hwe",
                              sep=""),head=T)
    aut_hwe$cohort<-cohorts[i]
  }
  if (i!=1) {
    tmp<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/hwe_eur_nonjewish_autosomal_",cohorts[i],".hwe",
                          sep=""),head=T)
    tmp$cohort<-cohorts[i]
    aut_hwe<-rbind(aut_hwe,tmp)
  }
}


aut_hwe<-aut_hwe[which(!aut_hwe$CHR %in% c(23,25)),]

dim(aut_hwe[which(aut_hwe$TEST=="AFF" & aut_hwe$P<=1E-12),])

remove_1<-aut_hwe[which(aut_hwe$TEST=="UNAFF" & aut_hwe$P<=1E-5),]
remove_2<-aut_hwe[which(aut_hwe$TEST=="AFF" & aut_hwe$P<=1E-12),]
remove_2<-remove_2[which(!remove_2$SNP %in% remove_1$SNP),]



for (i in 1:length(cohorts)){
  if(i==1){
    chr23_hwe<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/hwe_eur_nonjewish_chr23_females_",cohorts[i],".hwe",
                                sep=""),head=T)
    chr23_hwe$cohort<-cohorts[i]
  }
  if (i!=1) {
    tmp<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/hwe_eur_nonjewish_chr23_females_",cohorts[i],".hwe",
                          sep=""),head=T)
    tmp$cohort<-cohorts[i]
    chr23_hwe<-rbind(chr23_hwe,tmp)
  }
}

remove_3<-chr23_hwe[which(chr23_hwe$TEST=="UNAFF" & chr23_hwe$P<=1E-5),]
remove_4<-chr23_hwe[which(chr23_hwe$TEST=="AFF" & chr23_hwe$P<=1E-12),]
remove_4<-remove_4[which(!remove_4$SNP %in% remove_3$SNP),]

dim(remove_1)
dim(remove_2)
dim(remove_3)
dim(remove_4)

remove<-rbind(remove_1,remove_2,remove_3,remove_4)
remove<-remove[,c(1:5)]
remove<-remove[!duplicated(remove),]
dim(remove)
# [1] 296    5

table(dim(remove)==dim(remove[!duplicated(remove$SNP),]))

table(remove$TEST)
# 

table(remove$CHR)
#

write.table(remove[,"SNP",drop=F],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_var_exclude_hwe",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


########

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy \
--allow-no-sex \
--exclude ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/list_var_exclude_hwe \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_hwe

# 556704 variants and 769 people pass filters and QC.
# Among remaining phenotypes, 769 are cases and 0 are controls.



/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy \
--allow-no-sex \
--exclude ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/list_var_exclude_hwe \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe

# 556704 variants and 782 people pass filters and QC.
# Among remaining phenotypes, 782 are cases and 0 are controls.



############################################################################################################################################

#######################
# 15 - HETEROZYGOSITY #
#######################

# see script pca_iibdgc_1000GP_iibdgc_projection.R

##############################################################################
# 15.1 Use list non duplicated samples to estimate heterozigosity limits:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_hwe \
--allow-no-sex \
--chr 1-22 --het \
--out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/autosomal_het_nodup


#### R

pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/pca_1000gp_iibdgc_pca_pcaEur_withDuplicates.tsv",sep=""),head=T,sep="\t")

# set thresholds from non duplicated set:

het<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/autosomal_het_nodup.het",sep=""),head=T)
dim(het)
#

het$het<-(het$N.NM.-het$O.HOM)/het$N.NM.

het<-merge(het[,c("FID","O.HOM.","E.HOM.","N.NM.","F","het")],pca[,c("FID","cohort","inferred_population","self_jewish","pca_jewish","PC1_EUR")],
           by="FID",all.x=T)

table(het$inferred_population)
# AFR AMR EAS EUR SAS 
# 1  21   0 747   0 

het$inferred_population_2<-het$inferred_population
het$inferred_population_2<-as.character(het$inferred_population_2)
het$inferred_population_2[which(het$pca_jewish=="Jewish")]<-"EUR_Jewish"
het$inferred_population_2[which(het$pca_jewish=="Non-Jewish")]<-"EUR_Non-Jewish"

table(het$inferred_population_2)
# AFR            AMR     EUR_Jewish EUR_Non-Jewish 
# 1             21            283            464 


het$inferred_population_2<-as.factor(het$inferred_population_2)

#### CREATE AN HISTOGRAM

# number of bins in histogram
fd=function(x) {
  n=length(x)
  r=IQR(x)
  2*r/n^(1/3)
}

limits<-matrix(ncol=3,nrow=length(levels(het$inferred_population_2)) )
limits<-as.data.frame(limits)
colnames(limits)<-c("pop","lower","upper")

for (i in c(1:length(levels(het$inferred_population_2))) ) {
  
  limits$pop[i]<-levels(het$inferred_population_2)[i]
  
  limits$lower[i]<-mean(het$het[which(het$inferred_population_2==levels(het$inferred_population_2)[i])])-(4*sd(het$het[which(het$inferred_population_2==levels(het$inferred_population_2)[i])]))
  limits$upper[i]<-mean(het$het[which(het$inferred_population_2==levels(het$inferred_population_2)[i])])+(4*sd(het$het[which(het$inferred_population_2==levels(het$inferred_population_2)[i])]))
  
}

##############################################################################
# 15.1 Apply limits to all samples, including duplicated:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe \
--exclude ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/list_var_exclude_hwe \
--allow-no-sex \
--chr 1-22 --het \
--out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/autosomal_het_all
# 542625 variants and 782 people pass filters and QC.
# Among remaining phenotypes, 782 are cases and 0 are controls.


####

pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/pca_1000gp_iibdgc_pca_pcaEur_withDuplicates.tsv",sep=""),head=T,sep="\t")

# set thresholds from non duplicated set:

het<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/autosomal_het_all.het",sep=""),head=T)
dim(het)
# 

het$het<-(het$N.NM.-het$O.HOM)/het$N.NM.

het<-merge(het[,c("FID","O.HOM.","E.HOM.","N.NM.","F","het")],pca[,c("FID","cohort","inferred_population","self_jewish","pca_jewish","PC1_EUR")],
           by="FID",all.x=T)

table(het$inferred_population,useNA="ifany")
# AFR AMR EAS EUR SAS 
# 1  21   0 760   0

het$inferred_population_2<-het$inferred_population
het$inferred_population_2<-as.character(het$inferred_population_2)
het$inferred_population_2[which(het$pca_jewish=="Jewish")]<-"EUR_Jewish"
het$inferred_population_2[which(het$pca_jewish=="Non-Jewish")]<-"EUR_Non-Jewish"

table(het$inferred_population_2,useNA="ifany")
# AFR            AMR     EUR_Jewish EUR_Non-Jewish 
# 1             21            291            469 

het$inferred_population_2<-as.factor(het$inferred_population_2)

#### CREATE AN HISTOGRAM

# number of bins in histogram
fd=function(x) {
  n=length(x)
  r=IQR(x)
  2*r/n^(1/3)
}


# to emulate ggplot2 palette:
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = length(levels(het$inferred_population_2))
cols = gg_color_hue(n)

xmax<-max(het$het)+0.01
xmin<-min(het$het)

for (i in c(1:length(levels(het$inferred_population_2))) ) {
  print(i)
  fd_bin<-fd(het$het[which(het$inferred_population_2==levels(het$inferred_population_2)[i])])
  assign(paste("p",i,sep=""),ggplot(het[which(het$inferred_population_2==levels(het$inferred_population_2)[i]),],aes(x=het,color=inferred_population_2,fill=inferred_population_2)) +
           geom_histogram(position="identity", alpha=0.5, binwidth=fd_bin) +
           scale_fill_manual(values=cols[i]) + scale_color_manual(values=cols[i]) + 
           geom_vline(xintercept = limits$lower[i],linetype="dotted", color = "black") +  
           geom_vline(xintercept = limits$upper[i],linetype="dotted", color = "black") + ggtitle(levels(het$inferred_population_2)[i]) +
           xlim(xmin,xmax) )
}

pdf(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_het_histogram.pdf",sep=""),width=18,height=(length(levels(het$inferred_population_2)) ) )
ggarrange(p1,p2,p3,p4,nrow=length(levels(het$inferred_population_2))/2,ncol=2 )
dev.off()
system(paste("cp ",path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_het_histogram.pdf ~/tmp_plots/",sep=""))


for (i in c(1:nrow(limits)) ) {
  if(i==1) {
    remove<-het[which( (het$inferred_population_2==limits$pop[i]) &
                         ( (het$het<limits$lower[i]) |
                             (het$het>limits$upper[i]) ) ), ]
  } else {
    tmp<-het[which( (het$inferred_population_2==limits$pop[i]) &
                      ( (het$het<limits$lower[i]) |
                          (het$het>limits$upper[i]) ) ), ]
    remove<-rbind(remove,tmp)
  }
  
}

table(remove$inferred_population_2)
# AFR            AMR     EUR_Jewish EUR_Non-Jewish 
# 0              0              4              2 

dim(remove)
# [1] 6  12

write.table(remove[,c("FID","FID")],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_samples_het_outliers",sep=""),col.names=T,row.names=F,quote=F,sep="\t")



# no duplicated:

fam<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_hwe.fam",sep="")
                ,head=F)

dim(remove[which(remove$FID %in% fam$V1),])
# [1] 5 12


######################

### SAMPLES WITHOUT DUPLICATES

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_hwe \
--allow-no-sex \
--remove ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/list_samples_het_outliers \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_hwe_het
# 556704 variants and 764 people pass filters and QC.
# Among remaining phenotypes, 764 are cases and 0 are controls.



### SAMPLES WITH DUPLICATES

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe \
--exclude ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/list_var_exclude_hwe \
--remove ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/list_samples_het_outliers \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het
# 556704 variants and 776 people pass filters and QC.
# Among remaining phenotypes, 776 are cases and 0 are controls.




############################################################################################################################################

######################
# 16.- UPDATE TO B38 #
######################

# NOTE: data in build 37, needs to be liftover first:

### SAMPLES WITH DUPLICATES

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het \
--recode tab --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het

# ${path_gwas}previous_qced_b38/liftover/liftOverx
# usage:
# liftOver oldFile map.chain newFile unMapped

#### R convert map into bed:

map<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.map",sep=""),head=F)

map$chr<-paste("chr",map$V1,sep="")
map$chr[which(map$V1=="23")]<-"chrX"
map$chr[which(map$V1=="25")]<-"chrX"
map$chr[which(map$V1=="24")]<-"chrY"
map$chromStart<-map$V4
map$chromEnd<-map$V4+1

map$chromStart<-format(map$chromStart, scientific=F)
map$chromEnd<-format(map$chromEnd, scientific=F)

dim(map)
# [1] 556704      7

write.table(map[,c(5:7,2)],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_postqc.bed",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#### lift positions

${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_postqc.bed \
${path_gwas}previous_qced_b38/liftover/hg19ToHg38.over.chain.gz \
${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_postqc_lifted_hg38 \
${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_postqc_no_lifted_hg38

cut -f 4 ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_postqc_no_lifted_hg38 | sed "/^#/d" \
> ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/nonlifted_hg38_variants_to_exclude_tmp.dat

grep "alt" ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_postqc_lifted_hg38 | cut -f 4 \
| cat - ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/nonlifted_hg38_variants_to_exclude_tmp.dat > \
${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/nonlifted_hg38_variants_to_exclude.dat

wc -l ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/nonlifted_hg38_variants_to_exclude.dat
# 132

### exclude non lifted:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het \
--exclude ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/nonlifted_hg38_variants_to_exclude.dat \
--recode tab --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_postqc_lifted_hg38_liftedvariants
# 556572 variants and 776 people pass filters and QC.
# Among remaining phenotypes, 776 are cases and 0 are controls.



## R

bed_lifted<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_postqc_lifted_hg38",sep=""),head=F)
excluded<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/nonlifted_hg38_variants_to_exclude.dat",sep=""),head=F)
bed_lifted<-bed_lifted[which(!bed_lifted$V4 %in% excluded$V1),]
map<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_postqc_lifted_hg38_liftedvariants.map",sep=""),head=F)

table(as.character(map$V2)==as.character(bed_lifted$V4))
# TRUE 
# 556572 

map$pos<-bed_lifted$V2
table(as.character(map$V2)==as.character(bed_lifted$V4))
table(map$pos==bed_lifted$V2)
# TRUE 
# 556572 

write.table(map[,c(1:3,5)],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_postqc_lifted_hg38_liftedvariants_edited.map",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


########

# put together updated .ped plus lifted.map
/path/to/software/plink_linux_x86_64_20181202/./plink \
--ped ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_postqc_lifted_hg38_liftedvariants.ped \
--map ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_postqc_lifted_hg38_liftedvariants_edited.map \
--merge-x 'no-fail' \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38
# 556572 variants and 776 people pass filters and QC.
# Among remaining phenotypes, 776 are cases and 0 are controls.


######################################################
# 16.1 - REMOVE MONOMORPHIC VARIANTS

### R

bim<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38.bim",sep=""),head=F)

table(bim$V5,bim$V6)
#

table(bim[which(bim$V5=="0"),"V1"])
#

monom<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]
dim(monom)
#[1] 357

write.table(monom[,"V2",drop=F],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_indel_var_exclude_2",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

######################################################
# 16.2 - CREATE FILE FORCING A1 ALLELE TO BE REFERENCE

# force A1 and A2 to be ref and alt alleles:

zcat ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_posstrandaligned.vcf.gz | cut -f '1-5' | awk '{print $1":"$2"_"$4"_"$5,$4}' | \
sed '/^##/d' > ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_variants_rioux_igenomed_gsa_hg19_posstrandaligned_with_A1
wc -l ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_variants_rioux_igenomed_gsa_hg19_posstrandaligned_with_A1
# 571174 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/list_variants_gwas3_hg19_posstrandaligned_with_A1

sort ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_variants_rioux_igenomed_gsa_hg19_posstrandaligned_with_A1 \
| uniq > ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_variants_rioux_igenomed_gsa_hg19_posstrandaligned_with_A1_ed

wc -l ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_variants_rioux_igenomed_gsa_hg19_posstrandaligned_with_A1_ed
# 570880 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/list_variants_rioux_igenomed_gsa_hg19_posstrandaligned_with_A1_ed


/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38 \
--exclude ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/list_indel_var_exclude_2 \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom
# 556215 variants and 776 people pass filters and QC.
# Among remaining phenotypes, 776 are cases and 0 are controls.

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom \
--allow-no-sex \
--a2-allele ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_variants_rioux_igenomed_gsa_hg19_posstrandaligned_with_A1_ed \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt
# --a2-allele: 556215 assignments made.
# 556215 variants and 776 people pass filters and QC.
# Among remaining phenotypes, 776 are cases and 0 are controls.



##############################################################
# 16.1 DOUBLE CHECK REF ALT ALLELE WITH HG38 REF SEQUENCE:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt \
--allow-no-sex \
--keep-allele-order \
--output-chr M --recode vcf-iid --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt

# Note that most PLINK analyses treat the A1 (usually minor) allele as the reference allele, which makes sense when only biallelic variants are involved.
# However, since it is conventional for VCF files to set the major allele as the reference allele instead

# double check alleles, variants:
export BCFTOOLS_PLUGINS=/path/to/software/bcftools-1.9/plugins 

/path/to/software/bcftools-1.9/./bcftools +fixref ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt.vcf \
-Oz -o ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned.vcf.gz \
-- -f ${path_gwas}resources/hg38/hg38_edited.fa -m top

# # SC, guessed strand convention
# SC	TOP-compatible	0
# SC	BOT-compatible	0
# # ST, substitution types
# ST	A>C	23799	4.3%
# ST	A>G	96889	17.4%
# ST	A>T	733	0.1%
# ST	C>A	28020	5.0%
# ST	C>G	1182	0.2%
# ST	C>T	127732	23.0%
# ST	G>A	127573	22.9%
# ST	G>C	1153	0.2%
# ST	G>T	27753	5.0%
# ST	T>A	773	0.1%
# ST	T>C	96967	17.4%
# ST	T>G	23641	4.3%
# # NS, Number of sites:
# NS	total        	556215
# NS	ref match    	554577	99.7%
# NS	ref mismatch 	1638	0.3%
# NS	flipped      	139	0.0%
# NS	swapped      	1349	0.2%
# NS	flip+swap    	2092	0.4%
# NS	unresolved   	5	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0

# REMOVE MISSMATCH ALLELES, THEY  IMPUTATION SERVER TO CRASH:
/path/to/software/bcftools-1.9/./bcftools norm \
--check-ref x \
-f /${path_gwas}resources/hg38/hg38_edited.fa \
${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned.vcf.gz \
-Oz -o ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_2.vcf.gz
# Lines   total/split/realigned/skipped:	556215/0/0/5

# Double check:
/path/to/software/bcftools-1.9/./bcftools +fixref ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_2.vcf.gz \
-- -f ${path_gwas}resources/hg38/hg38_edited.fa

# # SC, guessed strand convention
# SC	TOP-compatible	0
# SC	BOT-compatible	0
# # ST, substitution types
# ST	A>C	23808	4.3%
# ST	A>G	96919	17.4%
# ST	A>T	736	0.1%
# ST	C>A	28021	5.0%
# ST	C>G	1188	0.2%
# ST	C>T	127681	23.0%
# ST	G>A	127555	22.9%
# ST	G>C	1145	0.2%
# ST	G>T	27733	5.0%
# ST	T>A	767	0.1%
# ST	T>C	97006	17.4%
# ST	T>G	23651	4.3%
# # NS, Number of sites:
# NS	total        	556210
# NS	ref match    	556210	100.0%
# NS	ref mismatch 	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0


#### VCF to BED

/path/to/software/plink_linux_x86_64_20181202/./plink \
--vcf ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_2.vcf.gz \
--keep-allele-order --allow-no-sex \
--double-id \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned
# 556210 variants and 776 people pass filters and QC.
# Note: No phenotypes present.

# zcat ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_2.vcf.gz | head -40 | cut -f 1-5
# head -10 ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned.bim


##############################################################
# 16.2 UPDATE NAME VARIANTS TO CHR:POSITION_REF_ALT in b38

zcat ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_2.vcf.gz | cut -f '1-5' | awk '{print $3,$1":"$2"_"$4"_"$5}' \
> ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_variants_rioux_igenomed_gsa_hg38_posstrandaligned


### R

bim<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned.bim",sep=""),sep="\t",head=F)
ids<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_variants_rioux_igenomed_gsa_hg38_posstrandaligned",sep=""),sep=" ",head=F,skip=34)

colnames(ids)[2]<-"ids"

table(bim$V2==ids$V1)
# TRUE 
# 556210 

#the major allele is set to A2 by default by Plink, keep ids with real ref/alt as in vcf using the ids file
bim.1<-cbind(bim,ids[,"ids",drop=F])

bim.1$ids<-gsub("X:","23:",bim.1$ids)

# two variants from GSA liftover to same position, b37 10:17808880_T_C and 10:18055858_C_T, remove second (= rs1556464660):
bim.1$ids<-as.character(bim.1$ids)
bim.1[which(bim.1$ids=="10:17766881_T_C"),]
bim.1$ids[which(bim.1$V2=="10:18055858_C_T")]<-"10:17766881_T_C_dup"
write.table(bim.1$ids[which(bim.1$V2=="10:18055858_C_T")],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_variants_exclude_dup2",sep="")
            ,col.names=F,row.names=F,quote=F,sep="\t")


write.table(bim.1[,c(1,7,3:6)],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_edited.bim",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


##############################################################
# 16.3 RENAME VARIANTS USING B38

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned.bed \
--bim ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_edited.bim \
--fam ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt.fam \
--exclude ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_variants_exclude_dup2 \
--keep-allele-order --allow-no-sex \
--freq counts \
--make-bed --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated
# 556209 variants and 776 people pass filters and QC.
# Among remaining phenotypes, 776 are cases and 0 are controls.



###################################################################################################################################################################################

##################################
# 17.- KEEP ONLY TOPMed VARIANTS #
##################################

path_gwas=/path/to/ibdgwas/IIBDGC/
  MEM=200

bsub -J"TOPMed" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/logs/stderr_TOPMed_subset \
-o ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/logs/stdout_TOPMed_subset \
"awk 'NR==FNR{vals[\$2];next} (\$9) in vals' ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated.bim \
${path_gwas}resources/TOPMed/PASS.Variantsbravo-dbsnp-all_edited_3.vcf > ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_TOPMed_variants"
# Job <695548> is submitted to queue <normal>.


cat ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_TOPMed_variants | cut -f 9 > ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_TOPMed_variants_ids

wc -l ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_TOPMed_variants_ids
# 538573 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_TOPMed_variants


### R

bim<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated.bim",sep=""),head=F)
tm<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_TOPMed_variants",sep=""),head=F)
frq<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated.frq.counts",sep=""),head=T)


frq$MAF<-frq$C1/(frq$C1+frq$C2)
frq$MAF[which(frq$MAF>0.5)]<-1-frq$MAF[which(frq$MAF>0.5)]


dim(tm)
#

dim(bim[which(!bim$V2 %in% tm$V9),])
# [1] 17636     6

table(bim[which(!bim$V2 %in% tm$V9),"V1"])
# 

# # save this list for Kyle:
# write.table(notm,paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsalist_variants_not_TOPMed",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
# system(paste("cp ",path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsalist_variants_not_TOPMed ~/",sep=""))

#######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated \
--keep-allele-order --allow-no-sex \
--extract ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_TOPMed_variants_ids \
--make-bed --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed
# 538573 variants and 776 people pass filters and QC.
# Among remaining phenotypes, 776 are cases and 0 are controls.

###################################################################################################################################################################################


####################################
# 18.- DOUBLE CHEK A/T C/G ALLELES #
####################################


awk -v OFS='\t' '{print $2,$6}' ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim > \
${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_variants_rioux_igenomed_gsa_hg38_posstrandaligned_with_A1


# NOTE: this cohort does not include ctr!

# # EDIT CONTROL LIST, SO FID AND IID ARE THE SAME
# awk -v OFS='\t' '{print $2,$2}' ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_controls > ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_controls2


### EUR NON-JEWISH CTR:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed \
--keep-allele-order --allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_nonjewish_ancestry_samples_noDuplicates \
--make-bed --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_eur


/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_eur \
--keep-allele-order --allow-no-sex \
# --keep ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_controls2 \
--freq counts \
--out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_eur_ctr


# GNOMAD:

path_gwas=/path/to/ibdgwas/IIBDGC/
  MEM=200

bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/logs/stderr_gnomad_subset \
-o ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/logs/stdout_gnomad_subset \
"awk 'NR==FNR{vals[\$2];next} (\$1) in vals' ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim \
${path_gwas}resources/gnomad/gnomad_freq_edited > ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_gnomad_variants"
# Job <697020> is submitted to queue <normal>.


### R

gnomad<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_gnomad_variants",sep=""),head=F)
colnames(gnomad)<-c("SNP","CHROM","POS","REF","ALT","AF","AF_nfe","AF_afr","AF_amr","AF_eas","AF_sas","AF_asj")
gnomad$AF_nfe<-as.numeric(as.character(gnomad$AF_nfe))

vfreq<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_eur_ctr.frq.counts",sep=""),head=T)
vfreq$freq_alt<-vfreq$C1/(vfreq$C1+vfreq$C2)
vfreq$freq_alt<-as.numeric(vfreq$freq_alt)

dim(vfreq[which(!vfreq$SNP %in% gnomad$SNP),])
# 
summary(vfreq[which(!vfreq$SNP %in% gnomad$SNP),"freq_alt"])
# 
dim(vfreq[which( (!vfreq$SNP %in% gnomad$SNP) & ((vfreq$freq_alt<0.001) | (vfreq$freq_alt>0.999)) ),])
# 

table(vfreq[which(!vfreq$SNP %in% gnomad$SNP),"CHR"])
# 

all<-merge(vfreq,gnomad,by="SNP",sort=F)

table(all$A2,all$REF)
#

table(all$A1,all$ALT)
#


remove1<-all[which( (all$freq_alt>=0.45 & all$freq_alt<=0.55) & 
                      ( (all$A1=="G" & all$A2=="C") | (all$A1=="C" & all$A2=="G") | 
                          (all$A1=="A" & all$A2=="T") | (all$A1=="T" & all$A2=="A") ) ),]
dim(remove1)
# 

# remove A/T G/C that cannot be evaluated:
remove2<-vfreq[which(!(vfreq$SNP %in% gnomad$SNP) & ( (vfreq$A1=="G" & vfreq$A2=="C") | (vfreq$A1=="C" & vfreq$A2=="G") | 
                                                        (vfreq$A1=="A" & vfreq$A2=="T") | (vfreq$A1=="T" & vfreq$A2=="A") ) ),]
dim(remove2)
# 
remove<-rbind(remove1[,"SNP",drop=F],remove2[,"SNP",drop=F])
dim(remove)
# [1] 135   1

flip<-all[which( (all$freq_alt<0.45 | all$freq_alt>0.55) & 
                   ( (all$A1=="G" & all$A2=="C") | (all$A1=="C" & all$A2=="G") | 
                       (all$A1=="A" & all$A2=="T") | (all$A1=="T" & all$A2=="A") ) ),]
flip<-flip[which( (flip$freq_alt>0.5 & flip$AF_nfe<0.5) | (flip$freq_alt<0.5 & flip$AF_nfe>0.5) ),]
dim(flip)
# [1] 1699   19

write.table(remove[,"SNP",drop=F],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_variants_to_remove_AT_CG_rioux_igenomed_gsa",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
write.table(flip[,"SNP",drop=F],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_variants_to_flip_AT_CG_rioux_igenomed_gsa",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


##########################################################################################
# 18.2 REMOVE VARIANTS WE CANNOT BE SURE ARE IN THE RIGHT STRAND, AND FLIP THE OTHERS

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed \
--allow-no-sex \
--exclude ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_variants_to_remove_AT_CG_rioux_igenomed_gsa \
--flip ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_variants_to_flip_AT_CG_rioux_igenomed_gsa \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip
# --flip: 1699 SNPs flipped.
# 538438 variants and 776 people pass filters and QC.
# Among remaining phenotypes, 776 are cases and 0 are controls.

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip \
--allow-no-sex \
--a2-allele ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_variants_rioux_igenomed_gsa_hg38_posstrandaligned_with_A1 \
--make-bed --out ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2
# --a2-allele: 538438 assignments made.
# 538438 variants and 776 people pass filters and QC.
# Among remaining phenotypes, 776 are cases and 0 are controls.



################################################################################################################################################################

#############################################################
# 19 - COMPARE ALLELE FREQUENCIES WITH REFERENCE POPULATION #
#############################################################

##################################################################################################
# 19.1 SPLIT SAMPLES INTO CONTINENTAL ANCESTRY GROUPS BASED ON PCA TO COMPARE FRQ - NON DUPLICATES

# EUR - keep Europeans 

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2 \
--keep-allele-order --allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_nonjewish_ancestry_samples_noDuplicates \
--freq counts \
--make-bed --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_eur

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_eur \
--keep-allele-order --allow-no-sex \
# --keep ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_controls2 \
--freq counts \
--out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_eur_ctr


### R

gnomad<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_gnomad_variants",sep=""),head=F)
colnames(gnomad)<-c("SNP","CHROM","POS","REF","ALT","AF","AF_nfe","AF_afr","AF_amr","AF_eas","AF_sas","AF_asj")
gnomad$AF_nfe<-as.numeric(as.character(gnomad$AF_nfe))

ancestry<-c("eur")

rm(freq_all)
for (i in 1:length(ancestry)) {
  
  if (ancestry[i]=="eur") {
    file_anc<-paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_",ancestry[i],"_ctr.frq.counts",sep="")
  } else {
    file_anc<-paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_",ancestry[i],".frq.counts",sep="")
  }
  
  
  if (file.exists(file_anc)) {
    
    vfreq<-read.table(file_anc,head=T)
    vfreq[,paste(ancestry[i],"freq_alt",sep="_")]<-vfreq$C1/(vfreq$C1+vfreq$C2)
    
    if(!exists("freq_all")) {
      freq_all<-vfreq[,c(1:4,8)]
    }else{
      freq_all<-merge(freq_all,vfreq[,c(2,ncol(vfreq))],by="SNP",sort=F)
    }
  }
  
}

all<-merge(freq_all,gnomad,by="SNP",sort=F,all.x=T)


table(all$A2,all$REF)

table(all$A1,all$ALT)

dim(vfreq)-dim(gnomad)
# 

table(vfreq[which(!vfreq$SNP %in% gnomad$SNP),"CHR"])
# 

ancestry<-c("eur")

for (i in 1:length(ancestry)) {
  
  all[paste(ancestry[i],"value",sep="_")]<-((all[,paste(ancestry[i],"freq_alt",sep="_"),drop=F]-all[,"AF_nfe",drop=F])^2)/
    ((all[,paste(ancestry[i],"freq_alt",sep="_"),drop=F]+all[,"AF_nfe",drop=F])*
       (2-all[,paste(ancestry[i],"freq_alt",sep="_"),drop=F]-all[,"AF_nfe",drop=F]))
  
}

summary(all$eur_value)
# 

summary(all$AF_nfe[which(is.na(all$eur_value))])
# 


tm<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_TOPMed_variants",sep=""),head=F)
tm<-tm[,c(9,8)]
colnames(tm)<-c("SNP","topmed_freq_alt")

all<-merge(all,tm,by="SNP",all.x=T,sort=F)


all$AF_nfe<-as.numeric(as.character(all$AF_nfe))

cohort<-"rioux_igenomed_gsa"

# use different thresholds by pop:

cbPalette <- c("#999999", "#E69F00", "#56B4E9")

for (i in 1:length(ancestry)) {
  
  tmp<-all[,c("CHR","POS","SNP",paste(ancestry[i],"freq_alt",sep="_"),"AF_nfe",paste(ancestry[i],"value",sep="_"),"topmed_freq_alt")]
  
  colnames(tmp)<-c("CHR","POS","SNP","freq_alt","AF_nfe","value","topmed_freq_alt")
  tmp$col<-NA
  
  if (i==1) {
    value<-0.025
  } 
  else if (i==2) {
    value<-0.025
  } else if (i==4) {
    value<-0.0375
  } else if (i==3) {
    value<-0.05
  } else {
    value<-0.075
  }
  
  tmp$col[which(tmp$value<=value)]<-paste("<=",value,sep="")
  tmp$col[which(tmp$value>value)]<-paste(">",value,sep="")
  tmp$col[which(is.na(tmp$value))]<-"monomorphic"
  list_variants_exclude<-tmp[which(tmp$value>value),"SNP",drop=F]
  
  print(ancestry[i])
  print(table(tmp$col))
  print(dim(list_variants_exclude))
  
  p3<-ggplot(tmp, aes(y=freq_alt, x=AF_nfe)) +
    geom_point(aes(colour = col)) + ylab(paste(cohort,"FRQ Alt")) + xlab("Gnomad FRQ Alt") + scale_colour_manual(values=cbPalette) +
    ggtitle(paste(cohort,toupper(ancestry[i]),sep=" "))
  
  pdf(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_plot_maf_list_variants_toexclude_maf_differs_Gnomad_value_",ancestry[i],"_",value,".pdf",sep=""),width = 6, height = 5)
  print(p3)
  dev.off()
  system(paste("cp ",path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_plot_maf_list_variants_toexclude_maf_differs_Gnomad_value_",ancestry[i],"_",value,".pdf ~/tmp_plots/",sep=""))
  
  write.table(list_variants_exclude,paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_list_variants_toexclude_maf_differs_",ancestry[i],"_gnomad",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  
}

# [1] "eur"
# <=0.025      >0.025 monomorphic 
#     536344        1195         899 


### double check with TOPMED variants that remain after using TOPMed to double check frequency:

all<-all[which(!all$SNP %in% list_variants_exclude$SNP),]

ancestry<-c("eur")

for (i in 1:length(ancestry)) {
  
  all[paste(ancestry[i],"value_tm",sep="_")]<-((all[,paste(ancestry[i],"freq_alt",sep="_"),drop=F]-all[,"topmed_freq_alt",drop=F])^2)/
    ((all[,paste(ancestry[i],"freq_alt",sep="_"),drop=F]+all[,"topmed_freq_alt",drop=F])*
       (2-all[,paste(ancestry[i],"freq_alt",sep="_"),drop=F]-all[,"topmed_freq_alt",drop=F]))
  
}

# use different thresholds by pop:

cbPalette <- c("#999999", "#E69F00", "#56B4E9")

for (i in 1:length(ancestry)) {
  
  tmp<-all[,c("CHR","POS","SNP",paste(ancestry[i],"freq_alt",sep="_"),"topmed_freq_alt",paste(ancestry[i],"value_tm",sep="_"))]
  
  colnames(tmp)<-c("CHR","POS","SNP","freq_alt","topmed_freq_alt","value")
  tmp$col<-NA
  
  if (i==1) {
    value<-0.125
  } 
  
  tmp$col[which(tmp$value<=value)]<-paste("<=",value,sep="")
  tmp$col[which(tmp$value>value)]<-paste(">",value,sep="")
  tmp$col[which(is.na(tmp$value))]<-"monomorphic"
  list_variants_exclude2<-tmp[which(tmp$value>value),"SNP",drop=F]
  
  print(ancestry[i])
  print(table(tmp$col))
  # print(dim(list_variants_exclude2))
  
  p3<-ggplot(tmp, aes(y=freq_alt, x=topmed_freq_alt)) +
    geom_point(aes(colour = col)) + ylab(paste(cohort,"FRQ Alt")) + xlab("TOPMed FRQ Alt") + scale_colour_manual(values=cbPalette) +
    ggtitle(paste(cohort,toupper(ancestry[i]),sep=" "))
  
  pdf(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_plot_maf_list_variants_toexclude_maf_differs_TOPMed_value_",ancestry[i],"_",value,".pdf",sep=""),width = 6, height = 5)
  print(p3)
  dev.off()
  system(paste("cp ",path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_plot_maf_list_variants_toexclude_maf_differs_TOPMed_value_",ancestry[i],"_",value,".pdf ~/tmp_plots/",sep=""))
  
  list_variants_exclude<-rbind(list_variants_exclude,list_variants_exclude2)
  print(paste("Final N to exclude:",nrow(list_variants_exclude)))
  write.table(list_variants_exclude,paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_list_variants_toexclude_maf_differs_",ancestry[i],"_gnomad",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  
}


# [1] "eur"
# <=0.125  >0.125 
#  537207      36

# [1] "Final N to exclude: 1231"


#########################
# 19.2 EXCLUDE VARIANTS

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2 \
--keep-allele-order --allow-no-sex \
--exclude ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_list_variants_toexclude_maf_differs_eur_gnomad \
--make-bed --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck
# 537207 variants and 776 people pass filters and QC.
# Among remaining phenotypes, 776 are cases and 0 are controls.




################################################################################################################################################################

######################################################
# 20 - EXCLUDE GSA VARIANTS FLAGGED BY GENOME STUDIO #
######################################################

# GSA array specific step

## R

bim<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck.bim",sep=""),head=F)

fl<-read.table(paste(path,"from_kyle/list_variants_flagged_genome_studio_b38.tsv",sep=""),head=F,sep="\t")
fl$ID1<-paste(fl$V2,":",fl$V3,"_",fl$V4,"_",fl$V5,sep="")
fl$ID2<-paste(fl$V2,":",fl$V3,"_",fl$V5,"_",fl$V4,sep="")
fl$X<-paste(fl$V2,":",fl$V3,sep="")

fl$a1<-""
fl$a1[which(fl$V4=="A")]<-"T"
fl$a1[which(fl$V4=="T")]<-"A"
fl$a1[which(fl$V4=="C")]<-"G"
fl$a1[which(fl$V4=="G")]<-"C"

fl$a2<-""
fl$a2[which(fl$V5=="A")]<-"T"
fl$a2[which(fl$V5=="T")]<-"A"
fl$a2[which(fl$V5=="C")]<-"G"
fl$a2[which(fl$V5=="G")]<-"C"

fl$ID3<-paste(fl$V2,":",fl$V3,"_",fl$a1,"_",fl$a2,sep="")
fl$ID4<-paste(fl$V2,":",fl$V3,"_",fl$a2,"_",fl$a1,sep="")


excl<-bim[which( (bim$V2 %in% fl$ID1) | (bim$V2 %in% fl$ID2) | (bim$V2 %in% fl$ID3) | (bim$V2 %in% fl$ID4)),]
dim(excl)
# 42541

write.table(excl[,"V2",drop=F],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/list_gsa_variants_flagged_to_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#########


/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck \
--keep-allele-order --allow-no-sex \
--exclude ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_gsa_variants_flagged_to_exclude \
--make-bed --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged
# 494666 variants and 776 people pass filters and QC.
# Among remaining phenotypes, 776 are cases and 0 are controls.



################################################################################################################################################################

########################################################
# 21 - SPLIT COHORTS INTO BROAD CONTINENTAL ANCESTRIES #
########################################################

#####################################
# 21.1 SPLIT (INCLUDING Duplicates)

ancestry=(eur)

for i in ${ancestry[@]}
do

# EUR
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged \
--keep-allele-order --allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_${i}_ancestry_samples_withDuplicates \
--freq \
--make-bed --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${i}

# non-EUR
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged \
--keep-allele-order --allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/pca_1000gp/list_${i}_ancestry_samples_withDuplicates \
--freq \
--make-bed --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_non${i}

done

wc -l ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_*.fam
#   754 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur.fam
# 22 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur.fam
# 776 total

#####################################################
# 21.2 SPLIT (Non Duplicates) - just to get exact N

ancestry=(eur amr afr sas eas)

for i in ${ancestry[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged \
--keep-allele-order --allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_${i}_ancestry_samples_noDuplicates \
--make-bed --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${i}_nodup
done

wc -l ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_*_nodup.fam
# 1 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_afr_nodup.fam
# 21 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_amr_nodup.fam
# 742 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nodup.fam
# 764 total



#####################################################
# 21.3 REMOVE MONOMORPHIC VARIANTS PER BATCH

### R

ancestry<-c("eur","noneur")

for (i in 1:length(ancestry)) {
  frq<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_",ancestry[i],".frq",sep=""),
                  head=T)
  mono<-frq[which(frq$MAF==0),]
  print(paste("N Monomorphic ",ancestry[i],": ",nrow(mono),sep=""))
  write.table(mono[,"SNP",drop=F],paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_list_variants_toexclude_",ancestry[i],"_monomorphic",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  rm(mono)
}

# [1] "N Monomorphic eur: 1687"
# [1] "N Monomorphic noneur: 102193"



#####################################
# 21.1 SPLIT (INCLUDING Duplicates)

ancestry=(eur noneur)

for i in ${ancestry[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${i} \
--keep-allele-order --allow-no-sex \
--exclude ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_list_variants_toexclude_${i}_monomorphic \
--freq \
--make-bed --out ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${i}_nomonom
done

wc -l ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_*_nomonom.fam
# 754 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam
# 22 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam
# 776 total

wc -l ${path_gwas}/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_*_nomonom.bim
# 492979 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim
# 392473 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim

##################################################################################################################################################################

######################
# 21.- CREATE  FILES #
######################

###########################################
# 21.1 Convert ped/map files to VCF files

mkdir ${path_gwas}imputation_ready/rioux_igenomed_gsa/
  
  ancestry=(eur noneur)

for i in ${ancestry[@]}
do
for chr in {1..23}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${i}_nomonom \
--allow-no-sex \
--keep-allele-order \
--chr ${chr} \
--output-chr chr26 --recode vcf-iid --out ${path_gwas}imputation_ready/rioux_igenomed_gsa/rioux_igenomed_gsa_${i}_chr${chr}
done
done


###########################################
# 21.2 Create a sorted vcf.gz

for i in ${ancestry[@]}
do
for chr in {1..23}
do
/path/to/software/bcftools-1.9/./bcftools sort ${path_gwas}imputation_ready/rioux_igenomed_gsa/rioux_igenomed_gsa_${i}_chr${chr}.vcf \
-Oz -o ${path_gwas}imputation_ready/rioux_igenomed_gsa/rioux_igenomed_gsa_${i}_chr${chr}.vcf.gz
done
done


# rm intermediate files:
rm ${path_gwas}imputation_ready/rioux_igenomed_gsa/rioux_igenomed_gsa_*_chr*.vcf
rm ${path_gwas}imputation_ready/rioux_igenomed_gsa/*.hh
rm ${path_gwas}imputation_ready/rioux_igenomed_gsa/*.log
rm ${path_gwas}imputation_ready/rioux_igenomed_gsa/*.nosex



##################################################################################################################################################################

##############################
# 22.- IMPUTATION SERVER QC  #
##############################

###########################################
# 22.1 Submit files:

python /path/to/user/scripts/IIBDGC/topmed_qc_eur.py rioux_igenomed_gsa
# QC  rioux_igenomed_gsa eur job-20201010-211949-685

python /path/to/user/scripts/IIBDGC/topmed_qc_noneur.py rioux_igenomed_gsa
# QC  rioux_igenomed_gsa noneur job-20201010-211954-583

###########################################
# 22.2 Download files:

path_gwas=/path/to/ibdgwas/IIBDGC/
  mkdir ${path_gwas}imputed/rioux_igenomed_gsa/ 
  mkdir ${path_gwas}imputed/rioux_igenomed_gsa/qc/
  cd ${path_gwas}imputed/rioux_igenomed_gsa/qc/
  
  eth=eur 
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/160248/1cf9a78ae9841467b7ed8885160bdb73 | bash
mv chunks-excluded.txt chunks-excluded_${eth}.txt
mv snps-excluded.txt snps-excluded_${eth}.txt
mv typed-only.txt typed-only_${eth}.txt

eth=noneur
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/160269/3fa1b885f1b16a0b3bcece69d6d90793 | bash
mv chunks-excluded.txt chunks-excluded_${eth}.txt
mv snps-excluded.txt snps-excluded_${eth}.txt
mv typed-only.txt typed-only_${eth}.txt

##### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

bim<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim",sep=""),head=F)
files_list<-list.files(paste(path,"imputed/rioux_igenomed_gsa/qc/",sep=""))

files_snp_exclude<-files_list[grep("snps-excluded",files_list)]

for (i in 1:length(files_snp_exclude)){
  tmp<-fread(paste(path,"imputed/rioux_igenomed_gsa/qc/",files_snp_exclude[i],sep=""),head=T,sep="\t",fill=TRUE)
  tmp<-as.data.frame(tmp)
  colnames(tmp)<-c("Position","FilterType","Info")
  if(i==1){
    snp_excl<-tmp
  } else {
    snp_excl<-rbind(snp_excl,tmp)
  }
}

snp_excl<-snp_excl[!duplicated(snp_excl$Position),]

table(snp_excl$FilterType)
# Allele mismatch   Low call rate 
# 109              23

snp_excl<-snp_excl[which(snp_excl$FilterType=="Allele mismatch"),]


files_typed_only<-files_list[grep("typed-only",files_list)]

for (i in 1:length(files_typed_only)){
  tmp<-fread(paste(path,"imputed/rioux_igenomed_gsa/qc/",files_typed_only[i],sep=""),head=T,sep="\t",fill=TRUE)
  tmp<-as.data.frame(tmp)
  colnames(tmp)<-c("Position")
  if(i==1){
    snp_tyonly<-tmp
  } else {
    snp_tyonly<-rbind(snp_tyonly,tmp)
  }
}

snp_tyonly<-snp_tyonly[!duplicated(snp_tyonly$Position),,drop=F]
dim(snp_tyonly)
# [1] 14905   1


list_exclude<-rbind(snp_tyonly,snp_excl[,1,drop=F])

list_exclude$SNP<-gsub(":","_",list_exclude$Position)
list_exclude$SNP<-sub("_",":",list_exclude$SNP)
list_exclude$SNP<-sub("chrX","chr23",list_exclude$SNP)
list_exclude$SNP<-sub("chr","",list_exclude$SNP)

eur<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim",sep=""),head=F)
noneur<-read.table(paste(path,"pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim",sep=""),head=F)

dim(list_exclude)
# [1] 15014     2

dim(list_exclude[which(list_exclude$SNP %in% eur$V2),])
# [1] 14478     2

dim(list_exclude[which(list_exclude$SNP %in% noneur$V2),])
# [1] 1386     2

write.table(list_exclude[,"SNP",drop=F],paste(path,"imputed/rioux_igenomed_gsa/qc/list_variants_to_exclude_eur",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


###########################################
# 22.2 Exclude variants from files and Convert ped/map files to VCF files

ancestry=(eur noneur)

for i in ${ancestry[@]}
do
for chr in {1..23}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${i}_nomonom \
--allow-no-sex \
--keep-allele-order \
--exclude ${path_gwas}imputed/rioux_igenomed_gsa/qc/list_variants_to_exclude_eur \
--chr ${chr} \
--output-chr chr26 --recode vcf-iid --out ${path_gwas}imputation_ready/rioux_igenomed_gsa/rioux_igenomed_gsa_${i}_chr${chr}
done
done


###########################################
# 22.3 Create a final sorted vcf.gz

for i in ${ancestry[@]}
do
for chr in {1..23}
do
/path/to/software/bcftools-1.9/./bcftools sort ${path_gwas}imputation_ready/rioux_igenomed_gsa/rioux_igenomed_gsa_${i}_chr${chr}.vcf \
-Oz -o ${path_gwas}imputation_ready/rioux_igenomed_gsa/rioux_igenomed_gsa_${i}_chr${chr}.vcf.gz
done
done

# rm intermediate files:
rm ${path_gwas}imputation_ready/rioux_igenomed_gsa/rioux_igenomed_gsa_*_chr*.vcf
rm ${path_gwas}imputation_ready/rioux_igenomed_gsa/*.hh
rm ${path_gwas}imputation_ready/rioux_igenomed_gsa/*.log
rm ${path_gwas}imputation_ready/rioux_igenomed_gsa/*.nosex



##################################################################################################################################################################

###################
# 23.- IMPUTATION #
###################

###########################################
# 23.1 Submit files:

python /path/to/user/scripts/IIBDGC/topmed_imputation_eur.py rioux_igenomed_gsa
# Imputation  rioux_igenomed_gsa eur job-20201010-215852-298

python /path/to/user/scripts/IIBDGC/topmed_imputation_noneur.py rioux_igenomed_gsa
# Imputation  rioux_igenomed_gsa noneur job-20201011-070139-322


###########################################
# 22.1 Download files:

path_gwas=/path/to/ibdgwas/IIBDGC/
  cd ${path_gwas}imputed/rioux_igenomed_gsa/
  
  mkdir eur
mkdir noneur

## EUR SAMPLES - download commands

cd ${path_gwas}imputed/rioux_igenomed_gsa/eur/
  curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/160374/b68d1c8c89d9869bae3cc987196b8510 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/160377/b2e0ac77ccfbe4a8fe5da7e9861e0121 | bash
# curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/160376/4b6ba7eef88ec8db875fada97e6a26d8 | bash

## NONEUR SAMPLES - download commands

cd ${path_gwas}imputed/rioux_igenomed_gsa/noneur/
  curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/160437/f0b4ce82821b5db8a7313a5d76a2d3e1 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/160440/8f2800a0c53f75d2d301a124f9d49c74 | bash
# curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/160439/8e7ebe0e27a169a5e170942fe676060b | bash

################################################################
# Create a script to download all imputation files in one go:

cd ${path_gwas}imputed/rioux_igenomed_gsa/
  
  ##########
# EUR:
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/160376/4b6ba7eef88ec8db875fada97e6a26d8 > tmp_eur

# add path to eur file:
sed -i '2 i\path=\/lustre\/scratch123/hgi\/projects\/ibdgwas\/IIBDGC\/imputed\/rioux_igenomed_gsa\/eur\/' tmp_eur

##########
# NON-EUR:
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/160439/8e7ebe0e27a169a5e170942fe676060b > tmp_noneur

# remove hearder:
sed -i -e '1,4d' tmp_noneur

# add path to noneur file:
sed -i '1 i\path=\/lustre\/scratch123/hgi\/projects\/ibdgwas\/IIBDGC\/imputed\/rioux_igenomed_gsa\/noneur\/' tmp_noneur


cat tmp_eur tmp_noneur > download_rioux_igenomed_gsa_eur_noneur
rm tmp_*
  
  sed -i -e 's/-o .*.zip --create-dirs/-P ${path}/g' download_rioux_igenomed_gsa_eur_noneur
sed -i -e 's/curl -L/wget --tries=75 -c/g' download_rioux_igenomed_gsa_eur_noneur
more download_rioux_igenomed_gsa_eur_noneur

nohup bash /path/to/ibdgwas/IIBDGC/imputed/rioux_igenomed_gsa/download_rioux_igenomed_gsa_eur_noneur &
  # downloaded
  
  
  ###########################################
# 22.3 Unzip files:

path_gwas="/path/to/ibdgwas/IIBDGC/"
cd ${path_gwas}imputed/rioux_igenomed_gsa/eur/
  for i in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${i}.zip;done
unzip -P PkSw6GlBg5VGuc chr_X.zip

# cd ${path_gwas}imputed/rioux_igenomed_gsa/noneur/
# for i in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${i}.zip;done
# unzip -P PkSw6GlBg5VGuc chr_X.zip


###########################################
# 22.4 create imputation summary files and plots:

path_gwas="/path/to/ibdgwas/IIBDGC/"

mkdir ${path_gwas}imputed/rioux_igenomed_gsa/eur/plots/
  mkdir ${path_gwas}imputed/rioux_igenomed_gsa/noneur/plots/
  
  for i in {1..23}; do bash ${path_gwas}scripts/post_imputation_summary_1.sh rioux_igenomed_gsa ${i};done
# Job <862340..862364> is submitted to queue <normal>.

bash ${path_gwas}scripts/post_imputation_summary_2.sh rioux_igenomed_gsa
# Job <862365> is submitted to queue <normal>.

cp ${path_gwas}imputed/rioux_igenomed_gsa/eur/plots/*_EmpRsq_mean.pdf ~/tmp_plots/
  
  
  
  ##########################################
## 23.5 Run post-imputation summary:

# see post_imputation_summary_2.R

wc -l ${path_gwas}imputed/rioux_igenomed_gsa/qc/list_variants_0.5_EmpRsq_toexclude
# 1199 /path/to/ibdgwas/IIBDGC/imputed/rioux_igenomed_gsa/qc/list_variants_low_EmpRsq_toexclude

cat ${path_gwas}imputed/rioux_igenomed_gsa/qc/list_variants_0.5_EmpRsq_toexclude ${path_gwas}imputed/rioux_igenomed_gsa/qc/list_variants_to_exclude_eur > \
${path_gwas}imputed/rioux_igenomed_gsa/qc/list_variants_toexclude

###########################################
## 23.6 Exclude variants EmpRsq < 0.25

ancestry=(eur noneur)

for i in ${ancestry[@]}
do
for chr in {1..23}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${i}_nomonom \
--allow-no-sex \
--a2-allele ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim 6 2 \
--flip ${path_gwas}imputed/rioux_igenomed_gsa/qc/list_variants_negative_EmpR_toflip \
--exclude ${path_gwas}imputed/rioux_igenomed_gsa/qc/list_variants_toexclude \
--chr ${chr} \
--output-chr chr26 --recode vcf-iid --out ${path_gwas}imputation_ready/rioux_igenomed_gsa/rioux_igenomed_gsa_${i}_chr${chr}
done
done


###########################################
# 23.7 Create a final sorted vcf.gz

for i in ${ancestry[@]}
do
for chr in {1..23}
do
/path/to/software/bcftools-1.9/./bcftools sort ${path_gwas}imputation_ready/rioux_igenomed_gsa/rioux_igenomed_gsa_${i}_chr${chr}.vcf \
-Oz -o ${path_gwas}imputation_ready/rioux_igenomed_gsa/rioux_igenomed_gsa_${i}_chr${chr}.vcf.gz
done
done


# rm intermediate files:
rm ${path_gwas}imputation_ready/rioux_igenomed_gsa/rioux_igenomed_gsa_*_chr*.vcf
rm ${path_gwas}imputation_ready/rioux_igenomed_gsa/*.hh
rm ${path_gwas}imputation_ready/rioux_igenomed_gsa/*.log
rm ${path_gwas}imputation_ready/rioux_igenomed_gsa/*.nosex
# rm ${path_gwas}imputed/rioux_igenomed_gsa/eur/*.info.gz
# rm ${path_gwas}imputed/rioux_igenomed_gsa/eur/*.dose.vcf.gz

###########################################
# 23.8 Submit imputation files:

python /path/to/user/scripts/IIBDGC/topmed_imputation_eur.py rioux_igenomed_gsa
# Imputation  rioux_igenomed_gsa eur job-20201105-164335-456

python /path/to/user/scripts/IIBDGC/topmed_imputation_noneur.py rioux_igenomed_gsa
# Imputation  rioux_igenomed_gsa noneur job-20201105-192039-642


###########################################
# 23.9 Remove previous files:

# temporarily move them (just test for Norway, Belgium 2, Italy):

# cd ${path_gwas}imputed/rioux_igenomed_gsa/eur/
#   mkdir old/
#   mv *.zip  old/
#   mv *.log  old/
#   mv *.info.gz old/
#   mv *.dose.vcf.gz old/
#   mv snps-excluded.txt old/
#   mv typed-only.txt old/
#   mv chunks-excluded.txt old/
#   mv *.info.gz old/
#   mv *.dose.vcf.gz old/
# 
#   cd ${path_gwas}imputed/rioux_igenomed_gsa/noneur/
#   mkdir old/
#   mv *.zip  old/
#   mv *.log  old/
#   mv snps-excluded.txt old/
#   mv typed-only.txt old/
#   mv chunks-excluded.txt old/


###########################################
# 23.10 Download new files:

path_gwas=/path/to/ibdgwas/IIBDGC/
  cd ${path_gwas}imputed/rioux_igenomed_gsa/
  
  
  ## EUR SAMPLES
  
  cd ${path_gwas}imputed/rioux_igenomed_gsa/eur/
  curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/177013/b1ca59a37c64f357de0443244da47631 | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/177016/5ddcf066680d5631adc0b426c9d247c8 | bash
# curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/177015/5a990695dd1413399c30c0915fb90b66 | bash

## NONEUR SAMPLES

cd ${path_gwas}imputed/rioux_igenomed_gsa/noneur/
  curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/177181/cd67f7302d72fa9d4fed47672eb87dfa | bash
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/177184/2bc3d4c9b2ba8d95e473cdc790fc77b0 | bash
# curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/177183/c500cc08c575822674be4cb86f604beb | bash


################################################################
# Create a script to download all imputation files in one go:

cd ${path_gwas}imputed/rioux_igenomed_gsa/
  rm nohup.out

##########
# EUR:
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/177015/5a990695dd1413399c30c0915fb90b66 > tmp_eur

# add path to eur file:
sed -i '2 i\path=\/lustre\/scratch123/hgi\/projects\/ibdgwas\/IIBDGC\/imputed\/rioux_igenomed_gsa\/eur\/' tmp_eur

##########
# NON-EUR:
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/177183/c500cc08c575822674be4cb86f604beb > tmp_noneur

# remove hearder:
sed -i -e '1,4d' tmp_noneur

# add path to noneur file:
sed -i '1 i\path=\/lustre\/scratch123/hgi\/projects\/ibdgwas\/IIBDGC\/imputed\/rioux_igenomed_gsa\/noneur\/' tmp_noneur


cat tmp_eur tmp_noneur > download_rioux_igenomed_gsa_eur_noneur
rm tmp_*
  
  sed -i -e 's/-o .*.zip --create-dirs/-P ${path}/g' download_rioux_igenomed_gsa_eur_noneur
sed -i -e 's/curl -L/wget --tries=75 -c/g' download_rioux_igenomed_gsa_eur_noneur
more download_rioux_igenomed_gsa_eur_noneur


nohup bash /path/to/ibdgwas/IIBDGC/imputed/rioux_igenomed_gsa/download_rioux_igenomed_gsa_eur_noneur &
  # downloaded
  
  
  ###########################################
# 23.11 Unzip files:

cd ${path_gwas}imputed/rioux_igenomed_gsa/eur/
  for i in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${i}.zip;done
unzip -P PkSw6GlBg5VGuc chr_X.zip

# cd ${path_gwas}imputed/rioux_igenomed_gsa/noneur/
# for i in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${i}.zip;done
# unzip -P PkSw6GlBg5VGuc chr_X.zip

# CONTINUE HERE

###########################################
# 23.12 create imputation summary files and plots:

path_gwas="/path/to/ibdgwas/IIBDGC/"

# final rainfall plots
for i in {1..23}; do bash ${path_gwas}scripts/post_imputation_summary_1.sh rioux_igenomed_gsa ${i};done
# Job <296456..296487> is submitted to queue <normal>.

bash ${path_gwas}scripts/post_imputation_summary_3.sh rioux_igenomed_gsa
# Job <296488> is submitted to queue <normal>.

less /path/to/ibdgwas/IIBDGC/scripts/logs/post_imputation_summary_3.R_output_rioux_igenomed_gsa_allchr 


################################################################
# 23.13 Create a list of variants that do not pass HWE in EUR

# for chr in {1..22}
# do
# /path/to/software/./plink2 \
# --vcf ${path_gwas}imputed/rioux_igenomed_gsa/eur/chr${chr}.dose.vcf.gz \
# --double-id \
# --keep ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_controls \
# --hardy \
# --out ${path_gwas}imputed/rioux_igenomed_gsa/eur/chr${chr}_hardy_ctr
# done



for chr in {1..22}
do
/path/to/software/./plink2 \
--vcf ${path_gwas}imputed/rioux_igenomed_gsa/eur/chr${chr}.dose.vcf.gz \
--double-id \
--remove ${path_gwas}pre_imputation/QC/rioux_igenomed_gsa/list_controls \
--hardy \
--out ${path_gwas}imputed/rioux_igenomed_gsa/eur/chr${chr}_hardy_cases
done

# merge files

# cat ${path_gwas}imputed/rioux_igenomed_gsa/eur/chr*_hardy_ctr.hardy > ${path_gwas}imputed/rioux_igenomed_gsa/eur/allchr_hardy_ctr.hardy
cat ${path_gwas}imputed/rioux_igenomed_gsa/eur/chr*_hardy_cases.hardy > ${path_gwas}imputed/rioux_igenomed_gsa/eur/allchr_hardy_cases.hardy

# remove intermediate files
# rm ${path_gwas}imputed/rioux_igenomed_gsa/eur/chr*_hardy_ctr.*
rm ${path_gwas}imputed/rioux_igenomed_gsa/eur/chr*_hardy_cases.*
  
  # # create list with variants to exclude from analysis
  # awk ' $10 < 1E-5 ' ${path_gwas}imputed/rioux_igenomed_gsa/eur/allchr_hardy_ctr.hardy \
  # > ${path_gwas}imputed/rioux_igenomed_gsa/eur/list_variants_fail_hwe_ctr
  # wc -l ${path_gwas}imputed/rioux_igenomed_gsa/eur/list_variants_fail_hwe_ctr
  # # 
  
  awk ' $10 < 1E-12 ' ${path_gwas}imputed/rioux_igenomed_gsa/eur/allchr_hardy_cases.hardy \
> ${path_gwas}imputed/rioux_igenomed_gsa/eur/list_variants_fail_hwe_cases
wc -l ${path_gwas}imputed/rioux_igenomed_gsa/eur/list_variants_fail_hwe_cases
# 7 /path/to/ibdgwas/IIBDGC/imputed/rioux_igenomed_gsa/eur/list_variants_fail_hwe_cases

cat ${path_gwas}imputed/rioux_igenomed_gsa/eur/list_variants_fail_hwe_ctr ${path_gwas}imputed/rioux_igenomed_gsa/eur/list_variants_fail_hwe_cases \
| cut -f 2 | awk '!seen[$0]++'> ${path_gwas}imputed/rioux_igenomed_gsa/eur/list_variants_post_imputation_hwe_toremove

wc -l ${path_gwas}imputed/rioux_igenomed_gsa/eur/list_variants_post_imputation_hwe_toremove
# 7 /path/to/ibdgwas/IIBDGC/imputed/rioux_igenomed_gsa/eur/list_variants_post_imputation_hwe_toremove

# continue here

# remove intermediate files
rm ${path_gwas}imputed/rioux_igenomed_gsa/eur/allchr_hardy_ctr.hardy
rm ${path_gwas}imputed/rioux_igenomed_gsa/eur/allchr_hardy_cases.hardy
rm ${path_gwas}imputed/rioux_igenomed_gsa/eur/list_variants_fail_hwe_cases
rm ${path_gwas}imputed/rioux_igenomed_gsa/eur/list_variants_fail_hwe_ctr

wc -l ${path_gwas}imputed/rioux_igenomed_gsa/eur/plots/rioux_igenomed_gsa_eur_list_variants_maf_0.001_rsq_0.04_updated.tsv
# 17799407 /path/to/ibdgwas/IIBDGC/imputed/rioux_igenomed_gsa/eur/plots/rioux_igenomed_gsa_eur_list_variants_maf_0.001_rsq_0.04_updated.tsv


grep -Fvf ${path_gwas}imputed/rioux_igenomed_gsa/eur/list_variants_post_imputation_hwe_toremove ${path_gwas}imputed/rioux_igenomed_gsa/eur/plots/rioux_igenomed_gsa_eur_list_variants_maf_0.001_rsq_0.04_updated.tsv \
> ${path_gwas}imputed/rioux_igenomed_gsa/eur/plots/rioux_igenomed_gsa_eur_list_variants_maf_0.001_rsq_0.04_updated_hwe.tsv
wc -l ${path_gwas}imputed/rioux_igenomed_gsa/eur/plots/rioux_igenomed_gsa_eur_list_variants_maf_0.001_rsq_0.04_updated_hwe.tsv
# 17799400 /path/to/ibdgwas/IIBDGC/imputed/rioux_igenomed_gsa/eur/plots/rioux_igenomed_gsa_eur_list_variants_maf_0.001_rsq_0.04_updated_hwe.tsv

###########################################
# 23.14 Remove previous files:

cd ${path_gwas}imputed/rioux_igenomed_gsa/eur/old/
  rm *.zip
rm *.log

cd ${path_gwas}imputed/rioux_igenomed_gsa/noneur/old/
  rm *.zip
rm *.log


# lfs quota -hg ibdgwas /path/to/project

########################################################################################################################################################################
########################################################################################################################################################################


# Clean qc folder, keep start and end files:

path_gwas=/path/to/ibdgwas/IIBDGC/
  study=rioux_igenomed_gsa

cd ${path_gwas}/pre_imputation/QC/${study}/
  rm -r logs/
  
  # keep :
  
  mkdir keep/
  mv ${study}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_*_nomonom* keep/
  mv *_hg19.* keep/
  mv *_hg19_edited.fam keep/
  cd keep/
  ls -la

cd ..

ls -la 1000GP_EUR_b37_${study}_variants*
  rm 1000GP_EUR_b37_${study}_variants*
  
  ls -la ${study}_hg17_noind*
  rm ${study}_hg17_noind*
  ls -la ${study}_hg18_noind*
  rm ${study}_hg18_noind*
  ls -la ${study}_hg19_*
  rm ${study}_hg19_*
  ls -la ${study}_hg19_postqc_lifted_hg38*
  rm ${study}_hg19_postqc_lifted_hg38*
  ls -la ${study}_postqc_lifted*
  rm ${study}_postqc_lifted*
  rm autosomal_het_*
  ls -la hwe_*
  rm hwe_*
  rm ${study}_hg19.vcf
ls -la *_postqc_preimp_*
  rm *_postqc_preimp_*
  ls -la ${study}_check_sex_*
  rm ${study}_check_sex_*
  rm ${study}_posstrandaligned.vcf.gz
rm ${study}_posstrandaligned_2.vcf.gz













