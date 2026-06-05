# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# harmonise nomenclature and liftover to b38:

# # singularity exec iibdgc_postprocess_10_singularity.sif

# MEM=25000
# bsub -Is -M"$MEM" -R"select[model==Intel_Platinum && mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group R \

library(data.table)
library(R.utils)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

args<-commandArgs(trailingOnly=TRUE)
study<-args[1]

# for testing purposes
# study<-"RE-LYMP(L)%"

df<-fread(paste0("/path/to/project",study,".tsv.gz"),head=T)

# replace some special characters in the name:
study_edited<-gsub("\\(","_",study)
study_edited<-gsub("\\)","_",study_edited)
study_edited<-gsub("\\(=","_",study_edited)
study_edited<-gsub("-","_",study_edited)

print(study_edited)

# keep only SNPs with imputation info>0.3 and MAF>=0.001
df<-df[which(df$INFO>=0.3 & df$ALT_FREQ>=0.001 & df$ALT_FREQ<=0.999),]

# keep only SNPs
df<-df[which(nchar(df$REF)==1 & nchar(df$ALT)==1),]

print(dim(df))
#[1] 17363385 6

df<-df[which(df$REF!=df$ALT),]
dim(df)
#[1][1] 17363385       10 

df$variant<-paste("chr",df$CHR,":",df$BP,":",df$REF,":",df$ALT,sep="")

bed<-df[,c("CHR","BP","variant","REF","ALT")]

bed<-bed[!duplicated(bed),]
dim(bed)
#[1] 17363383        5

bed$end<-bed$BP+1

bed$Pos<-format(bed$BP, scientific=F)
bed$end<-format(bed$end, scientific=F)

bed$chr<-paste("chr",bed$CHR,sep="")
setorder(bed,cols="chr","Pos")

file_out<-paste(path_gwas,"/resources/gwas_summary_statistics/",study_edited,"_buildGRCh37_b37.bed",sep="")
write.table(bed[,c("chr","Pos","end","variant")],file_out,col.names=F,row.names=F,quote=F)


#### lift positions
system(paste("/path/to/software/username/./liftOver ",file_out,
" /path/to/ibdgwas/IIBDGC/previous_qced_b38/liftover/hg19ToHg38.over.chain.gz ",file_out,
"_lifted_hg38 ",file_out,"_no_lifted_hg38",sep=""))

system(paste("cut -f 4 ",file_out,"_no_lifted_hg38 | sed '/^#/d' > ",file_out,"_no_lifted_hg38_variants_to_exclude_tmp",sep=""))
system(paste("grep alt ",file_out,"_lifted_hg38 | cut -f 4 | cat - ",file_out,"_no_lifted_hg38_variants_to_exclude_tmp > ",
             file_out,"_no_lifted_hg38_variants_to_exclude",sep=""))

system(paste("wc -l ",file_out,"_no_lifted_hg38_variants_to_exclude",sep=""))
# 6836 /path/to/ibdgwas/IIBDGC//resources/gwas_summary_statistics/RE_LYMP_L_%_buildGRCh37_b37.bed_no_lifted_hg38_variants_to_exclude

### exclude non lifted:

bed_up<-fread(paste(file_out,"_lifted_hg38",sep=""),head=F)
colnames(bed_up)[2]<-"pos_b38"
colnames(bed_up)[4]<-"variant"
dim(bed_up)
# [1] 17359722        4

dim(bed)
# [1] 17363383       7
bed<-merge(bed,bed_up[,c("variant","pos_b38")],by="variant",all.y=T,sort=F)
dim(bed)
# [1] 17359722       8


# create tfam and tped like file to convert to vcf file
bed$morgan<-0
bed$genotype<-paste(bed$REF,bed$ALT)

tped<-bed[,c("chr","variant","morgan","pos_b38","genotype")]

file_out_tped<-paste(path_gwas,"resources/gwas_summary_statistics/",study_edited,"_b38.tped",sep="")
write.table(tped,file_out_tped,
            col.names=F,row.names=F,sep="\t",quote=F)

tfam<-c("ID_1","ID_1","0","0","1","1")
tfam<-t(as.data.frame(tfam))
file_out_tfam<-paste(path_gwas,"resources/gwas_summary_statistics/",study_edited,"_b38.tfam",sep="")
write.table(tfam,file_out_tfam,
            col.names=F,row.names=F,sep="\t",quote=F)
  
# create A1|REF|REF allele
file_out_allele<-paste(path_gwas,"resources/gwas_summary_statistics/",study_edited,"_Ref",sep="")
write.table(bed[,c("variant","REF")],file_out_allele,
            col.names=F,row.names=F,sep="\t",quote=F)



# create vcf:
system(paste("/path/to/software/username/plink_linux_x86_64_20181202/./plink ",
             "--tfile ",path_gwas,"resources/gwas_summary_statistics/",study_edited,"_b38 ",
             "--allow-no-sex ",
             "--a2-allele ",path_gwas,"resources/gwas_summary_statistics/",study_edited,"_Ref ",
             "--keep-allele-order --output-chr M --recode vcf-iid ",
             "--out ",path_gwas,"resources/gwas_summary_statistics/",study_edited,"_b38",sep=""))

# # Note that most PLINK analyses treat the A1 (usually minor) allele as the reference allele, which makes sense when only biallelic variants are involved.
# # However, since it is conventional for VCF files to set the major allele as the reference allele instead

# # # double check alleles, variants:

system(paste0("bcftools +fixref ",path_gwas,
"resources/gwas_summary_statistics/",study_edited,"_b38.vcf -Oz -o ",
path_gwas,"resources/gwas_summary_statistics/",study_edited,"_b38.vcf.gz -- -f ",path_gwas,
"resources/hg38/hg38_edited.fa -m top"))

# # Double check conversion
system(paste("bcftools +fixref ",path_gwas,
"resources/gwas_summary_statistics/",study_edited,"_b38.vcf.gz -- -f ",path_gwas,
"resources/hg38/hg38_edited.fa",sep=""))


# # #### VCF to BED

system(paste("/path/to/software/username/plink_linux_x86_64_20181202/./plink --vcf ",path_gwas,
"resources/gwas_summary_statistics/",study_edited,"_b38.vcf.gz  --keep-allele-order --allow-no-sex --double-id --make-bed --out ",
path_gwas,"resources/gwas_summary_statistics/",study_edited,"_b38_2",sep=""))


# # ##############################################################
# # # 16.2 UPDATE NAME VARIANTS TO CHR:POSITION_REF_ALT in b38

system(paste("zcat ",
path_gwas,"resources/gwas_summary_statistics/",study_edited,
"_b38.vcf.gz | cut -f '1-5' | awk '{print $3,$1\":\"$2\"_\"$4\"_\"$5}' > ",path_gwas,
"resources/gwas_summary_statistics/list_variants_",study_edited,"_b38",sep=""))

########################

bim<-fread(paste(path_gwas,"/resources/gwas_summary_statistics/",study_edited,"_b38_2.bim",sep=""),head=F)
colnames(bim)[2]<-"variant"
ids<-fread(paste(path_gwas,"resources/gwas_summary_statistics/list_variants_",study_edited,"_b38",sep=""),head=F)
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
# [1] 17359722       14

rm(bim,tped,tfam)

### add freq from gnomad:

system(paste0("awk 'NR==FNR{vals[$2];next} $1 in vals' ",
path_gwas,"resources/gwas_summary_statistics/list_variants_",study_edited,"_b38  ",
path_gwas,"resources/gnomad/gnomad_freq_edited | gzip > ",
path_gwas,"resources/gwas_summary_statistics/list_variants_",study_edited,"_b38_gnomad.gz"))

##############

gnomad<-read.table(paste(path_gwas,"resources/gwas_summary_statistics/list_variants_",study_edited,"_b38_gnomad.gz",sep=""),head=F)
colnames(gnomad)<-c("SNP","CHROM","POS","REF","ALT","AF","AF_nfe","AF_afr","AF_amr","AF_eas","AF_sas","AF_asj")
gnomad$AF_nfe<-as.numeric(as.character(gnomad$AF_nfe))

bed<-merge(bed,gnomad[,c("SNP","AF_nfe")],by.x="variant_b38",by.y="SNP",all.x=T,sort=F)
rm(gnomad)

dim(df[which(df$variant %in% bed$variant),])
# [1] 17359724       10

dim(bed)
# [1] 17359722       13

head(bed)
head(df)


head(df)
df<-merge(df[,c("variant","EFFECT","SE","P","INFO","ALT_FREQ")],bed,by="variant",all.y=T)
# # rm(bed)

df<-as.data.frame(df)

df$N<-41515

# BETA IS FOR ALT allele as described in /path/to/project
df$beta_ed<-NA
df$beta_ed[which(df$REF==df$A0 & df$ALT==df$A1)]<-df$EFFECT[which(df$REF==df$A0 & df$ALT==df$A1)]
df$beta_ed[which(df$REF==df$A1 & df$ALT==df$A0)]<-df$EFFECT[which(df$REF==df$A1 & df$ALT==df$A0)]*-1

dim(df[which(is.na(df$beta_ed)),])
# [1] 19858    23

# VARIANTS WHERE REF AND ALT HAVE BEEN SWAPPED -  same effect:
df$beta_ed[which( is.na(df$beta_ed) & (df$A1=="A") & (df$REF=="T") & (df$A0=="C") & (df$ALT=="G"))]<-df$EFFECT[which( is.na(df$beta_ed) & (df$A1=="A") & (df$REF=="T") & (df$A0=="C") & (df$ALT=="G"))]
df$beta_ed[which( is.na(df$beta_ed) & (df$A1=="A") & (df$REF=="T") & (df$A0=="G") & (df$ALT=="C"))]<-df$EFFECT[which( is.na(df$beta_ed) & (df$A1=="A") & (df$REF=="T") & (df$A0=="G") & (df$ALT=="C"))]
df$beta_ed[which( is.na(df$beta_ed) & (df$A1=="C") & (df$REF=="G") & (df$A0=="A") & (df$ALT=="T"))]<-df$EFFECT[which( is.na(df$beta_ed) & (df$A1=="C") & (df$REF=="G") & (df$A0=="A") & (df$ALT=="T"))]
df$beta_ed[which( is.na(df$beta_ed) & (df$A1=="C") & (df$REF=="G") & (df$A0=="T") & (df$ALT=="A"))]<-df$EFFECT[which( is.na(df$beta_ed) & (df$A1=="C") & (df$REF=="G") & (df$A0=="T") & (df$ALT=="A"))]
df$beta_ed[which( is.na(df$beta_ed) & (df$A1=="G") & (df$REF=="C") & (df$A0=="T") & (df$ALT=="A"))]<-df$EFFECT[which( is.na(df$beta_ed) & (df$A1=="G") & (df$REF=="C") & (df$A0=="T") & (df$ALT=="A"))]
df$beta_ed[which( is.na(df$beta_ed) & (df$A1=="G") & (df$REF=="C") & (df$A0=="A") & (df$ALT=="T"))]<-df$EFFECT[which( is.na(df$beta_ed) & (df$A1=="G") & (df$REF=="C") & (df$A0=="A") & (df$ALT=="T"))]
df$beta_ed[which( is.na(df$beta_ed) & (df$A1=="T") & (df$REF=="A") & (df$A0=="C") & (df$ALT=="G"))]<-df$EFFECT[which( is.na(df$beta_ed) & (df$A1=="T") & (df$REF=="A") & (df$A0=="C") & (df$ALT=="G"))]
df$beta_ed[which( is.na(df$beta_ed) & (df$A1=="T") & (df$REF=="A") & (df$A0=="G") & (df$ALT=="C"))]<-df$EFFECT[which( is.na(df$beta_ed) & (df$A1=="T") & (df$REF=="A") & (df$A0=="G") & (df$ALT=="C"))]
  
dim(df[which(is.na(df$beta_ed)),])
# [1] 9961   23

# VARIANTS WHERE REF AND ALT HAVE BEEN SWAPPED AND FLIPPED - effect*-1:
df$beta_ed[which( is.na(df$beta_ed) & (df$A1=="A") & (df$REF=="G") & (df$A0=="C") & (df$ALT=="T"))]<-df$EFFECT[which( is.na(df$beta_ed) & (df$A1=="A") & (df$REF=="G") & (df$A0=="C") & (df$ALT=="T"))]*(-1)
df$beta_ed[which( is.na(df$beta_ed) & (df$A1=="A") & (df$REF=="C") & (df$A0=="G") & (df$ALT=="T"))]<-df$EFFECT[which( is.na(df$beta_ed) & (df$A1=="A") & (df$REF=="C") & (df$A0=="G") & (df$ALT=="T"))]*(-1)
df$beta_ed[which( is.na(df$beta_ed) & (df$A1=="C") & (df$REF=="T") & (df$A0=="A") & (df$ALT=="G"))]<-df$EFFECT[which( is.na(df$beta_ed) & (df$A1=="C") & (df$REF=="T") & (df$A0=="A") & (df$ALT=="G"))]*(-1)
df$beta_ed[which( is.na(df$beta_ed) & (df$A1=="C") & (df$REF=="A") & (df$A0=="T") & (df$ALT=="G"))]<-df$EFFECT[which( is.na(df$beta_ed) & (df$A1=="C") & (df$REF=="A") & (df$A0=="T") & (df$ALT=="G"))]*(-1)
df$beta_ed[which( is.na(df$beta_ed) & (df$A1=="G") & (df$REF=="A") & (df$A0=="T") & (df$ALT=="C"))]<-df$EFFECT[which( is.na(df$beta_ed) & (df$A1=="G") & (df$REF=="A") & (df$A0=="T") & (df$ALT=="C"))]*(-1)
df$beta_ed[which( is.na(df$beta_ed) & (df$A1=="G") & (df$REF=="T") & (df$A0=="A") & (df$ALT=="C"))]<-df$EFFECT[which( is.na(df$beta_ed) & (df$A1=="G") & (df$REF=="T") & (df$A0=="A") & (df$ALT=="C"))]*(-1)
df$beta_ed[which( is.na(df$beta_ed) & (df$A1=="T") & (df$REF=="G") & (df$A0=="C") & (df$ALT=="A"))]<-df$EFFECT[which( is.na(df$beta_ed) & (df$A1=="T") & (df$REF=="G") & (df$A0=="C") & (df$ALT=="A"))]*(-1)
df$beta_ed[which( is.na(df$beta_ed) & (df$A1=="T") & (df$REF=="C") & (df$A0=="G") & (df$ALT=="A"))]<-df$EFFECT[which( is.na(df$beta_ed) & (df$A1=="T") & (df$REF=="C") & (df$A0=="G") & (df$ALT=="A"))]*(-1)

dim(df[which(is.na(df$beta_ed)),])
# 0


# UPDATE ALT_FREQ
df$A1FREQ_ed<-NA
df$A1FREQ_ed[which(df$REF==df$A0 & df$ALT==df$A1)]<-df$ALT_FREQ[which(df$REF==df$A0 & df$ALT==df$A1)]
df$A1FREQ_ed[which(df$REF==df$A1 & df$ALT==df$A0)]<-1-(df$ALT_FREQ[which(df$REF==df$A1 & df$ALT==df$A0)])

dim(df[which(is.na(df$A1FREQ_ed)),])
# [1] 19858    23

# VARIANTS WHERE REF AND ALT HAVE BEEN SWAPPED -  same effect:
df$A1FREQ_ed[which( is.na(df$A1FREQ_ed) & (df$A1=="A") & (df$REF=="T") & (df$A0=="C") & (df$ALT=="G"))]<-df$ALT_FREQ[which( is.na(df$A1FREQ_ed) & (df$A1=="A") & (df$REF=="T") & (df$A0=="C") & (df$ALT=="G"))]
df$A1FREQ_ed[which( is.na(df$A1FREQ_ed) & (df$A1=="A") & (df$REF=="T") & (df$A0=="G") & (df$ALT=="C"))]<-df$ALT_FREQ[which( is.na(df$A1FREQ_ed) & (df$A1=="A") & (df$REF=="T") & (df$A0=="G") & (df$ALT=="C"))]
df$A1FREQ_ed[which( is.na(df$A1FREQ_ed) & (df$A1=="C") & (df$REF=="G") & (df$A0=="A") & (df$ALT=="T"))]<-df$ALT_FREQ[which( is.na(df$A1FREQ_ed) & (df$A1=="C") & (df$REF=="G") & (df$A0=="A") & (df$ALT=="T"))]
df$A1FREQ_ed[which( is.na(df$A1FREQ_ed) & (df$A1=="C") & (df$REF=="G") & (df$A0=="T") & (df$ALT=="A"))]<-df$ALT_FREQ[which( is.na(df$A1FREQ_ed) & (df$A1=="C") & (df$REF=="G") & (df$A0=="T") & (df$ALT=="A"))]
df$A1FREQ_ed[which( is.na(df$A1FREQ_ed) & (df$A1=="G") & (df$REF=="C") & (df$A0=="T") & (df$ALT=="A"))]<-df$ALT_FREQ[which( is.na(df$A1FREQ_ed) & (df$A1=="G") & (df$REF=="C") & (df$A0=="T") & (df$ALT=="A"))]
df$A1FREQ_ed[which( is.na(df$A1FREQ_ed) & (df$A1=="G") & (df$REF=="C") & (df$A0=="A") & (df$ALT=="T"))]<-df$ALT_FREQ[which( is.na(df$A1FREQ_ed) & (df$A1=="G") & (df$REF=="C") & (df$A0=="A") & (df$ALT=="T"))]
df$A1FREQ_ed[which( is.na(df$A1FREQ_ed) & (df$A1=="T") & (df$REF=="A") & (df$A0=="C") & (df$ALT=="G"))]<-df$ALT_FREQ[which( is.na(df$A1FREQ_ed) & (df$A1=="T") & (df$REF=="A") & (df$A0=="C") & (df$ALT=="G"))]
df$A1FREQ_ed[which( is.na(df$A1FREQ_ed) & (df$A1=="T") & (df$REF=="A") & (df$A0=="G") & (df$ALT=="C"))]<-df$ALT_FREQ[which( is.na(df$A1FREQ_ed) & (df$A1=="T") & (df$REF=="A") & (df$A0=="G") & (df$ALT=="C"))]
  
dim(df[which(is.na(df$A1FREQ_ed)),])
# [1] 9961   23

# VARIANTS WHERE REF AND ALT HAVE BEEN SWAPPED AND FLIPPED - effect*-1:
df$A1FREQ_ed[which( is.na(df$A1FREQ_ed) & (df$A1=="A") & (df$REF=="G") & (df$A0=="C") & (df$ALT=="T"))]<-1-(df$ALT_FREQ[which( is.na(df$A1FREQ_ed) & (df$A1=="A") & (df$REF=="G") & (df$A0=="C") & (df$ALT=="T"))])
df$A1FREQ_ed[which( is.na(df$A1FREQ_ed) & (df$A1=="A") & (df$REF=="C") & (df$A0=="G") & (df$ALT=="T"))]<-1-(df$ALT_FREQ[which( is.na(df$A1FREQ_ed) & (df$A1=="A") & (df$REF=="C") & (df$A0=="G") & (df$ALT=="T"))])
df$A1FREQ_ed[which( is.na(df$A1FREQ_ed) & (df$A1=="C") & (df$REF=="T") & (df$A0=="A") & (df$ALT=="G"))]<-1-(df$ALT_FREQ[which( is.na(df$A1FREQ_ed) & (df$A1=="C") & (df$REF=="T") & (df$A0=="A") & (df$ALT=="G"))])
df$A1FREQ_ed[which( is.na(df$A1FREQ_ed) & (df$A1=="C") & (df$REF=="A") & (df$A0=="T") & (df$ALT=="G"))]<-1-(df$ALT_FREQ[which( is.na(df$A1FREQ_ed) & (df$A1=="C") & (df$REF=="A") & (df$A0=="T") & (df$ALT=="G"))])
df$A1FREQ_ed[which( is.na(df$A1FREQ_ed) & (df$A1=="G") & (df$REF=="A") & (df$A0=="T") & (df$ALT=="C"))]<-1-(df$ALT_FREQ[which( is.na(df$A1FREQ_ed) & (df$A1=="G") & (df$REF=="A") & (df$A0=="T") & (df$ALT=="C"))])
df$A1FREQ_ed[which( is.na(df$A1FREQ_ed) & (df$A1=="G") & (df$REF=="T") & (df$A0=="A") & (df$ALT=="C"))]<-1-(df$ALT_FREQ[which( is.na(df$A1FREQ_ed) & (df$A1=="G") & (df$REF=="T") & (df$A0=="A") & (df$ALT=="C"))])
df$A1FREQ_ed[which( is.na(df$A1FREQ_ed) & (df$A1=="T") & (df$REF=="G") & (df$A0=="C") & (df$ALT=="A"))]<-1-(df$ALT_FREQ[which( is.na(df$A1FREQ_ed) & (df$A1=="T") & (df$REF=="G") & (df$A0=="C") & (df$ALT=="A"))])
df$A1FREQ_ed[which( is.na(df$A1FREQ_ed) & (df$A1=="T") & (df$REF=="C") & (df$A0=="G") & (df$ALT=="A"))]<-1-(df$ALT_FREQ[which( is.na(df$A1FREQ_ed) & (df$A1=="T") & (df$REF=="C") & (df$A0=="G") & (df$ALT=="A"))])

dim(df[which(is.na(df$A1FREQ_ed)),])
# [1]  0 23


df$value<-((df$A1FREQ_ed-df$AF_nfe)^2)/(df$A1FREQ_ed+df$AF_nfe)*(2-df$A1FREQ_ed-df$AF_nfe)

dim(df[which(df$value>0.05),])
# [1] 157815     24

df<-df[which(df$value<=0.05 | is.na(df$A1FREQ_ed)),]

dim(df)
# [1] 16380141       24

df<-df[,c("variant_b38","CHR","pos_b38","A0","A1","beta_ed","SE","P","INFO","N","A1FREQ_ed")]

colnames(df)<-c("SNP","CHR","BP","A0","A1","BETA","SE","PVALUE","INFO","N","A1FREQ")
df<-df[which(df$A1FREQ>=0.001 & df$A1FREQ<=0.999),]
dim(df)
# [1] 16380141       11

df$SNP<-gsub("_",":",df$SNP)
df$SNP<-paste("chr",df$SNP,sep="")

# remove duplicated ids:
dups<-df$SNP[duplicated(df$SNP)]
length(dups)
# [1] 94
df<-df[!df$SNP %in% dups,]
dim(df)
# [1] 16379953       11

fwrite(df,paste(path_gwas,"resources/gwas_summary_statistics/",study_edited,"_edited.tsv.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

q("no")