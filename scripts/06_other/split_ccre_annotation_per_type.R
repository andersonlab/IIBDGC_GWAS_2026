# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# test
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -n 2 R

# How to submit:

# cd /path/to/project
# files_ccre=($(ls | grep bed | sed 's/\.bed\.gz//g' ))

# MEM=1200
# for i in ${files_ccre[@]}
# do
# bsub -J"ccre" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 2 \
# -e ${path_gwas}post_imputation/log/${i}_ccre_stderr \
# -o ${path_gwas}post_imputation/log/${i}_ccre_stdout \
# "Rscript ~/git/IIBDGC_GWAS/scripts/other/split_ccre_annotation_per_type.R ${i} > \
# ${path_gwas}post_imputation/2022/log/split_ccre_annotation_per_type${i}.Rout"
# done

library(data.table)

# provide trait (in this case combination of database + trait) as input
args<-commandArgs(trailingOnly=TRUE)
ccre_id<-args[1]

# for testing purposes:
# ccre_id<-"ENCFF478PHY"
# ccre_id<-"ENCFF999DTZ"
# ccre_id<-"ENCFF047FTG"
# ccre_id<-"ENCFF271XXO"

path<-"/path/to/ibdgwas/IIBDGC/"
dat<-fread("/path/to/project",sep="\t",head=T)

dat<-dat[which(dat$'Audit WARNING'!='derived from revoked file'),]
dat<-dat[which(dat$'Biosample term name'!=''),]

dat$biosample<-gsub(", ",",",dat$'Biosample term name')
dat$biosample<-gsub(" ","_",dat$biosample)
dat$biosample<-gsub("-","_",dat$biosample)
dat$biosample<-gsub("/","_",dat$biosample)
dat$biosample<-gsub("'","",dat$biosample)

# treat as binomial data, create a tmp file:
dat<-dat[which(dat$"File accession"==ccre_id),]

print(ccre_id)

# CCREs are colored and labeled according to classification by regulatory signature:
# Color		UCSC label	ENCODE classification	ENCODE label
# red	prom	promoter-like signature	PLS
# orange	enhP	proximal enhancer-like signature	pELS
# yellow	enhD	distal enhancer-like signature	dELS
# pink	K4m3	DNase-H3K4me3	DNase-H3K4me3
# blue	CTCF	CTCF-only	CTCF-only


print(paste(dat$biosample,dat$"File accession"))

tmp<-fread(paste("/path/to/project",ccre_id,".bed.gz",sep=""))
tmp$V10<-gsub(",","_",tmp$V10)
tmp$V10<-gsub("-","_",tmp$V10)
tmp$V10<-gsub("/","_",tmp$V10)
conditions<-levels(as.factor(tmp$V10))
conditions<-conditions[!grepl(c("Unclassified"),conditions)]

print(conditions)
  
for (j in 1:length(conditions)) {
    
    print(conditions[j])
    tmp1<-tmp[which(tmp$V10 %in% conditions[j]),]
    tmp1$chr<-gsub("chr","",tmp1$V1)
    tmp1<-tmp1[which(tmp1$chr %in% c(seq(1:22),"X")),]
    tmp1$chr<-as.numeric(tmp1$chr)
    tmp1$chr[which(tmp1$V1=="chrX")]<-23
    tmp1<-tmp1[order(tmp1$chr,tmp1$V2,decreasing=F),]
    
    write.table(tmp1[,c(1:3)],paste(path,"post_imputation/analysis/metaanalysis/annotation/encode/cre_",dat$biosample,"_",dat$"File accession","_",conditions[j],"_tmp.bed",sep=""),
                col.names=F,row.names=F,quote=F,sep="\t")

    # system(paste("bedtools intersect -a ",path,"post_imputation/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed -b ",
    #              path,"post_imputation/analysis/metaanalysis/annotation/encode/cre_",dat$biosample[i],"_",conditions[j],"_tmp.bed -c > ",path,
    #              "post_imputation/analysis/metaanalysis/annotation/encode/allchr_",dat$biosample[i],"_",conditions[j],"_",dat$file[i],"_b38_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed",sep=""))
    # rm(tmp1)
    # system(paste("gzip ",path,
    #              "post_imputation/analysis/metaanalysis/annotation/encode/allchr_",dat$biosample[i],"_",conditions[j],"_",dat$file[i],"_b38_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed",sep=""))

    # system(paste("rm ",path,"post_imputation/analysis/metaanalysis/annotation/encode/cre_",dat$biosample[i],"_",conditions[j],"_tmp.bed",sep=""))
  
  rm(tmp1)

}

q("no")
