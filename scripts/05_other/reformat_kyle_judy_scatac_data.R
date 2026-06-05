# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# reformat_kyle_judy_scatac_data

# Kyle, Judy 
# Kyle provides a list of variants (out of 49283610) that overlaps with each atac-seq cell type

MEM=12500
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -n 2 R

library(data.table)

celltypes<-c("BcellOverlapVariants","EpithelialOverlapVariants",
"FibroblastOverlapVariants","MyeloidOverlapVariants",
"PlasmaCellOverlapVariants","StromalOverlapVariants",
"TcytotoxicOverlapVariants","TlymphoidOverlapVariants")
celltypes<-gsub("OverlapVariants","",celltypes)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

snpids<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed",sep=""),head=F)

for (i in 1:length(celltypes)) {

  tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/levantovsky_2024/ATACoverlapNoThreshold/",celltypes[i],"OverlapVariants.bed",sep=""),head=F)
  tmp$MarkerName<-paste(tmp$V1,tmp$V2,tmp$V4,tmp$V5,sep=":")

  tmp1<-snpids
  tmp1$tmp<-0
  tmp1$tmp[which(tmp$MarkerName %in% snpids$V4)]<-1

  print(nrow(tmp))
  table(tmp1$tmp)

  tmp1<-tmp1[,c(4,5)]
  colnames(tmp1)[2]<-paste(celltypes[i],"sc_atac_perianal_cd",sep="_")
  colnames(tmp1)[1]<-"variant"

  tmp1<-tmp1[match(snpids$V4,tmp1$variant),]
  print(table(tmp1$variant==snpids$V4))

  write.table(tmp1,paste(path_gwas,"/post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/allchr_",celltypes[i],"_sc_atac_perianal_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

  rm(tmp,tmp1)

}

q("no")