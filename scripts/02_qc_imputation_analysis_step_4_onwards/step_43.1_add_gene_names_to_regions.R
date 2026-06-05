# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# add genes in region to the region file:

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=8000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(ggplot2)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

reg<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_regions_420_cojo_supervised_plus_unsupervised.tsv.gz",sep=""))
reg<-as.data.frame(reg)

reg$chr<-as.numeric(gsub("_.*","",reg$updated_region))

reg$min<-gsub("^[0-9]{1,2}_","",reg$updated_region)
reg$max<-as.numeric(gsub("^[0-9]*_","",reg$min))
reg$min<-as.numeric(gsub("_[0-9]*$","",reg$min))

reg$genes_coding<-NA
reg$closest_gene<-NA

gtf<-rtracklayer::import(paste(path_gwas,'post_imputation/2022/analysis/metaanalysis/annotation/gencode.v44.annotation.gtf.gz',sep=""))
gene<-as.data.frame(gtf)
rm(gtf)


for (i in 1:nrow(reg)) {
  
  tmp<-gene[which( (gene$seqnames==paste("chr",reg$chr[i],sep="")) &
                     ((gene$start<=reg$min[i] & gene$end>=reg$min[i]) | 
                        (gene$start<=reg$max[i] & gene$end>=reg$max[i]) | 
                        (gene$start>=reg$min[i] & gene$end<=reg$max[i])) )
            ,c("gene_name","gene_type","start","end","type","strand")]
  tmp<-tmp[!duplicated(tmp$gene_name),]
  tmp<-tmp[which(tmp$gene_type=="protein_coding"),]
  
  reg$genes_coding[i]<-paste(tmp$gene_name[!duplicated(tmp$gene_name)],collapse="|")
    
  rm(tmp)

}

fwrite(reg,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_regions_420_cojo_supervised_plus_unsupervised_with_gene_names.tsv.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

q("no")