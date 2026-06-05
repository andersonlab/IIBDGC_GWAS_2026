# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# plot summary_number_variants_per_sample_size

library(data.table)
# library(qqman)
library(colorspace)

path<-"/path/to/ibdgwas/IIBDGC/"

pheno<-c("ibd","cd","uc")

cols<-c("#db7107","#004488","#BB5566")
col_ligth<-lighten(cols, 0.70, space = "HCL")

data<-as.data.frame(matrix(ncol=12,nrow=3))
colnames(data)<-c("pheno","N_variants","N_variants_0.50","N_variants_0.55","N_variants_0.60","N_variants_0.65",
                  "N_variants_0.70","N_variants_0.75","N_variants_0.80","N_variants_0.85","N_variants_0.90","N_variants_0.95")

for (j in 1:length(pheno)) {

  for (chr in c(1:22)) {
    
    tmp<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/",pheno[j],"/chr",chr,"_eur_all_",pheno[j],"_meta_noGC_PCs_firthse_1.txt.gz",sep=""),head=T)
    
    if(chr==1) {
      tmp1<-tmp
    } else {
      tmp1<-rbind(tmp,tmp1)
    }
    rm(tmp)
    
  }
  
  data$pheno[j]<-pheno[j]
  data$N_variants[j]<-nrow(tmp1)
  data$N_variants_0.50[j]<-nrow(tmp1[which(tmp1$TotalSampleSize>=(0.50*max(tmp1$TotalSampleSize))),])
  data$N_variants_0.55[j]<-nrow(tmp1[which(tmp1$TotalSampleSize>=(0.55*max(tmp1$TotalSampleSize))),])
  data$N_variants_0.60[j]<-nrow(tmp1[which(tmp1$TotalSampleSize>=(0.60*max(tmp1$TotalSampleSize))),])
  data$N_variants_0.65[j]<-nrow(tmp1[which(tmp1$TotalSampleSize>=(0.65*max(tmp1$TotalSampleSize))),])
  data$N_variants_0.70[j]<-nrow(tmp1[which(tmp1$TotalSampleSize>=(0.70*max(tmp1$TotalSampleSize))),])
  data$N_variants_0.75[j]<-nrow(tmp1[which(tmp1$TotalSampleSize>=(0.75*max(tmp1$TotalSampleSize))),])
  data$N_variants_0.80[j]<-nrow(tmp1[which(tmp1$TotalSampleSize>=(0.80*max(tmp1$TotalSampleSize))),])
  data$N_variants_0.85[j]<-nrow(tmp1[which(tmp1$TotalSampleSize>=(0.85*max(tmp1$TotalSampleSize))),])
  data$N_variants_0.90[j]<-nrow(tmp1[which(tmp1$TotalSampleSize>=(0.90*max(tmp1$TotalSampleSize))),])
  data$N_variants_0.95[j]<-nrow(tmp1[which(tmp1$TotalSampleSize>=(0.95*max(tmp1$TotalSampleSize))),])
  
  write.table(data,paste(path,"post_imputation/2022/analysis/list_variants_to_include/summary_number_variants_per_percentage_samples_with_data.txt.gz",sep=""),
              col.names=T,row.names=T,quote=F,sep="\t")
  
}

q("no")


# MEM=4000
# pheno=(ibd cd uc)
# for ph in ${pheno[@]}
# do
# bsub -J"plot_imp2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -m "modern_hardware" -G ibdgwas -n 2 \
# -e ${path_gwas}post_imputation/analysis/2022/metaanalysis/log/plot_summary_number_variants_per_sample_size_stderr \
# -o ${path_gwas}post_imputation/analysis/2022/metaanalysis/log/plot_summary_number_variants_per_sample_size_stdout \
# "/software/R-4.3.1/bin/Rscript ${path_gwas}scripts/plot_summary_number_variants_per_sample_size.R ${ph} \
# > ${path_gwas}scripts/logs/plot_summary_number_variants_per_sample_size.Rout"
# done
