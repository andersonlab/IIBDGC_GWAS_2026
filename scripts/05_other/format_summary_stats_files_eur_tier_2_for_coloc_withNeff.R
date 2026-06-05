# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
### HOW TO SUBMIT:
# pheno=(ibd cd uc)
# 
# for ph in ${pheno[@]}
# do
# for chr in {1..22} X
# do
# bsub -J"format" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
# -e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${chr}_format_files_eur_tier_2_for_coloc_stderr \
# -o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${chr}_format_files_eur_tier_2_for_colocC_stdout \
# "/software/R-4.3.1/bin/Rscript ~/git/IIBDGC_GWAS/scripts/other/format_summary_stats_files_eur_tier_2_for_coloc.R ${ph} ${chr} \
# > ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/format_summary_stats_files_eur_tier_2_for_coloc_${ph}_${chr}.Rout"
# done
# done

# MEM=6000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G ibdgwas -q "yesterday" R

library(data.table)

path<-"/path/to/project"
# pheno<-c("ibd","cd","uc")
# pheno<-c("ibd")

args <- commandArgs()

pheno<-args[6]
chr<-args[7]

print(pheno)
print(chr)

for (j in 1:length(pheno)){
  
  tmp<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/",pheno[j],"/",chr,"_",pheno[j],"_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz",sep=""),head=T)
  tmp<-tmp[!duplicated(tmp),]
  
  dup<-tmp$MarkerName[duplicated(tmp$MarkerName)]
  print(paste("Number duplicated ids: ",length(dup),sep=""))
  
  if (length(dup)>0) {
    
    for (i in 1:length(dup)) {
      
      # retain the one with higher info:
      tmp1<-tmp[which(tmp$MarkerName==dup[i]),]
      print(tmp1)
      tmp1<-tmp1$INFO[which(tmp1$INFO==max(tmp1$INFO))]
      
      tmp$MarkerName[which(tmp$MarkerName==dup[i] & tmp$INFO!=tmp1)]<-paste(dup[i],"_rm",sep="")
      
    }

    dup2<-paste(dup,"_rm",sep="")
    tmp<-tmp[which(!tmp$MarkerName %in% dup2),]
    
  }

  tmp$Disease<-pheno[j]
  
  if(pheno[j]=="ibd") {
    maxNeff<-240386/2
  } else if (pheno[j]=="cd") {
    maxNeff<-118922/2
  } else if (pheno[j]=="uc") {
    maxNeff<-145099/2
  }

  print(nrow(tmp))
  tmp<-tmp[which(tmp$Neff>=maxNeff),]
  print(nrow(tmp))
  
  # exclude variants with large HetPval:
  tmp<-tmp[which(tmp$HetPVal>1E-15),]
  print(nrow(tmp))
  
  tmp$MAF<-pmin(tmp$avgA2FREQ,1-tmp$avgA2FREQ)
  tmp$EAF<-tmp$avgA2FREQ
  
  tmp$Pos<-tmp$Position_b38
  tmp$Chr<-chr
  tmp$OR<-NA
  tmp$log_OR<-NA
  tmp$z.score<-NA
  tmp$PubmedID<-NA
  tmp$used_file<-NA
  
  tmp$Eff_allele<-tmp$A2
  
  tmp$RSid<-tmp$MarkerName
  
  colnames(tmp)[7]<-"pval"
  colnames(tmp)[5]<-"beta"
  colnames(tmp)[6]<-"se"
  colnames(tmp)[14]<-"imputation_rsqr"
  
  tmp<-tmp[,c("RSid","Chr","Pos","Eff_allele","MAF","pval","beta","OR","log_OR","se","z.score","Disease","PubmedID","used_file")]
  
  write.table(tmp,paste("/path/to/project",pheno[j],"_chr",chr,"_summary_stats_metal_nogc_Oct2023_eur_tier2_Neff.tsv",sep=""),
              col.names=T,row.names=F,quote=F,sep="\t")
}

q("no")



