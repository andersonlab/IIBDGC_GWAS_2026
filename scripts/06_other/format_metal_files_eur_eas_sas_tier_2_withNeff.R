# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
###########################################################################################

# MEM=25000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)
library(ggpubr)
library(ggplot2)
library(ggExtra)
library(stringr)
library(textclean)

rm(list=ls())

args <- commandArgs()

pheno<-args[6]
chr<-args[7]

print(pheno)
print(chr)

path<-"/path/to/ibdgwas/IIBDGC/"

if(chr=="X" && pheno!="cd") {
  array<-c("affymetrix6","humanomniexpress","gsa","ukbb","finngen","decode","danish","eas","ukbbsas","intervalibdbioresource","intervalibdbioresourcesas")
} else if (chr=="X" && pheno=="cd") {
  array<-c("affymetrix6","gsa","ukbb","finngen","decode","danish","eas","intervalibdbioresource","intervalibdbioresourcesas")
} else if (pheno=="cd") {
  array<-c("illumina370","affymetrix6","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome","ukbb","finngen","decode","danish","eas","intervalibdbioresource","intervalibdbioresourcesas")
} else {
  array<-c("illumina370","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome","ukbb","finngen","decode","danish","eas","ukbbsas","intervalibdbioresource","intervalibdbioresourcesas")
}


for (j in 1:length(pheno)) {
  
  tmp<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/",pheno[j],"/chr",chr,"_eur_eas_sas_",pheno[j],"_meta_noGC_PCs_firthse_1.txt.gz",sep=""),head=F)
  header<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/",pheno[j],"/header",sep=""),head=T)
  colnames(tmp)<-colnames(header)
  rm(header)
    
  if (chr=="X") {
    tmp$tmp<-gsub("chrX:[0-9]*:","",tmp$MarkerName)
  } else {
    tmp$tmp<-gsub("chr[0-9]{1,2}:[0-9]*:","",tmp$MarkerName)
  }

  tmp$A1<-gsub(":[A-Z]*$","",tmp$tmp)
  tmp$A2<-gsub("^[A-Z]*:","",tmp$tmp)
  
  tmp$BETA<-NA
  tmp$BETA[which( tmp$A1==toupper(tmp$Allele1) & tmp$A2==toupper(tmp$Allele2))]<-tmp$Effect[which( tmp$A1==toupper(tmp$Allele1) & tmp$A2==toupper(tmp$Allele2))]
  tmp$BETA[which( tmp$A2==toupper(tmp$Allele1) & tmp$A1==toupper(tmp$Allele2))]<-tmp$Effect[which( tmp$A2==toupper(tmp$Allele1) & tmp$A1==toupper(tmp$Allele2))]*-1
  
  tmp$SE<-tmp$StdErr
  
  tmp$tmp<-swap(tmp$Direction,"+","-")
  
  # update direction of effect accordingly
  tmp$Direction_ed<-NA
  tmp$Direction_ed[which( tmp$A1==toupper(tmp$Allele1) & tmp$A2==toupper(tmp$Allele2))]<-tmp$Direction[which( tmp$A1==toupper(tmp$Allele1) & tmp$A2==toupper(tmp$Allele2))]
  tmp$Direction_ed[which( tmp$A2==toupper(tmp$Allele1) & tmp$A1==toupper(tmp$Allele2))]<-tmp$tmp[which( tmp$A2==toupper(tmp$Allele1) & tmp$A1==toupper(tmp$Allele2))]
  

  # if chrX keep pseudoautosomal region only:
  
  if (chr=="X") {
    tmp$Position_b38<-gsub("chrX:","",tmp$MarkerName)
    tmp$Position_b38<-as.numeric(gsub(":[A-Z]*:[A-Z]*$","",tmp$Position_b38))
    tmp<-as.data.frame(tmp)
    tmp<-tmp[which(tmp$Position_b38<=2781479 | tmp$Position_b38>=155701383),]
  } else {
    tmp$Position_b38<-gsub("chr[0-9]{1,2}:","",tmp$MarkerName)
    tmp$Position_b38<-as.numeric(gsub(":[A-Z]*:[A-Z]*$","",tmp$Position_b38))
  }
  
  tmp$Position_b38<-format(tmp$Position_b38, scientific=FALSE)
  
  tmp<-tmp[,c("MarkerName","Position_b38","A1","A2","BETA","SE","P-value","Direction_ed","HetISq","HetChiSq","HetDf","HetPVal","TotalSampleSize")]
  
  
  # ADD AVG ALT (AND EFFECT) ALLELE FREQUENCY, PLUS INFO
  rm(tmp1)
  
  for (jj in 1:length(array)) {
    

    if (array[jj]=="ukbb") {
      file_tmp11<-paste(path,"post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/",pheno[j],"/chr",chr,"_",array[jj],"_eur_step2_",pheno[j],"_eur_sex_PCs_firthse_",pheno[j],".regenie.gz",sep="")
    } else if (array[jj]=="finngen") {
      file_tmp11<-paste(path,"post_imputation/2022/analysis/stage_2_summary_statistics/",array[jj],"/",pheno[j],"/chr",chr,"_",array[jj],"_r10_K11_",pheno[j],".regenie.gz",sep="")
    } else if (array[jj]=="decode") {
      file_tmp11<-paste(path,"post_imputation/2022/analysis/stage_2_summary_statistics/",array[jj],"/",pheno[j],"/allchr_",array[jj],"_",pheno[j],"_30062021_edited_noMult.regenie.gz",sep="")
    } else if (array[jj]=="danish") {
      file_tmp11<-paste(path,"post_imputation/2022/analysis/stage_2_summary_statistics/",array[jj],"/",pheno[j],"/allchr_",array[jj],"_gsa_",pheno[j],"_eur_sex_10PCs_saige_spa.regenie.gz",sep="")
    } else if (array[jj]=="danish") {
      file_tmp11<-paste(path,"post_imputation/2022/analysis/stage_2_summary_statistics/",array[jj],"/",pheno[j],"/allchr_",array[jj],"_gsa_",pheno[j],"_eur_sex_10PCs_saige_spa.regenie.gz",sep="")
    } else if (array[jj]=="ukbbsas") {
      file_tmp11<-paste(path,"post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/",pheno[j],"/chr",chr,"_ukbb_sas_step2_",pheno[j],"_sas_sex_PCs_firthse_",pheno[j],".regenie.gz",sep="")
    } else if (array[jj]=="intervalibdbioresource") {
      file_tmp11<-paste(path,"post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/",pheno[j],"/",chr,"_interval_ibdbioresource_eur_step2_",pheno[j],"_eur_sex_PCs_firthse_",pheno[j],".regenie.gz",sep="")
    } else if (array[jj]=="intervalibdbioresourcesas") {
      file_tmp11<-paste(path,"post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/",pheno[j],"/",chr,"_interval_ibdbioresource_noneur_step2_",pheno[j],"_sas_sex_PCs_firthse_",pheno[j],".regenie.gz",sep="")
    } else if (array[jj]=="eas") {
      file_tmp11<-paste(path,"post_imputation/2022/analysis/stage_2_summary_statistics/eas_iibdgc/",pheno[j],"/ibd_EAS_SiKJ_meta_",pheno[j],".regenie.gz",sep="")
    } else {
      file_tmp11<-paste(path,"post_imputation/2022/analysis/regenie/",array[jj],"/",pheno[j],"/chr",chr,"_",array[jj],"_eur_all_step2_",pheno[j],"_eur_sex_PCs_firthse_",pheno[j],".regenie.gz",sep="")
    }

    if(file.exists(file_tmp11)) {
      
      tmp11<-fread(file_tmp11,head=T)
      tmp11$CHROM<-gsub("chr","",tmp11$CHROM)
      tmp11<-tmp11[which(tmp11$CHROM==chr),]
      
      tmp11<-tmp11[,c("ID","A1FREQ","A1FREQ_CASES","A1FREQ_CONTROLS","INFO","N","N_CASES","N_CONTROLS")]
      tmp11$Neff<-round(4/((1/as.numeric(tmp11$N_CASES))+(1/as.numeric(tmp11$N_CONTROLS))))
      
      colnames(tmp11)[2:ncol(tmp11)]<-paste(colnames(tmp11)[2:ncol(tmp11)],array[jj],sep="_")
      
      if(!exists("tmp1")) {
        tmp1<-tmp11
      } else {
        tmp1<-merge(tmp1,tmp11,by="ID",all=T)
      }
      
      rm(tmp11)
      
    } else {
      print(array[jj])
    }
    
  }
  
  
  # note that regenie files use A0/A1 notation, and metal output files use A1/A2
  
  tmp1<-as.data.frame(tmp1)
  
  ### OVERALL:
  toMatch<-c("N_[a-z0-9]*$","N_[a-z]_*$")
  N<-tmp1[,colnames(tmp1)[grep(paste(toMatch,collapse="|"),colnames(tmp1))],drop=F]

  toMatch<-c("Neff_[a-z0-9]*$","Neff_[a-z]_*$")
  Neff<-tmp1[,colnames(tmp1)[grep(paste(toMatch,collapse="|"),colnames(tmp1))],drop=F]
  
  A1FREQ<-tmp1[,colnames(tmp1)[grep("A1FREQ_[a-z0-9]*$",colnames(tmp1))],drop=F]
  INFO<-tmp1[,colnames(tmp1)[grep("INFO_[a-z0-9]*$",colnames(tmp1))],drop=F]
  
  A1FREQ[INFO<0.4]<-NA
  INFO[A1FREQ<0.001 | A1FREQ>0.999]<-NA
  
  A1FREQ[A1FREQ<0.001 | A1FREQ>0.999]<-NA
  INFO[INFO<0.4]<-NA
  
  A1FREQ[is.na(INFO)]<-NA
  INFO[is.na(A1FREQ)]<-NA
  
  N[is.na(A1FREQ) | is.na(INFO)]<-NA
  Neff[is.na(A1FREQ) | is.na(INFO)]<-NA
  
  tmp1$N<-NA
  tmp1$N<-rowSums(N,na.rm=T)

  tmp1$Neff<-NA
  tmp1$Neff<-rowSums(Neff,na.rm=T)
  
  A1FREQ2<-A1FREQ*N
  tmp1$avgA2FREQ<-NA
  tmp1$avgA2FREQ<-rowSums(A1FREQ2,na.rm=T)/rowSums(N,na.rm=T)
  
  INFO2<-INFO*N
  tmp1$INFO<-NA
  tmp1$INFO<-rowSums(INFO2,na.rm=T)/rowSums(N,na.rm=T)
  
  
  ### IN CASES:
  
  N_CASES<-tmp1[,colnames(tmp1)[grep("N_CASES_",colnames(tmp1))],drop=F]
  N_CASES[is.na(A1FREQ) | is.na(INFO)]<-NA
  
  A1FREQ_CASES<-tmp1[,colnames(tmp1)[grep("A1FREQ_CASES",colnames(tmp1))],drop=F]
  A1FREQ_CASES2<-A1FREQ_CASES*N_CASES
  A1FREQ_CASES2[is.na(A1FREQ) | is.na(INFO)]<-NA
  
  
  tmp1$avgA2FREQ_CASES<-NA
  tmp1$avgA2FREQ_CASES<-rowSums(A1FREQ_CASES2,na.rm=T)/rowSums(N_CASES,na.rm=T)
  
  tmp1$N_CASES<-NA
  tmp1$N_CASES<-rowSums(N_CASES,na.rm=T)
  
  
  ### IN CONTROLS:
  
  N_CONTROLS<-tmp1[,colnames(tmp1)[grep("N_CONTROLS_",colnames(tmp1))],drop=F]
  N_CONTROLS[is.na(A1FREQ) | is.na(INFO)]<-NA
  
  A1FREQ_CONTROLS<-tmp1[,colnames(tmp1)[grep("A1FREQ_CONTROLS",colnames(tmp1))],drop=F]
  A1FREQ_CONTROLS2<-A1FREQ_CONTROLS*N_CONTROLS
  A1FREQ_CONTROLS2[is.na(A1FREQ) | is.na(INFO)]<-NA
  
  
  tmp1$avgA2FREQ_CONTROLS<-NA
  tmp1$avgA2FREQ_CONTROLS<-rowSums(A1FREQ_CONTROLS2,na.rm=T)/rowSums(N_CONTROLS,na.rm=T)
  
  tmp1$N_CONTROLS<-NA
  tmp1$N_CONTROLS<-rowSums(N_CONTROLS,na.rm=T)
  
  # combine with metal data
  
  tmp<-merge(tmp,tmp1[,c("ID","INFO","N","avgA2FREQ","N_CASES","avgA2FREQ_CASES","N_CONTROLS","avgA2FREQ_CONTROLS","Neff")],by.x="MarkerName",by.y="ID",sort=F,all.x=T)
  
  nmax<-max(tmp$N,na.rm=T)
  tmp$rate_total_sample<-tmp$N/nmax
  
  Neffmax<-max(tmp$Neff,na.rm=T)
  tmp$rate_Neff<-tmp$Neff/Neffmax
  
  # OK - recalculated N matches metal output
  print(table(tmp$TotalSampleSize-tmp$N))
  
  fwrite(tmp,paste(path,"post_imputation/2022/analysis/metaanalysis/",pheno[j],"/",chr,"_",pheno[j],"_meta_eur_eas_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz",sep=""),
         col.names=T,row.names=F,quote=F,sep="\t")
  
  rm(tmp,sample)
  
}

q("no")


### HOW TO SUBMIT:
# pheno=(ibd cd uc)
# 
# for j in ${ancestry[@]}
# do
# for ph in ${pheno[@]}
# do
# for chr in {1..22} X
# do
# sub -J"format" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
# -e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${j}_${chr}_format_metal_IIBDGC_stderr \
# -o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_${j}_${chr}_format_metal_IIBDGC_stdout \
# "/software/R-4.3.1/bin/Rscript ${path_gwas}scripts/format_metal_files_eur_tier.R ${ph} ${i} \
# > ${path_gwas}scripts/logs/format_metal_files_eur_tier_${ph}_${j}_${chr}.Rout"
# done
# done
