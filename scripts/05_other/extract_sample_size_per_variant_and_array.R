# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
##########################################################################################

# to test:
# # singularity exec iibdgc_postprocess_10_singularity.sif

# MEM=55000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R 

library(data.table)
library(R.utils)

args <- commandArgs()

pheno<-args[6]
chr<-args[7]

print(pheno)
print(chr)

path<-"/path/to/ibdgwas/IIBDGC/"

if(chr=="X" && pheno!="cd") {
  array<-c("affymetrix6","humanomniexpress","gsa")
} else if (chr=="X" && pheno=="cd") {
  array<-c("affymetrix6","gsa")
} else {
  array<-c("illumina370","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")
}

for (j in 1:length(pheno)) {
  
  tmp<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/",pheno[j],"/",chr,"_",pheno[j],"_meta_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),head=T)
  print(nrow(tmp))
    #   [1] 364182

  tmp<-tmp[,c("MarkerName","A1","A2")]

  for (jj in 1:length(array)) {
    
    file_tmp11<-paste(path,"post_imputation/2022/analysis/regenie/",array[jj],"/",pheno[j],"/chr",chr,"_",array[jj],"_eur_all_step2_",pheno[j],"_eur_sex_PCs_firthse_",pheno[j],".regenie.gz",sep="")
    
    if(file.exists(file_tmp11)) {
        
      print(array[jj])
      
      tmp11<-fread(file_tmp11,head=T)
      tmp11<-tmp11[,c("ID","ALLELE0","ALLELE1","A1FREQ","A1FREQ_CASES","A1FREQ_CONTROLS","INFO","N","N_CASES","N_CONTROLS")]

      # set to missing variants MAF<0.001 | INFO<0.4
      tmp11<-tmp11[which(tmp11$INFO>=0.4 & tmp11$A1FREQ>=0.001 & tmp11$A1FREQ<=0.999),]

      tmp11<-merge(tmp,tmp11,by.x="MarkerName",by.y="ID",all.x=T,sort=F)

      # double check alleles match
      print(table(tmp11$ALLELE1==tmp11$A2))
      
      tmp11<-tmp11[,c("MarkerName","A1FREQ","A1FREQ_CASES","A1FREQ_CONTROLS","INFO","N","N_CASES","N_CONTROLS")]
      colnames(tmp11)<-c("MarkerName","A2FREQ","A2FREQ_CASES","A2FREQ_CONTROLS","INFO","N","N_CASES","N_CONTROLS")

      colnames(tmp11)[2:ncol(tmp11)]<-paste(colnames(tmp11)[2:ncol(tmp11)],array[jj],sep="_")
      
      tmp<-merge(tmp,tmp11,by="MarkerName",all.x=T,sort=F)
      }
      
    rm(tmp11,file_tmp11)
    
    }
    
    print(nrow(tmp))
    #   [1] 364182
    
    fwrite(tmp,paste(path,"to_upload_to_iibdgc_globus/chr",chr,"_",pheno[j],"_meta_noGC_PCs_firthse_formatted_A2_effect_A2freq_info_and_N_by_array.tsv.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  
}

q("no")


# ### HOW TO SUBMIT:
# # singularity exec iibdgc_postprocess_10_singularity.sif

# MEM=55000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R 

# pheno=(ibd cd uc)

# # for ph in ${pheno[@]}
# # do
# # for chr in {1..22} X
# # do
# bsub -J"format" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
# -e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_eur_all_${chr}_extract_N_info_freq_IIBDGC_stderr \
# -o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_eur_all_${chr}_extract_N_info_freq_IIBDGC_stdout \
# "Rscript ~/git/IIBDGC_GWAS/scripts/other/extract_sample_size_per_variant_and_array.R ${ph} ${chr} \
# > ${path_gwas}scripts/logs/extract_sample_size_per_variant_and_array_${ph}_eur_all_${chr}.Rout"
# # done
# # done
