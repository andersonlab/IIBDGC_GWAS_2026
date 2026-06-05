# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# # how to submit
# # singularity exec iibdgc_postprocess_10_singularity.sif

# MEM=6000
# bsub -J"plot_imp1" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -m "modern_hardware" -G ibdgwas \
# -e ${path_gwas}post_imputation/analysis/metaanalysis/log/draw_manhattan_plots_stderr \
# -o ${path_gwas}post_imputation/analysis/metaanalysis/log/draw_manhattan_plots_stdout \
# "Rscript ${path_gwas}scripts/draw_manhattan_plots.R > \
# ${path_gwas}scripts/logs/draw_manhattan_plots.Rout"

# for testing purposes
# MEM=6000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(qqman)
library(colorspace)

path<-"/path/to/ibdgwas/IIBDGC/"
path_2<-"~/git/IIBDGC_GWAS/"

pheno<-c("ibd","cd","uc")
cols<-c("#db7107","#004488","#BB5566")
col_ligth<-lighten(cols, 0.70, space = "HCL")

args <- commandArgs()
release<-args[6]

# options for release are: 
#   eur_tier_1
#   eur_tier_2
#   sas_tier_2
#   eur_eas_sas_tier_2

for (i in 1:length(pheno)){
  
  print(pheno[i])
  
  for (chr in c(1:22,"X")) {
    
    if (release=="eur_tier_1") {
        tmp<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),)
        tmp<-tmp[which(tmp$N>=(0.8*max(tmp$N,na.rm=T))),]
    } else {
        tmp<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_",release,"_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),)
        tmp<-tmp[which(tmp$N>=(0.4*max(tmp$N,na.rm=T))),]
    } 

    tmp$'P-value'<-as.numeric(tmp$'P-value')
    tmp<-tmp[which(tmp$'P-value'<0.1),]

    tmp$HetPVal<-as.numeric(tmp$HetPVal)

    if(chr=="X") {
      tmp$pos<-gsub("chrX:","",tmp$MarkerName)
      tmp$pos<-as.numeric(gsub(":[A-Z]*:[A-Z]*","",tmp$pos))
      # keep only PAR:
      tmp<-tmp[which(tmp$pos<=2781479 | tmp$pos>=155701383),] 
    }else{
       tmp$pos<-gsub("chr[0-9]{1,2}:","",tmp$MarkerName)
       tmp$pos<-as.numeric(gsub(":[A-Z]*:[A-Z]*","",tmp$pos))
    }
      
    tmp$chr<-chr

    if(chr=="X") {
        tmp$chr<-23
    }

    # plot only signals with no significant het in effect size:
    if (release=="sas_tier_2") {
        tmp<-tmp[which(tmp$HetPVal>=1E-3),]
    } else {
        tmp<-tmp[which(tmp$HetPVal>=1E-10),]
    }

    tmp<-tmp[,c("chr","MarkerName","pos","P-value")]
    
    if(chr==1) {
      tmp1<-tmp
    } else {
      tmp1<-rbind(tmp,tmp1)
    }
    rm(tmp)
    
  }

  tmp1$chr<-as.numeric(tmp1$chr)


  if (release=="sas_tier_2") {
    ylim_max<-1E-30
    pdf(paste(path_2,"plots/metaanalysis/",pheno[i],"_",release,"_meta_noGC_PCs_firthse_ylim30.pdf",sep=""),width=14, height=5)
    manhattan(tmp1, genomewideline = FALSE, suggestiveline = FALSE, main = paste(toupper(pheno[i]),release),
    chr = "chr", bp = "pos", p = "P-value", snp = "MarkerName",ylim=c(0,-log10(ylim_max)),col = c(cols[i],col_ligth[i],"gray60"))
    dev.off()

  } else {
      ylim_max<-1E-300
      pdf(paste(path_2,"plots/metaanalysis/",pheno[i],"_",release,"_meta_noGC_PCs_firthse.pdf",sep=""),width=14, height=5)
      manhattan(tmp1, genomewideline = FALSE, suggestiveline = FALSE, main = paste(toupper(pheno[i]),release),
      chr = "chr", bp = "pos", p = "P-value", snp = "MarkerName",ylim=c(0,-log10(ylim_max)),col = c(cols[i],col_ligth[i],"gray60"))
      dev.off()

     ylim_max<-1E-30
     pdf(paste(path_2,"plots/metaanalysis/",pheno[i],"_",release,"_meta_noGC_PCs_firthse_ylim30.pdf",sep=""),width=14, height=5)
     manhattan(tmp1, genomewideline = FALSE, suggestiveline = FALSE, main = paste(toupper(pheno[i]),release),
     chr = "chr", bp = "pos", p = "P-value", snp = "MarkerName",ylim=c(0,-log10(ylim_max)),col = c(cols[i],col_ligth[i],"gray60"))
     dev.off()
     
  }
  
  rm(tmp1)

}

q("no")









