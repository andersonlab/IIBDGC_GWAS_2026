# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# REFORMAT FILES IN REGENIE FORMAT

# MEM=5000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

# # # # HOW TO SUBMIT:
# # # #
# path_gwas=/path/to/ibdgwas/IIBDGC/
# pheno=(cd ibd uc)
# MEM=5000
# 
# i=finngen
#
# for ph in ${pheno[@]}
# do
# for chr in {1..22}
# do
# bsub -J"form" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -m "modern_hardware" \
# -e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/stderr_format_meta_${i}_${ph}_${chr} \
# -o {path_gwas}post_imputation/2022/analysis/metaanalysis/log/stdout_format_meta_${i}_${ph}_${chr} \
# "/software/R-4.3.1/bin/Rscript ${path_ukbb}scripts/format_finngen_files_for_metal.R ${ph} ${chr} > \
# ${path_gwas}scripts/format_finngen_files_for_metal_${ph}_${chr}.Rout"
# done
# done

library(data.table)

args <- commandArgs()
pheno<-args[6]
chr<-args[7]

print(pheno)
print(chr)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

pheno_1<-toupper(pheno)
print(pheno_1)
  
print(chr)

dat<-fread(paste(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/finngen/edited_files/chr",chr,"_finngen_r10_K11_",pheno_1,".gz",sep=""),head=T)
colnames(dat)[1]<-"CHROM"

if(chr==23) {
  dat$CHROM<-as.character(dat$CHROM)
  dat$CHROM<-"X"
}

dat$ID<-paste("chr",dat$CHROM,":",dat$pos,":",dat$ref,":",dat$alt,sep="")

# force INFO to be set to 1 - retain for meta-analysis

dat$INFO<-1
dat$EXTRA<-"NA"
dat$TEST<-"ADD"
dat$CHISQ.Y1<-(dat$beta/dat$sebeta)^2

# updated N

if(pheno=="ibd") {
  dat$N<-sum(9083,392974)
  dat$N_CASES<-9083
  dat$N_CONTROLS<-392974
}
if(pheno=="cd") {
  dat$N<-sum(2033,409940)
  dat$N_CASES<-2033
  dat$N_CONTROLS<-409940
}
if(pheno=="uc") {
  dat$N<-sum(5931,405386)
  dat$N_CASES<-5931
  dat$N_CONTROLS<-405386
}

dat<-dat[,c("CHROM","pos","ID","ref","alt","af_alt","af_alt_cases","af_alt_controls","INFO","N","N_CASES","N_CONTROLS","TEST",
            "beta","sebeta","CHISQ.Y1","mlogp","EXTRA")]
colnames(dat)<-c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","A1FREQ_CASES","A1FREQ_CONTROLS","INFO","N","N_CASES","N_CONTROLS","TEST",
                 "BETA","SE","CHISQ","LOG10P","EXTRA")
if(chr==23) {
  fwrite(dat,paste(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/finngen/",pheno,"/chrX_finngen_r10_K11_",pheno,".regenie",sep=""),
         col.names=T,row.names=F,sep=" ",quote=F)
} else {
  fwrite(dat,paste(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/finngen/",pheno,"/chr",chr,"_finngen_r10_K11_",pheno,".regenie",sep=""),
         col.names=T,row.names=F,sep=" ",quote=F,na="NA")
}

q("no")

