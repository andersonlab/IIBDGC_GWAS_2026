# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#


cd /path/to/ibdgwas/IIBDGC/post_imputation/analysis/stage_2_summary_statistics/finngen/
gunzip finngen_R7_K11_*

wc -l finngen_R7_K11_*
# 16383316 finngen_R7_K11_CD_STRICT
# 16383316 finngen_R7_K11_IBD_STRICT
# 16383316 finngen_R7_K11_UC_STRICT
# 49149948 total

gzip finngen_R7_K11_*

#chrom	pos	ref	alt	rsids	nearest_genes	pval	mlogp	beta	sebeta	af_alt	af_alt_cases	af_alt_controls

# split files in chr:
  
path_gwas=/path/to/ibdgwas/IIBDGC/
pheno=(IBD CD UC)

MEM=1000

for ph in ${pheno[@]}
do
for chr in {1..22}
do
bsub -J"fing_split" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${ph}_${chr}_finngen_split_stderr \
-o ${path_gwas}post_imputation/log/${ph}_${chr}_finngen_split_stdout \
"zcat ${path_gwas}post_imputation/analysis/stage_2_summary_statistics/finngen/finngen_R7_K11_${ph}_STRICT.gz \
| grep -w ^${chr} | cat <(zcat ${path_gwas}post_imputation/analysis/stage_2_summary_statistics/finngen/finngen_R7_K11_${ph}_STRICT.gz | head -1) -  \
| gzip > ${path_gwas}post_imputation/analysis/stage_2_summary_statistics/finngen/chr${chr}_finngen_R7_K11_${ph}_STRICT.gz"
done
done

# REFORMAT FILES IN REGENIE AND METAL FORMAT

MEM=8000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)
path<-"/path/to/ibdgwas/IIBDGC/"

pheno_1<-c("CD","IBD","UC")
pheno_2<-c("cd","ibd","uc")

for (i in 1:length(pheno_1)) {
  
  print(pheno_1[i])
  
  for (chr in 1:22) {
    
    print(chr)
    dat<-fread(paste(path,"post_imputation/analysis/stage_2_summary_statistics/finngen/edited_files/chr",chr,"_finngen_R7_K11_",pheno_1[i],"_STRICT.gz",sep=""),head=T)
    colnames(dat)[1]<-"CHROM"
    
    dat$ID<-paste("chr",dat$CHROM,":",dat$pos,":",dat$ref,":",dat$alt,sep="")
    dat$INFO<-"NA"
    dat$TEST<-"ADD"
    dat$CHISQ.Y1<-(dat$beta/dat$sebeta)^2
    dat<-dat[,c("CHROM","pos","ID","ref","alt","af_alt","INFO","TEST","beta","sebeta","CHISQ.Y1","mlogp")]
    colnames(dat)<-c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","TEST","BETA.Y1","SE.Y1","CHISQ.Y1","LOG10P.Y1")
    fwrite(dat,paste(path,"post_imputation/analysis/stage_2_summary_statistics/finngen/",pheno_2[i],"/chr",chr,"_finngen_R7_K11_",pheno_2[i],"_STRICT.regenie",sep=""),
           col.names=T,row.names=F,sep=" ",quote=F)
    
    dat$OR<-exp(dat$BETA.Y1)
    dat$OR_95L<-exp(dat$BETA.Y1-(1.96*dat$SE.Y1))
    dat$OR_95U<-exp(dat$BETA.Y1+(1.96*dat$SE.Y1))

    # note N might not be updated
    if(pheno_2[i]=="ibd") {
      dat$N<-sum(5671,300000)
    }
    if(pheno_2[i]=="cd") {
      dat$N<-sum(1307,300000)
    }
    if(pheno_2[i]=="uc") {
      dat$N<-sum(4024,300000)
    }

    colnames(dat)<-c("CHROMOSOME","POSITION","MARKERNAME","NEA","EA","EAF","INFO","TEST","BETA.Y1","SE.Y1"
                     ,"CHISQ.Y1","LOG10P.Y1","OR","OR_95L","OR_95U","N")

    dat<-dat[,c("MARKERNAME","EA","NEA","OR","OR_95L","OR_95U","EAF","N","CHROMOSOME","POSITION")]
    fwrite(dat,paste(path,"post_imputation/analysis/stage_2_summary_statistics/finngen/",pheno_2[i],"/chr",chr,"_finngen_R7_K11_",pheno_2[i],"_STRICT_MRMEGA_format.txt.gz",sep=""),
           col.names=T,row.names=F,sep="\t",quote=F)

    rm(dat)
  }
}

###########################################################################################################################################
###########################################################################################################################################

# best definitions (case numbers in R10):
# K11_UC_STRICT2  (n=5931)
# K11_CD_STRICT2 (n=2033)
# K11_KELAIBD  (n=9083)


# mv data to specific folders:

cd /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/stage_2_summary_statistics/finngen/
pheno=(ibd cd uc)

for ph in ${pheno[@]}
do mkdir ${ph}
done

# cd finngen_ibd

wc -l Finngen_r10_K11_*


  
# split files in chr:
  
path_gwas=/path/to/ibdgwas/IIBDGC/
pheno=(CD UC)

MEM=1000

for ph in ${pheno[@]}
do
for chr in {1..23}
do
bsub -J"fing_split" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${ph}_${chr}_finngen_split_stderr \
-o ${path_gwas}post_imputation/log/${ph}_${chr}_finngen_split_stdout \
"zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/finngen/raw/Finngen_r10_K11_${ph}_STRICT2.gz \
| grep -w ^${chr} | cat <(zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/finngen/raw/Finngen_r10_K11_${ph}_STRICT2.gz | head -1) -  \
| gzip > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/finngen/edited_files/chr${chr}_finngen_r10_K11_${ph}.gz"
done
done


pheno=(IBD)

for ph in ${pheno[@]}
do
for chr in {1..23}
do
bsub -J"fing_split" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${ph}_${chr}_finngen_split_stderr \
-o ${path_gwas}post_imputation/log/${ph}_${chr}_finngen_split_stdout \
"zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/finngen/raw/Finngen_r10_K11_KELAIBD.gz \
| grep -w ^${chr} | cat <(zcat ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/finngen/raw/Finngen_r10_K11_KELAIBD.gz | head -1) -  \
| gzip > ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/finngen/edited_files/chr${chr}_finngen_r10_K11_${ph}.gz"
done
done

pheno=(IBD CD UC)

for ph in ${pheno[@]}
do echo ${ph} && for chr in {1..23}
do echo ${chr} && tail -50  ${path_gwas}post_imputation/log/${ph}_${chr}_finngen_split_stdout | grep "Successfully"
done
done


# # REFORMAT FILES IN REGENIE (2023) FORMAT

path_gwas=/path/to/ibdgwas/IIBDGC/
pheno=(cd ibd uc)
MEM=5000

i=finngen

for ph in ${pheno[@]}
do
for chr in {1..23}
do
bsub -J"form" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -m "modern_hardware" \
-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/stderr_format_meta_${i}_${ph}_${chr} \
-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/stdout_format_meta_${i}_${ph}_${chr} \
"/software/R-4.3.1/bin/Rscript ${path_gwas}scripts/format_finngen_files_for_metal.R ${ph} ${chr} > \
${path_gwas}scripts/log/format_finngen_files_for_metal_${ph}_${chr}.Rout"
done
done



##############################
# MEM=8000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1
# 
# library(data.table)
# path_gwas<-"/path/to/ibdgwas/IIBDGC/"
# 
# pheno_1<-c("CD","IBD","UC")
# pheno_2<-c("cd","ibd","uc")
# 
# for (i in 1:length(pheno_1)) {
#   
#   print(pheno_1[i])
#   
#   for (chr in c(1:23)) {
#     
#     print(chr)
#     dat<-fread(paste(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/finngen/edited_files/chr",chr,"_finngen_r10_K11_",pheno_1[i],".gz",sep=""),head=T)
#     colnames(dat)[1]<-"CHROM"
#     
#     if(chr==23) {
#       dat$CHROM<-as.character(dat$CHROM)
#       dat$CHROM<-"X"
#     }
#     
#     dat$ID<-paste("chr",dat$CHROM,":",dat$pos,":",dat$ref,":",dat$alt,sep="")
#     dat$INFO<-1
#     dat$EXTRA<-"NA"
#     dat$TEST<-"ADD"
#     dat$CHISQ.Y1<-(dat$beta/dat$sebeta)^2
#     
#     # updated N
#     
#     if(pheno_2[i]=="ibd") {
#       dat$N<-sum(9083,392974)
#       dat$N_CASES<-9083
#       dat$N_CONTROLS<-392974
#     }
#     if(pheno_2[i]=="cd") {
#       dat$N<-sum(2033,409940)
#       dat$N_CASES<-2033
#       dat$N_CONTROLS<-409940
#     }
#     if(pheno_2[i]=="uc") {
#       dat$N<-sum(5931,405386)
#       dat$N_CASES<-5931
#       dat$N_CONTROLS<-405386
#     }
#     
#     dat<-dat[,c("CHROM","pos","ID","ref","alt","af_alt","af_alt_cases","af_alt_controls","INFO","N","N_CASES","N_CONTROLS","TEST",
#                 "beta","sebeta","CHISQ.Y1","mlogp","EXTRA")]
#     colnames(dat)<-c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","A1FREQ_CASES","A1FREQ_CONTROLS","INFO","N","N_CASES","N_CONTROLS","TEST",
#                      "BETA","SE","CHISQ","LOG10P","EXTRA")
#     if(chr==23) {
#       fwrite(dat,paste(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/finngen/",pheno_2[i],"/chrX_finngen_r10_K11_",pheno_2[i],".regenie",sep=""),
#              col.names=T,row.names=F,sep=" ",quote=F)
#     } else {
#       fwrite(dat,paste(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/finngen/",pheno_2[i],"/chr",chr,"_finngen_r10_K11_",pheno_2[i],".regenie",sep=""),
#              col.names=T,row.names=F,sep=" ",quote=F)
#     }
# 
#     
#     
#     rm(dat)
#   }
# }
# 


# K11_UC_STRICT2  (n=5931)
# K11_CD_STRICT2 (n=2033)
# K11_KELAIBD  (n=9083)

# K11_UC_STRICT2 405,386
# K11_CD_STRICT2  409,940
# K11_KELAIBD  392,974

