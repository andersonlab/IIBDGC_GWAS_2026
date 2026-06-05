# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# HOW TO SUBMIT:
# path_gwas=/path/to/ibdgwas/IIBDGC/
# MEM="1500"
# 
#ancestry=(eur sas eas)
#pheno=(ibd cd uc)
#for j in ${release[@]}
#do
#for ph in ${pheno[@]}
#do
#bsub -J"for_MRMEGA" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
#-e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/format_summary_stats_to_MRMEGA_${ph}_${chr}_${j}_stderr \
#-o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/format_summary_stats_to_MRMEGA_${ph}_${chr}_${j}_stdout \
#"/software/R-4.3.1/bin/Rscript ~/git/IIBDGC_GWAS/scripts/other/format_summary_stats_files_for_MRMEGA.R ${ph} ${j} \
#> ${path_gwas}scripts/logs/format_summary_stats_IIBDGC_danish_ukbb_to_MRMEGA_${ph}_${chr}_${j}.Rout"
#done
#done


# MEM=45000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)
path<-"/path/to/ibdgwas/IIBDGC/"

args <- commandArgs()
pheno<-args[6]
ancestry<-args[7]

# options for ancestry are
# eur
# sas
# eas

if (ancestry=="sas") {
    if(pheno=="ibd") {
        array<-c("ukbb","interval_ibdbioresource")
    } else {
        array<-c("interval_ibdbioresource")
    }
} else if (ancestry=="eas") {
    array<-c("eas_iibdgc")
} else if (ancestry=="eur") {
    if (pheno=="cd") {
        array<-c("illumina370","affymetrix6","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome","ukbb","finngen","decode","danish","interval_ibdbioresource")
    } else {
        array<-c("illumina370","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome","ukbb","finngen","decode","danish","interval_ibdbioresource")
    }
} 

print(ancestry)

# take as input the regenie files, all chr combined, same input as for metal:

for (i in 1:length(array)) {


    if (array[i] %in% c("illumina370","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")) {
        file_tmp<-paste(path,"post_imputation/2022/analysis/regenie/",array[i],"/",pheno,"/allchr_",array[i],"_",ancestry,"_all_step2_",pheno,"_",ancestry,"_sex_PCs_firthse_",pheno,".regenie.gz",sep="")
    } else if (array[i]=="interval_ibdbioresource" & ancestry=="eur") { 
        file_tmp<-paste(path,"post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/",pheno,"/allchr_",array[i],"_",ancestry,"_step2_",pheno,"_",ancestry,"_sex_PCs_firthse_",pheno,".regenie.gz",sep="")
    } else if (array[i]=="interval_ibdbioresource" & ancestry=="sas") { 
        file_tmp<-paste(path,"post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/",pheno,"/allchr_",array[i],"_noneur_step2_",pheno,"_",ancestry,"_sex_PCs_firthse_",pheno,".regenie.gz",sep="")
    } else if (array[i]=="ukbb") {
        file_tmp<-paste(path,"post_imputation/2022/analysis/stage_2_summary_statistics/",array[i],"/",pheno,"/allchr_",array[i],"_",ancestry,"_step2_",pheno,"_",ancestry,"_sex_PCs_firthse_",pheno,".regenie.gz",sep="")
    } else if (array[i]=="finngen") {
        file_tmp<-paste(path,"post_imputation/2022/analysis/stage_2_summary_statistics/",array[i],"/",pheno,"/allchr_",array[i],"_r10_K11_",pheno,".regenie.gz",sep="")
    } else if (array[i]=="danish") {
        file_tmp<-paste(path,"post_imputation/2022/analysis/stage_2_summary_statistics/",array[i],"/",pheno,"/allchr_",array[i],"_gsa_",pheno,"_",ancestry,"_sex_10PCs_saige_spa.regenie.gz",sep="")
    } else if (array[i]=="decode") {
        file_tmp<-paste(path,"post_imputation/2022/analysis/stage_2_summary_statistics/",array[i],"/",pheno,"/allchr_",array[i],"_",pheno,"_30062021_edited_noMult.regenie.gz",sep="")
    } else if (array[i]=="eas_iibdgc") {
        file_tmp<-paste(path,"post_imputation/2022/analysis/stage_2_summary_statistics/",array[i],"/",pheno,"/ibd_EAS_SiKJ_meta_",pheno,".regenie.gz",sep="")
    }

    print(array[i])
    # debug:
    #system(paste("ls -la",file_tmp))
    all<-fread(file_tmp)
    all$CHROM<-gsub("chr","",all$CHROM)

    all$CHROM<-as.character(all$CHROM)
    all$CHROM[which(all$CHROM=="23")]<-"X"

    # subset here the input data - as done for metal:
    all<-all[which(all$INFO>=0.4 & all$A1FREQ>=0.001 & all$A1FREQ<= 0.999),]

    all<-all[,c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","N","BETA","SE")]       
    colnames(all)<-c("CHROMOSOME","POSITION","MarkerName","A1","A2","avgA2FREQ","N","BETA","SE")

    # MR-MEGA format
    all$OR<-exp(all$BETA)
    all$OR_95L<-exp(all$BETA-(1.96*all$SE))
    all$OR_95U<-exp(all$BETA+(1.96*all$SE))

    all<-all[,c("MarkerName","A2","A1","OR","OR_95L","OR_95U","avgA2FREQ","N","CHROMOSOME","POSITION")]
    colnames(all)<-c("MARKERNAME","EA","NEA","OR","OR_95L","OR_95U","EAF","N","CHROMOSOME","POSITION")

    file_out<-paste(path,"post_imputation/2022/analysis/metaanalysis/mrmega/input_files/",pheno,"/allchr_",array[i],"_",ancestry,"_step2_",pheno,"_",ancestry,"_sex_PCs_firthse_regenie_MRMEGA_format.txt.gz",sep="")
    fwrite(all,file_out,col.names=T,row.names=F,sep="\t",quote=F)

    rm(all,file_out)

}

q("no")

