# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752

# MEM=15000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q cpu-interactive R

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/sex_specific_analyses/"
path_gwas<-"/path/to/ibdgwas/IIBDGC/"


# provide pheno
args = commandArgs(trailingOnly=TRUE)

ph<-args[1]

# for troubleshooting purposes
# ph<-"ibd"

for (chr in c(1:22)) {

    tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr",chr,"_subset_included_in_ibd_analysis.bim"))
    if (chr==1) {
        bim<-tmp
    } else {
        bim<-rbind(tmp,bim)
    }
    rm(tmp)
}

for (chr in c(1:22)) {

    print(chr)

    tmp<-fread(paste0(path,ph,"_sexhet_cojo/chr",chr,"_sex_",toupper(ph),"_gwama_metal_withHCE_linked_out_InfoMafFilt.txt"))
    tmp$SNP<-tmp$rs_number
    tmp$SNP2<-gsub(":[A-Z]*:[A-Z]*$","",tmp$SNP)
    tmp$SNP2<-paste0(tmp$SNP2,":",tmp$reference_allele,":",tmp$other_allele)

    # gwama effect allele
    tmp$A1<-tmp$reference_allele
    tmp$A2<-tmp$other_allele

    tmp$beta_gwama_male_sex_specific<-log(as.numeric(tmp$male_OR))
    tmp$se_gwama_male_sex_specific<-tmp$beta_gwama_male_sex_specific/as.numeric(tmp$male_z)
    tmp$N_gwama_male_sex_specific<-as.numeric(tmp$male_n_samples)
    tmp$pvalue_gwama_male_sex_specific<-tmp$'male_p-value'


    tmp$beta_gwama_female_sex_specific<-log(as.numeric(tmp$female_OR))
    tmp$se_gwama_female_sex_specific<-tmp$beta_gwama_female_sex_specific/as.numeric(tmp$female_z)
    tmp$N_gwama_female_sex_specific<-as.numeric(tmp$female_n_samples)
    tmp$pvalue_gwama_female_sex_specific<-tmp$'female_p-value'

    tmp$beta_gwama_male_female<-log(as.numeric(tmp$OR))
    tmp$se_gwama_male_female<-tmp$beta_gwama_male_female/as.numeric(tmp$z)
    tmp$N_gwama_male_female<-as.numeric(tmp$n_samples)
    tmp$pvalue_gwama_male_female<-tmp$'p-value'


    tmp1<-tmp[which(tmp$SNP %in% bim$V2),]
    tmp2<-tmp[which(tmp$SNP2 %in% bim$V2),]

    print(nrow(tmp1))
    # print(nrow(tmp2))

    # if (nrow(tmp2)>0) {

    #     tmp2$beta_gwama_male_sex_specific<-(tmp2$beta_gwama_male_sex_specific)*-1
    #     tmp2$beta_gwama_female_sex_specific<-(tmp2$beta_gwama_male_sex_specific)*-1
    #     tmp2$beta_gwama_male_female<-(tmp2$beta_gwama_male_female)*-1

    #     tmp2$A2<-tmp2$reference_allele
    #     tmp2$A1<-tmp2$other_allele

    #     tmp2$male_eaf<-1-tmp2$male_eaf
    #     tmp2$female_eaf<-1-tmp2$female_eaf
    #     tmp2$eaf<-1-tmp2$eaf

    # }

    # print(nrow(tmp))
    # tmp<-rbind(tmp1,tmp2)
    print(nrow(tmp))


    tmp<-tmp[,c("SNP","A1","A2","male_eaf","beta_gwama_male_sex_specific","se_gwama_male_sex_specific","N_gwama_male_sex_specific","pvalue_gwama_male_sex_specific",
        "beta_gwama_female_sex_specific","female_eaf","se_gwama_female_sex_specific","N_gwama_female_sex_specific","pvalue_gwama_female_sex_specific",
        "beta_gwama_male_female","eaf","se_gwama_male_female","N_gwama_male_female","pvalue_gwama_male_female")]

    if (chr==1) {
        dat<-tmp
    } else {
        dat<-rbind(dat,tmp)
    }
    rm(tmp)

}


### create 3 sets of files, one per trait (female specific, male specifc, )

male<-dat[,c("SNP","A1","A2","male_eaf","beta_gwama_male_sex_specific","se_gwama_male_sex_specific","pvalue_gwama_male_sex_specific","N_gwama_male_sex_specific")]
colnames(male)<-c("SNP","A1","A2","freq","beta","se","p","N")
dim(male)
male<-male[which(male$N>=(0.5*max(male$N)) & male$freq>=0.001 & male$freq<=0.999),]
dim(male)
fwrite(male,paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/allchr_",ph,"_meta_eur_tier_1_male_sex_specific_gwama_cojo_format.ma"),col.names=T,row.names=F,quote=F,sep="\t")


female<-dat[,c("SNP","A1","A2","female_eaf","beta_gwama_female_sex_specific","se_gwama_female_sex_specific","pvalue_gwama_female_sex_specific","N_gwama_female_sex_specific")]
colnames(female)<-c("SNP","A1","A2","freq","beta","se","p","N")
dim(female)
female<-female[which(female$N>=(0.5*max(female$N)) & female$freq>=0.001 & female$freq<=0.999),]
dim(female)
fwrite(female,paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/allchr_",ph,"_meta_eur_tier_1_female_sex_specific_gwama_cojo_format.ma"),col.names=T,row.names=F,quote=F,sep="\t")


all<-dat[,c("SNP","A1","A2","eaf","beta_gwama_male_female","se_gwama_male_female","pvalue_gwama_male_female","N_gwama_male_female")]
colnames(all)<-c("SNP","A1","A2","freq","beta","se","p","N")
dim(all)
all<-all[which(all$N>=(0.5*max(all$N)) & all$freq>=0.001 & all$freq<=0.999),]
dim(all)
fwrite(all,paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/allchr_",ph,"_meta_eur_tier_1_female_male_sex_specific_gwama_cojo_format.ma"),col.names=T,row.names=F,quote=F,sep="\t")


q("no")

