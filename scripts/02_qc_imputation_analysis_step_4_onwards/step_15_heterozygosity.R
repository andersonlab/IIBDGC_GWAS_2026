# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
#######################
# 15 - HETEROZYGOSITY #
#######################

# see script pca_iibdgc_1000GP_iibdgc_projection.R

##############################################################################
# 15.1 Use list non duplicated samples to estimate heterozigosity limits:


studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

### define in noDuplicates
for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_48_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_48_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_hwe \
--allow-no-sex \
--chr 1-22 --het \
--out ${path_gwas}pre_imputation/QC/${i}/autosomal_het_nodup"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_48_${i} | grep -E "completed"
done



### apply to all
for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_49_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_49_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe \
--exclude ${path_gwas}/pre_imputation/QC/${i}/list_var_exclude_hwe \
--allow-no-sex \
--chr 1-22 --het \
--out ${path_gwas}pre_imputation/QC/${i}/autosomal_het_all"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_49_${i} | grep -E "completed"
done


## /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

length(cohorts)
# [1] 58


pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/pca_1000gp_iibdgc_pca_pcaEur_withDuplicates.tsv",sep=""),head=T,sep="\t")


for (j in 1:length(cohorts)) {
  
  print(cohorts[j])
  
  # 1.- set thresholds from non duplicated set:
  
  het<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/autosomal_het_nodup.het",sep=""),head=T)
  dim(het)
  
  het$het<-(het$N.NM.-het$O.HOM)/het$N.NM.
  
  het<-merge(het[,c("FID","O.HOM.","E.HOM.","N.NM.","F","het")],pca[,c("FID","cohort","inferred_population","self_jewish","pca_jewish","PC1_EUR")],
             by="FID",all.x=T)
  
  table(het$inferred_population)
 
  het$inferred_population_2<-het$inferred_population
  het$inferred_population_2<-as.character(het$inferred_population_2)
  het$inferred_population_2[which(het$pca_jewish=="Jewish")]<-"EUR_Jewish"
  het$inferred_population_2[which(het$pca_jewish=="Non-Jewish")]<-"EUR_Non-Jewish"
  
  table(het$inferred_population_2)
  
  het$inferred_population_2<-as.factor(het$inferred_population_2)
  
  #### CREATE AN HISTOGRAM
  
  # number of bins in histogram
  fd=function(x) {
    n=length(x)
    r=IQR(x)
    2*r/n^(1/3)
  }
  
  limits<-matrix(ncol=3,nrow=length(levels(het$inferred_population_2)) )
  limits<-as.data.frame(limits)
  colnames(limits)<-c("pop","lower","upper")
  
  for (i in c(1:length(levels(het$inferred_population_2))) ) {
    
    limits$pop[i]<-levels(het$inferred_population_2)[i]
    
    sds<-(4*sd(het$het[which(het$inferred_population_2==levels(het$inferred_population_2)[i])]))
    if(is.na(sds)) {
      sds<-0.01
    }
    limits$lower[i]<-mean(het$het[which(het$inferred_population_2==levels(het$inferred_population_2)[i])])-sds
    limits$upper[i]<-mean(het$het[which(het$inferred_population_2==levels(het$inferred_population_2)[i])])+sds
    
  }
  
  ## 2.- apply to the rest of the samples
  
  het<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/autosomal_het_all.het",sep=""),head=T)
  dim(het)
  
  het$het<-(het$N.NM.-het$O.HOM)/het$N.NM.
  
  het<-merge(het[,c("FID","O.HOM.","E.HOM.","N.NM.","F","het")],pca[,c("FID","cohort","inferred_population","self_jewish","pca_jewish","PC1_EUR")],
             by="FID",all.x=T)
  
  table(het$inferred_population,useNA="ifany")
  
  het$inferred_population_2<-het$inferred_population
  het$inferred_population_2<-as.character(het$inferred_population_2)
  het$inferred_population_2[which(het$pca_jewish=="Jewish")]<-"EUR_Jewish"
  het$inferred_population_2[which(het$pca_jewish=="Non-Jewish")]<-"EUR_Non-Jewish"
  
  table(het$inferred_population_2,useNA="ifany")
  het$inferred_population_2<-as.factor(het$inferred_population_2)
  
  #### CREATE AN HISTOGRAM
  
  # number of bins in histogram
  fd=function(x) {
    n=length(x)
    r=IQR(x)
    2*r/n^(1/3)
  }
  
  
  # to emulate ggplot2 palette:
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  n = length(levels(het$inferred_population_2))
  cols = gg_color_hue(n)
  
  xmax<-max(max(het$het),max(limits$upper))+0.01
  xmin<-min(min(het$het),min(limits$lower))-0.01
  
  for (i in c(1:length(levels(het$inferred_population_2))) ) {
    print(i)
    fd_bin<-fd(het$het[which(het$inferred_population_2==levels(het$inferred_population_2)[i])])
    assign(paste("p",i,sep=""),ggplot(het[which(het$inferred_population_2==levels(het$inferred_population_2)[i]),],aes(x=het,color=inferred_population_2,fill=inferred_population_2)) +
             geom_histogram(position="identity", alpha=0.5, binwidth=fd_bin) +
             scale_fill_manual(values=cols[i]) + scale_color_manual(values=cols[i]) + 
             geom_vline(xintercept = limits$lower[i],linetype="dotted", color = "black") +  
             geom_vline(xintercept = limits$upper[i],linetype="dotted", color = "black") + ggtitle(levels(het$inferred_population_2)[i]) +
             xlim(xmin,xmax) )
  }
  
  if(exists("p6")) {
    pdf(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_het_histogram.pdf",sep=""),width=18,height=(length(levels(het$inferred_population_2)) ) )
    print(ggarrange(p1,p2,p3,p4,p5,p6,nrow=ceiling(length(levels(het$inferred_population_2))/2),ncol=2 ))
    dev.off()
    system(paste("cp ",path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_het_histogram.pdf ~/tmp_plots/",sep=""))
  }else if(exists("p5")){
    pdf(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_het_histogram.pdf",sep=""),width=18,height=(length(levels(het$inferred_population_2)) ) )
    print(ggarrange(p1,p2,p3,p4,p5,nrow=ceiling(length(levels(het$inferred_population_2))/2),ncol=2 ))
    dev.off()
    system(paste("cp ",path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_het_histogram.pdf ~/tmp_plots/",sep=""))
  }else if(exists("p4")){
    pdf(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_het_histogram.pdf",sep=""),width=18,height=(length(levels(het$inferred_population_2)) ) )
    print(ggarrange(p1,p2,p3,p4,nrow=ceiling(length(levels(het$inferred_population_2))/2),ncol=2 ))
    dev.off()
    system(paste("cp ",path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_het_histogram.pdf ~/tmp_plots/",sep=""))
  }else if(exists("p3")){
    pdf(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_het_histogram.pdf",sep=""),width=18,height=(length(levels(het$inferred_population_2)) ) )
    print(ggarrange(p1,p2,p3,nrow=ceiling(length(levels(het$inferred_population_2))/2),ncol=2 ))
    dev.off()
    system(paste("cp ",path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_het_histogram.pdf ~/tmp_plots/",sep=""))
  }else if(exists("p2")){
    pdf(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_het_histogram.pdf",sep=""),width=18,height=(length(levels(het$inferred_population_2)) ) )
    print(ggarrange(p1,p2,nrow=ceiling(length(levels(het$inferred_population_2))/2),ncol=2 ))
    dev.off()
    system(paste("cp ",path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_het_histogram.pdf ~/tmp_plots/",sep=""))
  }else if(exists("p1")){
    pdf(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_het_histogram.pdf",sep=""),width=9,height=(length(levels(het$inferred_population_2))+1) )
    print(p1)
    dev.off()
    system(paste("cp ",path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_het_histogram.pdf ~/tmp_plots/",sep=""))
  }
  
  
  
  for (i in c(1:nrow(limits)) ) {
    if(i==1) {
      remove<-het[which( (het$inferred_population_2==limits$pop[i]) &
                           ( (het$het<limits$lower[i]) |
                               (het$het>limits$upper[i]) ) ), ]
    } else {
      tmp<-het[which( (het$inferred_population_2==limits$pop[i]) &
                        ( (het$het<limits$lower[i]) |
                            (het$het>limits$upper[i]) ) ), ]
      remove<-rbind(remove,tmp)
    }
    
  }
  
  print(table(remove$inferred_population_2))
  print(paste("N samples remove - with dups: ",nrow(remove),sep=""))
  write.table(remove[,c("FID","FID")],paste(path,"pre_imputation/QC/",cohorts[j],"/list_samples_het_outliers",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  
  # no duplicated:
  fam<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_hwe.fam",sep="")
                  ,head=F)
  print(paste("N samples remove - no dups: ",nrow(remove[which(remove$FID %in% fam$V1),]),sep=""))
  
  rm(list=ls()[!ls() %in% c("path","j","cohorts","pca")])
}

q("no")


######################


### SAMPLES WITHOUT DUPLICATES

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_50_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_50_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_hwe \
--allow-no-sex \
--remove ${path_gwas}/pre_imputation/QC/${i}/list_samples_het_outliers \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_hwe_het"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_50_${i} | grep -E "completed"
done


### SAMPLES WITH DUPLICATES

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_51_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_51_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe \
--allow-no-sex \
--remove ${path_gwas}/pre_imputation/QC/${i}/list_samples_het_outliers \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_51_${i} | grep -E "completed"
done


for i in ${studies[@]}
do echo ${i} && wc -l  ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam \
&& wc -l  ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
done

#########
# all_hce
# 23763 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/all_hce/all_hce_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 434598 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/all_hce/all_hce_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# iddk_old_gwas
# 2732 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 300166 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# australia_omniexome
# 1299 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/australia_omniexome_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 739312 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/australia_omniexome_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# gwas1
# 4659 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/gwas1_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 453416 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/gwas1_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# gwas2
# 7773 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas2/gwas2_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 773228 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas2/gwas2_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# pittsburgh_gsa
# 2709 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 897842 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# spain_gsa
# 3408 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/spain_gsa/spain_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 575804 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/spain_gsa/spain_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# italy_gsa
# 943 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/italy_gsa/italy_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 573850 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/italy_gsa/italy_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# kiel_austria_sibdcs_gsa
# 14329 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 564511 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# netherlands_gsa
# 4537 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 666191 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# slovenia_gsa
# 263 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 529608 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# sweden_gsa
# 1396 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/sweden_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 558607 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/sweden_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# niddk_broad_gsa
# 5426 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 590626 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# niddk_feinstein_gsa
# 8124 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 603214 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# basque_gsa
# 1498 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/basque_gsa/basque_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 567401 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/basque_gsa/basque_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# lithuania_gsa
# 2225 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lithuania_gsa/lithuania_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 554135 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lithuania_gsa/lithuania_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# belgium_louis_gsa
# 1511 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 591948 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# belgium_franchimont_gsa
# 1491 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 561994 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# belgium_vermeire_gsa
# 3974 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 653147 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# prism_nfe_gsa
# 462 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 559017 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# prism_nfe_gwas
# 810 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 242899 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# finland_illugwas
# 444 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/finland_illugwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 237055 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/finland_illugwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# german_affy6_old_gwas
# 2797 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 783439 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# norway_affy6_old_gwas
# 544 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 725899 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# belgium_inf1_old_gwas
# 1406 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 301714 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# belgium_inf2_old_gwas
# 271 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 290191 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# cedars_370k_old_gwas
# 597 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 338717 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# cedars_610k_old_gwas
# 881 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 583907 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# cedars_omni_old_gwas
# 1207 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 719139 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# swedish_uc_old_gwas
# 1258 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 296648 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# mccauley_gsa
# 776 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 556682 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# ccfa_gsa
# 2167 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 599981 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# cedars_gsa
# 3059 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_gsa/cedars_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 646302 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_gsa/cedars_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# bernstein_gsa
# 509 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/bernstein_gsa/bernstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 553944 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/bernstein_gsa/bernstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# farkkila_gsa
# 68 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 511512 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# franchimont_gsa
# 2769 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franchimont_gsa/franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 648211 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franchimont_gsa/franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# franke_gsa
# 864 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franke_gsa/franke_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 557496 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franke_gsa/franke_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# helmsley_prism_gsa
# 753 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 243103 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# helmsley_xavier_prism_gsa
# 1281 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 243082 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# hyams_protect_gsa
# 416 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 569468 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# lewis_sparc_gsa
# 2846 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 607715 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# mccauley_new_gsa
# 1610 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 578685 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# mcgovern_gsa
# 5978 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 630955 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# moayyedi_imagine_gsa
# 1120 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 608003 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# newberry_share_gsa
# 858 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 583694 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# niddk_cho_gsa
# 1754 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 593576 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# niddk_duerr_gsa
# 1937 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 581142 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# niddk_rioux_gsa
# 913 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 574971 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# niddk_silverberg_gsa
# 2347 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 611240 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# palotie_hus_gsa
# 869 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 553667 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# pekow_share_gsa
# 632 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 574060 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# rioux_igenomed_gsa
# 181 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 543342 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# sands_msccr_gsa
# 1419 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 595558 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# stampfer_gsa
# 1464 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 643686 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# vermeire_gsa
# 4659 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/vermeire_gsa/vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 653343 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/vermeire_gsa/vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# weersma_gsa
# 705 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/weersma_gsa/weersma_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 566401 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/weersma_gsa/weersma_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# xavier_prism_gsa
# 686 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 569740 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
# xavier_share_gsa
# 692 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.fam
# 568510 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het.bim
#########
