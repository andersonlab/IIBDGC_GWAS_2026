# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
#############################################################
# 19 - COMPARE ALLELE FREQUENCIES WITH REFERENCE POPULATION #
#############################################################

##################################################################################################
# 19.1 SPLIT SAMPLES INTO CONTINENTAL ANCESTRY GROUPS BASED ON PCA TO COMPARE FRQ - NON DUPLICATES

# EUR - keep Europeans 

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500
  
for i in ${studies[@]}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_compare_freq_1_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_compare_freq_1_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2 \
--keep-allele-order --allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_nonjewish_ancestry_samples_withDuplicates \
--freq counts \
--make-bed --out ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_eur"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_compare_freq_1_${i} | grep -E "completed"
done

  
for i in ${studies[@]}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_compare_freq_2_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_compare_freq_2_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_eur \
--keep-allele-order --allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/${i}/list_controls2 \
--freq counts \
--out ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_eur_ctr"
done

# resubmit for those with cases only - keeping all samples (all cases)
studies_noctr=(finland_illugwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas mccauley_gsa ccfa_gsa farkkila_gsa hyams_protect_gsa lewis_sparc_gsa newberry_share_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa)

for i in ${studies_noctr[@]}
do echo ${i} && wc -l ${path_gwas}pre_imputation/QC/${i}/list_controls2
done

for i in ${studies_noctr[@]}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_compare_freq_2_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_compare_freq_2_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_eur \
--keep-allele-order --allow-no-sex \
--freq counts \
--out ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_eur_ctr"
done


for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_compare_freq_2_${i} | grep -E "completed"
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

for (j in 1:length(cohorts)) {
  
  print(cohorts[j])
  
  gnomad<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_gnomad_variants",sep=""),head=F)
  colnames(gnomad)<-c("SNP","CHROM","POS","REF","ALT","AF","AF_nfe","AF_afr","AF_amr","AF_eas","AF_sas","AF_asj")
  gnomad$AF_nfe<-as.numeric(as.character(gnomad$AF_nfe))
  
  ancestry<-c("eur")
  
  rm(freq_all)
  for (i in 1:length(ancestry)) {
    
    if (ancestry[i]=="eur") {
      file_anc<-paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_",ancestry[i],"_ctr.frq.counts",sep="")
    } else {
      file_anc<-paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_",ancestry[i],".frq.counts",sep="")
    }
    
    
    if (file.exists(file_anc)) {
      
      vfreq<-read.table(file_anc,head=T)
      vfreq[,paste(ancestry[i],"freq_alt",sep="_")]<-vfreq$C1/(vfreq$C1+vfreq$C2)
      
      if(!exists("freq_all")) {
        freq_all<-vfreq[,c(1:4,8)]
      }else{
        freq_all<-merge(freq_all,vfreq[,c(2,ncol(vfreq))],by="SNP",sort=F)
      }
    }
    
  }
  
  all<-merge(freq_all,gnomad,by="SNP",sort=F,all.x=T)
  
  table(all$A2,all$REF)
  
  table(all$A1,all$ALT)
  
  dim(vfreq)-dim(gnomad)
  # 
  
  table(vfreq[which(!vfreq$SNP %in% gnomad$SNP),"CHR"])
  # 
  
  ancestry<-c("eur")
  
  for (i in 1:length(ancestry)) {
    
    all[paste(ancestry[i],"value",sep="_")]<-((all[,paste(ancestry[i],"freq_alt",sep="_"),drop=F]-all[,"AF_nfe",drop=F])^2)/
      ((all[,paste(ancestry[i],"freq_alt",sep="_"),drop=F]+all[,"AF_nfe",drop=F])*
         (2-all[,paste(ancestry[i],"freq_alt",sep="_"),drop=F]-all[,"AF_nfe",drop=F]))
    
  }
  
  summary(all$eur_value)
  # 
  
  summary(all$AF_nfe[which(is.na(all$eur_value))])
  # 
  
  
  tm<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_TOPMed_variants",sep=""),head=F)
  tm<-tm[,c(9,8)]
  colnames(tm)<-c("SNP","topmed_freq_alt")
  
  all<-merge(all,tm,by="SNP",all.x=T,sort=F)
  
  
  all$AF_nfe<-as.numeric(as.character(all$AF_nfe))
  
  cohort<-cohorts[j]
  
  # use different thresholds by pop:
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9")
  
  for (i in 1:length(ancestry)) {
    
    tmp<-all[,c("CHR","POS","SNP",paste(ancestry[i],"freq_alt",sep="_"),"AF_nfe",paste(ancestry[i],"value",sep="_"),"topmed_freq_alt")]
    
    colnames(tmp)<-c("CHR","POS","SNP","freq_alt","AF_nfe","value","topmed_freq_alt")
    tmp$col<-NA
    
    if (i==1) {
      value<-0.025
    } else if (i==2) {
      value<-0.025
    } else if (i==4) {
      value<-0.0375
    } else if (i==3) {
      value<-0.05
    } else {
      value<-0.075
    }
    
    tmp$col[which(tmp$value<=value)]<-paste("<=",value,sep="")
    tmp$col[which(tmp$value>value)]<-paste(">",value,sep="")
    tmp$col[which(is.na(tmp$value))]<-"monomorphic"
    list_variants_exclude<-tmp[which(tmp$value>value),"SNP",drop=F]
    
    print(ancestry[i])
    print(table(tmp$col))
    print(dim(list_variants_exclude))
    
    p3<-ggplot(tmp, aes(y=freq_alt, x=AF_nfe)) +
      geom_point(aes(colour = col)) + ylab(paste(cohort,"FRQ Alt")) + xlab("Gnomad FRQ Alt") + scale_colour_manual(values=cbPalette) +
      ggtitle(paste(cohort,toupper(ancestry[i]),sep=" "))
    
    pdf(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_plot_maf_list_variants_toexclude_maf_differs_Gnomad_value_",ancestry[i],"_",value,".pdf",sep=""),width = 6, height = 5)
    print(p3)
    dev.off()
    system(paste("cp ",path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_plot_maf_list_variants_toexclude_maf_differs_Gnomad_value_",ancestry[i],"_",value,".pdf ~/tmp_plots/",sep=""))
    
    write.table(list_variants_exclude,paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_list_variants_toexclude_maf_differs_",ancestry[i],"_gnomad",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    
  }

  
  ### double check with TOPMED variants that remain after using TOPMed to double check frequency:
  
  all<-all[which(!all$SNP %in% list_variants_exclude$SNP),]
  
  ancestry<-c("eur")
  
  for (i in 1:length(ancestry)) {
    
    all[paste(ancestry[i],"value_tm",sep="_")]<-((all[,paste(ancestry[i],"freq_alt",sep="_"),drop=F]-all[,"topmed_freq_alt",drop=F])^2)/
      ((all[,paste(ancestry[i],"freq_alt",sep="_"),drop=F]+all[,"topmed_freq_alt",drop=F])*
         (2-all[,paste(ancestry[i],"freq_alt",sep="_"),drop=F]-all[,"topmed_freq_alt",drop=F]))
    
  }
  
  # use different thresholds by pop:
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9")
  
  for (i in 1:length(ancestry)) {
    
    tmp<-all[,c("CHR","POS","SNP",paste(ancestry[i],"freq_alt",sep="_"),"topmed_freq_alt",paste(ancestry[i],"value_tm",sep="_"))]
    
    colnames(tmp)<-c("CHR","POS","SNP","freq_alt","topmed_freq_alt","value")
    tmp$col<-NA
    
    if (i==1) {
      value<-0.125
    } 
    
    tmp$col[which(tmp$value<=value)]<-paste("<=",value,sep="")
    tmp$col[which(tmp$value>value)]<-paste(">",value,sep="")
    tmp$col[which(is.na(tmp$value))]<-"monomorphic"
    list_variants_exclude2<-tmp[which(tmp$value>value),"SNP",drop=F]
    
    print(ancestry[i])
    print(table(tmp$col))
    # print(dim(list_variants_exclude2))
    
    p3<-ggplot(tmp, aes(y=freq_alt, x=topmed_freq_alt)) +
      geom_point(aes(colour = col)) + ylab(paste(cohort,"FRQ Alt")) + xlab("TOPMed FRQ Alt") + scale_colour_manual(values=cbPalette) +
      ggtitle(paste(cohort,toupper(ancestry[i]),sep=" "))
    
    pdf(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_plot_maf_list_variants_toexclude_maf_differs_TOPMed_value_",ancestry[i],"_",value,".pdf",sep=""),width = 6, height = 5)
    print(p3)
    dev.off()
    system(paste("cp ",path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_plot_maf_list_variants_toexclude_maf_differs_TOPMed_value_",ancestry[i],"_",value,".pdf ~/tmp_plots/",sep=""))
    
    list_variants_exclude<-rbind(list_variants_exclude,list_variants_exclude2)
    print(paste("Final N to exclude:",nrow(list_variants_exclude)))
    write.table(list_variants_exclude,paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_list_variants_toexclude_maf_differs_",ancestry[i],"_gnomad",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    
  }
  
  rm(list=ls()[!ls() %in% c("path","j","cohorts")]) 
}

##########
# [1] "xavier_share_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 543981        1235        4295 
# [1] 1235    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 548229      47 
# [1] "Final N to exclude: 1282"

# [1] "australia_omniexome"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 725202         244        1743 
# [1] 244   1
# [1] "eur"
# 
# <=0.125  >0.125 
# 726852      93 
# [1] "Final N to exclude: 337"
# [1] "gwas1"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 437021         158          25 
# [1] 158   1
# [1] "eur"
# 
# <=0.125  >0.125 
# 437006      40 
# [1] "Final N to exclude: 198"
# [1] "gwas2"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 757301         350          32 
# [1] 350   1
# [1] "eur"
# 
# <=0.125  >0.125 
# 757166     167 
# [1] "Final N to exclude: 517"
# [1] "all_hce"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 373093          91       29586 
# [1] 91  1
# [1] "eur"
# 
# <=0.125      >0.125 monomorphic 
# 402616          58           5 
# [1] "Final N to exclude: 149"
# [1] "pittsburgh_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 872950         369         443 
# [1] 369   1
# [1] "eur"
# 
# <=0.125  >0.125 
# 873299      94 
# [1] "Final N to exclude: 463"
# [1] "spain_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 567231         299          41 
# [1] 299   1
# [1] "eur"
# 
# <=0.125  >0.125 
# 567253      19 
# [1] "Final N to exclude: 318"
# [1] "italy_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 542847        1440         593 
# [1] 1440    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 543388      52 
# [1] "Final N to exclude: 1492"
# [1] "kiel_austria_sibdcs_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 518975         599        9273 
# [1] 599   1
# [1] "eur"
# 
# <=0.125  >0.125 
# 528193      55 
# [1] "Final N to exclude: 654"
# [1] "netherlands_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 556252        1002        9283 
# [1] 1002    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 565440      95 
# [1] "Final N to exclude: 1097"
# [1] "slovenia_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 507806         906         106 
# [1] 906   1
# [1] "eur"
# 
# <=0.125  >0.125 
# 507833      79 
# [1] "Final N to exclude: 985"
# [1] "sweden_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 535497         935        2102 
# [1] 935   1
# [1] "eur"
# 
# <=0.125  >0.125 
# 537523      76 
# [1] "Final N to exclude: 1011"
# [1] "niddk_broad_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 555997         948       13606 
# [1] 948   1
# [1] "eur"
# 
# <=0.125  >0.125 
# 569555      48 
# [1] "Final N to exclude: 996"
# [1] "niddk_feinstein_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 561248         943       19818 
# [1] 943   1
# [1] "eur"
# 
# <=0.125  >0.125 
# 581019      47 
# [1] "Final N to exclude: 990"
# [1] "basque_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 542993         822         446 
# [1] 822   1
# [1] "eur"
# 
# <=0.125  >0.125 
# 543395      44 
# [1] "Final N to exclude: 866"
# [1] "prism_nfe_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 526254       11123        3544 
# [1] 11123     1
# [1] "eur"
# 
# <=0.125  >0.125 
# 529520     278 
# [1] "Final N to exclude: 11401"
# [1] "lithuania_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 533135        1633        1357 
# [1] 1633    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 534408      84 
# [1] "Final N to exclude: 1717"
# [1] "belgium_louis_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 541030        1310        1677 
# [1] 1310    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 542660      47 
# [1] "Final N to exclude: 1357"
# [1] "belgium_franchimont_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 532934        1448        1555 
# [1] 1448    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 534439      50 
# [1] "Final N to exclude: 1498"
# [1] "belgium_vermeire_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 577719        1980       25061 
# [1] 1980    1
# [1] "eur"
# 
# <=0.125      >0.125 monomorphic 
# 602709          69           2 
# [1] "Final N to exclude: 2049"
# [1] "prism_nfe_gwas"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 239743         100           7 
# [1] 100   1
# [1] "eur"
# 
# <=0.125  >0.125 
# 239737      13 
# [1] "Final N to exclude: 113"
# [1] "finland_illugwas"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 231872        1222           7 
# [1] 1222    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 231862      17 
# [1] "Final N to exclude: 1239"
# [1] "german_affy6_old_gwas"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 766491         487          79 
# [1] 487   1
# [1] "eur"
# 
# <=0.125  >0.125 
# 766440     130 
# [1] "Final N to exclude: 617"
# [1] "norway_affy6_old_gwas"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 708120         405          92 
# [1] 405   1
# [1] "eur"
# 
# <=0.125  >0.125 
# 708069     143 
# [1] "Final N to exclude: 548"
# [1] "belgium_inf1_old_gwas"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 298227          95           8 
# [1] 95  1
# [1] "eur"
# 
# <=0.125  >0.125 
# 298229       6 
# [1] "Final N to exclude: 101"
# [1] "belgium_inf2_old_gwas"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 286805         111           9 
# [1] 111   1
# [1] "eur"
# 
# <=0.125  >0.125 
# 286780      34 
# [1] "Final N to exclude: 145"
# [1] "niddk_old_gwas"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 296706          89           9 
# [1] 89  1
# [1] "eur"
# 
# <=0.125  >0.125 
# 296707       8 
# [1] "Final N to exclude: 97"
# [1] "cedars_370k_old_gwas"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 333577         139          41 
# [1] 139   1
# [1] "eur"
# 
# <=0.125  >0.125 
# 333600      18 
# [1] "Final N to exclude: 157"
# [1] "cedars_610k_old_gwas"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 574750         227          87 
# [1] 227   1
# [1] "eur"
# 
# <=0.125  >0.125 
# 574805      32 
# [1] "Final N to exclude: 259"
# [1] "cedars_omni_old_gwas"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 601905       12938        3239 
# [1] 12938     1
# [1] "eur"
# 
# <=0.125  >0.125 
# 603028    2116 
# [1] "Final N to exclude: 15054"
# [1] "swedish_uc_old_gwas"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 293212          90           8 
# [1] 90  1
# [1] "eur"
# 
# <=0.125  >0.125 
# 293205      15 
# [1] "Final N to exclude: 105"
# [1] "mccauley_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 536394        1192         837 
# [1] 1192    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 537194      37 
# [1] "Final N to exclude: 1229"
# [1] "ccfa_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 568957        1498        8804 
# [1] 1498    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 577710      51 
# [1] "Final N to exclude: 1549"
# [1] "cedars_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 567927        2593       18140 
# [1] 2593    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 586021      46 
# [1] "Final N to exclude: 2639"
# [1] "bernstein_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 460237       74581        1394 
# [1] 74581     1
# [1] "eur"
# 
# <=0.125  >0.125 
# 461288     343 
# [1] "Final N to exclude: 74924"
# [1] "farkkila_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 471125        8530          43 
# [1] 8530    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 471036     132 
# [1] "Final N to exclude: 8662"
# [1] "franchimont_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 547916        1556        3971 
# [1] 1556    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 551840      47 
# [1] "Final N to exclude: 1603"
# [1] "franke_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 536022        1150        2409 
# [1] 1150    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 538359      72 
# [1] "Final N to exclude: 1222"
# [1] "helmsley_prism_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 239940          85           7 
# [1] 85  1
# [1] "eur"
# 
# <=0.125  >0.125 
# 239932      15 
# [1] "Final N to exclude: 100"
# [1] "helmsley_xavier_prism_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 239978          81           7 
# [1] 81  1
# [1] "eur"
# 
# <=0.125  >0.125 
# 239968      17 
# [1] "Final N to exclude: 98"
# [1] "hyams_protect_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 545655        1259        3824 
# [1] 1259    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 549437      42 
# [1] "Final N to exclude: 1301"
# [1] "lewis_sparc_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 573778        1487        9908 
# [1] 1487    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 583637      49 
# [1] "Final N to exclude: 1536"
# [1] "mccauley_new_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 552209        1711        5007 
# [1] 1711    1
# [1] "eur"
# 
# <=0.125      >0.125 monomorphic 
# 557206           9           1 
# [1] "Final N to exclude: 1720"
# [1] "mcgovern_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 577356        2710       23168 
# [1] 2710    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 600477      47 
# [1] "Final N to exclude: 2757"
# [1] "moayyedi_imagine_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 559680        1140        9209 
# [1] 1140    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 568807      82 
# [1] "Final N to exclude: 1222"
# [1] "newberry_share_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 551490        1166        4667 
# [1] 1166    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 556105      52 
# [1] "Final N to exclude: 1218"
# [1] "niddk_cho_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 561554        1444        9372 
# [1] 1444    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 570877      49 
# [1] "Final N to exclude: 1493"
# [1] "niddk_duerr_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 554598        1365        4860 
# [1] 1365    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 559403      55 
# [1] "Final N to exclude: 1420"
# [1] "niddk_rioux_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 545993        1301        3978 
# [1] 1301    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 549927      44 
# [1] "Final N to exclude: 1345"
# [1] "niddk_silverberg_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 565633        1521       12274 
# [1] 1521    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 577869      38 
# [1] "Final N to exclude: 1559"
# [1] "palotie_hus_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 531780        3177         395 
# [1] 3177    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 532079      96 
# [1] "Final N to exclude: 3273"
# [1] "pekow_share_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 549548        1184        3721 
# [1] 1184    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 553216      53 
# [1] "Final N to exclude: 1237"
# [1] "rioux_igenomed_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 524332        1122         780 
# [1] 1122    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 525049      63 
# [1] "Final N to exclude: 1185"
# [1] "sands_msccr_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 563674        1265       10376 
# [1] 1265    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 573991      59 
# [1] "Final N to exclude: 1324"
# [1] "stampfer_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 544893        2440        3671 
# [1] 2440    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 548501      63 
# [1] "Final N to exclude: 2503"
# [1] "vermeire_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 554772        1586        7396 
# [1] 1586    1
# [1] "eur"
# 
# <=0.125      >0.125 monomorphic 
# 562099          68           1 
# [1] "Final N to exclude: 1654"
# [1] "weersma_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 539026        1087        2567 
# [1] 1087    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 541502      91 
# [1] "Final N to exclude: 1178"
# [1] "xavier_prism_gsa"
# [1] "eur"
# 
# <=0.025      >0.025 monomorphic 
# 542545        2201        5240 
# [1] 2201    1
# [1] "eur"
# 
# <=0.125  >0.125 
# 547633     152 
# [1] "Final N to exclude: 2353"

##########

#########################
# 19.2 EXCLUDE VARIANTS


path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

for i in ${studies[@]}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_compare_freq_3_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_compare_freq_3_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2 \
--keep-allele-order --allow-no-sex \
--exclude ${path_gwas}pre_imputation/QC/${i}/${i}_list_variants_toexclude_maf_differs_eur_gnomad \
--make-bed --out ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_compare_freq_3_${i} | grep -E "completed"
done





