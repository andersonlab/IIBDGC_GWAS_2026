# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
#######################################
# 24.- INTRA-STUDY DUPLICATED SAMPLES #
#######################################

##################################################
# 24.1.- REMOVE INTRA-STUDY DUPLICATED SAMPLES - OUT OF SAMPLES REMAINING AT END OF QC / pre-liftover to b37

path_gwas=/path/to/ibdgwas/IIBDGC/

MEM=3500
studies=(all_hce)

MEM=2800
studies=(niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do 
bsub -J"kg1" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_king_25.1_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_king_25.1_${i} \
"/path/to/software/./king \
-b ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24.bed \
--related --cpus 4 --prefix ${path_gwas}pre_imputation/QC/${i}/${i}_king_step25"
done

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_king_25.1_${i} | grep -E "completed"
done



## /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"
studies<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
           "kiel_austria_sibdcs_gsa","netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa",
           "niddk_feinstein_gsa","basque_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa",
           "belgium_vermeire_gsa","prism_nfe_gsa","prism_nfe_gwas","finland_illugwas","german_affy6_old_gwas",
           "norway_affy6_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
           "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas",
           "swedish_uc_old_gwas",
           "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

studies_new<-c("bernstein_gsa","farkkila_gsa","franchimont_gsa",
               "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
               "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
               "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
               "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
               "xavier_share_gsa")



################################################################################################
# RE-do case control status for all new studies (batches from 2021 and 2022)

pheno<-fread(paste(path,"pheno/gwas-mega-2-core-phenotypes-28ceede.csv.gz",sep=""),head=T)
pheno<-as.data.frame(pheno)

colnames(pheno)[6]<-"gender"

pheno$sex<-0
pheno$sex[which(pheno$gender=="Female")]<-2
pheno$sex[which(pheno$gender=="Male")]<-1
table(pheno$gender,pheno$sex,useNA="ifany")
#             0     1     2
#          3261     0     0
# Female      0     0 55427
# Male        0 51830     0
# Unknown  8872     0     0


pheno$pheno<-NA
pheno$pheno[which(pheno$control=="1")]<-"1"
pheno$pheno[which(pheno$affection=="Affected")]<-"2"
table(pheno$pheno,pheno$affection,useNA="ifany")
#            Affected Unaffected Unknown
# 1        0        0      28970       0
# 2        0    85835          0       0
# <NA>  3967        0        143     475
table(pheno$pheno,pheno$diag,useNA="ifany")
#            Crohn's Disease Indeterminate Ulcerative Colitis Unknown
# 1    28970               0             0                  0       0
# 2       17           49394          1221              31608    3595
# <NA>  4345               0             0                  0     240


### add sex from latest update (see step24_define_final_inferred_sex)

sex<-fread(paste(path,"pheno/gwas-mega-2-core_updated_sex_gender_Dec22.txt.gz",sep=""),head=T)

################################################################
# collect all post-QC fam files, before being split into of ancestry

for (j in 1:length(studies)) {
  
  print(studies[j])
  fam_tmp<-fread(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24.fam",sep=""),head=F)
  fam_tmp<-as.data.frame(fam_tmp)
  fam_tmp$study<-studies[j]
  
  # add sex and pheno from updated files
  
  if (!studies[j] %in% c("all_hce","gwas1","gwas2","german_affy6_old_gwas",
                         "norway_affy6_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
                         "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas",
                         "swedish_uc_old_gwas","niddk_old_gwas")) {
    fam_tmp<-merge(fam_tmp,pheno[,c("sample_id","pheno")],by.x="V2",by.y="sample_id",all.x=T,sort=F)
  } else {
    fam_tmp$pheno<-fam_tmp$V6
  }
  
  if(!studies[j] %in% c("swedish_uc_old_gwas")) {
    fam_tmp<-merge(fam_tmp,sex[,c("IID","SNPSEX2")],by.x="V2",by.y="IID",all.x=T,sort=F)
  } else {
    fam_tmp$SNPSEX2<-fam_tmp$V5
  }
  
  
  if(j==1) {
    fam_all<-fam_tmp
  }else{
    fam_all<-rbind(fam_all,fam_tmp)
  }
  
  rm(fam_tmp)
  
}

colnames(fam_all)[9]<-"sex"
table(fam_all$study,fam_all$pheno,useNA="ifany")
table(fam_all$study,fam_all$sex,useNA="ifany")

fwrite(fam_all,paste(path,"pheno/gwas-mega-2-core_updated_sex_gender_phenotype_Dec22.txt.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

data_final<-as.data.frame(matrix(nrow=length(studies),ncol=4))
colnames(data_final)<-c("study","n_excluded_inconsist","n_excluded_eur","n_excluded_noneur")
data_final$study<-studies

for (j in 1:length(studies)) {
  
  print(studies[j])
  
  file.kin<-paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_king_step25.kin0",sep="")
  
  if(file.exists(file.kin)) {
    
    kin<-read.table(file.kin,head=T)
    print(table(kin$InfType))
    
    print(table(cut(kin$Kinship,breaks=c(1,0.354,0.177))))
    
    dup<-kin[which(kin$Kinship>0.354 & kin$InfType=="Dup/MZ"),c("FID1","FID2","Kinship")]
    
    # file_to_exclude<-paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_list_samples_sex_discrepancy_step24",sep="")
    # if(file.exists(file_to_exclude)) {
    #   to_exclude<-fread(file_to_exclude,head=T)
    #   dup<-dup[which(!dup$FID1 %in% to_exclude$FID),]
    # }
    
    
    if(nrow(dup)>0) {
      
      dup_ids<-c(as.character(dup$FID1),as.character(dup$FID2))
      length(dup_ids)
      table(length(dup_ids)==(nrow(dup)*2))
      
      dup_ids<-dup_ids[!duplicated(dup_ids)]
      length(dup_ids)
      
      # CHECK ALSO THAT BOTH SAMPLES HAVE SAME PHENOTYPE AND GENDER, AND REMOVE SAMPLES LOWER CALL RATES
      
      fam<-fam_all[which(fam_all$study==studies[j]),]
      sample_miss<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy.imiss",sep=""),head=T)
      
      all<-merge(fam[,c("V1","sex","pheno")],sample_miss[,c("FID","F_MISS")],by.x="V1",by.y="FID",all.x=T,sort=F)
      
      tmp<-all[which(all$V1 %in% dup_ids),]
      dup_ids<-dup_ids[match(tmp$V1,dup_ids)]
      rm(tmp)
      
      for (i in 1:length(dup_ids)) {
        
        tmp1<-dup[which(dup$FID1==dup_ids[i]),]
        tmp2<-dup[which(dup$FID2==dup_ids[i]),]
        colnames(tmp2)[1:2]<-colnames(tmp2)[2:1]
        
        tmp<-rbind(tmp1,tmp2)
        
        ids_tmp<-c(as.character(tmp$FID1),as.character(tmp$FID2))
        ids_tmp<-ids_tmp[!duplicated(ids_tmp)]
        
        #keep same order as in dup_ids and in all
        ids_tmp<-ids_tmp[match(dup_ids[which(dup_ids %in% ids_tmp)],ids_tmp)]
        
        # get number of possible combinations
        n_possible_combinations<-nrow(permutations(length(ids_tmp), 2, letters[1:length(ids_tmp)]))/2
        
        
        if ( nrow(dup[which((dup$FID1 %in% ids_tmp) | (dup$FID2 %in% ids_tmp)),])==n_possible_combinations ) {
          
          # for duplicated samples that have same phenotype and sex, remove the ones with smaller call rate
          
          data<-as.data.frame(matrix(ncol=1,nrow=length(ids_tmp)))
          data$V1<-ids_tmp
          data<-merge(data,all,by="V1",all.x=T,sort=F)
          
          if ( (dim(table(data$sex))==1 & dim(table(data$pheno))==1) ) {
            
            keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]
            
            if(!exists("data_remove")) {
              data_remove<-data[which(!data$V1 %in% keep_sample),]
            } else {
              data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
            }
            
          } 
          
          # to deal with sex = 0
          else if ( any(data$sex==0) & any(data$sex %in% c(1,2)) ) {
            
            keep_sample<-data$V1[which(data$sex!=0)]
            
            if(!exists("data_remove")) {
              data_remove<-data[which(!data$V1 %in% keep_sample),]
            } else {
              data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
            }
            
          }
          
          else {
            
            data<-as.data.frame(matrix(ncol=1,nrow=length(ids_tmp)))
            data$V1<-ids_tmp
            data<-merge(data,all,by="V1",all.x=T,sort=F)
            
            # for duplicated samples that do not have same pheno and sex, remove all
            
            # print("Samples with different sex/pheno:")
            # print(data)
            
            if(!exists("data_remove")) {
              
              data_remove<-data
              
            } else {
              
              data_remove<-rbind(data_remove,data)
              
            }
            
            if(!exists("data_inconsist")){
              
              jj<-1
              data_inconsist<-data
              data_inconsist$group<-jj
              
            } else {
              
              jj<-jj+1
              data$group<-jj
              data_inconsist<-rbind(data_inconsist,data)
              
            }
          }
          
        } else {
          
          # not all combinations of duplicated pairs found, show list of IDs, to manually inspect issues
          
          print(paste("Number of expected combinations: ",n_possible_combinations,sep=""))
          print(paste("Number of observed combinations: ",nrow(dup[which((dup$FID1 %in% ids_tmp) | (dup$FID2 %in% ids_tmp)),]),sep=""))
          print(ids_tmp)
          
        }
        
      }
      
      data_remove<-data_remove[!duplicated(data_remove$V1),]
      # print(dim(data_remove))
      
      data_remove<-data_remove[,c(1,1)]
      colnames(data_remove)<-c("FID","IID")
      write.table(data_remove,paste(path,"pre_imputation/QC/",studies[j],"/list_duplicated_samples_step25",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
      
      
      #### create a file to report inconsistend data:
      
      if(exists("data_inconsist")){
        
        data_inconsist<-data_inconsist[!duplicated(data_inconsist$V1),]
        print(paste("N samples with inconsistent sex|pheno:",nrow(data_inconsist)))
        
        data_final$n_excluded_inconsist[j]<-nrow(data_inconsist)
        
        data_inconsist$cohort<-NA
        
        if (studies[j] %in% c("all_hce","australia_omniexome","kiel_austria_sibdcs_gsa")) {
          if(studies[j]=="all_hce") {
            cohorts<-c("gwas3","new_wave")
          } else if (studies[j]=="australia_omniexome") {
            cohorts<-c("australia_2012_omniexome","australia_2014_omniexome")
          } else if (studies[j]=="kiel_austria_sibdcs_gsa") {
            cohorts<-c("kiel_foc_gsa","kiel_eze_gsa","kiel_bc_gsa","kiel_ibd_gsa"
                       ,"kiel_hlitm_gsa","austria_gsa","sibdcs_gsa")
          }
          
          for (i in 1:length(cohorts)){
            
            file_fam<-paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19.fam",sep="")
            
            if (studies[j]=="all_hce" & cohorts[i]=="gwas3") {
              file_fam<-paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_tmp.fam",sep="")
            }
            
            if (studies[j]=="kiel_austria_sibdcs_gsa") {
              file_fam<-paste(path,"pre_imputation/QC/",studies[j],"/",cohorts[i],"_hg19.fam",sep="")
            }
            
            fam<-read.table(file_fam,head=F)
            
            data_inconsist$cohort[which(data_inconsist$V1 %in% fam$V1)]<-cohorts[i]
            
          }
        } 
        
        data_inconsist
        write.table(data_inconsist,paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_list_duplicated_samples_inconsistent_phenotype_step25",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
      }
    } else {
      data_remove<-as.data.frame(matrix(ncol=2,nrow=0))
      colnames(data_remove)<-c("FID","IID")
      write.table(data_remove,paste(path,"pre_imputation/QC/",studies[j],"/list_duplicated_samples_step25",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    }
  } else {
    data_remove<-as.data.frame(matrix(ncol=2,nrow=0))
    colnames(data_remove)<-c("FID","IID")
    write.table(data_remove,paste(path,"pre_imputation/QC/",studies[j],"/list_duplicated_samples_step25",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  }
  
  
  print(paste("Number of samples to exclude: ",nrow(data_remove),sep=""))
  eur<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam",sep=""),head=F)
  
  print(paste("Number of EUR samples to exclude: ",nrow(data_remove[which(data_remove$IID %in% eur$V1),]),sep=""))
  data_final$n_excluded_eur[j]<-nrow(data_remove[which(data_remove$IID %in% eur$V1),])
  
  file_noneur<-paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam",sep="")
  
  if(file.exists(file_noneur)){
    noneur<-read.table(file_noneur,head=F)
    print(paste("Number of nonEUR samples to exclude: ",nrow(data_remove[which(data_remove$IID %in% noneur$V1),]),sep=""))
    data_final$n_excluded_noneur[j]<-nrow(data_remove[which(data_remove$IID %in% noneur$V1),])
    
  }
  
  
  rm(list=ls()[!ls() %in% c("path","j","studies","fam_all","data_final")])
}


write.table(data_final,"~/tmp_plots/summary_duplicates_excuded_step25.tsv",col.names=T,row.names=F,quote=F,sep="\t")

  
######

# [1] "all_hce"
# 
# 3rd Dup/MZ     FS     PO 
# 1   1389    186    915 
# 
# (0.177,0.354]     (0.354,1] 
# 1101          1390 
# [1] "N samples with inconsistent sex|pheno: 2"
# [1] "Number of samples to exclude: 1265"
# [1] "Number of EUR samples to exclude: 1161"
# [1] "Number of nonEUR samples to exclude: 104"
# [1] "niddk_old_gwas"
# 
# Dup/MZ     FS     PO 
# 8      8      6 
# 
# (0.177,0.354]     (0.354,1] 
# 14             8 
# [1] "Number of samples to exclude: 8"
# [1] "Number of EUR samples to exclude: 8"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "australia_omniexome"
# 
# Dup/MZ     FS     PO 
# 8      1      3 
# 
# (0.177,0.354]     (0.354,1] 
# 4             8 
# [1] "Number of samples to exclude: 8"
# [1] "Number of EUR samples to exclude: 7"
# [1] "Number of nonEUR samples to exclude: 1"
# [1] "gwas1"
# 
# Dup/MZ 
# 1 
# 
# (0.177,0.354]     (0.354,1] 
# 0             1 
# [1] "Number of samples to exclude: 1"
# [1] "Number of EUR samples to exclude: 1"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "gwas2"
# 
# Dup/MZ     FS 
# 1      1 
# 
# (0.177,0.354]     (0.354,1] 
# 1             1 
# [1] "Number of samples to exclude: 1"
# [1] "Number of EUR samples to exclude: 1"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "pittsburgh_gsa"
# 
# FS PO 
# 18  3 
# 
# (0.177,0.354]     (0.354,1] 
# 21             0 
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "spain_gsa"
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "italy_gsa"
# 
# Dup/MZ     FS     PO 
# 3      2      8 
# 
# (0.177,0.354]     (0.354,1] 
# 10             3 
# [1] "Number of samples to exclude: 3"
# [1] "Number of EUR samples to exclude: 3"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "kiel_austria_sibdcs_gsa"
# 
# Dup/MZ     FS     PO 
# 537    176    301 
# 
# (0.177,0.354]     (0.354,1] 
# 477           537 
# [1] "N samples with inconsistent sex|pheno: 70"
# [1] "Number of samples to exclude: 494"
# [1] "Number of EUR samples to exclude: 479"
# [1] "Number of nonEUR samples to exclude: 15"
# [1] "netherlands_gsa"
# 
# Dup/MZ     FS     PO 
# 20     35     36 
# 
# (0.177,0.354]     (0.354,1] 
# 71            20 
# [1] "N samples with inconsistent sex|pheno: 4"
# [1] "Number of samples to exclude: 22"
# [1] "Number of EUR samples to exclude: 22"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "slovenia_gsa"
# 
# Dup/MZ     PO 
# 2      1 
# 
# (0.177,0.354]     (0.354,1] 
# 1             2 
# [1] "Number of samples to exclude: 2"
# [1] "Number of EUR samples to exclude: 2"
# [1] "sweden_gsa"
# 
# Dup/MZ     FS     PO 
# 21      4     14 
# 
# (0.177,0.354]     (0.354,1] 
# 18            21 
# [1] "N samples with inconsistent sex|pheno: 2"
# [1] "Number of samples to exclude: 22"
# [1] "Number of EUR samples to exclude: 22"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "niddk_broad_gsa"
# 
# Dup/MZ     FS     PO 
# 52     29     41 
# 
# (0.177,0.354]     (0.354,1] 
# 70            52 
# [1] "N samples with inconsistent sex|pheno: 4"
# [1] "Number of samples to exclude: 53"
# [1] "Number of EUR samples to exclude: 47"
# [1] "Number of nonEUR samples to exclude: 6"
# [1] "niddk_feinstein_gsa"
# 
# Dup/MZ     FS     PO 
# 208     92    114 
# 
# (0.177,0.354]     (0.354,1] 
# 206           208 
# [1] "N samples with inconsistent sex|pheno: 5"
# [1] "Number of samples to exclude: 200"
# [1] "Number of EUR samples to exclude: 164"
# [1] "Number of nonEUR samples to exclude: 36"
# [1] "basque_gsa"
# 
# Dup/MZ     FS     PO 
# 17     15     10 
# 
# (0.177,0.354]     (0.354,1] 
# 25            17 
# [1] "Number of samples to exclude: 17"
# [1] "Number of EUR samples to exclude: 17"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "lithuania_gsa"
# 
# Dup/MZ     FS     PO 
# 21      8      9 
# 
# (0.177,0.354]     (0.354,1] 
# 17            21 
# [1] "Number of samples to exclude: 21"
# [1] "Number of EUR samples to exclude: 21"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "belgium_louis_gsa"
# 
# Dup/MZ     FS     PO 
# 6     10    117 
# 
# (0.177,0.354]     (0.354,1] 
# 127             6 
# [1] "Number of samples to exclude: 6"
# [1] "Number of EUR samples to exclude: 6"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "belgium_franchimont_gsa"
# 
# Dup/MZ     FS     PO 
# 7      9     16 
# 
# (0.177,0.354]     (0.354,1] 
# 25             7 
# [1] "Number of samples to exclude: 7"
# [1] "Number of EUR samples to exclude: 7"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "belgium_vermeire_gsa"
# 
# Dup/MZ     FS     PO 
# 9     12     13 
# 
# (0.177,0.354]     (0.354,1] 
# 25             9 
# [1] "Number of samples to exclude: 9"
# [1] "Number of EUR samples to exclude: 9"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "prism_nfe_gsa"
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "prism_nfe_gwas"
# 
# Dup/MZ     FS     PO 
# 1      6     15 
# 
# (0.177,0.354]     (0.354,1] 
# 15             1 
# [1] "Number of samples to exclude: 1"
# [1] "Number of EUR samples to exclude: 1"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "finland_illugwas"
# 
# Dup/MZ 
# 2 
# 
# (0.177,0.354]     (0.354,1] 
# 0             2 
# [1] "Number of samples to exclude: 2"
# [1] "Number of EUR samples to exclude: 2"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "german_affy6_old_gwas"
# 
# Dup/MZ     FS     PO 
# 56     11     65 
# 
# (0.177,0.354]     (0.354,1] 
# 76            56 
# [1] "Number of samples to exclude: 36"
# [1] "Number of EUR samples to exclude: 36"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "norway_affy6_old_gwas"
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "belgium_inf1_old_gwas"
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "belgium_inf2_old_gwas"
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "cedars_370k_old_gwas"
# 
# FS PO 
# 4  3 
# 
# (0.177,0.354]     (0.354,1] 
# 7             0 
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "cedars_610k_old_gwas"
# 
# FS PO 
# 6  6 
# 
# (0.177,0.354]     (0.354,1] 
# 12             0 
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "cedars_omni_old_gwas"
# 
# Dup/MZ     FS 
# 2      1 
# 
# (0.177,0.354]     (0.354,1] 
# 1             2 
# [1] "Number of samples to exclude: 2"
# [1] "Number of EUR samples to exclude: 1"
# [1] "Number of nonEUR samples to exclude: 1"
# [1] "swedish_uc_old_gwas"
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "mccauley_gsa"
# 
# FS PO 
# 5  3 
# 
# (0.177,0.354]     (0.354,1] 
# 8             0 
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "ccfa_gsa"
# 
# Dup/MZ     FS     PO 
# 1      4      4 
# 
# (0.177,0.354]     (0.354,1] 
# 8             1 
# [1] "Number of samples to exclude: 1"
# [1] "Number of EUR samples to exclude: 1"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "cedars_gsa"
# 
# Dup/MZ     FS     PO 
# 2     31     97 
# 
# (0.177,0.354]     (0.354,1] 
# 128             2 
# [1] "Number of samples to exclude: 2"
# [1] "Number of EUR samples to exclude: 2"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "bernstein_gsa"
# 
# Dup/MZ     FS     PO 
# 31      2      6 
# 
# (0.177,0.354]     (0.354,1] 
# 8            31 
# [1] "Number of samples to exclude: 30"
# [1] "Number of EUR samples to exclude: 29"
# [1] "Number of nonEUR samples to exclude: 1"
# [1] "farkkila_gsa"
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "franchimont_gsa"
# 
# Dup/MZ     FS     PO 
# 51     20     27 
# 
# (0.177,0.354]     (0.354,1] 
# 47            51 
# [1] "Number of samples to exclude: 51"
# [1] "Number of EUR samples to exclude: 45"
# [1] "Number of nonEUR samples to exclude: 6"
# [1] "franke_gsa"
# 
# Dup/MZ     FS 
# 7      1 
# 
# (0.177,0.354]     (0.354,1] 
# 1             7 
# [1] "Number of samples to exclude: 7"
# [1] "Number of EUR samples to exclude: 7"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "helmsley_prism_gsa"
# 
# Dup/MZ     FS     PO 
# 1      3     10 
# 
# (0.177,0.354]     (0.354,1] 
# 13             1 
# [1] "Number of samples to exclude: 1"
# [1] "Number of EUR samples to exclude: 1"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "helmsley_xavier_prism_gsa"
# 
# FS PO 
# 8  8 
# 
# (0.177,0.354]     (0.354,1] 
# 16             0 
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "hyams_protect_gsa"
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "lewis_sparc_gsa"
# 
# Dup/MZ     FS     PO 
# 2      4      4 
# 
# (0.177,0.354]     (0.354,1] 
# 8             2 
# [1] "Number of samples to exclude: 2"
# [1] "Number of EUR samples to exclude: 2"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "mccauley_new_gsa"
# 
# Dup/MZ     FS     PO 
# 7      7     22 
# 
# (0.177,0.354]     (0.354,1] 
# 29             7 
# [1] "Number of samples to exclude: 7"
# [1] "Number of EUR samples to exclude: 5"
# [1] "Number of nonEUR samples to exclude: 2"
# [1] "mcgovern_gsa"
# 
# Dup/MZ     FS     PO 
# 40    286    508 
# 
# (0.177,0.354]     (0.354,1] 
# 794            40 
# [1] "N samples with inconsistent sex|pheno: 42"
# [1] "Number of samples to exclude: 61"
# [1] "Number of EUR samples to exclude: 59"
# [1] "Number of nonEUR samples to exclude: 2"
# [1] "moayyedi_imagine_gsa"
# 
# FS PO 
# 3 18 
# 
# (0.177,0.354]     (0.354,1] 
# 21             0 
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "newberry_share_gsa"
# 
# Dup/MZ     FS     PO 
# 2      3      2 
# 
# (0.177,0.354]     (0.354,1] 
# 5             2 
# [1] "Number of samples to exclude: 2"
# [1] "Number of EUR samples to exclude: 2"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "niddk_cho_gsa"
# 
# Dup/MZ     FS     PO 
# 16      8      3 
# 
# (0.177,0.354]     (0.354,1] 
# 11            16 
# [1] "N samples with inconsistent sex|pheno: 6"
# [1] "Number of samples to exclude: 18"
# [1] "Number of EUR samples to exclude: 16"
# [1] "Number of nonEUR samples to exclude: 2"
# [1] "niddk_duerr_gsa"
# 
# Dup/MZ     FS     PO 
# 120     29     28 
# 
# (0.177,0.354]     (0.354,1] 
# 57           120 
# [1] "N samples with inconsistent sex|pheno: 12"
# [1] "Number of samples to exclude: 126"
# [1] "Number of EUR samples to exclude: 125"
# [1] "Number of nonEUR samples to exclude: 1"
# [1] "niddk_rioux_gsa"
# 
# Dup/MZ     FS     PO 
# 5      7     15 
# 
# (0.177,0.354]     (0.354,1] 
# 22             5 
# [1] "N samples with inconsistent sex|pheno: 2"
# [1] "Number of samples to exclude: 6"
# [1] "Number of EUR samples to exclude: 6"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "niddk_silverberg_gsa"
# 
# Dup/MZ     FS     PO 
# 39     25     50 
# 
# (0.177,0.354]     (0.354,1] 
# 75            39 
# [1] "N samples with inconsistent sex|pheno: 8"
# [1] "Number of samples to exclude: 40"
# [1] "Number of EUR samples to exclude: 36"
# [1] "Number of nonEUR samples to exclude: 4"
# [1] "palotie_hus_gsa"
# 
# FS PO 
# 2  1 
# 
# (0.177,0.354]     (0.354,1] 
# 3             0 
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "pekow_share_gsa"
# 
# FS PO 
# 1  1 
# 
# (0.177,0.354]     (0.354,1] 
# 2             0 
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "rioux_igenomed_gsa"
# 
# Dup/MZ 
# 6 
# 
# (0.177,0.354]     (0.354,1] 
# 0             6 
# [1] "Number of samples to exclude: 6"
# [1] "Number of EUR samples to exclude: 6"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "sands_msccr_gsa"
# 
# Dup/MZ     FS     PO 
# 7      2      8 
# 
# (0.177,0.354]     (0.354,1] 
# 10             7 
# [1] "N samples with inconsistent sex|pheno: 2"
# [1] "Number of samples to exclude: 7"
# [1] "Number of EUR samples to exclude: 4"
# [1] "Number of nonEUR samples to exclude: 3"
# [1] "stampfer_gsa"
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "vermeire_gsa"
# 
# Dup/MZ     FS     PO 
# 14    104     85 
# 
# (0.177,0.354]     (0.354,1] 
# 189            14 
# [1] "Number of samples to exclude: 14"
# [1] "Number of EUR samples to exclude: 14"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "weersma_gsa"
# 
# FS PO 
# 3  7 
# 
# (0.177,0.354]     (0.354,1] 
# 10             0 
# [1] "Number of samples to exclude: 0"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 0"
# [1] "xavier_prism_gsa"
# 
# Dup/MZ     PO 
# 2      1 
# 
# (0.177,0.354]     (0.354,1] 
# 1             2 
# [1] "Number of samples to exclude: 2"
# [1] "Number of EUR samples to exclude: 0"
# [1] "Number of nonEUR samples to exclude: 2"
# [1] "xavier_share_gsa"
# 
# Dup/MZ     FS     PO 
# 3      1      2 
# 
# (0.177,0.354]     (0.354,1] 
# 3             3 
# [1] "Number of samples to exclude: 3"
# [1] "Number of EUR samples to exclude: 3"
# [1] "Number of nonEUR samples to exclude: 0"

######

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_25.2_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_25.2_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24 \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/${i}/list_duplicated_samples_step25 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_25.2_${i} | grep -E "completed"
done



############################

################################################
# 25.2 - REMOVE INTRA-STUDY DUPLICATED SAMPLES #
################################################

# list variants set up in script - create_list_common_variants_among_cohorts
# evaluating relatedness per country / cohorts in:

# evaluate_relatedness_UKIBDGC_step25.R - COMPLETED
# evaluate_relatedness_belgium_step25.R - COMPLETED
# evaluate_relatedness_german_step25.R - COMPLETED
# evaluate_relatedness_spain_step25.R - COMPLETED
# evaluate_relatedness_niddk_step25.R - COMPLETED
# evaluate_relatedness_cedars_step25.R - COMPLETED 
# evaluate_relatedness_SWEDEN_step25.R - COMPLETED 
# evaluate_relatedness_netherlands_step25.R - COMPLETED
# evaluate_relatedness_prism_step25.R - COMPLETED
# evaluate_relatedness_mccauley.R - COMPLETED
# evaluate_relatedness_ccfa.R - COMPLETED

# updated in 2022
path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_25.3_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_25.3_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy"
done

# ukibdgc
studies=(gwas1 gwas2 all_hce)

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_25.3_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_25.3_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_ukiibdgc_data_remove_step25 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy"
done

# german

studies=(kiel_austria_sibdcs_gsa german_illu550_old_gwas german_affy6_old_gwas franke_gsa)
for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_25.3_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_25.3_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_german_data_remove_step25 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy"
done


# niddk
studies=(niddk_broad_gsa niddk_feinstein_gsa niddk_old_gwas niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa)

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_25.3_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_25.3_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_niddk_data_remove_step25 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy"
done

# spain - no related

# belgium
studies=(belgium_inf1_old_gwas belgium_inf2_old_gwas belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa 
         vermeire_gsa franchimont_gsa)
for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_25.3_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_25.3_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_belgium_data_remove_step25 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy"
done

# prism
studies=(prism_nfe_gsa prism_nfe_gwas helmsley_prism_gsa helmsley_xavier_prism_gsa xavier_prism_gsa xavier_share_gsa)
for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_34_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_25.3_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_prism_data_remove_step25 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy"
done

# sweden
studies=(swedish_uc_old_gwas sweden_gsa)
for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_25.3_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_25.3_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_sweden_data_remove_step25 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy"
done

# cedars
studies=(cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas cedars_gsa mcgovern_gsa)
for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_25.3_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_25.3_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_cedars_data_remove_step25 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy"
done

# netherlands
studies=(netherlands_gsa weersma_gsa)
for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_25.3_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_25.3_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_netherlands_data_remove_step25 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy"
done

# mccauley
studies=(mccauley_new_gsa mccauley_gsa)
for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_25.3_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_25.3_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_mccauley_data_remove_step25 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy"
done

# ccfa
studies=(ccfa_gsa lewis_sparc_gsa)
for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_25.3_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_25.3_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_ccfa_data_remove_step25 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy"
done


studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_25.3_${i} | grep -E "completed"
done


############################

# get summary per ancestry:

### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

studies<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
           "kiel_austria_sibdcs_gsa","netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa",
           "niddk_feinstein_gsa","basque_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa",
           "belgium_vermeire_gsa","prism_nfe_gsa","prism_nfe_gwas","finland_illugwas","german_affy6_old_gwas",
           "norway_affy6_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
           "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas",
           "swedish_uc_old_gwas",
           "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

for (j in 1:length(studies)) {
  
  eur_tmp<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.fam",sep=""),head=F)
  eur_tmp$study<-studies[j]
  
  file_noneur<-paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.fam",sep="")
  if(file.exists(file_noneur)){
    noneur_tmp<-read.table(file_noneur,head=F)
    noneur_tmp$study<-studies[j]
    if(j==1) {
      noneur<-noneur_tmp
    }else{
      noneur<-rbind(noneur,noneur_tmp)
    }
  }
  
  if(j==1) {
    eur<-eur_tmp
  }else{
    eur<-rbind(eur,eur_tmp)
  }
}


centers<-c("ukibdgc","german","niddk","belgium","prism","sweden","cedars","netherlands","mccauley","ccfa")


for (j in 1:length(centers)){
  
  print(centers[j])
  tmp<-read.table(paste(path,"pre_imputation/QC/relatedness/",centers[j],"_data_remove_step25",sep=""),head=T)
 
  print(paste("N samples exclude:",nrow(tmp)))
  
  print("EUR")
  print(table(tmp$cohort[which(tmp$V1 %in% eur$V1)]))
  
  print("nonEUR")
  print(table(tmp$cohort[which(tmp$V1 %in% noneur$V1)]))
  
  tmp_file<-paste(path,"pre_imputation/QC/relatedness/",centers[j],"_data_inconsist_step25",sep="")
  if(file.exists(tmp_file)) {
    tmp<-read.table(tmp_file,head=T)
    print("N samples different pheno")
    print(table(tmp$cohort))
  } else {
    print("No pheno discrepancies")
  }
  rm(tmp,tmp_file)
  
}

##################
# [1] "ukibdgc"
# [1] "N samples exclude: 2831"
# 
# [1] "EUR"
# all_hce   gwas1   gwas2 
# 184       3    2638 
# 
# [1] "nonEUR"
# all_hce   gwas1   gwas2 
# 2       0       4 
# 
# [1] "N samples different pheno"
# all_hce   gwas1   gwas2 
# 3       2       2 


# [1] "german"
# [1] "N samples exclude: 260"
# [1] "EUR"
# franke_gsa   german_affy6_old_gwas kiel_austria_sibdcs_gsa 
# 90                     165                       5
# 
# [1] "nonEUR"
# franke_gsa   german_affy6_old_gwas kiel_austria_sibdcs_gsa 
# 0                       0                       0 
# 
# [1] "N samples different pheno"
# franke_gsa   german_affy6_old_gwas kiel_austria_sibdcs_gsa 
# 4                       2                       4 


# [1] "niddk"
# [1] "N samples exclude: 8162"

# [1] "EUR"
# niddk_broad_gsa        niddk_cho_gsa      niddk_duerr_gsa 
#            4395                  293                  541 
# niddk_feinstein_gsa       niddk_old_gwas      niddk_rioux_gsa 
# 491                 1555                  267 
# niddk_silverberg_gsa 
# 388 

# [1] "nonEUR"
# niddk_broad_gsa        niddk_cho_gsa      niddk_duerr_gsa 
# 152                    1                    4 
# niddk_feinstein_gsa       niddk_old_gwas      niddk_rioux_gsa 
# 57                   12                    2 
# niddk_silverberg_gsa 
# 4 

# [1] "N samples different pheno"
# niddk_broad_gsa        niddk_cho_gsa      niddk_duerr_gsa 
# 1338                  284                  437 
# niddk_feinstein_gsa       niddk_old_gwas      niddk_rioux_gsa 
# 101                  446                  262 
# niddk_silverberg_gsa 
# 372 


# [1] "N samples different pheno"
# niddk_broad_gsa        niddk_cho_gsa      niddk_duerr_gsa 
# 1338                  284                  437 
# niddk_feinstein_gsa       niddk_old_gwas      niddk_rioux_gsa 
# 101                  446                  262 
# niddk_silverberg_gsa 
# 372 


# [1] "belgium"
# [1] "N samples exclude: 5799"
# [1] "EUR"
# belgium_franchimont_gsa   belgium_inf1_old_gwas   belgium_inf2_old_gwas 
# 1449                     349                     113 
# belgium_louis_gsa    belgium_vermeire_gsa         franchimont_gsa 
# 9                    3761                       5 
# vermeire_gsa 
# 4 
# [1] "nonEUR"
# belgium_franchimont_gsa   belgium_inf1_old_gwas   belgium_inf2_old_gwas 
# 31                       1                       0 
# belgium_louis_gsa    belgium_vermeire_gsa         franchimont_gsa 
# 1                      76                       0 
# vermeire_gsa 
# 0 
# [1] "No pheno discrepancies"


# [1] "prism"
# [1] "N samples exclude: 1763"

# [1] "EUR"
# helmsley_prism_gsa helmsley_xavier_prism_gsa             prism_nfe_gsa 
# 590                       292                        30 
# prism_nfe_gwas          xavier_prism_gsa          xavier_share_gsa 
# 290                       371                        28 

# [1] "nonEUR"
# helmsley_prism_gsa helmsley_xavier_prism_gsa             prism_nfe_gsa 
# 71                        28                         6 
# prism_nfe_gwas          xavier_prism_gsa          xavier_share_gsa 
# 27                        30                         0 

# [1] "N samples different pheno"
# helmsley_prism_gsa      prism_nfe_gsa     prism_nfe_gwas   xavier_prism_gsa 
# 1                  2                  2                  2 
# xavier_share_gsa 
# 3 


# [1] "sweden"
# [1] "N samples exclude: 16"
# [1] "EUR"
# swedish_uc_old_gwas 
# 16 
# [1] "nonEUR"
# swedish_uc_old_gwas 
# 0 
# [1] "No pheno discrepancies"


# [1] "cedars"
# [1] "N samples exclude: 4005"
# [1] "EUR"
# cedars_370k_old_gwas cedars_610k_old_gwas           cedars_gsa 
# 114                  339                 2654 
# cedars_omni_old_gwas         mcgovern_gsa 
# 426                    7 

# [1] "nonEUR"
# cedars_370k_old_gwas cedars_610k_old_gwas           cedars_gsa 
# 40                   23                  357 
# cedars_omni_old_gwas         mcgovern_gsa 
# 45                    0 

# [1] "N samples different pheno"
# cedars_gsa cedars_omni_old_gwas         mcgovern_gsa 
# 6                    1                    7 



# [1] "netherlands"
# [1] "N samples exclude: 693"
# [1] "EUR"
# netherlands_gsa     weersma_gsa 
# 1             672 
# [1] "nonEUR"
# netherlands_gsa     weersma_gsa 
# 0              20 
# [1] "N samples different pheno"
# netherlands_gsa     weersma_gsa 
# 1               1 


# [1] "mccauley"
# [1] "N samples exclude: 771"
# [1] "EUR"
# mccauley_gsa 
# 754 
# [1] "nonEUR"
# mccauley_gsa 
# 17 
# [1] "No pheno discrepancies"


# [1] "ccfa"
# [1] "N samples exclude: 2166"
# [1] "EUR"
# 
# ccfa_gsa 
# 1931 
# [1] "nonEUR"
# 
# ccfa_gsa 
# 235 
# [1] "No pheno discrepancies"

##################

############################


#################################################
# 25.3 - REMOVE INTER-COHORT DUPLICATED SAMPLES #
#################################################

# see script evaluate_relatedness_IIBDGC.R


path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa 
         netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa 
         belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas 
         german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas 
         cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa 
         ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa 
         helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa 
         moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa 
         palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa 
         xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_35_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_35_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_all_iibdgc_data_remove_step25 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_35_${i} | grep -E "completed"
done

for i in ${studies[@]}
do echo ${i} && wc -l ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
done

###############
# all_hce
# all_hce
# 22312 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/all_hce/all_hce_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# niddk_old_gwas
# 1017 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# australia_omniexome
# 1291 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/australia_omniexome_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# gwas1
# 4655 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/gwas1_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# gwas2
# 5130 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas2/gwas2_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# pittsburgh_gsa
# 2057 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# spain_gsa
# 3408 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/spain_gsa/spain_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# italy_gsa
# 940 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/italy_gsa/italy_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# kiel_austria_sibdcs_gsa
# 13795 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# netherlands_gsa
# 4510 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# slovenia_gsa
# 261 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# sweden_gsa
# 1374 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/sweden_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# niddk_broad_gsa
# 185 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# niddk_feinstein_gsa
# 6075 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# basque_gsa
# 1481 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/basque_gsa/basque_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# lithuania_gsa
# 2204 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lithuania_gsa/lithuania_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# belgium_louis_gsa
# 1495 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# belgium_franchimont_gsa
# 4 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# belgium_vermeire_gsa
# 128 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# prism_nfe_gsa
# 419 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# prism_nfe_gwas
# 478 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# finland_illugwas
# 423 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/finland_illugwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# german_affy6_old_gwas
# 2595 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# norway_affy6_old_gwas
# 544 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# belgium_inf1_old_gwas
# 1056 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# belgium_inf2_old_gwas
# 158 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# cedars_370k_old_gwas
# 382 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# cedars_610k_old_gwas
# 323 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# cedars_omni_old_gwas
# 468 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# swedish_uc_old_gwas
# 1242 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# mccauley_gsa
# 5 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# ccfa_gsa # cero samples
# NO SAMPLES /path/to/ibdgwas/IIBDGC/pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# cedars_gsa
# 29 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_gsa/cedars_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# bernstein_gsa
# 467 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/bernstein_gsa/bernstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# farkkila_gsa
# 68 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# franchimont_gsa
# 2593 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franchimont_gsa/franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# franke_gsa
# 767 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franke_gsa/franke_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# helmsley_prism_gsa
# 91 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# helmsley_xavier_prism_gsa
# 954 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# hyams_protect_gsa
# 410 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# lewis_sparc_gsa
# 2725 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# mccauley_new_gsa
# 1581 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# mcgovern_gsa
# 5882 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# moayyedi_imagine_gsa
# 1093 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# newberry_share_gsa
# 846 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# niddk_cho_gsa
# 1289 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# niddk_duerr_gsa
# 1216 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# niddk_rioux_gsa
# 631 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# niddk_silverberg_gsa
# 1905 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# palotie_hus_gsa
# 867 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# pekow_share_gsa
# 615 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# rioux_igenomed_gsa
# 161 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# sands_msccr_gsa
# 1162 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# stampfer_gsa
# 1464 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# vermeire_gsa
# 4636 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/vermeire_gsa/vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# weersma_gsa
# 13 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/weersma_gsa/weersma_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# xavier_prism_gsa
# 283 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam
# xavier_share_gsa
# 648 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam




##### update eur list of no-duplicates

#### /software/R-4.3.1/bin/R

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"


pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/pca_1000gp_iibdgc_pca_pcaEur_withDuplicates.tsv",sep=""),head=T,sep="\t")
dim(pca)

### some old studies keep only small N of samples - repeat with non duplicates per cohort


cohorts<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
           "kiel_austria_sibdcs_gsa","netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa",
           "niddk_feinstein_gsa","basque_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa",
           "belgium_vermeire_gsa","prism_nfe_gsa","prism_nfe_gwas","finland_illugwas","german_affy6_old_gwas",
           "norway_affy6_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
           "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas",
           "swedish_uc_old_gwas",
           "mccauley_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

length(cohorts)
# [1] 58
# "ccfa_gsa", exclude - no remainng samples that were not included in lewish_part

for (i in 1:length(cohorts)){
  a<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam",sep=""),head=F)
  a$cohort<-cohorts[i]
  if(i==1){
    fam<-a
  }else{
    fam<-rbind(a,fam)
  }
}
dim(fam)
# [1] 112811      7

dim(pca)
# [1] 149457     30

table(pca$cohort[which(pca$FID %in% fam$V1)])
# 1000GP                   all_hce       australia_omniexome 
# 0                     22312                      1291 
# basque_gsa   belgium_franchimont_gsa     belgium_inf1_old_gwas 
# 1481                         4                      1056 
# belgium_inf2_old_gwas         belgium_louis_gsa      belgium_vermeire_gsa 
# 158                      1495                       128 
# bernstein_gsa                  ccfa_gsa      cedars_370k_old_gwas 
# 467                         0                       382 
# cedars_610k_old_gwas                cedars_gsa      cedars_omni_old_gwas 
# 323                        29                       468 
# farkkila_gsa          finland_illugwas           franchimont_gsa 
# 68                       423                      2593 
# franke_gsa     german_affy6_old_gwas                     gwas1 
# 767                      2595                      4655 
# gwas2        helmsley_prism_gsa helmsley_xavier_prism_gsa 
# 5130                        91                       954 
# hyams_protect_gsa                 italy_gsa   kiel_austria_sibdcs_gsa 
# 410                       940                     13795 
# lewis_sparc_gsa             lithuania_gsa              mccauley_gsa 
# 2725                      2204                         5 
# mccauley_new_gsa              mcgovern_gsa      moayyedi_imagine_gsa 
# 1581                      5882                      1093 
# netherlands_gsa        newberry_share_gsa           niddk_broad_gsa 
# 4510                       846                       185 
# niddk_cho_gsa           niddk_duerr_gsa       niddk_feinstein_gsa 
# 1289                      1216                      6075 
# niddk_old_gwas           niddk_rioux_gsa      niddk_silverberg_gsa 
# 1017                       631                      1905 
# norway_affy6_old_gwas           palotie_hus_gsa           pekow_share_gsa 
# 544                       867                       615 
# pittsburgh_gsa             prism_nfe_gsa            prism_nfe_gwas 
# 2057                       419                       478 
# rioux_igenomed_gsa           sands_msccr_gsa              slovenia_gsa 
# 161                      1162                       261 
# spain_gsa              stampfer_gsa                sweden_gsa 
# 3408                      1464                      1374 
# swedish_uc_old_gwas              vermeire_gsa               weersma_gsa 
# 1242                      4636                        13 
# xavier_prism_gsa          xavier_share_gsa 
# 283                       648 




array<-c("illumina370","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

illumina370<-c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","swedish_uc_old_gwas","niddk_old_gwas")
affymetrix6<-c("german_affy6_old_gwas","norway_affy6_old_gwas","gwas2")
humanomniexpress<-c("australia_omniexome")
affymetrix500<-c("gwas1")
humancoreexome<-c("all_hce")
humanomni1<-c("pittsburgh_gsa")
quad610<-c("spain_gsa")
gsa<-c("italy_gsa","kiel_austria_sibdcs_gsa"
       ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
       ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
       ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa","franke_gsa","hyams_protect_gsa",
       "lewis_sparc_gsa","mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
       "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
       "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
       "xavier_share_gsa")
illuminaexome<-c("prism_nfe_gwas","helmsley_prism_gsa","helmsley_xavier_prism_gsa")

pca$array<-NA

for (i in 1:length(array)) {
  
  cohort<-get(array[i])
  print(cohort)
  
  pca$array[which(pca$cohort %in% cohort)]<-array[i]
  
}


pca_nondup_eur<-pca[which(pca$FID %in% fam$V1 & pca$inferred_population=="EUR"),]
dim(pca_nondup_eur)
# [1] [1] 106630     30

table(pca_nondup_eur$array)
# affymetrix500      affymetrix6              gsa   humancoreexome 
#          4648             8246            59509            21243 
# humanomni1 humanomniexpress      illumina370    illuminaexome 
#        2053             1238             3452             1363 
# quad610 
#    3396 

table(pca_nondup_eur$cohort)
# 1000GP                   all_hce       australia_omniexome 
# 0                     21243                      1238 
# basque_gsa   belgium_franchimont_gsa     belgium_inf1_old_gwas 
# 1470                         3                      1047 
# belgium_inf2_old_gwas         belgium_louis_gsa      belgium_vermeire_gsa 
# 157                      1484                       124 
# bernstein_gsa                  ccfa_gsa      cedars_370k_old_gwas 
# 442                         0                       339 
# cedars_610k_old_gwas                cedars_gsa      cedars_omni_old_gwas 
# 301                        26                       423 
# farkkila_gsa          finland_illugwas           franchimont_gsa 
# 68                       419                      2194 
# franke_gsa     german_affy6_old_gwas                     gwas1 
# 756                      2585                      4648 
# gwas2        helmsley_prism_gsa helmsley_xavier_prism_gsa 
# 5119                        78                       872 
# hyams_protect_gsa                 italy_gsa   kiel_austria_sibdcs_gsa 
# 318                       935                     13570 
# lewis_sparc_gsa             lithuania_gsa              mccauley_gsa 
# 2416                      2192                         4 
# mccauley_new_gsa              mcgovern_gsa      moayyedi_imagine_gsa 
# 1296                      5055                       985 
# netherlands_gsa        newberry_share_gsa           niddk_broad_gsa 
# 4288                       737                       174 
# niddk_cho_gsa           niddk_duerr_gsa       niddk_feinstein_gsa 
# 1167                      1152                      5157 
# niddk_old_gwas           niddk_rioux_gsa      niddk_silverberg_gsa 
# 1009                       591                      1682 
# norway_affy6_old_gwas           palotie_hus_gsa           pekow_share_gsa 
# 542                       851                       525 
# pittsburgh_gsa             prism_nfe_gsa            prism_nfe_gwas 
# 2053                       389                       413 
# rioux_igenomed_gsa           sands_msccr_gsa              slovenia_gsa 
# 151                       896                       261 
# spain_gsa              stampfer_gsa                sweden_gsa 
# 3396                      1439                      1342 
# swedish_uc_old_gwas              vermeire_gsa               weersma_gsa 
# 1239                      4517                        10 
# xavier_prism_gsa          xavier_share_gsa 
# 248                       594 

write.table(pca_nondup_eur[,c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eur_ancestry_samples_noDuplicates_perstudy_interstudy_2",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca_nondup_eur[which(pca_nondup_eur$pca_jewish=="Non-Jewish"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eur_nonjewish_ancestry_samples_noDuplicates_perstudy_interstudy_2",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca_nondup_eur[which(pca_nondup_eur$pca_jewish=="Jewish"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eur_jewish_ancestry_samples_noDuplicates_perstudy_interstudy_2",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")


