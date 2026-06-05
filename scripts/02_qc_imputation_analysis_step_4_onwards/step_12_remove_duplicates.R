# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
##############################################
# 12.- REMOVE INTRA-STUDY DUPLICATED SAMPLES #
##############################################

# The amount of computer memory required by KING analysis is modest, at ~N ✕ M / 4 (where N is the number of samples and M is the number of SNPs)
# plus a small percentage of overhead cost. E.g., for a dataset consisting of 100,000 samples each genotyped at 1,000,000 SNPs, 
# the required memory size is ~25GB.

# max required - 
# by kiel
# (566234*14190/4)/1E+6
# [1] 2008.715

# by gwas3
# (441981*23800/4)/1E+6
# [1] 2629.787


MEM=2800
MEM=3500 (gwas3)

path_gwas=/path/to/ibdgwas/IIBDGC/
studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
  
for i in ${studies[@]}
do 
bsub -J"kg1" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_king_1_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_king_1_${i} \
"/path/to/software/./king \
-b ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy.bed \
--related --cpus 4 --prefix ${path_gwas}pre_imputation/QC/${i}/${i}_king"
done


for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_king_1_${i} | grep -E "completed"
done


  
# The first algorithm estimates pair-wise kinship coefficients (through option --kinship), and the second algorithm infers pairwise IBD 
# (identical by descent) segments (through option --ibdseg). Both algorithms can also be integrated in a single inference procedure through 
# option --related.

# Duplicated or Monozygotic twin  Dup/MZ
# Parent–offspring                PO
# Full sib                        FS
# 2nd Degree                      2nd

#kinship coefficient range
# >0.354 duplicate/MZ twin 
# [0.177, 0.354] 1st-degree 
# [0.0884, 0.177] 2nd-degree


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

for (j in 1:length(studies)) {
  
  print(studies[j])
  
  file.kin<-paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_king.kin0",sep="")
  
  if(file.exists(file.kin)) {
    
    kin<-read.table(file.kin,head=T)
    print(table(kin$InfType))
    
    print(table(cut(kin$Kinship,breaks=c(1,0.354,0.177))))
    
    dup<-kin[which(kin$Kinship>0.354 & kin$InfType=="Dup/MZ"),c("FID1","FID2","Kinship")]
    
    if(nrow(dup)>0) {
      dup_ids<-c(as.character(dup$FID1),as.character(dup$FID2))
      length(dup_ids)
      table(length(dup_ids)==(nrow(dup)*2))
      
      dup_ids<-dup_ids[!duplicated(dup_ids)]
      length(dup_ids)
      
      # CHECK ALSO THAT BOTH SAMPLES HAVE SAME PHENOTYPE AND GENDER, AND REMOVE SAMPLES LOWER CALL RATES
      
      fam<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy.fam",sep=""),head=F)
      colnames(fam)[5:6]<-c("sex","pheno")
      sample_miss<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy.imiss",sep=""),head=T)
      
      all<-merge(fam[,c(1,5,6)],sample_miss[,c("FID","F_MISS")],by.x="V1",by.y="FID",all.x=T,sort=F)
      
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
            
            print("Samples with different sex/pheno:")
            print(data)
            
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
      print(dim(data_remove))
      
      data_remove<-data_remove[,c(1,1)]
      colnames(data_remove)<-c("FID","IID")
      write.table(data_remove,paste(path,"pre_imputation/QC/",studies[j],"/list_duplicated_samples",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

      
      #### create a file to report inconsitend data:
      
      if(exists("data_inconsist")){
        
        data_inconsist<-data_inconsist[!duplicated(data_inconsist$V1),]
        print(dim(data_inconsist))

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
        write.table(data_inconsist,paste(path,"pre_imputation/QC/",studies[j],"/_list_duplicated_samples_inconsistent_phenotype",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
      }
    } 
    else {
      data_remove<-as.data.frame(matrix(ncol=2,nrow=0))
      colnames(data_remove)<-c("FID","IID")
      write.table(data_remove,paste(path,"pre_imputation/QC/",studies[j],"/list_duplicated_samples",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    }
  } else {
    data_remove<-as.data.frame(matrix(ncol=2,nrow=0))
    colnames(data_remove)<-c("FID","IID")
    write.table(data_remove,paste(path,"pre_imputation/QC/",studies[j],"/list_duplicated_samples",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  }
  
  rm(list=ls()[!ls() %in% c("path","j","studies")])
}


######

# [1] "all_hce"
# 
# 2nd    3rd Dup/MZ     FS     PO 
# 128      1   1400    190   1202 
# 
# (0.177,0.354]     (0.354,1] 
# 1292          1434 
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 295082_A09_usgwas5494223   1     1 0.0001559
# 2  351593_F04_IBD_R5792538   1     2 0.0001834
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 295082_A09_usgwas5494223   1     1 0.0001559
# 2  351593_F04_IBD_R5792538   1     2 0.0001834
# 
# [1] 1276    4 # OK
# [1] 2 5


# [1] "niddk_old_gwas"
# Dup/MZ     FS     PO 
# 8      8      6 
# (0.177,0.354]     (0.354,1] 
# 14             8 
# [1] 8 4

# [1] "australia_omniexome"
# Dup/MZ     FS     PO 
# 8      1      4 
# (0.177,0.354]     (0.354,1] 
# 5             8 
# [1] 8 4


# [1] "gwas1"
# Dup/MZ 
# 1 
# (0.177,0.354]     (0.354,1] 
# 0             1 
# [1] 1 4

# [1] "gwas2"
# Dup/MZ     FS 
# 1      1 
# (0.177,0.354]     (0.354,1] 
# 1             1 
# [1] 1 4

# [1] "pittsburgh_gsa"
# FS PO 
# 18  3 
# (0.177,0.354]     (0.354,1] 
# 21             0 

# [1] "spain_gsa"

# [1] "italy_gsa"
# Dup/MZ     FS     PO 
# 4      2      8 
# (0.177,0.354]     (0.354,1] 
# 10             4 
# [1] 4 4

# [1] "kiel_austria_sibdcs_gsa"
# Dup/MZ     FS     PO 
# 539    176    303 
# (0.177,0.354]     (0.354,1] 
# 479           539 
#
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.004869
# 2 sample_id   1     2 0.001110
# 3 sample_id   1     2 0.005532
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.014090
# 2 sample_id   1     2 0.008092
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0021430
# 2 sample_id   2     2 0.0005919
# 3 sample_id   2     2 0.0035550
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0004700
# 2 sample_id   1     2 0.0007262
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.001709
# 2 sample_id   2     2 0.007468
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.005564
# 2 sample_id   2     2 0.002368
# 3 sample_id   2     2 0.004189
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0094010
# 2 sample_id   2     2 0.0067760
# 3 sample_id   2     2 0.0007898
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.009334
# 2 sample_id   1     2 0.001631
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0075120
# 2 sample_id   2     2 0.0079580
# 3 sample_id   2     2 0.0008498
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0012920
# 2 sample_id   2     2 0.0006608
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0008004
# 2 sample_id   1     2 0.0020550
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.001979
# 2 sample_id   1     2 0.008428
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0005141
# 2 sample_id   2     2 0.0016060
# 3 sample_id   2     2 0.0011470
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.003154
# 2 sample_id   1     2 0.004657
# 3 sample_id   1     2 0.001615
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.001239
# 2 sample_id   2     2 0.004514
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.001270
# 2 sample_id   2     2 0.005083
# 3 sample_id   2     2 0.007489
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.001150
# 2 sample_id   2     2 0.002991
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.003740
# 2 sample_id   1     2 0.001475
# 3 sample_id   1     2 0.005749
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0004682
# 2 sample_id   2     2 0.0010480
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.001686
# 2 sample_id   2     2 0.002656
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.002209
# 2 sample_id   2     2 0.008228
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0040950
# 2 sample_id   1     2 0.0008693
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.006562
# 2 sample_id   1     2 0.002046
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.004435
# 2 sample_id   2     2 0.004488
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.002615
# 2 sample_id   2     2 0.002807
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.002304
# 2 sample_id   2     2 0.002578
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.010480
# 2 sample_id   1     2 0.004233
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0028160
# 2 sample_id   2     2 0.0005654
# 3 sample_id   2     2 0.0011730
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.001866
# 2 sample_id   2     2 0.007654
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.001270
# 2 sample_id   2     2 0.005083
# 3 sample_id   2     2 0.007489
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.010480
# 2 sample_id   1     2 0.004233
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.005564
# 2 sample_id   2     2 0.002368
# 3 sample_id   2     2 0.004189
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.003154
# 2 sample_id   1     2 0.004657
# 3 sample_id   1     2 0.001615
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.002304
# 2 sample_id   2     2 0.002578
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.001709
# 2 sample_id   2     2 0.007468
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.002615
# 2 sample_id   2     2 0.002807
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0040950
# 2 sample_id   1     2 0.0008693
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0005141
# 2 sample_id   2     2 0.0016060
# 3 sample_id   2     2 0.0011470
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     2 0.001375
# 2 sample_id   1     1 0.001770
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0028160
# 2 sample_id   2     2 0.0005654
# 3 sample_id   2     2 0.0011730
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0075120
# 2 sample_id   2     2 0.0079580
# 3 sample_id   2     2 0.0008498
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.009334
# 2 sample_id   1     2 0.001631
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.003740
# 2 sample_id   1     2 0.001475
# 3 sample_id   1     2 0.005749
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0004682
# 2 sample_id   2     2 0.0010480
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0094010
# 2 sample_id   2     2 0.0067760
# 3 sample_id   2     2 0.0007898
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.001239
# 2 sample_id   2     2 0.004514
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0004700
# 2 sample_id   1     2 0.0007262
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.001150
# 2 sample_id   2     2 0.002991
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0021430
# 2 sample_id   2     2 0.0005919
# 3 sample_id   2     2 0.0035550
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.004869
# 2 sample_id   1     2 0.001110
# 3 sample_id   1     2 0.005532
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0021430
# 2 sample_id   2     2 0.0005919
# 3 sample_id   2     2 0.0035550
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0008004
# 2 sample_id   1     2 0.0020550
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.001866
# 2 sample_id   2     2 0.007654
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0094010
# 2 sample_id   2     2 0.0067760
# 3 sample_id   2     2 0.0007898
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     2 0.001375
# 2 sample_id   1     1 0.001770
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0012920
# 2 sample_id   2     2 0.0006608
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.005564
# 2 sample_id   2     2 0.002368
# 3 sample_id   2     2 0.004189
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.003154
# 2 sample_id   1     2 0.004657
# 3 sample_id   1     2 0.001615
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.004869
# 2 sample_id   1     2 0.001110
# 3 sample_id   1     2 0.005532
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.001979
# 2 sample_id   1     2 0.008428
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.014090
# 2 sample_id   1     2 0.008092
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.002209
# 2 sample_id   2     2 0.008228
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.004435
# 2 sample_id   2     2 0.004488
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0005141
# 2 sample_id   2     2 0.0016060
# 3 sample_id   2     2 0.0011470
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.001686
# 2 sample_id   2     2 0.002656
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0075120
# 2 sample_id   2     2 0.0079580
# 3 sample_id   2     2 0.0008498
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0028160
# 2 sample_id   2     2 0.0005654
# 3 sample_id   2     2 0.0011730
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.006562
# 2 sample_id   1     2 0.002046
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.003740
# 2 sample_id   1     2 0.001475
# 3 sample_id   1     2 0.005749
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     1 0.001270
# 2 sample_id   2     2 0.005083
# 3 sample_id   2     2 0.007489
# 
# [1] 496   4
# [1] 70   5


# [1] "netherlands_gsa"
# Dup/MZ     FS     PO 
# 20     35     36 
# (0.177,0.354]     (0.354,1] 
# 71            20 
##
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     2 0.0005544
# 2 sample_id   1     1 0.0008346
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     2 0.0010640
# 2 sample_id   2     1 0.0009724
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     2 0.0010640
# 2 sample_id   2     1 0.0009724
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     2 0.0005544
# 2 sample_id   1     1 0.0008346
# 
# [1] 22  4
# [1] 4 5  

# [1] "slovenia_gsa"
# Dup/MZ     PO 
# 2      1 
# (0.177,0.354]     (0.354,1] 
# 1             2 
# [1] 2 4

# [1] "sweden_gsa"
# Dup/MZ     FS     PO 
# 21      4     14 
# (0.177,0.354]     (0.354,1] 
# 18            21 
#
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     2 0.004135
# 2 sample_id   1     1 0.001151
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     2 0.004135
# 2 sample_id   1     1 0.001151
# 
# [1] 22  4
# [1] 2 5

# [1] "niddk_broad_gsa"
# Dup/MZ     FS     PO 
# 53     29     41 
# 
# (0.177,0.354]     (0.354,1] 
# 70            53 
#
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     2 0.001711
# 2 sample_id   2     1 0.002882
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   2     2 0.001711
# 2 sample_id   2     1 0.002882
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     2 0.001127
# 2 sample_id   1     1 0.001869
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     2 0.001127
# 2 sample_id   1     1 0.001869
# 
# [1] 54  4
# [1] 4 5

# [1] "niddk_feinstein_gsa"
# Dup/MZ     FS     PO 
# 208     92    116 
# 
# (0.177,0.354]     (0.354,1] 
# 208           208 
##
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0001605
# 2 sample_id   2     2 0.0003408
# 3 sample_id   2     1 0.0017920
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0001605
# 2 sample_id   2     2 0.0003408
# 3 sample_id   2     1 0.0017920
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0003855
# 2 sample_id   1     2 0.0017970
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0003855
# 2 sample_id   1     2 0.0017970
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0001605
# 2 sample_id   2     2 0.0003408
# 3 sample_id   2     1 0.0017920
# 
# [1] 200   4
# [1] 5  5


# [1] "basque_gsa"
# Dup/MZ     FS     PO 
# 17     16     10 
# (0.177,0.354]     (0.354,1] 
# 26            17 
# [1] 17  4

# [1] "lithuania_gsa"
# Dup/MZ     FS     PO 
# 21      8      9 
# (0.177,0.354]     (0.354,1] 
# 17            21 
# [1] 21  4

# [1] "belgium_louis_gsa"
# 2nd Dup/MZ     FS     PO 
# 1      6     10    119 
# (0.177,0.354]     (0.354,1] 
# 130             6 
# [1] 6 4


# [1] "belgium_franchimont_gsa"
# Dup/MZ     FS     PO 
# 7     10     16 
# (0.177,0.354]     (0.354,1] 
# 26             7 
# [1] 7 4

# [1] "belgium_vermeire_gsa"
# 2nd Dup/MZ     FS     PO 
# 1      9     12     19 
# (0.177,0.354]     (0.354,1] 
# 31            10 
# [1] 9 4

# [1] "prism_nfe_gsa"

# [1] "prism_nfe_gwas"
# Dup/MZ     FS     PO 
# 1      6     15 
# (0.177,0.354]     (0.354,1] 
# 15             1 
# [1] 1 4

# [1] "finland_illugwas"
# Dup/MZ 
# 2 
# (0.177,0.354]     (0.354,1] 
# 0             2 
# [1] 2 4

[1] "german_affy6_old_gwas"
# Dup/MZ     FS     PO 
# 56     12     66 
# (0.177,0.354]     (0.354,1] 
# 78            56 
# [1] 36  4

# [1] "norway_affy6_old_gwas"
# [1] "belgium_inf1_old_gwas"
# [1] "belgium_inf2_old_gwas"
# 
# [1] "cedars_370k_old_gwas"
# FS PO 
# 4  3 
# (0.177,0.354]     (0.354,1] 
# 7             0 

# [1] "cedars_610k_old_gwas"
# FS PO 
# 6  6 
# (0.177,0.354]     (0.354,1] 
# 12             0 

# [1] "cedars_omni_old_gwas"
# Dup/MZ     FS 
# 2      1 
# 
# (0.177,0.354]     (0.354,1] 
# 1             2 
# [1] 2 4

# [1] "swedish_uc_old_gwas"

# [1] "mccauley_gsa"
# Dup/MZ     FS     PO 
# 1      5      4 
# (0.177,0.354]     (0.354,1] 
# 9             1 
# [1] 1 4

# [1] "ccfa_gsa"
# Dup/MZ     FS     PO 
# 1      4      4 
# (0.177,0.354]     (0.354,1] 
# 8             1 
# [1] 1 4

# [1] "cedars_gsa"
# Dup/MZ     FS     PO 
# 2     32     99 
# (0.177,0.354]     (0.354,1] 
# 131             2 
# [1] 2 4

# [1] "bernstein_gsa"
# Dup/MZ     FS     PO 
# 31      2      6 
# (0.177,0.354]     (0.354,1] 
# 8            31 
# [1] 30  4

# [1] "farkkila_gsa"

# [1] "franchimont_gsa"
# Dup/MZ     FS     PO 
# 61     21     29 
# (0.177,0.354]     (0.354,1] 
# 50            61 
# [1] 61  4

# [1] "franke_gsa"
# Dup/MZ     FS 
# 7      1 
# (0.177,0.354]     (0.354,1] 
# 1             7 
# [1] 7 4

# [1] "helmsley_prism_gsa"
# Dup/MZ     FS     PO 
# 1      3     10 
# (0.177,0.354]     (0.354,1] 
# 13             1 
# [1] 1 4

# [1] "helmsley_xavier_prism_gsa"
# 2nd  FS  PO 
# 1   8   8 
# (0.177,0.354]     (0.354,1] 
# 17             0 

# [1] "hyams_protect_gsa"

# [1] "lewis_sparc_gsa"
# Dup/MZ     FS     PO 
# 2      4      5 
# (0.177,0.354]     (0.354,1] 
# 9             2 
# [1] 2 4

# [1] "mccauley_new_gsa"
# Dup/MZ     FS     PO 
# 7      7     23 
# (0.177,0.354]     (0.354,1] 
# 30             7 
# [1] 7 4

# [1] "mcgovern_gsa"
# Dup/MZ     FS     PO 
# 41    292    513 
# (0.177,0.354]     (0.354,1] 
# 805            41 
#
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0006269
# 2 sample_id   1     2 0.0003237
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0008321
# 2 sample_id   2     2 0.0003758
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0005005
# 2 sample_id   1     2 0.0004390
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0010960
# 2 sample_id   1     2 0.0005305
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0038500
# 2 sample_id   1     2 0.0003537
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0005921
# 2 sample_id   1     2 0.0003884
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0006063
# 2 sample_id   1     2 0.0005353
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.001587
# 2 sample_id   1     2 0.002016
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0008226
# 2 sample_id   2     2 0.0003837
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     2 0.0013340
# 2 sample_id   1     1 0.0006348
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0006521
# 2 sample_id   1     2 0.0009537
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0004595
# 2 sample_id   1     2 0.0005511
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0013670
# 2 sample_id   2     2 0.0005053
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.001985
# 2 sample_id   1     2 0.001669
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0010550
# 2 sample_id   1     2 0.0007832
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0004926
# 2 sample_id   1     2 0.0006205
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0004753
# 2 sample_id   2     2 0.0008684
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0008211
# 2 sample_id   1     2 0.0002574
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0006521
# 2 sample_id   1     2 0.0009537
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0005921
# 2 sample_id   2     2 0.0072110
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0005921
# 2 sample_id   1     2 0.0003884
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0005921
# 2 sample_id   2     2 0.0072110
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.001587
# 2 sample_id   1     2 0.002016
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0004753
# 2 sample_id   2     2 0.0008684
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0004926
# 2 sample_id   1     2 0.0006205
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.001985
# 2 sample_id   1     2 0.001669
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     2 0.0013340
# 2 sample_id   1     1 0.0006348
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0013820
# 2 sample_id   2     2 0.0004626
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0013670
# 2 sample_id   2     2 0.0005053
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0008321
# 2 sample_id   2     2 0.0003758
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0004595
# 2 sample_id   1     2 0.0005511
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0008226
# 2 sample_id   2     2 0.0003837
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0038500
# 2 sample_id   1     2 0.0003537
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0027570
# 2 sample_id   1     2 0.0004121
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0008211
# 2 sample_id   1     2 0.0002574
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0006063
# 2 sample_id   1     2 0.0005353
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0006269
# 2 sample_id   1     2 0.0003237
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0010960
# 2 sample_id   1     2 0.0005305
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0010550
# 2 sample_id   1     2 0.0007832
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0005005
# 2 sample_id   1     2 0.0004390
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0027570
# 2 sample_id   1     2 0.0004121
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0013820
# 2 sample_id   2     2 0.0004626
#
# [1] 62  4
# [1] 42  5

# [1] "moayyedi_imagine_gsa"
# FS PO 
# 3 18 
# (0.177,0.354]     (0.354,1] 
# 21             0 

# [1] "newberry_share_gsa"
# Dup/MZ     FS     PO 
# 2      3      2 
# (0.177,0.354]     (0.354,1] 
# 5             2 
# [1] 2 4

# [1] "niddk_cho_gsa"
# Dup/MZ     FS     PO 
# 16      8      3 
# (0.177,0.354]     (0.354,1] 
# 11            16 
#
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     2 0.0010310
# 2 sample_id   1     1 0.0009469
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.000693
# 2 sample_id   1     2 0.001171
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     2 0.0010310
# 2 sample_id   1     1 0.0009469
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0003061
# 2 sample_id   2     2 0.0006021
# [1] "Samples with different sex/pheno:"
# V1 sex pheno   F_MISS
# 1 sample_id   1     1 0.000693
# 2 sample_id   1     2 0.001171
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0003061
# 2 sample_id   2     2 0.0006021
#
# [1] 18  4
# [1] 6  5


# [1] "niddk_duerr_gsa"
# Dup/MZ     FS     PO 
# 121     29     29 
# (0.177,0.354]     (0.354,1] 
# 58           121 
###
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0003039
# 2 sample_id   1     2 0.0004533
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0004602
# 2 sample_id   2     2 0.0006010
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     2 0.0006164
# 2 sample_id   2     1 0.0005185
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     2 0.0005151
# 2 sample_id   2     1 0.0003674
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     2 0.0006164
# 2 sample_id   2     1 0.0005185
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0008001
# 2 sample_id   1     2 0.0007538
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     2 0.0005151
# 2 sample_id   2     1 0.0003674
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0004739
# 2 sample_id   1     2 0.0007898
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0004602
# 2 sample_id   2     2 0.0006010
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0003039
# 2 sample_id   1     2 0.0004533
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0004739
# 2 sample_id   1     2 0.0007898
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0008001
# 2 sample_id   1     2 0.0007538
# 
# [1] 127   4
# [1] 12  5


# [1] "niddk_rioux_gsa"
# Dup/MZ     FS     PO 
# 5      7     16 
# (0.177,0.354]     (0.354,1] 
# 23             5 
###
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0003786
# 2 sample_id   2     2 0.0005869
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0003786
# 2 sample_id   2     2 0.0005869
# 
# [1] 6 4
# [1] 2 5


# [1] "niddk_silverberg_gsa"
# Dup/MZ     FS     PO 
# 40     26     51 
# (0.177,0.354]     (0.354,1] 
# 77            40 
#
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     2 0.0003871
# 2 sample_id   2     1 0.0003495
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     2 0.0003871
# 2 sample_id   2     1 0.0003495
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0007023
# 2 sample_id   2     2 0.0006827
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0040410
# 2 sample_id   2     2 0.0005537
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0007023
# 2 sample_id   2     2 0.0006827
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0005276
# 2 sample_id   1     2 0.0008591
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   1     1 0.0005276
# 2 sample_id   1     2 0.0008591
# [1] "Samples with different sex/pheno:"
# V1 sex pheno    F_MISS
# 1 sample_id   2     1 0.0040410
# 2 sample_id   2     2 0.0005537
# 
# [1] 41  4
# [1] 8  5

# [1] "palotie_hus_gsa"
# FS PO 
# 2  1 
# (0.177,0.354]     (0.354,1] 
# 3             0 

# [1] "pekow_share_gsa"
# Dup/MZ     FS     PO 
# 1      1      1 
# (0.177,0.354]     (0.354,1] 
# 2             1 
# [1] 1 4

# [1] "rioux_igenomed_gsa"
# Dup/MZ 
# 6 
# (0.177,0.354]     (0.354,1] 
# 0             6 
# [1] 6 4

# [1] "sands_msccr_gsa"
# Dup/MZ     FS     PO 
# 10      2      8 
# (0.177,0.354]     (0.354,1] 
# 10            10 
# [1] 8 4

# [1] "stampfer_gsa"
# PO 
# 1 
# (0.177,0.354]     (0.354,1] 
# 1             0 

# [1] "vermeire_gsa"
# Dup/MZ     FS     PO 
# 14    105     86 
# (0.177,0.354]     (0.354,1] 
# 191            14 
# [1] 14  4

# [1] "weersma_gsa"
# FS PO 
# 3  7 
# (0.177,0.354]     (0.354,1] 
# 10             0 

# [1] "xavier_prism_gsa"
# Dup/MZ     PO 
# 2      1 
# (0.177,0.354]     (0.354,1] 
# 1             2 
# [1] 2 4

# [1] "xavier_share_gsa"
# Dup/MZ     FS     PO 
# 3      1      2 
# (0.177,0.354]     (0.354,1] 
# 3             3 
# [1] 3 4

######

######

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_33_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_33_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/${i}/list_duplicated_samples \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_33_${i} | grep -E "completed"
done

############################

############
# 12.2 - REMOVE INTRA-STUDY DUPLICATED SAMPLES #
############

# list variants set up in script - create_list_common_variants_among_cohorts
# evaluating relatedness per country / cohorts in:

# evaluate_relatedness_belgium.R
# evaluate_relatedness_cedars.R
# evaluate_relatedness_german.R
# evaluate_relatedness_niddk.R
# evaluate_relatedness_spain.R
# evaluate_relatedness_SWEDEN.R
# evaluate_relatedness_UKIBDGC.R

# updated in 2022
path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_34_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_34_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy"
done


# belgium
studies=(belgium_inf1_old_gwas belgium_inf2_old_gwas belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa)

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_34_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_34_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_belgium_data_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy"
done

# ukibdgc
studies=(gwas1 gwas2 all_hce)

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_34_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_34_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_ukiibdgc_data_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy"
done

studies=(kiel_austria_sibdcs_gsa german_illu550_old_gwas german_affy6_old_gwas)

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_34_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_34_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_german_affy_kiel_data_remove_2 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy"
done

studies=(niddk_broad_gsa niddk_feinstein_gsa niddk_old_gwas)
for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_34_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_34_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_niddk_data_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy"
done



studies=(spain_gsa basque_gsa) # no overalp
for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_34_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_34_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy"
done


studies=(prism_nfe_gsa prism_nfe_gwas)
for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_34_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_34_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_prism_data_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy"
done

studies=(swedish_uc_old_gwas sweden_gsa)
for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_34_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_34_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_sweden_data_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy"
done

studies=(cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas cedars_gsa)
for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_34_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_34_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_cedars_data_remove.tsv \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy"
done



studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_34_${i} | grep -E "completed"
done


############################

#################################################
# 12.3 - REMOVE INTER-COHORT DUPLICATED SAMPLES #
#################################################

# see script evaluate_relatedness_IIBDGC.R


path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_35_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_35_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_iibdgc_data_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy"
done


for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_35_${i} | grep -E "completed"
done


for i in ${studies[@]}
do echo ${i} && wc -l ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
done

###############
# all_hce
# 22426 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/all_hce/all_hce_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# niddk_old_gwas
# 1346 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# australia_omniexome
# 1298 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/australia_omniexome_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# gwas1
# 4676 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/gwas1_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# gwas2
# 5132 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas2/gwas2_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# pittsburgh_gsa
# 2030 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# spain_gsa
# 3443 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/spain_gsa/spain_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# italy_gsa
# 949 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/italy_gsa/italy_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# kiel_austria_sibdcs_gsa
# 13870 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# netherlands_gsa
# 4555 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# slovenia_gsa
# 262 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# sweden_gsa
# 1379 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/sweden_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# niddk_broad_gsa
# 5122 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# niddk_feinstein_gsa
# 7960 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# basque_gsa
# 1491 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/basque_gsa/basque_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# lithuania_gsa
# 2210 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lithuania_gsa/lithuania_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# belgium_louis_gsa
# 1506 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# belgium_franchimont_gsa
# 1488 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# belgium_vermeire_gsa
# 3979 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# prism_nfe_gsa
# 464 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# prism_nfe_gwas
# 735 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# finland_illugwas
# 444 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/finland_illugwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# german_affy6_old_gwas
# 2631 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# norway_affy6_old_gwas
# 550 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# belgium_inf1_old_gwas
# 1067 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# belgium_inf2_old_gwas
# 161 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# cedars_370k_old_gwas
# 453 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# cedars_610k_old_gwas
# 467 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# cedars_omni_old_gwas
# 673 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# swedish_uc_old_gwas
# 1248 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# mccauley_gsa
# 769 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# ccfa_gsa
# 2064 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# cedars_gsa
# 1803 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_gsa/cedars_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# bernstein_gsa
# 481 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/bernstein_gsa/bernstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# farkkila_gsa
# 68 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# franchimont_gsa
# 2727 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franchimont_gsa/franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# franke_gsa
# 860 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franke_gsa/franke_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# helmsley_prism_gsa
# 760 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# helmsley_xavier_prism_gsa
# 1293 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# hyams_protect_gsa
# 418 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# lewis_sparc_gsa
# 2857 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# mccauley_new_gsa
# 1612 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# mcgovern_gsa
# 5960 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# moayyedi_imagine_gsa
# 1128 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# newberry_share_gsa
# 861 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# niddk_cho_gsa
# 1746 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# niddk_duerr_gsa
# 1817 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# niddk_rioux_gsa
# 913 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# niddk_silverberg_gsa
# 2330 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# palotie_hus_gsa
# 878 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# pekow_share_gsa
# 633 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# rioux_igenomed_gsa
# 176 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# sands_msccr_gsa
# 1422 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# stampfer_gsa
# 1468 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# vermeire_gsa
# 4683 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/vermeire_gsa/vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# weersma_gsa
# 709 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/weersma_gsa/weersma_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# xavier_prism_gsa
# 687 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam
# xavier_share_gsa
# 693 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy.fam

#################

############################

#################################################
# 12.3 - REMOVE INTER-COHORT DUPLICATED SAMPLES #
#################################################

# see script evaluate_relatedness_samples_old_broad_batch_vs_new_batch_2022.R


path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_36_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_36_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy \
--allow-no-sex \
--remove ${path_gwas}pre_imputation/QC/relatedness/list_iibdgc_old_new_broad_batch_data_remove \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_36_${i} | grep -E "completed"
done

for i in ${studies[@]}
do echo ${i} && wc -l ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
done


#################
# all_hce
# 22426 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/all_hce/all_hce_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# niddk_old_gwas
# 1346 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# australia_omniexome
# 1298 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/australia_omniexome_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# gwas1
# 4676 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/gwas1_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# gwas2
# 5132 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas2/gwas2_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# pittsburgh_gsa
# 2030 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# spain_gsa
# 3443 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/spain_gsa/spain_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# italy_gsa
# 949 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/italy_gsa/italy_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# kiel_austria_sibdcs_gsa
# 13777 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# netherlands_gsa
# 4554 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# slovenia_gsa
# 262 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# sweden_gsa
# 1379 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/sweden_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# niddk_broad_gsa
# 185 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# niddk_feinstein_gsa
# 5911 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# basque_gsa
# 1491 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/basque_gsa/basque_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# lithuania_gsa
# 2210 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lithuania_gsa/lithuania_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# belgium_louis_gsa
# 1502 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# belgium_franchimont_gsa
# 2 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# belgium_vermeire_gsa
# 119 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# prism_nfe_gsa
# 42 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# prism_nfe_gwas
# 569 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# finland_illugwas
# 425 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/finland_illugwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# german_affy6_old_gwas
# 2631 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# norway_affy6_old_gwas
# 550 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# belgium_inf1_old_gwas
# 1067 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# belgium_inf2_old_gwas
# 161 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# cedars_370k_old_gwas
# 453 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# cedars_610k_old_gwas
# 467 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# cedars_omni_old_gwas
# 673 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# swedish_uc_old_gwas
# 1248 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# mccauley_gsa
# 5 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# ccfa_gsa
# 1 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# cedars_gsa
# 17 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_gsa/cedars_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# bernstein_gsa
# 464 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/bernstein_gsa/bernstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# farkkila_gsa
# 68 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# franchimont_gsa
# 2722 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franchimont_gsa/franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# franke_gsa
# 857 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franke_gsa/franke_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# helmsley_prism_gsa
# wc: /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam: No such file or directory
# helmsley_xavier_prism_gsa
# 936 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# hyams_protect_gsa
# 418 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# lewis_sparc_gsa
# 2780 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# mccauley_new_gsa
# 1602 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# mcgovern_gsa
# 5898 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# moayyedi_imagine_gsa
# 1107 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# newberry_share_gsa
# 860 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# niddk_cho_gsa
# 1448 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# niddk_duerr_gsa
# 1320 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# niddk_rioux_gsa
# 640 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# niddk_silverberg_gsa
# 1952 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# palotie_hus_gsa
# 877 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# pekow_share_gsa
# 600 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# rioux_igenomed_gsa
# 167 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# sands_msccr_gsa
# 1080 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# stampfer_gsa
# 1467 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# vermeire_gsa
# 4678 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/vermeire_gsa/vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# weersma_gsa
# 13 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/weersma_gsa/weersma_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# xavier_prism_gsa
# 663 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam
# xavier_share_gsa
# 651 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_2.fam



