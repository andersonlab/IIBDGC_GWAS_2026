# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# RUN PCA IN EUROPEAN ANCESTRY SAMPLES, ESTIMATE PCS IN NON DUPLICATED EUR SUBSET, AND PROJECT ALL


############################

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=5000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G ibdgwas R

library(data.table)
library(ggplot2)
library(ggpubr)


path<-"/path/to/ibdgwas/IIBDGC/"

array<-c("affymetrix6","affymetrix500","humancoreexome")

affymetrix6<-c("gwas2")
affymetrix500<-c("gwas1")
humancoreexome<-c("all_hce")

cohorts<-c(affymetrix6,affymetrix500,humancoreexome)


length(cohorts)
# [1]


for (i in 1:length(cohorts)){
  a<-read.table(paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_withDuplicates.fam",sep=""),head=F)
  a$cohort<-cohorts[i]
  if(i==1){
    fam<-a
  }else{
    fam<-rbind(a,fam)
  }
}

dim(fam)
# [1] 36348      7

# add country of origin

fam$country<-NA
fam$country[which(fam$cohort %in% c("australia_omniexome"))]<-"Australia"
fam$country[which(fam$cohort %in% c("gwas1","gwas2","all_hce"))]<-"UK"
fam$country[which(fam$cohort %in% c("spain_gsa","basque_gsa"))]<-"Spain"
fam$country[which(fam$cohort %in% c("kiel_austria_sibdcs_gsa","german_affy6_old_gwas"))]<-"Germany"                 
fam$country[which(fam$cohort %in% c("italy_gsa"))]<-"Italy"   
fam$country[which(fam$cohort %in% c("netherlands_gsa"))]<-"Netherlands"                     
fam$country[which(fam$cohort %in% c("slovenia_gsa"))]<-"Slovenia"                     
fam$country[which(fam$cohort %in% c("sweden_gsa","swedish_uc_old_gwas"))]<-"Sweden"                    
fam$country[which(fam$cohort %in% c("finland_illugwas","palotie_hus_gsa","farkkila_gsa"))]<-"Finland"
fam$country[which(fam$cohort %in% c("norway_affy6_old_gwas"))]<-"Norway"
fam$country[which(fam$cohort %in% c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","belgium_louis_gsa","belgium_franchimont_gsa"
                                    ,"belgium_vermeire_gsa"))]<-"Belgium"
fam$country[which(fam$cohort %in% c("lithuania_gsa"))]<-"Lithuania"
fam$country[which(fam$cohort %in% c("prism_nfe_gwas","prism_nfe_gsa","cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","cedars_gsa",
                                    "niddk_broad_gsa","niddk_feinstein_gsa","niddk_old_gwas","pittsburgh_gsa","mccauley_gsa","ccfa_gsa",
                                    "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","franchimont_gsa",
                                    "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
                                    "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
                                    "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","pekow_share_gsa",
                                    "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
                                    "xavier_share_gsa"))]<-"USA"


table(fam$country,useNA="ifany")
fam$country<-as.factor(fam$country)
fam$country<-factor(fam$country, levels=c("Spain","Italy","Slovenia","Belgium","Germany","Netherlands","UK","Lithuania",
                                          "Sweden","Norway","Finland","Australia","USA"))

table(fam$cohort,fam$country,useNA="ifany")

######################################################################################################################################################



ancestry<-c("eur_all")

  
for (ii in 1:length(array)) {
    
    print(array[ii])
    pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/",array[ii],"_all_studies_merged_",ancestry,"_pruned_withDuplicates_ref_pcs_new_projection.sscore",sep=""),head=F)
    colnames(pca)<-c("FID","IID","PHENO1","ALLELE_CT","NAMED_ALLELE_DOSAGE_SUM",paste("PC",seq(1:40),sep=""))
    pca<-merge(pca,fam[,c("V1","cohort")],by.x="FID",by.y="V1",all.x=T)
    
    if (ii==1) {
        pca_all<-pca
    } else {
        pca_all<-rbind(pca_all,pca)
    }
}

table(pca_all$cohort)
# all_hce   gwas1   gwas2 
#   22588    4652    7758 

fwrite(pca_all,"/path/to/project",
col.names=T,row.names=F,quote=F,sep="\t")

## final check
tmp<-fread("/path/to/project")

tmp<-merge(tmp,pca_all,by="IID",all.x=T)
table(tmp$PC1.x==tmp$PC1.y)
#  TRUE 
# 34998 