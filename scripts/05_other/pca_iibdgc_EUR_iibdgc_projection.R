# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# RUN PCA IN EUROEPAN ANCESTRY SAMPLES, ESTIMATE PCS IN NON DUPLICATED EUR SUBSET, AND PROJECT ALL:


#########################################################
## 1.1 - ESTIMATE PC, EXPORT PCS, AND ALLELE FREQUENCIES:

# The following command exports PCs to project onto, along with the allele frequencies needed to calibrate the 'variance-standardize' operation:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_noHighLD_pruned_noassoclogP4 \
--keep ${path}pre_imputation/QC/pca_1000gp/list_eur_ancestry_samples_noDuplicates \
--make-bed --out ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_noHighLD_pruned_noassoclogP4_EURonly_noDuplicates
# 2323 variants and 102134 people pass filters and QC.

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_withDuplicates_subset_iibdgc \
--keep ${path}pre_imputation/QC/pca_1000gp/list_eur_ancestry_samples_withDuplicates \
--make-bed --out ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_noHighLD_pruned_noassoclogP4_EURonly_withDuplicates
# 2323 variants and 113112 people pass filters and QC.


path=/path/to/ibdgwas/IIBDGC/
MEM=5000 # next time only this
  
bsub -J"bcftools" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 10 \
-e ${path}pre_imputation/QC/pca_1000gp/logs/stderr_pca_eur_nodup \
-o ${path}pre_imputation/QC/pca_1000gp/logs/stdout_pca_eur_nodup \
"/path/to/software/./plink2  \
--bfile ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_noHighLD_pruned_noassoclogP4_EURonly_noDuplicates \
--freq counts \
--pca approx allele-wts \
--out ${path}pre_imputation/QC/pca_1000gp/iibdgc_EUR_nodup_ref_pcs \
--threads 8 --allow-no-sex --memory $MEM"
# Job <993962> is submitted to queue <normal>.


#########################################################
## 1.2 - PROJECT ONTO THOSE PCS ALL IIBDGC

path=/path/to/ibdgwas/IIBDGC/
MEM=1000 # next time only this

bsub -J"bcftools" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 10 \
-e ${path}pre_imputation/QC/pca_1000gp/logs/stdout_pca_eur_nodup_project_withDuplicates \
-o ${path}pre_imputation/QC/pca_1000gp/logs/stdout_pca_eur_nodup_project_withDuplicates \
"/path/to/software/./plink2 \
--bfile ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_noHighLD_pruned_noassoclogP4_EURonly_withDuplicates \
--read-freq ${path}pre_imputation/QC/pca_1000gp/iibdgc_EUR_nodup_ref_pcs.acount \
--score ${path}pre_imputation/QC/pca_1000gp/iibdgc_EUR_nodup_ref_pcs.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
--score-col-nums 6-15 \
--out ${path}pre_imputation/QC/pca_1000gp/iibdgc_EUR_withDuplicates_ref_pcs_new_projection"
# Job <995336> is submitted to queue <normal>.



###############################
  
### /software/R-4.3.1/bin/R


library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(viridis)

path<-"/path/to/ibdgwas/IIBDGC/"

# no german_illu
cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"chop_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas"
           ,"swedish_uc_old_gwas","mccauley_gsa","ccfa_gsa","cedars_gsa")

length(cohorts)
# [1] 34


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
# [1] 117923      7


pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_EUR_withDuplicates_ref_pcs_new_projection.sscore",sep=""),head=F)
colnames(pca)<-c("FID","IID","PHENO1","ALLELE_CT","NAMED_ALLELE_DOSAGE_SUM","PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG"
                 ,"PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG")
dim(pca)
# [1]  113112     15
colnames(pca)<-gsub("_AVG","",colnames(pca))

pca<-merge(pca,fam[,c("V1","cohort")],by.x="FID",by.y="V1",all.x=T)
dim(pca)
# [1] 113112    13

table(pca$cohort)
# all_hce     australia_omniexome              basque_gsa 
# 22692                    1249                    1498 
# belgium_franchimont_gsa   belgium_inf1_old_gwas   belgium_inf2_old_gwas 
# 1464                    1403                     269 
# belgium_louis_gsa    belgium_vermeire_gsa                ccfa_gsa 
# 1489                    3905                    1930 
# cedars_370k_old_gwas    cedars_610k_old_gwas              cedars_gsa 
# 508                     830                    2708 
# cedars_omni_old_gwas           chop_old_gwas        finland_illugwas 
# 1107                    8467                     443 
# german_affy6_old_gwas                   gwas1                   gwas2 
# 2804                    4663                    7751 
# italy_gsa kiel_austria_sibdcs_gsa           lithuania_gsa 
# 947                   14126                    2195 
# mccauley_gsa         netherlands_gsa         niddk_broad_gsa 
# 758                    4347                    5146 
# niddk_feinstein_gsa          niddk_old_gwas   norway_affy6_old_gwas 
# 6956                    2731                     545 
# pittsburgh_gsa           prism_nfe_gsa          prism_nfe_gwas 
# 2713                     429                     719 
# slovenia_gsa               spain_gsa              sweden_gsa 
# 264                    3438                    1324 
# swedish_uc_old_gwas 
# 1260 


# TOTAL VARIANCE EXPLAINED BY EACH PC:


eigenval<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_EUR_nodup_ref_pcs.eigenval",sep=""),head=F)

eigenval$var_exp<-NA
for (i in 1:nrow(eigenval)){
  eigenval$var_exp[i]<-(eigenval$V1[i] / sum(eigenval$V1))*100
}

eigenval
#          V1   var_exp
# 1  214.6160 22.699497
# 2  115.4610 12.212074
# 3   92.6266  9.796927
# 4   78.9639  8.351851
# 5   76.2488  8.064680
# 6   76.0964  8.048561
# 7   74.5226  7.882103
# 8   73.9027  7.816538
# 9   72.5297  7.671318
# 10  70.4982  7.456451



cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"chop_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas"
           ,"swedish_uc_old_gwas","mccauley_gsa","ccfa_gsa","cedars_gsa")

pca$study<-NA

pca$study[which(pca$cohort=="australia_omniexome")]<-"Australia"
pca$study[which(pca$cohort %in% c("gwas1","gwas2","all_hce"))]<-"UK"
pca$study[which(pca$cohort %in% c("spain_gsa","basque_gsa"))]<-"Spain"
pca$study[which(pca$cohort %in% c("kiel_austria_sibdcs_gsa","german_affy6_old_gwas"))]<-"Germany"                 
pca$study[which(pca$cohort %in% c("italy_gsa"))]<-"Italy"   
pca$study[which(pca$cohort %in% c("netherlands_gsa"))]<-"Netherlands"                     
pca$study[which(pca$cohort %in% c("slovenia_gsa"))]<-"Slovenia"                     
pca$study[which(pca$cohort %in% c("sweden_gsa","swedish_uc_old_gwas"))]<-"Sweden"                    
pca$study[which(pca$cohort %in% c("finland_illugwas"))]<-"Finland"
pca$study[which(pca$cohort %in% c("chop_old_gwas"))]<-"CHOP"
pca$study[which(pca$cohort %in% c("norway_affy6_old_gwas"))]<-"Norway"
pca$study[which(pca$cohort %in% c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","belgium_louis_gsa","belgium_franchimont_gsa"
                                  ,"belgium_vermeire_gsa"))]<-"Belgium"
pca$study[which(pca$cohort %in% c("lithuania_gsa"))]<-"Lithuania"
pca$study[which(pca$cohort %in% c("prism_nfe_gwas","prism_nfe_gsa","cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","cedars_gsa",
                                  "niddk_broad_gsa","niddk_feinstein_gsa","niddk_old_gwas","pittsburgh_gsa","mccauley_gsa","ccfa_gsa"))]<-"USA"

table(pca$study)
# Australia     Belgium        CHOP     Finland     Germany       Italy 
# 1251        8534        8475         443       16937         947 
# Lithuania Netherlands      Norway    Slovenia       Spain      Sweden 
# 2198        4352         547         264        4935        2584 
# UK         USA 
# 35114       26531 

pca$study<-as.factor(pca$study)
pca$study<-factor(pca$study, levels=c("Spain","Italy","Slovenia","Belgium","Germany","Netherlands","UK","Lithuania",
                                      "Sweden","Norway","Finland","Australia","CHOP","USA"))


table(pca$study,useNA="ifany")
# Spain       Italy    Slovenia     Belgium     Germany Netherlands 
# 4935         947         264        8534       16937        4352 
# UK   Lithuania      Sweden      Norway     Finland   Australia 
# 35114        2198        2584         547         443        1251 
# CHOP         USA 
# 8475       26531



########################
# 6.- PLOT EUR RESULTS #
########################

######################################################
# 4.1 PLOT 1000GP ALONE, COLOUR BY COUNTRY:

# create plots PC1 to PC10

# pna<-qplot()+theme(
#   panel.background = element_rect(fill = "transparent") # bg of the panel
#   , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
#   , panel.grid.major = element_blank() # get rid of major grid
#   , panel.grid.minor = element_blank() # get rid of minor grid
#   , legend.background = element_rect(fill = "transparent") # get rid of legend bg
#   , axis.title=element_blank()
#   , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
# )
# 
# 
# p11<-qplot(PC1,PC1, data = pca, colour = study)
# p12<-qplot(PC2,PC1, data = pca, colour = study)
# p13<-qplot(PC3,PC1, data = pca, colour = study)
# p14<-qplot(PC4,PC1, data = pca, colour = study)
# p15<-qplot(PC5,PC1, data = pca, colour = study)
# p16<-qplot(PC6,PC1, data = pca, colour = study)
# p17<-qplot(PC7,PC1, data = pca, colour = study)
# p18<-qplot(PC8,PC1, data = pca, colour = study)
# p19<-qplot(PC9,PC1, data = pca, colour = study)
# p110<-qplot(PC10,PC1, data = pca, colour = study)
# 
# p21<-qplot(PC1,PC2, data = pca, colour = study)
# p22<-qplot(PC2,PC2, data = pca, colour = study)
# p23<-qplot(PC3,PC2, data = pca, colour = study)
# p24<-qplot(PC4,PC2, data = pca, colour = study)
# p25<-qplot(PC5,PC2, data = pca, colour = study)
# p26<-qplot(PC6,PC2, data = pca, colour = study)
# p27<-qplot(PC7,PC2, data = pca, colour = study)
# p28<-qplot(PC8,PC2, data = pca, colour = study)
# p29<-qplot(PC9,PC2, data = pca, colour = study)
# p210<-qplot(PC10,PC2, data = pca, colour = study)
# 
# p31<-qplot(PC1,PC3, data = pca, colour = study)
# p32<-qplot(PC2,PC3, data = pca, colour = study)
# p33<-qplot(PC3,PC3, data = pca, colour = study)
# p34<-qplot(PC4,PC3, data = pca, colour = study)
# p35<-qplot(PC5,PC3, data = pca, colour = study)
# p36<-qplot(PC6,PC3, data = pca, colour = study)
# p37<-qplot(PC7,PC3, data = pca, colour = study)
# p38<-qplot(PC8,PC3, data = pca, colour = study)
# p39<-qplot(PC9,PC3, data = pca, colour = study)
# p310<-qplot(PC10,PC3, data = pca, colour = study)
# 
# p41<-qplot(PC1,PC4, data = pca, colour = study)
# p42<-qplot(PC2,PC4, data = pca, colour = study)
# p43<-qplot(PC3,PC4, data = pca, colour = study)
# p44<-qplot(PC4,PC4, data = pca, colour = study)
# p45<-qplot(PC5,PC4, data = pca, colour = study)
# p46<-qplot(PC6,PC4, data = pca, colour = study)
# p47<-qplot(PC7,PC4, data = pca, colour = study)
# p48<-qplot(PC8,PC4, data = pca, colour = study)
# p49<-qplot(PC9,PC4, data = pca, colour = study)
# p410<-qplot(PC10,PC4, data = pca, colour = study)
# 
# p51<-qplot(PC1,PC5, data = pca, colour = study)
# p52<-qplot(PC2,PC5, data = pca, colour = study)
# p53<-qplot(PC3,PC5, data = pca, colour = study)
# p54<-qplot(PC4,PC5, data = pca, colour = study)
# p55<-qplot(PC5,PC5, data = pca, colour = study)
# p56<-qplot(PC6,PC5, data = pca, colour = study)
# p57<-qplot(PC7,PC5, data = pca, colour = study)
# p58<-qplot(PC8,PC5, data = pca, colour = study)
# p59<-qplot(PC9,PC5, data = pca, colour = study)
# p510<-qplot(PC10,PC5, data = pca, colour = study)
# 
# p61<-qplot(PC1,PC6, data = pca, colour = study)
# p62<-qplot(PC2,PC6, data = pca, colour = study)
# p63<-qplot(PC3,PC6, data = pca, colour = study)
# p64<-qplot(PC4,PC6, data = pca, colour = study)
# p65<-qplot(PC5,PC6, data = pca, colour = study)
# p66<-qplot(PC6,PC6, data = pca, colour = study)
# p67<-qplot(PC7,PC6, data = pca, colour = study)
# p68<-qplot(PC8,PC6, data = pca, colour = study)
# p69<-qplot(PC9,PC6, data = pca, colour = study)
# p610<-qplot(PC10,PC6, data = pca, colour = study)
# 
# p61<-qplot(PC1,PC6, data = pca, colour = study)
# p62<-qplot(PC2,PC6, data = pca, colour = study)
# p63<-qplot(PC3,PC6, data = pca, colour = study)
# p64<-qplot(PC4,PC6, data = pca, colour = study)
# p65<-qplot(PC5,PC6, data = pca, colour = study)
# p66<-qplot(PC6,PC6, data = pca, colour = study)
# p67<-qplot(PC7,PC6, data = pca, colour = study)
# p68<-qplot(PC8,PC6, data = pca, colour = study)
# p69<-qplot(PC9,PC6, data = pca, colour = study)
# p610<-qplot(PC10,PC6, data = pca, colour = study)
# 
# p71<-qplot(PC1,PC7, data = pca, colour = study)
# p72<-qplot(PC2,PC7, data = pca, colour = study)
# p73<-qplot(PC3,PC7, data = pca, colour = study)
# p74<-qplot(PC4,PC7, data = pca, colour = study)
# p75<-qplot(PC5,PC7, data = pca, colour = study)
# p76<-qplot(PC6,PC7, data = pca, colour = study)
# p77<-qplot(PC7,PC7, data = pca, colour = study)
# p78<-qplot(PC8,PC7, data = pca, colour = study)
# p79<-qplot(PC9,PC7, data = pca, colour = study)
# p710<-qplot(PC10,PC7, data = pca, colour = study)
# 
# p81<-qplot(PC1,PC8, data = pca, colour = study)
# p82<-qplot(PC2,PC8, data = pca, colour = study)
# p83<-qplot(PC3,PC8, data = pca, colour = study)
# p84<-qplot(PC4,PC8, data = pca, colour = study)
# p85<-qplot(PC5,PC8, data = pca, colour = study)
# p86<-qplot(PC6,PC8, data = pca, colour = study)
# p87<-qplot(PC7,PC8, data = pca, colour = study)
# p88<-qplot(PC8,PC8, data = pca, colour = study)
# p89<-qplot(PC9,PC8, data = pca, colour = study)
# p810<-qplot(PC10,PC8, data = pca, colour = study)
# 
# p91<-qplot(PC1,PC9, data = pca, colour = study)
# p92<-qplot(PC2,PC9, data = pca, colour = study)
# p93<-qplot(PC3,PC9, data = pca, colour = study)
# p94<-qplot(PC4,PC9, data = pca, colour = study)
# p95<-qplot(PC5,PC9, data = pca, colour = study)
# p96<-qplot(PC6,PC9, data = pca, colour = study)
# p97<-qplot(PC7,PC9, data = pca, colour = study)
# p98<-qplot(PC8,PC9, data = pca, colour = study)
# p99<-qplot(PC9,PC9, data = pca, colour = study)
# p910<-qplot(PC10,PC9, data = pca, colour = study)
# 
# p101<-qplot(PC1,PC10, data = pca, colour = study)
# p102<-qplot(PC2,PC10, data = pca, colour = study)
# p103<-qplot(PC3,PC10, data = pca, colour = study)
# p104<-qplot(PC4,PC10, data = pca, colour = study)
# p105<-qplot(PC5,PC10, data = pca, colour = study)
# p106<-qplot(PC6,PC10, data = pca, colour = study)
# p107<-qplot(PC7,PC10, data = pca, colour = study)
# p108<-qplot(PC8,PC10, data = pca, colour = study)
# p109<-qplot(PC9,PC10, data = pca, colour = study)
# p1010<-qplot(PC10,PC10, data = pca, colour = study)
# 
# 
# r1<-ggarrange(pna,p12,p13,p14,p15,p16,p17,p18,p19,p110,ncol=10,legend=c("none"))
# r2<-ggarrange(p21,pna,p23,p24,p25,p26,p27,p28,p29,p210,ncol=10,legend=c("none"))
# r3<-ggarrange(p31,p32,pna,p34,p35,p36,p37,p38,p39,p310,ncol=10,legend=c("none"))
# r4<-ggarrange(p41,p42,p43,pna,p45,p46,p47,p48,p49,p410,ncol=10,legend=c("none"))
# r5<-ggarrange(p51,p52,p53,p54,pna,p56,p57,p58,p59,p510,ncol=10,legend=c("none"))
# r6<-ggarrange(p61,p62,p63,p64,p65,pna,p67,p68,p69,p610,ncol=10,legend=c("none"))
# r7<-ggarrange(p71,p72,p73,p74,p75,p76,pna,p78,p79,p710,ncol=10,legend=c("none"))
# r8<-ggarrange(p81,p82,p83,p84,p85,p86,p87,pna,p89,p810,ncol=10,legend=c("none"))
# r9<-ggarrange(p91,p92,p93,p94,p95,p96,p97,p98,pna,p910,ncol=10,legend=c("none"))
# r10<-ggarrange(p101,p102,p103,p104,p105,p106,p107,p108,p109,pna,ncol=10,common.legend = TRUE,legend=c("bottom"))
# 
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_PC7_PC8_PC9_PC10_3kvariants_iibdgc_EURonly.pdf",sep=""),width=45,height=45)
# ggarrange(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,nrow=10,common.legend = TRUE,legend=c("bottom"),heights = c(1,1,1,1,1,1,1,1,1,1.2))
# dev.off()
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_PC7_PC8_PC9_PC10_3kvariants_iibdgc_EURonly.pdf ~/tmp_plots/",sep=""))


### SELF-REPORTED ASKENAZI:

# only new cohorts and cedars will have data for this:
pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
colnames(pheno)[14]<-"self_jewish"
pheno_cd<-fread(paste(path,"pheno/list_cedars_old_gwas_ids_TH.txt",sep=""),head=T,sep="\t")

pheno<-rbind(pheno[,c("sample_id","self_jewish")],pheno_cd[,c("sample_id","self_jewish")])
table(pheno$self_jewish)
#               No    Unknown        Yes     Jewish Non-Jewish 
# 23348      33975      16224       4546       1217       1513 

pheno$self_jewish<-as.character(pheno$self_jewish)
pheno$self_jewish[which(pheno$self_jewish %in% c("","Unknown"))]<-"Unknown"
pheno$self_jewish[which(pheno$self_jewish %in% c("No","Non-Jewish"))]<-"Non-Jewish"
pheno$self_jewish[which(pheno$self_jewish %in% c("Yes","Jewish"))]<-"Jewish"
table(pheno$self_jewish)
# Jewish Non-Jewish    Unknown 
# 5763      35488      39572 



pca$self_jewish<-"Unknown"
pca$self_jewish[which(pca$FID %in% pheno$sample_id[which(pheno$self_jewish=="Non-Jewish")])]<-"Non-Jewish"
pca$self_jewish[which(pca$FID %in% pheno$sample_id[which(pheno$self_jewish=="Jewish")])]<-"Jewish"
table(pca$self_jewish)
#     Jewish Non-Jewish    Unknown 
#       5473      30081      77558


# p11<-qplot(pca$PC1,pca$PC1, data = pca, colour = self_jewish)
# p12<-qplot(pca$PC2,pca$PC1, data = pca, colour = self_jewish)
# p13<-qplot(pca$PC3,pca$PC1, data = pca, colour = self_jewish)
# p14<-qplot(pca$PC4,pca$PC1, data = pca, colour = self_jewish)
# p15<-qplot(pca$PC5,pca$PC1, data = pca, colour = self_jewish)
# p16<-qplot(pca$PC6,pca$PC1, data = pca, colour = self_jewish)
# 
# p21<-qplot(pca$PC1,pca$PC2, data = pca, colour = self_jewish)
# p22<-qplot(pca$PC2,pca$PC2, data = pca, colour = self_jewish)
# p23<-qplot(pca$PC3,pca$PC2, data = pca, colour = self_jewish)
# p24<-qplot(pca$PC4,pca$PC2, data = pca, colour = self_jewish)
# p25<-qplot(pca$PC5,pca$PC2, data = pca, colour = self_jewish)
# p26<-qplot(pca$PC6,pca$PC2, data = pca, colour = self_jewish)
# 
# p31<-qplot(pca$PC1,pca$PC3, data = pca, colour = self_jewish)
# p32<-qplot(pca$PC2,pca$PC3, data = pca, colour = self_jewish)
# p33<-qplot(pca$PC3,pca$PC3, data = pca, colour = self_jewish)
# p34<-qplot(pca$PC4,pca$PC3, data = pca, colour = self_jewish)
# p35<-qplot(pca$PC5,pca$PC3, data = pca, colour = self_jewish)
# p36<-qplot(pca$PC6,pca$PC3, data = pca, colour = self_jewish)
# 
# p41<-qplot(pca$PC1,pca$PC4, data = pca, colour = self_jewish)
# p42<-qplot(pca$PC2,pca$PC4, data = pca, colour = self_jewish)
# p43<-qplot(pca$PC3,pca$PC4, data = pca, colour = self_jewish)
# p44<-qplot(pca$PC4,pca$PC4, data = pca, colour = self_jewish)
# p45<-qplot(pca$PC5,pca$PC4, data = pca, colour = self_jewish)
# p46<-qplot(pca$PC6,pca$PC4, data = pca, colour = self_jewish)
# 
# p51<-qplot(pca$PC1,pca$PC5, data = pca, colour = self_jewish)
# p52<-qplot(pca$PC2,pca$PC5, data = pca, colour = self_jewish)
# p53<-qplot(pca$PC3,pca$PC5, data = pca, colour = self_jewish)
# p54<-qplot(pca$PC4,pca$PC5, data = pca, colour = self_jewish)
# p55<-qplot(pca$PC5,pca$PC5, data = pca, colour = self_jewish)
# p56<-qplot(pca$PC6,pca$PC5, data = pca, colour = self_jewish)
# 
# p61<-qplot(pca$PC1,pca$PC6, data = pca, colour = self_jewish)
# p62<-qplot(pca$PC2,pca$PC6, data = pca, colour = self_jewish)
# p63<-qplot(pca$PC3,pca$PC6, data = pca, colour = self_jewish)
# p64<-qplot(pca$PC4,pca$PC6, data = pca, colour = self_jewish)
# p65<-qplot(pca$PC5,pca$PC6, data = pca, colour = self_jewish)
# p66<-qplot(pca$PC6,pca$PC6, data = pca, colour = self_jewish)
# 
# 
# r1<-ggarrange(pna,p12,p13,p14,p15,p16,ncol=6,legend=c("none"))
# r2<-ggarrange(p21,pna,p23,p24,p25,p26,ncol=6,legend=c("none"))
# r3<-ggarrange(p31,p32,pna,p34,p35,p36,ncol=6,legend=c("none"))
# r4<-ggarrange(p41,p42,p43,pna,p45,p36,ncol=6,legend=c("none"))
# r5<-ggarrange(p51,p52,p53,p54,pna,p56,ncol=6,legend=c("none"))
# r6<-ggarrange(p61,p62,p63,p64,p65,pna,ncol=6,common.legend = TRUE,legend=c("bottom"))
# dev.off()
# 
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_3kvariants_iibdgc_EURonly_sr_jewish.pdf",sep=""),width=30,height=30)
# ggarrange(r1,r2,r3,r4,r5,r6,nrow=6,common.legend = TRUE,legend=c("bottom"),heights = c(1,1,1,1,1,1.2))
# dev.off()
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_3kvariants_iibdgc_EURonly_sr_jewish.pdf ~/tmp_plots/",sep=""))

# ############################################
# # PC1 - JEWISH ANCESTRY:
# 
# ylims<-c(min(pca$PC1),max(pca$PC1))
# xlims<-c(min(pca$PC2),max(pca$PC2))
# 
# p1<-ggplot(pca,aes(PC2,PC1)) +
#   geom_point(aes(color = study)) + 
#   xlim(xlims) + ylim(ylims)
# 
# p2<-ggplot(pca[which(pca$self_jewish=="Jewish"),],aes(PC2,PC1,color = self_jewish)) +
#   geom_point() + scale_color_manual(values=c("#E69F00")) +
#   xlim(xlims) + ylim(ylims)
# 
# p3<-ggplot(pca[which(pca$self_jewish=="Non-Jewish"),],aes(PC2,PC1,color = self_jewish)) +
#   geom_point() + scale_color_manual(values=c("#56B4E9")) +
#   xlim(xlims) + ylim(ylims)
# 
# p4<-ggplot(pca[which(pca$self_jewish=="Unknown"),],aes(PC2,PC1,color = self_jewish)) +
#   geom_point() + scale_color_manual(values=c("#999999")) +
#   xlim(xlims) + ylim(ylims)
# 
# p100<-ggarrange(p1,pna,pna,ncol=3)
# p234<-ggarrange(p2,p3,p4,ncol=3)
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_PC2_EURonly_sr_jewish.pdf",sep=""),width=22,height=12)
# ggarrange(p100,p234,nrow=2)
# dev.off()
# 
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_PC2_EURonly_sr_jewish.pdf ~/tmp_plots/",sep=""))
# 
# 
# ########################
# # DENSITY PLOT FOR PC1:
# 
# xlims<-c(min(pca$PC1),max(pca$PC1))
# 
# p1<-ggplot(pca[which(pca$self_jewish %in% c("Jewish","Non-Jewish")),], aes(x=PC1, fill=self_jewish)) +
#   geom_density(alpha=0.4) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
#   xlim(xlims)
# 
# p2<-ggplot(pca[which(pca$self_jewish %in% c("Unknown")),], aes(x=PC1, fill=self_jewish)) +
#   geom_density(alpha=0.4) + scale_fill_manual(values=c("#999999")) +
#   xlim(xlims)
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_density_EURonly_sr_jewish.pdf",sep=""),width=14,height=10)
# ggarrange(p1,p2,nrow=2)
# dev.off()
# 
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_density_EURonly_sr_jewish.pdf ~/tmp_plots/",sep=""))
# 
# #### histogram
# 
# fd=function(x) {
#   n=length(x)
#   r=IQR(x)
#   2*r/n^(1/3)
# }
# 
# fd_bin<-fd(pca$PC1)
# 
# p1<-ggplot(pca[which(pca$self_jewish %in% c("Jewish","Non-Jewish")),], aes(x=PC1,color=self_jewish, fill=self_jewish)) +
#   geom_histogram(position="identity", alpha=0.5, binwidth=fd_bin) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) + 
#   scale_color_manual(values=c("#E69F00", "#56B4E9")) +
#   xlim(xlims)
# 
# p2<-ggplot(pca[which(pca$self_jewish %in% c("Unknown")),], aes(x=PC1,color=self_jewish, fill=self_jewish)) +
#   geom_histogram(position="identity", alpha=0.5, binwidth=fd_bin) + scale_fill_manual(values=c("#999999")) + 
#   scale_color_manual(values=c("#999999")) +
#   xlim(xlims)
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_histogram_EURonly_sr_jewish.pdf",sep=""),width=14,height=10)
# ggarrange(p1,p2,nrow=2)
# dev.off()
# 
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_histogram_EURonly_sr_jewish.pdf ~/tmp_plots/",sep=""))
# 
# 
# ### ESTIMATE CUTOFF POINT:
# 
# cutpoint<- 0.045
# 
# ### density:
# 
# p1<-ggplot(pca[which(pca$self_jewish %in% c("Jewish","Non-Jewish")),], aes(x=PC1, fill=self_jewish)) +
#   geom_density(alpha=0.4) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) + geom_vline(xintercept = cutpoint, linetype="dotted", color = "darkgrey") +
#   xlim(xlims)
# p2<-ggplot(pca[which(pca$self_jewish %in% c("Unknown")),], aes(x=PC1, fill=self_jewish)) +
#   geom_density(alpha=0.4) + scale_fill_manual(values=c("#999999")) + geom_vline(xintercept = cutpoint, linetype="dotted", color = "black") +
#   xlim(xlims)
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_density_EURonly_sr_jewish_cutoff.pdf",sep=""),width=14,height=10)
# ggarrange(p1,p2,nrow=2)
# dev.off()
# 
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_density_EURonly_sr_jewish_cutoff.pdf ~/tmp_plots/",sep=""))
# 
# 
# ### histogram:
# 
# p1<-ggplot(pca[which(pca$self_jewish %in% c("Jewish","Non-Jewish")),], aes(x=PC1,color=self_jewish, fill=self_jewish)) +
#   geom_histogram(position="identity", alpha=0.5, binwidth=fd_bin) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) + 
#   scale_color_manual(values=c("#E69F00", "#56B4E9")) + geom_vline(xintercept = cutpoint, linetype="dotted", color = "black") +
#   xlim(xlims)
# 
# p2<-ggplot(pca[which(pca$self_jewish %in% c("Unknown")),], aes(x=PC1,color=self_jewish, fill=self_jewish)) +
#   geom_histogram(position="identity", alpha=0.5, binwidth=fd_bin) + scale_fill_manual(values=c("#999999")) + 
#   scale_color_manual(values=c("#999999")) + geom_vline(xintercept = cutpoint, linetype="dotted", color = "black") +
#   xlim(xlims)
# 
# pca$tmp<-"All IIBDG EUR"
# 
# p3<-ggplot(pca, aes(x=PC1,color=tmp, fill=tmp)) +
#   geom_histogram(position="identity", alpha=0.5, binwidth=fd_bin) + scale_fill_manual(values=c("lightblue4")) + 
#   scale_color_manual(values=c("lightblue4")) + geom_vline(xintercept = cutpoint, linetype="dotted", color = "black") +
#   xlim(xlims)
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_histogram_EURonly_sr_jewish_cutoff.pdf",sep=""),width=14,height=15)
# ggarrange(p1,p2,p3,nrow=3)
# dev.off()
# 
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_PC1_histogram_EURonly_sr_jewish_cutoff.pdf ~/tmp_plots/",sep=""))
# 


# reclassify:

cutpoint<- 0.045

pca$pca_jewish<-NA
pca$pca_jewish[which(pca$PC1<cutpoint)]<-"Non-Jewish"
pca$pca_jewish[which(pca$PC1>=cutpoint)]<-"Jewish"
table(pca$pca_jewish)
# Jewish Non-Jewish 
# 7894     105218    

# table(pca$pca_jewish,pca$self_jewish)
# #            Jewish Non-Jewish Unknown
# # Jewish       4186       1001    2707
# # Non-Jewish   1287      29080   74851
# 
# 4186/sum(4186,1001)
# # [1] 0.8070175 # sensitibity
# 
# 29080/sum(29080,1287)
# # [1] 0.9576185 # specificity
# 
# 
# 4186/sum(4186,1287)
# # [1] 0.7648456 PPV
# 
# 29080/sum(29080,1303)
# # [1] 0.9571142 NPV

# 
# #######################################################
# # PLOT CASES AND CONTROLS PER STUDY, FOR JEWISH SR:
# 
# pca$pheno<-NA
# pca$pheno[which(pca$FID %in% fam$V1[which(fam$V6==1)])]<-"1"
# pca$pheno[which(pca$FID %in% fam$V1[which(fam$V6==2)])]<-"2"

# table(pca$cohort[which(pca$pheno=="2")],pca$self_jewish[which(pca$pheno=="2")],useNA="ifany")
# # Jewish Non-Jewish Unknown
# # all_hce                      0          0   12342
# # australia_omniexome          1        658       0
# # basque_gsa                   0          0     518
# # belgium_franchimont_gsa      0          0     880
# # belgium_inf1_old_gwas        0          0     513
# # belgium_inf2_old_gwas        0          0     159
# # belgium_louis_gsa            0          0     897
# # belgium_vermeire_gsa         0          0    3098
# # ccfa_gsa                     0          0    1931
# # cedars_370k_old_gwas       236        273       0
# # cedars_610k_old_gwas       413        417       0
# # cedars_gsa                 511       1248     125
# # cedars_omni_old_gwas       537        566       1
# # chop_old_gwas                0          0    2371
# # finland_illugwas             0          0     443
# # german_affy6_old_gwas        0          0    1045
# # gwas1                        0          0    1740
# # gwas2                        0          0    2353
# # italy_gsa                    0        576       0
# # kiel_austria_sibdcs_gsa      0       9597       0
# # lithuania_gsa                0          0    1072
# # mccauley_gsa               325        402      33
# # netherlands_gsa              4        402    3412
# # niddk_broad_gsa            665       2860     632
# # niddk_feinstein_gsa       2243       2752     191
# # niddk_old_gwas               0          0    1808
# # norway_affy6_old_gwas        0          0     266
# # pittsburgh_gsa               0       1251       0
# # prism_nfe_gsa                0          0     401
# # prism_nfe_gwas               0          0     510
# # slovenia_gsa                 0         87       0
# # spain_gsa                    0          0    1957
# # sweden_gsa                   0          0     390
# # swedish_uc_old_gwas          0          0     921
# 
# 
# table(pca$cohort[which(pca$pheno=="1")],pca$self_jewish[which(pca$pheno=="1")],useNA="ifany")
# # Jewish Non-Jewish Unknown
# # all_hce                      0          0   10356
# # australia_omniexome          0        592       0
# # basque_gsa                   0          0     981
# # belgium_franchimont_gsa      0          0     586
# # belgium_inf1_old_gwas        0          0     891
# # belgium_inf2_old_gwas        0          0     111
# # belgium_louis_gsa            0          0     592
# # belgium_vermeire_gsa         0          0     807
# # cedars_gsa                 255        533      13
# # cedars_omni_old_gwas         1          2       0
# # chop_old_gwas                0          0    6104
# # german_affy6_old_gwas        0          0    1759
# # gwas1                        0          0    2923
# # gwas2                        0          0    5400
# # italy_gsa                    0        371       0
# # kiel_austria_sibdcs_gsa      0       4536       0
# # lithuania_gsa                0          0    1126
# # netherlands_gsa              0          0     534
# # niddk_broad_gsa            151        832       8
# # niddk_feinstein_gsa        131       1475     175
# # niddk_old_gwas               0          0     922
# # norway_affy6_old_gwas        0          0     281
# # pittsburgh_gsa               0        474     990
# # prism_nfe_gsa                0          0      28
# # prism_nfe_gwas               0          0     210
# # slovenia_gsa                 0        177       0
# # spain_gsa                    0          0    1479
# # sweden_gsa                   0          0     934
# # swedish_uc_old_gwas          0          0     339
# 
# 
# #######################################################
# # PLOT CASES AND CONTROLS PER STUDY, FOR JEWISH PCA:
# 
# 
# table(pca$cohort[which(pca$pheno=="2")],pca$pca_jewish[which(pca$pheno=="2")],useNA="ifany")
# #                         Jewish Non-Jewish
# # all_hce                    144      12198
# # australia_omniexome         12        647
# # basque_gsa                   1        517
# # belgium_franchimont_gsa     46        834
# # belgium_inf1_old_gwas       15        498
# # belgium_inf2_old_gwas        4        155
# # belgium_louis_gsa           21        876
# # belgium_vermeire_gsa        73       3025
# # ccfa_gsa                   147       1784
# # cedars_370k_old_gwas       173        336
# # cedars_610k_old_gwas       290        540
# # cedars_gsa                 474       1410
# # cedars_omni_old_gwas       366        738
# # chop_old_gwas              225       2146
# # finland_illugwas             1        442
# # german_affy6_old_gwas        2       1043
# # gwas1                       44       1696
# # gwas2                        1       2352
# # italy_gsa                  118        458
# # kiel_austria_sibdcs_gsa    173       9424
# # lithuania_gsa                6       1066
# # mccauley_gsa               291        469
# # netherlands_gsa             42       3776
# # niddk_broad_gsa            823       3334
# # niddk_feinstein_gsa       1994       3192
# # niddk_old_gwas             421       1387
# # norway_affy6_old_gwas        0        266
# # pittsburgh_gsa              32       1219
# # prism_nfe_gsa               66        335
# # prism_nfe_gwas              72        438
# # slovenia_gsa                 0         87
# # spain_gsa                   20       1937
# # sweden_gsa                  14        376
# # swedish_uc_old_gwas          0        921
# 
# table(pca$cohort[which(pca$pheno=="1")],pca$pca_jewish[which(pca$pheno=="1")],useNA="ifany")
# # Jewish Non-Jewish
# # all_hce                     46      10310
# # australia_omniexome          2        590
# # basque_gsa                   6        975
# # belgium_franchimont_gsa     25        561
# # belgium_inf1_old_gwas       14        877
# # belgium_inf2_old_gwas        1        110
# # belgium_louis_gsa           13        579
# # belgium_vermeire_gsa         4        803
# # cedars_gsa                 247        554
# # cedars_omni_old_gwas         1          2
# # chop_old_gwas              472       5632
# # german_affy6_old_gwas        2       1757
# # gwas1                       15       2908
# # gwas2                        0       5400
# # italy_gsa                   81        290
# # kiel_austria_sibdcs_gsa      9       4527
# # lithuania_gsa                0       1126
# # netherlands_gsa              0        534
# # niddk_broad_gsa            151        840
# # niddk_feinstein_gsa        153       1628
# # niddk_old_gwas             410        512
# # norway_affy6_old_gwas        0        281
# # pittsburgh_gsa              87       1377
# # prism_nfe_gsa                4         24
# # prism_nfe_gwas              26        184
# # slovenia_gsa                 0        177
# # spain_gsa                   14       1465
# # sweden_gsa                   0        934
# # swedish_uc_old_gwas          0        339
# 
# 
# 
# 
# 

# # PC2 - PC3
# 
# xlims<-c(min(pca$PC3),max(pca$PC3))
# ylims<-c(min(pca$PC2),max(pca$PC2))
# 
# p1<-ggplot(pca,aes(PC3,PC2)) +
#   geom_point(aes(color = study)) + 
#   xlim(xlims) + ylim(ylims)
# 
# p2<-ggplot(pca,aes(PC3,PC2)) +
#   geom_point(aes(color = study)) + 
#   xlim(xlims) + ylim(ylims)
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_PC2_PC3_EURonly.pdf",sep=""),width=7,height=7)
# p23
# dev.off()
# 
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_PC2_PC3_EURonly.pdf ~/tmp_plots/",sep=""))



# # PC2 - PC3 -  all studies
# 
# table(pca$study)
# # Spain       Italy    Slovenia     Belgium     Germany Netherlands 
# # 4935         947         264        8534       16937        4352 
# # UK   Lithuania      Sweden      Norway     Finland   Australia 
# # 35114        2198        2584         547         443        1251 
# # CHOP         USA 
# # 8475       26531 
# 
# 
# # to emulate ggplot2 palette:
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# 
# n = length(levels(pca$study))
# cols = gg_color_hue(n)
# 
# pca$PC2_inv<-pca$PC2*-1
# 
# xlims<-c(min(pca$PC2_inv),max(pca$PC2_inv))
# ylims<-c(min(pca$PC3),max(pca$PC3))
# 
# p0<-ggplot(pca,aes(PC2_inv,PC3)) +
#   geom_point(aes(color = study)) + 
#   xlim(xlims) + ylim(ylims)
# 
# for (i in 1:length(levels(pca$study))) {
#   print(i)
#   assign(paste("p",i,sep=""),ggplot(pca[which(pca$study==levels(pca$study)[i]),],aes(PC2_inv,PC3,color=study)) +
#            geom_point() + scale_color_manual(values=cols[i]) + 
#            xlim(xlims) + ylim(ylims))
# }
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_PC2_PC3_EURonly_per_study.pdf",sep=""),width=30,height=27)
# print(ggarrange(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,legend=c("right"),ncol=4,nrow=5))
# dev.off()
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_PC2_PC3_EURonly_per_study.pdf ~/tmp_plots/",sep=""))
# 

##### create two additional lists of samples:
# - EUR samples all
# - EUR samples non-jewish
# - EUR samples Jewish


# write.table(pca[which(pca$pca_jewish=="Jewish"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eur_jewish_ancestry_samples_withDuplicates",sep=""),
#             col.names=F,row.names=F,quote=F,sep="\t")
# write.table(pca[which(pca$pca_jewish=="Non-Jewish"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eur_nonjewish_ancestry_samples_withDuplicates",sep=""),
#             col.names=F,row.names=F,quote=F,sep="\t")
# write.table(pca,paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_EUR_withDuplicates_ref_pcs_new_projection_edited.sscore",sep=""),
#             col.names=T,row.names=F,quote=F,sep="\t")
# 
#############################################################################################################################################################

# EXPLORE what the other PCs account for?


eigenval
#          V1   var_exp
# 1  214.6160 22.699497 # Jewish ancestry
# 2  115.4610 12.212074 # north south
# 3   92.6266  9.796927 # EUR east west
# 4   78.9639  8.351851 # Spain to Italy
# 5   76.2488  8.064680
# 6   76.0964  8.048561
# 7   74.5226  7.882103
# 8   73.9027  7.816538
# 9   72.5297  7.671318
# 10  70.4982  7.456451

array<-c("illumina370","illumina550","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

illumina370<-c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","swedish_uc_old_gwas","niddk_old_gwas")
illumina550<-c("chop_old_gwas")
affymetrix6<-c("german_affy6_old_gwas","norway_affy6_old_gwas","gwas2")
humanomniexpress<-c("australia_omniexome")
affymetrix500<-c("gwas1")
humancoreexome<-c("all_hce")
humanomni1<-c("pittsburgh_gsa")
quad610<-c("spain_gsa")
gsa<-c("italy_gsa","kiel_austria_sibdcs_gsa"
       ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
       ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
       ,"mccauley_gsa","ccfa_gsa","cedars_gsa")
illuminaexome<-c("prism_nfe_gwas")

pca$array<-NA
pca$array[which(pca$cohort %in% illumina370)]<-"illumina370"
pca$array[which(pca$cohort %in% illumina550)]<-"illumina550"
pca$array[which(pca$cohort %in% affymetrix6)]<-"affymetrix6"
pca$array[which(pca$cohort %in% humanomniexpress)]<-"humanomniexpress"
pca$array[which(pca$cohort %in% affymetrix500)]<-"affymetrix500"
pca$array[which(pca$cohort %in% humancoreexome)]<-"humancoreexome"
pca$array[which(pca$cohort %in% humanomni1)]<-"humanomni1"
pca$array[which(pca$cohort %in% quad610)]<-"quad610"
pca$array[which(pca$cohort %in% gsa)]<-"gsa"
pca$array[which(pca$cohort %in% illuminaexome)]<-"illuminaexome"

table(pca$cohort,pca$array) # NA are studies not included in analysis (like old CEDARS)

# add pheno CD/UC/IBDu/Ctr:
pheno<-read.table(paste(path,"pheno/phenotype_data_for_all_cohorts_Dec2020_analysis.tsv",sep=""),head=T)

pca$pheno<-NA
pca$pheno[which(pca$FID %in% pheno$FID[which(pheno$ibd==0)])]<-"Control"
pca$pheno[which(pca$FID %in% pheno$FID[which(pheno$cd==0)])]<-"Control"
pca$pheno[which(pca$FID %in% pheno$FID[which(pheno$uc==0)])]<-"Control"

pca$pheno[which(pca$FID %in% pheno$FID[which(pheno$ibd==1)])]<-"IBDu"
pca$pheno[which(pca$FID %in% pheno$FID[which(pheno$cd==1)])]<-"CD"
pca$pheno[which(pca$FID %in% pheno$FID[which(pheno$uc==1)])]<-"UC"

pca$pheno<-factor(pca$pheno, levels=c("CD","IBDu","UC","Control"))

### histogram per PC to show what they account for:

## PC1 - Jewish ancestry

sp <- split(pca$PC1, pca$self_jewish)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(self_jewish = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p1<-ggplot(pca, aes(x=PC1,color=self_jewish, fill=self_jewish)) + geom_density(alpha=0.4) + scale_fill_manual(values=c("#E69F00", "#56B4E9","#999999")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9","#999999")) + geom_text(data = a, aes(x = xmax, y = ymax,label=self_jewish, vjust = - 0.5)) + 
  ggtitle("Jewish ancestry") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())

pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_european_ancestry_samples_PC1.pdf",sep=""),width=14,height=5)
p1
dev.off()

system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_european_ancestry_samples_PC1.pdf ~/tmp_plots/",sep=""))


## PC2 - Southwest to northeast

sp <- split(pca$PC2, pca$study)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(study = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p2<-ggplot(pca, aes(x=PC2,color=study, fill=study)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = study, vjust = - 0.5)) + 
  ggtitle("Country") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())

pdf(paste(path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC2.pdf",sep=""),width=14,height=5)
p2
dev.off()

system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC2.pdf ~/tmp_plots/",sep=""))


## PC3 - ?

sp <- split(pca$PC3, pca$study)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(study = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p3<-ggplot(pca, aes(x=PC3,color=study, fill=study)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = study, vjust = - 0.5)) + 
  ggtitle("Country") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())

pdf(paste(path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC3.pdf",sep=""),width=14,height=5)
p3
dev.off()

system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC3.pdf ~/tmp_plots/",sep=""))


## PC4 East to west

sp <- split(pca$PC4, pca$study)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(study = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p4<-ggplot(pca, aes(x=PC4,color=study, fill=study)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = study, vjust = - 0.5))+ 
  ggtitle("Country") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())

pdf(paste(path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC4.pdf",sep=""),width=14,height=5)
p4
dev.off()

system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC4.pdf ~/tmp_plots/",sep=""))


## PC5  _ finish vs non finish?

# by Country:

sp <- split(pca$PC5, pca$study)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(study = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p1<-ggplot(pca, aes(x=PC5,color=study, fill=study)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = study, vjust = - 0.5)) + 
  ggtitle("Country") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())


# by Array:

sp <- split(pca$PC5, pca$array)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(array = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p2<-ggplot(pca, aes(x=PC5,color=array, fill=array)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = array, vjust = - 0.5)) + 
  ggtitle("Array") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())


# by Phenotype:
sp <- split(pca$PC5, pca$pheno)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(pheno = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p3<-ggplot(pca, aes(x=PC5,color=pheno, fill=pheno)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = pheno, vjust = - 0.5)) + 
  ggtitle("Phenotype") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())



# by Cohort:

sp <- split(pca$PC5, pca$cohort)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(cohort = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p4<-ggplot(pca, aes(x=PC5,color=cohort, fill=cohort)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = cohort, vjust = - 0.5)) + 
  ggtitle("Cohort") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())

pdf(paste(path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC5.pdf",sep=""),width=14,height=20)
ggarrange(p1,p2,p3,p4,nrow=4,heights=c(1,1,0.8,1.2))
dev.off()

system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC5.pdf ~/tmp_plots/",sep=""))



## PC6  Italy Finland

# by Country:

sp <- split(pca$PC6, pca$study)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(study = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p1<-ggplot(pca, aes(x=PC6,color=study, fill=study)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = study, vjust = - 0.5)) + 
  ggtitle("Country") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())


# by Array:

sp <- split(pca$PC6, pca$array)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(array = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p2<-ggplot(pca, aes(x=PC6,color=array, fill=array)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = array, vjust = - 0.5)) + 
  ggtitle("Array") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())


# by Phenotype:
sp <- split(pca$PC6, pca$pheno)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(pheno = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p3<-ggplot(pca, aes(x=PC6,color=pheno, fill=pheno)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = pheno, vjust = - 0.5)) + 
  ggtitle("Phenotype") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())



# by Cohort:

sp <- split(pca$PC6, pca$cohort)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(cohort = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p4<-ggplot(pca, aes(x=PC6,color=cohort, fill=cohort)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = cohort, vjust = - 0.5)) + 
  ggtitle("Cohort") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())

pdf(paste(path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC6.pdf",sep=""),width=14,height=20)
ggarrange(p1,p2,p3,p4,nrow=4,heights=c(1,1,0.8,1.2))
dev.off()

system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC6.pdf ~/tmp_plots/",sep=""))

## PC7  - Finland

# by Country:

sp <- split(pca$PC7, pca$study)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(study = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p1<-ggplot(pca, aes(x=PC7,color=study, fill=study)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = study, vjust = - 0.5)) + 
  ggtitle("Country") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())


# by Array:

sp <- split(pca$PC7, pca$array)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(array = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p2<-ggplot(pca, aes(x=PC7,color=array, fill=array)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = array, vjust = - 0.5)) + 
  ggtitle("Array") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())


# by Phenotype:
sp <- split(pca$PC7, pca$pheno)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(pheno = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p3<-ggplot(pca, aes(x=PC7,color=pheno, fill=pheno)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = pheno, vjust = - 0.5)) + 
  ggtitle("Phenotype") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())


# by Cohort:

sp <- split(pca$PC7, pca$cohort)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(cohort = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p4<-ggplot(pca, aes(x=PC7,color=cohort, fill=cohort)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = cohort, vjust = - 0.5)) + 
  ggtitle("Cohort") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())

pdf(paste(path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC7.pdf",sep=""),width=14,height=20)
ggarrange(p1,p2,p3,p4,nrow=4,heights=c(1,1,0.8,1.2))
dev.off()

system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC7.pdf ~/tmp_plots/",sep=""))

## PC8  - Basques

# by Country:

sp <- split(pca$PC8, pca$study)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(study = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p1<-ggplot(pca, aes(x=PC8,color=study, fill=study)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = study, vjust = - 0.5)) + 
  ggtitle("Country") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())


# by Array:

sp <- split(pca$PC8, pca$array)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(array = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p2<-ggplot(pca, aes(x=PC8,color=array, fill=array)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = array, vjust = - 0.5)) + 
  ggtitle("Array") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())


# by Phenotype:
sp <- split(pca$PC8, pca$pheno)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(pheno = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p3<-ggplot(pca, aes(x=PC8,color=pheno, fill=pheno)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = pheno, vjust = - 0.5)) + 
  ggtitle("Phenotype") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())


# by Cohort:

sp <- split(pca$PC8, pca$cohort)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(cohort = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p4<-ggplot(pca, aes(x=PC8,color=cohort, fill=cohort)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = cohort, vjust = - 0.5)) + 
  ggtitle("Cohort") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())

pdf(paste(path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC8.pdf",sep=""),width=14,height=20)
ggarrange(p1,p2,p3,p4,nrow=4,heights=c(1,1,0.8,1.2))
dev.off()

system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC8.pdf ~/tmp_plots/",sep=""))



## PC9  ?

# by Country:

sp <- split(pca$PC9, pca$study)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(study = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p1<-ggplot(pca, aes(x=PC9,color=study, fill=study)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = study, vjust = - 0.5)) + 
  ggtitle("Country") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())


# by Array:

sp <- split(pca$PC9, pca$array)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(array = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p2<-ggplot(pca, aes(x=PC9,color=array, fill=array)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = array, vjust = - 0.5)) + 
  ggtitle("Array") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())


# by Phenotype:
sp <- split(pca$PC9, pca$pheno)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(pheno = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p3<-ggplot(pca, aes(x=PC9,color=pheno, fill=pheno)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = pheno, vjust = - 0.5)) + 
  ggtitle("Phenotype") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())


# by Cohort:

sp <- split(pca$PC9, pca$cohort)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(cohort = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p4<-ggplot(pca, aes(x=PC9,color=cohort, fill=cohort)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = cohort, vjust = - 0.5)) + 
  ggtitle("Cohort") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())

pdf(paste(path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC9.pdf",sep=""),width=14,height=20)
ggarrange(p1,p2,p3,p4,nrow=4,heights=c(1,1,0.8,1.2))
dev.off()

system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC9.pdf ~/tmp_plots/",sep=""))


## PC10  

# by Country:

sp <- split(pca$PC10, pca$study)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(study = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p1<-ggplot(pca, aes(x=PC10,color=study, fill=study)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = study, vjust = - 0.5)) + 
  ggtitle("Country") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())


# by Array:

sp <- split(pca$PC10, pca$array)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(array = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p2<-ggplot(pca, aes(x=PC10,color=array, fill=array)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = array, vjust = - 0.5)) + 
  ggtitle("Array") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())


# by Phenotype:
sp <- split(pca$PC10, pca$pheno)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(pheno = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p3<-ggplot(pca, aes(x=PC10,color=pheno, fill=pheno)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = pheno, vjust = - 0.5)) + 
  ggtitle("Phenotype") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())


# by Cohort:

sp <- split(pca$PC10, pca$cohort)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]])
  k <- which.max(d$y)
  data.frame(cohort = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

p4<-ggplot(pca, aes(x=PC10,color=cohort, fill=cohort)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = cohort, vjust = - 0.5)) + 
  ggtitle("Cohort") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.1) + theme(legend.title = element_blank())

pdf(paste(path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC10.pdf",sep=""),width=14,height=20)
ggarrange(p1,p2,p3,p4,nrow=4,heights=c(1,1,0.8,1.2))
dev.off()

system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC10.pdf ~/tmp_plots/",sep=""))





###### what drives each PC?

wg<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_EUR_nodup_ref_pcs.eigenvec.allele",sep=""),head=F)

for (i in c(1:10)) {
  print(i)
  print(summary(wg[,i+5]))
}



fd=function(x) {
  n=length(x)
  r=IQR(x)
  2*r/n^(1/3)
}

fd_bin<-fd(wg[,6])

p1<-ggplot(wg, aes(x=V6)) +
  geom_histogram(position="identity",binwidth=fd_bin) + ggtitle("PC1")
p2<-ggplot(wg, aes(x=V7)) +
  geom_histogram(position="identity",binwidth=fd_bin) + ggtitle("PC2")
p3<-ggplot(wg, aes(x=V8)) +
  geom_histogram(position="identity",binwidth=fd_bin) + ggtitle("PC3")
p4<-ggplot(wg, aes(x=V9)) +
  geom_histogram(position="identity",binwidth=fd_bin) + ggtitle("PC4")
p5<-ggplot(wg, aes(x=V10)) +
  geom_histogram(position="identity",binwidth=fd_bin) + ggtitle("PC5")
p6<-ggplot(wg, aes(x=V11)) +
  geom_histogram(position="identity",binwidth=fd_bin) + ggtitle("PC6")
p7<-ggplot(wg, aes(x=V12)) +
  geom_histogram(position="identity",binwidth=fd_bin) + ggtitle("PC7")
p8<-ggplot(wg, aes(x=V13)) +
  geom_histogram(position="identity",binwidth=fd_bin) + ggtitle("PC8")
p9<-ggplot(wg, aes(x=V14)) +
  geom_histogram(position="identity",binwidth=fd_bin) + ggtitle("PC9")
p10<-ggplot(wg, aes(x=V15)) +
  geom_histogram(position="identity",binwidth=fd_bin) + ggtitle("PC10")

  
pdf(paste(path,"pre_imputation/QC/pca_1000gp/histogram_weights_pcs.pdf",sep=""),width=7,height=50)
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,nrow=10)
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/histogram_weights_pcs.pdf ~/tmp_plots/",sep=""))

# FROM PC5 ONWARDS, EFFECT DRIVEN BY A COUPLE OF VARIANTS OUTLIERS:
frq<-read.table("/path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/iibdgc_merged_noHighLD_pruned_noassoclogP4_EURonly_noDuplicates.afreq",head=F)
colnames(frq)<-c("CHROM","ID","REF","ALT","ALT_FREQS","OBS_CT")

wg<-merge(wg,frq,by.x="V2",by.y="ID")


for (i in c(1:10)) {
  print(i)
  print(summary(wg[,5+i]))
}

PC5
wg[which( wg$V10> 5),1]
11:56294275_T_C rs658845 GBR 0.18, FIN 0.32 IBD 0.21
11:56869386_C_T rs527174 GBR 0.38 FIN 0.54 IBS 0.42
11:56910484_C_A rs648409 GBR 0.25 FIN 0.14 IBS 0.25
11:57120047_A_G rs3781902

PC6
wg[which( wg$V11> 5),1]
1:4980692_G_A rs12035499 GBR 0.2 IBD 0.19 FIN 0.14
1:4987424_T_G 
1:4996185_A_G

PC7
wg[which( wg$V12> 5),1]
7:28494445_G_A
7:28523442_G_A 
7:28545611_C_T
8:104171886_C_T
8:104196713_T_C
8:104249839_A_G

PC8
wg[which( wg$V13> 5),1]
7:28494445_G_A
7:28523442_G_A
7:28545611_C_T
8:104171886_C_T
8:104196713_T_C
8:104249839_A_G

PC9
wg[which( wg$V14> 5),1]
6:38728142_A_C
6:38953292_A_G
6:38988026_A_G

PC10
wg[which( wg$V15> 5),1]
6:56268377_G_A rs9475690 EUR 0.33 EAS 0.09
6:56281321_T_C
6:56507489_A_G

# 11 56294275 56294275 T C
# 11 56869386 56869386 C T
# 11 56910484 56910484 C A
# 11 57120047 57120047 A G
# 1 4980692 4980692 G A 
# 1 4987424 4987424 T G
# 1 4996185 4996185 A G
# 7 28494445 28494445 G A
# 7 28523442 28523442 G A 
# 7 28545611 28545611 C T
# 8 104171886 104171886 C T
# 8 104196713 104196713 T C
# 8 104249839 104249839 A G
# 7 28494445 28494445 G A
# 7 28523442 28523442 G A
# 7 28545611 28545611 C T
# 8 104171886 104171886 C T
# 8 104196713 104196713 T C
# 8 104249839 104249839 A G
# 6 38728142 38728142 A C
# 6 38953292 38953292 A G
# 6 38988026 38988026 A G
# 6 56268377 56268377 G A
# 6 56281321 56281321 T C
# 6 56507489 56507489 A G

# where outliers for PC9 and PC10 plot in all ancestries eth

pca2<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/1000GP_all_b37_plus_iibdgc_subset_iibdgc_new_projection_edited_withDuplicates.sscore",sep=""),head=T,sep="\t")
dim(pca2)

pc9<-pca[which( (pca$PC9 >= (mean(pca$PC9)+4*sd(pca$PC9)) ) | (pca$PC9 <= (mean(pca$PC9)-4*sd(pca$PC9)) ) ),]

pca2.0<-pca2[which(pca2$super_pop!="IIBDGC"),]
pca2.1<-pca2[which(pca2$IID %in% pc9$IID),]
pca2.1<-rbind(pca2.0,pca2.1)


# PC1 - PC2

ylims<-c(min(pca2.1$PC1),max(pca2.1$PC1))
xlims<-c(min(pca2.1$PC2),max(pca2.1$PC2))

p0<-ggplot(pca2.1[which(pca2.1$cohort=="1000GP"),],aes(PC2,PC1)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims)

p1<-ggplot(pca2.1[which(pca2.1$cohort!="1000GP"),],aes(PC2,PC1)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims)


pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_EURPC9_outliers.pdf",sep=""),width=14,height=6)
ggarrange(p0,p1,legend=c("right"),widths = c(1,1.1))
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_EURPC9_outliers.pdf ~/tmp_plots/",sep=""))

# PC3 - PC4

ylims<-c(min(pca2.1$PC3),max(pca2.1$PC3))
xlims<-c(min(pca2.1$PC4),max(pca2.1$PC4))

p0<-ggplot(pca2.1[which(pca2.1$cohort=="1000GP"),],aes(PC4,PC3)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims)

p1<-ggplot(pca2.1[which(pca2.1$cohort!="1000GP"),],aes(PC4,PC3)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims)


pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC3_PC4_EURPC9_outliers.pdf",sep=""),width=14,height=6)
ggarrange(p0,p1,legend=c("right"),widths = c(1,1.1))
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC3_PC4_EURPC9_outliers.pdf ~/tmp_plots/",sep=""))


# PC5 - PC6

ylims<-c(min(pca2.1$PC5),max(pca2.1$PC5))
xlims<-c(min(pca2.1$PC6),max(pca2.1$PC6))

p0<-ggplot(pca2.1[which(pca2.1$cohort=="1000GP"),],aes(PC6,PC5)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims)

p1<-ggplot(pca2.1[which(pca2.1$cohort!="1000GP"),],aes(PC6,PC5)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims)


pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC5_PC6_EURPC9_outliers.pdf",sep=""),width=14,height=6)
ggarrange(p0,p1,legend=c("right"),widths = c(1,1.1))
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC5_PC6_EURPC9_outliers.pdf ~/tmp_plots/",sep=""))





# PLOT ALL PCS COLOUR BY PHENOTYPE

pna<-qplot()+theme(
  panel.background = element_rect(fill = "transparent") # bg of the panel
  , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
  , panel.grid.major = element_blank() # get rid of major grid
  , panel.grid.minor = element_blank() # get rid of minor grid
  , legend.background = element_rect(fill = "transparent") # get rid of legend bg
  , axis.title=element_blank()
  , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
)


p11<-qplot(PC1,PC1, data = pca, colour = pheno)
p12<-qplot(PC2,PC1, data = pca, colour = pheno)
p13<-qplot(PC3,PC1, data = pca, colour = pheno)
p14<-qplot(PC4,PC1, data = pca, colour = pheno)
p15<-qplot(PC5,PC1, data = pca, colour = pheno)
p16<-qplot(PC6,PC1, data = pca, colour = pheno)
p17<-qplot(PC7,PC1, data = pca, colour = pheno)
p18<-qplot(PC8,PC1, data = pca, colour = pheno)
p19<-qplot(PC9,PC1, data = pca, colour = pheno)
p110<-qplot(PC10,PC1, data = pca, colour = pheno)

p21<-qplot(PC1,PC2, data = pca, colour = pheno)
p22<-qplot(PC2,PC2, data = pca, colour = pheno)
p23<-qplot(PC3,PC2, data = pca, colour = pheno)
p24<-qplot(PC4,PC2, data = pca, colour = pheno)
p25<-qplot(PC5,PC2, data = pca, colour = pheno)
p26<-qplot(PC6,PC2, data = pca, colour = pheno)
p27<-qplot(PC7,PC2, data = pca, colour = pheno)
p28<-qplot(PC8,PC2, data = pca, colour = pheno)
p29<-qplot(PC9,PC2, data = pca, colour = pheno)
p210<-qplot(PC10,PC2, data = pca, colour = pheno)

p31<-qplot(PC1,PC3, data = pca, colour = pheno)
p32<-qplot(PC2,PC3, data = pca, colour = pheno)
p33<-qplot(PC3,PC3, data = pca, colour = pheno)
p34<-qplot(PC4,PC3, data = pca, colour = pheno)
p35<-qplot(PC5,PC3, data = pca, colour = pheno)
p36<-qplot(PC6,PC3, data = pca, colour = pheno)
p37<-qplot(PC7,PC3, data = pca, colour = pheno)
p38<-qplot(PC8,PC3, data = pca, colour = pheno)
p39<-qplot(PC9,PC3, data = pca, colour = pheno)
p310<-qplot(PC10,PC3, data = pca, colour = pheno)

p41<-qplot(PC1,PC4, data = pca, colour = pheno)
p42<-qplot(PC2,PC4, data = pca, colour = pheno)
p43<-qplot(PC3,PC4, data = pca, colour = pheno)
p44<-qplot(PC4,PC4, data = pca, colour = pheno)
p45<-qplot(PC5,PC4, data = pca, colour = pheno)
p46<-qplot(PC6,PC4, data = pca, colour = pheno)
p47<-qplot(PC7,PC4, data = pca, colour = pheno)
p48<-qplot(PC8,PC4, data = pca, colour = pheno)
p49<-qplot(PC9,PC4, data = pca, colour = pheno)
p410<-qplot(PC10,PC4, data = pca, colour = pheno)

p51<-qplot(PC1,PC5, data = pca, colour = pheno)
p52<-qplot(PC2,PC5, data = pca, colour = pheno)
p53<-qplot(PC3,PC5, data = pca, colour = pheno)
p54<-qplot(PC4,PC5, data = pca, colour = pheno)
p55<-qplot(PC5,PC5, data = pca, colour = pheno)
p56<-qplot(PC6,PC5, data = pca, colour = pheno)
p57<-qplot(PC7,PC5, data = pca, colour = pheno)
p58<-qplot(PC8,PC5, data = pca, colour = pheno)
p59<-qplot(PC9,PC5, data = pca, colour = pheno)
p510<-qplot(PC10,PC5, data = pca, colour = pheno)

p61<-qplot(PC1,PC6, data = pca, colour = pheno)
p62<-qplot(PC2,PC6, data = pca, colour = pheno)
p63<-qplot(PC3,PC6, data = pca, colour = pheno)
p64<-qplot(PC4,PC6, data = pca, colour = pheno)
p65<-qplot(PC5,PC6, data = pca, colour = pheno)
p66<-qplot(PC6,PC6, data = pca, colour = pheno)
p67<-qplot(PC7,PC6, data = pca, colour = pheno)
p68<-qplot(PC8,PC6, data = pca, colour = pheno)
p69<-qplot(PC9,PC6, data = pca, colour = pheno)
p610<-qplot(PC10,PC6, data = pca, colour = pheno)

p61<-qplot(PC1,PC6, data = pca, colour = pheno)
p62<-qplot(PC2,PC6, data = pca, colour = pheno)
p63<-qplot(PC3,PC6, data = pca, colour = pheno)
p64<-qplot(PC4,PC6, data = pca, colour = pheno)
p65<-qplot(PC5,PC6, data = pca, colour = pheno)
p66<-qplot(PC6,PC6, data = pca, colour = pheno)
p67<-qplot(PC7,PC6, data = pca, colour = pheno)
p68<-qplot(PC8,PC6, data = pca, colour = pheno)
p69<-qplot(PC9,PC6, data = pca, colour = pheno)
p610<-qplot(PC10,PC6, data = pca, colour = pheno)

p71<-qplot(PC1,PC7, data = pca, colour = pheno)
p72<-qplot(PC2,PC7, data = pca, colour = pheno)
p73<-qplot(PC3,PC7, data = pca, colour = pheno)
p74<-qplot(PC4,PC7, data = pca, colour = pheno)
p75<-qplot(PC5,PC7, data = pca, colour = pheno)
p76<-qplot(PC6,PC7, data = pca, colour = pheno)
p77<-qplot(PC7,PC7, data = pca, colour = pheno)
p78<-qplot(PC8,PC7, data = pca, colour = pheno)
p79<-qplot(PC9,PC7, data = pca, colour = pheno)
p710<-qplot(PC10,PC7, data = pca, colour = pheno)

p81<-qplot(PC1,PC8, data = pca, colour = pheno)
p82<-qplot(PC2,PC8, data = pca, colour = pheno)
p83<-qplot(PC3,PC8, data = pca, colour = pheno)
p84<-qplot(PC4,PC8, data = pca, colour = pheno)
p85<-qplot(PC5,PC8, data = pca, colour = pheno)
p86<-qplot(PC6,PC8, data = pca, colour = pheno)
p87<-qplot(PC7,PC8, data = pca, colour = pheno)
p88<-qplot(PC8,PC8, data = pca, colour = pheno)
p89<-qplot(PC9,PC8, data = pca, colour = pheno)
p810<-qplot(PC10,PC8, data = pca, colour = pheno)

p91<-qplot(PC1,PC9, data = pca, colour = pheno)
p92<-qplot(PC2,PC9, data = pca, colour = pheno)
p93<-qplot(PC3,PC9, data = pca, colour = pheno)
p94<-qplot(PC4,PC9, data = pca, colour = pheno)
p95<-qplot(PC5,PC9, data = pca, colour = pheno)
p96<-qplot(PC6,PC9, data = pca, colour = pheno)
p97<-qplot(PC7,PC9, data = pca, colour = pheno)
p98<-qplot(PC8,PC9, data = pca, colour = pheno)
p99<-qplot(PC9,PC9, data = pca, colour = pheno)
p910<-qplot(PC10,PC9, data = pca, colour = pheno)

p101<-qplot(PC1,PC10, data = pca, colour = pheno)
p102<-qplot(PC2,PC10, data = pca, colour = pheno)
p103<-qplot(PC3,PC10, data = pca, colour = pheno)
p104<-qplot(PC4,PC10, data = pca, colour = pheno)
p105<-qplot(PC5,PC10, data = pca, colour = pheno)
p106<-qplot(PC6,PC10, data = pca, colour = pheno)
p107<-qplot(PC7,PC10, data = pca, colour = pheno)
p108<-qplot(PC8,PC10, data = pca, colour = pheno)
p109<-qplot(PC9,PC10, data = pca, colour = pheno)
p1010<-qplot(PC10,PC10, data = pca, colour = pheno)


r1<-ggarrange(pna,p12,p13,p14,p15,p16,p17,p18,p19,p110,ncol=10,legend=c("none"))
r2<-ggarrange(p21,pna,p23,p24,p25,p26,p27,p28,p29,p210,ncol=10,legend=c("none"))
r3<-ggarrange(p31,p32,pna,p34,p35,p36,p37,p38,p39,p310,ncol=10,legend=c("none"))
r4<-ggarrange(p41,p42,p43,pna,p45,p46,p47,p48,p49,p410,ncol=10,legend=c("none"))
r5<-ggarrange(p51,p52,p53,p54,pna,p56,p57,p58,p59,p510,ncol=10,legend=c("none"))
r6<-ggarrange(p61,p62,p63,p64,p65,pna,p67,p68,p69,p610,ncol=10,legend=c("none"))
r7<-ggarrange(p71,p72,p73,p74,p75,p76,pna,p78,p79,p710,ncol=10,legend=c("none"))
r8<-ggarrange(p81,p82,p83,p84,p85,p86,p87,pna,p89,p810,ncol=10,legend=c("none"))
r9<-ggarrange(p91,p92,p93,p94,p95,p96,p97,p98,pna,p910,ncol=10,legend=c("none"))
r10<-ggarrange(p101,p102,p103,p104,p105,p106,p107,p108,p109,pna,ncol=10,common.legend = TRUE,legend=c("bottom"))


pdf(paste(path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_PC7_PC8_PC9_PC10_3kvariants_iibdgc_EURonly_colour_by_pheno.pdf",sep=""),width=45,height=45)
ggarrange(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,nrow=10,common.legend = TRUE,legend=c("bottom"),heights = c(1,1,1,1,1,1,1,1,1,1.2))
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_PC7_PC8_PC9_PC10_3kvariants_iibdgc_EURonly_colour_by_pheno.pdf ~/tmp_plots/",sep=""))


# PLOT ALL PCS COLOUR BY ARRAY:


# create plots PC1 to PC6

pna<-qplot()+theme(
  panel.background = element_rect(fill = "transparent") # bg of the panel
  , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
  , panel.grid.major = element_blank() # get rid of major grid
  , panel.grid.minor = element_blank() # get rid of minor grid
  , legend.background = element_rect(fill = "transparent") # get rid of legend bg
  , axis.title=element_blank()
  , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
)


p11<-qplot(PC1,PC1, data = pca, colour = array)
p12<-qplot(PC2,PC1, data = pca, colour = array)
p13<-qplot(PC3,PC1, data = pca, colour = array)
p14<-qplot(PC4,PC1, data = pca, colour = array)
p15<-qplot(PC5,PC1, data = pca, colour = array)
p16<-qplot(PC6,PC1, data = pca, colour = array)
p17<-qplot(PC7,PC1, data = pca, colour = array)
p18<-qplot(PC8,PC1, data = pca, colour = array)
p19<-qplot(PC9,PC1, data = pca, colour = array)
p110<-qplot(PC10,PC1, data = pca, colour = array)

p21<-qplot(PC1,PC2, data = pca, colour = array)
p22<-qplot(PC2,PC2, data = pca, colour = array)
p23<-qplot(PC3,PC2, data = pca, colour = array)
p24<-qplot(PC4,PC2, data = pca, colour = array)
p25<-qplot(PC5,PC2, data = pca, colour = array)
p26<-qplot(PC6,PC2, data = pca, colour = array)
p27<-qplot(PC7,PC2, data = pca, colour = array)
p28<-qplot(PC8,PC2, data = pca, colour = array)
p29<-qplot(PC9,PC2, data = pca, colour = array)
p210<-qplot(PC10,PC2, data = pca, colour = array)

p31<-qplot(PC1,PC3, data = pca, colour = array)
p32<-qplot(PC2,PC3, data = pca, colour = array)
p33<-qplot(PC3,PC3, data = pca, colour = array)
p34<-qplot(PC4,PC3, data = pca, colour = array)
p35<-qplot(PC5,PC3, data = pca, colour = array)
p36<-qplot(PC6,PC3, data = pca, colour = array)
p37<-qplot(PC7,PC3, data = pca, colour = array)
p38<-qplot(PC8,PC3, data = pca, colour = array)
p39<-qplot(PC9,PC3, data = pca, colour = array)
p310<-qplot(PC10,PC3, data = pca, colour = array)

p41<-qplot(PC1,PC4, data = pca, colour = array)
p42<-qplot(PC2,PC4, data = pca, colour = array)
p43<-qplot(PC3,PC4, data = pca, colour = array)
p44<-qplot(PC4,PC4, data = pca, colour = array)
p45<-qplot(PC5,PC4, data = pca, colour = array)
p46<-qplot(PC6,PC4, data = pca, colour = array)
p47<-qplot(PC7,PC4, data = pca, colour = array)
p48<-qplot(PC8,PC4, data = pca, colour = array)
p49<-qplot(PC9,PC4, data = pca, colour = array)
p410<-qplot(PC10,PC4, data = pca, colour = array)

p51<-qplot(PC1,PC5, data = pca, colour = array)
p52<-qplot(PC2,PC5, data = pca, colour = array)
p53<-qplot(PC3,PC5, data = pca, colour = array)
p54<-qplot(PC4,PC5, data = pca, colour = array)
p55<-qplot(PC5,PC5, data = pca, colour = array)
p56<-qplot(PC6,PC5, data = pca, colour = array)
p57<-qplot(PC7,PC5, data = pca, colour = array)
p58<-qplot(PC8,PC5, data = pca, colour = array)
p59<-qplot(PC9,PC5, data = pca, colour = array)
p510<-qplot(PC10,PC5, data = pca, colour = array)

p61<-qplot(PC1,PC6, data = pca, colour = array)
p62<-qplot(PC2,PC6, data = pca, colour = array)
p63<-qplot(PC3,PC6, data = pca, colour = array)
p64<-qplot(PC4,PC6, data = pca, colour = array)
p65<-qplot(PC5,PC6, data = pca, colour = array)
p66<-qplot(PC6,PC6, data = pca, colour = array)
p67<-qplot(PC7,PC6, data = pca, colour = array)
p68<-qplot(PC8,PC6, data = pca, colour = array)
p69<-qplot(PC9,PC6, data = pca, colour = array)
p610<-qplot(PC10,PC6, data = pca, colour = array)

p61<-qplot(PC1,PC6, data = pca, colour = array)
p62<-qplot(PC2,PC6, data = pca, colour = array)
p63<-qplot(PC3,PC6, data = pca, colour = array)
p64<-qplot(PC4,PC6, data = pca, colour = array)
p65<-qplot(PC5,PC6, data = pca, colour = array)
p66<-qplot(PC6,PC6, data = pca, colour = array)
p67<-qplot(PC7,PC6, data = pca, colour = array)
p68<-qplot(PC8,PC6, data = pca, colour = array)
p69<-qplot(PC9,PC6, data = pca, colour = array)
p610<-qplot(PC10,PC6, data = pca, colour = array)

p71<-qplot(PC1,PC7, data = pca, colour = array)
p72<-qplot(PC2,PC7, data = pca, colour = array)
p73<-qplot(PC3,PC7, data = pca, colour = array)
p74<-qplot(PC4,PC7, data = pca, colour = array)
p75<-qplot(PC5,PC7, data = pca, colour = array)
p76<-qplot(PC6,PC7, data = pca, colour = array)
p77<-qplot(PC7,PC7, data = pca, colour = array)
p78<-qplot(PC8,PC7, data = pca, colour = array)
p79<-qplot(PC9,PC7, data = pca, colour = array)
p710<-qplot(PC10,PC7, data = pca, colour = array)

p81<-qplot(PC1,PC8, data = pca, colour = array)
p82<-qplot(PC2,PC8, data = pca, colour = array)
p83<-qplot(PC3,PC8, data = pca, colour = array)
p84<-qplot(PC4,PC8, data = pca, colour = array)
p85<-qplot(PC5,PC8, data = pca, colour = array)
p86<-qplot(PC6,PC8, data = pca, colour = array)
p87<-qplot(PC7,PC8, data = pca, colour = array)
p88<-qplot(PC8,PC8, data = pca, colour = array)
p89<-qplot(PC9,PC8, data = pca, colour = array)
p810<-qplot(PC10,PC8, data = pca, colour = array)

p91<-qplot(PC1,PC9, data = pca, colour = array)
p92<-qplot(PC2,PC9, data = pca, colour = array)
p93<-qplot(PC3,PC9, data = pca, colour = array)
p94<-qplot(PC4,PC9, data = pca, colour = array)
p95<-qplot(PC5,PC9, data = pca, colour = array)
p96<-qplot(PC6,PC9, data = pca, colour = array)
p97<-qplot(PC7,PC9, data = pca, colour = array)
p98<-qplot(PC8,PC9, data = pca, colour = array)
p99<-qplot(PC9,PC9, data = pca, colour = array)
p910<-qplot(PC10,PC9, data = pca, colour = array)

p101<-qplot(PC1,PC10, data = pca, colour = array)
p102<-qplot(PC2,PC10, data = pca, colour = array)
p103<-qplot(PC3,PC10, data = pca, colour = array)
p104<-qplot(PC4,PC10, data = pca, colour = array)
p105<-qplot(PC5,PC10, data = pca, colour = array)
p106<-qplot(PC6,PC10, data = pca, colour = array)
p107<-qplot(PC7,PC10, data = pca, colour = array)
p108<-qplot(PC8,PC10, data = pca, colour = array)
p109<-qplot(PC9,PC10, data = pca, colour = array)
p1010<-qplot(PC10,PC10, data = pca, colour = array)


r1<-ggarrange(pna,p12,p13,p14,p15,p16,p17,p18,p19,p110,ncol=10,legend=c("none"))
r2<-ggarrange(p21,pna,p23,p24,p25,p26,p27,p28,p29,p210,ncol=10,legend=c("none"))
r3<-ggarrange(p31,p32,pna,p34,p35,p36,p37,p38,p39,p310,ncol=10,legend=c("none"))
r4<-ggarrange(p41,p42,p43,pna,p45,p46,p47,p48,p49,p410,ncol=10,legend=c("none"))
r5<-ggarrange(p51,p52,p53,p54,pna,p56,p57,p58,p59,p510,ncol=10,legend=c("none"))
r6<-ggarrange(p61,p62,p63,p64,p65,pna,p67,p68,p69,p610,ncol=10,legend=c("none"))
r7<-ggarrange(p71,p72,p73,p74,p75,p76,pna,p78,p79,p710,ncol=10,legend=c("none"))
r8<-ggarrange(p81,p82,p83,p84,p85,p86,p87,pna,p89,p810,ncol=10,legend=c("none"))
r9<-ggarrange(p91,p92,p93,p94,p95,p96,p97,p98,pna,p910,ncol=10,legend=c("none"))
r10<-ggarrange(p101,p102,p103,p104,p105,p106,p107,p108,p109,pna,ncol=10,common.legend = TRUE,legend=c("bottom"))


pdf(paste(path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_PC7_PC8_PC9_PC10_3kvariants_iibdgc_EURonly_colour_by_array.pdf",sep=""),width=45,height=45)
ggarrange(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,nrow=10,common.legend = TRUE,legend=c("bottom"),heights = c(1,1,1,1,1,1,1,1,1,1.2))
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_PC7_PC8_PC9_PC10_3kvariants_iibdgc_EURonly_colour_by_array.pdf ~/tmp_plots/",sep=""))





#############################################################################################################################################################
# CREATE AN UNIQUE PCA FILE:

pca_eur<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_EUR_withDuplicates_ref_pcs_new_projection_edited.sscore",sep=""),head=T,sep="\t")
dim(pca_eur)
# [1] 113112    19

pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/1000GP_all_b37_plus_iibdgc_subset_iibdgc_new_projection_edited_withDuplicates.sscore",sep=""),head=T,sep="\t")
dim(pca)
# [1] 120427     21

colnames(pca_eur)[6:15]<-paste(colnames(pca_eur)[6:15],"_EUR",sep="")

pca<-merge(pca,pca_eur[,c(1,6:15,18,20)],by="FID",all.x=T)
pca<-pca[,c(1:3,6:32)]

table(pca$super_pop)
# AFR    AMR    EAS    EUR IIBDGC    SAS 
# 661    347    504    503 117923    489

write.table(pca,paste(path,"pre_imputation/QC/pca_1000gp/pca_1000gp_iibdgc_pca_pcaEur_withDuplicates.tsv",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")


#### create lists with no duplicates:

fam<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_nodup.fam",sep=""),head=F)

pca_nodup<-pca[which(pca$IID %in% fam$V2),]
pca_eur_nodup<-pca_nodup[which(pca_nodup$inferred_population=="EUR"),]
dim(pca_eur_nodup)
# [1] 102134     30

write.table(pca_eur_nodup[which(pca_eur_nodup$pca_jewish=="Jewish"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eur_jewish_ancestry_samples_noDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")
write.table(pca_eur_nodup[which(pca_eur_nodup$pca_jewish=="Non-Jewish"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eur_nonjewish_ancestry_samples_noDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")


table(pca_eur_nodup$cohort[which(pca_eur_nodup$PHENO1=="2")],pca_eur_nodup$pca_jewish[which(pca_eur_nodup$PHENO1=="2")],useNA="ifany")
#                         Jewish Non-Jewish
# 1000GP                       0          0
# all_hce                    130      10866
# australia_omniexome         12        642
# basque_gsa                   1        517
# belgium_franchimont_gsa     45        827
# belgium_inf1_old_gwas        5        223
# belgium_inf2_old_gwas        4         78
# belgium_louis_gsa           20        866
# belgium_vermeire_gsa        72       3013
# ccfa_gsa                   140       1681
# cedars_370k_old_gwas       140        255
# cedars_610k_old_gwas       171        269
# cedars_gsa                 183        744
# cedars_omni_old_gwas       221        379
# chop_old_gwas              158       1549
# finland_illugwas             1        440
# german_affy6_old_gwas        2        951
# gwas1                       44       1696
# gwas2                        1       2347
# italy_gsa                  118        455
# kiel_austria_sibdcs_gsa    167       9043
# lithuania_gsa                6       1049
# mccauley_gsa               283        464
# netherlands_gsa             42       3752
# niddk_broad_gsa            716       3143
# niddk_feinstein_gsa       1917       3137
# niddk_old_gwas             251        355
# norway_affy6_old_gwas        0        266
# pittsburgh_gsa              13        637
# prism_nfe_gsa               65        334
# prism_nfe_gwas              54        381
# slovenia_gsa                 0         87
# spain_gsa                   20       1937
# sweden_gsa                  13        363
# swedish_uc_old_gwas          0        919

table(pca_eur_nodup$cohort[which(pca_eur_nodup$PHENO1=="1")],pca_eur_nodup$pca_jewish[which(pca_eur_nodup$PHENO1=="1")],useNA="ifany")
#                         Jewish Non-Jewish
# 1000GP                       0          0
# all_hce                     46      10299
# australia_omniexome          2        587
# basque_gsa                   6        958
# belgium_franchimont_gsa     25        558
# belgium_inf1_old_gwas       13        818
# belgium_inf2_old_gwas        1         76
# belgium_louis_gsa           13        578
# belgium_vermeire_gsa         4        802
# ccfa_gsa                     0          0
# cedars_370k_old_gwas         0          0
# cedars_610k_old_gwas         0          0
# cedars_gsa                 193        448
# cedars_omni_old_gwas         1          0
# chop_old_gwas              469       5597
# finland_illugwas             0          0
# german_affy6_old_gwas        2       1657
# gwas1                       15       2904
# gwas2                        0       2767
# italy_gsa                   81        289
# kiel_austria_sibdcs_gsa      9       4400
# lithuania_gsa                0       1122
# mccauley_gsa                 0          0
# netherlands_gsa              0        532
# niddk_broad_gsa            151        835
# niddk_feinstein_gsa        151       1584
# niddk_old_gwas             403        263
# norway_affy6_old_gwas        0        281
# pittsburgh_gsa              85       1290
# prism_nfe_gsa                4         24
# prism_nfe_gwas              24        181
# slovenia_gsa                 0        175
# spain_gsa                   14       1465
# sweden_gsa                   0        927
# swedish_uc_old_gwas          0        325





