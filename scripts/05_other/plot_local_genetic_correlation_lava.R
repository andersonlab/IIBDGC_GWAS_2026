# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# singularity exec iibdgc_postprocess_10_singularity.sif

path_gwas="/path/to/ibdgwas/IIBDGC/"

MEM=800
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q normal R 

library(data.table)
library(R.utils)
library(stringr)
library(qvalue)
library(ggplot2)
library(ggpubr)

path_gwas="/path/to/ibdgwas/IIBDGC/"

# retrieve all results, and explore negative ones

list<-list.files(paste0(path_gwas,"post_imputation/2022/analysis/local_genetic_correlation/results/"))


# EXPLORE LOCAL GENETIC HERITABILITY:

lgc<-list[grep("local_heritability_per_phenotype",list)]

length(lgc)
# [1] 419

lgc_reg<-gsub("_local_heritability_per_phenotype.tsv","",lgc)

# reg<-fread(paste0(path_gwas,"post_imputation/2022/analysis/local_genetic_correlation/input_files/regions.txt"),head=F)
# # resubmitted those that fail:
# reg[which(!reg$V1 %in% lgc_reg),] 

for (i in 1:length(lgc)) {
    tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/local_genetic_correlation/results/",lgc[i]),head=T)
    tmp$region<-gsub("_local_heritability_per_phenotype.tsv","",lgc[i])

    if(i==1) {
        lgh_dat<-tmp
    } else {
        lgh_dat<-rbind(lgh_dat,tmp)
    }

    rm(tmp)

}

# keep only pairs region_disease for those with significant heritability:
lgh_dat<-lgh_dat[which(lgh_dat$phen %in% c("cd","uc","ibd")),]
lgh_dat<-lgh_dat[which(lgh_dat$p<0.05),]

table(lgh_dat$region,lgh_dat$phen)                  
#                          cd ibd uc
#   1_11296321_12587358     1   1  1
#   1_113334946_114334947   1   0  0
#   1_119380078_121468248   1   1  0
#   1_149759332_152333344   1   1  1
....

dim(table(lgh_dat$region,lgh_dat$phen)  )
# [1] 290   3


# EXPLORE LOCAL GENETIC CORRELATION:

lgc<-list[grep("genetic_correlation_between_phenotype_pairs.tsv",list)]

length(lgc)
# [1] 408

lgc_reg<-gsub("genetic_correlation_between_phenotype_pairs.tsv","",lgc)

# reg<-fread(paste0(path_gwas,"post_imputation/2022/analysis/local_genetic_correlation/input_files/regions.txt"),head=F)
# # resubmitted those that fail:
# reg[which(!reg$V1 %in% lgc_reg),] 

for (i in 1:length(lgc)) {
    tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/local_genetic_correlation/results/",lgc[i]),head=T)
    tmp$region<-gsub("_local_genetic_correlation_between_phenotype_pairs.tsv","",lgc[i])

    if(i==1) {
        lgc_dat<-tmp
    } else {
        lgc_dat<-rbind(lgc_dat,tmp)
    }

    rm(tmp)

}

table(lgc_dat$phen2)

#    cd   ibd    uc 
# 21817 22924 21316

lgc_dat<-lgc_dat[which(!is.na(lgc_dat$p)),]
lgc_dat<-lgc_dat[which(!lgc_dat$phen1 %in% c("ibd","cd","uc")),]


# relabel to phenotypes the study accessions
dat<-fread(paste(path_gwas,"resources/gwas_summary_statistics/gwas-catalog-v1.0.3.1-studies-r2025-05-13.tsv.gz",sep=""),quote="")
dat$study_accession<-dat$'STUDY ACCESSION'
dat$disease_trait<-dat$'DISEASE/TRAIT'

dat$group<-""
dat$group[which(dat$'PUBMED ID'==38448586)]<-"Metabolic biomarkers"
dat$group[which(dat$'PUBMED ID'==32888493)]<-"Blood cell counts"
dat$group[which(dat$'STUDY ACCESSION' %in% c("GCST90476007","GCST90475936","GCST90475967","GCST90475929","GCST90435254","GCST90475214","GCST90476007"))]<-"Cardiovascular"
dat$group[which(dat$'STUDY ACCESSION' %in% c("GCST004988","GCST90274714","GCST006464","GCST004748","GCST90018921","GCST90129505"))]<-"Cancer"
dat$group[which(dat$'STUDY ACCESSION' %in% c("GCST90027161","GCST005523",
"GCST90044763","GCST90018847","GCST005531","GCST90061440","GCST004030","GCST002318","GCST003156","GCST003566",
"GCST90014023","GCST90018926","GCST005529","GCST005527",
"GCST90132223","GCST011096","GCST90044158","GCST90472771","GCST90480502","GCST003566","GCST90013445","GCST90029017"))]<-"Immune Mediated"
dat$group[which(dat$'STUDY ACCESSION' %in% c("GCST90027158","GCST90018959","GCST90018947","GCST90018926","GCST90165267","GCST90029070"))]<-"Other"

dat1<-dat[,c("study_accession","disease_trait","group","MAPPED_TRAIT")]



# create a similar map for the interval data:

dat2<-fread("/path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/Akbari_2023_supplementary_data_2.txt")

dat2$study_accession<-dat2$user_friendly
# replace some special characters in the name:
dat2$study_accession<-gsub("\\(","_",dat2$study_accession)
dat2$study_accession<-gsub("\\)","_",dat2$study_accession)
dat2$study_accession<-gsub("\\(=","_",dat2$study_accession)
dat2$study_accession<-gsub("-","_",dat2$study_accession)
dat2$study_accession[which(dat2$study_accession=="DELTA_HGB")]<-"Delta_HGB"

dat2$disease_trait<-dat2$user_friendly
dat2$MAPPED_TRAIT<-dat2$'Long name'
dat2$group<-"Blood cell morphology"

dat2<-dat2[,c("study_accession","disease_trait","group","MAPPED_TRAIT","cell")]
dim(dat2[which(dat2$study_accession %in% all$p2),])
# [1] 63  4

dat<-rbind(dat1,dat2[,c("study_accession","disease_trait","group","MAPPED_TRAIT")])
rm(dat1,dat2)


# combine all

lgc_dat<-merge(lgc_dat,dat,by.x="phen1",by.y="study_accession",all.x=T)


## CREATE PLOTS PER GROUPS - similar to rg:

lgc_dat$group[grep("levels",lgc_dat$disease_trait)]<-"Metabolic biomarkers - Other"

lgc_dat$group[grep("levels",lgc_dat$disease_trait)]<-"Metabolic biomarkers - Fluid balance"

lgc_dat$group[grep("Glycoprotein acetylation|C-reactive",lgc_dat$disease_trait)]<-"Metabolic biomarkers - Inflammation"

lgc_dat$group[grep("cholesterol|VLDL particles|IDL particles|LDL particles|HDL particles",lgc_dat$disease_trait)]<-"Metabolic biomarkers - Cholesterol"
lgc_dat$group[grep("Cholesterol",lgc_dat$disease_trait)]<-"Metabolic biomarkers - Cholesterol"
lgc_dat$group[grep("Cholesteryl",lgc_dat$disease_trait)]<-"Metabolic biomarkers - Cholesterol"
lgc_dat$group[grep("Triglyceride",lgc_dat$disease_trait)]<-"Metabolic biomarkers - Lipids"
lgc_dat$group[grep("Total lipid|diacylglycerol",lgc_dat$disease_trait)]<-"Metabolic biomarkers - Lipids"
lgc_dat$group[grep("Phospholipid",lgc_dat$disease_trait)]<-"Metabolic biomarkers - Lipids"
lgc_dat$group[grep("phosphoglycerides",lgc_dat$disease_trait)]<-"Metabolic biomarkers - Lipids"

lgc_dat$group[grep("fatty|unsaturation|docosahexaenoic|linoleic",lgc_dat$disease_trait)]<-"Metabolic biomarkers - Fatty acids"

lgc_dat$group[grep("Citrate|Glucose|Glycerol|Lactate|Pyruvate",lgc_dat$disease_trait)]<-"Metabolic biomarkers - Glycolysis related metabolites"

lgc_dat$group[grep("Alanine|Glutamine|Glycine|Histidine|Isoleucine|Leucine|Phenylalanine|Tyrosine|Valine",lgc_dat$disease_trait)]<-"Metabolic biomarkers - Amino acids"

lgc_dat$group[grep("apolipoprotein|Apolipo",lgc_dat$disease_trait)]<-"Metabolic biomarkers - Apolipoproteins"

lgc_dat$group[grep("INTERVAL",lgc_dat$phen1)]<-"Blood cell morphology"

table(lgc_dat$group,useNA="ifany")
#                      Blood cell counts                  Blood cell morphology 
#                                   6340                                   3317 
#                                 Cancer                         Cardiovascular 
#                                   1921                                   3619 
#                        Immune Mediated     Metabolic biomarkers - Amino acids 
#                                   5349                                    557 
# Metabolic biomarkers - Apolipoproteins     Metabolic biomarkers - Cholesterol 
#                                   1244                                  14855 
#     Metabolic biomarkers - Fatty acids   Metabolic biomarkers - Fluid balance 
#                                   3779                                    655 
#    Metabolic biomarkers - Inflammation          Metabolic biomarkers - Lipids 
#                                   1292                                   7732 
#                                  Other 
#                                   2652


lgc_dat$qval<-qvalue(lgc_dat$p)$qvalues
lgc_dat<-lgc_dat[which(!is.na(lgc_dat$qval)),]


lgc_dat$direction_signif<-NA
lgc_dat$direction_signif[which(lgc_dat$qval>0.05)]<-"no signif"
lgc_dat$direction_signif[which(lgc_dat$qval<=0.05 & lgc_dat$rho>0)]<-"same direction"
lgc_dat$direction_signif[which(lgc_dat$qval<=0.05 & lgc_dat$rho<0)]<-"opposite direction"

lgc_dat<-lgc_dat[which(!is.na(lgc_dat$qval)),]
lgc_dat$direction_signif<-factor(lgc_dat$direction_signif,levels=c("opposite direction","same direction","no signif"))


head(lgc_dat)

# write.table(lgc_dat,
# "~/git/IIBDGC_GWAS/plots/paper_tables/Source_data_local_genetic_correlation.tsv",col.names=T,row.names=F,quote=F,sep="\t")



# Plot count bars per significant regions - split by direction of effect

p1<-ggplot(lgc_dat[which(lgc_dat$direction_signif %in% c("same direction","opposite direction")),], 
aes(y=disease_trait, fill=direction_signif)) +
   geom_bar(stat="count") + 
    facet_grid( group ~ toupper(phen2), scales = "free",space = "free") + xlim(0,40) +
   xlab("Number of IBD regions with significant genetic correlation") + scale_fill_manual(values=c("#B2182B","#2166AC")) + 
   theme(legend.title = element_blank(),axis.title.y=element_blank(),legend.position = "bottom",plot.margin = unit(c(1,1,1,1),"cm"))
# p1


# plot mean r2 per trait:


for (ph in c("ibd","cd","uc")) {

    tmp1<-lgc_dat[which(lgc_dat$direction_signif %in% c("same direction","opposite direction") & lgc_dat$phen2==ph),]
    print(dim(tmp1))
    tmp2<-tmp1[,c("disease_trait","phen2","group")]
    tmp2<-tmp2[!duplicated(tmp2),]

    tmp2$mean_rho<-NA
    tmp2$n_loci<-NA

    for (i in 1:nrow(tmp2)) {
        tmp2$mean_rho[i]<-mean(tmp1$rho[which(tmp1$disease_trait==tmp2$disease_trait[i])])
        tmp2$n_loci[i]<-length(tmp1$rho[which(tmp1$disease_trait==tmp2$disease_trait[i])])
            
    }
    
    if(ph=="ibd") {
        dat<-tmp2
    } else {
        dat<-rbind(dat,tmp2)
    }

    rm(tmp2,tmp1)

}

p2<-ggplot(dat, aes(y = disease_trait,x=phen2, fill = mean_rho)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#B2182B",
                       mid = "white",
                       high = "#2166AC") + 
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
    legend.title = element_blank(),axis.title.y=element_blank(),legend.position = "bottom",plot.margin = unit(c(1,1,1,2),"cm")) + 
    facet_grid(group ~ toupper(phen2), scales = "free",space = "free") + 
    geom_text(aes(label = n_loci), size = 8 / .pt)
# p2


table(dat$group)
#                      Blood cell counts                  Blood cell morphology 
#                                     24                                      6 
#                                 Cancer                         Cardiovascular 
#                                      9                                     12 
#                        Immune Mediated     Metabolic biomarkers - Amino acids 
#                                     30                                      3 
# Metabolic biomarkers - Apolipoproteins     Metabolic biomarkers - Cholesterol 
#                                      5                                     62 
#     Metabolic biomarkers - Fatty acids   Metabolic biomarkers - Fluid balance 
#                                     14                                      3 
#    Metabolic biomarkers - Inflammation          Metabolic biomarkers - Lipids 
#                                      5                                     29 
#                                  Other 
#                                      9

head(dat[order(dat$mean_rho,decreasing=T),],20)
#                                            disease_trait  phen2
#                                                   <char> <char>
#  1:                                               MicroR    ibd
#  2:                                               MicroR     cd
#  3:                                 Phenylalanine levels     uc
#  4:                                    Creatinine levels     uc
#  5:                                               MicroR     uc
#  6:                                    Creatinine levels    ibd
#  7:                                 Phenylalanine levels     cd
#  8:                       Primary sclerosing cholangitis     uc
#  9:                       Primary sclerosing cholangitis    ibd
# 10:                                 Phenylalanine levels    ibd
# 11:                                 Rheumatoid arthritis     uc
# 12:                                    Atopic dermatitis     uc
# 13: Total cholesterol to total lipids ratio in large HDL     uc
# 14:  Free cholesterol to total lipids ratio in large HDL     uc
# 15:                               Ankylosing spondylitis     cd
# 16:                               Ankylosing spondylitis     uc
# 17:                       Primary sclerosing cholangitis     cd
# 18:                          Primary biliary cholangitis    ibd
# 19:                               Ankylosing spondylitis    ibd
# 20:                                            Psoriasis     cd
#                                            disease_trait  phen2
#                                    group  mean_rho n_loci
#                                   <char>     <num>  <int>
#  1:                Blood cell morphology 1.0000000      1
#  2:                Blood cell morphology 1.0000000      1
#  3:   Metabolic biomarkers - Amino acids 1.0000000      2
#  4: Metabolic biomarkers - Fluid balance 0.9862365      2
#  5:                Blood cell morphology 0.9648660      1
#  6: Metabolic biomarkers - Fluid balance 0.9221530      3 = HIGHLIGHT
#  7:   Metabolic biomarkers - Amino acids 0.9156313      3 = HIGHLIGHT
#  8:                      Immune Mediated 0.9036405     11
#  9:                      Immune Mediated 0.8880160      7
# 10:   Metabolic biomarkers - Amino acids 0.8696017      3
# 11:                      Immune Mediated 0.8539958     10
# 12:                      Immune Mediated 0.8538025      2
# 13:   Metabolic biomarkers - Cholesterol 0.8533070      1
# 14:   Metabolic biomarkers - Cholesterol 0.8455910      1
# 15:                      Immune Mediated 0.8385815     17
# 16:                      Immune Mediated 0.8235301     11
# 17:                      Immune Mediated 0.8010427      3
# 18:                      Immune Mediated 0.8007168     15
# 19:                      Immune Mediated 0.7854479     18
# 20:                      Immune Mediated 0.7800877     26
#                                    group  mean_rho n_loci


head(dat[order(dat$mean_rho,decreasing=F),],20)
#                                               disease_trait  phen2
#                                                      <char> <char>
#  1:                                             Lung cancer    ibd
#  2:                                             Lung cancer     cd
#  3:         Atrial fibrillation and flutter (PheCode 427.2)     cd
#  4:                                             Lung cancer     uc
#  5:                                                 RBC-FSC     uc
#  6:                                                  RBC-He     uc
#  7:                Coronary atherosclerosis (PheCode 411.4)     cd
#  8:               Mean corpuscular hemoglobin concentration    ibd
#  9:                    Ischemic heart disease (PheCode 411)     uc
# 10:        Triglycerides to total lipids ratio in large HDL     cd
# 11:                              Triglycerides in large HDL     cd
# 12:        Triglycerides to total lipids ratio in large HDL    ibd
# 13:    Free cholesterol to total lipids ratio in medium HDL     cd
# 14:                                                 RBC-SSC     cd
# 15: Ratio of 22:6 docosahexaenoic acid to total fatty acids     cd
# 16:       Ratio of omega-3 fatty acids to total fatty acids     cd
# 17:                Coronary atherosclerosis (PheCode 411.4)     uc
# 18:                                       Colorectal cancer     cd
# 19:                    Ischemic heart disease (PheCode 411)     cd
# 20:                    Ischemic heart disease (PheCode 411)    ibd
#                                               disease_trait  phen2
#                                  group   mean_rho n_loci
#                                 <char>      <num>  <int>
#  1:                             Cancer -1.0000000      1
#  2:                             Cancer -1.0000000      1
#  3:                     Cardiovascular -1.0000000      1
#  4:                             Cancer -1.0000000      1
#  5:              Blood cell morphology -0.9527005      2
#  6:              Blood cell morphology -0.9360520      2
#  7:                     Cardiovascular -0.9290260      2
#  8:                  Blood cell counts -0.9212525      4 = HIGHLIGHT
#  9:                     Cardiovascular -0.9104730      1
# 10:      Metabolic biomarkers - Lipids -0.9025770      3 = HIGHLIGHT
# 11:      Metabolic biomarkers - Lipids -0.8869435      2
# 12:      Metabolic biomarkers - Lipids -0.8820160      2
# 13: Metabolic biomarkers - Cholesterol -0.8357690      1
# 14:              Blood cell morphology -0.8205395      2
# 15: Metabolic biomarkers - Fatty acids -0.8107080      1
# 16: Metabolic biomarkers - Fatty acids -0.8002260      1
# 17:                     Cardiovascular -0.7903880      2
# 18:                             Cancer -0.7879932      5
# 19:                     Cardiovascular -0.7620300      2
# 20:                     Cardiovascular -0.7529905      2
#                                  group   mean_rho n_loci



### which is the region/s with larger number of traits linked?

signif<-lgc_dat[which(lgc_dat$qval<0.05),c("region","disease_trait","rho","group","phen2")]

# number of pleiotropic IBD regions
length(names(table(signif$region)))
# [1] 146

ibd<-signif[which(signif$phen2=="ibd"),]

p<-ggplot(ibd, aes(region, disease_trait, fill= rho)) + 
  geom_tile() +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_uc<-ggplot(cd, aes(region, disease_trait, fill= rho)) + 
  geom_tile() +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") + theme(axis.text.x = element_text(angle = 45, hjust = 1)





p<-ggarrange(p1,p2,ncol=2,widths=c(4,2),align=c("h"),labels=c("a","b"))
# p

ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_figure_local_genetic_correlation.png",
  p,
  width = 220,
  height = ((nrow(dat)/3)*3)+5,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)


