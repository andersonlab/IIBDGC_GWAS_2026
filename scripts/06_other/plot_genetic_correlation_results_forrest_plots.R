# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# singularity exec iibdgc_postprocess_10_singularity.sif
path_gwas="/path/to/ibdgwas/IIBDGC/"

MEM=2000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(rcartocolor)
library(qvalue)

rm(list=ls())

path_gwas="/path/to/ibdgwas/IIBDGC/"

cols<-c("#db7107","#004488","#BB5566","darkgrey")

files<-list.files(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/genetic_correlation/",sep=""))
files<-files[grep(".log",files)]
files<-files[which(files!="ldsc.log")]
# files<-allfiles[!grepl("test",allfiles)]



# # for testing purposes
# files_hm3<-files[grep("hm3",files)]
# files<-files[!grepl("hm3",files)]

for (i in 1:length(files)) {
# for (i in c(1,262:264,525:527,788,789)) {

    file_tmp<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/genetic_correlation/",files[i],sep="")

    if (file.info(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/genetic_correlation/",files[i],sep=""))$size>0) {
        tmp<-fread(file_tmp,fill=T,skip=60)
        colnames(tmp)<-c("p1","p2","rg","se","z","p","h2_obs","h2_obs_se","h2_int","h2_int_se" ,"gcov_int","gcov_int_se")
        tmp<-tmp[grep("/lustre/",tmp$p1),]

        tmp$file<-file_tmp
    
        if(i==1) {
            all<-tmp
        } else {
            all<-rbind(all,tmp)
        }
        rm(tmp)
    }

    rm(file_tmp)
}

all<-as.data.frame(all)

all$p1<-gsub("/path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/","",all$p1)
all$p1<-gsub("_eur_tier2_list_variants_rate_0.5_info_0.9_only_SNPs_nohla_sumstats_munged.sumstats","",all$p1)
all$p1<-gsub("_only_SNPs_sumstats_munged.sumstats","",all$p1)
all$p1<-gsub("_sumstats_munged.sumstats","",all$p1)

all$p2<-gsub("/path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/","",all$p2)
all$p2<-gsub("_eur_tier2_list_variants_rate_0.5_info_0.9_only_SNPs_nohla_sumstats_munged.sumstats","",all$p2)
all$p2<-gsub("_only_SNPs_sumstats_munged.sumstats","",all$p2)
all$p2<-gsub("_sumstats_munged.sumstats","",all$p2)


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

all<-merge(all,dat,by.x="p2",by.y="study_accession",all.x=T)

# no data from mapp:

table(all$disease_trait,all$group)

threshold<-0.05/(nrow(all))
threshold

all$rg<-as.numeric(all$rg)


# CREATE PLOTS PER GROUPS - OTHERWISE TOO LARGE TO BE HANDLED:

all$group[grep("levels",all$disease_trait)]<-"Metabolic biomarkers - Other"

all$group[grep("levels",all$disease_trait)]<-"Metabolic biomarkers - Fluid balance"

all$group[grep("Glycoprotein acetylation|C-reactive",all$disease_trait)]<-"Metabolic biomarkers - Inflammation"

all$group[grep("cholesterol|VLDL particles|IDL particles|LDL particles|HDL particles",all$disease_trait)]<-"Metabolic biomarkers - Cholesterol"
all$group[grep("Cholesterol",all$disease_trait)]<-"Metabolic biomarkers - Cholesterol"
all$group[grep("Cholesteryl",all$disease_trait)]<-"Metabolic biomarkers - Cholesterol"
all$group[grep("Triglyceride",all$disease_trait)]<-"Metabolic biomarkers - Lipids"
all$group[grep("Total lipid|diacylglycerol",all$disease_trait)]<-"Metabolic biomarkers - Lipids"
all$group[grep("Phospholipid",all$disease_trait)]<-"Metabolic biomarkers - Lipids"
all$group[grep("phosphoglycerides",all$disease_trait)]<-"Metabolic biomarkers - Lipids"

all$group[grep("fatty|unsaturation|docosahexaenoic|linoleic",all$disease_trait)]<-"Metabolic biomarkers - Fatty acids"

all$group[grep("Citrate|Glucose|Glycerol|Lactate|Pyruvate",all$disease_trait)]<-"Metabolic biomarkers - Glycolysis related metabolites"

all$group[grep("Alanine|Glutamine|Glycine|Histidine|Isoleucine|Leucine|Phenylalanine|Tyrosine|Valine",all$disease_trait)]<-"Metabolic biomarkers - Amino acids"

all$group[grep("apolipoprotein|Apolipo",all$disease_trait)]<-"Metabolic biomarkers - Apolipoproteins"

all$group[grep("INTERVAL",all$p2)]<-"Blood cell morphology"

all$disease_trait<-as.factor(all$disease_trait)
all$disease_trait<-factor(all$disease_trait,levels=rev(levels(all$disease_trait)))


table(all$group,useNA="ifany")
#                                     Blood cell counts 
#                                                    45 
#                                 Blood cell morphology 
#                                                   189 
#                                                Cancer 
#                                                    18 
#                                        Cardiovascular 
#                                                    18 
#                                       Immune Mediated 
#                                                    36 
#                    Metabolic biomarkers - Amino acids 
#                                                    27 
#                Metabolic biomarkers - Apolipoproteins 
#                                                     9 
#                    Metabolic biomarkers - Cholesterol 
#                                                   330 
#                    Metabolic biomarkers - Fatty acids 
#                                                    51 
#                  Metabolic biomarkers - Fluid balance 
#                                                    36 
# Metabolic biomarkers - Glycolysis related metabolites 
#                                                    15 
#                   Metabolic biomarkers - Inflammation 
#                                                     6 
#                         Metabolic biomarkers - Lipids 
#                                                   228 
#                                                 Other 
#                                                    15 
#                                                  <NA> 
#                                                     9 

table(all$disease_trait,all$group,useNA="ifany")

all$rg<-as.numeric(all$rg)
all$se<-as.numeric(all$se)

all$ci_lower<-all$rg-(1.96*all$se)
all$ci_upper<-all$rg+(1.96*all$se)


#################################

all<-all[!grepl("hm3",all$p1),]
all<-all[!grepl("_metabolic_biomarkers_test_gc",all$file),]

all[which(all$p2 %in% c("cd","uc","ibd")),]
#       p2  p1     rg         se          z p h2_obs h2_obs_se h2_int h2_int_se
# 13    cd  cd 1.0002 4.3087e-05 23212.8936 0 0.3516    0.0292 1.0284    0.0056
# 14    cd  uc 0.6798 1.6900e-02    40.1284 0 0.3576    0.0281 1.0227    0.0042
# 15    cd ibd 0.9279 6.0000e-03    154.997 0 0.3541    0.0289 1.0263    0.0046
# 883  ibd ibd 1.0003 3.8741e-05 25820.8382 0 0.2308    0.0151 1.0484    0.0041
# 884  ibd  uc 0.9256 5.8000e-03   160.6781 0 0.2311    0.0151 1.0479    0.0041
# 885  ibd  cd 0.9279 6.0000e-03    154.997 0 0.2317    0.0154 1.0468    0.0038
# 1030  uc ibd 0.9256 5.8000e-03   160.6781 0  0.262    0.0165 1.0383    0.0038
# 1031  uc  cd 0.6798 1.6900e-02    40.1284 0 0.2643    0.0165 1.0352    0.0037
# 1032  uc  uc 1.0005 5.4010e-05 18523.5103 0  0.263    0.0167 1.0369    0.0037
#      gcov_int gcov_int_se
# 13     1.0283      0.0056
# 14     0.2414      0.0021
# 15     0.7042      0.0034
# 883    1.0484       0.004
# 884    0.8308      0.0035
# 885    0.7042      0.0034
# 1030   0.8308      0.0035
# 1031   0.2414      0.0021
# 1032   1.0368      0.0037

all<-all[which(!all$p2 %in% c("cd","uc","ibd")),]

# test usage of qval to determine significance:
all$p<-as.numeric(all$p)
all$pvalue<-all$p
all$qval<-qvalue(all$p)$qvalues

threshold<-0.05/(nrow(all))
threshold

all$pheno<-toupper(all$p1)

all$pheno_signif<-as.character(all$pheno)
all$pheno_signif[which(all$pvalue>threshold)]<-"NA"
all$pheno_signif<-factor(all$pheno_signif,levels=rev(c("CD","IBD","UC","NA")))

all$pheno_signif_qval<-as.character(all$pheno)
all$pheno_signif_qval[which(all$qval>=0.01)]<-"NA"
all$pheno_signif_qval<-factor(all$pheno_signif_qval,levels=rev(c("CD","IBD","UC","NA")))

all$pheno<-as.factor(all$pheno)
all$pheno<-factor(all$pheno,levels=c("CD","IBD","UC"))
table(all$pheno)

cols<-c("CD" = "#db7107","IBD" = "#004488","UC" = "#BB5566","NA" = "darkgrey")


max(all$ci_upper)
min(all$ci_lower)

all$disease_trait<-as.character(all$disease_trait)
all$disease_trait[which(all$disease_trait=="Clonal haematopoiesis of indeterminate potential (with or without mosaic chromosomal alteration)")]<-"Clonal haematopoiesis"
all$disease_trait<-as.factor(all$disease_trait)
all$disease_trait<-factor(all$disease_trait,levels=rev(levels(all$disease_trait)))



#######################################
# main figure keeping only significant:

signif<-names(table(as.character(all$disease_trait[which(all$pheno_signif_qval!="NA")])))
tmp<-all[which(all$disease_trait %in% signif),]

tmp$disease_trait<-as.character(tmp$disease_trait)
tmp$disease_trait[grep("INTERVAL",tmp$p2)]<-tmp$MAPPED_TRAIT[grep("INTERVAL",tmp$p2)]

dim(tmp)
# [1] 261  23
dim(tmp[!duplicated(tmp$disease_trait),])
# [1] 87 23


p<-ggplot(tmp, aes(x = rg, y = disease_trait, color = pheno_signif_qval)) +
    geom_point(shape = 16, size = 2) + xlim(-0.40,0.75) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
    # facet_nested( group + disease_trait ~ . , scales = "free",space = "free") +
    facet_grid(group ~ pheno, scales = "free",space = "free") +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) + 
    scale_color_manual(values=cols) + guides(fill = "none", color="none") + 
    geom_vline(xintercept = 0,linetype="dotted",colour="black")

ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_2_genetic_correlation_qval_forrest_plot.pdf",
  p,
  width = 140,
  height = ((nrow(tmp)/8)*6)+5,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)


###################################################
#### create panel a for a main figure:

panel_name<-"Figure_2a_genetic_correlation_qval_forrest_plot"

list_traits_tokeep<-tmp$disease_trait[which(tmp$group %in% c("Metabolic biomarkers - Apolipoproteins",
"Metabolic biomarkers - Amino acids","Metabolic biomarkers - Inflammation"))]

list_traits_tokeep<-c(list_traits_tokeep,"Total cholesterol in large HDL","Red blood cell count",
"Hematocrit","Mean corpuscular hemoglobin","Mean corpuscular hemoglobin concentration","Mean corpuscular volume",
"RET-He","RBC-He","MicroR","Hyper-He")

list_traits_tokeep<-list_traits_tokeep[!duplicated(list_traits_tokeep)]

tmp1<-tmp[which(tmp$disease_trait %in% c(list_traits_tokeep)),]

tmp1$rg_direction<-NA
tmp1$rg_direction[which(tmp1$rg>0)]<-"+"
tmp1$rg_direction[which(tmp1$rg<0)]<-"-"
tmp1$rg_direction[which(tmp1$pheno_signif_qval=="NA")]<-NA


write.table(tmp1,paste0("~/git/IIBDGC_GWAS/plots/paper_figures/",panel_name,".tsv"),col.names=T,row.names=F,sep="\t")


rm(tmp1)

####################################

# save results:

# save the list of singificant pairs:
list_signif_pairs<-all[which(all$pheno_signif_qval!="NA"),c("p1","p2","disease_trait","MAPPED_TRAIT","group","rg","se","qval")]
list_signif_pairs<-list_signif_pairs[!duplicated(list_signif_pairs),]

list_signif_pairs$disease_trait<-as.character(list_signif_pairs$disease_trait)

write.table(list_signif_pairs,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/genetic_correlation/list_pairs_traits_significantly_correlated_with_ibd_cd_uc_with_qval",sep="")
,col.names=T,row.names=F,quote=F,sep="\t")

############

all<-all[,c("p2","disease_trait","p1","rg","ci_lower","ci_upper","se","z","p","qval","h2_obs","h2_obs_se",
"h2_int","h2_int_se","gcov_int","gcov_int_se","group")]

# relabel to phenotypes the study accessions
dat<-fread(paste(path_gwas,"resources/gwas_summary_statistics/gwas-catalog-v1.0.3.1-studies-r2025-05-13.tsv.gz",sep=""),quote="")
dat$study_accession<-dat$'STUDY ACCESSION'
dat$disease_trait<-dat$'DISEASE/TRAIT'

dat<-dat[which(dat$study_accession %in% all$p2)]

dat$map_trait<-gsub("http://www.ebi.ac.uk/efo/","",dat$MAPPED_TRAIT_URI)

dat$path<-NA

# two formats of IDs
new<-c("GCST90132223","GCST011096","GCST90472771","GCST90480502","GCST90435254","GCST90476007","GCST90475929","GCST90475936","GCST90475967","GCST90475214","GCST90013445","GCST90274714")
rm(list_tmp)
for (i in 1:nrow(dat)) {

    if (!dat$'STUDY ACCESSION'[i] %in% new) {
        dat$path[i]<-paste("wget ",dat$'SUMMARY STATS LOCATION'[i],"/harmonised/",dat$'PUBMED ID'[i],"-",dat$'STUDY ACCESSION'[i],"-",dat$map_trait[i],".h.tsv.gz -O ",path_gwas,"/resources/gwas_summary_statistics/",dat$'STUDY ACCESSION'[i],".h.tsv.gz",sep="")
    } else {
        dat$path[i]<-paste("wget ",dat$'SUMMARY STATS LOCATION'[i],"/harmonised/",dat$'STUDY ACCESSION'[i],".h.tsv.gz -O ",path_gwas,"/resources/gwas_summary_statistics/",dat$'STUDY ACCESSION'[i],".h.tsv.gz",sep="")
    }
        
}


# manually add:
dat$path[which(dat$'STUDY ACCESSION'=="GCST90129505")]<-"http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90129001-GCST90130000/GCST90129505/GCST90129505_buildGRCh37.tsv"
dat$path[which(dat$'STUDY ACCESSION'=="GCST90029070")]<-"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90029001-GCST90030000/GCST90029070/GCST90029070_buildGRCh37.tsv.gz"


all<-merge(all,dat[,c("STUDY ACCESSION","path")],by.x="p2",by.y="STUDY ACCESSION",all.x=T)

all$path[which(is.na(all$path))]<-paste0("ftp://ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/sysmex_blood_cell_genetics/",all$p2[grep("INTERVAL",all$p2)],"_gwas_normalised_chr[1..22]_imputed.assoc.gz")

write.table(all,
"~/git/IIBDGC_GWAS/plots/paper_tables/Source_data_genetic_correlation.tsv",col.names=T,row.names=F,quote=F,sep="\t")

# save the list of IDs for the significant results:
list_signif<-all[which(all$disease_trait %in% signif),"p2"]
list_signif<-list_signif[!duplicated(list_signif)]

write.table(list_signif,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/genetic_correlation/list_traits_significantly_correlated_with_ibd_cd_uc",sep="")
,col.names=F,row.names=F,quote=F,sep="\t")

q("no")






##########################################################
#### test with cholesterol pair by pair comparisons:

# cho<-all[grep("_metabolic_biomarkers_test_gc",all$file),]
# cho$group<-"Cholesterol"

# # get a random list of IDs 20:
# ids<-names(table(as.character(cho$p2)))[sample(20)]
# cho$p<-as.numeric(cho$p)

# cho<-cho[which(cho$p1 %in% ids & cho$p2 %in% ids),]
# cho$text<-paste("rg=",round(cho$rg,digits=2),"\n","p=E-",round(-log10(cho$p),digits=0),sep="")

# p<-ggplot(cho, aes(x = toupper(p1), y = p2, fill = rg)) +
#   geom_tile(color = "white",lwd = 1.5,linetype = 1) +
#   scale_fill_gradient2(low = "#075AFF",
#                        mid = "white",
#                        high = "#FF0000",midpoint = 0,limits= c(-1.1,1.1)) +
# geom_text(aes(label = text), color = "black", size = 4) + theme(axis.title.x=element_blank(),axis.title.y=element_blank())

                       
# ggsave(
#   "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_cholesterol_traits_genetic_correlation.pdf",
#   p,
#   width = 180,
#   height = 180,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=3
# )


#### IBD - test 

# ibd<-all[which(all$p1 %in% c("ibd","cd","uc") & all$p2 %in% c("ibd","cd","uc","GCST004030")),]
# ibd$group<-"IBD"

# ibd2<-all[grep("hm3",all$p1),]
# ibd2$N_variants<-"944844"
# ibd2$ldscore_reference<-"NA"
# ibd2$ldscore_reference[grep("chrposid_b38",ibd2$p2)]<-"gsa_weights_and_ld_reference"
# ibd2$ldscore_reference[grep("rsid_b38",ibd2$p2)]<-"1000GP_weights_weights_estimated_with_only_hm3_and_ld_reference"
# ibd2$ldscore_reference[grep("rsid_b37",ibd2$p2)]<-"UKB_weights_and_ld_reference"
# ibd2$ldscore_reference[grep("gsa_weights_estimated_with_only_hm3_SNPs",ibd2$file)]<-"GSA_weights_estimated_with_only_hm3_and_ld_reference"


# # add original to the table:
# tmp<-ibd[which(ibd$p2=="GCST004030"),]
# tmp$N_variants<-"12123320"
# tmp$ldscore_reference<-"gsa_weights_and_ld_reference"

# ibd2<-rbind(ibd2,tmp)
# rm(tmp)
# ibd2$p2<-"GCST004030"
# ibd2$p1<-gsub("_eur_.*","",ibd2$p1)

# ibd2<-ibd2[,c("p1","p2","rg","se","z","p","h2_obs","h2_obs_se","h2_int","h2_int_se","gcov_int","gcov_int_se","p","N_variants","ldscore_reference")]
# write.table(ibd2,"~/git/IIBDGC_GWAS/plots/paper_tables/Supplementary_table_genetic_correlation_under_different_ld_references.tsv",col.names=T,row.names=F,quote=F,sep="\t")


# find duplicated traits:

# tmp<-all[which(all$pheno=="CD" & all$p2!="cd"),]

# dups<-tmp$MAPPED_TRAIT[duplicated(tmp$MAPPED_TRAIT)]
# dups<-dups[!duplicated(dups)]

# for (i in 1:length(dups)) {
#   print(dups[i])

#   print(tmp$p2[which(tmp$MAPPED_TRAIT==dups[i])])
#   # tmp[which(tmp$MAPPED_TRAIT==dups[i]),]
# }



#################################
# non metabolic biomarkers:

# tmp<-all[which(all$group %in% c("Blood cell counts","Cancer","Cardiovascular","Immune Mediated","Other")),]

# p1<-ggplot(tmp, aes(x = rg, y = disease_trait, color = pheno_signif)) +
#     geom_point(shape = 16, size = 2) + xlim(-0.40,0.75) +
#     geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
#     # facet_nested( group + disease_trait ~ . , scales = "free",space = "free") +
#     facet_grid(group ~ pheno, scales = "free",space = "free") +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) + 
#     scale_color_manual(values=cols) + guides(fill = "none", color="none") + 
#     geom_vline(xintercept = 0,linetype="dotted",colour="black")

# ggsave(
#   "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_9_non_met_genetic_correlation_bonferroni_forrest_plot.pdf",
#   p1,
#   width = 180,
#   height = ((nrow(tmp)/3)*4)+5,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=2
# )


# p2<-ggplot(tmp, aes(x = rg, y = disease_trait, color = pheno_signif_qval)) +
#     geom_point(shape = 16, size = 2) + xlim(-0.4,0.75) +
#     geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
#     # facet_nested( group + disease_trait ~ . , scales = "free",space = "free") +
#     facet_grid(group ~ pheno, scales = "free",space = "free") +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) + 
#     scale_color_manual(values=cols) + guides(fill = "none", color="none") + 
#     geom_vline(xintercept = 0,linetype="dotted",colour="black")

# ggsave(
#   "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_9_non_met_genetic_correlation_qval_forrest_plot.pdf",
#   p2,
#   width = 180,
#   height = ((nrow(tmp)/3)*4)+5,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=2
# )

#################################
# Blood morphology:

# tmp<-all[which(all$group %in% c("Blood cell morphology")),]

# tmp<-merge(tmp,dat2[,c("cell","study_accession")],by.x="p2",by.y="study_accession",all.x=T)

# p1<-ggplot(tmp, aes(x = rg, y = disease_trait, color = pheno_signif)) +
#     geom_point(shape = 16, size = 2) + xlim(-0.40,0.75) +
#     geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
#     # facet_nested( group + disease_trait ~ . , scales = "free",space = "free") +
#     facet_grid(group ~ pheno, scales = "free",space = "free") +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) + 
#     scale_color_manual(values=cols) + guides(fill = "none", color="none") + 
#     geom_vline(xintercept = 0,linetype="dotted",colour="black")

# ggsave(
#   "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_9_non_met_genetic_correlation_bonferroni_forrest_plot.pdf",
#   p1,
#   width = 180,
#   height = ((nrow(tmp)/3)*4)+5,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=2
# )

# p2<-ggplot(tmp, aes(x = rg, y = MAPPED_TRAIT, color = pheno_signif_qval)) +
#     geom_point(shape = 16, size = 2) + xlim(-0.40,0.75) +
#     geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
#     # facet_nested( group + disease_trait ~ . , scales = "free",space = "free") +
#     facet_grid(cell ~ pheno, scales = "free",space = "free") +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) + 
#     scale_color_manual(values=cols) + guides(fill = "none", color="none") + 
#     geom_vline(xintercept = 0,linetype="dotted",colour="black")

# ggsave(
#   "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_9_cell_morphology_genetic_correlation_qval_forrest_plot.pdf",
#   p2,
#   width = 180,
#   height = ((nrow(tmp)/3)*4)+5,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=2
# )



#################################
### Cholesterol

# tmp<-all[which(all$group %in% c("Metabolic biomarkers - Cholesterol")),]

# # subgroup:
# tmp$group<-"Cholesterol"
# tmp$group[grep("LDL",tmp$disease_trait)]<-"Cholesterol LDL"
# tmp$group[grep("HDL",tmp$disease_trait)]<-"Cholesterol HDL"
# tmp$group[grep("VLDL",tmp$disease_trait)]<-"Cholesterol VLDL"
# tmp$group[grep("IDL",tmp$disease_trait)]<-"Cholesterol IDL"
# tmp$group[grep("Cholesteryl",tmp$disease_trait)]<-"Cholesteryl"

# p1<-ggplot(tmp[which(tmp$group %in% c("Cholesterol","Cholesterol HDL","Cholesterol IDL","Cholesterol LDL")),], aes(x = rg, y = disease_trait, color = pheno_signif)) +
#     geom_point(shape = 16, size = 2) + xlim(-0.25,0.25) +
#     geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
#     facet_grid(group ~ pheno, scales = "free",space = "free") +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) + 
#     scale_color_manual(values=cols) + guides(fill = "none", color="none") + 
#     geom_vline(xintercept = 0,linetype="dotted",colour="black")

# p2<-ggplot(tmp[which(tmp$group %in% c("Cholesterol VLDL","Cholesteryl")),], aes(x = rg, y = disease_trait, color = pheno_signif)) +
#     geom_point(shape = 16, size = 2) + xlim(-0.25,0.25) +
#     geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
#     facet_grid(group ~ pheno, scales = "free",space = "free") +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) + 
#     scale_color_manual(values=cols) + guides(fill = "none", color="none") + 
#     geom_vline(xintercept = 0,linetype="dotted",colour="black")


# p<-ggarrange(p1,p2,ncol=2,common.legend =T,legend="right")

# ggsave(
#   "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_10_cholest_genetic_correlation_bonferroni_forrest_plot.pdf",
#   p,
#   width = 360,
#   height = ((nrow(tmp)/6)*4)+5,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=2
# )


# p1<-ggplot(tmp[which(tmp$group %in% c("Cholesterol","Cholesterol HDL","Cholesterol IDL","Cholesterol LDL")),], aes(x = rg, y = disease_trait, color = pheno_signif_qval)) +
#     geom_point(shape = 16, size = 2) + xlim(-0.25,0.25) +
#     geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
#     facet_grid(group ~ pheno, scales = "free",space = "free") +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) + 
#     scale_color_manual(values=cols) + guides(fill = "none", color="none") + 
#     geom_vline(xintercept = 0,linetype="dotted",colour="black")

# p2<-ggplot(tmp[which(tmp$group %in% c("Cholesterol VLDL","Cholesteryl")),], aes(x = rg, y = disease_trait, color = pheno_signif_qval)) +
#     geom_point(shape = 16, size = 2) + xlim(-0.25,0.25) +
#     geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
#     facet_grid(group ~ pheno, scales = "free",space = "free") +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) + 
#     scale_color_manual(values=cols) + guides(fill = "none", color="none") + 
#     geom_vline(xintercept = 0,linetype="dotted",colour="black")


# p<-ggarrange(p1,p2,ncol=2,common.legend =T,legend="right")

# ggsave(
#   "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_10_cholest_genetic_correlation_qval_forrest_plot.pdf",
#   p,
#   width = 360,
#   height = ((nrow(tmp)/6)*4)+5,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=2
# )



#################################
### Lipids

# tmp<-all[which(all$group %in% c("Metabolic biomarkers - Lipids")),]

# # subgroup:
# tmp$group<-"Lipids"
# tmp$group[grep("Triglyceride",tmp$disease_trait)]<-"Lipids - Triglyceride"
# tmp$group[grep("Phospholipid",tmp$disease_trait)]<-"Lipids - Phospholipid"

# min(tmp$ci_lower)
# max(tmp$ci_upper)

# p1<-ggplot(tmp[which(tmp$group %in% c("Lipids")),], aes(x = rg, y = disease_trait, color = pheno_signif)) +
#     geom_point(shape = 16, size = 2) + xlim(-0.3,0.25) +
#     geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
#     facet_grid(group ~ pheno, scales = "free",space = "free") +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) + 
#     scale_color_manual(values=cols) + guides(fill = "none", color="none") + 
#     geom_vline(xintercept = 0,linetype="dotted",colour="black")

# p2<-ggplot(tmp[which(tmp$group %in% c("Lipids - Triglyceride")),], aes(x = rg, y = disease_trait, color = pheno_signif)) +
#     geom_point(shape = 16, size = 2) + xlim(-0.3,0.25) +
#     geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
#     facet_grid(group ~ pheno, scales = "free",space = "free") +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) + 
#     scale_color_manual(values=cols) + guides(fill = "none", color="none") + 
#     geom_vline(xintercept = 0,linetype="dotted",colour="black")


# p3<-ggplot(tmp[which(tmp$group %in% c("Lipids - Phospholipid")),], aes(x = rg, y = disease_trait, color = pheno_signif)) +
#     geom_point(shape = 16, size = 2) + xlim(-0.3,0.25) +
#     geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
#     facet_grid(group ~ pheno, scales = "free",space = "free") +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) + 
#     scale_color_manual(values=cols) + guides(fill = "none", color="none") + 
#     geom_vline(xintercept = 0,linetype="dotted",colour="black")


# p<-ggarrange(p1,p2,p3,ncol=2,nrow=2,common.legend =T,legend="right",align = "v")


# ggsave(
#   "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_11_Lipids_genetic_correlation_bonferroni_forrest_plot.pdf",
#   p,
#   width = 360,
#   height = ((nrow(tmp)/6)*4)+5,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=2
# )




# p1<-ggplot(tmp[which(tmp$group %in% c("Lipids")),], aes(x = rg, y = disease_trait, color = pheno_signif_qval)) +
#     geom_point(shape = 16, size = 2) + xlim(-0.3,0.25) +
#     geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
#     facet_grid(group ~ pheno, scales = "free",space = "free") +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) + 
#     scale_color_manual(values=cols) + guides(fill = "none", color="none") + 
#     geom_vline(xintercept = 0,linetype="dotted",colour="black")

# p2<-ggplot(tmp[which(tmp$group %in% c("Lipids - Triglyceride")),], aes(x = rg, y = disease_trait, color = pheno_signif_qval)) +
#     geom_point(shape = 16, size = 2) + xlim(-0.3,0.25) +
#     geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
#     facet_grid(group ~ pheno, scales = "free",space = "free") +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) + 
#     scale_color_manual(values=cols) + guides(fill = "none", color="none") + 
#     geom_vline(xintercept = 0,linetype="dotted",colour="black")


# p3<-ggplot(tmp[which(tmp$group %in% c("Lipids - Phospholipid")),], aes(x = rg, y = disease_trait, color = pheno_signif_qval)) +
#     geom_point(shape = 16, size = 2) + xlim(-0.3,0.25) +
#     geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
#     facet_grid(group ~ pheno, scales = "free",space = "free") +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) + 
#     scale_color_manual(values=cols) + guides(fill = "none", color="none") + 
#     geom_vline(xintercept = 0,linetype="dotted",colour="black")


# p<-ggarrange(p1,p2,p3,ncol=2,nrow=2,common.legend =T,legend="right",align = "v")


# ggsave(
#   "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_11_Lipids_genetic_correlation_qval_forrest_plot.pdf",
#   p,
#   width = 360,
#   height = ((nrow(tmp)/6)*4)+5,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=2
# )


#################################
### Other levels

# tmp<-all[which(all$group %in% c("Metabolic biomarkers - Other levels")),]


# min(tmp$ci_lower)
# max(tmp$ci_upper)


# p<-ggplot(tmp, aes(x = rg, y = disease_trait, color = pheno_signif)) +
#     geom_point(shape = 16, size = 2) + xlim(-0.40,0.75) +
#     geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
#     # facet_nested( group + disease_trait ~ . , scales = "free",space = "free") +
#     facet_grid(group ~ pheno, scales = "free",space = "free") +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) + 
#     scale_color_manual(values=cols) + guides(fill = "none", color="none") + 
#     geom_vline(xintercept = 0,linetype="dotted",colour="black")


# ggsave(
#   "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_12_other_metab_genetic_correlation_bonferroni_forrest_plot.pdf",
#   p1,
#   width = 180,
#   height = ((nrow(tmp)/3)*4)+5,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=2
# )

# p<-ggplot(tmp, aes(x = rg, y = disease_trait, color = pheno_signif_qval)) +
#     geom_point(shape = 16, size = 2) + xlim(-0.25,0.25) +
#     geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25) +
#     # facet_nested( group + disease_trait ~ . , scales = "free",space = "free") +
#     facet_grid(group ~ pheno, scales = "free",space = "free") +
#     theme(axis.title.x=element_blank(),axis.title.y=element_blank(),strip.text.y = element_text(angle = 360,margin = margin(0,0,0,0))) + 
#     scale_color_manual(values=cols) + guides(fill = "none", color="none") + 
#     geom_vline(xintercept = 0,linetype="dotted",colour="black")


# ggsave(
#   "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_12_other_metab_genetic_correlation_qval_forrest_plot.pdf",
#   p,
#   width = 180,
#   height = ((nrow(tmp)/3)*4)+5,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=2
# )


