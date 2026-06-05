# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM=6000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q cpu-interactive R

library(data.table)
library(rtracklayer)
library(stringr)

rm(list=ls())

pheno<-c("ibd","cd","uc")

print(pheno)

path_gwas  <- "/path/to/ibdgwas/IIBDGC/"
path_eqtl  <- "/path/to/project"
R2_THRESHOLD <- 0.6

# Get gene names
gtf<-rtracklayer::import(paste(path_gwas,'post_imputation/2022/analysis/metaanalysis/annotation/gencode.v39.annotation.gtf.gz',sep=""))
gene<-as.data.frame(gtf)
gene<-gene[,c("gene_id","gene_name")]
gene$gene<-gsub("\\.[0-9]*","",gene$gene_id)
gene<-gene[!duplicated(gene$gene),]
rm(gtf)

# map to rename conditions:
map<-fread(paste0(path_eqtl,"map_files.txt.gz"),head=T)
map<-as.data.frame(map)

# GET LD DATA

ld_file <- paste0(path_gwas, "post_imputation/2022/analysis/conditional_analysis/eur/chr%d_list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_for_ld.ld")
ld <- rbindlist(lapply(1:22, function(chr) fread(sprintf(ld_file, chr))))
ld <- as.data.frame(ld)
dim(ld)
# [1] 1330481       7

# merge coloc data by ld
ld$ID1 <- paste0(ld$SNP_A, "_", ld$SNP_B)
ld$ID2 <- paste0(ld$SNP_B, "_", ld$SNP_A)

for (i in 1:length(pheno)) {

    print(pheno[i])


    tmp1<-fread(paste("/path/to/project",toupper(pheno[i]),".gz",sep=""),head=F)
    tmp2<-fread(paste("/path/to/project",toupper(pheno[i]),".gz",sep=""),head=F)
    # do not include edQTL data
    # tmp3<-fread(paste("/path/to/project",toupper(pheno[i]),".gz",sep=""),head=F)
    tmp4<-fread(paste("/path/to/project",toupper(pheno[i]),".gz",sep=""),head=F)
    tmp5<-fread(paste("/path/to/project",(pheno[i]),".gz",sep=""),head=F,skip=1)
    tmp6<-fread(paste("/path/to/project",toupper(pheno[i]),".gz",sep=""),head=F)
    tmp7<-fread(paste("/path/to/project",toupper(pheno[i]),".gz",sep=""),head=F)

    tmp7$V3<-paste0("chr",tmp7$V3)
    tmp7$V13<-paste0("chr",tmp7$V13)
    tmp7$V14<-paste0("chr",tmp7$V14)

    tmp<-rbind(tmp1,tmp2,tmp4,tmp5,tmp6,tmp7)
    rm(tmp1,tmp2,tmp4,tmp5,tmp6,tmp7)

    tmp<-as.data.frame(tmp)

    colnames(tmp)<-c("condition_name","phenotype_id","snp_id",".row","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","qtl_pval","gwas_pval","qtl_lead","gwas_lead","chr","gwas_lead_pos","qtl_lead_pos","gwas_trait")
    tmp$condition_name<-gsub("-","_",tmp$condition_name)
    tmp<-merge(tmp,map,by.x="condition_name",by.y="id_map")
    
    print(nrow(tmp))
    tmp$PP.H4.abf<-as.numeric(tmp$PP.H4.abf)
    tmp<-tmp[which(tmp$PP.H4.abf>=0.8),]
    print(nrow(tmp))
    
    tmp$gene<-gsub("\\.[0-9]{1,2}$","",tmp$phenotype_id)
    tmp<-merge(tmp,gene,by="gene",all.x=T,sort=F)

    # create a merge ID by combining gene and variant
    tmp$X<-paste(tmp$condition_name,tmp$phenotype_id,tmp$qtl_lead,sep="_")

    tmp$pheno_coloc<-pheno[i]
    
    if (i==1) {
        coloc<-tmp
    } else{
        coloc<-rbind(coloc,tmp)
    }
    rm(tmp)
}

table(coloc$pheno_coloc)
#   cd  ibd   uc 
# 8206 9138 7332

coloc<-as.data.frame(coloc)

table(coloc$study_label)
#           Alasoo_2018   Allegbe_Harris_2025             BLUEPRINT 
#                   195                 12173                   545 
# Bossini-Castillo_2019                   CAP                 CEDAR 
#                   140                   113                   364 
#            Cytoimmgen          Eldjarn_2023          Fairfax_2012 
#                   996                    40                    51 
#          Fairfax_2014               GENCORD              GEUVADIS 
#                   341                   237                   133 
#        Gilchrist_2021                  GTEx               Hu_2021 
#                    50                  1866                    48 
#           Kasela_2017     Kim-Hellmuth_2017            Lepik_2017 
#                    75                   162                   325 
#        Naranbhai_2015           Nathan_2022          Nedelec_2016 
#                    29                   553                   214 
#                OneK1K         Panousis_2025            Perez_2022 
#                   745                  3102                   133 
#            Quach_2016         Randolph_2021        Schmiedel_2018 
#                   613                    99                   894 
#                 sparc              Sun_2018               TwinsUK 
#                    18                    23                   399 

coloc[which(coloc$cohort=="pQTL_sparc"),]

table(coloc$cohort)
# eQTL_catalogue        hu_2021       IBDverse       macromap    pQTL_decode 
#           9295             48          12173           3102             40 
#     pQTL_sparc 
#             18 

coloc$MarkerName_gwas_local_lead<-gsub("_",":",coloc$gwas_lead)
coloc$MarkerName_eqtl_lead<-gsub("_",":",coloc$qtl_lead)

coloc$qtl_pval<-as.numeric(coloc$qtl_pval)
coloc$gwas_pval<-as.numeric(coloc$gwas_pval)


coloc$gwas_eqtl_lead_pairs<-paste0(coloc$MarkerName_gwas_local_lead,"_",coloc$MarkerName_eqtl_lead)

ld_tmp<-ld[which( (ld$ID1 %in% coloc$gwas_eqtl_lead_pairs) | (ld$ID2 %in% coloc$gwas_eqtl_lead_pairs)),]
dim(ld_tmp)
# [1] 7583    9

# ld[which(ld$SNP_A %in% c("chr5:159387525:A:C","chr5:159356215:G:T") & ld$SNP_B %in% c("chr5:159387525:A:C","chr5:159356215:G:T")),]


# merge list of colocs with ld data:

tmp1<-merge(coloc,ld_tmp[,c("ID1","R2")],by.x="gwas_eqtl_lead_pairs",by.y="ID1")
tmp2<-merge(coloc,ld_tmp[,c("ID2","R2")],by.x="gwas_eqtl_lead_pairs",by.y="ID2")

tmp3<-coloc[which(!coloc$gwas_eqtl_lead_pairs %in% ld_tmp$ID1 & !coloc$gwas_eqtl_lead_pairs %in% ld_tmp$ID2),]
tmp3$R2<-NA

coloc_results<-rbind(tmp1,tmp2,tmp3)
print(nrow(coloc_results)==nrow(coloc))
rm(tmp1,tmp2,tmp3)
nrow(coloc_results)
# [1] 24676

# fill in the cases where MarkerName_gwas_local_lead==MarkerName_eqtl_lead
coloc_results$R2[which(coloc_results$MarkerName_gwas_local_lead==coloc_results$MarkerName_eqtl_lead)]<-1
colnames(coloc_results)[colnames(coloc_results)=="R2"]<-"R2_MarkerName_gwas_local_lead_MarkerName_eqtl_lead"

rm(ld_tmp)

######################################################################################
## merge any GWAS lead variant, as defined by consensus leads, with eQTL variants

# LOAD THE INDEPENDENT SIGNALS - see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_38_2_define_independent_signals_using_as_ref_cojo_agnostic_results.R
all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 631 182


markernames<-c(all$MarkerName,coloc$MarkerName_eqtl_lead)
markernames<-markernames[!duplicated(markernames)]
length(markernames)
# [1] 7387

ld<-as.data.frame(ld)
ld_tmp<-ld[which( (ld$SNP_A %in% markernames) | (ld$SNP_B %in% markernames)),]
dim(ld_tmp)
# 1324328

ld_tmp1<-ld_tmp[which(ld_tmp$SNP_A %in% coloc$MarkerName_eqtl_lead & ld_tmp$SNP_B %in% all$MarkerName),]
ld_tmp1<-ld_tmp1[,c("SNP_A","SNP_B","R2")]
dim(ld_tmp1)
# [1] 104093       3

ld_tmp2<-ld_tmp[which(ld_tmp$SNP_B %in% coloc$MarkerName_eqtl_lead & ld_tmp$SNP_A %in% all$MarkerName),]
ld_tmp2<-ld_tmp2[,c("SNP_B","SNP_A","R2")]
dim(ld_tmp2)
# [1] 100010       3

rm(ld_tmp)

colnames(ld_tmp1)<-c("MarkerName_eqtl_lead","IIBDGC_GWAS_index_variant","R2")
colnames(ld_tmp2)<-c("MarkerName_eqtl_lead","IIBDGC_GWAS_index_variant","R2")

ld_tmp<-rbind(ld_tmp1,ld_tmp2)
rm(ld_tmp1,ld_tmp2)

ld_tmp<-ld_tmp[order(ld_tmp$R2,decreasing=T),]
ld_tmp_1<-ld_tmp[!duplicated(ld_tmp$MarkerName_eqtl_lead),]

dim(ld_tmp_1)
# [1] 5743    3

dim(coloc_results)
# [1] 24676   34


coloc_results<-merge(coloc_results,ld_tmp_1,by="MarkerName_eqtl_lead",all.x=T)
rm(ld_tmp_1)

coloc_results<-as.data.frame(coloc_results)
dim(coloc_results[which(coloc_results$MarkerName_eqtl_lead %in% all$MarkerName),])
# [1] 1253  36

# fill in those where lead eQTL is the same as the consensus index variant:
coloc_results$IIBDGC_GWAS_index_variant[which(coloc_results$MarkerName_eqtl_lead %in% all$MarkerName)]<-coloc_results$MarkerName_eqtl_lead[which(coloc_results$MarkerName_eqtl_lead %in% all$MarkerName)]
coloc_results$R2[which(coloc_results$MarkerName_eqtl_lead %in% all$MarkerName)]<-1

table(coloc_results$IIBDGC_GWAS_index_variant[which(coloc_results$MarkerName_eqtl_lead %in% all$MarkerName)]==coloc_results$MarkerName_eqtl_lead[which(coloc_results$MarkerName_eqtl_lead %in% all$MarkerName)])
# TRUE 
# 1253

colnames(coloc_results)[colnames(coloc_results)=="R2"]<-"R2_IIBDGC_GWAS_index_variant_MarkerName_eqtl_lead"

dim(coloc_results[which(!is.na(coloc_results$IIBDGC_GWAS_index_variant)),])
# [1] 21994    36

dim(coloc_results[which(!coloc_results$IIBDGC_GWAS_index_variant %in% all$MarkerName),])
# 2682   36
table(coloc_results$IIBDGC_GWAS_index_variant[which(!coloc_results$IIBDGC_GWAS_index_variant %in% all$MarkerName)],useNA="ifany")
# <NA> 
# 2682 

######################################################################################
## merge any GWAS lead variant, as defined by consensus leads, with local GWAS lead variants

markernames<-c(all$MarkerName,all$MarkerName_gwas_local_lead)
markernames<-markernames[!duplicated(markernames)]

ld<-as.data.frame(ld)
ld_tmp<-ld[which( (ld$SNP_A %in% markernames) | (ld$SNP_B %in% markernames)),]
dim(ld_tmp)
# 233238

ld_tmp1<-ld_tmp[which(ld_tmp$SNP_A %in% coloc$MarkerName_gwas_local_lead),]
ld_tmp1<-ld_tmp1[,c("SNP_A","SNP_B","R2")]
dim(ld_tmp1)
# [1] 64432       3
ld_tmp2<-ld_tmp[which(ld_tmp$SNP_B %in% coloc$MarkerName_gwas_local_lead),]
ld_tmp2<-ld_tmp2[,c("SNP_B","SNP_A","R2")]
dim(ld_tmp2)
# [1] 66558       3

rm(ld_tmp)

colnames(ld_tmp1)<-c("MarkerName_gwas_local_lead","IIBDGC_GWAS_index_variant","R2")
colnames(ld_tmp2)<-c("MarkerName_gwas_local_lead","IIBDGC_GWAS_index_variant","R2")

ld_tmp<-rbind(ld_tmp1,ld_tmp2)
rm(ld_tmp1,ld_tmp2)

ld_tmp<-ld_tmp[order(ld_tmp$R2,decreasing=T),]
ld_tmp_1<-ld_tmp[!duplicated(ld_tmp$MarkerName_gwas_local_lead),]

dim(ld_tmp_1)
# [1] 1133   3

dim(coloc_results)
# [1] 24676   38

coloc_results<-merge(coloc_results,ld_tmp_1,by="MarkerName_gwas_local_lead",all.x=T)

coloc_results<-as.data.frame(coloc_results)
dim(coloc_results[which(coloc_results$MarkerName_gwas_local_lead %in% all$MarkerName),])
# [1] 7177   38

coloc_results$IIBDGC_GWAS_index_variant.y[which(coloc_results$MarkerName_gwas_local_lead %in% all$MarkerName)]<-coloc_results$MarkerName_gwas_local_lead[which(coloc_results$MarkerName_gwas_local_lead %in% all$MarkerName)]
coloc_results$R2[which(coloc_results$MarkerName_gwas_local_lead %in% all$MarkerName)]<-1

table(coloc_results$IIBDGC_GWAS_index_variant.y[which(coloc_results$MarkerName_gwas_local_lead %in% all$MarkerName)]==coloc_results$MarkerName_gwas_local_lead[which(coloc_results$MarkerName_gwas_local_lead %in% all$MarkerName)])
# TRUE 
# 7177 

colnames(coloc_results)[colnames(coloc_results)=="R2"]<-"R2_IIBDGC_GWAS_index_variant_MarkerName_gwas_local_lead"

table(coloc_results$IIBDGC_GWAS_index_variant.y[which(!coloc_results$IIBDGC_GWAS_index_variant %in% all$MarkerName)],useNA="ifany")
# 0

# harmonise the IIBDG_GWAS_index variant
dim(coloc_results[which(!is.na(coloc_results$IIBDGC_GWAS_index_variant.x) & !is.na(coloc_results$IIBDGC_GWAS_index_variant.y)),])
# [1] 21993    38


dim(all[which(all$MarkerName %in% coloc_results$IIBDGC_GWAS_index_variant.x | all$MarkerName %in% coloc_results$IIBDGC_GWAS_index_variant.y),])
# [1] 534 192
dim(all)
# [1] 631 192

dim(all[which(!all$MarkerName %in% coloc_results$IIBDGC_GWAS_index_variant.x & !all$MarkerName %in% coloc_results$IIBDGC_GWAS_index_variant.y),])
# [1] 97 192

table(all$class_signal_final[which(!all$MarkerName %in% coloc_results$IIBDGC_GWAS_index_variant.x & !all$MarkerName %in% coloc_results$IIBDGC_GWAS_index_variant.y)])
# new_cojo_supervised_gw_significant_multiancestry_new_signal 
#                                                           3 
#                          new_cojo_unsupervised_known_signal 
#                                                          22 
#                            new_cojo_unsupervised_new_signal 
#                                                          72



########################

dim(coloc_results)
# [1] 24676   38

table(coloc_results$cohort)
# eQTL_catalogue        hu_2021       IBDverse       macromap    pQTL_decode 
#           9295             48          12173           3102             40 
#     pQTL_sparc 
#             18 


# retain only colocs from signals with leads (IIBDGC and eQTL) in LD:

# EXAMPLES FOR CD
dim(coloc_results[which(coloc_results$R2_IIBDGC_GWAS_index_variant_MarkerName_eqtl_lead>=R2_THRESHOLD),])
dim(coloc_results[which(coloc_results$R2_IIBDGC_GWAS_index_variant_MarkerName_eqtl_lead>=0.5),])
dim(coloc_results[which(coloc_results$R2_IIBDGC_GWAS_index_variant_MarkerName_eqtl_lead>=0.4),])
dim(coloc_results[which(coloc_results$R2_IIBDGC_GWAS_index_variant_MarkerName_eqtl_lead>=0.3),])

# > dim(coloc_results[which(coloc_results$R2_IIBDGC_GWAS_index_variant_MarkerName_eqtl_lead>=0.6),])
# [1] 10340    38
# > dim(coloc_results[which(coloc_results$R2_IIBDGC_GWAS_index_variant_MarkerName_eqtl_lead>=0.5),])
# [1] 11507    38
# > dim(coloc_results[which(coloc_results$R2_IIBDGC_GWAS_index_variant_MarkerName_eqtl_lead>=0.4),])
# [1] 12224    38
# > dim(coloc_results[which(coloc_results$R2_IIBDGC_GWAS_index_variant_MarkerName_eqtl_lead>=0.3),])
# [1] 12823    38

tmp1<-coloc_results[which(coloc_results$R2_IIBDGC_GWAS_index_variant_MarkerName_eqtl_lead>=R2_THRESHOLD),]
dim(tmp1)
# [1] 10340  38
tmp1<-tmp1[which(tmp1$IIBDGC_GWAS_index_variant.x %in% all$MarkerName),]
dim(tmp1)
# [1] 10340    38

summary(tmp1$R2_IIBDGC_GWAS_index_variant_MarkerName_eqtl_lead)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.6001  0.8604  0.9509  0.9071  0.9907  1.0000

tmp2<-coloc_results[which(coloc_results$R2_MarkerName_gwas_local_lead_MarkerName_eqtl_lead>=R2_THRESHOLD & coloc_results$R2_IIBDGC_GWAS_index_variant_MarkerName_gwas_local_lead>=R2_THRESHOLD),]
dim(tmp2)
# [1] 8296  38
tmp2<-tmp2[which(tmp2$IIBDGC_GWAS_index_variant.y %in% all$MarkerName),]
dim(tmp2)
# [1] 8296  38

summary(tmp2$R2_MarkerName_gwas_local_lead_MarkerName_eqtl_lead)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.6001  0.8451  0.9487  0.9083  0.9944  1.0000
summary(tmp2$R2_IIBDGC_GWAS_index_variant_MarkerName_gwas_local_lead)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.6103  0.9612  0.9996  0.9672  1.0000  1.0000



coloc_results<-rbind(tmp1,tmp2)
rm(tmp1,tmp2)
coloc_results<-coloc_results[!duplicated(coloc_results),]
dim(coloc_results)
# [1] 10420    38



summary(coloc_results$PP.H4.abf)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.8000  0.8858  0.9449  0.9312  0.9851  1.0000 


table(coloc_results$study_label)
#           Alasoo_2018   Allegbe_Harris_2025             BLUEPRINT 
#                    89                  5783                   200 
# Bossini-Castillo_2019                   CAP                 CEDAR 
#                    42                    41                   100 
#            Cytoimmgen          Eldjarn_2023          Fairfax_2012 
#                   316                    15                    15 
#          Fairfax_2014               GENCORD              GEUVADIS 
#                   122                    94                    50 
#        Gilchrist_2021                  GTEx               Hu_2021 
#                    15                   845                    25 
#           Kasela_2017     Kim-Hellmuth_2017            Lepik_2017 
#                    16                    42                   148 
#        Naranbhai_2015           Nathan_2022          Nedelec_2016 
#                     9                    83                    93 
#                OneK1K         Panousis_2025            Perez_2022 
#                   225                  1286                    45 
#            Quach_2016         Randolph_2021        Schmiedel_2018 
#                   205                    57                   294 
#              Sun_2018               TwinsUK 
#                    10                   155 

table(coloc_results$cohort)
# eQTL_catalogue        hu_2021       IBDverse       macromap    pQTL_decode 
#           3311             25           5783           1286             15

write.table(coloc_results,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_filtered_by_gwas_r2_0.6_allpheno_noFM.tsv"),
       col.names=T,row.names=F,quote=F,sep="\t")

q("no")



##############################################################################################################################
##############################################################################################################################

