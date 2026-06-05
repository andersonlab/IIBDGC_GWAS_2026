# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
## ADD COLOC RESULTS FROM SNAKEMAKE_COLOC PIPELINE
# See ~/git/snakemake_colocalisation

# cd /path/to/ibdgwas/IIBDGC/resources/otar/
# wget --recursive --no-parent --no-host-directories --cut-dirs 6 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output/association_overall_direct .
# wget --recursive --no-parent --no-host-directories --cut-dirs 6 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output/disease .
# wget --recursive --no-parent --no-host-directories --cut-dirs 6 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output/association_by_datasource_direct .

# singularity exec iibdgc_postprocess_10_singularity.sif

###############################################################

# 1.-  COMBINE THE RESULTS:

MEM=4000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -q cpu-interactive -G your_hpc_group R

library(data.table)
library(rtracklayer)
library(stringr)
library(arrow)


rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"


# LOAD THE INDEPENDENT SIGNALS - see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_38_2_define_independent_signals_using_as_ref_cojo_agnostic_results.R
all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 619 192

length(table(all$updated_region))
# [1] 362

table(all$phenotype,useNA="ifany")
            #  CD   IBD_saturated IBD_unsaturated              UC 
            # 128             104             269             118  

# LOAD THE MAP FOR COLOCALIZATION RESULTS - SEE ~/git/snakemake_colocalisation/scripts/create_maps_for_IIBDGC.R
map<-fread("/path/to/project",head=T)
map<-as.data.frame(map)
map<-map[which(map$cohort!="edQTL_GTEX"),]
map<-map[which(map$quant_method %in% c("aptamer","ge","microarray")),]

# LOAD THE COLOCALIZATION RESULTS - SEE ~/git/IIBDGC_GWAS/scripts/other/coloc_subset_results_by_ld_leads.R
coloc_results<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_filtered_by_gwas_r2_0.6_allpheno_noFM_with_beta.tsv",sep=""),head=T)
dim(coloc_results)
# [1] 10420    41

# LOAD THE NEAREST GENE - see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_44_add_closest_gene_names_to_variants.R
ng<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_closest_gene.tsv.gz",sep=""))
ng<-as.data.frame(ng)
ng<-merge(ng,all[,c("MarkerName","phenotype")],all.x=T)


# COMBINE BOTH
dat <- rbindlist(lapply(seq_len(nrow(all)), function(i) {
    # IIBDGC_GWAS_index_variant in LD with local GWAS lead:
    tmp1 <- coloc_results[which(coloc_results$IIBDGC_GWAS_index_variant.x==all$MarkerName[i] & coloc_results$R2_IIBDGC_GWAS_index_variant_MarkerName_eqtl_lead>=0.6),c("IIBDGC_GWAS_index_variant.x","pheno_coloc","qtl_lead","gene","gene_name","Beta_eqtl","SE_eqtl","Pvalue_eqtl","PP.H4.abf","cohort","study_label","condition_name","R2_IIBDGC_GWAS_index_variant_MarkerName_eqtl_lead")]
    # IIBDGC_GWAS_index_variant in LD with eQTL lead:
    tmp2 <- coloc_results[which(coloc_results$IIBDGC_GWAS_index_variant.y==all$MarkerName[i] & coloc_results$R2_IIBDGC_GWAS_index_variant_MarkerName_gwas_local_lead>=0.6),c("IIBDGC_GWAS_index_variant.y","pheno_coloc","qtl_lead","PP.H4.abf","Beta_eqtl","SE_eqtl","Pvalue_eqtl","cohort","gene","gene_name","study_label","condition_name","R2_IIBDGC_GWAS_index_variant_MarkerName_gwas_local_lead")]
    colnames(tmp1)[1] <- "IIBDGC_GWAS_index_variant"
    colnames(tmp2)[1] <- "IIBDGC_GWAS_index_variant"
    colnames(tmp1)[13] <- "R2"
    colnames(tmp2)[13] <- "R2"
    unique(rbind(tmp1, tmp2))
}))
dat <- unique(dat)

dat<-as.data.frame(dat)
dat$Y<-paste(dat$condition_name,dat$pheno_coloc,sep="_")

## add the score per gene:

# LOAD THE ENRICHMENT RESULTS - SEE ~/git/IIBDGC_GWAS/scripts/other/plot_enrichment_analyses.R
enrich<-fread("~/git/IIBDGC_GWAS/plots/paper_tables/Source_data_enrichment_analyses_eqtl_pqtl.tsv",head=T)
enrich<-as.data.frame(enrich)
enrich<-enrich[which(enrich$Category %in% map$id_map),]

enrich$pheno<-tolower(enrich$pheno)
enrich<-enrich[which(enrich$pheno %in% c("ibd","cd","uc")),]
dim(enrich)
enrich<-enrich[which(enrich$Enrichment_p<=(0.05)),]
dim(enrich)
# [1] 1472   20

enrich$Y<-paste(enrich$Category,enrich$pheno,sep="_")

dim(dat[which(!dat$condition_name %in% enrich$Category),])
# [1]  4 14
dim(dat[which(dat$condition_name %in% enrich$Category),])
# [1] 18472    14

dat<-merge(dat,enrich[,c("Y","Coefficient_zscore")],by="Y",all.x=T)
colnames(dat)[which(colnames(dat)=="Coefficient_zscore")]<-"tissue_coefficient_zscore"
dat<-as.data.frame(dat)

# relabel pqtl id:
dat$gene_name[which(is.na(dat$gene_name))]<-dat$gene[which(is.na(dat$gene_name))]

all<-as.data.frame(all)

all$eqtl_ibd<-""
all$eqtl_cd<-""
all$eqtl_uc<-""

all$eqtl_gene_ibd<-""
all$eqtl_gene_cd<-""
all$eqtl_gene_uc<-""

pheno<-c("ibd","cd","uc")

all_tmp<-all[,"MarkerName",drop=F]

for (i in 1:nrow(all_tmp)) {

    tmp_region<-dat[which(dat$IIBDGC_GWAS_index_variant==all_tmp$MarkerName[i]),]

    if (nrow(tmp_region)>0) {

        for (ph in pheno) {

            tmp<-tmp_region[which(tmp_region$pheno_coloc==ph),]

            if (nrow(tmp)>0) {
                tmp$eqtl<-paste(tmp$cohort,tmp$study_label,tmp$condition_name,tmp$tissue_cell_condition_label,tmp$gene_name,tmp$qtl_lead,formatC(tmp$Beta_eqtl,format="f",digits=2),formatC(tmp$SE_eqtl,format="f",digits=2),formatC(tmp$Pvalue_eqtl,format="E",digits=2),formatC(tmp$PP.H4.abf,format="f",digits=2),formatC(tmp$tissue_coefficient_zscore,format="f",digits=2),sep=";")
                tmp$eqtl_gene<-paste(tmp$gene_name[!duplicated(tmp$gene_name)],collapse="|")
                all[i,paste0("eqtl_",ph)]<-paste(tmp$eqtl,collapse="|")
                all[i,paste0("eqtl_gene_",ph)]<-paste(tmp$gene_name[!duplicated(tmp$gene_name)],collapse="|")
                rm(tmp)
            }

        }
    }

}

dim(all)
# [1] 619 198

dim(all[which(all$eqtl_ibd!="" | all$eqtl_cd!="" | all$eqtl_uc!=""),])
# [1] 305 188

all<-all[,c("MarkerName","phenotype","updated_region","BETA_cd_eur_tier_2","SE_cd_eur_tier_2","P-value_cd_eur_tier_2","BETA_ibd_eur_tier_2","P-value_ibd_eur_tier_2","BETA_uc_eur_tier_2","SE_uc_eur_tier_2","P-value_uc_eur_tier_2","eqtl_ibd","eqtl_gene_ibd","eqtl_cd","eqtl_uc","eqtl_gene_cd","eqtl_gene_uc")]

exo<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz",sep=""),
head=T)
exo<-exo[,c("MarkerName","exonic_variant_in_ld",
"exonic_variant_in_ld_gene_aac","exonic_variant_in_ld_status","exonic_variant_in_ld_associated_phenotype","class_signal_final_exome")]


all<-merge(all,exo,by="MarkerName",all.x=T)

mr<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_mr_results.tsv.gz",sep=""),
head=T)
colnames(mr)[which(colnames(mr) %in% colnames(all))]
mr<-mr[,c("MarkerName","mr_cd","mr_ibd","mr_uc","mr_gene_cd","mr_gene_ibd","mr_gene_uc")]

all<-merge(all,mr,by="MarkerName",all.x=T)

fwrite(all,paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants_effector_genes_coloc_mr.tsv.gz",sep=""),
col.names=T,row.names=F,quote=F,sep="\t")

# for how many signals we find a coloc, a MR, or coding

dim(all)
# [1] 619  28

## through MR:

dim(all[which((all$mr_gene_cd!="") | (all$mr_gene_ibd!="") | (all$mr_gene_uc!="")),])
#[1] 12 32

mr_genes<-c(all$mr_gene_cd,all$mr_gene_uc,all$mr_gene_ibd)
mr_genes<-mr_genes[mr_genes!=""]
mr_genes<-mr_genes[which(!is.na(mr_genes))]
mr_genes<-unlist(strsplit(mr_genes,"\\|"))
mr_genes<-mr_genes[!duplicated(mr_genes)]
length(mr_genes)
# [1] 18

## through ExWAS:

dim(all[which(all$exonic_variant_in_ld_gene_aac!=""),c("exonic_variant_in_ld_gene_aac"),drop=F])
# [1] 51  1

exo_genes<-c(all$exonic_variant_in_ld_gene_aac)
exo_genes<-exo_genes[exo_genes!=""]
length(exo_genes)
# [1] 51
exo_genes[which(exo_genes=="TSPO;TTLL12:TSPO:Thr147Ala")]<-"TSPO:Thr147Ala"
exo_genes[which(exo_genes=="MFNG:MFNG:Arg302Cys")]<-"MFNG:Arg302Cys"

exo_genes<-exo_genes[!duplicated(exo_genes)]
length(exo_genes)
exo_genes
#  [1] "PTPN22:W620R"           "SLAMF8:G99S"            "FCGR2A:H167R"          
#  [4] "PLEKHG5:A414T"          "IL23R:G149R"            "IL23R:V362I"           
#  [7] "IL23R:R381Q"            "NXPE1:G353R"            "IL10RA:P295L"          
# [10] "LACC1:I254V"            "PTGER2:C83G"            "GPR65:I231L"           
# [13] "ADCY7:D439E"            "NOD2:V793M"             "NOD2:R702W"            
# [16] "NOD2:R703C"             "NOD2:A755V"             "NOD2:G908R"            
# [19] "NOD2:p.Leu980ProfsTer2" "CARMIL2:V181M"          "PLCG2:H244R"           
# [22] "PLCG2:R268W"            "CCR7:M7V"               "TYK2:P1104A"           
# [25] "TYK2:A928V"             "TYK2:I684S"             "FUT2:p.Trp154Ter"      
# [28] "IFIH1:I923V"            "IFIH1:."                "GPR35:T108M"           
# [31] "GCKR:L446P"             "ARHGAP25:E254K"         "ADAM17:."              
# [34] "MFNG:Arg302Cys"         "TSPO:TSPO:Thr147Ala"    "MST1:p.Arg651Ter"      
# [37] "MST1:R703C"             "SLC39A8:A391T"          "DUSP1:A56T"            
# [40] "TRAF3IP2:R83W"          "ZC3H12D:K106R"          "ORC3:."                
# [43] "PTCD1:R113W"            "RCC1L:G30R"             "SHARPIN:S17F"          
# [46] "DOK2:L138S"             "DOK2:P274L"             "RIPK2:I259T"           
# [49] "CARD9:."                "CARD9:S12N"             "SPATA31F1:R1083S"    

length(exo_genes)
# [1] 51

exo_genes<-gsub(":.*","",exo_genes)
table(exo_genes)[order(table(exo_genes))]
exo_genes
#  [1] "PTPN22"    "SLAMF8"    "FCGR2A"    "PLEKHG5"   "IL23R"     "IL23R"    
#  [7] "IL23R"     "NXPE1"     "IL10RA"    "LACC1"     "PTGER2"    "GPR65"    
# [13] "ADCY7"     "NOD2"      "NOD2"      "NOD2"      "NOD2"      "NOD2"     
# [19] "NOD2"      "CARMIL2"   "PLCG2"     "PLCG2"     "CCR7"      "TYK2"     
# [25] "TYK2"      "TYK2"      "FUT2"      "IFIH1"     "IFIH1"     "GPR35"    
# [31] "GCKR"      "ARHGAP25"  "ADAM17"    "MFNG"      "TSPO"      "MST1"     
# [37] "MST1"      "SLC39A8"   "DUSP1"     "TRAF3IP2"  "ZC3H12D"   "ORC3"     
# [43] "ZNF655"    "RCC1L"     "SHARPIN"   "DOK2"      "DOK2"      "RIPK2"    
# [49] "CARD9"     "CARD9"     "SPATA31F1"


exo_genes<-exo_genes[!duplicated(exo_genes)]
length(exo_genes)
# [1] 37


### through coloc:
dim(all[which((all$eqtl_ibd!="") | (all$eqtl_cd!="") | (all$eqtl_uc!="")),])
# [1] 305  32
coloc_genes<-c(all$eqtl_gene_ibd,all$eqtl_gene_cd,all$eqtl_gene_uc)
coloc_genes<-coloc_genes[coloc_genes!=""]
coloc_genes<-coloc_genes[which(!is.na(coloc_genes))]
coloc_genes<-unlist(strsplit(coloc_genes,"\\|"))
coloc_genes<-coloc_genes[!duplicated(coloc_genes)]
length(coloc_genes)
# [1] 625

list_genes<-c(coloc_genes,mr_genes,exo_genes)
list_genes<-list_genes[!duplicated(list_genes)]
list_genes<-as.data.frame(list_genes)
list_genes<-list_genes[order(list_genes$list_genes),,drop=F]

list_genes$mr<-"0"
list_genes$mr[which(list_genes$list_genes %in% mr_genes)]<-"1"
table(list_genes$mr)
#   0   1 
# 646  18

list_genes$exonic<-"0"
list_genes$exonic[which(list_genes$list_genes %in% exo_genes)]<-"1"
table(list_genes$exonic)
#   0   1 
# 627  37 

list_genes$coloc<-0
list_genes$coloc[which(list_genes$list_genes %in% coloc_genes)]<-"1"
table(list_genes$coloc)
#   0   1 
#  39 625 

dim(list_genes)
# [1] 664   4

dim(list_genes[which( 
    (list_genes$mr==1 & list_genes$exonic==0 & list_genes$coloc==0) |
    (list_genes$mr==0 & list_genes$exonic==1 & list_genes$coloc==0) | 
    (list_genes$mr==0 & list_genes$exonic==0 & list_genes$coloc==1)),])
# [1] 648   4

664-648
# [1] 16

table(list_genes$mr,list_genes$coloc)
#       0   1
#   0  24 622
#   1  15   3

table(list_genes$exonic,list_genes$coloc)
#       0   1
#   0  15 612
#   1  24  13

table(list_genes$exonic,list_genes$mr)
#       0   1
#   0 609  18
#   1  37   0

list_genes[which(list_genes$mr==1 & list_genes$coloc==1),]
#     list_genes mr exonic coloc
# 273       CD28  1      0     1
# 283      PDCD1  1      0     1
# 370      TIMD4  1      0     1

list_genes[which(list_genes$exonic==1 & list_genes$coloc==1),]
#     list_genes mr exonic coloc
# 617      ADCY7  0      1     1
# 444      CARD9  0      1     1 - two signals, coloc exonic
# 250       FUT2  0      1     1
# 50       IL23R  0      1     1 - one signal, pqtl and eQTL, coding
# 132      LACC1  0      1     1
# 313       MFNG  0      1     1 - one signal, pqtl and eQTL, coding
# 344       MST1  0      1     1
# 84       NXPE1  0      1     1
# 578       ORC3  0      1     1
# 428      RCC1L  0      1     1
# 585      RIPK2  0      1     1
# 14      SLAMF8  0      1     1
# 319       TSPO  0      1     1

fwrite(list_genes,paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_noFM.tsv.gz",sep=""),
col.names=T,row.names=F,quote=F,sep="\t")


# for how many loci we find a coloc:
length(table(all$updated_region))
# [1] 362

length(table(all$updated_region[which((all$eqtl_ibd!="") | (all$eqtl_cd!="") | (all$eqtl_uc!=""))]))
# [1] 253

253/362
# 0.698895


## split per signal:

# exonic:
dim(all[which(all$exonic_variant_in_ld_gene_aac!=""),])
# [1] 51 32

# coloc:
dim(all[which((all$eqtl_ibd!="") | (all$eqtl_cd!="") | (all$eqtl_uc!="")),])
# [1] 305 32

# mr:
dim(all[which((all$mr_gene_cd!="") | (all$mr_gene_ibd!="") | (all$mr_gene_uc!="")),])
# [1] 12 32


dim(all[which( (all$exonic_variant_in_ld_gene_aac!="") |
(all$eqtl_ibd!="") | (all$eqtl_cd!="") | (all$eqtl_uc!="") | 
(all$mr_gene_cd!="") | (all$mr_gene_ibd!="") | (all$mr_gene_uc!="")) ,])
# [1] 341  32


## signals only exonic:
dim(all[which( ( (all$exonic_variant_in_ld_gene_aac!="")) 
& ((all$eqtl_ibd=="") & (all$eqtl_cd=="") & (all$eqtl_uc=="")) 
& ((all$mr_gene_cd=="") & (all$mr_gene_ibd=="") & (all$mr_gene_uc=="")) ),])
# [1] 31 32

## signals only coloc:
dim(all[which( ( (all$exonic_variant_in_ld_gene_aac=="")) 
& ((all$eqtl_ibd!="") | (all$eqtl_cd!="") | (all$eqtl_uc!="")) 
& ((all$mr_gene_cd=="") & (all$mr_gene_ibd=="") & (all$mr_gene_uc=="")) ),])
# [1] 279  32

## signals only MR:
dim(all[which( ((all$exonic_variant_in_ld_gene_aac=="")) 
& ((all$eqtl_ibd=="") & (all$eqtl_cd=="") & (all$eqtl_uc=="")) 
& ((all$mr_gene_cd!="") | (all$mr_gene_ibd!="") | (all$mr_gene_uc!="")) ),])
# [1] 4  32


# number of singals with evidence from > 1 method
341-31-279-4
# 27

# number of singals with evidence from only 1 method
341-27
# [1] 314

# rate of singals with evidence from only 1 method
314/341
# [1] 0.9215686


#######################################################################
## Option 1: create an unique list of genes_index variant

dat <- unique(rbindlist(lapply(pheno, function(ph) {
    tmp <- fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_matrices/matrix_",ph,"_zscore_genes_tissue_noFM.tsv.gz"),head=T)
    unique(as.data.frame(tmp[,c("IIBDGC_GWAS_index_variant","gene_name")]))
})))

# add phenotype as defined in file:
dat<-merge(dat,all[,c("MarkerName","phenotype")],by.x="IIBDGC_GWAS_index_variant",by.y="MarkerName",all.x=T)

# add nearest gene:
dat<-merge(dat,ng[,c("MarkerName","closest_gene")],by.x="IIBDGC_GWAS_index_variant",by.y="MarkerName",all.x=T)

dat_all_colocs<-dat
dim(dat_all_colocs)
# [1] 651   4

####################################################################################################
# Option 2: from coloc, only retain the gene with the largest tissue a score (rank1 in the tests)

# normalize score (values 0-1) per phenotype:
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

dat <- rbindlist(lapply(pheno, function(ph) {
    dat_final <- fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_matrices/matrix_",ph,"_zscore_genes_tissue_noFM.tsv.gz"),head=T)
    dat_final <- as.data.frame(dat_final[,c("IIBDGC_GWAS_index_variant","gene_name","score")])
    dat_final$score_norm <- scale_values(dat_final$score)
    dat_final
}))

dim(dat)
# 1956

# add phenotype as defined in file:
dat<-merge(dat,all[,c("MarkerName","phenotype")],by.x="IIBDGC_GWAS_index_variant",by.y="MarkerName",all.x=T)

# add nearest gene:
dat<-merge(dat,ng[,c("MarkerName","closest_gene")],by.x="IIBDGC_GWAS_index_variant",by.y="MarkerName",all.x=T)

length(table(dat$gene_name))
# [1] 625
length(table(dat$closest_gene))
# [1] 299



dat_colocs_rank1 <- as.data.table(dat)[, .SD[score == max(score)], by=IIBDGC_GWAS_index_variant]
dat_colocs_rank1 <- as.data.frame(dat_colocs_rank1)
dim(dat_colocs_rank1)
# [1] 321   6

# coloc_signals vs non coloc signals 
signals_coloc<-dat_all_colocs$IIBDGC_GWAS_index_variant[!duplicated(dat_all_colocs$IIBDGC_GWAS_index_variant)]
signals_no_coloc<-all$MarkerName[which(!all$MarkerName %in% signals_coloc)]

length(signals_coloc)
# [1] 305
length(signals_no_coloc)
# [1] 314

305+314
# [1] 619


# Codign signals and genes:
exo<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz",sep=""),
head=T)
exo<-as.data.frame(exo)

exo<-exo[which(exo$exonic_variant_in_ld_gene_aac!=""),]
exo$gene<-gsub(":.*","",exo$exonic_variant_in_ld_gene_aac)

exo<-exo[which(exo$gene!=""),c("MarkerName","gene")]

exo<-merge(exo,all[,c("MarkerName","phenotype")],by="MarkerName")


exo_genes<-exo$gene
# exo_genes<-exo_genes[!duplicated(exo_genes)]
length(exo_genes)
# [1] 51

############################
# MR signals and genes:
mr1<-mr[which(mr$mr_gene_cd!=""),c("MarkerName","mr_gene_cd")]
colnames(mr1)[2]<-"gene"
mr1_1<-mr1[!grepl("\\|",mr1$gene),]
mr1_2<-mr1[grepl("\\|",mr1$gene),]
mr1_2$gene1<-NA
mr1_2$gene2<-NA
mr1_2$gene1<-unlist(strsplit(mr1_2$gene,"\\|"))[c(1,3)]
mr1_2$gene2<-unlist(strsplit(mr1_2$gene,"\\|"))[c(2,4)]
mr1<-rbind(mr1_1,setNames(mr1_2[,c("MarkerName","gene1")],colnames(mr1_1)),setNames(mr1_2[,c("MarkerName","gene2")],colnames(mr1_1)))
rm(mr1_1,mr1_2)

mr2<-mr[which(mr$mr_gene_uc!=""),c("MarkerName","mr_gene_uc")]
colnames(mr2)[2]<-"gene"

mr3<-mr[which(mr$mr_gene_ibd!=""),c("MarkerName","mr_gene_ibd")]
colnames(mr3)[2]<-"gene"

mr<-rbind(mr1,mr2,mr3)
mr<-merge(mr,all[,c("MarkerName","phenotype")],by="MarkerName")

###############################

dim(dat_all_colocs)
# [1] 651   4

# dat_all_colocs<-dat_all_colocs[which(!dat_all_colocs$IIBDGC_GWAS_index_variant %in% exo$MarkerName),]
# dim(dat_all_colocs)
# # [1] 607   4

dim(dat_colocs_rank1)
# 1] 321   6
# dat_colocs_rank1<-dat_colocs_rank1[which(!dat_colocs_rank1$IIBDGC_GWAS_index_variant %in% exo$MarkerName),]
# dim(dat_colocs_rank1)
# # [1] 302   6

dim(dat)
# [1] 1953    6
# dat<-dat[which(!dat$IIBDGC_GWAS_index_variant %in% exo$MarkerName),]
# dim(dat)
# # [1] 1821    6

dim(ng)
# [1] 619   3

# coloc_signals vs non coloc signals 
signals_coloc<-dat_all_colocs$IIBDGC_GWAS_index_variant[!duplicated(dat_all_colocs$IIBDGC_GWAS_index_variant)]
signals_no_coloc<-all$MarkerName[which(!all$MarkerName %in% signals_coloc)]

length(signals_coloc)
# [1] 305
length(signals_no_coloc)
# [1] 314

305+314
# [1] 619


# ### GET THE NUMBER OF SIGNALS WITH ONLY ONE EFFECTOR GENE NOMINAGED:

# list_pheno_signals<-c("IBD_unsaturated","IBD_saturated","CD","UC")
# genes1<-dat[which(dat$phenotype %in% list_pheno_signals),c("IIBDGC_GWAS_index_variant","gene_name")]
# genes1$category<-"coloc"
# colnames(genes1)<-c("MarkerName","gene","category")
# genes2<-exo[which(exo$phenotype %in% list_pheno_signals),c("MarkerName","gene")]
# genes2$category<-"exonic"
# genes3<-mr[which(mr$phenotype %in% list_pheno_signals),c("MarkerName","gene")]
# genes3$category<-"mr"
# genes3<-as.data.frame(genes3)

# genes1<-genes1[!duplicated(genes1),]
# genes2<-genes2[!duplicated(genes2),]
# genes3<-genes3[!duplicated(genes3),]

# all_genes<-rbind(genes1,genes2,genes3)
# all_genes<-all_genes[!duplicated(all_genes),]

# length(all_genes$gene[!duplicated(all_genes$gene)])
# # [1] 664

# tmp<-all_genes$gene[!duplicated(all_genes$gene)]
# tmp[which(!tmp %in% list_genes$list_genes)]
# # 0


# all_genes<-all_genes[which(all_genes$gene!="TSPO;TTLL12"),]
# length(all_genes$gene[!duplicated(all_genes$gene)])
# # [1] 664 - OK

# dim(all_genes)
# # [1] 720   3

# # more than one gene:
# tmp<-all_genes$MarkerName[which(duplicated(all_genes$MarkerName))]
# tmp<-tmp[!duplicated(tmp)]
# length(tmp)
# # 169

# for (i in 1:length(tmp)) {

#     tmp_keep<-as.data.frame(table(all_genes$category[which(all_genes$MarkerName==tmp[i])]))
#     tmp_keep$N<-i

#     if (i==1) {
#         keep<-tmp_keep
#     } else {
#         keep<-rbind(keep,tmp_keep)
#     }

# }

# x<-as.data.frame.matrix(table(keep$N,keep$Var1))

# dim(x[which(x$coloc==1 & x$exonic==0 & x$mr==0),])
# # [1] 141   3

# dim(x[which(x$coloc==1 & (x$exonic==1 | x$mr==1)),])
# # [1] 26  3

# dim(x[which(x$coloc==0 & (x$exonic==1 | x$mr==1)),])
# # [1] 2 3

# 26+2+141
# # [1] 169


# x[which(x$coloc==1 & (x$exonic==1 | x$mr==1)),]

# # only one gene:
# tmp2<-all_genes$MarkerName[which(!all_genes$MarkerName %in% tmp)]
# length(tmp2)
# # [1] 172



for (ph in pheno) {

    print(ph)

    if (ph=="ibd") {
        list_pheno_signals<-c("IBD_unsaturated","IBD_saturated","CD","UC")
    } else if (ph=="cd") {
        list_pheno_signals<-c("IBD_unsaturated","IBD_saturated","CD")
    } else if (ph=="uc") {
        list_pheno_signals<-c("IBD_unsaturated","IBD_saturated","UC")
    }

    ######################### 

    print("list of all genes")

    genes1<-dat[which(dat$phenotype %in% list_pheno_signals),c("IIBDGC_GWAS_index_variant","gene_name")]
    genes2<-exo[which(exo$phenotype %in% list_pheno_signals),c("MarkerName","gene")]
    genes3<-mr[which(mr$phenotype %in% list_pheno_signals),c("MarkerName","gene")]

    genes<-rbind(genes1,setNames(genes2,colnames(genes1)),setNames(genes3,colnames(genes1)))
    rm(genes1,genes2,genes3)

    genes$uniq_pair<-paste(genes$IIBDGC_GWAS_index_variant,genes$gene_name,sep="_")
    genes<-genes$gene_name[!duplicated(genes$uniq_pair)]

    print(length(genes))

    write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_coloc_mr_exo_",ph),col.names=F,row.names=F,sep="\t")


    ######################### 

    print("list of all genes, coloc rank1")

    genes1<-dat_colocs_rank1[which(dat_colocs_rank1$phenotype %in% list_pheno_signals),c("IIBDGC_GWAS_index_variant","gene_name")]
    genes2<-exo[which(exo$phenotype %in% list_pheno_signals),c("MarkerName","gene")]
    genes3<-mr[which(mr$phenotype %in% list_pheno_signals),c("MarkerName","gene")]

    genes<-rbind(genes1,setNames(genes2,colnames(genes1)),setNames(genes3,colnames(genes1)))
    rm(genes1,genes2,genes3)

    genes$uniq_pair<-paste(genes$IIBDGC_GWAS_index_variant,genes$gene_name,sep="_")
    genes<-genes$gene_name[!duplicated(genes$uniq_pair)]

    print(length(genes))

    write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_coloc_rank1_mr_exo_",ph),col.names=F,row.names=F,sep="\t")

    ######################### 

    print("list of all coloc genes")

    genes<-dat_all_colocs[which(dat_all_colocs$phenotype %in% list_pheno_signals),]
    genes$uniq_pair<-paste(genes$IIBDGC_GWAS_index_variant,genes$gene_name,sep="_")
    genes<-genes$gene_name[!duplicated(genes$uniq_pair)]

    print(length(genes))

    write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_coloc_",ph),col.names=F,row.names=F,sep="\t")

    ######################### 

    print("list of all coloc genes, no coding")

    genes<-dat_all_colocs[which(dat_all_colocs$phenotype %in% list_pheno_signals & !dat_all_colocs$IIBDGC_GWAS_index_variant %in% exo$MarkerName),]
    genes$uniq_pair<-paste(genes$IIBDGC_GWAS_index_variant,genes$gene_name,sep="_")
    genes<-genes$gene_name[!duplicated(genes$uniq_pair)]

    print(length(genes))

    write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_coloc_no_coding_",ph),col.names=F,row.names=F,sep="\t")


    # ######################### 

    # print("list of all coloc genes by 0.5 score")

    # genes<-dat[which(dat$phenotype %in% list_pheno_signals),]
    # genes<-genes[which(genes$score_norm>=0.5),]
    # genes$uniq_pair<-paste(genes$IIBDGC_GWAS_index_variant,genes$gene_name,sep="_")
    # genes<-genes$gene_name[!duplicated(genes$uniq_pair)]
    # print(length(genes))

    # write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_coloc_score_0.5_",ph),col.names=F,row.names=F,sep="\t")

    # ######################### 

    # print("list of all coloc genes by 0.6 score")

    # genes<-dat[which(dat$phenotype %in% list_pheno_signals),]
    # genes<-genes[which(genes$score_norm>=0.6),]
    # genes$uniq_pair<-paste(genes$IIBDGC_GWAS_index_variant,genes$gene_name,sep="_")
    # genes<-genes$gene_name[!duplicated(genes$uniq_pair)]
    # print(length(genes))

    # write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_coloc_score_0.6_",ph),col.names=F,row.names=F,sep="\t")

    # ######################### 

    # print("list of all coloc genes by 0.7 score")

    # genes<-dat[which(dat$phenotype %in% list_pheno_signals),]
    # genes<-genes[which(genes$score_norm>=0.7),]
    # genes$uniq_pair<-paste(genes$IIBDGC_GWAS_index_variant,genes$gene_name,sep="_")
    # genes<-genes$gene_name[!duplicated(genes$uniq_pair)]
    # print(length(genes))

    # write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_coloc_score_0.7_",ph),col.names=F,row.names=F,sep="\t")


    # ######################### 

    # print("list of all coloc genes by 0.8 score")

    # genes<-dat[which(dat$phenotype %in% list_pheno_signals),]
    # genes<-genes[which(genes$score_norm>=0.8),]
    # genes$uniq_pair<-paste(genes$IIBDGC_GWAS_index_variant,genes$gene_name,sep="_")
    # genes<-genes$gene_name[!duplicated(genes$uniq_pair)]
    # print(length(genes))

    # write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_coloc_score_0.8_",ph),col.names=F,row.names=F,sep="\t")

    # ######################### 

    # print("list of all coloc genes by 0.9 score")

    # genes<-dat[which(dat$phenotype %in% list_pheno_signals),]
    # genes<-genes[which(genes$score_norm>=0.9),]
    # genes$uniq_pair<-paste(genes$IIBDGC_GWAS_index_variant,genes$gene_name,sep="_")
    # genes<-genes$gene_name[!duplicated(genes$uniq_pair)]
    # print(length(genes))

    # write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_coloc_score_0.9_",ph),col.names=F,row.names=F,sep="\t")

    # ######################### 

    # print("list of all coloc genes by 1 score")

    # genes<-dat[which(dat$phenotype %in% list_pheno_signals),]
    # genes<-genes[which(genes$score_norm>=1),]
    # genes$uniq_pair<-paste(genes$IIBDGC_GWAS_index_variant,genes$gene_name,sep="_")
    # genes<-genes$gene_name[!duplicated(genes$uniq_pair)]
    # print(length(genes))

    # write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_coloc_score_1_",ph),col.names=F,row.names=F,sep="\t")


    ######################### 
    
    print("list of rank1 genes")

    genes<-dat_colocs_rank1$gene_name[which(dat_colocs_rank1$phenotype %in% list_pheno_signals)]
    print(length(genes))

    write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_coloc_rank1_",ph),col.names=F,row.names=F,sep="\t")

    ######################### 
    
    print("list of rank1 genes, no coding")

    genes<-dat_colocs_rank1$gene_name[which(dat_colocs_rank1$phenotype %in% list_pheno_signals & !dat_colocs_rank1$IIBDGC_GWAS_index_variant %in% exo$MarkerName)]
    print(length(genes))

    write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_coloc_rank1_no_coding_",ph),col.names=F,row.names=F,sep="\t")

    ######################### 

    print("list of all closest genes, no coding")

    genes<-ng$closest_gene[which(ng$phenotype %in% list_pheno_signals & !ng$MarkerName %in% exo$MarkerName)]
    print(length(genes))

    write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_closest_no_coding_",ph),col.names=F,row.names=F,sep="\t")

    ######################### 

    print("list of all closest genes")

    genes<-ng$closest_gene[which(ng$phenotype %in% list_pheno_signals)]
    print(length(genes))

    write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_closest_",ph),col.names=F,row.names=F,sep="\t")

    ######################### 

    print("list of closest genes from coloc signals ")

    genes<-ng[which(ng$MarkerName %in% signals_coloc),]
    genes<-genes$closest_gene[which(genes$phenotype %in% list_pheno_signals)]
    print(length(genes))

    write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_closest_coloc_signals_",ph),col.names=F,row.names=F,sep="\t")

    ######################### 

    print("list of closest genes from non-coloc signals ")

    genes<-ng[which(ng$MarkerName %in% signals_no_coloc),]
    genes<-genes$closest_gene[which(genes$phenotype %in% list_pheno_signals)]
    print(length(genes))

    write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_closest_no_coloc_signals_",ph),col.names=F,row.names=F,sep="\t")

    ######################### 

    print("list of closest genes from coloc signals, no coding")

    genes<-ng[which(ng$MarkerName %in% signals_coloc & !ng$MarkerName %in% exo$MarkerName),]
    genes<-genes$closest_gene[which(genes$phenotype %in% list_pheno_signals)]
    print(length(genes))

    write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_closest_coloc_signals_no_coding_",ph),col.names=F,row.names=F,sep="\t")

    ######################### 

    print("list of closest genes from non-coloc signals, no coding")

    genes<-ng[which(ng$MarkerName %in% signals_no_coloc & !ng$MarkerName %in% exo$MarkerName),]
    genes<-genes$closest_gene[which(genes$phenotype %in% list_pheno_signals)]
    print(length(genes))

    write.table(genes,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/genes_closest_no_coloc_signals_no_coding_",ph),col.names=F,row.names=F,sep="\t")


    rm(list_pheno_signals)

}

# # [1] "ibd"
# [1] "list of all genes"
# [1] 708
# [1] "list of all genes, coloc rank1"
# [1] 383
# [1] "list of all coloc genes"
# [1] 651
# [1] "list of all coloc genes, no coding"
# [1] 606
# [1] "list of rank1 genes"
# [1] 321
# [1] "list of rank1 genes, no coding"
# [1] 301
# [1] "list of all closest genes, no coding"
# [1] 568
# [1] "list of all closest genes"
# [1] 619
# [1] "list of closest genes from coloc signals "
# [1] 305
# [1] "list of closest genes from non-coloc signals "
# [1] 314
# [1] "list of closest genes from coloc signals, no coding"
# [1] 286
# [1] "list of closest genes from non-coloc signals, no coding"
# [1] 282
# [1] "cd"
# [1] "list of all genes"
# [1] 587
# [1] "list of all genes, coloc rank1"
# [1] 317
# [1] "list of all coloc genes"
# [1] 538
# [1] "list of all coloc genes, no coding"
# [1] 494
# [1] "list of rank1 genes"
# [1] 263
# [1] "list of rank1 genes, no coding"
# [1] 244
# [1] "list of all closest genes, no coding"
# [1] 454
# [1] "list of all closest genes"
# [1] 501
# [1] "list of closest genes from coloc signals "
# [1] 248
# [1] "list of closest genes from non-coloc signals "
# [1] 253
# [1] "list of closest genes from coloc signals, no coding"
# [1] 230
# [1] "list of closest genes from non-coloc signals, no coding"
# [1] 224
# [1] "uc"
# [1] "list of all genes"
# [1] 564
# [1] "list of all genes, coloc rank1"
# [1] 301
# [1] "list of all coloc genes"
# [1] 521
# [1] "list of all coloc genes, no coding"
# [1] 481
# [1] "list of rank1 genes"
# [1] 254
# [1] "list of rank1 genes, no coding"
# [1] 238
# [1] "list of all closest genes, no coding"
# [1] 453
# [1] "list of all closest genes"
# [1] 491
# [1] "list of closest genes from coloc signals "
# [1] 247
# [1] "list of closest genes from non-coloc signals "
# [1] 244
# [1] "list of closest genes from coloc signals, no coding"
# [1] 232
# [1] "list of closest genes from non-coloc signals, no coding"
# [1] 221



## add the rank to the list of genes:

list_genes<-fread(paste0(path_gwas,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_noFM.tsv.gz"))
list_genes<-as.data.frame(list_genes)
dim(list_genes)
# [1] 664   4

# list_genes<-fread(paste0(path_gwas,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence_with_ligand_receptor_with_protein_complexes_with_monogenic_with_gps.tsv.gz"))
# list_genes<-as.data.frame(list_genes)
# dim(list_genes)

list_genes$coloc_rank1<-0
list_genes$coloc_rank1[which(list_genes$list_genes %in% dat_colocs_rank1$gene_name)]<-1
table(list_genes$coloc_rank1)
#   0   1 
# 352 312 

# list of genes:
dim(list_genes)
# [1] 664   5

dim(list_genes[which( (list_genes$mr!=0 | list_genes$exonic!=0) & list_genes$coloc!=0),])
# [1] 16 5


## add the signals:
dat_all_colocs<-dat_all_colocs[,1:3]
colnames(dat_all_colocs)<-colnames(exo)
dat_all<-rbind(exo,dat_all_colocs,mr)
dim(dat_all)
# [1] 722   3

variants_by_gene <- tapply(dat_all$MarkerName, dat_all$gene, function(x) paste(unique(x), collapse=";"))
list_genes$independent_index_variants <- variants_by_gene[list_genes$list_genes]

# genes nominated by more than one signal:
dim(list_genes[grep(";",list_genes$independent_index_variants),])
# [1] 37  6
list_genes[grep(";",list_genes$independent_index_variants),]
#          list_genes mr exonic coloc coloc_rank1
# 7             ABTB2  0      0     1           1
# 15            ADCY7  0      1     1           1
# 39            ASCL2  0      0     1           1
# 43             ATL3  0      0     1           1
# 56         C11orf21  0      0     1           0
# 66            CARD9  0      1     1           1 ####
# 84            CD244  0      0     1           1
# 85             CD28  1      0     1           1
# 89             CD81  0      0     1           0
# 125          DCBLD1  0      0     1           1
# 135            DOK2  0      1     0           0
# 143             EHF  0      0     1           1
# 205 ENSG00000272432  0      0     1           0
# 206 ENSG00000272449  0      0     1           0
# 236            ETS1  0      0     1           1 ####
# 253           FOXO1  0      0     1           1
# 273          GPRIN3  0      0     1           1
# 293           IFIH1  0      1     0           0
# 309           IL23R  0      1     1           1 ####
# 324            IRF8  0      0     1           1
# 333           ITLN1  0      0     1           0
# 357     LINC00266-1  0      0     1           0
# 376           MACO1  0      0     1           0
# 405            MST1  0      1     1           0 ####
# 425            NEU4  0      0     1           1
# 432            NOD2  0      1     0           0
# 456           PDCD1  1      0     1           1 ####
# 471          PLAAT3  0      0     1           0
# 475           PLCG2  0      1     0           0
# 504            PVT1  0      0     1           1
# 518             RHD  0      0     1           0
# 577           STK10  0      0     1           1
# 591           TIMD4  1      0     1           1 ####
# 596         TMEM50A  0      0     1           1
# 621         TSPAN32  0      0     1           0
# 628            TYK2  0      1     0           0
# 651           ZFP90  0      0     1           1

fwrite(list_genes,paste0(path_gwas,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM.tsv.gz"),col.names=T,row.names=F,quote=F,sep="\t")


q("no")

################################################################################################################################

path_gwas="/path/to/ibdgwas/IIBDGC/"

number_of_permutation=1000000
MEM=20000

tests=("genes_coloc" "genes_coloc_rank1" "genes_closest" "genes_closest_coloc_signals" "genes_closest_no_coloc_signals" "genes_coloc_mr_exo" "genes_coloc_rank1_mr_exo" "genes_closest_no_coding" "genes_coloc_rank1_no_coding" "genes_coloc_no_coding" "genes_closest_no_coloc_signals_no_coding" "genes_closest_coloc_signals_no_coding")

pheno=(ibd cd uc)


for ph in ${pheno[@]}
do
echo ${ph} && \
for test in ${tests[@]}
do
echo ${test} && wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/coloc_results/test/${test}_${ph}
done
done


# without NOD2

for ph in ${pheno[@]}
do
for test in ${tests[@]}
do
bsub -J"permut" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}post_imputation/2022/log/test_heritability_${test}_${ph}_nod2_exclussion_list_stderr \
-o ${path_gwas}post_imputation/2022/log/test_heritability_${test}_${ph}_nod2_exclussion_list_stdout \
"Rscript /path/to/user/git/IIBDGC_GWAS/scripts/other/burden_heritability_enrichment_no_filter_with_exclussion_list.R \
${path_gwas}post_imputation/2022/analysis/metaanalysis/coloc_results/test/${test}_${ph} \
${path_gwas}post_imputation/2022/analysis/metaanalysis/coloc_results/test/${test}_${ph}_nod2_exclussion_list_output \
${number_of_permutation} \
${ph} \
${path_gwas}post_imputation/2022/analysis/metaanalysis/coloc_results/test/nod2_exclussion_list > \
${path_gwas}post_imputation/2022/log/burden_heritability_enrichment_${test}_${ph}_nod2_exclussion_list.Rout"
done
done

# CONTINUE HERE

for ph in ${pheno[@]}
do
echo ${ph} && for test in ${tests[@]}
do
echo ${test} && \
tail -50 ${path_gwas}post_imputation/2022/log/test_heritability_${test}_${ph}_nod2_exclussion_list_stdout | grep "Successfully"
done
done

for ph in ${pheno[@]}
do
for test in ${tests[@]}
do
ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/coloc_results/test/${test}_${ph}_nod2_exclussion_list_output
done
done





####################################################################################
####################################################################################

# 2.-  USE OTAR TO LABEL :

MEM=4000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -q cpu-interactive -G your_hpc_group R

library(data.table)
library(rtracklayer)
library(stringr)
library(arrow)

rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

list_genes<-fread(paste0(path_gwas,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM.tsv.gz"))
list_genes<-as.data.frame(list_genes)


dir_list_files<-"/path/to/ibdgwas/IIBDGC/resources/otar/"

disease<-read_parquet(paste0(dir_list_files,"disease/disease.parquet"))

(disease[grep("bowel disease|Crohn|ulcerative colitis",disease$name),"name"])
#  [1] "Crohn jejunoileitis"                                                              
#  [2] "Crohn ileitis"                                                                    
#  [3] "gastroduodenal Crohn disease"                                                     
#  [4] "inflammatory bowel disease 1"                                                     
#  [5] "inflammatory bowel disease 5"                                                     
#  [6] "inflammatory bowel disease 10"                                                    
#  [7] "inflammatory bowel disease 13"                                                    
#  [8] "inflammatory bowel disease 14"                                                    
#  [9] "inflammatory bowel disease 17"                                                    
# [10] "inflammatory bowel disease 19"                                                    
# [11] "inflammatory bowel disease 25"                                                    
# [12] "inflammatory skin and bowel disease, neonatal, 2"                                 
# [13] "IL10-related early-onset inflammatory bowel disease"                              
# [14] "neonatal inflammatory skin and bowel disease"                                     
# [15] "Crohn jejunitis"                                                                  
# [16] "Crohn disease of the esophagus"                                                   
# [17] "inflammatory bowel disease (infantile ulcerative colitis) 31, autosomal recessive"
# [18] "inflammatory bowel disease 30"                                                    
# [19] "inflammatory bowel disease 29"                                                    
# [20] "Crohn's disease"                                                                  
# [21] "ulcerative colitis"                                                               
# [22] "inflammatory bowel disease"                                                       
# [23] "Crohn's colitis"                                                                  
# [24] "oral Crohn's disease"                                                             
# [25] "perianal Crohn's disease"                                                         
# [26] "small bowel Crohn's disease"                                                      
# [27] "inflammatory bowel disease, immunodeficiency, and encephalopathy"                 
# [28] "Autosomal recessive early-onset inflammatory bowel disease"

map_disease_ids<-disease[grep("bowel disease|Crohn|ulcerative colitis",disease$name),c("name","id")]
dim(map_disease_ids)
# [1] 28  2


read_otar_parquets <- function(subdir, disease_ids) {
    list_files <- setdiff(list.files(paste0(dir_list_files, subdir)), "_SUCCESS")
    rbindlist(lapply(list_files, function(f) {
        tmp <- read_parquet(paste0(dir_list_files, subdir, f))
        tmp[tmp$diseaseId %in% disease_ids, ]
    }))
}

distar <- read_otar_parquets("association_overall_direct/", map_disease_ids$id)

dim(distar)
# [1] 23403     4


# merge with gene ID:
gtf<-rtracklayer::import(paste(path_gwas,'post_imputation/2022/analysis/metaanalysis/annotation/gencode.v44.annotation.gtf.gz',sep=""))
gene<-as.data.frame(gtf)
gene<-gene[,c("gene_id","gene_name")]
gene<-gene[!duplicated(gene),]
rm(gtf)

gene$targetId<-gsub("\\.[0-9]{1,2}","",gene$gene_id)

distar<-merge(distar,gene,by="targetId",all.x=T)
dim(distar)
# [1] 23403     6

head(distar[duplicated(distar$targetId),])

# when duplicated entrances, retain the one with the larger score:
distar<-distar[order(distar$score,decreasing=T),]
distar<-distar[!duplicated(distar$gene_name),]

list_genes$otar_evidence<-0
list_genes$otar_evidence[which(list_genes$list_genes %in% distar$gene_name)]<-1

table(list_genes$otar_evidence)
#   0   1 
# 209 455 

table(list_genes$otar_evidence,list_genes$coloc_rank1)
#       0   1
#   0 123  86
#   1 229 226


distar <- read_otar_parquets("association_by_datasource_direct/", map_disease_ids$id)

dim(distar)
# [1] 28608     6

table(distar$datatypeId)
    #    animal_model genetic_association  genetic_literature          known_drug 
    #            2454                2745                  55                 435 
    #      literature      rna_expression    somatic_mutation 
    #           13161                9672                  86 

distar<-merge(distar,gene,by="targetId",all.x=T)
dim(distar)
# [1] 28608     6

distar<-distar[which(distar$datatypeId=="genetic_association"),]
dim(distar)
# [1] 2745    8

list_genes$otar_evidence_genetic_association<-0
list_genes$otar_evidence_genetic_association[which(list_genes$list_genes %in% distar$gene_name)]<-1


table(list_genes$otar_evidence_genetic_association)
#   0   1 
# 390 274 

table(list_genes$otar_evidence_genetic_association,list_genes$coloc_rank1)
#       0   1
#   0 205 185
#   1 147 127

dim(list_genes[which(list_genes$otar_evidence_genetic_association==0 & (list_genes$coloc_rank1==1 | list_genes$exonic==1 | list_genes$mr==1)),])
# [1] 203   8


# no otar evidence and high priority evidence
tmp0<-list_genes[which(list_genes$otar_evidence_genetic_association==1),]
dim(tmp0)
# [1] 274   8

tmp1<-list_genes[which(list_genes$otar_evidence_genetic_association==0 & (list_genes$coloc_rank1==1
 | list_genes$exonic==1 | list_genes$mr==1)),]
dim(tmp1)
# [1] 203   8

tmp2<-list_genes[which(list_genes$otar_evidence_genetic_association==0 & (list_genes$coloc_rank1==0
 & list_genes$exonic==0 & list_genes$mr==0)),]
dim(tmp2)
# [1] 187   8

vec<-c(tmp0$list_genes,tmp1$list_genes,tmp2$list_genes)

table(tmp1$otar_evidence)
#   0   1 
#  89 114 
table(tmp2$otar_evidence)
#   0   1 
# 120  67 

p <- c(114,89,67,120)
                                # Expected count in category 5
                                # is 1.86 < 5 ==> chi square approx.
chisq.test(as.table(p))   
#         Chi-squared test for given probabilities

# data:  as.table(p)
# X-squared = 18.267, df = 3, p-value = 0.0003875
as.table(p)
#   A   B   C   D 
# 114  89  67 120  


# OR = A*D/B*C
(112*120)/(89*67)
# [1] 2.253899

## reload again to keep track of the source of the evidence:
distar <- read_otar_parquets("association_by_datasource_direct/", map_disease_ids$id)

distar<-merge(distar,gene,by="targetId",all.x=T)
dim(distar)
# [1] 28608     6


otar_source_by_gene <- tapply(distar$datatypeId, distar$gene_name, function(x) paste(unique(x), collapse="|"))
list_genes$otar_evidence_source <- otar_source_by_gene[list_genes$list_genes]

list_genes[which(list_genes$list_genes %in% c("PDGFA","PDGFB")),]
#     list_genes mr exonic coloc coloc_rank1 otar_evidence
# 499      PDGFA  0      0     1           1             1
# 500      PDGFB  0      0     1           1             1
#     otar_evidence_genetic_association
# 499                                 0
# 500                                 1
#                              otar_evidence_source
# 499                       animal_model|literature
# 500 genetic_association|rna_expression|literature


list_genes[which(list_genes$list_genes %in% c("RGS14")),]
#     list_genes mr exonic coloc coloc_rank1 otar_evidence
# 564      RGS14  0      0     1           0             1
#     otar_evidence_genetic_association
# 564                                 1
#                              otar_evidence_source
# 564 genetic_association|literature|rna_expression

list_genes[which(list_genes$list_genes %in% c("KPNA1")),]
#     list_genes mr exonic coloc coloc_rank1 otar_evidence
# 375      KPNA1  0      0     1           1             1
#     otar_evidence_genetic_association otar_evidence_source
# 375                                 0           literature

fwrite(list_genes,paste0(path_gwas,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence.tsv.gz"),col.names=T,row.names=F,quote=F,sep="\t")

q("no")


###################################################################################################################

## showcase the results without NOD2, better background estimations not biased by the selection of NOD2; and the results including duplicated genes (nominated by different signals) - larger h2 and smaller pvalues; the overall results are similar otherwise



### CREATE AN UNIQUE FIGURE WITH TWO PANELS, INCLUDING THESE RESULTS AND THE RANK COMPARISON
/path/to/user/git/IIBDGC_GWAS/scripts/other/plot_comparisons_mean_burden_heritability_and_coloc_ranks_noFM.R


