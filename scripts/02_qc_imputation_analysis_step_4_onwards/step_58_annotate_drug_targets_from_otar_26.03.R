# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
### annotate drug targets:


# Files uploaded by Kyle:

# Retrieve OTAR data (run in shell) -  this is now in v25.12
# cd /path/to/ibdgwas/IIBDGC/resources/otar/v25.12/
# rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/25.12/output/known_drug . # deprecated in favour of clinical_report
# rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/25.12/output/drug_mechanism_of_action .
# rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/25.12/output/drug_molecule .
# rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/25.12/output/openfda_significant_adverse_drug_reactions .


# cd /path/to/ibdgwas/IIBDGC/resources/otar/26.03/
#  rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/26.03/output/disease .
# rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/26.03/output/clinical_report .
# rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/26.03/output/drug_mechanism_of_action .
# rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/26.03/output/drug_molecule .
# rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/26.03/output/openfda_significant_adverse_drug_reactions .
# rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/26.03/output/drug_warning .
# rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/26.03/output/clinical_target .
# rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/26.03/output/clinical_indication .

# Setup (run in shell before launching R):
# singularity exec iibdgc_postprocess_10_singularity.sif
MEM=8000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -q cpu-interactive -G your_hpc_group R


library(data.table)
library(stringr)
library(arrow)
library(ggplot2)
library(ggpubr)
library(gtools)
library(rtracklayer)
library(tidyr)
library(dplyr)

rm(list=ls())

# --- Paths ---
path         <- "/path/to/ibdgwas/IIBDGC/"
path_tables  <- paste0(path, "post_imputation/2022/analysis/final_tables/")
path_otar    <- "/path/to/ibdgwas/IIBDGC/resources/otar/26.03/"

# --- Helper: read all parquet files in a directory ---
read_parquet_dir <- function(dir_path) {
    files <- list.files(dir_path)
    files <- files[files != "_SUCCESS"]
    rbindlist(lapply(files, function(f) read_parquet(file.path(dir_path, f))))
}

list_genes <- fread(paste0(path_tables, "list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence_with_coloc_mr_direction_with_ligand_receptor_with_protein_complexes_with_monogenic.tsv.gz"))
dim(list_genes)
# [1] 664  20

otar <- fread(paste0(path_tables, "TherapyGroupGeneList_OTAR.csv"))
otar <- otar[otar$gene %in% list_genes$list_genes]
dim(otar)
# dim(otar): [1] 60 16; no duplicated genes

list_genes <- merge(list_genes, otar, by.x="list_genes", by.y="gene", all.x=TRUE)
dim(list_genes)

# link additional genes by protein complex:


prc <- fread("~/git/IIBDGC_GWAS/plots/paper_tables/Source_data_human_protein_complexes_corum.txt.gz")
table(prc$n_ibd_genes_complex)
#    0    1    2    3    4
# 4743  664   96    4    1

# extract the list of genes that contribute to those protein complexes
lgprc <- fread(paste0(path_tables, "protein_complexes_ibd_genes_corum.txt.gz"), head=TRUE)
lgprc <- lgprc[lgprc$protein_complex %in% prc$protein_complex]
lgprc <- lgprc[lgprc$ibd_gene == 1]
lgprc <- unique(lgprc)
lgprc <- lgprc[lgprc$protein_complex_drug != ""]
dim(lgprc)
# dim(lgprc): [1] 324   7
lgprc <- lgprc[lgprc$n_genes_complex < 5]
dim(lgprc)
# dim(lgprc): [1] 260   7

lgprc_agg <- lgprc[, .(
    protein_complex     = paste(protein_complex, collapse="|"),
    protein_complex_drug = paste(protein_complex_drug, collapse="|")
), by = list_genes]

list_genes <- merge(list_genes, lgprc_agg, by="list_genes", all.x=TRUE)
list_genes[is.na(list_genes$protein_complex),     "protein_complex"]      <- ""
list_genes[is.na(list_genes$protein_complex_drug), "protein_complex_drug"] <- ""

# add the gpr by Ron:

gps <- fread(paste0(path_tables, "IBD_Gene_GPS_machine_readable.csv"))
gps <- gps[gps$Phecode != "-99" & gps$Phecode != ""]

ids <- names(table(gps$Phecode))

# Reshape GPS and GPS_DOE to wide format, one row per gene
gps_wide_gps <- dcast(gps, genes ~ paste0("gps_",     Phecode), value.var="GPS")
gps_wide_doe <- dcast(gps, genes ~ paste0("gps_doe_", Phecode), value.var="GPS_DOE")
gps_wide     <- merge(gps_wide_gps, gps_wide_doe, by="genes")

list_genes <- merge(list_genes, gps_wide, by.x="list_genes", by.y="genes", all.x=TRUE)

fwrite(list_genes, paste0(path_tables, "list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence_with_coloc_mr_direction_with_ligand_receptor_with_protein_complexes_with_monogenic_with_gps.tsv.gz"),
       col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

list_genes<-as.data.frame(list_genes)

###################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################


#### EXPLORE DIRECTION SUPPORT FROM OTAR


###########################
# load OTAR evidence:

disease <- read_parquet(file.path(path_otar, "disease/disease.parquet"))

map_disease_ids <- disease[grep("bowel disease|Crohn|ulcerative colitis|crohn", disease$name), c("name", "id")]
dim(map_disease_ids)
# dim(map_disease_ids): [1] 28  2

names(table(map_disease_ids$name))
#  [1] "Autosomal recessive early-onset inflammatory bowel disease"                       
#  [2] "Crohn disease of the esophagus"                                                   
#  [3] "Crohn ileitis"                                                                    
#  [4] "Crohn jejunitis"                                                                  
#  [5] "Crohn jejunoileitis"                                                              
#  [6] "Crohn's colitis"                                                                  
#  [7] "Crohn's disease"                                                                  
#  [8] "gastroduodenal Crohn disease"                                                     
#  [9] "IL10-related early-onset inflammatory bowel disease"                              
# [10] "inflammatory bowel disease"                                                       
# [11] "inflammatory bowel disease (infantile ulcerative colitis) 31, autosomal recessive"
# [12] "inflammatory bowel disease 1"                                                     
# [13] "inflammatory bowel disease 10"                                                    
# [14] "inflammatory bowel disease 13"                                                    
# [15] "inflammatory bowel disease 14"                                                    
# [16] "inflammatory bowel disease 17"                                                    
# [17] "inflammatory bowel disease 19"                                                    
# [18] "inflammatory bowel disease 25"                                                    
# [19] "inflammatory bowel disease 29"                                                    
# [20] "inflammatory bowel disease 30"                                                    
# [21] "inflammatory bowel disease 5"                                                     
# [22] "inflammatory bowel disease, immunodeficiency, and encephalopathy"                 
# [23] "inflammatory skin and bowel disease, neonatal, 2"                                 
# [24] "neonatal inflammatory skin and bowel disease"                                     
# [25] "oral Crohn's disease"                                                             
# [26] "perianal Crohn's disease"                                                         
# [27] "small bowel Crohn's disease"                                                      
# [28] "ulcerative colitis

#####################
# CLINICAL TARGET (DRUG) OBJECT

ct <- read_parquet_dir(file.path(path_otar, "clinical_target"))
ct <- tidyr::unnest_longer(as.data.frame(ct), clinicalReportIds)
dim(ct)
# [1] 662733     6

colnames(ct)
# [1] "id"                "drugId"            "targetId"         
# [4] "diseases"          "clinicalReportIds" "maxClinicalStage" 

drugs<-ct[,c("drugId","targetId")]
head(drugs)
        #   drugId        targetId
#           <char>          <char>
# 1: CHEMBL2105741 ENSG00000126934

drugs<-drugs[!duplicated(drugs),]

#####################
# 2. OpenFDA significant adverse drug reactions:

ae <- read_parquet_dir(file.path(path_otar, "openfda_significant_adverse_drug_reactions"))

drugs$adverse_event<-0
drugs$adverse_event[which(drugs$drugId %in% ae$chembl_id)]<-1

table(ae$event[grep("colitis",ae$event)])

vec<-names(table(ae[grepl("colitis|gastrointestinal bleeding",ae$event),"event"]))
# vec<-vec[!grepl("cytomegalovirus|enterocolitis|cytomegalovirus|infectious|campylobacter|clostridium|antibiotic|herpes|microscopic|collagenous|ischaemic|necrotising|diverticular|pseudomembranous",vec)]
vec<-vec[grepl("autoimmune colitis|immune-mediated enterocolitis|eosinophilic colitis|neutropenic colitis",vec)]
vec
# [1] "autoimmune colitis" "immune-mediated enterocolitis"

drugs$adverse_event_colitis<-0
drugs$adverse_event_colitis[which(drugs$drugId %in% ae$chembl_id[grepl(paste(vec,collapse="|"),ae$event)])]<-1


table(drugs$adverse_event,drugs$adverse_event_colitis)
#        0    1
#   0 6783    0
#   1 6520  104

drugs[which(drugs$adverse_event_colitis==1),]
#           drugId        targetId adverse_event adverse_event_colitis

dim(drugs)
# [1] 13407     4
drugs<-drugs[!duplicated(drugs),]
dim(drugs)
[1] 13407     4

# map gene symbol:
gtf  <- import(paste0(path, "post_imputation/2022/analysis/metaanalysis/annotation/gencode.v44.annotation.gtf.gz"))
gene <- as.data.frame(gtf)
gene <- gene[!gene$seqnames %in% c("chrX", "chrY", "chrM"), ]
gene <- gene[gene$type == "exon", c("gene_name", "gene_type", "gene_id")]
gene <- gene[!duplicated(gene), ]
gene <- gene[gene$gene_type == "protein_coding", ]

gene$targetId <- gsub("\\.[0-9]*$","",gene$gene_id)

#########################################################################################################
#### CLINICAL REPORT OBJECT

cr <- read_parquet_dir(file.path(path_otar, "clinical_report"))
colnames(cr)
#  [1] "id"                        "clinicalStage"            
#  [3] "type"                      "source"                   
#  [5] "drugs"                     "trialDescription"         
#  [7] "countries"                 "phaseFromSource"          
#  [9] "trialOfficialTitle"        "diseases"                 
# [11] "trialNumberOfArms"         "hasExpertReview"          
# [13] "trialOverallStatus"        "year"                     
# [15] "trialStopReasonCategories" "trialLiterature"          
# [17] "sideEffects"               "trialPhase"               
# [19] "trialWhyStopped"           "url"                      
# [21] "trialStudyType"            "qualityControls"          
# [23] "trialPrimaryPurpose"       "trialStartDate"           
# [25] "title" 

table(cr$trialOverallStatus)
#     ACTIVE_NOT_RECRUITING    APPROVED_FOR_MARKETING                 AVAILABLE 
#                      9139                       203                       218 
#                 COMPLETED   ENROLLING_BY_INVITATION       NO_LONGER_AVAILABLE 
#                    128409                       892                       457 
#        NOT_YET_RECRUITING                RECRUITING                 SUSPENDED 
#                      7995                     22384                       682 
# TEMPORARILY_NOT_AVAILABLE                TERMINATED                   UNKNOWN 
#                        23                     19674                     28175 
#                 WITHDRAWN 
#                      8034


cr_tox<-cr[,c("id","sideEffects")]
cr_tox<-cr_tox[which(!is.na(cr_tox$sideEffects)),]
cr_tox_expanded <- tidyr::unnest(cr_tox, sideEffects)
dim(cr_tox_expanded)
# [1] 4904    3
cr_tox_expanded<-as.data.frame(cr_tox_expanded)

cr_drug<-cr[,c("id","drugs")]
cr_drug_expanded<-tidyr::unnest(as.data.frame(cr_drug), drugs)

cr_drug_expanded<-as.data.frame(cr_drug_expanded)


cr_tox<-merge(cr_tox_expanded,cr_drug_expanded,by="id",all.x=T)
dim(cr_tox)
# [1] 4904    4

rm(cr_tox_expanded,cr_drug_expanded,cr_drug)


## extract clinical trial status:
cr_st<-cr[,c("id","trialOverallStatus")]
cr_st<-cr_st[!is.na(cr_st$trialOverallStatus),]

table(cr_st$trialOverallStatus)
#     ACTIVE_NOT_RECRUITING    APPROVED_FOR_MARKETING                 AVAILABLE 
#                      9139                       203                       218 
#                 COMPLETED   ENROLLING_BY_INVITATION       NO_LONGER_AVAILABLE 
#                    128409                       892                       457 
#        NOT_YET_RECRUITING                RECRUITING                 SUSPENDED 
#                      7995                     22384                       682 
# TEMPORARILY_NOT_AVAILABLE                TERMINATED                   UNKNOWN 
#                        23                     19674                     28175 
#                 WITHDRAWN 
#                      8034 

#########################################################################################################
#### CLINICAL INDICATION OBJECT

ci <- read_parquet_dir(file.path(path_otar, "clinical_indication"))
ci <- tidyr::unnest_longer(as.data.frame(ci), clinicalReportIds)

dim(ci)
# [1] 228632      5

# merge with drugs:
ci<-merge(ci,drugs,by="drugId",all.x=T)
dim(ci)
# [1] 662829      8

# # Exclude drugs not linked to drug targets by OTAR
# ci<-ci[!is.na(ci$targetId),]
# dim(ci)
# # 599344      8

# merge with genes:
ci<-merge(ci,gene[,c("gene_name","targetId")],by="targetId",all.x=T)
dim(ci)
# [1] 662829      9



ibd<-ci[which(ci$diseaseId %in% map_disease_ids$id),]
dim(ibd)
# [1] 3340    9

dim(ibd[!duplicated(ibd$targetId) & !is.na(ibd$targetId),])
# [1] 229   9 - uniq genes

dim(ibd[!duplicated(ibd$drugId),])
# [1] 283   9 - uniq drugs

ibd$pair<-paste(ibd$drugId,ibd$targetId,sep="_")
dim(ibd[!duplicated(ibd$pair) & !is.na(ibd$targetId),])
# [1] 369   9 - uniq pairs drug-target


table(ibd$maxClinicalStage)
    #  APPROVAL EARLY_PHASE_1       PHASE_1     PHASE_1_2       PHASE_2 
    #      1906            27            77            43           546 
    # PHASE_2_3       PHASE_3   PREAPPROVAL   PRECLINICAL       UNKNOWN 
    #        98           595             1            11            36 

ibd$maxClinicalStage_rec<-NA
ibd$maxClinicalStage_rec[which(ibd$maxClinicalStage=="APPROVAL")]<-4
ibd$maxClinicalStage_rec[which(ibd$maxClinicalStage %in% c("PHASE_3","PHASE_2_3"))]<-3
ibd$maxClinicalStage_rec[which(ibd$maxClinicalStage %in% c("PHASE_2","PHASE_1_2"))]<-2
ibd$maxClinicalStage_rec[which(ibd$maxClinicalStage %in% c("PHASE_1","EARLY_PHASE_1"))]<-1
# ibd$maxClinicalStage_rec[which(ibd$maxClinicalStage %in% c("UNKNOWN"))]<-0


table(ibd$maxClinicalStage_rec,useNA="ifany")
#    1    2    3    4 <NA> 
#  104  589  693 1906   48 

# retain only one drug_pair,pair with highest phase:
tmp<-ibd[which(!is.na(ibd$maxClinicalStage_rec)),c("drugId","gene_name","targetId","pair","maxClinicalStage_rec")]
tmp<-tmp[order(tmp$maxClinicalStage_rec,tmp$pair,decreasing=T),]

tmp<-tmp[!duplicated(tmp$pair),]

table(tmp$maxClinicalStage_rec)
#   1   2   3   4 
#  33 281  92  73 

dim(tmp)
# [1] 479   4

tmp$phase_3_4<-0
tmp$phase_3_4[which(tmp$maxClinicalStage_rec>=3)]<-1

tmp$study_support<-0
tmp$study_support[which(tmp$gene_name %in% list_genes$list_genes)]<-1
table(tmp$study_support)
#   0   1 
# 435  44 

tmp$study_support<-0
tmp$study_support[which(tmp$gene_name %in% list_genes$list_genes[which(list_genes$coloc_rank1==1 | list_genes$exonic==1 | list_genes$mr==1)])]<-1
table(tmp$study_support)
#   0   1 
# 440  39


table(tmp$phase_3_4,tmp$study_support)
#       0   1
#   0 216  20
#   1 110  19

fisher.test(table(tmp$phase_3_4,tmp$study_support),alternative = "greater")

#         Fisher's Exact Test for Count Data

# data:  table(tmp$phase_3_4, tmp$study_support)
# p-value = 0.03949
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  1.038871      Inf
# sample estimates:
# odds ratio 
#   1.910307  

table(tmp$gene_name[which(tmp$study_support==1)],tmp$maxClinicalStage_rec[which(tmp$study_support==1)])
#            1 2 3 4
#   CXCR1    0 2 0 0
#   IL2RA    0 2 0 0
#   IL2RB    0 2 0 0
#   ITGA4    0 1 0 2
#   JAK2     0 1 4 1
#   MAP3K8   0 1 0 0
#   NDUFA4L2 0 1 0 0
#   PDCD1    0 1 0 0
#   PRKCB    0 1 0 0
#   PTGS2    1 0 1 6
#   TEC      0 1 0 0
#   TNFSF11  0 1 0 0
#   TXK      0 1 0 0
#   TYK2     0 4 4 1

dim(list_genes[which(list_genes$list_genes %in% ibd$gene_name),])
# [1] 17 59

list_genes$list_genes[which(list_genes$list_genes %in% ibd$gene_name)]
#  [1] "CCR5"     "CXCR1"    "CXCR2"    "IL2RA"    "IL2RB"    "ITGA4"   
#  [7] "JAK2"     "MAP3K8"   "NDUFA4L2" "PDCD1"    "PRKCB"    "PTGER4"  
# [13] "PTGS2"    "TEC"      "TNFSF11"  "TXK"      "TYK2" 

# missing TNFSF11 TXK from Kyle's merge, why???


names(table(ibd$gene_name))
#   [1] "ADORA2A"  "ADORA2B"  "ALOX5"    "CACNA1C"  "CACNA1D"  "CACNA1S" 
#   [7] "CCL11"    "CCR5"     "CCR9"     "CD14"     "CD3D"     "CD3E"    
#  [13] "CD3G"     "CD80"     "CD86"     "CNR1"     "CNR2"     "CRBN"    
#  [19] "CSF2RB"   "CSF3R"    "CUL4A"    "CX3CL1"   "CXCL10"   "CXCR1"   
#  [25] "CXCR2"    "DDB1"     "DHFR"     "DHODH"    "ELANE"    "ESR2"    
#  [31] "FFAR2"    "FKBP1A"   "GHR"      "GLP2R"    "GPD2"     "GPR84"   
#  [37] "GUCY2C"   "HDAC1"    "HDAC10"   "HDAC11"   "HDAC2"    "HDAC3"   
#  [43] "HDAC4"    "HDAC5"    "HDAC7"    "HDAC9"    "HMGCR"    "ICAM1"   
#  [49] "IFNAR1"   "IFNAR2"   "IFNG"     "IGF1R"    "IL12A"    "IL12B"   
#  [55] "IL13"     "IL17A"    "IL17RA"   "IL18"     "IL1A"     "IL1B"    
#  [61] "IL1RAP"   "IL1RL2"   "IL23A"    "IL2RA"    "IL2RB"    "IL4R"    
#  [67] "IL6"      "IMPDH1"   "IMPDH2"   "ITGA2"    "ITGA4"    "ITGB1"   
#  [73] "ITGB7"    "ITK"      "JAK1"     "JAK2"     "JAK3"     "MADCAM1" 
#  [79] "MAP3K8"   "MAPK14"   "MC2R"     "MMP1"     "MMP13"    "MMP7"    
#  [85] "MMP8"     "MMP9"     "MS4A1"    "MTNR1A"   "MTNR1B"   "NCBP1"   
#  [91] "NCBP2"    "NDUFA10"  "NDUFA11"  "NDUFA12"  "NDUFA13"  "NDUFA2"  
#  [97] "NDUFA3"   "NDUFA4"   "NDUFA4L2" "NDUFA5"   "NDUFA6"   "NDUFA7"  
# [103] "NDUFA8"   "NDUFA9"   "NDUFAB1"  "NDUFAF1"  "NDUFAF2"  "NDUFAF3" 
# [109] "NDUFAF4"  "NDUFB1"   "NDUFB10"  "NDUFB2"   "NDUFB3"   "NDUFB4"  
# [115] "NDUFB5"   "NDUFB6"   "NDUFB7"   "NDUFB8"   "NDUFB9"   "NDUFC1"  
# [121] "NDUFC2"   "NDUFS1"   "NDUFS2"   "NDUFS3"   "NDUFS4"   "NDUFS5"  
# [127] "NDUFS6"   "NDUFS7"   "NDUFS8"   "NDUFV1"   "NDUFV2"   "NDUFV3"  
# [133] "NR3C1"    "NR3C2"    "NTSR1"    "OPRD1"    "OPRK1"    "OPRM1"   
# [139] "OSMR"     "PDCD1"    "PDE10A"   "PDE1A"    "PDE1B"    "PDE1C"   
# [145] "PDE2A"    "PDE3A"    "PDE3B"    "PDE4A"    "PDE4B"    "PDE4C"   
# [151] "PDE4D"    "PDE5A"    "PDE6A"    "PDE6B"    "PDE6C"    "PDE6D"   
# [157] "PDE6G"    "PDE6H"    "PDE7A"    "PDE7B"    "PDE8A"    "PDE8B"   
# [163] "PDE9A"    "PIKFYVE"  "PPARA"    "PPARG"    "PPAT"     "PRKCA"   
# [169] "PRKCB"    "PRKCD"    "PRKCE"    "PRKCG"    "PRKCH"    "PRKCI"   
# [175] "PRKCQ"    "PRKCZ"    "PRKD1"    "PRKD3"    "PTGER4"   "PTGS1"   
# [181] "PTGS2"    "RARA"     "RARB"     "RBX1"     "REN"      "RIPK1"   
# [187] "S1PR1"    "S1PR2"    "S1PR3"    "S1PR4"    "S1PR5"    "SIRT1"   
# [193] "SLC5A2"   "SLC6A2"   "SLC6A3"   "SLC6A4"   "SMAD7"    "TEC"     
# [199] "TLR2"     "TLR9"     "TNF"      "TNFRSF4"  "TNFSF11"  "TNFSF14" 
# [205] "TPH1"     "TPH2"     "TXK"      "TYK2"     "VDR"      "XDH"   

ibd[which(ibd$gene_name %in% c("TEC","TXK")),]
#                    label approvedSymbol                approvedName targetClass
# 12360    Crohn's disease            TEC tec protein tyrosine kinase      Enzyme
# 31238 ulcerative colitis            TEC tec protein tyrosine kinase      Enzyme
#           prefName tradeNames                                 synonyms
# 12360 RITLECITINIB            PF-06651600 , Pf-06651600 , Ritlecitinib
# 31238 RITLECITINIB            PF-06651600 , Pf-06651600 , Ritlecitinib
#             drugType           mechanismOfAction        targetName
# 12360 Small molecule TEC family kinase inhibitor TEC family kinase
# 31238 Small molecule TEC family kinase inhibitor TEC family kinase
#                     name adverse_event adverse_event_colitis
# 12360    Crohn's disease             0                     0
# 31238 ulcerative colitis             0                     0


dim(ibd[!duplicated(ibd$gene_name),])
# [1] 211  11


dim(list_genes[which(list_genes$list_genes %in% ibd$gene_name),])
# [1] 17 59

dim(list_genes[which(list_genes$list_genes %in% ci$gene_name),])
# [1] 62 59

list_genes$list_genes[which(list_genes$list_genes %in% ci$gene_name)]

head(tmp)
tmp[duplicated(tmp$pair),]
#0

# save this table:
fwrite(tmp, paste0(path_tables, "list_drugs_genes_pairs_ibd_clinical_trials.tsv.gz"),
       col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

ibd<-tmp

#########################################################################################################
#### CLINICAL INDICATION OBJECT

ma <- read_parquet_dir(file.path(path_otar, "drug_mechanism_of_action"))
ma <- ma[,c("actionType","mechanismOfAction","chemblIds","targetType")]

ma_expanded <- tidyr::unnest(ma, chemblIds)
dim(ma_expanded)
# [1] 7921    3
rm(ma)


# contains duplicated mechanismOfAction per drug
ma_expanded[which(ma_expanded$chemblIds=="CHEMBL206815"),]
# # A tibble: 2 × 4
#   actionType mechanismOfAction                     chemblIds    targetType    
#   <chr>      <chr>                                 <chr>        <chr>         
# 1 INHIBITOR  ADAM17 inhibitor                      CHEMBL206815 single protein
# 2 INHIBITOR  Matrix metalloproteinase 13 inhibitor CHEMBL206815 single protein

############# compile a summarised version of drugs targetting the gene 

otar<-ci[which(ci$gene_name %in% list_genes$list_genes[which(list_genes$coloc_rank1==1 | list_genes$exonic==1 | list_genes$mr==1)]),]
dim(otar)

otar<-merge(otar,ma_expanded,by.x="drugId",by.y="chemblIds",all.x=T)
table(otar$actionType,useNA="ifany")

otar$mechanism_action_hm<-"other"

names(otar$mechanismOfAction[which(otar$actionType=="EXOGENOUS PROTEIN")])
# [1] "Urokinase-type plasminogen activator exogenous protein"
names(table(otar$mechanismOfAction[which(otar$actionType=="OPENER")]))
# [1] "Vanilloid receptor opener"

otar$mechanism_action_hm[which(otar$actionType %in% c("ACTIVATOR","AGONIST","PARTIAL AGONIST","EXOGENOUS PROTEIN","POSITIVE ALLOSTERIC MODULATOR"))]<-"activator"
otar$mechanism_action_hm[which(otar$actionType %in% c("ANTAGONIST","ANTISENSE INHIBITOR","BLOCKER","INHIBITOR","NEGATIVE ALLOSTERIC MODULATOR"))]<-"inhibitor"

table(otar$mechanism_action_hm)
# activator inhibitor     other 
#      2154     28203      1417  

dim(otar)
# [1] 31774    13

vec<-c("drugId","targetId","adverse_event","adverse_event_colitis","gene_name","actionType","targetType","mechanism_action_hm")
tmp<-otar[which(!otar$diseaseId %in% map_disease_ids$id),vec]
dim(tmp)
# [1] 30660    10

tmp<-tmp[!duplicated(tmp),]
dim(tmp)
# [1] 331   10

tmp<-tmp[which(!tmp$drugId %in% ibd$drugId),]
dim(tmp)
# [1] 291   10

dim(tmp[!duplicated(tmp$gene_name),])
# [1] 44  20 - uniq genes

tmp$gene_name[!duplicated(tmp$gene_name)]
#  [1] "PTGS2"    "JAK2"     "KCNH2"    "GUCY1A1"  "GUCY1B1"  "IL2RB"   
#  [7] "IL2RA"    "IFNGR2"   "ITGAL"    "ITGAV"    "SELL"     "BRD3"    
# [13] "PRKAB1"   "NDUFA4L2" "CCL2"     "NRP1"     "TYK2"     "KIF11"   
# [19] "CXCR1"    "PIM3"     "CCR2"     "ADAM17"   "LNPEP"    "ITGA4"   
# [25] "PTGER2"   "NOD2"     "TRAF3IP2" "PDCD1"    "CD28"     "LTBR"    
# [31] "ITGB8"    "PRKCB"    "PTK2B"    "EZH2"     "PDGFB"    "SLC9A3"  
# [37] "NFKB1"    "IL6ST"    "IL10RA"   "PPP5C"    "TSPO"     "STAT3"   
# [43] "TXK"      "TEC" 

dim(tmp[!duplicated(tmp$drugId),])
# [1] 242  20 - uniq drugs


tmp$pair<-paste(tmp$drugId,tmp$targetId,sep="_")
dim(tmp[!duplicated(tmp$pair),])
# [1] 269  10

dim(tmp)
# [1] 291  10

# add directio of effect from our analyses:

tmp$direction_ibd<-NA
tmp$direction_ibd[which(tmp$gene_name %in% list_genes$list_genes[which(list_genes$direction_effect_colocalization=="same" | list_genes$direction_effect_gsmr=="same")])]<-"same"
tmp$direction_ibd[which(tmp$gene_name %in% list_genes$list_genes[which(list_genes$direction_effect_colocalization=="opposite" | list_genes$direction_effect_gsmr=="opposite")])]<-"opposite"
table(tmp$direction_ibd,useNA="ifany")
# opposite     same     <NA> 
#      182      164       52 

tmp$drug_ibd_direction_concordance<-NA
tmp$drug_ibd_direction_concordance[which(tmp$mechanism_action_hm=="activator" & tmp$direction_ibd=="opposite")]<-1
tmp$drug_ibd_direction_concordance[which(tmp$mechanism_action_hm=="activator" & tmp$direction_ibd=="same")]<-0
tmp$drug_ibd_direction_concordance[which(tmp$mechanism_action_hm=="inhibitor" & tmp$direction_ibd=="opposite")]<-0
tmp$drug_ibd_direction_concordance[which(tmp$mechanism_action_hm=="inhibitor" & tmp$direction_ibd=="same")]<-1
table(tmp$drug_ibd_direction_concordance,useNA="ifany")
#    0    1 <NA> 
#  150  183   65 


# number of drugs targetting the gene/protein in the right direction:
drugs_right_direction<-tmp[which(tmp$drug_ibd_direction_concordance==1),c("drugId","gene_name","drug_ibd_direction_concordance")]
drugs_right_direction<-drugs_right_direction[!duplicated(drugs_right_direction),]
drugs_right_direction
#              drugId gene_name drug_ibd_direction_concordance
# 286   CHEMBL1078178      JAK2                              1
# 4807  CHEMBL1200689   GUCY1A1                              1
# 4808  CHEMBL1200689   GUCY1B1                              1
# 6363  CHEMBL1201550     IL2RB                              1
# 6367  CHEMBL1201550     IL2RA                              1
# 6658  CHEMBL1201605     IL2RA                              1
# 6661  CHEMBL1201605     IL2RB                              1
# 6950  CHEMBL1215923      SELL                              1
# 6982  CHEMBL1231124      JAK2                              1
# 6995  CHEMBL1236936   GUCY1B1                              1
# 6996  CHEMBL1236936   GUCY1A1                              1
# 7285  CHEMBL1287853      JAK2                              1
# 7379     CHEMBL1311   GUCY1B1                              1
# 7383     CHEMBL1311   GUCY1A1                              1
# 7522   CHEMBL136478   GUCY1A1                              1
# 7524   CHEMBL136478   GUCY1B1                              1
# 12640 CHEMBL1743089      NRP1                              1
# 12651 CHEMBL1789941      JAK2                              1
# 15393 CHEMBL1795071      JAK2                              1
# 15427 CHEMBL1829433     KIF11                              1
# 15504 CHEMBL1944698      JAK2                              1
# 15543 CHEMBL2029422      CCR2                              1
# 15546 CHEMBL2035187      JAK2                              1
# 15734 CHEMBL2097081   GUCY1A1                              1
# 15735 CHEMBL2097081   GUCY1B1                              1
# 15751 CHEMBL2103743      JAK2                              1
# 15782 CHEMBL2103847     LNPEP                              1
# 15812 CHEMBL2104967     ITGA4                              1
# 15822 CHEMBL2105661     KIF11                              1
# 15829 CHEMBL2105759      JAK2                              1
# 16216 CHEMBL2107729   GUCY1A1                              1
# 16217 CHEMBL2107729   GUCY1B1                              1
# 16232 CHEMBL2107823      JAK2                              1
# 16236 CHEMBL2107834   GUCY1B1                              1
# 16241 CHEMBL2107834   GUCY1A1                              1
# 16353 CHEMBL2108370     IL2RA                              1
# 16355 CHEMBL2108680     PDCD1                              1
# 16366 CHEMBL2108738     PDCD1                              1
# 17592 CHEMBL2109250      CD28                              1
# 17596 CHEMBL2109326     IL2RB                              1
# 17600 CHEMBL2109367     IL2RA                              1
# 17615 CHEMBL2109509      LTBR                              1
# 17616 CHEMBL2109610     IL2RB                              1
# 17620 CHEMBL2109621     ITGB8                              1
# 17715 CHEMBL2110727      CCR2                              1
# 17797 CHEMBL2178573      CCR2                              1
# 17802 CHEMBL2204263      CCR2                              1
# 18130  CHEMBL228814     KIF11                              1
# 18144 CHEMBL2347655     KIF11                              1
# 18925  CHEMBL300138     PRKCB                              1
# 19033 CHEMBL3137343     PDCD1                              1
# 21083 CHEMBL3287735      EZH2                              1
# 21843 CHEMBL3305901      CCR2                              1
# 21846 CHEMBL3414621      EZH2                              1
# 21943 CHEMBL3545215      JAK2                              1
# 21945 CHEMBL3545221      CCR2                              1
# 21952 CHEMBL3545241      JAK2                              1
# 21954 CHEMBL3545328      JAK2                              1
# 22630 CHEMBL3707228      CCR2                              1
# 22641 CHEMBL3707446      SELL                              1
# 22751 CHEMBL4066936   GUCY1A1                              1
# 22752 CHEMBL4066936   GUCY1B1                              1
# 23166 CHEMBL4297214   GUCY1B1                              1
# 23167 CHEMBL4297214   GUCY1A1                              1
# 23172 CHEMBL4297216      JAK2                              1
# 23180 CHEMBL4297300     PPP5C                              1
# 23201 CHEMBL4297363     ITGA4                              1
# 23258 CHEMBL4297507      JAK2                              1
# 23302 CHEMBL4297571     PDCD1                              1
# 23354 CHEMBL4297616   GUCY1A1                              1
# 23355 CHEMBL4297616   GUCY1B1                              1
# 23382 CHEMBL4297715     PDCD1                              1
# 23758 CHEMBL4297723     PDCD1                              1
# 23878 CHEMBL4297739     STAT3                              1
# 23883 CHEMBL4297829     PDCD1                              1
# 24228 CHEMBL4297831     PDCD1                              1
# 24260 CHEMBL4297840     PDCD1                              1
# 24632 CHEMBL4297843     PDCD1                              1
# 24973 CHEMBL4297856     PDCD1                              1
# 25003 CHEMBL4297858     PDCD1                              1
# 25027 CHEMBL4297877     PDCD1                              1
# 25069 CHEMBL4298037     PDCD1                              1
# 25116 CHEMBL4298124     PDCD1                              1
# 25173 CHEMBL4298167      JAK2                              1
# 25205 CHEMBL4298191     PDCD1                              1
# 25247 CHEMBL4303389      JAK2                              1
# 25334 CHEMBL4457723      CCR2                              1
# 25343 CHEMBL4594260      EZH2                              1
# 25348 CHEMBL4594275      JAK2                              1
# 25396 CHEMBL4594355     PDCD1                              1
# 25400 CHEMBL4594382      JAK2                              1
# 25416 CHEMBL4594419      CCR2                              1
# 25424 CHEMBL4594456     PDCD1                              1
# 25740 CHEMBL4594536     PDCD1                              1
# 25747 CHEMBL4594550     PDCD1                              1
# 25822 CHEMBL4594578     PDCD1                              1
# 25826 CHEMBL4594602     PDCD1                              1
# 25833 CHEMBL4650322   GUCY1A1                              1
# 25834 CHEMBL4650322   GUCY1B1                              1
# 25864 CHEMBL4650472     PDCD1                              1
# 25876 CHEMBL4650482     PDCD1                              1
# 25881 CHEMBL4650516     PDCD1                              1
# 25886  CHEMBL466659   GUCY1B1                              1
# 25889  CHEMBL466659   GUCY1A1                              1
# 26016 CHEMBL4802163      JAK2                              1
# 26030  CHEMBL494089     PRKCB                              1
# 26098  CHEMBL495727      JAK2                              1
# 26138 CHEMBL5095029      CCR2                              1
# 26144 CHEMBL5095049      JAK2                              1
# 29465  CHEMBL574737     PRKCB                              1
# 30400  CHEMBL603469      JAK2                              1
# 30467 CHEMBL6068336      JAK2                              1
# 30470  CHEMBL608533     PRKCB                              1
# 30806    CHEMBL6622   GUCY1A1                              1
# 30808    CHEMBL6622   GUCY1B1                              1
# 31342     CHEMBL730   GUCY1B1                              1
# 31343     CHEMBL730   GUCY1A1                              1
# 31749   CHEMBL91829     PRKCB                              1

# NUMBER DRUGS:
dim(drugs_right_direction[!duplicated(drugs_right_direction$drugId),])
# [1] 102   4

# NUMBER GENES:
dim(drugs_right_direction[!duplicated(drugs_right_direction$gene_name),])
# [1] 19  3


drugs_right_direction[!duplicated(drugs_right_direction$gene_name),"gene_name"]
#  [1] "JAK2"    "GUCY1A1" "GUCY1B1" "IL2RB"   "IL2RA"   "SELL"    "NRP1"   
#  [8] "KIF11"   "CCR2"    "LNPEP"   "ITGA4"   "PDCD1"   "CD28"    "LTBR"   
# [15] "ITGB8"   "PRKCB"   "EZH2"    "PPP5C"   "STAT3"  

# "CD40","ERAP2","INPP5D","FGFR4","MUC1","IL10RB","NQO1","TNNI2" -lost


dim(tmp)
# [1] 291  11

# save this table:
fwrite(tmp, paste0(path_tables, "list_drugs_genes_pairs_noibd_clinical_trials.tsv.gz"),
       col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#####################################################################################################################################


ids<-drugs_right_direction[!duplicated(drugs_right_direction$gene_name),"gene_name"]

rep<-matrix(ncol=3,nrow=length(ids))
colnames(rep)<-c("approvedSymbol","drugIds","mechanism_action")
rep<-as.data.frame(rep)

rep$gene_name<-ids


for (i in 1:nrow(rep)) {
    
    drugIDs<-drugs_right_direction[which(drugs_right_direction$gene_name==rep$gene_name[i]),"drugId"]
    drugIDs<-drugIDs[!duplicated(drugIDs)]

    rep$drugIds[i]<-paste(drugIDs,collapse="|")

    tmp_mechanism<-tmp$mechanism_action[which(tmp$drugId %in% drugIDs)]
    tmp_mechanism<-tmp_mechanism[!duplicated(tmp_mechanism)]

    rep$mechanism_action[i]<-paste(tmp_mechanism,collapse="|")

}



#####################################################################################################################################

####### adverse events:
tmp1<-tmp[,c("drugId","maxClinicalStage","mechanism_action_hm","drug_ibd_direction_concordance","adverse_event","adverse_event_colitis")]
tmp1<-tmp1[!duplicated(tmp1),]
dim(tmp1[duplicated(tmp1$drugId),])
# [1] 461     3

table(tmp1$adverse_event_colitis,tmp1$drug_ibd_direction_concordance)
#       0   1
#   0 307 257
#   1   7  37

tmp1[which(tmp1$adverse_event_colitis==1),]

table(otar$gene_name[which(otar$drugId %in% tmp1$drugId[which(tmp1$adverse_event_colitis==1)])])
# IFNGR2  PDCD1  PRKCB 
#     69   3998    275 






