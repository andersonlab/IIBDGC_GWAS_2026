# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# run enrichment analysis with list of genes from coloc + coding genes

# download lastest release of pathway gene set database dated 1 September 2021

# cd /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/pathways/
# mkdir November_03_2025

# last version:
# wget http://download.baderlab.org/EM_Genesets/November_03_2025/Human/Summary_GeneSet_Counts_symbol.txt
# wget http://download.baderlab.org/EM_Genesets/November_03_2025/Human/symbol/Human_GO_AllPathways_noPFOCR_no_GO_iea_November_03_2025_symbol.gmt 
# wget http://download.baderlab.org/EM_Genesets/November_03_2025/Human/symbol/Human_AllPathways_noPFOCR_November_03_2025_symbol.gmt 
# wget http://download.baderlab.org/EM_Genesets/November_03_2025/Human/symbol/symbol_translation_summary.log


# mkdir March_02_2026
# wget http://download.baderlab.org/EM_Genesets/March_02_2026/Human/Summary_GeneSet_Counts_symbol.txt
# wget http://download.baderlab.org/EM_Genesets/March_02_2026/Human/symbol/Human_GO_AllPathways_noPFOCR_no_GO_iea_March_02_2026_symbol.gmt 
# wget http://download.baderlab.org/EM_Genesets/March_02_2026/Human/symbol/Human_AllPathways_noPFOCR_March_02_2026_symbol.gmt 
# wget http://download.baderlab.org/EM_Genesets/March_02_2026/Human/symbol/symbol_translation_summary.log




# list of genes from WES (by 22/01/2025)
# list of annotations by James (by cell type and direction of effect)
# /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/pathways/directional_celltypes_MegaGWAS_coloc.txt
# /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/pathways/directional_celltypes_MegaGWAS_coloc_20260111.txt


######################################################################################################################################
# 1.- Get list with all genes identified by colocalization analysis:

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=4000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R


library(stringr)
library(qusage)
library(ggplot2)
library(reshape2)
library(grid)
# library(ggdendro)
library(ggpubr)
library(data.table)

rm(list=ls())
path<-"/path/to/ibdgwas/IIBDGC/"

library(qvalue)

# LOAD THE INDEPENDENT SIGNALS - see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_38_2_define_independent_signals_using_as_ref_cojo_agnostic_results.R
all<-fread(paste(path,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 619 192


# eqtl/coloc genes
all_genes<-fread(paste0(path,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence_with_coloc_mr_direction_with_ligand_receptor_with_protein_complexes_with_monogenic.tsv.gz"))

all_genes<-as.data.frame(all_genes[,1,drop=F])
dim(all_genes)
# [1] 664

all_genes$all_genes<-1

# list of genes by group as defined by James:
list_genes_updown<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/pathways/directional_celltypes_MegaGWAS_coloc_20260111.txt")
list_genes_updown<-as.data.frame(list_genes_updown)

colnames(list_genes_updown)

for (i in 1:length(colnames(list_genes_updown))) {

  tmp<-list_genes_updown[,i]
  tmp<-tmp[which(!tmp %in% c(""))]

  all_genes[,colnames(list_genes_updown)[i]]<-0
  all_genes[which(all_genes$list_genes %in% tmp),colnames(list_genes_updown)[i]]<-1

  rm(tmp)

}
# rm(list_genes_updown)


etstmp<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/pathways/ets2_genesetsSYMBOL_MASTER_forLF.csv")

ets2_all<-c(etstmp$ETS2_g1_g4_DN,etstmp$OE_ETS2_UP_combined)
ets2_all<-ets2_all[!duplicated(ets2_all)]
ets2_all<-ets2_all[!ets2_all %in% c("")]
length(ets2_all)
# [1] 1732

ets2_dn<-c(etstmp$ETS2_g1_g4_DN)
ets2_dn<-ets2_dn[!duplicated(ets2_dn)]
ets2_dn<-ets2_dn[!ets2_dn %in% c("")]
length(ets2_dn)
# [1] 672

ets2_up<-c(etstmp$OE_ETS2_UP_combined)
ets2_up<-ets2_up[!duplicated(ets2_up)]
ets2_up<-ets2_up[!ets2_up %in% c("")]
length(ets2_up)
# [1] 1217

ets2_up[which(ets2_up %in% ets2_dn)]
#   [1] "ACP3"        "ADA2"        "ADAM8"       "ADCY7"       "ADGRE2"     
#   [6] "ALOX5AP"     "ANKH"        "AP1B1"       "APLP2"       "ARHGAP26"   
#  [11] "ARHGAP4"     "ARHGAP9"     "ARHGDIB"     "ATP8B4"      "BST1"       
#  [16] "C1orf162"    "CCL3"        "CD163"       "CD209"       "CD33"       
#  [21] "CD36"        "CD4"         "CD93"        "CD99"        "CEMIP2"     
#  [26] "CHI3L2"      "CHSY1"       "CISH"        "COL23A1"     "COLEC12"    
#  [31] "CORO1A"      "CR1"         "CTNND1"      "CTSS"        "CYBB"       
#  [36] "DCUN1D3"     "DEF6"        "DENND4C"     "DHRS7"       "DNAJC10"    
#  [41] "DOCK2"       "DPYD"        "ECM1"        "EEF2K"       "EMILIN2"    
#  [46] "EMP1"        "EPAS1"       "ERCC1"       "ETS2"        "FCGR2A"     
#  [51] "FCGRT"       "FHL3"        "FMN1"        "FOXRED2"     "FUCA1"      
#  [56] "GASK1B"      "GBA"         "GBAP1"       "GBP5"        "GCH1"       
#  [61] "GGTA1"       "GMFG"        "GNA11"       "GNS"         "GPR84"      
#  [66] "GPX3"        "GRAMD4"      "HEG1"        "HEXA"        "HEXB"       
#  [71] "HFE"         "HK1"         "HLA-DMB"     "HM13"        "IL1B"       
#  [76] "IL1RAP"      "IRF8"        "ITGAM"       "ITGB2"       "KCNAB2"     
#  [81] "KCNJ2"       "LAPTM5"      "LDLRAD4"     "LILRB2"      "LILRB5"     
#  [86] "LIMK2"       "LIPA"        "LRMDA"       "LY86"        "MACF1"      
#  [91] "MAN2B1"      "MAN2B2"      "MARCO"       "MIR3945HG"   "MMP9"       
#  [96] "MOB3B"       "MPP1"        "MPRIP"       "MPV17"       "MRC1"       
# [101] "MSRA"        "MYO1E"       "NAAA"        "NAGA"        "NPC2"       
# [106] "NPL"         "PARVG"       "PIEZO1"      "PLAUR"       "PLCB2"      
# [111] "PLCG2"       "PLP2"        "PLXDC2"      "PNP"         "PNPO"       
# [116] "PRKCB"       "PTGS1"       "PTPRM"       "PYGL"        "RAC2"       
# [121] "RASGRP4"     "RGL1"        "RHOU"        "RNASET2"     "S100A8"     
# [126] "S100A9"      "SCPEP1"      "SEMA3A"      "SERPINA1"    "SERPINB9"   
# [131] "SGMS1"       "SIGLEC15"    "SIGLEC7"     "SIRT2"       "SLC12A5-AS1"
# [136] "SLC2A9"      "SLC39A11"    "SLCO2B1"     "SPARC"       "SUCNR1"     
# [141] "SULT1A1"     "SVIL"        "TCF4"        "TEP1"        "TIMP2"      
# [146] "TLR4"        "TM7SF3"      "TM9SF2"      "TMC8"        "TMEM106A"   
# [151] "TPP1"        "TRAF3IP3"    "TRAPPC9"     "TREM2"       "TRPV4"      
# [156] "VMP1"        "ZC3H12A"    

ets2_dn[which(ets2_dn %in% ets2_up)]
#   [1] "ACP3"        "ADA2"        "ADAM8"       "ADCY7"       "ADGRE2"     
#   [6] "ALOX5AP"     "ANKH"        "AP1B1"       "APLP2"       "ARHGAP26"   
#  [11] "ARHGAP4"     "ARHGAP9"     "ARHGDIB"     "ATP8B4"      "BST1"       
#  [16] "C1orf162"    "CCL3"        "CD163"       "CD209"       "CD33"       
#  [21] "CD36"        "CD4"         "CD93"        "CD99"        "CEMIP2"     
#  [26] "CHI3L2"      "CHSY1"       "CISH"        "COL23A1"     "COLEC12"    
#  [31] "CORO1A"      "CR1"         "CTNND1"      "CTSS"        "CYBB"       
#  [36] "DCUN1D3"     "DEF6"        "DENND4C"     "DHRS7"       "DNAJC10"    
#  [41] "DOCK2"       "DPYD"        "ECM1"        "EEF2K"       "EMILIN2"    
#  [46] "EMP1"        "EPAS1"       "ERCC1"       "ETS2"        "FCGR2A"     
#  [51] "FCGRT"       "FHL3"        "FMN1"        "FOXRED2"     "FUCA1"      
#  [56] "GASK1B"      "GBA"         "GBAP1"       "GBP5"        "GCH1"       
#  [61] "GGTA1"       "GMFG"        "GNA11"       "GNS"         "GPR84"      
#  [66] "GPX3"        "GRAMD4"      "HEG1"        "HEXA"        "HEXB"       
#  [71] "HFE"         "HK1"         "HLA-DMB"     "HM13"        "IL1B"       
#  [76] "IL1RAP"      "IRF8"        "ITGAM"       "ITGB2"       "KCNAB2"     
#  [81] "KCNJ2"       "LAPTM5"      "LDLRAD4"     "LILRB2"      "LILRB5"     
#  [86] "LIMK2"       "LIPA"        "LRMDA"       "LY86"        "MACF1"      
#  [91] "MAN2B1"      "MAN2B2"      "MARCO"       "MIR3945HG"   "MMP9"       
#  [96] "MOB3B"       "MPP1"        "MPRIP"       "MPV17"       "MRC1"       
# [101] "MSRA"        "MYO1E"       "NAAA"        "NAGA"        "NPC2"       
# [106] "NPL"         "PARVG"       "PIEZO1"      "PLAUR"       "PLCB2"      
# [111] "PLCG2"       "PLP2"        "PLXDC2"      "PNP"         "PNPO"       
# [116] "PRKCB"       "PTGS1"       "PTPRM"       "PYGL"        "RAC2"       
# [121] "RASGRP4"     "RGL1"        "RHOU"        "RNASET2"     "S100A8"     
# [126] "S100A9"      "SCPEP1"      "SEMA3A"      "SERPINA1"    "SERPINB9"   
# [131] "SGMS1"       "SIGLEC15"    "SIGLEC7"     "SIRT2"       "SLC12A5-AS1"
# [136] "SLC2A9"      "SLC39A11"    "SLCO2B1"     "SPARC"       "SUCNR1"     
# [141] "SULT1A1"     "SVIL"        "TCF4"        "TEP1"        "TIMP2"      
# [146] "TLR4"        "TM7SF3"      "TM9SF2"      "TMC8"        "TMEM106A"   
# [151] "TPP1"        "TRAF3IP3"    "TRAPPC9"     "TREM2"       "TRPV4"      
# [156] "VMP1"        "ZC3H12A"   

ets2 <- list(ETS2 = ets2_all,
        ETS2_DN = ets2_dn,
        ETS2_UP = ets2_up)


#####################
# A.- PATHWAYS ONLY #
#####################

p<-read.gmt(paste(path,"post_imputation/2022/analysis/metaanalysis/annotation/pathways/March_02_2026/Human_AllPathways_noPFOCR_March_02_2026_symbol.gmt",sep=""))
#list
length(names(p))
#[1] 24411

p<-append(p,ets2)
length(names(p))
#[1] 24411
 

tmp<-data.frame(pathway=names(p), n_genes_pathway=sapply(p, length))

tmp<-tmp[!grepl("PATHWHIZ|NETPATH|SMPDB|ECOCYC|BIOCYC|ARACYC|HUMANCYC",tmp$pathway),]
dim(tmp)
# [1] 18895     2

tmp<-tmp[!duplicated(tmp),]
dim(tmp)
# [1] 4726    2
summary(tmp$n_genes_pathway)

# tmp<-tmp[which( (tmp$n_genes_pathway<=600 & tmp$n_genes_pathway>=2) | (tmp$pathway %in% c("ETS2","ETS2_DN","ETS2_UP"))),]

p<-p[!duplicated(p)]

p<-p[which(names(p) %in% tmp$pathway)]
length(names(p))
#[1] 4138

rm(tmp)

##########################################################
# A.1.- generate the list of N variants, background level:

p_genes<-unique(unlist(p))
length(p_genes)
#[1] 14253

# number of genes in the union of all annotation sets of the selected typ an identifiable in the user-provided identifier (symbol) namespace (background)
Np<-length(p_genes)

################
# B.- GO TERMS #
################

po<-read.gmt(paste(path,"post_imputation/2022/analysis/metaanalysis/annotation/pathways/March_02_2026/Human_GO_AllPathways_noPFOCR_no_GO_iea_March_02_2026_symbol.gmt",sep=""))

#list
length(names(po))
# [1] 45545

#keep only GO terms, no pathways
po<-po[grep(paste(c("GOMF","GOBP","GOCC"),collapse="|"), names(po), value = TRUE)]
length(names(po))
# [1] 21137

#keep only GO terms with less than 600 genes:
tmp<-data.frame(pathway=names(po), n_genes_pathway=sapply(po, length))
summary(tmp$n_genes_pathway)

dim(tmp)
# [1] 21137     2
tmp<-tmp[!duplicated(tmp),]
dim(tmp)
# [1] 21137     2

# tmp<-tmp[which( (tmp$n_genes_pathway<=600 & tmp$n_genes_pathway>=2) | (tmp$pathway %in% c("ETS2","ETS2_DN","ETS2_UP"))),]

po<-po[!duplicated(po)]
po<-po[which(names(po) %in% tmp$pathway)]
length(names(po))
#[1] 16694

rm(tmp)

##########################################################
# B.1.- generate the list of N variants, background level:

po_genes<-unique(unlist(po))
length(po_genes)
#[1] 17015

# number of genes in the union of all annotation sets of the selected typ an identifiable in the user-provided identifier (symbol) namespace (background)
No<-length(po_genes)


#############################################
# create an output per list of gene:

# colnames(all_genes)[1]<-"all_genes"

for (jj in c(2:ncol(all_genes))) {

  print(jj)
  print(colnames(all_genes)[jj])

  # if (jj==1) {
  #   list_genes<-all_genes[,jj,drop=F]
  #   list_genes$tmp<-1
  # } else {
    print(table(all_genes[,jj]))
    list_genes<-all_genes[,c(1,jj)]
    colnames(list_genes)[2]<-"tmp"
    print(table(list_genes$tmp))
  # }

  colnames(list_genes)[1]<-"all_genes"
  



  list_genes<-list_genes[which(list_genes$tmp==1),"all_genes"]
  print(length(list_genes))

  print(list_genes)

  ###############################################################################
  # A.2.- identify pathways that include at least two IBD genes

  # number of genes that the user uploads an that are part of a least one annoation set of the selected type:
  mp<-length(list_genes[which(list_genes %in% p_genes)])
  print(paste("N IBD genes in pathways",mp))
  # [1] "N IBD genes in pathways 456"

  # list genes not in pathways:
  # list_genes[which(!list_genes %in% p_genes)]

  dat_p<-data.frame(
    pathway=names(p),
    n_genes_pathway=sapply(p, length),
    n_ibd_genes_in_pathway=sapply(p, function(x) sum(x %in% list_genes)),
    ibd=sapply(p, function(x) paste(x[x %in% list_genes], collapse=";"))
  )

  dat_p$p_value<-ifelse(
    dat_p$n_ibd_genes_in_pathway >= 2,
    1 - phyper(dat_p$n_ibd_genes_in_pathway - 1, dat_p$n_genes_pathway, Np - dat_p$n_genes_pathway, mp, lower.tail=TRUE, log.p=FALSE),
    1
  )


  # print(paste("N pathways:",nrow(dat_p)))

  # use fdr only for those sets where minimum overlap with input list (n=2 default) is fulfill, as in consensuspathdb. Use BH as in the consensuspathdb
  dat_p<-dat_p[which(dat_p$n_ibd_genes_in_pathway>=2),]

  dat_p<-dat_p[order(dat_p$p_value,decreasing=F),]
  # dat_p$q_value<-p.adjust(dat_p$p_value, method = c("BH"), n = nrow(dat_p))
  q_values<-try(qvalue(dat_p$p_value)$qval, silent = TRUE)

  if ("try-error" %in% class(q_values)) {
      # if qvalue fails this is probably because there are no p-values greater than
      # 0.95 (the highest lambda value
      # if so add a single p-value of 1 to try to combat the problem
      q_values = qvalue::qvalue(c(dat_p$p_value, 1))$qval
      dat_p$q_value = q_values[-length(q_values)]
  } else {
    dat_p$q_value = q_values
  }

  # #keep only significant at fdr 5%
  # dat_p<-dat_p[which(dat_p$q_value<=0.05),]


  ###############################################################################
  # B.2.- identify GO terms that include at least two IBD genes

  # number of genes that the user uploads an that are part of a least one annoation set of the selected type:
  mo<-length(list_genes[which(list_genes %in% po_genes)])
  print(paste("N IBD genes in GO",mo))
  # [1] "N IBD genes in GO 551"

  dat_po<-data.frame(
    pathway=names(po),
    n_genes_pathway=sapply(po, length),
    n_ibd_genes_in_pathway=sapply(po, function(x) sum(x %in% list_genes)),
    ibd=sapply(po, function(x) paste(x[x %in% list_genes], collapse=";"))
  )

  dat_po$p_value<-ifelse(
    dat_po$n_ibd_genes_in_pathway >= 2,
    1 - phyper(dat_po$n_ibd_genes_in_pathway - 1, dat_po$n_genes_pathway, No - dat_po$n_genes_pathway, mo, lower.tail=TRUE, log.p=FALSE),
    1
  )

  # print(paste("N GO:",nrow(dat_po)))

  # use fdr only for those sets where minimum overlap with input list (n=2 default) is fulfill, as in consensuspathdb. Use BH as in the consensuspathdb
  dat_po<-dat_po[which(dat_po$n_ibd_genes_in_pathway>=2),]

  dat_po<-dat_po[order(dat_po$p_value,decreasing=F),]

  q_values<-try(qvalue(dat_po$p_value)$qval, silent = TRUE)

  if ("try-error" %in% class(q_values)) {
      # if qvalue fails this is probably because there are no p-values greater than
      # 0.95 (the highest lambda value
      # if so add a single p-value of 1 to try to combat the problem
      q_values = qvalue::qvalue(c(dat_po$p_value, 1))$qval
      dat_po$q_value = q_values[-length(q_values)]
  } else {
    dat_po$q_value = q_values
  }

  # # plot top results:

  # dat_p$hitsPerc<-100*(dat_p$n_ibd_genes_in_pathway/dat_p$n_genes_pathway)
  # dat_p$path<-gsub("\\%.*","",dat_p$pathway)
  # dat_p$path[which(dat_p$path %in% dat_p$path[duplicated(dat_p$path)])]<-paste(dat_p$path[which(dat_p$path %in% dat_p$path[duplicated(dat_p$path)])],seq(1:length(dat_p$path[which(dat_p$path %in% dat_p$path[duplicated(dat_p$path)])])),sep="_")

  # dat_p_top<-dat_p[1:50,]
  # dat_p_top<-dat_p_top[which(dat_p_top$p_value!=""),]
  # dat_p_top$path<-factor(dat_p_top$path,levels=rev(dat_p_top$path))

  # p1<-ggplot(data =dat_p_top, aes(x=hitsPerc, 
  #                                 y=path, 
  #                                 colour=q_value, 
  #                                 size=n_ibd_genes_in_pathway)) +
  #   geom_point() +
  #   expand_limits(x=0) +
  #   labs(x="IBD Genes in Pathway (%)", y="", colour="q value", size="N IBD Genes")


  # # plot top results:

  # dat_po$hitsPerc<-100*(dat_po$n_ibd_genes_in_pathway/dat_po$n_genes_pathway)
  # dat_po$path<-gsub("\\%.*","",dat_po$pathway)

  # dat_po_top<-dat_po[1:50,]
  # dat_po_top<-dat_po_top[which(dat_po_top$p_value!=""),]
  # dat_po_top$path<-factor(dat_po_top$pathway,levels=rev(dat_po_top$pathway))

  # p2<-ggplot(data =dat_po_top, aes(x=hitsPerc, 
  #                             y=path, 
  #                             colour=q_value, 
  #                             size=n_ibd_genes_in_pathway)) +
  #   geom_point() +
  #   expand_limits(x=0) +
  #   labs(x="IBD Genes in GO term (%)", y="", colour="q value", size="N IBD Genes")

  # pdf(paste0("~/git/IIBDGC_GWAS/plots/pathways/figures/pathways_",colnames(all_genes)[jj],".pdf"),width=20, height=12, onefile=FALSE)
  # print(ggarrange(p1,p2,ncol=2,common.legend = T,legend="bottom"))
  # dev.off()

  if (nrow(dat_po) > 0 && nrow(dat_p) > 0) {

    dat_po$data<-"GeneOntology"
    dat_p$data<-"Pathways"

    all_dat<-rbind(dat_po,dat_p)

  } else if (nrow(dat_po) > 0 && nrow(dat_p) == 0) {

    print(paste(colnames(all_genes)[jj],": No significant pathways"))
    dat_po$data<-"GeneOntology"
    all_dat<-dat_po

  } else if (nrow(dat_po) == 0 && nrow(dat_p) > 0) {

    print(paste(colnames(all_genes)[jj],": No significant GO terms"))
    dat_p$data<-"Pathways"
    all_dat<-dat_p

  } else {

    print(paste(colnames(all_genes)[jj],": No significant pathways or GO terms"))
    rm(dat_po,dat_p)

  }


  if (exists("all_dat")) {

    all_dat$path<-gsub("\\%.*","",all_dat$pathway)

    write.table(all_dat,paste0("~/git/IIBDGC_GWAS/plots/pathways/tables/gene_enrichment_analysis_pathways_go_",colnames(all_genes)[jj],".csv"),
    col.names=T,row.names=F,quote=F,sep=",")

    rm(all_dat,dat_po,dat_p)

  }

}



q("no")




