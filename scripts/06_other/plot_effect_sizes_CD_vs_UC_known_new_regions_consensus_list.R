# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=6000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R


library(rcartocolor)
library(qvalue)
library(data.table)
library(ggplot2)
library(colorspace)
library(dplyr)
library(tibble)
library(ggrepel)
library(ggpubr)
library(tidyverse)
library(scales)
library(png)
library(grid)


rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz",sep=""),)
all<-as.data.frame(all)
dim(all)
# [1] 653 184

table(all$class_signal_final,useNA="ifany")
# new_cojo_supervised_gw_significant_multiancestry_known_signal 
#                                                             5 
#   new_cojo_supervised_gw_significant_multiancestry_new_signal 
#                                                            17 
#                            new_cojo_unsupervised_known_signal 
#                                                           234 
#                              new_cojo_unsupervised_new_signal 
#                                                           399 

# get the average frequency in EUR tier2 in controls:

chr<-as.numeric(unlist(labels(table(all$chr))))
chr<-chr[which(!is.na(chr))]

pheno<-c("ibd","cd","uc")

for (i in 1:length(pheno)) {
  
  print(pheno[i])
  
  for (j in 1:length(chr)) {
    
    tmp1<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr[j],"_",pheno[i],"_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz",sep=""),head=T)
    tmp1<-tmp1[which(tmp1$MarkerName %in% all$MarkerName),]

    tmp1<-tmp1[,c("MarkerName","avgA2FREQ","avgA2FREQ_CASES","avgA2FREQ_CONTROLS")]
    
    if(j==1) {
      tmp<-tmp1
    } else {
      tmp<-rbind(tmp,tmp1)
    }
    
    rm(tmp1)
    
  }

  colnames(tmp)[2:4]<-paste(colnames(tmp)[2:4],pheno[i],"eur_tier_2",sep="_")
  
  if(i==1) {
    dat<-tmp
  } else {
    dat<-merge(dat,tmp,by="MarkerName",sort=F)
  }
  rm(tmp)
}



all<-merge(all,dat,by="MarkerName",all.x=T)

all$MAF<-pmin(all$avgA2FREQ_CONTROLS_ibd_eur_tier_2,1-all$avgA2FREQ_CONTROLS_ibd_eur_tier_2)
all$MAF[which(is.na(all$MAF))]<-pmin(all$avgA2FREQ_CONTROLS_ibd_eur_tier_2[which(is.na(all$MAF))],1-all$avgA2FREQ_CONTROLS_ibd_eur_tier_2[which(is.na(all$MAF))])


table(all$class_signal_final_exome,useNA="ifany")
#                  new_cojo_supervised_gw_significant_multiancestry_known_signal 
#                                                                              5 
#                    new_cojo_supervised_gw_significant_multiancestry_new_signal  # new
#                                                                             16 
# new_cojo_supervised_gw_significant_multiancestry_new_signal_new_exonic_variant  # new
#                                                                              1 
#                                             new_cojo_unsupervised_known_signal 
#                                                                            232 
#                        new_cojo_unsupervised_known_signal_known_exonic_variant 
#                                                                              1 
#                          new_cojo_unsupervised_known_signal_new_exonic_variant
#                                                                              1 
#                                               new_cojo_unsupervised_new_signal # new
#                                                                            380 
#                          new_cojo_unsupervised_new_signal_known_exonic_variant 
#                                                                              5 
#                            new_cojo_unsupervised_new_signal_new_exonic_variant # new
#                                                                             14 

all$class_region<-NA
all$class_region[which(all$class_signal_final_exome %in% c("new_cojo_supervised_gw_significant_multiancestry_known_signal",
"new_cojo_unsupervised_known_signal",
"new_cojo_unsupervised_known_signal_known_exonic_variant",
"new_cojo_unsupervised_known_signal_new_exonic_variant",
"new_cojo_unsupervised_new_signal_known_exonic_variant")) ]<-"Known signals"

all$class_region[which(all$class_signal_final_exome %in% c("new_cojo_supervised_gw_significant_multiancestry_new_signal",
"new_cojo_supervised_gw_significant_multiancestry_new_signal_new_exonic_variant",
"new_cojo_unsupervised_new_signal",
"new_cojo_unsupervised_new_signal_new_exonic_variant"))]<-"New signals"

table(all$class_region,useNA="ifany")
# Known signals   New signals 
#           244           411 


table(all$phenotype,useNA="ifany")
            #  CD   IBD_saturated IBD_unsaturated              UC 
            # 148             114             289             104 



all$Effect<-NA
all$StdErr<-NA
all$pvalue<-NA


all$Effect[which(all$phenotype=="IBD_unsaturated")]<-all$BETA_ibd_eur_tier_2[which(all$phenotype=="IBD_unsaturated")]
all$StdErr[which(all$phenotype=="IBD_unsaturated")]<-all$SE_ibd_eur_tier_2[which(all$phenotype=="IBD_unsaturated")]
all$pvalue[which(all$phenotype=="IBD_unsaturated")]<-all$"P-value_ibd_eur_tier_2"[which(all$phenotype=="IBD_unsaturated")]

all$Effect[which(all$phenotype=="CD")]<-all$BETA_cd_eur_tier_2[which(all$phenotype=="CD")]
all$StdErr[which(all$phenotype=="CD")]<-all$SE_cd_eur_tier_2[which(all$phenotype=="CD")]
all$pvalue[which(all$phenotype=="CD")]<-all$"P-value_ibd_eur_tier_2"[which(all$phenotype=="CD")]

all$Effect[which(all$phenotype=="UC")]<-all$BETA_uc_eur_tier_2[which(all$phenotype=="UC")]
all$StdErr[which(all$phenotype=="UC")]<-all$SE_uc_eur_tier_2[which(all$phenotype=="UC")]
all$pvalue[which(all$phenotype=="UC")]<-all$"P-value_uc_eur_tier_2"[which(all$phenotype=="UC")]


test<-all[which(all$phenotype=="IBD_saturated"),c("MarkerName","P-value_ibd_eur_tier_2","P-value_cd_eur_tier_2","P-value_uc_eur_tier_2")]
test$model<-NA
test$model[which(test$"P-value_cd_eur_tier_2"<test$"P-value_uc_eur_tier_2")]<-"CD"
test$model[which(test$"P-value_cd_eur_tier_2">test$"P-value_uc_eur_tier_2")]<-"UC"
table(test$model)
# CD UC 
# 87 42  
cd_model<-test$MarkerName[which(test$model=="CD")]
uc_model<-test$MarkerName[which(test$model=="UC")]


# list of IBD-sat variants more significant for CD
all$Effect[which(all$MarkerName %in% cd_model)]<-all$BETA_cd_eur_tier_2[which(all$MarkerName %in% cd_model)]
all$StdErr[which(all$MarkerName %in% cd_model)]<-all$SE_cd_eur_tier_2[which(all$MarkerName %in% cd_model)]
all$pvalue[which(all$MarkerName %in% cd_model)]<-all$"P-value_cd_eur_tier_2"[which(all$MarkerName %in% cd_model)]

# list of IBD-sat variants more significant for UC
all$Effect[which(all$MarkerName %in% uc_model)]<-all$BETA_uc_eur_tier_2[which(all$MarkerName %in% uc_model)]
all$StdErr[which(all$MarkerName %in% uc_model)]<-all$SE_uc_eur_tier_2[which(all$MarkerName %in% uc_model)]
all$pvalue[which(all$MarkerName %in% uc_model)]<-all$"P-value_uc_eur_tier_2"[which(all$MarkerName %in% uc_model)]


table(all$class_region)
# Known signals   New signals 
#           244           411 


table(all$class_region,all$phenotype,useNA="ifany")
#                CD IBD_saturated IBD_unsaturated  UC
#   Known signals  51            71              95  27
#   New signals    97            43             194  77

x<-table(all$class_region,all$phenotype,useNA="ifany")

all$pheno_class<-"NA"
all$pheno_class<-as.character(all$pheno_class)
all$pheno_class[which(all$class_region=="New signals" & all$phenotype=="IBD_saturated")]<-paste("New IBD saturated signal (N=",x[2,2],")",sep="")
all$pheno_class[which(all$class_region=="New signals" & all$phenotype=="IBD_unsaturated")]<-paste("New IBD unsaturated signal (N=",x[2,3],")",sep="")
all$pheno_class[which(all$class_region=="New signals" & all$phenotype=="CD")]<-paste("New CD signal (N=",x[2,1],")",sep="")
all$pheno_class[which(all$class_region=="New signals" & all$phenotype=="UC")]<-paste("New UC signal (N=",x[2,4],")",sep="")
all$pheno_class[which(all$class_region=="Known signals" & all$phenotype=="IBD_saturated")]<-paste("Known IBD saturated signal (N=",x[1,2],")",sep="")
all$pheno_class[which(all$class_region=="Known signals" & all$phenotype=="IBD_unsaturated")]<-paste("Known IBD unsaturated signal (N=",x[1,3],")",sep="")
all$pheno_class[which(all$class_region=="Known signals" & all$phenotype=="CD")]<-paste("Known CD signal (N=",x[1,1],")",sep="")
all$pheno_class[which(all$class_region=="Known signals" & all$phenotype=="UC")]<-paste("Known UC signal (N=",x[1,4],")",sep="")

table(all$pheno_class,all$phenotype,useNA="ifany")


all$pheno_class<-factor(all$pheno_class,levels = (c(paste("New CD signal (N=",x[2,1],")",sep=""),paste("Known CD signal (N=",x[1,1],")",sep=""),
                                                    paste("New UC signal (N=",x[2,4],")",sep=""),paste("Known UC signal (N=",x[1,4],")",sep="")
                                                    ,paste("New IBD saturated signal (N=",x[2,2],")",sep=""),paste("Known IBD saturated signal (N=",x[1,2],")",sep=""),
                                                    paste("New IBD unsaturated signal (N=",x[2,3],")",sep=""),paste("Known IBD unsaturated signal (N=",x[1,3],")",sep=""))))
table(all$pheno_class)
#                New CD signal (N=97)              Known CD signal (N=51) 
#                                  97                                  51 
#                New UC signal (N=77)              Known UC signal (N=27) 
#                                  77                                  27 
#     New IBD saturated signal (N=43)   Known IBD saturated signal (N=71) 
#                                  43                                  71 
#  New IBD unsaturated signal (N=194) Known IBD unsaturated signal (N=95) 
#                                 194                                  95

cols<-c('#004488','#BB5566','#DDAA33',"#db7107")
col_ligth<-lighten(cols, 0.70, space = "HCL")
# col_dark<-darken(cols,0.2,space = "combined")
col_all<-c(cols[1],col_ligth[1],cols[2],col_ligth[2],cols[3],col_ligth[3],cols[4],col_ligth[4])

all$OR<-exp((all$Effect))
all$confint_low<-exp((all$Effect)-(1.96*all$StdErr))
all$confint_up<-exp((all$Effect)+(1.96*all$StdErr))

all<-all[order(all$MAF,decreasing=T),]
# all<-all[order(as.character(all$pheno_class),decreasing=F),]

write.table(all,"")
range(all$OR,na.rm=T)
# [1] 0.3998363 2.2443142
range(all$confint_low,na.rm=T)
# [1] 0.3511127 2.1041593
range(all$confint_up,na.rm=T)
# [1] [1] 0.4553213 2.9051586


all$gene<-NA
all$gene[which(all$class_signal_final_exome %in% c("new_cojo_unsupervised_new_signal_known_exonic_variant"))]<-all$exonic_variant_in_ld_gene_aac[which(all$class_signal_final_exome %in% c("new_cojo_unsupervised_new_signal_known_exonic_variant"))]


p1<-ggplot(all, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
  geom_errorbar(aes(ymin=confint_low, ymax=confint_up)) +
  geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
  scale_color_manual(values=col_all) +
  scale_shape_manual(values=c(19,20,19,20,19,20,19,20),guide=F)  +
  theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
  geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=all %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
  scale_x_continuous(trans='reverse',limits=c(-log10(0.001),-log10(0.6)),
  breaks = -log10(as.numeric(c("0.001","0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
  label = c("0.001","0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
  scale_y_continuous(limits=c(0.2,3),breaks = c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3))

pdf("~/git/IIBDGC_GWAS/plots/paper_figures/Figure_1_signals_effect_phenotype_KNOWN_EXONIC.pdf",width=10, height=5)
p1
dev.off()

all$gene<-NA

all$gene[which(all$class_signal_final_exome %in% c("new_cojo_unsupervised_new_signal_new_exonic_variant","new_cojo_supervised_gw_significant_multiancestry_new_signal_new_exonic_variant"))]<-all$exonic_variant_in_ld_gene_aac[which(all$class_signal_final_exome %in% c("new_cojo_unsupervised_new_signal_new_exonic_variant","new_cojo_supervised_gw_significant_multiancestry_new_signal_new_exonic_variant"))]
all$gene[which(all$gene %in% c("MFNG:MFNG:Arg302Cys"))]<-"MFNG:Arg302Cys"

all$gene[which(all$gene %in% c("IFIH1:.","ADAM17:.","ORC3:.","TSPO;TTLL12:TSPO:Thr147Ala"))]<-NA

all$gene[which(all$MarkerName %in% c("chr7:592301:G:A"))]<-"PDGFA"
all$gene[which(all$MarkerName %in% c("chr22:39264824:T:C"))]<-"PDGFB"
all$gene[which(all$MarkerName %in% c("chr2:218120526:G:C"))]<-"CXCR1"
all$gene[which(all$MarkerName %in% c("chr3:122386581:A:G"))]<-"KPNA1"


table(all$gene)


p1<-ggplot(all, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
  geom_errorbar(aes(ymin=confint_low, ymax=confint_up)) +
  geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
  scale_color_manual(values=col_all) +
  scale_shape_manual(values=c(19,20,19,20,19,20,19,20),guide=F)  +
  theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
  geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=all %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
  scale_x_continuous(trans='reverse',limits=c(-log10(0.001),-log10(0.6)),
  breaks = -log10(as.numeric(c("0.001","0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
  label = c("0.001","0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
  scale_y_continuous(limits=c(0.2,3),breaks = c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3))

pdf("~/git/IIBDGC_GWAS/plots/paper_figures/Figure_1_signals_effect_phenotype_NEW_EXONIC.pdf",width=10, height=5)
p1
dev.off()


p1<-ggplot(all, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
  geom_errorbar(aes(ymin=confint_low, ymax=confint_up)) +
  geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
  scale_color_manual(values=col_all) +
  scale_shape_manual(values=c(19,20,19,20,19,20,19,20),guide=F)  +
  theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
  scale_x_continuous(trans='reverse',limits=c(-log10(0.001),-log10(0.6)),
  breaks = -log10(as.numeric(c("0.001","0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
  label = c("0.001","0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
  scale_y_continuous(limits=c(0.2,3),breaks = c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3))

pdf("~/git/IIBDGC_GWAS/plots/paper_figures/Figure_1_signals_effect_phenotype_no_labels.pdf",width=10, height=5)
p1
dev.off()


# q("no")



all[which(all$MAF<0.002),]
#              MarkerName old_region_signal_name FM study BETA_cond_cd SE_cond_cd
# 192  chr16:50537867:C:T                                           NA         NA
# 99  chr11:117998788:C:T                                           NA         NA
#     Pvalue_cond_cd BETA_cond_ibd SE_cond_ibd Pvalue_cond_ibd BETA_cond_uc
# 192             NA            NA          NA              NA           NA
# 99              NA      0.714744    0.117208     1.07345e-09           NA
#     SE_cond_uc Pvalue_cond_uc chr dbsnp154_eur_tier_1 Position_b38_eur_tier_1
# 192         NA             NA  16         rs145005832                50537867
# 99          NA             NA  11          rs56143179               117998788
#     A1_eur_tier_1 A2_cd_eur_tier_1 BETA_cd_eur_tier_1 SE_cd_eur_tier_1
# 192             C                T             0.7130           0.1254
# 99              C                T             0.7816           0.1864
#     P-value_cd_eur_tier_1 Direction_ed_cd_eur_tier_1 HetISq_cd_eur_tier_1
# 192             1.288e-08                   +?+-++++                 37.7
# 99              2.767e-05                   --++-++-                 15.7
#     HetChiSq_cd_eur_tier_1 HetDf_cd_eur_tier_1 HetPVal_cd_eur_tier_1
# 192                  9.627                   6                0.1413
# 99                   8.299                   7                0.3070
#     rate_Neff_cd_eur_tier_1 BETA_ibd_eur_tier_1 SE_ibd_eur_tier_1
# 192               0.9964531              0.5486            0.1184
# 99                0.9909979              0.6740            0.1524
#     P-value_ibd_eur_tier_1 Direction_ed_ibd_eur_tier_1 HetISq_ibd_eur_tier_1
# 192              3.614e-06                   +??+-++++                  17.2
# 99               9.721e-06                   -+-?++++-                   7.2
#     HetChiSq_ibd_eur_tier_1 HetDf_ibd_eur_tier_1 HetPVal_ibd_eur_tier_1
# 192                   7.243                    6                 0.2989
# 99                    7.539                    7                 0.3750
#     rate_Neff_ibd_eur_tier_1 BETA_uc_eur_tier_1 SE_uc_eur_tier_1
# 192                0.8975930             0.0649           0.1648
# 99                 0.9464498             0.6460           0.1927
#     P-value_uc_eur_tier_1 Direction_ed_uc_eur_tier_1 HetISq_uc_eur_tier_1
# 192             0.6939000                  +???-+-++                  0.0
# 99              0.0007994                  -+-??+++-                 29.9
#     HetChiSq_uc_eur_tier_1 HetDf_uc_eur_tier_1 HetPVal_uc_eur_tier_1
# 192                  1.190                   5                0.9459
# 99                   8.565                   6                0.1996
#     rate_Neff_uc_eur_tier_1 Position_b38_eur_tier_2 A1_eur_tier_2 A2_eur_tier_2
# 192               0.8378729                50537867             C             T
# 99                0.7752097               117998788             C             T

#             updated_region class             class_signal variants_in_ld_r2_0.1
# 192   16_49599118_51605296   old new_fm_credible_set_lead                      
# 99  11_117488494_119415274   old    new_cojo_unsupervised                      
#      region_mb         region_label trait trait_CS trait_variant     variant
# 192 49.6-51.61 16_49599118_51605296    CD       CD            CD rs145005832
# 99                                                                          
#     position CS_index nSNP_CS purity_min alpha_lead BETA_lead SE_lead
# 192 50537867        9       4   0.841036  0.4888095     0.713  0.1254
# 99        NA       NA      NA         NA         NA        NA      NA
#     p_value_lead Cond_p_lead    MAF_lead rate_total_sample_lead
# 192    1.288e-08 5.32401e-12 0.003578782              0.9610245
# 99            NA          NA          NA                     NA
#     rate_total_effective_sample_lead    MAF_mean rate_total_sample_mean
# 192                        0.9900254 0.003621184              0.9501945
# 99                                NA          NA                     NA
#     rate_total_effective_sample_mean lead     alpha  BETA     SE   p_value
# 192                        0.9737459 TRUE 0.4888095 0.713 0.1254 1.288e-08
# 99                                NA   NA        NA    NA     NA        NA
#     novel pre_Pro pre_p_multi pre_tier2 pre_trait pre_trait_reassigned
# 192  TRUE       0           0         .         .                    .
# 99     NA      NA          NA                                         
#           Func           Gene             GeneDetail Exonic_Func AAChange
# 192 intergenic LINC02178;NKD1 dist=142734;dist=10529           .       .:
# 99                                                                       
#     Canonical AAChange_VEP      A2freq         MAF N_cases A2freq_cases
# 192         -            - 0.003578782 0.001691903   35369  0.005395106
# 99                                  NA 0.001065124      NA           NA
#     N_controls A2freq_controls rate_total_sample N_effective
# 192      34312     0.001706508         0.9610245       64615
# 99          NA              NA                NA          NA
#     rate_total_effective_sample       region_fm_2025 phenotype
# 192                   0.9900254                    Y        UC
# 99                           NA N_filtered_out_preFM        CD
#                     class_signal_postfm eqtl_ibd eqtl_cd eqtl_uc eqtl_gene_ibd
# 192 new_fm_credible_set_lead_new_signal                                       
# 99     new_cojo_unsupervised_new_signal                                       
#     eqtl_gene_cd eqtl_gene_uc exonic_variant_in_FM_credible_set
# 192                                                            
# 99                                                             
#     exonic_variant_in_FM_credible_set_gene_aac
# 192                                           
# 99                                            
#     exonic_variant_in_FM_credible_set_status
# 192                                         
# 99                                          
#     exonic_variant_in_FM_credible_set_associated_phenotype
# 192                                                       
# 99                                                        
#                          exonic_variant_in_ld exonic_variant_in_ld_gene_aac
# 192                                                                        
# 99  chr11:117998788:C:T|chr11:117998788:C:T|1                  IL10RA:P295L
#     exonic_variant_in_ld_status exonic_variant_in_ld_associated_phenotype
# 192                                                                      
# 99             Reported_variant                                       IBD
#                                 class_signal_postfm_exome
# 192                   new_fm_credible_set_lead_new_signal
# 99  new_cojo_unsupervised_new_signal_known_exonic_variant
#     avgA2FREQ_ibd_eur_tier_2 avgA2FREQ_CASES_ibd_eur_tier_2
# 192              0.003426183                    0.004461089
# 99               0.001396386                    0.002121627
#     avgA2FREQ_CONTROLS_ibd_eur_tier_2 avgA2FREQ_cd_eur_tier_2
# 192                       0.001691903             0.003578782
# 99                        0.001065124             0.001331707
#     avgA2FREQ_CASES_cd_eur_tier_2 avgA2FREQ_CONTROLS_cd_eur_tier_2
# 192                   0.005395106                      0.001706508
# 99                    0.002095357                      0.001006346
#     avgA2FREQ_uc_eur_tier_2 avgA2FREQ_CASES_uc_eur_tier_2
# 192             0.002299560                   0.003019799
# 99              0.001339885                   0.002227425
#     avgA2FREQ_CONTROLS_uc_eur_tier_2  class_region Effect StdErr   pvalue
# 192                      0.001814994   New signals 0.0649 0.1648   0.6939
# 99                       0.001137848 Known signals 0.7768 0.1478 5.79e-10
#                pheno_class       OR confint_low confint_up gene
# 192   New UC signal (N=52) 1.067052   0.7725118   1.473894   NA
# 99  Known CD signal (N=47) 2.174503   1.6276089   2.905159   NA


############################################################################################################################################################
############################################################################################################################################################

# PLOTS FOR ASHG'25 - IIBDGC calll


#### SIGNALS RECLASSIFIED AS KNOWN:
all$gene<-""


table(all$class_signal_postfm_exome,useNA="ifany")

#                         new_cojo_unsupervised_known_signal 
#                                                         77 # known no-FM
#                           new_cojo_unsupervised_new_signal 
#                                                        195 # new no-FM
#      new_cojo_unsupervised_new_signal_known_exonic_variant 
#                                                          2 # known no-FM
#        new_cojo_unsupervised_new_signal_new_exonic_variant 
#                                                          6 # new no-FM
#                      new_fm_credible_set_lead_known_signal 
#                                                        173 # known FM
  #                      new_fm_credible_set_lead_new_signal 
  #                                                       72  # new FM
  # new_fm_credible_set_lead_new_signal_known_exonic_variant 
  #                                                        2  # known FM
  #   new_fm_credible_set_lead_new_signal_new_exonic_variant 
  #                                                        2  # new FM


all[which(all$class_signal_postfm_exome %in% c("new_cojo_unsupervised_new_signal_known_exonic_variant")),c("MarkerName","class_signal_postfm_exome","exonic_variant_in_ld_gene_aac")]
#              MarkerName                             class_signal_postfm_exome
# 519   chr8:21909729:G:A new_cojo_unsupervised_new_signal_known_exonic_variant
# 99  chr11:117998788:C:T new_cojo_unsupervised_new_signal_known_exonic_variant
#     exonic_variant_in_ld_gene_aac
# 519                    DOK2:P274L
# 99                   IL10RA:P295L



all$gene[which(all$MarkerName=="chr11:117998788:C:T")]<-"IL10RA:P295L"
all$gene[which(all$MarkerName=="chr8:21909729:G:A")]<-"DOK2:P274L"


# reclassified as known
all[which(all$class_signal_postfm_exome %in% c("new_fm_credible_set_lead_new_signal_known_exonic_variant")),c("MarkerName","class_signal_postfm_exome","exonic_variant_in_ld_gene_aac")]
#             MarkerName                                class_signal_postfm_exome
# 252 chr19:10354167:G:A new_fm_credible_set_lead_new_signal_known_exonic_variant
# 195 chr16:50712018:C:T new_fm_credible_set_lead_new_signal_known_exonic_variant
#     exonic_variant_in_ld_gene_aac
# 252                    TYK2:A928V
# 195                    NOD2:R703C

all$gene[which(all$MarkerName=="chr19:10354167:G:A")]<-"TYK2:A928V"
all$gene[which(all$MarkerName=="chr16:50712018:C:T")]<-"NOD2:R703C"


p1<-ggplot(all, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
  geom_errorbar(aes(ymin=confint_low, ymax=confint_up)) +
  geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
  scale_color_manual(values=col_all) +
  scale_shape_manual(values=c(19,20,19,20,19,20,19,20),guide=F)  +
  theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
  geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=all %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
  scale_x_continuous(trans='reverse',limits=c(-log10(0.001),-log10(0.6)),
  breaks = -log10(as.numeric(c("0.001","0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
  label = c("0.001","0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
  scale_y_continuous(limits=c(0.2,3),breaks = c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3))

pdf("~/git/IIBDGC_GWAS/plots/paper_figures/Figure_signals_effect_phenotype_with_gene_names_exonic_reclassified_as_known.pdf",width=10, height=5)
p1
dev.off()


#### NEW SIGNALS / NEW VARIANTS FORM WES:


all[which(all$class_signal_postfm_exome %in% c("new_fm_credible_set_lead_new_signal_new_exonic_variant",
"new_cojo_unsupervised_new_signal_new_exonic_variant")),c("MarkerName","class_signal_postfm_exome","exonic_variant_in_ld_gene_aac")]
#             MarkerName                              class_signal_postfm_exome
# 253 chr19:10359299:A:C new_fm_credible_set_lead_new_signal_new_exonic_variant
# 201 chr16:67627037:C:T new_fm_credible_set_lead_new_signal_new_exonic_variant
# 456 chr6:149413099:A:G    new_cojo_unsupervised_new_signal_new_exonic_variant
# 518  chr8:21853230:T:G    new_cojo_unsupervised_new_signal_new_exonic_variant
# 351 chr22:37470025:G:A    new_cojo_unsupervised_new_signal_new_exonic_variant
# 324   chr2:9521321:A:G    new_cojo_unsupervised_new_signal_new_exonic_variant
# 321  chr2:68813372:G:A    new_cojo_unsupervised_new_signal_new_exonic_variant
# 157 chr14:52314795:T:G    new_cojo_unsupervised_new_signal_new_exonic_variant
#     exonic_variant_in_ld_gene_aac
# 253                    TYK2:I684S
# 201                 CARMIL2:V181M
# 456                 ZC3H12D:K106R
# 518                    DOK2:L138S
# 351           MFNG:MFNG:Arg302Cys
# 324                      ADAM17:.
# 321                ARHGAP25:E254K
# 157                   PTGER2:C83G

all$gene<-""

all$gene[which(all$MarkerName=="chr6:149413099:A:G")]<-"ZC3H12D:K106R"
all$gene[which(all$MarkerName=="chr8:21853230:T:G")]<-"DOK2:L138S"
all$gene[which(all$MarkerName=="chr22:37470025:G:A")]<-"MFNG:R302C"
all$gene[which(all$MarkerName=="chr2:68813372:G:A")]<-"ARHGAP25:E254K"
all$gene[which(all$MarkerName=="chr14:52314795:T:G")]<-"PTGER2:C83G"
all$gene[which(all$MarkerName=="chr2:9521321:A:G")]<-"ADAM17"

all$gene[which(all$MarkerName=="chr16:67627037:C:T")]<-"CARMIL2:V181M"
all$gene[which(all$MarkerName=="chr19:10359299:A:C")]<-"TYK2:I684S"


p1<-ggplot(all, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
  geom_errorbar(aes(ymin=confint_low, ymax=confint_up)) +
  geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
  scale_color_manual(values=col_all) +
  scale_shape_manual(values=c(19,20,19,20,19,20,19,20),guide=F)  +
  theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
  geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=all %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
  scale_x_continuous(trans='reverse',limits=c(-log10(0.001),-log10(0.6)),
  breaks = -log10(as.numeric(c("0.001","0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
  label = c("0.001","0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
  scale_y_continuous(limits=c(0.2,3),breaks = c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3))

pdf("~/git/IIBDGC_GWAS/plots/paper_figures/Figure_signals_effect_phenotype_with_gene_names_exonic_new.pdf",width=10, height=5)
p1
dev.off()


###### all signals likely driven by coding variants

all$gene<-all$exonic_variant_in_ld_gene_aac


p1<-ggplot(all, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
  geom_errorbar(aes(ymin=confint_low, ymax=confint_up)) +
  geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
  scale_color_manual(values=col_all) +
  scale_shape_manual(values=c(19,20,19,20,19,20,19,20),guide=F)  +
  theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
  geom_label_repel(size=2,point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=all %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
  scale_x_continuous(trans='reverse',limits=c(-log10(0.001),-log10(0.6)),
  breaks = -log10(as.numeric(c("0.001","0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
  label = c("0.001","0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
  scale_y_continuous(limits=c(0.2,3),breaks = c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3))

pdf("~/git/IIBDGC_GWAS/plots/paper_figures/Figure_signals_effect_phenotype_with_gene_names_exonic_all.pdf",width=10, height=5)
p1
dev.off()

############################################################################################################################################################
############################################################################################################################################################

# PLOTS FOR ASHG'25

all1<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_with_fm_signals_no_overlap_with_coloc_eqtl.tsv.gz",sep=""))
all2<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_regions_no_conditional_significant_eur_tier2.tsv.gz",sep=""))

all<-rbind(all1,all2,fill=TRUE)
all<-as.data.frame(all)
dim(all)
# [1] 618 217

table(all$class_signal_postfm,useNA="ifany")
#    new_cojo_unsupervised_known_signal      new_cojo_unsupervised_new_signal 
#                                    77                                   209 
# new_fm_credible_set_lead_known_signal   new_fm_credible_set_lead_new_signal 
#                                   192                                    76 
#                                  <NA> 
#                                    64  - the ones not signif in EUR_tier2


# get the average frequency in EUR tier2 in controls:

chr<-as.numeric(unlist(labels(table(all$chr))))
chr<-chr[which(!is.na(chr))]

pheno<-c("ibd","cd","uc")

for (i in 1:length(pheno)) {
  
  print(pheno[i])
  
  for (j in 1:length(chr)) {
    
    tmp1<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr[j],"_",pheno[i],"_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz",sep=""),head=T)
    tmp1<-tmp1[which(tmp1$MarkerName %in% all$MarkerName),]

    tmp1<-tmp1[,c("MarkerName","avgA2FREQ","avgA2FREQ_CASES","avgA2FREQ_CONTROLS")]
    
    if(j==1) {
      tmp<-tmp1
    } else {
      tmp<-rbind(tmp,tmp1)
    }
    
    rm(tmp1)
    
  }

  colnames(tmp)[2:4]<-paste(colnames(tmp)[2:4],pheno[i],"eur_tier_2",sep="_")
  
  if(i==1) {
    dat<-tmp
  } else {
    dat<-merge(dat,tmp,by="MarkerName",sort=F)
  }
  rm(tmp)
}



all<-merge(all,dat,by="MarkerName")

all$MAF<-pmin(all$avgA2FREQ_CONTROLS_ibd_eur_tier_2,1-all$avgA2FREQ_CONTROLS_ibd_eur_tier_2)
all$MAF[which(is.na(all$MAF))]<-pmin(all$avgA2FREQ_CONTROLS_ibd_eur_tier_2[which(is.na(all$MAF))],1-all$avgA2FREQ_CONTROLS_ibd_eur_tier_2[which(is.na(all$MAF))])



# remove Japanese and Korean loci:
all<-all[which(!all$MarkerName %in% c("chr2:5523876:G:T","chr4:38323415:T:C","chr9:114783386:C:T","chr10:110426390:C:T","chr13:40983974:A:G")),]

# remove Kenny 2012, not sign in de lange:
all<-all[which(!all$MarkerName %in% c("chr5:102611094:A:G")),]


#    new_cojo_unsupervised_known_signal      new_cojo_unsupervised_new_signal 
#                                    77                                   209 
# new_fm_credible_set_lead_known_signal   new_fm_credible_set_lead_new_signal 
#                                   192                                    76 
#                                  <NA> 
#                                    64  - the ones not signif in EUR_tier2

table(all$class_signal_postfm,useNA="ifany")

all$class_region<-NA
all$class_region[which(is.na(all$class_signal_postfm))]<-"Known signals"
all$class_region[which(all$class_signal_postfm %in% c("new_cojo_unsupervised_known_signal","new_fm_credible_set_lead_known_signal")) ]<-"Known signals"
all$class_region[which(all$class_signal_postfm %in% c("new_cojo_unsupervised_new_signal","new_fm_credible_set_lead_new_signal"))]<-"New signals"

table(all$class_signal_postfm,all$class_region,useNA="ifany")
#                                         Known signals New signals
#   new_cojo_unsupervised_known_signal               77           0
#   new_cojo_unsupervised_new_signal                  0         209
#   new_fm_credible_set_lead_known_signal           191           0
#   new_fm_credible_set_lead_new_signal               0          76
#   <NA>                                             59           0

table(all$phenotype,useNA="ifany")
#                  CD   IBD_saturated IBD_unsaturated              UC 
#  17             115             144             248              88 

list_ids<-all$MarkerName[which(all$phenotype=="")]

for (i in 1:length(list_ids)) {

    tmp<-as.data.frame(t(all[which(all$MarkerName==list_ids[i]),c("P-value_uc_eur_tier_2","P-value_cd_eur_tier_2","P-value_ibd_eur_tier_2")]))
    colnames(tmp)<-"ID"
    tmp<-rownames(tmp[which(tmp$ID==min(tmp$ID)),,drop=F])
    
    if (tmp=="P-value_uc_eur_tier_2") {
        all$phenotype[which(all$MarkerName==list_ids[i])]<-"UC"
    } else if (tmp=="P-value_cd_eur_tier_2"){
        all$phenotype[which(all$MarkerName==list_ids[i])]<-"CD"
    } else if (tmp=="P-value_ibd_eur_tier_2"){
        all$phenotype[which(all$MarkerName==list_ids[i])]<-"IBD_unsaturated"
    }
    
}

table(all$phenotype,useNA="ifany")
#  CD   IBD_saturated IBD_unsaturated              UC 
# 121             144             255              92 



all$Effect<-NA
all$StdErr<-NA
all$pvalue<-NA


all$Effect[which(all$phenotype=="IBD_unsaturated")]<-all$BETA_ibd_eur_tier_2[which(all$phenotype=="IBD_unsaturated")]
all$StdErr[which(all$phenotype=="IBD_unsaturated")]<-all$SE_ibd_eur_tier_2[which(all$phenotype=="IBD_unsaturated")]
all$pvalue[which(all$phenotype=="IBD_unsaturated")]<-all$"P-value_ibd_eur_tier_2"[which(all$phenotype=="IBD_unsaturated")]

all$Effect[which(all$phenotype=="CD")]<-all$BETA_cd_eur_tier_2[which(all$phenotype=="CD")]
all$StdErr[which(all$phenotype=="CD")]<-all$SE_cd_eur_tier_2[which(all$phenotype=="CD")]
all$pvalue[which(all$phenotype=="CD")]<-all$"P-value_ibd_eur_tier_2"[which(all$phenotype=="CD")]

all$Effect[which(all$phenotype=="UC")]<-all$BETA_uc_eur_tier_2[which(all$phenotype=="UC")]
all$StdErr[which(all$phenotype=="UC")]<-all$SE_uc_eur_tier_2[which(all$phenotype=="UC")]
all$pvalue[which(all$phenotype=="UC")]<-all$"P-value_uc_eur_tier_2"[which(all$phenotype=="UC")]


test<-all[which(all$phenotype=="IBD_saturated"),c("MarkerName","P-value_ibd_eur_tier_2","P-value_cd_eur_tier_2","P-value_uc_eur_tier_2")]
test$model<-NA
test$model[which(test$"P-value_cd_eur_tier_2"<test$"P-value_uc_eur_tier_2")]<-"CD"
test$model[which(test$"P-value_cd_eur_tier_2">test$"P-value_uc_eur_tier_2")]<-"UC"
table(test$model)
# CD UC 
# 72 72  
cd_model<-test$MarkerName[which(test$model=="CD")]
uc_model<-test$MarkerName[which(test$model=="UC")]


# list of IBD-sat variants more significant for CD
all$Effect[which(all$MarkerName %in% cd_model)]<-all$BETA_cd_eur_tier_2[which(all$MarkerName %in% cd_model)]
all$StdErr[which(all$MarkerName %in% cd_model)]<-all$SE_cd_eur_tier_2[which(all$MarkerName %in% cd_model)]
all$pvalue[which(all$MarkerName %in% cd_model)]<-all$"P-value_cd_eur_tier_2"[which(all$MarkerName %in% cd_model)]

# list of IBD-sat variants more significant for UC
all$Effect[which(all$MarkerName %in% uc_model)]<-all$BETA_uc_eur_tier_2[which(all$MarkerName %in% uc_model)]
all$StdErr[which(all$MarkerName %in% uc_model)]<-all$SE_uc_eur_tier_2[which(all$MarkerName %in% uc_model)]
all$pvalue[which(all$MarkerName %in% uc_model)]<-all$"P-value_uc_eur_tier_2"[which(all$MarkerName %in% uc_model)]


table(all$class_region)
# Known signals   New signals 
#           327           285


table(all$class_region,all$phenotype,useNA="ifany")
#                CD IBD_saturated IBD_unsaturated  UC
#   Known signals  61           107             120  39
#   New signals    60            37             135  53

x<-table(all$class_region,all$phenotype,useNA="ifany")

all$pheno_class<-"NA"
all$pheno_class<-as.character(all$pheno_class)
all$pheno_class[which(all$class_region=="New signals" & all$phenotype=="IBD_saturated")]<-paste("New IBD saturated signal (N=",x[2,2],")",sep="")
all$pheno_class[which(all$class_region=="New signals" & all$phenotype=="IBD_unsaturated")]<-paste("New IBD unsaturated signal (N=",x[2,3],")",sep="")
all$pheno_class[which(all$class_region=="New signals" & all$phenotype=="CD")]<-paste("New CD signal (N=",x[2,1],")",sep="")
all$pheno_class[which(all$class_region=="New signals" & all$phenotype=="UC")]<-paste("New UC signal (N=",x[2,4],")",sep="")
all$pheno_class[which(all$class_region=="Known signals" & all$phenotype=="IBD_saturated")]<-paste("Known IBD saturated signal (N=",x[1,2],")",sep="")
all$pheno_class[which(all$class_region=="Known signals" & all$phenotype=="IBD_unsaturated")]<-paste("Known IBD unsaturated signal (N=",x[1,3],")",sep="")
all$pheno_class[which(all$class_region=="Known signals" & all$phenotype=="CD")]<-paste("Known CD signal (N=",x[1,1],")",sep="")
all$pheno_class[which(all$class_region=="Known signals" & all$phenotype=="UC")]<-paste("Known UC signal (N=",x[1,4],")",sep="")

table(all$pheno_class,all$phenotype,useNA="ifany")
#                                         CD IBD_saturated IBD_unsaturated  UC
#   Known CD signal (N=61)                61             0               0   0
#   Known IBD saturated signal (N=107)     0           107               0   0
#   Known IBD unsaturated signal (N=120)   0             0             120   0
#   Known UC signal (N=39)                 0             0               0  39
#   New CD signal (N=60)                  60             0               0   0
#   New IBD saturated signal (N=37)        0            37               0   0
#   New IBD unsaturated signal (N=135)     0             0             135   0
#   New UC signal (N=53)                   0             0               0  53

all$pheno_class<-factor(all$pheno_class,levels = (c(paste("New CD signal (N=",x[2,1],")",sep=""),paste("Known CD signal (N=",x[1,1],")",sep=""),
                                                    paste("New UC signal (N=",x[2,4],")",sep=""),paste("Known UC signal (N=",x[1,4],")",sep="")
                                                    ,paste("New IBD saturated signal (N=",x[2,2],")",sep=""),paste("Known IBD saturated signal (N=",x[1,2],")",sep=""),
                                                    paste("New IBD unsaturated signal (N=",x[2,3],")",sep=""),paste("Known IBD unsaturated signal (N=",x[1,3],")",sep=""))))
table(all$pheno_class)
#                 New CD signal (N=60)               Known CD signal (N=61) 
#                                   60                                   61 
#                 New UC signal (N=53)               Known UC signal (N=39) 
#                                   53                                   39 
#      New IBD saturated signal (N=37)   Known IBD saturated signal (N=107) 
#                                   37                                  107 
#   New IBD unsaturated signal (N=135) Known IBD unsaturated signal (N=120) 
#                                  135                                  120 

cols<-c('#004488','#BB5566','#DDAA33',"#db7107")
col_ligth<-lighten(cols, 0.70, space = "HCL")
# col_dark<-darken(cols,0.2,space = "combined")
col_all<-c(cols[1],col_ligth[1],cols[2],col_ligth[2],cols[3],col_ligth[3],cols[4],col_ligth[4])

all$OR<-exp((all$Effect))
all$confint_low<-exp((all$Effect)-(1.96*all$StdErr))
all$confint_up<-exp((all$Effect)+(1.96*all$StdErr))

all<-all[order(all$MAF,decreasing=T),]
# all<-all[order(as.character(all$pheno_class),decreasing=F),]

all$gene<-NA

all[which(all$OR>4),]


p1<-ggplot(all, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
  geom_errorbar(aes(ymin=confint_low, ymax=confint_up)) +
  geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
  scale_color_manual(values=col_all) +
  scale_shape_manual(values=c(19,20,19,20,19,20,19,20),guide=F)  +
  theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
  geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=all %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
  scale_x_continuous(trans='reverse',limits=c(-log10(0.002),-log10(0.6)),
  breaks = -log10(as.numeric(c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
  label = c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
  scale_y_continuous(limits=c(0.2,2.6),breaks = c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6))

pdf("~/git/IIBDGC_GWAS/plots/metaanalysis/plot_effect_size_known_and_new_signals_no_gene_names.pdf",width=10, height=5)
p1
dev.off()

all$gene[which(all$MarkerName=="chr7:100701361:A:G")]<-"SERPINE1"
all$gene[which(all$MarkerName=="chr19:43648948:A:G")]<-"PLAUR"
all$gene[which(all$MarkerName=="chr10:73934322:G:GCT")]<-"PLAU"





p1<-ggplot(all, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
  geom_errorbar(aes(ymin=all$confint_low, ymax=all$confint_up)) +
  geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
  scale_color_manual(values=col_all) +
  scale_shape_manual(values=c(19,20,19,20,19,20,19,20),guide=F)  +
  theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
  geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=all %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
  scale_x_continuous(trans='reverse',limits=c(-log10(0.002),-log10(0.6)),
  breaks = -log10(as.numeric(c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
  label = c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
  scale_y_continuous(limits=c(0.2,2.6),breaks = c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6))

pdf("~/git/IIBDGC_GWAS/plots/metaanalysis/plot_effect_size_known_and_new_signals_with_plasminogen_gene_names.pdf",width=10, height=5)
p1
dev.off()




# # plot first the known signals - light colours

old<-all[which(all$class_region=="Known signals"),]

p1<-ggplot(old, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
  geom_errorbar(aes(ymin=confint_low, ymax=confint_up)) +
  geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
  scale_color_manual(values=cols) +
  scale_shape_manual(values=c(19,20,19,20,19,20,19,20),guide=F)  +
  theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
  geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=all %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
  scale_x_continuous(trans='reverse',limits=c(-log10(0.002),-log10(0.6)),
  breaks = -log10(as.numeric(c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
  label = c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
  scale_y_continuous(limits=c(0.2,2.6),breaks = c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6))

pdf("~/git/IIBDGC_GWAS/plots/metaanalysis/plot_effect_size_known_signals_no_gene_names_dark.pdf",width=10, height=4.5)
p1
dev.off()


p1<-ggplot(old, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
  geom_errorbar(aes(ymin=confint_low, ymax=confint_up)) +
  geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
  scale_color_manual(values=col_ligth) +
  scale_shape_manual(values=c(19,20,19,20,19,20,19,20),guide=F)  +
  theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
  geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=all %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
  scale_x_continuous(trans='reverse',limits=c(-log10(0.002),-log10(0.6)),
  breaks = -log10(as.numeric(c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
  label = c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
  scale_y_continuous(limits=c(0.2,2.6),breaks = c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6))

pdf("~/git/IIBDGC_GWAS/plots/metaanalysis/plot_effect_size_known_signals_no_gene_names_light.pdf",width=10, height=4.5)
p1
dev.off()










# all$gene<-NA
# # all$gene[which(all$MarkerName=="chr11:308059:T:C")]<-"IFITM2"
# all$gene[which(all$MarkerName=="chr10:71705481:G:A")]<-"VSIR"
# # all$gene[which(all$MarkerName=="chr10:14011254:G:T")]<-"FRMD4A"
# # all$gene[which(all$MarkerName=="chr12:807144:G:A")]<-"WNK1"
# # all$gene[which(all$MarkerName=="chr17:83088537:C:T")]<-"METRNL"
# # all$gene[which(all$MarkerName=="chr22:43164424:C:A")]<-"TSPO"
# # all$gene[which(all$MarkerName=="chr4:155676868:C:T")]<-"GUCY1A1|GUCY1B1"
# # all$gene[which(all$MarkerName=="chr6:109055444:C:G")]<-"SESN1"
# # all$gene[which(all$MarkerName=="chr1:169741528:A:G")]<-"SELL"
# # all$gene[which(all$MarkerName=="chr18:63286830:T:C")]<-"BCL2"
# # all$gene[which(all$MarkerName=="chr2:111175424:G:T")]<-"BCL2L11"
# # # all$gene[which(all$MarkerName=="chr19:54203491:G:A")]<-"LILRA6"

# # all$gene[which(all$MarkerName=="chr2:241859600:CG:C")]<-"PDCD1"
# # all$gene[which(all$MarkerName=="chr11:102788454:G:A")]<-"MMP10"


# all[which(!is.na(all$gene)),]

# p1<-ggplot(all, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
#   geom_errorbar(aes(ymin=all$confint_low, ymax=all$confint_up)) +
#   geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
#   scale_color_manual(values=col_all) +
#   scale_shape_manual(values=c(19,20,19,20,19,20,19,20),guide=F)  +
#   theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
#   geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=all %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
#   scale_x_continuous(trans='reverse',limits=c(-log10(0.002),-log10(0.6)),
#   breaks = -log10(as.numeric(c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
#   label = c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
#   scale_y_continuous(limits=c(0.8,3),breaks = c(0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3))

# pdf("~/git/IIBDGC_GWAS/plots/metaanalysis/plot_effect_size_known_regions_and_new_signals_with_gene_names.pdf",width=9, height=4)
# p1
# dev.off()


# # simplified version
# table(all$class_region,all$phenotype,useNA="ifany")
# #                CD IBD_saturated IBD_unsaturated  UC
#   # Known regions  81            81             150  58
#   # New regions    67            43             137  45

# all$phenotype_simplified<-all$phenotype
# all$phenotype_simplified[which(all$phenotype %in% c("IBD_saturated","IBD_unsaturated"))]<-"IBD"
# table(all$class_region,all$phenotype_simplified,useNA="ifany")
#   #                CD IBD  UC
#   # Known regions  81 231  58
#   # New regions    67 180  45



# x<-table(all$class_region,all$phenotype_simplified,useNA="ifany")

# all$pheno_class<-"NA"
# all$pheno_class<-as.character(all$pheno_class)
# all$pheno_class[which(all$class_region=="New regions" & all$phenotype_simplified=="IBD")]<-paste("New IBD signal (N=",x[2,2],")",sep="")
# all$pheno_class[which(all$class_region=="New regions" & all$phenotype_simplified=="CD")]<-paste("New CD signal (N=",x[2,1],")",sep="")
# all$pheno_class[which(all$class_region=="New regions" & all$phenotype_simplified=="UC")]<-paste("New UC signal (N=",x[2,3],")",sep="")
# all$pheno_class[which(all$class_region=="Known regions" & all$phenotype_simplified=="IBD")]<-paste("Known IBD signal (N=",x[1,2],")",sep="")
# all$pheno_class[which(all$class_region=="Known regions" & all$phenotype_simplified=="CD")]<-paste("Known CD signal (N=",x[1,1],")",sep="")
# all$pheno_class[which(all$class_region=="Known regions" & all$phenotype_simplified=="UC")]<-paste("Known UC signal (N=",x[1,3],")",sep="")


# all$pheno_class<-factor(all$pheno_class,levels = (c(paste("New CD signal (N=",x[2,1],")",sep=""),paste("Known CD signal (N=",x[1,1],")",sep=""),
#                                                     paste("New IBD signal (N=",x[2,2],")",sep=""),paste("Known IBD signal (N=",x[1,2],")",sep=""),
#                                                     paste("New UC signal (N=",x[2,3],")",sep=""),paste("Known UC signal (N=",x[1,3],")",sep=""))))
# table(all$pheno_class)

# table(all$pheno_class,all$phenotype_simplified,useNA="ifany")
#   #                           CD IBD  UC
#   # Known CD signal (N=81)    81   0   0
#   # Known IBD signal (N=231)   0 231   0
#   # Known UC signal (N=58)     0   0  58
#   # New CD signal (N=67)      67   0   0
#   # New IBD signal (N=180)     0 180   0
#   # New UC signal (N=45)       0   0  45


# # cols<-c('#BB5566','#004488',"#db7107")
# cols<-c('#1b9e77','#d95f02','#7570b3')
# col_ligth<-lighten(cols, 0.70, space = "HCL")
# # col_dark<-darken(cols,0.2,space = "combined")
# col_all<-c(cols[1],col_ligth[1],cols[2],col_ligth[2],cols[3],col_ligth[3])


# p1<-ggplot(all, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
#   geom_errorbar(aes(ymin=all$confint_low, ymax=all$confint_up)) +
#   geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
#   scale_color_manual(values=col_all) +
#   scale_shape_manual(values=c(19,20,19,20,19,20,19,20),guide=F)  +
#   theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
#   geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=all %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
#   scale_x_continuous(trans='reverse',limits=c(-log10(0.002),-log10(0.6)),
#   breaks = -log10(as.numeric(c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
#   label = c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
#   scale_y_continuous(limits=c(0.8,3),breaks = c(0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3))


# pdf("~/git/IIBDGC_GWAS/plots/metaanalysis/plot_effect_size_known_regions_and_new_signals_with_gene_names_pheno_simpl.pdf",width=9, height=4)
# p1
# dev.off()

# # simplified version with no labels:


# p1<-ggplot(all, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
#   geom_errorbar(aes(ymin=all$confint_low, ymax=all$confint_up)) +
#   geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
#   scale_color_manual(values=col_all) +
#   scale_shape_manual(values=c(19,20,19,20,19,20,19,20),guide=F)  +
#   theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
#   # geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=all %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
#   scale_x_continuous(trans='reverse',limits=c(-log10(0.002),-log10(0.6)),
#   breaks = -log10(as.numeric(c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
#   label = c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
#   scale_y_continuous(limits=c(0.8,3),breaks = c(0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3))


# pdf("~/git/IIBDGC_GWAS/plots/metaanalysis/plot_effect_size_known_regions_and_new_signals_no_gene_names_pheno_simpl.pdf",width=9, height=4)
# p1
# dev.off()


# # simplified version with just known signals:


# old<-all[which(all$class_signal=="old"),]
# table(old$pheno_class,old$phenotype,useNA="ifany")

# old$gene<-NA

# old$gene[which(old$MarkerName=="chr19:10352442:G:C")]<-"TYK2"
# old$gene[which(old$MarkerName=="chr1:67242007:G:A")]<-"IL23R"
# old$gene[which(old$MarkerName=="chr1:67214307:G:A")]<-"IL23R"
# # old$gene[which(old$MarkerName=="chr10:6052734:C:T")]<-"IL2RA"
# # old$gene[which(old$MarkerName=="chr2:181443625:A:G")]<-"ITGA4"
# # old$gene[which(old$MarkerName=="chr7:20537675:G:A")]<-"ITGB8"
# # old$gene[which(old$MarkerName=="chr16:30471173:T:C")]<-"ITGAL"
# # old$gene[which(old$MarkerName=="chr2:102322576:C:T")]<-"IL18RAP"
# old$gene[which(old$MarkerName %in% c("chr16:50660100:A:G","chr16:50712288:G:A","chr16:50722629:G:C","chr16:50728860:T:C","chr16:50729867:G:GC"))]<-"NOD2"




# p1<-ggplot(old, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
#   geom_errorbar(aes(ymin=old$confint_low, ymax=old$confint_up)) +
#   geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
#   scale_color_manual(values=col_ligth) +
#   scale_shape_manual(values=c(20,20,20,20),guide=F)  +
#   theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
#   # geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=old %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
#   scale_x_continuous(trans='reverse',limits=c(-log10(0.002),-log10(0.6)),
#   breaks = -log10(as.numeric(c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
#   label = c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
#   scale_y_continuous(limits=c(0.8,3),breaks = c(0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3))

# pdf("~/git/IIBDGC_GWAS/plots/metaanalysis/plot_effect_size_known_regions_no_gene_names_pheno_simpl.pdf",width=9, height=3.75)
# p1
# dev.off()

# p1<-ggplot(old, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
#   geom_errorbar(aes(ymin=old$confint_low, ymax=old$confint_up)) +
#   geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
#   scale_color_manual(values=cols) +
#   scale_shape_manual(values=c(20,20,20,20),guide=F)  +
#   theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
#   geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=old %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
#     scale_x_continuous(trans='reverse',limits=c(-log10(0.002),-log10(0.6)),
#   breaks = -log10(as.numeric(c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
#   label = c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+  
#   scale_y_continuous(limits=c(0.8,3),breaks = c(0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3))

# pdf("~/git/IIBDGC_GWAS/plots/metaanalysis/plot_effect_size_known_regions_with_gene_names_pheno_simpl.pdf",width=9, height=3.75)
# p1
# dev.off()

# p1<-ggplot(old, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
#   geom_errorbar(aes(ymin=old$confint_low, ymax=old$confint_up)) +
#   geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
#   scale_color_manual(values=cols) +
#   scale_shape_manual(values=c(20,20,20,20),guide=F)  +
#   theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
#   # geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=old %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
#     scale_x_continuous(trans='reverse',limits=c(-log10(0.002),-log10(0.6)),
#   breaks = -log10(as.numeric(c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
#   label = c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+  
#   scale_y_continuous(limits=c(0.8,3),breaks = c(0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3))

# pdf("~/git/IIBDGC_GWAS/plots/metaanalysis/plot_effect_size_known_regions_no_gene_names_bold_colors_pheno_simpl.pdf",width=9, height=3.75)
# p1
# dev.off()



# # # all loci:

# # p1<-ggplot(all, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
# #   geom_errorbar(aes(ymin=all$confint_low, ymax=all$confint_up)) +
# #   geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
# #   scale_color_manual(values=col_all) +
# #   scale_shape_manual(values=c(19,20,19,20,19,20,19,20),guide=F)  +
# #   theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
# #         panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
# #   # geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=all %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
# #   scale_x_continuous(trans='reverse',limits=c(2.50,0.30103),breaks = c(2.49,2,1.69897,1.30103,1,0.69897,0.5228787,0.39794,0.30103),label = c("0.003","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
# #   scale_y_continuous(limits=c(0.9,3),breaks = c(1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3))

# # pdf(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/plots/plot_effect_size_known_regions_and_new_signals.pdf",sep=""),width=9, height=4)
# # p1
# # dev.off()


# # all loci - with names

# p1<-ggplot(all, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
#   geom_errorbar(aes(ymin=all$confint_low, ymax=all$confint_up)) +
#   geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
#   scale_color_manual(values=col_all) +
#   scale_shape_manual(values=c(19,20,19,20,19,20,19,20),guide=F)  +
#   theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
#   geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=all %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
#   scale_x_continuous(trans='reverse',limits=c(-log10(0.002),-log10(0.6)),
#   breaks = -log10(as.numeric(c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
#   label = c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
#   scale_y_continuous(limits=c(0.8,3),breaks = c(0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3))

# pdf("~/git/IIBDGC_GWAS/plots/metaanalysis/plot_effect_size_known_regions_and_new_signals_with_gene_names.pdf",width=9, height=4)
# p1
# dev.off()


# all[which(all$OR>1.6 & all$MAF>0.01),]
# #  MarkerName chr      pos dbsnp154_eur_tier_1 Position_b38_eur_tier_1
# # 1:  chr16:50722629:G:C  16 50722629           rs2066845                50722629
# # 2:   chr1:67242007:G:A   1 67242007          rs11581607                67242007
# # 3:  chr16:50712015:C:T  16 50712015           rs2066844                50712015
# # 4: chr16:50729867:G:GC  16 50729867                                    50729867
# #          rsid V5  old_region_signal_name FM study Region Signal
# # 1:  rs2066845 ok Huang_region196_signal3    Huang    196      3
# # 2: rs11581607 ok   Huang_region7_signal1    Huang      7      1
# # 3:  rs2066844 ok Huang_region196_signal2    Huang    196      2
# # 4:  rs2066847 OK Huang_region196_signal1    Huang    196      1

# # all coding signals in NOD2 and IL23




# # # plot first the known signals - light colours

# old<-all[which(all$class_signal=="old"),]
# table(old$pheno_class,old$phenotype,useNA="ifany")



# old$gene<-NA

# old$gene[which(old$MarkerName=="chr19:10352442:G:C")]<-"TYK2"
# old$gene[which(old$MarkerName=="chr1:67242007:G:A")]<-"IL23R"
# old$gene[which(old$MarkerName=="chr1:67214307:G:A")]<-"IL23R"
# old$gene[which(old$MarkerName=="chr10:6052734:C:T")]<-"IL2RA"
# old$gene[which(old$MarkerName=="chr2:181443625:A:G")]<-"ITGA4"
# old$gene[which(old$MarkerName=="chr7:20537675:G:A")]<-"ITGB8"
# old$gene[which(old$MarkerName=="chr16:30471173:T:C")]<-"ITGAL"
# old$gene[which(old$MarkerName=="chr2:102322576:C:T")]<-"IL18RAP"





# # p1<-ggplot(old, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
# #   geom_errorbar(aes(ymin=old$confint_low, ymax=old$confint_up)) +
# #   geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
# #   scale_color_manual(values=col_ligth) +
# #   scale_shape_manual(values=c(20,20,20,20),guide=F)  +
# #   theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
# #         panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
# #   geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=old %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
# #   scale_x_continuous(trans='reverse',limits=c(2.50,0.30103),breaks = c(2.49,2,1.69897,1.30103,1,0.69897,0.5228787,0.39794,0.30103),label = c("0.003","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
# #   scale_y_continuous(limits=c(0.9,3),breaks = c(1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3))

# # pdf(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/plots/plot_effect_size_known_regions_light_colours.pdf",sep=""),width=9, height=3.75)
# # p1
# # dev.off()



# # # known signals - strong colours


# table(old$pheno_class,old$phenotype,useNA="ifany")
#   #                                       CD IBD_saturated IBD_unsaturated  UC
#   # Known CD signal (N=81)                81             0               0   0
#   # Known IBD saturated signal (N=81)      0            81               0   0
#   # Known IBD unsaturated signal (N=150)   0             0             150   0
#   # Known UC signal (N=58)                 0             0               0  58
#   # New CD signal (N=67)                  67             0               0   0
#   # New IBD saturated signal (N=43)        0            43               0   0
#   # New IBD unsaturated signal (N=137)     0             0             137   0
#   # New UC signal (N=45)                   0             0               0  45




# p1<-ggplot(old, aes(x=-log10(MAF), y=OR,color=pheno_class)) +
#   geom_errorbar(aes(ymin=old$confint_low, ymax=old$confint_up)) +
#   geom_point(aes(shape=pheno_class)) + ylab("Odds Ratio") + xlab("Minor Allele Frequency") +
#   scale_color_manual(values=cols) +
#   scale_shape_manual(values=c(19,20,19,20,19,20,19,20),guide=F)  +
#   theme(legend.position="bottom",legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"),legend.key=element_blank()) +
#   geom_label_repel(point.padding=0.45,min.segment.length=8,nudge_x=0.07,nudge_y=0.07,data=old %>% filter(gene!=""),aes(label=gene),show.legend = FALSE,fontface = "italic") +
#   scale_x_continuous(trans='reverse',limits=c(-log10(0.002),-log10(0.6)),
#   breaks = -log10(as.numeric(c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))),
#   label = c("0.002","0.003","0.005","0.01","0.02","0.05","0.1","0.2","0.3","0.4","0.5"))+
#   scale_y_continuous(limits=c(0.8,3),breaks = c(0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3))



# pdf("~/git/IIBDGC_GWAS/plots/metaanalysis/plot_effect_size_known_regions_strong_colours.pdf",width=9, height=3.75)
# p1
# dev.off()