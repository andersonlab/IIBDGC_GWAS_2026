# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# singularity exec iibdgc_postprocess_10_singularity.sif

# MEM=15000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q normal R 

library(data.table)
library(ggplot2)
library(ggpubr)
library(viridis)
library(R.utils)

path_gwas="/path/to/ibdgwas/IIBDGC/"

release<-"eur_tier_1"

print(release)

# load huang results:

hg<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/summary_results/list_credible_sets_huang_2017.csv")
table(hg$tier2)
        #    No single_method 
        #   139            15 

# keep only resutls from tier1 - whether this associaiton is reported in > one method
hg<-hg[which(hg$tier2=="No"),]
summary(hg$Prob.lead)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.0415  0.1824  0.3297  0.4269  0.6454  0.9999 

# for this list of variants, plot pvalue and compare with eur tier 1, with eur tier 2:

summary(hg$p_single)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.240   8.665  14.520     Inf  24.885     Inf

# single are non conditional results
hg$pvalue<-10^(-hg$p_single)

table(hg$trait.reassigned)
#  CD IBD  UC 
#  79  37  23 

# re-name rs5743293 - NOD2(fs1007insC) to its new ID
hg$variant.lead[which(hg$variant.lead=="rs5743293")]<-"rs2066847 (rs5743293)"

# number of variants to extract the same lead variant:
dim(dat[which(dat$variant.lead==dat$dbsnp154_eur_tier_1),])
# [1] 126  60

# re-map to the summary table:

all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_conditional_on_known_index_gwas_fm_2023_with_iibdgc_tier1_tier2_eas_sas_summary_stats_plus_unsupervised_conditional.tsv",sep=""))
all<-all[which(all$BestSNP_FM %in% hg$variant.lead),c("BestSNP_FM","MarkerName","dbsnp154_eur_tier_1")]

hg<-merge(hg,all,by.x="variant.lead",by.y="BestSNP_FM",all.x=T)

pheno<-c("cd","uc","ibd")

for (chr in 1:22) {
  
  print(chr)

  for(i in 1:length(pheno)) {

    print(pheno[i])

    if (release=="eur_tier_1") {
        tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/chr",chr,"_",pheno[i],"_meta_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_with_rsid_dbsnp154.txt.gz",sep=""))
        tmp<-tmp[which(tmp$MarkerName %in% hg$MarkerName[which(hg$trait.reassigned==toupper(pheno[i]))]),]

    } else {

        tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_",release,"_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""))
        tmp<-tmp[which(tmp$MarkerName %in% hg$MarkerName[which(hg$trait.reassigned==toupper(pheno[i]))]),]

        rm(tmp1)

    }

    tmp<-tmp[,c("MarkerName","P-value","BETA","SE","INFO","rate_total_sample","avgA2FREQ_CASES","avgA2FREQ_CONTROLS","N_CASES","N_CONTROLS")]
    tmp$pheno<-pheno[i]

    if (chr==1 & i==1) {
        dat<-tmp
    } else {
        dat<-rbind(dat,tmp)
    }
  }
}

dat<-merge(dat,hg,by="MarkerName")

table(dat$pheno,dat$trait.reassigned)
  #     CD IBD UC
  # cd  78   0  0
  # ibd  0  37  0
  # uc   0   0 23

dat$pvalue_new<-as.numeric(dat$'P-value')
dat$pvalue_huang<-10^(-dat$p_single)
dat$pvalue_huang[which(dat$pvalue_huang==0)]<-1E-320

maxlim<-max(-log10(as.numeric(dat$pvalue_huang)),-log10(as.numeric(dat$pvalue_new)))+10

dat$pp_50<-"no"
dat$pp_50[which(dat$Prob.lead>=0.50)]<-"yes"


dat[which(dat$variant.lead!=dat$dbsnp154_eur_tier_1),c("variant.lead","dbsnp154_eur_tier_1","MarkerName")]
#              variant.lead dbsnp154_eur_tier_1             MarkerName
#                    <char>              <char>                 <char>
#  1:           rs398097167          rs76908010    chr10:99528354:CA:C - in credible set
#  2:           rs200349593            rs661054    chr11:114559688:A:G - in credible set
#  3: rs2066847 (rs5743293)                        chr16:50729867:G:GC # SAME
#  4:             rs3883338          rs12942547     chr17:42375526:A:G - in credible set
#  5:            rs60542850          rs11879191     chr19:10402235:G:A - in credible set
#  6:        chr20:43258079          rs79493594     chr20:44561731:C:T - in credible set
#  7:             rs4456788           rs2838517     chr21:44193942:T:C - in credible set
#  8:            rs33993564          rs11741861     chr5:150898347:A:G - alternative lead variant from old GWAS
#  9:             rs7711427          rs11742570      chr5:40410482:T:C - alternative lead variant from old GWAS
# 10:             rs7723899           rs1363907      chr5:96917099:G:A - in credible set 
# 11:           rs796819847                         chr6:20640188:CT:C - rs11339738, not same
# 12:            rs56691573                       chr7:28151111:ATTT:A - not same!
# 13:            rs34841270          rs16903065     chr8:128528218:C:A - in credible set
# 14:        chr9:117571294           rs6478106     chr9:114783386:C:T - alternative lead variant from old GWAS
# 15:           rs146029108                     chr9:136435514:GTTAT:G - SAME

# manually add pvalue for NOD2fs:
dat$dbsnp154_eur_tier_1[which(dat$variant.lead %in% c("rs2066847 (rs5743293)","rs146029108"))]<-dat$variant.lead[which(dat$variant.lead %in% c("rs2066847 (rs5743293)","rs146029108"))]




p1<-ggplot(d, aes(x=AF_HC, y=avgA2FREQ_CONTROLS)) + 
              geom_point(shape=20) + xlim(0,1) + ylim(0,1) + 
              geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
              ggtitle(paste(toupper(pheno[i]),"Allele frequency in Controls, Huang vs IIBDGC 2022",release)) +
              geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8)

p2<-ggplot(dat[which(dat$pheno=="cd" & dat$variant.lead==dat$dbsnp154_eur_tier_1),], aes(x=AF_CD, y=avgA2FREQ_CASES)) + 
              geom_point(shape=20) + xlim(0,1) + ylim(0,1) + 
              geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
              ggtitle(paste(toupper(pheno[i]),"Allele frequency in CD, Huang vs IIBDGC 2022",release)) +
              geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8)

p3<-ggplot(dat[which(dat$pheno=="uc" & dat$variant.lead==dat$dbsnp154_eur_tier_1),], aes(x=AF_UC, y=avgA2FREQ_CASES)) + 
              geom_point(shape=20) + xlim(0,1) + ylim(0,1) + 
              geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
              ggtitle(paste(toupper(pheno[i]),"Allele frequency in UC, Huang vs IIBDGC 2022",release)) +
              geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8)


p<-ggarrange(p1,p2,p3,ncol=3)


ggsave(paste0("~/git/IIBDGC_GWAS/plots/metaanalysis/allele_frequency_huang_vs_iibdgc_",release,".pdf",sep=""),
    p,
    width = 15,
    height = 6,
    dpi = 200,
    units = c("in"),
    limitsize = T
)


# split by phenotype:

dat<-as.data.frame(dat)

for (i in 1:length(pheno)) {

  p1<-ggplot(dat[which(dat$pheno==pheno[i] & dat$variant.lead==dat$dbsnp154_eur_tier_1),], aes(x=-log10(pvalue_huang), y=(-log10(pvalue_new)), color=Prob.lead)) + 
              geom_point(shape=20) + xlim(0,maxlim) + ylim(0,maxlim) + 
              geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
              ggtitle(paste(toupper(pheno[i]),"P-value Huang vs IIBDGC 2022",release)) + 
              geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8) + scale_color_viridis(limits = c(0, 1),direction = -1,option = "inferno")
    
  p2<-ggplot(dat[which(dat$pheno==pheno[i] & dat$variant.lead==dat$dbsnp154_eur_tier_1),], aes(x=-log10(pvalue_huang), y=(-log10(pvalue_new)), color=INFO)) + 
              geom_point(shape=20) + xlim(0,maxlim) + ylim(0,maxlim) + 
              geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
              ggtitle(paste(toupper(pheno[i]),"P-value Huang vs IIBDGC 2022",release)) +
              geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8) + scale_color_viridis(limits = c(0, 1),direction = -1)

  p3<-ggplot(dat[which(dat$pheno==pheno[i] & dat$variant.lead==dat$dbsnp154_eur_tier_1),], aes(x=-log10(pvalue_huang), y=(-log10(pvalue_new)), color=rate_total_sample)) + 
              geom_point(shape=20) + xlim(0,maxlim) + ylim(0,maxlim) + 
              geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
              ggtitle(paste(toupper(pheno[i]),"P-value Huang vs IIBDGC 2022",release)) +
              geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8) + scale_color_viridis(limits = c(0, 1),direction = -1,option = "plasma")

  p<-ggarrange(p1,p2,p3,ncol=3)

  ggsave(paste0("~/git/IIBDGC_GWAS/plots/comparison_previous_datasets/pvalues_huang_vs_iibdgc_",release,"_",pheno[i],".pdf",sep=""),
      p,
    width = 15,
    height = 4,
      dpi = 200,
      units = c("in"),
      limitsize = T) 
  
}


### plot a more zoom in view:

top<-40

for (i in 1:length(pheno)) {

  p1<-ggplot(dat[which(dat$pheno==pheno[i] & dat$variant.lead==dat$dbsnp154_eur_tier_1),], aes(x=-log10(pvalue_huang), y=(-log10(pvalue_new)), color=Prob.lead)) + 
              geom_point(shape=20) + xlim(0,50) + 
              geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
              ggtitle(paste(toupper(pheno[i]),"P-value Huang vs IIBDGC 2022",release)) + 
              geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8) + scale_color_viridis(limits = c(0, 1),direction = -1,option = "inferno") + 
              xlim(0,top) + ylim(0,top)
    
  p2<-ggplot(dat[which(dat$pheno==pheno[i] & dat$variant.lead==dat$dbsnp154_eur_tier_1),], aes(x=-log10(pvalue_huang), y=(-log10(pvalue_new)), color=INFO)) + 
              geom_point(shape=20) + xlim(0,50) +  
              geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
              ggtitle(paste(toupper(pheno[i]),"P-value Huang vs IIBDGC 2022",release)) +
              geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8) + scale_color_viridis(limits = c(0, 1),direction = -1) + 
              xlim(0,top) + ylim(0,top)

  p3<-ggplot(dat[which(dat$pheno==pheno[i] & dat$variant.lead==dat$dbsnp154_eur_tier_1),], aes(x=-log10(pvalue_huang), y=(-log10(pvalue_new)), color=rate_total_sample)) + 
              geom_point(shape=20) + 
              geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
              ggtitle(paste(toupper(pheno[i]),"P-value Huang vs IIBDGC 2022",release)) +
              geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8) + scale_color_viridis(limits = c(0, 1),direction = -1,option = "plasma") + 
              xlim(0,top) + ylim(0,top)

  p<-ggarrange(p1,p2,p3,ncol=3)

  ggsave(paste0("~/git/IIBDGC_GWAS/plots/comparison_previous_datasets/pvalues_huang_vs_iibdgc_",release,"_",pheno[i],"_pvalue_10-",top,".pdf",sep=""),
      p,
    width = 15,
    height = 4,
      dpi = 200,
      units = c("in"),
      limitsize = T) 
  
}



head(dat)
#          pvalue dbsnp154_eur_tier_1 pvalue_new pvalue_huang  pp_50
#           <num>              <char>      <num>        <num> <char>
# 1: 7.943282e-08           rs2265189  7.915e-08 7.943282e-08     no
# 2: 2.137962e-13           rs2790159  1.360e-04 2.137962e-13     no
# 3: 1.513561e-09          rs12722504  9.120e-01 1.513561e-09    yes
# 4: 6.456542e-24          rs61839660  1.244e-18 6.456542e-24    yes
# 5: 5.248075e-13           rs7915475  1.863e-05 5.248075e-13    yes
# 6: 1.202264e-20           rs1250563  3.855e-23 1.202264e-20     no




dat[which(dat$pp_50=="yes" & dat$pvalue_new>5e-8),c("dbsnp154_eur_tier_1","variant.lead","pvalue_new","pvalue_huang","Prob.lead","INFO","pheno")]
#     dbsnp154_eur_tier_1   variant.lead pvalue_new pvalue_huang Prob.lead
#                  <char>         <char>      <num>        <num>     <num>
#  1:          rs12722504     rs12722504  9.120e-01 1.513561e-09    0.6146
#  2:           rs7915475      rs7915475  1.863e-05 5.248075e-13    0.5284
#  3:            rs630923       rs630923  8.597e-08 1.659587e-08    0.8195
#  4:           rs7307562      rs7307562  5.664e-06 2.511886e-14    0.9999
#  5:          rs35874463     rs35874463  1.952e-02 2.884032e-11    0.9893
#  6:           rs5743271      rs5743271  5.822e-02 9.332543e-09    0.9930
#  7:         rs104895444    rs104895444  3.359e-04 6.760830e-09    0.8647
#  8:         rs104895467    rs104895467  3.318e-06 1.230269e-06    0.8328
#  9:          rs72796367     rs72796367  1.023e-05 9.772372e-22    0.9829
# 10:         rs145530718    rs145530718  8.413e-07 5.011872e-08    0.7625
# 11:          rs76418789     rs76418789  1.169e-03 6.606934e-09    0.9368
# 12:          rs79493594 chr20:43258079  1.668e-02 6.456542e-08    0.7365
# 13:          rs74465132     rs74465132  9.236e-08 6.165950e-10    0.9941
#          INFO  pheno
#         <num> <char>
#  1: 0.9891912     cd
#  2: 0.9920654     cd
#  3: 0.9857865     cd
#  4: 0.9876365     cd
#  5: 0.8935222    ibd
#  6: 0.9434051     cd
#  7: 0.8617293     cd
#  8: 0.9286784     cd
#  9: 0.9559303     cd
# 10: 0.8980263     cd
# 11: 0.9842066     cd
# 12: 0.9376758     cd
# 13: 0.9653756    ibd

dim(dat[which(dat$pvalue_new>dat$pvalue_huang),])
# [1] 76 58
dim(dat[which(dat$pvalue_new<dat$pvalue_huang),])
# [1] 63 58

fwrite(dat[,c("variant.lead","HD","signal","p_multi","p_single",
"trait.reassigned","Prob.lead","dbsnp154_eur_tier_1","MarkerName","BETA","SE","INFO","rate_total_sample","pheno","avgA2FREQ_CASES","avgA2FREQ_CONTROLS","N_CASES","N_CONTROLS",
"pvalue_new","pvalue_huang")],
paste0("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/summary_results/list_credible_sets_huang_2017_with_iibdgc_",release,".tsv.gz"),
col.names=T,row.names=F,quote=F,sep="\t")


# add eur tier 2

for (chr in 1:22) {
  
  print(chr)

  for(i in 1:length(pheno)) {

    print(pheno[i])

    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""))
    tmp<-tmp[which(tmp$MarkerName %in% hg$MarkerName[which(hg$trait.reassigned==toupper(pheno[i]))]),]

    rm(tmp1)

    tmp<-tmp[,c("MarkerName","P-value","BETA","SE","INFO","rate_total_sample","avgA2FREQ_CASES","avgA2FREQ_CONTROLS","N_CASES","N_CONTROLS")]
    tmp$pheno<-pheno[i]

    if (chr==1 & i==1) {
        dat_tier2<-tmp
    } else {
        dat_tier2<-rbind(dat_tier2,tmp)
    }
  }
}



colnames(dat_tier2)[2:ncol(dat_tier2)]<-paste0(colnames(dat_tier2)[2:ncol(dat_tier2)],"_eur_tier_2")
dat<-merge(dat,dat_tier2,by="MarkerName")


fwrite(dat[,c("variant.lead","HD","signal","p_multi","p_single",
"trait.reassigned","Prob.lead","pheno","dbsnp154_eur_tier_1","MarkerName",
"BETA","SE","INFO","rate_total_sample","avgA2FREQ_CASES","avgA2FREQ_CONTROLS","N_CASES","N_CONTROLS",
"BETA_eur_tier_2","SE_eur_tier_2","INFO_eur_tier_2","rate_total_sample_eur_tier_2","avgA2FREQ_CASES_eur_tier_2","avgA2FREQ_CONTROLS_eur_tier_2","N_CASES_eur_tier_2","N_CONTROLS_eur_tier_2",
"pvalue_new","pvalue_huang")],
paste0("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/summary_results/list_credible_sets_huang_2017_with_iibdgc_",release,"_eur_tier_2.tsv.gz"),
col.names=T,row.names=F,quote=F,sep="\t")

# for each one of these variants, plot beta and SE per study:

dat$study<-paste0("Meta-analysis ",release)
dat$p.value<-format(dat$pvalue_new, scientific=TRUE,digits=2)


studies_ibd_uc<-c("gsa","humancoreexome","illuminaexome","quad610","humanomni1","affymetrix500","humanomniexpress","affymetrix6","illumina370","ibdbioresource","finngen","ukbb","decode","danish")
studies_cd<-c("gsa","humancoreexome","illuminaexome","quad610","humanomni1","affymetrix500","affymetrix6","illumina370","ibdbioresource","finngen","ukbb","decode","danish")

studies_tier_2<-c("ibdbioresource","finngen","ukbb","decode","danish")

dat$Neff<-NA

dat$A1FREQ_CASES<-""
dat$A1FREQ_CONTROLS<-""

for (j in 1:nrow(dat)) {

  if(dat$pheno[j]=="cd") {
    studies<-studies_cd
  } else {
    studies<-studies_ibd_uc
  }

  for (k in 1:length(studies)) {

    if(studies[k]=="ukbb") {
      tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/ukbb/",dat$pheno[j],"/chr",dat$chr[j],"_ukbb_eur_step2_",dat$pheno[j],"_eur_sex_PCs_firthse_",dat$pheno[j],".regenie.gz"))
    } else if (studies[k]=="finngen") {
      tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/finngen/",dat$pheno[j],"/chr",dat$chr[j],"_finngen_r10_K11_",dat$pheno[j],".regenie.gz"))
    } else if (studies[k]=="decode") {
      tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/decode/",dat$pheno[j],"/chr",dat$chr[j],"_DECODE_",dat$pheno[j],"_30062021_edited_noMult.regenie"))
    } else if (studies[k]=="danish") {
      tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/danish/",dat$pheno[j],"/chr",dat$chr[j],"_danish_gsa_",dat$pheno[j],"_eur_sex_10PCs_saige_spa.regenie"))
    } else if (studies[k]=="ibdbioresource") {
      tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/ibdbioresource/",dat$pheno[j],"/",dat$chr[j],"_interval_ibdbioresource_eur_step2_",dat$pheno[j],"_eur_sex_PCs_firthse_",dat$pheno[j],".regenie.gz"))
    } else {
      tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/regenie/",studies[k],"/",dat$pheno[j],"/chr",dat$chr[j],"_",studies[k],"_eur_all_step2_",dat$pheno[j],"_eur_sex_PCs_firthse_",dat$pheno[j],".regenie.gz"))
    }


    tmp<-tmp[which(tmp$ID==dat$MarkerName[j]),]

    if (k==1) {
      tmp$study<-studies[k]
      datfp<-tmp
    } else if (nrow(tmp)>0) {
      if(is.na(tmp$LOG10P)) {
        tmp$LOG10P<-pchisq(tmp$BETA/tmp$SE,1,lower.tail=T)
      }
      tmp$study<-studies[k]
      datfp<-rbind(datfp,tmp)
    }
    
  }

  datfp$p.value<-format(10^(-datfp$LOG10P), scientific=TRUE,digits=2)
  datfp$Neff<-round(4/((1/as.numeric(datfp$N_CASES))+(1/as.numeric(datfp$N_CONTROLS))))

  datfp<-datfp[,c("ID","BETA","SE","p.value","study","INFO","Neff","A1FREQ_CASES","A1FREQ_CONTROLS")]
  colnames(datfp)[1]<-"MarkerName"

  tmp<-dat[j,c("MarkerName","BETA_eur_tier_2","SE_eur_tier_2","P-value_eur_tier_2","Neff","INFO_eur_tier_2")]
  tmp$study<-"Meta-analysis eur_tier_2"
  tmp$Neff<-sum(datfp$Neff)
  tmp$A1FREQ_CASES<-""
  tmp$A1FREQ_CONTROLS<-""


  dat$Neff[j]<-round(sum(datfp$Neff[which(!datfp$study %in% studies_tier_2)]))

  datfp<-rbind(datfp,dat[j,c("MarkerName","BETA","SE","p.value","study","INFO","Neff","A1FREQ_CASES","A1FREQ_CONTROLS")])

  tmp<-tmp[,c("MarkerName","BETA_eur_tier_2","SE_eur_tier_2","P-value_eur_tier_2",
  "study","INFO_eur_tier_2","Neff","A1FREQ_CASES","A1FREQ_CONTROLS")]

  colnames(tmp)<-colnames(datfp)
  tmp$p.value<-format(as.numeric(tmp$p.value), scientific=TRUE,digits=2)
  datfp<-rbind(datfp,tmp)

  datfp$OR<-exp(datfp$BETA)
  datfp$low_ci<-exp(datfp$BETA-1.96*(datfp$SE))
  datfp$up_ci<-exp(datfp$BETA+1.96*(datfp$SE))

  # create truncated lines:
  datfp$low_ci[which(datfp$low_ci<0.5)]<-0.5
  datfp$up_ci[which(datfp$up_ci>2)]<-2

  datfp<-datfp[order(datfp$Neff,decreasing=T),]

  datfp<-datfp[order(datfp$Neff,decreasing=T),]
  datfp$tier<-""
  datfp$tier[which(datfp$study %in% studies)]<-"tier_1"
  datfp$tier[which(datfp$study %in% studies_tier_2)]<-"tier_2"

  # datfp$study<-factor(datfp$study,levels=c(studies,"Meta-analysis eur_tier_2","Meta-analysis eur_tier_1"))
  datfp$study<-factor(datfp$study,levels=c(datfp$study))

  datfp$INFO<-format(round(datfp$INFO,4),nsmall=4)

  p1<-ggplot(datfp,aes(y = study)) + 
    theme_classic() +
    geom_point(aes(x=OR), shape=15, size=3) +
    geom_linerange(aes(xmin=low_ci, xmax=up_ci)) + xlim(0.5,2) +
    geom_vline(xintercept = 1, linetype="dashed",color="grey") 


  p2<-ggplot(datfp) + geom_text(aes(x = 0.5, y = study, label = p.value),
      hjust = 0) + theme_void() + theme(plot.margin = margin(r = -0.5, l=-0.5))


  p3<-ggplot(datfp) + geom_text(aes(x = 0.5, y = study, label = Neff),
      hjust = 0) + theme_void() + theme(plot.margin = margin(l = -0.5))

  p4<-ggplot(datfp) + geom_text(aes(x = 0.5, y = study, label = INFO),
      hjust = 0) + theme_void() + theme(plot.margin = margin(l = -0.5))

  p5<-ggplot(datfp) + geom_text(aes(x = 0.5, y = study, label = tier,color=tier),
      hjust = 0) + theme_void() + theme(plot.margin = margin(l = -0.5)) + 
      scale_color_manual(guide="none",values = c( "red","black","blue"))

  p<-ggarrange(p1,p2,p3,p4,p5,ncol=5,widths=c(2.8,0.5,0.3,0.3,0.3),align = "h")

  p<-annotate_figure(p, top = text_grob(paste(dat$MarkerName[j]," - ",dat$pheno[j]), 
               color = "black", face = "bold", size = 14)) 

  ggsave(paste0("~/git/IIBDGC_GWAS/plots/forrest_plots/forrest_plot_",dat$pheno[j],"_",dat$MarkerName[j],"_huang_cs_lead.pdf",sep=""),
      p,
      width = 15,
      height = 5,
      dpi = 200,
      units = c("in"),
      limitsize = T
  )

  rm(datfp,tmp)

}

q("no")


# compress into one tar file:cd 

# cd ~/git/IIBDGC_GWAS/plots/
# tar -czvf forrest_plot_lead_variant_huang_2017.tar.gz forrest_plots/



