# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# # how to submit - see step_40
# # singularity exec iibdgc_postprocess_10_singularity.sif

# MEM=8000
# release="eur_eas_sas_tier_2"

# bsub -J"plot_imp1" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -m "modern_hardware" -G ibdgwas \
# -e ${path_gwas}scripts/logs/draw_manhattan_plots_stderr_${release} \
# -o ${path_gwas}scripts/logs/draw_manhattan_plots_stdout_${release} \
# "Rscript ~/git/IIBDGC_GWAS/scripts/other/draw_manhattan_plots_label_known_regions.R ${release} > \
# ${path_gwas}scripts/logs/draw_manhattan_plots_label_known_regions_${release}.Rout"

# for testing purposes
# MEM=6000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(ggplot2)
library(colorspace)
library(dplyr)
library(ggrepel)
library(ggpubr)

path<-"/path/to/ibdgwas/IIBDGC/"
path_2<-"~/git/IIBDGC_GWAS/"

pheno<-c("ibd","cd","uc")
cols<-c("#db7107","#004488","#BB5566")
col_ligth<-lighten(cols, 0.70, space = "HCL")

args <- commandArgs()
release<-args[6]

# options for release are: 
#   eur_tier_1
#   eur_tier_2
#   sas_tier_2
#   eur_eas_sas_tier_2

# define the genomic regions including the signals
all<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_conditional_on_known_index_gwas_fm_2023_with_iibdgc_tier1_tier2_eas_sas_summary_stats_plus_unsupervised_conditional_plus_new_fm_plus_pheno.tsv",sep=""),head=T)

# create a list of SNPs per region:
regions<-all[,c("chr","updated_region","min","max","class")]
regions<-regions[which(regions$updated_region!=""),]
regions<-regions[!duplicated(regions),]

table(regions$class)
# (old)   new old 
#    30   138 283


for (i in 2:length(pheno)){
  
  print(pheno[i])
  
  for (chr in c(1:22,"X")) {
    
    if (release=="eur_tier_1") {
        tmp<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),)
        tmp<-tmp[which(tmp$N>=(0.8*max(tmp$N,na.rm=T))),]
    } else {
        tmp<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_",release,"_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),)
        tmp<-tmp[which(tmp$N>=(0.4*max(tmp$N,na.rm=T))),]
    } 

    tmp$'P-value'<-as.numeric(tmp$'P-value')
    tmp<-tmp[which(tmp$'P-value'<0.05),]

    tmp$HetPVal<-as.numeric(tmp$HetPVal)

    if(chr=="X") {
      tmp$pos<-gsub("chrX:","",tmp$MarkerName)
      tmp$pos<-as.numeric(gsub(":[A-Z]*:[A-Z]*","",tmp$pos))
      # keep only PAR:
      tmp<-tmp[which(tmp$pos<=2781479 | tmp$pos>=155701383),] 
    }else{
       tmp$pos<-gsub("chr[0-9]{1,2}:","",tmp$MarkerName)
       tmp$pos<-as.numeric(gsub(":[A-Z]*:[A-Z]*","",tmp$pos))
    }
      
    tmp$chr<-chr

    if(chr=="X") {
        tmp$chr<-23
    }

    # plot only signals with no significant het in effect size:
    if (release=="sas_tier_2") {
        tmp<-tmp[which(tmp$HetPVal>=1E-3),]
    } else {
        tmp<-tmp[which(tmp$HetPVal>=1E-10),]
    }

    tmp<-tmp[,c("chr","MarkerName","pos","P-value")]
    
    if(chr==1) {
      tmp1<-tmp
    } else {
      tmp1<-rbind(tmp,tmp1)
    }
    rm(tmp)
    
  }

  tmp1$chr<-as.numeric(tmp1$chr)

  don <- tmp1 %>% 
  
    # Compute chromosome size
    group_by(chr) %>% 
    summarise(chr_len=max(pos)) %>% 
  
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
  
    # Add this info to the initial dataset
    left_join(tmp1, ., by=c("chr"="chr")) %>%
  
    # Add a cumulative position of each SNP
    arrange(chr, pos) %>%
    mutate( poscum=pos+tot)

    axisdf = don %>%
    group_by(chr) %>%
    summarize(center=( max(poscum) + min(poscum) ) / 2 )

  rm(tmp1)

  don$P<-(-log10(don$'P-value'))

  don$is_highlight<-"no"
  for (jj in 1:nrow(regions)) {

    if (regions$class[jj] %in% c("old","(old)")) {
        don$is_highlight[which(don$chr==regions$chr[jj] & don$pos>=regions$min[jj]-250000 & don$pos<=regions$max[jj]+250000)]<-"yes_old"
    }

    if (regions$class[jj]=="new") {
        don$is_highlight[which(don$chr==regions$chr[jj] & don$pos>=regions$min[jj] & don$pos<=regions$max[jj] & don$'P-value'<5E-8)]<-"yes_new"
    }

  }

  # label anything in or nearby HLA as known:
  don$is_highlight[which(don$chr==6 & don$pos>=28510120-500000 & don$pos<=33480577+500000 & don$'P-value'<5E-8)]<-"yes_old"
  don$P[which(don$P=="Inf")]<-330

  # plot all results - not truncated:
  p1<-ggplot(don, aes(x=poscum, y=P)) +
    
    # Show all points
    geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c(col_ligth[i],"gray60"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0),limits=c(0,340)) +     # remove space between plot area and x axis
  
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) + xlab("Chromosome") + ylab(bquote(-log[10](Pval)))

  # plot all results - excluding known regions:
  p2<-ggplot(don[which(don$is_highlight!="yes_old")], aes(x=poscum, y=P)) +
    
    # Show all points
    geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c(col_ligth[i],"gray60"), 22 )) +

    geom_point(data=subset(don, is_highlight=="yes_new"), color=cols[i]) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0),limits=c(1,20)) +     # remove space between plot area and x axis
  
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) + xlab("Chromosome") + ylab(bquote(-log[10](Pval)))

  if (i==1) {
    vec<-c("a", "b")
  } else if (i==2) {
    vec<-c("c", "d")
  } else {
    vec<-c("e", "f")
  }

  assign(pheno[i],ggarrange(p1,p2,ncol=2,labels = vec))

  # rm(don)

}

# p<-ggarrange(p1,p2)
p<-ggarrange(ibd,cd,uc,nrow=3)

ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Figure_1_manhattan_plots.pdf",
  p,
  width = 180,
  height = 150,
  dpi = 400,
  units = c("mm"),
  scale = 2.7
)


q("no")









