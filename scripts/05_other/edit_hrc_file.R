# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# retrieve frequency from hrc file

path_gwas=/path/to/ibdgwas/IIBDGC/
  
# cd /path/to/ibdgwas/IIBDGC/resources/HRC
# wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz
  
  
  /path/to/software/bcftools-1.9/./bcftools query -f '%CHROM %POS %REF %ALT %AF\n' \
/path/to/ibdgwas/IIBDGC/resources/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz -o ${path_gwas}resources/HRC/HRC_freq

awk -v OFS='' '{print $1,":",$2,"_",$3,"_",$4}' ${path_gwas}resources/HRC/HRC_freq | sed -e 's/X:/23:/g' -e 's/Y:/24:/g' \
| paste -d'\t' - ${path_gwas}resources/HRC/HRC_freq | sed -e 's/ /\t/g' > ${path_gwas}resources/HRC/HRC_freq_edited

rm ${path_gwas}resources/HRC/HRC_freq

gzip ${path_gwas}resources/HRC/HRC_freq_edited