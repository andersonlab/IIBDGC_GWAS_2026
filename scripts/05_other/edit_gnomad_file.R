# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# retrieve frequency from gnomad file

path_gwas=/path/to/ibdgwas/IIBDGC/

# /path/to/project


bcftools query -f '%CHROM %POS %REF %ALT\n' file.bcf | head -3

zcat /path/to/project | grep "##INFO="
##INFO=<ID=AF_nfe,Number=A,Type=Float,Description="Alternate allele frequency in samples of Non-Finnish European ancestry">
##INFO=<ID=AF_afr,Number=A,Type=Float,Description="Alternate allele frequency in samples of African-American/African ancestry">
##INFO=<ID=AF_eas,Number=A,Type=Float,Description="Alternate allele frequency in samples of East Asian ancestry">
##INFO=<ID=AF_sas,Number=A,Type=Float,Description="Alternate allele frequency in samples of South Asian ancestry">
##INFO=<ID=AF_asj,Number=A,Type=Float,Description="Alternate allele frequency in samples of Ashkenazi Jewish ancestry">

##INFO=<ID=AF_ami,Number=A,Type=Float,Description="Alternate allele frequency in samples of Amish ancestry">
##INFO=<ID=AF_oth,Number=A,Type=Float,Description="Alternate allele frequency in samples of Other ancestry">
##INFO=<ID=AF_sas,Number=A,Type=Float,Description="Alternate allele frequency in samples of South Asian ancestry">
##INFO=<ID=AF_fin,Number=A,Type=Float,Description="Alternate allele frequency in samples of Finnish ancestry">
##INFO=<ID=AF_amr,Number=A,Type=Float,Description="Alternate allele frequency in samples of Latino ancestry">



path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=200

bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q long \
-e ${path_gwas}pre_imputation/QC/log/stderr_gnomad \
-o ${path_gwas}pre_imputation/QC/log/stdout_gnomad \
"/path/to/software/bcftools-1.9/./bcftools query -i 'TYPE=\"snp\"' -f '%CHROM %POS %REF %ALT %AF %AF_nfe %AF_afr %AF_amr %AF_eas %AF_sas %AF_asj\n' \
/path/to/project -o ${path_gwas}resources/gnomad/gnomad_freq"
# Job <885660> is submitted to queue <long>.

wc -l ${path_gwas}resources/gnomad/gnomad_freq
# 602821631 /path/to/ibdgwas/IIBDGC/resources/gnomad/gnomad_freq


########

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=1500

bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/log/stderr_gnomad_2 \
-o ${path_gwas}pre_imputation/QC/log/stdout_gnomad_2 \
"awk -v OFS='' '{print \$1,\":\",\$2,\"_\",\$3,\"_\",\$4}' ${path_gwas}resources/gnomad/gnomad_freq | sed 's/chr//g' | sed -e 's/X:/23:/g' -e 's/Y:/24:/g' \
| paste -d'\t' - ${path_gwas}resources/gnomad/gnomad_freq > ${path_gwas}resources/gnomad/gnomad_freq_edited"
# Job <963129> is submitted to queue <normal>.

wc -l ${path_gwas}resources/gnomad/gnomad_freq_edited
# 602821631 /path/to/ibdgwas/IIBDGC/resources/gnomad/gnomad_freq_edited



########

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=1500

bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/log/stderr_gnomad_3 \
-o ${path_gwas}pre_imputation/QC/log/stdout_gnomad_3 \
"sed '1 i\SNP\tCHROM\tPOS\tREF\tALT\tAF\tAF_nfe\tAF_afr\tAF_amr\tAF_eas\tAF_sas\tAF_asj' ${path_gwas}resources/gnomad/gnomad_freq_edited \
> ${path_gwas}resources/gnomad/gnomad_freq_edited2d"
# Job <964240> is submitted to queue <normal>.


wc -l ${path_gwas}resources/gnomad/gnomad_freq_edited2d
# 602821632 /path/to/ibdgwas/IIBDGC/resources/gnomad/gnomad_freq_edited2d

mv ${path_gwas}resources/gnomad/gnomad_freq_edited2d ${path_gwas}resources/gnomad/gnomad_freq_edited
rm ${path_gwas}resources/gnomad/gnomad_freq
ls -lh
# -rw-r--r-- 1 username ibdgwas 40G Jul  3 09:09 gnomad_freq_edited



########################################################################################################################################################################
########################################################################################################################################################################

# create similar resource but for b37

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=200

for chr in {1..22} X
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/log/stderr_gnomad_chr${chr} \
-o ${path_gwas}pre_imputation/QC/log/stdout_gnomad_chr${chr} \
"/path/to/software/bcftools-1.16/./bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' /path/to/project \
| /path/to/software/bcftools-1.16/./bcftools \
query -i 'TYPE=\"snp\"' -f '%ID %CHROM %POS %REF %ALT %AF %AF_nfe %AF_afr %AF_amr %AF_eas %AF_asj\n' | gzip > ${path_gwas}resources/gnomad/gnomad_freq_edited_b37_chr${chr}.tsv.gz"
done

for chr in {1..22} X
do tail -50 ${path_gwas}pre_imputation/QC/log/stdout_gnomad_chr${chr} | grep "Successfully"
done



