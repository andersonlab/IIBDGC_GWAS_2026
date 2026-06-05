#!/bin/bash
# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# Author: https://orcid.org/0000-0002-7256-9752
#

# create fam files from final call report style files created by L. Philip Schumm <pschumm@uchicago.edu>

study=$1
path_gwas=$2
array=$3

# execute command
zcat ${path_gwas}pre_imputation/raw/${study}/${study}_${array}_sorted.txt.gz | awk '{if (NR!=1) {print $1}}' | awk '!x[$0]++' | awk -F"\t" -v OFS="\t" '{$2="0" OFS $2; $3="-9" OFS $3; print $1,$1,$2,$2,$3,$3}' > ${path_gwas}pre_imputation/raw/${study}/${study}_${array}.fam

exit $?
  
