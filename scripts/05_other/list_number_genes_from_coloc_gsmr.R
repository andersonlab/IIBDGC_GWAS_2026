# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
#!/bin/bash

study=$1
imputation_panel=$2

path_gwas=/path/to/ibdgwas/IIBDGC/

# execute command

/path/to/user/software/bcftools-1.9/./bcftools \
convert ${path_gwas}imputed/${study}/${study}${imputation_panel}.vcfs/${LSB_JOBINDEX}.vcf.gz \
-g ${path_gwas}imputed/${study}/${study}${imputation_panel}.vcfs/${LSB_JOBINDEX}_imputed_${study}${imputation_panel}

exit $?
