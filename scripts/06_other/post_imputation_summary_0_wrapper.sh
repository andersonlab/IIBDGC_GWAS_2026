#!/bin/bash
# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# Author: https://orcid.org/0000-0002-7256-9752
#

study=$1
chr=$2
r_script=$3

r_exec=/software/R-4.3.1/bin/Rscript

# create R output
output_r_script=/path/to/ibdgwas/IIBDGC/scripts/logs/post_imputation_summary_0.R_output_${study}_${chr}

# execute command
${r_exec} ${r_script} ${study} ${chr} > ${output_r_script}

exit $?