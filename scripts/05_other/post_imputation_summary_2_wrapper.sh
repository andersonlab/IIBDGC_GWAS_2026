#!/bin/bash
# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# Author: https://orcid.org/0000-0002-7256-9752
#

study=$1
r_script=$2

r_exec=/software/R-4.3.1/bin/Rscript

# create R output
output_r_script=/path/to/ibdgwas/IIBDGC/scripts/logs/post_imputation_summary_2.R_output_${study}_allchr

# execute command
${r_exec} ${r_script} ${study} > ${output_r_script}

exit $?