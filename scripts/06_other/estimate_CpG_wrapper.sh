#!/bin/bash
# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# Author: https://orcid.org/0000-0002-7256-9752
#

path=/path/to/ibdgwas/IIBDGC/

r_script=$1
chr=$2

r_exec=Rscript

# execute command
${r_exec} ${r_script} ${chr} > ${path}scripts/logs/${chr}_CpG.output

exit $?