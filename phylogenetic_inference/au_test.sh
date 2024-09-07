#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# iqtree2 needs to be included in $PATH (v2.2.0; http://www.iqtree.org/)

## Command-line args:
nt=$1
alignment=$2
trees=$3
out_file=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### au_test.sh: Starting script."
echo -e "#### au_test.sh: Number of threads: $nt"
echo -e "#### au_test.sh: Input alignment: $alignment"
echo -e "#### au_test.sh: File with trees to compare: $trees"
echo -e "#### au_test.sh: Output file: $out_file \n\n"

################################################################################
#### APPROXIMATELY UNBIASED TEST ####
################################################################################
echo -e "#### au_test.sh: Conducting approximately unbiased test ...\n"

iqtree2 -T $nt -s $alignment -m GTR+G --seqtype DNA -n 0 -z $trees -zb 10000 -au --prefix $out_file

## Report:
echo -e "\n#### au_test.sh: Done with script."
date




