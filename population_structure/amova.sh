#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:
scripts_dir=$1
vcf_in=$2
sample_file=$3
out=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### amova.sh: Starting script."
echo -e "#### amova.sh: Directory with scripts: $scripts_dir"
echo -e "#### amova.sh: Input VCF file: $vcf_in"
echo -e "#### amova.sh: File with individuals and associated populations: $sample_file"
echo -e "#### amova.sh: Output file: $out \n\n"

################################################################################
#### CONDUCT AMOVA ####
################################################################################
echo -e "#### amova.sh: Conducting AMOVA ...\n"
Rscript $scripts_dir/amova.R $vcf_in $sample_file $out

## Report:
echo -e "\n#### amova.sh: Done with script."
date

