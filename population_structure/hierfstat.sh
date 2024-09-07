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
out_file=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### hierfstat.sh: Starting script."
echo -e "#### hierfstat.sh: Directory with scripts: $scripts_dir"
echo -e "#### hierfstat.sh: Input VCF file: $vcf_in"
echo -e "#### hierfstat.sh: File with individuals and associated populations: $sample_file"
echo -e "#### hierfstat.sh: Output file: $out_file \n\n"

################################################################################
#### ESTIMATE PAIRWISE F_ST BETWEEN POPULATIONS ####
################################################################################
echo -e "#### hierfstat.sh: Estimating pairwise F_ST between populations ...\n"
Rscript $scripts_dir/hierfstat.R $vcf_in $sample_file $out_file

## Report:
echo -e "\n#### hierfstat.sh: Done with script."
date

