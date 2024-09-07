#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)

## Command-line args:
vcf_in=$1
min_dp=$2
mean_dp=$3
vcf_out=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 01_filter_min-dp.sh: Starting script."
echo -e "#### 01_filter_min-dp.sh: Input VCF: $vcf_in"
echo -e "#### 01_filter_min-dp.sh: Minimum depth: $min_dp"
echo -e "#### 01_filter_min-dp.sh: Minimum mean depth: $mean_dp"
echo -e "#### 01_filter_min-dp.sh: Output VCF: $vcf_out \n\n"

################################################################################
#### FILTER FOR MINIMUM DEPTHS ####
################################################################################
## Filter with VCFtools for minimum depth
echo -e "#### 01_filter_min-dp.sh: Filtering for minimum depth ...\n"
vcftools --vcf $vcf_in --minDP $min_dp --min-meanDP $mean_dp --recode --recode-INFO-all --stdout > $vcf_out

## Report:
nvar_in=$(grep -cv "^#" $vcf_in)
nvar_out=$(grep -cv "^#" $vcf_out)
nvar_filt=$(( $nvar_in - $nvar_out ))

echo -e "\n\n"
echo -e "#### 01_filter_min-dp.sh: Number of SNPs in input VCF: $nvar_in"
echo -e "#### 01_filter_min-dp.sh: Number of SNPs in output VCF: $nvar_out"
echo -e "#### 01_filter_min-dp.sh: Number of SNPs filtered: $nvar_filt"
echo
echo -e "#### 01_filter_min-dp.sh: Listing output VCF:"
ls -lh $vcf_out
[[ $(grep -cv "^#" $vcf_out) = 0 ]] && echo -e "\n\n#### 01_filter_min-dp.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 01_filter_min-dp.sh: Done with script."
date
