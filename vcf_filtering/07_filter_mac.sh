#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)

## Command-line args:
vcf_in=$1
vcf_out=$2
mac=$3

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 07_filter_mac.sh: Starting script."
echo -e "#### 07_filter_mac.sh: Input VCF: $vcf_in"
echo -e "#### 07_filter_mac.sh: Output VCF: $vcf_out"
echo -e "#### 07_filter_mac.sh: Minor allele count: $mac \n\n"

################################################################################
#### FILTER WITH MINOR ALLELE COUNT ####
################################################################################
echo -e "#### 07_filter_mac.sh: Filtering for minor allele count $mac ...\n"
vcftools --vcf $vcf_in --mac $mac --recode --recode-INFO-all --stdout > $vcf_out

## Report:
nvar_in=$(grep -cv "^#" $vcf_in || true)
nvar_out=$(grep -cv "^#" $vcf_out || true)
nvar_filt=$(( $nvar_in - $nvar_out ))

echo -e "\n#### 07_filter_mac.sh: Number of SNPs before MAC filtering: $nvar_in"
echo -e "#### 07_filter_mac.sh: Number of SNPs filtered: $nvar_filt"
echo -e "#### 07_filter_mac.sh: Number of SNPs after MAC filtering: $nvar_out"

echo -e "\n#### 07_filter_mac.sh: Listing output VCF:"
ls -lh $vcf_out
[[ $(grep -cv "^#" $vcf_out) = 0 ]] && echo -e "\n\n#### 07_filter_mac.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 07_filter_mac.sh: Done with script."
date

