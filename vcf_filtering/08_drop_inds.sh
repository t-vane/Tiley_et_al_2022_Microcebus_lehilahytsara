#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# vcftools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)
# bcftools needs to be included in $PATH (v1.11; http://www.htslib.org/)

## Command-line args:
vcf_in=$1
drop_inds=$2
vcf_out=$3

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 08_drop_inds.sh: Starting script."
echo -e "#### 08_drop_inds.sh: Input VCF: $vcf_in"
echo -e "#### 08_drop_inds.sh: List of individuals to remove: $drop_inds"
echo -e "#### 08_drop_inds.sh: Output VCF: $vcf_out \n\n"

################################################################################
#### REMOVE UNDESIRED INDIVIDUALS ####
################################################################################
echo -e "#### 08_drop_inds.sh: Removing undesired individuals ...\n"
vcftools --remove $drop_inds --vcf $vcf_in --recode --recode-INFO-all --stdout > $vcf_out

## Report:
nind_in=$(bcftools query -l | wc -l $vcf_in)
nind_out=$(bcftools query -l | wc -l $vcf_out)
nind_filt=$(( $nind_in - $nind_out ))

echo -e "\n#### 08_drop_inds.sh: Number of individuals before filtering: $nind_in"
echo -e "#### 08_drop_inds.sh: Number of individuals filtered: $nind_filt"
echo -e "#### 08_drop_inds.sh: Number of individuals after filtering: $nind_out"

echo -e "\n#### 08_drop_inds.sh: Listing output VCF:"
ls -lh $vcf_out
[[ $(grep -cv "^#" $vcf_out) = 0 ]] && echo -e "\n\n#### 08_drop_inds.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 08_drop_inds.sh: Done with script."
date

