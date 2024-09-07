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
vcf_out=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 05_filter_max-dp.sh: Starting script."
echo -e "#### 05_filter_max-dp.sh: Input VCF: $vcf_in"
echo -e "#### 05_filter_max-dp.sh: Output VCF: $vcf_out \n\n"

################################################################################
#### FILTER FOR MAXIMUM DEPTH ####
################################################################################
## Create a file with original site depths per locus
mkdir -p $(dirname $vcf_out)/tmp
dp_file=$(dirname $vcf_out)/tmp/vcf_depth.$RANDOM
echo -e "#### 05_filter_max-dp.sh: Extracting original site depths per locus ...\n"
vcftools --vcf $vcf_in --site-depth --stdout | tail -n +2 | cut -f 3 > $dp_file

## Calculate the mean maximum depth: (mean depth + 2 * standard deviation) / number of individuals
echo -e "#### 05_filter_max-dp.sh: Estimating mean maximum depth ...\n"
n_ind=$(bcftools query -l $vcf_in | wc -l || true)
dp_mean_all=$(awk '{ total += $1; count++ } END { print total/count }' $dp_file)
dp_sd_all=$(awk '{ sum+=$1; sumsq+=$1*$1 } END { print sqrt(sumsq/NR - (sum/NR)^2) }' $dp_file)
dp_hi_all=$(python -c "print($dp_mean_all + (2* $dp_sd_all))")
dp_mean_ind=$(python -c "print($dp_mean_all / $n_ind)")
dp_max_ind=$(python -c "print($dp_hi_all / $n_ind)")

## Filter with VCFtools for maximum depth
echo -e "#### 05_filter_max-dp.sh: Filtering for maximum depth ...\n"
vcftools --vcf $vcf_in --max-meanDP $dp_max_ind --recode --recode-INFO-all --stdout > $vcf_out

## Save a separate file with SNPs with too high depth
echo -e "#### 05_filter_max-dp.sh: Creating file with SNPs with too high depth ...\n"
vcftools --vcf $vcf_in --min-meanDP $dp_max_ind --recode --recode-INFO-all --stdout > ${vcf_out//.vcf/_too-high-DP.vcf}

## Remove temporary file
rm $dp_file

## Report:
nvar_in=$(grep -cv "^#" $vcf_in)
nvar_out=$(grep -cv "^#" $vcf_out)
nvar_filt=$(( $nvar_in - $nvar_out ))

echo -e "\n\n"
echo -e "#### 05_filter_max-dp.sh: Number of SNPs in input VCF: $nvar_in"
echo -e "#### 05_filter_max-dp.sh: Number of SNPs in output VCF: $nvar_out"
echo -e "#### 05_filter_max-dp.sh: Number of SNPs filtered: $nvar_filt"
echo
echo -e "#### 05_filter_max-dp.sh: Listing output VCF:"
ls -lh $vcf_out
[[ $(grep -cv "^#" $vcf_out) = 0 ]] && echo -e "\n\n#### 05_filter_max-dp.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 05_filter_max-dp.sh: Done with script."
date

