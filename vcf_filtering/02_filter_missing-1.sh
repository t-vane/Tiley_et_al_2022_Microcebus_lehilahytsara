#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)
# BCFtools needs to be included in $PATH (v1.11; http://www.htslib.org/)

## Command-line args:
vcf_in=$1
vcf_out=$2
maxmiss_geno1=$3
maxmiss_geno2=$4
maxmiss_geno3=$5
filter_inds=$6
maxmiss_ind1=$7
maxmiss_ind2=$8
maxmiss_ind3=$9

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 02_filter_missing-1.sh: Starting script."
echo -e "#### 02_filter_missing-1.sh: Input VCF: $vcf_in"
echo -e "#### 02_filter_missing-1.sh: Output VCF: $vcf_out"
echo -e "#### 02_filter_missing-1.sh: Maximum missingness value for genotypes (round 1): $maxmiss_geno1"
echo -e "#### 02_filter_missing-1.sh: Maximum missingness value for genotypes (round 2): $maxmiss_geno2"
echo -e "#### 02_filter_missing-1.sh: Maximum missingness value for genotypes (round 3): $maxmiss_geno3"
echo -e "#### 02_filter_missing-1.sh: NOTE: Maximum missingness is inverted for genotypes, with a value of 1 indicating no missing data (this is how the '--max-missing' option in VCFtools works)."
echo -e "#### 02_filter_missing-1.sh: Filter individuals: $filter_inds"
[[ $filter_inds ]] && echo -e "#### 02_filter_missing-1.sh: Maximum missingness value for individuals (round 1): $maxmiss_ind1"
[[ $filter_inds ]] && echo -e "#### 02_filter_missing-1.sh: Maximum missingness value for individuals (round 2): $maxmiss_ind2"
[[ $filter_inds ]] && echo -e "#### 02_filter_missing-1.sh: Maximum missingness value for individuals (round 3): $maxmiss_ind3 \n\n"

################################################################################
#### FILTER FOR MISSINGNESS PER GENOTYPE (AND INDIVIDUAL) ####
################################################################################
## Estimate number of indidviduals in input
echo -e "#### 02_filter_missing-1.sh: Estimating number of individuals in $vcf_in ...\n"
nind_in=$(bcftools query -l $vcf_in | wc -l || true)

## Create temporary filter files
mkdir -p $(dirname $vcf_out)/tmp
filterfile_prefix=$(dirname $vcf_out)/tmp/filterfile.$RANDOM
vcf_tmp_prefix=$(dirname $vcf_out)/tmp/vcf.$RANDOM

## Run filtering
if [[ $filter_inds ]]
then
	echo -e "#### 02_filter_missing-1.sh: Filtering by missing data at genotype and individual level in three rounds ...\n"
	for i in 1 2 3 
	do
		[[ $i == 1 ]] && maxmiss_geno=$maxmiss_geno1 && maxmiss_ind=$maxmiss_ind1
		[[ $i == 2 ]] && maxmiss_geno=$maxmiss_geno2 && maxmiss_ind=$maxmiss_ind2
		[[ $i == 3 ]] && maxmiss_geno=$maxmiss_geno3 && maxmiss_ind=$maxmiss_ind3
		
		echo -e "## 02_filter_missing-1.sh: Filtering genotypes by missing data - round $i ..."
		vcftools --vcf $vcf_in --max-non-ref-af 0.99 --min-alleles 2 --max-missing $maxmiss_geno --recode --recode-INFO-all --stdout > $vcf_tmp_prefix.R${i}a.vcf
		
		echo -e "## 02_filter_missing-1.sh: Filtering individuals by missing data - round $i ..."
		# Get amount of missing data per individual
		vcftools --vcf $vcf_tmp_prefix.R${i}a.vcf --missing-indv --stdout > $filterfile_prefix.round$i.imiss
		# Get list of individuals with too much missing data
		tail -n +2 $filterfile_prefix.round$i.imiss | awk -v var=$maxmiss_ind '$5 > var' | cut -f1 > $filterfile_prefix.HiMissInds$i
		# Remove individuals with too much missing data
		vcftools --vcf $vcf_tmp_prefix.R${i}a.vcf --remove $filterfile_prefix.HiMissInds$i --recode --recode-INFO-all --stdout > $vcf_tmp_prefix.R${i}b.vcf
		[[ $i == 3 ]] && mv $vcf_tmp_prefix.R${i}b.vcf $vcf_out		
	done
	
	##Report:
	nvar_in=$(grep -cv "^#" $vcf_in || true)
	nvar_out=$(grep -cv "^#" $vcf_out || true)
	
	nvar_r1=$(grep -cv "^#" $vcf_tmp_prefix.R1a.vcf || true)
	nvar_r2=$(grep -cv "^#" $vcf_tmp_prefix.R2a.vcf || true)
	nvar_r3=$(grep -cv "^#" $vcf_tmp_prefix.R3a.vcf || true)
	
	nvar_filt_r1=$(( $nvar_in - $nvar_r1 ))
	nvar_filt_r2=$(( $nvar_in - $nvar_r2 ))
	nvar_filt_r3=$(( $nvar_in - $nvar_r3 ))
	
	nind_filt_r1=$(wc -l < $filterfile_prefix.HiMissInds1 || true)
	nind_filt_r2=$(wc -l < $filterfile_prefix.HiMissInds2 || true)
	nind_filt_r3=$(wc -l < $filterfile_prefix.HiMissInds3 || true)
	
	nind_out=$(bcftools query -l $vcf_out | wc -l || true)
	
	echo -e "\n#### 02_filter_missing-1.sh: Number of indidviduals before individual filtering: $nind_in"
	echo -e "#### 02_filter_missing-1.sh: Number of indidviduals filtered in round 1: $nind_filt_r1"
	echo -e "#### 02_filter_missing-1.sh: Number of indidviduals filtered in round 2: $nind_filt_r2"
	echo -e "#### 02_filter_missing-1.sh: Number of indidviduals filtered in round 3: $nind_filt_r3"
	echo -e "#### 02_filter_missing-1.sh: Number of indidviduals left after individual filtering: $nind_out \n"
	
	echo -e "#### 02_filter_missing-1.sh: Number of SNPs before genotype filtering: $nvar_in"
	echo -e "#### 02_filter_missing-1.sh: Number of SNPs filtered in round 1: $nvar_filt_r1"
	echo -e "#### 02_filter_missing-1.sh: Number of SNPs filtered in round 2: $nvar_filt_r2"
	echo -e "#### 02_filter_missing-1.sh: Number of SNPs filtered in round 3: $nvar_filt_r3"
	echo -e "#### 02_filter_missing-1.sh: Number of SNPs left after genotype filtering: $nvar_out"	
else
	echo -e "#### 02_filter_missing-1.sh: Only filtering by missing data at genotype level in a single round (no individuals will be removed) ..."
	vcftools --vcf $vcf_in --max-non-ref-af 0.99 --min-alleles 2 --max-missing $maxmiss_geno3 --recode --recode-INFO-all --stdout > $vcf_out
	
	## Report:
	nvar_in=$(grep -cv "^#" $vcf_in || true)
	nvar_out=$(grep -cv "^#" $vcf_out || true)
	nvar_filt=$(( $nvar_in - $nvar_out ))
	
	echo -e "\n#### 02_filter_missing-1.sh: Number of SNPs before genotype filtering: $nvar_in"
	echo -e "#### 02_filter_missing-1.sh: Number of SNPs after genotype filtering: $nvar_out"
	echo -e "#### 02_filter_missing-1.sh: Number of SNPs filtered: $nvar_filt"
fi

## Remove temporary files
rm $filterfile_prefix*
rm $vcf_tmp_prefix*

## Report:
echo -e "\n#### 02_filter_missing-1.sh: Listing output VCF:"
ls -lh $vcf_out
[[ $(grep -cv "^#" $vcf_out) = 0 ]] && echo -e "\n\n#### 02_filter_missing-1.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 02_filter_missing-1.sh: Done with script."
date
