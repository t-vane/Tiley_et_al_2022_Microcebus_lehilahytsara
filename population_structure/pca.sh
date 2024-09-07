#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
pcangsd=/home/nibtve93/software/pcangsd/pcangsd.py # (v1.01; https://github.com/Rosemeis/pcangsd)

## Command-line args:
nt=$1
beagle=$2
out_dir=$3
scripts_dir=$4
ind_file=$5
set_id=$6

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### pca.sh: Starting script."
echo -e "#### pca.sh: Number of threads: $nt"
echo -e "#### pca.sh: Genotype likelihood file in beagle format: $beagle"
echo -e "#### pca.sh: Output directory: $out_dir"
echo -e "#### pca.sh: Directory with scripts: $scripts_dir"
echo -e "#### pca.sh: File that maps individuals to populations: $ind_file"
echo -e "#### pca.sh: Set ID: $set_id \n\n"

################################################################################
#### CONDUCT PRINCIPAL COMPONENT ANALYSIS ####
################################################################################
echo -e "#### pca.sh: Estimating covariance matrix ...\n"
python $pcangsd -threads $nt -beagle $beagle -out $out_dir/cov_matrix.txt -tree

echo -e "#### pca.sh: Obtaining principal components and plotting ...\n"
Rscript $scripts_dir/pca.R $out_dir $out_dir/cov_matrix.txt $ind_file $set_id

echo -e "\n#### pca.sh: Done with script."
date

