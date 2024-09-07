#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args
scripts_dir=$1
geo_dist=$2
gen_dist=$3
out=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### ibd.sh: Starting script."
echo -e "#### ibd.sh: Directory with scripts: $scripts_dir"
echo -e "#### ibd.sh: Geographic distance matrix: $geo_dist"
echo -e "#### ibd.sh: Genetic distance matrix: $gen_dist"
echo -e "#### ibd.sh: Output prefix: $out \n\n"

################################################################################
#### CONDUCT MANTEL TEST AND PLOT ISOLATION BY DISTANCE ####
################################################################################
echo -e "#### ibd.sh: Conducting Mantel test and plotting ...\n"
Rscript $scripts_dir/ibd.R $geo_dist $gen_dist $out

## Report:
echo -e "\n#### ibd.sh: Done with script."
date

