#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:
scripts_dir=$1
out_dir=$2
like_values=$3
ind_file=$4
set_id=$5

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### plot_ngsadmix.sh: Starting script."
echo -e "#### plot_ngsadmix.sh: Directory with scripts: $scripts_dir"
echo -e "#### plot_ngsadmix.sh: Output directory: $out_dir"
echo -e "#### plot_ngsadmix.sh: File with likelihood values: $like_values"
echo -e "#### plot_ngsadmix.sh: File that maps individuals to populations: $ind_file"
echo -e "#### plot_ngsadmix.sh: Set ID: $set_id \n\n"

################################################################################
#### PLOT ADMIXTURE RESULTS ####
################################################################################
echo -e "#### plot_ngsadmix.sh: Running script to plot admixture results ...\n"
Rscript $scripts_dir/plot_ngsadmix.R $out_dir $like_values $ind_file $set_id

echo -e "\n#### plot_ngsadmix.sh: Done with script."
date

