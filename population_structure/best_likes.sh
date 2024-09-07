#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:
log=$1
like_file=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### print_likes.sh: Starting script."
echo -e "#### print_likes.sh: Respective log file: $log"
echo -e "#### print_likes.sh: Likelihoods summary file: $like_file \n\n"

################################################################################
#### EXTRACT LIKELIHOOD ####
################################################################################
echo -e "#### print_likes.sh: Extracting likelihood ...\n"
grep "best" $log | awk '{print $K}' | cut -d'=' -f2- | sort -g | sed "s/after/$K/g" | sed "s/iterations/$SEED/g" >> $like_file

echo -e "#### print_likes.sh: Removing fopt.gz files ...\n"
rm $(dirname $log)/$(basename $log .log).fopt.gz

echo -e "\n#### print_likes.sh: Done with script."
date

