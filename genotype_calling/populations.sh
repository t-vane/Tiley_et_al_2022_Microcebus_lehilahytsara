#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# populations of Stacks needs to be included in $PATH (v2.53; http://catchenlab.life.illinois.edu/stacks/)

## Command-line args:
nt=$1
in_dir=$2
pop_dir=$3
popmap=$4
filters=$5
output=$6

## Activate conda environment:
conda activate stacks

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### populations.sh: Starting script."
echo -e "#### populations.sh: Number of threads: $nt"
echo -e "#### populations.sh: Directory with created stacks: $in_dir"
echo -e "#### populations.sh: Output directory: $pop_dir"
echo -e "#### populations.sh: Population map: $popmap"
echo -e "#### populations.sh: Desired filters: $filters"
echo -e "#### populations.sh: Desired outputs: $output \n\n"

################################################################################
#### EXTRACT VCF WITH POPULATIONS ####
################################################################################
echo -e "#### populations.sh: Extracting VCF file from created stacks with populations ...\n"
populations -t $nt -P $in_dir -O $pop_dir -M $popmap $filters $output

## Report:
echo -e "\n#### populations.sh: Done with script."
date







