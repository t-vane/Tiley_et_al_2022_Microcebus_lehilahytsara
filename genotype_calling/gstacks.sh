#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# gstacks of Stacks needs to be included in $PATH (v2.53; http://catchenlab.life.illinois.edu/stacks/)

## Command-line args:
nt=$1
in_dir=$2
out_dir=$3
popmap=$4
suffix=$5

## Activate conda environment
conda activate stacks

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### gstacks.sh: Starting script."
echo -e "#### gstacks.sh: Number of threads: $nt"
echo -e "#### gstacks.sh: Directory with BAM files: $in_dir"
echo -e "#### gstacks.sh: Output directory: $out_dir"
echo -e "#### gstacks.sh: Population map: $popmap"
echo -e "#### gstacks.sh: Suffix for BAM files: $suffix \n\n"

################################################################################
#### CREATE STACKS WITH GSTACKS ####
################################################################################
echo -e "#### gstacks.sh: Creating stacks with gstacks for individuals in $popmap ...\n"
gstacks -t $nt -I $in_dir -O $out_dir -M $popmap -S $suffix.bam

## Report:
echo -e "\n#### gstacks.sh: Done with script."
date