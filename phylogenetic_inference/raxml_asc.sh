#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# RAxML-NG needs to be included in $PATH (v1.0.2; https://github.com/amkozlov/raxml-ng)

## Command-line args:
nt=$1
bootstrap=$2
outgroup=$3
in_file=$4
out_file=$5
command=$6

stamatakis1=$(awk '{print $1}' $in_file.stamatakis)
stamatakis2=$(awk '{print $2}' $in_file.stamatakis)
stamatakis3=$(awk '{print $3}' $in_file.stamatakis)
stamatakis4=$(awk '{print $4}' $in_file.stamatakis)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### raxml_asc.sh: Starting script."
echo -e "#### raxml_asc.sh: Number of threads: $nt"
echo -e "#### raxml_asc.sh: Number of bootstrap replicates: $bootstrap"
echo -e "#### raxml_asc.sh: Outgroup: $outgroup"
echo -e "#### raxml_asc.sh: Input alignment file: $in_file"
echo -e "#### raxml_asc.sh: Output prefix: $out_file"
echo -e "#### raxml_asc.sh: Additional commmands: $command \n\n"

################################################################################
#### INFER ML PHYLOGENY WITH ASCERTAINMENT BIAS CORRECTION IN RAXML-NG ####
################################################################################
echo -e "#### raxml_asc.sh: Maximum likelihood phylogenetic inference ...\n"
raxml-ng --all --threads $nt --bs-trees $bootstrap --model GTR+G+ASC_STAM{$stamatakis1/$stamatakis2/$stamatakis3/$stamatakis4} --outgroup $outgroup --seed 12345 --msa $in_file \
	--msa-format PHYLIP --data-type DNA --prefix $out_file $command

## Report:
echo -e "\n#### raxml_asc.sh: Done with script."
date

