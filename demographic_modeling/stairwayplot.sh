#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
stairwayplot=/home/nibtve93/software/stairway_plot_v2.1.1/stairway_plot_es # (v2.0; https://github.com/xiaoming-liu/stairway-plot-v2)

## Command-line args:
blueprint=$1

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### stairwayplot.sh: Starting script."
echo -e "#### stairwayplot.sh: Additional options: $options \n\n"

################################################################################
#### ESTIMATE POPULATION SIZE CHANGES ####
################################################################################
echo -e "#### stairwayplot.sh: Creating shell script ...\n"
java -cp $stairwayplot Stairbuilder $blueprint

echo -e "#### stairwayplot.sh: Running Stairway Plot ...\n"
bash $blueprint.sh 

## Report:
echo -e "\n#### stairwayplot.sh: Done with script."
date
