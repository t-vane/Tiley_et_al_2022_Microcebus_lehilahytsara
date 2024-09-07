#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# realSFS needs to be included in $PATH (http://www.popgen.dk/angsd/index.php/RealSFS)

## Command-line args:
nt=$1
saf_in=$2
out_file=$3
options=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### realsfs.sh: Starting script."
echo -e "#### realsfs.sh: Number of threads: $nt"
echo -e "#### realsfs.sh: Site allele frequency file(s) for population(s): $saf_in"
echo -e "#### realsfs.sh: Output file: $out_file"
echo -e "#### realsfs.sh: Additional options: $options \n\n"

################################################################################
#### ESTIMATE MINOR ALLELE FREQUENCY SPECTRUM ####
################################################################################
echo -e "#### realsfs.sh: Estimating minor allele frequency spectrum ...\n"
realSFS $saf_in -P $nt -fold 1 $options > $out_file

## Report:
echo -e "\n#### realsfs.sh: Done with script."
date

