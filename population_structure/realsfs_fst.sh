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
sfs=$3
out=$4
options=$5

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### realsfs_fst.sh: Starting script."
echo -e "#### realsfs_fst.sh: Number of threads: $nt"
echo -e "#### realsfs_fst.sh: Site allele frequency file(s) for population(s): $saf_in"
echo -e "#### realsfs_fst.sh: Joint minor allele frequency spectrum: $sfs"
echo -e "#### realsfs_fst.sh: Output prefix: $out"
echo -e "#### realsfs_fst.sh: Additional options: $options \n\n"

################################################################################
#### ESTIMATE PAIRWISE F_ST BETWEEN POPULATIONS ####
################################################################################
echo -e "#### realsfs_fst.sh: Estimating pairwise F_ST between populations ...\n"
realSFS fst index $saf_in -P $nt -sfs $sfs -fstout $out
realSFS fst stats $out.fst.idx > $out.fst

## Report:
echo -e "\n#### realsfs_fst.sh: Done with script."
date

