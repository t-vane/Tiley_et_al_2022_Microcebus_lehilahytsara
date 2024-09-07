#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
trimmomatic=/home/nibtve93/software/Trimmomatic-0.39/trimmomatic-0.39.jar (v0.39; https://github.com/usadellab/Trimmomatic)

## Command-line args:
mode=$1
in_file=$2
nt=$3
in_dir=$4
out_dir=$5
adapters=$6
indv=$(sed -n "$SLURM_ARRAY_TASK_ID"p $in_file)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### read_trimming.sh: Starting script."
echo -e "#### read_trimming.sh: Paired-end or single-end: $mode"
echo -e "#### read_trimming.sh: File with individuals: $in_file"
echo -e "#### read_trimming.sh: Number of threads: $nt"
echo -e "#### read_trimming.sh: Directory with raw reads: $in_dir"
echo -e "#### read_trimming.sh: Directory for trimmed reads: $out_dir"
echo -e "#### read_trimming.sh: Adapter file: $adapters"
echo -e "#### reference_mapping.sh: Individual: $indv \n\n"

################################################################################
#### TRIM RAW READS ####
################################################################################
if [[ $mode == "PE" ]]
then
	echo -e "#### read_trimming.sh: Trimming paired-end reads for individual $indv ...\n"
	java -jar $trimmomatic PE -threads $nt -phred33 $in_dir/$indv.1.fq.gz $in_dir/$indv.2.fq.gz $out_dir/$indv.trimmed.1.fq.gz $out_dir/$indv.1.rem.fq.gz $out_dir/$indv.trimmed.2.fq.gz $out_dir/$indv.2.rem.fq.gz \
		ILLUMINACLIP:$adapters:2:30:10 AVGQUAL:20 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:60 
elif [[ $mode == "SE" ]]
then
	echo -e "#### read_trimming.sh: Trimming single-end reads for individual $indv ...\n"
	java -jar $trimmomatic SE -threads $nt -phred33 $in_dir/$indv.1.fq.gz $out_dir/$indv.trimmed.1.fq.gz ILLUMINACLIP:$adapters:2:30:10 AVGQUAL:20 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:60
else
	echo -e "#### read_trimming.sh: Invalid sequencing mode provided - only PE and SE allowed. ...\n" && exit 1
fi

## Report:
echo -e "\n#### read_trimming.sh: Done with script."
date
