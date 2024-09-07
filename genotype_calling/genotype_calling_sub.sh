################################################################################
#### GENOTYPE CALLING WITH STACKS####
################################################################################
scripts_dir=/home/nibtve93/scripts/stacks

nt=16

set_id=lehilahytsara
bam_dir=$PWORK/bamFiles/$set_id
out_dir=$PWORK/$set_id/stacks
popmap=$out_dir/$set_id.popmap # List of samples included in genotyping
suffix=auto	# Suffix for final BAM files (see scripts for reference mapping)

mkdir -p $out_dir/logFiles
mkdir -p $out_dir/populations

## Create stacks with gstacks
sbatch --job-name=stacks --output=$out_dir/logFiles/gstacks.$set_id.oe $scripts_dir/gstacks.sh $nt $bam_dir $out_dir $popmap 

## Extract VCF file with populations
filters="-p 1 -r 0.75"
output="--vcf"
sbatch --job-name=stacks --dependency=singleton --output=$out_dir/logFiles/populations.$set_id.oe $scripts_dir/populations.sh $out_dir $out_dir/populations $popmap $nt "$filters" "$output"
