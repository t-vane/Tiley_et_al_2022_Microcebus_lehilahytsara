################################################################################
#### READ TRIMMING ####
################################################################################
scripts_dir=/home/nibtve93/scripts/readTrimming

nt=6

set_id=lehilahytsara
in_dir=$PWORK/rawReads/$set_id
out_dir=$PWORK/trimmedReads/$set_id
mkdir $out_dir/logFiles

## Submit read trimming script for paired- and single-end samples (PE and SE, respectively)
for i in PE SE
do
	in_file=$in_dir/trim_$i.txt # List of samples for which reads shall be trimmed
	no_inds=$(cat $in_file | wc -l)
	adapters=$PWORK/trimmomaticAdapters_$i.fa # FASTA file with adapter sequences
	
	sbatch --array=1-$no_inds --output=$out_dir/logFiles/read_trimming_$i.%A_%a.$set_id.oe $scripts_dir/read_trimming.sh $i $in_file $nt $in_dir $out_dir $adapters
done

