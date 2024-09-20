################################################################################
#### GENOTYPE LIKELIHOOD INFERENCE ####
################################################################################
scripts_dir=/home/nibtve93/scripts/genotypeLikelihoods

set_id=lehilahytsara
bam_dir=$PWORK/bamFiles/$set_id
reference_dir=$PWORK/references/mmur3
reference=$reference_dir/GCF_000165445.2_Mmur_3.0_genomic.fna # Reference genome in fasta format
angsd_dir=$PWORK/$set_id/angsd
nt=80

mkdir -p $angsd_dir/logFiles


#################################################################
#### 1 ESTIMATE GENOTYPE LIKELIHOODS ACROSS ALL SAMPLES ####
#################################################################
in_file_all=$angsd_dir/angsd_all.txt # List of all samples (without file extensions) that shall be included for global genotype likelihood estimation
no_inds_all=$(cat $in_file_all | wc -l)

mkdir -p $angsd_dir/bamHits

## Create bamHits file with locus coverages for each individual
for i in PE SE
do
	in_file=$angsd_dir/angsd_$i.txt # List of samples (without file extensions) that shall be included for global genotype likelihood estimation, separated by sequencing mode (PE or SE)
	no_inds=$(cat $in_file | wc -l)

	sbatch --wait --array=1-$no_inds --output=$angsd_dir/logFiles/bamHits.$set_id.%A_%a.oe $scripts_dir/coverage.sh $i $in_file $bam_dir $angsd_dir/bamHits
done

## Estimate and plot coverage distributions for each individual
sbatch --wait --output=$angsd_dir/logFiles/cov_plot.$set_id.oe $scripts_dir/cov_plot.sh $scripts_dir $in_file_all $angsd_dir/bamHits $set_id

## Set thresholds for angsd
mindepthind=$(cat $angsd_dir/bamHits/statistics/$set_id.minmax.txt | cut -d " " -f2 | sort -n | head -1) # Minimum depth per individual
maxdepthind=$(cat $angsd_dir/bamHits/statistics/$set_id.minmax.txt | cut -d " " -f3 | sort -n | tail -1) # Maximum depth per individual
gmin=$(cat $angsd_dir/bamHits/statistics/$set_id.minmax.txt | cut -d " " -f2 | paste -sd+ | bc) # Minimum depth across all individuals
gmax=$(cat $angsd_dir/bamHits/statistics/$set_id.minmax.txt | cut -d " " -f3 | paste -sd+ | bc) # Maximum depth across all individuals
percentage="75/100" # Minimum percentage of represented individuals
minind=$(($NO_INDS_ALL * $percentage )) # Minimum number of represented individuals

## Create BAM list as input to angsd
rm $angsd_dir/$set_id.bamlist
while read indv
do
echo $bam_dir/$indv.auto.bam >> $angsd_dir/$set_id.bamlist
done < $in_file_all

## Run angsd
filters="-setMinDepth $gmin -setMaxDepth $gmax -setMaxDepthInd $maxdepthind -setMinDepthInd $mindepthind -minInd $minind -SNP_pval 1e-5 -minQ 20 -minMapQ 20 -minMaf 0.05 -uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -only_proper_pairs 1 -baq 1 -C 50"
todo="-GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1"
sbatch --output=$angsd_dir/logFiles/angsd.$set_id.oe $scripts_dir/angsd.sh $nt $reference $angsd_dir/$set_id.bamlist "$todo" "$filters" $angsd_dir/$set_id


#################################################################
#### 2 ESTIMATE SITE ALLELE FREQUENCY LIKELIHOODS PER POPULATION ####
#################################################################
pops="riamalandy ambavala ambatovy ambohitantely anjanaharibe anjiahely ankafobe marojejy tsinjoarivo metapopulation"

## Estimate site allele frequency likelihoods per population
for i in $pops
do
	echo -e "#### Processing population $i ...\n"
	in_file=$angsd_dir/angsd_$i.txt # List of samples in the population
	no_inds=$(cat $in_file | wc -l)

	## Estimate and plot coverage distributions for each individual
	sbatch --wait --output=$angsd_dir/logFiles/cov_plot.$i.oe $scripts_dir/cov_plot.sh $scripts_dir $in_file $angsd_dir/bamHits $i

	## Set thresholds for angsd
	mindepthind=$(cat $angsd_dir/bamHits/statistics/$i.minmax.txt | cut -d " " -f2 | sort -n | head -1) # Minimum depth per individual
	maxdepthind=$(cat $angsd_dir/bamHits/statistics/$i.minmax.txt | cut -d " " -f3 | sort -n | tail -1) # Maximum depth per individual
	gmin=$(cat $angsd_dir/bamHits/statistics/$i.minmax.txt | cut -d " " -f2 | paste -sd+ | bc) # Minimum depth across all individuals
	gmax=$(cat $angsd_dir/bamHits/statistics/$i.minmax.txt | cut -d " " -f3 | paste -sd+ | bc) # Maximum depth across all individuals
	percentage="75/100" # Minimum percentage of represented individuals
	minind=$(( $no_inds * $percentage )) # Minimum number of represented individuals

	## Create BAM list as input to angsd
	rm $angsd_dir/$i.bamlist
	while read indv
	do
	echo $bam_dir/$indv.auto.bam >> $angsd_dir/$i.bamlist
	done < $in_file

	## Run angsd
	todo="-GL 1 -doSaf 1 -doCounts 1 -anc $reference"
	filters="-setMinDepth $gmin -setMaxDepth $gmax -setMaxDepthInd $maxdepthind -setMinDepthInd $mindepthind -minInd $minind -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 -C 50 -baq 1 -only_proper_pairs 1"
	sbatch --output=$angsd_dir/logFiles/angsd.saf.$i.oe $scripts_dir/angsd.sh $nt $reference $angsd_dir/$i.bamlist "$todo" "$filters" $angsd_dir/$i
done
