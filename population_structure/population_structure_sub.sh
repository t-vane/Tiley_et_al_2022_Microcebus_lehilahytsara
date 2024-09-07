################################################################################
#### POPULATION STRUCTURE ANALYSES ####
################################################################################
scripts_dir=/home/nibtve93/scripts/populationStructure

set_id=lehilahytsara
pop_dir=$PWORK/$set_id/populationStructure
beagle=$PWORK/$set_id/angsd/$set_id.beagle.gz # Genotype likelihood file created in genotype_likelihoods_sub.sh

mkdir -p $pop_dir/logFiles


#################################################################
#### 1 PRINCIPAL COMPONENT ANALYSIS ####
#################################################################

##################################################
#### 1.1 GENOTYPE CALL BASED ####
##################################################

## See Dryad digital repository (https://doi.org/10.5061/dryad.dncjsxkvz).

##################################################
#### 1.2 GENOTYPE LIKELIHOOD BASED ####
##################################################
mkdir -p $pop_dir/pca/genotypeLikelihoods

ind_file=$pop_dir/$set_id.txt # File with individual IDs in first column and population assignments in second column
nt=20

sbatch --output=$pop_dir/logFiles/pca.$set_id.oe $scripts_dir/pca.sh $nt $beagle $pop_dir/pca/genotypeLikelihoods $scripts_dir $ind_file $set_id


#################################################################
#### 2 ADMIXTURE ####
#################################################################

##################################################
#### 2.1 GENOTYPE CALL BASED ####
##################################################

## See Dryad digital repository (https://doi.org/10.5061/dryad.dncjsxkvz).

##################################################
#### 2.2 GENOTYPE LIKELIHOOD BASED ####
##################################################
mkdir -p $pop_dir/ngsadmix/genotypeLikelihoods

clusters=12 # Maximum number of clusters to assume in admixture analysis
repeats=10 # Number of independent runs
percentage="75/100" # Minimum percentage of represented individuals
minind=$(( ($(zcat $beagle | head -1 | wc -w)/3-1) * $percentage )) # Minimum number of represented individuals
nt=80

## Submit array job to infer individual ancestries for each number of clusters (ranging from 1 to $clusters), using $repeats repetitions 
for k in $(seq 1 $clusters)
do
	jid=$(sbatch --array=1-$repeats --output=$pop_dir/logFiles/ngsadmix$k.$set_id.%A_%a.oe $scripts_dir/ngsadmix.sh $nt $k $beagle $pop_dir/ngsadmix/genotypeLikelihoods $minind $set_id)
	declare RUNID_$k=${jid##* }
done


## Print likelihood values file
like_file=$pop_dir/ngsadmix/genotypeLikelihoods/likevalues.$set_id.txt # File for likelihoods summary
rm $like_file; touch $like_file
for k in $(seq 1 $clusters); 
do
	for seed in $(seq 1 $repeats)
	do
		[[ $k == 1 ]] && [[ $seed == 1 ]] && varname=RUNID_$k && jid=$(sbatch --account=nib00015 --dependency=afterok:${!varname} --output=$pop_dir/logFiles/print_likes.$set_id.oe $scripts_dir/print_likes.sh $pop_dir/ngsadmix/genotypeLikelihoods/$set_id.K$k.seed$seed.log $like_file)
		[[ $k != 1 ]] || [[ $seed != 1 ]] && varname=RUNID_$k && jid=$(sbatch --account=nib00015 --dependency=afterok:${!varname}:${jid##* } --output=$pop_dir/logFiles/print_likes.$set_id.oe $scripts_dir/print_likes.sh $pop_dir/ngsadmix/genotypeLikelihoods/$set_id.K$k.seed$seed.log $like_file)
	done
done

## Plot results
ind_file=$pop_dir/$set_id.txt # File with individual IDs in first columns and population assignments in second column

until [[ $(cat $like_file | wc -l) == $(( $clusters*$repeats )) ]]
do
	sleep 5
done

sbatch --output=$pop_dir/logFiles/plot_ngsadmix.$set_id.oe $scripts_dir/plot_ngsadmix.sh $scripts_dir $pop_dir/ngsadmix/genotypeLikelihoods $like_file $ind_file $set_id


#################################################################
#### 3 ISOLATION-BY-DISTANCE ####
#################################################################

##################################################
#### 3.1 GENOTYPE CALL BASED ####
##################################################
mkdir -p $pop_dir/ibd/genotypeCalls

vcf_in=$PWORK/$set_id/vcf/$set_id.populations.snps.08filt.vcf.gz # Filtered input VCF without outgroups and undesired individuals created in vcf_filtering_sub.sh
sample_file=$pop_dir/ibd/genotypeCalls/$set_id.samples.pops.txt # List with indviduals in first column and associated populations in second column, without row or column names
out=$pop_dir/ibd/genotypeCalls/$set_id

## Submit script to infer pairwise F_ST between populations
sbatch --output=$pop_dir/logFiles/$set_id.hierfstat.oe $scripts_dir/hierfstat.sh $scripts_dir $vcf_in $sample_file $out

## Conduct Mantel tests and plot IBD
geo_dist=$pop_dir/ibd/genotypeCalls/geo_dist.txt # Distance matrix with mean geographic distances between population pairs as estimated with Geographic Distance Matrix Generator v1.2.3 (https://biodiversityinformatics.amnh.org/open_source/gdmg/), with row and column names
gen_dist=$pop_dir/ibd/genotypeCalls/gen_dist.txt # Distance matrix with pairwise F_ST between populations as estimated with hierfstat, with row and column names
sbatch --output=$pop_dir/logFiles/$set_id.genotypeCalls.ibd.oe $scripts_dir/ibd.sh $geo_dist $gen_dist $pop_dir/ibd/genotypeCalls/$set_id

##################################################
#### 3.2 GENOTYPE LIKELIHOOD BASED ####
##################################################
mkdir -p $pop_dir/ibd/genotypeLikelihoods

saf_dir=$PWORK/$set_id/angsd # Directory with site allele frequency likelihoods for each population inferred in genotype_likelihoods_sub.sh
pops="ambavala ambatovy ambohitantely anjanaharibe anjiahely ankafobe marojejy tsinjoarivo"
nt=24

## Get pairwise population comparisons
awk -f $scripts_dir/combinations.awk <<< $pops > $pop_dir/ibd/genotypeLikelihoods/pop_combinations.txt

## Initialize summary file
echo "pair unweighted weighted" > $pop_dir/ibd/genotypeLikelihoods/$set_id.fst_sumstats.txt

## Get pairwise F_ST between populations
while read combination
do
	first=$(awk '{print $1}' <<< $combination)
	second=$(awk '{print $2}' <<< $combination)
	
	echo -e "#### Processing combination $first $second ...\n"
	# Estimate joint minor allele frequency spectrum
	out_file=$pop_dir/ibd/genotypeLikelihoods/$set_id.$first.$second.ml
	sbatch --job-name=fst --output=$pop_dir/logFiles/$set_id.realsfs.2d.$first.$second.oe $scripts_dir/realsfs.sh $nt "$saf_dir/$first.saf.idx $saf_dir/$second.saf.idx" $out_file
	
	# Estimate F_ST values
	out=$pop_dir/ibd/genotypeLikelihoods/$set_id.$first.$second
	sbatch --wait --job-name=fst --dependency=singleton --output=$pop_dir/logFiles/$set_id.realsfs_fst.$first.$second.oe $scripts_dir/realsfs_fst.sh $nt "$saf_dir/$first.saf.idx $saf_dir/$second.saf.idx" $out
	
	# Write to summary file 
	echo ${first}_$second $(cat $out.fst) >> $pop_dir/ibd/genotypeLikelihoods/$set_id.fst_sumstats.txt
done < $pop_dir/ibd/genotypeLikelihoods/pop_combinations.txt

## Conduct Mantel tests and plot IBD
geo_dist=$pop_dir/ibd/genotypeLikelihoods/geo_dist.txt # Distance matrix with mean geographic distances between population pairs as estimated with Geographic Distance Matrix Generator v1.2.3 (https://biodiversityinformatics.amnh.org/open_source/gdmg/), with row and column names
gen_dist=$pop_dir/ibd/genotypeLikelihoods/gen_dist.txt # Distance matrix with weighted pairwise F_ST between populations as estimated with realSFS, with row and column names
sbatch --output=$pop_dir/logFiles/$set_id.genotypeLikelihoods.ibd.oe $scripts_dir/ibd.sh $geo_dist $gen_dist $pop_dir/ibd/genotypeLikelihoods/$set_id


#################################################################
#### AMOVA ####
#################################################################
mkdir -p $pop_dir/amova
vcf_in=$PWORK/$set_id/vcf/$set_id.populations.snps.08filt.vcf.gz # Filtered input VCF without outgroups and undesired individuals created in vcf_filtering_sub.sh
sample_file=$pop_dir/amova/$set_id.samples.pops.txt # List with indviduals in first column, populations in second column and major clusters (north/south/central) in third column, without row or column names
OUT=$pop_dir/amova/$set_id

## Submit script for AMOVA
sbatch --output=$pop_dir/logFiles/$set_id.amova.oe $scripts_dir/amova.sh $scripts_dir $vcf_in $sample_file $out
