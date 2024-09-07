################################################################################
#### DEMOGRAPHIC MODELLING ####
################################################################################
scripts_dir=/home/nibtve93/scripts/demographicModeling

set_id=lehilahytsara
demo_dir=$PWORK/$set_id/demographicModeling
saf_dir=$PWORK/$set_id/angsd # Directory with site allele frequency likelihoods for each population inferred in genotype_likelihoods_sub.sh

mkdir -p $demo_dir/logFiles

#################################################################
#### 1 MODEL DEMOGRAPHIC HISTORY ####
#################################################################
pops="riamalandy ambavala ambatovy ambohitantely anjanaharibe anjiahely ankafobe marojejy tsinjoarivo metapopulation"
nt=24

## Get pairwise population comparisons
awk -f $scripts_dir/combinations.awk <<< $pops > $demo_dir/pop_combinations.txt

## Estimate joint minor allele frequency (MAF) spectra for selected population comparisons
while read combination
do
	first=$(awk '{print $1}' <<< $combination)
	second=$(awk '{print $2}' <<< $combination)
	
	echo -e "#### Processing combination $first $second ...\n"
	
	# Across all sites for demographic modeling
	out_file=$demo_dir/$set_id.$first.$second.ml
	sbatch --output=$demo_dir/logFiles/$set_id.realsfs.2d.$first.$second.oe $scripts_dir/realsfs.sh $nt "$saf_dir/$first.saf.idx $saf_dir/$second.saf.idx" $out_file
	# Over blocks of 10000 bp (with bootstrapping) to generate confidence intervals in demographic modelling
	options="-nSites 10000 -bootstrap 100"
	out_file=$demo_dir/$set_id.$first.$second.block.ml
	sbatch --output=$demo_dir/logFiles/$set_id.realsfs.2d.block.$first.$second.oe $scripts_dir/realsfs.sh $nt "$saf_dir/$FIRST.saf.idx $saf_dir/$SECOND.saf.idx" $out_file "$options"
done < $demo_dir/pop_combinations.txt

## See Dryad digital repository (https://doi.org/10.5061/dryad.dncjsxkvz) for scripts to assemble block bootstraps and run fastsimcoal v2.6.0.3 (http://cmpg.unibe.ch/software/fastsimcoal27/).


#################################################################
#### 2 CHARACTERIZE POPULATION SIZE CHANGES THROUGH TIME ####
#################################################################
## Software:
stairwayplot=/home/nibtve93/software/stairway_plot_v2.1.1/stairway_plot_es # (v2.0; https://github.com/xiaoming-liu/stairway-plot-v2)

pops="ambavala ambatovy ankafobe metapopulation"
nt=24
mu="1.2e-8" # Mutation rate for study system
gen_time="2.5" # Generation time for study system
seed=$RANDOM

for i in $pops
do
	echo -e "#### Processing population $i ...\n"
	
	## Estimate minor allele frequency (MAF) spectra for selected population
	out_file=$demo_dir/$i.sfs
	sbatch --wait --output=$demo_dir/logFiles/realsfs.1d.$i.oe $scripts_dir/realsfs.sh $nt $saf_dir/$i.saf.idx $out_file
	
	## Print blueprint file for Stairway Plot
	in_file=$PWORK/$set_id/angsd/angsd_$i.txt # List of samples in the population
	no_seq=$(( $(cat $IN_FILE | wc -l) * 2 )) # Multiplied by 2 because of diploid genome
	saf_out=$PWORK/$set_id/angsd/logFiles/angsd.saf.$i.oe
	
	rm $demo_dir/populationSize/$i.blueprint; touch $demo_dir/populationSize/$i.blueprint
	echo -e "popid: $i" >> $demo_dir/populationSize/$i.blueprint
	echo -e "nseq: $no_seq" >> $demo_dir/populationSize/$i.blueprint
	echo -e "L: $(cat $saf_out | grep nSites | sed -n 1p | awk '{ print $3}')" >> $demo_dir/populationSize/$i.blueprint
	echo -e "whether_folded: true" >> $demo_dir/populationSize/$i.blueprint
	echo -e "SFS: $(cut -d' ' -f2-$(( 1 + $no_seq/2)) $out_file)" >> $demo_dir/populationSize/$i.blueprint
	echo -e "smallest_size_of_SFS_bin_used_for_estimation: 1" >> $demo_dir/populationSize/$i.blueprint
	echo -e "largest_size_of_SFS_bin_used_for_estimation: $(cat $IN_FILE | wc -l)" >> $demo_dir/populationSize/$i.blueprint
	echo -e "pct_training: 0.67" >> $demo_dir/populationSize/$i.blueprint
	echo -e "nrand: $(( ($no_seq-2)/4 )) $(( ($no_seq-2)/2 )) $(( ($no_seq-2)*3/4 )) $(( $no_seq-2 ))" >> $demo_dir/populationSize/$i.blueprint
	echo -e "project_dir: $demo_dir/populationSize" >> $demo_dir/populationSize/$i.blueprint
	echo -e "stairway_plot_dir: $stairwayplot" >> $demo_dir/populationSize/$i.blueprint
	echo -e "ninput: 200" >> $demo_dir/populationSize/$i.blueprint
	echo -e "random_seed: $seed" >> $demo_dir/populationSize/$i.blueprint
	echo -e "mu: $mu" >> $demo_dir/populationSize/$i.blueprint
	echo -e "year_per_generation: $gen_time" >> $demo_dir/populationSize/$i.blueprint
	echo -e "plot_title: $i" >> $demo_dir/populationSize/$i.blueprint
	echo -e "xrange: 0,0" >> $demo_dir/populationSize/$i.blueprint
	echo -e "yrange: 0,0" >> $demo_dir/populationSize/$i.blueprint
	echo -e "xspacing: 2" >> $demo_dir/populationSize/$i.blueprint
	echo -e "yspacing: 2" >> $demo_dir/populationSize/$i.blueprint
	echo -e "fontsize: 12" >> $demo_dir/populationSize/$i.blueprint
	
	## Submit script to run Stairway Plot
	sbatch --output=$demo_dir/logFiles/stairwayplot.$i.oe $scripts_dir/stairwayplot.sh $demo_dir/populationSize/$i.blueprint
done

## Plot Stairway Plot output
Rscript $scripts_dir/plot_stairwayplot.R $demo_dir/populationSize/$(echo $pops | awk '{print $1}').final.summary $demo_dir/populationSize/$(echo $pops | awk '{print $2}').final.summary $demo_dir/populationSize/$(echo $pops | awk '{print $3}').final.summary $demo_dir/populationSize/$(echo $pops | awk '{print $4}').final.summary

