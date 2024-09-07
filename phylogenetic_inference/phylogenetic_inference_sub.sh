################################################################################
#### PHYLOGENETIC INFERENCE ####
################################################################################

scripts_dir=/home/nibtve93/scripts/phylogeneticInference

set_id=lehilahytsara
phyl_dir=$PWORK/$set_id/phylogeneticInference
vcf_in=$PWORK/$set_id/vcf/populations.snps.07filt.vcf 

mkdir -p $phyl_dir/logFiles

#################################################################
#### 1 MAXIMUM LIKELIHOOD INFERENCE WITH ASCERTAINMENT BIAS CORRECTION ####
#################################################################
mkdir -p $phyl_dir/ml

## Convert VCF file to PHYLIP format and calculate basic alignment statistics
format=phylip # Output format of alignment (phylip or nexus)
sbatch --job-name=ml_inference --output=$phyl_dir/logFiles/vcf_convert.$format.$set_id.oe $scripts_dir/vcf_convert.sh $scripts_dir $vcf_in $phyl_dir/ml $format

## Run phylogenetic inference with ascertainment bias correction in RAxML-NG
nt=80
bootstrap=100 # Number of bootstrap replicates
outgroup="Mmur_RMR44,Mmur_RMR45,Mmur_RMR49" # Outgroup individuals in alignment file
sbatch --job-name=ml_inference --dependency=singleton --output=$phyl_dir/logFiles/raxml_asc.$set_id.oe $scripts_dir/raxml_asc.sh $nt $bootstrap $outgroup $phyl_dir/ml/populations.snps.07filt.noinv.phy $phyl_dir/ml/populations.snps.07filt.noinv.tre

#################################################################
#### 2 QUARTET-BASED INFERENCE FOR INDIVIDUAL AND POPULATION ASSIGNMENT ####
#################################################################
mkdir -p $phyl_dir/quartet

## Thin VCF file by 10,000 bp intervals to ensure independence of SNPs
vcftools --vcf $vcf_in --thin 10000 --recode --recode-INFO-all --stdout > $(dirname $vcf_in)/$(basename $vcf_in .vcf).thin10k.vcf

## Convert VCF file to NEXUS format and calculate basic alignment statistics
format=nexus # Output format of alignment (phylip or nexus)
sbatch --wait --output=$phyl_dir/logFiles/vcf_convert.$format.$set_id.oe $scripts_dir/vcf_convert.sh $scripts_dir $vcf_in $phyl_dir/quartet $format

## Create taxon partitions block files
# For population assignment, the file has to be created manually
# For individual assignment:
echo "BEGIN SETS;" > $phyl_dir/quartet/$set_id.taxPartitions.individual.nex
echo -e "\t TAXPARTITION SPECIES =" >> echo "BEGIN SETS;" > $phyl_dir/quartet/$set_id.taxPartitions.individual.nex
for i in $(awk '{ print $1 }' $phyl_dir/ml/populations.snps.07filt.noinv.phy | tail -n+2)
do
	echo -e "\t\t${i}:${i}," >> echo "BEGIN SETS;" > $phyl_dir/quartet/$set_id.taxPartitions.individual.nex
done
sed -i '$ s/,$//' echo "BEGIN SETS;" > $phyl_dir/quartet/$set_id.taxPartitions.individual.nex # Remove last comma
echo -e "\t\t;" >> echo "BEGIN SETS;" > $phyl_dir/quartet/$set_id.taxPartitions.individual.nex
echo "END;" >> echo "BEGIN SETS;" > $phyl_dir/quartet/$set_id.taxPartitions.individual.nex
echo "" >> echo "BEGIN SETS;" > $phyl_dir/quartet/$set_id.taxPartitions.individual.nex

## Create PAUP block files
nt=80
seed=$RANDOM
# For population assignment:
echo "BEGIN PAUP;" > $phyl_dir/quartet/$set_id.paup.population.nex
echo -e "\toutgroup SPECIES.Mmurinus;" >> $phyl_dir/quartet/$set_id.paup.population.nex
echo -e "\tset root=outgroup outroot=monophyl;" >> $phyl_dir/quartet/$set_id.paup.population.nex
echo -e "\tsvdq nthreads=$nt evalQuartets=all taxpartition=SPECIES bootstrap=standard seed=${SEED};" >> $phyl_dir/quartet/$set_id.paup.population.nex
echo -e "\tsavetrees format=Newick file=$phyl_dir/quartet/$set_id.populations.snps.07filt.population.svdq.tre savebootp=nodelabels;" >> $phyl_dir/quartet/$set_id.paup.population.nex
echo -e "\tquit;" >> $phyl_dir/quartet/$set_id.paup.population.nex
echo "END;" >> $phyl_dir/quartet/$set_id.paup.population.nex
# For individual assignment
echo "BEGIN PAUP;" > $phyl_dir/quartet/$set_id.paup.individual.nex
echo -e "\toutgroup SPECIES.Mmur_RMR44 SPECIES.Mmur_RMR45 SPECIES.Mmur_RMR49;" >> $phyl_dir/quartet/$set_id.paup.individual.nex
echo -e "\tset root=outgroup outroot=monophyl;" >> $phyl_dir/quartet/$set_id.paup.individual.nex
echo -e "\tsvdq nthreads=$nt evalQuartets=all taxpartition=SPECIES bootstrap=standard seed=${SEED};" >> $phyl_dir/quartet/$set_id.paup.individual.nex
echo -e "\tsavetrees format=Newick file=$phyl_dir/quartet/$set_id.populations.snps.07filt.individual.svdq.tre savebootp=nodelabels;" >> $phyl_dir/quartet/$set_id.paup.individual.nex
echo -e "\tquit;" >> $phyl_dir/quartet/$set_id.paup.individual.nex
echo "END;" >> $phyl_dir/quartet/$set_id.paup.individual.nex

## Concatenate files and submit SVDquartets job
for i in population individual
do
	cat $phyl_dir/quartet/populations.snps.07filt.thin10k.noinv.nex $phyl_dir/quartet/$set_id.taxPartitions.$i.nex $phyl_dir/quartet/$set_id.paup.$i.nex > $phyl_dir/quartet/$set_id.paup.$i.concat.nex
	sbatch --job-name=quartet_inference --output=$phyl_dir/logFiles/svdq.$set_id.oe $scripts_dir/svdq.sh $phyl_dir/quartet/$set_id.paup.$i.concat.nex $phyl_dir/quartet/$set_id.paup.$i.concat.nex.log
done

#################################################################
#### 3 APPROXIMATELY UNBIASED (AU) TEST ####
#################################################################
mkdir -p $phyl_dir/au_test

## Infer phylogenies under specified constraints
nt=80
bootstrap=0 # Number of bootstrap replicates
outgroup="Mmur_RMR44,Mmur_RMR45,Mmur_RMR49" # Outgroup individuals in alignment file
constraints="ambatovy.anjz.monophyl ambatovy.monophyl east.west north.south snapp snapp.polytomy svdq svdq.polytomy" # String with prefixes of files containing constrained trees for AU test (format $PREFIX.constraint.nwk)
for i in $constraints
do
	constraint_file=$i.constraint.nwk
	command="-g $constraint_file"
	sbatch --dependency=afterok:$(squeue --noheader --format %i --name ml_inference):$(squeue --noheader --format %i --name quartet_inference) --output=$phyl_dir/logFiles/raxml_asc.$set_id.constraint.$i.oe $scripts_dir/raxml_asc.sh \
		$nt $bootstrap $outgroup $phyl_dir/ml/populations.snps.07filt.noinv.phy $phyl_dir/au_test/populations.snps.07filt.noinv.$i.tre "$command"
done

## Concatenate tree constrained and unconstrained tree files
cat $phyl_dir/au_test/populations.snps.07filt.noinv.*.tree.*support $phyl_dir/ml/populations.snps.07filt.noinv.tre.*support $phyl_dir/quartet/$set_id.populations.snps.07filt.individual.svdq.tre > $phyl_dir/au_test/au_test.trees

## Perform AU test
nt=40
sbatch --output=$phyl_dir/logFiles/au_test.$set_id.oe $scripts_dir/au_test.sh $nt $phyl_dir/ml/populations.snps.07filt.noinv.phy $phyl_dir/au_test/au_test.trees $phyl_dir/au_test/$set_id.au

