################################################################################
#### VCF FILTERING ####
################################################################################
## Pipeline adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# bgzip needs to be included in $PATH (v1.11; http://www.htslib.org/)

scripts_dir=/home/nibtve93/scripts/vcfFiltering

set_id=lehilahytsara
vcf_raw=$PWORK/$set_id/stacks/populations/populations.snps.vcf
vcf_dir=$PWORK/$set_id/vcf
bam_dir=$PWORK/bamFiles/$set_id
reference=$PWORK/references/mmur3/GCF_000165445.2_Mmur_3.0_genomic_reduced.fna # Reference genome in fasta format, not containing unlocalized chromosomal scaffolds between chromosomal ones

mkdir -p $vcf_dir/logFiles

## Softlink raw VCF to VCF directory
ln -s $vcf_raw $vcf_dir

#################################################################
#### 1 MAIN FILTERING PIPELINE ####
#################################################################
## Filter for minimum depth
min_dp=5
mean_dp=5
sbatch --job-name=vcf_filter_pip --output=$vcf_dir/logFiles/01_filter_min-dp.$set_id.oe $scripts_dir/01_filter_min-dp.sh $vcf_raw $min_dp $mean_dp $vcf_dir/$set_id.populations.snps.01filt.vcf

## Apply three rounds of filtering for missing data across individuals and genotypes
maxmiss_geno1=0.5 # One minus maxmimum missingness across genotypes (filtering round 1), i.e., maximum missingness = 1-$maxmiss_geno1
maxmiss_geno2=0.6 # One minus maxmimum missingness across genotypes (filtering round 2), i.e., maximum missingness = 1-$maxmiss_geno2
maxmiss_geno3=0.7 # One minus maxmimum missingness across genotypes (filtering round 3), i.e., maximum missingness = 1-$maxmiss_geno3
filter_inds=TRUE # Boolean specifying whether to filter for missingness across individuals
maxmiss_ind1=0.875 # Maxmimum missingness across individuals (filtering round 1)
maxmiss_ind2=0.7 # Maxmimum missingness across individuals (filtering round 2)
maxmiss_ind3=0.5 # Maxmimum missingness across individuals (filtering round 3)
sbatch --job-name=vcf_filter_pip --dependency=singleton --output=$vcf_dir/logFiles/02_filter_missing-1.$set_id.oe $scripts_dir/02_filter_missing-1.sh \
	$vcf_dir/$set_id.populations.snps.01filt.vcf $vcf_dir/$set_id.populations.snps.02filt.vcf $maxmiss_geno1 $maxmiss_geno2 $maxmiss_geno3 $filter_inds $maxmiss_ind1 $maxmiss_ind2 $maxmiss_ind3

## Annotate with INFO fields FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance
suffix=auto
sbatch --job-name=vcf_filter_pip --dependency=singleton --output=$vcf_dir/logFiles/03_annot_gatk.$set_id.oe $scripts_dir/03_annot_gatk.sh \
	$vcf_dir/$set_id.populations.snps.02filt.vcf $vcf_dir/$set_id.populations.snps.03filt.vcf $bam_dir $reference $suffix

## Filter for INFO fields FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance
sbatch --job-name=vcf_filter_pip --dependency=singleton --output=$vcf_dir/logFiles/04_filter_gatk.$set_id.oe $scripts_dir/04_filter_gatk.sh \
	$vcf_dir/$set_id.populations.snps.03filt.vcf $vcf_dir/$set_id.populations.snps.04filt-soft.vcf $vcf_dir/$set_id.populations.snps.04filt-hard.vcf $reference

## Filter for maximum depth as (mean depth + 2 * standard deviation) / number of individuals
sbatch --job-name=vcf_filter_pip --dependency=singleton --output=$vcf_dir/logFiles/05_filter_max-dp.$set_id.oe $scripts_dir/05_filter_max-dp.sh \
	$vcf_dir/$set_id.populations.snps.04filt-hard.vcf $vcf_dir/$set_id.populations.snps.05filt.vcf 

## Apply final round of filtering for missing data across individuals and genotypes
maxmiss_geno=0.9 # One minus maxmimum missingness across genotypes, i.e., maximum missingness = 1-$maxmiss_geno
filter_inds=TRUE # Boolean specifying whether to filter for missingness across individuals
maxmiss_ind=0.5 # Maxmimum missingness across individuals
sbatch --job-name=vcf_filter_pip --dependency=singleton --output=$vcf_dir/logFiles/06_filter_missing-2.$set_id.oe $scripts_dir/06_filter_missing-2.sh \
	$vcf_dir/$set_id.populations.snps.05filt.vcf $vcf_dir/$set_id.populations.snps.06filt.vcf $maxmiss_geno $filter_inds $maxmiss_ind

## Apply minor allele count filter
mac=3
sbatch --job-name=vcf_filter_pip --dependency=singleton --output=$vcf_dir/logFiles/07_filter_mac.$set_id.oe $scripts_dir/07_filter_mac.sh \
	$vcf_dir/$set_id.populations.snps.06filt.vcf $vcf_dir/$set_id.populations.snps.07filt.vcf $mac

## Remove outgroup individuals and those not needed for population structure analyses
drop_inds=$vcf_dir/$set_id.samples.drop.txt # List with individuals to remove from VCF (i.e., outgroups and undesired populations), without row or column names
sbatch --job-name=vcf_filter_pip --dependency=singleton --output=$vcf_dir/logFiles/08_drop_inds.$set_id.oe $scripts_dir/08_drop_inds.sh \
	$vcf_dir/$set_id.populations.snps.07filt.vcf $drop_inds $vcf_dir/$set_id.populations.snps.08filt.vcf 
# bgzip final file
bgzip -c $vcf_dir/$set_id.populations.snps.08filt.vcf  > $vcf_dir/$set_id.populations.snps.08filt.vcf.gz
