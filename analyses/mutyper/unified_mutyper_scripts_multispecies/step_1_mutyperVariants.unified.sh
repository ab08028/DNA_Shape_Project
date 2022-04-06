#! /bin/bash
#$ -l h_rt=50:00:00,mfree=6G
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/mutyper
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/mutyper
#$ -m bea
#$ -M annabel.beichman@gmail.com


module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.12 # for filtering out fixed sites (higher version of bcftools)
module load bedtools/2.29.2

######## generic script #########
# need either/or for a few things
# removing individuals due to quality or other issues (gorilla, chimp, vaquita, fin whale)
# using strict or not (humans only)
# whether to filter on pass or not (fin whales dont just get pass sites)

# for apes, need to go to the interval in question using -R because vcf files aren't per chromosome. 

set -exou pipefail

configfile=$1

todaysdate=`date +%Y%m%d` # for the log file 



############### NOW SOURCE CONFIG FILE ONCE INTERVALS HAVE BEEN SET ###############
source $configfile

###### set up outdir/outfiles #######
wd=/net/harris/vol1/home/beichman/DNAShape/analyses/mutyper/unified_mutyper_results/$label/$species
variantdir=$wd/mutyper_variant_files/
spectrumdir=$wd/mutyper_spectrum_files/
ksfsdir=$wd/mutyper_ksfs_files/

mkdir -p $wd
mkdir -p $variantdir
mkdir -p $spectrumdir
mkdir -p $ksfsdir

#make a copy of config file as it was used:
cp $configfile $wd/COPYOFCONGIGFILEUSEDON.${todaysdate}.txt


######## set up log file #########
log=$wd/${species}.${todaysdate}.mutyper_variants.log
> $log

echo "vcffile: $vcfdir/$vcffilename" >> $log
echo "ancestral fasta: $ancestralFastafilename" >> $log
echo "negative mask: $NEGATIVEMASK" >> $log

######### set up output name ###########
mutypervariantsoutputname=${species}.int_or_chr_${interval}.mutyper.variants.SomeRevComped.SeeLogForFilters.${kmersize}mer.vcf.gz

################# restricting to only PASS sites ##################
if [ $passOption = "TRUE" ]
then
	echo "select only PASS sites" >> $log
	pass_snippet='-f PASS' 
elif [ $passOption = "FALSE" ]
then
	echo "not restricting to PASS sites" >> $log
	pass_snippet=''
else
	echo "NOT A VALID passOption option"
	exit 1
fi

#################### using --strict ######################
if [ $strictOption = "TRUE" ]
then
	echo "using --strict for $species - are you sure you want to do this? should be humans only" >> $log
	strict_snippet='--strict'
elif [ $strictOption = "FALSE" ]
then
	echo "not using --strict" >> $log
	strict_snippet=''
else
	echo "NOT A VALID strictOption option"
	exit 1
fi

############ if you need to subset the input vcf by chromosome before processing (apes need this) ########
if [ $vcfNeedsToBeSubsetByChr = "TRUE" ]
then
	echo "need to subset the vcf by chromosome/interval first (needed for apes)" >> $log
	subset_vcf_snippet='-R $intervalLabel'
elif [ $vcfNeedsToBeSubsetByChr = "FALSE" ]
then
	echo "don't need to subset vcf by chr/interval" >> $log
	subset_vcf_snippet=''
else
	echo "not a valid vcfNeedsToBeSubsetByChr option"
fi

##################### removing individuals ###################
# this checks if variable is empty: 
if [ -e $individualsToExclude ]
then
	echo "no individuals to remove" >> $log
	rm_inds_snippet=''
else
	echo "removing the following individuals : $individualsToExclude " >> $log
	rm_inds_snippet='-s $individualsToExclude'
fi


############ build mutyper variants code #########

# code snippets: 
# MUST use double quotes! 
initialize_subsetifneeded_snippet="bcftools view $subset_vcf_snippet $vcfdir/$vcffilename -Ou" # initialize with this. if no subsetting needed it'll just start reading the vcf. saves issues with needing to supply vcf at diff parts of pipe
filter_snippet="bcftools view $rm_inds_snippet -T ^$NEGATIVEMASK -m2 -M2 -v snps $pass_snippet -Ou" # if passOption=False then $pass_snippet will be ''; if no inds to remove it will be blank
no_fixed_sites_snippet="bcftools view -c 1:minor -Ou" # this will removed 0/0 sites but keep in fixed 1/1 sites
missing_data_snippet="bcftools view -g ^miss -Ou" # removes missing data 
mutyper_variants_snippet="mutyper variants --k $kmersize --chrom_pos 0 $strict_snippet $ancestralFastafilename -  | bcftools convert -Oz -o $variantdir/${mutypervariantsoutputname}" # if strictoption=FALSE then it will be '' and not used

######### need to subset vcf prior to processing? ############
# requires double quotes (single don't work)
lineOfCode="${initialize_subsetifneeded_snippet} | ${filter_snippet} | ${no_fixed_sites_snippet} | ${missing_data_snippet} | ${mutyper_variants_snippet}" 
echo "FINAL CODE LINE: \n\n"
echo $lineOfCode >> $log

# make a function and run it:
myMutyperVariantsFunc() {

	"$lineOfCode" # needs double quotes around whole thing otherwise get '|' around pipes.

}
# always start with bcftools view then add in other snippets ; if you don't need to pre-subset the vcf by chr then $subset_vcf_snippet will be empty and it'll just read the vcf

# run the code: (not yet )
myMutyperVariantsFunc
####### need to deal with intervals/chromosomes --- especially for vaquita. wrapper script? how to pull -t from config file?  