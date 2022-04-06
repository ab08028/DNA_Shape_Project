#! /bin/bash
#$ -l h_rt=50:00:00,mfree=6G
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/mutyper
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/mutyper
#$ -m bea
#$ -M annabel.beichman@gmail.com

######## generic script #########
# need either/or for a few things
# removing individuals due to quality or other issues (gorilla, chimp, vaquita, fin whale)
# using strict or not (humans only)
# whether to filter on pass or not (fin whales dont just get pass sites)

# for apes, need to go to the interval in question using -R because vcf files aren't per chromosome. 

set -eou pipefail

configfile=$1
interval_or_chr_or_all=$2
prepend0=$3
todaysdate=`date +%Y%m%d` # for the log file 


echo " these are the config settings "
cat $configfile
echo "end of config settings"

########## get interval (or chr) number ############
if [ $interval_or_chr_or_all = "interval" ]
then
	# prepend 01 02 if needed 
	if [ $prepend0 = "TRUE" ]
	then
		interval=`printf %02d ${SGE_TASK_ID}`
	elif [ $prepend0 = "FALSE" ]
	then
		interval=${SGE_TASK_ID}
	else
		echo "invalid prepend0 option"
		exit 1
	fi
	# set interval label: 
	intervalLabel=interval_${interval}
elif [ $interval_or_chr_or_all = "chr" ]
then
	interval=${SGE_TASK_ID}
	intervalLabel=chr${interval} # label for output files 
elif [ $interval_or_chr_or_all= "allautos" ]
then
	interval='' # no intervals 
	intervalLabel="allautos"
fi

############### NOW SOURCE CONFIG FILE ONCE INTERVALS HAVE BEEN SET ###############
source $configfile

echo " these are the config settings "
cat $configfile
echo "end of config settings"

ls $vcfdir/$vcffile
################# restricting to only PASS sites ##################
if [ $passOption = "TRUE" ]
then
	echo "select only PASS sites"
	pass_snippet='-f PASS' 
elif [ $passOption = "FALSE" ]
then
	echo "not restricting to PASS sites"
	pass_snippet=''
else
	echo "NOT A VALID passOption option"
	exit 1
fi

#################### using --strict ######################
if [ $strictOption = "TRUE" ]
then
	echo "using strict for $species - are you sure you want to do this? should be humans only"
	strict_snippet='--strict'
elif [ $strictOption = "FALSE" ]
then
	echo "not using strict"
	strict_snippet=''
else
	echo "NOT A VALID strictOption option"
	exit 1
fi

############ if you need to subset the input vcf by chromosome before processing (apes need this) ########
if [ $vcfNeedsToBeSubsetByChr = "TRUE" ]
then
	echo "need to subset the vcf by chromosome/interval first (needed for apes)"
	subset_vcf_snippet='-R $intervalLabel'
elif [ $vcfNeedsToBeSubsetByChr = "FALSE" ]
then
	echo "don't need to subset vcf by chr/interval"
	subset_vcf_snippet=''
else
	echo "not a valid vcfNeedsToBeSubsetByChr option"
fi

##################### removing individuals ###################
# this checks if variable is empty: 
if [ -e $individualsToExclude ]
then
	echo "no individuals to remove"
	rm_inds_snippet=''
else
	rm_inds_snippet='-s $individualsToExclude'
fi


############ build mutyper variants code #########

# code snippets: 
# MUST use double quotes! 
initialize_subsetifneeded_snippet="bcftools view $subset_vcf_snippet $vcfdir/$vcffilename -Ou" # initialize with this. if no subsetting needed it'll just start reading the vcf. saves issues with needing to supply vcf at diff parts of pipe
filter_snippet="bcftools view $rm_inds_snippet -T ^$NEGATIVEMASK -m2 -M2 -v snps $pass_snippet -Ou" # if passOption=False then $pass_snippet will be ''; if no inds to remove it will be blank
no_fixed_sites_snippet="bcftools view -c 1:minor -Ou" # this will removed 0/0 sites but keep in fixed 1/1 sites
missing_data_snippet="bcftools view -g ^miss -Ou" # removes missing data 
mutyper_variants_snippet="mutyper variants --k $kmersize --chrom_pos 0 $strict_snippet $ancestralFastafilename -  | bcftools convert -Oz -o ${mutypervariantsoutputname}" # if strictoption=FALSE then it will be '' and not used

# could make and run a piece of code instead? with a header and -t and everything? maybe that's the way. and then submit it?
# but still would need something for fin whale intervals that are 01 02 03 : think about on wednesday! 
# assemble code:
 # need to do this part at some point: 
wd=/net/harris/vol1/home/beichman/DNAShape/analyses/mutyper/unified_mutyper_results/$label/$species
variantdir=$wd/mutyper_variant_files/
spectrumdir=$wd/mutyper_spectrum_files/
ksfsdir=$wd/mutyper_ksfs_files/

#mkdir -p $wd
#mkdir -p $variantdir
#mkdir -p $spectrumdir
#mkdir -p $ksfsdir

######### need to subset vcf prior to processing? ############
# requires double quotes (single don't work)
lineOfCode="${initialize_subsetifneeded_snippet} | ${filter_snippet} | ${no_fixed_sites_snippet} | ${missing_data_snippet} | ${mutyper_variants_snippet}" 
echo $lineOfCode > ${species}.${todaysdate}.log
# always start with bcftools view then add in other snippets ; if you don't need to pre-subset the vcf by chr then $subset_vcf_snippet will be empty and it'll just read the vcf

# run the code: (not yet )
#### $lineOfCode
####### need to deal with intervals/chromosomes --- especially for vaquita. wrapper script? how to pull -t from config file?  