#! /bin/bash
#$ -l h_rt=50:00:00,mfree=6G
#$ -o /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified/reports.nobackup/mutyper_ksfs
#$ -e /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified/reports.nobackup/mutyper_ksfs
#$ -m bea
#$ -M annabel.beichman@gmail.com

############## mutyper ksfs #################
# must be run after variants because you need the masked fasta file (can be run concurrently with targets/spectrum)

module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.12 # for filtering out fixed sites (higher version of bcftools)
module load bedtools/2.29.2

######### set up script ########

set -exou pipefail

configfile=$1

todaysdate=`date +%Y%m%d` # for the log file 



############### NOW SOURCE CONFIG FILE ONCE INTERVALS HAVE BEEN SET ###############
source $configfile
# note that this loads all the variables you need for the species, including kmer size etc ^^^


###### set up outdir/outfiles #######
wd=/net/harris/vol1/home/beichman/allspecies_mutyper_results_unified/$label/$species/mutyper_results_masked_${maskLabel}
variantdir=$wd/mutyper_variant_files
spectrumdir=$wd/mutyper_spectrum_files
ksfsdir=$wd/mutyper_ksfs_files
targetdir=$wd/mutyper_target_files

mkdir -p $wd
mkdir -p $wd/logs
mkdir -p $variantdir
mkdir -p $spectrumdir
mkdir -p $ksfsdir
mkdir -p $targetdir


log=$wd/logs/${species}.${intervalLabel}.${todaysdate}.mutyper_ksfs.log
> $log

################### get mutyper variants file #############
# note this must match step_1 output exactly! trying to think of a way to make that more efficient.
mutypervariantsoutputname=${species}.int_or_chr_${interval}.mutyper.variants.SomeRevComped.SeeLogForFilters.${maskLabel}.${kmersize}mer.vcf.gz


######### iterate through populations if pops are provided #######
# if no pops provided, all inds are part of the population 
# if there aren't any pops defined then run without subsetting inds out: 
echo "pops:" $pops
if [ -z $pops ] # if no pops defined
then
	echo "not splitting into pops" >> $log
	bcftools view -S $popfile $ALLFREQSvariantsOutfile |  bcftools view -c 1:minor  | mutyper ksfs - > $ksfsdir/${species}.int_or_chr_${interval}.mutyper.ksfs.SeeLogForFilters.${maskLabel}.${kmersize}mer.txt
else
	echo "splitting into pops: $pops" >> $log
	for pop in $pops
		do
		echo -e "starting $pop"
		popfile=$poplistdir/${pop}.${popfilesuffix} 

		ksfspopdir=$ksfsdir/$pop
		mkdir -p $ksfspopdir

		######## NOTE: it is *extremely* important for mutyper that the AN and AC fields are correct. BCFTOOLS correctly updates them through filtering (I checked). ####
		# select individuals in file | exclude fixed sites within that population | mutyper variants --> 
		# NOTE this REMOVES in fixed sites per population (sites fixed across all inds have already been removed)

		# per population:
		bcftools view -S $popfile $mutypervariantsoutputname |  bcftools view -c 1:minor  | mutyper ksfs - > $ksfspopdir/${species}.${pop}.int_or_chr_${interval}.mutyper.ksfs.SeeLogForFilters.${maskLabel}.${kmersize}mer.txt
		done

fi


exitVal=$?

if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper variants"
	exit 1
else
	echo "finished script"
fi


######## TO FIGURE OUT: HAVEN'T FIGURED OUT RANDOMIZE ACROSS POPULATIONS WHILE STILL REMOVING FIXED SITES BETWEEN POPULATIONS (NEED TO FIGURE THIS OUT FOR KSFS TOO)
