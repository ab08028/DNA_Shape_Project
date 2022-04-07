#! /bin/bash
#$ -l h_rt=50:00:00,mfree=6G
#$ -o /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified/reports.nobackup/mutyper_targets
#$ -e /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified/reports.nobackup/mutyper_targets
#$ -m bea
#$ -M annabel.beichman@gmail.com

############## mutyper targets ###########
###### NOTE this step can be run concurrently as mutyper variants 
# get from config file: intervals; ancestral fasta; negative mask
# also need to restrict fasta to just those regions that I'm considering (some sort of positive mask)
# eg bears its only the scaffolds that are in my contigs
# for humans its chrs 1-22 (autosomes)
# for vaquita it's autosomes (1-??)

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


log=$wd/logs/${species}.${intervalLabel}.${todaysdate}.mutyper_targets.log
> $log
echo "ancestral fasta: $ancestralFastafilename" >> $log
echo "" >> $log
echo "negative mask: $NEGATIVEMASK" >> $log
echo "" >> $log
echo "################## COPY OF CONTIG FILE: ##########################" >> $log
echo "" >> $log
cat $configfile >> $log
echo "#########################################################" >> $log

#################### using --strict ######################
# need strict for humans only
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

#################### mask the ancestral fasta using the negative mask file *regions you DONT want * ###########
#### note that these fastas don't contain any non-autosome chroms (made sure for each species. but you should make sure before running)
# note this is a HARD MASK (NNNNN) (soft mask wouldn't work without --strict in mutyper variants).

bedtools maskfasta -fi $ancestralFastafilename -bed $NEGATIVEMASK -fo $targetdir/${ancestralFastafilename%.fasta}.NEGATIVEMASKED.temp.fasta

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in bedtools maskfasta"
	exit 1
else
	echo "finished"
fi



######## mutyper targets: input the masked fasta ###########
targetsoutfilename=${species}.int_or_chr_${interval}.mutyper.targets.SeeLogForFilters.${maskLabel}.${kmersize}mer.txt

mutyper targets ${strict_snippet} --chrom_pos $chrom_pos --k $kmersize $targetdir/${ancestralFastafilename%.fasta}.NEGATIVEMASKED.temp.fasta > $targetdir/$targetsoutfilename

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper targets"
	exit 1
else
	# if it worked then get rid of temp fasta file : 
	rm $targetdir/${ancestralFastafilename%.fasta}.NEGATIVEMASKED.temp.fasta
	echo "finished"
fi


######## note about non-autosome chroms/scaffs not included in genotype calling: made sue they aren't in these anc genomes : #######
# bears: ancestral fasta already in intervals so doesn't contain scaffs <1mb (done)
# fin whale: already in intervals (done)
# human ancestor: already in chromosomes (done)
# apes ancestors: already in chromosomes (done)
# vaquita : aha! this one has non autosomes. have updated anc fasta to just be autosomes 1-21. note for past target calc I was using a positve bed mask so ti was fine. but now need ot deal with it. restrict to 1-21
# mice : already in chromosomes (done)


# note one issue with targets: dont have any sort of callability filter for the all-called sites, so targets is an overestimate. but proportions should be about right.
# could try to improve this at some point. Used to use callability from my vcf files for bear, vaquita; had callability mask for humans/apes. don't have anything for mouse. so not using for any species to try to keep things equal
# more of an issue for mushi when rates matter. less of an issue for spectra.



