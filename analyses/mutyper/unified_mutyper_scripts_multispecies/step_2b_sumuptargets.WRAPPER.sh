#! /bin/bash
#$ -l h_rt=50:00:00,mfree=6G
#$ -o /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified/reports.nobackup/mutyper_spectrum
#$ -e /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified/reports.nobackup/mutyper_spectrum
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N sumup_targets
###################### sum up ksfs ###################
module load modules modules-init modules-gs # initialize modules 
module load pcre2/10.35 hdf5/1.10.1 R/4.0.4
module load gcc/10.2.0 # necessary for reshape2


scriptdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/DNA_Shape_Project/analyses/mutyper/unified_mutyper_scripts_multispecies
script=$scriptdir/step_2b_sumuptargets.R

configdir=$scriptdir/config_files_per_species


#speciesList='mice bears fin_whale vaquita Gorilla_gorilla Pan_troglodytes Pan_paniscus Pongo_abelii Pongo_pygmaeus humans'
speciesList="bears vaquita"
for species in $speciesList
do
echo "starting $species"

configfile=$configdir/config_${species}.sh
source $configfile # load species info l 
# these are coming from config file: $species $interval_count $prepend0; but pops are coming from the indir 

indir=/net/harris/vol1/home/beichman/allspecies_mutyper_results_unified/${label}/${species}/mutyper_results_masked_${maskLabel}/mutyper_target_files/
outdir=/net/harris/vol1/home/beichman/allspecies_mutyper_results_unified/${label}/allspecies_summed_up_over_intervals_forTransfer/${species}/mutyper_results_masked_${maskLabel}/mutyper_target_files
inputfilesuffix=".mutyper.targets.SeeLogForFilters.maskALL.7mer.txt"

mkdir -p $outdir
#echo $indir
#echo $outdir
Rscript $script $indir $outdir $species $interval_count $prepend0 $inputfilesuffix


done