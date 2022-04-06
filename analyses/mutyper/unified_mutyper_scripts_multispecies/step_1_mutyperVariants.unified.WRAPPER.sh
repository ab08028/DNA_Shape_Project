######### Wrapper script to submit mutyper variants #############
scriptdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/DNA_Shape_Project/analyses/mutyper/unified_mutyper_scripts_multispecies
configdir=$scriptdir/config_files_per_species

script=$scriptdir/step_1_mutyperVariants.unified.sh

speciesList='humans mice bears fin_whale vaquita Gorilla_gorilla Pan_troglodytes Pan_paniscus Pongo_abelii Pongo_pygmaeus'
speciesList='Gorilla_gorilla Pan_troglodytes vaquita'
# for vaquita need to submit with no intervals 

for species in $speciesList
do

configfile=$configdir/config_${species}.sh
source $configfile # load species info 


if [ $interval_or_chr_or_all = "allautos" ]
then
	# if is all autosomes you don't need intervals: 
	qsub -N ${species} $script $configfile

else
	qsub -N ${species} -t 1-${interval_count} $script $configfile
fi

done 

# want to be able to rerun a single species or interval eventually ? 
