######### Wrapper script to submit mutyper variants #############
# must be run after step_1 because you need the masked ancestral fasta


scriptdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/DNA_Shape_Project/analyses/mutyper/unified_mutyper_scripts_multispecies
configdir=$scriptdir/config_files_per_species

script=$scriptdir/step_4a_mutyper_spectrum.INPROGRESS.sh

# speciesList='humans mice bears fin_whale vaquita Gorilla_gorilla Pan_troglodytes Pan_paniscus Pongo_abelii Pongo_pygmaeus'
speciesList="bears vaquita"
for species in $speciesList
do

configfile=$configdir/config_${species}.sh
source $configfile # load species info 

#make error file:
mkdir -p /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified/reports.nobackup/mutyper_spectrum
# for vaquita need to submit with no intervals 
if [ $interval_or_chr_or_all = "allautos" ]
then
	# if is all autosomes you don't need intervals: 
	qsub -N ${species}_spectrum $script $configfile

else
	qsub -N ${species}_spectrum -t 1-1 $script $configfile ### -t 1-${interval_count}
fi

done 

###### still need to figure out the shared variation issue ######
##### maybe not as much of an issue if I'm only choosing one pop per species for the work? more of an issue for dogs and hamming distance stuff.
