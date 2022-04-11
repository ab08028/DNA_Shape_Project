######### Wrapper script to submit mutyper variants #############
# must be run after step_1 because you need the masked ancestral fasta


scriptdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/DNA_Shape_Project/analyses/mutyper/unified_mutyper_scripts_multispecies
configdir=$scriptdir/config_files_per_species

script=$scriptdir/step_3a_mutyper_ksfs.sh 

mkdir -p /net/harris/vol1/home/beichman/allspecies_mutyper_results_unified/reports.nobackup/mutyper_ksfs


speciesList='humans mice bears fin_whale vaquita Gorilla_gorilla Pan_troglodytes Pan_paniscus Pongo_abelii Pongo_pygmaeus'
#speciesList="bears" # vaquita"
for species in $speciesList
do

configfile=$configdir/config_${species}.sh
source $configfile # load species info 

#make error file:


# for vaquita need to submit with no intervals 
if [ $interval_or_chr_or_all = "allautos" ]
then
	# if is all autosomes you don't need intervals: 
	qsub -N ${species}_ksfs $script $configfile

else
	qsub -N ${species}_ksfs -t 1-${interval_count} $script $configfile ### -t 1-${interval_count}
fi

done 
