######### Wrapper script to submit mutyper variants #############
scriptdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/DNA_Shape_Project/analyses/mutyper/unified_mutyper_scripts_multispecies
configdir=$scriptdir/config_files_per_species

script=$scriptdir/step_1_mutyperVariants.unified.sh

speciesList='bears,mice,humans,fin_whale,vaquita,Gorilla_gorilla,Pan_troglodytes,Pan_paniscus,Pongo_abelii,Pongo_pygmaeus'

# for vaquita need to submit with no intervals 

for species in $speciesList
do

configfile=$configdir/config_${species}.sh
source $configfile # load species info 

if [ $interval_or_chr_or_all != "allautos" ]
then
	qsub -N ${species} -t 1-1 $script $configfile # -t 1-${interval_count}
elif [ $interval_or_chr_or_all = "allautos"]
then
	# if is all autosomes you don't need intervals: 
	qsub -N ${species} $script $configfile

else
	echo "not a valid interval_or_chr_or_all option"
	exit 1
fi

done 

# want to be able to rerun a single species or interval eventually ? 
