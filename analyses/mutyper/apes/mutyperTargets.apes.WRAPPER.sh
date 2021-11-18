# wrap over species (chromosomes are sge task id )

gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/DNA_Shape_Project
scriptdir=$gitdir/analyses/mutyper/apes/
script=$scriptdir/mutyperTargets.apes.sh


speciesList='Gorilla Pan_paniscus Pan_troglodytes Pongo_abelii Pongo_pygmaeus'
for species in $speciesList
do

qsub $script $species

done
