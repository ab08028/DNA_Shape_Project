######## trying to run modeling script in hoffman #######


module load modules modules-init modules-gs # initialize modules 
module load pcre2/10.35 R/4.0.4
# need to load newer gcc for installations! 
module load gcc/10.2.0
# need to set up R environment: 
#R 
#install.packages("tidymodels")
# install.packages("tidyverse")
# install.packages("workflows")
# install.packages("tune")
# install.packages("vip")
# install.packages("ranger")
# install.packages("ggplot2")
# install.packages("ggrepel")
# install.packages("ggbeeswarm")
# install.packages("reshape2")
# install.packages("devtools")

# need a data dir with shapes
# need a data dir with species data 
# these are just inside the R script
#spectrumdir=/net/harris/vol1/home/beichman/DNAShape/spectrumDataForModeling/mouse
#shapedir=/net/harris/vol1/home/beichman/DNAShape/shapeDataForModeling

gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/DNA_Shape_Project/analyses/modeling
scriptdir=$gitdir/
script=20210616.TryRunningOnSage.RF.multispecies.R # don't have dates on these
# but date this running script each time

todaysdate=`date +%Y%m%d`
##### CHANGE THESE APPROPRIATELY TO BE ABOUT SCRIPT YOU"RE RUNNING #########
outcomeLabel="MutationRate"
modelLabel="RF"
otherInfo="Multispecies.DummyPopVar"

description=${modelLabel}.${outcomeLabel}.${otherInfo}

outdir="/net/harris/vol1/home/beichman/DNAShape/analyses/modeling/experiments/".${todaysdate}"_"${description}"/" #  date specific 

mkdir -p $outdir
cp $scriptdir/$script $outdir/$script.COPYRUNON${todaysdate} # copy it to the outdir 
# this Rscript has !#/usr/bin/env Rscript as setting

qsub -V -o $outdir -e $outdir -N $description \
-l h_rt=20:00:00,mfree=20G -pe serial 10 \
$scriptdir/$script $description $outdir

# would like output files and the script to be copied to the same outdir as a record. 

