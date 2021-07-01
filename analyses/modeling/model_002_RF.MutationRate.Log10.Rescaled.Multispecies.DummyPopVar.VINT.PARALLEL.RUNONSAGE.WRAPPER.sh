#! /bin/bash
#$ -l h_rt=40:00:00,mfree=2G
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/modeling/vint
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/modeling/vint
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N model_002_vint_parallel
#$ -pe serial 5
#$ -t 1-10

#### NEED TO CHANGE THIS BASED ON HOW MANY TOP PARAMETERS YOU'RE USING
# THIS ISN'T WELL AUTOMATED --- JUST GET IT WORKING WELL ENOUGH FOR NOW.
# NEED TO MANUALLY CALCULATE 5 CHOOSE 2 = 10 TO GET THE SIZE OF THE ARRAY
######## trying to run modeling script in hoffman #######

###### SAVE NEW WRAPPER SCRIPT EACH TIME ######
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

# but date this running script each time

todaysdate=`date +%Y%m%d`
##### CHANGE THESE APPROPRIATELY TO BE ABOUT SCRIPT YOU"RE RUNNING #########
#outcomeLabel="FracSegSites.Log10"
#modelLabel="RF"
#otherInfo="Multispecies.DummyPopVar"

#description=${modelLabel}.${outcomeLabel}.${otherInfo}



description=model_002_RF.MutationRate.Log10.Rescaled.Multispecies.DummyPopVar


indir="/net/harris/vol1/home/beichman/DNAShape/analyses/modeling/experiments/20210618"_${description}"/" # must end in "/"
#script=${description}.RUNONSAGE.R # this may not be the best way to do this
#script=RF.FracSegSites.Log10.Multispecies.DummyPopVar.2.VINT.RUNONSAGE.R
script=VINT.RUNONSAGE.parallelPairsOfParams.R

# need to specify top number of parameters (5) to compare interactions of and also the index will be SGE_TASK_ID


topXParams=5 ####### IF YOU CHANGE THIS YOU MUST ALSO CHANGE THE -t setting up above to topXparams choose 2 to get all pairs as an array ######### 
Rscript $scriptdir/$script $indir $topXParams $SGE_TASK_ID

# would like output files and the script to be copied to the same outdir as a record. 

