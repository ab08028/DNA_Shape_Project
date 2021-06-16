#! /bin/bash
#$ -l h_rt=20:00:00,h_data=50G
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N sandboxMouseModel
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/modeling
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/modeling


######## trying to run modeling script in hoffman #######


module load modules modules-init modules-gs # initialize modules 
module load pcre2/10.35 R/4.0.4

# need to set up R environment: 
R 
install.packages("tidymodels","tidyverse","workflows","tune","vip","ranger","ggplot2","ggrepel","ggbeeswarm","reshape2","devtools",verbose=T)
# need a data dir with shapes
# need a data dir with species data 
# these are just inside the R script
#spectrumdir=/net/harris/vol1/home/beichman/DNAShape/spectrumDataForModeling/mouse
#shapedir=/net/harris/vol1/home/beichman/DNAShape/shapeDataForModeling

gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/DNA_Shape_Project/analyses/modeling
scriptdir=$gitdir/analyses/modeling
script=20210616.TryRunningOnHoffman.R

Rscript $scriptdir/$script