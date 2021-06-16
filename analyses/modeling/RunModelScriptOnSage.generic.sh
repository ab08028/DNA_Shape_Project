#! /bin/bash
#$ -l h_rt=20:00:00,h_data=10G
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N sandboxMouseModel
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/modeling
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/modeling
#$ -pe serial 10

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
script=20210616.TryRunningOnSage.RF.multispecies.R

Rscript $scriptdir/$script