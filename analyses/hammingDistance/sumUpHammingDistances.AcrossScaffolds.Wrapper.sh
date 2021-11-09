#! /bin/bash
#$ -l h_rt=10:00:00,h_data=5G
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/hamming
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/hamming
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N sumUpHamming

####### sum up human ksfses ########
module load modules modules-init modules-gs # initialize modules 
module load pcre2/10.35 hdf5/1.10.1 R/4.0.4
module load gcc/10.2.0 # necessary for reshape2

script=sumUpHammingDistances.AcrossScaffolds.R
scriptdir=/Users/annabelbeichman/Documents/UW/DNAShapeProject/scripts/DNA_Shape_Project/analyses/hammingDistance

Rscript $scriptdir/$script

