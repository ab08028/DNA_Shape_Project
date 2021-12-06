#! /bin/bash
#$ -l h_rt=60:00:00,h_data=8G
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/mutyper
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/mutyper
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N lemurs_targets
######## human targets #########
module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
#module load samtools/1.9 # contains bgzip
#module load htslib/1.9 bcftools/1.9 # for filtering out fixed sites

k=7
species=gray_mouse_lemur

# NOT POLARIZED YET FOR LEMURS:

refgenome=/net/harris/vol1/home/beichman/reference_genomes/gray_mouse_lemur/GCF_000165445.2_Mmur_3.0_genomic.fna.gz # NOTE : NOT YET POLARIZED

outdir=/net/harris/vol1/home/beichman/lemurs/analyses/mutyper/mutyperResults_20211206/mutyper_targets_files
mkdir -p $outdir
mutyper targets $refgenome --k ${k} > $outdir/${species}.mutyper.targets.NOSTRICT.NOTPOLARIZED.NOMASKING.${k}mers.txt
