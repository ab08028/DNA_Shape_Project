#! /bin/bash
#$ -l h_rt=50:00:00,mfree=6G
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/repeats
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/repeats
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N TRF_RepMasker



# script to run TRF (tandem repeat finder) https://github.com/Benson-Genomics-Lab/TRF
# and repmasker 
# on genomes that don't yet have a TRF file. (just brown bear for now; polar bear also)
# have trf for: humans (covers apes as well), mice, minke whale, vaquita,
# need it for: brown bear (and polar bear for completeness I suppose)
# lets do it for dog as well 

module load modules modules-init modules-gs # initialize modules 
module load trf/4.09 # load trf (slightly older version)
module load RepeatMasker/4.0.8 

wd=/net/harris/vol1/home/beichman/reference_genomes/brown_bear


genomeFasta=brown_bear.fasta



########### run repeat masker ############
# needs a fair bit of memory (min 10G or gets killed)

# mkdir  -p $wd/RepeatMaskerOutput_parallel
# cd $wd/RepeatMaskerOutput_parallel
# RepeatMasker -species carnivore -pa 10 -dir $wd/RepeatMaskerOutput_parallel $wd/$genomeFasta
# 
# exitVal=$?
# if [ ${exitVal} -ne 0 ]; then
# 	echo "error in repmasker"
# 	exit 1
# else
# 	echo "finished"
# fi


########## run trf ##########
#running trf:
# recommended settings:
##trf yourfile.fa 2 5 7 [keep as is] 80 10 [keep as is] 50 2000 <-- change to 12 to match UCSC
# only thing changing from trf defaults is last parameter: period size (size of repeat itself, not number of copies)
# setting to 12 to match UCSC 
#1. File: The sequence file to be analyzed in FASTA format (see for details). Multiple sequence in the same file are allowed.
#2. Match, Mismatch, and Delta: Weights for match, mismatch and indels. These parameters are for Smith-Waterman style local alignment using wraparound dynamic programming. Lower weights allow alignments with more mismatches and indels. A match weight of 2 has proven effective with mismatch and indel penalties in the range of 3 to 7. Mismatch and indel weights are interpreted as negative numbers. A 3 is more permissive and a 7 less permissive. The recomended values for Match Mismatch and Delta are 2, 7, and 7 respectively.
#3. PM and PI: Probabilistic data is available for PM values of 80 and 75 and PI values of 10 and 20. The best performance can be achieved with values of PM=80 and PI=10. Values of PM=75 and PI=20 give results which are very similar, but often require as much as ten times the processing time when compared with values of PM=80 and PI=10.
#4. Minscore: The alignment of a tandem repeat must meet or exceed this alignment score to be reported. For example, if we set the matching weight to 2 and the minimun score to 50, assuming perfect alignment, we will need to align at least 25 characters to meet the minimum score (for example 5 copies with a period of size 5). [not sure if should keep as is ?] 
#5. Maxperiod: Period size is the program's best guess at the pattern size of the tandem repeat. The program will find all repeats with period size between 1 and 2000, but the output can be limited to a smaller range. [[[ In UCSC this is set to 12 ]]] 

mkdir -p $wd/trf
cd $wd/trf # maybe this will work
trf -h -d 2 5 7 80 10 50 12 $wd/$genomeFasta
# h: suppresss html
# d: produce data file 


exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in trf"
	exit 1
else
	echo "finished"
fi


