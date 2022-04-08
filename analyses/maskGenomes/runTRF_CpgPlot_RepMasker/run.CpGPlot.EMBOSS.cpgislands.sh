#! /bin/bash
#$ -l h_rt=50:00:00,mfree=6G
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/repeats
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/repeats
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N emboss_cpgplot_brownbear
############################ brown bear #########################
###### detect cpg islands #######
set -euo pipefail

module load modules modules-init modules-gs # initialize modules 

module load EMBOSS/6.6.0
#http://emboss.open-bio.org/rel/rel6/apps/cpgplot.html
module load bedtools/2.29.2 

# fin whale : already has them masked from vcf. could rerun it on them anway; or could get from Meixi
# vaquita: JAR provided file
# bear : needs to be run
# dog/mouse/human can I find from UCSC? 
# so maybe just run on fin whale and bear 

# to download on UCSC: table browser --> choose version of genome --> all tracks --> CpG islands --> bed format --> name file/gzip --> download and xfer to sage
# doing for: mouse, minke, human, dog
# running for: bear
# have from JAR: vaquita
# NOTE FOR MINKE WHALE: cpg islands were wrong chrom names from ucsc, but meixi converted them so use her version of cpg islands (updated excel file)

# so just need to get cpg islands for the brown bear:
# going to run with default params: seems similar to what ucsc does
# http://emboss.sourceforge.net/apps/cvs/emboss/apps/cpgplot.html
wd=/net/harris/vol1/home/beichman/reference_genomes/brown_bear
genome=$wd/brown_bear.fasta
cpgplot $genome -outfeat $wd/CpGIslands.fromcpgplot.emboss.output.gff -auto -graph none
# auto disables prompts
# dfaults windows of 100, minoe 0.6, min len of island is 200
exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in cpgplot"
	exit 1
else
	echo "finished"
fi

# convert from gff3 to bed file (skip header) and sort/merge:
grep -v "#" $wd/CpGIslands.fromcpgplot.emboss.output.gff |  awk 'BEGIN {OFS="\t"} {print $1,$4-1,$5}' | bedtools sort -i stdin | bedtools merge -i stdin > $wd/CpGIslands.fromcpgplot.emboss.output.bed