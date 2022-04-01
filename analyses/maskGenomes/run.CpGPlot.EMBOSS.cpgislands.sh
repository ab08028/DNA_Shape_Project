#! /bin/bash
#$ -l h_rt=50:00:00,mfree=6G
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/repeats
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/repeats
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N emboss_cpgplot_brownbear
############################ brown bear #########################
###### detect cpg islands #######
module load modules modules-init modules-gs # initialize modules 

module load EMBOSS/6.6.0
#http://emboss.open-bio.org/rel/rel6/apps/cpgplot.html

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
cpgplot $genome -outfeat $wd/cpgplot.emboss.output.cpgislands.gff -auto
# auto disables prompts
# dfaults windows of 100, minoe 0.6, min len of island is 200