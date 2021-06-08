#! /bin/bash
#$ -l h_rt=40:00:00,h_data=1G
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -t 1-19

######## getting even/odd targets ##########

# chromosome lengths:
wd=/net/harris/vol1/home/beichman/mice/analyses/ancestralReferenceFastas/separate_even_odd_bp_fromRefGenome

# bedtools fmt genome file (chrname\tlength)

lenfile=/net/harris/vol1/home/beichman/reference_genomes/mouse/mm10_aka_mm38/mm10.fa.bedtoolsGenomeFile.ChrLengths.txt # file of chr lengths


#ancestralFastaDir=/net/harris/vol1/home/beichman/mice/analyses/ancestralReferenceFastas/20210225_focalMmd_out1Mmm_out2Ms
# MAKE sure you're doing targets over on ancestral fastas!


# okay want to make a bed file of every other position -- will be big
# but is per chromosome so isn't going to be much bigger than a fasta

chr=${SGE_TASK_ID}
len=`awk -v chr=$chr '{if($1==chr)print $2}' $lenfile`
#cat $lenfile | while read record
#do
#chr=`echo $record | awk '{print $1}'`
#len=`echo $record | awk '{print $2}'`
echo "starting $chr, len $len"
# these positions will be ODD in a 1-based coordinate sstem (chr1 pos 0 in 1-based is pos 1 = odd), even though start coordinate is even 
oddOutfile=$wd/chr${chr}.oddPositionsIfConvertedTo1Based.0basedCoords.AllSites.NoCallabilityFilter.bed
# these positions will be EVEN in a 1-based coordinate sstem (chr1 pos 1 in 1-based is pos 2 = even), even though start coordinate is odd
evenOutfile=$wd/chr${chr}.evenPositionsIfConvertedTo1Based.0basedCoords.AllSites.NoCallabilityFilter.bed 
> $oddOutfile
> $evenOutfile
# subtract one from len for the counter otherwise it adds an extra site at the end of each file beyond the length
for (( COUNTER=0; COUNTER<=(("$len"-1)); COUNTER+=2 )); do
    echo -e "$chr\t$COUNTER\t$((COUNTER+1))" >> $oddOutfile # eventually maybe need some sort of callability mask 
    echo -e "$chr\t$((COUNTER+1))\t$((COUNTER+2))" >> $evenOutfile # eventually maybe need some sort of callability mask 

done



#done

# do I want odd / even to be one based or not? tricky be careful here. 
# be odd in the 1-based land 
# or I guess just be *consistent* and check carefully
#chr1 2 3   # this is 1-based site 3 
#chr1 4 5  #this is 1-based site 5
# chr1 0 1 # this is position 1 in 1-based coords so is odd
# chr1 1 2 # this is position 2 in 1-based coords so is even 

# this makes the bed files 

