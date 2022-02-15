# get ape species lists for hamming

vcfdir=/net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs

sampleListDir=/net/harris/vol1/home/beichman/apes/sampleLists
> $sampleListDir/apes.txt
for spp in Gorilla Pan_paniscus Pan_troglodytes Pongo_abelii Pongo_pygmaeus
do
zcat ${spp}.vcf.gz | grep -v "##" | grep -m1 "#" | sed 's/#//g' | sed 's/\t/\n/g' | grep $spp | awk -v spp=$spp 'BEGIN {OFS="\t"}; {print spp,$1}' >> $sampleListDir/apes.txt
done