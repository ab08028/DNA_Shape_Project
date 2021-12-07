#! /bin/bash
#$ -l h_rt=20:00:00,h_data=8G
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/distance
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/distance
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N mergeApesForHamming

# want to merge ape vcfs
######### NOTE: this doesn't yet use any targets or bed files
module load modules modules-init modules-gs # initialize modules 
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.9 # for filtering out fixed sites

vcfdir=/net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs

outdir=/net/harris/vol1/home/beichman/apes/merged_vcf

# merge them with bcftools
# all lets all things be merged so no dup records are created (hopefully)
# use the bcf of PanTro 

whole_callability_mask=/net/harris/vol1/home/beichman/apes/callability_mask/Intersect_filtered_cov8.bed.gz # same for all species 


#bcftools merge -m all $vcfdir/Gorilla.vcf.gz $vcfdir/Pan_paniscus.vcf.gz $vcfdir/Pan_troglodytes.bcf $vcfdir/Pongo_abelii.vcf.gz $vcfdir/Pongo_pygmaeus.vcf.gz -Ou | bcftools view -c 1:minor -m2 -M2 -v snps -Ou | bcftools convert -Oz -o $outdir/ALLAPES.mergedVCFs.ForHammingDistance.vcf.gz
# then restrict to just snps again 
# needs index and bgzip
#
# note this is going to be a bit different from other species 
# want  to exclude sex chromosomes and only be chrs 1-22
# needs to be indexed
bcftools index $outdir/ALLAPES.mergedVCFs.ForHammingDistance.vcf.gz
bcftools view -R ${whole_callability_mask} $outdir/ALLAPES.mergedVCFs.ForHammingDistance.vcf.gz > $outdir/ALLAPES.mergedVCFs.ForHammingDistance.CALLABILITYMASKED.vcf.gz