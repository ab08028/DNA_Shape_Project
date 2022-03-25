#! /bin/bash
#$ -l h_rt=100:00:00,h_data=8G
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

#whole_callability_mask=/net/harris/vol1/home/beichman/apes/callability_mask/Intersect_filtered_cov8.bed.gz # same for all species 

# chr1
# chr10
# chr11
# chr12
# chr13
# chr14
# chr15
# chr16
# chr17
# chr18
# chr19
# chr2
# chr20
# chr21
# chr22
# chr3
# chr4
# chr5
# chr6
# chr7
# chr8
# chr9
# need to remove:
# chrX
# chrY

# restrict this to just 1-22
# only once: zcat $whole_callability_mask | grep -v "chrX" | grep -v "chrY" | bgzip > ${whole_callability_mask%.bed.gz}.chr1-22.bed.gz

whole_callability_mask_1_22=/net/harris/vol1/home/beichman/apes/callability_mask/Intersect_filtered_cov8.chr1-22.bed.gz
#### 20220215: making a change. Since these are SNP sets, they introduce a LOT of missing data that seems to mess up my hamming calculations
# instead want to assume that if a site is seen in Gorilla.vcf but NOT in Pan-pan that doesn't mean that all of Pan_pan inds had missing data there, but rather that they are monomorphic there
# but mono for reference or mono for alternate!?!?! would have to be for reference otherwise would be a fixed site (?)
# so going to use the flag https://samtools.github.io/bcftools/bcftools.html#merge
# --missing-to-ref so that a site that is missing from pan_pan gets set to all 0/0 
# not totally sure if this will work but let's try it

#bcftools merge -m all $vcfdir/Gorilla.vcf.gz $vcfdir/Pan_paniscus.vcf.gz $vcfdir/Pan_troglodytes.bcf $vcfdir/Pongo_abelii.vcf.gz $vcfdir/Pongo_pygmaeus.vcf.gz -Ou | bcftools view -c 1:minor -m2 -M2 -v snps -Ou | bcftools convert -Oz -o $outdir/ALLAPES.mergedVCFs.ForHammingDistance.vcf.gz
#bcftools merge -m all --missing-to-ref $vcfdir/Gorilla.vcf.gz $vcfdir/Pan_paniscus.vcf.gz $vcfdir/Pan_troglodytes.bcf $vcfdir/Pongo_abelii.vcf.gz $vcfdir/Pongo_pygmaeus.vcf.gz -Ou | bcftools view -c 1:minor -m2 -M2 -v snps -Oz -o $outdir/ALLAPES.missingToREF.mergedVCFs.ForHammingDistance.vcf.gz

# then restrict to just snps again 
# needs index and bgzip
#
# note this is going to be a bit different from other species 
# want  to exclude sex chromosomes and only be chrs 1-22
# needs to be indexed

#bcftools index $outdir/ALLAPES.missingToREF.mergedVCFs.ForHammingDistance.vcf.gz
# aha! forgot to bgzip -- whoops! in future will use -Oz 
bcftools view -R ${whole_callability_mask_1_22} $outdir/ALLAPES.missingToREF.mergedVCFs.ForHammingDistance.vcf.gz -Oz -o $outdir/ALLAPES.missingToREF.mergedVCFs.ForHammingDistance.CALLABILITYMASKED.chr1-22only.vcf.gz

# note I then deleted the non callabilit masked vcf to save space 
# still is plenty of missing data 
