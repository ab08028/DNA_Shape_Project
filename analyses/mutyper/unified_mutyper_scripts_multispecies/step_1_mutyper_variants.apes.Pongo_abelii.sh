#! /bin/bash
#$ -l h_rt=50:00:00,mfree=6G
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N Pongo_abelii_mutypervariants
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/unified_mutyper
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/unified_mutyper
#$ -t 1-22

# Pongo_abelii doesn't need anything special (not removing any inds), except as all apes: need to jump to right chr
######################## Pongo_abelii ################

intervallabel=chr${SGE_TASK_ID} # skipping xy un m for now; to label chr or interval

module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.12 # for filtering out fixed sites (higher version of bcftools)

module load bedtools/2.29.2

set -euo pipefail

kmersize=7 # just want the basic categories of A--> G etc.
label=${kmersize}'mer'

species=Pongo_abelii
speciesLabel=apes_${species}

vcfdir=/net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs

vcffilename=${species}.vcf.gz

# apes have own anc ref genomes made by Michael Goldberg 
ancestralFastafilename=/net/harris/vol1/home/beichman/apes/polarized_ref_fastas/hg18_references_with_gagp_ancestral_alleles_exclude_recurrent/${species}_${intervallabel}.fa


# apes use same negative mask as humans (no longer using ape callability mask)
NEGATIVEMASK=/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/humans_GRCh38/perInterval/${intervallabel}.humans_GRCh38.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed ## ONE NEGATIVE MASK FILE PER INTERVAL
# note this is a *NEGATIVE* mask (regions you DONT want )

########################## run mutyper #############################

wd=/net/harris/vol1/home/beichman/DNAShape/analyses/mutyper/unified_mutyper_results/$label/$species
variantdir=$wd/mutyper_variant_files/
spectrumdir=$wd/mutyper_spectrum_files/
ksfsdir=$wd/mutyper_ksfs_files/

mkdir -p $wd
mkdir -p $variantdir
mkdir -p $spectrumdir
mkdir -p $ksfsdir


mutypervariantsoutfile=$variantdir/${species}.${intervallabel}.mutyper.variants.SomeRevComped.noMissingData.NegativeMasked.PASS.${kmersize}mer.vcf.gz # no longer specifying strict; simplifying name 

# first jump with -R to the correct chromosome since the vcf has all chrs, then apply chr-specific negative mask: 

bcftools view -R $intervalLabel $vcfdir/$vcffilename -Ou | bcftools view -T ^$NEGATIVEMASK -m2 -M2 -v snps -f PASS -Ou | bcftools view -c 1:minor -Ou | bcftools view -g ^miss -Ou  | mutyper variants --k $kmersize --chrom_pos 0 $ancestralFastafilename - | bcftools convert -Oz -o ${mutypervariantsoutfile}
### this used strict for humans (don't use for any others) #### 
exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper variants"
	exit 1
else
	echo "finished"
fi
