#! /bin/bash
#$ -l h_rt=100:00:00,mfree=6G
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N bears_mutypervariants
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/unified_mutyper
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/unified_mutyper
#$ -t 1-31

# what's different about bears: 31 intervals (intervals not chromosomes)
######################## bears ################

interval=${SGE_TASK_ID} # skipping xy un m for now; to label chr or interval

module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.12 # for filtering out fixed sites (higher version of bcftools)

module load bedtools/2.29.2

set -euo pipefail

kmersize=7 # just want the basic categories of A--> G etc.
label=${kmersize}'mer'

species=bears
reference=brown_bear
vcfdir=/net/harris/vol1/home/beichman/bears/variant_calling/mapped_to_${reference}/vcfs/vcf_20200916_${reference}/interval_${interval}/SNPsOnly/phased
# this vcf file has had sites that couldn't be polarized removed. 
vcffilename=PHASED.mergedSamples.POLARIZED.SNPs.filt_variants_4b.LowCovIndsREMOVED.nomissfiltapplied.PASS.WARN.Only.CustomFilt.HardFilt.AnnotVar.TrimAlt.ref_${reference}.int_${interval}.withINFO.vcf.gz


NEGATIVEMASK=/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/brown_bear_GCF_003584765.1/perInterval/interval${interval}.brown_bear_GCF_003584765.1.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed
# note this is a *NEGATIVE* mask (regions you DONT want )

ancestralFastafilename=/net/harris/vol1/home/beichman/bears/analyses/ancestralReferenceFastas/$reference/interval_${interval}/MODIFIED.ANCESTRAL.${reference}.interval_${interval}.fasta


wd=/net/harris/vol1/home/beichman/DNAShape/analyses/mutyper/unified_mutyper_results/$label/$species
variantdir=$wd/mutyper_variant_files/
spectrumdir=$wd/mutyper_spectrum_files/
ksfsdir=$wd/mutyper_ksfs_files/

mkdir -p $wd
mkdir -p $variantdir
mkdir -p $spectrumdir
mkdir -p $ksfsdir
 

mutypervariantsoutfile=$variantdir/${species}.interval${interval}.mutyper.variants.SomeRevComped.noMissingData.NegativeMasked.PASS.${kmersize}mer.vcf.gz # no longer specifying strict; simplifying name 


bcftools view -c 1:minor -T ^$NEGATIVEMASK -m2 -M2 -v snps -f PASS -Ou $vcfdir/$vcffilename | bcftools view -g ^miss -Ou  | mutyper variants --k $kmersize --chrom_pos 0  $ancestralFastafilename - | bcftools convert -Oz -o ${mutypervariantsoutfile}

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper variants"
	exit 1
else
	echo "finished"
fi


