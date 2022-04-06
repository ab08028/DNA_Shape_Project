#! /bin/bash
#$ -l h_rt=50:00:00,mfree=6G
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N vaquita_mutypervariants
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/unified_mutyper
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/unified_mutyper


# what's different for vaquita: need to exclude relatives from a file; only one interval; strict doesn't make a difference because genome isn't soft masked (but don't use it)
######################## vaquita ################

# only one interval intervallabel=chr${SGE_TASK_ID} # skipping xy un m for now; to label chr or interval
intervallabel=allintervals

module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.12 # for filtering out fixed sites (higher version of bcftools)

module load bedtools/2.29.2

set -euo pipefail

kmersize=7 # just want the basic categories of A--> G etc.
label=${kmersize}'mer'

species=vaquita

# need to use polarized vcfs: 
vcfdir=/net/harris/vol1/home/beichman/vaquita/vcfs/SNPsOnly/polarized/20210216
vcffilename=$vcfdir/vaquita_20_simple_PASS_autos_variants_outgroupAlleles_ancAllele.vcf.gz


ancestralFastaDir=/net/harris/vol1/home/beichman/vaquita/analyses/ancestralReferenceFasta/

ancestralFastafilename=$ancestralFastaDir/MODIFIED.ANCESTRAL.vaquitaAncestral.fasta

 
NEGATIVEMASK=/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/vaquita_mPhoSin1/vaquita_mPhoSin1.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed 
# note this is a *NEGATIVE* mask (regions you DONT want )

excludeFile=/net/harris/vol1/home/beichman/scriptsAndGitDirs/VaquitaSpectrum/samplesToExclude.relatives.txt # need to exclude these 5 individuals bc are relatives
# excluding z0001663, z0004380, z0004393, z0004394, z0185383. 


########################## run mutyper #############################

wd=/net/harris/vol1/home/beichman/DNAShape/analyses/mutyper/unified_mutyper_results/$label/$species
variantdir=$wd/mutyper_variant_files/
spectrumdir=$wd/mutyper_spectrum_files/
ksfsdir=$wd/mutyper_ksfs_files/

mkdir -p $wd
mkdir -p $variantdir
mkdir -p $spectrumdir
mkdir -p $ksfsdir



mutypervariantsoutfile=$variantdir/${species}.${intervallabel}.mutyper.variants.SomeRevComped.noMissingData.NegativeMasked.PASS.rmRelatives.${kmersize}mer.vcf.gz # no longer specifying strict; simplifying name 

# restrict to mask first
# note here we are removing relatives !!!! (prior to doing anything else)
bcftools view -T ^$NEGATIVEMASK -S ^$excludeFile -m2 -M2 -v snps -f PASS -Ou $vcfdir/$vcffilename | bcftools view -c 1:minor -Ou | bcftools view -g ^miss -Ou  | mutyper variants --k $kmersize --chrom_pos 0  $ancestralFastafilename - | bcftools convert -Oz -o ${mutypervariantsoutfile}


exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper variants"
	exit 1
else
	echo "finished"
fi


