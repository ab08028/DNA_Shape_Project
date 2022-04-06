#! /bin/bash
#$ -l h_rt=50:00:00,mfree=6G
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N humans_mutypervariants
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/unified_mutyper
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/unified_mutyper
#$ -t 1-22


######################## humans ################

intervallabel=chr${SGE_TASK_ID} # skipping xy un m for now; to label chr or interval

module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.12 # for filtering out fixed sites (higher version of bcftools)

module load bedtools/2.29.2

set -euo pipefail

kmersize=7 # just want the basic categories of A--> G etc.
label=${kmersize}'mer'

species=humans

vcfdir=/net/harris/vol1/data/30x1KG/
# this vcf file has had sites that couldn't be polarized removed. 
vcffilename=CCDG_13607_B01_GRM_WGS_2019-02-19_${intervallabel}.recalibrated_variants.vcf.gz


ancestralFastafilename=/net/harris/vol1/home/beichman/reference_genomes/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${intervallabel}.fa

 
NEGATIVEMASK=/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/humans_GRCh38/perInterval/${intervallabel}.humans_GRCh38.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed ## ONE NEGATIVE MASSK FILE PER INTERVAL
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



mutypervariantsoutfile=$variantdir/${species}.${intervallabel}.mutyper.variants.SomeRevComped.noMissingData.NegativeMasked.PASS.strictforHumansOnly.${kmersize}mer.vcf.gz # no longer specifying strict; simplifying name 

#### NOTE FOR HUMANS ONLY *** using STRICT because it is meaningful -- don't use for any others 
########### USES STRICT ###################
bcftools view -c 1:minor -T ^$NEGATIVEMASK -m2 -M2 -v snps -f PASS -Ou $vcfdir/$vcffilename | bcftools view -g ^miss -Ou  | mutyper variants --k $kmersize --chrom_pos 0 --strict $ancestralFastafilename - | bcftools convert -Oz -o ${mutypervariantsoutfile}
### this used strict for humans (don't use for any others) #### 
exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper variants"
	exit 1
else
	echo "finished"
fi
