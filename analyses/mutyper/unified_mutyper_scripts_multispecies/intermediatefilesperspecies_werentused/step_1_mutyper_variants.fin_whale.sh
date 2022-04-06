#! /bin/bash
#$ -l h_rt=100:00:00,mfree=6G
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N fin_whale_mutypervariants
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/unified_mutyper
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/unified_mutyper
#$ -t 1-96

# what's different about fin whales: 96 intervals (intervals not chromosomes); not all sites are "PASS": some are 
# failCpG_rep which were filtered by Meixi/Sergio. I want to be a bit less conservative on filtering those and so am allowing them to pass and am filtering with my own mask files.
# also the intervals are numbered as 01, 02, etc.
# need to remove bad qual individuals !!!
######################## fin whale ################
# note these are 01..96 so need to prepend a 0 if lt 10
interval=`printf %02d ${SGE_TASK_ID}` # skipping xy un m for now; to label chr or interval

module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.12 # for filtering out fixed sites (higher version of bcftools)

module load bedtools/2.29.2

set -euo pipefail

kmersize=7 # just want the basic categories of A--> G etc.
label=${kmersize}'mer'

species=fin_whale

######vcfdir=/net/harris/vol1/home/beichman/bears/variant_calling/mapped_to_${reference}/vcfs/vcf_${genotypeDate}_${reference}/interval_${interval}/SNPsOnly/phased
# this vcf file has had sites that couldn't be polarized removed. 
########vcffilename=PHASED.mergedSamples.POLARIZED.SNPs.filt_variants_4b.LowCovIndsREMOVED.nomissfiltapplied.PASS.WARN.Only.CustomFilt.HardFilt.AnnotVar.TrimAlt.ref_${reference}.int_${interval}.withINFO.vcf.gz

ancestralFastaDir=/net/harris/vol1/home/beichman/fin_whale/analyses/ancestralReferenceFastas/20211020_focalENP_out1MegNov_out2BalMus

# proper polarized:
ancestralFastafilename=$ancestralFastaDir/FINAL.MODIFIED.ANCESTRAL.interval_${interval}_finWhale_Ancestral_extSFS_gte90-lte10Range_NOTEthisisMinkeWhaleRefGenome.butAncCallsAreforfinWhaleHumpbackBlueWhale_notforMinke.fasta

# NOTE: THESE ARE POLARIZED VCFS WITH UNPOLARIZABLE SITES REMOVED, AND FILTERED TO KEEP IN PASS+FAILCPGREP sites (so I can refilter less conservatively for repeat regions)
# these are new as of 20220405
vcfdir=/net/harris/vol1/home/beichman/fin_whale/analyses/keightley_polarization/20211020_focalENP_out1MegNov_out2BalMus/vcfs_with_ancestral_alleles/perInterval

vcffilename=$vcfdir/OnlySitesWithConfidentAncAllelesFromExtSFS.90-10Range.PASSSITESPlusWarnCpGRefSitesONLY.lowqualremoved.JointCalls_f50b4_08_B_VariantFiltration_${interval}.vcf.gz

indsToRemove=BalAcu02,BalMus01,MegNov01,EubGla01,ENPOR12,ENPAK28 

NEGATIVEMASK=/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/minke_whale_GCF_000493695.1_BalAcu1.0/perInterval/interval${interval}.minke_whale_GCF_000493695.1_BalAcu1.0.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed
# note this is a *NEGATIVE* mask (regions you DONT want )


wd=/net/harris/vol1/home/beichman/DNAShape/analyses/mutyper/unified_mutyper_results/$label/$species
variantdir=$wd/mutyper_variant_files/
spectrumdir=$wd/mutyper_spectrum_files/
ksfsdir=$wd/mutyper_ksfs_files/

mkdir -p $wd
mkdir -p $variantdir
mkdir -p $spectrumdir
mkdir -p $ksfsdir


mutypervariantsoutfile=$variantdir/${species}.interval${interval}.mutyper.variants.SomeRevComped.noMissingData.NegativeMasked.PASSplusoldfailCpgRep.rmBadInds.${kmersize}mer.vcf.gz # no longer specifying strict; simplifying name 

# taking out: -f PASS because don't just want to restrict to PASS (fin whale only!!!) . have already taken out things failing other filters, only things left are PASS and failCpG_Rep!
# need to restrict individuals ! 
bcftools view -T ^$NEGATIVEMASK -s ^BalAcu02,BalMus01,MegNov01,EubGla01,ENPOR12,ENPAK28 -m2 -M2 -v snps -Ou $vcfdir/$vcffilename | bcftools view -c 1:minor -Ou | bcftools view -g ^miss -Ou  | mutyper variants --k $kmersize --chrom_pos 0  $ancestralFastafilename - | bcftools convert -Oz -o ${mutypervariantsoutfile}

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper variants"
	exit 1
else
	echo "finished"
fi
