#! /bin/bash
#$ -l h_rt=50:00:00,mfree=6G
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/mutyper
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/mutyper
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N lemurs_mutyper
######### script to run mutyper on any vcf/ref pair

######### NOTE: this doesn't yet use any targets or bed files
module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.9 # for filtering out fixed sites

set -o pipefail


todaysdate=`date +%Y%m%d` # date you run mutyper 
#todaysdate=20211101

# had to bgzip this
vcffilename=/net/harris/vol1/home/beichman/lemurs/VCF/L35_d6_wholeGenome.bam.vcf.gz
# need to index
# ONLY ONCE: bcftools index $vcffilename
refgenome=/net/harris/vol1/home/beichman/reference_genomes/gray_mouse_lemur/GCF_000165445.2_Mmur_3.0_genomic.fna # NOTE : NOT YET POLARIZED
species=gray_mouse_lemur
IndsToINCLUDE='GL/CRC-L2,GL/CRC-L6'



# can add as additional params
#sep="\s"
kmersize=7

outdir=/net/harris/vol1/home/beichman/lemurs/analyses/mutyper/mutyperResults_${todaysdate}
mkdir -p $outdir
ksfsdir=$outdir/mutyper_ksfs_files
variantsdir=$outdir/mutyper_variant_files
mkdir -p $variantsdir

mkdir -p $ksfsdir


variantsoutfile=$variantsdir/${species}.mutyper.variants.mutationTypes.noMissingData.noFixedSites.${kmersize}mers.nostrict.unrelatedIndsOnly.NOTPOLARIZED.noXChrom.vcf.gz

# no callability mask.
# no rep mask (yet)

########### mutyper variants : all pops together
# this removes missing data and sites that are fixed: 
# 20211101 modified from will's code 
# adding bed mask 
# need to remove sepcific inds for chimp and gorilla -- how? 


# need to exclude: NC_033692.1 (x chromosome) 
# can only use ^ with -t not with -r in bcftools fyi
bcftools view -t ^'NC_033692.1' $vcffilename -Ou | bcftools view -s $IndsToINCLUDE -Ou | bcftools view -c 1:minor -m2 -M2 -v snps -Ou | bcftools view -g ^miss -Ou | mutyper variants --k $kmersize --sep "\s" $refgenome - | bcftools convert -Oz -o ${variantsoutfile} # from Will's code, different way to output bcftools output

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper variants"
	exit 1
else
	echo "finished"
fi


######## skipping spectrum for now; just need 7mer ksfs #########

############ ksfs #################


ksfsoutfile=$ksfsdir/${species}.mutyper.ksfs.nostrict.unrelatedIndsOnly.NOTPOLARIZED.txt
# restrict to just one population and sites that are variant for that population: 
bcftools view -c 1:minor $variantsoutfile -Oz | mutyper ksfs - > $ksfsoutfile

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper ksfs"
	exit 1
else
	echo "finished ksfs "
fi

