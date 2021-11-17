#! /bin/bash
#$ -l h_rt=50:00:00,h_data=1G
#$ -o /net/harris/vol1/home/beichman/vaquita/reports/mutyper/
#$ -e /net/harris/vol1/home/beichman/vaquita/reports/mutyper/
#$ -m bea
#$ -M annabel.beichman@gmail.com

######### script to run mutyper on any vcf/ref pair

######### NOTE: this doesn't yet use any targets or bed files
module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.9 # for filtering out fixed sites

set -o pipefail


todaysdate=`date +%Y%m%d` # date you run mutyper 
#todaysdate=20211101


vcffilename=$1
refgenome=$2
species=$3
label=$4
IndividualsToExclude=$5


# can add as additional params
#sep="\s"
kmersize=7

outdir=/net/harris/vol1/home/beichman/apes/analyses/mutyper/mutyperResults_${todaysdate}
mkdir $outdir
ksfsdir=$outdir/mutyper_ksfs_files
variantsdir=$outdir/mutyper_variant_files

mkdir -p $ksfsdir


variantsoutfile=$variantsdir/${species}.${label}.mutyper.variants.mutationTypes.noMissingData.noFixedSites.${kmersize}mers.nostrict.PASSONLY.BEDMASKED.vcf.gz

# ONLY DO ONCE 
whole_callability_mask=/net/harris/vol1/home/beichman/apes/callability_mask/Intersect_filtered_cov8.bed.gz # same for all species 
# is originally from: /net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs/Callable_regions # is positive mask

# separate mask by chromosomes:
chr_callability_mask=${whole_callability_mask%.bed}.$label.bed
grep -w $label ${whole_callability_mask} | awk 'BEGIN {OFS="\t"}; {print $1,$2,$3}'> ${chr_callability_mask}



########### mutyper variants : all pops together
# this removes missing data and sites that are fixed: 
# 20211101 modified from will's code 
# adding bed mask 
# need to remove sepcific inds for chimp and gorilla -- how? 

if individualsToExclude=="none"
then

bcftools view -c 1:minor -R ${chr_callability_mask} -m2 -M2 -v snps -f PASS -Ou  $vcffilename | bcftools view -g ^miss -Ou |  mutyper variants --k $kmersize $refgenome - | bcftools convert -Oz -o ${variantsoutfile} # from WIll's code, different way to output bcftools output

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper variants"
	exit 1
else
	echo "finished"
fi

elif individualsToExclude!="none"

then
bcftools view -s ^$individualsToExclude -c 1:minor -R ${chr_callability_mask} -m2 -M2 -v snps -f PASS -Ou  $vcffilename | bcftools view -g ^miss -Ou |  mutyper variants --k $kmersize $refgenome - | bcftools convert -Oz -o ${variantsoutfile} # from WIll's code, different way to output bcftools output
exitVal=$?

if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper variants"
	exit 1
else
	echo "finished"
fi


fi



############ ksfs #################


ksfsoutfile=$ksfsoutdir/${species}.${label}.mutyper.ksfs.nostrict.PASSONLY.BEDMASKED.txt
# restrict to just one population and sites that are variant for that population: 
bcftools view -c 1:minor $variantsoutfile | mutyper ksfs - > $ksfsoutfile

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper ksfs"
	exit 1
else
	echo "finished ksfs "
fi

