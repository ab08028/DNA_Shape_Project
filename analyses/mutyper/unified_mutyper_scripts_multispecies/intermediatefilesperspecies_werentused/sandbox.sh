
vcffilename=$1
refgenome=$2
species=$3
label=$4
IndividualsToExclude=$5
sep=$6


kmersize=7

outdir=/net/harris/vol1/home/beichman/apes/analyses/mutyper/mutyperResults_${todaysdate}
mkdir -p $outdir
ksfsdir=$outdir/mutyper_ksfs_files
variantsdir=$outdir/mutyper_variant_files
mkdir -p $variantsdir

mkdir -p $ksfsdir


variantsoutfile=$variantsdir/${species}.${label}.mutyper.variants.${kmersize}mers.noMissingData.noFixedSites.NEGATIVEMASKED.vcf.gz


if individualsToExclude=="none"
then

bcftools view -c 1:minor -^T ${negativeBedMask} -m2 -M2 -v snps -f PASS -Ou $vcffilename | bcftools view -g ^miss -Ou |  mutyper variants --k $kmersize --sep $sep $refgenome - | bcftools convert -Oz -o ${variantsoutfile} # from WIll's code, different way to output bcftools output

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper variants"
	exit 1
else
	echo "finished"
fi

elif individualsToExclude!="none"

then
bcftools view -s ^$individualsToExclude -c 1:minor -^T ${negativeBedMask}  -m2 -M2 -v snps -f PASS -Ou  $vcffilename | bcftools view -g ^miss -Ou |  mutyper variants --k $kmersize --sep $sep $refgenome - | bcftools convert -Oz -o ${variantsoutfile} # from WIll's code, different way to output bcftools output
exitVal=$?

if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper variants"
	exit 1
else
	echo "finished"
fi


fi


if $species=="humans"
then
# need to use strict
echo "using --strict for $label = humans because in the human anc ref sequence, softmasked sites indicate non-polarized. not true for any other species"



fi


if $species=="fin_whale"
then
echo "keeping in both PASS sites and FAIL_CpGRep sites because I want to filter the repeat regions less conservatively than the original vcf. Most should still get masked away"




fi

