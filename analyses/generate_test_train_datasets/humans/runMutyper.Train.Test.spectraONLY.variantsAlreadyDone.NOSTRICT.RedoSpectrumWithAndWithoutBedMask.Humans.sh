#! /bin/bash
#$ -l h_rt=40:00:00,h_data=1G
#$ -o /net/harris/vol1/home/beichman/vaquita/reports/mutyper/
#$ -e /net/harris/vol1/home/beichman/vaquita/reports/mutyper/
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -t 1-22
######### script to run mutyper on any vcf/ref pair

######### NOTE: this doesn't yet use any targets or bed files
module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.9 # for filtering out fixed sites

set -o pipefail

chromosome=chr${SGE_TASK_ID}

#todaysdate=`date +%Y%m%d` # date you run mutyper 
todaysdate=20210329


sampleList=/net/harris/vol1/home/beichman/humans/sample_information/${pop}.sampleList.txt
# can add as additional params
#sep="\s"
kmersize=7
bedmask=/net/harris/vol1/home/beichman/reference_genomes/human_positive_mask/20160622.allChr.mask.bed # note this mask is based on low coverage 1K genomes
# data -- so may not be as appropriate for high cov data. but hopefully will help us avoid low complexity regions etc
# but may need to do more systematic filtering in future.
# represents callable regions at low cov, so will be callable at high cov, but may be overly stringent compared to other species and make cross species comparisons hard. 
outdir=/net/harris/vol1/home/beichman/humans/analyses/mutyper/mutyperResults_${todaysdate}
variantsdir=$outdir/mutyper_variant_files 
spectrumdir_withMask=$outdir/mutyper_spectrum_files_WITHBEDMASK_PASSONLY
spectrumdir_withoutMask=$outdir/mutyper_spectrum_files_noBEDMASK_PASSONLY

testtraindir_withMask=$outdir/training_and_test_spectra_perChr_WITHBEDMASK_PASSONLY
testtraindir_withoutMask=$outdir/training_and_test_spectra_perChr_noBEDMASK_PASSONLY

mkdir -p $outdir
mkdir -p $variantsdir
mkdir -p $spectrumdir
mkdir -p $testtraindir

# refdir=/net/harris/vol1/home/beichman/reference_genomes/homo_sapiens_ancestor_GRCh38
# reference=homo_sapiens_ancestor_${chr}.fa
# vcf=CCDG_13607_B01_GRM_WGS_2019-02-19_chr${chr}.recalibrated_variants.vcf.gz


# once: get populations
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
#sampleInfo=/net/harris/vol1/home/beichman/humans/sample_information/integrated_call_samples_v3.20130502.ALL.panel
# grep to get rid of header
#superpops=`grep -v sample $sampleInfo | awk '{print $3}' | sort | uniq`
# only need to do once: 
#for superpop in $superpops
#do
# write out sample set:
#awk -v superpop=$superpop '{if($3==superpop)print $1}' $sampleInfo > /net/harris/vol1/home/beichman/humans/sample_information/${superpop}.sampleList.txt

#done
pops='AFR AMR EAS EUR SAS'



########### mutyper variants | mutyper spectrum because don't want to have a 1TB of variants files #######
# gets rid of missing data, and selects only variant sites and then does mutyper variants
# WITHOUT STRICT and without Randomize (do I want population? Yes! )
# NOT using the callability mask for exploration of the particular 7mer
# need to select biallelic snps only

## ACTUALLY: forgot about separating by population
# makes more sense to do variants once and then sep by pop
#### make this more like old script actually! 
# need to restrict to biallelic snps tho : -m2 -M2 -v snps (restricts to biallelic snps not indels)
# indels are present


variantsoutfile=$variantsdir/${chromosome}.mutyper.variants.mutationTypes.noMissingData.noFixedSites.${kmersize}mers.NOSTRICT.vcf.gz


########### mutyper variants : all pops together
# this removes missing data and sites that are fixed: 
#### consider using strict at some point later which in this case means sites that are not well polarized. for now not using. 
# this makes 1TB of files!!! be careful. 
# don't do again -- already did once: zcat $vcffilename | grep -v "\.\/\." | bcftools view -c 1:minor -m2 -M2 -v snps  | mutyper variants --k $kmersize $refgenome - | bgzip -cf > ${variantsoutfile}

#exitVal=$?
#if [ ${exitVal} -ne 0 ]; then
#	echo "error in mutyper variants"
#	exit 1
#else
#	echo "finished"
#fi


###### mutyper spectrum: loop over pops 

for pop in $pops
do
sampleList=/net/harris/vol1/home/beichman/humans/sample_information/${pop}.sampleList.txt
spectrumoutfile_withMask=${spectrumdir_withMask}/${chromosome}_${pop}_mutyper.spectrum.PERPOPULATION.ALLFREQS.NOSTRICT.WITHBEDMASK.PASSONLY.txt
spectrumoutfile_withoutMask=${spectrumdir_withMask}/${chromosome}_${pop}_mutyper.spectrum.PERPOPULATION.ALLFREQS.NOSTRICT.noBEDMASK.PASSONLY.txt

# restrict to just one population and sites that are variant for that population: 
# this is positive mask 
### NEED TO ONLY SELECT PASS SITES!! not lower tranche sites 
bcftools view -S $sampleList -R $bedmask -f PASS $variantsoutfile | bcftools view -c 1:minor | mutyper spectra --population - > ${spectrumoutfile_withMask}

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper spectra"
	exit 1
else
	echo "finished $pop"
fi

bcftools view -S $sampleList -f PASS $variantsoutfile | bcftools view -c 1:minor | mutyper spectra --population - > ${spectrumoutfile_withoutMask}

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper spectra"
	exit 1
else
	echo "finished $pop"
fi


##### make training and test data sets #########

######### training set (odd bp) ##########
# make SURE it's !=0 here for ODD
bcftools view -S $sampleList -R $bedmask -f PASS $variantsoutfile | bcftools view -c 1:minor |  awk '{if(/#/ || $2%2!=0)print}' | mutyper spectra --population - > $testtraindir_withMask/${chromosome}_${pop}_TRAINING.odd_bpOnly.mutyper.spectra.PERPOPULATION.ALLFREQS.NOSTRICT.WITHBEDMASK.PASSONLY.txt

# going to also make without the mask


exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper train spectra"
	exit 1
else
	echo "finished $pop"
fi



bcftools view -S $sampleList -R $bedmask -f PASS $variantsoutfile | bcftools view -c 1:minor |  awk '{if(/#/ || $2%2==0)print}' | mutyper spectra --population - > $testtraindir_withMask/${chromosome}_${pop}_TESTING.even_bpOnly.mutyper.spectra.PERPOPULATION.ALLFREQS.NOSTRICT.WITHBEDMASK.PASSONLY.txt



exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper test spectra"
	exit 1
else
	echo "finished $pop"
fi




##### make training and test data sets WITHOUT mask #########
######### training set (odd bp) ##########
# make SURE it's !=0 here for ODD
bcftools view -S $sampleList -f PASS $variantsoutfile | bcftools view -c 1:minor |  awk '{if(/#/ || $2%2!=0)print}' | mutyper spectra --population - > ${testtraindir_withoutMask}/${chromosome}_${pop}_TRAINING.odd_bpOnly.mutyper.spectra.PERPOPULATION.ALLFREQS.NOSTRICT.noBEDMASK.PASSONLY.txt

# going to also make without the mask


exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper train spectra"
	exit 1
else
	echo "finished $pop"
fi



bcftools view -S $sampleList -f PASS $variantsoutfile | bcftools view -c 1:minor |  awk '{if(/#/ || $2%2==0)print}' | mutyper spectra --population - > ${testtraindir_withoutMask}/${chromosome}_${pop}_TESTING.even_bpOnly.mutyper.spectra.PERPOPULATION.ALLFREQS.NOSTRICT.noBEDMASK.PASSONLY.txt



exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in mutyper test spectra"
	exit 1
else
	echo "finished $pop"
fi





done
