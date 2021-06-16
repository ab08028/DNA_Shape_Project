#! /bin/bash
#$ -l h_rt=20:00:00,h_data=6G
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N mouseMutyper7merSpectrum
#$ -o /net/harris/vol1/home/beichman/mice/reports/mutyper
#$ -e /net/harris/vol1/home/beichman/mice/reports/mutyper
#$ -t 1-19
########## script to generate test/train dataset based on odd/even bps ##################

###### assumes you've already run mutyper variants ! ############ 
module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.9 # for filtering out fixed sites


set -o pipefail

chromosome=chr${SGE_TASK_ID} # skipping xy un m for now

 
mutyperdir=/net/harris/vol1/home/beichman/mice/analyses/mutyper/mutyperResults_20210317_NOSTRICT_7mer # prexisting path to 7mer mutyper results


variantdir=$mutyperdir/mutyper_variant_files # variant files (assumes you've already run mutyper variants)


poplistdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/BearAdmixtureProject/samples/mouseSamples/populationFilesForMutyper



outdir=$mutyperdir/training_and_test_spectra_perChr
mkdir -p $outdir

vcf=$variantdir/${chromosome}.mutyper.variants.mutationTypes.noFixedSites.AncestralDerivedNotRefAlt.SomeRevComped.NOSTRICT.ALLFREQS.7mer.vcf.gz

# testing even/odd separation method:
# total non-header lines in vcf chr 19 = 2515610

# this will print a line if it's a header line (contains #) *OR* if it's an odd-numbered position (using $2 modulo 2 != 0 )
# odd = train
# even = test
# make sure adds up to total 
# checks: 
# Odd sites:
# bcftools view $vcf | awk '{if(/#/ || $2%2!=0)print}' | grep -v -c "#" # how many total odd sites
# 1258474 total sites -- is 50% 

# Even sites: this will print a line if it's a header line (contains #) *OR* if it's an even-numbered position (using $2 modulo 2 == 0 )
# bcftools view $vcf | awk '{if(/#/ || $2%2==0)print}' | grep -v -c "#" # how many total even sites 
# 1257136 -- 50% but slightly different number which is good 


######### population list:
populations='Mmd Mmc Mmm Ms' # note allNonABC contains EUR, MT and AK and PB contains EGR and WGR
# could also do Mmd_FRA.txt
#Mmd_GER.txt
#(Mmd_HEL.SKIP.txt)
#Mmd_IRA.txt
#Mmm_AFG.txt
#Mmm_CZE.txt
#Mmm_KAZ.txt
# if I want to do sub populations.

# for now keeping singletons and low complexity regions in
for pop in $populations
do
echo -e "starting $pop"
popfile=$poplistdir/${pop}.txt


# don't need to remake variants file for each pop; instead just parse it right into mutyper spectrum. 
######## NOTE: it is *extremely* important for mutyper that the AN and AC fields are correct. BCFTOOLS correctly updates them through filtering (I checked). ####
# select individuals in file | exclude fixed sites within that population | mutyper variants --> 
# NOTE this REMOVES in fixed sites per population (sites fixed across all inds have already been removed)


# selects population, removes fixed sites per population, then separates into even/odd 
# note that for each population the train/test datasets will be at mostly (but not exactly due to removign fixed sites) the same sites
# if mapped to same ref genome
# that seems not ideal but maybe can alternate even/odd when testing across species or something

######### training set (odd bp) ##########
# make SURE it's !=0 here for ODD
bcftools view -S $popfile $vcf |  bcftools view -c 1:minor  | awk '{if(/#/ || $2%2!=0)print}' |  mutyper spectra --population - > $outdir/${chromosome}_${pop}_TRAINING.mutyper.spectra.odd_bpOnly.PERPOPULATION.ALLFREQS.NOSTRICT.txt


########### test set (even bp) ######################
# make SURE it's ==0 here for EVEN

bcftools view -S $popfile $vcf |  bcftools view -c 1:minor  |  awk '{if(/#/ || $2%2==0)print}' |  mutyper spectra --population - > $outdir/${chromosome}_${pop}_TESTING.mutyper.spectra.even_bpOnly.PERPOPULATION.ALLFREQS.NOSTRICT.txt



done

# then in R will combine these across chromosomes