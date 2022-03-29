#! /bin/bash
#$ -l h_rt=10:00:00,h_data=8G
#$ -o /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/distance
#$ -e /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/distance
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N distCanids
#$ -t 1-38
###### NOTE: for other species, I used mutyper variants files
# but for apes I need to merge their vcfs. am worried about merging post mutyper variants
# because of teh ref allele recoding etc
# so I am merging the pre-mutyper vcfs and bcfs
# but then I am *CALLABILITY MASKING* them in the same way that I did the mutyper variants files
# but note that if you do more masking in future mutyper runs you should probably remask these
# but won't make a huge difference. 
# so the vcf you're using here is merged across all apes and has been callability masked (see merge.ape.vcfs.ForHamming.sh script)
# so can use the targets with the callability mask as approximate denominator 
# experiment with plink
# 
module load modules modules-init modules-gs # initialize modules 
module load  plink/1.90b6.18 # stick to 1.9, plink2 seems v different
# https://www.cog-genomics.org/plink/1.9/distance#read_dists
# this is super fast
# this does allele counts (then maybe could divide myself?)
interval=chr${SGE_TASK_ID}
label=canids
# running it on the same vcf I used for mutyper variants (phased SNPs -- so maybe missing some ; might undercall distance but should be ballpark right)
# maybe I should rerun bears without phasing? tbd 
#vcfdir=/net/harris/vol1/home/beichman/bears/variant_calling/mapped_to_brown_bear/vcfs/vcf_20200916_brown_bear/interval_${interval}/SNPsOnly/phased
# mutyper variants! 
vcfdir=/net/harris/vol1/home/beichman/dogs/analyses/mutyper/mutyperResults_20220103/mutyper_variant_files




### using 7mer mutyper variant files (so will be missing some SNPs that can't be phased/pol (was true when using vcf) or that don't have 7mer context)
# logic is that I want same masking/filtering that is going into the mutation spectrum to be used for the hamming distance. as similar as possible to spectrum calcs
# but that does mean distances are probably slightly undercalled  (though hopefully similar masking is being used for target files which will be denominator)
# so I think this is best -case

# NOTE: unlike other spp this is NOT the vcf I used for mutyper variants
# but I used the callability mask on it the same way so muytper targets is okay to use with it
# this has been callablity masked (not rep masked)
# have to index for excluding regions

#vcf=$vcfdir/ALLAPES.mergedVCFs.ForHammingDistance.CALLABILITYMASKED.vcf.gz
vcf=$vcfdir/${interval}.canid.mutyper.variants.mutationTypes.MissingDataAllowed.noFixedSites.3mers.nostrict.PASSONLY.vcf.gz



outdir=/net/harris/vol1/home/beichman/DNAShape/analyses/hamming_distance/${label}/perInterval
outroot=$outdir/plink.${label}.interval.${interval}.FromMutyperVariantsVCF
mkdir -p $outdir

plink --vcf $vcf --distance square0 flat-missing --const-fid --allow-extra-chr --out $outroot
# ^^ use this! 
# square0 gives you allele counts in a grid with upper triangle 0'd out. ids are in .ids file to combine in R.
# this does the full on hamming distance (but denominator may be too small since this is snps only?)

 #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4342193/
 
 #### note with flat-missing you scale up the distance by 1-missing frequency (so gets a bit bigger)
 # there's another version (default) that scales up by avg contribution of missing to distance --
# i'm going to go with flat missing for now. 

# don't want to use this because vcf doesn't contain monomorphic sites
# plink --vcf SNPs.filt_variants_4b.LowCovIndsREMOVED.nomissfiltapplied.PASS.WARN.Only.CustomFilt.HardFilt.AnnotVar.TrimAlt.ref_brown_bear.int_19.vcf.gz --distance 1-ibs --const-fid --allow-extra-chr
 

 # okay so I think i just want the allele count or 1-ibs (don't want ibs because that is similarity not distance)
 
 # so plink 1-ibs is just the allele count divided by total sites in the vcf (with maybe some missingness correction)
 # so I can just get the allele counts and then divide by total genome size in some fashion 
 # easy peasy! 
 # make sure I understand how the allele count dist is done. I think it's just hamming distance with some missingness correction
 # yes it's basically hamming but then divided by 1-missingness frequency (not totally sure how calc'd) to make distance slightly bigger
 # maybe need to deal with missingness in targets but not yet. for now just use 2*targets as denominator in R
 
