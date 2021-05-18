module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
module load samtools/1.9 # contains bgzip
module load htslib/1.9 bcftools/1.9 # for filtering out fixed sites


####### ksfs of TTTAAAA>TTTTAAA motif
wd=/net/harris/vol1/home/beichman/vaquita/analyses/mutyper/mutyperResults_20210225_7mer/mutyper_variant_files
# get just the TTTAAAA > TTTTAAA motif for now:
outvcf=$wd/TTTAAAA.TTTTAAA.Variants.vcf
# get header:
zcat $wd/allInds.exceptRelatives.mutyper.variants.mutationTypes.noFixedSites.AncestralDerivedNotRefAlt.SomeRevComped.Strict.7mer.ALLFREQS.vcf.gz | grep "#" > $outvcf
# get sites with motif:
zcat $wd/allInds.exceptRelatives.mutyper.variants.mutationTypes.noFixedSites.AncestralDerivedNotRefAlt.SomeRevComped.Strict.7mer.ALLFREQS.vcf.gz | grep "TTTAAAA>TTTTAAA" >> $outvcf

outdir=/net/harris/vol1/home/beichman/vaquita/analyses/mutyper/mutyperResults_20210225_7mer/mutyper_ksfs_files
ksfsoutfile=$outdir/TTTAAAA.TTTTAAA.ksfs.txt
mutyper ksfs $outvcf > $ksfsoutfile


# compare to ksfs for just A>T and C>G ? 