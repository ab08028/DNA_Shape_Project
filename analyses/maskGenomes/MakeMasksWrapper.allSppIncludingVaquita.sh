######## Wrapper script for all species #############
#### this will make/combine negative bed mask files for exons +-10kb from gtf/gff files, rep masker , tandem repeat finder
# these will get sorted and merged into one big bed file that is a *negative mask* (regions you DONT want ) 
# want to mask a) ksfs b) spectra c) targets
# note : no callability masks used (?) have i already used those for mutyper variants (?) I may have. 
# note that vaquita uses a slightly diff script because JAR combined repmask+trf already
# there's an excel file that keeps track of these files

# note that vaquita uses a different script #

# format for submitting (except vaquita script submits combo file for repmask+trf)
#label=$1 # species label
#faiFile=$2 # full path to species fai file 
#gff_or_gtf=$3 # full path to species annotation file (gff or gtf is fine); must be gzipped
#repeatMaskerBed=$4 # full path to species rep masker file *in bed format *
#trfBed=$5 # full path to species trf file *in bed format*  # note for vaquita that trf and rep mask are alrready combined so need a slightly different script.

scriptdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/DNA_Shape_Project/analyses/maskGenomes

############ humans ###########
label=humans_GRCh38
faiFile=/net/harris/vol1/home/beichman/reference_genomes/human_GRCh38_annotation/hg38.fa.gz.fai # note is ancestral fai, which is fine - still is GrCh38
gff_or_gtf=/net/harris/vol1/home/beichman/reference_genomes/human_GRCh38_annotation/hg38.ensGene.gtf.gz # usin ensembl genes
repeatMaskerBed=/net/harris/vol1/home/beichman/reference_genomes/human_GRCh38_annotation/repeatmask.sorted.bed
trfBed=/net/harris/vol1/home/beichman/reference_genomes/human_GRCh38_annotation/hg38.trf.bed.gz
cpgIslands=/net/harris/vol1/home/beichman/reference_genomes/human_GRCh38_annotation/humanGRCh38.CpGIslandsTrack.UCSC.bed.gz
qsub -N $label $scriptdir/makeMaskBedFiles.negativemasks.sh $label $faiFile $gff_or_gtf $repeatMaskerBed $trfBed $CpGIslands


######## minke whale ########
label=minke_whale_GCF_000493695.1_BalAcu1.0
faiFile=/net/harris/vol1/home/beichman/reference_genomes/minke_whale/minke_whale_genome/GCF_000493695.1_BalAcu1.0/BalAcu1.0.fa.fai
gff_or_gtf=/net/harris/vol1/home/beichman/reference_genomes/minke_whale/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.gff.gz
repeatMaskerBed=/net/harris/vol1/home/beichman/reference_genomes/minke_whale/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_rm.out.bed
trfBed=/Users/annabelbeichman/Documents/UW/WhaleProject/fin_whales/fin_whale_scripts/fin_whale_spectrum/KeightleyPolarizationMethod_fin_whale/trf.output.sorted.merged.0based.bed
# note this is version of trf that I ran for minke whale, not the version from ucsc which had weird chr names
cpgIslands=/net/harris/vol1/home/beichman/reference_genomes/minke_whale/minke_whale_genome/GCF_000493695.1_BalAcu1.0/ucsc/CpG_ExtUnmasked.bed
qsub -N $label $scriptdir/makeMaskBedFiles.negativemasks.sh $label $faiFile $gff_or_gtf $repeatMaskerBed $trfBed

###### mouse ###########

# once for mice: need to combine fai across all chrs:
#cd /net/harris/vol1/home/beichman/reference_genomes/mouse/mm10_aka_mm38/REPEATS/trfMaskChrom
#> mm10.trf.chr1-19.ABcombined.bed # combined by me 
#for i in {1..19}
#do
#cat chr${i}.bed >> mm10.trf.chr1-19_ONLY.ABcombined.bed
#done 

label=mouse_mm10
faiFile=/net/harris/vol1/home/beichman/reference_genomes/mouse/mm10_aka_mm38/mm10.fa.ADDEDCHRTONUMBERS.fai
gff_or_gtf=/net/harris/vol1/home/beichman/reference_genomes/mouse/mm10_aka_mm38/annotation/mm10.ensGene.gtf.gz
repeatMaskerBed=/net/harris/vol1/home/beichman/reference_genomes/mouse/mm10_aka_mm38/REPEATS/mm10.rmsk.bed
trfBed=/net/harris/vol1/home/beichman/reference_genomes/mouse/mm10_aka_mm38/REPEATS/trfMaskChrom/mm10.trf.chr1-19_ONLY.ABcombined.bed
cpgIslands=/net/harris/vol1/home/beichman/reference_genomes/mouse/mm10_aka_mm38/CpGIslands_fromUCSC/mm10.CpGIslandsTrack.UCSC.bed.gz

qsub -N $label $scriptdir/makeMaskBedFiles.negativemasks.sh $label $faiFile $gff_or_gtf $repeatMaskerBed $trfBed

########## brown_bear #######
# convert using bedops:
# #gunzip [ BROWN BEAR REP MASK OUTPUT]
# from bedops:
#cd /net/harris/vol1/home/beichman/reference_genomes/brown_bear/RepeatMaskerOutput_parallel
#module load modules modules-init modules-gs # initialize modules 
#module load bedops/2.4.35
#rmsk2bed < brown_bear.fasta.out > brown_bear.fasta.RepeatMasker.output.sorted.bed 
# do I need to convert trf output to bed as well? 

label=brown_bear_GCF_003584765.1
faiFile=/net/harris/vol1/home/beichman/reference_genomes/brown_bear/brown_bear.fasta.fai
gff_or_gtf=/net/harris/vol1/home/beichman/reference_genomes/brown_bear/GCF_003584765.1_ASM358476v1_genomic.gff.gz
repeatMaskerBed=/net/harris/vol1/home/beichman/reference_genomes/brown_bear/RepeatMaskerOutput_parallel/brown_bear.fasta.RepeatMasker.output.sorted.bed
trfBed=/net/harris/vol1/home/beichman/reference_genomes/brown_bear/trf/trf.output.sorted.merged.0based.bed
cpgIslands=/net/harris/vol1/home/beichman/reference_genomes/brown_bear/CpGIslands.fromcpgplot.emboss.output.bed
# generated rep mask, trf and cpg myself for bears
 
qsub -N $label $scriptdir/makeMaskBedFiles.negativemasks.sh $label $faiFile $gff_or_gtf $repeatMaskerBed $trfBed


########## dog ##########


# need to convert dog repmask output with bedops: 
#gunzip canFam3.fa.out.gz
# from bedops:
#rmsk2bed < canFam3.fa.out > canFam3.fa.RepeatMasker.sorted.bed 
# seems like bedtools sort can't handle >3 columns -- just leads to blank! so can use sort-bed in bedops instead -- I think this was actually a variable name issue that has been resolved so either should be ok


label=dog_canFam3
faiFile=/net/harris/vol1/home/beichman/reference_genomes/canFam3/canFam3.fa.fai
gff_or_gtf=/net/harris/vol1/home/beichman/reference_genomes/canFam3/annotation/canFam3.ensGene.gtf.gz
repeatMaskerBed=/net/harris/vol1/home/beichman/reference_genomes/canFam3/canFam3.fa.RepeatMasker.sorted.bed
trfBed=/net/harris/vol1/home/beichman/reference_genomes/canFam3/canFam3.trf.bed.gz
cpgIslands=/net/harris/vol1/home/beichman/reference_genomes/canFam3/canFam3.CpGIslandsTrack.UCSC.bed.gz

qsub -N $label $scriptdir/makeMaskBedFiles.negativemasks.sh $label $faiFile $gff_or_gtf $repeatMaskerBed $trfBed

######## vaquita: NOTE HAS A DIFF SCRIPT BECAUSE OF COMBO RM+TRF FILE ############

label=vaquita_mPhoSin1
faiFile=/net/harris/vol1/home/beichman/reference_genomes/vaquita/mPhoSin1.pri.cur.20190723_rename/mPhoSin1.pri.cur.20190723_rename.fasta.fai
gff_or_gtf=/net/harris/vol1/home/beichman/reference_genomes/vaquita/mPhoSin1.pri.cur.20190723_rename/GCF_008692025.1_mPhoSin1.pri_genomic_rename.gtf.gz
repeatMaskerPlusTrfBed=/net/harris/vol1/home/beichman/reference_genomes/vaquita/mPhoSin1.pri.cur.20190723_rename/mPhoSin1.pri.cur.20190723_rename_repeats_TRF_RM.bed
cpgIslands=/net/harris/vol1/home/beichman/reference_genomes/vaquita/mPhoSin1.pri.cur.20190723_rename/mPhoSin1.pri.cur.20190723_rename_cpgIslands.bed

qsub -N $label $scriptdir/makeMaskBedFiles.negativemasks.vaquitaOnly.sh $label $faiFile $gff_or_gtf $repeatMaskerPlusTrfBed

