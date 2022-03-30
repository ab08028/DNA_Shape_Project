######## Wrapper script for all species #############


# note that vaquita uses a different script #

# format for submitting (except vaquita script submits combo file for repmask+trf)
#label=$1 # species label
#faiFile=$2 # full path to species fai file 
#gtf_or_gff=$3 # full path to species annotation file (gff or gtf is fine); must be gzipped
#repeatMaskerBed=$4 # full path to species rep masker file *in bed format *
#trfBed=$5 # full path to species trf file *in bed format*  # note for vaquita that trf and rep mask are alrready combined so need a slightly different script.

scriptdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/DNA_Shape_Project/analyses/maskGenomes

############ humans ###########
label=humans
faiFile=/net/harris/vol1/home/beichman/reference_genomes/homo_sapiens_ancestor_GRCh38/allChrs.fai
gtf_or_gff=/net/harris/vol1/home/beichman/reference_genomes/human_GRCh38_annotation/GRCh38_latest_genomic.gff.gz
repeatMaskerBed=/net/harris/vol1/home/beichman/reference_genomes/human_GRCh38_annotation/repeatmask.sorted.bed
trfBed=/net/harris/vol1/home/beichman/reference_genomes/human_GRCh38_annotation/hg38.trf.bed.gzÊ

qsub -N $label $scriptdir/makeMaskBedFiles.negativemasks.sh $label $faiFile $gtf_or_gff $repeatMaskerBed $trfBed


######## minke whale ########
label=minke_whale
faiFile=/net/harris/vol1/home/beichman/reference_genomes/minke_whale/minke_whale_genome/GCF_000493695.1_BalAcu1.0/BalAcu1.0.fa.fai
gtf_or_gff=/net/harris/vol1/home/beichman/reference_genomes/minke_whale/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.gff.gz
repeatMaskerBed=/net/harris/vol1/home/beichman/reference_genomes/minke_whale/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_rm.out.bed
trfBed=/net/harris/vol1/home/beichman/reference_genomes/minke_whale/minke_whale_genome/GCF_000493695.1_BalAcu1.0/AnnabelDownloads_notusedbyMeixi/balAcu1.trf.bed.gz

qsub -N $label $scriptdir/makeMaskBedFiles.negativemasks.sh $label $faiFile $gtf_or_gff $repeatMaskerBed $trfBed

###### mouse ###########
label=mouse
faiFile=/net/harris/vol1/home/beichman/reference_genomes/mouse/mm10_aka_mm38/mm10.fa.fai
gtf_or_gff=/net/harris/vol1/home/beichman/reference_genomes/mouse/mm10_aka_mm38/annotation/mm10.ensGene.gtf.gz
repeatMaskerBed=/net/harris/vol1/home/beichman/reference_genomes/mouse/mm10_aka_mm38/REPEATS/mm10.rmsk.bed.gz
trfBed=/net/harris/vol1/home/beichman/reference_genomes/mouse/mm10_aka_mm38/REPEATS/trfMaskChrom/mm10.trf.chr1-19_ONLY.ABcombined.bed

qsub -N $label $scriptdir/makeMaskBedFiles.negativemasks.sh $label $faiFile $gtf_or_gff $repeatMaskerBed $trfBed

########## brown_bear #######
label=brown_bear
faiFile=/net/harris/vol1/home/beichman/reference_genomes/brown_bear/brown_bear.fasta.fai
gtf_or_gff=/net/harris/vol1/home/beichman/reference_genomes/brown_bear/GCF_003584765.1_ASM358476v1_genomic.gff.gz
repeatMaskerBed= # RUNNING ON SAGE --> then need to convert to bed with bedops.
trfBed= # RUNNING ON SAGE >> need to convert to bed?
# STILL RUNNING 
# qsub -N $label $scriptdir/makeMaskBedFiles.negativemasks.sh $label $faiFile $gtf_or_gff $repeatMaskerBed $trfBed


########## dog ##########
label=dog
faiFile=/net/harris/vol1/home/beichman/reference_genomes/canFam3/canFam3.fa.fai
gtf_or_gff=/net/harris/vol1/home/beichman/reference_genomes/canFam3/annotation/canFam3.ensGene.gtf.gz
repeatMaskerBed=/net/harris/vol1/home/beichman/reference_genomes/canFam3/canFam3.fa.RepeatMasker.sorted.bed
trfBed=/net/harris/vol1/home/beichman/reference_genomes/canFam3/canFam3.trf.bed.gz

qsub -N $label $scriptdir/makeMaskBedFiles.negativemasks.sh $label $faiFile $gtf_or_gff $repeatMaskerBed $trfBed

######## vaquita: NOTE HAS A DIFF SCRIPT BECAUSE OF COMBO RM+TRF FILE ############

label=vaquita
faiFile=/net/harris/vol1/home/beichman/reference_genomes/vaquita/mPhoSin1.pri.cur.20190723_rename/mPhoSin1.pri.cur.20190723_rename.fasta.fai
gtf_or_gff=/net/harris/vol1/home/beichman/reference_genomes/vaquita/mPhoSin1.pri.cur.20190723_rename/GCF_008692025.1_mPhoSin1.pri_genomic_rename.gtf.gz
repeatMaskerPlusTrfBed=/net/harris/vol1/home/beichman/reference_genomes/vaquita/mPhoSin1.pri.cur.20190723_rename/mPhoSin1.pri.cur.20190723_rename_repeats_TRF_RM.bed

qsub -N $label $scriptdir/makeMaskBedFiles.negativemasks.vaquitaOnly.sh $label $faiFile $gtf_or_gff $repeatMaskerPlusTrfBed

