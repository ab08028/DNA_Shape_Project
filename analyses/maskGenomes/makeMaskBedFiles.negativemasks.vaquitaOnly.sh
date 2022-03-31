#! /bin/bash
#$ -l h_rt=10:00:00,mfree=2G
#$ -o  /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/masks
#$ -e  /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/masks
#$ -m bea
#$ -M annabel.beichman@gmail.com


########## script to generate bed masks for ref genomes in consistent way #########
set -euxo pipefail


module load modules modules-init modules-gs # initialize modules 
module load bedtools/2.29.2 
module load bedops/2.4.35

####### note for vaquita there is no separate trf bed file so gets a special script ########

label=$1 # species label
faiFile=$2 # full path to species fai file 
gff_or_gtf=$3 # full path to species annotation file (gff or gtf is fine); must be gzipped
repeatMaskerPlusTrfBed=$4 # full path to combined species rep masker + trf output *in bed format *


outdir=/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/$label
mkdir -p $outdir


# then a script where you feed in all the gffs and beds and make a set of masks and one combo mask?
# want to keep them separate so you can mask on just one at a time
# don't redo mutyper variants, just redo mutyper ksfs (!) will be much faster
# mask one at a time and then mask all 
# kind of a pain, but will be useful for SI
# also need to do downsampling and pca on 5 inds per species  with individual spectra
# need to figure out shared variation across species 
# one thing at a time

########## make negative bed mask files ####

###### mask exons +-1kb from GTF files ########

# make genome file for bedtools

chrLenFile=${faiFile%.fai}.CHRLENGTHSFORBEDTOOLS # make this in line below 
awk 'BEGIN{OFS="\t"} {print $1,$2}' $faiFile > $chrLenFile

# so this will skip "#" lines, then for exons only extract the chrom, start-1 and stop positions. it will then add +-10kb to either end of each exon
# choosing exons because include UTRs and some gffs don't have 'gene' annotations (eg mouse doesn't) so masking anything with exon+-10kb.
# format:
# skip "#" lines

# $3 is gene, cds or exon # note that exons include UTRs, which I want to also mask

# https://bedtools.readthedocs.io/en/latest/content/tools/slop.html # bedtools slop is great: it will not make things negative and will 'clip' to length of chr
#name output file:
exonfinal=$outdir/${label}.exonMask.fromGFF_or_GTF.plusminus10kb.0based.sorted.merged.bed
zcat $gff_or_gtf | grep -v "#" | awk 'BEGIN{OFS="\t"} {if($3=="exon") print $1,$4-1,$5}' | bedtools slop -i stdin -g $chrLenFile -b 10000 | bed-sort - | bedtools merge -i stdin > $exonfinal

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in gff conversion"
	exit 1
else
	echo "finished"
fi


# some exons are repeated/overlapping, but with sort/merge that doesn't matter

######### repeat masker + trf : need to make sure are sorted and bed formatted (not all are )#######
repeatsfinal=$outdir/${label}.repeatsOnly.repmask.trf.NEGATIVEMASK.merged.bed # name outfile
sort-bed $repeatMaskerPlusTrfBed | bedtools merge -i stdin > $repeatsfinal

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in repmask conversion"
	exit 1
else
	echo "finished"
fi


########### combine into one giant negative mask #############
bedops --merge $exonfinal $repeatsfinal > $outdir/${label}.exon10kb.repmask.trf.NEGATIVEMASK.merged.USETHIS.bed

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in merging"
	exit 1
else
	echo "finished"
fi
