#! /bin/bash
#$ -l h_rt=10:00:00,mfree=2G
#$ -o  /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/masks
#$ -e  /net/harris/vol1/home/beichman/DNAShape/reports.nobackup/masks
#$ -m bea
#$ -M annabel.beichman@gmail.com
####### need to figure out: splitting into intervals/chromosomes. a new script? fin whale will be most annoying. bear too. 
# figure out next. 
########## script to generate bed masks for ref genomes in consistent way #########
set -euxo pipefail
# super useful: https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425
# e cancels whole script if any error occurs
# u exits if there's an undefined variable (saves you from passing empty variables) ** super useful ** 
# o: useful! 
module load modules modules-init modules-gs # initialize modules 
module load bedtools/2.29.2 
module load bedops/2.4.35


label=$1 # species label
faiFile=$2 # full path to species fai file 
gff_or_gtf=$3 # full path to species annotation file (gff or gtf is fine); must be gzipped
repeatMaskerBed=$4 # full path to species rep masker file *in bed format *
trfBed=$5 # full path to species trf file *in bed format*  # note for vaquita that trf and rep mask are alrready combined so need a slightly different script.
cpgIslandsBed=$6 # path to cpg islands bed file 
# just once for humans: combine .fai files:
#> allChrs.fai
#for i in {1..22}
#do
#cat homo_sapiens_ancestor_${i}.fa.fai >> allChrs.fai
#done

outdir=/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/$label
mkdir -p $outdir

echo "files used:" > $outdir/readme
echo "fai: " $faiFile >> $outdir/readme
echo "gff/gtf: " $gff_or_gtf >> $outdir/readme
echo "rep masker: " $repeatMaskerBed >> $outdir/readme
echo "trf: " $trfBed >> $outdir/readme
echo "cpgislands:" $cpgIslandsBed >> $outdir/readme

########## make negative bed mask files ####

###### mask exons +-1kb from GTF files ########

# make genome file for bedtools
echo "making chr file for bedtools"
chrLenFile=${faiFile%.fai}.CHRLENGTHSFORBEDTOOLS # make this in line below 
awk 'BEGIN{OFS="\t"} {print $1,$2}' $faiFile > $chrLenFile

# so this will skip "#" lines, then for exons only extract the chrom, start-1 and stop positions. it will then add +-10kb to either end of each exon
# choosing exons because include UTRs and some gffs don't have 'gene' annotations (eg mouse doesn't) so masking anything with exon+-10kb.
# format:
# skip "#" lines

# $3 is gene, cds or exon # note that exons include UTRs, which I want to also mask

# https://bedtools.readthedocs.io/en/latest/content/tools/slop.html # bedtools slop is great: it will not make things negative and will 'clip' to length of chr
#name output file:
echo "converting gff/gtf file to bed and adding +-10kb buffer to exons"
exonfinal=$outdir/${label}.exonMask.fromGFF_or_GTF.plusminus10kb.0based.sorted.merged.bed
zcat $gff_or_gtf | grep -v "#" | awk 'BEGIN{OFS="\t"} {if($3=="exon") print $1,$4-1,$5}' | bedtools slop -i stdin -g $chrLenFile -b 10000 | sort-bed - | bedtools merge -i stdin > $exonfinal

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in gff conversion"
	exit 1
else
	echo "finished"
fi

awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'  $exonfinal > ${exonfinal%.bed}.TOTALAMOUNTOFSEQUENCE.txt
# some exons are repeated/overlapping, but with sort/merge that doesn't matter

######### repeat masker: need to make sure are sorted and bed formatted (not all are )#######
echo "sorting and merging rep mask bed"
repmaskfinal=$outdir/${label}.repeatMasker.0based.sorted.merged.bed # name outfile
awk 'BEGIN {OFS="\t"} {print $1,$2,$3}' $repeatMaskerBed | sort-bed - | bedtools merge -i stdin > $repmaskfinal
# pulling 3 columns because there seem to be some formatting issues
#sort-bed $repeatMaskerBed | bedtools merge -i stdin > $repmaskfinal

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in repmask conversion"
	exit 1
else
	echo "finished"
fi

awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'  $repmaskfinal > ${repmaskfinal%.bed}.TOTALAMOUNTOFSEQUENCE.txt


############## trf output ################
# already in bed format ; merge and sort it over 
echo "sorting and merging trf bed"
trffinal=$outdir/${label}.trf.0based.sorted.merged.bed
bedtools sort -i $trfBed | bedtools merge -i stdin > $trffinal
# using bedtools sort here beause trf is gzipped
exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in trf conversion"
	exit 1
else
	echo "finished"
fi

awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'  $trffinal > ${trffinal%.bed}.TOTALAMOUNTOFSEQUENCE.txt


# output want to have a bed file for the 3 masks: coding sequence +-10kb, rep masker, and trf . need to use these for targets and for 
# what about 'callability' ? mapping depth? 

########## CpG Islands ############
echo "CpG Islands"
cpgIslandsFinal=$outdir/${label}.cpgIslands.0based.sorted.merged.bed
bedtools sort -i $cpgIslandsBed | bedtools merge -i stdin > $cpgIslandsFinal

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in cpg islands conversion"
	exit 1
else
	echo "finished"
fi

awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'  $cpgIslandsFinal > ${cpgIslandsFinal%.bed}.TOTALAMOUNTOFSEQUENCE.txt


########### combine into one giant negative mask #############
echo "combining exons repmasker and trf and cpgislands"
bedops --merge $exonfinal $repmaskfinal $trffinal $cpgIslandsFinal > $outdir/${label}.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed
exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in merging0"
	exit 1
else
	echo "finished"
fi

awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $outdir/${label}.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed > $outdir/${label}.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.TOTALAMOUNTOFSEQUENCE.txt


####### leave out cpg islands ############
echo "combining exons rm and trf"
bedops --merge $exonfinal $repmaskfinal $trffinal > $outdir/${label}.exon10kb.repmask.trf.NEGATIVEMASK.merged.bed
exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in merging1"
	exit 1
else
	echo "finished"
fi

awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $outdir/${label}.exon10kb.repmask.trf.NEGATIVEMASK.merged.bed > $outdir/${label}.exon10kb.repmask.trf.NEGATIVEMASK.merged.TOTALAMOUNTOFSEQUENCE.txt

####### combine just the repeats (to match JAR's RM+trf) ###########
echo "combining repeats only (rm and trf)"
bedops --merge $repmaskfinal $trffinal > $outdir/${label}.repeatsOnly.repmask.trf.NEGATIVEMASK.merged.bed

exitVal=$?
if [ ${exitVal} -ne 0 ]; then
	echo "error in merging2"
	exit 1
else
	echo "finished"
fi

awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $outdir/${label}.repeatsOnly.repmask.trf.NEGATIVEMASK.merged.bed > $outdir/${label}.repeatsOnly.repmask.trf.NEGATIVEMASK.merged.TOTALAMOUNTOFSEQUENCE.txt