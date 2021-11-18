#! /bin/bash
#$ -l h_rt=60:00:00,h_data=8G
#$ -o /net/harris/vol1/home/beichman/vaquita/reports/mutyper/
#$ -e /net/harris/vol1/home/beichman/vaquita/reports/mutyper/
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -t 1-22
######## human targets #########
module load modules modules-init modules-gs # initialize modules 
module load python/3.7.7 # need python >3.7
#module load samtools/1.9 # contains bgzip
#module load htslib/1.9 bcftools/1.9 # for filtering out fixed sites

species=$1
# from will
# https://github.com/harrispopgen/mushi-pipelines/blob/main/1KG/1KG.nf

#params.mask = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/StrictMask/20160622.allChr.mask.bed"
# need to sep by chroms
#whole_callability_mask=/net/harris/vol1/home/beichman/apes/callability_mask/Intersect_filtered_cov8.bed.gz # same for all species 
# mask isn't perfect for targets but will do for now.
# maybe should separate mask by chr. or maybe it's okay. 

# separate mask by chromosomes:
# this was done when running mutyper variants so don't redo. 
chr_callability_mask=/net/harris/vol1/home/beichman/apes/callability_mask/Intersect_filtered_cov8.bed.gz.chr${SGE_TASK_ID}.bed
#zcat ${whole_callability_mask} | grep -w $label | awk 'BEGIN {OFS="\t"}; {print $1,$2,$3}'> ${chr_callability_mask}

# using awk to get rid of 4th column that says "strict" -- (will did as well); don't need to change coordinates though since are already bed fmt
# need -w otherwise chr1 and chr 10 chr11 would be grabbed for chr1. 

k=7

ancestralFasta=/net/harris/vol1/home/beichman/apes/polarized_ref_fastas/hg18_references_with_gagp_ancestral_alleles_exclude_recurrent/${species}_chr${SGE_TASK_ID}.fa

outdir=/net/harris/vol1/home/beichman/apes/analyses/mutyper/mutyperResults_2021118/mutyper_targets_files
mkdir -p $outdir
mutyper targets $ancestralFasta --k ${k} --bed ${chr_callability_mask} > $outdir/${species}.mutyper.targets.chr.${SGE_TASK_ID}.NOSTRICT.BedMasked.${k}mers.txt
# okay so here strict is acceptable because IN THE ANCESTRAL FASTA it refers to regions that couldn't be well polarized. only place that strict makes sense.