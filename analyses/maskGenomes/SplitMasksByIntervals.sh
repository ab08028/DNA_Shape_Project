######### need to split up masks ##############
# note that sum of resulting interval files may not equal original bed files because some chrs are getting excluding (like M and X and Y and unplaced scaffs)
 
suffixForMask=exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed # the mask to use for each species in format SPECIES.$suffix
###### mice:
label=mouse_mm10
wd=/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/$label
mkdir $wd/perInterval
negativemask_allchr=${label}.${suffixForMask}
for i in {1..19} # 19 for mice 
do
awk -v chr=chr${i} 'BEGIN {OFS="\t"} {if($1==chr)print}' $wd/$negativemask_allchr > $wd/perInterval/chr${i}.${negativemask_allchr}


done


###### humans:
label=humans_GRCh38
wd=/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/$label
mkdir $wd/perInterval
negativemask_allchr=${label}.${suffixForMask}
for i in {1..22} # 22 for humans 
do
awk -v chr=chr${i} 'BEGIN {OFS="\t"} {if($1==chr)print}' $wd/$negativemask_allchr > $wd/perInterval/chr${i}.${negativemask_allchr}


done


##### apes mapped to hg18 #########
label=apes_mapped_to_hg18
wd=/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/apes_mapped_to_hg18
mkdir $wd/perInterval
negativemask_allchr=${label}.${suffixForMask}
for i in {1..22} # 22 for humans/apes 
do
awk -v chr=chr${i} 'BEGIN {OFS="\t"} {if($1==chr)print}' $wd/$negativemask_allchr > $wd/perInterval/chr${i}.${negativemask_allchr}


done

### dogs:
label=dog_canFam3
wd=/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/$label
mkdir $wd/perInterval
negativemask_allchr=${label}.${suffixForMask}
for i in {1..38} # 38 for dogs 
do
awk -v chr=chr${i} 'BEGIN {OFS="\t"} {if($1==chr)print}' $wd/$negativemask_allchr > $wd/perInterval/chr${i}.${negativemask_allchr}


done


### 

############ more tricky: brown bear and minke whale need to be separated into intervals #############
# brown bear: scaffolds gt 1Mb in intervals: 
label=brown_bear_GCF_003584765.1
wd=/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/$label
mkdir $wd/perInterval
negativemask_allchr=${label}.${suffixForMask}
intervaldir=/net/harris/vol1/home/beichman/reference_genomes/brown_bear/genomeIntervalBedFiles_ScaffsGT1Mb # for bears: tells what scaffolds are in each interval
for i in {1..31}   #31 for bb
do
scaffoldbed=${intervaldir}/brown_bear.interval_${i}.bed
scaffs=`awk '{print $1}' $scaffoldbed` # get a list of scaffolds that are in the interval's bed file 
> $wd/perInterval/interval${i}.${negativemask_allchr} # initialize empty file
# go through the scaffolds: 
for scaff in $scaffs
do
awk -v scaff=$scaff 'BEGIN {OFS="\t"} {if($1==scaff)print}' $wd/$negativemask_allchr >> $wd/perInterval/interval${i}.${negativemask_allchr}
done


done

######## minke whale (fin whale) : 
label=minke_whale_GCF_000493695.1_BalAcu1.0
wd=/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/$label
mkdir $wd/perInterval
negativemask_allchr=${label}.${suffixForMask}
intervaldir=/net/harris/vol1/home/beichman/reference_genomes/minke_whale/minke_whale_genome/GCF_000493695.1_BalAcu1.0/contiglist # tells which scaffs are in each interval
for i in {01..96}   #31 for bb
do
contiglist=${intervaldir}/BalAcu1.0_genomic.contiglist_${i}.list
scaffs=`awk '{print $1}' $contiglist` # get a list of scaffolds that are in the interval's bed file 
> $wd/perInterval/interval${i}.${negativemask_allchr} # initialize empty file
# go through the scaffolds: 
for scaff in $scaffs
do
awk -v scaff=$scaff 'BEGIN {OFS="\t"} {if($1==scaff)print}' $wd/$negativemask_allchr >> $wd/perInterval/interval${i}.${negativemask_allchr}
done


done
# vaquita is fine as is. 