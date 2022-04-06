#### run config file as a script to set variables ###



####### species specific files and parameters: #####
species=bears


##### flags for dealing with genome splitting by chr or interval ########
interval_or_chr_or_all=interval  # are genomes split into "interval"", "chr" or "allautos" (allautosomes) ; note that chimp vcfs are whole genome but am splitting them up during processing
prepend0=FALSE # TRUE if intervals are counted as 01 02 03 etc. (fin whale); only valid with "interval" (not with chr)
interval_count=31


########## set interval value from SGE_TASK_ID (*code is the same for all species*) ########
if [ $interval_or_chr_or_all = "interval" ]
then
	# prepend 01 02 if needed 
	if [ $prepend0 = "TRUE" ]
	then
		interval=`printf %02d ${SGE_TASK_ID}`
	elif [ $prepend0 = "FALSE" ]
	then
		interval=${SGE_TASK_ID}
	else
		echo "invalid prepend0 option"
		exit 1
	fi
	# set interval label: 
	intervalLabel=interval_${interval}
elif [ $interval_or_chr_or_all = "chr" ]
then
	interval=${SGE_TASK_ID}
	intervalLabel=chr${interval} # label for output files 
elif [ $interval_or_chr_or_all= "allautos" ]
then
	interval='' # no intervals 
	intervalLabel="allautos"
fi


######### kmer size  ########

kmersize=7
label=${kmersize}_mer

########## input vcf files (should be polarized and have non -pol sites removed (except for humans/apes)) ############
vcfdir="/net/harris/vol1/home/beichman/bears/variant_calling/mapped_to_brown_bear/vcfs/vcf_20200916_brown_bear/interval_${interval}/SNPsOnly/phased"
vcffilename="PHASED.mergedSamples.POLARIZED.SNPs.filt_variants_4b.LowCovIndsREMOVED.nomissfiltapplied.PASS.WARN.Only.CustomFilt.HardFilt.AnnotVar.TrimAlt.ref_brown_bear.int_${interval}.withINFO.vcf.gz"
vcfNeedsToBeSubsetByChr=FALSE


############# *negative* mask (regions you DON'T want to use in spectrum)  ################
NEGATIVEMASK="/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/brown_bear_GCF_003584765.1/perInterval/interval${interval}.brown_bear_GCF_003584765.1.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed"


############ ancestral fasta info ####################
ancestralFastafilename="/net/harris/vol1/home/beichman/bears/analyses/ancestralReferenceFastas/brown_bear/interval_${interval}/MODIFIED.ANCESTRAL.brown_bear.interval_${interval}.fasta"
chrom_pos=0 # 0-based position of chr name in fasta (e.g. >chr1 blahlbah blahblah)

########### individuals to exclude ############
individualsToExclude='' # keep this empty if you don't wasnt to exclude any individuals

######## divide inds into pops : ##########
#pops=''
#popFileDir=


####### mutyper variants options ######
passOption=TRUE  # options: TRUE or FALSE. do you want to only select sites that are PASS? TRUE for most species, except for fin whales

strictOption=FALSE # options: TRUE OR FALSE. 
# do you want to ignore softmasked sites in ancestral ref fasta? 
# If TRUE you will ignore those sites, if FALSE you will not ignore them. 
# Softmasked sites only meaningful in human ancestral ref fasta (indicate sites that couldn't be polarized). 
# so only want this true for humans
