#### run config file as a script to set variables ###



####### species specific files and parameters: #####
species=mice


##### flags for dealing with genome splitting by chr or interval ########
interval_or_chr_or_all=chr  # are genomes split into "interval"", "chr" or "allautos" (allautosomes) ; note that chimp vcfs are whole genome but am splitting them up during processing
prepend0=FALSE # TRUE if intervals are counted as 01 02 03 etc. (fin whale); only valid with "interval" (not with chr)
interval_count=19


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
vcfdir="/net/harris/vol1/home/beichman/mice/analyses/keightley_polarization/20210225_focalMmd_out1Mmm_out2Ms/vcfs_with_ancestral_alleles/perChromosome"
vcffilename="chr${interval}.OnlySitesWithConfidentAncAllelesFromExtSFS.90-10Range.AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz"
vcfNeedsToBeSubsetByChr=FALSE


############# *negative* mask (regions you DON'T want to use in spectrum)  ################
NEGATIVEMASK="/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/mouse_mm10/perInterval/chr${interval}.mouse_mm10.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed" 
maskLabel=maskALL


############ ancestral fasta info ####################
# make sure doesn't contain non-autosome chrs/scaffs that weren't used in genotype calling (i have confirmed for all my species)
ancestralFastaDir="/net/harris/vol1/home/beichman/mice/analyses/ancestralReferenceFastas/20210225_focalMmd_out1Mmm_out2Ms/"
ancestralFastafilename="FINAL.MODIFIED.ANCESTRAL.chr${interval}_mouseAncestral_extSFS_gte90-lte10Range.fasta"
chrom_pos=0 # 0-based position of chr name in fasta (e.g. >chr1 blahlbah blahblah)

########### individuals to exclude ############
individualsToExclude='' # keep this empty if you don't wasnt to exclude any individuals

######## divide inds into pops : ##########
pops='Mmd Mmc Mmm Ms'
poplistdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/BearAdmixtureProject/samples/mouseSamples/populationFilesForMutyper
popfilesuffix=".txt" # format pop.txt


####### mutyper variants options ######
passOption=TRUE  # options: TRUE or FALSE. do you want to only select sites that are PASS? TRUE for most species, except for fin whales

strictOption=FALSE # options: TRUE OR FALSE. should only be true for humans.
# do you want to ignore softmasked sites in ancestral ref fasta? 
# If TRUE you will ignore those sites, if FALSE you will not ignore them. 
# Softmasked sites only meaningful in human ancestral ref fasta (indicate sites that couldn't be polarized). 
# so only want this true for humans
