#### run config file as a script to set variables ###



####### species specific files and parameters: #####
species=fin_whale


##### flags for dealing with genome splitting by chr or interval ########
interval_or_chr_or_all=interval  # are genomes split into "interval"", "chr" or "allautos" (allautosomes) ; note that chimp vcfs are whole genome but am splitting them up during processing
prepend0=TRUE # TRUE if intervals are counted as 01 02 03 etc. (fin whale); only valid with "interval" (not with chr)
interval_count=96


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
label=${kmersize}_mer_masked # can add other notes here 

########## input vcf files (should be polarized and have non -pol sites removed (except for humans/apes)) ############
vcfdir="/net/harris/vol1/home/beichman/fin_whale/analyses/keightley_polarization/20211020_focalENP_out1MegNov_out2BalMus/vcfs_with_ancestral_alleles/perInterval"
vcffilename="$vcfdir/OnlySitesWithConfidentAncAllelesFromExtSFS.90-10Range.PASSSITESPlusWarnCpGRefSitesONLY.lowqualremoved.JointCalls_f50b4_08_B_VariantFiltration_${interval}.vcf.gz"
vcfNeedsToBeSubsetByChr=FALSE


############# *negative* mask (regions you DON'T want to use in spectrum)  ################
NEGATIVEMASK="/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/minke_whale_GCF_000493695.1_BalAcu1.0/perInterval/interval${interval}.minke_whale_GCF_000493695.1_BalAcu1.0.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed"


############ ancestral fasta info ####################
ancestralFastafilename="/net/harris/vol1/home/beichman/fin_whale/analyses/ancestralReferenceFastas/20211020_focalENP_out1MegNov_out2BalMus/FINAL.MODIFIED.ANCESTRAL.interval_${interval}_finWhale_Ancestral_extSFS_gte90-lte10Range_NOTEthisisMinkeWhaleRefGenome.butAncCallsAreforfinWhaleHumpbackBlueWhale_notforMinke.includesCpGandRepeatSites.fasta"
chrom_pos=0 # 0-based position of chr name in fasta (e.g. >chr1 blahlbah blahblah)

########### individuals to exclude ############
individualsToExclude='BalAcu02,BalMus01,MegNov01,EubGla01,ENPOR12,ENPAK28' # keep this empty if you don't wasnt to exclude any individuals

######## divide inds into pops : ##########
#pops=''
#popFileDir=


####### mutyper variants options ######
passOption=FALSE  # options: TRUE or FALSE. do you want to only select sites that are PASS? TRUE for most species, except for fin whales where I want to keep in the sites previously filtered as CPGs and Repeats (already removed all otehr bad sites)
strictOption=FALSE # options: TRUE OR FALSE. 
# do you want to ignore softmasked sites in ancestral ref fasta? 
# If TRUE you will ignore those sites, if FALSE you will not ignore them. 
# Softmasked sites only meaningful in human ancestral ref fasta (indicate sites that couldn't be polarized). 
# so only want this true for humans
