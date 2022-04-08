#### run config file as a script to set variables ###



####### species specific files and parameters: #####
species=Gorilla_gorilla


##### flags for dealing with genome splitting by chr or interval ########
interval_or_chr_or_all=chr  # are genomes split into "interval"", "chr" or "allautos" (allautosomes) ; note that chimp vcfs are whole genome but am splitting them up during processing
prepend0=FALSE # TRUE if intervals are counted as 01 02 03 etc. (fin whale); only valid with "interval" (not with chr)
interval_count=22


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
vcfdir="/net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs"
vcffilename="Gorilla.vcf.gz" # note this whole genome; but Michael's anc fasta are per chrom so i am subsetting this vcf into chrs when it's processed
vcfNeedsToBeSubsetByChr=TRUE # because it needs to be subset during processing at mutyper variants step (apes only)


############# *negative* mask (regions you DON'T want to use in spectrum)  ################
NEGATIVEMASK="/net/harris/vol1/home/beichman/reference_genomes/unifiedBedMasksForAllGenomes/apes_mapped_to_hg18/perInterval/chr${interval}.apes_mapped_to_hg18.exon10kb.repmask.trf.cpgIslands.NEGATIVEMASK.merged.USETHIS.bed" # ONE MASK FILE PER INTERVAL
# note: apes cannot use same negative mask as humans because are mapped to diff ref genomes (humans: hg38; apes: hg18 old), so have made a different hg18 specific mask for the apes
maskLabel=maskALL


############ ancestral fasta info ####################
# make sure doesn't contain non-autosome chrs/scaffs that weren't used in genotype calling (i have confirmed for all my species)
ancestralFastaDir="/net/harris/vol1/home/beichman/apes/polarized_ref_fastas/hg18_references_with_gagp_ancestral_alleles_exclude_recurrent/"
ancestralFastafilename="Gorilla_chr${interval}.fa" # made by Michael
chrom_pos=0 # 0-based position of chr name in fasta (e.g. >chr1 blahlbah blahblah)

########### individuals to exclude ############
individualsToExclude='Gorilla_gorilla_gorilla-X00108_Abe,Gorilla_gorilla_gorilla-KB7973_Porta,Gorilla_beringei_graueri-9732_Mkubwa,Gorilla_beringei_graueri-A929_Kaisi,Gorilla_beringei_graueri-Victoria,Gorilla_gorilla_dielhi-B646_Nyango'  
# excluded due to bad qual (first 2) and due to being berigei or dielhi 

######## divide inds into pops : ##########
pops='' # have it blank if no pops to split into and don't need to define other terms



####### mutyper variants options ######
passOption=TRUE  # options: TRUE or FALSE. do you want to only select sites that are PASS? TRUE for most species, except for fin whales

strictOption=FALSE # options: TRUE OR FALSE. ; should only be true for humans
# do you want to ignore softmasked sites in ancestral ref fasta? 
# If TRUE you will ignore those sites, if FALSE you will not ignore them. 
# Softmasked sites only meaningful in human ancestral ref fasta (indicate sites that couldn't be polarized). 
# so only want this true for humans

