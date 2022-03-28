############# gather all spectra and calculate rescaled mutation rates  ###########
require(tidyverse)
require(reshape2)
require(dplyr)
require(ggplot2)
require(tidyr)
require(compositions)  # for clr
require(broom) # for tidy
require(ggrepel)
require(ggfortify) # for pca autoplot
require(gridExtra)
require(sjPlot) # for grids of plots
require(ggpmisc) # for nicer plotting of linear models
require(ape)
require(sigfit)
require(spgs)  # for reverse complements
require(gtools) # for mixedsort
require(RColorBrewer)
set.seed(42) # setting a seed from now on (20220309+)
########### NOTE FOR SIGFIT: 
# COUNTS, SIGNATURES, TARGETS MUST BE IN SAME ORDER. NOT SUFFICIENT TO HAVE SAME NAMES
# SO IF WORKING WITH TRIPLETS YOU MUST GET INTO SAME ORDER AS COSMIC IF YOU WANT TO COMPARE TO COSMIC 
niter=10000 #default is 2000 ; but want more 
model="multinomial" # (NMF)
mantelTestPermutations=99999 # up this on 20220207 

separateCpGs="yes" # yes or no
todaysdate=format(Sys.Date(),"%Y%m%d")
flagOfAnalysis=paste0(todaysdate,".sigfitResults.RescaleToHumanGenomeTargets.highAdaptDelta.YoungParent13.plusSoloHumans.",model) # a label that will be in title of plots and outdir 
plotdir=paste0("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/sigfit/prepSpectraForSigfit/",flagOfAnalysis,"/")
dir.create(plotdir,showWarnings = F,recursive = T)

####### signatures ##########
# OLD ; used through 20220209: used non-exponentiated exponents so young parent signature had v few c>gs
# agingSignatures <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/jonsson_parental_aging/jonsson.agingsignatures.tables9.convertedtoproportionsinexcel.REVCOMP.REORDEREDTOMATCHMYDATA.txt",header=T,sep="\t",row.names = "signature") # USE VERSION THAT HAS BEEN REVCOMPED 

##### NEW: young_parent_13 signature calculated from age 13 instead of intercepts:
agingSignatures = read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/jonsson_parental_aging/jonsson.agingsignatures.tables9.convertedtoproportionsinexcel.REVCOMP.REORDEREDTOMATCHMYDATA.YOUNGPARENTCALCAT13.notintercept.txt",header=T,sep="\t",row.names="signature")

agingSignatures ####### START HERE TOMORROW: 
########## test:adding a made up cpg signature ##########
CpGSignature = data.frame(row.names = "CpG.TpG signature",A.C=0,A.G=0,	A.T=0,	C.A=0	,C.G=0,	C.T=0,	C.T_CpG=1.0) # this signature is just CpG>T (100% of mutations from this signature)
# Can I use it to sop up CpGs? and then look more closely at other signatures?
agingPlusCpGSignatures <- bind_rows(CpGSignature,agingSignatures)

####### cosmic signatures #########
data(cosmic_signatures_v3.2)
cosmic_signatures_v3.2
# keeping these as human targets since I'm rescaling everything to human. but note if I return to using 'opportunities' then I need to rescale the cosmic signatures away from human targets. are fine for now though. (see github sigfit notes for details)
########## table of target locations ########

tableOfTargetsToInclude=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/spectra/tableOfSpectraToIncludeInComparisons.SpeciesWithMultipleIndividualsOnly.PlusApes.MaskedMouse.Vaquita.20211123.txt",header=T,sep="\t")  # UPDATED this file; it now has bedmasked vaquita targets and repmasked mouse targets

############## all projected spectra #############
# okay have all the spectra together already in:
# THESE have been projected down to k haploids:
#allProjectedSpectra <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/allKSFSes_and_projected_Spectra/All_SpectrumCounts_ProjectedDownTo.16.Haploids.txt",header=T)
allProjectedSpectra <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/allKSFSes_and_projected_Spectra/20211123_All_SpectrumCounts_ProjectedDownTo.10.Haploids.includesApes.RepMaskedMice.txt",header=T) # updated this file 20211123 to include repmasked mouse 
# this now includes bears_EUR and all apes (projectd down to 5 diploids, 10 haps)
speciesCodes=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/scripts/DNA_Shape_Project/analyses/PhylogeneticDistances/speciesCodes.txt",header=T) # for now just using AFR as human ; added apes
head(speciesCodes) # for getting phylo distances 

######  get phylo distances from upham et al Raxml tree:  #########
# used to read these in, now I process it her : 
#phyloDistances = read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/raxml_tree/raxml.cophenetic.pairwisedistances.fromUphametal.MYSPECIES.usethis.txt",header=T) # added apes
#head(phyloDistances) # not all spp are present
# used to get these distances in a separate script. now doing it here. 

# process raxml tree here instead of in previous 
raxmlTree=ape::read.tree("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/raxml_tree/RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick")

# want to get rid of the family/order info at end of tip. just keep species : 
raxmlTree_renamedTips <- raxmlTree
raxmlTree_renamedTips$tip.label <- paste0(unlist(lapply(strsplit(raxmlTree_renamedTips$tip.label,"_"),"[",1)),"_",unlist(lapply(strsplit(raxmlTree_renamedTips$tip.label,"_"),"[",2)))

speciesToInclude_in_raxmlTree=speciesList=c("Balaenoptera_physalus" ,"Phocoena_sinus" ,"Ursus_maritimus" ,"Ursus_arctos" ,"Mus_musculus" ,"Mus_spretus","Homo_sapiens","Pongo_pygmaeus","Pongo_abelii","Pan_troglodytes","Pan_paniscus","Gorilla_gorilla")  # raxml name fmt
# adding in apes 
# restrict to just the species to include
raxmlTree_renamedTips_subset <- keep.tip(raxmlTree_renamedTips,speciesToInclude_in_raxmlTree)

raxml_cophenetic_dist <- cophenetic(raxmlTree_renamedTips_subset)
# save this as an object

melt_cophenetic_func <- function(cophenetic_dist){
  cophenetic_dist_melt <- melt(cophenetic_dist)
  colnames(cophenetic_dist_melt) <- c("Sp1","Sp2","cophenetic_distance")
  # get species names # ASSUMES ARE SEPARATED BY "_" eg mus_musculus
  cophenetic_dist_melt$Sp1_species <- paste0(unlist(lapply(strsplit(as.character(cophenetic_dist_melt$Sp1),"_"),"[",1)),"_",unlist(lapply(strsplit(as.character(cophenetic_dist_melt$Sp1),"_"),"[",2)))
  
  cophenetic_dist_melt$Sp2_species <- paste0(unlist(lapply(strsplit(as.character(cophenetic_dist_melt$Sp2),"_"),"[",1)),"_",unlist(lapply(strsplit(as.character(cophenetic_dist_melt$Sp2),"_"),"[",2)))
  
  # need to deal with reciprocal duplicates
  # put in alphabetical order so I can get rid of dups:
  cophenetic_dist_melt$comparisonLabel <- paste0(pmin(cophenetic_dist_melt$Sp1_species,cophenetic_dist_melt$Sp2_species),".",pmax(cophenetic_dist_melt$Sp1_species,cophenetic_dist_melt$Sp2_species))  # # need to use pmin and pmax to get row by row min max. otherwise gives min max of whoel vector
  # get rid of self to self comparisons:
  cophenetic_dist_melt_distinct <- cophenetic_dist_melt[cophenetic_dist_melt$Sp1!=cophenetic_dist_melt$Sp2,]
  # and get rid of reciprocal dups (sp A - B and sp B- A)
  cophenetic_dist_melt_distinct <- cophenetic_dist_melt_distinct %>%
    ungroup() %>% # adding ungroup as a safety measure on 20211201 not sure if ottally necessary though
    distinct(comparisonLabel,.keep_all=T)  # restrict to just distinct comparison labels. 
  dim(cophenetic_dist_melt_distinct)
  return(cophenetic_dist_melt_distinct)
  
  
}

phyloDistances = melt_cophenetic_func(raxml_cophenetic_dist)


######### 1mers are dealt with in the if statement below this (dealing with CpGs) ########
# get ancestral kmers:
allProjectedSpectra$ancestral7mer <- substr(allProjectedSpectra$variable,1,7)
allProjectedSpectra$ancestral5mer <- substr(allProjectedSpectra$variable,2,6)
allProjectedSpectra$ancestral3mer <- substr(allProjectedSpectra$variable,3,5)

# get central mutation types:
allProjectedSpectra$mutation_7mer <- allProjectedSpectra$variable
allProjectedSpectra$mutation_5mer <- paste0(substr(allProjectedSpectra$variable,2,6),".",substr(allProjectedSpectra$variable,10,14))
allProjectedSpectra$mutation_3mer <- paste0(substr(allProjectedSpectra$variable,3,5),".",substr(allProjectedSpectra$variable,11,13))

# want to get ancestral 1mers for all data frames so I can plot things with central mutattion type separated later:

######## SEPARATE OUT  -- OPTIONAL ############
# just want to separate out C>Ts not other CpG types here (different in other script)
if(separateCpGs=="yes"){
  #### label CpGs and specify them as ancestral target as well
  print("separating out CpG sites!")
  allProjectedSpectra$mutation_1mer_CpGNotLabeled <- paste0(substr(allProjectedSpectra$variable,4,4),".",substr(allProjectedSpectra$variable,12,12))
  allProjectedSpectra$ancestral1mer_CpGNotLabeled <- paste0(substr(allProjectedSpectra$variable,4,4)) # adding this on 20220112 
  #### label CpGs and specify them as ancestral target as well
  allProjectedSpectra$CpGLabel <- ""
  # find all CpG sites:
  allProjectedSpectra[grepl("CG",allProjectedSpectra$ancestral3mer) & allProjectedSpectra$mutation_1mer_CpGNotLabeled=="C.T",]$CpGLabel <- "_CpG" # on 20220112 am changing this to just do C.T mtuations 
  allProjectedSpectra$mutation_1mer <- paste0(allProjectedSpectra$mutation_1mer_CpGNotLabeled,allProjectedSpectra$CpGLabel)
  # label ancestral state with CpG as well
  allProjectedSpectra$ancestral1mer <- allProjectedSpectra$ancestral1mer_CpGNotLabeled
  allProjectedSpectra[allProjectedSpectra$CpGLabel=="_CpG",]$ancestral1mer <- "C_CpG"
  
  allProjectedSpectra$mutation_1mer <- paste0(allProjectedSpectra$mutation_1mer_CpGNotLabeled,allProjectedSpectra$CpGLabel)
  
  allProjectedSpectra$ancestral1mer_CpGNotLabeled <- substr(allProjectedSpectra$variable,4,4) # not yet CpG identified 
  
  
} else if(separateCpGs=="no"){
  print("note that CpGs are not separated out!")
  allProjectedSpectra$mutation_1mer <- paste0(substr(allProjectedSpectra$variable,4,4),".",substr(allProjectedSpectra$variable,12,12)) # if you don't want to separate out CpGs then this is just the regular center bp
  # still get ancestral1mer (but don't specify if CpG or not)
  allProjectedSpectra$ancestral1mer <- substr(allProjectedSpectra$variable,4,4) # if you don't want to separate out CpGs then this is just the regular center bp
  
  
  
} else {
  print("invalid choice for separateCpGs")
  break
}



# get spectrum summed up per type
all7merSpectraOnly <- allProjectedSpectra %>%
  # adding groupbing by central 1mer but that doesn't actually change the sums since 7mers already group by those 
  # I just want to keep it in the dataframe 
  group_by(species,population,label,ancestral7mer,mutation_7mer,ancestral1mer,mutation_1mer) %>%
  summarise(total_mutations_projected = sum(totalSites_projectedDown_allFreqsSummed)) # this will be identical to original spectrum (7mer), just doing for consistency of formatting

all5merSpectraOnly <- allProjectedSpectra %>%
  group_by(species,population,label,ancestral5mer,mutation_5mer,ancestral1mer,mutation_1mer) %>%
  summarise(total_mutations_projected = sum(totalSites_projectedDown_allFreqsSummed)) 

all3merSpectraOnly <- allProjectedSpectra %>%
  group_by(species,population,label,ancestral3mer,mutation_3mer,ancestral1mer,mutation_1mer) %>%
  summarise(total_mutations_projected = sum(totalSites_projectedDown_allFreqsSummed)) 

all1merSpectraOnly <- allProjectedSpectra %>%
  group_by(species,population,label,ancestral1mer,mutation_1mer) %>%
  summarise(total_mutations_projected = sum(totalSites_projectedDown_allFreqsSummed)) 
# if you don't want to partition over CpGs you can use mutation_1mer_CpGNotLabeled instead
# okay NOTE HERE: "C" counts are nonCpG and then C_CpG are CpG ancestral targets


########## read in targets ##############
all7merTargets=data.frame()
all5merTargets=data.frame()
all3merTargets = data.frame()
all1merTargets = data.frame()

for(sample in tableOfTargetsToInclude$sample){
  # some 'samples' may contain multipops like mice and bears
  sampleinfo = tableOfTargetsToInclude[tableOfTargetsToInclude$sample==sample,]
  targetsfilename = paste0(sampleinfo$targetsindir,"/",sampleinfo$targetsFile)
  
  targets=read.table(targetsfilename,header=sampleinfo$targetsHeader) # assigns whether has header or not
  colnames(targets) <- c("target_7mer","countOverAllIntervals")
  
  # get central 5mer, 3mer and 1mer:
  targets <- targets %>%
    mutate(target_5mer=substr(target_7mer,2,6),target_3mer=substr(target_7mer,3,5),target_1mer= substr(target_7mer,4,4)) # get substrings 
  
  ######### relabel target_1mer with CpGs if you are separating out CpGs ###########
  if(separateCpGs=="yes"){
    targets$CpGLabel <- ""
    targets[grepl("CG",targets$target_3mer),]$CpGLabel <- "_CpG"
    targets$target_1mer_CpGsNotLabeled <- targets$target_1mer # update with this 
    targets$target_1mer <- paste0(targets$target_1mer_CpGsNotLabeled,targets$CpGLabel)
    
  } # don't need an else statement, just keep going with targets target$target_1mer as is.
  
  
  
  # TARGETS per 3mer type #
  targets_7mer <- targets %>% group_by(target_7mer) %>% summarise(total_target_count=sum(countOverAllIntervals)) # original targets file, just doing to be consistent but it's identical to original targets
  targets_3mer <- targets %>% group_by(target_3mer) %>% summarise(total_target_count=sum(countOverAllIntervals))
  targets_5mer <- targets %>% group_by(target_5mer) %>% summarise(total_target_count=sum(countOverAllIntervals))
  targets_1mer <- targets %>% group_by(target_1mer) %>% summarise(total_target_count=sum(countOverAllIntervals))
  targets_7mer$species <- sample
  targets_5mer$species <- sample
  targets_3mer$species <- sample
  targets_1mer$species <- sample
  # make all species targets:
  
  
  # combine  all targets:
  all7merTargets = bind_rows(all7merTargets,targets_7mer)
  all5merTargets = bind_rows(all5merTargets,targets_5mer)
  all3merTargets = bind_rows(all3merTargets,targets_3mer)
  all1merTargets = bind_rows(all1merTargets,targets_1mer)
  
  
}
# plot targets:
target1merPlot <- ggplot(all1merTargets,aes(y=species,x=total_target_count,fill=target_1mer))+
  geom_col(position="dodge")

########## get target proportions (will be added to merged spectra ) ########
all1merTargets <- all1merTargets %>%
  group_by(species) %>%
  mutate(target_proportion = total_target_count / sum(total_target_count)) %>%
  ungroup()


all3merTargets <- all3merTargets %>%
  group_by(species) %>%
  mutate(target_proportion = total_target_count / sum(total_target_count)) %>%
  ungroup()

all5merTargets <- all5merTargets %>%
  group_by(species) %>%
  mutate(target_proportion = total_target_count / sum(total_target_count)) %>%
  ungroup()


all7merTargets <- all7merTargets %>%
  group_by(species) %>%
  mutate(target_proportion = total_target_count / sum(total_target_count)) %>%
  ungroup()

target1merPlot2 <-   ggplot(all1merTargets,aes(y=species,x=target_proportion,fill=target_1mer))+geom_col()+facet_wrap(~target_1mer,scales="free")+ggtitle("target proportions per genome")
target1merPlot2

ggsave(paste0(plotdir,"targetPlot.1mers.png"),target1merPlot,height=9,width=5)
ggsave(paste0(plotdir,"targetPlot.1mers.proportions.png"),target1merPlot2,height=9,width=12)

#ggplot(all1merSpectraOnly,aes(y=species,x=total_mutations_projected,fill=mutation_1mer))+
#  geom_col(position="dodge")
# complete so that missing mutation types are filled in 
# okay this is kind of sinister, depending on how it was grouped above this may or may not work (but doesn't throw an error. Frustrating!)
# so put in a check
# keep an eye on this in other code. Don't think it's made things go wrong before but be alert.
all7merSpectraOnly_filledin <- all7merSpectraOnly %>%
  ungroup() %>%
  complete(nesting(mutation_7mer,ancestral7mer,ancestral1mer,mutation_1mer),nesting(species,population,label),fill=list(total_mutations_projected=0))%>%
  mutate(mutation_label=mutation_7mer)

if(dim(all7merSpectraOnly_filledin)[1]==dim(all7merSpectraOnly)[1]){
  print("filling in didn't work!")
}
all5merSpectraOnly_filledin <- all5merSpectraOnly %>%
  ungroup() %>% # need to pre-ungroup ; got grouped up above
  complete(nesting(mutation_5mer,ancestral5mer,ancestral1mer,mutation_1mer),nesting(species,population,label),fill=list(total_mutations_projected=0))%>%
  mutate(mutation_label=mutation_5mer) # adding mutation label for use in functions

all3merSpectraOnly_filledin <- all3merSpectraOnly %>%
  ungroup() %>%
  complete(nesting(mutation_3mer,ancestral3mer,ancestral1mer,mutation_1mer),nesting(species,population,label),fill=list(total_mutations_projected=0))%>%
  mutate(mutation_label=mutation_3mer)

all1merSpectraOnly_filledin <- all1merSpectraOnly %>%
  ungroup() %>%
  complete(nesting(mutation_1mer,ancestral1mer),nesting(species,population,label),fill=list(total_mutations_projected=0)) %>%
  mutate(mutation_label=mutation_1mer)
# fill in should only be needed for 7mer (won't change dim of others)

# so when I was merging before all species were combined there was a bug where missing mtuation types wouldn't get targets filled in which led to them being NA downstream 
# okay now fill in missing mutation types : as 0s PRIOR TO MERGING WITH TARGEST
#  note that targets are the same within species (mouse1 and mouse2 have same targets; bear 1 and bear2 have targets)

#### to merge with ancestral targets, need apes species not to be 'apes' (but still want label to be apes_Species)
# so let's just update them here: # making the species equal to the popualtion (which is where I stored ape species. not ideal from a gneralizing point of view but works fine)
############### fill in ape species info
all7merSpectraOnly_filledin[all7merSpectraOnly_filledin$species=="apes",]$species <- all7merSpectraOnly_filledin[all7merSpectraOnly_filledin$species=="apes",]$population


all5merSpectraOnly_filledin[all5merSpectraOnly_filledin$species=="apes",]$species <- all5merSpectraOnly_filledin[all5merSpectraOnly_filledin$species=="apes",]$population


all3merSpectraOnly_filledin[all3merSpectraOnly_filledin$species=="apes",]$species <- all3merSpectraOnly_filledin[all3merSpectraOnly_filledin$species=="apes",]$population

all1merSpectraOnly_filledin[all1merSpectraOnly_filledin$species=="apes",]$species <- all1merSpectraOnly_filledin[all1merSpectraOnly_filledin$species=="apes",]$population

############## merge with targets ###############

all7merSpectra <- merge(all7merSpectraOnly_filledin, all7merTargets,by.x=c("species","ancestral7mer"),by.y=c("species","target_7mer"))

if(dim(all7merSpectra)[1] != dim(all7merSpectraOnly_filledin)[1]){
  print("something went wrong with merge!")
  
}


all5merSpectra <- merge(all5merSpectraOnly_filledin, all5merTargets,by.x=c("species","ancestral5mer"),by.y=c("species","target_5mer"))
all3merSpectra <- merge(all3merSpectraOnly_filledin, all3merTargets,by.x=c("species","ancestral3mer"),by.y=c("species","target_3mer"))
all1merSpectra <- merge(all1merSpectraOnly_filledin, all1merTargets,by.x=c("species","ancestral1mer"),by.y=c("species","target_1mer"))
# these now include target proportion


##################### rescale by human target proportions ############
humanTargetProportions_1mer <- all1merTargets[all1merTargets$species=="humans",]
humanTargetProportions_3mer <- all3merTargets[all3merTargets$species=="humans",]
humanTargetProportions_5mer <- all5merTargets[all5merTargets$species=="humans",]
humanTargetProportions_7mer <- all7merTargets[all7merTargets$species=="humans",]

# how to do this : count(vaquita) * [ (target proportion in human)/(target proportion in vaquita)], where target proportion is target i / total targets. So for instance, if humans have proportionally more of target i, then human/vaqu ratio will be >1 and then your count will go up.
# whereas if they have proportionally less, count will go down.
# can also do this with proportion(vaquita) instead of counts -- will work the same way,where prop(vaq) is frac seg sites
# I'm doing it with projected counts. Then still going to multinomial downsample it (??)
# maybe not? 

# so need to get target ratio for every species (done)
# and divide human by each species target ratio <- how?
# then multiply across
# okay I think I can do that

######## >>>> use these for spectrum distances at bottom of script <<<< ##############
# note these don't yet have the human counts; need to run processSpectra_AdjustForHumanTargets on them (is done in combo distance script at bottom fo the script )
all1merSpectra_MERGEDWITHUMANTARGETPROPORTIONS <- merge(all1merSpectra,humanTargetProportions_1mer[,c("target_1mer","target_proportion")],by.x=c("ancestral1mer"),by.y=c("target_1mer"),suffixes=c("_thisSpecies","_HUMAN"))

# during processing will rescale  
all3merSpectra_MERGEDWITHUMANTARGETPROPORTIONS <- merge(all3merSpectra,humanTargetProportions_3mer[,c("target_3mer","target_proportion")],by.x=c("ancestral3mer"),by.y=c("target_3mer"),suffixes=c("_thisSpecies","_HUMAN"))

all5merSpectra_MERGEDWITHUMANTARGETPROPORTIONS <- merge(all5merSpectra,humanTargetProportions_5mer[,c("target_5mer","target_proportion")],by.x=c("ancestral5mer"),by.y=c("target_5mer"),suffixes=c("_thisSpecies","_HUMAN"))

all7merSpectra_MERGEDWITHUMANTARGETPROPORTIONS <- merge(all7merSpectra,humanTargetProportions_7mer[,c("target_7mer","target_proportion")],by.x=c("ancestral7mer"),by.y=c("target_7mer"),suffixes=c("_thisSpecies","_HUMAN"))


################ multinomial downsample during process spectra ######################


# add multinomial variables to this function.
# new version that deals with new counts that have been rescaled for human targets 
# doesn't do all the epsilon stuff (don't need for sigfit )

processSpectra_AdjustForHumanTargets <- function(spectradf,wd){
  
  spectradf <- spectradf %>%
    group_by(species,population,label) %>%
    mutate(total_mutations_projected_RESCALED_BY_HUMAN_TARGETS = total_mutations_projected* (target_proportion_HUMAN /target_proportion_thisSpecies),frac_seg_sites_RESCALED_BY_HUMAN_TARGETS= total_mutations_projected_RESCALED_BY_HUMAN_TARGETS/sum(total_mutations_projected_RESCALED_BY_HUMAN_TARGETS)) # doing this based on non-rescaled by human version 
  # note that this leads to differing numbers of mutations due to rescaling of CpG targets etc.
  
  # add multinomial sampling:
  # calculate min number of seg sites:
  totalSegSites <- spectradf %>%
    group_by(species,population,label) %>%
    summarise(totalSegSites_RESCALED_BY_HUMAN_TARGETS=sum(total_mutations_projected_RESCALED_BY_HUMAN_TARGETS))
  
  subsampleValue=min(totalSegSites$totalSegSites_RESCALED_BY_HUMAN_TARGETS)
  
  subsampleSpecies=totalSegSites[totalSegSites$totalSegSites_RESCALED_BY_HUMAN_TARGETS==min(totalSegSites$totalSegSites_RESCALED_BY_HUMAN_TARGETS),]$species
  
  multinomlogfile <- file(paste0(wd,"multinomialLogFile.txt"))
  
  writeLines(paste0("downsampling to ",subsampleValue," sites (",subsampleSpecies,")"),paste0(plotdir,"multinomialDownsampling.info.txt"))
  close(multinomlogfile)
  
  print(paste0("downsampling to ",subsampleValue," sites (",subsampleSpecies,")"))
        # carry out multinomial sampling:
        
  spectradf <- spectradf %>%
        group_by(species,label,population) %>%
        mutate(total_mutations_projected_RESCALED_BY_HUMAN_TARGETS_multinom_downsampled=as.numeric(rmultinom(n=1,size=subsampleValue,prob=frac_seg_sites_RESCALED_BY_HUMAN_TARGETS))) %>%  # n=1 is how many times to sample (reps); size is how many sites (~230,000), probs is vector of probabilities aka the raw fraction (no +epsilon, no target correction)
          ungroup()
  
  return(spectradf)
        
}

pivotSpectra_perVariable <- function(processedspectradf,variableToPivot) {
  processedspectradf_wider <-  pivot_wider(processedspectradf,id_cols=c(mutation_label),names_from=c("label"),values_from =all_of(variableToPivot ))
  processedspectradf_wider$variable <- variableToPivot
  # pivots on a variable 
  return(processedspectradf_wider)
}


################# process spectra #################
# for 1mers this doesn't matter, but for 3mers need to rev comp ala SPE for sigfit to compare to cosmic signatures!!! (pull code from ProjSpectra...script for writing out for SPE)
# also need to write out processed spectra and oppos here
# need to get opportunities (oppos) to work -- how do they get matched? is ">" essential? does it work for non triplets? 
spectrum_1mer_processed = processSpectra_AdjustForHumanTargets(all1merSpectra_MERGEDWITHUMANTARGETPROPORTIONS,plotdir)

spectrum_3mer_processed = processSpectra_AdjustForHumanTargets(all3merSpectra_MERGEDWITHUMANTARGETPROPORTIONS,plotdir)

######### need to add rev comp information to 3mers (don't need for 1mers) to relate to cosmic signatures ##########
# if ancestral1mer is C it's good to go. but if it's not then you need to revcomp
spectrum_3mer_processed$mutation_label_for_sigfit <-  spectrum_3mer_processed$mutation_label
# get derived 3mer
spectrum_3mer_processed$derived3mer <- substr(spectrum_3mer_processed$mutation_label,5,7)
# need to revcomp any that
spectrum_3mer_processed[spectrum_3mer_processed$ancestral1mer=="A",]$mutation_label_for_sigfit <- paste0(toupper(reverseComplement(spectrum_3mer_processed[spectrum_3mer_processed$ancestral1mer=="A",]$ancestral3mer)),".",toupper(reverseComplement(spectrum_3mer_processed[spectrum_3mer_processed$ancestral1mer=="A",]$derived3mer)))
# SO note that this relabels 3mers that have "A" as central bp to their rev comps so it's in the same format as cosmic 
# you can check with:
sum(!gsub("\\.",">",spectrum_3mer_processed$mutation_label_for_sigfit) %in% colnames(cosmic_signatures_v3.2)) # should be 0 -- checks for any names that aren't in the cosmic names once periods subbed for >s. Yay they are all in there. 
# note that rescaling by human targets changes the total number of SNPs 
# doesn't just redistribute because of different rates of different types. CpGs make the biggest change due to their high rate
# so for instance if you have more of a high rate target in one species, but human has less
# proportionally get redistrubted so that you observe less of that but that has a big impact on the number of mutations observed (fewer)

# okay now want to format for sigfit:
# want sample column and then want mtuaiton labels as column names

######## get counts for sigfit:  ###########

spectrum_1mer_forSigfit_downsampledCounts <- pivot_wider(spectrum_1mer_processed[,c("label","total_mutations_projected_RESCALED_BY_HUMAN_TARGETS_multinom_downsampled","mutation_label")],id_cols=c("label"),names_from = "mutation_label",values_from ="total_mutations_projected_RESCALED_BY_HUMAN_TARGETS_multinom_downsampled" )
########## plot spectra #########
plot1merspectra <- ggplot(melt(spectrum_1mer_forSigfit_downsampledCounts),aes(y=label,x=value,fill=variable))+
  geom_col()+
  facet_wrap(~variable,scales="free")+
  ggtitle("rescaled to human targets")
ggsave(paste0(plotdir,"plot1merspectra.rescaledToHumanTargets.png"),plot1merspectra,height=7,width=12)

# note for 3mers  need to use sigfit labelled 3mers so they match cosmic signatures 
spectrum_3mer_forSigfit_downsampledCounts <-  pivot_wider(spectrum_3mer_processed[,c("label","total_mutations_projected_RESCALED_BY_HUMAN_TARGETS_multinom_downsampled","mutation_label_for_sigfit")],id_cols=c("label"),names_from = "mutation_label_for_sigfit",values_from ="total_mutations_projected_RESCALED_BY_HUMAN_TARGETS_multinom_downsampled" )

# 220921 projected down to this (note fewer than 230K from before due to human target rescaling)
######### okay so a bug was introduced here: for some reason after this rescaling, the types were reordered so they no longer match signatures: fixed now.

# change 'label' column name to 'sample' for sigfit format 
# rename C.T_CpG to CpG.TpG
############ move label column to rownames (important ) #############
spectrum_1mer_forSigfit_downsampledCounts <- spectrum_1mer_forSigfit_downsampledCounts %>%
  remove_rownames %>% # have to remove any preexisting rownames 
  column_to_rownames(var="label") # adds row names

spectrum_3mer_forSigfit_downsampledCounts <- spectrum_3mer_forSigfit_downsampledCounts %>%
  remove_rownames %>% # have to remove any preexisting rownames 
  column_to_rownames(var="label") # adds row names

# reorder to match signatures (can do this with cosmic as well once I'm doing triplets)
spectrum_1mer_forSigfit_downsampledCounts <- spectrum_1mer_forSigfit_downsampledCounts[,names(agingSignatures)]

write.table(spectrum_1mer_forSigfit_downsampledCounts,paste0(plotdir,"spectrum_1mer.forsigfit.downsampledCounts.RESCALEDTOHUMANTARGETS.txt"),row.names=T,sep="\t",quote=F)

# reorder 3mers to be in same order as cosmic signatures: 

# need to add ">" to my data intead of . so i can order it like cosmic
# also need to rev comp (didn't need to revcomp 1mers because wasn't dealing with cosmic signatures and already rev comped the aging signatures)
colnames(spectrum_3mer_forSigfit_downsampledCounts) <- gsub("\\.",">",colnames(spectrum_3mer_forSigfit_downsampledCounts))
head(spectrum_3mer_forSigfit_downsampledCounts)

spectrum_3mer_forSigfit_downsampledCounts <- spectrum_3mer_forSigfit_downsampledCounts[,colnames(cosmic_signatures_v3.2)]

write.table(spectrum_3mer_forSigfit_downsampledCounts,paste0(plotdir,"spectrum_3mer.forsigfit.downsampledCounts.RESCALEDTOHUMANTARGETS.revcompedlabelstomatchcosmic.txt"),row.names=T,sep="\t",quote=F)

############### get human-only with and without downsampling #########
# just need 1mer and 3mer for sigfit
humans=c("humans_AFR","humans_AMR","humans_EAS","humans_SAS","humans_EUR")

# with downsampling:
spectrum_1mer_forSigfit_downsampledCounts_humansOnly <- spectrum_1mer_forSigfit_downsampledCounts[humans,]

spectrum_3mer_forSigfit_downsampledCounts_humansOnly <- spectrum_3mer_forSigfit_downsampledCounts[humans,]


# withOUT downsampling (doesn't need human target rescaling bc is already mapped to same genome)
spectrum_1mer_forSigfit_NONdownsampledCounts <- pivot_wider(spectrum_1mer_processed[,c("label","total_mutations_projected","mutation_label")],id_cols=c("label"),names_from = "mutation_label",values_from ="total_mutations_projected" )

spectrum_1mer_forSigfit_NONdownsampledCounts_humansOnly <- spectrum_1mer_forSigfit_NONdownsampledCounts %>%
  filter(label %in% humans) %>%
  remove_rownames %>% # have to remove any preexisting rownames 
  column_to_rownames(var="label") # adds row names
# reorder to match aging :
spectrum_1mer_forSigfit_NONdownsampledCounts_humansOnly <- spectrum_1mer_forSigfit_NONdownsampledCounts_humansOnly[,names(agingSignatures)]
##### 3mer:
spectrum_3mer_forSigfit_NONdownsampledCounts <- pivot_wider(spectrum_3mer_processed[,c("label","total_mutations_projected","mutation_label")],id_cols=c("label"),names_from = "mutation_label",values_from ="total_mutations_projected" )

spectrum_3mer_forSigfit_NONdownsampledCounts_humansOnly <- spectrum_3mer_forSigfit_NONdownsampledCounts %>%
  filter(label %in% humans) %>%
  remove_rownames %>% # have to remove any preexisting rownames 
  column_to_rownames(var="label") # adds row names

# reorder to match cosmic :
spectrum_3mer_forSigfit_NONdownsampledCounts_humansOnly <- spectrum_3mer_forSigfit_NONdownsampledCounts_humansOnly[,names(cosmic_signatures_v3.2)]
############## NO LONGER DOING SIGFIT VERSION OF OPPOTURNITIES: makes CpGs too dominant #############
# instead am scaling by human oppos above
# need unique because only care about 1 appearance of each target type: 
#spectrum_1mer_forSigfit_opportunities <- pivot_wider(unique(spectrum_1mer_forSigfit[,c("label","total_target_count","mutation_label")]),id_cols=c("label"),names_from = "mutation_label",values_from ="total_target_count" )

# okay actually I want this to be for every mutation type 
# so C>A and C>T would have same oppos; C>T_CpG would have CpG oppos; A>T and A>C woudl have same oppos etc. 
# okay the https://github.com/kgori/sigfit/issues/58 they say to make opportunities normalize to 1 -- that's intriguing... let's try it! 
#spectrum_1mer_forSigfit_opportunities <- spectrum_1mer_forSigfit_opportunities %>%
#  remove_rownames %>% # have to remove any preexisting rownames 
#  column_to_rownames(var="label") # adds row names




#if(sum(names(spectrum_1mer_forSigfit_opportunities)!=names#(spectrum_1mer_forSigfit_downsampledCounts))!=0){
#  print("opportunities and counts are not in same order! this is fatally bad!")
#  exit()
#}

#write.table(spectrum_1mer_forSigfit_opportunities,paste0(plotdir,"opportunities_1mer.forsigfit.noteThatTargetsAreRepeatedPerMutType.txt"),row.names=T,sep="\t",quote=F)

########### normalize opportunities -- this gets done internally with sigfit don't have to do it #############
#spectrum_1mer_forSigfit_opportunities_RowSums <- rowSums(spectrum_1mer_forSigfit_opportunities)
# normalize (note these are not GENOME FRACTIONS -- they are just normalized counts; will be the same for A>X and C>X because it's just the A, C counts. but have to add them up to get genome fractions -- weird. not sure about this.)
#spectrum_1mer_forSigfit_opportunities_normalized <- spectrum_1mer_forSigfit_opportunities %>%
#  mutate_all(~./spectrum_1mer_forSigfit_opportunities_RowSums)
# this is waht the sigfit dev wants but note that fractions DO NOT REPRESENT GENOME FRACTIONS. BECAUSE TARGETS ARE MULTI-COUNTED. I'M NOT SURE ABOUT THIS. maintains right relationships between them (e.g. A / C targets is still correct ratio ) ; but doesn't rep genome fraction. 
# does let you see that fin whales are weird. 



######## think about contextualizing these

# to try:
# figure out how to plot cosines per species and extract more info! 
# try sigfit with normalized oppos
# try triplets
# try other models 
############################# SIGFIT STUFF STARTS HERE ########################################


######### function to get cosine similarities in a dataframe instead of as plot titles ######
GetCosineSimilarity_allSamples <- function(mcmc_results,runlabel){
  counts=mcmc_results$data$counts_int
  samples=rownames(counts)
  reconstructions=retrieve_pars(mcmc_results,"reconstructions")$mean # get mean
  NSAMP=length(samples)
  allCos=data.frame()
  for(i in seq(1,NSAMP)){
    samplename=samples[i]
    cossimilarity=cosine_sim(counts[i,],reconstructions[i,]) # get cosine similarity (this is how sigfit does it in plotting code -- not very efficient, but oh well )
    cosdf = data.frame(species=samplename,cosine_similarity=cossimilarity)
    allCos=bind_rows(allCos,cosdf)
  }
  allCos$runlabel <- runlabel
  return(allCos)
}


makeHeatmapPlot <- function(sigfitResults,runlabel,outdir)
{
  reconstructions <- retrieve_pars(sigfitResults,par="reconstructions")$mean
  data <- sigfitResults$data$counts_real
  reconstructions_over_data <- reconstructions / data # divide tables (works bc are same dimensions)
  reconstructions_over_data_melt <- melt(reconstructions_over_data)
  colnames(reconstructions_over_data_melt) <- c("pop","mutationLabel","ratio_recon_over_data")
  # detect if it's 1mer or 3mer:
  if(nchar(as.character(reconstructions_over_data_melt[1,]$mutationLabel))==3){ # 1mer 
    reconstructions_over_data_melt$ancestralbp <- substr(reconstructions_over_data_melt$mutationLabel,1,1)
    reconstructions_over_data_melt$derivedbp <- substr(reconstructions_over_data_melt$mutationLabel,3,3)
    # deal with CpG:
    reconstructions_over_data_melt[reconstructions_over_data_melt$mutationLabel=="C.T_CpG",]$ancestralbp <- "CpG"
    
    heatmap1 <- ggplot(reconstructions_over_data_melt,aes(x=pop,y=mutationLabel,fill=ratio_recon_over_data))+
      geom_tile()+
      scale_fill_gradient2(low = "blue",high="red",mid="white",midpoint = 1,limits=c(0.5,2.5))+
      ggtitle(paste0("ratio of reconstruction/data\n",runlabel))+
      geom_text(aes(label=round(ratio_recon_over_data,2)))
    heatmap1
    ggsave(paste0(outdir,"heatmapOfReconstructions_vs_data.",runlabel,".png"),width=dim(data)[1]*1.5,height=5)
    return(heatmap1)
    
  } else if(nchar(as.character(reconstructions_over_data_melt[1,]$mutationLabel))==7) {
    reconstructions_over_data_melt$centralBP <- paste0(substr(reconstructions_over_data_melt$mutationLabel,2,2),".",substr(reconstructions_over_data_melt$mutationLabel,6,6))
    reconstructions_over_data_melt$ancestralbp <- substr(reconstructions_over_data_melt$centralBP,1,1)
    reconstructions_over_data_melt$derivedbp <- substr(reconstructions_over_data_melt$centralBP,3,3)    
    
    ####### NEED TO CODE UP HEATMAP FOR 3mers ##### 
    
  } else {
    print('something has gone wrong -- are your data not 1mers or 3mers?')
  }

  
}

makeSigfitResidualsPlot <- function(sigfitResults,runlabel,outdir){
  reconstructions <- retrieve_pars(sigfitResults,par="reconstructions")$mean
  data <- sigfitResults$data$counts_real
  residuals <-  data - reconstructions # subtract recon from data (works bc are same dimensions)
  residuals_melt <- melt(residuals)
  colnames(residuals_melt) <- c("pop","mutationLabel","residual_data_minus_model")
  # want to plot residual vs fitted value https://statisticsbyjim.com/regression/check-residual-plots-regression-analysis/
  reconstructions_melt <- melt(reconstructions)
  colnames(reconstructions_melt) <- c("pop","mutationLabel","reconstructionValue")
  # merge:
  combo_residuals_reconstructions <- merge(residuals_melt,reconstructions_melt,by=c("pop","mutationLabel"))
  

  residualsPlot <- ggplot(combo_residuals_reconstructions,aes(x=mutationLabel,y=residual_data_minus_model))+
    geom_violin()+
    geom_point()+
    geom_hline(yintercept = 0)+
    ggtitle(paste0("residuals (data - reconstruction)\n",runlabel))
  residualsPlot
  ggsave(paste0(outdir,"residualsDataMinusReconstruction.",runlabel,".png"),residualsPlot,width=dim(data)[1]*.75,height=5)
  return(residualsPlot)
}
  
# if you don't want to plot, you just awnt the residuals (for example if comibing across runs in one plot)
calculateSigfitResiduals <- function(sigfitResults){
  reconstructions <- retrieve_pars(sigfitResults,par="reconstructions")$mean
  data <- sigfitResults$data$counts_real
  residuals <-  data - reconstructions # subtract recon from data (works bc are same dimensions)
  residuals_melt <- melt(residuals)
  colnames(residuals_melt) <- c("pop","mutationLabel","residual_data_minus_model")
  reconstructions_melt <- melt(reconstructions)
  colnames(reconstructions_melt) <- c("pop","mutationLabel","reconstructionValue")
  # merge:
  # also put data in too:
  data_melt <- melt(data)
  colnames(data_melt) <- c("pop","mutationLabel","counts_real")
  combo_residuals_reconstructions <- merge(residuals_melt,reconstructions_melt,by=c("pop","mutationLabel"))
  combo_residuals_reconstructions_data <- merge(combo_residuals_reconstructions,data_melt,by=c("pop","mutationLabel"))
  

  return(combo_residuals_reconstructions_data)
}

########## for extracting de novo signatures ###############
run_plot_sigfit_EXTRACTING <- function(counts,opportunities=NULL,nsig_min,nsig_max,niter,model,wd,note){
  
  runlabel=paste("sigfit.extract_signatures",note,model,"iter",niter,sep=".")
  outdir=paste0(wd,runlabel,"/")
  
  dir.create(outdir,showWarnings = F,recursive=T)
  # first need to check that counts and oppos have same names in order (same species and same signatures)
  if(!is.null(opportunities)){
  if(sum(rownames(counts)!=rownames(opportunities))!=0){
    print("different rownames for counts and opportunities")
    break
  }
  if(sum(names(counts)!=names(opportunities))!=0){
    print("different column names/order for counts and opportunities")
    break
  }
  }
  
  # first fit with a small number of iterations (2000) ; then redo 'best' with more! 
  sigfit_results <- extract_signatures(counts = counts,opportunities = opportunities,nsignatures=nsig_min:nsig_max,iter=2000,model=model) # noet this iter is NOT niter; is permanently 2000; note adding adapt_delta 0.99 just for final best fit to make sure no divergent results. not adding it here as its too slow
  saveRDS(sigfit_results,paste0(outdir,"sigfit.extract_signatures.",note,".",model,".allSignatureCounts.lowiter.1000.rds"))
  
  # what is going on with R?
  best=plot_gof(sigfit_results) # see what fits best -- still trying to save this plot 
  # need to save this plot somehow.
  #plot_all(sigfit_results[[best]],out_path = outdir,prefix=runlabel)
  
  # rerun with higher iterations (niter) and the best set:
  sigfit_results_best <- extract_signatures(counts = counts,opportunities = opportunities,nsignatures=best,iter=niter,model=model,control=list(adapt_delta=0.99)) # adding high adapt delta (smaller step sizes) to avoid divergent results (see stan docs)
  saveRDS(sigfit_results_best,paste0(outdir,runlabel,".BEST.rds"))
  plot_all(sigfit_results_best,out_path = outdir,prefix=paste0(runlabel,".BEST"))
  # pull out and plot cosine similarities
  cosine_sims = GetCosineSimilarity_allSamples(sigfit_results_best,runlabel)
  write.table(cosine_sims,file=paste0(outdir,runlabel,".BEST.cosineSimilarities.txt"),row.names = F,quote=F,sep="\t")
  
  # save plot
  png(file=paste0(outdir,runlabel,".gof.png"))
  print(plot_gof(sigfit_results)) # see what fits best -- still trying to save this plot 
  dev.off()
  
  makeHeatmapPlot(sigfit_results_best,paste0(runlabel,".BEST"),outdir) # makes and saves a heatmap plot of reconstructions vs data.
  makeSigfitResidualsPlot(sigfit_results_best,paste0(runlabel,".BEST"),outdir) # makes residuals plot
  
  
  return(cosine_sims) # maybe return this so you can save it 
  # pull out and plot exposures (? maybe later)
  
}


####### for fitting a prior signatures ###############
# make a heatmap plot of reconstructions:


run_plot_sigfit_FITTING <- function(counts,opportunities=NULL,signaturesToFit,signatureslabel,niter,model,wd,note){
  runlabel=paste("sigfit.fit",signatureslabel,"signatures",note,model,"iter",niter,sep=".")
  
  outdir=paste0(wd,runlabel,"/")
  dir.create(outdir,showWarnings = F,recursive=T)
  # check if counts and sigs are in same order;:
  if(sum(names(counts)!=names(signaturesToFit))!=0){
    print("counts and signatures are not in same order!")
    break
  }
  if(!is.null(opportunities)){
    if(sum(rownames(counts)!=rownames(opportunities))!=0){
      print("different rownames for counts and opportunities")
      break
    }
    if(sum(names(counts)!=names(opportunities))!=0){
      print("different column names/order for counts and opportunities")
      break
    }
  }
  
  sigfit_results <- fit_signatures(counts,signatures = signaturesToFit,opportunities = opportunities,iter=niter,model=model,control=list(adapt_delta=0.99))
  saveRDS(sigfit_results,paste0(outdir,runlabel,".rds"))
  plot_all(sigfit_results,out_path = outdir,prefix=runlabel)
  # pull out and plot cosine similarities
  cosine_sims = GetCosineSimilarity_allSamples(sigfit_results,runlabel)
  write.table(cosine_sims,file=paste0(outdir,runlabel,".cosineSimilarities.txt"),row.names = F,quote=F,sep="\t")
  # make a heatmap plot:
  makeHeatmapPlot(sigfit_results,runlabel,outdir) # makes and saves a heatmap plot of reconstructions vs data.
  makeSigfitResidualsPlot(sigfit_results,runlabel,outdir) # makes and saves a residuals plot
  return(cosine_sims) # maybe return this so you can save it 
  # pull out and plot exposures (? maybe later)
  
}

# if you just want cosine sims and nothing written out:
run_plot_sigfit_FITTING_NOPLOTTING <- function(counts,opportunities=NULL,signaturesToFit,signatureslabel,niter,model,wd,note){
  runlabel=paste("sigfit.fit",signatureslabel,"signatures",note,model,"iter",niter,sep=".")
  
  #outdir=paste0(wd,runlabel,"/")
  #dir.create(outdir,showWarnings = F,recursive=T)
  # check if counts and sigs are in same order;:
  if(sum(names(counts)!=names(signaturesToFit))!=0){
    print("counts and signatures are not in same order!")
    break
  }
  if(!is.null(opportunities)){
    if(sum(rownames(counts)!=rownames(opportunities))!=0){
      print("different rownames for counts and opportunities")
      break
    }
    if(sum(names(counts)!=names(opportunities))!=0){
      print("different column names/order for counts and opportunities")
      break
    }
  }
  
  sigfit_results <- fit_signatures(counts,signatures = signaturesToFit,opportunities = opportunities,iter=niter,model=model,control=list(adapt_delta=0.99))
  #saveRDS(sigfit_results,paste0(outdir,runlabel,".rds"))
  #plot_all(sigfit_results,out_path = outdir,prefix=runlabel)
  # pull out and plot cosine similarities
  cosine_sims = GetCosineSimilarity_allSamples(sigfit_results,runlabel)
  #write.table(cosine_sims,file=paste0(outdir,runlabel,".cosineSimilarities.txt"),row.names = F,quote=F,sep="\t")
  makeHeatmapPlot(sigfit_results,runlabel,outdir) # makes and saves a heatmap plot of reconstructions vs data.
  makeSigfitResidualsPlot(sigfit_results,runlabel,outdir) # makes and saves a residuals plot
  
  return(cosine_sims) # maybe return this so you can save it 
  # pull out and plot exposures (? maybe later)
  
}
##### for fitting a priori signatures then extracting additional ones: 
# this will extract the 'best' and rerun it with more iterations
run_plot_sigfit_FIT_EXT  <- function(counts,opportunities=NULL,signaturesToFit,signatureslabel,ndenovo, niter,model,wd,note) {
  if(sum(names(counts)!=names(signaturesToFit))!=0){
    print("counts and signatures are not in same order!")
    break
  }
  
  runlabel=paste("sigfit.fit",signatureslabel,"signatures",note,model,"iter",niter,sep=".")
  
  outdir=paste0(wd,runlabel,"/")
  dir.create(outdir,showWarnings =F,recursive=T)
  
  sigfit_results <-  fit_extract_signatures(counts,signatures = signaturesToFit,num_extra_sigs = ndenovo,opportunities = opportunities,iter=niter,model=model,control=list(adapt_delta=0.99))
  # adding adapt delta to improve robustness; slower
  saveRDS(sigfit_results,paste0(outdir,runlabel,".rds"))
  plot_all(sigfit_results,out_path = outdir,prefix=runlabel)
  # pull out and plot cosine similarities
  cosine_sims = GetCosineSimilarity_allSamples(sigfit_results,runlabel)
  write.table(cosine_sims,file=paste0(outdir,runlabel,".cosineSimilarities.txt"),row.names = F,quote=F,sep="\t")
  makeHeatmapPlot(sigfit_results,runlabel,outdir) # makes and saves a heatmap plot of reconstructions vs data.
  
  return(cosine_sims) # maybe return this so you can save it 
  # pull out and plot exposures (? maybe later)
  
}



#################### RUN SIGFIT : 1mers ##########
# eventually loop this better 
allcosineSims = data.frame()
###### extract sigs #########
cosinesims1 = run_plot_sigfit_EXTRACTING(counts=spectrum_1mer_forSigfit_downsampledCounts,opportunities = NULL,nsig_min = 2,nsig_max = 8,model = model,niter = niter,wd=plotdir,note = "extract_sigs") 


allcosineSims = bind_rows(allcosineSims,cosinesims1)
# read back in the signatures:
# Regress exposure of puberty signature against gen time # 
# get cosine similarity between de novo signatures and aging signatures 
deNovoSigs_best <- readRDS(paste0(plotdir,paste0("sigfit.extract_signatures.extract_sigs.",model,".iter.",niter,"/sigfit.extract_signatures.extract_sigs.",model,".iter.",niter,".BEST.rds")))

makeHeatmapPlot(deNovoSigs_best,"extract_signatures_best",outdir = plotdir)
makeSigfitResidualsPlot(deNovoSigs_best,"extract_signatures_best",outdir = plotdir)

deNovoSigs_best_signatures <- retrieve_pars(deNovoSigs_best,par="signatures")$mean
deNovoSigs_best_exposures <- retrieve_pars(deNovoSigs_best,par="exposures")$mean
###### fit aging sigs: ##########
cosinesims2 <- run_plot_sigfit_FITTING(counts=spectrum_1mer_forSigfit_downsampledCounts,signaturesToFit = agingSignatures,signatureslabel = "agingSignatures",opportunities = NULL,wd=plotdir,niter=niter,model=model,note="fittingAgingSignatures")
allcosineSims = bind_rows(allcosineSims,cosinesims2)

cosinesims2_results <- readRDS(paste0(plotdir,"sigfit.fit.agingSignatures.signatures.fittingAgingSignatures.",model,".iter.",niter,"/sigfit.fit.agingSignatures.signatures.fittingAgingSignatures.",model,".iter.",niter,".rds"))
makeHeatmapPlot(cosinesims2_results,"agingSignatures",plotdir)
makeSigfitResidualsPlot(cosinesims2_results,"agingSignatures",plotdir)
####### idea: plot residuals and choose model that has non-correlated residuals (stripes)

########## fit aging + CpG (not needed-- can skip) ########## 
#cosinesims3 <- run_plot_sigfit_FITTING(counts=spectrum_1mer_forSigfit_downsampledCounts,signaturesToFit = agingPlusCpGSignatures,signatureslabel = "agingPlusCpGSignatures",opportunities = NULL,wd=plotdir,niter=niter,model=model,note="fittingAgingSignatures")
#allcosineSims = bind_rows(allcosineSims,cosinesims3)

######## FIT+Extract: aging +1 denovo #######
cosinesims4 <- run_plot_sigfit_FIT_EXT(counts=spectrum_1mer_forSigfit_downsampledCounts,opportunities =NULL,signaturesToFit = agingSignatures,signatureslabel = "agingSignatures",ndenovo = 1,niter = niter,model = model,note = "agingPlus1DeNovo" ,wd=plotdir)
allcosineSims = bind_rows(allcosineSims,cosinesims4)

######## FIT+Extract: aging +2 denovo #######
#cosinesims5 <- run_plot_sigfit_FIT_EXT(counts=spectrum_1mer_forSigfit_downsampledCounts,opportunities =NULL,signaturesToFit = agingSignatures,signatureslabel = "agingSignatures",ndenovo = 2,niter = niter,model = model,note = "agingPlus2DeNovo",wd=plotdir )
#allcosineSims = bind_rows(allcosineSims,cosinesims5)

####### SKIP: FIT+ extract: aging+Cpg+1 de novo #############
cosinesims6 <- run_plot_sigfit_FIT_EXT(counts=spectrum_1mer_forSigfit_downsampledCounts,opportunities =NULL,signaturesToFit = agingPlusCpGSignatures,signatureslabel = "agingPlusCpGSignatures",ndenovo = 1,niter = niter,model = model,note = "agingPlusCpGPlus1DeNovo",wd=plotdir )
allcosineSims = bind_rows(allcosineSims,cosinesims6)


########### Extract: HUMANS ONLY : downsampled ##############
cosinesims7 <- run_plot_sigfit_EXTRACTING(counts=spectrum_1mer_forSigfit_downsampledCounts_humansOnly,opportunities = NULL,nsig_min = 2,nsig_max = 8,model = model,niter = niter,wd=plotdir,note = "extract_sigs_humansOnly") 

#deNovoSigs_best_HUMANONLY <-  readRDS(paste0(plotdir,paste0("sigfit.extract_signatures.extract_sigs.",model,".iter.",niter,"/sigfit.extract_signatures.extract_sigs_humansOnly.",model,".iter.",niter,".BEST.rds")))


#deNovoSigs_best_HUMANONLY_signatures <- retrieve_pars(deNovoSigs_best_HUMANONLY,par="signatures")$mean
#deNovoSigs_best_HUMANONLY_exposures <- retrieve_pars(deNovoSigs_best_HUMANONLY,par="exposures")$mean

####### Extract: humans only : non-downsampled ######

########### FIT AGING: HUMANS ONLY : downsampled ##############
cosinesims8 <- run_plot_sigfit_FITTING(counts=spectrum_1mer_forSigfit_downsampledCounts_humansOnly,signaturesToFit = agingSignatures,signatureslabel = "agingSignatures",opportunities = NULL,wd=plotdir,niter=niter,model=model,note="fittingAgingSignatures_humansonly")
# maybe dont use downsampled humans? 
# maybe dont use downsampled humans? 
cosinesims8_results <- readRDS(paste0(plotdir,"/sigfit.fit.agingSignatures.signatures.fittingAgingSignatures_humansonly.",model,".iter.",niter,"/sigfit.fit.agingSignatures.signatures.fittingAgingSignatures_humansonly.",model,".iter.",niter,".rds"))
cosinesims8_exposures <- retrieve_pars(cosinesims8_results,par="exposures")$mean
cosinesims8_exposures$pop  <- rownames(cosinesims8_exposures)

cosinesims8_exposures_plot <- ggplot(melt(cosinesims8_exposures),aes(x=pop,y=value))+
  facet_wrap(~variable,scales="free_y")+
  geom_point()
cosinesims8_exposures_plot
ggsave(paste0(plotdir,"exposuresToAgingSigs.HumansOnly.closeup.png"),cosinesims8_exposures_plot,height=4,width=14)

makeHeatmapPlot(cosinesims8_results,"test",plotdir)


########### FIT AGING: HUMANS ONLY : not downsampled ##############
cosinesims9 <- run_plot_sigfit_FITTING(counts=spectrum_1mer_forSigfit_NONdownsampledCounts_humansOnly,signaturesToFit = agingSignatures,signatureslabel = "agingSignatures",opportunities = NULL,wd=plotdir,niter=niter,model=model,note="fittingAgingSignatures_humansonly_notdownsampled")



########## plot all cossime similarities ########
# not including the humans only cosinesims
write.table(allcosineSims,paste0(plotdir,"AllcosineSims.txt"),row.names = F,quote=F,sep="\t")



cosineSimPlot <- ggplot(allcosineSims,aes(y=species,x=cosine_similarity,color=runlabel))+
  geom_point()
cosineSimPlot
ggsave(paste0(plotdir,"cosineSimilarityPlot.comparingRuns.png"),cosineSimPlot,height=7,width=18)

cosineSimPlot2 <- ggplot(allcosineSims,aes(x=cosine_similarity,y=reorder(runlabel,desc(runlabel))))+
  geom_boxplot()+
  ylab("")
cosineSimPlot2 
ggsave(paste0(plotdir,"cosineSimilarityPlot.comparingRuns2.png"),cosineSimPlot2,height=7,width=15)
  
# pick only most relevant: (don't need CpG signature stuff)
runsToPlot = c(paste0("sigfit.extract_signatures.extract_sigs.",model,".iter.",niter),paste0("sigfit.fit.agingSignatures.signatures.fittingAgingSignatures.",model,".iter.",niter),paste0("sigfit.fit.agingSignatures.signatures.agingPlus1DeNovo.",model,".iter.",niter))

cosineSimPlot3 <- ggplot(allcosineSims[allcosineSims$runlabel %in% runsToPlot,],aes(y=species,x=cosine_similarity,color=runlabel))+
  geom_point()
cosineSimPlot3
ggsave(paste0(plotdir,"cosineSimilarityPlot.subsetOfRuns.png"),cosineSimPlot3,height=7,width=18)

cosineSimPlot4 <- ggplot(allcosineSims[allcosineSims$runlabel %in% runsToPlot,],aes(x=cosine_similarity,y=reorder(runlabel,desc(runlabel))))+
  geom_boxplot()+
  theme_bw()+
  ylab("")
cosineSimPlot4
ggsave(paste0(plotdir,"cosineSimilarityPlot.subsetOfRuns.boxplot.png"),cosineSimPlot4,height=7,width=12)

# TO DO: make another plot here is that is just aging, de novo / agin , de novo and sbs1+5

####### show most improved species from aging --> aging +1 extra df#########
# question: but shouldn't it always improve when you add an extra
allcosineSims_CompareAgingPlusOne <- allcosineSims[allcosineSims$runlabel %in% runsToPlot[c(2,3)],]

allcosineSims_CompareAgingPlusOne$runLabel_nicer <- ""
allcosineSims_CompareAgingPlusOne[allcosineSims_CompareAgingPlusOne$runlabel == paste0("sigfit.fit.agingSignatures.signatures.fittingAgingSignatures.",model,".iter.",niter),]$runLabel_nicer <- "aging"
allcosineSims_CompareAgingPlusOne[allcosineSims_CompareAgingPlusOne$runlabel == paste0("sigfit.fit.agingSignatures.signatures.agingPlus1DeNovo.",model,".iter.",niter),]$runLabel_nicer <- "aging_plus_de_novo"

allcosineSims_CompareAgingPlusOne_wide <- pivot_wider(allcosineSims_CompareAgingPlusOne,names_from = runLabel_nicer,values_from = cossine_similarity,id_cols = species)

head(allcosineSims_CompareAgingPlusOne_wide)
allcosineSims_CompareAgingPlusOne_wide$delta_improvement_adding_denovo <- allcosineSims_CompareAgingPlusOne_wide$aging_plus_de_novo - allcosineSims_CompareAgingPlusOne_wide$aging

improvementAdding1DenovoToAging <- ggplot(allcosineSims_CompareAgingPlusOne_wide,aes(y=species,x=delta_improvement_adding_denovo))+
  geom_point()+
  theme_bw()+
  ggtitle("improvement in cosine similarity when adding a de novo signature")+
  xlab("improvement in cosine similarity (aging + de novo - aging alone)")
ggsave(paste0(plotdir,"improvement_adding_1denovo_toAging.png"),improvementAdding1DenovoToAging)
####### just plot fit of aging signatures per species ###########
cosineSimPlot4 <- ggplot(allcosineSims[allcosineSims$runlabel==paste0("sigfit.fit.agingSignatures.signatures.fittingAgingSignatures.",model,".iter.",niter),],aes(y=species,x=cosine_similarity,color=runlabel))+
  geom_point(size=3)+
  theme(text=element_text(size=14))
cosineSimPlot4
ggsave(paste0(plotdir,"cosineSimilarityPlot.AgingOnly.png"),cosineSimPlot4,height=7,width=14)

############# run sigfit on humans by themselves ############
# without downsampling (!) <-----
# with downsampling 
################ compare de novo signatures to aging signatures ##############


GetCosineSimilarity_betweenSignatures <- function(de_novo_signature_set,a_priori_signature_set){
  # make sure names are in same order
  if(sum(colnames(de_novo_signature_set)!=colnames(a_priori_signature_set))!=0){
    print("different column names for signature sets")
    break
  }
  nsig_denovo=dim(de_novo_signature_set)[1] # number of sigs
  nsig_apriori=dim(a_priori_signature_set)[1] # number of sigs
  allCos=data.frame()
  for(i in seq(1,nsig_denovo)){
    for(j in seq(1,nsig_apriori)){
    cossimilarity=cosine_sim(as.numeric(de_novo_signature_set[i,]),as.numeric(a_priori_signature_set[j,])) # get cosine similarity (this is how sigfit does it in plotting code -- not very efficient, but oh well )
    cosdf = data.frame(de_novo_signature=rownames(de_novo_signature_set)[i],a_priori_signature=rownames(a_priori_signature_set)[j],cosine_similarity=cossimilarity)
    allCos=bind_rows(allCos,cosdf)
  }}
  return(allCos)
}

cosineSimBetweenSignatures <- GetCosineSimilarity_betweenSignatures(deNovoSigs_best_signatures,agingSignatures)
write.table(cosineSimBetweenSignatures,paste0(plotdir,"cosineSimilarityBetweenDeNovoAndAPrioriSignatures.txt"),row.names = F,quote=F,sep="\t")

cosine_similarity_among_signaturesPlot1 <- ggplot(cosineSimBetweenSignatures,aes(y=de_novo_signature,fill=cosine_similarity,x=a_priori_signature))+
  geom_tile()+
  scale_fill_viridis_c()+
  geom_text(aes(label=round(cosine_similarity,2)))
cosine_similarity_among_signaturesPlot1
ggsave(paste0(plotdir,"cosine_similarity_among_signaturesPlot1.png"),cosine_similarity_among_signaturesPlot1)

# this does the same thing as the in house sigfit function match_signatures, but is more useful output
match_signatures(agingSignatures,deNovoSigs_best_signatures)
# Optimal assignment:
# multinomial: 1 => 3, 2 => 2, 3 => 4
# poisson: 1 => 2, 2 => 1, 3 => 4

# get clr-aitchison distance between signatures :: try this!
GetAitchisonDistance_betweenSignatures <- function(signature_set_1,signature_set_2){
  # make sure names are in same order
  if(sum(colnames(signature_set_1)!=colnames(signature_set_2))!=0){
    print("different column names for signature sets")
    break
  }
  nsig1=dim(signature_set_1)[1] # number of sigs
  nsig2=dim(signature_set_2)[1] # number of sigs
  allDist=data.frame()
  for(i in seq(1,nsig1)){
    for(j in seq(1,nsig2)){
      # need to make a matrix: 
      dist=tidy(dist(rbind(clr(as.numeric(signature_set_1[i,])),clr(as.numeric(signature_set_2[j,]))))) # get cosine similarity (this is how sigfit does it in plotting code -- not very efficient, but oh well )
      distdf = data.frame(signature_set_1=rownames(signature_set_1)[i],signature_set_2=rownames(signature_set_2)[j],aitchison_clr_dist=dist$distance)
      allDist=bind_rows(allDist,distdf)
    }}
  return(allDist)
  
  
}



########### how similar are aging signatures to each other? ##########
aging_aging <- GetCosineSimilarity_betweenSignatures(agingSignatures,agingSignatures)

aging_agingPlot <- ggplot(aging_aging,aes(x=de_novo_signature,y=a_priori_signature,fill=cosine_similarity))+
  geom_tile()+
  scale_fill_viridis_c()+
  ylab("")+
  xlab("")+
  geom_text(aes(label=round(cosine_similarity,2)))
aging_agingPlot
ggsave(paste0(plotdir,"cosine_similarity_among_Aging_ToThemselves.png"),aging_agingPlot)

aging_aging_Aitchison <- GetAitchisonDistance_betweenSignatures(agingSignatures,agingSignatures)
ggplot(aging_aging_Aitchison,aes(x=signature_set_1,y=signature_set_2,fill=aitchison_clr_dist))+
  scale_fill_viridis_c(direction = -1)+
  geom_tile()+
  geom_text(aes(label=signif(aitchison_clr_dist,2)),color="lightgray")
######### compare aging plus one de novo to de novo signatures: #############
# need to make this read-in more generic! if I change niter etc this will break.
# now best to do ? 
aging_plusOneDenovo <- readRDS(paste0(plotdir,paste0("sigfit.fit.agingSignatures.signatures.agingPlus1DeNovo.",model,".iter.",niter,"/sigfit.fit.agingSignatures.signatures.agingPlus1DeNovo.",model,".iter.",niter,".rds")))
aging_plusOneDenovo_signatures <- retrieve_pars(aging_plusOneDenovo,par="signatures")$mean


# rename so that it has aging names:
rownames(aging_plusOneDenovo_signatures) <- c(rownames(agingSignatures),"de_novo")

compareDeNovo_toAgingPlusOneDenovo <- GetCosineSimilarity_betweenSignatures(de_novo_signature_set = deNovoSigs_best_signatures,a_priori_signature_set = aging_plusOneDenovo_signatures)

compareDeNovo_toAgingPlusOneDenovo <- ggplot(compareDeNovo_toAgingPlusOneDenovo,aes(x=a_priori_signature,y=de_novo_signature,fill=cosine_similarity))+
  geom_tile()+
  scale_fill_viridis_c()+
  geom_text(aes(label=round(cosine_similarity,2)))
compareDeNovo_toAgingPlusOneDenovo
ggsave(paste0(plotdir,"comparingAgingPlusOneDeNovo.ToDenovoSignatures.png"),compareDeNovo_toAgingPlusOneDenovo)
# so matches aging sig 1 (maternal) to 3 of de novo
# matches aging sig 2 (paternal) to 2 of de novo 
# matches aging sig 3 (puberty) to 4 of de novo
# so 1 of de novo doesn't get matched (?)

# generation times <-
reproductiveAges=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/mammal_generation_times/subset.small.txt",header=T)
# this is just a subset of species
########## get exposures to aging signatures #########
fitToAging <- readRDS(paste0(plotdir,"/sigfit.fit.agingSignatures.signatures.fittingAgingSignatures.",model,".iter.",niter,"/sigfit.fit.agingSignatures.signatures.fittingAgingSignatures.",model,".iter.",niter,".rds"))
fitToAging_exposures <- retrieve_pars(fitToAging,par="exposures")$mean
# now regress somehow
fitToAging_exposures
fitToAging_exposures$label <- rownames(fitToAging_exposures)
fitToAging_exposures_mergedWithreproages <- merge(fitToAging_exposures,reproductiveAges,by="label")

# plotting exposures vs puberty etc
fitToAging_exposures_mergedWithreproages_melt <- melt(select(fitToAging_exposures_mergedWithreproages,-source),id.vars=c("label","Rspan_d","AFR_d"),value.name = "exposure")
head(fitToAging_exposures_mergedWithreproages_melt)

############### reproductive span ##########################
rspan_exposureToAgingSignatures <- ggplot(fitToAging_exposures_mergedWithreproages_melt,aes(x=Rspan_d,y=exposure,color=variable))+
  geom_point()+
  facet_wrap(~variable,scales="free")+
  geom_smooth(method="lm",se = F)+
  geom_text(aes(label=label),size=2)+
  ggtitle("Reproductive Span\nnote that without phylogenetically independent contrasts, lm significance isn't reliable")+
  theme_bw()+
  xlab("Reproductive span (days)")
rspan_exposureToAgingSignatures
ggsave(paste0(plotdir,"rspan.vs.agingSignatureExposures.png"),rspan_exposureToAgingSignatures,height=7,width=12)
       
############### age at first reproduction #############
AFR_exposureToAgingSignatures <- ggplot(fitToAging_exposures_mergedWithreproages_melt,aes(x=AFR_d,y=exposure,color=variable))+
  geom_point()+
  facet_wrap(~variable,scales="free")+
  geom_smooth(method="lm",se = F)+
  geom_text(aes(label=label),size=2)+
  ggtitle("Reproductive Span\nnote that without phylogenetically independent contrasts, lm significance isn't reliable")+
  theme_bw()+
  xlab("Age at first reproduction (days)")
AFR_exposureToAgingSignatures
ggsave(paste0(plotdir,"afr.vs.agingSignatureExposures.png"),AFR_exposureToAgingSignatures,height=7,width=12)

################################ run sigfit: 3mers  #############################################
outdir_3mer=paste0(plotdir,"3mers/")
dir.create(outdir_3mer,showWarnings = F)

allcosineSims_3mers = data.frame()
###### extract denovo triplet sigs #########
triplet_sigfit1 = run_plot_sigfit_EXTRACTING(counts=spectrum_3mer_forSigfit_downsampledCounts,opportunities = NULL,nsig_min = 2,nsig_max = 8,model = model,niter = niter,wd=outdir_3mer,note = "extract_de_novo_sigs") 

allcosineSims_3mers = bind_rows(allcosineSims,triplet_sigfit1)
extractTripletsResult <- readRDS(paste0(outdir_3mer,"sigfit.extract_signatures.extract_de_novo_sigs.",model,".iter.",niter,"/sigfit.extract_signatures.extract_de_novo_sigs.",model,".iter.",niter,".BEST.rds"))

triplet_sigfit1_sigs <- retrieve_pars(extractTripletsResult,par="signatures")$mean
triplet_sigfit1_exposures <-  retrieve_pars(extractTripletsResult,par="exposures")$mean
triplet_sigfit1_exposures$label <- rownames(triplet_sigfit1_exposures)

# match to cosmic signatures
compareDenovoTOCosmic_3mers <- GetCosineSimilarity_betweenSignatures(triplet_sigfit1_sigs,cosmic_signatures_v3.2)
write.table(compareDenovoTOCosmic_3mers,paste0(outdir_3mer,"cosineSimilarityBetweenDeNovo.cosmic.txt"),row.names=F,quote=F,sep="\t")

# reorder sbs levels to be in order
compareDenovoTOCosmic_3mers$a_priori_signature <- factor(compareDenovoTOCosmic_3mers$a_priori_signature,levels=mixedsort(unique(compareDenovoTOCosmic_3mers$a_priori_signature))) 

similarityToCosmicSignaturesPlot <- ggplot(compareDenovoTOCosmic_3mers,aes(y=a_priori_signature,x=de_novo_signature,fill=cosine_similarity))+
  geom_tile()+
  scale_fill_viridis_c()+
  geom_text(aes(label=round(cosine_similarity,2)))+
  ggtitle("comparing de novo triplet signatures to SBS signatures")
similarityToCosmicSignaturesPlot
ggsave(paste0(outdir_3mer,"similarityOFDenovo.ToCosmicSignaturesPlot.png"),height=16,width=8)
# need to sort cosmic
match_signatures(triplet_sigfit1_sigs,cosmic_signatures_v3.2)
rownames(cosmic_signatures_v3.2)[76]
# de novo signature matches: 1 => 5, 2 => 76, 3 => 37, 4 => 1 ; note these numbers aren't SBS numbers, but the index of the SBS names. since there are 7a 7b 7c etc. 76 actually is sbs92!

######### extract de novo signatures: humans only! downsampled #######
spectrum_1mer_forSigfit_downsampledCounts_humansOnly

triplet_sigfit1_humansOnly = run_plot_sigfit_EXTRACTING(counts=spectrum_3mer_forSigfit_downsampledCounts_humansOnly,opportunities = NULL,nsig_min = 2,nsig_max = 8,model = model,niter = niter,wd=outdir_3mer,note = "extract_de_novo_sigs_humansOnly") 


############ make signatures from Hamid biorxiv paper #########
hamid_SIG1 <- 0.13*cosmic_signatures_v3.2["SBS1",]+0.87*cosmic_signatures_v3.2["SBS5",]
hamid_SIG2 <- 0.36*cosmic_signatures_v3.2["SBS7a",]+0.64*cosmic_signatures_v3.2["SBS11",]
# 3 and 4 are likely artifacts and aren't well reconstructed by cosmics (cosine sim 0.83 and 0.75 respectively)
hamid_SIG3 <- 0.26*cosmic_signatures_v3.2["SBS43",]+0.58*cosmic_signatures_v3.2["SBS47",]+0.16*cosmic_signatures_v3.2["SBS51",]
hamid_SIG4 <- 0.29*cosmic_signatures_v3.2["SBS16",]+0.23*cosmic_signatures_v3.2["SBS37",]+0.48*cosmic_signatures_v3.2["SBS39",]

hamid_df <- bind_rows(hamid_SIG1,hamid_SIG2,hamid_SIG3,hamid_SIG4)
rownames(hamid_df) <- c("hamid_SIG1","hamid_SIG2","hamid_SIG3","hamid_SIG4")

triplet_sigfit_hamid_humansOnly = run_plot_sigfit_FITTING(counts=spectrum_3mer_forSigfit_downsampledCounts_humansOnly,signaturesToFit =hamid_df, opportunities = NULL,model = model,niter = niter,wd=outdir_3mer,signatureslabel="hamid4sigs",note = "fit_hamid_4sigs") 
# note that these are imperfect reconstructions using cosmic signatures (cosine sim 0.9-1 to the hamid denovo signatures)
# what to do with these? fit them?
triplet_sigfit_hamid_humansOnly_results <- readRDS(paste0(outdir_3mer,"sigfit.fit.hamid4sigs.signatures.fit_hamid_4sigs.",model,".iter.",niter,"/sigfit.fit.hamid4sigs.signatures.fit_hamid_4sigs.",model,".iter.",niter,".rds"))
triplet_sigfit_hamid_humansOnly_Exposures <- retrieve_pars(triplet_sigfit_hamid_humansOnly_results,par="exposures")$mean
triplet_sigfit_hamid_humansOnly_Exposures$population <- rownames(triplet_sigfit_hamid_humansOnly_Exposures)

triplet_sigfit_hamid_humansOnly_Exposures_plot <- ggplot(melt(triplet_sigfit_hamid_humansOnly_Exposures),aes(x=population,y=value))+
  geom_point()+
  facet_wrap(~variable,scales="free_y")
triplet_sigfit_hamid_humansOnly_Exposures_plot
ggsave(paste0(plotdir,"exposuresToHamidSignatures.Closeup.png"),triplet_sigfit_hamid_humansOnly_Exposures_plot)
###################### project de novo 3mer signtaures down to 1mers and compare to de novo 1mer signatures *IMPORANt RESULT* ##################
triplet_sigfit1_sigs$signature <- rownames(triplet_sigfit1_sigs)
triplet_sigfit1_sigs_melt <- melt(triplet_sigfit1_sigs)
head(triplet_sigfit1_sigs_melt)
triplet_sigfit1_sigs_melt$centralBP <- paste0(substr(triplet_sigfit1_sigs_melt$variable,2,2),".",substr(triplet_sigfit1_sigs_melt$variable,6,6))
triplet_sigfit1_sigs_melt$centralBP_KHFormat <- triplet_sigfit1_sigs_melt$centralBP
triplet_sigfit1_sigs_melt[triplet_sigfit1_sigs_melt$centralBP=="T.G",]$centralBP_KHFormat <- "A.C"
triplet_sigfit1_sigs_melt[triplet_sigfit1_sigs_melt$centralBP=="T.A",]$centralBP_KHFormat <- "A.T"
triplet_sigfit1_sigs_melt[triplet_sigfit1_sigs_melt$centralBP=="T.C",]$centralBP_KHFormat <- "A.G"
triplet_sigfit1_sigs_melt$ancestral3mer <- substr(triplet_sigfit1_sigs_melt$variable,1,3)
# note that this works because central BP has to be C so you won't grep CGX you'll only grep XCG
# but also need to specify that it's C>T 
triplet_sigfit1_sigs_melt[grepl("CG",triplet_sigfit1_sigs_melt$ancestral3mer) & triplet_sigfit1_sigs_melt$centralBP=="C.T",]$centralBP_KHFormat <- "C.T_CpG"

triplet_sigfit1_sigs_melt_SUMMEDUPTO1MERS <- triplet_sigfit1_sigs_melt %>%
  group_by(centralBP_KHFormat,signature) %>%
  summarise(summedUpFraction=sum(value))
# widen and make in same order as de novo signatures:
triplet_sigfit1_sigs_melt_SUMMEDUPTO1MERS_wide <- pivot_wider(triplet_sigfit1_sigs_melt_SUMMEDUPTO1MERS,names_from = centralBP_KHFormat,values_from = summedUpFraction)

triplet_sigfit1_sigs_melt_SUMMEDUPTO1MERS_wide <- triplet_sigfit1_sigs_melt_SUMMEDUPTO1MERS_wide %>%
  remove_rownames %>% # have to remove any preexisting rownames 
  column_to_rownames(var="signature") # adds row names
rownames(triplet_sigfit1_sigs_melt_SUMMEDUPTO1MERS_wide) <- paste0(rownames(triplet_sigfit1_sigs_melt_SUMMEDUPTO1MERS_wide)," from 3mers")
# make in same order as de novo 1mer signatures:
triplet_sigfit1_sigs_melt_SUMMEDUPTO1MERS_wide <- data.frame(triplet_sigfit1_sigs_melt_SUMMEDUPTO1MERS_wide[,names(deNovoSigs_best_signatures)])

# okay now compare
cosineSimlarity_comparingDeNovo3merToDenovo1mer <- GetCosineSimilarity_betweenSignatures(triplet_sigfit1_sigs_melt_SUMMEDUPTO1MERS_wide,deNovoSigs_best_signatures)

cosineSimlarity_comparingDeNovo3merToDenovo1mer_PLOT <- ggplot(cosineSimlarity_comparingDeNovo3merToDenovo1mer,aes(x=de_novo_signature,y=a_priori_signature,fill=cosine_similarity))+
  geom_tile()+
  scale_fill_viridis_c()+
  xlab("de novo 3mer signatures projected down to 1mers")+
  ylab("de novo 1mer signatures")+
  geom_text(aes(label=round(cosine_similarity,2)))
cosineSimlarity_comparingDeNovo3merToDenovo1mer_PLOT
ggsave(paste0(plotdir,"cosineSimlarity_comparingDeNovo3merToDenovo1mer_PLOT.png"),cosineSimlarity_comparingDeNovo3merToDenovo1mer_PLOT)




######## fit cosmic signatures ########
triplet_sigfit2 <- run_plot_sigfit_FITTING(counts=spectrum_3mer_forSigfit_downsampledCounts,signaturesToFit = cosmic_signatures_v3.2,signatureslabel = "cosmicSignatures",opportunities = NULL,wd=outdir_3mer,niter=niter,model=model,note="fitting_cosmic_signatures")
allcosineSims_3mers = bind_rows(allcosineSims,triplet_sigfit2)

triplet_sigfit2_results <- readRDS(paste0(outdir_3mer,"sigfit.fit.cosmicSignatures.signatures.fitting_cosmic_signatures.",model,".iter.",niter,"/sigfit.fit.cosmicSignatures.signatures.fitting_cosmic_signatures.",model,".iter.",niter,".rds"))

triplet_sigfit2_results_reconstructions <- retrieve_pars(triplet_sigfit2_results,par="reconstructions")
# want to correlate exposures with AFR:
triplet_sigfit2_results_exposures <- retrieve_pars(triplet_sigfit2_results,par="exposures")$mean


# combine with reproductive information
triplet_sigfit2_results_exposures$label <- rownames(triplet_sigfit2_results_exposures)

triplet_sigfit2_results_exposures_PlusReproductiveAges <- merge(triplet_sigfit2_results_exposures,reproductiveAges,by="label")

triplet_sigfit2_results_exposures_melt <- melt(select(triplet_sigfit2_results_exposures,c(starts_with("SBS"),"label")))
# plot for each signature to look for species who stand out
#facet wrap by signature?
triplet_sigfit2_results_exposures_melt$label <- factor(triplet_sigfit2_results_exposures_melt$label,levels=c("bears_EUR","bears_ABC","bears_PB","fin_whale_GOC","fin_whale_ENP","vaquita","humans_EUR"  ,"humans_EAS","humans_SAS","humans_AMR","humans_AFR" , "apes_Pan_paniscus", "apes_Pan_troglodytes" ,"apes_Gorilla"  , "apes_Pongo_pygmaeus"  , "apes_Pongo_abelii","mice_Mmd","mice_Mmm","mice_Mmc","mice_Ms" ))

cosmicExposuresPerSpeciesPlot <- ggplot(triplet_sigfit2_results_exposures_melt,aes(x=label,y=variable,fill=log10(value)))+
  geom_tile()+
  theme_bw()+
  scale_fill_viridis_c(limits=c(-3,0))+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  ggtitle("Exposure of COSMIC signatures across species from sigfit\nexposures below 0.001 are grayed out ")
cosmicExposuresPerSpeciesPlot
ggsave(paste0(outdir_3mer,"cosmicexposuresperspecies.png"),cosmicExposuresPerSpeciesPlot,height=16,width=16)

# interested in mouse sbs30 -- does it explain the excess C>T
triplet_sigfit2_results_reconstructions <- data.frame(triplet_sigfit2_results_reconstructions)

triplet_sigfit2_results_reconstructions$label <- rownames(triplet_sigfit2_results_reconstructions)

triplet_sigfit2_results_reconstructions_melt <- melt(triplet_sigfit2_results_reconstructions)

ggplot(triplet_sigfit2_results_reconstructions_melt,aes(x=variable,))

########## plot SBS5 exposure with aging ##########
sbs5exposureplot_AFR <- ggplot(triplet_sigfit2_results_exposures_PlusReproductiveAges,aes(x=AFR_d,y=SBS5,label=label))+
  geom_point()+
  geom_text(size=2)+
  ggtitle('exposure to SBS5 (3mer) vs age at first reproduction\nnot phylogenetically corrected')
sbs5exposureplot_AFR
ggsave(paste0(outdir_3mer,"expsosureToSBS5.AgeAtFirstReproduction.png"),sbs5exposureplot_AFR)

sbs1exposureplot_Rspan <- ggplot(triplet_sigfit2_results_exposures_PlusReproductiveAges,aes(x=Rspan_d,y=SBS1,label=label))+
  geom_point()+
  geom_text(size=2)+
  ggtitle('exposure to SBS1 (3mer) vs reproductive span\nnot phylogentically corrected')
sbs1exposureplot_Rspan
ggsave(paste0(outdir_3mer,"expsosureToSBS1.ReproductiveSpan.png"),sbs1exposureplot_Rspan)


sbs1exposureplot_AFR <- ggplot(triplet_sigfit2_results_exposures_PlusReproductiveAges,aes(x=AFR_d,y=SBS1,label=label))+
  geom_point()+
  geom_text(size=2)+
  ggtitle('exposure to SBS1 (3mer) vs age at first reproduction\nnot phylogenetically corrected')
sbs1exposureplot_AFR
ggsave(paste0(outdir_3mer,"expsosureToSBS1.AgeAtFirstReproduction.png"),sbs1exposureplot_AFR)

sbs5exposureplot_Rspan <- ggplot(triplet_sigfit2_results_exposures_PlusReproductiveAges,aes(x=Rspan_d,y=SBS5,label=label))+
  geom_point()+
  geom_text(size=2)+
  ggtitle('exposure to SBS5 (3mer) vs reproductive span\nnot phylogentically corrected')
sbs5exposureplot_Rspan
ggsave(paste0(outdir_3mer,"expsosureToSBS5.ReproductiveSpan.png"),sbs5exposureplot_Rspan)

############ fit SBS5 and SBS1 then extract more de novo signatures ############
triplet_sigfit_fitSBS15_extract <- run_plot_sigfit_FIT_EXT(counts=spectrum_3mer_forSigfit_downsampledCounts,signaturesToFit = cosmic_signatures_v3.2[c("SBS1","SBS5"),],ndenovo = 2, signatureslabel = "cosmic_SBS1_SBS5",opportunities = NULL,wd=outdir_3mer,niter=niter,model=model,note="fitSBS1_SBS5_extract2DeNovo")

allcosineSims_3mers = bind_rows(allcosineSims,triplet_sigfit_fitSBS15_extract)

triplet_sigfit_fitSBS15_extract_RESULTS <- readRDS(paste0(outdir_3mer,"/sigfit.fit.cosmic_SBS1_SBS5.signatures.fitSBS1_SBS5_extract2DeNovo.",model,".iter.",niter,"/sigfit.fit.cosmic_SBS1_SBS5.signatures.fitSBS1_SBS5_extract2DeNovo.",model,".iter.",niter,".rds"))

triplet_sigfit_fitSBS15_extract_RESULTS_sigs <- retrieve_pars(triplet_sigfit_fitSBS15_extract_RESULTS,par='signatures')$mean

FittingSBS1_SBS5_Denovo_CosineSimToCosmic <- GetCosineSimilarity_betweenSignatures(triplet_sigfit_fitSBS15_extract_RESULTS_sigs,cosmic_signatures_v3.2)

FittingSBS1_SBS5_Denovo_CosineSimToCosmic$a_priori_signature <- factor(FittingSBS1_SBS5_Denovo_CosineSimToCosmic$a_priori_signature, levels=mixedsort(unique(FittingSBS1_SBS5_Denovo_CosineSimToCosmic$a_priori_signature)))
# sig a = sbs 1 ; sig b = sbs 5; so only want to plot c and d 
denovo3mersignatures <- ggplot(FittingSBS1_SBS5_Denovo_CosineSimToCosmic[FittingSBS1_SBS5_Denovo_CosineSimToCosmic$de_novo_signature %in% c("Signature C","Signature D"),],aes(x=de_novo_signature,y=a_priori_signature,fill=cosine_similarity))+
  geom_tile()+
  scale_fill_viridis_c()+
  theme_bw()+
  geom_text(aes(label=round(cosine_similarity,2)),size=2)
ggsave(paste0(outdir_3mer,"comparingDeNovoSignaturesExtractedAfterSBS1.SBS5fit.CompareTOCosmic.png"),denovo3mersignatures,height=12,width=9)

############ fit *JUST* SBS5 and SBS1 ############
justSBS1_SBS5_3mer <- run_plot_sigfit_FITTING(counts=spectrum_3mer_forSigfit_downsampledCounts,signaturesToFit = cosmic_signatures_v3.2[c("SBS1","SBS5"),],signatureslabel = "SBS1and5",opportunities = NULL,wd=outdir_3mer,niter=niter,model=model,note="fittingJustSBS1_and_SBS5")

justSBS1_SBS5_3merPlot <- ggplot(justSBS1_SBS5_3mer,aes(y=species,x=cosine_similarity))+
  geom_point()+
  ggtitle("Just sbs5, sbs1 -- 3mers")
ggsave(paste0(outdir_3mer,"justSBS5.SBS1.3mersVersion.png"),justSBS1_SBS5_3merPlot)
############### TROUBLESHOOTING: DIVERGENT WARMUPS WARNING -- fixed ##############
# TRYING TO get rid of warnigns
# try making adapt_delta (step size) -- smaller
# will make it take longer
#test1 <- sigfit::extract_signatures(counts=spectrum_3mer_forSigfit_downsampledCounts,nsignatures =5,control=list(adapt_delta=0.85)) # https://mc-stan.org/users/documentation/case-studies/divergences_and_bias.html # this is with default adapt delta/ did  not have divergent warmup warning. that's interesting
# see what warnings look like
#test2<- sigfit::extract_signatures(counts=spectrum_3mer_forSigfit_downsampledCounts,nsignatures =5,control=list(adapt_delta=0.99)) # this is with higher 
# try with a huge number of iterations?

#test4 <- sigfit::fit_extract_signatures(counts=spectrum_3mer_forSigfit_downsampledCounts,signatures = cosmic_signatures_v3.2[c("SBS1","SBS5"),],num_extra_sigs = 2,control=list(adapt_delta=0.99),iter=niter) # this doesn't have divergent warmup warning (good); but needs more iterations 
#plot_all(test4,out_path = paste0(outdir_3mer,"TEST4.tryingFitExtractSBS1.5.WithAdaptDelta99"))
#test3<- sigfit::fit_extract_signatures(counts=spectrum_3mer_forSigfit_downsampledCounts,signatures = cosmic_signatures_v3.2[c("SBS1","SBS5"),],num_extra_sigs = 2,control=list(adapt_delta=0.85),iter=niter) # okay this did have the divergent warmup issue! so we do need the adapt delta 
#plot_all(test3,out_path = paste0(outdir_3mer,"TEST3.tryingFitExtractSBS1.5.WithAdaptDelta85"))
# okay based on these tests it looks like upping adapt_delta 
# does this have bad warnings?
# need to rerun everything with 10K perhaps? more? 
# maybe my snp numbers are too high?
######## is sbs5 similar to aging signatures if I reduce it to 1mers? ##########
sbs5 <- cosmic_signatures_v3.2["SBS5",]
# 
sbs5_melt <- melt(sbs5)
head(sbs5_melt)
sbs5_melt$triplet <- rownames(sbs5_melt)
sbs5_melt$centralBP <- paste0(substr(sbs5_melt$triplet,2,2),".",substr(sbs5_melt$triplet,6,6))
sbs5_melt$centralBP_KHFormat <- sbs5_melt$centralBP
sbs5_melt[sbs5_melt$centralBP=="T.G",]$centralBP_KHFormat <- "A.C"
sbs5_melt[sbs5_melt$centralBP=="T.A",]$centralBP_KHFormat <- "A.T"
sbs5_melt[sbs5_melt$centralBP=="T.C",]$centralBP_KHFormat <- "A.G"
sbs5_melt$ancestral3mer <- substr(sbs5_melt$triplet,1,3)
# note that this works because central BP has to be C so you won't grep CGX you'll only grep XCG
# but also need to specify that it's C>T 
sbs5_melt[grepl("CG",sbs5_melt$ancestral3mer) & sbs5_melt$centralBP=="C.T",]$centralBP_KHFormat <- "C.T_CpG"




sbs5_SUMMEDUPTO1MERS <- sbs5_melt %>%
  group_by(centralBP_KHFormat) %>%
  summarise(summedUpFraction=sum(value))
sbs5_SUMMEDUPTO1MERS_wide <- pivot_wider(sbs5_SUMMEDUPTO1MERS,names_from = centralBP_KHFormat,values_from = summedUpFraction)
# make in same order as aging:
sbs5_SUMMEDUPTO1MERS_wide <- data.frame(sbs5_SUMMEDUPTO1MERS_wide[,names(agingSignatures)])
rownames(sbs5_SUMMEDUPTO1MERS_wide) <- "sbs5_summed_over1mers"

sbs5cosineSim <- GetCosineSimilarity_betweenSignatures(sbs5_SUMMEDUPTO1MERS_wide,agingSignatures)

sbs5cosineSimPlot <- ggplot(sbs5cosineSim,aes(x=a_priori_signature,de_novo_signature,fill=cosine_similarity))+
  geom_tile()+
  scale_fill_viridis_c()+
  geom_text(aes(label=round(cosine_similarity,2)))+
  ylab("")+
  xlab("")
sbs5cosineSimPlot
ggsave(paste0(plotdir,"sbs5.summedUpOver1mers.ComparedToAgingSigs.png"),sbs5cosineSimPlot)

sbs5_SUMMEDUPTO1MERS_wide$label <- "SBS5_summedOver1mers"
agingSignatures_addLabel <- agingSignatures
agingSignatures_addLabel$label <- rownames(agingSignatures)
aging_plus_sbs5 <- bind_rows(agingSignatures_addLabel,sbs5_SUMMEDUPTO1MERS_wide)

aging_plus_sbs5_melt <- melt(aging_plus_sbs5)

sbs5colPlot <- ggplot(aging_plus_sbs5_melt,aes(x=mixedsort(variable),y=value,fill=label))+
  geom_col(position="dodge")
sbs5colPlot
ggsave(paste0(plotdir,"sbs5.summedUpOver1mers.ComparedToAgingSigs.columns.png"),sbs5colPlot)

########### do this for all  cosmic signatures ###############
# 
cosmic_melt <- melt(cosmic_signatures_v3.2)
head(cosmic_melt)
colnames(cosmic_melt) <- c("signature","triplet","value")
cosmic_melt$centralBP <- paste0(substr(cosmic_melt$triplet,2,2),".",substr(cosmic_melt$triplet,6,6))
cosmic_melt$centralBP_KHFormat <- cosmic_melt$centralBP
cosmic_melt[cosmic_melt$centralBP=="T.G",]$centralBP_KHFormat <- "A.C"
cosmic_melt[cosmic_melt$centralBP=="T.A",]$centralBP_KHFormat <- "A.T"
cosmic_melt[cosmic_melt$centralBP=="T.C",]$centralBP_KHFormat <- "A.G"
cosmic_melt$ancestral3mer <- substr(cosmic_melt$triplet,1,3)
# note that this works because central BP has to be C so you won't grep CGX you'll only grep XCG
# but also need to specify that it's C>T 
cosmic_melt[grepl("CG",cosmic_melt$ancestral3mer) & cosmic_melt$centralBP=="C.T",]$centralBP_KHFormat <- "C.T_CpG"

### sum up over central bp:
cosmic_SUMMEDUPTO1MERS <- cosmic_melt %>%
  group_by(signature,centralBP_KHFormat) %>%
  summarise(summedUpFraction=sum(value))


cosmic_SUMMEDUPTO1MERS_wide <- data.frame(pivot_wider(cosmic_SUMMEDUPTO1MERS,names_from = centralBP_KHFormat,values_from = summedUpFraction,id_cols = signature))
# make in same order as aging:
rownames(cosmic_SUMMEDUPTO1MERS_wide) <- cosmic_SUMMEDUPTO1MERS_wide$signature

# this sets columns in same order as aging signatures and gets rid of 'signature ' column
cosmic_SUMMEDUPTO1MERS_wide <- data.frame(cosmic_SUMMEDUPTO1MERS_wide[,names(agingSignatures)])

cosmic1mer_vs_aging_cosineSim <- GetCosineSimilarity_betweenSignatures(agingSignatures,cosmic_SUMMEDUPTO1MERS_wide)

cosmic1mer_vs_aging_cosineSim$a_priori_signature <- factor(cosmic1mer_vs_aging_cosineSim$a_priori_signature,levels=mixedsort(unique(cosmic1mer_vs_aging_cosineSim$a_priori_signature)))

cosmic1mer_vs_aging_cosineSimPlot <- ggplot(cosmic1mer_vs_aging_cosineSim,aes(y=a_priori_signature,x=de_novo_signature,fill=cosine_similarity))+
  geom_tile()+
  scale_fill_viridis_c()+
  geom_text(aes(label=round(cosine_similarity,2)),size=1.5)+
  ylab("")+
  xlab("")+
  ggtitle("similarity of sbs signatures summed over central BP to parental age signatures")
cosmic1mer_vs_aging_cosineSimPlot
ggsave(paste0(plotdir,"comparisongOfCosmicSignaturesReducedTo1mers.ToAgingSignatures.png"),cosmic1mer_vs_aging_cosineSimPlot,height=12,width=8)

############## fit JUST sbs5-1mer and sbs1-mer ###########
justSBS1_SBS5_1mer <- run_plot_sigfit_FITTING(counts=spectrum_1mer_forSigfit_downsampledCounts,signaturesToFit = cosmic_SUMMEDUPTO1MERS_wide[c("SBS1","SBS5"),],signatureslabel = "SBS1and5_1merVersion",opportunities = NULL,wd=plotdir,niter=niter,model=model,note="fittingJustSBS1_and_SBS5_1merVersion")

justSBS1_SBS5_1mer_results <- readRDS(paste0(plotdir,"sigfit.fit.SBS1and5_1merVersion.signatures.fittingJustSBS1_and_SBS5_1merVersion.",model,".iter.",niter,"/sigfit.fit.SBS1and5_1merVersion.signatures.fittingJustSBS1_and_SBS5_1merVersion.",model,".iter.",niter,".rds"))

justSBS1_SBS5_1merPlot <- ggplot(justSBS1_SBS5_1mer,aes(y=species,x=cosine_similarity))+
  geom_point()+
  ggtitle("Just sbs5 + sbs1 --- 1mer version")
ggsave(paste0(plotdir,"JustSBS1_and_SBS5_1merversion.png"),justSBS1_SBS5_1merPlot)
############# try fitting every possible pair of 1mer version cosmic signatures -- VERY SLOW> did once.#############

#outdir_PairsOfCosmic = paste0(plotdir,"comparingAllPairsOfCosmicSignatures/")
#dir.create(outdir_PairsOfCosmic,showWarnings = F)
###### THIS TAKES A LONG TIME (>18 hours )
#allPairsCosineSimilarities = data.frame()
#cosmic_possible_pairs <- combn(rownames(cosmic_SUMMEDUPTO1MERS_wide),2,simplify = F)
# don't want duplicates 
#length(cosmic_possible_pairs) # 3003
# for(i in seq(1,length(cosmic_possible_pairs))){
#   pairOfSigs=cosmic_possible_pairs[[i]]
#   sig1Name=pairOfSigs[1]
#   sig2Name=pairOfSigs[2]
#   print(paste0("starting ",sig1Name,"+ ",sig2Name))
#   
#   cos_intermediate <- run_plot_sigfit_FITTING_NOPLOTTING(counts=spectrum_1mer_forSigfit_downsampledCounts,signaturesToFit = cosmic_SUMMEDUPTO1MERS_wide[c(sig1Name,sig2Name),],signatureslabel = paste0(sig1Name,"_",sig2Name),opportunities = NULL,wd=outdir_PairsOfCosmic,niter=2000,model=model,note="pairsOfSigs") # starting with 2000 iterations here so that it's fast
#   # update funciton ; doesn't need all the outdir stuff actually
#   cos_intermediate$sig1 <- sig1Name
#   cos_intermediate$sig2 <- sig2Name
#   cos_intermediate$label <- paste0(pmin(sig1Name,sig2Name),"_",pmax(sig1Name,sig2Name))
#   allPairsCosineSimilarities = bind_rows(allPairsCosineSimilarities,cos_intermediate)
# }
# # write.table(allPairsCosineSimilarities,paste0(outdir_PairsOfCosmic,"pairsOfCosmic1merSignatures.CosineSimilarities.allpairs.txt"),row.names = F,quote=F,sep="\t")
# 
# allPairsCosineSimilarities <- read.table(paste0(outdir_PairsOfCosmic,"pairsOfCosmic1merSignatures.CosineSimilarities.allpairs.txt"),header=T)
# 
# # label if it involves sbs5:
# allPairsCosineSimilarities$colorLabel1 <- ""
# allPairsCosineSimilarities[allPairsCosineSimilarities$sig1=="SBS5" | allPairsCosineSimilarities$sig2=="SBS5",]$colorLabel1 <- "contains_SBS5"
# 
# allPairsCosineSimilarities[allPairsCosineSimilarities$label=="SBS1_SBS5",]$colorLabel1 <- "SBS1_SBS5"
# 
# 
# # find the top for each species:
# 
# allPairsCosineSimilaritiesPlot <- ggplot(allPairsCosineSimilarities,aes(x=cosine_similarity,fill=colorLabel1))+
#   geom_histogram(bins=100)+
#   #geom_point(data=allPairsCosineSimilarities[allPairsCosineSimilarities$colorLabel1=="contains_SBS5",],aes(x=cosine_similarity,y=1))+
#   facet_wrap(~species)+
#   ggtitle("as in Rahbari, looking at explanatory value of every pair of SBS signatures (but in their 1mer form!)")
# ggsave(paste0(outdir_PairsOfCosmic,"PairsOfCosmicSignatures.1merVersion.ExplanatoryPower.png"))
# 
# allPairsCosineSimilarities %>%
#   group_by(species) %>%
#   filter(cosine_similarity==max(cosine_similarity)) %>%
#   ggplot(.,aes(x=label))+
#   geom_histogram(stat="count",aes(fill=species))
# 
# allPairsCosineSimilarities %>%
#   filter(cosine_similarity>=0.95) %>%
#   ggplot(.,aes(x=cosine_similarity,fill=colorLabel1))+
#   geom_histogram(bins=100)+
#   facet_wrap(~species)
# 
# allPairsCosineSimilarities %>%
#   filter(cosine_similarity>=0.99) %>%
#   ggplot(.,aes(x=species,y=label))+
#   geom_point()+
#   facet_wrap(~species,scales="free")
# 
# 
# allPairsCosineSimilarities %>%
#   group_by(species) %>%
#   filter(cosine_similarity==max(cosine_similarity))  %>%
#   select(c(species,cosine_similarity,label)) %>%
#   write.table(.,paste0(outdir_PairsOfCosmic,"MaxCosineSimBasedOnPairsOFSignatures.txt"),row.names=F,quote=F,sep="\t")
# 

########### !!! sandbox !!!hihgly experimental: subtract sbs5 from parental signatures ############
# kelley suggests rescaling sbs5 so that only one mutaiton type goes to zero, not multiple.
# how best to do that... 
# rescale sbs5 so that A.G matches maternal age -- get a scaling factor that is maternal A.G / young parent A.G and then multiply sbs5 that to scale things down so that only A>G cancels out
### STUCK ON THIS -- talk to kelley
# then set anything negative to 0 and renormalize
# need to repeat sbs5 wide 3x so taht there's a row each time 
agingSignaturesMinusSBS51mer <- agingSignatures - sbs5_SUMMEDUPTO1MERS_wide[rep(1,3),names(agingSignatures)] # this worked! nice! https://stackoverflow.com/questions/35892111/subtract-one-row-from-another-row-in-df
agingSignaturesMinusSBS51mer$label <- rownames(agingSignaturesMinusSBS51mer)
# want to make any entry that is negative into 0
agingSignaturesMinusSBS51mer <- melt(agingSignaturesMinusSBS51mer)
head(agingSignaturesMinusSBS51mer)
agingSignaturesMinusSBS51mer$value_new_negativeTo0_CpGTo0 <- agingSignaturesMinusSBS51mer$value
agingSignaturesMinusSBS51mer[agingSignaturesMinusSBS51mer$value<0,]$value_new_negativeTo0_CpGTo0 <- 0
# really not sure about this. 
agingSignaturesMinusSBS51mer[agingSignaturesMinusSBS51mer$variable=="C.T_CpG",]$value_new_negativeTo0_CpGTo0 <- 0

### FOR NOW ALSO SETTING CpG RATE to 1 (like subtracting off SBS1)
# now rescale to 1
agingSignaturesMinusSBS51mer <- agingSignaturesMinusSBS51mer %>%
  group_by(label) %>%
  mutate(RENORMALIZED_value_new_negativeTo0_CpGTo0 =value_new_negativeTo0_CpGTo0 / sum(value_new_negativeTo0_CpGTo0) )
# hmm so they end up all cpgs. so also want to subtract of SBS1 in some manner ? for now setting CpG to 0 before renormalization. so then paternal age signature just becomes A.C and maternal is mostly C.G and puberty is all C.T and A.T 
agingSignaturesMinusSBS51mer
ggplot(agingSignaturesMinusSBS51mer,aes(x=variable,y=RENORMALIZED_value_new_negativeTo0_CpGTo0,fill=label))+
  geom_col(position='dodge')+
  ggtitle("subtracted sbs5 proportions from each aging signature. set any negative values to 0. set CpG>TpG to 0 (like subtracting off SBS1? not sure)\nrenormalized")

agingSignaturesMinusSBS51mer_wide <- pivot_wider(agingSignaturesMinusSBS51mer,id_cols=c(label),names_from = variable,values_from =RENORMALIZED_value_new_negativeTo0_CpGTo0 )
# could try to fit these signatures? 
agingSignaturesMinusSBS51mer_wide_plusSBS5 <- data.frame(bind_rows(agingSignaturesMinusSBS51mer_wide,sbs5_SUMMEDUPTO1MERS_wide))
rownames(agingSignaturesMinusSBS51mer_wide_plusSBS5) <- agingSignaturesMinusSBS51mer_wide_plusSBS5$label
agingSignaturesMinusSBS51mer_wide_plusSBS5 <- select(agingSignaturesMinusSBS51mer_wide_plusSBS5,-label)

TEST_signaturesSubtractSBS5_fit <- run_plot_sigfit_FITTING(counts = spectrum_1mer_forSigfit_downsampledCounts,signaturesToFit = agingSignaturesMinusSBS51mer_wide_plusSBS5,signatureslabel = "agingMinusSBS5",niter=niter,model = model,wd = plotdir,note = "agingMinusSBS5") # interesting -- so this has the sbs5 dominate

############# fit SBS5-1mer version + aging signatures (not subtracted) ########
# keeping sbs5 and aging as separate signatures (not subtracting off sbs5)
sbs5_1merVersion_plusAging <- run_plot_sigfit_FITTING(counts = spectrum_1mer_forSigfit_downsampledCounts,signaturesToFit = select(aging_plus_sbs5,-label),signatureslabel = "agingPlusSBS5",niter=niter,model = model,wd = plotdir,note = "agingplusSBS5")

sbs5_1merVersion_plusAging_results <- readRDS(pste0(plotdir,""))######## ENTER IT HERE 
############### maybe: fit SBS5-1mer + SBS1-1mer by themselves ###########



######### could do like rahbari and see distribution of fitting all pairs of sbs -1mer signatres ? ##########


######## maybe: fit SBS51mer version + subtracted mat/pat sigs + original  young parent signature #########
#sbs1mer_subtractedMatPat_origYoungParent <- bind_rows(agingSignaturesMinusSBS51mer_wide[agingSignaturesMinusSBS51mer_wide$label %in% c("maternal_age","paternal_age"),],data.frame(label="puberty_notSubSBS5_butWithCpG0",agingSignatures["puberty",]),sbs5_SUMMEDUPTO1MERS_wide) # need to change pubery to 'young_Parent'
#sbs1mer_subtractedMatPat_origYoungParent
######## YOU ARE HERE: you are trying this with puberty signature : trying to decide whether to remove CpGs from it. I think yes. ###########
######## FORMALIZE? subtract SBS1 as well? CODE THAT UP? try to do reconstructions with JUST SBS5 next! 
# fit just sbs5 and sbs1 and see how good reconstructions are

sbs5_sbs1_1merVersion_only <- run_plot_sigfit_FITTING(counts = spectrum_1mer_forSigfit_downsampledCounts,signaturesToFit =  cosmic_SUMMEDUPTO1MERS_wide[c("SBS1","SBS5"),],signatureslabel = "sbs5_sbs1_1merOnly",niter=niter,model = model,wd = plotdir,note = "sbs5_sbs1_1merVersions_only")


############ >>>>>>>>>>>>>> NICE PLOTS FOR MANUSCRIPT: <<<<<<<<<<<< #################
# ones I want to plot: aging1mer, 1merdenovo, aging SBS1+5, 3mer: sbs1+5, de novo

######## make a list of 1mer results you want to plot ####################

#deNovoSigs_best # de novo signatures (best # of sigs) # it's cosine similarities are in:
#cosinesims2_results # aging signtures (1mer)
#justSBS1_SBS5_1mer_results # just 1mer version of sbs1 and sbs5


# this is a list of the sigfit results 
list_1merResultsToPlotTogether=list("de novo 1-mer signatures"=deNovoSigs_best,"aging signatures"=cosinesims2_results,"SBS1 plus SBS5 (1-mer)"=justSBS1_SBS5_1mer_results)
# here are the cosine similarities (at some point make this all more consistent in namign conventions)
list_1merPlus3merResultsToPlotTogether=list("de novo 1-mer signatures"=deNovoSigs_best,"aging signatures"=cosinesims2_results,"SBS1 plus SBS5 (1-mer)"=justSBS1_SBS5_1mer_results,"de novo 3-mer signatures"=extractTripletsResult)
############### super short species names ##############3
speciesOrder = c("mice_Ms","mice_Mmd","mice_Mmm","mice_Mmc","apes_Pongo_abelii","apes_Pongo_pygmaeus","apes_Gorilla","apes_Pan_paniscus","apes_Pan_troglodytes","humans_AFR","humans_EUR","humans_SAS","humans_EAS","humans_AMR","bears_ABC","bears_EUR","bears_PB","fin_whale_GOC","fin_whale_ENP","vaquita")
speciesShortNames = data.frame(speciesOriginal=speciesOrder,speciesSuperShort=c("Ms","Mmd","Mmm","Mmc","SumO","BorO","Gor","Bon","Chi","H_AFR","H_EUR","H_SAS","H_EAS","H_AMR","BB_ABC","BB_EUR","PB","FW_GOC","FW_ENP","Vaq"))

modelOrder_1mer = c("SBS1 plus SBS5 (1-mer)","aging signatures","de novo 1-mer signatures")
modelOrder_1merPlus3mer = c("SBS1 plus SBS5 (1-mer)","aging signatures","de novo 1-mer signatures","de novo 3-mer signatures")

list_1merResultsToPlotTogether_signatures=list("de novo 1-mer signatures"=retrieve_pars(deNovoSigs_best,par="signatures")$mean,"aging signatures"=data.frame(cosinesims2_results$data$signatures),"SBS1 plus SBS5 (1-mer)"=data.frame(justSBS1_SBS5_1mer_results$data$signatures))

# add rownames
list_1merResultsToPlotTogether_signatures = lapply(list_1merResultsToPlotTogether_signatures,rownames_to_column,var="signature")

# with 3mer sigs included: 
list_1merPlus3merResultsToPlotTogether_signatures=list("de novo 3-mer signatures"=retrieve_pars(extractTripletsResult,par='signatures')$mean,"de novo 1-mer signatures"=retrieve_pars(deNovoSigs_best,par="signatures")$mean,"aging signatures"=data.frame(cosinesims2_results$data$signatures),"SBS1 plus SBS5 (1-mer)"=data.frame(justSBS1_SBS5_1mer_results$data$signatures))

# add rownames
list_1merPlus3merResultsToPlotTogether_signatures = lapply(list_1merPlus3merResultsToPlotTogether_signatures,rownames_to_column,var="signature")

# set order of signatures based on order of models that you want: 
signatureOrder_1mer = c(list_1merResultsToPlotTogether_signatures$`de novo signatures`$signature,list_1merResultsToPlotTogether_signatures$`aging signatures`$signature,list_1merResultsToPlotTogether_signatures$`SBS1 plus SBS5 (1-mer)`$signature)
  
# will need to deal with de novo signatures in 3mer and 1mer being named the same but NOT being the same signatures.

signatureOrder_plus3mer = c(paste0(list_1merPlus3merResultsToPlotTogether_signatures$`de novo 3-mer signatures`$signature," (3-mer)"),paste0(list_1merPlus3merResultsToPlotTogether_signatures$`de novo 1-mer signatures`$signature," (1-mer)"),list_1merPlus3merResultsToPlotTogether_signatures$`aging signatures`$signature,list_1merPlus3merResultsToPlotTogether_signatures$`SBS1 plus SBS5 (1-mer)`$signature)
#  RColorBrewer::brewer.pal(name="Set1",n=8)
# "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF"
# dont include red and green (one or other)
modelColors=c("de novo 1-mer signatures"="#FF7F00","SBS1 plus SBS5 (1-mer)"= "#4DAF4A" ,"aging signatures"="#984EA3" ,"de novo 3-mer signatures"="#377EB8")

signatureColors_incl3mer=c(brewer.pal(name="Dark2",n=4),brewer.pal(name="Set3",n=8))
# reorder species in phylo order 
# 1mer:


########## 1mer cosine comparison plot ##############

# calculate cosine similarities for every sample between reconstruction and data:
list_1merResultsToPlotTogether_cosineSimilarities <- lapply(list_1merResultsToPlotTogether,GetCosineSimilarity_allSamples,runlabel=NA)

df_1merResultsToPlotTogether_cosineSimilarities <- bind_rows(list_1merResultsToPlotTogether_cosineSimilarities,.id="id")

df_1merResultsToPlotTogether_cosineSimilarities$id <- factor(df_1merResultsToPlotTogether_cosineSimilarities$id,levels=modelOrder_1mer)

df_1merResultsToPlotTogether_cosineSimilarities$species <- factor(df_1merResultsToPlotTogether_cosineSimilarities$species,levels=speciesOrder)

cosineComparisonPlot <- ggplot(df_1merResultsToPlotTogether_cosineSimilarities,aes(x=species,y=cosine_similarity,color=id))+
  geom_point(size=4)+ # maybe try to add
  theme_bw()+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle=45,hjust=1),text=element_text(size=14))+
  ggtitle("Cosine similarity between reconstruction and data comparison")+
  scale_color_manual(values=modelColors[c(1,2,3)])+
  xlab("")
cosineComparisonPlot
ggsave(paste0(plotdir,"comparisons.cosineSimilarityComparisonPlot.png"),cosineComparisonPlot,height=4,width=8)

########## 1mer plus 3mer cosine comparison plot ##############

# calculate cosine similarities for every sample between reconstruction and data:
list_1merPlus3merResultsToPlotTogether_cosineSimilarities <- lapply(list_1merPlus3merResultsToPlotTogether,GetCosineSimilarity_allSamples,runlabel=NA)

df_1merPlus3merResultsToPlotTogether_cosineSimilarities <- bind_rows(list_1merPlus3merResultsToPlotTogether_cosineSimilarities,.id="id")

df_1merPlus3merResultsToPlotTogether_cosineSimilarities$id <- factor(df_1merPlus3merResultsToPlotTogether_cosineSimilarities$id,levels=modelOrder_1merPlus3mer)

df_1merPlus3merResultsToPlotTogether_cosineSimilarities$species <- factor(df_1merPlus3merResultsToPlotTogether_cosineSimilarities$species,levels=speciesOrder)

cosineComparisonPlot_1merplus3mer <- ggplot(df_1merPlus3merResultsToPlotTogether_cosineSimilarities,aes(x=species,y=cosine_similarity,color=id))+
  geom_point(size=4,alpha=0.8)+ # maybe try to add
  theme_bw()+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle=45,hjust=1),text=element_text(size=14))+
  ggtitle("Cosine similarity between reconstruction and data comparison")+
  scale_color_manual(values=modelColors)+
  xlab("")
cosineComparisonPlot_1merplus3mer
ggsave(paste0(plotdir,"comparisons.cosineSimilarityComparisonPlot.plus3mer.png"),cosineComparisonPlot_1merplus3mer,height=4,width=8)

########## 1mer signatures comparison plot ##############
# aha so for these ones the signatures that were just fit they are in data not in results; annoying 

# get residuals (my custom function)
# set levels:

df_1merResultsToPlotTogether_signatures =bind_rows(list_1merResultsToPlotTogether_signatures, .id = "id")


df_1merResultsToPlotTogether_signatures_melt <- melt(df_1merResultsToPlotTogether_signatures)
df_1merResultsToPlotTogether_signatures_melt$variable <- gsub("\\.",">",df_1merResultsToPlotTogether_signatures_melt$variable)

df_1merResultsToPlotTogether_signatures_melt$signature <- factor(df_1merResultsToPlotTogether_signatures_melt$signature,levels=signatureOrder)

df_1merResultsToPlotTogether_signatures_melt$id <- factor(df_1merResultsToPlotTogether_signatures_melt$id,levels=modelOrder_1mer)

signaturesComparisonPlot <- ggplot(df_1merResultsToPlotTogether_signatures_melt,aes(x=variable,y=value,fill=signature))+
  geom_col(position="dodge")+
  facet_wrap(~id,ncol=1,scales="free_y")+
  scale_fill_brewer(palette = "Set3")+
  theme_bw()+
  ylab("mutation fraction")+
  xlab("")+
  ggtitle("Signatures")+
  theme(text=element_text(size=14))
signaturesComparisonPlot
ggsave(paste0(plotdir,"comparisons.signaturesComparisonPlot.png"),signaturesComparisonPlot,height=8,width=6)

########### 1mer exposures comparison plot ############

# this  makes a list of each ones' exposures (can I select mean? )
list_1merResultsToPlotTogether_exposures = lapply(list_1merResultsToPlotTogether,retrieve_pars,par="exposures")
# get means: 
list_1merResultsToPlotTogether_exposures_mean = lapply(list_1merResultsToPlotTogether_exposures,"[[","mean")
# make rownames a column: 

list_1merResultsToPlotTogether_exposures_mean= lapply(list_1merResultsToPlotTogether_exposures_mean,rownames_to_column,var="species")


# need to melt and combine:
list_1merResultsToPlotTogether_exposures_mean_melt <- lapply(list_1merResultsToPlotTogether_exposures_mean,melt)

# combine:
df_1merResultsToPlotTogether_exposures_mean_melt <- bind_rows(list_1merResultsToPlotTogether_exposures_mean_melt,.id="id")

# set levels:


df_1merResultsToPlotTogether_exposures_mean_melt$variable <- factor(df_1merResultsToPlotTogether_exposures_mean_melt$variable,levels=c("Signature A","Signature B","Signature C","maternal_age","paternal_age","young_parent_13","SBS1","SBS5"))

df_1merResultsToPlotTogether_exposures_mean_melt$id <- factor(df_1merResultsToPlotTogether_exposures_mean_melt$id,levels=c("SBS1 plus SBS5 (1-mer)","aging signatures","de novo signatures"))

df_1merResultsToPlotTogether_exposures_mean_melt$species <- factor(df_1merResultsToPlotTogether_exposures_mean_melt$species,levels=speciesOrder)

exposuresComparisonPlot <- ggplot(df_1merResultsToPlotTogether_exposures_mean_melt,aes(x=species,y=value,fill=variable))+
  geom_col(position="stack")+
  theme_bw()+
  facet_wrap(~id,ncol=1)+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  scale_fill_brewer(palette = "Set3")+
  theme(legend.title=element_blank(),text=element_text(size=14))+
  ylab("exposure fraction")+
  xlab("")+
  ggtitle("Exposures comparison")
exposuresComparisonPlot
ggsave(paste0(plotdir,"comparisons.exposuresComparisonPlot.png"),exposuresComparisonPlot,height=8,width=8)


######## 1mer + 3mer exposures plot #############
# this  makes a list of each ones' exposures (can I select mean? )
list_1merPlus3merResultsToPlotTogether_exposures = lapply(list_1merPlus3merResultsToPlotTogether,retrieve_pars,par="exposures")
# get means: 
list_1merPlus3merResultsToPlotTogether_exposures_mean = lapply(list_1merPlus3merResultsToPlotTogether_exposures,"[[","mean")
# make rownames a column: 

list_1merPlus3merResultsToPlotTogether_exposures_mean= lapply(list_1merPlus3merResultsToPlotTogether_exposures_mean,rownames_to_column,var="species")


# need to melt and combine:
list_1merPlus3merResultsToPlotTogether_exposures_mean_melt <- lapply(list_1merPlus3merResultsToPlotTogether_exposures_mean,melt)

# combine:
df_1merPlus3merResultsToPlotTogether_exposures_mean_melt <- bind_rows(list_1merPlus3merResultsToPlotTogether_exposures_mean_melt,.id="id")

# make broad comparison categories to facet over
df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$broadCategory <- ""
df_1merPlus3merResultsToPlotTogether_exposures_mean_melt[grepl("mice_",df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$species),]$broadCategory <- "mice"

df_1merPlus3merResultsToPlotTogether_exposures_mean_melt[grepl("apes_",df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$species) | grepl("humans_",df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$species),]$broadCategory <- "great apes (incl. humans)"

df_1merPlus3merResultsToPlotTogether_exposures_mean_melt[grepl("bears_",df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$species),]$broadCategory <- "bears"

df_1merPlus3merResultsToPlotTogether_exposures_mean_melt[grepl("fin_whale",df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$species) | grepl("vaquita",df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$species),]$broadCategory <- "whales"


# set levels:

# need to rename signatures
df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$variable_modelSpecified <- as.character(df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$variable)
# specify if sig A is from 1mer or 3mer model:
df_1merPlus3merResultsToPlotTogether_exposures_mean_melt[df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$id =="de novo 1-mer signatures",]$variable_modelSpecified <- as.character(paste0(df_1merPlus3merResultsToPlotTogether_exposures_mean_melt[df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$id=="de novo 1-mer signatures",]$variable_modelSpecified," (1-mer)"))

##### add in short species labels
df_1merPlus3merResultsToPlotTogether_exposures_mean_melt <- merge(df_1merPlus3merResultsToPlotTogether_exposures_mean_melt,speciesShortNames,by.x="species",by.y="speciesOriginal")

df_1merPlus3merResultsToPlotTogether_exposures_mean_melt[df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$id=="de novo 3-mer signatures",]$variable_modelSpecified <- paste0(df_1merPlus3merResultsToPlotTogether_exposures_mean_melt[df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$id=="de novo 3-mer signatures",]$variable_modelSpecified," (3-mer)")

df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$variable_modelSpecified <- factor(df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$variable_modelSpecified,levels=signatureOrder_plus3mer)

df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$id <- factor(df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$id,levels=modelOrder_1merPlus3mer)

df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$species <- factor(df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$species,levels=speciesOrder)

df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$speciesSuperShort <- factor(df_1merPlus3merResultsToPlotTogether_exposures_mean_melt$speciesSuperShort,levels=speciesShortNames$speciesSuperShort)

exposuresComparisonPlot_1merPlus3mer <- ggplot(df_1merPlus3merResultsToPlotTogether_exposures_mean_melt,aes(x=speciesSuperShort,y=value,fill=variable_modelSpecified))+
  geom_col(position="stack")+
  theme_bw()+
  facet_wrap(~id,ncol=1)+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))+
  scale_fill_manual(values=signatureColors_incl3mer)+
  theme(legend.title=element_blank(),text=element_text(size=14))+
  ylab("exposure fraction")+
  xlab("")+
  ggtitle("Exposures comparison")
exposuresComparisonPlot_1merPlus3mer
ggsave(paste0(plotdir,"comparisons.exposuresComparisonPlot.Includes3mers.png"),exposuresComparisonPlot_1merPlus3mer
,height=11,width=10)
######## want to find some way to demarcate clades
######## 1mer residuals comparison plot ########
list_1merResultsToPlotTogether_residuals = lapply(list_1merResultsToPlotTogether,calculateSigfitResiduals)
# collapse it:
df_1merResultsToPlotTogether_residuals = bind_rows(list_1merResultsToPlotTogether_residuals, .id = "id")
df_1merResultsToPlotTogether_residuals

df_1merResultsToPlotTogether_residuals$id <- factor(df_1merResultsToPlotTogether_residuals$id,levels=c("de novo 1-mer signatures","aging signatures","SBS1 plus SBS5 (1-mer)"))

df_1merResultsToPlotTogether_residuals$mutationLabel <- gsub("\\.",">",df_1merResultsToPlotTogether_residuals$mutationLabel)
# make a plot
residualsComparisonPlot <- ggplot(df_1merResultsToPlotTogether_residuals,aes(x=id,y=residual_data_minus_model,fill=id))+
  geom_hline(yintercept=0,linetype="dashed",color="darkgray")+
  #geom_point()+
  geom_boxplot()+
  facet_wrap(~mutationLabel,scales="free_x",nrow=1)+
  ylab("residual (data - reconstruction)")+
  xlab("")+
  theme_bw()+
  theme(legend.title = element_blank(),text=element_text(size=10))+
  ggtitle("Residuals comparison")+
  scale_fill_manual(values=modelColors[1:3])+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.background = element_rect(fill="white"))
residualsComparisonPlot
ggsave(paste0(plotdir,"comparisons.ResidualsComparisonPlot.png"),residualsComparisonPlot,height=4,width=10)


# 3 mer:
#extractTripletsResult

####################### Distances between Exposures: Mantel test ###########################
# want to take exposures to 3 aging signatures and calculate pairwise distances between them
# need to do ilr/clr transformation
speciesToInclude_phylo=c("mice_Mmd","mice_Ms","bears_ABC","bears_PB","vaquita","fin_whale_GOC","humans_AFR","apes_Gorilla","apes_Pan_paniscus","apes_Pan_troglodytes","apes_Pongo_abelii","apes_Pongo_pygmaeus") 

# making this now take in full sigfit results:
clrTransformExposures_getDistances_forMantelTest <- function(sigfitResultsFull,speciesToInclude,speciesCodes){
  exposures <- retrieve_pars(sigfitResultsFull,par="exposures")$mean # mean exposures
  exposures$label <- rownames(exposures)
  # subset to phylo species: 
  exposures_phyloSpp <- exposures[exposures$label %in% speciesToInclude,]
  exposures_phyloSpp_fullSppNames <- merge(exposures_phyloSpp,speciesCodes,by.x="label",by.y="code")
  rownames(exposures_phyloSpp_fullSppNames) <- exposures_phyloSpp_fullSppNames$species
  exposures_clr <- clr(select(exposures_phyloSpp_fullSppNames,names(select(exposures,-label)))) # need to switch to young_parent -- will throw an error
  print("clr transform (head):")
  print(head(exposures_clr))
  
  # for Mantel test, keep as a full matrix with full species names:
  exposures_clr_dist_MatrixForMantel <- as.matrix(dist(exposures_clr,diag=T,upper=T))
  return(exposures_clr_dist_MatrixForMantel)
  
}

# note for this funciton it will reorder the phylo dist matrix so you don't have to
# and it will make sure it is in right order 
# how to get common names in here 
carryOutMantelTest <- function(phylo_dist_matrix,exposure_dist_matrix_forMantel,npermut,label){
  # reorder phylo dist matrix to match col and row names of exposures:
  phylo_dist_matrix <- as.matrix(phylo_dist_matrix)[rownames(exposure_dist_matrix_forMantel),colnames(exposure_dist_matrix_forMantel)]
  if(sum(colnames(phylo_dist_matrix)!=colnames(exposure_dist_matrix_forMantel))!=0){
    print("something is wrong with ordering of colnames")
    break
  }
  if(sum(rownames(phylo_dist_matrix)!=rownames(exposure_dist_matrix_forMantel))!=0){
    print("something is wrong with ordering of rownames")
    break
  }

  # need to reorder the phylo dist to make sure is in same order as distances 
  mtr = vegan::mantel(xdis=phylo_dist_matrix,ydis=exposure_dist_matrix_forMantel,permutations = npermut,method = "pearson")
  # make a data frame
  mantel_ExposureToAgingSignatures_df <- data.frame(call="vegan_mantel",method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif),label=label)
  
  return(mantel_ExposureToAgingSignatures_df)
}


meltDistanceCombineWithPhyloToPlot <- function(phylo_dist_matrix,exposure_dist_matrix_forMantel,speciesCodes){
  # melt phylo distances and add comparison label
  phylo_dist_melted <- melt_cophenetic_func(phylo_dist_matrix) # melts and adds comparison label 
  exposures_melt <- melt(exposure_dist_matrix_forMantel) # contains symmetrical dups
  dim(exposures_melt)
  colnames(exposures_melt) <- c("sp1","sp2","distance")
  exposures_melt$comparisonLabel <-  paste0(pmin(as.character(exposures_melt$sp1),as.character(exposures_melt$sp1)),".",pmax(as.character(exposures_melt$sp1),as.character(exposures_melt$sp2)))
  # make unique:
  exposures_melt_unique <- unique(exposures_melt[,c("comparisonLabel","distance")])
  dim(exposures_melt_unique)
  # merge with phylo distances
  exposures_and_phylo_dist <- merge(exposures_melt_unique,phylo_dist_melted,by="comparisonLabel") 
  return(exposures_and_phylo_dist)
  
}

# this function returns all the different tables as a list (nice!)
# gives you the distances of expsoures and phylo distances for makign plots
# as well as a matrix of exposure distances used for mantel test (probably not needed)
# table of mantel test results 
combo_MantelTestPlottingScript <- function(sigfitResultsFull,speciesToInclude,speciesCodes,phylo_dist_matrix,npermut,label,outdir){
  # get clr distances based on exposures
  exposureDistancesForMantelTest <- clrTransformExposures_getDistances_forMantelTest(sigfitResultsFull,speciesToInclude,speciesCodes)
  # do mantel test:
  mantelTestResults <- carryOutMantelTest(phylo_dist_matrix,exposureDistancesForMantelTest,npermut,label=label)
  write.table(mantelTestResults,paste0(outdir,label,".clrDistances.mantelTestResults.txt"),row.names=F,quote=F,sep="\t")
  # get distances in a format that's easy to plot: 
  exposureDistancesForPlotting <- meltDistanceCombineWithPhyloToPlot(phylo_dist_matrix,exposureDistancesForMantelTest,speciesCodes)
  write.table(exposureDistancesForPlotting,paste0(outdir,label,".clrDistances.forplotting.txt"),row.names=F,quote=F,sep="\t")
  # make a plot
  DistancePlot <- ggplot(exposureDistancesForPlotting,aes(x=cophenetic_distance,y=distance))+
    geom_point(shape=16,alpha=0.75,color="gray")+
    geom_text(aes(label=comparisonLabel),size=1)+
    theme_bw()+
    ggtitle(label)+
    geom_text(data=mantelTestResults,aes(x=.5,y=Inf,vjust=2,label=paste0("Pearson's r:  ",round(statistic,2),"\np-value: ",as.character(signif)," (Mantel Test; ",permCount," perm.)")),size=5)
  
  print(DistancePlot)
  # save plot:
  ggsave(paste0(outdir,label,".clrDistances.phylodistance.mantel.png"),DistancePlot)
  returnList <- list(all_distancesForPlotting=exposureDistancesForPlotting,mantelResults=mantelTestResults,exposureDistancesForMantelTest=exposureDistancesForMantelTest,plot=DistancePlot)
  return(returnList) #returnst he distances in case you want to make another plot (but would need to read in mantel test?)

}
############### mantel test on aging signatures ###############

exp_distancesOutdir=paste0(plotdir,"exposure_distances/")
dir.create(exp_distancesOutdir)

agingExposureMantelComboResults = combo_MantelTestPlottingScript(sigfitResultsFull=fitToAging,speciesToInclude = speciesToInclude_phylo,speciesCodes = speciesCodes,phylo_dist_matrix = raxml_cophenetic_dist,npermut = mantelTestPermutations,label = "aging1merSignatures",outdir=exp_distancesOutdir)


################## Mantel test: exposures to de novo 1mer signatures ##############
denovo1merSignaturesExposureMantelComboResults <- combo_MantelTestPlottingScript(sigfitResultsFull=deNovoSigs_best,speciesToInclude = speciesToInclude_phylo,speciesCodes = speciesCodes,phylo_dist_matrix = raxml_cophenetic_dist,npermut = mantelTestPermutations,label = "denovo1mersignatures",outdir=exp_distancesOutdir)

################## Mantel test: exposures to de novo 3mer signatures ###############

denovo3merSignaturesExposureMantelComboResults <- combo_MantelTestPlottingScript(sigfitResultsFull=extractTripletsResult,speciesToInclude = speciesToInclude_phylo,speciesCodes = speciesCodes,phylo_dist_matrix = raxml_cophenetic_dist,npermut = mantelTestPermutations,label = "denovo3mersignatures",outdir=exp_distancesOutdir)

# To do: with 1mers,subtract off SBS5 from mat/pat signatures if possible?
########## try to subtract SBS5 from the mat/pat signatures (what if something goes negatie? set to 0?)
######
################ >>>>>>>>>>>> TROUBLESHOOTING: SEPARATE COUNTS BY BGC CATEGORY <<<<<<<<< ###################


spectrum_1mer_forSigfit_downsampledCounts_forBGC <- spectrum_1mer_forSigfit_downsampledCounts
spectrum_1mer_forSigfit_downsampledCounts_forBGC$label <- rownames(spectrum_1mer_forSigfit_downsampledCounts_forBGC)

ggplot(melt(spectrum_1mer_forSigfit_downsampledCounts_forBGC[spectrum_1mer_forSigfit_downsampledCounts_forBGC$label %in% c("mice_Mmd","mice_Mmm","mice_Mmc","mice_Ms","humans_AFR"),]),aes(x=variable,y=value,fill=label))+
  geom_col(position="dodge")

spectrum_1mer_forSigfit_downsampledCounts_forBGC <- melt(spectrum_1mer_forSigfit_downsampledCounts_forBGC)


spectrum_1mer_forSigfit_downsampledCounts_forBGC$BGC_label <- ""

spectrum_1mer_forSigfit_downsampledCounts_forBGC[spectrum_1mer_forSigfit_downsampledCounts_forBGC$variable %in% c("A.C","A.G"),]$BGC_label <- "WtS"

spectrum_1mer_forSigfit_downsampledCounts_forBGC[spectrum_1mer_forSigfit_downsampledCounts_forBGC$variable %in% c("C.T","C.A"),]$BGC_label <- "StW_exclCpG"


spectrum_1mer_forSigfit_downsampledCounts_forBGC[spectrum_1mer_forSigfit_downsampledCounts_forBGC$variable %in% c("C.T_CpG"),]$BGC_label <- "C.T_CpG"

spectrum_1mer_forSigfit_downsampledCounts_forBGC[spectrum_1mer_forSigfit_downsampledCounts_forBGC$variable %in% c("C.G","A.T"),]$BGC_label <- "BGC_conserved"

spectrum_1mer_forSigfit_downsampledCounts_forBGC_summarized <-  spectrum_1mer_forSigfit_downsampledCounts_forBGC %>%
  group_by(label,BGC_label) %>%
  summarise(BGC_count_humanTargets =sum(value)) %>%
  ungroup() %>%
  group_by(label) %>%
  mutate(BGC_proportion_humanTargets=BGC_count_humanTargets/sum(BGC_count_humanTargets))

ggplot(spectrum_1mer_forSigfit_downsampledCounts_forBGC_summarized,aes(x=BGC_label,y=BGC_proportion_humanTargets,fill=label))+
  geom_col(position="dodge")

ggplot(spectrum_1mer_forSigfit_downsampledCounts_forBGC_summarized[spectrum_1mer_forSigfit_downsampledCounts_forBGC_summarized$label %in% c("mice_Mmd","mice_Mmm","mice_Mmc","mice_Ms","humans_AFR"),],aes(x=BGC_label,y=BGC_proportion_humanTargets,fill=label))+
  geom_col(position="dodge")

###### what about mispol? #########
# if you make it a matrix then rownames get included in melt
spectrum_1mer_forSigfit_downsampledCounts_forMisPol <- melt(as.matrix(spectrum_1mer_forSigfit_downsampledCounts))
spectrum_1mer_forSigfit_downsampledCounts_forMisPol$misPolLabel <- ""

spectrum_1mer_forSigfit_downsampledCounts_forMisPol[spectrum_1mer_forSigfit_downsampledCounts_forMisPol$Var2 %in% c("A.C","C.A"),]$misPolLabel <- "A.C_and_C.A"

spectrum_1mer_forSigfit_downsampledCounts_forMisPol[spectrum_1mer_forSigfit_downsampledCounts_forMisPol$Var2 %in% c("C.G"),]$misPolLabel <- "C.G"

spectrum_1mer_forSigfit_downsampledCounts_forMisPol[spectrum_1mer_forSigfit_downsampledCounts_forMisPol$Var2 %in% c("A.T"),]$misPolLabel <- "A.T"

spectrum_1mer_forSigfit_downsampledCounts_forMisPol[spectrum_1mer_forSigfit_downsampledCounts_forMisPol$Var2 %in% c("A.G","C.T","C.T_CpG"),]$misPolLabel <- "A.G_and_C.T_inclCpG"

spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized <- spectrum_1mer_forSigfit_downsampledCounts_forMisPol %>%
  group_by(Var1,misPolLabel) %>%
  summarise(misPolLabel_count=sum(value))

ggplot(spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized[spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized$Var1 %in% c("mice_Mmd","mice_Mmm","mice_Mmc","mice_Ms","humans_AFR"),],aes(x=misPolLabel,y=misPolLabel_count,fill=Var1))+
  geom_col(position="dodge")


ggplot(spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized,aes(x=misPolLabel,y=misPolLabel_count,fill=Var1))+
  geom_col(position="dodge")

#### now get proportions and distances :
spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized <- spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized %>%
  group_by(Var1) %>%
  mutate(misPolLabel_fraction=misPolLabel_count/sum(misPolLabel_count))

spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide <- data.frame(pivot_wider(spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized,id_cols = c("Var1"),names_from = misPolLabel,values_from=misPolLabel_fraction))

rownames(spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide) <- spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide$Var1

spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide_clr <- clr(select(spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide,-Var1))

heatmap(as.matrix(dist(spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide_clr)),symm=T)

spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide_clr_dist <- tidy(dist(spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide_clr))

spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide_clr_dist <- addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide_clr_dist,speciesCodes)
#NAs are ok here
spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide_clr_dist <- addPhyloDistances_excludeNas(spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide_clr_dist,phylodistancesdf = phyloDistances)
head(spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide_clr_dist)

# set to species in phylo list
spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide_clr_dist_subset <- spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide_clr_dist[spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide_clr_dist$item1 %in% speciesToInclude_phylo & spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide_clr_dist$item2 %in% speciesToInclude_phylo,]

ggplot(spectrum_1mer_forSigfit_downsampledCounts_forMisPol_summarized_wide_clr_dist_subset,aes(y=distance,x=cophenetic_distance))+
  geom_point(alpha=0.75,color="gray")+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=1.5)+
  theme_bw()+
  ggtitle("collapsed 1mer spectrum experiment")
############## OkAY 
################# redoing distances with human-targets ###############
######## functions from other script, modified to work with the humantarget data
clr_or_ilr_orNoTransform_calculations <- function(pivotedspectradf_plusEpsilonVersionOfVariables,ilr_or_clr_or_none){
  # remove character vectors:
  tableWithoutIDVars <- select(pivotedspectradf_plusEpsilonVersionOfVariables, -c("variable","mutation_label"))
  # AHA so if you run clr just on a matrix it goes ROW-WISE not column wise! annoying! so need to transpose
  #clr_matrix <- clr(t(tableWithoutIDVars)) # if you do it this way it matches when I do it just on on un-tranposed column by itself. if you do it without transposing then you get weird row-wise results that are WRONG. so be cautious here. 
  # so alternative without transposing is to use sapply -- this is a lot nicer I think
  if(ilr_or_clr_or_none=="clr"){
    print("carrying out clr transformation")
    df <- data.frame(sapply(tableWithoutIDVars,clr)) # this way you don't have to transpose and you get a much nicer shaped output
    return(df)
  } else if(ilr_or_clr_or_none=="ilr"){
    print("carrying out ilr transformation")
    
    df <- data.frame(sapply(tableWithoutIDVars,ilr)) # this way you don't have to transpose and you get a much nicer shaped output
    return(df)
    
  } else if(ilr_or_clr_or_none=="none") {
    print("carrying out no transformation")
    df <- tableWithoutIDVars # just return the table without id vars and no transformation
    return(df)
  } else {
    print("invalid option")
    break
  }
  
}
euclideanDistance <- function(pivotedspectradf_noIDVars){
  # make it a matrix with each row as a species (transposed and all numeric)
  # variable in question:
  # distances are row-wise nto column wise so need to transpose 
  distances <- tidy(dist(t(pivotedspectradf_noIDVars))) # transpose for dist otherwise it goes along rows and is wrong
  colnames(distances) <- c("item1","item2","spectrum_distance")
  return(distances)
}

# match with species codes for raxml tree distances
addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf <- function(distancesdf,speciesCodesdf){
  merge1 <- merge(distancesdf,speciesCodesdf,by.x="item1",by.y="code",all.x=TRUE,all.y=FALSE) # only want comparisons which have phylo distances (?)
  merge2 <- merge(merge1,speciesCodesdf,by.x="item2",by.y="code",all.x=TRUE,all.y=FALSE,suffixes=c(".item1",".item2"))
  print("if there are NAs in species name, it's because I don't have it in raxml tree")
  return(merge2)
}

# add in phylo distances
addPhyloDistances_excludeNas <- function(dfWithSpeciesCodesAdded,phylodistancesdf) {
  # remove NAs (no phylo distance present)
  dfWithSpeciesCodesAdded <- na.omit(dfWithSpeciesCodesAdded)
  # add alphabetical comparison label:
  dfWithSpeciesCodesAdded$comparisonLabel <- paste0(pmin(dfWithSpeciesCodesAdded$species.item1,dfWithSpeciesCodesAdded$species.item2),".",pmax(dfWithSpeciesCodesAdded$species.item1,dfWithSpeciesCodesAdded$species.item2))
  # add a common name comparison label (note order may be diff from latin name because alphabetical)
  dfWithSpeciesCodesAdded$comparisonLabel_common_alphabetical <-  paste0(pmin(dfWithSpeciesCodesAdded$common_name.item1,dfWithSpeciesCodesAdded$common_name.item2),".",pmax(dfWithSpeciesCodesAdded$common_name.item1,dfWithSpeciesCodesAdded$common_name.item2))
  
  dfWithSpeciesCodesAdded$comparisonLabel_broad_alphabetical <-  paste0(pmin(dfWithSpeciesCodesAdded$broad_classification.item1,dfWithSpeciesCodesAdded$broad_classification.item2),".",pmax(dfWithSpeciesCodesAdded$broad_classification.item1,dfWithSpeciesCodesAdded$broad_classification.item2))
  
  merge1 <- merge(dfWithSpeciesCodesAdded,phylodistancesdf,by="comparisonLabel")
  return(merge1)
}

combo_function_phylo_humanTargets <- function(spectrumdf,variable,ilr_or_clr_or_none,speciesToInclude,speciesCodes,phyloDistances,wd){
  distance_dataframe <- 
    processSpectra_AdjustForHumanTargets(spectrumdf,wd) %>% # note making this the new process and it needs wd 
    pivotSpectra_perVariable(.,variable)  %>%
    select(c(variable,mutation_label,speciesToInclude)) %>% # just get the cols for the species you want 
    clr_or_ilr_orNoTransform_calculations(.,ilr_or_clr_or_none) %>%
    euclideanDistance(.) %>%
    addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
    addPhyloDistances_excludeNas(.,phyloDistances) %>%
    select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
    mutate(variable=variable,transform_label=ilr_or_clr_or_none)
  
  
  return(distance_dataframe)
  
  
}
############ distance functions for old mantel test ###########
euclideanDistance_keepMatrix_UPPERANDLOWER <- function(pivotedspectradf_noIDVars){
  # make it a matrix with each row as a species (transposed and all numeric)
  # variable in question:
  # distances are row-wise nto column wise so need to transpose 
  distances <- dist(t(pivotedspectradf_noIDVars),diag=T,upper=T) # transpose for dist otherwise it goes along rows and is wrong
  #colnames(distances) <- c("item1","item2","spectrum_distance")
  return(distances)
}
pivotSpectra_perVariable_KEEPFULLSPECIESNAMES <- function(processedspectradf,variableToPivot) {
  processedspectradf_wider <-  pivot_wider(processedspectradf,id_cols=c(mutation_label),names_from=c("species.sppCodes"),values_from =all_of(variableToPivot ))
  processedspectradf_wider$variable <- variableToPivot
  # pivots on a variable 
  return(processedspectradf_wider)
}

comboFunction_mantelTest_forSpectrumDistances_includesPlots <- function(spectrum, phylo_dist,variable,ilr_or_clr_or_none,speciesToInclude,outdir,npermut,label){
  spectrum_dist <- spectrum %>% 
    processSpectra_AdjustForHumanTargets(.,outdir) %>% 
    merge(.,speciesCodes,by.x="label",by.y="code",suffixes=c(".spectrum",".sppCodes")) %>% 
    filter(label %in% speciesToInclude) %>%
    pivotSpectra_perVariable_KEEPFULLSPECIESNAMES(variable) %>%
    clr_or_ilr_orNoTransform_calculations(.,ilr_or_clr_or_none) %>% # note in previous versino of function this was fixed as clr
    euclideanDistance_keepMatrix_UPPERANDLOWER(.)  
  # the matrices need to be ordered in the same way 
  # reorder the time matrix in the same way as the spectrum matrix:
  phylo_dist_reordered <- as.matrix(phylo_dist)[rownames(as.matrix(spectrum_dist)),colnames(as.matrix(spectrum_dist))]
  head(phylo_dist_reordered)
  
  mtr = vegan::mantel(xdis=phylo_dist_reordered,ydis=spectrum_dist,permutations = npermut,method = "pearson")
  # make a data frame
  mantelTestResult_df <- data.frame(call="vegan_mantel",variable=variable,method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif))
  # mantel p value is calc as (ngreater+1)/(length(rsample)+1 so if there are none that are as well correlated then you'd get 0+1 / 999+1 = 1/ 10000 <- 0.001 which is what I get. which means in 999 permutations there are none that are as correlated. am I doing this right?
  DistancesForPlotting <- meltDistanceCombineWithPhyloToPlot(phylo_dist,as.matrix(spectrum_dist),speciesCodes) # can use original or reordered phylo dist here 
  write.table(DistancesForPlotting,paste0(outdir,label,".Distances.forplotting.txt"),row.names=F,quote=F,sep="\t")
  # make a plot
  DistancePlot <- ggplot(DistancesForPlotting,aes(x=cophenetic_distance,y=distance))+
    geom_point(shape=16,alpha=0.75,color="gray")+
    geom_text(aes(label=comparisonLabel),size=1)+
    theme_bw()+
    ggtitle(label)+
    ylab("spectrum distance")+
    geom_text(data=mantelTestResult_df,aes(x=.5,y=Inf,vjust=2,label=paste0("Pearson's r:  ",round(statistic,2),"\np-value: ",as.character(signif)," (Mantel Test; ",permCount," perm.)")),size=5)
  
  print(DistancePlot)
  # save plot:
  ggsave(paste0(outdir,label,".Distances.phylodistance.HUMANTARGETSCORRECTED.mantel.png"),DistancePlot)
  returnList <- list(all_distancesForPlotting=DistancesForPlotting,mantelResults=mantelTestResult_df,spectrumDistancesForMantelTest_matrix=spectrum_dist,plot=DistancePlot)
  return(returnList)
  
}


################ rescale 3mer 5 mer 7mer by human targets ###############

############### >> Get spectrum distances using human corrected proportions << ##############

####### outdir:
spectrum_DistancesOutdir = paste0(plotdir,"/spectrum_distances/")
dir.create(spectrum_DistancesOutdir)
# # regular way of processing:
# trial_distanceDF_basedOnHumanTargetRescaling_1mer <- combo_function_phylo_humanTargets(spectrumdf=all1merSpectra_MERGEDWITHUMANTARGETPROPORTIONS,variable ="frac_seg_sites_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none = "clr",speciesToInclude = speciesToInclude_phylo[!speciesToInclude_phylo %in% c("mice_Mmd","mice_Ms","fin_whale_GOC","bears_PB")],phyloDistances = phyloDistances,speciesCodes = speciesCodes,wd=spectrum_distancesOutdir)
# 
# ggplot(trial_distanceDF_basedOnHumanTargetRescaling_1mer,aes(x=cophenetic_distance,y=spectrum_distance,label=comparisonLabel_common_alphabetical))+
#   geom_point()+
#   geom_text(size=1.5) # okay so this is still good 
# cor(trial_distanceDF_basedOnHumanTargetRescaling_1mer$cophenetic_distance,trial_distanceDF_basedOnHumanTargetRescaling_1mer$spectrum_distance)

# mantel test version:
trial_distanceDF_basedOnHumanTargetRescaling_1mer_MANTEL <- comboFunction_mantelTest_forSpectrumDistances_includesPlots(spectrum =all1merSpectra_MERGEDWITHUMANTARGETPROPORTIONS,variable="frac_seg_sites_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none="clr",phylo_dist =raxml_cophenetic_dist,npermut = mantelTestPermutations,outdir = spectrum_DistancesOutdir,speciesToInclude = speciesToInclude_phylo,label="1merDistance.humanscaledtargets.allSpecies"  ) # note transform is fixed as clr but could make an option in the funciton if desired # fixing as clr -- could set as an option in the function tho; note it needs matrix of distance
# okay so it's still good . 
trial_distanceDF_basedOnHumanTargetRescaling_1mer_MANTEL_exludeMs <-  comboFunction_mantelTest_forSpectrumDistances_includesPlots(spectrum =all1merSpectra_MERGEDWITHUMANTARGETPROPORTIONS,variable="frac_seg_sites_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none="clr",phylo_dist =raxml_cophenetic_dist,npermut = mantelTestPermutations,outdir = spectrum_DistancesOutdir,speciesToInclude = speciesToInclude_phylo[speciesToInclude_phylo!="mice_Ms"],label="1merDistance.humanscaledtargets.noMs"  )


trial_distanceDF_basedOnHumanTargetRescaling_1mer_MANTEL_exludeMs_bears_fw <-  comboFunction_mantelTest_forSpectrumDistances_includesPlots(spectrum =all1merSpectra_MERGEDWITHUMANTARGETPROPORTIONS,variable="frac_seg_sites_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none="clr",phylo_dist =raxml_cophenetic_dist,npermut = mantelTestPermutations,outdir = spectrum_DistancesOutdir,speciesToInclude = speciesToInclude_phylo[!speciesToInclude_phylo %in% c("mice_Ms","bears_ABC","fin_whale_GOC","bears_PB")],label="1merDistance.humanscaledtargets.noMsBearsFW"  )


trial_distanceDF_basedOnHumanTargetRescaling_1mer_MANTEL_exludeMs_bears_vaq <-  comboFunction_mantelTest_forSpectrumDistances_includesPlots(spectrum =all1merSpectra_MERGEDWITHUMANTARGETPROPORTIONS,variable="frac_seg_sites_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none="clr",phylo_dist =raxml_cophenetic_dist,npermut = mantelTestPermutations,outdir = spectrum_DistancesOutdir,speciesToInclude = speciesToInclude_phylo[!speciesToInclude_phylo %in% c("mice_Ms","bears_ABC","vaquita","bears_PB")],label="1merDistance.humanscaledtargets.noMsBearsVaquita"  )


trial_distanceDF_basedOnHumanTargetRescaling_1mer_MANTEL_exludeMs_bears <-  comboFunction_mantelTest_forSpectrumDistances_includesPlots(spectrum =all1merSpectra_MERGEDWITHUMANTARGETPROPORTIONS,variable="frac_seg_sites_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none="clr",phylo_dist =raxml_cophenetic_dist,npermut = mantelTestPermutations,outdir = spectrum_DistancesOutdir,speciesToInclude = speciesToInclude_phylo[!speciesToInclude_phylo %in% c("mice_Ms","bears_ABC","bears_PB")],label="1merDistance.humanscaledtargets.noMsBears"  )

# try just having one species per dataset (deal with batch effects)
trial_distanceDF_basedOnHumanTargetRescaling_1mer_MANTEL_OneSpeciesPerDataset <-  comboFunction_mantelTest_forSpectrumDistances_includesPlots(spectrum =all1merSpectra_MERGEDWITHUMANTARGETPROPORTIONS,variable="frac_seg_sites_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none="clr",phylo_dist =raxml_cophenetic_dist,npermut = mantelTestPermutations,outdir = spectrum_DistancesOutdir,speciesToInclude = speciesToInclude_phylo[speciesToInclude_phylo %in% c("mice_Mmd","bears_ABC","vaquita","fin_whale_GOC","humans_AFR","apes_Pan_paniscus")],label="1merDistance.humanscaledtargets.OneSpeciesPerDataset"  )

############# 3mer distances ####################

trial_distanceDF_basedOnHumanTargetRescaling_3mer_MANTEL <- comboFunction_mantelTest_forSpectrumDistances_includesPlots(spectrum =all3merSpectra_MERGEDWITHUMANTARGETPROPORTIONS,variable="frac_seg_sites_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none="clr",phylo_dist =raxml_cophenetic_dist,npermut = mantelTestPermutations,outdir = spectrum_DistancesOutdir,speciesToInclude = speciesToInclude_phylo,label="3merDistance.humanscaledtargets.allSpecies"  )

trial_distanceDF_basedOnHumanTargetRescaling_3mer_MANTEL_exludeMs_bears <-  comboFunction_mantelTest_forSpectrumDistances_includesPlots(spectrum =all3merSpectra_MERGEDWITHUMANTARGETPROPORTIONS,variable="frac_seg_sites_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none="clr",phylo_dist =raxml_cophenetic_dist,npermut = mantelTestPermutations,outdir = spectrum_DistancesOutdir,speciesToInclude = speciesToInclude_phylo[!speciesToInclude_phylo %in% c("mice_Ms","bears_ABC","bears_PB")],label="3merDistance.humanscaledtargets.noMsBears"  )

trial_distanceDF_basedOnHumanTargetRescaling_3mer_MANTEL_OneSpeciesPerDataset <-  comboFunction_mantelTest_forSpectrumDistances_includesPlots(spectrum =all3merSpectra_MERGEDWITHUMANTARGETPROPORTIONS,variable="frac_seg_sites_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none="clr",phylo_dist =raxml_cophenetic_dist,npermut = mantelTestPermutations,outdir = spectrum_DistancesOutdir,speciesToInclude = speciesToInclude_phylo[speciesToInclude_phylo %in% c("mice_Mmd","bears_ABC","vaquita","fin_whale_GOC","humans_AFR","apes_Pan_paniscus")],label="3merDistance.humanscaledtargets.OneSpeciesPerDataset"  )
############# 5mer distances ##############

trial_distanceDF_basedOnHumanTargetRescaling_5mer_MANTEL <- comboFunction_mantelTest_forSpectrumDistances_includesPlots(spectrum =all5merSpectra_MERGEDWITHUMANTARGETPROPORTIONS,variable="frac_seg_sites_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none="clr",phylo_dist =raxml_cophenetic_dist,npermut = mantelTestPermutations,outdir = spectrum_DistancesOutdir,speciesToInclude = speciesToInclude_phylo,label="5merDistance.humanscaledtargets.allSpecies"  )

############# 7mer distances ##############

trial_distanceDF_basedOnHumanTargetRescaling_7mer_MANTEL <- comboFunction_mantelTest_forSpectrumDistances_includesPlots(spectrum =all7merSpectra_MERGEDWITHUMANTARGETPROPORTIONS,variable="frac_seg_sites_RESCALED_BY_HUMAN_TARGETS",ilr_or_clr_or_none="clr",phylo_dist =raxml_cophenetic_dist,npermut = mantelTestPermutations,outdir = spectrum_DistancesOutdir,speciesToInclude = speciesToInclude_phylo,label="7merDistance.humanscaledtargets.allSpecies"  )
# wow is no longer significant! interesting! try with clock time?



###########>>>  try distances with collapsed types to avoid issues with polarization? <<< ########### YOU ARE HERE : 
# want to collapse 3mer types : 
# want to do as part of my processing framework. 
# after processing but before pivoting? 
# first let's just see what happens. 
spectrum_1mer_forSigfit_downsampledCounts_ExperimentingWithCollapsing <- spectrum_1mer_forSigfit_downsampledCounts
# set label
spectrum_1mer_forSigfit_downsampledCounts_ExperimentingWithCollapsing$label <- rownames(spectrum_1mer_forSigfit_downsampledCounts_ExperimentingWithCollapsing)
spectrum_1mer_forSigfit_downsampledCounts_ExperimentingWithCollapsing <- melt(spectrum_1mer_forSigfit_downsampledCounts_ExperimentingWithCollapsing)
#spectrum_1mer_forSigfit_downsampledCounts_ExperimentingWithCollapsing$collapsedLabel <- ######## YOU ARE HERE!!!!!!!!

################## log likelihood ##############
# sigfit as a get_loglik function in docs; why isn't it loading? 
# try it:
# MESSAGED GITHUB -- this isn't available yet! get_loglik(deNovoSigs_best) # doesn't work !
# want to collapse mis pol types to TTC --> TCC and TCC --> TTC get collapsed
# but some types then also need to be rev comped. 
# YOU ARE HERE: COLLAPSE possibly mis polarized mutation types and see if distances still good
# also could subjset to just A.T 
# try with 1mers (and 3mers?)

# how is this going to work target size? 
################## try doing distances based on collpsed spectra? stilla phylo signal? waht about if exclude ms? # 


####### run on dogs at some point 