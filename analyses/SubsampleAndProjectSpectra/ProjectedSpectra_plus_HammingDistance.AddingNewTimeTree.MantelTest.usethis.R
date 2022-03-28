############# gather all spectra and calculate rescaled mutation rates  ###########
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
########### not e: dplyr is a bit dangerous when dividing by sum(X) -- if the df has ever been grouped by things, it maintains that grouping. This manifests as sneaky bugs for fraction of segregating sites.
# Be very cautious ! !!! need to go back and check some resutls 
epsilon=1 # amount to add to 0 bins
separateCpGs="no" # yes or no
todaysdate=format(Sys.Date(),"%Y%m%d")
flagOfAnalysis=paste0(todaysdate,".plots.projectedCounts.epsilon.",epsilon,".sepCpGs.",separateCpGs,".addingApes_maskedMouseVaq.addingMultinomialDownsampling.AllAnalyses.plusSqrtPatristic") # a label that will be in title of plots and outdir 
plotdir=paste0("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/comparing_spectrum_and_phylo_divergence_WITHPROJECTIONS/",flagOfAnalysis,"/")
dir.create(plotdir,showWarnings = F)
mantelTestPermutationCount=99999
########## table of target locations ########
# restricting to just species with >1 individual:
#tableOfTargetsToInclude=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/spectra/tableOfSpectraToIncludeInComparisons.SpeciesWithMultipleIndividualsOnly.20211026.txt",header=T,sep="\t") # this also has info on spectra but don't use that. 
# adding apes:
tableOfTargetsToInclude=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/spectra/tableOfSpectraToIncludeInComparisons.SpeciesWithMultipleIndividualsOnly.PlusApes.MaskedMouse.Vaquita.20211123.txt",header=T,sep="\t")  # UPDATED this file; it now has bedmasked vaquita targets and repmasked mouse targets
# still waiting on vaquita targets (not done yet --- download when done!)
############## all projected spectra #############
# okay have all the spectra together already in:
# THESE have been projected down to k haploids:
#allProjectedSpectra <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/allKSFSes_and_projected_Spectra/All_SpectrumCounts_ProjectedDownTo.16.Haploids.txt",header=T)
allProjectedSpectra <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/allKSFSes_and_projected_Spectra/20211123_All_SpectrumCounts_ProjectedDownTo.10.Haploids.includesApes.RepMaskedMice.txt",header=T) # updated this file 20211123 to include repmasked mouse 
# this now includes bears_EUR and all apes (projectd down to 5 diploids, 10 haps)
speciesCodes=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/scripts/DNA_Shape_Project/analyses/PhylogeneticDistances/speciesCodes.txt",header=T) # for now just using AFR as human ; added apes
head(speciesCodes) # for getting phylo distances 

########### species to include ######################

speciesToInclude_in_raxmlTree=speciesList=c("Balaenoptera_physalus" ,"Phocoena_sinus" ,"Ursus_maritimus" ,"Ursus_arctos" ,"Mus_musculus" ,"Mus_spretus","Homo_sapiens","Pongo_pygmaeus","Pongo_abelii","Pan_troglodytes","Pan_paniscus","Gorilla_gorilla")  # raxml name fmt

speciesToInclude_phylo=c("mice_Mmd","mice_Ms","bears_ABC","bears_PB","vaquita","fin_whale_GOC","humans_AFR","apes_Gorilla","apes_Pan_paniscus","apes_Pan_troglodytes","apes_Pongo_abelii","apes_Pongo_pygmaeus") # human_AFR# switching to GOC fin whale for now due to ksfs weirdness
# adding in apes 
# 
#for hamming distances:
speciesToInclude_hamming = c("bears_ABC"  , "bears_PB","bears_EUR", "fin_whale_ENP" ,"fin_whale_GOC", "humans_AFR"   , "humans_AMR" ,"humans_EAS"  ,  "humans_EUR"  ,  "humans_SAS" ,   "mice_Mmc"      ,"mice_Mmd"   ,   "mice_Mmm" ,   "mice_Ms","apes_Gorilla","apes_Pan_paniscus","apes_Pan_troglodytes","apes_Pongo_abelii","apes_Pongo_pygmaeus")  # 20220214 -- adding apes to hamming

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

# restrict to just the species to include
raxmlTree_renamedTips_subset <- keep.tip(raxmlTree_renamedTips,speciesToInclude_in_raxmlTree)

raxml_cophenetic_dist <- cophenetic(raxmlTree_renamedTips_subset)
# save this as an object
saveRDS(raxml_cophenetic_dist,file=paste0(plotdir,"raxml_cophenetic_dist.dist"))

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

########## read in hamming distances:  ##########
hammingDistances = read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/hamming_distance/ALLSPECIESCOMBINED.totalAlleleCounts.HammingDistance.allIntervals.WITHPOPINFO.AVERAGEDPERCOMPARISONTYPE.DIVIDEDBY2xGenomeSIZE.INCLUDESAPES.usethis.txt",header=T)
### 20220214: added apes

##
######### read in time-calibrated tree cophenetic distances -- this is from UpHam et al. may want to use TimeTree Instead in future!  ########3
#timeDistances = read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/time_calibrated_trees/20211122_addingApes/All_Species/tip_dated/cophenetic.Pairwise.BranchLengthSums.AveragedOverTrees.avg.sd.txt",header=T) # figure out how I want to incorproate these. 
# am labeling the avg branch len as cophenetic distance (which it is) so that it works with existing code
#timeDistances$cophenetic_distance <- timeDistances$avgBranchLen # thisis avg cophenetic distance of tiem cal tree (across 1000 trees)
# switching to use timeTree isntead of Upham et al tree because I don't trust the Upham times for whales that much. 
# is this good? maybe it's wrong to not just use the Upham... well let's see.
timeTree=ape::read.tree("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/timeTreePhyloTree/listOfSpecies.ForTimeTree.20211210(1).nwk")
# get cophenetic distances:
timeTree_cophenetic_dist <- cophenetic(timeTree)
timeDistances = melt_cophenetic_func(timeTree_cophenetic_dist)

# now do same processing:
############ read in gneration times in days ##############
# from https://go.gale.com/ps/i.do?id=GALE%7CA541334609&sid=googleScholar&v=2.1&it=r&linkaccess=abs&issn=13146947&p=AONE&sw=w&userGroupName=anon%7Ec3d53501

reproductiveAges=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/mammal_generation_times/doi_10.5061_dryad.gd0m3__v1/GenerationLengthMammals.txt",header=T,sep="\t")
# relevant columsn are AFR_d (age at first reproduction (days)) and Rspan_d (reproductive lifespan) (days)
reproductiveAges$species <- gsub(" ","_",reproductiveAges$Scientific_name)
reproductiveAges[reproductiveAges$species %in% speciesCodes$species,]
##### 1mers are dealt with in the if statement below this (dealing with CpGs) ########
# get ancestral kmers:
allProjectedSpectra$ancestral7mer <- substr(allProjectedSpectra$variable,1,7)
allProjectedSpectra$ancestral5mer <- substr(allProjectedSpectra$variable,2,6)
allProjectedSpectra$ancestral3mer <- substr(allProjectedSpectra$variable,3,5)

# get central mutation types:
allProjectedSpectra$mutation_7mer <- allProjectedSpectra$variable
allProjectedSpectra$mutation_5mer <- paste0(substr(allProjectedSpectra$variable,2,6),".",substr(allProjectedSpectra$variable,10,14))
allProjectedSpectra$mutation_3mer <- paste0(substr(allProjectedSpectra$variable,3,5),".",substr(allProjectedSpectra$variable,11,13))

# want to get ancestral 1mers for all data frames so I can plot things with central mutattion type separated later:

#### SEPARATE OUT  -- OPTIONAL ############

if(separateCpGs=="yes"){
  #### label CpGs and specify them as ancestral target as well
  print("separating out CpG sites!")
  allProjectedSpectra$mutation_1mer_CpGNotLabeled <- paste0(substr(allProjectedSpectra$variable,4,4),".",substr(allProjectedSpectra$variable,12,12))
  allProjectedSpectra$ancestral1mer_CpGNotLabeled <- paste0(substr(allProjectedSpectra$variable,4,4)) # adding this on 20220112 
  #### label CpGs and specify them as ancestral target as well
  allProjectedSpectra$CpGLabel <- ""
  # find all CpG sites:
  allProjectedSpectra[grepl("CG",allProjectedSpectra$ancestral3mer) & allProjectedSpectra$ancestral1mer_CpGNotLabeled=="C",]$CpGLabel <- "_CpG" 
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

# okay do need to process the targets #
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
  # 20220112 -- note that grepping central CpG is fine because there are no CGX ancestral 3mers because they would be revcomped to YCG so you're good here!  
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
ggsave(paste0(plotdir,"targetPlot.1mers.png"),target1merPlot,height=9,width=5)
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



####### define some functions ###########
# okay you want to get spectrum and the mutation rate and want to add 1 to all counts to deal with 0 entries (assuming every mutation type would happen if you sampled enough -- though note this will have a bigger i)
totalSegSites <- all1merSpectra %>%
  group_by(species,label) %>%
  summarise(totalSegSites=sum(total_mutations_projected))

multinomlogfile <- file(paste0(plotdir,"multinomialLogFile.txt"))

writeLines(paste0("downsampling to ",min(totalSegSites$totalSegSites)," sites (",totalSegSites[totalSegSites$totalSegSites==min(totalSegSites$totalSegSites),]$species,")"),paste0(plotdir,"multinomialDownsampling.info.txt"))
close(multinomlogfile)

# add multinomial variables to this function.
processSpectra <- function(spectradf){
  # get mutation rates and fraction of segregating sites: also adding a column in where I add +epsilon (1) to every count to deal with 0 entries ; taking out adding epsilon to target counts since I fixed the empty target bug
  spectradf <- spectradf %>%
    mutate(total_mutations_projected_plusEpsilon = (total_mutations_projected+epsilon),mutation_rate_projected = (total_mutations_projected/total_target_count), mutation_rate_projected_plusEpsilon=(total_mutations_projected_plusEpsilon/total_target_count)) %>%
    group_by(species,population,label) %>%
    mutate(mutation_rate_projected_normalized = (mutation_rate_projected/sum(mutation_rate_projected)),fractionOfSegregatingSites_projected=(total_mutations_projected/sum(total_mutations_projected)), mutation_rate_projected_normalized_plusEpsilon=(mutation_rate_projected_plusEpsilon/sum(mutation_rate_projected_plusEpsilon)),fractionOfSegregatingSites_projected_plusEpsilon=(total_mutations_projected_plusEpsilon/sum(total_mutations_projected_plusEpsilon))) # adding plus1 versions.
  
  # add multinomial sampling:
  # calculate min number of seg sites:
  totalSegSites <- spectradf %>%
    group_by(species,population,label) %>%
    summarise(totalSegSites=sum(total_mutations_projected))
  
  subsampleValue=min(totalSegSites$totalSegSites)
  
  subsampleSpecies=totalSegSites[totalSegSites$totalSegSites==min(totalSegSites$totalSegSites),]$species
  
  print(paste0("downsampling to ",subsampleValue," sites (",subsampleSpecies,")"))
  # carry out multinomial sampling:
  
  spectradf <- spectradf %>%
    group_by(species,label,population) %>%
    mutate(total_mutations_projected_multinom_downsampled=as.numeric(rmultinom(n=1,size=subsampleValue,prob=fractionOfSegregatingSites_projected))) %>%  # n=1 is how many times to sample (reps); size is how many sites (~230,000), probs is vector of probabilities aka the raw fraction (no +epsilon, no target correction)
    ungroup() %>% # ungroup just in case here 
    # now am adding epsilon and doing all transforms to the downsampled counts 
    mutate(total_mutations_projected_multinom_downsampled_plusEpsilon = (total_mutations_projected_multinom_downsampled+epsilon),mutation_rate_projected_multinom_downsampled = (total_mutations_projected_multinom_downsampled/total_target_count), mutation_rate_projected_multinom_downsampled_plusEpsilon=(total_mutations_projected_multinom_downsampled_plusEpsilon/total_target_count)) %>%
    group_by(species,population,label) %>%
    mutate(mutation_rate_projected_multinom_downsampled_normalized = (mutation_rate_projected_multinom_downsampled/sum(mutation_rate_projected_multinom_downsampled)), mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon=(mutation_rate_projected_multinom_downsampled_plusEpsilon/sum(mutation_rate_projected_multinom_downsampled_plusEpsilon)),fractionOfSegregatingSites_projected_multinom_downsampled_plusEpsilon=(total_mutations_projected_multinom_downsampled/sum(total_mutations_projected_multinom_downsampled))) # adding plus1 versions.
  # get mutation rates and fraction of segregating sites: also adding a column in where I add +epsilon (1) to every count to deal with 0 entries ; taking out adding epsilon to target counts since I fixed the empty target bug
  # adding write out.
  return(spectradf)

}
# test function:
test = processSpectra(all1merSpectra)
head(test)
# check that sums up add up to 1 (if not would indicate a grouping error in dplyr)
test %>% group_by(species,label,population) %>%
  summarise(sum=sum(fractionOfSegregatingSites_projected_plusEpsilon),checkSum = sum(fractionOfSegregatingSites_projected_plusEpsilon)==1) # floating pt issues don't quite equal one but very close. 

# need to pivot it wider? 
pivotSpectra_perVariable <- function(processedspectradf,variableToPivot) {
  processedspectradf_wider <-  pivot_wider(processedspectradf,id_cols=c(mutation_label),names_from=c("label"),values_from =all_of(variableToPivot ))
  processedspectradf_wider$variable <- variableToPivot
  # pivots on a variable 
  return(processedspectradf_wider)
}

# this is just for mantel test:
pivotSpectra_perVariable_KEEPFULLSPECIESNAMES <- function(processedspectradf,variableToPivot) {
  processedspectradf_wider <-  pivot_wider(processedspectradf,id_cols=c(mutation_label),names_from=c("species.sppCodes"),values_from =all_of(variableToPivot ))
  processedspectradf_wider$variable <- variableToPivot
  # pivots on a variable 
  return(processedspectradf_wider)
}

# test function
testpivoted = pivotSpectra_perVariable(test,"mutation_rate_projected_normalized_plusEpsilon")
head(testpivoted)
# want to speed up with only a subset of species.
# ilr/clr calcs require no 0 entries
# can also specify no transform 
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
# note for manual clr to check things:
# exp(mean(log(testpivoted$bears_EUR))) gets geometric mean; then  you divide by geo meean
# and then you take log()
# it works! (checked with bear_EUR)
# some values will be negative and some will be pos (neg if less than geo mean; pos if gt than geo mean) -- then distances all end up positive once you take euclidean distance.

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

addHammingDistances <- function(spectrumdistancesdf,hammingdistancesdf) {
  
  # add alphabetical comparison label:
  # switch to characters not factors:
  spectrumdistancesdf$item1 <- as.character(spectrumdistancesdf$item1)
  spectrumdistancesdf$item2 <- as.character(spectrumdistancesdf$item2)
  
  spectrumdistancesdf$comparisonLabel <- paste0(pmin(spectrumdistancesdf$item1,spectrumdistancesdf$item2),".",pmax(spectrumdistancesdf$item1,spectrumdistancesdf$item2))
  # note hamming distances already has a dumb comp label that is just ABC.ABC = not useful. need ot use labels
  hammingdistancesdf$comparisonLabel <- paste0(pmin(hammingdistancesdf$pop1,hammingdistancesdf$pop2),".",pmax(hammingdistancesdf$pop1,hammingdistancesdf$pop2))
  
  merge1 <- merge(spectrumdistancesdf,hammingdistancesdf,by="comparisonLabel")
  return(merge1)
}



###### Get ILR and CLR distances ########
# this is slow so restrict just to species you want
##### this subsets just to the species you want above ^^ 
############# clr mutation rate ##############

####### combine functions ##########
# input: spectrum, variable (frac seg sites or mut rate), label (1mer , 5mer etc, and choice of ilr or clr)
combo_function_phylo <- function(spectrumdf,variable,ilr_or_clr_or_none,speciesToInclude,speciesCodes,phyloDistances){
  distance_dataframe <- 
    processSpectra(spectrumdf) %>%
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

######### hamming distance funciton ###########
# input: spectrum, variable (frac seg sites or mut rate), label (1mer , 5mer etc, and choice of ilr or clr)
combo_function_hamming <- function(spectrumdf,variable,ilr_or_clr_or_none,speciesToInclude,hammingDistances){
  distance_dataframe <- 
    processSpectra(spectrumdf) %>%
    pivotSpectra_perVariable(.,variable)  %>%
    select(c(variable,mutation_label,speciesToInclude)) %>% # just get the cols for the species you want 
    clr_or_ilr_orNoTransform_calculations(.,ilr_or_clr_or_none) %>%
    euclideanDistance(.) %>%
    addHammingDistances(.,hammingDistances) %>%
    mutate(variable=variable,transform_label=ilr_or_clr_or_none)
  
  
  return(distance_dataframe)
  
  
}


# formalize here and make some plots! 
# sweet this works!! and is very similar to original values. woooo.
############## run combo functions ######################
# do some loops
# all Distances
listOfDFs = list(spectra_1mer=all1merSpectra,spectra_3mer=all3merSpectra,spectra_5mer=all5merSpectra,spectra_7mer = all7merSpectra)
# lapply applies combo_function_phylo to each of the above spectra ^^ 
# if you set the names of a list, they get preserved in the output of lapply! cool!

listOfVars = c("mutation_rate_projected_normalized_plusEpsilon","fractionOfSegregatingSites_projected_plusEpsilon","mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon","fractionOfSegregatingSites_projected_multinom_downsampled_plusEpsilon")
listOfTransforms=c("none","clr") # skipping ilr for now 

########### run combo function on ALL data frames with all conditions; do for phylo and hamming ############
all_distances_all_conditions_phylo = data.frame()
all_distances_all_conditions_hamming = data.frame()
all_distances_all_conditions_time = data.frame()
for(variable in listOfVars){
  for(transform in listOfTransforms){
    # get for all dfs: 
    ########## phylo distances ###############
    list_of_distance_dfs_phylo <- lapply(listOfDFs,combo_function_phylo, variable =variable, ilr_or_clr_or_none = transform,speciesToInclude = speciesToInclude_phylo,speciesCodes = speciesCodes,phyloDistances = phyloDistances)
    # this will return a list of named results. How to deal with that? 
    # can bind them together and add some ids:
    distance_dfs_phylo <- bind_rows(list_of_distance_dfs_phylo, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
    all_distances_all_conditions_phylo <- bind_rows(all_distances_all_conditions_phylo,distance_dfs_phylo)
    rm(distance_dfs_phylo,list_of_distance_dfs_phylo) # to keep memory open
    
    ############# hamming distances ###############
    list_of_distance_dfs_hamming <- lapply(listOfDFs,combo_function_hamming, variable =variable, ilr_or_clr_or_none = transform,speciesToInclude = speciesToInclude_hamming,hammingDistances=hammingDistances)
    # this will return a list of named results. 
    # can bind them together and add some ids:
    distance_dfs_hamming <- bind_rows(list_of_distance_dfs_hamming, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
    all_distances_all_conditions_hamming <- bind_rows(all_distances_all_conditions_hamming,distance_dfs_hamming)
    rm(distance_dfs_hamming,list_of_distance_dfs_hamming) # to keep memory open
    
    
    ########### time distances (from time calibrated phylogeny )#########
    list_of_distance_dfs_time <- lapply(listOfDFs,combo_function_phylo, variable =variable, ilr_or_clr_or_none = transform,speciesToInclude = speciesToInclude_phylo,speciesCodes = speciesCodes,phyloDistances = timeDistances) # just put in time distances instead of raw phylo distances 
    # this will return a list of named results. How to deal with that? 
    # can bind them together and add some ids:
    distance_dfs_time <- bind_rows(list_of_distance_dfs_time, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
    all_distances_all_conditions_time <- bind_rows(all_distances_all_conditions_time,distance_dfs_time)
    rm(distance_dfs_time,list_of_distance_dfs_time) # to keep memory open
    
  }
}

########### an oddity here: after clr, frac seg sites and mut rate have EXACT same spectrum distance. don't expect that to be the case, right? what's going on with that? okay so I think this is just for species that are maped to same genome -- so somehow targets are canceling out during transform. ###########
write.table(all_distances_all_conditions_phylo,paste0(plotdir,"all_distances_all_conditions_phylo.txt"),row.names=F,quote=F,sep="\t")

write.table(all_distances_all_conditions_hamming,paste0(plotdir,"all_distances_all_conditions_hamming.txt"),row.names=F,quote=F,sep="\t")


write.table(all_distances_all_conditions_time,paste0(plotdir,"all_distances_all_conditions_timeCalibrated.txt"),row.names=F,quote=F,sep="\t")

########### want to write out processed spectra ##############
# need to keep in labels 
combo_function_JustGetTransformedSpectra_notDistances <- function(spectrumdf,variable,clr_or_none,speciesToInclude){
  procSpec <- processSpectra(spectrumdf) %>%
    pivotSpectra_perVariable(.,variable)  %>%
    select(c(variable,mutation_label,speciesToInclude)) # just get the cols for the species you want

  # do clr transform
  if(clr_or_none=="clr"){
    procSpec_transform <- data.frame(sapply(select(procSpec,-c(variable,mutation_label)),clr)) # this way you don't have to transpose and you get a much nicer shaped output
    procSpec_transform$mutation_label <- procSpec$mutation_label
    procSpec <- procSpec_transform 
  } else {
    
  }
  # want mutation labels to stay in here 
  return(procSpec)
}

all_processed_spectra_to_write_out = data.frame()

for(variable in listOfVars){
  for(transform in c("clr","none")){
    list_of_processed_spectra_notdistances <- lapply(listOfDFs,combo_function_JustGetTransformedSpectra_notDistances, variable =variable, clr_or_none = transform,speciesToInclude=unique(all1merSpectra$label)) # unique(all1merSpectra$label) is just to get all the species names. isn't 1mer specific 
    
    all_spectra <- bind_rows(list_of_processed_spectra_notdistances, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
    all_spectra$variable <- variable
    all_spectra$transform_id <- transform 
    all_processed_spectra_to_write_out <- bind_rows(all_processed_spectra_to_write_out,all_spectra)
    rm(all_spectra,list_of_processed_spectra_notdistances) # to keep memory open
    
  }}

# then writet his out.
write.table(all_processed_spectra_to_write_out,paste0(plotdir,'AllProcessedSpectra.AllConditions.AllVariables.AllSpecies.txt'),quote=F,sep="\t",row.names = F)
################ phylo plots #################
##### be careful that you aren't falsely combining things 

phylo_plot1 <- ggplot(all_distances_all_conditions_phylo,aes(x=cophenetic_distance,y=spectrum_distance,color=variable))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=1.8)+
  ggtitle(paste0("phylo distances\n",flagOfAnalysis))

phylo_plot1
ggsave(paste0(plotdir,"phylo_plot1_transforms.allconditions.png"),height=12,width=24,phylo_plot1)

# plot each individual statistic:
my.formula=y~x
### plot individual stats:
for(variable in listOfVars){
  phylo_plot1b <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variable,],aes(x=cophenetic_distance,y=spectrum_distance))+
    facet_wrap(~transform_label~id,scales="free",ncol=4)+
    geom_point()+
    theme_bw()+
    geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=.8)+
    ggtitle(paste0(variable,"\n",flagOfAnalysis))+
    geom_smooth(method="lm",se=F)+
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 label.x.npc = "left", label.y.npc = 0.95,
                 formula = my.formula, parse = TRUE, size = 4,color='red')
  ggsave(paste0(plotdir,"phylo_plot1b_transforms.",variable,".png"),height=8,width=16,phylo_plot1b)
}
########### phylo plot: apes only ########
phylo_plot2_apes <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$comparisonLabel_broad_alphabetical=="apes.apes",],aes(x=cophenetic_distance,y=spectrum_distance))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=1)+
  ggtitle(paste0(variable,"\n",flagOfAnalysis))
phylo_plot2_apes
ggsave(paste0(plotdir,"phylo_plot2_apes_transforms.",variable,".png"),height=8,width=16,phylo_plot2_apes)
####### plot both frac seg sites and mut rate with and without downsamp ######
pairsOfVars=list(mutrate=c("mutation_rate_projected_normalized_plusEpsilon","mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon"),fracSegSites=c("fractionOfSegregatingSites_projected_plusEpsilon","fractionOfSegregatingSites_projected_multinom_downsampled_plusEpsilon"))
for(index in c("mutrate","fracSegSites")){
  variablepair=pairsOfVars[index]
  phylo_plot1c <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable %in% unlist(variablepair),],aes(x=cophenetic_distance,y=spectrum_distance,color=variable))+
    facet_wrap(~transform_label~id,scales="free",ncol=4)+
    geom_point()+
    theme_bw()+
    geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=.8)+
    ggtitle(paste0(index,"\n",flagOfAnalysis))+
    geom_smooth(method="lm",se=F)
  ggsave(paste0(plotdir,"phylo_plot1c_transforms.comparingDownsampling.",index,".png"),height=10,width=24,phylo_plot1c)
}

############## hamming plots ###############
# 20220214 adding apes
hamming_plot1 <- ggplot(all_distances_all_conditions_hamming,aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,color=variable))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  geom_text_repel(aes(label=comparisonLabel),size=1.8)+
  ggtitle(paste0("hamming distances\n",flagOfAnalysis))

hamming_plot1
ggsave(paste0(plotdir,"hamming_plot1_transforms.allconditions.png"),height=12,width=24,hamming_plot1)
# note that in hamming plot, looks like when transformed differences betwn frac seg sites and mut rate disappear for genomes with same targets (eg PB adn ABC or human1 and human2 -- the ones in the hamming plot) ; must be some algebra that cancels out targets in some way. does'nt happen in phylo plot bc are mapped to diff genomes

# log 10 scale hamming plot:
hamming_plot1b <- hamming_plot1 + scale_x_log10() + xlab("hamming distance (log10 scaled)")
hamming_plot1b
ggsave(paste0(plotdir,"hamming_plot1b_log10_transforms.allconditions.png"),height=12,width=24,hamming_plot1b)

hamming_plot1c <-  ggplot(all_distances_all_conditions_hamming[all_distances_all_conditions_hamming$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon",],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  geom_text(aes(label=comparisonLabel),size=1)+
  ggtitle(paste0("hamming distances\n",flagOfAnalysis)) #+
  #scale_x_log10()

hamming_plot1c
ggsave(paste0(plotdir,"hamming_plot1_transforms.subsetConditions.png"),height=12,width=24,hamming_plot1c)

######## hamming: apes only ############
# just mutation rate variable with (or wi/out ?) downsampling
# pan paniscus distances all seem too small. 

hamming_plot2 <- ggplot(all_distances_all_conditions_hamming[all_distances_all_conditions_hamming$transform_label=="clr" & all_distances_all_conditions_hamming$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_hamming$label=="apes",],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  geom_text(aes(label=comparisonLabel),size=1.8)+
  ggtitle(paste0("hamming distances\n",flagOfAnalysis,"\nmutation_rate_projected_multinom_downsampled_normalized_plusEpsilon"))

hamming_plot2
ggsave(paste0(plotdir,"hamming_plot2_APESONLY.subsetConditions.png"),height=5,width=16,hamming_plot2)

hamming_plot2b_log10 <- hamming_plot2+xlab("hamming distance -log10 scaled axis")
ggsave(paste0(plotdir,"hamming_plot2b_APESONLY.subsetConditions.LOG10.png"),height=5,width=16,hamming_plot2b_log10)
# just mutation rate variable:
######## hamming: humans only #########
hamming_plot3 <- ggplot(all_distances_all_conditions_hamming[ all_distances_all_conditions_hamming$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_hamming$label=="humans" & all_distances_all_conditions_hamming$transform_label=="clr",],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  geom_text(aes(label=comparisonLabel),size=1.8)+
  ggtitle(paste0("hamming distances\n",flagOfAnalysis,"\nmutation_rate_projected_multinom_downsampled_normalized_plusEpsilon"))

hamming_plot3

ggsave(paste0(plotdir,"hamming_plot3_HUMANSONLY.subsetConditions.png"),height=5,width=16,hamming_plot3)

hamming_plot3b_log10 <- hamming_plot3+scale_x_log10()+xlab("hamming distance -log10 scaled axi")

ggsave(paste0(plotdir,"hamming_plot3b_HUMANSONLY.subsetConditions.LOG10.png"),height=5,width=16,hamming_plot3b_log10)

######### hamming: humans, non downsampled ###########
hamming_plot3c <- ggplot(all_distances_all_conditions_hamming[ all_distances_all_conditions_hamming$variable=="mutation_rate_projected_normalized_plusEpsilon" & all_distances_all_conditions_hamming$label=="humans" & all_distances_all_conditions_hamming$transform_label=="clr",],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  geom_text(aes(label=comparisonLabel),size=1.8)+
  ggtitle(paste0("hamming distances\n",flagOfAnalysis,"\nmutation_rate_projected_normalized_plusEpsilon"))

hamming_plot3c

ggsave(paste0(plotdir,"hamming_plot3c_HUMANSONLY.subsetConditions.NOTDOWNSAMPLED.png"),height=5,width=16,hamming_plot3c)
######### hamming: apes + humans only #########
hamming_plot4 <- ggplot(all_distances_all_conditions_hamming[ all_distances_all_conditions_hamming$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_hamming$label %in% c("apes","humans") & all_distances_all_conditions_hamming$transform_label=="clr",],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,color=label))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  geom_text(aes(label=comparisonLabel),size=1.8)+
  ggtitle(paste0("hamming distances\n",flagOfAnalysis,"\nmutation_rate_projected_multinom_downsampled_normalized_plusEpsilon"))

hamming_plot4

ggsave(paste0(plotdir,"hamming_plot4_APES_PLUS_HUMANSONLY.subsetConditions.png"),height=5,width=16,hamming_plot4)

hamming_plot4b_log10 <- hamming_plot4+scale_x_log10()+xlab("hamming distance -- log10 scaled axis")
hamming_plot4b_log10

ggsave(paste0(plotdir,"hamming_plot4b_APES_PLUS_HUMANSONLY.subsetConditions.LOG10.png"),height=5,width=16,hamming_plot4b_log10)


############### hamming: all species ###########
hamming_plot5 <- ggplot(all_distances_all_conditions_hamming[ all_distances_all_conditions_hamming$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_hamming$transform_label=="clr",],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,color=label))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  geom_text(aes(label=comparisonLabel),size=1.8)+
  ggtitle(paste0("hamming distances\n",flagOfAnalysis,"\nmutation_rate_projected_multinom_downsampled_normalized_plusEpsilon"))

hamming_plot5

ggsave(paste0(plotdir,"hamming_plot5_ALLSPECIES.subsetConditions.png"),height=5,width=16,hamming_plot5)

hamming_plot5b_log10 <- hamming_plot5+scale_x_log10()+xlab("hamming distance -- log10 scaled axis")
hamming_plot5b_log10

ggsave(paste0(plotdir,"hamming_plot5b_ALLSPECIES.subsetConditions.LOG10.png"),height=5,width=16,hamming_plot5b_log10)

############ hamming + phylo plots #############
hamming_phylo_combo_plot1 <- ggplot(all_distances_all_conditions_hamming,aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,color=variable,shape="hamming"))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  geom_point(data=all_distances_all_conditions_phylo,aes(x=cophenetic_distance,y=spectrum_distance,color=variable,shape="phylo"))+
  theme_bw()+
  #geom_text_repel(aes(label=comparisonLabel),size=1.8)+
  ggtitle(paste0("hamming distances and phylo distances \n",flagOfAnalysis))+
  xlab("hamming or phylo distance")
hamming_phylo_combo_plot1
ggsave(paste0(plotdir,"hamming_phylo_combo_plot1_transforms.allconditions.png"),height=6,width=16,hamming_phylo_combo_plot1)

# log10 scale it:
hamming_phylo_combo_plot1b <- hamming_phylo_combo_plot1+
  scale_x_log10()+
  xlab("hamming or phylo distance (log10 scaled)")
hamming_phylo_combo_plot1b
ggsave(paste0(plotdir,"hamming_phylo_combo_plot1b_log10scaled_transforms.allconditions.png"),height=6,width=16,hamming_phylo_combo_plot1b)


# just plot mutaiton rate: 
# plot each variable separately:
for(variable in listOfVars){
hamming_phylo_combo_plot1c <- ggplot(all_distances_all_conditions_hamming[all_distances_all_conditions_hamming$variable==variable,],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,shape="hamming"))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  geom_point(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable==variable,],aes(x=cophenetic_distance,y=spectrum_distance,shape="phylo"))+
  theme_bw()+
  #geom_text_repel(aes(label=comparisonLabel),size=1.8)+
  ggtitle(paste0("hamming distances and phylo distances\n",variable,"\n",flagOfAnalysis))+
  xlab("hamming or phylo distance")
hamming_phylo_combo_plot1c
ggsave(paste0(plotdir,"hamming_phylo_combo_plot1_transforms.",variable,".png"),height=12,width=24,hamming_phylo_combo_plot1c)

}



############ humans and mice only ##########
# humansMice=c("humans_AFR","humans_EUR","humans_EAS","humans_SAS","humans_AMR","mice_Mmd","mice_Ms","mouse_Ms","mouse_Mmd","mice_Mmm","mice_Mmc")
# 
# hamming_phylo_combo_plot2 <- ggplot(all_distances_all_conditions_hamming[all_distances_all_conditions_hamming$item1 %in% humansMice & all_distances_all_conditions_hamming$item2 %in% humansMice,],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,color=variable,shape="hamming"))+
#   facet_wrap(~transform_label~id,scales="free",ncol=4)+
#   geom_point()+
#   geom_text_repel(aes(label=comparisonLabel))+
#   geom_point(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$item1 %in% humansMice & all_distances_all_conditions_phylo$item2 %in% humansMice,],aes(x=cophenetic_distance,y=spectrum_distance,color=variable,shape="phylo"))+
#   theme_bw()+
#   #geom_text_repel(aes(label=comparisonLabel),size=1.8)+
#   ggtitle(paste0("hamming distances and phylo distances \n",flagOfAnalysis))+
#   xlab("hamming or phylo distance")+
#   geom_text_repel(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$item1 %in% humansMice & all_distances_all_conditions_phylo$item2 %in% humansMice,],aes(x=cophenetic_distance,y=spectrum_distance,label=comparisonLabel_common_alphabetical))
# hamming_phylo_combo_plot2
# ggsave(paste0(plotdir,"hamming_phylo_combo_plot2_humans_mice_only_transforms.allconditions.png"),height=6,width=16,hamming_phylo_combo_plot2)
# 
# ########## humans mice apes only  ##############
# humansMiceApes=c("humans_AFR","humans_EUR","humans_EAS","humans_SAS","humans_AMR","mice_Mmd","mice_Ms","mouse_Ms","mouse_Mmd","mice_Mmm","mice_Mmc","apes_Pan_paniscus","apes_Pan_troglodytes","apes_Pongo_abelii","apes_Gorilla","apes_Pongo_pygmaeus")
# 
# hamming_phylo_combo_plot3 <- ggplot(all_distances_all_conditions_hamming[all_distances_all_conditions_hamming$item1 %in% humansMiceApes & all_distances_all_conditions_hamming$item2 %in% humansMiceApes,],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,color=variable,shape="hamming"))+
#   facet_wrap(~transform_label~id,scales="free",ncol=4)+
#   geom_point()+
#   geom_text_repel(aes(label=comparisonLabel))+
#   geom_point(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$item1 %in% humansMiceApes & all_distances_all_conditions_phylo$item2 %in% humansMiceApes,],aes(x=cophenetic_distance,y=spectrum_distance,color=variable,shape="phylo"))+
#   theme_bw()+
#   #geom_text_repel(aes(label=comparisonLabel),size=1.8)+
#   ggtitle(paste0("hamming distances and phylo distances \n",flagOfAnalysis))+
#   xlab("hamming or phylo distance")+
#   geom_text_repel(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$item1 %in% humansMiceApes & all_distances_all_conditions_phylo$item2 %in% humansMiceApes,],aes(x=cophenetic_distance,y=spectrum_distance,label=comparisonLabel_common_alphabetical),size=0.8)
# hamming_phylo_combo_plot3
# ggsave(paste0(plotdir,"hamming_phylo_combo_plot3_humans_mice_apes_only_transforms.allconditions.png"),height=6,width=16,hamming_phylo_combo_plot3)



#############  sandbox: restrict to semi-independent comparisons ###################

# want 1 each: human-human, ape-ape, human-ape, mouse-mouse, human-mouse, whale-mouse, bear-bear, whale-whale ; don't want multiple: 
# can compare within a clade but don't do human to both chimps etc.
# both chimps, both orangs, human pops,homo:1chimp; gorilla:homo; bear:bear, whale:whale; 1bear:1whale; mouse;mouse; mouse: one ape; mouse: one whale or bear. 
#comparisonsIWant=c("Pongo_abelii.Pongo_pygmaeus","Pan_paniscus.Pan_troglodytes" , "Mus_musculus.Mus_spretus" ,"Ursus_arctos.Ursus_maritimus" ,"Homo_sapiens.Pan_paniscus" ,"Gorilla_gorilla.Homo_sapiens" ,"Phocoena_sinus.Ursus_arctos" , "Homo_sapiens.Mus_musculus" )
#comparisonsIWant_some_indepdendence1_phylo = c("mouse_Mmd.mouse_Ms","brown_bear_ABC.polar_bear" ,"bornean_orangutan.sumatran_orangutan","bonobo.chimpanzee" ,"fin_whale_GOC.vaquita","chimpanzee.humans_AFR" ,"gorilla.humans_AFR" ,"brown_bear_ABC.vaquita" , "humans_AFR.mouse_Mmd","humans_AFR.vaquita" ,    "mouse_Mmd.vaquita") # not sure if it's fair to have mouse:vaquita and human:vaquita; hm. note AFR-EUR comparison is hamming distance only.
# comparisonsIWant_some_indepdendence1_phylo = c("mouse_Mmd.mouse_Ms","brown_bear_ABC.polar_bear" ,"bornean_orangutan.sumatran_orangutan","bonobo.chimpanzee" ,"fin_whale_GOC.vaquita","chimpanzee.humans_AFR" ,"gorilla.humans_AFR" ,"brown_bear_ABC.fin_whale_GOC" , "humans_AFR.mouse_Mmd","fin_whale_GOC.humans_AFR" ,    "mouse_Mmd.fin_whale_GOC")
# comparisonsIWant_some_indepdendence1_hamming = c("bears_ABC.bears_EUR","fin_whale_ENP.fin_whale_GOC","humans_AFR.humans_EUR","mice_Mmc.mice_Mmd")
# # plot all points in light gray, then add in points of interest as color
# all_distances_all_conditions_hamming$ptOfInterest <- "grayout"
# all_distances_all_conditions_hamming[all_distances_all_conditions_hamming$comparisonLabel %in% comparisonsIWant_some_indepdendence1_hamming,]$ptOfInterest <- "mark"
# all_distances_all_conditions_phylo$ptOfInterest <- "grayout"
# all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$comparisonLabel_common_alphabetical %in% comparisonsIWant_some_indepdendence1_phylo,]$ptOfInterest <- "mark"
# 
# hamming_phylo_combo_plot4 <- ggplot(all_distances_all_conditions_hamming[all_distances_all_conditions_hamming$comparisonLabel %in% comparisonsIWant_some_indepdendence1_hamming & all_distances_all_conditions_hamming$variable=="mutation_rate_projected_normalized_plusEpsilon",],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,shape="hamming"))+
#   facet_wrap(~transform_label~id,scales="free",ncol=4)+
#   geom_point()+
#   geom_text_repel(aes(label=comparisonLabel),size=1.8)+
#   geom_point(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$comparisonLabel_common_alphabetical %in% comparisonsIWant_some_indepdendence1_phylo & all_distances_all_conditions_phylo$variable=="mutation_rate_projected_normalized_plusEpsilon",],aes(x=cophenetic_distance,y=spectrum_distance,shape="phylo"))+
#   theme_bw()+
#   #geom_text_repel(aes(label=comparisonLabel),size=1.8)+
#   ggtitle(paste0("hamming distances and phylo distances \n",flagOfAnalysis))+
#   xlab("hamming or phylo distance")+
#   geom_text_repel(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$comparisonLabel_common_alphabetical %in% comparisonsIWant_some_indepdendence1_phylo & all_distances_all_conditions_phylo$variable=="mutation_rate_projected_normalized_plusEpsilon",],aes(x=cophenetic_distance,y=spectrum_distance,label=comparisonLabel_common_alphabetical),size=1.8)
# 
# hamming_phylo_combo_plot4
# 
# ggsave(paste0(plotdir,"hamming_phylo_combo_plot4_TRYINGFORSOMEINDEPENDENCE_only_transforms.allconditions.png"),height=6,width=16,hamming_phylo_combo_plot4)
# 
# ########## gray out other points ##########
# hamming_phylo_combo_plot4b <- ggplot(all_distances_all_conditions_hamming[all_distances_all_conditions_hamming$variable=="mutation_rate_projected_normalized_plusEpsilon",],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,shape="hamming",color=ptOfInterest))+
#   facet_wrap(~transform_label~id,scales="free",ncol=4)+
#   geom_point()+
#   #geom_text_repel(aes(label=comparisonLabel),size=1)+
#   geom_point(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_normalized_plusEpsilon",],aes(x=cophenetic_distance,y=spectrum_distance,shape="phylo",color=ptOfInterest))+
#   scale_color_manual(values=c("grey50","red"))+
#   theme_bw()+
#   #geom_text_repel(aes(label=comparisonLabel),size=1)+
#   ggtitle(paste0("hamming distances and phylo distances \n",flagOfAnalysis))+
#   xlab("hamming or phylo distance") #+
#   #geom_text_repel(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$comparisonLabel_common_alphabetical %in% comparisonsIWant_some_indepdendence1_phylo & all_distances_all_conditions_phylo$variable=="mutation_rate_projected_normalized_plusEpsilon",],aes(x=cophenetic_distance,y=spectrum_distance,label=comparisonLabel_common_alphabetical),size=1.8)
# 
# hamming_phylo_combo_plot4b
# 
# 
# ggsave(paste0(plotdir,"hamming_phylo_combo_plot4b_TRYINGFORSOMEINDEPENDENCE_only_transforms.allconditions.GRAYOUT.png"),height=6,width=16,hamming_phylo_combo_plot4b)
# ####### restricted plot:
############# time plots: time calibrated (should make mice an outlier ) ############
time_plot1 <- ggplot(all_distances_all_conditions_time,aes(x=cophenetic_distance,y=spectrum_distance,color=variable))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  xlab("cophenetic distance (time-calibrated tree)")+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=0.8)+
  ggtitle(paste0("time distances\n",flagOfAnalysis))

time_plot1
ggsave(paste0(plotdir,"time_plot1_transforms.allconditions.png"),height=6,width=16,time_plot1)

# plot each individual statistic:
time_plot1b <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$variable=="mutation_rate_projected_normalized_plusEpsilon",],aes(x=cophenetic_distance,y=spectrum_distance))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  xlab("cophenetic distance (time-calibrated tree)")+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=.8,max.overlaps = 10)+
  ggtitle(paste0("time distances\nmutation_rate_projected_normalized_plusEpsilon\n",flagOfAnalysis))+
  geom_smooth(method="lm")

time_plot1b
ggsave(paste0(plotdir,"time_plot1b_transforms.mutationrate.png"),height=6,width=16,time_plot1b)

# color by comparison type:
time_plot1c <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$variable=="mutation_rate_projected_normalized_plusEpsilon",],aes(x=cophenetic_distance,y=spectrum_distance,color=comparisonLabel_broad_alphabetical))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  xlab("cophenetic distance (time-calibrated tree)")+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=.8)+
  ggtitle(paste0("time distances\nmutation_rate_projected_normalized_plusEpsilon\n",flagOfAnalysis))# +
  #geom_smooth(method="lm")

time_plot1c
ggsave(paste0(plotdir,"time_plot1c_transforms.mutationrate.ColoredByType.png"),height=6,width=16,time_plot1c)



################# PCA ######################
pca_function <- function(processeddf,kmersizelabel,variable,ilr_or_clr_or_none,loadings=T,plotdir){
  outdir=paste0(plotdir,"pcaPlots/")
  dir.create(outdir,showWarnings = F)
  # specify if you want loadings (don't want for 5mer or 7mer -- just turns into a cloud) 
  # this will do pca and plot it from a processed --NOT PIVOTED-- df
  # for pca need a different pivot -- need df to be 'wide' with mutation types along the top (columns)
  # and species as the rows. 
  # need to make it so that mutaiton types are columns and species are rows: 
  df_WIDE <-  spread(data.frame(processeddf[,c("label","species","population","mutation_label",variable)]),key = mutation_label,variable) # spread the variable you want. 
  

  ### apply transform if it is specified: how that I've pivoted teh df to have mutation types along columns
  # I need to transform ROWS not colums (different from above)
  if(ilr_or_clr_or_none=="clr"){
    transformed_df <- data.frame(clr(select(df_WIDE,-c("label","species","population")))) # clr works ROWWISE which here is appropriate since species are ROWs. above the species were columns so instead I had to use sapply.
    transformed_df$species <- df_WIDE$species
    transformed_df$population <- df_WIDE$population
    transformed_df$label <- df_WIDE$label
    
  } else if(ilr_or_clr_or_none=="ilr") {
    
    transformed_df <- data.frame(ilr(select(df_WIDE,-c("label","species","population")))) # clr works ROWWISE which here is appropriate since species are ROWs. above the species were columns so instead I had to use sapply.
    transformed_df$species <- df_WIDE$species
    transformed_df$population <- df_WIDE$population
    transformed_df$label <- df_WIDE$label
    # note with ilr you won't get mutation types as columns any more 
  } else if(ilr_or_clr_or_none=="none"){
    transformed_df = df_WIDE # don't transform it. 
  } else {
    print("invalid transform choice")
    break
  }
  
  pca <- prcomp(subset(transformed_df, select=-c(label,species,population)),scale=T,center=T) # Michael says to scale/center; similar with and without
  if(loadings==T){
  pcaPlot_withLoadings <- autoplot(pca,data=transformed_df,colour="species",size=4,label.colour="black",loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black")+theme_bw()+geom_text(aes(color=species,label=label))+ggtitle(paste0(kmersizelabel,"\n",flagOfAnalysis,"\n",variable,"\ntransform = ",ilr_or_clr_or_none))+theme(legend.position = "none")
  ggsave(paste0(outdir,"pca.",kmersizelabel,".",variable,".transform.",transform,".Plot_withLoadings.png"),pcaPlot_withLoadings,height=7,width=10)
  }
  # plot without loadings as well: 
  pcaPlot_noLoadings <- autoplot(pca,data=transformed_df,colour="species",size=4,label.colour="black")+theme_bw()+geom_text(aes(color=species,label=label))+ggtitle(paste0(kmersizelabel,"\n",flagOfAnalysis,"\n",variable,"\ntransform = ",ilr_or_clr_or_none))+theme(legend.position = "none")
  ggsave(paste0(outdir,"pca.",kmersizelabel,".",variable,".transform.",transform,".Plot_noLoadings.png"),pcaPlot_noLoadings,height=7,width=10)
  
  # don't need to return anything
  # return(pca)
  return(pcaPlot_noLoadings)
} 

combo_function_pca <- function(spectrumdf,kmsersizelabel,variable,ilr_or_clr_or_none,loadings=T,plotdir,speciesToExclude) {
  processSpectra(spectrumdf) %>%
  filter(!label %in% speciesToExclude) %>%
  pca_function(.,kmsersizelabel,variable,ilr_or_clr_or_none,loadings=loadings,plotdir)
    
}


########### apply pca function to all DFs with all transforms and variables ###########

#listOfDFs = list(spectra_1mer=all1merSpectra,spectra_3mer=all3merSpectra,spectra_5mer=all5merSpectra,spectra_7mer = all7merSpectra)
# lapply applies combo_function_phylo to each of the above spectra ^^ 
# if you set the names of a list, they get preserved in the output of lapply! cool!

#listOfVars = c("mutation_rate_projected_normalized_plusEpsilon","fractionOfSegregatingSites_projected_plusEpsilon") # use same list of vars
#listOfTransforms=c("none","clr") # skipping ilr for now 
speciesToExclude=c("fin_whale_GOC","fin_whale_ENP","mice_Mmc","mice_Mmm","humans_EAS","humans_SAS","humans_AMR") # excluding these to reduce sample sizes to make more comparable; excluding fin whale because of data filtering issues. 
# want to also downsample humans mice to two each 
for(variable in listOfVars){
  listOfPCAPlots = list()
  for(transform in listOfTransforms){
    # get for all dfs: 
    ########## pca ###############
    # hm okay need to deal with labeling each one 
    # maybe just run separately because 5mer/7mer loadings are crazy
    # not going to run on list of DFs -- too annoying. 
    # 1mer:
    pcaplot1mer <- combo_function_pca(all1merSpectra,"1mer_spectrum",variable,transform,loadings=T,plotdir,speciesToExclude)
    # 3mer:
    pcaplot3mer <- combo_function_pca(all3merSpectra,"3mer_spectrum",variable,transform,loadings=T,plotdir,speciesToExclude)
    # 5mer: (turn off loadings)
    pcaplot5mer <- combo_function_pca(all5merSpectra,"5mer_spectrum",variable,transform,loadings=F,plotdir,speciesToExclude)
    # 7mer: (turn off loadings)
    pcaplot7mer <- combo_function_pca(all7merSpectra,"7mer_spectrum",variable,transform,loadings=F,plotdir,speciesToExclude)
    listOfPCAPlots=append(listOfPCAPlots,list(pcaplot1mer,pcaplot3mer,pcaplot5mer,pcaplot7mer))

  }
  grid <- grid.arrange(grobs=listOfPCAPlots,ncol=4)
  ggsave(paste0(plotdir,"pcaPlots/gridOFPCAPlots.",variable,".png"),grid,height=10,width=20)
  
}


########### troubleshooting: why are bear and vaquita so similar ?? ####################
trblshootdir=paste0(plotdir,"troubleshooting/")
dir.create(trblshootdir,showWarnings = F)
listOfDFs_forTrblshooting=list(spectrum_1mer=all1merSpectra,spectrum_3mer=all3merSpectra)
for(i in seq(1,length(listOfDFs_forTrblshooting))){
  label=names(listOfDFs_forTrblshooting)[i]
  df=data.frame(listOfDFs_forTrblshooting[[i]])
  for(variable in listOfVars){
    tblshootdf <- processSpectra(df) %>%
      pivotSpectra_perVariable(.,variable)
    
    trblshootPlot1 <-  ggplot(tblshootdf,aes(x=vaquita,y=bears_ABC))+
      geom_point(aes(color="bears_ABC"))+
      geom_abline()+
      geom_point(data=tblshootdf,aes(x=vaquita,y=bears_PB,color="bears_PB"))+
      geom_point(data=tblshootdf,aes(x=vaquita,y=fin_whale_ENP,color="fin_whale_ENP"))+
      geom_point(data=tblshootdf,aes(x=vaquita,y=humans_AFR,color="humans_AFR"))+
      geom_point(data=tblshootdf,aes(x=vaquita,y=mice_Mmd,color="mice_Mmd"))+
      ylab("other species")+
      ggtitle(variable)
    ggsave(paste0(trblshootdir,label,".",variable,".trblshootPlot1.png"),trblshootPlot1,height=6,width=8)
    
  }
}

############## add in generation times ########################

#merge_WithGenerationTimes <- merge(all_distances_all_conditions_phylo,reproductiveAges,by.x)
# what do I want to plot -- reproductive timing?
# some have no info 

genTimePlot1 <- ggplot(reproductiveAges[reproductiveAges$species %in% speciesCodes$species & reproductiveAges$AFR_d!="no information",],aes(y=species,x=as.numeric(AFR_d)))+
  geom_point()+
  geom_point(aes(y="Phocoena_sinus",x=3*365,color="Annabel rough estimate"))+
  geom_point(aes(y="Homo_sapiens",x=21*365,color="Annabel rough estimate"))+
  ggtitle("estimates of age at first reproduction (days)\n from Pacifici et al")
genTimePlot1
ggsave(paste0(plotdir,"generationTimes.plot1.png"),height=8,width=8,genTimePlot1)


genTimePlot2 <- ggplot(reproductiveAges[reproductiveAges$species %in% speciesCodes$species & reproductiveAges$Rspan_d!="no information",],aes(y=species,x=as.numeric(Rspan_d)))+
  geom_point()+
  geom_point(aes(y="Phocoena_sinus",x=(25-3)*365,color="Annabel rough estimate"))+
  geom_point(aes(y="Homo_sapiens",x=(40-21)*365,color="Annabel rough estimate"))+
  ggtitle("estimates of reproductive lifespan (days)\n from Pacifici et al")
genTimePlot2
ggsave(paste0(plotdir,"generationTimes.plot2.reproductiveLifespan.png"),height=8,width=8,genTimePlot2)


############ try to make nj tree based on spectrum distances instead of genetic dsitance #######
# from ape, need distance matrix 
require(ape)
ape::nj()
# keep it as distance mat
euclideanDistance_keepMatrix <- function(pivotedspectradf_noIDVars){
  # make it a matrix with each row as a species (transposed and all numeric)
  # variable in question:
  # distances are row-wise nto column wise so need to transpose 
  distances <- dist(t(pivotedspectradf_noIDVars)) # transpose for dist otherwise it goes along rows and is wrong
  #colnames(distances) <- c("item1","item2","spectrum_distance")
  return(distances)
}


euclideanDistance_keepMatrix_UPPERANDLOWER <- function(pivotedspectradf_noIDVars){
  # make it a matrix with each row as a species (transposed and all numeric)
  # variable in question:
  # distances are row-wise nto column wise so need to transpose 
  distances <- dist(t(pivotedspectradf_noIDVars),diag=T,upper=T) # transpose for dist otherwise it goes along rows and is wrong
  #colnames(distances) <- c("item1","item2","spectrum_distance")
  return(distances)
}

combo_function_njTree <- function(spectrumdf,variable,ilr_or_clr_or_none,speciesToInclude){
  distance_dataframe <- 
    processSpectra(spectrumdf) %>%
    pivotSpectra_perVariable(.,variable)  %>%
    select(c(variable,mutation_label,speciesToInclude)) %>% # just get the cols for the species you want 
    clr_or_ilr_orNoTransform_calculations(.,ilr_or_clr_or_none) %>%
    euclideanDistance_keepMatrix(.)
  
  return(distance_dataframe)
  
  
}
distanceMatsFromSpectrumDistancesForMakingTrees_mutRate <- lapply(listOfDFs,combo_function_njTree,variable="mutation_rate_projected_normalized_plusEpsilon",ilr_or_clr_or_none="clr",speciesToInclude=c(speciesToInclude_hamming,speciesToInclude_phylo))
# I think something is wrong here. Why are Mmd and Ms so different in the distance matrix. 

distanceMatsFromSpectrumDistancesForMakingTrees_mutRate_multinomDownsamp <- lapply(listOfDFs,combo_function_njTree,variable="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon",ilr_or_clr_or_none="clr",speciesToInclude=c(speciesToInclude_hamming,speciesToInclude_phylo))


trees_based_on_spectra <- lapply(distanceMatsFromSpectrumDistancesForMakingTrees_mutRate,ape::nj)
trees_based_on_spectra_multinomDownsamp <- lapply(distanceMatsFromSpectrumDistancesForMakingTrees_mutRate_multinomDownsamp,ape::nj)

treesBasedOnSpectraPlotDir=paste0(plotdir,"/treesBasedOnSpectrumDistance/")
dir.create(treesBasedOnSpectraPlotDir,showWarnings = F)
# make tree plots:
for(spectrum in names(trees_based_on_spectra)){
  png(paste0(treesBasedOnSpectraPlotDir,"nj.tree.basedon.distances.clr.mutationrate.",spectrum,".png"),height=8,width=8,units="in",res=300)
  plot.phylo(trees_based_on_spectra[[spectrum]],type="unrooted")
  title(paste0("neighbor joining tree based on clr distance (norm mut rate)\n",spectrum))
  dev.off()
  
}

for(spectrum in names(trees_based_on_spectra_multinomDownsamp)){
  png(paste0(treesBasedOnSpectraPlotDir,"nj.tree.basedon.distances.clr.mutationrate.multinomDownsampled.",spectrum,".png"),height=8,width=8,units="in",res=300)
  plot.phylo(trees_based_on_spectra_multinomDownsamp[[spectrum]],type="unrooted",)
  title(paste0("neighbor joining tree based on clr distance (norm mut rate with multinom downsampling)\n",spectrum))
  dev.off()
  
}

# do this based on downsampled mutation rate as well ! 

####### TO DO: need to get diffs between rep times #######

###### plot some correlations between confounders ##########
# difference in callable sites?
# ratio of callable sites?
# het
# gen time
head(all1merSpectra)
totalTargetSizes_ballpark <- unique(all1merSpectra[,c("species","label","total_target_count","ancestral1mer")]) %>% # need to unique() here because otherwise it would sum over multiple instances of target sizes. working from this instead of all1merTargets so I have nice label
  group_by(species,label) %>%
  summarise(totalTargets=sum(total_target_count))
totalTargetSizes_ballpark

totalSegSites_ballpark <- all1merSpectra %>%
  group_by(species,label) %>%
  summarise(totalSegSites=sum(total_mutations_projected))
totalSegSites_ballpark

reproductiveAges_codes <- merge(reproductiveAges,speciesCodes,by="species")
# ugly code for now:
# okay want to merge all these in 
# dim shouldn't change
# add targets:
merge_withMeta1_phylo <- merge(all_distances_all_conditions_phylo,totalTargetSizes_ballpark,by.x="item1",by.y="label") # targets for item 1

merge_withMeta2_phylo <- merge(merge_withMeta1_phylo,totalTargetSizes_ballpark,by.x="item2",by.y="label",suffixes=c(".item1",".item2")) # targets for item 1

# add total seg sites:
merge_withMeta3_phylo <- merge(merge_withMeta2_phylo,totalSegSites_ballpark,by.x=c("item1","species.item1"),by.y=c("label","species"))

merge_withMeta4_phylo <- merge(merge_withMeta3_phylo,totalSegSites_ballpark,by.x=c("item2","species.item2"),by.y=c("label","species"),suffixes=c(".item1",".item2"))

# add repro ages:
# add in humans
reproductiveAges_codes_human_vaquita_ballpark <- data.frame(code=c("humans_AFR","vaquita"),AFR_d=c(as.character(22*365),as.character(3*365)),Max_longevity_d=c(as.character(90*365),as.character(25*365)),Rspan_d=c(as.character((40-22)*365),as.character((25-3)*365))) # ball park making up human and vaquita parameters here
reproductiveAges_codes_human_vaquita_ballpark
# there was a 'no info' vaquita row; get rid of that
reproductiveAges_codes_addInHumans_Vaquita <- bind_rows(reproductiveAges_codes[reproductiveAges_codes$code!="vaquita",c("Max_longevity_d","Rspan_d","AFR_d","code")],reproductiveAges_codes_human_vaquita_ballpark)

merge_withMeta5_phylo <- merge(merge_withMeta4_phylo,reproductiveAges_codes_addInHumans_Vaquita,by.x="item1",by.y="code") # 


merge_withMeta6_phylo <- merge(merge_withMeta5_phylo,reproductiveAges_codes_addInHumans_Vaquita,by.x="item2",by.y="code",suffixes = c(".item1",".item2")) # 

####### okay so this is all merged ######
merge_withMeta_final_phylo <- merge_withMeta6_phylo

####### get diffs 
# age at first reproduction:
merge_withMeta_final_phylo$absDiff.AFR_d.Item1.Item2 <- abs(as.numeric(merge_withMeta_final_phylo$AFR_d.item1) - as.numeric(merge_withMeta_final_phylo$AFR_d.item2)) # some have no info. could add in ms and others as needed

merge_withMeta_final_phylo$absDiff.targetSize <- abs(merge_withMeta_final_phylo$totalTargets.item1 - merge_withMeta_final_phylo$totalTargets.item2)


merge_withMeta_final_phylo$absDiff.totalSegSites <- abs(merge_withMeta_final_phylo$totalSegSites.item1 - merge_withMeta_final_phylo$totalSegSites.item2)

# also get seg sites / target size (don't need to use harm mean because already samp size corrected?)
merge_withMeta_final_phylo$segSites.divBy.Targets.item1 <- merge_withMeta_final_phylo$totalSegSites.item1/merge_withMeta_final_phylo$totalTargets.item1

merge_withMeta_final_phylo$segSites.divBy.Targets.item2 <- merge_withMeta_final_phylo$totalSegSites.item2/merge_withMeta_final_phylo$totalTargets.item2

merge_withMeta_final_phylo$absDiff.segSites.divBy.Targets.item1.item2 <- abs(merge_withMeta_final_phylo$segSites.divBy.Targets.item1 - merge_withMeta_final_phylo$segSites.divBy.Targets.item2 )


merge_withMeta_final_phylo$absDiff.Longevity_d <- abs(as.numeric(merge_withMeta_final_phylo$Max_longevity_d.item1) - as.numeric(merge_withMeta_final_phylo$Max_longevity_d.item2 )) # some species missing info

merge_withMeta_final_phylo$absDiff.Rspan_d <- abs(as.numeric(merge_withMeta_final_phylo$Rspan_d.item1) - as.numeric(merge_withMeta_final_phylo$Rspan_d.item2 )) # some species missing info




confounderPlotsDir=paste0(plotdir,"/possibleConfounders/")
dir.create(confounderPlotsDir,showWarnings = F)

my.formula <- y~x

segSitesVsCopheneticDistance_plot <- ggplot(merge_withMeta_final_phylo[merge_withMeta_final_phylo$variable=="mutation_rate_projected_normalized_plusEpsilon" & merge_withMeta_final_phylo$transform_label=="clr" & merge_withMeta_final_phylo$id=="spectra_3mer",],aes(x=cophenetic_distance,y=absDiff.totalSegSites))+
  geom_point()+
  geom_text(aes(label=comparisonLabel_common_alphabetical))+
  geom_smooth(method="lm",se = F)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.95,
               formula = my.formula, parse = TRUE, size = 4,color='red')+
  ggtitle("difference in total segregating sites between species vs cophenetic distance")
segSitesVsCopheneticDistance_plot

ggsave(paste0(confounderPlotsDir,"totalSegSites.vs.Cophenetic.png"),segSitesVsCopheneticDistance_plot)

segSitesDivByTargetsVsCophenetic_plot <- ggplot(merge_withMeta_final_phylo[merge_withMeta_final_phylo$variable=="mutation_rate_projected_normalized_plusEpsilon" & merge_withMeta_final_phylo$transform_label=="clr" & merge_withMeta_final_phylo$id=="spectra_3mer",],aes(x=cophenetic_distance,y=absDiff.segSites.divBy.Targets.item1.item2))+
  geom_point()+
  geom_text(aes(label=comparisonLabel_common_alphabetical))+
  geom_smooth(method="lm",se = F)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.95,
               formula = my.formula, parse = TRUE, size = 4,color='red')+
  ggtitle("difference in total segregating sites divided by total target size (counts are already sample size corrected so didn't divide by harm number for now) vs. cophenetic distance")
segSitesDivByTargetsVsCophenetic_plot
ggsave(paste0(confounderPlotsDir,"totalSegSites.divBy.TotalTargets.vs.Cophenetic.png"),segSitesVsCopheneticDistance_plot)


segSitesDivByTargetsVsSpectrumDistance_plot <- ggplot(merge_withMeta_final_phylo[merge_withMeta_final_phylo$variable=="mutation_rate_projected_normalized_plusEpsilon" & merge_withMeta_final_phylo$transform_label=="clr",],aes(y=spectrum_distance,x=absDiff.segSites.divBy.Targets.item1.item2))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  #geom_text(aes(label=comparisonLabel_common_alphabetical))+
  geom_smooth(method="lm",se = F)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = my.formula, parse = TRUE, size = 4,color='red')+
  ggtitle("difference in total segregating sites divided by total target size\n(counts are already sample size corrected so didn't divide by harm number for now)\nvs. spectrum (clr) distance")

segSitesDivByTargetsVsSpectrumDistance_plot
ggsave(paste0(confounderPlotsDir,"totalSegSites.divBy.TotalTargets.vs.SpectrumDistance.nolabel.png"),segSitesDivByTargetsVsSpectrumDistance_plot)

segSitesDivByTargetsVsSpectrumDistance_plot_label + geom_text(aes(label=comparisonLabel_common_alphabetical))
ggsave(paste0(confounderPlotsDir,"totalSegSites.divBy.TotalTargets.vs.SpectrumDistance.withlabel.png"),segSitesDivByTargetsVsSpectrumDistance_plot_label)
# target size: 
targetsVsCopheneticDistance_plot <- ggplot(merge_withMeta_final_phylo[merge_withMeta_final_phylo$variable=="mutation_rate_projected_normalized_plusEpsilon" & merge_withMeta_final_phylo$transform_label=="clr" & merge_withMeta_final_phylo$id=="spectra_3mer",],aes(x=cophenetic_distance,y=absDiff.targetSize))+
  geom_point()+
  geom_smooth(method="lm",se = F)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = my.formula, parse = TRUE, size = 4,color='red')+
  ggtitle("difference in target size (callable sites of genome) and cophenetic distance")#+
  #geom_text(aes(label=comparisonLabel_common_alphabetical),size=3)
targetsVsCopheneticDistance_plot
ggsave(paste0(confounderPlotsDir,"totalGenomeTargetSize.vs.CopheneticDistance.png"),targetsVsCopheneticDistance_plot)


targetsVsSpectrumDistance_plot <- ggplot(merge_withMeta_final_phylo[merge_withMeta_final_phylo$variable=="mutation_rate_projected_normalized_plusEpsilon" & merge_withMeta_final_phylo$transform_label=="clr",],aes(y=spectrum_distance,x=absDiff.targetSize))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  ggtitle("difference in target size (callable sites of genome) and spectrum distance")+
  #geom_text(aes(label=comparisonLabel_common_alphabetical),size=3) +
  geom_smooth(method="lm",se = F)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = my.formula, parse = TRUE, size = 4,color='red')
targetsVsSpectrumDistance_plot
ggsave(paste0(confounderPlotsDir,"totalGenomeTargetSize.vs.SpectrumDistance.nolabel.png"),targetsVsSpectrumDistance_plot)


targetsVsSpectrumDistance_plot_label <- targetsVsSpectrumDistance_plot+geom_text(aes(label=comparisonLabel_common_alphabetical),size=3)
ggsave(paste0(confounderPlotsDir,"totalGenomeTargetSize.vs.SpectrumDistance.withlabel.png"),targetsVsSpectrumDistance_plot_label)


genTimeVsCopheneticDistance_plot <- ggplot(merge_withMeta_final_phylo[merge_withMeta_final_phylo$variable=="mutation_rate_projected_normalized_plusEpsilon" & merge_withMeta_final_phylo$transform_label=="clr" & merge_withMeta_final_phylo$id=="spectra_3mer",],aes(x=cophenetic_distance,y=absDiff.AFR_d.Item1.Item2))+
  geom_point()+
  #geom_smooth(method="lm")+
  ggtitle("difference in age at first reproduction (days) and cophenetic distance")+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=3)+
  geom_smooth(method="lm",se = F)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = my.formula, parse = TRUE, size = 4,color='red')

genTimeVsCopheneticDistance_plot
ggsave(paste0(confounderPlotsDir,"diffInAgeAtFirstRepro.vs.CopheneticDistance.MissingSomeSpeciesInfo.png"),genTimeVsCopheneticDistance_plot)

genTimeVsSpectrumDistance_plot <- ggplot(merge_withMeta_final_phylo[merge_withMeta_final_phylo$variable=="mutation_rate_projected_normalized_plusEpsilon" & merge_withMeta_final_phylo$transform_label=="clr",],aes(y=spectrum_distance,x=absDiff.AFR_d.Item1.Item2))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  #geom_smooth(method="lm")+
  ggtitle("difference in age at first reproduction (days) and spectrum distance [ some species missing info ]")+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=3)+
  geom_smooth(method="lm",se = F)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = my.formula, parse = TRUE, size = 4,color='red')

genTimeVsSpectrumDistance_plot
ggsave(paste0(confounderPlotsDir,"diffInAgeAtFirstRepro.vs.SpectrumDistance.MissingSomeSpeciesInfo.png"),genTimeVsSpectrumDistance_plot)


longevityVsSpectrumDistance_plot <- ggplot(merge_withMeta_final_phylo[merge_withMeta_final_phylo$variable=="mutation_rate_projected_normalized_plusEpsilon" & merge_withMeta_final_phylo$transform_label=="clr",],aes(y=spectrum_distance,x=absDiff.Longevity_d))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  #geom_smooth(method="lm")+
  ggtitle("difference in max longevity (days) and spectrum distance [ some species missing info ]")+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=3)+
  geom_smooth(method="lm",se = F)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = my.formula, parse = TRUE, size = 4,color='red')
longevityVsSpectrumDistance_plot
ggsave(paste0(confounderPlotsDir,"diffInMaxLongevity.vs.SpectrumDistance.MissingSomeSpeciesInfo.png"),longevityVsSpectrumDistance_plot)


rspanVsSpectrumDistance_plot <- ggplot(merge_withMeta_final_phylo[merge_withMeta_final_phylo$variable=="mutation_rate_projected_normalized_plusEpsilon" & merge_withMeta_final_phylo$transform_label=="clr",],aes(y=spectrum_distance,x=absDiff.Rspan_d))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  #geom_smooth(method="lm")+
  ggtitle("difference in repdroductive span (days) and spectrum distance [ some species missing info ]")+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=3)+
  geom_smooth(method="lm",se = F)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = my.formula, parse = TRUE, size = 4,color='red')
rspanVsSpectrumDistance_plot
ggsave(paste0(confounderPlotsDir,"diffInReproSpan.vs.SpectrumDistance.MissingSomeSpeciesInfo.png"),rspanVsSpectrumDistance_plot)
##################### LINEAR MODELING:  linear modeling of phylo and time-scaled times and spectrum distance #################
# want to lm 
# want to compare fits 
# and want to try with and without mice
# includes mice and phylo distance
# want to group by kmer type and transform 
require(broom)
# https://stackoverflow.com/questions/22713325/fitting-several-regression-models-with-dplyr - scroll down to 2020 update of solution
lmplotdir=paste0(plotdir,'/lmPlots/')
dir.create(lmplotdir)
########## get fit: #################
lm_phylo_1_summaryStats <- all_distances_all_conditions_phylo %>% group_by(transform_label,id,variable) %>%   
  do(m1_tidy = glance(lm(spectrum_distance ~ cophenetic_distance, data = .))) %>% 
  # use broom::glance to get summary stats
  unnest(m1_tidy) %>%
  mutate(distance_label="substitution phylogentic distance",mouse_label="mice included")
lm_phylo_1_summaryStats



# get neg log10 pvalue for comparisons
lm_phylo_1_summaryStats$negLog10Pval <- -log10(lm_phylo_1_summaryStats$p.value)
lm_phylo_1_summaryStats_melt <- melt(lm_phylo_1_summaryStats,variable.name  = "stat",id.vars = c("transform_label","id","variable"))


lm_time_1_summaryStats <- all_distances_all_conditions_time %>% unique() %>% group_by(transform_label,id,variable) %>%   
  do(m1_tidy = glance(lm(spectrum_distance ~ cophenetic_distance, data = .))) %>% 
  # use broom::glance to get summary stats
  unnest(m1_tidy) %>%
  mutate(distance_label="time scaled distance",mouse_label="mice included")
lm_time_1_summaryStats

lm_hamming_1_summaryStats <- 
  all_distances_all_conditions_hamming %>% group_by(transform_label,id,variable) %>%   
  do(m1_tidy = glance(lm(spectrum_distance ~ average_HammingDistance_divBy2xGenomeSize, data = .))) %>% 
  # use broom::glance to get summary stats
  unnest(m1_tidy) %>%
  mutate(distance_label="hamming distance",mouse_label="mice included")

lm_hamming_1_summaryStats

######### try without mice ########
# do without mice:
mice=c("mice_Mmd","mice_Ms","mice_Mmm","mice_Mmc") # note phylo nly has Mmd and Ms
lm_phylo_2_noMice_summaryStats <- all_distances_all_conditions_phylo %>% 
  filter(!item1 %in% mice, !item2 %in% mice) %>% # checked that this worked
  group_by(transform_label,id,variable) %>%   
  do(m1_tidy = glance(lm(spectrum_distance ~ cophenetic_distance, data = .))) %>% 
  # use broom::glance to get summary stats
  unnest(m1_tidy) %>%
  mutate(distance_label="substitution phylogentic distance",mouse_label="mice excluded")
lm_phylo_2_noMice_summaryStats

lm_time_2_noMice_summaryStats <- all_distances_all_conditions_time %>% 
  filter(!item1 %in% mice, !item2 %in% mice) %>% # checked that this worked
  group_by(transform_label,id,variable) %>%   
  do(m1_tidy = glance(lm(spectrum_distance ~ cophenetic_distance, data = .))) %>% 
  # use broom::glance to get summary stats
  unnest(m1_tidy) %>%
  mutate(distance_label="time scaled distance",mouse_label="mice excluded")
lm_time_2_noMice_summaryStats


lm_hamming_2_noMice_summaryStats <- all_distances_all_conditions_hamming %>% 
  filter(!item1 %in% mice, !item2 %in% mice) %>% # checked that this worked
  group_by(transform_label,id,variable) %>%   
  do(m1_tidy = glance(lm(spectrum_distance ~ average_HammingDistance_divBy2xGenomeSize, data = .))) %>% 
  # use broom::glance to get summary stats
  unnest(m1_tidy)  %>%
  mutate(distance_label="hamming distance",mouse_label="mice excluded")
lm_hamming_2_noMice_summaryStats

########## bind the model stats all together ############
lm_summaryStats_ALL <- bind_rows(lm_phylo_1_summaryStats,lm_time_1_summaryStats,lm_hamming_1_summaryStats,lm_phylo_2_noMice_summaryStats,lm_time_2_noMice_summaryStats,lm_hamming_2_noMice_summaryStats)

# write out the stats
write.table(lm_summaryStats_ALL,paste0(lmplotdir,"all.lmModels.SummaryStats.txt"),quote=F,sep = "\t",row.names = F)

########### plot lm stats ###################

lm_summaryStats_ALL$variable <- factor(lm_summaryStats_ALL$variable,levels=c("mutation_rate_projected_normalized_plusEpsilon","mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon"  ,"fractionOfSegregatingSites_projected_plusEpsilon","fractionOfSegregatingSites_projected_multinom_downsampled_plusEpsilon"))
lmPlot1_alt1 <- ggplot(lm_summaryStats_ALL[lm_summaryStats_ALL$mouse_label=="mice included",],aes(x=id,y=adj.r.squared,fill=variable))+
  geom_col(position="dodge")+
  facet_wrap(distance_label~transform_label,ncol=2)+
  theme_bw() +
  scale_fill_brewer(palette = "Paired",direction = -1)

lmPlot1_alt1
ggsave(paste0(lmplotdir,"RegressionPlot1-alt.Columns.png"),lmPlot1_alt1,height=12,width=18)
########### nice pag plot of rsq #############
# make a version for PAG talk that is just downsampled mut rate and is just phylo dist:
# with time and phylo:

lmPlot1_alt2a <- ggplot(lm_summaryStats_ALL[lm_summaryStats_ALL$mouse_label=="mice included" & lm_summaryStats_ALL$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & lm_summaryStats_ALL$distance_label =="substitution phylogentic distance" & lm_summaryStats_ALL$transform_label=="clr",],aes(x=id,y=adj.r.squared,fill=distance_label))+
  geom_col(position="dodge")+
  #facet_wrap(distance_label~transform_label,ncol=2)+
  theme_bw() +
  scale_fill_manual(values=c("#1F78B4"))+
  xlab("")+
  theme(text=element_text(size=14))

lmPlot1_alt2a
ggsave(paste0(lmplotdir,"RegressionPlot1-alt2.Columns.FORPAG.JustMultinomDownsampMutRate.JUSTphylo.clrOnly.png"),lmPlot1_alt2a,height=5,width=9)


lmPlot1_alt2b <- ggplot(lm_summaryStats_ALL[lm_summaryStats_ALL$mouse_label=="mice included" & lm_summaryStats_ALL$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & lm_summaryStats_ALL$distance_label %in% c("substitution phylogentic distance","time scaled distance") & lm_summaryStats_ALL$transform_label=="clr",],aes(x=id,y=adj.r.squared,fill=distance_label))+
  geom_col(position="dodge")+
  #facet_wrap(distance_label~transform_label,ncol=2)+
  theme_bw() +
  scale_fill_brewer(palette = "Paired",direction = -1)+
  xlab("")+
  theme(text=element_text(size=14),legend.title = element_blank())

lmPlot1_alt2b
ggsave(paste0(lmplotdir,"RegressionPlot1-alt2.Columns.FORPAG.JustMultinomDownsampMutRatePhyoAndTime.clrOnly.png"),lmPlot1_alt2b,height=5,width=9)



####### plot all lm metrics to gain intuition ########


lmPlot1b <- ggplot(lm_phylo_1_summaryStats_melt[lm_phylo_1_summaryStats_melt$stat %in% c("r.squared","adj.r.squared","negLog10Pval","logLik","AIC"),],aes(x=id,y=as.numeric(value),color=variable))+
  facet_wrap(~transform_label~stat,scales="free",ncol=5)+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))+
  xlab("")+
  ylab("")
lmPlot1b
ggsave(paste0(lmplotdir,"RegressionPlot1b.ComparingMetricsForPhyloModelOnly.png"),lmPlot1b,height=7,width=16)


######### plot with and without mice ##########
# paired colors
lm_summaryStats_ALL$mouse_label <- factor(lm_summaryStats_ALL$mouse_label,levels=c("mice included","mice excluded"))

lmPlot1_alt2_noMice <- ggplot(lm_summaryStats_ALL,aes(x=id,y=adj.r.squared,fill=interaction(mouse_label,variable)))+
  geom_col(position="dodge")+
  facet_wrap(distance_label~transform_label,ncol=2)+
  theme_bw()+
  scale_fill_brewer(palette = "Paired")+
  ggtitle(paste0("fit of linear model of spectrum dist ~ divergence\n",flagOfAnalysis))
lmPlot1_alt2_noMice

ggsave(paste0(lmplotdir,"RegressionPlot1-alt2.Columns.MOUSEEffect.png"),lmPlot1_alt2_noMice,height=12,width=18)
## make a nice version of this plot for pag:


lmPlot1_alt3_noMice_subsetOfConditions <- ggplot(lm_summaryStats_ALL[lm_summaryStats_ALL$transform_label=="clr" & lm_summaryStats_ALL$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & lm_summaryStats_ALL$distance_label!="hamming distance",],aes(x=id,y=adj.r.squared,fill=distance_label))+
  geom_col(position="dodge")+
  theme_bw()+
  facet_wrap(~mouse_label)+
  scale_fill_brewer(palette = "Paired",direction=-1)+
  ggtitle(paste0("fit of linear model of spectrum dist ~ divergence\n",flagOfAnalysis))+
  theme(text=element_text(size=14),legend.title = )+
  xlab("")
lmPlot1_alt3_noMice_subsetOfConditions

ggsave(paste0(lmplotdir,"RegressionPlot1-alt2.Columns.MOUSEEffect.NICEPLOT.FORPAG.png"),lmPlot1_alt3_noMice_subsetOfConditions,height=6,width=12)
### MAKE SOME FUNCTIONS 
######## get residuals:  ########
lm_phylo_1_fitAndResiduals <- all_distances_all_conditions_phylo %>% group_by(transform_label,id,variable) %>%   
  do(m1_aug = augment(lm(spectrum_distance ~ cophenetic_distance, data = .))) %>% 
  # use broom::augment to get fit of indivivual points
  unnest(m1_aug)  %>%
  mutate(label="phylo")

lm_time_1_fitAndResiduals <- all_distances_all_conditions_time %>% group_by(transform_label,id,variable) %>%   
  do(m1_aug = augment(lm(spectrum_distance ~ cophenetic_distance, data = .))) %>% 
  # use broom::augment to get fit of indivivual points
  unnest(m1_aug)  %>%
  mutate(label="time")

lm_hamming_1_fitAndResiduals <- all_distances_all_conditions_hamming %>% group_by(transform_label,id,variable) %>% 
  do(m1_aug = augment(lm(spectrum_distance ~ average_HammingDistance_divBy2xGenomeSize, data = .))) %>% 
  # use broom::augment to get fit of indivivual points
  unnest(m1_aug) %>%
  mutate(label="hamming")

head(lm_phylo_1_fitAndResiduals)

########### PLOT REGRESSIONS  #########
# https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
my.formula <- y ~ x # set formula
# set order or variables:

all_distances_all_conditions_phylo$variable <- factor(all_distances_all_conditions_phylo$variable,levels=c("mutation_rate_projected_normalized_plusEpsilon","mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" ,"fractionOfSegregatingSites_projected_plusEpsilon","fractionOfSegregatingSites_projected_multinom_downsampled_plusEpsilon"))
all_distances_all_conditions_time$variable <- factor(all_distances_all_conditions_time$variable,levels=c("mutation_rate_projected_normalized_plusEpsilon","mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" ,"fractionOfSegregatingSites_projected_plusEpsilon","fractionOfSegregatingSites_projected_multinom_downsampled_plusEpsilon"))
all_distances_all_conditions_hamming$variable <- factor(all_distances_all_conditions_hamming$variable,levels=c("mutation_rate_projected_normalized_plusEpsilon","mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" ,"fractionOfSegregatingSites_projected_plusEpsilon","fractionOfSegregatingSites_projected_multinom_downsampled_plusEpsilon"))

lmPlot2_phylo_FitToData <- ggplot(all_distances_all_conditions_phylo,aes(x=cophenetic_distance,y=spectrum_distance,color="data"))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE,show.legend=T,aes(color="model"))+
  #geom_point(data=lm_phylo_1_fitAndResiduals,aes(x=cophenetic_distance,y=.fitted,color="model"))+
  facet_wrap(~variable~transform_label~id,scales="free",ncol=4)+
  scale_color_manual(values=c("black","gray","gray"))+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = my.formula, parse = TRUE, size = 4)+
  theme_bw()+
  ggtitle(paste0("substitution phylogenetic distance regressions\n",flagOfAnalysis))+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=1,color="blue")+
  theme(legend.position="none")
lmPlot2_phylo_FitToData
ggsave(paste0(lmplotdir,"RegressionPlot2.PhyloLM.WithEquation.png"),lmPlot2_phylo_FitToData,height=20,width=22)

lmPlot2_time_FitToData <- ggplot(all_distances_all_conditions_time,aes(x=cophenetic_distance,y=spectrum_distance,color="data"))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE,show.legend=T,aes(color="lm"))+
  #geom_point(data=lm_time_1_fitAndResiduals,aes(x=cophenetic_distance,y=.fitted,color="model"))+
  facet_wrap(~variable~transform_label~id,scales="free",ncol=4)+
  scale_color_manual(values=c("black","gray","gray"))+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = my.formula, parse = TRUE, size = 4)+
  theme_bw()+
  ggtitle(paste0("time scaled distance regressions\n",flagOfAnalysis))+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=1,color="blue")+
  theme(legend.position="none")
lmPlot2_time_FitToData
ggsave(paste0(lmplotdir,"RegressionPlot2.TimeLM.WithEquation.png"),lmPlot2_time_FitToData,height=20,width=22)

lmPlot2_hamming_FitToData <- ggplot(all_distances_all_conditions_hamming,aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,color="data"))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE,show.legend=T,aes(color="lm"))+
  #geom_point(data=lm_hamming_1_fitAndResiduals,aes(x=cophenetic_distance,y=.fitted,color="model"))+
  facet_wrap(~variable~transform_label~id,scales="free",ncol=4)+
  scale_color_manual(values=c("black","gray","gray"))+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = my.formula, parse = TRUE, size = 4)+
  theme_bw()+
  ggtitle(paste0("hamming distance regressions\n",flagOfAnalysis))+
  geom_text(aes(label=comparisonLabel),size=1,color="blue")+
  theme(legend.position="none")
lmPlot2_hamming_FitToData
ggsave(paste0(lmplotdir,"RegressionPlot2.HammingLM.WithEquation.png"),lmPlot2_hamming_FitToData,height=20,width=22)

########### plot REGRESSIONS without mice #############
lmPlot3_phylo_FitToData_noMice <- ggplot(all_distances_all_conditions_phylo[!all_distances_all_conditions_phylo$item1 %in% mice & !all_distances_all_conditions_phylo$item2 %in% mice,],aes(x=cophenetic_distance,y=spectrum_distance,color="data"))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE,show.legend=T,aes(color="model"))+
  facet_wrap(~variable~transform_label~id,scales="free",ncol=4)+
  scale_color_manual(values=c("black","red","gray"))+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = my.formula, parse = TRUE, size = 4)+
  geom_point(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$item1 %in% mice | all_distances_all_conditions_phylo$item2 %in% mice,],aes(x=cophenetic_distance,y=spectrum_distance,color="excluded mice"))+ # adding inmouse points that don't go into the lm
  theme_bw()+
  ggtitle(paste0("substitution phylogenetic distance regressions\n",flagOfAnalysis,"\nComparisons involving mice (red) excluded from lm"))+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=1,color="blue")+
  theme(legend.position="none")

lmPlot3_phylo_FitToData_noMice
ggsave(paste0(lmplotdir,"RegressionPlot3.PhyloLM.WithEquation.NOMICE.png"),lmPlot3_phylo_FitToData_noMice,height=20,width=22)

lmPlot3_time_FitToData_noMice <- ggplot(all_distances_all_conditions_time[!all_distances_all_conditions_time$item1 %in% mice & !all_distances_all_conditions_time$item2 %in% mice,],aes(x=cophenetic_distance,y=spectrum_distance,color="data"))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE,show.legend=T,aes(color="model"))+
  facet_wrap(~variable~transform_label~id,scales="free",ncol=4)+
  scale_color_manual(values=c("black","red","gray"))+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = my.formula, parse = TRUE, size = 4)+
  geom_point(data=all_distances_all_conditions_time[all_distances_all_conditions_time$item1 %in% mice | all_distances_all_conditions_time$item2 %in% mice,],aes(x=cophenetic_distance,y=spectrum_distance,color="excluded mice"))+ # adding inmouse points that don't go into the lm
  theme_bw()+
  ggtitle(paste0("time scaled distance regressions\n",flagOfAnalysis,"\nComparisons involving mice (red) excluded from lm"))+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=1,color="blue")+
  theme(legend.position="none")

lmPlot3_time_FitToData_noMice
ggsave(paste0(lmplotdir,"RegressionPlot3.TimeLM.WithEquation.NOMICE.png"),lmPlot3_time_FitToData_noMice,height=20,width=22)

lmPlot3_hamming_FitToData_noMice <-  ggplot(all_distances_all_conditions_hamming[!all_distances_all_conditions_hamming$item1 %in% mice & !all_distances_all_conditions_hamming$item2 %in% mice,],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,color="data"))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE,show.legend=T,aes(color="model"))+
  facet_wrap(~variable~transform_label~id,scales="free",ncol=4)+
  scale_color_manual(values=c("black","red","gray"))+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = my.formula, parse = TRUE, size = 4)+
  geom_point(data=all_distances_all_conditions_hamming[all_distances_all_conditions_hamming$item1 %in% mice | all_distances_all_conditions_hamming$item2 %in% mice,],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,color="excluded mice"))+ # adding inmouse points that don't go into the lm
  theme_bw()+
  ggtitle(paste0("Hamming distance regressions\n",flagOfAnalysis,"\nComparisons involving mice (red) excluded from lm"))+
  geom_text(aes(label=comparisonLabel),size=1,color="blue")+
  theme(legend.position="none")

lmPlot3_hamming_FitToData_noMice
ggsave(paste0(lmplotdir,"RegressionPlot3.HammingLM.WithEquation.NOMICE.png"),lmPlot3_hamming_FitToData_noMice,height=20,width=22)

######## get coefficients ###########
lm_phylo_1_coefficients <-  all_distances_all_conditions_phylo %>% group_by(transform_label,id) %>%   
  do(m1_tidy = tidy(lm(spectrum_distance ~ cophenetic_distance, data = .))) %>% 
  # use broom::glance to get coefficients
  unnest(m1_tidy) 

lm_time_1_coefficients <-  all_distances_all_conditions_time %>% group_by(transform_label,id) %>%   
  do(m1_tidy = tidy(lm(spectrum_distance ~ cophenetic_distance, data = .))) %>% 
  # use broom::glance to get coefficients
  unnest(m1_tidy) 

lm_hamming_1_coefficients <-  all_distances_all_conditions_hamming %>% group_by(transform_label,id) %>%   
  do(m1_tidy = tidy(lm(spectrum_distance ~ average_HammingDistance_divBy2xGenomeSize, data = .))) %>% 
  # use broom::glance to get coefficients
  unnest(m1_tidy) 
# write out
write.table(lm_phylo_1_coefficients,paste0(lmplotdir,"phyloModels.Coefficients.txt"),quote=F,sep = "\t",row.names = F)
write.table(lm_time_1_coefficients,paste0(lmplotdir,"timeModels.Coefficients.txt"),quote=F,sep = "\t",row.names = F)
write.table(lm_hamming_1_coefficients,paste0(lmplotdir,"hammingModels.Coefficients.txt"),quote=F,sep ="\t",row.names = F)


############# make nice plot for a talk #################
# want to exclude fin whale because it's overfiltered (?) for now
# want to only plot clr and mutation rate (maybe just 3mer?)

# version with fin whale and with downsamp
lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$transform_label=="clr",],aes(x=cophenetic_distance,y=spectrum_distance,color="data"))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE,show.legend=T,aes(color="model"))+
  facet_wrap(~id,scales="free",ncol=2)+
  scale_color_manual(values=c("black","red","gray"))+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.07,
               formula = my.formula, parse = TRUE, size = 6)+
  theme_bw()+
  ggtitle(paste0("substitution phylogenetic distance regressions\n",flagOfAnalysis,"\nkeeping in fw; multinom downsampled mut rate"))+
  #geom_text(aes(label=comparisonLabel_common_alphabetical),size=2,color="blue")+
  theme(legend.position="none",text=element_text(size=14))
lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels
ggsave(paste0(lmplotdir,"RegressionPlot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.png"),lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels,height=8,width=12)
ggsave(paste0(lmplotdir,"RegressionPlot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.pdf"),lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels,height=8,width=12)

# make a version without the regression (Will suggestion)

lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_NOREGRESSION <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$transform_label=="clr",],aes(x=cophenetic_distance,y=spectrum_distance,color="data"))+
  geom_point()+
  facet_wrap(~id,scales="free",ncol=2)+
  scale_color_manual(values=c("black","red","gray"))+
  theme_bw()+
  ggtitle(paste0("substitution phylogenetic distance regressions\n",flagOfAnalysis,"\nkeeping in fw; multinom downsampled mut rate"))+
  #geom_text(aes(label=comparisonLabel_common_alphabetical),size=2,color="blue")+
  theme(legend.position="none",text=element_text(size=14))

lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_NOREGRESSION
ggsave(paste0(lmplotdir,"RegressionPlot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.NOREGRESSION.png"),lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_NOREGRESSION,height=8,width=9)
ggsave(paste0(lmplotdir,"RegressionPlot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.NOREGRESSION.pdf"),lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_NOREGRESSION,height=8,width=9)


########### make a version with cor.test information ########
require(ggpubr)
cortestPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_CorTestLabels <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$transform_label=="clr",],aes(x=cophenetic_distance,y=spectrum_distance,color="data"))+
  geom_point()+
  facet_wrap(~id,scales="free",ncol=2)+
  scale_color_manual(values=c("black","red","gray"))+
  stat_cor(method="pearson",cor.coef.name = "R",label.y.npc = 0.9)+
  stat_cor(method="spearman",cor.coef.name = "rho",label.y.npc=0.78)+
  theme_bw()+
  ggtitle(paste0("substitution phylogenetic distance regressions\n",flagOfAnalysis,"\nkeeping in fw; multinom downsampled mut rate"))+
  #geom_text(aes(label=comparisonLabel_common_alphabetical),size=2,color="blue")+
  theme(legend.position="none",text=element_text(size=14))
cortestPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_CorTestLabels

ggsave(paste0(lmplotdir,"CorTestPlot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.png"),cortestPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_CorTestLabels,height=8,width=9)
ggsave(paste0(lmplotdir,"CorTestPlot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.pdf"),cortestPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_CorTestLabels,height=8,width=9)
# good explanation of pearson v spearman https://towardsdatascience.com/clearly-explained-pearson-v-s-spearman-correlation-coefficient-ada2f473b8
# spearman rank (rho) tests for whether relationship is explained by monotonic function; doesn't have to be linear
# pearson tests for linear relationship with constant rate 

# can thin to reduce shared branch lengths
comparisonChoices1=c("Balaenoptera_physalus.Gorilla_gorilla" ,"Balaenoptera_physalus.Mus_musculus","Balaenoptera_physalus.Phocoena_sinus","Balaenoptera_physalus.Ursus_arctos" ,"Gorilla_gorilla.Homo_sapiens"  , "Gorilla_gorilla.Mus_musculus"     , "Gorilla_gorilla.Ursus_arctos" ,"Mus_musculus.Mus_spretus"  ,  "Pan_paniscus.Pan_troglodytes"  ,"Pongo_abelii.Pongo_pygmaeus" ,  "Ursus_arctos.Ursus_maritimus" ,"Mus_musculus.Ursus_arctos"   )

lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_reduced <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$transform_label=="clr" & all_distances_all_conditions_phylo$comparisonLabel %in% comparisonChoices1,],aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE,show.legend=T,aes(color="model"))+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.07,
               formula = my.formula, parse = TRUE, size = 6)+
  facet_wrap(~id,scales="free",ncol=2)+
  #scale_color_manual(values=c("black","red","gray"))+
  theme_bw()+
  ggtitle(paste0("substitution phylogenetic distance regressions\n",flagOfAnalysis,"\nkeeping in fw; multinom downsampled mut rate"))+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=2,color="blue")+
  theme(text=element_text(size=14))

lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_reduced

ggsave(paste0(lmplotdir,"RegressionPlot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.ReducedSpeciesSet.png"),lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_reduced,height=8,width=9)
ggsave(paste0(lmplotdir,"RegressionPlot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.ReducedSpeciesSet.pdf"),lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_reduced,height=8,width=9)
# make a version with pearson/spearman:

cortestPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_reduced <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$transform_label=="clr" & all_distances_all_conditions_phylo$comparisonLabel %in% comparisonChoices1,],aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point()+
  stat_cor(method="pearson",cor.coef.name = "R",label.y.npc = 0.9)+
  stat_cor(method="spearman",cor.coef.name = "rho",label.y.npc=0.78)+
  facet_wrap(~id,scales="free",ncol=2)+
  #scale_color_manual(values=c("black","red","gray"))+
  theme_bw()+
  ggtitle(paste0("substitution phylogenetic distance regressions\n",flagOfAnalysis,"\nkeeping in fw; multinom downsampled mut rate"))+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=2,color="blue")+
  theme(text=element_text(size=14))

cortestPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_reduced

ggsave(paste0(lmplotdir,"CorTest.Plot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.ReducedSpeciesSet.png"),cortestPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_reduced,height=8,width=9)
ggsave(paste0(lmplotdir,"CorTest.Plot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.ReducedSpeciesSet.pdf"),cortestPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_reduced,height=8,width=9)
# make a version with pearson/spearman:


  
  
# reduce even ore
# still some slight overlap of branch lengths though 
comparisonChoices2=c("Balaenoptera_physalus.Phocoena_sinus","Balaenoptera_physalus.Ursus_arctos" , "Gorilla_gorilla.Mus_musculus"     , "Mus_musculus.Mus_spretus"  ,  "Pan_paniscus.Pan_troglodytes"  ,"Pongo_abelii.Pongo_pygmaeus" ,  "Ursus_arctos.Ursus_maritimus")


lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_FURTHERreduced <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$transform_label=="clr" & all_distances_all_conditions_phylo$comparisonLabel %in% comparisonChoices2,],aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point()+
  stat_cor(method="pearson",cor.coef.name = "R",label.y.npc = 0.9)+
  stat_cor(method="spearman",cor.coef.name = "rho",label.y.npc=0.78)+
  facet_wrap(~id,scales="free",ncol=2)+
  #scale_color_manual(values=c("black","red","gray"))+
  theme_bw()+
  ggtitle(paste0("substitution phylogenetic distance regressions\n",flagOfAnalysis,"\nkeeping in fw; multinom downsampled mut rate"))+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=2,color="blue")+
  theme(text=element_text(size=14))

lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_FURTHERreduced

ggsave(paste0(lmplotdir,"RegressionPlot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.EVENMOREReducedSpeciesSet.png"),lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_FURTHERreduced,height=8,width=9)
ggsave(paste0(lmplotdir,"RegressionPlot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.EVENMOREReducedSpeciesSet.pdf"),lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_FURTHERreduced,height=8,width=9)

# version with spearman/pearson
CorTestPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_FURTHERreduced <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$transform_label=="clr" & all_distances_all_conditions_phylo$comparisonLabel %in% comparisonChoices2,],aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE,show.legend=T,aes(color="model"))+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.07,
               formula = my.formula, parse = TRUE, size = 6)+
  facet_wrap(~id,scales="free",ncol=2)+
  #scale_color_manual(values=c("black","red","gray"))+
  theme_bw()+
  ggtitle(paste0("substitution phylogenetic distance regressions\n",flagOfAnalysis,"\nkeeping in fw; multinom downsampled mut rate"))+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=2,color="blue")+
  theme(text=element_text(size=14))

CorTestPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_FURTHERreduced

ggsave(paste0(lmplotdir,"CorTest.Plot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.EVENMOREReducedSpeciesSet.png"),CorTestPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_FURTHERreduced,height=8,width=9)
ggsave(paste0(lmplotdir,"CorTest.Plot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.EVENMOREReducedSpeciesSet.pdf"),CorTestPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_FURTHERreduced,height=8,width=9)


# absolutely no shared branches: try to have diversity of branch lengths 
# gorilla-mouse longest; then whale bear ; then human -chimp, then orangorang
comparisonChoices3=c("Gorilla_gorilla.Mus_musculus", "Homo_sapiens.Pan_troglodytes" ,  "Pongo_abelii.Pongo_pygmaeus", "Balaenoptera_physalus.Ursus_arctos")

lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_NOSHAREDBRANCHES <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$transform_label=="clr" & all_distances_all_conditions_phylo$comparisonLabel %in% comparisonChoices3,],aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE,show.legend=T,aes(color="model"))+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.07,
               formula = my.formula, parse = TRUE, size = 6)+
  facet_wrap(~id,scales="free",ncol=2)+
  #scale_color_manual(values=c("black","red","gray"))+
  theme_bw()+
  ggtitle(paste0("substitution phylogenetic distance regressions\n",flagOfAnalysis,"\nkeeping in fw; multinom downsampled mut rate"))+
  theme(text=element_text(size=14))

lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_NOSHAREDBRANCHES

ggsave(paste0(lmplotdir,"RegressionPlot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.NOSHAREDBRANCHES.png"),lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_NOSHAREDBRANCHES,height=8,width=9)
ggsave(paste0(lmplotdir,"RegressionPlot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.NOSHAREDBRANCHES.pdf"),lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_NOSHAREDBRANCHES,height=8,width=9)

# version with labels:
lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_withLabels_NOSHAREDBRANCHES <- lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_NOSHAREDBRANCHES+geom_text(aes(label=comparisonLabel_common_alphabetical),size=2,color="blue")
ggsave(paste0(lmplotdir,"RegressionPlot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.WITHLABELS.NOSHAREDBRANCHES.png"),lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_withLabels_NOSHAREDBRANCHES,height=8,width=9)
ggsave(paste0(lmplotdir,"RegressionPlot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.WITHLABELS.NOSHAREDBRANCHES.pdf"),lmPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_withLabels_NOSHAREDBRANCHES,height=8,width=9)

# version with cor.test:
CorTestPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_NOSHAREDBRANCHES <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$transform_label=="clr" & all_distances_all_conditions_phylo$comparisonLabel %in% comparisonChoices3,],aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point()+
  stat_cor(method="pearson",cor.coef.name = "R",label.y.npc = 0.9)+
  stat_cor(method="spearman",cor.coef.name = "rho",label.y.npc=0.78)+
  facet_wrap(~id,scales="free",ncol=2)+
  #scale_color_manual(values=c("black","red","gray"))+
  theme_bw()+
  ggtitle(paste0("substitution phylogenetic distance regressions\n",flagOfAnalysis,"\nkeeping in fw; multinom downsampled mut rate"))+
  theme(text=element_text(size=14))

CorTestPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_NOSHAREDBRANCHES

ggsave(paste0(lmplotdir,"CorTest.Plot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.NOSHAREDBRANCHES.png"),CorTestPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_NOSHAREDBRANCHES,height=8,width=9)
ggsave(paste0(lmplotdir,"CorTest.Plot4a.NICESUBSETFORTALKS.PhyloLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.NOSHAREDBRANCHES.pdf"),CorTestPlot4a_phylo_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_NOSHAREDBRANCHES,height=8,width=9)



############## NEED TO DO COR TESTS EN MASSE HERE AND GET P VALUES ###########
### do plots with geom_cor from ggpubr ! 
# not significant bc too few points?
# maybe add spearmans rho and pearsons r
# testing: 
#cor.test(x=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$transform_label=="clr"  & all_distances_all_conditions_phylo$id=="spectra_1mer",]$cophenetic_distance,y=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$transform_label=="clr" & all_distances_all_conditions_phylo$id=="spectra_1mer",]$spectrum_distance,method="spearman")
# okay so that'sc ool



# trying something: scale up x-axis by 10 to be more similar to time-scaled and see if Rsq changes
# see if there's an effect of diff x axis scale
# expt_xTimes100 <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$transform_label=="clr",],aes(x=100*cophenetic_distance,y=spectrum_distance,color="data"))+
#   geom_point()+
#   geom_smooth(method="lm",se=FALSE,show.legend=T,aes(color="model"))+
#   facet_wrap(~id,scales="free",ncol=2)+
#   scale_color_manual(values=c("black","red","gray"))+
#   stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#                label.x.npc = "right", label.y.npc = 0.07,
#                formula = my.formula, parse = TRUE, size = 6)+
#   theme_bw()+
#   ggtitle(paste0("substitution phylogenetic distance regressions\n",flagOfAnalysis,"\nkeeping in fw; multinom downsampled mut rate"))+
#   #geom_text(aes(label=comparisonLabel_common_alphabetical),size=2,color="blue")+
#   theme(legend.position="none",text=element_text(size=14))
# expt_xTimes100

# make a version with hamming included:

lmPlot4a_hamming_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels <- ggplot(all_distances_all_conditions_hamming[all_distances_all_conditions_hamming$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_hamming$transform_label=="clr",],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,color="data"))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE,show.legend=T,aes(color="model"))+
  facet_wrap(~id,scales="free",ncol=2)+
  scale_color_manual(values=c("black","red","gray"))+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.1,
               formula = my.formula, parse = TRUE, size = 6)+
  theme_bw()+
  ggtitle(paste0("substitution phylogenetic distance regressions\n",flagOfAnalysis,"\nkeeping in fw; multinom downsampled mut rate"))+
  #geom_text(aes(label=comparisonLabel_common_alphabetical),size=2,color="blue")+
  theme(legend.position="none",text=element_text(size=14))
lmPlot4a_hamming_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels
ggsave(paste0(lmplotdir,"RegressionPlot4a.NICESUBSETFORTALKS.HammingLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.png"),lmPlot4a_hamming_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels,height=8,width=12)
# and save as pdf
ggsave(paste0(lmplotdir,"RegressionPlot4a.NICESUBSETFORTALKS.HammingLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.pdf"),lmPlot4a_hamming_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels,height=8,width=12)


############### with time scaled  mouse labels ####################3
all_distances_all_conditions_time$mouseLabel <- "no mouse"
all_distances_all_conditions_time[all_distances_all_conditions_time$item1 %in% mice,]$mouseLabel <- "mouse"
all_distances_all_conditions_time[all_distances_all_conditions_time$item2 %in% mice,]$mouseLabel <- "mouse"


lmPlot4b_time_NICEFORPLOTS_keepfw_labelmice_noLabels <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_time$transform_label=="clr",],aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point(aes(shape=mouseLabel,color=mouseLabel))+
  geom_smooth(method="lm",se=FALSE,show.legend=T,color="red")+
  facet_wrap(~id,scales="free",ncol=2)+
  scale_color_manual(values=c("blue","black"))+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.85,
               formula = my.formula, parse = TRUE, size = 6)+
  theme_bw()+
  ggtitle(paste0("time scaled distance regressions\n",flagOfAnalysis,"\nKeeping Fin Whale comparisons"))+
  scale_shape_manual(values=c(8,16))#+
  #geom_text(aes(label=comparisonLabel_common_alphabetical),size=2,color="blue")+
 # theme(legend.position="none")

lmPlot4b_time_NICEFORPLOTS_keepfw_labelmice_noLabels
ggsave(paste0(lmplotdir,"RegressionPlot4b.NICESUBSETFORTALKS.TimeLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.png"),lmPlot4b_time_NICEFORPLOTS_keepfw_labelmice_noLabels,height=8,width=12)
lmPlot4b_time_NICEFORPLOTS_keepfw_labelmice_noLabels
ggsave(paste0(lmplotdir,"RegressionPlot4b.NICESUBSETFORTALKS.TimeLM.WithEquation.MULTINOMDOWNSAMPLE.KeptFINWHALE.NOLABELS.pdf"),lmPlot4b_time_NICEFORPLOTS_keepfw_labelmice_noLabels,height=8,width=12)


lmPlot4b_time_NICEFORPLOTS_keepfw_DONOTlabelmice_noLabels <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_time$transform_label=="clr",],aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE,show.legend=T,color="red")+
  facet_wrap(~id,scales="free",ncol=2)+
  scale_color_manual(values=c("blue","black"))+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.85,
               formula = my.formula, parse = TRUE, size = 6)+
  theme_bw()+
  ggtitle(paste0("time scaled distance regressions\n",flagOfAnalysis,"\nKeeping Fin Whale comparisons"))+
  scale_shape_manual(values=c(8,16))#+


lmPlot4b_time_NICEFORPLOTS_keepfw_DONOTlabelmice_noLabels
ggsave(paste0(lmplotdir,"RegressionPlot4b.NICESUBSETFORTALKS.TimeLM.WithEquation.KeptFINWHALE.DONOTLABELMICe.NOLABELS.png"),lmPlot4b_time_NICEFORPLOTS_keepfw_DONOTlabelmice_noLabels,height=8,width=11)
ggsave(paste0(lmplotdir,"RegressionPlot4b.NICESUBSETFORTALKS.TimeLM.WithEquation.KeptFINWHALE.DONOTLABELMICe.NOLABELS.pdf"),lmPlot4b_time_NICEFORPLOTS_keepfw_DONOTlabelmice_noLabels,height=8,width=11)

########### try a hamming plot cutting species differences, just populations ##########
# no apes (not yet at least)
# so humans, mice, bears, fin whales
reducedHammingComparisons=c("bears_ABC","bears_EUR","fin_whale_GOC","fin_whale_ENP","humans_AFR","humans_AMR","humans_EAS","humans_SAS","humans_EUR","mice_Mmd","mice_Mmc","mice_Mmm")

lmPlot4a_hamming_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_WITHINSPECIESONLY <- ggplot(all_distances_all_conditions_hamming[all_distances_all_conditions_hamming$item1 %in%reducedHammingComparisons &  all_distances_all_conditions_hamming$item2 %in% reducedHammingComparisons & all_distances_all_conditions_hamming$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_hamming$transform_label=="clr",],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,color="data"))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE,show.legend=T,aes(color="model"))+
  facet_wrap(~id,scales="free",ncol=2)+
  scale_color_manual(values=c("black","red","gray"))+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.1,
               formula = my.formula, parse = TRUE, size = 6)+
  theme_bw()+
  geom_text(aes(label=comparisonLabel))
  ggtitle(paste0("substitution phylogenetic distance regressions\n",flagOfAnalysis,"\nkeeping in fw; multinom downsampled mut rate"))+
  #geom_text(aes(label=comparisonLabel_common_alphabetical),size=2,color="blue")+
  theme(legend.position="none",text=element_text(size=14))
lmPlot4a_hamming_NICEFORPLOTS_keepFinWhale_multinomDownsamp_noLabels_WITHINSPECIESONLY


################### MANTEL TEST ########################
# need dist matrix and phylo dist matrix of phylo distances
# may want to sq euclidean distances (bc expected to sale linearly with brownian motion apparently)
############ do for hamming as well! 
####### need as.dist objects; easiest to calc distances directly ? write out as dist mat? #####
#
manteldir=paste0(plotdir,"mantelTest/")
dir.create(manteldir)
comboFunction_mantelTest <- function(spectrum, phylo_dist,variable,speciesToInclude,mantelTestPermutationCount){
  spectrum_dist <- spectrum %>% 
    processSpectra() %>% 
    merge(.,speciesCodes,by.x="label",by.y="code",suffixes=c(".spectrum",".sppCodes")) %>% 
    filter(label %in% speciesToInclude) %>%
    pivotSpectra_perVariable_KEEPFULLSPECIESNAMES(variable) %>%
    clr_or_ilr_orNoTransform_calculations(.,"clr") %>% # fixing as clr -- could set as an option in the function tho
    euclideanDistance_keepMatrix_UPPERANDLOWER(.)  
  # the matrices need to be ordered in the same way 
  # reorder the time matrix in the same way as the spectrum matrix:
  phylo_dist_reordered <- as.matrix(phylo_dist)[rownames(as.matrix(spectrum_dist)),colnames(as.matrix(spectrum_dist))]
  head(phylo_dist_reordered)

  mtr = vegan::mantel(xdis=phylo_dist_reordered,ydis=spectrum_dist,permutations = mantelTestPermutationCount,method = "pearson")
  # make a data frame
  mantelTestResult_df <- data.frame(call="vegan_mantel",variable=variable,method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif))
  # mantel p value is calc as (ngreater+1)/(length(rsample)+1 so if there are none that are as well correlated then you'd get 0+1 / 999+1 = 1/ 10000 <- 0.001 which is what I get. which means in 999 permutations there are none that are as correlated. am I doing this right? 
  
  return(mantelTestResult_df)
  
}
# according to Hardy and Pavoine you get most power from euclidean dist vs sqrt(patristic)
# harmon and glor do euclidean squared vs patristic. either way you want one term to be squared or rooted because that is expectation for linearity under brownian motion
comboFunction_mantelTest_SQRTDISTANCE <- function(spectrum, phylo_dist,variable,speciesToInclude,mantelTestPermutationCount){
  spectrum_dist <- spectrum %>% 
    processSpectra() %>% 
    merge(.,speciesCodes,by.x="label",by.y="code",suffixes=c(".spectrum",".sppCodes")) %>% 
    filter(label %in% speciesToInclude) %>%
    pivotSpectra_perVariable_KEEPFULLSPECIESNAMES(variable) %>%
    clr_or_ilr_orNoTransform_calculations(.,"clr") %>% # fixing as clr -- could set as an option in the function tho
    euclideanDistance_keepMatrix_UPPERANDLOWER(.)  
  # the matrices need to be ordered in the same way 
  # reorder the time matrix in the same way as the spectrum matrix:
  phylo_dist_reordered <- as.matrix(phylo_dist)[rownames(as.matrix(spectrum_dist)),colnames(as.matrix(spectrum_dist))]
  head(phylo_dist_reordered)
  phylo_dist_reordered_sqrt <- sqrt(phylo_dist_reordered)
  mtr = vegan::mantel(xdis=phylo_dist_reordered_sqrt,ydis=spectrum_dist,permutations = mantelTestPermutationCount,method = "pearson")
  # make a data frame
  mantelTestResult_df <- data.frame(call="vegan_mantel_SQRTPHYLODISTANCE",variable=variable,method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif))
  # mantel p value is calc as (ngreater+1)/(length(rsample)+1 so if there are none that are as well correlated then you'd get 0+1 / 999+1 = 1/ 10000 <- 0.001 which is what I get. which means in 999 permutations there are none that are as correlated. am I doing this right? 
  
  return(mantelTestResult_df)
  
}

mantel_test_results_raxml_list = lapply(listOfDFs,comboFunction_mantelTest,phylo_dist= raxml_cophenetic_dist,variable="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon",speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount)
# bind into a df:
mantel_test_results_raxml_df = bind_rows(mantel_test_results_raxml_list,.id="id")
mantel_test_results_raxml_df$distance_metric <- "raxml_tree"

mantel_test_results_raxml_SQRTDISTANCE_list =  lapply(listOfDFs,comboFunction_mantelTest_SQRTDISTANCE,phylo_dist= raxml_cophenetic_dist,variable="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon",speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount)
mantel_test_results_raxml_SQRTDISTANCE_df = bind_rows(mantel_test_results_raxml_SQRTDISTANCE_list,.id="id")
mantel_test_results_raxml_SQRTDISTANCE_df$distance_metric <- "raxml_tree_SQRT"


mantel_test_results_time_list = lapply(listOfDFs,comboFunction_mantelTest,phylo_dist= timeTree_cophenetic_dist,variable="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon",speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount)
mantel_test_results_time_df = bind_rows(mantel_test_results_time_list,.id="id")
mantel_test_results_time_df$distance_metric <- "time_tree"

mantel_test_results_time_SQRTDISTANCE_list = lapply(listOfDFs,comboFunction_mantelTest_SQRTDISTANCE,phylo_dist= timeTree_cophenetic_dist,variable="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon",speciesToInclude_phylo,mantelTestPermutationCount=mantelTestPermutationCount)
mantel_test_results_time_SQRTDISTANCE_df = bind_rows(mantel_test_results_time_SQRTDISTANCE_list,.id="id")
mantel_test_results_time_SQRTDISTANCE_df$distance_metric <- "time_tree"


mantel_results_all <- bind_rows(mantel_test_results_raxml_SQRTDISTANCE_df,mantel_test_results_time_SQRTDISTANCE_df,mantel_test_results_raxml_df,mantel_test_results_time_df)
write.table(mantel_results_all,paste0(manteldir,"mantelTest.Results.txt"),row.names=F,quote=F,sep="\t")
############ plot mantel r values:  ############
# just raxml: 
mantel_r_plot_raxml <- 
  ggplot(mantel_test_results_raxml_df,aes(x=id,y=statistic,fill=distance_metric))+
  geom_col(position="dodge")+
  #facet_wrap(distance_label~transform_label,ncol=2)+
  theme_bw() +
  scale_fill_manual(values=c("#1F78B4"))+
  xlab("")+
  theme(text=element_text(size=14),legend.title = element_blank())+
  ggtitle("Mantel Test: r statistic (Pearson)")+
  ylab("correlation")
mantel_r_plot_raxml
ggsave(paste0(manteldir,"mantelTest.r.stat.PerSpectrum.raxmlTreeOnly.pdf"),mantel_r_plot_raxml,height=4,width=7)
# want to add mantel test information to plots: 

######## raxml and time :
mantel_r_plot_raxml_both <-   ggplot(mantel_results_all,aes(x=id,y=statistic,fill=distance_metric))+
  geom_col(position="dodge")+
  #facet_wrap(distance_label~transform_label,ncol=2)+
  theme_bw() +
  scale_fill_brewer(palette = "Paired",direction = -1)+
  xlab("")+
  theme(text=element_text(size=14),legend.title = element_blank())+
  ggtitle("Mantel Test: r statistic (Pearson)")+
  ylab("correlation")
mantel_r_plot_raxml_both
ggsave(paste0(manteldir,"mantelTest.r.stat.PerSpectrum.raxml.and.time.pdf"),mantel_r_plot_raxml_both,height=4,width=7)
# want to add mantel test information to plots: 




##################### plot correlations with mantel coefficients #############
####### PLOTTING BOTh SQRT AND NON SQRT DISTANCE ala Hardy and Pavoine. SQRT cophenetic is best!
mantelTestCorrelationsPlot_phylo <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$transform_label=="clr",],aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_df,aes(x=1.25,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",round(statistic,2),"\np-value: ",as.character(signif)," (Mantel Test)")),size=5)+
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_df$permCount)," permutations"))+
  theme_bw()

mantelTestCorrelationsPlot_phylo
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.pdf"),mantelTestCorrelationsPlot_phylo,height=8,width=11)


mantelTestCorrelationsPlot_phylo_SQRTDISTANCE <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_df,aes(x=1,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",round(statistic,2),"\np-value: ",as.character(signif)," (Mantel Test)")),size=5)+
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance"))+
  theme_bw()

mantelTestCorrelationsPlot_phylo_SQRTDISTANCE
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCE,height=8,width=11)



mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  geom_text(data=mantel_test_results_raxml_SQRTDISTANCE_df,aes(x=1,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",round(statistic,2),"\np-value: ",as.character(signif)," (Mantel Test)")),size=5)+
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance"))+
  theme_bw()+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=1)

mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.labels.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCEb,height=8,width=11)

# compare with and without sqrt patterns:

mantelTestCorrelationsPlot_phylo_SQRTDISTANCEc <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(aes(color="sqrt"))+
  geom_point(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_phylo$transform_label=="clr",],aes(x=cophenetic_distance,y=spectrum_distance,color="nonsqrt"))+
  facet_wrap(~id,scales="free")+
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_SQRTDISTANCE_df$permCount)," permutations-- SQRT of phylo distance"))+
  theme_bw()+
  xlab("sqrt or non-sqrt cophenetic distance")

mantelTestCorrelationsPlot_phylo_SQRTDISTANCEc
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.phylo.raxml.SQRTDISTANCE.COMPARISON.pdf"),mantelTestCorrelationsPlot_phylo_SQRTDISTANCEc,height=8,width=11)

########### mantel plot: time tree ###################
mantelTestCorrelationsPlot_time <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_time$transform_label=="clr",],aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  #geom_text(data=mantel_test_results_raxml_df,aes(x=1.25,y=-Inf,hjust=.5,vjust=-.5,label=paste0("Mantel test r (Pearson):  ",round(statistic,2),"\np-value: ",as.character(signif)," (",permCount," permutations)")))+
  geom_text(data=mantel_test_results_time_df,aes(x=150,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",round(statistic,2),"\np-value: ",as.character(signif)," (Mantel Test)")),size=5)+
  ggtitle(paste0("Distances based on time tree\nCorrelation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_df$permCount)," permutations"))+
  theme_bw()
mantelTestCorrelationsPlot_time
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.timetree.pdf"),mantelTestCorrelationsPlot_time,height=8,width=11)


mantelTestCorrelationsPlot_time_SQRTDISTANCE <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_time$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  #geom_text(data=mantel_test_results_raxml_df,aes(x=1.25,y=-Inf,hjust=.5,vjust=-.5,label=paste0("Mantel test r (Pearson):  ",round(statistic,2),"\np-value: ",as.character(signif)," (",permCount," permutations)")))+
  geom_text(data=mantel_test_results_time_SQRTDISTANCE_df,aes(x=12,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",round(statistic,2),"\np-value: ",as.character(signif)," (Mantel Test)")),size=5)+
  ggtitle(paste0("Distances based on time tree\nCorrelation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_time_SQRTDISTANCE_df$permCount)," permutations"))+
  theme_bw()
mantelTestCorrelationsPlot_time_SQRTDISTANCE
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.timetree.SQRTDISTANCE.pdf"),mantelTestCorrelationsPlot_time_SQRTDISTANCE,height=8,width=11)

########### mantel test with mouse labels ##########
mantelTestCorrelationsPlot_time_mouselabel <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_time$transform_label=="clr",],aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point(aes(color=mouseLabel,shape=mouseLabel))+
  facet_wrap(~id,scales="free")+
  #geom_text(data=mantel_test_results_raxml_df,aes(x=1.25,y=-Inf,hjust=.5,vjust=-.5,label=paste0("Mantel test r (Pearson):  ",round(statistic,2),"\np-value: ",as.character(signif)," (",permCount," permutations)")))+
  #geom_text(data=mantel_test_results_time_df,aes(x=150,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",round(statistic,2),"\np-value: ",as.character(signif)," (Mantel Test)")),size=5)+
  ggtitle(paste0("Distances based on time tree\nCorrelation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_df$permCount)," permutations"))+
  theme_bw()+
  scale_color_manual(values=c("blue","black"))+
  scale_shape_manual(values=c(5,16))+
  theme(legend.position = "none")
mantelTestCorrelationsPlot_time_mouselabel 
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.timetree.mouseLabel.pdf"),mantelTestCorrelationsPlot_time_mouselabel,height=8,width=11)


mantelTestCorrelationsPlot_time_mouselabel_SQRTDISTANCE <- ggplot(all_distances_all_conditions_time[all_distances_all_conditions_time$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_time$transform_label=="clr",],aes(x=sqrt(cophenetic_distance),y=spectrum_distance))+
  geom_point(aes(color=mouseLabel,shape=mouseLabel))+
  facet_wrap(~id,scales="free")+
  ggtitle(paste0("Distances based on time tree\nCorrelation coefficients from Mantel Test (Pearson's r)\n",unique(mantel_test_results_raxml_df$permCount)," permutations"))+
  theme_bw()+
  scale_color_manual(values=c("blue","black"))+
  scale_shape_manual(values=c(5,16))+
  theme(legend.position = "none")
mantelTestCorrelationsPlot_time_mouselabel_SQRTDISTANCE
ggsave(paste0(manteldir,"mantelTest.CorrelationPlot.timetree.mouseLabel.SQRTDISTANCE.pdf"),mantelTestCorrelationsPlot_time_mouselabel_SQRTDISTANCE,height=8,width=11)

############ mantel test on hamming distances ###########
#issue: don't have a full hamming distance matrix amongst all pairs; instead
# need to just use mantel test on apes alone, humans alone, mice alone #
# so need to get a matrix of hamming distances:
# adding labeltoinclude so it's like 'apes' or 'humans'
comboFunction_mantelTest_hammingWithinSpecies <- function(spectrum, hammingDistances,variable,labelToInclude,mantelTestPermutationCount){
  
  
  spectrum_dist <- spectrum %>% 
    processSpectra() %>% 
    filter(grepl(labelToInclude,label)) %>%
    pivotSpectra_perVariable(.,variable)  %>%
    clr_or_ilr_orNoTransform_calculations(.,"clr") %>% # fixing as clr -- could set as an option in the function tho
    euclideanDistance_keepMatrix_UPPERANDLOWER(.)  
  # the matrices need to be ordered in the same way 
  # reorder the time matrix in the same way as the spectrum matrix:
  # need to make a hamming distance matrix. (need it to be symmetrical)
  # hmmmm
  # pivot it wide: 
  hamming_wide <- data.frame(pivot_wider(hammingDistances[hammingDistances$label==labelToInclude,c("pop1","pop2","average_HammingDistance_divBy2xGenomeSize"),],names_from = c("pop2"),values_from = "average_HammingDistance_divBy2xGenomeSize"))
  rownames(hamming_wide) <- hamming_wide$pop1
  hamming_wide_mat <- as.matrix(select(hamming_wide,-"pop1"))
  # here the diag isn't zero and it's not a sq matrix:
  # note: diagonal isn't zero because there is internal hamming distances within a pop
  # but I'm setting diagonal to zero here so that it is similar to phylo distance
  # so set lower tri to 0 and then add together:
  hamming_wide_mat[!upper.tri(hamming_wide_mat)] <- 0
  hamming_wide_mat_SQUARE <- hamming_wide_mat + t(hamming_wide_mat) # now add to its own transpose to make it sq. code from: https://stackoverflow.com/questions/25911067/convert-triangular-parwise-matrix-to-symmetric-matrix-with-0-diagonal-r
  # checked it! yes this worked! 
  
  hamming_wide_mat_SQUARE_reordered <- as.matrix(hamming_wide_mat_SQUARE)[rownames(as.matrix(spectrum_dist)),colnames(as.matrix(spectrum_dist))]
  head(hamming_wide_mat_SQUARE_reordered)

  mtr = vegan::mantel(xdis=hamming_wide_mat_SQUARE_reordered,ydis=spectrum_dist,permutations = mantelTestPermutationCount,method = "pearson")
  # make a data frame
  mantelTestResult_df <- data.frame(call="vegan_mantel",variable=variable,method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif))
  # mantel p value is calc as (ngreater+1)/(length(rsample)+1 so if there are none that are as well correlated then you'd get 0+1 / 999+1 = 1/ 10000 <- 0.001 which is what I get. which means in 999 permutations there are none that are as correlated. am I doing this right? 
  
  return(mantelTestResult_df)
  
}

# do apes, humans and mice

all_hamming_mantel_results <- data.frame()
for(label in c("apes","humans","mice")){
hamming_mantel_perspecies_list <- lapply(listOfDFs,comboFunction_mantelTest_hammingWithinSpecies,hammingDistances= hammingDistances,variable="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon",label,mantelTestPermutations)

hamming_mantel_perspecies <- bind_rows(hamming_mantel_perspecies_list, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 

hamming_mantel_perspecies$label <- label 

# plot:
# positioning for geom_text is mean of hamming distance
geom_text_xval = mean(all_distances_all_conditions_hamming[all_distances_all_conditions_hamming$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_hamming$transform_label=="clr" & all_distances_all_conditions_hamming$label==label,]$average_HammingDistance_divBy2xGenomeSize)

mantelHamming_plot <- ggplot(all_distances_all_conditions_hamming[all_distances_all_conditions_hamming$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & all_distances_all_conditions_hamming$transform_label=="clr" & all_distances_all_conditions_hamming$label==label,],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance))+
  geom_point()+
  facet_wrap(~id,scales="free")+
  geom_text(aes(label=populationComparisonLabel_alphabetical),size=2)+
  geom_text(data=hamming_mantel_perspecies,aes(x=geom_text_xval,y=-Inf,hjust=.3,vjust=-.5,label=paste0("Pearson's r:  ",round(statistic,2),"\np-value: ",as.character(round(signif,3))," (Mantel Test)\n",permCount," permutations")),size=3)+
  ggtitle(paste0("Correlation coefficients from Mantel Test (Pearson's r)"))+
  theme_bw()
mantelHamming_plot
ggsave(paste0(manteldir,"hamming.",label,".MantelTest.png"),mantelHamming_plot,height=7,width=9)
all_hamming_mantel_results <- bind_rows(all_hamming_mantel_results,hamming_mantel_perspecies)

}


write.table(all_hamming_mantel_results,paste0(manteldir,"hammingtest.APES.HUMANS.MICE.txt"),row.names = F,quote=F,sep="\t")



# so here's the issue -- there aren't enough points to be significant!
#### to do: save these results. think about a different way to test significance.
# maybe PICs? but too few? 
# maybe apes+humans but not mantel test? stick with phylogeny? need human-ape distances?
# THINK !
# maybe pics? maybe something else? need distance matrix?
############### want to try to use phylosignal to look at changes in mutation types across phylogeny #########
#require(phylosignal) # https://cran.r-project.org/web/packages/phylosignal/vignettes/Basics.html
#require(adephylo)
#require(phylobase)
#dataForPhylosignal <- all1merSpectra %>%
#  processSpectra(.) %>%
#    merge(.,speciesCodes,by.x="label",by.y="code",suffixes=c(".spectrum",".sppCodes")) %>% 
#    filter(label %in% speciesToInclude_phylo) %>%
#    select(mutation_1mer,species.sppCodes,mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon) %>%
#  pivot_wider(., id_cols=species.sppCodes, values_from = #mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon,names_from = mutation_1mer) 

# set row names
#rownames(dataForPhylosignal) <- dataForPhylosignal$species.sppCodes # this odesn't work ; it's not just row names; you need it to be in order of tips of tree! annoying!
#dataForPhylosignal <- dataForPhylosignal[order(match(dataForPhylosignal$species.sppCodes,timeTree$tip.label)),] # reorder

# okay cool


# wants a time tree not a phylo tree (?)
#p4d <- phylo4d(timeTree,dataForPhylosignal) #oh these somehow have to be in same order! doesn't match on spp label. that is lame! 
#p4d
#barplot.phylo4d(p4d, tree.type = "phylo", tree.ladderize = TRUE,center = F,scale = F,label)
# THIS ISN'T WORKING AND IS WEIRD
# THIS IS WRONG UNTIL YOU CAN GET dATA IN SAME ORDER 
# how annoying 
# OKAY SOMETHING IS WRONG HERE ; it's not adding the data correctly. 
####### make a heatmap ##########
# use the downsampled metric
pivot1merForHeatmap <- all1merSpectra %>%
  processSpectra() %>%
  pivotSpectra_perVariable("mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon")
pivot3merForHeatmap <- all3merSpectra %>%
  processSpectra() %>%
  pivotSpectra_perVariable("mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon")
pivot5merForHeatmap <- all5merSpectra %>%
  processSpectra() %>%
  pivotSpectra_perVariable("mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon")
pivot7merForHeatmap <- all7merSpectra %>%
  processSpectra() %>%
  pivotSpectra_perVariable("mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon")


##### make a heatmap comparison of a few species
# ape/human
pivot3merForHeatmap$apes_Pan_troglodytes_over_humans_AFR <- pivot3merForHeatmap$apes_Pan_troglodytes / pivot3merForHeatmap$humans_AFR
# ape / mouse 
pivot3merForHeatmap$apes_Pan_troglodytes_over_mice_Mmd <- pivot3merForHeatmap$apes_Pan_troglodytes / pivot3merForHeatmap$mice_Mmd
# ape / whale
pivot3merForHeatmap$apes_Pan_troglodytes_over_vaquita <- pivot3merForHeatmap$apes_Pan_troglodytes / pivot3merForHeatmap$vaquita
#human/ mouse
pivot3merForHeatmap$apes_humans_AFR_over_mice_Mmd<- pivot3merForHeatmap$humans_AFR / pivot3merForHeatmap$mice_Mmd
# vaquita/mouse
pivot3merForHeatmap$vaquita_over_mice_Mmd<- pivot3merForHeatmap$vaquita / pivot3merForHeatmap$mice_Mmd
# polar bear / mouse 
pivot3merForHeatmap$bears_PB_over_mice_Mmd<- pivot3merForHeatmap$bears_PB / pivot3merForHeatmap$mice_Mmd
# polar bear / bb
pivot3merForHeatmap$bears_PB_over_bears_ABC <- pivot3merForHeatmap$bears_PB / pivot3merForHeatmap$bears_ABC

# try other mouse
pivot3merForHeatmap$mice_Ms_over_mice_Mmd <- pivot3merForHeatmap$mice_Ms / pivot3merForHeatmap$mice_Mmd

comparisons=c("apes_Pan_troglodytes_over_humans_AFR","apes_Pan_troglodytes_over_mice_Mmd","apes_Pan_troglodytes_over_vaquita","apes_humans_AFR_over_mice_Mmd","vaquita_over_mice_Mmd","bears_PB_over_mice_Mmd","bears_PB_over_bears_ABC","mice_Ms_over_mice_Mmd")

# make heatmap:
pivot3merForHeatmap$centralBP <- paste0(substr(pivot3merForHeatmap$mutation_label,2,2),".",substr(pivot3merForHeatmap$mutation_label,6,6))

pivot3merForHeatmap$five_prime_flanking_base = substr(pivot3merForHeatmap$mutation_label,1,1)
pivot3merForHeatmap$three_prime_flanking_base = substr(pivot3merForHeatmap$mutation_label,3,3)

pivot3merForHeatmap_melt_subset <- melt(pivot3merForHeatmap[,c("mutation_label","centralBP","five_prime_flanking_base","three_prime_flanking_base",comparisons)])

pivot3merForHeatmap_melt_subset$nicerCentralBPLabel <- gsub("\\.",">",pivot3merForHeatmap_melt_subset$centralBP,)

heatmap_lotsOfComparisons <- ggplot(pivot3merForHeatmap_melt_subset,aes(x=three_prime_flanking_base,y=five_prime_flanking_base,fill=value))+
  geom_tile()+
  facet_grid(~nicerCentralBPLabel~variable,switch="y")+
  scale_fill_gradient2(high="red",mid="white",low="blue",midpoint = 1)+
  theme_bw()+
  geom_text(aes(label=mutation_label),size=1)
heatmap_lotsOfComparisons
ggsave(paste0(plotdir,"heatmap.LotsOfComparisons.pdf"),heatmap_lotsOfComparisons,height=10,width=16)

# smaller subet:
pivot3merForHeatmap_melt_subset_smaller <- pivot3merForHeatmap_melt_subset[pivot3merForHeatmap_melt_subset$variable %in% c("bears_PB_over_bears_ABC","apes_humans_AFR_over_mice_Mmd","apes_Pan_troglodytes_over_vaquita"),]

pivot3merForHeatmap_melt_subset_smaller$variable <- factor(pivot3merForHeatmap_melt_subset_smaller$variable,levels=c("bears_PB_over_bears_ABC","apes_humans_AFR_over_mice_Mmd","apes_Pan_troglodytes_over_vaquita"))
# nicer version for for a talk
heatmap_niceVersion <- ggplot(pivot3merForHeatmap_melt_subset_smaller,aes(x=three_prime_flanking_base,y=five_prime_flanking_base,fill=value))+
  geom_tile()+
  facet_grid(~nicerCentralBPLabel~variable,switch = "y")+
  scale_fill_gradient2(high="red",mid="white",low="blue",midpoint = 1)+
  theme_bw()+
  geom_text(aes(label=mutation_label),size=2)+
  theme(strip.placement = "outside",strip.background = element_rect(fill="white"))
heatmap_niceVersion
ggsave(paste0(plotdir,"heatmap.3Comparisons.pdf"),heatmap_niceVersion,height=7,width=12)

###### try to show a different way? not by ratios?
head(pivot1merForHeatmap)
pivot1merForHeatmap_melt <- melt(pivot1merForHeatmap,variable.name = "spp")
head(pivot1merForHeatmap_melt)
looking_at_1mer_rates <- ggplot(pivot1merForHeatmap_melt,aes(y=spp,x=value,fill=mutation_label))+
  geom_col()+
  facet_wrap(~mutation_label,scales="free_x")
looking_at_1mer_rates
ggsave(paste0(plotdir,"1mer.mutrate.mulitnomdownsamp.perspcies.LookatC.G.png"),looking_at_1mer_rates,height=8,width=12)

# maybe plot C>G fraction vs Rspan?

############### more trying to troubleshoot fin whale /vaquita ##########
# something must be wrong with fw? so filtered, so little data? cpg islands? mapped to outgroup genome? lots of questions. polarization weird? I don't know!! mask it? 
ggplot(pivot1merForHeatmap,aes(x=vaquita,y=fin_whale_ENP))+
  geom_point(aes(color="1-mer"))+
  geom_text(aes(label=mutation_label))+
  geom_point(data=pivot3merForHeatmap,aes(x=vaquita,y=fin_whale_ENP,color="3-mer"))+
  geom_point(data=pivot5merForHeatmap,aes(x=vaquita,y=fin_whale_ENP,color="5-mer"))+
  #geom_point(data=pivot7merForHeatmap,aes(x=vaquita,y=fin_whale_GOC,color="7-mer"))+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline()+
  theme_bw()


pivot5merForHeatmap_fracSeg <- all5merSpectra %>%
  processSpectra() %>%
  pivotSpectra_perVariable("fractionOfSegregatingSites_projected_multinom_downsampled_plusEpsilon")
# plot a heatmap


pivot3merForHeatmap_fracSeg <- all3merSpectra %>%
  processSpectra() %>%
  pivotSpectra_perVariable("fractionOfSegregatingSites_projected_multinom_downsampled_plusEpsilon")
# plot a heatmap


ggplot(pivot5merForHeatmap_fracSeg,aes(x=humans_AFR,y=fin_whale_ENP))+
  geom_point(aes(color="5-mer"))+
  #scale_x_log10()+
  #scale_y_log10()+
  geom_abline()+
  geom_text(aes(label=mutation_label))+
  theme_bw()

########### write out downsampled counts in SigProfilerExtractor format ###########
# need to also do this in sigfit format
#require(spgs)  # for reverse complements

#Mutation Types	PD4199a	PD4005a	PD3851a	PD4116a	PD4086a	PD4194a	PD4248a	PD4120a	PD4198a	PD3904a	PD3945a	PD4107a	PD3905a	PD4192a	PD4109a	PD4103a	PD4115a	PD4085a	PD3890a	PD4006a	PD4088a
#A[C>A]A	58	74	31	128	31	28	64	165	58	122	243	210	94	48	122	112	228	61	110	198	34
#A[C>A]C	36	66	34	138	19	24	43	55	72	112	163	176	69	42	133	50	169	51	91	173	18
#A[C>A]G	13	12	9	15	8	6	7	20	12	13	24	33	11	16	12	7	21	7	33	3
#A[C>A]T	37	64	21	142	23	13	35	65	45	107	155	176	65	40	96	44	158	49	87	192	26
# # have to rev comp those with central -A- ancestral state to be ancestral T state
# have to fully rev comp both mutations 
# so A[A>G]C (AAC > AGC) is going to be come GTT > GCT so that the center position is T
# but don't want to change direction of mutation (ancestral is still ancestral, derived is still derived)
# example of rev comping:
# separate ancestral and derived 3mers:
#all3merSpectra_processed <- processSpectra(all3merSpectra)
## already has anc3mer in there
#all3merSpectra_processed$derived3mer <- substr(all3merSpectra_processed$mutation_label,5,7)
#head(all3merSpectra_processed)

# get ancestral central BP: 
#all3merSpectra_processed$ancestralCentralBP <- substr(all3merSpectra_processed$mutation_label,2,2)
#head(all3merSpectra_processed)
#all3merSpectra_processed$derivedCentralBP <- substr(all3merSpectra_processed$mutation_label,6,6)
#head(all3merSpectra_processed)
# get the reverse complement of all of them (then going to put together)

#all3merSpectra_processed$ancestral3mer.RevComplement <- toupper(reverseComplement(all3merSpectra_processed$ancestral3mer))

#all3merSpectra_processed$derived3mer.RevComplement <- toupper(reverseComplement(all3merSpectra_processed$derived3mer))

# so to be in sigprofiler format, if the ancestral central BP is an "A" then it needs to be rev complemented so that central BP is a "T"
# if central bp is a C then it's already good to go for sigprofiler
# note here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6717374/
# 
#all3merSpectra_processed$`Mutation Types` <- ""

# so if the center bp is a "C" it is good to go (no rev comping needed)
#all3merSpectra_processed[all3merSpectra_processed$ancestralCentralBP=="C",]$`Mutation Types` <-  paste0(substr(all3merSpectra_processed[all3merSpectra_processed$ancestralCentralBP=="C",]$mutation_label,1,1),"[",substr(all3merSpectra_processed[all3merSpectra_processed$ancestralCentralBP=="C",]$mutation_label,2,2),">",substr(all3merSpectra_processed[all3merSpectra_processed$ancestralCentralBP=="C",]$mutation_label,6,6),"]",substr(all3merSpectra_processed[all3merSpectra_processed$ancestralCentralBP=="C",]$mutation_label,3,3))

# but if the center bp is an "A" then you need to form the rev complement:
#all3merSpectra_processed[all3merSpectra_processed$ancestralCentralBP=="A",]$`Mutation Types` <-  paste0(substr(all3merSpectra_processed[all3merSpectra_processed$ancestralCentralBP=="A",]$ancestral3mer.RevComplement,1,1),"[",substr(all3merSpectra_processed[all3merSpectra_processed$ancestralCentralBP=="A",]$ancestral3mer.RevComplement,2,2),">",substr(all3merSpectra_processed[all3merSpectra_processed$ancestralCentralBP=="A",]$derived3mer.RevComplement,2,2),"]",substr(all3merSpectra_processed[all3merSpectra_processed$ancestralCentralBP=="A",]$ancestral3mer.RevComplement,3,3))


#head(all3merSpectra_processed)
# pivot: 
# this is the downsampled version:
#all3merSpectra_processed_wider_multiNomDownsampled <- pivot_wider(select(all3merSpectra_processed,c(`Mutation Types`,label,total_mutations_projected_multinom_downsampled)),id_cols = `Mutation Types`,names_from =  label,values_from = total_mutations_projected_multinom_downsampled)
#head(all3merSpectra_processed_wider_multiNomDownsampled)                    
 
#dir.create(paste0(plotdir,"/SigProfilerExtractorFiles/"),showWarnings = F)
#write.table(all3merSpectra_processed_wider_multiNomDownsampled,paste0(plotdir,"/SigProfilerExtractorFiles/AllSpecies.All_Intervals_Counts.perPopulation.SAMPLESIZECORRECTED.MULTINOMIALDOWNSAMPLED.COUNTS.CentralAsRevCompedToTs.SIGPROFILEREXTRACTORFORMAT.txt"),row.names=F,quote=F,sep="\t")

# this is the non-downsampled version *but is sample size corrected*
#all3merSpectra_processed_wider_notDownsampled <- pivot_wider(select(all3merSpectra_processed,c(`Mutation Types`,label,total_mutations_projected)),id_cols = `Mutation Types`,names_from =  label,values_from = total_mutations_projected)
#head(all3merSpectra_processed_wider_notDownsampled)                    

#write.table(all3merSpectra_processed_wider_notDownsampled,paste0(plotdir,"/SigProfilerExtractorFiles/AllSpecies.All_Intervals_Counts.perPopulation.SAMPLESIZECORRECTED.notdownsampled.COUNTS.CentralAsRevCompedToTs.SIGPROFILEREXTRACTORFORMAT.txt"),row.names=F,quote=F,sep="\t")

######## SEPARATE DISTANCES BY CENTRAL BP #####################
combo_function_separateByMutationType <- function(spectrumdf,variable,ilr_or_clr_or_none,speciesToInclude,speciesCodes,phyloDistances,centralMutationType){
  distance_dataframe <- 
    processSpectra(spectrumdf) %>%
    filter(mutation_1mer==centralMutationType) %>%
    pivotSpectra_perVariable(.,variable)  %>%
    select(c(variable,mutation_label,speciesToInclude)) %>% # just get the cols for the species you want 
    clr_or_ilr_orNoTransform_calculations(.,ilr_or_clr_or_none) %>%
    euclideanDistance(.) %>%
    addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
    addPhyloDistances_excludeNas(.,phyloDistances) %>%
    select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
    mutate(variable=variable,transform_label=ilr_or_clr_or_none,centralMutationType=centralMutationType)
  
  
  return(distance_dataframe)
  
  
}
######## mantel test faceted by central bp: ###############
combo_function_separateByMutationType_MantelTest <- function(spectrum, phylo_dist,variable,speciesToInclude,centralMutationType,mantelTestPermutationCount){
  spectrum_dist <- spectrum %>% 
    processSpectra() %>% 
    filter(mutation_1mer==centralMutationType) %>%
    merge(.,speciesCodes,by.x="label",by.y="code",suffixes=c(".spectrum",".sppCodes")) %>% 
    filter(label %in% speciesToInclude) %>%
    pivotSpectra_perVariable_KEEPFULLSPECIESNAMES(variable) %>%
    clr_or_ilr_orNoTransform_calculations(.,"clr") %>% # fixing as clr -- could set as an option in the function tho
    euclideanDistance_keepMatrix_UPPERANDLOWER(.)  
  # the matrices need to be ordered in the same way 
  # reorder the time matrix in the same way as the spectrum matrix:
  phylo_dist_reordered <- as.matrix(phylo_dist)[rownames(as.matrix(spectrum_dist)),colnames(as.matrix(spectrum_dist))]
  head(phylo_dist_reordered)
  print(paste0("carrying out mantel test with",mantelTestPermutationCount," permutations"))
  
  mtr = vegan::mantel(xdis=phylo_dist_reordered,ydis=spectrum_dist,permutations = mantelTestPermutationCount,method = "pearson")
  # make a data frame
  mantelTestResult_df <- data.frame(call="vegan_mantel",variable=variable,method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif),centralMutationType=centralMutationType)
  # mantel p value is calc as (ngreater+1)/(length(rsample)+1 so if there are none that are as well correlated then you'd get 0+1 / 999+1 = 1/ 10000 <- 0.001 which is what I get. which means in 999 permutations there are none that are as correlated. am I doing this right? 
  
  return(mantelTestResult_df)
  
}

#centralMutationTypes=c("C.T","C.A","C.G","A.C","A.T","A.G")
# for now haven't been separating out CpGs. 
centralMutationTypes=unique(all1merSpectra$mutation_1mer) # forn ow doesn't sep CpGs out

##### separate 3mer, 5mer and 7mer distances by central bp: ##########
all_distances_phylo_downsampMutRate_clr_SEPBYCENTRALBP = data.frame()
all_mantelresults_phylo_downsampMutRate_clr_SEPBYCENTRALBP = data.frame()
for(centralMutationType in centralMutationTypes){
  
  list_of_distance_dfs_phylo_SEPBYCENTRALBP <- lapply(listOfDFs[c("spectra_3mer","spectra_5mer","spectra_7mer")]
,combo_function_separateByMutationType, variable ="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon", ilr_or_clr_or_none = "clr",speciesToInclude = speciesToInclude_phylo,speciesCodes = speciesCodes,phyloDistances = phyloDistances,centralMutationType=centralMutationType)
  
  distance_dfs_phylo_sepByCentralBP <- bind_rows(list_of_distance_dfs_phylo_SEPBYCENTRALBP, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
  
  all_distances_phylo_downsampMutRate_clr_SEPBYCENTRALBP <- bind_rows(all_distances_phylo_downsampMutRate_clr_SEPBYCENTRALBP,distance_dfs_phylo_sepByCentralBP)
  
  rm(list_of_distance_dfs_phylo_SEPBYCENTRALBP)
  
  # carry out mantel test:
  list_of_mantel_dfs_phylo_SEPBYCENTRALBP <- lapply(listOfDFs[c("spectra_3mer","spectra_5mer","spectra_7mer")],combo_function_separateByMutationType_MantelTest, variable ="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon", speciesToInclude = speciesToInclude_phylo,phylo_dist = raxml_cophenetic_dist,centralMutationType=centralMutationType,mantelTestPermutationCount=mantelTestPermutationCount)
  
  mantel_dfs_phylo_sepByCentralBP <- bind_rows(list_of_mantel_dfs_phylo_SEPBYCENTRALBP, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
  
  all_mantelresults_phylo_downsampMutRate_clr_SEPBYCENTRALBP <- bind_rows(all_mantelresults_phylo_downsampMutRate_clr_SEPBYCENTRALBP,mantel_dfs_phylo_sepByCentralBP)
  
  rm(list_of_mantel_dfs_phylo_SEPBYCENTRALBP)
  
}

allDistances_facetedByCentralBP_plot <- ggplot(all_distances_phylo_downsampMutRate_clr_SEPBYCENTRALBP,aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point()+
  facet_grid(id~centralMutationType,scales="free_y")+
  theme_bw()+
  ggtitle("mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon -- clr")+
  geom_text(data=all_mantelresults_phylo_downsampMutRate_clr_SEPBYCENTRALBP,aes(x=1,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",round(statistic,2),"\np-value: ",as.character(signif)," (Mantel Test)")),size=2)

allDistances_facetedByCentralBP_plot

ggsave(paste0(manteldir,"FacetedByCentralBP.downsampmutrate.clr.png"),allDistances_facetedByCentralBP_plot,height=8,width=16)

all_mantelresults_phylo_downsampMutRate_clr_SEPBYCENTRALBP$negLog10Pvalue <- -log10(all_mantelresults_phylo_downsampMutRate_clr_SEPBYCENTRALBP$signif)

all_mantelresults_phylo_downsampMutRate_clr_SEPBYCENTRALBP_PlotPearsonR <- ggplot(all_mantelresults_phylo_downsampMutRate_clr_SEPBYCENTRALBP,aes(x=centralMutationType,y=statistic,fill=id))+
  geom_col(position="dodge")+
  theme_bw()+
  ggtitle("correlation between spectrum distance (downsampmutrate; clr) and cophenetic distance ")+
  ylab("Pearson's r")
all_mantelresults_phylo_downsampMutRate_clr_SEPBYCENTRALBP_PlotPearsonR

ggsave(paste0(manteldir,"FacetedByCentralBP.downsampmutrate.clr.mantelTestPearsonValues.png"),all_mantelresults_phylo_downsampMutRate_clr_SEPBYCENTRALBP_PlotPearsonR,height=8,width=16)


all_mantelresults_phylo_downsampMutRate_clr_SEPBYCENTRALBP_Plotlog10pvalue <- ggplot(all_mantelresults_phylo_downsampMutRate_clr_SEPBYCENTRALBP,aes(x=centralMutationType,y=negLog10Pvalue,fill=id))+
  geom_col(position="dodge")+
  theme_bw()+
  ggtitle("correlation between spectrum distance (downsampmutrate; clr) and cophenetic distance\ndashed line is 0.05/18 to correct for multiple testing (6 mutation types * 3 types of 3mer)")+
  ylab("log10 p-value")+
  geom_hline(yintercept = -log10(0.05/(6*3)),linetype="dashed")

all_mantelresults_phylo_downsampMutRate_clr_SEPBYCENTRALBP_Plotlog10pvalue

ggsave(paste0(manteldir,"FacetedByCentralBP.downsampmutrate.clr.mantelTestPValues.png"),all_mantelresults_phylo_downsampMutRate_clr_SEPBYCENTRALBP_Plotlog10pvalue,height=8,width=16)


######### BGC categories ##########
BGC_WS = c("A.C","A.G")
BGC_conserved = c("A.T","C.G")
BGC_SW = c("C.A","C.T")
BGC_Categories = list(BGC_WS=BGC_WS,BGC_conserved=BGC_conserved,BGC_SW=BGC_SW)

combo_function_separateByBGCCategory <- function(spectrumdf,variable,ilr_or_clr_or_none,speciesToInclude,speciesCodes,phyloDistances,BGC_Categories,BGC_Category){
  BGC_category_mutationTypes=BGC_Categories[[BGC_Category]]
  distance_dataframe <- 
    processSpectra(spectrumdf) %>%
    filter(mutation_1mer %in% BGC_category_mutationTypes) %>%
    pivotSpectra_perVariable(.,variable)  %>%
    select(c(variable,mutation_label,speciesToInclude)) %>% # just get the cols for the species you want 
    clr_or_ilr_orNoTransform_calculations(.,ilr_or_clr_or_none) %>%
    euclideanDistance(.) %>%
    addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
    addPhyloDistances_excludeNas(.,phyloDistances) %>%
    select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
    mutate(variable=variable,transform_label=ilr_or_clr_or_none,BGC_Category=BGC_Category)
  
  
  return(distance_dataframe)
  
  
}

######## mantel test faceted by central bp: ###############
combo_function_separateByBGCCategory_MantelTest <- function(spectrum, phylo_dist,variable,speciesToInclude,BGC_Categories,BGC_Category,mantelTestPermutationCount){
  BGC_category_mutationTypes=BGC_Categories[[BGC_Category]]
  
  spectrum_dist <- spectrum %>% 
    processSpectra() %>% 
    filter(mutation_1mer %in% BGC_category_mutationTypes) %>%
    merge(.,speciesCodes,by.x="label",by.y="code",suffixes=c(".spectrum",".sppCodes")) %>% 
    filter(label %in% speciesToInclude) %>%
    pivotSpectra_perVariable_KEEPFULLSPECIESNAMES(variable) %>%
    clr_or_ilr_orNoTransform_calculations(.,"clr") %>% # fixing as clr -- could set as an option in the function tho
    euclideanDistance_keepMatrix_UPPERANDLOWER(.)  
  # the matrices need to be ordered in the same way 
  # reorder the time matrix in the same way as the spectrum matrix:
  phylo_dist_reordered <- as.matrix(phylo_dist)[rownames(as.matrix(spectrum_dist)),colnames(as.matrix(spectrum_dist))]
  head(phylo_dist_reordered)
  print(paste0("carrying out mantel test with",mantelTestPermutationCount," permutations"))
  
  mtr = vegan::mantel(xdis=phylo_dist_reordered,ydis=spectrum_dist,permutations = mantelTestPermutationCount,method = "pearson")
  # make a data frame
  mantelTestResult_df <- data.frame(call="vegan_mantel",variable=variable,method=as.character(mtr$method),statistic=as.numeric(mtr$statistic),permCount=as.numeric(mtr$permutations),signif=as.numeric(mtr$signif),centralMutationType=centralMutationType,BGC_Category=BGC_Category)
  # mantel p value is calc as (ngreater+1)/(length(rsample)+1 so if there are none that are as well correlated then you'd get 0+1 / 999+1 = 1/ 10000 <- 0.001 which is what I get. which means in 999 permutations there are none that are as correlated. am I doing this right? 
  
  return(mantelTestResult_df)
  
}

##### separate 3mer, 5mer and 7mer distances by bgc category: ##########
all_distances_phylo_downsampMutRate_clr_SEPBYBGC = data.frame()
all_mantelresults_phylo_downsampMutRate_clr_SEPBYBGC = data.frame()
for(BGC_Category in c("BGC_WS","BGC_conserved","BGC_SW")){
  print(BGC_Category)  
  list_of_distance_dfs_phylo_SEPBYBGC<- lapply(listOfDFs[c("spectra_3mer","spectra_5mer","spectra_7mer")],combo_function_separateByBGCCategory, variable ="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon", ilr_or_clr_or_none = "clr",speciesToInclude = speciesToInclude_phylo,speciesCodes = speciesCodes,phyloDistances = phyloDistances,BGC_Categories,BGC_Category)
  
  distance_dfs_phylo_sepByBGC <- bind_rows(list_of_distance_dfs_phylo_SEPBYBGC, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
  
  all_distances_phylo_downsampMutRate_clr_SEPBYBGC <- bind_rows(all_distances_phylo_downsampMutRate_clr_SEPBYBGC,distance_dfs_phylo_sepByBGC)
  
  rm(list_of_distance_dfs_phylo_SEPBYBGC)
  
  # carry out mantel test:
  list_of_mantel_dfs_phylo_SEPBYBGC <- lapply(listOfDFs[c("spectra_3mer","spectra_5mer","spectra_7mer")],combo_function_separateByBGCCategory_MantelTest, variable ="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon", speciesToInclude = speciesToInclude_phylo,phylo_dist = raxml_cophenetic_dist,BGC_Categories=BGC_Categories,BGC_Category=BGC_Category,mantelTestPermutationCount=mantelTestPermutationCount)
  
  mantel_dfs_phylo_sepByBGC <- bind_rows(list_of_mantel_dfs_phylo_SEPBYBGC, .id = "id") # this keeps the 'name' from the list of dfs as the kmer label! nice!! 
  
  all_mantelresults_phylo_downsampMutRate_clr_SEPBYBGC <- bind_rows(all_mantelresults_phylo_downsampMutRate_clr_SEPBYBGC,mantel_dfs_phylo_sepByBGC)
  
  rm(list_of_mantel_dfs_phylo_SEPBYBGC)
  
}

allDistances_facetedByBGC <- ggplot(all_distances_phylo_downsampMutRate_clr_SEPBYBGC,aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point()+
  facet_grid(id~BGC_Category,scales="free_y")+
  theme_bw()+
  ggtitle("mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon -- clr")+
  geom_text(data=all_mantelresults_phylo_downsampMutRate_clr_SEPBYBGC,aes(x=1,y=-Inf,hjust=.8,vjust=-.5,label=paste0("Pearson's r:  ",round(statistic,2),"\np-value: ",as.character(signif)," (Mantel Test)")),size=2)

allDistances_facetedByBGC

ggsave(paste0(manteldir,"FacetedByBGC.downsampmutrate.clr.png"),allDistances_facetedByBGC,height=8,width=12)

all_mantelresults_phylo_downsampMutRate_clr_SEPBYBGC$negLog10Pvalue <- -log10(all_mantelresults_phylo_downsampMutRate_clr_SEPBYBGC$signif)

all_mantelresults_phylo_downsampMutRate_clr_SEPBYBGC_PlotPearsonR <- ggplot(all_mantelresults_phylo_downsampMutRate_clr_SEPBYBGC,aes(x=BGC_Category,y=statistic,fill=id))+
  geom_col(position="dodge")+
  theme_bw()+
  ggtitle("correlation between spectrum distance (downsampmutrate; clr) and cophenetic distance ")+
  ylab("Pearson's r")
all_mantelresults_phylo_downsampMutRate_clr_SEPBYBGC_PlotPearsonR

ggsave(paste0(manteldir,"FacetedByBGC.downsampmutrate.clr.mantelTestPearsonValues.png"),all_mantelresults_phylo_downsampMutRate_clr_SEPBYBGC_PlotPearsonR,height=5,width=8)


all_mantelresults_phylo_downsampMutRate_clr_SEPBYBGC_Plotlog10pvalue <- ggplot(all_mantelresults_phylo_downsampMutRate_clr_SEPBYBGC,aes(x=BGC_Category,y=negLog10Pvalue,fill=id))+
  geom_col(position="dodge")+
  theme_bw()+
  ggtitle("correlation between spectrum distance (downsampmutrate; clr) and cophenetic distance\ndashed line is 0.05/9 to correct for multiple testing (3 mutation types * 3 types of 3mer)")+
  ylab("log10 p-value")+
  geom_hline(yintercept = -log10(0.05/(3*3)),linetype="dashed")

all_mantelresults_phylo_downsampMutRate_clr_SEPBYBGC_Plotlog10pvalue

ggsave(paste0(manteldir,"FacetedByBGC.downsampmutrate.clr.mantelTestPValues.png"),all_mantelresults_phylo_downsampMutRate_clr_SEPBYBGC_Plotlog10pvalue,height=5,width=8)


