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
########### note: dplyr is a bit dangerous when dividing by sum(X) -- if the df has ever been grouped by things, it maintains that grouping. This manifests as sneaky bugs for fraction of segregating sites.
# Be very cautious ! !!! need to go back and check some resutls 
epsilon=1 # amount to add to 0 bins
separateCpGs="no" # yes or no
todaysdate=format(Sys.Date(),"%Y%m%d")
flagOfAnalysis=paste0(todaysdate,".plots.projectedCounts.epsilon.",epsilon,".sepCpGs.",separateCpGs,".addingApes") # a label that will be in title of plots and outdir 
plotdir=paste0("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/comparing_spectrum_and_phylo_divergence_WITHPROJECTIONS_CpGSitesSeparatedOut/",flagOfAnalysis,"/")
dir.create(plotdir,showWarnings = F)
########## table of target locations ########
# restricting to just species with >1 individual:
#tableOfTargetsToInclude=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/spectra/tableOfSpectraToIncludeInComparisons.SpeciesWithMultipleIndividualsOnly.20211026.txt",header=T,sep="\t") # this also has info on spectra but don't use that. 
# adding apes:
tableOfTargetsToInclude=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/spectra/tableOfSpectraToIncludeInComparisons.SpeciesWithMultipleIndividualsOnly.PlusApes.20211119.txt",header=T,sep="\t") 
############## all projected spectra #############
# okay have all the spectra together already in:
# THESE have been projected down to k haploids:
#allProjectedSpectra <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/allKSFSes_and_projected_Spectra/All_SpectrumCounts_ProjectedDownTo.16.Haploids.txt",header=T)
allProjectedSpectra <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/allKSFSes_and_projected_Spectra/All_SpectrumCounts_ProjectedDownTo.10.Haploids.includesAPes.txt",header=T) # where is vaquita! hmmm weird.
# this now includes bears_EUR and all apes (projectd down to 5 diploids, 10 haps)
speciesCodes=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/scripts/DNA_Shape_Project/analyses/PhylogeneticDistances/speciesCodes.txt",header=T) # for now just using AFR as human ; added apes
head(speciesCodes) # for getting phylo distances 

###### read in phylo distances:  #########
phyloDistances = read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/raxml_tree/raxml.cophenetic.pairwisedistances.fromUphametal.MYSPECIES.usethis.txt",header=T) # added apes
head(phyloDistances) # not all spp are present


########## read in hamming distances:  ##########
hammingDistances = read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/hamming_distance/ALLSPECIESCOMBINED.totalAlleleCounts.HammingDistance.allIntervals.WITHPOPINFO.AVERAGEDPERCOMPARISONTYPE.DIVIDEDBY2xGenomeSIZE.usethis.txt",header=T)
### STILL NEED APES ^^ (actually don't need apes. would need to merge vcfs to get hamming distances and in fact have their distances from phylogeny. only need hamming ot comapre popualtions)

##
######### read in time-calibrated tree cophenetic distances ########3
timeDistances = read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/time_calibrated_trees/20211122_addingApes/All_Species/tip_dated/cophenetic.Pairwise.BranchLengthSums.AveragedOverTrees.avg.sd.txt",header=T) # figure out how I want to incorproate these. 
# am labeling the avg branch len as cophenetic distance (which it is) so that it works with existing code
timeDistances$cophenetic_distance <- timeDistances$avgBranchLen # thisis avg cophenetic distance of tiem cal tree (across 1000 trees)
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
########## YOU ARE HERE -- deal with cpg targets ##########
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
processSpectra <- function(spectradf){
  # get mutation rates and fraction of segregating sites: also adding a column in where I add +epsilon (1) to every count to deal with 0 entries ; taking out adding epsilon to target counts since I fixed the empty target bug
  spectradf <- spectradf %>%
    mutate(total_mutations_projected_plusEpsilon = (total_mutations_projected+epsilon),mutation_rate_projected = (total_mutations_projected/total_target_count), mutation_rate_projected_plusEpsilon=(total_mutations_projected_plusEpsilon/total_target_count)) %>%
    group_by(species,population,label) %>%
    mutate(mutation_rate_projected_normalized = (mutation_rate_projected/sum(mutation_rate_projected)),fractionOfSegregatingSites_projected=(total_mutations_projected/sum(total_mutations_projected)), mutation_rate_projected_normalized_plusEpsilon=(mutation_rate_projected_plusEpsilon/sum(mutation_rate_projected_plusEpsilon)),fractionOfSegregatingSites_projected_plusEpsilon=(total_mutations_projected_plusEpsilon/sum(total_mutations_projected_plusEpsilon))) # adding plus1 versions.
  
  return(spectradf)
  
}
# test function:
test = processSpectra(all1merSpectra)
head(test)
# check that sums up add up to 1 (if not would indicate a grouping error in dplyr)
test %>% group_by(species,label,population) %>%
  summarise(checkSum = sum(fractionOfSegregatingSites_projected_plusEpsilon)==1)
########## need to do sample size correction/projection before further processing ############
# want to add one before or after sample size correction? come back to this. 
#all7merSpectra_proc <- processSpectra(all7merSpectra)
#head(all7merSpectra_proc)
# need to pivot it wider? 
pivotSpectra_perVariable <- function(processedspectradf,variableToPivot) {
  processedspectradf_wider <-  pivot_wider(processedspectradf,id_cols=c(mutation_label),names_from=c("label"),values_from =all_of(variableToPivot ))
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

speciesToInclude_phylo=c("mice_Mmd","mice_Ms","bears_ABC","bears_PB","vaquita","fin_whale_GOC","humans_AFR","apes_Gorilla","apes_Pan_paniscus","apes_Pan_troglodytes","apes_Pongo_abelii","apes_Pongo_pygmaeus") # human_AFR# switching to GOC fin whale for now due to ksfs weirdness
# adding in apes 
#for hamming distances:
speciesToInclude_hamming = c("bears_ABC"  , "bears_PB","bears_EUR", "fin_whale_ENP" ,"fin_whale_GOC", "humans_AFR"   , "humans_AMR" ,"humans_EAS"  ,  "humans_EUR"  ,  "humans_SAS" ,   "mice_Mmc"      ,"mice_Mmd"   ,   "mice_Mmm" ,   "mice_Ms") # excluding bears_EUR because would be too low a projection value



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


# do some loops
# all Distances
listOfDFs = list(spectra_1mer=all1merSpectra,spectra_3mer=all3merSpectra,spectra_5mer=all5merSpectra,spectra_7mer = all7merSpectra)
# lapply applies combo_function_phylo to each of the above spectra ^^ 
# if you set the names of a list, they get preserved in the output of lapply! cool!

listOfVars = c("mutation_rate_projected_normalized_plusEpsilon","fractionOfSegregatingSites_projected_plusEpsilon")
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


################ phylo plots #################
##### be careful that you aren't falsely combining things 

phylo_plot1 <- ggplot(all_distances_all_conditions_phylo,aes(x=cophenetic_distance,y=spectrum_distance,color=variable))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=1.8)+
  ggtitle(paste0("phylo distances\n",flagOfAnalysis))

phylo_plot1
ggsave(paste0(plotdir,"phylo_plot1_transforms.allconditions.png"),height=6,width=16,phylo_plot1)

# plot each individual statistic:
phylo_plot1b <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_normalized_plusEpsilon",],aes(x=cophenetic_distance,y=spectrum_distance))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=.8)+
  ggtitle(paste0("phylo distances\nmutation_rate_projected_normalized_plusEpsilon\n",flagOfAnalysis))+
  geom_smooth(method="lm")

phylo_plot1b
ggsave(paste0(plotdir,"phylo_plot1b_transforms.mutationrate.png"),height=6,width=16,phylo_plot1b)


# plot each individual statistic:
phylo_plot1c <- ggplot(all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="fractionOfSegregatingSites_projected_plusEpsilon",],aes(x=cophenetic_distance,y=spectrum_distance))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=.8)+
  ggtitle(paste0("note: apes likely outliers due to very few sites relative to other species (a lot of masking, so diffs not corrected for genome size may be more different?)\nphylo distances\nfractionOfSegregatingSites_projected_plusEpsilon\n",flagOfAnalysis))+
  geom_smooth(method="lm")

phylo_plot1c
ggsave(paste0(plotdir,"phylo_plot1c_transforms.fracsegsites.png"),height=6,width=16,phylo_plot1c)



############## hamming plots ###############
hamming_plot1 <- ggplot(all_distances_all_conditions_hamming,aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,color=variable))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  theme_bw()+
  geom_text_repel(aes(label=comparisonLabel),size=1.8)+
  ggtitle(paste0("hamming distances\n",flagOfAnalysis))

hamming_plot1
ggsave(paste0(plotdir,"hamming_plot1_transforms.allconditions.png"),height=6,width=16,hamming_plot1)
# note that in hamming plot, looks like when transformed differences betwn frac seg sites and mut rate disappear for genomes with same targets (eg PB adn ABC or human1 and human2 -- the ones in the hamming plot) ; must be some algebra that cancels out targets in some way. does'nt happen in phylo plot bc are mapped to diff genomes

# log 10 scale hamming plot:
hamming_plot1b <- hamming_plot1 + scale_x_log10() + xlab("hamming distance (log10 scaled)")
hamming_plot1b
ggsave(paste0(plotdir,"hamming_plot1b_log10_transforms.allconditions.png"),height=6,width=16,hamming_plot1b)

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
hamming_phylo_combo_plot1c <- ggplot(all_distances_all_conditions_hamming[all_distances_all_conditions_hamming$variable=="mutation_rate_projected_normalized_plusEpsilon",],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,shape="hamming"))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  geom_point(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="mutation_rate_projected_normalized_plusEpsilon",],aes(x=cophenetic_distance,y=spectrum_distance,shape="phylo"))+
  theme_bw()+
  #geom_text_repel(aes(label=comparisonLabel),size=1.8)+
  ggtitle(paste0("hamming distances and phylo distances \nmutation_rate_projected_normalized_plusEpsilon\n",flagOfAnalysis))+
  xlab("hamming or phylo distance")
hamming_phylo_combo_plot1c
ggsave(paste0(plotdir,"hamming_phylo_combo_plot1_transforms.mutationrate.png"),height=6,width=16,hamming_phylo_combo_plot1c)

# just plot frac seg sites: 
hamming_phylo_combo_plot1d <- ggplot(all_distances_all_conditions_hamming[all_distances_all_conditions_hamming$variable=="fractionOfSegregatingSites_projected_plusEpsilon",],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,shape="hamming"))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  geom_point(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$variable=="fractionOfSegregatingSites_projected_plusEpsilon",],aes(x=cophenetic_distance,y=spectrum_distance,shape="phylo"))+
  theme_bw()+
  #geom_text_repel(aes(label=comparisonLabel),size=1.8)+
  ggtitle(paste0("hamming distances and phylo distances \nfractionOfSegregatingSites_projected_plusEpsilon\n",flagOfAnalysis))+
  xlab("hamming or phylo distance")
hamming_phylo_combo_plot1d
ggsave(paste0(plotdir,"hamming_phylo_combo_plot1_transforms.fracsegsites.png"),height=6,width=16,hamming_phylo_combo_plot1d
)


############ humans and mice only ##########
humansMice=c("humans_AFR","humans_EUR","humans_EAS","humans_SAS","humans_AMR","mice_Mmd","mice_Ms","mouse_Ms","mouse_Mmd","mice_Mmm","mice_Mmc")

hamming_phylo_combo_plot2 <- ggplot(all_distances_all_conditions_hamming[all_distances_all_conditions_hamming$item1 %in% humansMice & all_distances_all_conditions_hamming$item2 %in% humansMice,],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,color=variable,shape="hamming"))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  geom_text_repel(aes(label=comparisonLabel))+
  geom_point(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$item1 %in% humansMice & all_distances_all_conditions_phylo$item2 %in% humansMice,],aes(x=cophenetic_distance,y=spectrum_distance,color=variable,shape="phylo"))+
  theme_bw()+
  #geom_text_repel(aes(label=comparisonLabel),size=1.8)+
  ggtitle(paste0("hamming distances and phylo distances \n",flagOfAnalysis))+
  xlab("hamming or phylo distance")+
  geom_text_repel(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$item1 %in% humansMice & all_distances_all_conditions_phylo$item2 %in% humansMice,],aes(x=cophenetic_distance,y=spectrum_distance,label=comparisonLabel_common_alphabetical))
hamming_phylo_combo_plot2
ggsave(paste0(plotdir,"hamming_phylo_combo_plot2_humans_mice_only_transforms.allconditions.png"),height=6,width=16,hamming_phylo_combo_plot2)

########## humans mice apes only  ##############
humansMiceApes=c("humans_AFR","humans_EUR","humans_EAS","humans_SAS","humans_AMR","mice_Mmd","mice_Ms","mouse_Ms","mouse_Mmd","mice_Mmm","mice_Mmc","apes_Pan_paniscus","apes_Pan_troglodytes","apes_Pongo_abelii","apes_Gorilla","apes_Pongo_pygmaeus")

hamming_phylo_combo_plot3 <- ggplot(all_distances_all_conditions_hamming[all_distances_all_conditions_hamming$item1 %in% humansMiceApes & all_distances_all_conditions_hamming$item2 %in% humansMiceApes,],aes(x=average_HammingDistance_divBy2xGenomeSize,y=spectrum_distance,color=variable,shape="hamming"))+
  facet_wrap(~transform_label~id,scales="free",ncol=4)+
  geom_point()+
  geom_text_repel(aes(label=comparisonLabel))+
  geom_point(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$item1 %in% humansMiceApes & all_distances_all_conditions_phylo$item2 %in% humansMiceApes,],aes(x=cophenetic_distance,y=spectrum_distance,color=variable,shape="phylo"))+
  theme_bw()+
  #geom_text_repel(aes(label=comparisonLabel),size=1.8)+
  ggtitle(paste0("hamming distances and phylo distances \n",flagOfAnalysis))+
  xlab("hamming or phylo distance")+
  geom_text_repel(data=all_distances_all_conditions_phylo[all_distances_all_conditions_phylo$item1 %in% humansMiceApes & all_distances_all_conditions_phylo$item2 %in% humansMiceApes,],aes(x=cophenetic_distance,y=spectrum_distance,label=comparisonLabel_common_alphabetical),size=0.8)
hamming_phylo_combo_plot3
ggsave(paste0(plotdir,"hamming_phylo_combo_plot3_humans_mice_apes_only_transforms.allconditions.png"),height=6,width=16,hamming_phylo_combo_plot3)



############# restrict to semi-independent comparisons ###################

# want 1 each: human-human, ape-ape, human-ape, mouse-mouse, human-mouse, whale-mouse, bear-bear, whale-whale ; don't want multiple:


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

listOfVars = c("mutation_rate_projected_normalized_plusEpsilon","fractionOfSegregatingSites_projected_plusEpsilon")
listOfTransforms=c("none","clr") # skipping ilr for now 
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


##################### YOU'VE MODIFIED UP TO HERE #############
#################### want to separate distances by central mutation type ##############

combo_function_separateByMutationType <- function(spectrumdf,variable,label,ilr_or_clr,centralMutationType){
  distance_dataframe <- 
    processSpectra(spectrumdf) %>%
    filter(mutation_1mer==centralMutationType) %>%
    pivotSpectra_perVariable(.,variable)  %>%
    select(c(variable,mutation_label,speciesToInclude_phylo)) %>% # just get the cols for the species you want 
    clr_or_ilr_calculations(.,ilr_or_clr) %>%
    euclideanDistance_ilr_or_clr(.) %>%
    addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
    addPhyloDistances_excludeNas(.,phyloDistances) %>%
    select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
    mutate(label=label,variable=variable,ilr_or_clr=ilr_or_clr,centralMutationType=centralMutationType)
  
  
  return(distance_dataframe)
  
  
}

centralMutationTypes=c("C.T","C.T_CpG","C.A","C.A_CpG","C.G","C.G_CpG","A.C","A.T","A.G")


########### 3mer separated by central type -- THI SIS COOL  ########


clr_distance_rescaledMutRate_3mer_perMutationType <- data.frame()
df <- data.frame()
for(centralMutationType in centralMutationTypes){
  df <- combo_function_separateByMutationType(all3merSpectra,"mutation_rate_projected_normalized_plusEpsilon","3mer_clr","clr",centralMutationType)
  clr_distance_rescaledMutRate_3mer_perMutationType <- bind_rows(clr_distance_rescaledMutRate_3mer_perMutationType,df)
}


perMutTypePlot_3mer <- ggplot(clr_distance_rescaledMutRate_3mer_perMutationType,aes(x=cophenetic_distance,y=spectrum_distance,color=centralMutationType))+
  #geom_line()+
  geom_point()+
  geom_smooth(se=F)+
  ggtitle(unique(clr_distance_rescaledMutRate_3mer_perMutationType$label))
perMutTypePlot_3mer

ggsave(paste0(plotdir,"3mer.clr.PerMutationType.epsilon.",epsilon,".png"),perMutTypePlot_3mer,height=8,width=12)

########### 5mer separated by central type -- THI SIS COOL  ########

clr_distance_rescaledMutRate_5mer_perMutationType <- data.frame()
df <- data.frame()
for(centralMutationType in centralMutationTypes){
  df <- combo_function_separateByMutationType(all5merSpectra,"mutation_rate_projected_normalized_plusEpsilon","5mer_clr","clr",centralMutationType)
  clr_distance_rescaledMutRate_5mer_perMutationType <- bind_rows(clr_distance_rescaledMutRate_5mer_perMutationType,df)
}

perMutTypePlot_5mer <- ggplot(clr_distance_rescaledMutRate_5mer_perMutationType,aes(x=cophenetic_distance,y=spectrum_distance,color=centralMutationType))+
  #geom_line()+
  geom_point()+
  geom_smooth(se=F)+
  ggtitle(unique(clr_distance_rescaledMutRate_5mer_perMutationType$label))
perMutTypePlot_5mer
ggsave(paste0(plotdir,"5mer.clr.PerMutationType.epsilon.",epsilon,".png"),perMutTypePlot_5mer,height=8,width=12)


########### 7mer separated by central type -- THI SIS COOL  ########


clr_distance_rescaledMutRate_7mer_perMutationType <- data.frame()
df <- data.frame()
for(centralMutationType in centralMutationTypes){
  df <- combo_function_separateByMutationType(all7merSpectra,"mutation_rate_projected_normalized_plusEpsilon","7mer_clr","clr",centralMutationType)
  clr_distance_rescaledMutRate_7mer_perMutationType <- bind_rows(clr_distance_rescaledMutRate_7mer_perMutationType,df)
}

perMutTypePlot_7mer <- ggplot(clr_distance_rescaledMutRate_7mer_perMutationType,aes(x=cophenetic_distance,y=spectrum_distance,color=centralMutationType))+
  #geom_line()+
  geom_point()+
  geom_smooth(se=F)+
  ggtitle(unique(clr_distance_rescaledMutRate_7mer_perMutationType$label))
perMutTypePlot_7mer
ggsave(paste0(plotdir,"7mer.clr.PerMutationType.epsilon.",epsilon,".png"),perMutTypePlot_7mer,height=8,width=12)

########### explore 1mer spectra ##########

all1mer_processed <- processSpectra(all1merSpectra)
all1mer_processed_pivot <- pivotSpectra_perVariable(all1mer_processed,"mutation_rate_projected_normalized_plusEpsilon")
all7mer_processed <- processSpectra(all7merSpectra)
head(all1mer_processed)
spectrum1merPlot1 <- ggplot(all1mer_processed,aes(x=mutation_1mer,y=mutation_rate_projected_normalized_plusEpsilon,fill=label))+
  geom_col(position="dodge")
spectrum1merPlot1
ggsave(paste0(plotdir,"1mer.mutationRate.PerMutationType.epsilon.",epsilon,".png"),spectrum1merPlot1,height=8,width=12)

spectrum1merPlot2 <- ggplot(all1mer_processed,aes(x=mutation_1mer,y=fractionOfSegregatingSites_projected_plusEpsilon,fill=label))+
  geom_col(position="dodge")
spectrum1merPlot2
ggsave(paste0(plotdir,"1mer.fracSegSites.PerMutationType.epsilon.",epsilon,".png"),spectrum1merPlot1,height=8,width=12)

# play with dist a bit:
tidy(dist(t(select(all1mer_processed_pivot,c(vaquita,mice_Mmd)))))
tidy(dist(t(select(all1mer_processed_pivot,c(vaquita,fin_whale_ENP)))))
tidy(dist(t(select(all1mer_processed_pivot,c(vaquita,bears_ABC)))))
tidy(dist(t(select(all1mer_processed_pivot,c(vaquita,bears_PB)))))


head(all7mer_processed)
View(all7mer_processed[all7mer_processed$mutation_1mer=="C.T_CpG" & all7mer_processed$label=="bears_ABC",])
