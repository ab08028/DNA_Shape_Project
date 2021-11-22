############# gather all spectra and calculate rescaled mutation rates  ###########
require(reshape2)
require(dplyr)
require(ggplot2)
require(tidyr)
require(compositions)  # for clr
require(broom) # for tidy
require(ggrepel)

epsilon=1 # amount to add to 0 bins
plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/comparing_spectrum_and_phylo_divergence_WITHPROJECTIONS_CpGSitesSeparatedOut/"
dir.create(plotdir,showWarnings = F)

#tableOfSpectraToInclude=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/spectra/tableOfSpectraToIncludeInComparisons.20211013.txt",header=T,sep="\t")  # you can add to this if you want more includd ; currenlty human and fin whale missing because need targets and genome ide measures for fin whale
# restricting to just species with >1 individual:
tableOfTargetsToInclude=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/spectra/tableOfSpectraToIncludeInComparisons.SpeciesWithMultipleIndividualsOnly.20211026.txt",header=T,sep="\t") # this also has info on spectra but don't use that. 
# okay have all the spectra together already in:
# THESE have been projected down to k haploids:
allProjectedSpectra <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/allKSFSes_and_projected_Spectra/All_SpectrumCounts_ProjectedDownTo.16.Haploids.txt",header=T)

speciesCodes=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/scripts/DNA_Shape_Project/analyses/PhylogeneticDistances/speciesCodes.txt",header=T) # for now just using AFR as human
head(speciesCodes)

###### read in phylo distances: 
phyloDistances = read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/raxml_tree/raxml.cophenetic.pairwisedistances.fromUphametal.MYSPECIES.usethis.txt",header=T)
head(phyloDistances) # not all spp are present

# dleu2 has this one weird non7mer in its spectrum-- GAAAG.GCAAG
# that messes stuff up (why didn't that come up before?)
# okay have to work ith table 

# no longer have to deal with exceptions -- all spectra are in same format because are from projected ksfs


# get ancestral kmers:
allProjectedSpectra$ancestral7mer <- substr(allProjectedSpectra$variable,1,7)
allProjectedSpectra$ancestral5mer <- substr(allProjectedSpectra$variable,2,6)
allProjectedSpectra$ancestral3mer <- substr(allProjectedSpectra$variable,3,5)
allProjectedSpectra$ancestral1mer_CpGNotLabeled <- substr(allProjectedSpectra$variable,4,4) # not yet CpG identified 

# get central mutation types:
allProjectedSpectra$mutation_7mer <- allProjectedSpectra$variable
allProjectedSpectra$mutation_5mer <- paste0(substr(allProjectedSpectra$variable,2,6),".",substr(allProjectedSpectra$variable,10,14))
allProjectedSpectra$mutation_3mer <- paste0(substr(allProjectedSpectra$variable,3,5),".",substr(allProjectedSpectra$variable,11,13))

#### ADD IN CPG LABEL 
allProjectedSpectra$mutation_1mer_CpGNotLabeled <- paste0(substr(allProjectedSpectra$variable,4,4),".",substr(allProjectedSpectra$variable,12,12))
# need to add CpG Label 

#### label CpGs and specify them as ancestral target as well
allProjectedSpectra$CpGLabel <- ""
# find all CpG sites:
allProjectedSpectra[grepl("CG",allProjectedSpectra$ancestral3mer) & allProjectedSpectra$ancestral1mer_CpGNotLabeled=="C",]$CpGLabel <- "_CpG" 
allProjectedSpectra$mutation_1mer <- paste0(allProjectedSpectra$mutation_1mer_CpGNotLabeled,allProjectedSpectra$CpGLabel)
# label ancestral state with CpG as well
allProjectedSpectra$ancestral1mer <- allProjectedSpectra$ancestral1mer_CpGNotLabeled
allProjectedSpectra[allProjectedSpectra$CpGLabel=="_CpG",]$ancestral1mer <- "C_CpG"



allProjectedSpectra$mutation_1mer <- paste0(allProjectedSpectra$mutation_1mer_CpGNotLabeled,allProjectedSpectra$CpGLabel)


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
    mutate(target_5mer=substr(target_7mer,2,6),target_3mer=substr(target_7mer,3,5),target_1mer_CpGsNotLabeled= substr(target_7mer,4,4)) # get substrings  # need CpG targets 
  # need to label CpGs 
  targets$CpGLabel <- ""
  targets[grepl("CG",targets$target_3mer),]$CpGLabel <- "_CpG"
  targets$target_1mer <- paste0(targets$target_1mer_CpGsNotLabeled,targets$CpGLabel)
  
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
test = processSpectra(all7merSpectra)
head(test)

########## need to do sample size correction/projection before further processing ############
# want to add one before or after sample size correction? come back to this. 
#all7merSpectra_proc <- processSpectra(all7merSpectra)
#head(all7merSpectra_proc)
# need to pivot it wider? 
pivotSpectra_perVariable <- function(processedspectradf,variableToPivot) {
  processedspectradf_wider <-  pivot_wider(processedspectradf,id_cols=mutation_label,names_from=c("label"),values_from =all_of(variableToPivot ))
  processedspectradf_wider$variable <- variableToPivot
  # pivots on a variable 
  return(processedspectradf_wider)
}

# test function
testpivoted = pivotSpectra_perVariable(test,"mutation_rate_projected_normalized_plusEpsilon")
head(testpivoted)
# want to speed up with only a subset of species.
# ilr/clr calcs require no 0 entries
clr_or_ilr_calculations <- function(pivotedspectradf_plusEpsilonVersionOfVariables,ilr_or_clr){
  # remove character vectors:
  tableWithoutIDVars <- select(pivotedspectradf_plusEpsilonVersionOfVariables, -c("variable","mutation_label"))
  # AHA so if you run clr just on a matrix it goes ROW-WISE not column wise! annoying! so need to transpose
  #clr_matrix <- clr(t(tableWithoutIDVars)) # if you do it this way it matches when I do it just on on un-tranposed column by itself. if you do it without transposing then you get weird row-wise results that are WRONG. so be cautious here. 
  # so alternative without tranposing is to use sapply -- this is a lot nicer I think
  if(ilr_or_clr=="clr"){
  df <- data.frame(sapply(tableWithoutIDVars,clr)) # this way you don't have to transpose and you get a much nicer shaped output
  } else if(ilr_or_clr=="ilr"){
    df <- data.frame(sapply(tableWithoutIDVars,ilr)) # this way you don't have to transpose and you get a much nicer shaped output
    
  } else {
    print("invalid option")
    break
  }
  return(df)
  
}
#vars=c("fractionOfSegregatingSites","mutation_rate_normalized")
#all7merSpectra_pivot <- pivotSpectra_perVariable(all7merSpectra_proc,"fractionOfSegregatingSites")

# ilr_calculations <- function(pivotedspectradf_plusEpsilonVersionOfVariables){
#   # remove character vectors:
#   tableWithoutIDVars <- select(pivotedspectradf_plusEpsilonVersionOfVariables, -c("variable","mutation_label"))
#   # AHA so if you run clr just on a matrix it goes ROW-WISE not column wise! annoying! so need to transpose
#   #clr_matrix <- clr(t(tableWithoutIDVars)) # if you do it this way it matches when I do it just on on un-tranposed column by itself. if you do it without transposing then you get weird row-wise results that are WRONG. so be cautious here. 
#   # so alternative without tranposing is to use sapply -- this is a lot nicer I think
#   ilr_df <- data.frame(sapply(tableWithoutIDVars,ilr)) # this way you don't have to transpose and you get a much nicer shaped output
#   return(ilr_df)
#   
# }

# get distances
# for regular distance. 
euclideanDistance <- function(pivotedspectradf){
  # make it a matrix with each row as a species (transposed and all numeric)
  # variable in question:
  variable = unique(pivotedspectradf$variable)
  matrix <- t(select(pivotedspectradf,-c("mutation_label","variable")))
  distances <- tidy(dist(matrix))
  colnames(distances) <- c("item1","item2","spectrum_distance")
  distances$variable <- variable # label
  return(distances)
}

euclideanDistance_ilr_or_clr <- function(pivotedspectradf){
# make it a matrix with each row as a species (transposed and all numeric)
# variable in question:
  # distances are row-wise nto column wise so need to transpose 
  distances <- tidy(dist(t(pivotedspectradf))) # transpose for dist otherwise it goes along rows and is wrong
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
######### get euclidean distances for fraction of seg sites ########## :
# aha it's a problem that bears are also EUR! and humans are EUR! aha! that's funny. Okay need to deal with t
# updating to do clr? 
#DO CLR and ILR here# HERE: 
speciesToInclude_1=c("mice_Mmd","mice_Ms","bears_ABC","bears_PB","vaquita","fin_whale_ENP","humans_AFR") # human_AFR
# comparisons to include : 
# mouse-mouse, bear-bear, whale-whale, human-mouse, human-bear, human-whale 


###### Get ILR and CLR distances ########
# this is slow so restrict just to species you want
##### this subsets just to the species you want above ^^ 
############# clr mutation rate ##############

####### combine functions ##########
# input: spectrum, variable (frac seg sites or mut rate), label (1mer , 5mer etc, and choice of ilr or clr)
combo_function <- function(spectrumdf,variable,label,ilr_or_clr){
  distance_dataframe <- 
    processSpectra(spectrumdf) %>%
    pivotSpectra_perVariable(.,variable)  %>%
    select(c(variable,mutation_label,speciesToInclude_1)) %>% # just get the cols for the species you want 
    clr_or_ilr_calculations(.,ilr_or_clr) %>%
    euclideanDistance_ilr_or_clr(.) %>%
    addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
    addPhyloDistances_excludeNas(.,phyloDistances) %>%
    select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
    mutate(label=label,variable=variable,ilr_or_clr=ilr_or_clr)
  
  
  return(distance_dataframe)
  
  
}

# get all types and stack them:
######## MUTATION RATE ##########
####### CLR: ########
clr_distance_rescaledMutRate_1mer = combo_function(all1merSpectra,"mutation_rate_projected_normalized_plusEpsilon","1mer_clr","clr")

clr_distance_rescaledMutRate_3mer = combo_function(all3merSpectra,"mutation_rate_projected_normalized_plusEpsilon","3mer_clr","clr")

clr_distance_rescaledMutRate_5mer = combo_function(all5merSpectra,"mutation_rate_projected_normalized_plusEpsilon","5mer_clr","clr")

clr_distance_rescaledMutRate_7mer = combo_function(all7merSpectra,"mutation_rate_projected_normalized_plusEpsilon","7mer_clr","clr")


######### ILR  -- SLOW !!! ##############
ilr_distance_rescaledMutRate_1mer = combo_function(all1merSpectra,"mutation_rate_projected_normalized_plusEpsilon","1mer_ilr","ilr")

ilr_distance_rescaledMutRate_3mer = combo_function(all3merSpectra,"mutation_rate_projected_normalized_plusEpsilon","3mer_ilr","ilr")

ilr_distance_rescaledMutRate_5mer = combo_function(all5merSpectra,"mutation_rate_projected_normalized_plusEpsilon","5mer_ilr","ilr")

# TOO SLOW exhausts vector memory: ilr_distance_rescaledMutRate_7mer = combo_function(all7merSpectra,"mutation_rate_projected_normalized_plusEpsilon","7mer_ilr","ilr")

############# combine mut rate dists ################
all_clr_ilr_dists_mutrate <- bind_rows(ilr_distance_rescaledMutRate_1mer,ilr_distance_rescaledMutRate_3mer,ilr_distance_rescaledMutRate_5mer,clr_distance_rescaledMutRate_1mer,clr_distance_rescaledMutRate_3mer,clr_distance_rescaledMutRate_5mer,clr_distance_rescaledMutRate_7mer)

#all_clr_ilr_dists_mutrate <- bind_rows(clr_distance_rescaledMutRate_1mer,clr_distance_rescaledMutRate_3mer,clr_distance_rescaledMutRate_5mer,clr_distance_rescaledMutRate_7mer)

write.table(all_clr_ilr_dists_mutrate,paste0(plotdir,"ilr.clr.distances.RescaledMutationRate.epsilon.",as.character(epsilon),".PROJECTEDCOUNTS.txt"),row.names=F,quote=F,sep="\t")

all_clr_ilr_dists_mutrate$label <- factor(all_clr_ilr_dists_mutrate$label,levels=c("1mer_clr","3mer_clr","5mer_clr","7mer_clr","1mer_ilr","3mer_ilr","5mer_ilr","7mer_ilr"))
ilrclr_plot1 <- ggplot(all_clr_ilr_dists_mutrate,aes(x=cophenetic_distance,y=spectrum_distance))+
  facet_wrap(~label,scales="free",ncol=4)+
  geom_point()+
  geom_smooth()+
  ggtitle("mut rate")+
  geom_text(aes(label=comparisonLabel_common_alphabetical),size=1.8)
ilrclr_plot1
ggsave(paste0(plotdir,"ilr.clr.distances.mutationRate.plusEpsilon.",as.character(epsilon),".PROJECTEDCOUNTS.png"),ilrclr_plot1,width=13,height=7)

######## Frac seg sites ##########
####### CLR: ########
clr_distance_FracSegSites_1mer = combo_function(all1merSpectra,"fractionOfSegregatingSites_plusEpsilon","1mer_clr","clr")

clr_distance_FracSegSites_3mer = combo_function(all3merSpectra,"fractionOfSegregatingSites_plusEpsilon","3mer_clr","clr")

clr_distance_FracSegSites_5mer = combo_function(all5merSpectra,"fractionOfSegregatingSites_plusEpsilon","5mer_clr","clr")

clr_distance_FracSegSites_7mer = combo_function(all7merSpectra,"fractionOfSegregatingSites_plusEpsilon","7mer_clr","clr")


######### ILR  -- SLOW !!! ##############
ilr_distance_FracSegSites_1mer = combo_function(all1merSpectra,"fractionOfSegregatingSites_plusEpsilon","1mer_ilr","ilr")

ilr_distance_FracSegSites_3mer = combo_function(all3merSpectra,"fractionOfSegregatingSites_plusEpsilon","3mer_ilr","ilr")

ilr_distance_FracSegSites_5mer = combo_function(all5merSpectra,"fractionOfSegregatingSites_plusEpsilon","5mer_ilr","ilr")

# TOO SLOW exhausts vector memory: ilr_distance_FracSegSites_7mer = combo_function(all7merSpectra,"mutation_rate_projected_normalized_plusEpsilon","7mer_ilr","ilr")
####### combine and plot ######


all_clr_ilr_dists_fracSegSites <- bind_rows(ilr_distance_fracSegSites_1mer,ilr_distance_fracSegSites_3mer,ilr_distance_fracSegSites_5mer,clr_distance_fracSegSites_1mer,clr_distance_fracSegSites_3mer,clr_distance_fracSegSites_5mer,clr_distance_fracSegSites_7mer) # keeping out for now ilr_distance_fracSegSites_7mer


write.table(all_clr_ilr_dists_fracSegSites,paste0(plotdir,"ilr.clr.distances.FracSegSites.epsilon.",as.character(epsilon),".PROJECTEDCOUNTS.txt"),row.names=F,quote=F,sep="\t")

all_clr_ilr_dists_fracSegSites$label <- factor(all_clr_ilr_dists_fracSegSites$label,levels=c("1mer_clr","3mer_clr","5mer_clr","7mer_clr","1mer_ilr","3mer_ilr","5mer_ilr","7mer_ilr"))


ilrclr_plot2 <- ggplot(all_clr_ilr_dists_fracSegSites,aes(x=cophenetic_distance,y=spectrum_distance))+
  facet_wrap(~label,scales="free",ncol=4)+
  geom_point()+
  geom_smooth()+
  ggtitle("frac seg sites")
ilrclr_plot2
ggsave(paste0(plotdir,"ilr.clr.distances.fracSegSites.plusEpsilon.",as.character(epsilon),".PROJECTEDCOUNTS.png"),ilrclr_plot2,width=13,height=7)



############ regular stuff non clr ilr ##########
spectrum_and_phylo_dist_1mer_fracSegSites <- 
  processSpectra(all1merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites_projected_plusEpsilon")  %>%
  select(c(variable,mutation_label,all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="1mer_noTransform_projected_plusEpsilon")

spectrum_and_phylo_dist_3mer_fracSegSites <- 
  processSpectra(all3merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites_projected_plusEpsilon")  %>%
  select(c(variable,mutation_label,all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="3mer_noTransform_projected_plusEpsilon")

spectrum_and_phylo_dist_5mer_fracSegSites <- 
  processSpectra(all5merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites_projected_plusEpsilon")  %>%
  select(c(variable,mutation_label,all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="5mer_noTransform_projected_plusEpsilon")

spectrum_and_phylo_dist_7mer_fracSegSites <- 
  processSpectra(all7merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites_projected_plusEpsilon")  %>%
  select(c(variable,mutation_label,all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="7mer_noTransform_projected_plusEpsilon")

# bind them all together:
spectrum_and_phylo_dist_allKmers_fracSegSites <- bind_rows(spectrum_and_phylo_dist_1mer_fracSegSites,spectrum_and_phylo_dist_3mer_fracSegSites,spectrum_and_phylo_dist_5mer_fracSegSites,spectrum_and_phylo_dist_7mer_fracSegSites)


######### get euclidean distances for normalized mut rates (use the +epsilon version!) #######
spectrum_and_phylo_dist_1mer_normMutRate <- 
  processSpectra(all1merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_projected_normalized_plusEpsilon")  %>%
  select(c(variable,mutation_label,all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="1mer_noTransform_projected_plusEpsilon")

spectrum_and_phylo_dist_3mer_normMutRate <- 
  processSpectra(all3merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_projected_normalized_plusEpsilon")  %>%
  select(c(variable,mutation_label,all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="3mer_noTransform_projected_plusEpsilon")

spectrum_and_phylo_dist_5mer_normMutRate <- 
  processSpectra(all5merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_projected_normalized_plusEpsilon")  %>%
  select(c(variable,mutation_label,all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="5mer_noTransform_projected_plusEpsilon")

spectrum_and_phylo_dist_7mer_normMutRate <- 
  processSpectra(all7merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_projected_normalized_plusEpsilon")  %>%
  select(c(variable,mutation_label,all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="7mer_noTransform_projected_plusEpsilon")

spectrum_and_phylo_dist_allKmers_normMutRate <- bind_rows(spectrum_and_phylo_dist_1mer_normMutRate,spectrum_and_phylo_dist_3mer_normMutRate,spectrum_and_phylo_dist_5mer_normMutRate,spectrum_and_phylo_dist_7mer_normMutRate)


########## combine raw euclidean with transforms to plot ########

allDistanceTypes_mutRate <- bind_rows(all_clr_ilr_dists_mutrate,spectrum_and_phylo_dist_allKmers_normMutRate)

write.table(allDistanceTypes_mutRate,paste0(plotdir,"ALLDISTANCETYPES.includesNonTransformed.distances.NormMutRate.epsilon.",as.character(epsilon),".PROJECTEDCOUNTS.txt"),row.names=F,quote=F,sep="\t")


allDistanceTypes_fracSegSites <- bind_rows(all_clr_ilr_dists_fracSegSites,spectrum_and_phylo_dist_allKmers_fracSegSites)
write.table(allDistanceTypes_fracSegSites,paste0(plotdir,"ALLDISTANCETYPES.includesNonTransformed.distances.FracSegSites.epsilon.",as.character(epsilon),".PROJECTEDCOUNTS.txt"),row.names=F,quote=F,sep="\t")


allDistanceTypes_mutRate$label <- factor(allDistanceTypes_mutRate$label,  levels=c("1mer_noTransform_projected_plusEpsilon","3mer_noTransform_projected_plusEpsilon","5mer_noTransform_projected_plusEpsilon","7mer_noTransform_projected_plusEpsilon","1mer_clr","3mer_clr","5mer_clr","7mer_clr","1mer_ilr","3mer_ilr","5mer_ilr"))

allDistanceTypes_fracSegSites$label <- factor(allDistanceTypes_fracSegSites$label, levels=c("1mer_noTransform_projected_plusEpsilon","3mer_noTransform_projected_plusEpsilon","5mer_noTransform_projected_plusEpsilon","7mer_noTransform_projected_plusEpsilon","1mer_clr","3mer_clr","5mer_clr","7mer_clr","1mer_ilr","3mer_ilr","5mer_ilr"))

allDistanceTypesPlot1 <- ggplot(allDistanceTypes_mutRate,aes(x=cophenetic_distance,y=spectrum_distance))+
  facet_wrap(~label,scales="free_y",ncol=4)+
  geom_point()+
  ggtitle(paste0("mutation rate; projected; epsilon = ",epsilon))+
  theme_bw()+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=2)+
  theme(text=element_text(size=14))



allDistanceTypesPlot1
ggsave(paste0(plotdir,"allDistanceTypes.AllWithPlus1Transform.MutationRate.epsilon.",as.character(epsilon),".png"),allDistanceTypesPlot1,height=8,width=16)

allDistanceTypesPlot2 <- ggplot(allDistanceTypes_fracSegSites,aes(x=cophenetic_distance,y=spectrum_distance))+
  facet_wrap(~label,scales="free_y",ncol=4)+
  geom_point()+
  ggtitle(paste0("fraction of segregating sites; projected; epsilon = ",epsilon))+
  theme_bw()+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=2)+
  theme(text=element_text(size=14))

allDistanceTypesPlot2
ggsave(paste0(plotdir,"allDistanceTypes.AllWithPlus1Transform.FracSegSites.epsilon.",as.character(epsilon),".png"),allDistanceTypesPlot2,height=8,width=16)

######### let's make some more plots ###############

fracsegsites_plot1_withlabels <- ggplot(spectrum_and_phylo_dist_allKmers_fracSegSites,aes(x=cophenetic_distance,y=spectrum_distance,color=label))+
  geom_point()+
  #geom_smooth()+
  facet_wrap(~label,scales="free")+
  theme_bw()+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=3)+
  ggtitle("fraction of segregating sites")
fracsegsites_plot1_withlabels
ggsave(paste0(plotdir,"fracsegsitesplot1_withlabels.png"),fracsegsitesplot1_withlabels,height=15,width=20)


mutRateNorm_plot1_withlabels <- ggplot(spectrum_and_phylo_dist_allKmers_normMutRate,aes(x=cophenetic_distance,y=spectrum_distance,color=label))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~label,scales="free")+
  theme_bw()+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=3)
mutRateNorm_plot1_withlabels
ggsave(paste0(plotdir,"mutRateNorm_plot1_withlabels.png"),mutRateNorm_plot1_withlabels,height=15,width=20)

##### color plots by broad classificaiton type #######

fracsegsites_plot2_withlabels <- ggplot(spectrum_and_phylo_dist_allKmers_fracSegSites,aes(x=cophenetic_distance,y=spectrum_distance,color=label))+
  geom_point()+
  #geom_smooth()+
  facet_wrap(~label,scales="free")+
  theme_bw()+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=4)+
  ggtitle("fraction of segregating sites")
fracsegsites_plot1_withlabels
ggsave(paste0(plotdir,"fracsegsitesplot1_withlabels.png"),fracsegsitesplot1_withlabels,height=15,width=20)


mutRateNorm_plot2_withlabels <- ggplot(spectrum_and_phylo_dist_allKmers_normMutRate,aes(x=cophenetic_distance,y=spectrum_distance,color=comparisonLabel_broad_alphabetical))+
  geom_point()+
  facet_wrap(~label,scales="free")+
  theme_bw()+
  geom_text_repel(data=spectrum_and_phylo_dist_allKmers_normMutRate[spectrum_and_phylo_dist_allKmers_normMutRate$comparisonLabel_broad_alphabetical=="baleen_whale.toothed_whale",],aes(x=cophenetic_distance,y=spectrum_distance,label=comparisonLabel_common_alphabetical),size=1,color="black")
mutRateNorm_plot2_withlabels
ggsave(paste0(plotdir,"mutRateNorm_plot1_withlabels.colorByBroad.png"),mutRateNorm_plot2_withlabels,height=15,width=20)

########## try without CpG Sites?? ###############
###### simplify down to just species for which I have multi individual data ######
# species to include
speciesToInclude=c("mice_Mmd","mice_Ms","bears_ABC","bears_PB","vaquita") # human_AFR

mutRateNorm_plot3_subsetOfSpecies <- ggplot(spectrum_and_phylo_dist_allKmers_normMutRate[spectrum_and_phylo_dist_allKmers_normMutRate$item1 %in% speciesToInclude & spectrum_and_phylo_dist_allKmers_normMutRate$item2 %in% speciesToInclude, ],aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point()+
  facet_wrap(~label,scales="free")+
  theme_bw()+
  geom_text_repel(aes(x=cophenetic_distance,y=spectrum_distance,label=comparisonLabel_common_alphabetical),size=4,color="black")#+
  #geom_smooth()
mutRateNorm_plot3_subsetOfSpecies
ggsave(paste0(plotdir,"mutRateNorm_SUBSETOFSPECIES.png"),mutRateNorm_plot3_subsetOfSpecies,height=15,width=20)


fracsegsites_plot3_subsetOfSpecies <- ggplot(spectrum_and_phylo_dist_allKmers_fracSegSites[spectrum_and_phylo_dist_allKmers_fracSegSites$item1 %in% speciesToInclude & spectrum_and_phylo_dist_allKmers_fracSegSites$item2 %in% speciesToInclude,],aes(x=cophenetic_distance,y=spectrum_distance))+
  geom_point()+
  #geom_smooth()+
  facet_wrap(~label,scales="free")+
  theme_bw()+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=4)+
  ggtitle("fraction of segregating sites") #+
  #geom_smooth()
fracsegsites_plot3_subsetOfSpecies
ggsave(paste0(plotdir,"fracsegsitesplot.SUBSETOFSPECIES.png"),fracsegsites_plot3_subsetOfSpecies,height=15,width=20)

#################### want to separate distances by central mutation type ##############

combo_function_separateByMutationType <- function(spectrumdf,variable,label,ilr_or_clr,centralMutationType){
  distance_dataframe <- 
    processSpectra(spectrumdf) %>%
    filter(mutation_1mer==centralMutationType) %>%
    pivotSpectra_perVariable(.,variable)  %>%
    select(c(variable,mutation_label,speciesToInclude_1)) %>% # just get the cols for the species you want 
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
ggsave(paste0(plotdir,"1mer.fracSegSites.PerMutationType.epsilon.",epsilon,".png"),spectrum1merPlot2,height=8,width=12)

# play with dist a bit:
tidy(dist(t(select(all1mer_processed_pivot,c(vaquita,mice_Mmd)))))
tidy(dist(t(select(all1mer_processed_pivot,c(vaquita,fin_whale_ENP)))))
tidy(dist(t(select(all1mer_processed_pivot,c(vaquita,bears_ABC)))))
tidy(dist(t(select(all1mer_processed_pivot,c(vaquita,bears_PB)))))


head(all7mer_processed)
all7mer_processed[all7mer_processed$mutation_1mer=="C.T_CpG" & all7mer_processed$label=="bears_ABC",]
all1mer_processed[all1mer_processed$mutation_1mer=="C.T_CpG" & all1mer_processed$label %in% c("vaquita","bears_ABC"),]
 # target coutns? 

# proportion of targets? 
all1mer_processed <- all1mer_processed %>% 
  group_by(species,label) %>%
  mutate(total_target_proportion=total_target_count/sum(total_target_count))


targetPlot1 <- ggplot(all1mer_processed[all1mer_processed$label %in% c("mice_Mmd","vaquita","humans_AFR","fin_whale_ENP","bears_ABC"),],aes(y=label,x=total_target_count,fill=ancestral1mer))+
  geom_col()+
  facet_wrap(~ancestral1mer,scales="free_x",ncol=1)+
  ggtitle("target proportions (ancestral1mers) in callable genome")
targetPlot1
ggsave(paste0(plotdir,"targetCounts.1mer.png"),targetPlot1,width=5,height=10)

targetPlot2 <- ggplot(all1mer_processed[all1mer_processed$label %in% c("mice_Mmd","vaquita","humans_AFR","fin_whale_ENP","bears_ABC"),],aes(y=label,x=total_target_proportion,fill=ancestral1mer))+
  geom_col()+
  facet_wrap(~ancestral1mer,scales="free_x",ncol=1)+
  ggtitle("target proportions (ancestral1mers) in callable genome")
targetPlot2
ggsave(paste0(plotdir,"targetProportion.1mer.png"),targetPlot2,width=5,height=10)

### try putting CpGs back in with C>Ts
all1mer_processed$mutation_1mer_noCpGLabel <- all1mer_processed$mutation_1mer

all1mer_processed[all1mer_processed$mutation_1mer=="C.T_CpG",]$mutation_1mer_noCpGLabel <- "C.T"
all1mer_processed[all1mer_processed$mutation_1mer=="C.A_CpG",]$mutation_1mer_noCpGLabel <- "C.A"
all1mer_processed[all1mer_processed$mutation_1mer=="C.G_CpG",]$mutation_1mer_noCpGLabel <- "C.G"

all1mer_processed_noCpG <- all1mer_processed %>% group_by(species,population,label,mutation_1mer_noCpGLabel) %>%
  summarise(total_mutations_projected2=sum(total_mutations_projected),total_targets2=sum(total_target_count))

head(all1mer_processed_noCpG)
all1mer_processed_noCpG$mutation_rate <- all1mer_processed_noCpG$total_mutations_projected2/all1mer_processed_noCpG$total_targets2

all1mer_processed_noCpG <- all1mer_processed_noCpG %>%
  group_by(species,population,label) %>%
  mutate(mutation_rate_normalized=mutation_rate/sum(mutation_rate),frac_seg_sites=total_mutations_projected2/sum(total_mutations_projected2))
# check
sum(all1mer_processed_noCpG[all1mer_processed_noCpG$label=="fin_whale_GOC",]$mutation_rate_normalized)
# good
# hmmm so the C>T rate is just so high that it makes all else very small and dominates? 

head(all1mer_processed_noCpG)

spectrum1merPlot3_noCpG <- ggplot(all1mer_processed_noCpG,aes(x=mutation_1mer_noCpGLabel,y=mutation_rate_normalized,fill=label))+
  geom_col(position="dodge")
spectrum1merPlot3_noCpG
ggsave(paste0(plotdir,"1mer.mutRate.PerMutationType.epsilon.",epsilon,".NOCpG.png"),spectrum1merPlot3_noCpG,height=8,width=12)

spectrum1merPlot4_noCpG <- ggplot(all1mer_processed_noCpG,aes(x=mutation_1mer_noCpGLabel,y=frac_seg_sites,fill=label))+
  geom_col(position="dodge")
spectrum1merPlot4_noCpG
ggsave(paste0(plotdir,"1mer.fracSegSites.PerMutationType.epsilon.",epsilon,".NOCpG.png"),spectrum1merPlot4_noCpG,height=8,width=12)
