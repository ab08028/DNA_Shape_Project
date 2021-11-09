############# gather all spectra and calculate rescaled mutation rates  ###########
require(reshape2)
require(dplyr)
require(ggplot2)
require(tidyr)
require(compositions)  # for clr
require(broom) # for tidy
require(ggrepel)

epsilon=0.1 # amount to add to 0 bins
plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/comparing_spectrum_and_phylo_divergence/"

#tableOfSpectraToInclude=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/spectra/tableOfSpectraToIncludeInComparisons.20211013.txt",header=T,sep="\t")  # you can add to this if you want more includd ; currenlty human and fin whale missing because need targets and genome ide measures for fin whale
# restricting to just species with >1 individual:
tableOfSpectraToInclude=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/spectra/tableOfSpectraToIncludeInComparisons.SpeciesWithMultipleIndividualsOnly.20211026.txt",header=T,sep="\t")

speciesCodes=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/scripts/DNA_Shape_Project/analyses/PhylogeneticDistances/speciesCodes.txt",header=T) # for now just using AFR as human
head(speciesCodes)

###### read in phylo distances: 
phyloDistances = read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/raxml_tree/raxml.cophenetic.pairwisedistances.fromUphametal.MYSPECIES.usethis.txt",header=T)
head(phyloDistances) # not all spp are present

# dleu2 has this one weird non7mer in its spectrum-- GAAAG.GCAAG
# that messes stuff up (why didn't that come up before?)
# okay have to work ith table 
all7merSpectraOnly = data.frame()
all5merSpectraOnly = data.frame()
all3merSpectraOnly = data.frame()
all1merSpectraOnly = data.frame()

all7merTargets=data.frame()
all5merTargets=data.frame()
all3merTargets = data.frame()
all1merTargets = data.frame()
for(sample in tableOfSpectraToInclude$sample){
  # some 'samples' may contain multipops like mice and bears
  sampleinfo = tableOfSpectraToInclude[tableOfSpectraToInclude$sample==sample,]
  sampleinfo
  spectrumfilename = paste0(sampleinfo$spectrumindir,"/",sampleinfo$spectrumFile)

  spectrum = read.table(spectrumfilename,header=T,sep="\t") # alwyas has header
  # check for weird shorter than 7mer mutation types (dleu1 has one) and remove
  # SOME WEIRD EXCEPTIONS:
  # if it hasn't been melted yet , melt it: 
  if(dim(spectrum)[2]>dim(spectrum)[1]){
    spectrum <- melt(spectrum)
  } # means it hasn't been melted 
  # for the bears I need to do something special:
  # for bears and mice select the relevant columns (ignore the extra columns)
  # these are different because were summed up over intervals:
  if(sample %in% c("bears","mice","fin_whale")){
    spectrum <- select(spectrum,c(population,variable,totalSitesAllIntervals))
    colnames(spectrum) <- c("population","variable","value")
    spectrum$sample <- sample # give it a generic name
    # relabel populations to be bears_ or mice_ to make clearer what htey are 
    spectrum$population <- paste0(sample,"_",spectrum$population)
    
  } else if (sample %in% c("human")){
      spectrum <- select(spectrum,c(population,variable,countOverAllIntervals)) # has slightly diff header names
      colnames(spectrum) <- c("population","variable","value")
      spectrum$sample <- sample # give it a generic name
      # label each pop as human_ to disambiguate with bears
      spectrum$population <- paste0(sample,"_",spectrum$population)
  } else if(sample == "vaquita"){
      spectrum$sample <- "vaquita"
      spectrum$population <- "vaquita"
  
  } else if(!"population" %in% names(spectrum)){
    spectrum$population <- sample
  }
  # check for length of mutation type (dleu2) has one weird 5mer in the 7mer spectrum 
  numberOfNon7mers <- sum(nchar(as.character(spectrum$variable))!=15)
  
  if(numberOfNon7mers!=0){
    print(paste0("there are",numberOfNon7mers," non-7mer mutation(s) presenet for",sample,"; they are being removed from spectrum:"))
    # print them out:
    spectrum[nchar(as.character(spectrum$variable))!=15,]
    spectrum <- spectrum[nchar(as.character(spectrum$variable))==15,]
        
  }  
  # get ancestral kmers:
  spectrum$ancestral7mer <- substr(spectrum$variable,1,7)
  spectrum$ancestral5mer <- substr(spectrum$variable,2,6)
  spectrum$ancestral3mer <- substr(spectrum$variable,3,5)
  spectrum$ancestral1mer <- substr(spectrum$variable,4,4)
  
  # get central mutation types:
  spectrum$mutation_7mer <- spectrum$variable
  spectrum$mutation_5mer <- paste0(substr(spectrum$variable,2,6),".",substr(spectrum$variable,10,14))
  spectrum$mutation_3mer <- paste0(substr(spectrum$variable,3,5),".",substr(spectrum$variable,11,13))
  spectrum$mutation_1mer <- paste0(substr(spectrum$variable,4,4),".",substr(spectrum$variable,12,12))
  
  # get spectrum summed up her type
  spectrum_7mer <- spectrum %>%
    group_by(sample,population,ancestral7mer,mutation_7mer) %>%
    summarise(total_mutations = sum(value)) # this will be identical to original spectrum (7mer), just doing for consistency of formatting
  
  spectrum_5mer <- spectrum %>%
    group_by(sample,population,ancestral5mer,mutation_5mer) %>%
    summarise(total_mutations = sum(value)) 
  
  spectrum_3mer <- spectrum %>%
    group_by(sample,population,ancestral3mer,mutation_3mer) %>%
    summarise(total_mutations = sum(value)) 
  
  spectrum_1mer <- spectrum %>%
    group_by(sample,population,ancestral1mer,mutation_1mer) %>%
    summarise(total_mutations = sum(value)) 
  
  # targets:
  targetsfilename = paste0(sampleinfo$targetsindir,"/",sampleinfo$targetsFile)
  
  targets=read.table(targetsfilename,header=sampleinfo$targetsHeader) # assigns whether has header or not
  colnames(targets) <- c("target_7mer","countOverAllIntervals")
  
  # get central 5mer, 3mer and 1mer:
  targets <- targets %>%
    mutate(target_5mer=substr(target_7mer,2,6),target_3mer=substr(target_7mer,3,5),target_1mer= substr(target_7mer,4,4)) # get substrings 

  # TARGETS per 3mer type #
  targets_7mer <- targets %>% group_by(target_7mer) %>% summarise(total_target_count=sum(countOverAllIntervals)) # original targets file, just doing to be consistent but it's identical to original targets
  targets_3mer <- targets %>% group_by(target_3mer) %>% summarise(total_target_count=sum(countOverAllIntervals))
  targets_5mer <- targets %>% group_by(target_5mer) %>% summarise(total_target_count=sum(countOverAllIntervals))
  targets_1mer <- targets %>% group_by(target_1mer) %>% summarise(total_target_count=sum(countOverAllIntervals))
  targets_7mer$sample <- sample
  targets_5mer$sample <- sample
  targets_3mer$sample <- sample
  targets_1mer$sample <- sample
  # make all species targets:
  

  # let's combine all spectra
  all7merSpectraOnly = bind_rows(all7merSpectraOnly,spectrum_7mer)
  all5merSpectraOnly = bind_rows(all5merSpectraOnly,spectrum_5mer)
  all3merSpectraOnly = bind_rows(all3merSpectraOnly,spectrum_3mer)
  all1merSpectraOnly = bind_rows(all1merSpectraOnly,spectrum_1mer)
  
  # and all targets:
  all7merTargets = bind_rows(all7merTargets,targets_7mer)
  all5merTargets = bind_rows(all5merTargets,targets_5mer)
  all3merTargets = bind_rows(all3merTargets,targets_3mer)
  all1merTargets = bind_rows(all1merTargets,targets_1mer)
  

  
}

# complete so that missing mutation types are filled in 
all7merSpectraOnly_filledin <- all7merSpectraOnly %>%
  complete(nesting(mutation_7mer,ancestral7mer),nesting(sample,population),fill=list(total_mutations=0))
all5merSpectraOnly_filledin <- all5merSpectraOnly %>%
  complete(nesting(mutation_5mer,ancestral5mer),nesting(sample,population),fill=list(total_mutations=0))
all3merSpectraOnly_filledin <- all3merSpectraOnly %>%
  complete(nesting(mutation_3mer,ancestral3mer),nesting(sample,population),fill=list(total_mutations=0))
all1merSpectraOnly_filledin <- all1merSpectraOnly %>%
  complete(nesting(mutation_1mer,ancestral1mer),nesting(sample,population),fill=list(total_mutations=0))


# so when I was merging before all species were combined there was a bug where missing mtuation types wouldn't get targets filled in which led to them being NA downstream 
# okay now fill in missing mutation types : as 0s PRIOR TO MERGING WITH TARGEST

all7merSpectra <- merge(all7merSpectraOnly_filledin, all7merTargets,by.x=c("sample","ancestral7mer"),by.y=c("sample","target_7mer"))
all5merSpectra <- merge(all5merSpectraOnly_filledin, all5merTargets,by.x=c("sample","ancestral5mer"),by.y=c("sample","target_5mer"))
all3merSpectra <- merge(all3merSpectraOnly_filledin, all3merTargets,by.x=c("sample","ancestral3mer"),by.y=c("sample","target_3mer"))
all1merSpectra <- merge(all1merSpectraOnly_filledin, all1merTargets,by.x=c("sample","ancestral1mer"),by.y=c("sample","target_1mer"))

#### okay bug fixed--- these now have targets filled in

####### define some functions ###########
# okay you want to get spectrum and the mutation rate and want to add 1 to all counts to deal with 0 entries (assuming every mutation type would happen if you sampled enough -- though note this will have a bigger i)
processSpectra <- function(spectradf){
  # get mutation rates and fraction of segregating sites: also adding a column in where I add +epsilon (1) to every count to deal with 0 entries ; adding it to targets as well since FWs missing some targets (is that interesting?)
  spectradf <- spectradf %>%
    mutate(total_mutations_plusEpsilon = (total_mutations+epsilon),total_target_count_plusEpsilon=(total_target_count+epsilon),mutation_rate = (total_mutations/total_target_count), mutation_rate_plusEpsilon=(total_mutations_plusEpsilon/total_target_count_plusEpsilon)) %>%
    group_by(sample,population) %>%
    mutate(mutation_rate_normalized = (mutation_rate/sum(mutation_rate)),fractionOfSegregatingSites=(total_mutations/sum(total_mutations)), mutation_rate_normalized_plusEpsilon=(mutation_rate_plusEpsilon/sum(mutation_rate_plusEpsilon)),fractionOfSegregatingSites_plusEpsilon=(total_mutations_plusEpsilon/sum(total_mutations_plusEpsilon))) # adding plus1 versions.
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
  processedspectradf_wider <-  pivot_wider(processedspectradf,id_cols=c(ends_with("mer")),names_from=c("population"),values_from =all_of(variableToPivot ))
  processedspectradf_wider$variable <- variableToPivot
  # pivots on a variable 
  return(processedspectradf_wider)
}
####### HOW TO deal with missing targets? (fin whale) -- need to deal with. TALK TO KH about nullomers. Hm.
# maybe be chiller about how I call targets in FW and vaq? maybe collapse down to 5mers? 
# test function
testpivoted = pivotSpectra_perVariable(test,"mutation_rate_normalized_plusEpsilon")
head(testpivoted)
# want to speed up with only a subset of species.
# ilr/clr calcs require no 0 entries
clr_calculations <- function(pivotedspectradf_plusEpsilonVersionOfVariables){
  # remove character vectors:
  tableWithoutIDVars <- select(pivotedspectradf_plusEpsilonVersionOfVariables, -c("variable",ends_with("mer")))
  # AHA so if you run clr just on a matrix it goes ROW-WISE not column wise! annoying! so need to transpose
  #clr_matrix <- clr(t(tableWithoutIDVars)) # if you do it this way it matches when I do it just on on un-tranposed column by itself. if you do it without transposing then you get weird row-wise results that are WRONG. so be cautious here. 
  # so alternative without tranposing is to use sapply -- this is a lot nicer I think
  clr_df <- data.frame(sapply(tableWithoutIDVars,clr)) # this way you don't have to transpose and you get a much nicer shaped output
  return(clr_df)
  
}
#vars=c("fractionOfSegregatingSites","mutation_rate_normalized")
#all7merSpectra_pivot <- pivotSpectra_perVariable(all7merSpectra_proc,"fractionOfSegregatingSites")

ilr_calculations <- function(pivotedspectradf_plusEpsilonVersionOfVariables){
  # remove character vectors:
  tableWithoutIDVars <- select(pivotedspectradf_plusEpsilonVersionOfVariables, -c("variable",ends_with("mer")))
  # AHA so if you run clr just on a matrix it goes ROW-WISE not column wise! annoying! so need to transpose
  #clr_matrix <- clr(t(tableWithoutIDVars)) # if you do it this way it matches when I do it just on on un-tranposed column by itself. if you do it without transposing then you get weird row-wise results that are WRONG. so be cautious here. 
  # so alternative without tranposing is to use sapply -- this is a lot nicer I think
  ilr_df <- data.frame(sapply(tableWithoutIDVars,ilr)) # this way you don't have to transpose and you get a much nicer shaped output
  return(ilr_df)
  
}

# get distances
# for regular distance. 
euclideanDistance <- function(pivotedspectradf){
  # make it a matrix with each row as a species (transposed and all numeric)
  # variable in question:
  variable = unique(pivotedspectradf$variable)
  matrix <- t(select(pivotedspectradf,-c(ends_with("mer"),"variable")))
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
speciesToInclude_1=c("human_AFR","mice_Mmd","mice_Ms","bears_ABC","bears_PB","vaquita","fin_whale_ENP")
# comparisons to include : 
# mouse-mouse, bear-bear, whale-whale, human-mouse, human-bear, human-whale 

########## investigating missing targets #######

targetCounts7mer <- 
  processSpectra(all7merSpectra) %>%
  pivotSpectra_perVariable(.,"total_target_count")
write.table(targetCounts7mer,paste(plotdir,"allSpecies.TargetCounts.7mer.txt"),quote=F,row.names = F,sep="\t")

targetCounts5mer <- 
  processSpectra(all5merSpectra) %>%
  pivotSpectra_perVariable(.,"total_target_count")
write.table(targetCounts5mer,paste(plotdir,"allSpecies.TargetCounts.5mer.txt"),quote=F,row.names = F,sep="\t")

targetCounts3mer <- 
  processSpectra(all3merSpectra) %>%
  pivotSpectra_perVariable(.,"total_target_count")
write.table(targetCounts3mer,paste(plotdir,"allSpecies.TargetCounts.3mer.txt"),quote=F,row.names = F,sep="\t")


###### Get ILR and CLR distances ########
# this is slow so restrict just to species you want
##### this subsets just to the species you want above ^^ 
############# clr mutation rate ##############

# Make a function:
#comboFunction_FullyProcessOutputDistances <- function(inputSpectrum,variableToPivotOver,label) later 
clr_distance_rescaledMutRate_1mer <- 
  processSpectra(all1merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_normalized_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),speciesToInclude_1)) %>% # just get the cols for the species you want 
  clr_calculations(.) %>%
  euclideanDistance_ilr_or_clr(.) %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="1mer_clr",variable="mutation_rate_normalized_plusEpsilon")


clr_distance_rescaledMutRate_3mer <- 
  
  processSpectra(all3merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_normalized_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  clr_calculations(.) %>%
  euclideanDistance_ilr_or_clr(.) %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="3mer_clr",variable="mutation_rate_normalized_plusEpsilon")

clr_distance_rescaledMutRate_5mer <- 
  processSpectra(all5merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_normalized_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  clr_calculations(.) %>%
  euclideanDistance_ilr_or_clr(.) %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="5mer_clr",variable="mutation_rate_normalized_plusEpsilon")

clr_distance_rescaledMutRate_7mer <- 
  processSpectra(all7merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_normalized_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  clr_calculations(.) %>%
  euclideanDistance_ilr_or_clr(.) %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="7mer_clr",variable="mutation_rate_normalized_plusEpsilon")

############### ilr mutation rate #############
###### note ILR will result in n-1 bins of vector 
ilr_distance_rescaledMutRate_1mer <- 
  processSpectra(all1merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_normalized_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  ilr_calculations(.) %>%
  euclideanDistance_ilr_or_clr(.) %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="1mer_ilr",variable="mutation_rate_normalized_plusEpsilon")

ilr_distance_rescaledMutRate_3mer <- 
  processSpectra(all3merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_normalized_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  ilr_calculations(.) %>%
  euclideanDistance_ilr_or_clr(.) %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="3mer_ilr",variable="mutation_rate_normalized_plusEpsilon")



ilr_distance_rescaledMutRate_5mer <- 
  processSpectra(all5merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_normalized_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  ilr_calculations(.) %>%
  euclideanDistance_ilr_or_clr(.) %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="5mer_ilr",variable="mutation_rate_normalized_plusEpsilon")

# ilr_distance_rescaledMutRate_7mer <- 
#   processSpectra(all7merSpectra) %>%
#   pivotSpectra_perVariable(.,"mutation_rate_normalized_plusEpsilon")  %>%
#   select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
#   ilr_calculations(.) %>%
#   euclideanDistance_ilr_or_clr(.) %>%
#   addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
#   addPhyloDistances_excludeNas(.,phyloDistances) %>%
#   select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
#   mutate(label="7mer_ilr",variable="mutation_rate_normalized_plusEpsilon")

############# combine mut rate dists ################
all_clr_ilr_dists_mutrate <- bind_rows(ilr_distance_rescaledMutRate_1mer,ilr_distance_rescaledMutRate_3mer,ilr_distance_rescaledMutRate_5mer,clr_distance_rescaledMutRate_1mer,clr_distance_rescaledMutRate_3mer,clr_distance_rescaledMutRate_5mer,clr_distance_rescaledMutRate_7mer)

write.table(all_clr_ilr_dists_mutrate,paste0(plotdir,"ilr.clr.distances.RescaledMutationRate.epsilon.",as.character(epsilon),".txt"),row.names=F,quote=F,sep="\t")

all_clr_ilr_dists_mutrate$label <- factor(all_clr_ilr_dists_mutrate$label,levels=c("1mer_clr","3mer_clr","5mer_clr","7mer_clr","1mer_ilr","3mer_ilr","5mer_ilr","7mer_ilr"))
ilrclr_plot1 <- ggplot(all_clr_ilr_dists_mutrate,aes(x=cophenetic_distance,y=spectrum_distance))+
  facet_wrap(~label,scales="free",ncol=4)+
  geom_point()+
  geom_smooth()+
  ggtitle("mut rate")
ilrclr_plot1
ggsave(paste0(plotdir,"ilr.clr.distances.mutationRate.plusEpsilon.",as.character(epsilon),".png"),ilrclr_plot1,width=13,height=7)
############ CLR fraction of seg sites ##########
clr_distance_fracSegSites_1mer <- 
  processSpectra(all1merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  clr_calculations(.) %>%
  euclideanDistance_ilr_or_clr(.) %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="1mer_clr",variable="fractionOfSegregatingSites_plusEpsilon")

clr_distance_fracSegSites_3mer <- 
  processSpectra(all3merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  clr_calculations(.) %>%
  euclideanDistance_ilr_or_clr(.) %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="3mer_clr",variable="fractionOfSegregatingSites_plusEpsilon")


clr_distance_fracSegSites_5mer <- 
  processSpectra(all5merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  clr_calculations(.) %>%
  euclideanDistance_ilr_or_clr(.) %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="5mer_clr",variable="fractionOfSegregatingSites_plusEpsilon")

clr_distance_fracSegSites_7mer <- 
  processSpectra(all7merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  clr_calculations(.) %>%
  euclideanDistance_ilr_or_clr(.) %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="7mer_clr",variable="fractionOfSegregatingSites_plusEpsilon")
######### ILR frac seg sites ##############
ilr_distance_fracSegSites_1mer <- 
  processSpectra(all1merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  ilr_calculations(.) %>%
  euclideanDistance_ilr_or_clr(.) %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="1mer_ilr",variable="fractionOfSegregatingSites_plusEpsilon")

ilr_distance_fracSegSites_3mer <- 
  processSpectra(all3merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  ilr_calculations(.) %>%
  euclideanDistance_ilr_or_clr(.) %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="3mer_ilr",variable="fractionOfSegregatingSites_plusEpsilon")


ilr_distance_fracSegSites_5mer <- 
  processSpectra(all5merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  ilr_calculations(.) %>%
  euclideanDistance_ilr_or_clr(.) %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="5mer_ilr",variable="fractionOfSegregatingSites_plusEpsilon")

# ilr_distance_fracSegSites_7mer <- 
#   processSpectra(all7merSpectra) %>%
#   pivotSpectra_perVariable(.,"fractionOfSegregatingSites_plusEpsilon")  %>%
#   select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
#   ilr_calculations(.) %>%
#   euclideanDistance_ilr_or_clr(.) %>%
#   addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
#   addPhyloDistances_excludeNas(.,phyloDistances) %>%
#   select(comparisonLabel,item1,item2,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
#   mutate(label="7mer_ilr",variable="fractionOfSegregatingSites_plusEpsilon")


####### combine and plot ######


all_clr_ilr_dists_fracSegSites <- bind_rows(ilr_distance_fracSegSites_1mer,ilr_distance_fracSegSites_3mer,ilr_distance_fracSegSites_5mer,clr_distance_fracSegSites_1mer,clr_distance_fracSegSites_3mer,clr_distance_fracSegSites_5mer,clr_distance_fracSegSites_7mer) # keeping out for now ilr_distance_fracSegSites_7mer


write.table(all_clr_ilr_dists_fracSegSites,paste0(plotdir,"ilr.clr.distances.FracSegSites.epsilon.",as.character(epsilon),".txt"),row.names=F,quote=F,sep="\t")

all_clr_ilr_dists_fracSegSites$label <- factor(all_clr_ilr_dists_fracSegSites$label,levels=c("1mer_clr","3mer_clr","5mer_clr","7mer_clr","1mer_ilr","3mer_ilr","5mer_ilr","7mer_ilr"))


ilrclr_plot2 <- ggplot(all_clr_ilr_dists_fracSegSites,aes(x=cophenetic_distance,y=spectrum_distance))+
  facet_wrap(~label,scales="free",ncol=4)+
  geom_point()+
  geom_smooth()+
  ggtitle("frac seg sites")
ilrclr_plot2
ggsave(paste0(plotdir,"ilr.clr.distances.fracSegSites.plusEpsilon.",as.character(epsilon),".png"),ilrclr_plot2,width=13,height=7)





####### YOU ARE HERE - play with these ! #########
# why does fin whale have slightly larger cophenetic distance than vaquita from mouse? shouldn't they be equidistant? hm. oh but it's substitutions so includes changes in gen time etc?

############ regular stuff non clr ilr ##########
spectrum_and_phylo_dist_1mer_fracSegSites <- 
  processSpectra(all1merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="1mer_noTransform_plusEpsilon")

spectrum_and_phylo_dist_3mer_fracSegSites <- 
  processSpectra(all3merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="3mer_noTransform_plusEpsilon")

spectrum_and_phylo_dist_5mer_fracSegSites <- 
  processSpectra(all5merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="5mer_noTransform_plusEpsilon")

spectrum_and_phylo_dist_7mer_fracSegSites <- 
  processSpectra(all7merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="7mer_noTransform_plusEpsilon")

# bind them all together:
spectrum_and_phylo_dist_allKmers_fracSegSites <- bind_rows(spectrum_and_phylo_dist_1mer_fracSegSites,spectrum_and_phylo_dist_3mer_fracSegSites,spectrum_and_phylo_dist_5mer_fracSegSites,spectrum_and_phylo_dist_7mer_fracSegSites)


######### get euclidean distances for normalized mut rates (use the +epsilon version!) #######
spectrum_and_phylo_dist_1mer_normMutRate <- 
  processSpectra(all1merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_normalized_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="1mer_noTransform_plusEpsilon")

spectrum_and_phylo_dist_3mer_normMutRate <- 
  processSpectra(all3merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_normalized_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="3mer_noTransform_plusEpsilon")

spectrum_and_phylo_dist_5mer_normMutRate <- 
  processSpectra(all5merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_normalized_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="5mer_noTransform_plusEpsilon")

spectrum_and_phylo_dist_7mer_normMutRate <- 
  processSpectra(all7merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_normalized_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1))) %>% # just get the cols for the species you want
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="7mer_noTransform_plusEpsilon")

spectrum_and_phylo_dist_allKmers_normMutRate <- bind_rows(spectrum_and_phylo_dist_1mer_normMutRate,spectrum_and_phylo_dist_3mer_normMutRate,spectrum_and_phylo_dist_5mer_normMutRate,spectrum_and_phylo_dist_7mer_normMutRate)


########## combine raw euclidean with transforms to plot ########

allDistanceTypes_mutRate <- bind_rows(all_clr_ilr_dists_mutrate,spectrum_and_phylo_dist_allKmers_normMutRate)

allDistanceTypes_fracSegSites <- bind_rows(all_clr_ilr_dists_fracSegSites,spectrum_and_phylo_dist_allKmers_fracSegSites)


allDistanceTypes_mutRate$label <- factor(allDistanceTypes_mutRate$label, levels=c("1mer_noTransform_plusEpsilon","3mer_noTransform_plusEpsilon","5mer_noTransform_plusEpsilon","7mer_noTransform_plusEpsilon","1mer_clr","3mer_clr","5mer_clr","7mer_clr","1mer_ilr","3mer_ilr","5mer_ilr","7mer_ilr"))

allDistanceTypes_fracSegSites$label <- factor(allDistanceTypes_fracSegSites$label, levels=c("1mer_noTransform_plusEpsilon","3mer_noTransform_plusEpsilon","5mer_noTransform_plusEpsilon","7mer_noTransform_plusEpsilon","1mer_clr","3mer_clr","5mer_clr","7mer_clr","1mer_ilr","3mer_ilr","5mer_ilr","7mer_ilr"))

allDistanceTypesPlot1 <- ggplot(allDistanceTypes_mutRate,aes(x=cophenetic_distance,y=spectrum_distance))+
  facet_wrap(~label,scales="free_y",ncol=4)+
  geom_point()+
  ggtitle("mutation rate")+
  theme_bw()+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=2)+
  theme(text=element_text(size=14))



allDistanceTypesPlot1
ggsave(paste0(plotdir,"allDistanceTypes.AllWithPlus1Transform.MutationRate.epsilon.",as.character(epsilon),".png"),allDistanceTypesPlot1,height=8,width=12)

allDistanceTypesPlot2 <- ggplot(allDistanceTypes_mutRate,aes(x=cophenetic_distance,y=spectrum_distance))+
  facet_wrap(~label,scales="free_y",ncol=4)+
  geom_point()+
  ggtitle("fraction of segregating sites")+
  theme_bw()+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=2)+
  theme(text=element_text(size=14))

allDistanceTypesPlot2
ggsave(paste0(plotdir,"allDistanceTypes.AllWithPlus1Transform.FracSegSites.epsilon.",as.character(epsilon),".png"),allDistanceTypesPlot2,height=8,width=12)

######### let's make some plots ###############

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
############ YOU ARE HERE ###########

########## try without CpG Sites?? ###############
###### simplify down to just species for which I have multi individual data ######
# species to include
speciesToInclude=c("human_AFR","mice_Mmd","mice_Ms","bears_ABC","bears_PB","vaquita")

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

########## exploring 1mer spectra ########
test1mer <-    processSpectra(all1merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites_plusEpsilon")  %>%
  select(c(variable,ends_with("mer"),all_of(speciesToInclude_1)))
ggplot(test1mer,aes(x=vaquita,y=fin_whale_ENP))+
  geom_point()+
  geom_label_repel(aes(label=mutation_1mer))+
  geom_abline()


ggplot(test1mer,aes(x=mice_Ms,y=fin_whale_ENP))+
  geom_point()+
  geom_label_repel(aes(label=mutation_1mer))+
  geom_abline()


ggplot(test1mer,aes(x=vaquita,y=mice_Ms))+
  geom_point()+
  geom_label_repel(aes(label=mutation_1mer))+
  geom_abline()
