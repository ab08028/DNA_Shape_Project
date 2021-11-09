############# gather all spectra and calculate rescaled mutation rates  ###########
require(reshape2)
require(dplyr)
require(ggplot2)
require(tidyr)
require(compositions)  # for clr
require(broom) # for tidy
require(ggrepel)
plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/comparing_spectrum_and_phylo_divergence/"

tableOfSpectraToInclude=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/spectra/tableOfSpectraToIncludeInComparisons.20211013.txt",header=T,sep="\t")  # you can add to this if you want more includd ; currenlty human and fin whale missing because need targets and genome ide measures for fin whale
# FOR NOW JUST USING THE FIRST INDIVIDUAL FOR EACH COMPARISON - maybe?
speciesCodes=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/scripts/DNA_Shape_Project/analyses/PhylogeneticDistances/speciesCodes.txt",header=T)
head(speciesCodes)

###### read in phylo distances: 
phyloDistances = read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/raxml_tree/raxml.cophenetic.pairwisedistances.fromUphametal.MYSPECIES.usethis.txt",header=T)
head(phyloDistances) # not all spp are present

# okay have to work ith table 
all7merSpectra = data.frame()
all5merSpectra = data.frame()
all3merSpectra = data.frame()
all1merSpectra = data.frame()
for(sample in tableOfSpectraToInclude$sample){
  # some 'samples' may contain multipops like mice and bears
  sampleinfo = tableOfSpectraToInclude[tableOfSpectraToInclude$sample==sample,]
  sampleinfo
  spectrumfilename = paste0(sampleinfo$spectrumindir,"/",sampleinfo$spectrumFile)

  spectrum = read.table(spectrumfilename,header=T,sep="\t") # alwyas has header
  # SOME WEIRD EXCEPTIONS:
  # if it hasn't been melted yet , melt it: 
  if(dim(spectrum)[2]>dim(spectrum)[1]){
    spectrum <- melt(spectrum)
  } # means it hasn't been melted 
  # for the bears I need to do something special:
  # for bears and mice select the relevant columns (ignore the extra columns)
  # these are different because were summed up over intervals:
  if(sample %in% c("bears","mice")){
    spectrum <- select(spectrum,c(population,variable,totalSitesAllIntervals))
    colnames(spectrum) <- c("population","variable","value")
    spectrum$sample <- sample # give it a generic name
  } else if(sample == "vaquita"){
      spectrum$sample <- "vaquita"
      spectrum$population <- "vaquita"
  
  } else if(!"population" %in% names(spectrum)){
    spectrum$population <- sample
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
  
  targets=read.table(targetsfilename,header=sampleinfo$targetsHeader) # assigns whetehr has header or not
  colnames(targets) <- c("target_7mer","countOverAllIntervals")
  
  # get central 5mer, 3mer and 1mer:
  targets <- targets %>%
    mutate(target_5mer=substr(target_7mer,2,6),target_3mer=substr(target_7mer,3,5),target_1mer= substr(target_7mer,4,4)) # get substrings 

  # TARGETS per 3mer type #
  targets_7mer <- targets %>% group_by(target_7mer) %>% summarise(total_target_count=sum(countOverAllIntervals)) # original targets file, just doing to be consistent but it's identical to original targets
  targets_3mer <- targets %>% group_by(target_3mer) %>% summarise(total_target_count=sum(countOverAllIntervals))
  targets_5mer <- targets %>% group_by(target_5mer) %>% summarise(total_target_count=sum(countOverAllIntervals))
  targets_1mer <- targets %>% group_by(target_1mer) %>% summarise(total_target_count=sum(countOverAllIntervals))

  # okay time to merge:
  merged_7mer <- merge(spectrum_7mer, targets_7mer,by.x="ancestral7mer",by.y="target_7mer")
  merged_5mer <- merge(spectrum_5mer, targets_5mer,by.x="ancestral5mer",by.y="target_5mer")
  merged_3mer <- merge(spectrum_3mer, targets_3mer,by.x="ancestral3mer",by.y="target_3mer")
  merged_1mer <- merge(spectrum_1mer, targets_1mer,by.x="ancestral1mer",by.y="target_1mer")
  
  # let's combine:
  all7merSpectra = bind_rows(all7merSpectra,merged_7mer)
  all5merSpectra = bind_rows(all5merSpectra,merged_5mer)
  all3merSpectra = bind_rows(all3merSpectra,merged_3mer)
  all1merSpectra = bind_rows(all1merSpectra,merged_1mer)
  
  
}


####### define some functions ###########
# okay you want to get spectrum and the mutation rate
processSpectra <- function(spectradf){
  # get mutation rates and fraction of segregating sites:
  spectradf <- spectradf %>%
    mutate(mutation_rate = total_mutations/total_target_count) %>%
    group_by(sample,population) %>%
    mutate(mutation_rate_normalized = mutation_rate/sum(mutation_rate),fractionOfSegregatingSites=total_mutations/sum(total_mutations))
  return(spectradf)

}


#all7merSpectra_proc <- processSpectra(all7merSpectra)
#head(all7merSpectra_proc)
# need to pivot it wider? 
pivotSpectra_perVariable <- function(processedspectradf,variableToPivot) {
  processedspectradf_wider <-  pivot_wider(processedspectradf,id_cols=c(ends_with("mer")),names_from=c("population"),values_from =all_of(variableToPivot ))
  processedspectradf_wider$variable <- variableToPivot
  # pivots on a variable 
  return(processedspectradf_wider)
}

vars=c("fractionOfSegregatingSites","mutation_rate_normalized")
all7merSpectra_pivot <- pivotSpectra_perVariable(all7merSpectra_proc,"fractionOfSegregatingSites")
# get distances
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

#clr_Distance <- function(pivotedspectradf){
  # make it a matrix with each row as a species (transposed and all numeric)
  # variable in question:
  # do I CLR over whole matrix? or just over the pairs?
  # need to figure this out.
#}

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
spectrum_and_phylo_dist_1mer_fracSegSites <- 
  processSpectra(all1merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites")  %>%
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="1mer")

spectrum_and_phylo_dist_3mer_fracSegSites <- 
  processSpectra(all3merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites")  %>%
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="3mer")

spectrum_and_phylo_dist_5mer_fracSegSites <- 
  processSpectra(all5merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites")  %>%
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="5mer")

spectrum_and_phylo_dist_7mer_fracSegSites <- 
  processSpectra(all7merSpectra) %>%
  pivotSpectra_perVariable(.,"fractionOfSegregatingSites")  %>%
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="7mer")

# bind them all together:
spectrum_and_phylo_dist_allKmers_fracSegSites <- bind_rows(spectrum_and_phylo_dist_1mer_fracSegSites,spectrum_and_phylo_dist_3mer_fracSegSites,spectrum_and_phylo_dist_5mer_fracSegSites,spectrum_and_phylo_dist_7mer_fracSegSites)


######### get euclidean distances for normalized mut rates #######
spectrum_and_phylo_dist_1mer_normMutRate <- 
  processSpectra(all1merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_normalized")  %>%
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="1mer")

spectrum_and_phylo_dist_3mer_normMutRate <- 
  processSpectra(all3merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_normalized")  %>%
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="3mer")

spectrum_and_phylo_dist_5mer_normMutRate <- 
  processSpectra(all5merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_normalized")  %>%
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="5mer")

spectrum_and_phylo_dist_7mer_normMutRate <- 
  processSpectra(all7merSpectra) %>%
  pivotSpectra_perVariable(.,"mutation_rate_normalized")  %>%
  euclideanDistance() %>%
  addFullSpeciesNamesIfPresentInRaxmlTreeToDistancedf(.,speciesCodes)  %>%
  addPhyloDistances_excludeNas(.,phyloDistances) %>%
  select(comparisonLabel,item1,item2,variable,common_name.item1,common_name.item2,comparisonLabel_common_alphabetical,comparisonLabel_broad_alphabetical,spectrum_distance,cophenetic_distance)  %>%
  mutate(label="7mer")

spectrum_and_phylo_dist_allKmers_normMutRate <- bind_rows(spectrum_and_phylo_dist_1mer_normMutRate,spectrum_and_phylo_dist_3mer_normMutRate,spectrum_and_phylo_dist_5mer_normMutRate,spectrum_and_phylo_dist_7mer_normMutRate)

######### let's make some plots ###############

fracsegsites_plot1_withlabels <- ggplot(spectrum_and_phylo_dist_allKmers_fracSegSites,aes(x=cophenetic_distance,y=spectrum_distance,color=label))+
  geom_point()+
  #geom_smooth()+
  facet_wrap(~label,scales="free")+
  theme_bw()+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=4)+
  ggtitle("fraction of segregating sites")
fracsegsites_plot1_withlabels
ggsave(paste0(plotdir,"fracsegsitesplot1_withlabels.png"),fracsegsitesplot1_withlabels,height=15,width=20)


mutRateNorm_plot1_withlabels <- ggplot(spectrum_and_phylo_dist_allKmers_normMutRate,aes(x=cophenetic_distance,y=spectrum_distance,color=label))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~label,scales="free")+
  theme_bw()+
  geom_text_repel(aes(label=comparisonLabel_common_alphabetical),size=4)
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
