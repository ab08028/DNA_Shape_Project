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
speciesCodes=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/scripts/DNA_Shape_Project/analyses/PhylogeneticDistances/speciesCodes.txt",header=T) # for now just using AFR as human
head(speciesCodes)

# exclude humans
tableOfSpectraToInclude <- tableOfSpectraToInclude[tableOfSpectraToInclude$sample!="human",]
###### read in phylo distances: 
phyloDistances = read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/raxml_tree/raxml.cophenetic.pairwisedistances.fromUphametal.MYSPECIES.usethis.txt",header=T)
head(phyloDistances) # not all spp are present

# dleu2 has this one weird non7mer in its spectrum-- GAAAG.GCAAG
# that messes stuff up (why didn't that come up before?)
# okay have to work ith table 
all7merSpectra = data.frame()
#all5merSpectra = data.frame()
#all3merSpectra = data.frame()
#all1merSpectra = data.frame()
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
  #spectrum$ancestral5mer <- substr(spectrum$variable,2,6)
  #spectrum$ancestral3mer <- substr(spectrum$variable,3,5)
  #spectrum$ancestral1mer <- substr(spectrum$variable,4,4)
  
  # get central mutation types:
  spectrum$mutation_7mer <- spectrum$variable
  #spectrum$mutation_5mer <- paste0(substr(spectrum$variable,2,6),".",substr(spectrum$variable,10,14))
  #spectrum$mutation_3mer <- paste0(substr(spectrum$variable,3,5),".",substr(spectrum$variable,11,13))
  #spectrum$mutation_1mer <- paste0(substr(spectrum$variable,4,4),".",substr(spectrum$variable,12,12))
  
  # get spectrum summed up her type
  spectrum_7mer <- spectrum %>%
    group_by(sample,population,ancestral7mer,mutation_7mer) %>%
    summarise(total_mutations = sum(value)) # this will be identical to original spectrum (7mer), just doing for consistency of formatting
  
 # spectrum_5mer <- spectrum %>%
  #  group_by(sample,population,ancestral5mer,mutation_5mer) %>%
  #  summarise(total_mutations = sum(value)) 
  
 # spectrum_3mer <- spectrum %>%
 #   group_by(sample,population,ancestral3mer,mutation_3mer) %>%
 #   summarise(total_mutations = sum(value)) 
  
#  spectrum_1mer <- spectrum %>%
  #  group_by(sample,population,ancestral1mer,mutation_1mer) %>%
  #  summarise(total_mutations = sum(value)) 
  
  
  # targets:
  targetsfilename = paste0(sampleinfo$targetsindir,"/",sampleinfo$targetsFile)
  
  targets=read.table(targetsfilename,header=sampleinfo$targetsHeader) # assigns whether has header or not
  colnames(targets) <- c("target_7mer","countOverAllIntervals")
  
  # get central 5mer, 3mer and 1mer:
  targets <- targets %>%
    mutate(target_5mer=substr(target_7mer,2,6),target_3mer=substr(target_7mer,3,5),target_1mer= substr(target_7mer,4,4)) # get substrings 
  
  # TARGETS per 3mer type #
  targets_7mer <- targets %>% group_by(target_7mer) %>% summarise(total_target_count=sum(countOverAllIntervals)) # original targets file, just doing to be consistent but it's identical to original targets
 # targets_3mer <- targets %>% group_by(target_3mer) %>% summarise(total_target_count=sum(countOverAllIntervals))
 # targets_5mer <- targets %>% group_by(target_5mer) %>% summarise(total_target_count=sum(countOverAllIntervals))
 # targets_1mer <- targets %>% group_by(target_1mer) %>% summarise(total_target_count=sum(countOverAllIntervals))
  
  # okay time to merge:
  merged_7mer <- merge(spectrum_7mer, targets_7mer,by.x="ancestral7mer",by.y="target_7mer")
  #merged_5mer <- merge(spectrum_5mer, targets_5mer,by.x="ancestral5mer",by.y="target_5mer")
 # merged_3mer <- merge(spectrum_3mer, targets_3mer,by.x="ancestral3mer",by.y="target_3mer")
 # merged_1mer <- merge(spectrum_1mer, targets_1mer,by.x="ancestral1mer",by.y="target_1mer")
  
  # let's combine:
  all7merSpectra = bind_rows(all7merSpectra,merged_7mer)
  #all5merSpectra = bind_rows(all5merSpectra,merged_5mer)
  #all3merSpectra = bind_rows(all3merSpectra,merged_3mer)
  #all1merSpectra = bind_rows(all1merSpectra,merged_1mer)
  
  
}

# just want rates of one 7mer so not filling in for now 


########## carry out CLR and ILR transformations ##############

####### define some functions ###########
# okay you want to get spectrum and the mutation rate and want to add 1 to all counts to deal with 0 entries (assuming every mutation type would happen if you sampled enough -- though note this will have a bigger i)
processSpectra <- function(spectradf){
  # get mutation rates and fraction of segregating sites: also adding a column in where I add +1 to every count to deal with 0 entries 
  spectradf <- spectradf %>%
    mutate(total_mutations_plus1 = (total_mutations+1),mutation_rate = (total_mutations/total_target_count), mutation_rate_plus1=(total_mutations_plus1/total_target_count)) %>%
    group_by(sample,population) %>%
    mutate(mutation_rate_normalized = (mutation_rate/sum(mutation_rate)),fractionOfSegregatingSites=(total_mutations/sum(total_mutations)), mutation_rate_normalized_plus1=(mutation_rate_plus1/sum(mutation_rate_plus1)),fractionOfSegregatingSites_plus1=(total_mutations_plus1/sum(total_mutations_plus1))) # adding plus1 versions.
  return(spectradf)
  
}
# test function:
all7merSpectra_processed = processSpectra(all7merSpectra)

# pull out one of interest
# 7mer of interest
desired7mer="TTTAAAA.TTTTAAA"
# relative mutation rate per species
# 
head(all7merSpectra_processed)
desired7merdf <- all7merSpectra_processed[all7merSpectra_processed$mutation_7mer==desired7mer,]

desired7merdf_withname <- merge(desired7merdf,speciesCodes,by.x="population",by.y="code")

plot1 <- ggplot(desired7merdf_withname,aes(y=common_name,x=mutation_rate_normalized))+
  geom_point()+
  ggtitle(paste0(desired7mer," with rescaled mutation rates. note how low fin whale is... filtering issue?"))
plot1
ggsave("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/diagonal_comparison_plots/TTTAAAA.TTTTAAA.RescaledMutationRate.noteLowFinWhale.png",plot1,height=8,width=8)
