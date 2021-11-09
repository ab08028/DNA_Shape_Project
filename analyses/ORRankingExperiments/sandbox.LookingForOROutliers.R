require(dplyr)
require(reshape2)
require(tidyr)
require(GGally)
require(SARP.compo) # for getting all ratios
plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/diagonal_comparison_plots/"
# want 3mer 5mer 7mer and want to show correlations 
# want it to be corrected for targets ideally
# proportion of segregating sites ?
# or rescaled mut rate <-- do this way

# whales, mice, bears, (condor leave out for now), humans 

#allCetaceans # /Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_spectrum_files/noSTRICTflag/AllCetaceans_melted.IncludesTwoConditions.HetsOnly.and.HetsPlusHoms.NOSTRICTFLAG.txt

# cetacean targets /Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_targets_files

# maybe just pick a couple cetaceans?
# let's do

# 4 wild mice, (fin whale -- has weirdness with being neutral so let's skip for now), minke whale, vaquita, polar bear, brown bear (eur) [condor? human?]

## maybe for now we do a case study -- just show that two mouse species and two whale species and mouse + whale -- I think that might be good.

##################### mice ################
mouseSpectrum <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/summaryTables/Mouse7merSpectrum.SummedUpOverAllChrs.NOSTRICT.txt",header=T)

head(mouseSpectrum)
# get anc 7mer:
mouseSpectrum$ancestral7mer <- unlist(lapply(strsplit(mouseSpectrum$variable,"\\."),"[",1))
mouseSpectrum$derived7mer <- unlist(lapply(strsplit(mouseSpectrum$variable,"\\."),"[",2))
head(mouseSpectrum)
tail(mouseSpectrum)
mouseTargets <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/summaryTables/ALLINTERVALS.mouse.mutyper.targets.7mer.nostrict.txt",header=T)

# merge them:
head(mouseTargets)
colnames(mouseTargets) <- c("ancestral7mer","targetCount_allIntervals")

mouseAll <- merge(mouseSpectrum,mouseTargets,by="ancestral7mer")



mouseAll$mutationRate_unscaled <- mouseAll$totalSitesAllIntervals / mouseAll$targetCount_allIntervals

mouseAll <- mouseAll %>%
  group_by(population) %>%
  mutate(mutationRate_rescaled = mutationRate_unscaled/sum(mutationRate_unscaled))

mouseAll_wide <- pivot_wider(mouseAll[,c("population","variable","mutationRate_rescaled")],id_cols = c("variable"),values_from = "mutationRate_rescaled",names_from = "population")

head(mouseAll_wide)
############# pick a whale #############
# ok let's just do minke what the hell ; for now? maybe human would be better but don't have targets on hand. maybe want to keep vaquita private?
whale1Choice="bacu2" # minke whale 
whale1ChoiceNiceLabel="minke whale"
allCetaceansSpectrum <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_spectrum_files/noSTRICTflag/AllCetaceans_melted.IncludesTwoConditions.HetsOnly.and.HetsPlusHoms.NOSTRICTFLAG.txt",sep="\t",header=T)

# restrict to just one whale for now
whale1Spectrum <- allCetaceansSpectrum[allCetaceansSpectrum$sample==whale1Choice & allCetaceansSpectrum$condition=="hetsonly",]

whale1Spectrum$ancestral7mer <- unlist(lapply(strsplit(whale1Spectrum$variable,"\\."),"[",1))


# species specific targets: 
whale1targetfilename=list.files("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_targets_files", pattern=whale1Choice,full.names = T)
whale1targetfilename #make sure this is right
whale1targets <- read.table(whale1targetfilename ,header=F)
# careful here; make sure it's the right file 
head(whale1targets)
colnames(whale1targets) <- c("ancestral7mer","targetCount_allIntervals")

whale1All <- merge(whale1Spectrum,whale1targets,by="ancestral7mer")

whale1All$mutationRate_unscaled <- whale1All$value / whale1All$targetCount_allIntervals

whale1All <- whale1All %>%
  group_by(sample) %>%
  mutate(mutationRate_rescaled = mutationRate_unscaled/sum(mutationRate_unscaled))

whale1All_wide = pivot_wider(whale1All[,c("sample","variable","mutationRate_rescaled")],id_cols = c("variable"),values_from = "mutationRate_rescaled",names_from = "sample")
head(whale1All_wide)

######## pick a second whale #######
# bottlenose dolphin
whale2Choice="ttru1" # bottlenosedolphin 
whale2ChoiceNiceLabel="dolphin" # bottlenosedolphin 

# restrict to just one whale for now
whale2Spectrum <- allCetaceansSpectrum[allCetaceansSpectrum$sample==whale2Choice & allCetaceansSpectrum$condition=="hetsonly",]

whale2Spectrum$ancestral7mer <- unlist(lapply(strsplit(whale2Spectrum$variable,"\\."),"[",1))


# species specific targets: 
whale2targetfilename=list.files("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_targets_files", pattern=whale2Choice,full.names = T)
whale2targetfilename #make sure this is right
whale2targets <- read.table(whale2targetfilename ,header=F)
# careful here; make sure it's the right file 
head(whale2targets)
colnames(whale2targets) <- c("ancestral7mer","targetCount_allIntervals")

whale2All <- merge(whale2Spectrum,whale2targets,by="ancestral7mer")

whale2All$mutationRate_unscaled <- whale2All$value / whale2All$targetCount_allIntervals

whale2All <- whale2All %>%
  group_by(sample) %>%
  mutate(mutationRate_rescaled = mutationRate_unscaled/sum(mutationRate_unscaled))

whale2All_wide = pivot_wider(whale2All[,c("sample","variable","mutationRate_rescaled")],id_cols = c("variable"),values_from = "mutationRate_rescaled",names_from = "sample")
head(whale2All_wide)
############# add in bears ###########
bearSpectrum <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210331_7mer_NoStrict/mapped_to_brown_bear/summaryTables/mapped_to_brown_bear.AllPops.bears.7merspectrum.SummedOverAllIntervals.ALLFREQS.NOSTRICT.txt",header=T,sep="\t")
head(bearSpectrum)
bearSpectrum$ancestral7mer <- unlist(lapply(strsplit(bearSpectrum$variable,"\\."),"[",1))
head(bearSpectrum)
bearTargets <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210331_7mer_NoStrict/mapped_to_brown_bear/summaryTables/ALLINTERVALS.brown_bear.mutyper.targets.7mer.nostrict.txt",header=T)
head(bearTargets)
colnames(bearTargets) <- c("ancestral7mer","targetCount_allIntervals")

bearAll <- merge(bearSpectrum,bearTargets,by="ancestral7mer")

bearAll$mutationRate_unscaled <- bearAll$totalSitesAllIntervals / bearAll$targetCount_allIntervals


bearAll <- bearAll %>%
  group_by(population) %>%
  mutate(mutationRate_rescaled = mutationRate_unscaled/sum(mutationRate_unscaled))

head(bearAll)


bearAll_wide = pivot_wider(bearAll[,c("population","variable","mutationRate_rescaled")],id_cols = c("variable"),values_from = "mutationRate_rescaled",names_from = "population")
head(bearAll_wide)

########## okay let's do some plotting ###########
merge1 <- merge(mouseAll_wide,whale1All_wide,by="variable")
head(merge1)

merge2 <- merge(merge1,whale2All_wide,by="variable")

merge3 <- merge(merge2,bearAll_wide,by="variable")
# can add more merges as needed
head(merge3)

mergeFinal <- merge3

#species=c("Mmd","Ms","EUR","PB","bacu2","ttru1") # don't restrict to species because want all mice in there, why not. 


########### want to get all pairwise ratios ########
# https://rdrr.io/cran/SARP.compo/man/rapports.html
??rapports
#species list is the list of columns to ratio 

# noms = col names to be used for ratios 
# don't use log = T (only if values are already log scaled in which case you'd take diff rather than ratio ) but that's not what I'm doing
# "These variables have their name constructed as the concatenation of the names of the two variables used, the first one being at the numerator, separated with a dot and with the additional suffix .r (or .r.log is working on difference of logarithms).

#Their order is determined by the order given in noms: the first variable of the list, V1, is used to compute ratios with all others (V1/V2, V1/V3 and so on). Then the second one is used for ratios further ones (V2/V3 and so on), and so on until the last one."
# note new df keeps original columns too, that's hand

oddsRatios_all <- calc.rapports(mergeFinal,noms=names(select(mergeFinal,-variable)))
# get rid of original columns
names(oddsRatios_all)
# 
oddsRatios_all <- select(oddsRatios_all,variable,ends_with(".r")) # just keep the ratio columns
names(oddsRatios_all)
# so Ms.bacu2.r is OR of Ms/bacu

# okay now I need rankings (or can I just compare odds ratios and look for outliers?)
# want to come up with a stat where if an OR differs in one comparison then it's notable

oddsRatios_all_melt <- melt(oddsRatios_all,variable.name = "speciesComparison")
head(oddsRatios_all_melt)
# need ranks https://rdrr.io/cran/matrixStats/man/rowRanks.html

# look for OR outliers (that are OR outliers in all comparisons or just in one?) 
# how to detect 
### figure out how to deal with ORs where one species has mutation rate of 0 (!!)
ggplot(oddsRatios_all_melt[oddsRatios_all_melt$value!=Inf,],aes(x=speciesComparison,y=value,group=variable))+
  geom_line()


####### rough stat: set an epsilon, if a row has all ORs with epsilon of 1, then it's labeled stable, but if at least one goes +- episolon then it's unstable in at least one species ... let's try it ########
epsilon = .25
stable <- oddsRatios_all_melt %>%
  group_by(variable) %>%
  summarise(countOfOutliers=(sum((value - 1) > epsilon)))  %>% 
  filter(countOfOutliers==0)


# tehse 'stable mutaiton types don't have any pairwise comparisons that have an outlier
stable$centralMutationType <- paste0(substr(stable$variable,4,4),".",substr(stable$variable,12,12))

ggplot(stable,aes(x=centralMutationType))+
  geom_histogram(stat="count")+
  ggtitle("Mutation Types with 0 OR outliers across pairwise comparisons")

oddsRatios_all[oddsRatios_all$variable %in% stable$variable,]
# counting up number of outliers per comparison 



######### need to get ranks and then look for rank outliers ###########
OR_rankings <- data.frame(sapply(select(oddsRatios_all,-variable),rank)) # ranking default way of dealing with ties is to give ties their average rank (I think that's good?) or you can replace their rank with their min rank or their max rank for ties .... 
# for now going to keep ties as average of tied rankings 
# OKAY note that rank orders from low --> high by default!
# this keeps order 
head(OR_rankings) # note some rankings end in .0 because averages are taken for ties
# need to think through rankings a bit more because rank 1 -2 is highest OR but 
# lowest ranks are ALSO interesting because those are low ORs.
# want a way to rank things close to one.

### SO I think I want to either work with ORs -or- with ranks, but not both. ####

# but isn't there another issue with rankings ORs that something could switch from low to high due to what's in the numerator?? hm. need to thnk this through more.

# so what If I didn't do ORs but I did rank. Rank on mutation rate?
# okay I think I either need to do ranks or OR.

# okay KH suggests a ranking where ORs close to 1/50th percentile as boring


# so ORs that are ranked far apart are interesting, but those near the 50th percentile (of all ORs? or the ones in the row?) are not? 

# idea: rank within central mut type categories?

############ rank mutations within central mutation type categories ########
mergeFinal$centralMutationType <- paste0(substr(mergeFinal$variable,4,4),".",substr(mergeFinal$variable,12,12))
head(mergeFinal)

mutationRate_rankings <- sapply(select(mergeFinal,-c(variable,centralMutationType)),rank) # ranking default way of dealing with ties is to give ties their average rank (I think that's good?) or you can replace their rank with their min rank or their max rank for ties .... 
head(mutationRate_rankings)
colnames(mutationRate_rankings) <- paste0(colnames(mutationRate_rankings),".rank") # which direction is this ranking? Is 1 high or low?
mutationRate_rankings <- cbind(mutationRate_rankings,mergeFinal) # CAREFUL HERE
head(mutationRate_rankings)
# okay so this is ranked from low to high 

### OKAY I AM JUST STUCK AND NEED TO THINK THIS THROUGH BETTER###

## make some plots? 
#### figure out some code to rank per central mut type
################ percentiles ############### #
testQuants <- quantile(oddsRatios_all$Mmc.Mmm.r  ,probs=c(.05,.95),names = F) # get 25th and 75th percentile cutoff values (so interested in <the 25 cutoff and >75, everything in between is boring)
# can set to names = F to just get quants back
testQuants

testOUtlierORs_MmcMmm <- oddsRatios_all[oddsRatios_all$Mmc.Mmm.r < testQuants[1] | oddsRatios_all$Mmc.Mmm.r >testQuants[2] ,c("variable","Mmc.Mmm.r")]


testQuants2 <- quantile(oddsRatios_all$EUR.PB.r,probs=c(.05,.95),names = F) # get 25th and 75th percentile cutoff values (so interested in <the 25 cutoff and >75, everything in between is boring)
testQuants2
# AHA this is the problem: need to use different cutoffs for bears: testQuants2 (BUG!)
testOUtlierORs_EurPB <- oddsRatios_all[oddsRatios_all$EUR.PB.r<testQuants2[1] | oddsRatios_all$EUR.PB.r >testQuants2[2] ,c("variable","EUR.PB.r")]
testOUtlierORs_EurPB

length(intersect(testOUtlierORs_MmdMs$variable,testOUtlierORs_EurPB$variable))
# okay 475 are the same 1839 are different length

(setdiff(testOUtlierORs_MmdMs$variable,testOUtlierORs_EurPB$variable))

#### okay there was a bug here where I wasn't using the correct cutoffs for the bears. (I was doing this during lab mtg so was distracted)
# once I use correct cutofs, it's more that are different than are overlapping! 
# So is NOT enriched for overlap. 


###########