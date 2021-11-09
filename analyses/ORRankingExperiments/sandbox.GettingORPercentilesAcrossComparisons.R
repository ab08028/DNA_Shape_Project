require(dplyr)
require(reshape2)
require(tidyr)
require(GGally)
require(SARP.compo) # for getting all ratios
plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/rankingORs/sandbox/"
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
#ggplot(oddsRatios_all_melt[oddsRatios_all_melt$value!=Inf,],aes(x=speciesComparison,y=value,group=variable))+
#  geom_line()

################ percentiles ###############
# focus on disjoint comparisons mouse-mouse bear-bear whale-whale
disjointComparisons=c("Mmd.Ms.r","EUR.PB.r","bacu2.ttru1.r")
# want to restrict just to these


oddsRatios_disjointCompsOnly <- oddsRatios_all_melt[oddsRatios_all_melt$speciesComparison %in% disjointComparisons,]

head(oddsRatios_disjointCompsOnly)

lowPercentile=.05
highPercentile=.95

# adding in quant low high values (same for all 7mers within each species comparison)
oddsRatios_disjointCompsOnly <- oddsRatios_disjointCompsOnly %>%
  group_by(speciesComparison) %>%
  mutate(quantLow = quantile(value,probs=lowPercentile),quantHigh=quantile(value,probs=highPercentile))

head(oddsRatios_disjointCompsOnly)

# want to label 7mers that are outside quant range
oddsRatios_disjointCompsOnly$quantLabel <- "withinIQR"
# label those that are below or above the percentils
oddsRatios_disjointCompsOnly[oddsRatios_disjointCompsOnly$value < oddsRatios_disjointCompsOnly$quantLow | oddsRatios_disjointCompsOnly$value > oddsRatios_disjointCompsOnly$quantHigh, ]$quantLabel <- "outsideIQR"


head(oddsRatios_disjointCompsOnly)

# pivot it wider:
oddsRatios_disjointCompsOnly_wider <- pivot_wider(oddsRatios_disjointCompsOnly,id_cols = c(variable),names_from=speciesComparison,values_from=quantLabel)

tail(oddsRatios_disjointCompsOnly_wider)



OutsideIQRCount <- oddsRatios_disjointCompsOnly %>%
  group_by(variable) %>% 
  # condition on being an outlier in at least one 
  summarise(outlierInXComparisons=sum(quantLabel=="outsideIQR"))

head(OutsideIQRCount)
head(OutsideIQRCount[OutsideIQRCount$outlierInXComparisons>0,])

ggplot(OutsideIQRCount[OutsideIQRCount$outlierInXComparisons>0,],aes(x=outlierInXComparisons))+
  geom_histogram(stat="count")+
  ggtitle(paste0("(outliers are <",round(lowPercentile*100),"th percentile and > ",round(highPercentile*100),"th percentile)"))

# want to be able to generalize to more comparisons

# want to look at overlap

# want to make a heatmap with shapes

# before when I was doing Mmc and Mmm I got way more overlap with bears whats going on there 
######### make a heatmap showing shape features ########
#shape features

# for testing on home computer: 
#shapes <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/METHYLATEDfirstOrder_featureTypes_allPossible7mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.METHYLATED.txt",header=T,sep="\t") # instead of using type with units want the min-max normalized for the heatmap

shapes <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/firstOrder_featureTypes_allPossible7mers.FeatureValues.min-maxNormalized.WillWorkForAnySpecies.goodForPCAOrAnythingWhereDontWantUnits.txt",header=T,sep="\t")
rownames(shapes) <- shapes$motif

### add in shapes ####
# all 7mers just appear once for eah species. no train/test chromosomes or even/odd bp.
oddsRatios_disjointCompsOnly_wider$ancestral7mer <- substr(oddsRatios_disjointCompsOnly_wider$variable,1,7)
oddsRatios_disjointCompsOnly_wider$derived7mer <- substr(oddsRatios_disjointCompsOnly_wider$variable,9,15)

head(oddsRatios_disjointCompsOnly_wider)

oddsRatios_disjointCompsOnly_wider_combinedWithShapes_intermediate <-merge(oddsRatios_disjointCompsOnly_wider,shapes,by.x="ancestral7mer",by.y="motif")

oddsRatios_disjointCompsOnly_wider_combinedWithShapes <- merge(oddsRatios_disjointCompsOnly_wider_combinedWithShapes_intermediate,shapes,by.x="derived7mer",by.y="motif",suffixes=c(".ancestral",".derived"))


head(oddsRatios_disjointCompsOnly_wider_combinedWithShapes)



######## add seq info (not yet 1-hot encoded, will encode with tidy models )##########
# just want ancestral 7mer 
oddsRatios_disjointCompsOnly_wider_combinedWithShapes <- oddsRatios_disjointCompsOnly_wider_combinedWithShapes %>% separate(ancestral7mer,into=c("Pos_1","Pos_2","Pos_3","Pos_4.ancestral","Pos_5","Pos_6","Pos_7"),remove=F,sep=c(1,2,3,4,5,6))
# also get derived bp:
oddsRatios_disjointCompsOnly_wider_combinedWithShapes$Pos_4.derived <- substr(oddsRatios_disjointCompsOnly_wider_combinedWithShapes$derived7mer,4,4)


head(oddsRatios_disjointCompsOnly_wider_combinedWithShapes)

#########  let's get ORs that are in at least 1, 2 or 3 outlier groups ########

# okay let's make a heatmap

# need a matrix

# ?heatmap
#3# add in count of outlier in X comparisons:
oddsRatios_disjointCompsOnly_wider_combinedWithShapes <- merge(oddsRatios_disjointCompsOnly_wider_combinedWithShapes,OutsideIQRCount,by="variable")

#ORoutliersWithFeatures_mat_outlierin3 <- oddsRatios_disjointCompsOnly_wider_combinedWithShapes %>%
#  filter(variable %in% OutsideIQRCount[OutsideIQRCount$outlierInXComparisons==3,]$variable) %>%
#  select(c(starts_with("feature_"),"variable"))  


#ORoutliersWithFeatures_mat_outlierin3
#rownames(ORoutliersWithFeatures_mat_outlierin3) <- ORoutliersWithFeatures_mat_outlierin3$variable
# need to 1-hot encode pos variables before I can use them
#heatmap(data.matrix(select(ORoutliersWithFeatures_mat_outlierin3,-variable)))

# use ggplot:
oddsRatios_disjointCompsOnly_wider_combinedWithShapes_melt <- melt(select(oddsRatios_disjointCompsOnly_wider_combinedWithShapes,c(starts_with("feature_"),"variable","outlierInXComparisons")),variable.name  ="feature",id.vars = c("variable","outlierInXComparisons"))

# just do features rather than other things 
heatmap1 <- ggplot(oddsRatios_disjointCompsOnly_wider_combinedWithShapes_melt,aes(x=variable,y=feature,fill=as.numeric(value)))+
  geom_tile()+
  facet_wrap(~as.factor(outlierInXComparisons),scales="free_x",nrow=1)+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  scale_fill_viridis_c()+
  ggtitle(paste0("Separating 7mers by whether they appear as OR outliers in 0, 1, 2, or 3 disjoint species comparisons\n(outliers are <",round(lowPercentile*100),"th percentile and > ",round(highPercentile*100),"th percentile)"))
heatmap1
ggsave(paste0(plotdir,"ORs.outliers.",round(lowPercentile*100),"thPercentile.",round(highPercentile*100),"thPercentile.heatmap.png"),heatmap1,width=12,height=7)
# note the melting turned thingsinto factors and that was slowing down
# another way to plot

# want to plot 7mers along the x axis and changes in features along y? 
# ggplot(oddsRatios_disjointCompsOnly_wider_combinedWithShapes_melt,aes(x=variable,y=value,color=feature,fill=feature,group=feature))+
#   geom_line()+
#   facet_wrap(~outlierInXComparisons,scales="free_x",nrow=1)+
#   theme(axis.ticks.x = element_blank(),
#         axis.text.x = element_blank())+
#   theme(legend.position="none")

#  geom_tile() #+
  #theme(axis.ticks = "none",axis.text = "none")+
  #facet_wrap(~as.factor(outlierInXComparisons),scales ="free_x")


######### see what is unique about the outliers ;; abundance? #####
# get raw counts in one species (?) or average abudnance (?) 

head(oddsRatios_disjointCompsOnly_wider)
head(mouseSpectrum)
oddsRatios_disjointCompsOnly_wider_withAbundance1 <- merge(oddsRatios_disjointCompsOnly_wider,mouseSpectrum[mouseSpectrum$population=="Mmd",c("variable","totalSitesAllIntervals")],by='variable')

oddsRatios_disjointCompsOnly_wider_withAbundance2 <- merge(oddsRatios_disjointCompsOnly_wider_withAbundance1,bearSpectrum[bearSpectrum$population=="EUR", c("variable","totalSitesAllIntervals")],by="variable",suffixes = c(".MmdTotalSites",".EURTotalSites"))

head(oddsRatios_disjointCompsOnly_wider_withAbundance2)
### leaving out whales for now just want to see if abundant

ggplot(oddsRatios_disjointCompsOnly_wider_withAbundance2,aes(x=totalSitesAllIntervals.MmdTotalSites,fill=Mmd.Ms.r))+
  geom_histogram(position="dodge",bins=100)+
  scale_x_log10()
######## YOU ARE HERE ########
