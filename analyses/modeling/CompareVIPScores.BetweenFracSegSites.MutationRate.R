######## compare the VIP scores of frac seg sites and mutation rate #######33
require(ggplot2)
require(dplyr)
require(reshape2)
fracSegSitesVIPScores=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210618_RF.FracSegSites.Log10.Multispecies.DummyPopVar/modelTrainedOnOneFold.VIPScores.PerMutatationType.AssessedOnChr7.txt",header=T)
head(fracSegSitesVIPScores)
# assign ranking (note are already in order but ranking to be careful)
fracSegSitesVIPScores <- fracSegSitesVIPScores %>% 
  mutate(rank = dense_rank(desc(Importance)))
head(fracSegSitesVIPScores)
#fracSegSitesVIPScores$label <- "FracSegSites"

mutationrateVIPScores=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210618_RF.MutationRate.Log10.Rescaled.Multispecies.DummyPopVar/modelTrainedOnOneFold.VIPScores.PerMutatationType.AssessedOnChr7.txt",header=T)
head(mutationrateVIPScores)

mutationrateVIPScores <- mutationrateVIPScores %>% 
  mutate(rank = dense_rank(desc(Importance)))
head(mutationrateVIPScores)
#mutationrateVIPScores$label <- "MutationRate"


merged <- merge(mutationrateVIPScores,fracSegSitesVIPScores,by="Variable",suffixes = c(".mutationRate",".fracSegSites"))
head(merged)

merged$variable_group <- unlist(lapply(strsplit(merged$Variable,"_"),"[",3))
merged$variable_shape <- unlist(lapply(strsplit(merged$Variable,"\\."),"[",2))

merged[merged$Variable=="population_Ms",]$variable_group <- "population"
merged[merged$Variable=="population_Ms",]$variable_shape <- "population"

# note that importance are in diff units because are diff measures of mut rate so only focus on rank
vipcomparisonplot <- ggplot(merged,aes(x=rank.fracSegSites,y=rank.mutationRate,color=variable_group,shape=variable_shape))+
  geom_point(size=2)+
  geom_abline(slope=1,intercept=0)+
  theme_bw()+
  ggtitle("Comparing Rank of frac seg sites and mutation rate outcomes\nVI Scores calculated based on permutation importance")+
  scale_shape_manual(values=c(1,16,8))
vipcomparisonplot
ggsave("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/comparisonPlots/comparingVIPScores.MutationRatelog10Rescaled.FracSegSites.png",vipcomparisonplot,height=5,width=8)

# top 10 for either one:
topScores <- subset(merged, rank.fracSegSites<=10 | rank.mutationRate <=10)$Variable
topScores 
# and need to get rid of .derived and .ancestral and
topScores <- gsub(".derived","",topScores)
topScores <- gsub(".ancestral","",topScores)
topScores <- unique(topScores)
topScores
###### let's look at shear values for each mutation type #########
shapes <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/firstOrder_featureTypes_allPossible7mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.txt",header=T)
head(shapes)
shapes$centralBP <- substr(shapes$motif,4,4)
head(shapes)
tail(shapes)

ggplot(shapes,aes(x=centralBP,y=feature_1_Shear_2))+
  geom_boxplot()

######## plot the top scores:
shapes_melt <- melt(shapes)
head(shapes_melt)

shapes_melt_subset_top <- subset(shapes_melt, variable %in% topScores)
dim(shapes_melt)
dim(shapes_melt_subset_top)

topScoresPlot <- ggplot(shapes_melt_subset_top,aes(x=centralBP,y=value))+
  facet_wrap(~variable,scales="free")+
  geom_boxplot()+
  theme_bw()+
  ggtitle("Top scores from frac seg sites and mutation rate models\nboxplot of 7mer shape values separated by central mutation type")
topScoresPlot
ggsave("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/comparisonPlots/comparingShapes.AcrossTopFeatures.png",topScoresPlot,height=5,width=8)

##### plot the bottom scores:
# bottom 10 for either one:
bottomScores <- subset(merged, rank.fracSegSites>=87 | rank.mutationRate >=87)$Variable
bottomScores 
# and need to get rid of .derived and .ancestral and
bottomScores <- gsub(".derived","",bottomScores)
bottomScores <- gsub(".ancestral","",bottomScores)
bottomScores <- unique(bottomScores)
bottomScores

shapes_melt_subset_bottom <- subset(shapes_melt, variable %in% bottomScores)

bottomScoresPlot <- ggplot(shapes_melt_subset_bottom,aes(x=centralBP,y=value))+
  facet_wrap(~variable,scales="free")+
  geom_boxplot()+
  theme_bw()+
  ggtitle("Bottom scores from frac seg sites and mutation rate models\nboxplot of 7mer shape values separated by central mutation type")
bottomScoresPlot
ggsave("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/comparisonPlots/comparingShapes.AcrossBottomFeatures.png",bottomScoresPlot,height=5,width=8)

#### plot all features:
allFeatures <- ggplot(shapes_melt,aes(x=centralBP,y=value))+
  facet_wrap(~variable,scales="free")+
  geom_boxplot()+
  theme_bw()
allFeatures
ggsave("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/comparisonPlots/comparingAllFeatures.PerCentralMutation.png",allFeatures,height=12,width=12)

