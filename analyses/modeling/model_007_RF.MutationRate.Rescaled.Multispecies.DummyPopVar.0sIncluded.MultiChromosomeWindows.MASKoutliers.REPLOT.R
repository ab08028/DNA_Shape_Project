####### replotting from:
wd="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210622_RF.MutationRate.Rescaled.Multispecies.DummyPopVar.0sIncluded.MultiChromosomeWindows.MASKoutliers/"

predictiondf <- readRDS(paste0(wd,"modelTrainedOnOneFold.PREDICTIONS.onChr7.rds"))
predictiondf

head(predictiondf)
predictiondf$centralMutationType <- paste0(substr(predictiondf$mutationType,4,4),".",substr(predictiondf$mutationType,12,12))

# add in tiny epsilon AFTER the model fact 
epsilon=1e-6
predictiondf$.pred_PLUSEPSILON <- predictiondf$.pred+epsilon
predictiondf$outcome_PLUSEPSILON <- predictiondf$outcome + epsilon

outcomePlusEpsilonlog10plot <- ggplot(predictiondf,aes(x=outcome_PLUSEPSILON,y=.pred_PLUSEPSILON,color=centralMutationType))+
  geom_point()+
  geom_abline(slope=1,intercept=0)+
  facet_wrap(~centralMutationType~population,scales="free")+
  scale_x_log10()+
  scale_y_log10()+
  ggtitle(paste0("added epsilon = ",epsilon, " after the fact to plot on log10 scale to all points"))
outcomePlusEpsilonlog10plot
ggsave(paste0(wd,"logScale.FacetedByCentralMutationType.PLUSepsilon",epsilon,"AfterthefactForPlottingOnly.png"),outcomePlusEpsilonlog10plot,height=12,width=12)
