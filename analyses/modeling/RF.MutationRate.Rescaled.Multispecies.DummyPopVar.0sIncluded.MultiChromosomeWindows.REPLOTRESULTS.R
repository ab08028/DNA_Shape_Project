########## Replotting some stuff from this model:
# 20210622_RF.MutationRate.Rescaled.Multispecies.DummyPopVar.0sIncluded.MultiChromosomeWindows
wd="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210622_RF.MutationRate.Rescaled.Multispecies.DummyPopVar.0sIncluded.MultiChromosomeWindows/"

model <- readRDS(paste0(wd,"modelTrainedOnOneFold.rds"))

require(scales)
require(tidymodels)
tidymodels_prefer()
# https://stackoverflow.com/questions/40219639/how-to-deal-with-zero-in-log-plot 
# ggplot has a pseudo log transform that transitions to linear scale around 0 
truth_prediction_df <- readRDS(paste0(wd,"modelTrainedOnOneFold.PREDICTIONS.onChr7.rds"))
truth_prediction_df$centralMutationType <- paste0(substr(truth_prediction_df$mutationType,4,4),".",substr(truth_prediction_df$mutationType,12,12))

rand_forest_Fold01_fit_predictions_plot <-  ggplot(truth_prediction_df, aes(y=.pred,x=outcome,color=centralMutationType,shape=population))+
  geom_point()+
  geom_abline()+
  facet_wrap(~population~centralMutationType,scales="free")+
  scale_x_log10()+
  scale_y_log10()+ 
  #scale_y_continuous(breaks=c(1e-6,1e-5,1e-4,1e-3,1e-2))+
  #scale_x_continuous(breaks=c(1e-6,1e-5,1e-4,1e-3,1e-2))+
  #ggtitle(paste0(description,"\ntrained on Fold01\n(all odd chrs but one, tested on just chr",toString(unique(truth_prediction_df$newGroup)),")"))+
  theme_bw()
rand_forest_Fold01_fit_predictions_plot



rsqsPerSpeciesAndMutationType <- truth_prediction_df %>%
  group_by(centralMutationType,population) %>%
  rsq(truth=outcome,estimate=.pred)
rsqsPerSpeciesAndMutationType

####  add rsq values to plot if possible ###
epsilon=1e-6
truth_prediction_df$.pred_ADDINEPSILON <- truth_prediction_df$.pred + epsilon
truth_prediction_df$outcome_ADDINEPSILON <- truth_prediction_df$outcome + epsilon

rand_forest_Fold01_fit_predictions_plot_faceted <-  ggplot(truth_prediction_df, aes(y=.pred_ADDINEPSILON,x=outcome_ADDINEPSILON,color=centralMutationType))+
  geom_point()+
  geom_abline()+
  geom_text(data=rsqsPerSpeciesAndMutationType,aes(x=1e-5,y=1e-4,label=round(.estimate,4)),color="black")+
  facet_grid(~centralMutationType~population,scales="free")+
  #ggtitle(paste0(description,"\ntrained on Fold01\n(all odd chrs but one, tested on just chr",toString(unique(truth_prediction_df$newGroup)),")"))+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  ggtitle(paste0("add in epsilon of ",epsilon," to enable log scaling (post-modeling)"))
rand_forest_Fold01_fit_predictions_plot_faceted
ggsave(paste0(wd,"addingInEpsilon",epsilon,"AfterTheFactOfTheModel.Faceted.LogScaled.png"),rand_forest_Fold01_fit_predictions_plot_faceted,height=12,width=7)



rand_forest_Fold01_fit_predictions_plot_faceted_noLog <-  ggplot(truth_prediction_df, aes(y=.pred,x=outcome,color=centralMutationType))+
  geom_point()+
  geom_abline()+
  geom_text(data=rsqsPerSpeciesAndMutationType,aes(x=1e-5,y=1e-4,label=round(.estimate,4)),color="black")+
  facet_grid(~centralMutationType~population,scales="free")+
  #ggtitle(paste0(description,"\ntrained on Fold01\n(all odd chrs but one, tested on just chr",toString(unique(truth_prediction_df$newGroup)),")"))+
  #scale_y_log10()+
  #scale_x_log10()+
  theme_bw()
rand_forest_Fold01_fit_predictions_plot_faceted_noLog
ggsave(paste0(wd,"plottedAfterTheFact.facetedByCentralMutation.NOTLOGSCALED.NOEPSILON.png"),rand_forest_Fold01_fit_predictions_plot_faceted_noLog,height=12,width=7)
