require(tidyr)
require(ggplot2)
require(tidymodels)
######## want to plot species comparision faceted :
outdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210625_model_009_RF.plusSequenceFeats.MutationRate.Rescaled.Multispecies.DummyPopVar.0sIncluded.MultiChromosomeWindows/"
truth_prediction_df <- readRDS(paste0(outdir,"modelTrainedOnOneFold.PREDICTIONS.onWindow7.rds"))

truth_prediction_df$centralMutationType <- paste0(substr(truth_prediction_df$mutationType,4,4),".",substr(truth_prediction_df$mutationType,12,12))

truth_prediction_df$ancestral_central3mer <- substr(truth_prediction_df$mutationType,start=3,stop=5)
truth_prediction_df$centralCpGLabel <- "NoCentralCpG"
truth_prediction_df[grep('CG',truth_prediction_df$ancestral_central3mer),]$centralCpGLabel <- "YesCentralCpG"

truth_prediction_df$centralMutationType <- paste0(substr(truth_prediction_df$mutationType,4,4),".",substr(truth_prediction_df$mutationType,12,12))

# make labels: 
truth_prediction_df$mutationLabel <- paste0(truth_prediction_df$centralMutationType,"_",truth_prediction_df$centralCpGLabel)

truth_prediction_df_spread <- pivot_wider(truth_prediction_df[,c("mutationType","population","outcome",".pred","centralMutationType","mutationCount","mutationLabel")],id_cols = c(mutationType,"mutationLabel",population,centralMutationType),names_from=population,values_from=c(outcome,.pred,mutationCount)) 
# adding mutation count ^


head(truth_prediction_df_spread)
# this is really useful: plot the difference between the two species and how well the predictions do (species are off of y=x line because of diff muttion rates; model picks that up! super cool)
speciesComparisonPlot <- ggplot(truth_prediction_df_spread,aes(x=outcome_Mmd,y=outcome_Ms))+
  geom_point()+
  geom_point(aes(x=.pred_Mmd,y=.pred_Ms),shape=1,color="red")+
  geom_abline()+
  facet_wrap(~mutationLabel,scales="free")+
  theme_bw()
speciesComparisonPlot

ggsave(paste0(outdir,"modelTrainedOnOneFold.SpeciesXYComparison.AssessedOnWindow",toString(unique(truth_prediction_df$newGroup)),"FACETED.mutationLabel.png"),speciesComparisonPlot,height=5,width=7)
# facet by mutation type:
truth_prediction_df_spread$meanMutationCountBetweenSpeciesWithinWindow <- (truth_prediction_df_spread$mutationCount_Ms+ truth_prediction_df_spread$mutationCount_Mmd)/2

speciesComparisonPlot2 <- ggplot(truth_prediction_df_spread,aes(x=outcome_Mmd,y=outcome_Ms,color=log10(meanMutationCountBetweenSpeciesWithinWindow)))+
  geom_point()+
  geom_point(aes(x=.pred_Mmd,y=.pred_Ms),shape=1,color="red")+
  geom_abline()+
  facet_wrap(~mutationLabel,scales="free")+
  scale_color_viridis_c()+
  theme_bw()
speciesComparisonPlot2
ggsave(paste0(outdir,"modelTrainedOnOneFold.SpeciesXYComparison.AssessedOnWindow",toString(unique(truth_prediction_df$newGroup)),"COLOREDBYABUNDANCE.mutationLabel.png"),speciesComparisonPlot2,height=5,width=12)


########3 want to replot model fit with free scales ##########

require(tidymodels)
rsqsPerSpeciesAndMutationType <- truth_prediction_df %>%
  group_by(mutationLabel,population) %>%
  rsq(truth=outcome,estimate=.pred)
rsqsPerSpeciesAndMutationType 
### need to add a false 'outcome' to the rsq to serve as its position

rsqsPerSpeciesAndMutationType$xPosition = 1e-4
rsqsPerSpeciesAndMutationType$yPosition = 1e-4
rsqsPerSpeciesAndMutationType[rsqsPerSpeciesAndMutationType$mutationLabel=="C.G_NoCentralCpG",]$xPosition = 1.5e-5
rsqsPerSpeciesAndMutationType[rsqsPerSpeciesAndMutationType$mutationLabel=="C.G_NoCentralCpG",]$yPosition = 2e-5
rsqsPerSpeciesAndMutationType[rsqsPerSpeciesAndMutationType$mutationLabel=="A.C_NoCentralCpG",]$xPosition = 7e-5
# add a position for each rsq?

rand_forest_Fold01_fit_predictions_plot_faceted <-  ggplot(truth_prediction_df, aes(y=.pred,x=outcome,color=mutationLabel))+
  geom_point()+
  geom_abline()+
  geom_text(data=rsqsPerSpeciesAndMutationType,aes(x=xPosition,y=yPosition,label=round(.estimate,4)),color="black",position = position_stack(vjust = 0.5))+
  facet_wrap(~population~mutationLabel,scales="free",ncol = 2,dir = "v")+
  #ggtitle(paste0(description,"\ntrained on Fold01\n(all odd Windows but one, tested on just Window",toString(unique(truth_prediction_df$newGroup)),")"))+
  theme_bw()+
  theme(legend.position = "none")
rand_forest_Fold01_fit_predictions_plot_faceted
ggsave(paste0(outdir,"modelTrainedOnOneFold.PredictionsPlot.FacetedPerMutationType.AssessedOnWindow",toString(unique(truth_prediction_df$newGroup)),"FACETWRAPFREESCALES.png"),rand_forest_Fold01_fit_predictions_plot_faceted,height=16,width=7)






# log10 -scale (it keeps in 0 entries but places them at far left - ok?)
rand_forest_Fold01_fit_predictions_plot_faceted_LOG10 <-  ggplot(truth_prediction_df, aes(y=.pred,x=outcome,color=mutationLabel))+
  geom_point()+
  geom_abline()+
  geom_text(data=rsqsPerSpeciesAndMutationType,aes(x=xPosition,y=yPosition,label=round(.estimate,4)),color="black",position = position_stack(vjust = 0.5))+
  facet_wrap(~population~mutationLabel,scales="free",ncol = 2,dir = "v")+
  #ggtitle(paste0(description,"\ntrained on Fold01\n(all odd Windows but one, tested on just Window",toString(unique(truth_prediction_df$newGroup)),")"))+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_log10()+
  scale_y_log10()
rand_forest_Fold01_fit_predictions_plot_faceted_LOG10
ggsave(paste0(outdir,"modelTrainedOnOneFold.PredictionsPlot.FacetedPerMutationType.AssessedOnWindow",toString(unique(truth_prediction_df$newGroup)),"FACETWRAPFREESCALES.LOG10.png"),rand_forest_Fold01_fit_predictions_plot_faceted_LOG10,height=16,width=7)

#plot without facets:
rand_forest_Fold01_fit_predictions_plot <-  ggplot(truth_prediction_df, aes(y=.pred,x=outcome,color=mutationLabel))+
  geom_point()+
  geom_abline()+
  theme_bw()+
  facet_wrap(~population)
rand_forest_Fold01_fit_predictions_plot
ggsave(paste0(outdir,"modelTrainedOnOneFold.PredictionsPlot.NotFaceted.AssessedOnWindow",toString(unique(truth_prediction_df$newGroup)),".png"),rand_forest_Fold01_fit_predictions_plot,height=8,width=12)


# log scaled -- lose the 0 entry mutations though
rand_forest_Fold01_fit_predictions_plot_log10 <-  ggplot(truth_prediction_df, aes(y=.pred,x=outcome,color=mutationLabel))+
  geom_point()+
  geom_abline()+
  theme_bw()+
  facet_wrap(~population)+
  scale_y_log10()+
  scale_x_log10()
rand_forest_Fold01_fit_predictions_plot_log10
ggsave(paste0(outdir,"modelTrainedOnOneFold.PredictionsPlot.NotFaceted.AssessedOnWindow",toString(unique(truth_prediction_df$newGroup)),"log10Scaled.KeepsNAsIn.png"),rand_forest_Fold01_fit_predictions_plot_log10,height=8,width=12)

########### try separataing by 5mer and looking at results #########
# want to see a cleaner model for 5mers without noise of 7mers
truth_prediction_df$central5mer <-  paste0(substr(truth_prediction_df$mutationType,2,6),".",substr(truth_prediction_df$mutationType,10,14))

# want to group by central 5mer and recalculate things 
central5mer_prediction_df <- truth_prediction_df %>%
  group_by(population,newGroup,central5mer,mutationLabel) %>% 
  summarise(.pred_SumOVerCentral5mer=sum(.pred),outcome_SumOverCentral5mer=sum(outcome))
central5mer_prediction_df
  #summarise(.pred_AvgOverCentral5mer=mean(.pred),outcome_AvgOverCentral5mer=mean(outcome))
# don't want average, want sum of 7mers within 5mer! this is wrong. 
### get rsq of 5mers:
rsqsPerSpeciesAndMutationType_central5mer <- central5mer_prediction_df %>%
  group_by(mutationLabel,population) %>%
  rsq(truth=outcome_SumOverCentral5mer,estimate=.pred_SumOVerCentral5mer)
rsqsPerSpeciesAndMutationType_central5mer 

## add in locations 

rsqsPerSpeciesAndMutationType_central5mer$xPosition = 5e-4
rsqsPerSpeciesAndMutationType_central5mer$yPosition = 5e-4
rsqsPerSpeciesAndMutationType_central5mer[rsqsPerSpeciesAndMutationType_central5mer$mutationLabel=="C.T_YesCentralCpG",]$xPosition = 5e-3
rsqsPerSpeciesAndMutationType_central5mer[rsqsPerSpeciesAndMutationType_central5mer$mutationLabel=="C.T_YesCentralCpG",]$yPosition = 5e-3
rsqsPerSpeciesAndMutationType_central5mer[rsqsPerSpeciesAndMutationType_central5mer$mutationLabel=="C.T_NoCentralCpG",]$xPosition = .001
rsqsPerSpeciesAndMutationType_central5mer[rsqsPerSpeciesAndMutationType_central5mer$mutationLabel=="C.T_NoCentralCpG",]$yPosition= 0.001
rsqsPerSpeciesAndMutationType_central5mer[rsqsPerSpeciesAndMutationType_central5mer$mutationLabel=="C.G_NoCentralCpG",]$xPosition = 2e-4
rsqsPerSpeciesAndMutationType_central5mer[rsqsPerSpeciesAndMutationType_central5mer$mutationLabel=="A.C_NoCentralCpG",]$xPosition = 3e-4
# add a position for each rsq?

dim(central5mer_prediction_df)

central5merplot <- ggplot(central5mer_prediction_df,aes(x=outcome_SumOverCentral5mer,y=.pred_SumOVerCentral5mer,color=mutationLabel))+
  facet_wrap(~population~mutationLabel,scales="free",ncol=2,dir="v")+
  geom_point()+
  geom_abline()+
  #scale_y_log10()+
  #scale_x_log10()+
  ggtitle("summed up 7mers over central 5mer")+
  geom_text(data=rsqsPerSpeciesAndMutationType_central5mer,aes(x=xPosition,y=yPosition,label=round(.estimate,4)),color="black",position = position_stack(vjust = 0.5))+
  theme_bw()
central5merplot
ggsave(paste0(outdir,"modelTrainedOnOneFold.PredictionsPlot.NotFaceted.AssessedOnWindow",toString(unique(central5mer_prediction_df$newGroup)),".SUMMEDOVERCENTRAL5mer.png"),central5merplot,height=12,width=10)

### compare between species #########
truth_prediction_df_spread <- pivot_wider(central5mer_prediction_df,id_cols = c(central5mer,mutationLabel,newGroup,population),names_from=population,values_from=c(outcome_SumOverCentral5mer,.pred_SumOVerCentral5mer)) 
head(truth_prediction_df_spread)

speciesComparison_5mers <- ggplot(truth_prediction_df_spread,aes(x=outcome_SumOverCentral5mer_Mmd,y=outcome_SumOverCentral5mer_Ms))+
  geom_point()+
  geom_point(aes(x=.pred_SumOVerCentral5mer_Mmd,y=.pred_SumOVerCentral5mer_Ms),shape=1,color="red")+
  facet_wrap(~mutationLabel,scales="free")+
  geom_abline()+
  theme_bw()
speciesComparison_5mers
ggsave(paste0(outdir,"modelTrainedOnOneFold.speciesComparison_5mers.AssessedOnWindow",toString(unique(central5mer_prediction_df$newGroup)),".SUMMEDOVERCENTRAL5mer.png"),speciesComparison_5mers,height=12,width=12)

#### compare 3mers ##########
truth_prediction_df$central3mer <-  paste0(substr(truth_prediction_df$mutationType,3,5),".",substr(truth_prediction_df$mutationType,11,13))
tail(truth_prediction_df[,c('mutationType',"central3mer","central5mer")])

  
# sum up 7mer predictions per central 3mer: 
central3mer_prediction_df <- truth_prediction_df %>%
  group_by(population,newGroup,central3mer,mutationLabel) %>% 
  summarise(.pred_SumOVerCentral3mer=sum(.pred),outcome_SumOverCentral3mer=sum(outcome))
central3mer_prediction_df
# check:
sum(central3mer_prediction_df[central3mer_prediction_df$population=="Mmd",]$outcome_SumOverCentral3mer) # cool sums to 1

rsqsPerSpeciesAndMutationType_central3mer <- central3mer_prediction_df %>%
  group_by(mutationLabel,population) %>%
  rsq(truth=outcome_SumOverCentral3mer,estimate=.pred_SumOVerCentral3mer)
rsqsPerSpeciesAndMutationType_central3mer 


central3merplot <- ggplot(central3mer_prediction_df,aes(x=outcome_SumOverCentral3mer,y=.pred_SumOVerCentral3mer,color=mutationLabel))+
  facet_wrap(~population~mutationLabel,scales="free",ncol=2,dir="v")+
  geom_point()+
  geom_abline()+
  #scale_y_log10()+
  #scale_x_log10()+
  ggtitle("summed up 7mers over central 5mer")+
  geom_text(data=rsqsPerSpeciesAndMutationType_central3mer,aes(x=1e-3,y=1e-3,label=round(.estimate,4)),color="black")+
  theme_bw()
central3merplot
ggsave(paste0(outdir,"modelTrainedOnOneFold.PredictionsPlot.NotFaceted.AssessedOnWindow",toString(unique(central5mer_prediction_df$newGroup)),".SUMMEDOVERCENTRAL3mer.png"),central3merplot,height=12,width=10)

### compare between species #########
truth_prediction_df_spread <- pivot_wider(central3mer_prediction_df,id_cols = c(central3mer,mutationLabel,newGroup,population),names_from=population,values_from=c(outcome_SumOverCentral3mer,.pred_SumOVerCentral3mer)) 
head(truth_prediction_df_spread)

speciesComparison_3mers <- ggplot(truth_prediction_df_spread,aes(x=outcome_SumOverCentral3mer_Mmd,y=outcome_SumOverCentral3mer_Ms))+
  geom_point()+
  geom_point(aes(x=.pred_SumOVerCentral3mer_Mmd,y=.pred_SumOVerCentral3mer_Ms),shape=1,color="red")+
  facet_wrap(~mutationLabel,scales="free")+
  geom_abline()+
  theme_bw()
speciesComparison_3mers
ggsave(paste0(outdir,"modelTrainedOnOneFold.speciesComparison_3mers.AssessedOnWindow",toString(unique(central3mer_prediction_df$newGroup)),".SUMMEDOVERCENTRAL3mer.png"),speciesComparison_3mers,height=12,width=12)
