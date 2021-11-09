require(tidyr)
require(ggplot2)
require(tidymodels)

######## want to plot species comparision faceted :
outdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210719_model_020_RF.plusSequenceFeats.MutationRate.Rescaled.Multispecies.DummyPopVar.SingleChrWindows.FIXED0ENTRYPROBLEM/"
truth_prediction_df <- readRDS(paste0(outdir,"modelTrainedOnOneFold.PREDICTIONS.onWindow7.rds"))

truth_prediction_df$centralMutationType <- paste0(substr(truth_prediction_df$mutationType,4,4),".",substr(truth_prediction_df$mutationType,12,12))

truth_prediction_df$ancestral_central3mer <- substr(truth_prediction_df$mutationType,start=3,stop=5)
truth_prediction_df$centralCpGLabel <- ""

# okay I want to change how I do this from how I was doing it in model 009
# I want to do it with dimers instead 
truth_prediction_df$ancestralDimer <- substr(truth_prediction_df$ancestral7mer,4,5)
  
# only want to label a CpG if the center bp is a C
# so CCG would be a _CpG but CGG would not be. 
# previously I was just labeling any CpG containing central 3mer as a CpG.
# but for now want to do it this way

# old way I did it: truth_prediction_df[grep('CG',truth_prediction_df$ancestral_central3mer),]$centralCpGLabel <- "_CpG"

truth_prediction_df[truth_prediction_df$ancestralDimer=="CG",]$centralCpGLabel <- "CpG_"

head(truth_prediction_df[,c("mutationType","centralCpGLabel")])

truth_prediction_df$MutationTypePlusCpGLabel <- paste0(truth_prediction_df$centralCpGLabel,gsub("\\.",">",truth_prediction_df$centralMutationType))

#View(truth_prediction_df[,c("mutationType","MutationTypePlusCpGLabel")])
#tail(truth_prediction_df[,c("mutationType","MutationTypePlusCpGLabel")])

truth_prediction_df_spread <- pivot_wider(truth_prediction_df[,c("mutationType","population","outcome",".pred","centralMutationType","mutationCount","MutationTypePlusCpGLabel")],id_cols = c(mutationType,MutationTypePlusCpGLabel,population,centralMutationType),names_from=population,values_from=c(outcome,.pred,mutationCount)) 

head(truth_prediction_df_spread)

speciesComparisonPlot <- ggplot(truth_prediction_df_spread,aes(x=outcome_Mmd,y=outcome_Ms))+
  geom_point()+
  geom_point(aes(x=.pred_Mmd,y=.pred_Ms),shape=1,color="red")+
  geom_abline()+
  facet_wrap(~MutationTypePlusCpGLabel,scales="free")+
  theme_bw()+
  xlab("Mmd rescaled mutation rate")+
  ylab("Ms rescaled mutation rate")
speciesComparisonPlot

ggsave(paste0(outdir,"modelTrainedOnOneFold.SpeciesXYComparison.AssessedOnWindow",toString(unique(truth_prediction_df$newGroup)),"FACETED.mutationLabel.png"),speciesComparisonPlot,height=5,width=7)

overallModelPlot1 <- ggplot(truth_prediction_df,aes(x=outcome,y=.pred,color=MutationTypePlusCpGLabel))+
  geom_point(alpha=0.5)+
  theme_bw()+
  geom_abline(slope=1,intercept=0)+
  facet_wrap(~population,ncol=1)+
  xlab("truth")+
  ylab("model prediction")
overallModelPlot1

ggsave(paste0(outdir,"ModelFit",toString(unique(truth_prediction_df$newGroup)),".allMutTypes.png"),overallModelPlot1,height=3,width=5)

overallModelPlot2 <- ggplot(truth_prediction_df,aes(x=outcome,y=.pred,color=MutationTypePlusCpGLabel))+
  geom_point(alpha=0.5)+
  theme_bw()+
  geom_abline(slope=1,intercept=0)+
  facet_wrap(~population,ncol=1)+
  scale_x_log10()+
  scale_y_log10()+
  xlab("truth (log10 scaled)")+
  ylab("model prediction (log10 scaled)")
overallModelPlot2

ggsave(paste0(outdir,"ModelFit",toString(unique(truth_prediction_df$newGroup)),".allMutTypes.log10.png"),overallModelPlot2,height=3,width=5)
