########## comparison plots between model 006 and 009 (adding in seq features ) ############
  # and label outliers 
# both are tested on chr 7
plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/comparisonPlots/"
model_006_predictions <- readRDS("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210622_model_006_RF.MutationRate.Rescaled.Multispecies.DummyPopVar.0sIncluded.MultiChromosomeWindows/modelTrainedOnOneFold.PREDICTIONS.onChr7.rds")
head(model_006_predictions)
model_006_predictions$model <- "model_006_ShapesOnly"

model_009_predictions <- readRDS("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210625_model_009_RF.plusSequenceFeats.MutationRate.Rescaled.Multispecies.DummyPopVar.0sIncluded.MultiChromosomeWindows/modelTrainedOnOneFold.PREDICTIONS.onWindow7.rds")
head(model_009_predictions)
model_009_predictions$model <- "model_009_SeqPlusShapes"


head(model_009_predictions)

combo <- rbind(model_009_predictions[,c("mutationType",".pred","outcome","population","model")],model_006_predictions[,c("mutationType",".pred","outcome","population","model")])

############ label CpGs ###################
combo$ancestral_central3mer <- substr(combo$mutationType,start=3,stop=5)
combo$centralCpGLabel <- "NoCentralCpG"
combo[grep('CG',combo$ancestral_central3mer),]$centralCpGLabel <- "YesCentralCpG"

combo$centralMutationType <- paste0(substr(combo$mutationType,4,4),".",substr(combo$mutationType,12,12))

# make labels: 
combo$mutationLabel <- paste0(combo$centralMutationType,"_",combo$centralCpGLabel)
# note this doesn't have the 1-coded pops unless you juice() it but rows are still in same order 


########### facet by central mutation type and get individual rsqs ##############
rsqsPerModelSpeciesAndMutationType <- combo %>%
  group_by(mutationLabel,population,model) %>%
  rsq(truth=outcome,estimate=.pred)
rsqsPerModelSpeciesAndMutationType


combo_predictionsPlot <-  ggplot(combo, aes(y=.pred,x=outcome,color=model,shape=population))+
  geom_point(alpha=0.6)+
  geom_abline()+
  facet_wrap(~mutationLabel~population,scales="free")+
  #ggtitle(paste0(description,"\ntrained on Fold01\n(all odd windows but one, tested on just window ",toString(unique(truth_prediction_df$newGroup)),")"))+
  theme_bw()+
  ggtitle("the sequence+shape and shape only models fit equivalently")+
  geom_text(data=rsqsPerModelSpeciesAndMutationType[rsqsPerModelSpeciesAndMutationType$model=="model_009_SeqPlusShapes",],aes(x=2e-5,y=1.2e-4,label=round(.estimate,4)))+
  geom_text(data=rsqsPerModelSpeciesAndMutationType[rsqsPerModelSpeciesAndMutationType$model=="model_006_ShapesOnly",],aes(x=2e-5,y=1e-4,label=round(.estimate,4)))

combo_predictionsPlot

ggsave(paste0(plotdir,"comparingFitsOFModel006_009_shapesOnlyvsSeqPlusShape.png"),combo_predictionsPlot,height=12,width=12)

########### Plot labeling outliers ###############

combo$residual_abs <- abs(combo$.pred - combo$outcome)

meanResidualsPerGroup <- combo %>%
  group_by(model,population,mutationLabel) %>%
  summarise(mean_residual=mean(residual_abs))

combo <- merge(combo,meanResidualsPerGroup,by=c("model","population","mutationLabel"))
head(combo)

combo$fracOfMeanResidual <- abs((combo$residual - combo$mean_residual) / combo$mean_residual)

#combo %>%
#  subset((residual - mean_residual)/mean_residual > 0.5)

combo_predictionsPlot_labelOutliers <-  ggplot(combo, aes(y=.pred,x=outcome,color=model,shape=population))+
  geom_point(alpha=0.6)+
  geom_abline()+
  facet_wrap(~mutationLabel~population,scales="free")+
  #ggtitle(paste0(description,"\ntrained on Fold01\n(all odd windows but one, tested on just window ",toString(unique(truth_prediction_df$newGroup)),")"))+
  theme_bw()+
  geom_text_repel(data=combo[combo$fracOfMeanResidual>10,],aes(label=mutationType),size=4)+
  ggtitle("the sequence+shape and shape only models fit equivalently")
  

combo_predictionsPlot_labelOutliers
######## seems like it's related to cpg presence in the 7mer but not central? ########
# label any that contains a cpg in the ancestral motif (doesn't have to be central)
combo$ancestralMotif <- substr(combo$mutationType,1,7)
head(combo)
combo$CpGAnywhereinAncestral7Mer <- "No"
combo[grep("CG",combo$ancestralMotif),]$CpGAnywhereinAncestral7Mer <- "Yes"

#### just one model and pop for now: 
combo_predictionsPlot_labelCpGAnywhere <-  ggplot(combo[combo$model=="model_009_SeqPlusShapes" & combo$population=="Mmd",], aes(y=.pred,x=outcome,color=CpGAnywhereinAncestral7Mer,shape=CpGAnywhereinAncestral7Mer))+
  geom_point(alpha=0.6,size=.8)+
  geom_abline()+
  facet_wrap(~mutationLabel,scales="free")+
  theme_bw()
combo_predictionsPlot_labelCpGAnywhere
ggsave(paste0(plotdir,"model_009_Mmd_only_IMPORTANTPLOT_SHOWINGLACKOFFITTO7mersContainingCpg.png"),width=10,height=8)


########### figure out if CpG containing central 3mer mutation types also have additional CpGs elsewhere ################
# use str_count() to get CpG count per motif
combo$ancestral_CpG_Count <- str_count(combo$ancestralMotif,"CG")
# want not abs value residual
combo$residual <- combo$.pred - combo$outcome
residualPlot <- ggplot(combo[combo$model=="model_009_SeqPlusShapes",],aes(x=as.factor(ancestral_CpG_Count),y=residual,group=ancestral_CpG_Count))+
  geom_boxplot()+
  facet_wrap(~mutationLabel~model,scales="free")+
  theme_bw()
residualPlot
ggsave(paste0(plotdir,"model_009_ResidualsVsCpGCount.png"),residualPlot,height=12,width=14)

residualPlot2 <- ggplot(combo[combo$model=="model_009_SeqPlusShapes",],aes(fill=as.factor(ancestral_CpG_Count),x=residual,group=ancestral_CpG_Count))+
  geom_histogram()+
  facet_wrap(~mutationLabel,scales="free")+
  theme_bw()
residualPlot2
ggsave(paste0(plotdir,"model_009_ResidualsVsCpGCount.histrogram.png"),residualPlot2,height=12,width=14)

View(combo)

### plot based on number of CpGs
#### just one model and pop for now: 
combo_predictionsPlot_labelCpGCount <-  ggplot(combo[combo$model=="model_009_SeqPlusShapes" & combo$population=="Mmd",], aes(y=.pred,x=outcome,color=as.factor(ancestral_CpG_Count)))+
  geom_point(size=1.2,alpha=1)+
  geom_abline()+
  facet_wrap(~mutationLabel,scales="free")+
  theme_bw()+
  scale_color_viridis_d(option = "inferno")
combo_predictionsPlot_labelCpGCount
ggsave(paste0(plotdir,"model_009_Mmd_only_IMPORTANTPLOT_SHOWINGLACKOFFITTO7mersContainingCpg.CpGCount.png"),combo_predictionsPlot_labelCpGCount,width=14,height=12)

####### try comparing species #########

combo_spread <- pivot_wider(combo[,c("mutationType","population","model","outcome",".pred","centralMutationType","CpGAnywhereinAncestral7Mer","mutationLabel","ancestral_CpG_Count")],id_cols = c(mutationType,population,centralMutationType,model,CpGAnywhereinAncestral7Mer,mutationLabel,ancestral_CpG_Count),names_from=population,values_from=c(outcome,.pred)) 

head(combo_spread)
# this is really useful: plot the difference between the two species and how well the predictions do (species are off of y=x line because of diff muttion rates; model picks that up! super cool)

# restrict to just one model
speciesComparisonPlot <- ggplot(combo_spread[combo_spread$model=="model_009_SeqPlusShapes",],aes(x=outcome_Mmd,y=outcome_Ms,color=as.factor(ancestral_CpG_Count)))+
  geom_point()+
  geom_point(aes(x=.pred_Mmd,y=.pred_Ms),shape=1,color="darkgray")+
  geom_abline()+
  facet_wrap(~mutationLabel,scales="free")+
  scale_color_viridis_d(option = "inferno")+
  theme_bw()
speciesComparisonPlot
ggsave(paste0(plotdir,"model_009_IMPORTANTPLOT_SHOWINGLACKOFFITTO7mersContainingCpgWhenComparingSpecies.png"),speciesComparisonPlot,height=12,width=14)

