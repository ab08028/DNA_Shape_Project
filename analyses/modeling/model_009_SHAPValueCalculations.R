require(tidymodels) # installing on hoffman 
require(tidyverse) # instead of caret going to use tidymodels
require(workflows)
require(tune)
require(vip) # for importance and shap values
require(ranger) # for random forest model
#require(glmnet) # for lasso and linear modeling (maybe)
require(ggplot2)
require(ggrepel)
require(ggbeeswarm)
require(reshape2)
require(devtools)
library(doParallel)
registerDoParallel()

######### try to get SHAP values to work ########
#dir on home laptop:
#outdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210625_model_009_RF.plusSequenceFeats.MutationRate.Rescaled.Multispecies.DummyPopVar.0sIncluded.MultiChromosomeWindows/"
# outdir on sage:
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
sink(paste0(outdir,"SHAP.logfile.sink.txt"),type="output") # this will only sink output not errors (errors will still go into errors dir)

#outdir="/net/harris/vol1/home/beichman/DNAShape/analyses/modeling/experiments/20210625_model_009_RF.plusSequenceFeats.MutationRate.Rescaled.Multispecies.DummyPopVar.0sIncluded.MultiChromosomeWindows/"
split <- readRDS(paste0(outdir,"split.rds"))
########## RECIPE #########
####### should I sum up each non held-out fold somehow? skip for now. ##########
rand_forest_processing_recipe <- 
  # which consists of the formula (outcome ~ predictors) (don't want to include the 'variable' column)
  
  recipe(outcome ~ .,data=training(split)) %>% # 
  update_role(mutationType, new_role="7mer mutation type label") %>%
  step_rm(derived7mer,ancestral7mer, mutationCount,mutationCount_divByTargetCount,newGroup,label,ancestral7merCount) %>%
  step_dummy(all_nominal_predictors())
rand_forest_processing_recipe %>% summary()
rand_forest_processing_recipe
# can extract data like this try it :
#prepped <- prep(rand_forest_processing_recipe,training(split)) %>%
#  juice()

######### Need ranger object ####### 

model <- readRDS(paste0(outdir,"modelTrainedOnOneFold.rds"))
# need ranger obj
ranger_obj <- pull_workflow_fit(model)$fit


######### Need to juice the training and assessment data and exclude anything that isn't features ###### 

# juice and get rid of outcome variables and anything that isn't predictions (like mutaiton type)
Xtrain <- prep(rand_forest_processing_recipe, training(split)) %>% 
  juice() %>% 
  select(-c(mutationType,outcome)) %>% 
  as.matrix()
head(Xtrain)


Xtest <- prep(rand_forest_processing_recipe, testing(split)) %>% 
  juice() %>% 
  select(-c(mutationType,outcome)) %>% 
  as.matrix()
head(Xtest)

# function describing how to get predictions from the model 
pfun <- function(object, newdata) {
  predict(object, data = newdata)$predictions
}


# a bit slow: 

shap <- fastshap::explain(ranger_obj, X = Xtrain,pred_wrapper=pfun,newdata=Xtest) ##
saveRDS(shap,file=paste0(outdir,"SHAPResults.rds"))
#read back in if you need it
#shap <- readRDS(paste0(outdir,"SHAPResults.rds"))
#eventually want: nsim=5,newdata=Xtest,adjust=T)
# oooh is it working? try to increase nsim (number of MC simulations ); adjust ot satisfy additivity principle
# VERY SLOW
#shap

#hist(shap$feature_1_Slide_4.ancestral)
#hist(shap$feature_1_Shift_4.derived)

# Aggregate Shapley values
shap_imp <- data.frame(
  Variable = names(shap),
  Importance = apply(shap, MARGIN = 2, FUN = function(x) sum(abs(x)))
)
shap_imp

shapplot1 <- ggplot(shap_imp, aes(reorder(Variable, Importance), Importance)) +
  geom_col() +
  coord_flip() +
  xlab("") +
  ylab("mean(|Shapley value|)")
shapplot1
ggsave(paste0(outdir,"random_forest.shapImportance.png"),shapplot1,height=10,width=6)
######## shapley dependence plot ########
# need dataframe of training data:
# note X above was a matrix but this has to be a df
Xdftrain <- prep(rand_forest_processing_recipe, training(split)) %>% 
  juice()  #%>%

head(Xdftrain)

# note X above was a matrix but this has to be a df
Xdftest <- prep(rand_forest_processing_recipe, testing(split)) %>% 
  juice()  #%>%

head(Xdftest)



featuresToPlot=names(shap)
dir.create(paste0(outdir,"perFeatureDependencePlots/"),showWarnings = F)
for(feature in featuresToPlot){
  
  shapdependenceplot1 <- autoplot(shap, 
                                  type = "dependence", 
                                  feature = feature, 
                                  X = Xdftest, # this hsould be test data if you ran shap with newdata=testing data! sholud be train if you just ran it without any new data
                                  smooth = TRUE, color_by = "outcome")+
    scale_color_viridis_c(trans="log10")+
    ggtitle(feature)+
    theme_bw()+
    labs(color="outcome")
  shapdependenceplot1
  ggsave(paste0(outdir,"/perFeatureDependencePlots/",mutationType,".random_forest.shapDependence.",feature,".png"),shapdependenceplot1,height=4,width=6)
}
# can use these to try and detect interactions by coloring by other variables -- is there a systematic way to detect though?
# color_by = "fractionOfAllSegSites"
# for a single motif of interest
motifOfInterest="TTTAAAA.TTTTAAA"
rowNumberForForce=which(Xdftest$variable==motifOfInterest)
rowNumberForForce

shapplot2 <- autoplot(shap, type = "contribution", row_num = rowNumberForForce) +
  ggtitle(motifOfInterest)
shapplot2
ggsave(paste0(outdir,mutationType,".random_forest.shapImportance.",motifOfInterest,".png"),shapplot2,height=10,width=6)



# requires python and shap to be correctly pathed: let's wait on this because don't want to deal with python on sage ? 
#forceplot1 <- force_plot(object = shap[rowNumberForForce,], 
#                         feature_values = select(Xdftest[rowNumberForForce,],-c#(mutationType,outcome)), 
#                         display = "html")
#write.table(forceplot1,paste0(outdir,mutationType,".random_forest.shapForcePlot.",motifOfInterest,".html"),row.names=F,quote=F,col.names=F)

# multiple obs at once:
#forceplot2 <- force_plot(object = shap, 
#                         feature_values = select(Xdftest,-c(mutationType,outcome)), 
#                         display = "html") # , 
#link = "logit") 
#write.table(forceplot2,paste0(outdir,mutationType,".random_forest.shapForcePlot.AllSamples.html"),row.names=F,quote=F,col.names=F)
# this saves an html plot file! 


####### SHAP summary plot (really useful) ######
shap_labeled <- shap
shap_labeled$mutationtype <- Xdftest$variable
head(shap_labeled)
shap_melt <- melt(shap_labeled)
head(shap_melt)
colnames(shap_melt) <- c("mutationtype","feature","SHAP.value")
# combine with xdf:
Xdftrain_melt <- melt(Xdftrain,varnames=c(mutationtype)) # are these based on train or test?
colnames(Xdftrain_melt) <- c("mutationtype","feature","feature.value")
Xdftest_melt <- melt(Xdftest) # are these based on train or test?
colnames(Xdftest_melt) <- c("mutationtype","feature","feature.value")


shap_plusFeatures <- merge(shap_melt,Xdftest_melt,by=c("mutationtype","feature"))
###### do min max normalization: 
shap_plusFeatures <- shap_plusFeatures %>%
  group_by(feature) %>%
  mutate(minmaxnormalized.feature.value=(feature.value-min(feature.value))/(max(feature.value)-min(feature.value)))
head(shap_plusFeatures)
# try to merge shap and Xdf 
# want to order by importance (mean abs shap )
featuresInOrderOfImportance <- shap_imp %>%
  arrange(desc(Importance))
# order features in order of importance:
shap_plusFeatures$feature <- factor(shap_plusFeatures$feature,levels=rev(featuresInOrderOfImportance$Variable))
shapSummaryPlot <- ggplot(shap_plusFeatures,aes(y=feature,x=SHAP.value,color=minmaxnormalized.feature.value))+
  #geom_point()+
  geom_quasirandom(groupOnX = F) + # bee swarm so ti's like a violin plot
  theme_bw()+
  ggtitle(paste0(mutationType ,"SHAP summary plot"))+
  scale_color_viridis_c(labels=c("low","high"),breaks=c(0,1),limits=c(0,1),option = "plasma")
shapSummaryPlot
ggsave(paste0(outdir,mutationType,".random_forest.shapSummaryPlot.USEFUL.png"),shapSummaryPlot,height=12,width=10)

# just do top 10 most important:
top10 <- head(featuresInOrderOfImportance$Variable,10)
shapSummaryPlot_top10 <- ggplot(shap_plusFeatures[shap_plusFeatures$feature %in% top10,],aes(y=feature,x=SHAP.value,color=minmaxnormalized.feature.value))+
  #geom_point()+
  geom_quasirandom(groupOnX = F) + # bee swarm so ti's like a violin plot; want to cluster on y though
  
  theme_bw()+
  ggtitle(paste0(mutationType ,"SHAP summary plot"))+
  scale_color_viridis_c(labels=c("low","high"),breaks=c(0,1),limits=c(0,1),option = "plasma")
shapSummaryPlot_top10
ggsave(paste0(outdir,mutationType,".random_forest.shapSummaryPlot.top10.USEFUL.png"),shapSummaryPlot_top10,height=4,width=10)

### just do all the population ones: 
shapSummaryPlot_populationOnly <- ggplot(shap_plusFeatures[grepl("population",shap_plusFeatures$feature),],aes(y=feature,x=SHAP.value,color=minmaxnormalized.feature.value))+
  #geom_point()+
  geom_quasirandom(groupOnX = F) + # bee swarm so ti's like a violin plot; want to cluster on y though
  
  theme_bw()+
  ggtitle(paste0(mutationType ,"SHAP summary plot"))+
  scale_color_viridis_c(labels=c("low","high"),breaks=c(0,1),limits=c(0,1),option = "plasma")
shapSummaryPlot_top10
ggsave(paste0(outdir,mutationType,".random_forest.shapSummaryPlot.top10.USEFUL.png"),shapSummaryPlot_top10,height=4,width=10)

