#!/usr/bin/env Rscript
#$ -m bea
#$ -M annabel.beichman@gmail.com

############## try a new chromosome -based model to practice ###########
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
#devtools::install_github("MI2DataLab/randomForestExplainer")
#require(randomForestExplainer) # too slow
#devtools::install_github('ModelOriented/treeshap')
#require(treeshap)
#require(shapper)

#install_shap()
tidymodels_prefer() # use this 'tidymodels_prefer()' uses the 'conflicted' package to handle common conflicts with tidymodels and other packages. <-- should add this to my script
set.seed(42) # so results are reproducible 
# doing this in the wrapper script instead:


args <- commandArgs(trailingOnly = TRUE)
description <- args[1]
outdir <- args[2]
# trying to use sink() to catch all output: (errors will go to a different file)
#sink(paste0(outdir,"logfile.sink.txt"),type="output") # this will only sink output not errors (errors will still go into errors dir)

#todaysdate=format(Sys.Date(), "%Y%m%d")
#outcomeLabel="MutationRate"
#modelLabel="RF"
#description=paste0(modelLabel,".",outcomeLabel,".Multispecies.DummyPopVar") # short description of these experiments
#description
chromCount=19
# per pop spectrum files:
populations=c("Mmd","Ms")
# modeling outdir:
#outdir=paste0("/net/harris/vol1/home/beichman/DNAShape/analyses/modeling/experiments/",todaysdate,"_",description,"/") #  date specific 
#dir.create(outdir,recursive = T,showWarnings = F)

############## read in shapes ##############
print('reading in shapes')
shapedir="/net/harris/vol1/home/beichman/DNAShape/shapeDataForModeling/"
shapes <- read.table(paste0(shapedir,"firstOrder_featureTypes_allPossible7mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.txt"),header=T,sep="\t")
rownames(shapes) <- shapes$motif


print('reading in spectrum')
spectrumdir="/net/harris/vol1/home/beichman/DNAShape/spectrumDataForModeling/mouse/"
allData_multipop <- read.table(paste0(spectrumdir,"MULTIPOPULATION_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.txt"),header=T) # 

# and now want to split train is sp A + B spectra across odd chroms, train is sp A + B spectra across even chroms. want to have species membership as a feature! see if it's VIP or not

# merge with shapes

allData_multipop$derived7mer <- substr(allData_multipop$mutationType,9,15)
allData_multipop_intermediate <-merge(allData_multipop,shapes,by.x="ancestral7mer",by.y="motif")

allData_withShapes_unprocessed <- merge(allData_multipop_intermediate,shapes,by.x="derived7mer",by.y="motif",suffixes=c(".ancestral",".derived"))




####### make your splits into train/test ##########
# not splitting by population so there will be 2x as many train and test observations, one from each species
indices <-
  list(analysis   = which(allData_withShapes_unprocessed$label=="TRAIN"), 
       assessment = which(allData_withShapes_unprocessed$label=="TEST"))

split <- make_splits(indices,allData_withShapes_unprocessed)

# if you want to see what's what:
#head(training(split),4)
#head(testing(split),4)

train_data_cv <- group_vfold_cv(training(split),group=window) # okay so I can make windows based on a grouping variable like so.

# to see a fold
# one fold: train_data_cv[[1]][[1]]
#unique(analysis(train_data_cv[[1]][[1]])$window)# the inside windows (chromosomes in this case): 
#unique(assessment(train_data_cv[[1]][[1]])$window) # the held out window: 
# okay this works great! 

########## RECIPE #########
####### should I sum up each non held-out fold somehow? skip for now. ##########
print('setting up recipe')
rand_forest_processing_recipe_OutcomeFracSegSites_withPopAsPredictor <- 
  # which consists of the formula (outcome ~ predictors) (don't want to include the 'variable' column)
  recipe(mutationCount_divByTargetCount ~ .,data=training(split)) %>% # 
  update_role(mutationType, new_role="7mer mutation type label") %>%
  step_rm(derived7mer,ancestral7mer, mutationCount,window,label,ancestral7merCount) %>%
  step_dummy(all_nominal_predictors()) # KEEPING population in here as a predictor careful here that nothing else slips in! but dummy encoding it;; which RF doesn't need but xgboost and SHAP values does so just doing it 
rand_forest_processing_recipe_OutcomeFracSegSites_withPopAsPredictor %>% summary()
rand_forest_processing_recipe_OutcomeFracSegSites_withPopAsPredictor


######### MODEL SPECIFICATION #########
print('setting up model specs')
rand_forest_ranger_model_specs <-
  rand_forest(trees = 1000, mtry = 32, min_n = 5) %>% # I added in tree number = 1000
  set_engine('ranger',importance="permutation",respect.unordered.factors="order",verbose=TRUE,num.threads=10) %>%
  set_mode('regression')
rand_forest_ranger_model_specs

############ WORKFLOW #############
######## make a workflow with recipe and model specs ########
print('setting up model workflow')
rand_forest_workflow <- workflow() %>%
  # add the recipe
  add_recipe(rand_forest_processing_recipe_OutcomeFracSegSites_withPopAsPredictor) %>%
  # add the model
  add_model(rand_forest_ranger_model_specs)
rand_forest_workflow

########### TRAIN/TEST MODEL JUST USING ONE FOLD SET #########  
print('pulling out one fold to train/assess with for now')
oneFoldSetToTrainAndAssessOn <- train_data_cv[[1]][[1]]
saveRDS(oneFoldSetToTrainAndAssessOn, file = paste0(outdir,"oneFoldSetToTrainAndAssessOn.rds"))
# load back in:
#oneFoldSetToTrainAndAssessOn <- readRDS(paste0(outdir,"oneFoldSetToTrainAndAssessOn.rds"))
# rand_forest_Fold01_fit <- rand_forest_workflow %>%
#   last_fit(oneFoldSetToTrainAndAssessOn) # is very slow (a couple hours because it's based on all 7mer types)
# 
# ########### assess fit to held-out fold #########
# rand_forest_Fold01_fit %>% collect_metrics()
# 
# rand_forest_Fold01_fit_predictions <- rand_forest_Fold01_fit %>% 
#   collect_predictions()
# 
# head(rand_forest_Fold01_fit_predictions)
# 
# rand_forest_Fold01_fit_predictions$MutationType <- assessment(oneFoldSetToTrainAndAssessOn)$mutationType
# rand_forest_Fold01_fit_predictions$population <- assessment(oneFoldSetToTrainAndAssessOn)$population
# rand_forest_Fold01_fit_predictions$centralMutationType <- paste0(substr(rand_forest_Fold01_fit_predictions$MutationType,4,4),".",substr(rand_forest_Fold01_fit_predictions$MutationType,12,12))
# rand_forest_Fold01_fit_predictions$window <- assessment(oneFoldSetToTrainAndAssessOn)$window
# 
# rand_forest_Fold01_fit_predictions_plot <-  ggplot(rand_forest_Fold01_fit_predictions, aes(y=.pred,x=fractionOfSegSites_perWindow,color=centralMutationType,shape=population))+
#   geom_point()+
#   geom_abline()+
#   #scale_x_log10()+
#   #scale_y_log10()+
#   facet_wrap(~population)+
#   ggtitle("OUTCOME: frac of seg sites; ALL 7mer mutation types; random forest, trained on Fold01 ( all odd chrs but one, tested on just one)")
# rand_forest_Fold01_fit_predictions_plot


######## hmm need to use fit() not last_fit() that's irritating ##########
print('starting to fit model')
rand_forest_Fold01_fit_notlastfit <- fit(rand_forest_workflow,data = analysis(oneFoldSetToTrainAndAssessOn))
rand_forest_Fold01_fit_notlastfit

saveRDS(rand_forest_Fold01_fit_notlastfit, file = paste0(outdir,"modelTrainedOnOneFold.rds"))
# want to load it back in
#rand_forest_Fold01_fit_notlastfit <- readRDS(paste0(outdir,"modelTrainedOnOneFold.rds"))
# OOB prediction error (MSE):       3.130142e-06
# R squared (OOB):                  0.9729254 
# growing trees takes 1/2 hour; permutation importance takes another 1/2 hour.
# then predict based on held out assessment set of fold split:
print(paste(Sys.time(), 'starting predictions')
rand_forest_Fold01_predictions <- predict(object =rand_forest_Fold01_fit_notlastfit, new_data=assessment(oneFoldSetToTrainAndAssessOn))
# rand_forest_Fold01_predictions

truth_prediction_df <- cbind(assessment(oneFoldSetToTrainAndAssessOn),rand_forest_Fold01_predictions)
View(truth_prediction_df) # VIEW doesn't show .pred column for some reason

windowOfAssessment=toString(unique(truth_prediction_df$window))
saveRDS(truth_prediction_df,file=paste0(outdir,"modelTrainedOnOneFold.PREDICTIONS.onChr",windowOfAssessment,".rds")) # ah that's a problem for loading it back in -- need the window ID 

# load it back in: 
#truth_prediction_df <- readRDS(paste0(outdir,"modelTrainedOnOneFold.PREDICTIONS.onChr",windowOfAssessment,".rds"))
# get rsq
rsq(truth_prediction_df,truth=mutationCount_divByTargetCount,estimate=.pred)
# 1 rsq     standard       0.980
rmse(truth_prediction_df,truth=mutationCount_divByTargetCount,estimate=.pred)

#  rmse    standard     0.00155 ; kind of big? 
truth_prediction_df$centralMutationType <- paste0(substr(truth_prediction_df$mutationType,4,4),".",substr(truth_prediction_df$mutationType,12,12))
# note this doesn't have the 1-coded pops unless you juice() it but rows are still in same order 
rand_forest_Fold01_fit_predictions_plot <-  ggplot(truth_prediction_df, aes(y=.pred,x=mutationCount_divByTargetCount,color=centralMutationType,shape=population))+
  geom_point()+
  geom_abline()+
  scale_x_log10()+
  scale_y_log10()+
  facet_wrap(~population)+
  ggtitle(paste0("OUTCOME: ",outcomeLabel,"\nALL 7mer mutation types;\n",modelLabel," trained on Fold01\n(all odd chrs but one, tested on just chr",toString(unique(truth_prediction_df$window)),")"))+
  theme_bw()
rand_forest_Fold01_fit_predictions_plot

ggsave(paste0(outdir,"modelTrainedOnOneFold.PredictionsPlot.AssessedOnChr",windowOfAssessment,".png"),rand_forest_Fold01_fit_predictions_plot,height=6,width=9)


########### facet by central mutation type and get individual rsqs ##############
rsqsPerSpeciesAndMutationType <- truth_prediction_df %>%
  group_by(centralMutationType,population) %>%
  rsq(truth=mutationCount_divByTargetCount,estimate=.pred)
rsqsPerSpeciesAndMutationType
write.table(rsqsPerSpeciesAndMutationType,paste0(outdir,"modelTrainedOnOneFold.Rsq.PerMutatationType.AssessedOnChr",toString(unique(truth_prediction_df$window)),".txt"),quote = F,row.names=F,sep="\t")
rsqPerMutplot <- ggplot(rsqsPerSpeciesAndMutationType,aes(x=centralMutationType,y=.estimate,fill=population))+
  geom_col(position="dodge")+
  theme_bw()+
  ylab("r-squared")
rsqPerMutplot

ggsave(paste0(outdir,"modelTrainedOnOneFold.Rsq.PerMutationType.AssessedOnChr",windowOfAssessment,".png"),rsqPerMutplot,height=3,width=5)

####  add rsq values to plot if possible ###
rand_forest_Fold01_fit_predictions_plot_faceted <-  ggplot(truth_prediction_df, aes(y=.pred,x=mutationCount_divByTargetCount,color=centralMutationType))+
  geom_point()+
  geom_abline()+
  geom_text(data=rsqsPerSpeciesAndMutationType,aes(x=1e-4,y=0.1,label=round(.estimate,4)),color="black")+
  scale_x_log10()+
  scale_y_log10()+
  facet_grid(~centralMutationType~population)+
  ggtitle(paste0("OUTCOME: ",outcomeLabel,"\nALL 7mer mutation types;\n",modelLabel," trained on Fold01\n(all odd chrs but one, tested on just chr",windowOfAssessment,")"))+
  theme_bw()
rand_forest_Fold01_fit_predictions_plot_faceted

ggsave(paste0(outdir,"modelTrainedOnOneFold.PredictionsPlot.FacetedPerMutationType.AssessedOnChr",windowOfAssessment,".png"),rand_forest_Fold01_fit_predictions_plot_faceted,height=12,width=9)

############ VIP: variable importance ###########
ranger_obj <- pull_workflow_fit(rand_forest_Fold01_fit_notlastfit)$fit
ranger_obj
#OOB prediction error (MSE):       3.130142e-06 
# R squared (OOB):                  0.9729254
vi_scores <- vip::vi(ranger_obj)
write.table(vi_scores,paste0(outdir,"modelTrainedOnOneFold.VIPScores.PerMutatationType.AssessedOnChr",windowOfAssessment,".txt"),quote = F,row.names=F,sep="\t")

vip_plot <- vip(vi_scores,include_type = T,num_features = 100)
vip_plot

ggsave(paste0(outdir,"modelTrainedOnOneFold.VIP.Plot.AssessedOnChr",windowOfAssessment,".png"),vip_plot,height=12,width=5)
#vip(ranger_obj,include_type = T,num_features = 20)

########## try to find interactions? ############
#vint(ranger_obj,feature_names=c("feature_1_Shear_2.derived","feature_1_Shear_2.ancestral",progress=T),train=analysis(oneFoldSetToTrainAndAssessOn)) # CRASHED

# can't get it to work with categorical variables ^ focus on this! try vint again now that I have dummy vars??

######## want to try to find sites that are at very different rates between the mouse species and focus on those and what makes them tick #############
head(truth_prediction_df)
truth_prediction_df_spread <- pivot_wider(truth_prediction_df[,c("mutationType","population","mutationCount_divByTargetCount",".pred","centralMutationType")],id_cols = c(mutationType,population,centralMutationType),names_from=population,values_from=c(mutationCount_divByTargetCount,.pred)) 

head(truth_prediction_df_spread)
# this is really useful: plot the difference between the two species and how well the predictions do (species are off of y=x line because of diff muttion rates; model picks that up! super cool)
speciesComparisonPlot <- ggplot(truth_prediction_df_spread,aes(x=mutationCount_divByTargetCount_Mmd,y=mutationCount_divByTargetCount_Ms))+
  geom_point()+
  geom_point(aes(x=.pred_Mmd,y=.pred_Ms),shape=1,color="red")+
  geom_abline()
speciesComparisonPlot
ggsave(paste0(outdir,"modelTrainedOnOneFold.SpeciesXYComparison.AssessedOnChr",windowOfAssessment,".png"),speciesComparisonPlot,height=5,width=7)

# plot both species together
plotBothSpeciesTogether <- ggplot(truth_prediction_df,aes(y=.pred,x=mutationCount_divByTargetCount,color=population))+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  geom_point(size=0.1)+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline()
plotBothSpeciesTogether
ggsave(paste0(outdir,"modelTrainedOnOneFold.ObsExpected.SpeciesTogether.AssessedOnChr",windowOfAssessment,".png"),plotBothSpeciesTogether,height=5,width=7)


############# try random forest explainer -- TOO SLOW #######################
#require(randomForestExplainer)
# try dist of min depth (slow)
#min_depth_frame <- min_depth_distribution(ranger_obj)
#head(min_depth_frame, n = 10)


# measure importance  -- too slow! 
#importance_frame <- measure_importance(ranger_obj) # incredibly slow!!!! 
# started around 3pm on 6/14
#saveRDS(importance_frame,file=paste0(outdir,"modelTrainedOnOneFold.IMPORATANCEFRAME.RFExplainer.Assessed.onChr",windowOfAssessment,".rds"))



######## do some shap value stuff ##############
Xtrain <- prep(rand_forest_processing_recipe_OutcomeFracSegSites_withPopAsPredictor, analysis(oneFoldSetToTrainAndAssessOn)) %>% 
  juice() %>% 
  select(-c(mutationType)) %>% 
  as.matrix()
head(Xtrain)

analysisdf_processed <- prep(rand_forest_processing_recipe_OutcomeFracSegSites_withPopAsPredictor, analysis(oneFoldSetToTrainAndAssessOn)) %>% 
  juice() %>% 
  select(-c(mutationType)) 
head(analysisdf_processed)
# now has population_Ms
# problem: population needs to be encoded not as categorcial for shapley to work -- annoying 

# matrix version:
Xtest <- prep(rand_forest_processing_recipe_OutcomeFracSegSites_withPopAsPredictor, assessment(oneFoldSetToTrainAndAssessOn)) %>% 
  juice() %>% 
  select(-c(mutationType)) %>% 
  as.matrix()
head(Xtest)

# df version:
assessmentdf_processed <- prep(rand_forest_processing_recipe_OutcomeFracSegSites_withPopAsPredictor, assessment(oneFoldSetToTrainAndAssessOn)) %>% 
  juice() %>% 
  select(-c(mutationType)) 
head(assessmentdf_processed)


pfun <- function(object, newdata) {
  predict(object, data = newdata)$predictions
}


# a bit slow: 

shap <- fastshap::explain(ranger_obj, X = Xtrain,pred_wrapper=pfun,newdata=Xtest) ##
#eventually want: nsim=5,newdata=Xtest,adjust=T)
# oooh is it working? try to increase nsim (number of MC simulations ); adjust ot satisfy additivity principle
# VERY SLOW
shap




saveRDS(shap, file = paste0(outdir,"modelTrainedOnOneFold.SHAPValues.rds"))


hist(shap$feature_1_Slide_4.ancestral)
hist(shap$feature_1_Shift_4.derived)

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

ggsave(paste0(outdir,"modelTrainedOnOneFold.SHAPVALUES.AssessedOnChr",windowOfAssessment,".png"),shapplot1,height=12,width=6)


#sink() # end sink