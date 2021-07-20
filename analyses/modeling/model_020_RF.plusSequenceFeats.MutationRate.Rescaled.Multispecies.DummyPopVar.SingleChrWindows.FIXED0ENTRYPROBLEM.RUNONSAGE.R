############## try a new window (summed up chrs) -based model to practice ###########
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
sink(paste0(outdir,"logfile.sink.txt"),type="output") # this will only sink output not errors (errors will still go into errors dir)
#outdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210719_model_018_model_018_RF.plusSequenceFeats.MutationRate.Rescaled.Multispecies.DummyPopVar.SingleChrWindows_ranOnLaptop/"
#dir.create(outdir)
#todaysdate=format(Sys.Date(), "%Y%m%d")
#outcomeLabel="MutationRate"
#modelLabel="RF"
#description=paste0(modelLabel,".",outcomeLabel,".Multispecies.DummyPopVar") # short description of these experiments
#description
#chromCount=19
# per pop spectrum files:
#populations=c("Mmd","Ms")
# modeling outdir:
#outdir=paste0("/net/harris/vol1/home/beichman/DNAShape/analyses/modeling/experiments/",todaysdate,"_",description,"/") #  date specific 
#dir.create(outdir,recursive = T,showWarnings = F)

############## read in shapes ##############
print('reading in shapes')
shapedir="/net/harris/vol1/home/beichman/DNAShape/shapeDataForModeling/"
shapes <- read.table(paste0(shapedir,"firstOrder_featureTypes_allPossible7mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.txt"),header=T,sep="\t")
rownames(shapes) <- shapes$motif
# for testing on home computer: 
#shapes <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/firstOrder_featureTypes_allPossible7mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.txt",header=T,sep="\t")

print('reading in spectrum')
spectrumdir="/net/harris/vol1/home/beichman/DNAShape/spectrumDataForModeling/mouse/"

###### 20210719: fixed 0-entries that were incorrect before. the issue was:
# I had tried to put missing 7mers and missing mutation types back in, but I was not filling in target counts correctly for 0-count mutations; were ending up as NA and then being converted to 0 so seeming as missing targets, instead of targets being present but the mutation count = 0
# The problem was that some were actually 0 and were working, but then if mutations based on the ancestral 7mer were missing *entirely* from the chromosome my code  was treating that like the ancestral7mer was missing! When in fact it was present but just had no mutations of any type. So if an ancestral 7mer had at least one mutation of any type, then the 0-entries were filled in correctly. But if hadn't mutated at all on a chromosome then it was treated like the target count was 0 (incorrect!!!) --> 0/0 = NA and removed from the data. So I have fixed it that so that those chromosomes
allData_multipop <- read.table(paste0(spectrumdir,"TEMPORARY.NEWGROUPS.MULTIPOPULATION_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.INCLUDESCorrected0Entries.txt"),header=T) # this is the 7mer spectrum but down below you are going to sum up over 5mers  # 20210719 this file has corrected 0 entries
# note this is NOT summed up over chrs
# for testing on home computer:
#allData_multipop <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_spectrum_target_merged_PerChr/TEMPORARY.NEWGROUPS.MULTIPOPULATION_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.INCLUDESCorrected0Entries.txt",header=T) # 20210719 this file has corrected 0 entries 

dim(allData_multipop)
# need to get rid of NA entries
print("getting rid of NA mutation rates due to no ancestral targets observed") ## there wern't actually as many of these as I thought there were -- it was the issue with the 0 entries described above!!! 
allData_multipop <- na.omit(allData_multipop) 
dim(allData_multipop) # doesn't change any more
# and now want to split train is sp A + B spectra across odd chroms, train is sp A + B spectra across even chroms. want to have species membership as a feature! see if it's VIP or not

########## merge with shapes ##########

allData_multipop$derived7mer <- substr(allData_multipop$mutationType,9,15)
allData_multipop_intermediate <-merge(allData_multipop,shapes,by.x="ancestral7mer",by.y="motif")

allData_withShapes_unprocessed <- merge(allData_multipop_intermediate,shapes,by.x="derived7mer",by.y="motif",suffixes=c(".ancestral",".derived"))


######## add seq info (not yet 1-hot encoded, will encode with tidy models )##########
# just want ancestral 7mer 
allData_withShapes_unprocessed <- allData_withShapes_unprocessed %>% separate(ancestral7mer,into=c("Pos_1","Pos_2","Pos_3","Pos_4.ancestral","Pos_5","Pos_6","Pos_7"),remove=F,sep=c(1,2,3,4,5,6))
# also get derived bp:
allData_withShapes_unprocessed$Pos_4.derived <- substr(allData_withShapes_unprocessed$derived7mer,4,4)

########### want to rescale outcome variable rate so it's relative #######
# not adding epsilon because 
# not log scaling # call whatever you want your final outcome to be 'outcome' so that it's the same in all plots
# picked this because it's smaller than 1/genome size 
allData_withShapes_unprocessed <- allData_withShapes_unprocessed %>%
  group_by(population,newGroup,label) %>%
  mutate(outcome=mutationCount_divByTargetCount/sum(mutationCount_divByTargetCount)) 

# this maintains the ranking but rescales the outcome ; try it! 

####### make your splits into train/test ##########
# not splitting by population so there will be 2x as many train and test observations, one from each species
indices <-
  list(analysis   = which(allData_withShapes_unprocessed$label=="TRAIN"), 
       assessment = which(allData_withShapes_unprocessed$label=="TEST"))

split <- make_splits(indices,allData_withShapes_unprocessed)
saveRDS(split, file = paste0(outdir,"split.rds"))

# if you want to see what's what:
head(training(split),4)
head(testing(split),4)

train_data_cv <- group_vfold_cv(training(split),group=newGroup) # okay so I can make windows based on a grouping variable like so.

# to see a fold
# one fold: train_data_cv[[1]][[1]]
unique(analysis(train_data_cv[[1]][[1]])$newGroup)# the inside newGroups (chromosomes in this case): 
unique(assessment(train_data_cv[[1]][[1]])$newGroup) # the held out newGroup: 
# okay this works great! 

########### TRAIN/TEST MODEL JUST USING ONE FOLD SET #########  
oneFoldSetToTrainAndAssessOn <- train_data_cv[[1]][[1]]
saveRDS(oneFoldSetToTrainAndAssessOn, file = paste0(outdir,"oneFoldSetToTrainAndAssessOn.rds"))
# load back in:
#oneFoldSetToTrainAndAssessOn <- readRDS(paste0(outdir,"oneFoldSetToTrainAndAssessOn.rds"))
# rand_forest_Fold01_fit <- rand_forest_workflow %>%
#   last_fit(oneFoldSetToTrainAndAssessOn) # is very slow (a couple hours because it's based on all 7mer types)

########## RECIPE #########
####### should I sum up each non held-out fold somehow? skip for now. ##########
# Note about recipes: You initialize the recipe with some data so that it gets the column names, but then when I do fit, I specify which data to use for fitting the model and it's not data that is in the recipe (so the recipe was initialized with the full trainign split, but then model is fit with the fold I am using, I got worried that somehow there could be some data leakage with that, but I don't think there is because the model used to initialize the recipe isn't used for fitting if I don't say to. But just to be safe in future models maybe initialize the recipe with the fold you're using # 
rand_forest_processing_recipe <- 
  # which consists of the formula (outcome ~ predictors) (don't want to include the 'variable' column)
  
  recipe(outcome ~ .,data=training(oneFoldSetToTrainAndAssessOn)) %>% # 
  update_role(mutationType, new_role="7mer mutation type label") %>%
  step_rm(derived7mer,ancestral7mer, mutationCount,mutationCount_divByTargetCount,newGroup,label,ancestral7merCount) %>%
  step_dummy(all_nominal_predictors()) # KEEPING population in here as a predictor careful here that nothing else slips in! but dummy encoding it;; which RF doesn't need but xgboost and SHAP values does so just doing it ; this also dummy encodes the sequence with A as 0 0 0 , 1 0 0 = c etc. 
#Encoding each ancestral bp as a feature Pos_1, Pos_2, Pos_3 etc.
#Pos_4.ancestral and Pos_4.derived
#Then convert to one-hot binary variables where A is 0 0 0 , C is 1 0 0 , G is 0 1 0, T is 0 0 1 (ala L&S)
#But unlike L&S I’m also keeping in the derived position of the central bp

rand_forest_processing_recipe %>% summary()
rand_forest_processing_recipe
# can extract data like this try it :
#prepped <- prep(rand_forest_processing_recipe,training(split)) %>%
#  juice()

######### MODEL SPECIFICATION #########
rand_forest_ranger_model_specs <-
  rand_forest(trees = 1000, mtry = 32, min_n = 5) %>% # I added in tree number = 1000
  set_engine('ranger',importance="permutation",respect.unordered.factors="order",verbose=TRUE,num.threads=5) %>%
  set_mode('regression')
rand_forest_ranger_model_specs

############ WORKFLOW #############
######## make a workflow with recipe and model specs ########
rand_forest_workflow <- workflow() %>%
  # add the recipe
  add_recipe(rand_forest_processing_recipe) %>%
  # add the model
  add_model(rand_forest_ranger_model_specs)
rand_forest_workflow


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
# rand_forest_Fold01_fit_predictions$newGroup <- assessment(oneFoldSetToTrainAndAssessOn)$newGroup
# 
# rand_forest_Fold01_fit_predictions_plot <-  ggplot(rand_forest_Fold01_fit_predictions, aes(y=.pred,x=fractionOfSegSites_pernewGroup,color=centralMutationType,shape=population))+
#   geom_point()+
#   geom_abline()+
#   #scale_x_log10()+
#   #scale_y_log10()+
#   facet_wrap(~population)+
#   ggtitle("OUTCOME: frac of seg sites; ALL 7mer mutation types; random forest, trained on Fold01 ( all odd chrs but one, tested on just one)")
# rand_forest_Fold01_fit_predictions_plot


######## hmm need to use fit() not last_fit() that's irritating ##########
rand_forest_Fold01_fit_notlastfit <- fit(rand_forest_workflow,data = analysis(oneFoldSetToTrainAndAssessOn))
rand_forest_Fold01_fit_notlastfit

saveRDS(rand_forest_Fold01_fit_notlastfit, file = paste0(outdir,"modelTrainedOnOneFold.rds"))
# want to load it back in
#rand_forest_Fold01_fit_notlastfit <- readRDS(paste0(outdir,"modelTrainedOnOneFold.rds"))
# OOB prediction error (MSE):       3.130142e-06
# R squared (OOB):                  0.9729254 
# growing trees takes 1/2 hour; permutation importance takes another 1/2 hour.
# then predict based on held out assessment set of fold split:
rand_forest_Fold01_predictions <- predict(object =rand_forest_Fold01_fit_notlastfit, new_data=assessment(oneFoldSetToTrainAndAssessOn))
rand_forest_Fold01_predictions

truth_prediction_df <- cbind(assessment(oneFoldSetToTrainAndAssessOn),rand_forest_Fold01_predictions)
#View(truth_prediction_df) # VIEW doesn't show .pred column for some reason

windowOfAssessment=toString(unique(truth_prediction_df$newGroup))
windowOfAssessment
#windowOfAssessment=9
saveRDS(truth_prediction_df,file=paste0(outdir,"modelTrainedOnOneFold.PREDICTIONS.onWindow",windowOfAssessment,".rds")) # ah that's a problem for loading it back in -- need the window ID 

# load it back in: 
#truth_prediction_df <- readRDS(paste0(outdir,"modelTrainedOnOneFold.PREDICTIONS.onChr",windowOfAssessment,".rds"))
# get rsq
rsq(truth_prediction_df,truth=outcome,estimate=.pred)
# 1 rsq     standard       0.980
rmse(truth_prediction_df,truth=outcome,estimate=.pred)

#  rmse    standard     0.00155 ; kind of big? 
truth_prediction_df$centralMutationType <- paste0(substr(truth_prediction_df$mutationType,4,4),".",substr(truth_prediction_df$mutationType,12,12))
# note this doesn't have the 1-coded pops unless you juice() it but rows are still in same order 
rand_forest_Fold01_fit_predictions_plot <-  ggplot(truth_prediction_df, aes(y=.pred,x=outcome,color=centralMutationType,shape=population))+
  geom_point()+
  geom_abline()+
  facet_wrap(~population)+
  #ggtitle(paste0(description,"\ntrained on Fold01\n(all odd windows but one, tested on just window ",toString(unique(truth_prediction_df$window)),")"))+
  theme_bw()
rand_forest_Fold01_fit_predictions_plot

ggsave(paste0(outdir,"modelTrainedOnOneFold.PredictionsPlot.AssessedOnWindow",toString(unique(truth_prediction_df$newGroup)),".png"),rand_forest_Fold01_fit_predictions_plot,height=6,width=9)


########### facet by central mutation type and get individual rsqs ##############
rsqsPerSpeciesAndMutationType <- truth_prediction_df %>%
  group_by(centralMutationType,population) %>%
  rsq(truth=outcome,estimate=.pred)
rsqsPerSpeciesAndMutationType
write.table(rsqsPerSpeciesAndMutationType,paste0(outdir,"modelTrainedOnOneFold.Rsq.PerMutatationType.AssessedOnWindow",toString(unique(truth_prediction_df$newGroup)),".txt"),quote = F,row.names=F,sep="\t")
rsqPerMutplot <- ggplot(rsqsPerSpeciesAndMutationType,aes(x=centralMutationType,y=.estimate,fill=population))+
  geom_col(position="dodge")+
  theme_bw()+
  ylab("r-squared")
rsqPerMutplot

ggsave(paste0(outdir,"modelTrainedOnOneFold.Rsq.PerMutationType.AssessedOnWindow",toString(unique(truth_prediction_df$newGroup)),".png"),rsqPerMutplot,height=3,width=5)

####  add rsq values to plot if possible ###
rand_forest_Fold01_fit_predictions_plot_faceted <-  ggplot(truth_prediction_df, aes(y=.pred,x=outcome,color=centralMutationType))+
  geom_point()+
  geom_abline()+
  geom_text(data=rsqsPerSpeciesAndMutationType,aes(x=3e-6,y=1e-4,label=round(.estimate,4)),color="black",size=8)+
  facet_wrap(~centralMutationType~population,scales="free")+
  scale_y_log10()+
  scale_x_log10()+
  #ggtitle(paste0(description,"\ntrained on Fold01\n(all odd Windows but one, tested on just Window",toString(unique(truth_prediction_df$window)),")"))+
  theme_bw()+
  theme(text=element_text(size=14))
rand_forest_Fold01_fit_predictions_plot_faceted

ggsave(paste0(outdir,"modelTrainedOnOneFold.PredictionsPlot.FacetedPerMutationType.AssessedOnWindow",toString(unique(truth_prediction_df$newGroup)),".png"),rand_forest_Fold01_fit_predictions_plot_faceted,height=12,width=15)

############ VIP: variable importance ###########
ranger_obj <- pull_workflow_fit(rand_forest_Fold01_fit_notlastfit)$fit
ranger_obj
#OOB prediction error (MSE):       3.130142e-06 
# R squared (OOB):                  0.9729254
vi_scores <- vip::vi(ranger_obj)
write.table(vi_scores,paste0(outdir,"modelTrainedOnOneFold.VIPScores.PerMutatationType.AssessedOnWindow",toString(unique(truth_prediction_df$newGroup)),".txt"),quote = F,row.names=F,sep="\t")

vip_plot <- vip(vi_scores,include_type = T,num_features = 100)
vip_plot

ggsave(paste0(outdir,"modelTrainedOnOneFold.VIP.Plot.AssessedOnWindow",toString(unique(truth_prediction_df$newGroup)),".png"),vip_plot,height=12,width=5)
#vip(ranger_obj,include_type = T,num_features = 20)

########## try to find interactions? ############
#vint(ranger_obj,feature_names=c("feature_1_Shear_2.derived","feature_1_Shear_2.ancestral"),progress=T),train=analysis(oneFoldSetToTrainAndAssessOn)) # CRASHED

# can't get it to work with categorical variables 

######## want to try to find sites that are at very different rates between the mouse species and focus on those and what makes them tick #############
head(truth_prediction_df)
truth_prediction_df_spread <- pivot_wider(truth_prediction_df[,c("mutationType","population","outcome",".pred","centralMutationType")],id_cols = c(mutationType,population,centralMutationType),names_from=population,values_from=c(outcome,.pred)) 

head(truth_prediction_df_spread)
# this is really useful: plot the difference between the two species and how well the predictions do (species are off of y=x line because of diff muttion rates; model picks that up! super cool)
speciesComparisonPlot <- ggplot(truth_prediction_df_spread,aes(x=outcome_Mmd,y=outcome_Ms))+
  geom_point()+
  geom_point(aes(x=.pred_Mmd,y=.pred_Ms),shape=1,color="red")+
  geom_abline()+
  facet_wrap(~centralMutationType,scales="free")
speciesComparisonPlot
ggsave(paste0(outdir,"modelTrainedOnOneFold.SpeciesXYComparison.AssessedOnWindow",toString(unique(truth_prediction_df$newGroup)),".png"),speciesComparisonPlot,height=5,width=7)

# plot both species together
plotBothSpeciesTogether <- ggplot(truth_prediction_df,aes(y=.pred,x=outcome,color=population))+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  geom_point(size=0.1)+
  #scale_x_log10()+
  #scale_y_log10()+
  theme_bw()+
  geom_abline()
plotBothSpeciesTogether
ggsave(paste0(outdir,"modelTrainedOnOneFold.ObsExpected.SpeciesTogether.AssessedOnWindow",toString(unique(truth_prediction_df$newGroup)),".png"),plotBothSpeciesTogether,height=5,width=7)

sink()
