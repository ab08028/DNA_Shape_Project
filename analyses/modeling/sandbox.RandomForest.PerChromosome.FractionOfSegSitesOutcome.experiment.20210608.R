
############## try a new chromosome -based model to practice ###########


require(tidymodels)
require(tidyverse) # instead of caret going to use tidymodels
require(workflows)
require(tune)
require(vip) # for importance and shap values
require(ranger) # for random forest model
require(glmnet) # for lasso and linear modeling (maybe)
require(ggplot2)
require(ggrepel)
require(ggbeeswarm)
require(reshape2)




# want counts divided by targets
targetdir="/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_target_files"
# per pop spectrum files:
population="Mmd"
spectrumdir=paste0("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_spectrum_files/",population,"/")
spectrumoutdir="/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_spectrum_target_merged_PerChr/"
dir.create(spectrumoutdir,showWarnings=F)
# merge per chromosome # 
# need to separate odd and even
chromCount=19
#oddChr=seq(1,chromCount,by=2)
#evenChr=seq(2,chromCount-1,by=2)
####### okay let's try this ###########

tidymodels_prefer() # use this 'tidymodels_prefer()' uses the 'conflicted' package to handle common conflicts with tidymodels and other packages. <-- should add this to my script
set.seed(42) # so results are reproducible 
outdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/sandbox/"

####### only once: then read it in below instead
allData <- data.frame()
for(chr in seq(1,chromCount)){
  print(chr)
  targets=read.table(paste0(targetdir,"//mutyper.targets.7mer.chr",chr,".nostrict.txt"),header=T)
  colnames(targets) <- c("ancestral7mer","ancestral7merCount")
  spectrum=read.table(paste0(spectrumdir,"chr",chr,"_",population,"_samples.mutyper.spectra.PERPOPULATION.ALLFREQS.NOSTRICT.txt"),header=T)
  spectrum_melt <- melt(spectrum)
  colnames(spectrum_melt) <- c("mutationType","mutationCount")
  spectrum_melt$window <- chr
  spectrum_melt$population <- population
  if(chr %% 2==0){
    spectrum_melt$label <- "TEST"
  } else if(chr %% 2 !=0){
    spectrum_melt$label <- "TRAIN"
  }
  spectrum_melt$ancestral7mer <- substr(spectrum_melt$mutationType,1,7)
  specrum_targets_merged <- merge(spectrum_melt,targets,by="ancestral7mer")
  specrum_targets_merged$mutationCount_divByTargetCount <- specrum_targets_merged$mutationCount/specrum_targets_merged$ancestral7merCount
  allData <- bind_rows(allData,specrum_targets_merged)
}
write.table(allData,paste0(spectrumoutdir,population,"_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.txt"),row.names=F,quote=F,sep="\t")



###### just one mutation type for now: #########
#mutationType="A.T.no_ancestralCpG"
#CVcount=8 # total folds for k-fold cross validation 
######## chosen species/population: ###########
#chosenPop="Mmd"
# (note that random forest doesn't need you to transform the data)
shapes <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/firstOrder_featureTypes_allPossible7mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.txt",header=T,sep="\t")
dim(shapes)
# need to assign row names:
rownames(shapes) <- shapes$motif

# read in:
allData <- read.table(paste0(spectrumoutdir,population,"_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.txt"),header=T)


###### need to sum up TEST windows into one spectrum ###########
############# also want to calculate the fraction of seg sites per window as an alternate outcome #########

# this gets frac of seg sites per window/population/label 
# so it doesn't have data leakage from the other folds within training set .
# but different chroms have different amounts of targets, should keep that in mind
allData <- allData %>%
  group_by(label,window,population) %>%
  mutate(fractionOfSegSites_perWindow = mutationCount/sum(mutationCount))
# merge with shapes

allData$derived7mer <- substr(allData$mutationType,9,15)
allData_intermediate <-merge(allData,shapes,by.x="ancestral7mer",by.y="motif")

allData_withShapes_unprocessed <- merge(allData_intermediate,shapes,by.x="derived7mer",by.y="motif",suffixes=c(".ancestral",".derived"))




####### make your splits into train/test ##########

indices <-
  list(analysis   = which(allData_withShapes_unprocessed$label=="TRAIN"), 
       assessment = which(allData_withShapes_unprocessed$label=="TEST"))

split <- make_splits(indices,allData_withShapes_unprocessed)

# if you want to see what's what:
head(training(split),4)
head(testing(split),4)

########### could sum up training data ###########3
#train_data_unprocessed_allWindowsSummedUp <- training(split)  %>%
#  group_by(mutationType,ancestral7mer) %>%
#  summarise(mutationCount_allTrainingWindows=sum(mutationCount),ancestral7merCount_allTrainingWindows=sum(ancestral7merCount))

#train_data_unprocessed_allWindowsSummedUp$mutationRate_allTrainingWindows <- train_data_unprocessed_allWindowsSummedUp$mutationCount_allTrainingWindows/train_data_unprocessed_allWindowsSummedUp$ancestral7merCount_allTrainingWindows



train_data_cv <- group_vfold_cv(training(split),group=window) # okay so I can make windows based on a grouping variable like so.


# to see a fold
# one fold: train_data_cv[[1]][[1]]
unique(analysis(train_data_cv[[1]][[1]])$window)# the inside windows (chromosomes in this case): all but 7
unique(assessment(train_data_cv[[1]][[1]])$window) # the held out window: 7 
# okay this works great! 
length(unique(analysis(train_data_cv[[1]][[1]])$mutationType))
length((analysis(train_data_cv[[1]][[1]])$mutationType))


######### want to sum up all in-folds and leave out each out fold ? that's hard... not sure if I want to do that. maybe keep folds? ##############

#analysis(train_data_cv[[1]][[1]]) %>%
#  group_by(derived7mer,ancestral7mer,mutationType,label) %>%
#  summarise(mutationCount_allWindows=sum(mutationCount),ancestral7merCount_allWindows=sum(ancestral7merCount),mutationCount_divByTargetCount=mutationCount_allWindows/ancestral7merCount_allWindows)



####### should I sum up each non held-out fold somehow? skip for now. ##########
rand_forest_processing_recipe_OutcomeFracSegSites <- 
  # which consists of the formula (outcome ~ predictors) (don't want to include the 'variable' column)
  recipe(fractionOfSegSites_perWindow ~ .,data=training(split)) %>% # 
  update_role(mutationType, new_role="7mer mutation type label") %>%
  step_rm(mutationCount_divByTargetCount,derived7mer,ancestral7mer, mutationCount,window,label,ancestral7merCount,population) # careful here that nothing else slips in!


rand_forest_processing_recipe_OutcomeFracSegSites %>% summary() 
# should remove the excess features but make sure it does
# add a role for the variable (motif) label that isn't part of the model but keeps track of the label (so it was included in the "." above but now we make it not an outcome or predictor variable so it's not going to be part of the model but will be a label. you see that this works )
rand_forest_processing_recipe_OutcomeFracSegSites

####### 
rand_forest_ranger_model_specs <-
  rand_forest(trees = 1000, mtry = tune(), min_n = tune()) %>% # I added in tree number = 1000
  set_engine('ranger',importance="permutation",respect.unordered.factors="order",verbose=TRUE,num.threads=3) %>%
  set_mode('regression')

######## make a workflow with recipe and model specs ########
rand_forest_workflow <- workflow() %>%
  # add the recipe
  add_recipe(rand_forest_processing_recipe_OutcomeFracSegSites) %>%
  # add the model
  add_model(rand_forest_ranger_model_specs)
rand_forest_workflow

mtrygoal=round(predictorCountdf$totalPredictors/3)
mtrygoal
mtrygrid=c(2,mtrygoal/2,mtrygoal,mtrygoal*2,mtrygoal*3) # range from 2-total features with mtrygoal (feat/3) in the middle
mtrygrid
min_ngrid = c(5,10,20) # not sure about these values but have to start somewhere. default is 5 for regression (minimum node size before gets converted to leaf)
rand_forest_tuning_grid <- expand.grid(mtry=mtrygrid,min_n=min_ngrid)

####### tune: SLOW ############

####### TUNE: **VERY ** slow! 
rf_tune_results <- rand_forest_workflow %>%
  tune_grid(resamples = train_data_cv, #CV object
            grid = rand_forest_tuning_grid, # grid of values to try
            metrics = metric_set(rmse, rsq) # metrics we care about: RMSE and R^2
  )
# started at 3:15 6/4 with different CV folds ; not specifying how bootstrapping is occurring though.
####### could manualy set instead
#paramdf = tibble(mtry=32,min_n=5)
#rand_forest_workflow <- rand_forest_workflow %>%
#  finalize_workflow(paramdf)
rand_forest_workflow
######## try training/testing on just one fold combo (all chr but one in train, tested on one chr)
oneFoldSetToTrainAndAssessOn <- train_data_cv[[1]][[1]]
rand_forest_Fold01_fit <- rand_forest_workflow %>%
  last_fit(oneFoldSetToTrainAndAssessOn) # is very slow (a couple hours because it's based on all 7mer types)

#Export option
# be careful here~! 
############ save model fit for later ##############
save(rand_forest_Fold01_fit,file=paste0(outdir,"rand_forest_Fold01_fit.RData"))
############# load model fit back in ######################
load(paste0(outdir,"rand_forest_Fold01_fit.RData")) # loads it back in using its original name of 
# rand_forest_Fold01_fit 
# and now it just exists again as rand_forest_Fold01_fit! but you had to know what it was called 


rand_forest_Fold01_fit %>% collect_metrics()

rand_forest_Fold01_fit_predictions <- rand_forest_Fold01_fit %>% 
  collect_predictions()

head(rand_forest_Fold01_fit_predictions)

rand_forest_Fold01_fit_predictions$MutationType <- assessment(oneFoldSetToTrainAndAssessOn)$mutationType
rand_forest_Fold01_fit_predictions$population <- assessment(oneFoldSetToTrainAndAssessOn)$population
rand_forest_Fold01_fit_predictions$centralMutationType <- paste0(substr(rand_forest_Fold01_fit_predictions$MutationType,4,4),".",substr(rand_forest_Fold01_fit_predictions$MutationType,12,12))
rand_forest_Fold01_fit_predictions$window <- assessment(oneFoldSetToTrainAndAssessOn)$window

rand_forest_Fold01_fit_predictions_plot <-  ggplot(rand_forest_Fold01_fit_predictions, aes(y=.pred,x=fractionOfSegSites_perWindow,color=centralMutationType,shape=population))+
  geom_point()+
  geom_abline()+
  #scale_x_log10()+
  #scale_y_log10()+
  facet_wrap(~population)+
  ggtitle("OUTCOME: frac of seg sites; ALL 7mer mutation types; random forest, trained on Fold01 ( all odd chrs but one, tested on just one)")
rand_forest_Fold01_fit_predictions_plot

rand_forest_Fold01_fit


# could run on train/test
#rand_forest_TrainTest_fit <- rand_forest_workflow %>%
#  last_fit(split)
##########3 after tuning maybe sum up training windows into one spectrum and train on full spectrum (?) but what about RF bootstraps! #############
# try to train on full training set? 





