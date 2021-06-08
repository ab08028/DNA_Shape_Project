
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




# want counts divided by targets
targetdir="/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_target_files"
# per pop spectrum files:
population="Mmd"
spectrumdir=paste0("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_spectrum_files/",population,"/")
outdir="/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_spectrum_target_merged_PerChr/"
dir.create(outdir,showWarnings=F)
# merge per chromosome # 
# need to separate odd and even
chromCount=19
#oddChr=seq(1,chromCount,by=2)
#evenChr=seq(2,chromCount-1,by=2)
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
write.table(allData,paste0(outdir,population,"_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.txt"),row.names=F,quote=F,sep="\t")


####### okay let's try this ###########

tidymodels_prefer() # use this 'tidymodels_prefer()' uses the 'conflicted' package to handle common conflicts with tidymodels and other packages. <-- should add this to my script
set.seed(42) # so results are reproducible 
outdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/sandbox/"

######## chosen species/population: ###########

# (note that random forest doesn't need you to transform the data)
shapes <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/firstOrder_featureTypes_allPossible7mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.txt",header=T,sep="\t")
dim(shapes)
# need to assign row names:
rownames(shapes) <- shapes$motif

# read in:
allData <- read.table(paste0(outdir,population,"_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.txt"),header=T)

# merge with shapes
allData$centralMutationType <- paste0(substr(allData$mutationType,4,4),".",substr(allData$mutationType,12,12))

allData$derived7mer <- substr(allData$mutationType,9,15)
allData_intermediate <-merge(allData,shapes,by.x="ancestral7mer",by.y="motif")


allData_withShapes_unprocessed <- merge(allData_intermediate,shapes,by.x="derived7mer",by.y="motif",suffixes=c(".ancestral",".derived"))


######### ########## TEMPORARY: RESTRICT TO ONE TYPE OF SHAPE ###########
allData_withShapes_unprocessed <- allData_withShapes_unprocessed[allData_withShapes_unprocessed$centralMutationType %in% "A.T",]
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
#   group_by(mutationType,ancestral7mer) %>%
#   summarise(mutationCount_allTrainingWindows=sum(mutationCount),ancestral7merCount_allTrainingWindows=sum(ancestral7merCount))
# 
# train_data_unprocessed_allWindowsSummedUp$mutationRate_allTrainingWindows <- train_data_unprocessed_allWindowsSummedUp$mutationCount_allTrainingWindows/train_data_unprocessed_allWindowsSummedUp$ancestral7merCount_allTrainingWindows



train_data_cv <- group_vfold_cv(training(split),group=window) # okay so I can make windows based on a grouping variable like so.


# to see a fold
# one fold: train_data_cv[[1]][[1]]
unique(analysis(train_data_cv[[1]][[1]])$window)# the inside windows (chromosomes in this case): all but 7
unique(assessment(train_data_cv[[1]][[1]])$window) # the held out window: 7 
# okay this works great! 
length(unique(analysis(train_data_cv[[1]][[1]])$mutationType))
length((analysis(train_data_cv[[1]][[1]])$mutationType))


xgboost_processing_recipe <- 
  # which consists of the formula (outcome ~ predictors) (don't want to include the 'variable' column)
  recipe(mutationCount_divByTargetCount ~ .,data=training(split)) %>% 
  step_rm(derived7mer,ancestral7mer, mutationCount,window,label,ancestral7merCount,centralMutationType)  %>%
  update_role(mutationType, new_role="7mer mutation type label")  %>%
  step_dummy(all_nominal_predictors())  # don't currently have any categorial vars but need to transform if I get some
xgboost_processing_recipe
### need to recalculate things actually based on the totalSitesAllIntervals rather than frac of seg sites - because should be recalc'd for each fold ###

# add a role for the variable (motif) label that isn't part of the model but keeps track of the label (so it was included in the "." above but now we make it not an outcome or predictor variable so it's not going to be part of the model but will be a label. you see that this works )

boost_tree_xgboost_spec_NOTUNING <-
  #boost_tree(tree_depth = tune(), trees = 1000, learn_rate = tune(), min_n = tune(), loss_reduction = tune()) %>%
  boost_tree() %>% # trying with defaults
  set_engine('xgboost',verbose=T,lambda=1,nthread=3) %>%
  set_mode('regression')
boost_tree_xgboost_spec
translate(boost_tree_xgboost_spec_NOTUNING)

boost_tree_xgboost_spec_TUNE <-
  boost_tree(tree_depth = tune(), trees = 1000, learn_rate = tune(), min_n = tune(), loss_reduction = tune()) %>%
  #boost_tree() %>% # trying with defaults
  set_engine('xgboost',verbose=T,lambda=1,nthread=3) %>%
  set_mode('regression')
boost_tree_xgboost_spec
translate(boost_tree_xgboost_spec_TUNE)
# other things you could tune: sample_size = tune()

######## make a workflow with recipe and model specs ########
xgboost_workflow_noTUNE <- workflow() %>%
  # add the recipe
  add_recipe(xgboost_processing_recipe) %>%
  # add the model
  add_model(boost_tree_xgboost_spec_NOTUNING)
xgboost_workflow_noTUNE

######### set up parameter tuning grid #######
# grid specification
xgboost_params <- 
  dials::parameters(
    min_n(),
    tree_depth(),
    learn_rate(),
    loss_reduction()
  )

# Next we set up the grid space. The dails::grid_* functions support several methods for defining the grid space. We are using the dails::grid_max_entropy() function which covers the hyperparameter space such that any portion of the space has an observed combination that is not too far from it.
xgboost_grid <- 
  dials::grid_max_entropy(
    xgboost_params, 
    size = 60
  )
#knitr::kable(head(xgboost_grid))
head(xgboost_grid)

xgboost_tune_results <- xgboost_workflow %>%
  tune_grid(resamples = train_data_cv, #CV object
            grid = xgboost_grid, # grid of values to try
            metrics = metric_set(rmse, rsq) # metrics we care about: RMSE and R^2
  )
# started at 1:00 on 6/7
xgboost_tune_results 


# trying on the first fold (10 chr in, 1 out)
xgboost_defaults_results_Fold01 <- xgboost_workflow %>%
  last_fit(train_data_cv[[1]][[1]])
xgboost_defaults_results_performance_Fold01 <- xgboost_defaults_results_Fold01 %>% collect_metrics()
xgboost_defaults_results_performance_Fold01

# You can also extract the test set predictions themselves using the collect_predictions() function. Note that there are 192 rows in the predictions object below which matches the number of test set observations (just to give you some evidence that these are based on the test set rather than the training set).

xgboost_defaults_results_predictions_Fold01 <- xgboost_defaults_results_Fold01 %>% collect_predictions()

xgboost_defaults_results_predictions_Fold01$mutationType <- assessment(train_data_cv[[1]][[1]])$mutationType
xgboost_defaults_results_predictions_Fold01$window <- assessment(train_data_cv[[1]][[1]])$window
xgboost_defaults_results_predictions_Fold01$centralMutationType <- assessment(train_data_cv[[1]][[1]])$centralMutationType

head(xgboost_defaults_results_predictions_Fold01)

ggplot(xgboost_defaults_results_predictions_Fold01,aes(y=.pred,x=mutationCount_divByTargetCount,color=centralMutationType,shape=as.factor(window)))+geom_point()+geom_abline()+scale_y_log10()+scale_x_log10()



xgboost_defaults_results_Fold02 <- xgboost_workflow %>%
  last_fit(train_data_cv[[1]][[2]])
xgboost_defaults_results_performance_Fold02 <- xgboost_defaults_results_Fold02 %>% collect_metrics()
xgboost_defaults_results_performance_Fold02


xgboost_defaults_results_predictions_Fold02 <- xgboost_defaults_results_Fold02 %>% collect_predictions()
ggplot(xgboost_defaults_results_predictions_Fold02,aes(y=.pred,x=mutationCount_divByTargetCount))+geom_point()+geom_abline()+scale_y_log10()+scale_x_log10()

