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
require(doParallel)

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
# for testing on home computer: 
#shapes <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/firstOrder_featureTypes_allPossible7mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.txt",header=T,sep="\t")
rownames(shapes) <- shapes$motif
print('reading in spectrum')
spectrumdir="/net/harris/vol1/home/beichman/DNAShape/spectrumDataForModeling/mouse/"

# note: don't need to read in sequence encoded data -- am doing it below with tidymodels
####### Now including data that has 0 entries and for now has large multi-chromosome windows ############
## eventually make better windows 
allData_multipop <- read.table(paste0(spectrumdir,"TEMPORARY.NEWGROUPS.MULTIPOPULATION_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.INCLUDES0Entries.txt"),header=T) # 
# for testing on home computer:
#allData_multipop <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_spectrum_target_merged_PerChr/TEMPORARY.NEWGROUPS.MULTIPOPULATION_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.INCLUDES0Entries.txt",header=T)

dim(allData_multipop)
# need to get rid of NA entries
print("getting rid of NA mutation rates due to no ancestral targets observed")
allData_multipop <- na.omit(allData_multipop) 
dim(allData_multipop)
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

train_data_cv <- group_vfold_cv(training(split),group=newGroup) # okay so I can make newGroups based on a grouping variable like so.
saveRDS(train_data_cv, file = paste0(outdir,"train_data_cv.rds"))

# to see a fold
# one fold: train_data_cv[[1]][[1]]
unique(analysis(train_data_cv[[1]][[1]])$newGroup)# the inside newGroups (chromosomes in this case): 
unique(assessment(train_data_cv[[1]][[1]])$newGroup) # the held out newGroup: 
# okay this works great! 

###################### RECIPE ###############
# note this recipe is the same as my random forest recipe (just need ot make sure nominal predictors are dummy-fied which I was already doing)
recipe <- 
  # which consists of the formula (outcome ~ predictors) (don't want to include the 'variable' column)
  
  recipe(outcome ~ .,data=training(split)) %>% # 
  update_role(mutationType, new_role="7mer mutation type label") %>%
  step_rm(derived7mer,ancestral7mer, mutationCount,mutationCount_divByTargetCount,newGroup,label,ancestral7merCount) %>%
  step_dummy(all_nominal_predictors()) # KEEPING population in here as a predictor careful here that nothing else slips in! but dummy encoding it;; which RF doesn't need but xgboost and SHAP values does so just doing it ; this also dummy encodes the sequence with A as 0 0 0 , 1 0 0 = c etc. 
recipe
saveRDS(recipe, file = paste0(outdir,"recipe.rds"))

############ MODEL SPECS ##############
# tuning based on: https://juliasilge.com/blog/xgboost-tune-volleyball/
# more info: https://towardsdatascience.com/xgboost-fine-tune-and-optimize-your-model-23d996fab663
boost_tree_xgboost_spec_TUNE <-
  boost_tree(
    trees = 1000, 
    tree_depth = tune(),  # model complexity
    min_n = tune(), # model complexity (number of obs per leaf)
    loss_reduction = tune(), # model compelxity
    sample_size = tune(),  ## randomness (fraction of data to sample for each tree)
    mtry = tune(),         ## randomness
    learn_rate = tune(),  ## step size (how much it learns from previous tree)
  ) %>% 
  #boost_tree() %>% # trying with defaults
  set_engine('xgboost',verbose=T,nthread=25) %>%
  set_mode('regression')
boost_tree_xgboost_spec_TUNE
translate(boost_tree_xgboost_spec_TUNE)
# not yet setting lambda = 1 -- wait on that!
############### TUNE GRID ###############


#"Notice that we had to treat mtry() differently because it depends on the actual number of predictors in the data."
xgb_grid <- grid_latin_hypercube(
  tree_depth(),
  min_n(),
  loss_reduction(),
  sample_size = sample_prop(),
  finalize(mtry(),training(split)), # this sets range of 1,total predictors for mtry
  learn_rate(),
  size = 30
)
# ranges have model specific defaults - what are they?
print('grid of tuning parameters (latin hypercube)')
xgb_grid
dim(xgb_grid)
saveRDS(xgb_grid, file = paste0(outdir,"tuningParameterGrid.rds"))

########### WORKFLOW ##############
xgboost_workflow_TUNE <- workflow() %>%
  # add the recipe
  add_recipe(recipe) %>%
  # add the model
  add_model(boost_tree_xgboost_spec_TUNE)
xgboost_workflow_TUNE

########### TUNE -- SLOW (do on sage) #############
doParallel::registerDoParallel()
print(paste0(Sys.time(),' starting tuning'))
xgb_tuning_results <- tune_grid(
  xgboost_workflow_TUNE,
  resamples = train_data_cv,
  grid = xgb_grid,
  control = control_grid(save_pred = TRUE) # saving predictions
)

xgb_tuning_results
saveRDS(xgb_tuning_results, file = paste0(outdir,"xgb_tuning_results.rds"))

print(paste0(Sys.time(),' finished tuning'))

metrics <- collect_metrics(xgb_tuning_results)
write.table(metrics,paste0(outdir,"tuningMetrics.txt"),sep="\t",rownames=F,quote=F)

best_parameters <- show_best(xgb_tuning_results, "rmse") # try for rmse 
write.table(best_parameters,paste0(outdir,"bestParameters.fromtuning.txt"),sep="\t",rownames=F,quote=F)


sink()
