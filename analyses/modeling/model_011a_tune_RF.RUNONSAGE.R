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
# for testing on home computer: 
#shapes <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/firstOrder_featureTypes_allPossible7mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.txt",header=T,sep="\t")

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

# to see a fold
# one fold: train_data_cv[[1]][[1]]
unique(analysis(train_data_cv[[1]][[1]])$newGroup)# the inside newGroups (chromosomes in this case): 
unique(assessment(train_data_cv[[1]][[1]])$newGroup) # the held out newGroup: 
# okay this works great! 

########## RECIPE #########
####### should I sum up each non held-out fold somehow? skip for now. ##########
rand_forest_processing_recipe <- 
  # which consists of the formula (outcome ~ predictors) (don't want to include the 'variable' column)
  
  recipe(outcome ~ .,data=training(split)) %>% # 
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
  #rand_forest(trees = 1000, mtry = 32, min_n = 5) %>% # I added in tree number = 1000
  rand_forest(trees = 1000, mtry = tune(), min_n = tune()) %>% # I added in tree number = 1000
  set_engine('ranger',importance="permutation",respect.unordered.factors="order",verbose=TRUE,num.threads=15) %>%
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

########### TUNE ##########
doParallel::registerDoParallel()

tune_results <- tune_grid(
  rand_forest_workflow,
  resamples = train_data_cv,
  grid = 20
)

tunePlot <- tune_results %>%
  collect_metrics() %>%
  filter(.metric == "rmse") %>%
  select(mean, min_n, mtry) %>%
  pivot_longer(min_n:mtry,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "rmse")
tunePlot
ggsave(paste0(outdir,"tuningPlot.png"),tunePlot,height=6, width=9)


tune_metrics <- collect_metrics(tune_results)
write.table(tune_metrics,paste0(outdir,"tuningMetrics.txt"),sep="\t",row.names=F,quote=F)


best_tuning_params_rmse <- select_best(tune_results,"rmse")
write.table(best_tuning_params_rmse,paste0(outdir,"best_tuning_params_rmse.txt"),sep="\t",row.names=F,quote=F)


best_tuning_params_rsq <- select_best(tune_results,"rsq")
write.table(best_tuning_params_rsq,paste0(outdir,"best_tuning_params_rsq.txt"),sep="\t",row.names=F,quote=F)
####### take a look at results and possibly narrow them down:
# https://juliasilge.com/blog/sf-trees-random-tuning/


sink()
