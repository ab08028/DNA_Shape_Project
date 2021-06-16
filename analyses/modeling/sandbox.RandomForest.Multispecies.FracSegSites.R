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
tidymodels_prefer() # use this 'tidymodels_prefer()' uses the 'conflicted' package to handle common conflicts with tidymodels and other packages. <-- should add this to my script
set.seed(42) # so results are reproducible 
chromCount=19
# want counts divided by targets
targetdir="/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_target_files"
# per pop spectrum files:
populations=c("Mmd","Ms")
spectrumoutdir="/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_spectrum_target_merged_PerChr/"
dir.create(spectrumoutdir,showWarnings=F)

# modeling outdir:
outdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/sandbox/"
shapes <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/firstOrder_featureTypes_allPossible7mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.txt",header=T,sep="\t")
rownames(shapes) <- shapes$motif

######### only need to do once to make the dataframes for the multipops #######
#allData_multipop <- data.frame() # store multiple species in here
# for(population in populations){
#   spectrumdir=paste0("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_spectrum_files/",population,"/")
#   for(chr in seq(1,chromCount)){
#     print(chr)
#     targets=read.table(paste0(targetdir,"//mutyper.targets.7mer.chr",chr,".nostrict.txt"),header=T)
#     colnames(targets) <- c("ancestral7mer","ancestral7merCount")
#     spectrum=read.table(paste0(spectrumdir,"chr",chr,"_",population,"_samples.mutyper.spectra.PERPOPULATION.ALLFREQS.NOSTRICT.txt"),header=T)
#     spectrum_melt <- melt(spectrum)
#     colnames(spectrum_melt) <- c("mutationType","mutationCount")
#     spectrum_melt$window <- chr
#     spectrum_melt$population <- population
#     if(chr %% 2==0){
#       spectrum_melt$label <- "TEST"
#     } else if(chr %% 2 !=0){
#       spectrum_melt$label <- "TRAIN"
#     }
#     spectrum_melt$ancestral7mer <- substr(spectrum_melt$mutationType,1,7)
#     specrum_targets_merged <- merge(spectrum_melt,targets,by="ancestral7mer")
#     specrum_targets_merged$mutationCount_divByTargetCount <- specrum_targets_merged$mutationCount/specrum_targets_merged$ancestral7merCount
#     allData_multipop <- bind_rows(allData_multipop,specrum_targets_merged)
#   }
# }

#write.table(allData_multipop,paste0(spectrumoutdir,"MULTIPOPULATION_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.txt"),row.names=F,quote=F,sep="\t")
# read in multiple popualtions:
allData_multipop <- read.table(paste0(spectrumoutdir,"MULTIPOPULATION_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.txt"),header=T) # 

# and now want to split train is sp A + B spectra across odd chroms, train is sp A + B spectra across even chroms. want to have species membership as a feature! see if it's VIP or not

allData_multipop <- allData_multipop %>%
  group_by(label,window,population) %>%
  mutate(fractionOfSegSites_perWindow = mutationCount/sum(mutationCount))
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
head(training(split),4)
head(testing(split),4)

train_data_cv <- group_vfold_cv(training(split),group=window) # okay so I can make windows based on a grouping variable like so.

# to see a fold
# one fold: train_data_cv[[1]][[1]]
unique(analysis(train_data_cv[[1]][[1]])$window)# the inside windows (chromosomes in this case): 
unique(assessment(train_data_cv[[1]][[1]])$window) # the held out window: 
# okay this works great! 

########## RECIPE #########
####### should I sum up each non held-out fold somehow? skip for now. ##########
rand_forest_processing_recipe_OutcomeFracSegSites_withPopAsPredictor <- 
  # which consists of the formula (outcome ~ predictors) (don't want to include the 'variable' column)
  recipe(fractionOfSegSites_perWindow ~ .,data=training(split)) %>% # 
  update_role(mutationType, new_role="7mer mutation type label") %>%
  step_rm(mutationCount_divByTargetCount,derived7mer,ancestral7mer, mutationCount,window,label,ancestral7merCount) # KEEPING population in here as a predictor careful here that nothing else slips in!
rand_forest_processing_recipe_OutcomeFracSegSites_withPopAsPredictor %>% summary()



######### MODEL SPECIFICATION #########
rand_forest_ranger_model_specs <-
  rand_forest(trees = 1000, mtry = 32, min_n = 5) %>% # I added in tree number = 1000
  set_engine('ranger',importance="permutation",respect.unordered.factors="order",verbose=TRUE,num.threads=3) %>%
  set_mode('regression')
rand_forest_ranger_model_specs

############ WORKFLOW #############
######## make a workflow with recipe and model specs ########
rand_forest_workflow <- workflow() %>%
  # add the recipe
  add_recipe(rand_forest_processing_recipe_OutcomeFracSegSites_withPopAsPredictor) %>%
  # add the model
  add_model(rand_forest_ranger_model_specs)
rand_forest_workflow

########### TRAIN/TEST MODEL JUST USING ONE FOLD SET #########  
oneFoldSetToTrainAndAssessOn <- train_data_cv[[1]][[1]]

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
rand_forest_Fold01_fit_notlastfit <- fit(rand_forest_workflow,data = analysis(oneFoldSetToTrainAndAssessOn))
# growing trees takes 1/2 hour; permutation importance takes another 1/2 hour.
# then predict based on held out assessment set of fold split:
rand_forest_Fold01_predictions <- predict(object =rand_forest_Fold01_fit_notlastfit, new_data= assessment(oneFoldSetToTrainAndAssessOn))
rand_forest_Fold01_predictions

assessment(oneFoldSetToTrainAndAssessOn) <- rand_forest_Fold01_predictions
############ VIP ###########
ranger_obj <- pull_workflow_fit(rand_forest_Fold01_fit_notlastfit)$fit
ranger_obj

vip_plot <- vip(ranger_obj,include_type = T,num_features = 100)
vip_plot

vip(ranger_obj,include_type = T,num_features = 10)

###### try to find interactions? ############
#vint(ranger_obj,feature_names=c("feature_1_Shear_2.derived","feature_1_Shear_2.ancestral",progress=T),train=analysis(oneFoldSetToTrainAndAssessOn)) # CRASHED

# can't get it to work with categorical variables 
