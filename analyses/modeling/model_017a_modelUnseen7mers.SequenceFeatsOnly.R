############## want the full spectrum of each species ################### 
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

outdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210714_model_017a_modelUnseen7mers.SequenceFeatsOnly_ranOnLaptop/"
dir.create(outdir,showWarnings = F)
############## read in shapes ##############
# for testing on home computer: 
shapes <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/METHYLATEDfirstOrder_featureTypes_allPossible7mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.METHYLATED.txt",header=T,sep="\t")
rownames(shapes) <- shapes$motif


##### read in spectra and sum up fully over all train and test chromosomes #######
allData_multipop <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_spectrum_target_merged_PerChr/TEMPORARY.NEWGROUPS.MULTIPOPULATION_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.INCLUDES0Entries.txt",header=T)

################ condense spectrum: #####################
allData_multipop <- allData_multipop %>%
  group_by(mutationType,population,ancestral7mer) %>%
  summarise(mutationCount_AllChrs=sum(mutationCount),ancestral7merCount_AllChrs=sum(ancestral7merCount)) # this is ok for summing up ancestral 7mers because you grouped by mutation type so even though they appear mult times in df they are getting summed approporiately 
dim(allData_multipop)
# all 7mers just appear once for eah species. no train/test chromosomes or even/odd bp.
allData_multipop$derived7mer <- substr(allData_multipop$mutationType,9,15)
allData_multipop_intermediate <-merge(allData_multipop,shapes,by.x="ancestral7mer",by.y="motif")

allData_withShapes_unprocessed <- merge(allData_multipop_intermediate,shapes,by.x="derived7mer",by.y="motif",suffixes=c(".ancestral",".derived"))


######## add seq info (not yet 1-hot encoded, will encode with tidy models )##########
# just want ancestral 7mer 
allData_withShapes_unprocessed <- allData_withShapes_unprocessed %>% separate(ancestral7mer,into=c("Pos_1","Pos_2","Pos_3","Pos_4.ancestral","Pos_5","Pos_6","Pos_7"),remove=F,sep=c(1,2,3,4,5,6))
# also get derived bp:
allData_withShapes_unprocessed$Pos_4.derived <- substr(allData_withShapes_unprocessed$derived7mer,4,4)


###### DO NOT RESCALE OUTCOME >>> WOULD RESULT IN DATA LEAKAGE #######
# because can't rescale per fold (per chrom)
# just divide to get mutation rate:
allData_withShapes_unprocessed$outcome <- allData_withShapes_unprocessed$mutationCount_AllChrs / allData_withShapes_unprocessed$ancestral7merCount_AllChrs
head(allData_withShapes_unprocessed) ####### outcome is just raw mutation rate
################# split the data ##############
# keep 3/4 of 7mers in model

#### *want* to randomly pick 1/2 of 7mers but have it be the *same* 7mers in both populations
length(unique(allData_withShapes_unprocessed$mutationType)) # just the 7mers
# pick out random set of 7mers in one species and assign the mutation types to train/test for both species

training1Species <- allData_withShapes_unprocessed %>%
  subset(population=="Mmd") %>%
  sample_frac(0.75) %>% ### NOTE 3/4 HERE
  mutate(label="TRAINING")

training7merlist = training1Species$mutationType
length(training7merlist)

# and now can label these 7mers as training and rest as test and willb e the same for both species!

allData_withShapes_unprocessed$label <- "TEST" 
allData_withShapes_unprocessed[allData_withShapes_unprocessed$mutationType %in% training7merlist,]$label <- "TRAIN"

# make sure same for both pops:
length(allData_withShapes_unprocessed[allData_withShapes_unprocessed$population=="Mmd" & allData_withShapes_unprocessed$label=="TRAIN",]$mutationType)

setequal(allData_withShapes_unprocessed[allData_withShapes_unprocessed$population=="Mmd" & allData_withShapes_unprocessed$label=="TRAIN",]$mutationType , allData_withShapes_unprocessed[allData_withShapes_unprocessed$population=="Ms" & allData_withShapes_unprocessed$label=="TRAIN",]$mutationType) # checks if are equal 


setequal(allData_withShapes_unprocessed[allData_withShapes_unprocessed$population=="Mmd" & allData_withShapes_unprocessed$label=="TEST",]$mutationType , allData_withShapes_unprocessed[allData_withShapes_unprocessed$population=="Ms" & allData_withShapes_unprocessed$label=="TEST",]$mutationType) # checks if are equal 

sum(allData_withShapes_unprocessed[allData_withShapes_unprocessed$label=="TEST",]$mutationType %in% allData_withShapes_unprocessed[allData_withShapes_unprocessed$label=="TRAIN",]$mutationType) # 0 good
  
head(allData_withShapes_unprocessed)
# now to make the split, I am going to split on TRAIN/TEST
indices <-
  list(analysis   = which(allData_withShapes_unprocessed$label=="TRAIN"), 
       assessment = which(allData_withShapes_unprocessed$label=="TEST"))

split <- make_splits(indices,allData_withShapes_unprocessed)
saveRDS(split, file = paste0(outdir,"split.rds"))

# if you want to see what's what:
head(training(split),4)
head(testing(split),4)

# not using yet:
train_data_cv <- vfold_cv(training(split),v=5) 
train_data_cv

########## recipe -- for this model want SEQUENCE FEATS ONLY #########
rand_forest_processing_recipe <- 
  # which consists of the formula (outcome ~ predictors) (don't want to include the 'variable' column)
  
  recipe(outcome ~ .,data=training(split)) %>% # 
  update_role(mutationType, new_role="7mer mutation type label") %>%
  step_rm(derived7mer,ancestral7mer, mutationCount_AllChrs,ancestral7merCount_AllChrs,label) %>%
  step_rm(starts_with("feature_")) %>%
  step_dummy(all_nominal_predictors())

rand_forest_processing_recipe
rand_forest_processing_recipe %>% summary()

# these are what's left after the recipe (note that mutation type isn't a predictor jsut a label so isn't used for prediction in the model)

prep(rand_forest_processing_recipe) %>%
  juice() %>%
  names() 
# make sure it removed it:


saveRDS(rand_forest_processing_recipe,file=paste0(outdir,"recipe.rds"))
######### MODEL SPECIFICATION #########
rand_forest_ranger_model_specs <-
  ### NOTE CHANGE HERE OF Min_N to 5
  # NOTE: only 25 vars here so mtry sould be ~8
  rand_forest(trees = 1000, mtry = 8, min_n = 5) %>% # I added in tree number = 1000
  set_engine('ranger',importance="permutation",respect.unordered.factors="order",verbose=TRUE) %>%
  set_mode('regression')
rand_forest_ranger_model_specs

########## WORKFLOW #############
######## make a workflow with recipe and model specs ########
rand_forest_workflow <- workflow() %>%
  # add the recipe
  add_recipe(rand_forest_processing_recipe) %>%
  # add the model
  add_model(rand_forest_ranger_model_specs)
rand_forest_workflow

########## tuning grid ##########

# rf_grid <- grid_regular(
#   mtry(range = c(10, 50)),
#   min_n(range = c(1, 8)),
#   levels = 5
# )
# dim(rf_grid)
# rf_grid
# ######### try to tune #########
# doParallel::registerDoParallel()
# 
# 
# tune_results_regular_grid <- tune_grid(
#   rand_forest_workflow_tune,
#   resamples = train_data_cv,
#   grid = rf_grid
# )
# 
# tune_results_regular_grid
####### LAST FIT ON SPLIT ############ 
last_rf_fit <- 
  rand_forest_workflow %>% 
  last_fit(split)
last_rf_fit 
saveRDS(last_rf_fit, file = paste0(outdir,"FitOfFullModelToHeldOut7mers.rds"))


predictions <- last_rf_fit %>%
  collect_predictions()

metrics <- last_rf_fit %>%
  collect_metrics()

metrics
last_rf_fit$.workflow # oob metrics

write.table(metrics,paste0(outdir,'metrics.txt'),row.names=F,quote=F,sep="\t")
  
predictions_withData <- bind_cols(predictions,testing(split)[,c("population","mutationType")])
predictions_withData$centralMutationType <- paste0(substr(predictions_withData$mutationType,4,4),".",substr(predictions_withData$mutationType,12,12))

saveRDS(predictions_withData, file = paste0(outdir,"ModelPredictions.WithData.rds"))

dim(predictions_withData)
head(predictions_withData)

#### get rsq per type ######
rsqsPerSpeciesAndMutationType <- predictions_withData %>%
  group_by(centralMutationType,population) %>%
  rsq(truth=outcome,estimate=.pred)
rsqsPerSpeciesAndMutationType
write.table(rsqsPerSpeciesAndMutationType,paste0(outdir,"rsq.percentraltype.txt"),quote = F,row.names=F,sep="\t")


# make a label that has rsq in the central mutation type
predictions_withData_WithPerMutationTypeMetrics <- merge(predictions_withData,rsqsPerSpeciesAndMutationType,by=c("population","centralMutationType"))


########## combine pop and mutaiton type and rsq into one label to facet over 
predictions_withData_WithPerMutationTypeMetrics$comboLabel <- paste0(predictions_withData_WithPerMutationTypeMetrics$population,"\n",predictions_withData_WithPerMutationTypeMetrics$centralMutationType,"\nrsq = ",round(predictions_withData_WithPerMutationTypeMetrics$.estimate,3))

head(predictions_withData_WithPerMutationTypeMetrics)

plotWithRsq <-  ggplot(predictions_withData_WithPerMutationTypeMetrics, aes(y=.pred,x=outcome,color=centralMutationType))+
  geom_point()+
  geom_abline()+
  #geom_text(data=rsqsPerSpeciesAndMutationType,aes(x=5e-3,y=4e-3,label=round(.estimate,4)),color="black")+
  facet_wrap(~comboLabel,scales="free",ncol=2,dir = "v")+
  #ggtitle(paste0(description,"\ntrained on Fold01\n(all odd Windows but one, tested on just Window",toString(unique(truth_prediction_df$newGroup)),")"))+
  theme_bw() +
  theme(text=element_text(size=14))+
  ggtitle("trained on 1/3 of 7mers across species, tested on 1/2 of 7mers\n min_n =5 to try to replicate bad behavior\nof OOB from my early models")
  #scale_x_log10()+
  #scale_y_log10()
plotWithRsq

ggsave(paste0(outdir,"plotWithRsq.png"),plotWithRsq,height=13,width=9)


#### want to plot species differences ######### 
predictions_withData_spread <- pivot_wider(predictions_withData[,c("mutationType","population","outcome",".pred","centralMutationType")],id_cols = c(mutationType,population,centralMutationType),names_from=population,values_from=c(outcome,.pred)) 

head(predictions_withData_spread)
# this is really useful: plot the difference between the two species and how well the predictions do (species are off of y=x line because of diff muttion rates; model picks that up! super cool)
speciesComparisonPlot <- ggplot(predictions_withData_spread,aes(x=outcome_Mmd,y=outcome_Ms))+
  geom_point()+
  geom_point(aes(x=.pred_Mmd,y=.pred_Ms),shape=1,color="red")+
  geom_abline()+
  facet_wrap(~centralMutationType,scales="free")
speciesComparisonPlot
ggsave(paste0(outdir,"SpeciesComparison.png"),speciesComparisonPlot,height=5,width=7)
