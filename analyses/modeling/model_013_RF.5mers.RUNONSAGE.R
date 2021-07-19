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

# running on laptop for now
args <- commandArgs(trailingOnly = TRUE)
description <- args[1]
outdir <- args[2]
#outdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210714_model_013_RF.5mers_runOnLaptop/"
dir.create(outdir,showWarnings = F)
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
shapes <- read.table(paste0(shapedir,"firstOrder_featureTypes_allPossible5mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.5mer.txt"),header=T,sep="\t")

# for testing on home computer: 
#shapes <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/firstOrder_featureTypes_allPossible5mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.5mer.txt",header=T,sep="\t")
rownames(shapes) <- shapes$motif

print('reading in spectrum')
spectrumdir="/net/harris/vol1/home/beichman/DNAShape/spectrumDataForModeling/mouse/"

# note: don't need to read in sequence encoded data -- am doing it below with tidymodels
####### Now including data that has 0 entries and for now has large multi-chromosome windows ############
## eventually make better windows : this is the 7mer spectrum but am going to sum up over 5mers below so don't worry :) 
allData_multipop <- read.table(paste0(spectrumdir,"TEMPORARY.NEWGROUPS.MULTIPOPULATION_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.INCLUDES0Entries.txt"),header=T) # this is the 7mer spectrum but down below you are going to sum up over 5mers 
# for testing on home computer:
#allData_multipop <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_spectrum_target_merged_PerChr/TEMPORARY.NEWGROUPS.MULTIPOPULATION_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.INCLUDES0Entries.txt",header=T)

dim(allData_multipop)
# need to get rid of NA entries
print("getting rid of NA mutation rates due to no ancestral targets observed")
allData_multipop <- na.omit(allData_multipop) 
dim(allData_multipop)
# and now want to split train is sp A + B spectra across odd chroms, train is sp A + B spectra across even chroms. want to have species membership as a feature! see if it's VIP or not

########## !!!!!! condense 7mer spectrum down to 5mer spectrum  !!!!!!! #########
allData_multipop$mutationType_5mer <- paste0(substr(allData_multipop$mutationType,2,6),".",substr(allData_multipop$mutationType,10,14))
# doing some debuggin:
#allData_multipop$ancestral5mer <- substr(allData_multipop$mutationType,2,6)
#allData_multipop[allData_multipop$ancestral5mer=="CGACG",]
#allData_multipop %>%
#  group_by(population,newGroup,mutationType) %>%
#  subset(ancestral5mer=="CGACG") %>%
#  summarise(test_Ancestral5merCount=sum(ancestral7merCount)) %>%
#  View()
  
head(allData_multipop)
tail(allData_multipop)

allData_multipop <- allData_multipop %>%
  group_by(population,newGroup,label,mutationType_5mer) %>%
  summarise(mutationCount_5mer = sum(mutationCount),ancestral5merCount=sum(ancestral7merCount),mutationCount_divByTargetCount_5mer = mutationCount_5mer/ancestral5merCount) ##### is this double counting 5mer targets?


dim(allData_multipop)
dim(allData_multipop[allData_multipop$population=="Mmd" & allData_multipop$newGroup==1,]) # there are 1536 types of central 5mer so this checks out.

# note sum of all ancestral targest across mutations is larger than the available space because each mutation type has some of the same ancestral targets, so that's why it's important to group by a particular mutation type 
# make sure these sums work 
  # newGroup here not window; eventually go back to window?

########## merge with shapes ##########
# need to get ancestral and derived 5mers
# be careful here to do derived/ancestral properly
allData_multipop$ancestral5mer <- substr(allData_multipop$mutationType_5mer,1,5)
allData_multipop$derived5mer <- substr(allData_multipop$mutationType_5mer,7,11)
head(allData_multipop[,c("mutationType_5mer","ancestral5mer","derived5mer")])


allData_multipop_intermediate <-merge(allData_multipop,shapes,by.x="ancestral5mer",by.y="motif")

allData_withShapes_unprocessed <- merge(allData_multipop_intermediate,shapes,by.x="derived5mer",by.y="motif",suffixes=c(".ancestral",".derived"))

head(allData_withShapes_unprocessed)

######## add seq info (not yet 1-hot encoded, will encode with tidy models )##########
# just want ancestral 5mer : so now it's Pos3 that is central position (not pos 4 the way it was for 7mers)
allData_withShapes_unprocessed <- allData_withShapes_unprocessed %>% separate(ancestral5mer,into=c("Pos_1","Pos_2","Pos_3.ancestral","Pos_4","Pos_5"),remove=F,sep=c(1,2,3,4))
head(allData_withShapes_unprocessed)

# also get derived bp:
allData_withShapes_unprocessed$Pos_3.derived <- substr(allData_withShapes_unprocessed$derived5mer,3,3) # it's position 3 not 4 now

head(allData_withShapes_unprocessed)
# check random rows; allData_withShapes_unprocessed[c(45,600,30),]
########### want to rescale outcome variable rate so it's relative #######
# not adding epsilon because 
# not log scaling # call whatever you want your final outcome to be 'outcome' so that it's the same in all plots
# picked this because it's smaller than 1/genome size 
allData_withShapes_unprocessed <- allData_withShapes_unprocessed %>%
  group_by(population,newGroup,label) %>%
  mutate(outcome=mutationCount_divByTargetCount_5mer/sum(mutationCount_divByTargetCount_5mer)) # didn't change the name of this to reflect 5mers so this code is good to go 

# this maintains the ranking but rescales the outcome ; try it! 
sum(allData_withShapes_unprocessed[allData_withShapes_unprocessed$population=="Mmd" & allData_withShapes_unprocessed$newGroup==1,]$outcome) # should be 1
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

########### TRAIN/TEST MODEL JUST USING ONE FOLD SET #########  
oneFoldSetToTrainAndAssessOn <- train_data_cv[[1]][[1]]
saveRDS(oneFoldSetToTrainAndAssessOn, file = paste0(outdir,"oneFoldSetToTrainAndAssessOn.rds"))
# load back in:
#oneFoldSetToTrainAndAssessOn <- readRDS(paste0(outdir,"oneFoldSetToTrainAndAssessOn.rds"))
# rand_forest_Fold01_fit <- rand_forest_workflow %>%
#   last_fit(oneFoldSetToTrainAndAssessOn) # is very slow (a couple hours because it's based on all 7mer types)

########## RECIPE #########
####### should I sum up each non held-out fold somehow? skip for now. ##########
# Note about recipes: You initialize the recipe with some data so that it gets the column names, but then when I do fit, I specify which data to use for fitting the model and it's not data that is in the recipe (so the recipe was initialized with the full trainign split, but then model is fit with the fold I am using, I got worried that somehow there could be some data leakage with that, but I don't think there is because the model used to initialize the recipe isn't used for fitting if I don't say to. But just to be safe in future models maybe initialize the recipe with the fold you're using # So for this model on 
# see variables:
names(analysis(oneFoldSetToTrainAndAssessOn))
rand_forest_processing_recipe <- 
  # which consists of the formula (outcome ~ predictors) (don't want to include the 'variable' column)
  
  recipe(outcome ~ .,data=analysis(oneFoldSetToTrainAndAssessOn)) %>% # 
  update_role(mutationType_5mer, new_role="5mer mutation type label") %>%
  step_rm(derived5mer,ancestral5mer,mutationCount_5mer,mutationCount_divByTargetCount_5mer,newGroup,label,ancestral5merCount) %>%
  step_dummy(all_nominal_predictors()) # KEEPING population in here as a predictor careful here that nothing else slips in! but dummy encoding it;; which RF doesn't need but xgboost and SHAP values does so just doing it ; this also dummy encodes the sequence with A as 0 0 0 , 1 0 0 = c etc. 
#Encoding each ancestral bp as a feature Pos_1, Pos_2, Pos_3 etc.
#Pos_4.ancestral and Pos_4.derived
#Then convert to one-hot binary variables where A is 0 0 0 , C is 1 0 0 , G is 0 1 0, T is 0 0 1 (ala L&S)
#But unlike L&S I’m also keeping in the derived position of the central bp

rand_forest_processing_recipe %>% summary()
rand_forest_processing_recipe

# extract it just to see:
print('processing training data to make sure variables I want are there')
juicedData <- prep(rand_forest_processing_recipe,training(oneFoldSetToTrainAndAssessOn)) %>%
  juice()
dim(juicedData)   # so there are ~57 features now so should adjust mtry
names(juicedData)
# can extract data like this try it :
#prepped <- prep(rand_forest_processing_recipe,training(split)) %>%
#  juice()

######### MODEL SPECIFICATION #########
rand_forest_ranger_model_specs <-
  rand_forest(trees = 1000, mtry = 19, min_n = 1) %>% # I added in tree number = 1000
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
# CAREFUL HERE: need to make sure subsetting correclty now that it's 5mers
# (you are)
truth_prediction_df$centralMutationType <- paste0(substr(truth_prediction_df$mutationType_5mer,3,3),".",substr(truth_prediction_df$mutationType_5mer,9,9))
# check it:
tail(truth_prediction_df[,c("mutationType_5mer","centralMutationType")])
# note this doesn't have the 1-coded pops unless you juice() it but rows are still in same order 
rand_forest_Fold01_fit_predictions_plot <-  ggplot(truth_prediction_df, aes(y=.pred,x=outcome,color=centralMutationType,shape=population))+
  geom_point()+
  geom_abline()+
  facet_wrap(~population)+
  #ggtitle(paste0(description,"\ntrained on Fold01\n(all odd windows but one, tested on just window ",toString(unique(truth_prediction_df$newGroup)),")"))+
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


# make a label that has rsq in the central mutation type
predictions_withData_WithPerMutationTypeMetrics <- merge(truth_prediction_df,rsqsPerSpeciesAndMutationType,by=c("population","centralMutationType"))


########## combine pop and mutaiton type and rsq into one label to facet over 
predictions_withData_WithPerMutationTypeMetrics$comboLabel <- paste0(predictions_withData_WithPerMutationTypeMetrics$population,"\n",predictions_withData_WithPerMutationTypeMetrics$centralMutationType,"\nrsq = ",round(predictions_withData_WithPerMutationTypeMetrics$.estimate,3))

rand_forest_Fold01_fit_predictions_plot_faceted <-  ggplot(predictions_withData_WithPerMutationTypeMetrics, aes(y=.pred,x=outcome,color=centralMutationType))+
  geom_point()+
  geom_abline()+
  facet_wrap(~comboLabel,ncol=2,dir="v",scales="free")+
  #ggtitle(paste0(description,"\ntrained on Fold01\n(all odd Windows but one, tested on just Window",toString(unique(truth_prediction_df$newGroup)),")"))+
  theme_bw()
rand_forest_Fold01_fit_predictions_plot_faceted

ggsave(paste0(outdir,"modelTrainedOnOneFold.PredictionsPlot.FacetedPerMutationType.AssessedOnWindow",toString(unique(truth_prediction_df$newGroup)),".png"),rand_forest_Fold01_fit_predictions_plot_faceted,height=12,width=9)

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
truth_prediction_df_spread <- pivot_wider(truth_prediction_df[,c("mutationType_5mer","population","outcome",".pred","centralMutationType")],id_cols = c(mutationType_5mer,population,centralMutationType),names_from=population,values_from=c(outcome,.pred)) 

head(truth_prediction_df_spread)
# this is really useful: plot the difference between the two species and how well the predictions do (species are off of y=x line because of diff muttion rates; model picks that up! super cool)
speciesComparisonPlot <- ggplot(truth_prediction_df_spread,aes(x=outcome_Mmd,y=outcome_Ms))+
  geom_point()+
  geom_point(aes(x=.pred_Mmd,y=.pred_Ms),shape=1,color="red")+
  geom_abline()+
  facet_wrap(~centralMutationType,scales="free")+
  theme_bw()
speciesComparisonPlot
ggsave(paste0(outdir,"modelTrainedOnOneFold.SpeciesXYComparison.AssessedOnWindow",toString(unique(truth_prediction_df$newGroup)),".png"),speciesComparisonPlot,height=5,width=7)

######## SHAP VALUES ###########
library(doParallel)
registerDoParallel()

######### Need to juice the training and assessment data and exclude anything that isn't features ###### 

# juice and get rid of outcome variables and anything that isn't predictions (like mutaiton type)
Xtrain <- prep(rand_forest_processing_recipe, training(oneFoldSetToTrainAndAssessOn)) %>% 
  juice() %>% 
  select(-c(mutationType_5mer,outcome)) %>% 
  as.matrix()
head(Xtrain)


Xtest <- prep(rand_forest_processing_recipe, testing(oneFoldSetToTrainAndAssessOn)) %>% 
  juice() %>% 
  select(-c(mutationType_5mer,outcome)) %>% 
  as.matrix()
head(Xtest)

# function describing how to get predictions from the model 
pfun <- function(object, newdata) {
  predict(object, data = newdata)$predictions
}


# a bit slow: 

shap <- fastshap::explain(ranger_obj, X = Xtrain,pred_wrapper=pfun,newdata=Xtest) ##
saveRDS(shap,file=paste0(outdir,"SHAPResults.rds")) # 1536 *2 because one obs per species
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
Xdftrain <- prep(rand_forest_processing_recipe, training(oneFoldSetToTrainAndAssessOn)) %>% 
  juice()  #%>%


# note X above was a matrix but this has to be a df
Xdftest <- prep(rand_forest_processing_recipe, testing(oneFoldSetToTrainAndAssessOn)) %>% 
  juice()  #%>%

# label categorically for nicer plots ; jsut make sure you get it right 
Xdftest$populationLabel <- ""
Xdftest[Xdftest$population_Ms==0,]$populationLabel <- "Mmd"
Xdftest[Xdftest$population_Ms==1,]$populationLabel <- "Ms"
Xdftest_melt <- melt(Xdftest,id.vars=c("mutationType_5mer","populationLabel")) # are these based on train or test?

dim(Xdftest)

head(Xdftest)



featuresToPlot=names(shap)
dir.create(paste0(outdir,"perFeatureDependencePlots/"),showWarnings = F)
for(feature in featuresToPlot){
  
  shapdependenceplot1 <- autoplot(shap, 
                                  type = "dependence", 
                                  feature = feature, 
                                  X = Xdftest, # this hsould be test data if you ran shap with newdata=testing data! sholud be train if you just ran it without any new data
                                  smooth = TRUE, color_by = "populationLabel")+ # going to try coloring by species membership
    scale_color_manual(values=c("orange","steelblue"))+
    ggtitle(feature)+
    theme_bw()+
    labs(color="population")
  shapdependenceplot1
  ggsave(paste0(outdir,"/perFeatureDependencePlots/random_forest.shapDependence.",feature,".png"),shapdependenceplot1,height=4,width=6)
}


####### SHAP summary plot (really useful) ######
shap_labeled <- shap
# need to label rows by their mutation type and pop observed in
# note that shap$population_Ms is the shap values for that feature
shap_labeled$mutationType_5mer <- Xdftest$mutationType_5mer
shap_labeled$populationLabel <- Xdftest$populationLabel




shap_melt <- melt(shap_labeled,id.vars = c("mutationType_5mer","populationLabel")) # dims are 1536*2 pops * 57 features = 175104
head(shap_melt)
colnames(shap_melt) <- c("mutationType_5mer","populationLabel","feature","SHAP.value")
# combine with xdf:
#Xdftrain_melt <- melt(Xdftrain,id.vars =c("mutationType_5mer","population_Ms","newGroup")) # are these based on train or test?
#head(Xdftrain_melt)
#colnames(Xdftrain_melt) <- c("mutationType_5mer","feature","feature.value")

dim(Xdftest_melt)
head(Xdftest_melt)
colnames(Xdftest_melt) <- c("mutationType_5mer","populationLabel","feature","feature.value")
# still includes the 'outcome column fyi -- doesn't mess anything up tho

shap_plusFeatures <- merge(shap_melt,Xdftest_melt[Xdftest_melt$feature!="outcome",],by=c("mutationType_5mer","feature","populationLabel")) # this merge created duplicates because I wasn't labelling points by what population they came from; fixed it
dim(shap_plusFeatures) # yes now this is right with no dups
###### do min max normalization of originalf feature values (for ease of plotting): 
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
  ggtitle(paste0("SHAP summary plot"))+
  scale_color_viridis_c(labels=c("low","high"),breaks=c(0,1),limits=c(0,1),option = "plasma")
shapSummaryPlot
ggsave(paste0(outdir,"random_forest.shapSummaryPlot.USEFUL.png"),shapSummaryPlot,height=12,width=10)

# just do top 10 most important:
top10 <- head(featuresInOrderOfImportance$Variable,10)
shapSummaryPlot_top10 <- ggplot(shap_plusFeatures[shap_plusFeatures$feature %in% top10,],aes(y=feature,x=SHAP.value,color=minmaxnormalized.feature.value))+
  #geom_point()+
  geom_quasirandom(groupOnX = F) + # bee swarm so ti's like a violin plot; want to cluster on y though
  
  theme_bw()+
  ggtitle(paste0("SHAP summary plot"))+
  scale_color_viridis_c(labels=c("low","high"),breaks=c(0,1),limits=c(0,1),option = "plasma")
shapSummaryPlot_top10
ggsave(paste0(outdir,"random_forest.shapSummaryPlot.top10.USEFUL.png"),shapSummaryPlot_top10,height=4,width=10)

### just do all the population ones: 
justPopulationRelatedFeatures=shap_plusFeatures[grepl("population",shap_plusFeatures$feature),]

shapSummaryPlot_populationOnly <- ggplot(justPopulationRelatedFeatures,aes(y=feature,x=SHAP.value,color=minmaxnormalized.feature.value))+
  #geom_point()+
  geom_quasirandom(groupOnX = F) + # bee swarm so ti's like a violin plot; want to cluster on y though
  
  theme_bw()+
  ggtitle(paste0("SHAP summary plot"))+
  scale_color_viridis_c(labels=c("low","high"),breaks=c(0,1),limits=c(0,1),option = "plasma")
shapSummaryPlot_populationOnly
ggsave(paste0(outdir,"random_forest.shapSummaryPlot.top10.USEFUL.png"),shapSummaryPlot_top10,height=4,width=10)

### label outliers #########
shapSummaryPlot_populationOnly_labelOutliers <- ggplot(justPopulationRelatedFeatures,aes(y=feature,x=SHAP.value,color=minmaxnormalized.feature.value))+
  #geom_point()+
  geom_quasirandom(groupOnX = F) + # bee swarm so ti's like a violin plot; want to cluster on y though
  
  theme_bw()+
  ggtitle(paste0("SHAP summary plot"))+
  scale_color_viridis_c(labels=c("low","high"),breaks=c(0,1),limits=c(0,1),option = "plasma")+
  geom_label_repel(data=subset(justPopulationRelatedFeatures,abs(SHAP.value)>=2.5e-4),aes(label=mutationType_5mer),color="black",  size = 2,
                   box.padding = unit(0.25, "lines"),
                   point.padding = unit(0.5, "lines"),force = T,nudge_y = 0.1)+
  facet_wrap(~populationLabel)

shapSummaryPlot_populationOnly_labelOutliers
ggsave(paste0(outdir,"random_forest.shapSummaryPlot.top10.USEFUL.png"),shapSummaryPlot_top10,height=4,width=10)


##### want to get 5mers that are outliers for population SHAP values ######
# not sure about this cutoff -- should automate 
populationSHAPoutlier_5mers_list <- subset(justPopulationRelatedFeatures,abs(SHAP.value)>=2.5e-4)$mutationType_5mer
head(populationSHAPoutlier_5mers_list)

# want to get heatmap: 
truth_prediction_df_spread$truth_oddsRatio_mmd_ms_rescaledMutationRate <- truth_prediction_df_spread$outcome_Mmd / truth_prediction_df_spread$outcome_Ms


truth_prediction_df_spread$modelPrediction_oddsRatio_mmd_ms_rescaledMutationRate <- truth_prediction_df_spread$.pred_Mmd / truth_prediction_df_spread$.pred_Ms

truth_prediction_df_spread$outlierLabel <- "doesn't have a high SHAP value for population status"
truth_prediction_df_spread[truth_prediction_df_spread$mutationType_5mer %in% populationSHAPoutlier_5mers_list,]$outlierLabel <- "population outlier from SHAP values"

sandbox_predictionsOutliersPlot <- ggplot(truth_prediction_df_spread,aes(x=mutationType_5mer,y=truth_oddsRatio_mmd_ms_rescaledMutationRate,color=outlierLabel,shape="truth"))+
  geom_point()+
  geom_point(aes(x=mutationType_5mer,y=modelPrediction_oddsRatio_mmd_ms_rescaledMutationRate,color=outlierLabel,shape="model prediction"))+
  scale_shape_manual(values=c(8,1))+
  ggtitle("SHAP outliers for population don't have extreme odds ratios")+
  facet_wrap(~centralMutationType)+
  geom_hline(yintercept = 1)+
  theme(axis.text.x = element_blank())
sandbox_predictionsOutliersPlot
ggsave(paste0(outdir,"sandbox.OutliersForPopulationSHap.NotOddsRatioOutliers.png"),sandbox_predictionsOutliersPlot,height=4,width=10)


########## look at original data of 5mers with large odds ratios of the rescaled mutation rates #####
## get shap values for the odds ratio outliers and see what's up
listOfOutlier5mers_OddsRatio <- subset(truth_prediction_df_spread,truth_oddsRatio_mmd_ms_rescaledMutationRate>1.5 | truth_oddsRatio_mmd_ms_rescaledMutationRate<0.75)$mutationType_5mer

# in truth dataset
truth_prediction_df[truth_prediction_df$mutationType_5mer %in% listOfOutlier5mers_OddsRatio,]

# do they have an overall weird odds ratio in seg sites?
# sum up over groups/labels
fracSegSites <- truth_prediction_df %>%
  group_by(population,newGroup) %>%
  mutate(fracOfSegSites = mutationCount_5mer/sum(mutationCount_5mer))
head(fracSegSites)
sum(fracSegSites[fracSegSites$population=="Ms" & truth_prediction_df$newGroup==1,]$fracOfSegSites) #  1 ;good

fracSegSites_spread <- pivot_wider(fracSegSites[,c("mutationType_5mer","population","fracOfSegSites")],id_cols = c(mutationType_5mer,population),names_from=population,values_from=c(fracOfSegSites)) 

head(fracSegSites_spread)
fracSegSites_spread$fracSegSites_oddsRatio_Mmd_ms <- fracSegSites_spread$Mmd/fracSegSites_spread$Ms
head(fracSegSites_spread)
# outliers:
listOfOutlier5mers_OddsRatio_FracSegSites <- subset(fracSegSites_spread,oddsRatio_Mmd_ms>1.5 | oddsRatio_Mmd_ms<0.75)$mutationType_5mer
listOfOutlier5mers_OddsRatio_FracSegSites

# let's compare odds ratios:
mergeOddsRatios <- merge(fracSegSites_spread[,c("mutationType_5mer","fracSegSites_oddsRatio_Mmd_ms")],truth_prediction_df_spread[,c("mutationType_5mer","truth_oddsRatio_mmd_ms_rescaledMutationRate")],by="mutationType_5mer")

comparingFracSegSites_MutationRateOddsRatios <- ggplot(mergeOddsRatios,aes(x=fracSegSites_oddsRatio_Mmd_ms,y=truth_oddsRatio_mmd_ms_rescaledMutationRate))+
  geom_point()+
  geom_abline()+
  geom_text(data=subset(mergeOddsRatios,(fracSegSites_oddsRatio_Mmd_ms-truth_oddsRatio_mmd_ms_rescaledMutationRate)>0.2),aes(label=mutationType_5mer))
comparingFracSegSites_MutationRateOddsRatios
ggsave(paste0(outdir,"sandbox.comparingFracsegSitestoMutationRateOddsRatios.FracSegSitesAlwaysBigger.png"),comparingFracSegSites_MutationRateOddsRatios,height=4,width=7)

# all the outliers are CGACG> ??? what's up with that.
truth_prediction_df[truth_prediction_df$ancestral5mer=="CGACG",c("mutationType_5mer","mutationCount_5mer","ancestral5merCount","mutationCount_divByTargetCount_5mer","population")]
# Okay *why* aren't then ancestral5mer counts the same HERE?
#sink()

######### want to try tree shap interactions ########
# https://www.r-bloggers.com/2021/01/treeshap%E2%80%8A-%E2%80%8Aexplain-tree-based-models-with-shap-values/
require(treeshap)
ranger_obj
model_unified <- ranger.unify(ranger_obj, select(Xdftrain,-c(outcome,mutationType_5mer))) # unify model with the training data 
# ranger.unify(ranger_obj,training data (must have same columns as training data so no outcome or type))
# Reference dataset. A data.frame or matrix with the same columns as in the training set of the model. Usually dataset used to train model.  https://rdrr.io/github/ModelOriented/treeshap/man/randomForest.unify.html
# note you can change ref dataset: unified2 <- set_reference_dataset(model_unified, aps_data[1:2000, ])
# compute shap values:

print(paste0(Sys.time()," staring treeshap with interactions"))
treeshap_res <- treeshap(model_unified, select(Xdftest,-c(outcome,mutationType_5mer,populationLabel)),interactions = T) ## Interactions TRUE! 

saveRDS(treeshap_res , paste0(outdir,"treeshapInteractions.rds"))
# head(treeshap_res$shaps)
# you get this error if a column is present that sholudn't be
#Error in treeshap_cpp(x2, is_na, roots, yes, no, missing, feature, split,  : 
#                        Not compatible with requested type: [type=character; target=double].

sink()
                              