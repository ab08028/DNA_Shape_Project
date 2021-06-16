################################## Trying a random forest model for mutation types ###############################
# single mutation type to start with (?)
# this is after going through the tidymodels tutorial and book

########## packages ###################
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
# Install the latest development version from GitHub:
#if (!requireNamespace("remotes")) {
#  install.packages("remotes")
#}
#remotes::install_github("bgreenwell/fastshap") # note: make sure ~/.R/Makevars file is empty 
require(fastshap)
require(reshape2)
Sys.setenv(RETICULATE_PYTHON="/opt/anaconda3/bin/python3.8")
require(reticulate) # for force plots
# to make reticulate python path work:
# cd ~
# nano .Renviron
# RETICULATE_PYTHON="enter your desired path here" 
# check with Sys.which("python") that it's right version
tidymodels_prefer() # use this 'tidymodels_prefer()' uses the 'conflicted' package to handle common conflicts with tidymodels and other packages. <-- should add this to my script
set.seed(42) # so results are reproducible 
outdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/sandbox/"

###### just one mutation type for now: #########
mutationType="A.T.no_ancestralCpG"
CVcount=8 # total folds for k-fold cross validation 
######## chosen species/population: ###########
chosenPop="Mmd"

###### read in shape data (works for all species; all possible 7mer motifs) : NOT normalized data! ###########
# (note that random forest doesn't need you to transform the data)
shapes <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/firstOrder_featureTypes_allPossible7mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.txt",header=T,sep="\t")
dim(shapes)
# need to assign row names:
rownames(shapes) <- shapes$motif


############ read in train and test data (NOTE MAY BE MULTIPLE POPS PRESENT IN DATA) ###################
########## train data:
train_spectra <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/training_test_datasets/mouse/mutyperResults_20210317_NOSTRICT_7mer/training_and_test_spectra_SUMMEDOVERCHROMOSOMES/mouse.TRAINING.odd_bpOnly.allPops.summedOverChr.spectra.NOSTRICT.txt",header=T)
# restrict to one pop and one mutation type (for now; maybe eventually will do more of both)
train_spectra <- train_spectra %>%
  filter(population == chosenPop & mutationClassLabel == mutationType) 
# this is a little finicky -- make sure it worked with this test:
if(length(unique(train_spectra$population))>1){
  stop("MORE THAN ONE SPECIES!")
}

######### test data:
test_spectra <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/training_test_datasets/mouse/mutyperResults_20210317_NOSTRICT_7mer/training_and_test_spectra_SUMMEDOVERCHROMOSOMES/mouse.TESTING.even_bpOnly.allPops.summedOverChr.spectra.NOSTRICT.txt",header=T)
test_spectra <- test_spectra %>%
  filter(population == chosenPop & mutationClassLabel == mutationType) 
if(length(unique(test_spectra$population))>1){
  stop("MORE THAN ONE SPECIES!")
}

##### check that first entry of test/train data is not equal #######
test_spectra[1,"fractionOfAllSegSites"] != train_spectra[1,"fractionOfAllSegSites"] # should be TRUE that these are not equal (may be *very* similar though)

######### merge shape data for ancestral and derived mutation motifs with the outcome values from the spectra #############
# get ancestral derived:
train_intermediate <- merge(train_spectra,shapes,by.x="ancestralkmer",by.y="motif")
# get the derived features (and set suffixes as ancestral and .derived):
train_data_unprocessed <- merge(train_intermediate,shapes,by.x="derivedkmer",by.y="motif",suffixes=c(".ancestral",".derived")) # sweet, I checked that this works and it does! 
# note that this data is UNPROCESSED (don't just use raw in a model-- needs to go through a recipe first)
rownames(train_data_unprocessed) <- train_data_unprocessed$variable # important so you can keep track of what motif is in what row without a column for it 

# note that the data includes columns that we don't want in the model (metadata)
# # get rid of all extra columns except variable and 

train_data_unprocessed <- train_data_unprocessed %>%
  select("fractionOfAllSegSites","variable",starts_with("feature_"))

######### process test data set #########
# get ancestral derived shape features:
test_intermediate <- merge(test_spectra,shapes,by.x="ancestralkmer",by.y="motif")
# get the derived features (and set suffixes as ancestral and .derived):
test_data_unprocessed <- merge(test_intermediate,shapes,by.x="derivedkmer",by.y="motif",suffixes=c(".ancestral",".derived")) # sweet, I checked that this works and it does! 
rownames(test_data_unprocessed) <- test_data_unprocessed$variable # important so you can keep track of what motif is in what row without a column for it 
# select relevant columns (exclude dataset)
test_data_unprocessed <- test_data_unprocessed %>%
  select("fractionOfAllSegSites","variable",starts_with("feature_")) 
####### check that dimensions are the same and check heads (make sure no extra populations are present, etc)
head(train_data_unprocessed)
head(test_data_unprocessed)


####################### set up recipe for data processing: random forest processing ###################
# see here for recommended processing for diff models: https://www.tmwr.org/pre-proc-table.html#pre-proc-table
# based on the table at that link, doing the non-zero variance removal and the decorrelation *may* help but also aren't necessary. Can compare with another recipe -- not doing for now.
# to get the basic code you need 
rand_forest_processing_recipe <- 
  # which consists of the formula (outcome ~ predictors) (don't want to include the 'variable' column)
  recipe(fractionOfAllSegSites ~ .,data=train_data_unprocessed) %>% # 
  update_role(variable, new_role="7mer mutation type label") # add a role for the variable (motif) label that isn't part of the model but keeps track of the label (so it was included in the "." above but now we make it not an outcome or predictor variable so it's not going to be part of the model but will be a label. you see that this works )
rand_forest_processing_recipe
# 
#Inputs:
  
#role #variables
#7mer mutation type label          1
#outcome          1
#predictor         96
# so this is the correct number ^^ 96 shape features and 1 outcome, and 1 label 
summary_rand_forest_processing_recipe <- summary(rand_forest_processing_recipe) # summary
predictorCountdf  <- summary_rand_forest_processing_recipe %>%
  filter(role=='predictor') %>%
  summarise(totalPredictors=n()) # count up total predictors
# be SUPER careful with *subset* if you forget = it will just push all through without throwing an error. 
predictorCountdf # should be 96 unless you removed correlated variables 
########### optional: could 'juice' the recipe to get out the processed data. not doign that here though (see my earlier scripts) because I don't need it #############

####### set up model specifications (not fitting the model yet) ###############
# you can use to get the right code for ranger RF with regression and to tune some parameters:
#==> to get model specific code: parsnip_addin() # gives you the basic code, and then you can add to it  
# to get the basic code you need ; going to tune two parameters 
rand_forest_ranger_model_specs <-
  rand_forest(trees = 1000, mtry = tune(), min_n = tune()) %>% # I added in tree number = 1000
  set_engine('ranger',importance="permutation",respect.unordered.factors="order",verbose=TRUE) %>%
  set_mode('regression')
# note about respect.unordered.factors : this is for categorical variables but is recommended by multiple tutorails
# https://win-vector.com/2016/05/30/on-ranger-respect-unordered-factors/ <-- doesn't really apply to my model, but adding it to be safe for the future so it get used if I copy code over.
rand_forest_ranger_model_specs

######## make a workflow with recipe and model specs ########
rand_forest_workflow <- workflow() %>%
  # add the recipe
  add_recipe(rand_forest_processing_recipe) %>%
  # add the model
  add_model(rand_forest_ranger_model_specs)
rand_forest_workflow


######### set up cross-validation folds from training dataset ##########
train_data_CVfolds_unprocessed <- vfold_cv(train_data_unprocessed,v = CVcount) # split into folds (unprocessed data because processing is going to occur during the workflow per-split)
train_data_CVfolds_unprocessed # has been split into 8 folds 

###### set up parameter tuning grid ##########
# want to get all combinations of parameters
# mtry: want to be ~ # features/3 for regression (sqrt(number of features) for classification) 
# one tutorial suggested even grid from 2 --> total number of parameters, with midpoint
mtrygoal=round(predictorCountdf$totalPredictors/3)
mtrygoal
mtrygrid=c(2,mtrygoal/2,mtrygoal,mtrygoal*2,mtrygoal*3) # range from 2-total features with mtrygoal (feat/3) in the middle
mtrygrid
min_ngrid = c(5,10,20) # not sure about these values but have to start somewhere. default is 5 for regression (minimum node size before gets converted to leaf)
rand_forest_tuning_grid <- expand.grid(mtry=mtrygrid,min_n=min_ngrid)

####### tune: SLOW ############

####### TUNE: **VERY ** slow! 
rf_tune_results <- rand_forest_workflow %>%
  tune_grid(resamples = train_data_CVfolds_unprocessed, #CV object
            grid = rand_forest_tuning_grid, # grid of values to try
            metrics = metric_set(rmse, rsq) # metrics we care about: RMSE and R^2
  )
# started at 3:45 on weds 5/26
rf_tune_results # finished 6 pm (or earlier)

# if you've already tuned you can just read it in and skip lower steps: 
#tune_metrics <- read.table(paste0(outdir,mutationType,".random_forest.ParameterTuning.txt",sep=""),header=T)

tune_metrics <- rf_tune_results %>%
  collect_metrics()
write.table(tune_metrics,paste0(outdir,mutationType,".random_forest.ParameterTuningMetrics.txt",sep=""),row.names = F,quote=F)
write.table(rf_tune_results,paste0(outdir,mutationType,".random_forest.ParameterTuningResults.txt",sep=""),row.names = F,quote=F)


 # see which metric gives best performance and set that value
# We can extract the best value for the accuracy metric by applying the select_best() function to the tune object.

param_final_rsq <- rf_tune_results %>%
  select_best(metric = "rsq")
param_final_rsq   

# ok so defaults are the best ! maybe earch more? 

param_final_rmse <- rf_tune_results %>%
  select_best(metric = "rmse")
param_final_rmse   
# hhmm are pretty different for rmse! 16 and 10... 

metricTuningPlot <- ggplot(tune_metrics,aes(x=mtry,y=mean,color=as.factor(min_n)))+
  facet_wrap(~.metric,scales="free")+
  geom_line()+
  geom_point()+
  theme_bw()+
  ggtitle(paste0("note values are averaged over all folds (some folds much better than others)\n"))
metricTuningPlot
ggsave(paste0(outdir,mutationType,".random_forest.ParameterTuning.plot.png"),metricTuningPlot,height=5,width=9)
# get predictions per fold:
tune_metrics_perFold <- collect_metrics(rf_tune_results,summarize = F) # to get predictions per fold with summarise=F

# interesting that rsq so variable between folds -- plot that as well! Not sure how to extract right now but can look in tutorial R script  
metricTuningPlot_perfold <- ggplot(tune_metrics_perFold,aes(x=mtry,y=.estimate,color=as.factor(id),group=interaction(min_n,id),shape=as.factor(min_n)))+
  facet_wrap(~.metric,scales="free")+
  geom_line()+
  geom_point()+
  theme_bw()+
  ggtitle(paste0("different folds have very different rsq/rmse values"))
metricTuningPlot_perfold
ggsave(paste0(outdir,mutationType,".random_forest.ParameterTuning.PERFOLD.plot.png"),metricTuningPlot_perfold,height=5,width=9)

  
#Optimal parameters based on rsq: mtry = 32, min_n = 5 (the defaults)
#Optimal parameters based on rmse: mtry = 16, min_n =10
# But functionally very similar results (mean rmse with 16,10 =1.896978e-5 and mean rmse with 32,5= 1.900329e-5)
#-- so going with the defaults for now: param_final_rsq

param_final <- param_final_rsq

######## FINALIZE workflow: add tuned parameters to workflow ############
rand_forest_workflow <- rand_forest_workflow %>%
  finalize_workflow(param_final)
rand_forest_workflow

# if you are setting best tuned aprameters manually you can do : 
#paramdf = tibble(mtry=32,min_n=5)
#rand_forest_workflow <- rand_forest_workflow %>%
#  finalize_workflow(paramdf)
#rand_forest_workflow

############# make a train/test split object for training the model and testing it ############
test_data_unprocessed$label <- "TEST"
train_data_unprocessed$label <- "TRAIN"
ALLDATA_dontusetotrain <- bind_rows(train_data_unprocessed,test_data_unprocessed)

# make an indicator of what is what: have to be labeled analysis (training) and assessment (test)
indices <-
  list(analysis   = which(ALLDATA_dontusetotrain$label=="TRAIN"), 
       assessment = which(ALLDATA_dontusetotrain$label=="TEST"))

split <- make_splits(indices,select(ALLDATA_dontusetotrain,-label)) # get rid of label column ; this separates train and test by index
# if you want to see what's what:
head(training(split),4)
head(testing(split),4)
# NOTE THAT THESE ARE UNPROCESSED 
# ^ how you get test/train data out of split



############ fit the model using the split (train/test) object #########
#Now we’ve defined our recipe, our model, and tuned the model’s parameters, we’re ready to actually fit the final model. Since all of this information is contained within the workflow object, we will apply the last_fit() function to our workflow and our train/test split object. This will automatically train the model specified by the workflow using the training data, and produce evaluations based on the test set.

# or training on the training data and fitting to the test data : (better)
rand_forest_fit <- rand_forest_workflow %>%
  # fit on the training set and evaluate on test set
  last_fit(split) # this is the split object we created above. it has the training and test spectra (they were bound together and then split up by index into a custom split object. represents spectrum measured from even/odd bp)

rand_forest_performanceOnTest <- rand_forest_fit %>% collect_metrics()
rand_forest_performanceOnTest
write.table(rand_forest_performanceOnTest,paste0(outdir,mutationType,".random_forest.PerformanceOnTestDataset.txt",sep=""),row.names = F,quote=F)

# title for plot:
rand_forest_performanceTitle=paste0("rsq: ",signif(filter(rand_forest_performanceOnTest,.metric=="rsq")$.estimate,digits=3),"\n","rmse: ",signif(filter(rand_forest_performanceOnTest,.metric=="rmse")$.estimate,digits=3))
rand_forest_performanceTitle

# pull out predictions of test data:
rand_forest_predictionsOnTest <- rand_forest_fit %>% collect_predictions()
rand_forest_predictionsOnTest

rand_forest_predictionsOnTest$residual_abs <- abs(rand_forest_predictionsOnTest$fractionOfAllSegSites - rand_forest_predictionsOnTest$.pred)

# cbind with original testing data to get mutation type labels:
######### IF you did preprocessing this might not work well..... figure out a better way to label!!!! 
rand_forest_predictionsOnTest_labeled <- cbind(rand_forest_predictionsOnTest,select(testing(split),variable)) # 

write.table(rand_forest_predictionsOnTest_labeled,paste0(outdir,mutationType,".random_forest.PredictionsOfTestDataset.txt",sep=""),row.names = F,quote=F)

# % of residual of total?
rand_forest_predictionsOnTest_labeled$residual_abs_FractionOfObserved <- abs(rand_forest_predictionsOnTest_labeled$fractionOfAllSegSites - rand_forest_predictionsOnTest$.pred) / rand_forest_predictionsOnTest_labeled$fractionOfAllSegSites
head(rand_forest_predictionsOnTest_labeled)

fitToTestData <- ggplot(rand_forest_predictionsOnTest_labeled,aes(x=fractionOfAllSegSites,y=.pred))+
  geom_point(alpha=0.6)+
  geom_abline(slope=1,intercept=0)+
  theme_bw()+
  geom_text_repel(data=filter(rand_forest_predictionsOnTest_labeled,rand_forest_predictionsOnTest_labeled$residual_abs>0.00005),aes(label=variable))+
  ggtitle(rand_forest_performanceTitle)
fitToTestData
ggsave(paste0(outdir,mutationType,".random_forest.FitToTestData.plot.png"),fitToTestData,height=5,width=9)


fitToTestData_log10 <- ggplot(rand_forest_predictionsOnTest_labeled,aes(x=fractionOfAllSegSites,y=.pred))+
  geom_point(alpha=0.6)+
  geom_abline(slope=1,intercept=0)+
  theme_bw()+
  geom_text_repel(data=filter(rand_forest_predictionsOnTest_labeled,rand_forest_predictionsOnTest_labeled$residual_abs>0.00005),aes(label=variable))+
  ggtitle(rand_forest_performanceTitle)+
  geom_text_repel(data=filter(rand_forest_predictionsOnTest_labeled,residual_abs_FractionOfObserved>0.3),aes(label=variable))+
  scale_x_log10()+
  scale_y_log10()
fitToTestData_log10
ggsave(paste0(outdir,mutationType,".random_forest.FitToTestData.logscale.plot.png"),fitToTestData_log10,height=5,width=9)


####### explore results using the .workflow section ########
workflow_fit_info <- rand_forest_fit %>% 
  pluck(".workflow", 1) %>%
  pull_workflow_fit()
# contains all the ranger info and oob stats:
ranger_object <- workflow_fit_info$fit
# to get OOB rsq:
oobrsq <- ranger_object$r.squared
oobrsq
#### how do I make this a nice table?
######### vip:  variable importance #############

vipPlot <- rand_forest_fit %>% 
  pluck(".workflow", 1) %>%   
  pull_workflow_fit() %>% 
  vip(num_features=96,include_type=T)
vipPlot
ggsave(paste0(outdir,mutationType,".random_forest.FitToTestData.logscale.plot.png"),vipPlot,height=10,width=6)


######### try to get SHAP values to work ########
# https://www.hfshr.xyz/posts/2020-06-07-variable-importance-with-fastshap/
#### get the processed train data: ####### 
# Apply the preprocessing steps with prep and juice to the training data to make sure it's processed the same way it was in the model fitting: 
 # try just fitting to train dataset (not doing test predictions)
rf_wflow <- workflow() %>%
  add_recipe(rand_forest_processing_recipe) %>%
  add_model(rand_forest_ranger_model_specs) %>% 
  finalize_workflow(param_final) %>% # need to add tune params # can put in paramdf if no param_final is present 
  fit(training(split))
rf_wflow

Xtrain <- prep(rand_forest_processing_recipe, training(split)) %>% 
  juice() %>% 
  select(-c(variable,fractionOfAllSegSites)) %>% 
  as.matrix()
head(Xtrain)


Xtest <- prep(rand_forest_processing_recipe, testing(split)) %>% 
  juice() %>% 
  select(-c(variable,fractionOfAllSegSites)) %>% 
  as.matrix()
head(Xtest)

pfun <- function(object, newdata) {
  predict(object, data = newdata)$predictions
}


# a bit slow: 

shap <- fastshap::explain(ranger_object, X = Xtrain,pred_wrapper=pfun,newdata=Xtest) ##
                          #eventually want: nsim=5,newdata=Xtest,adjust=T)
# oooh is it working? try to increase nsim (number of MC simulations ); adjust ot satisfy additivity principle
# VERY SLOW
shap

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
ggsave(paste0(outdir,mutationType,".random_forest.shapImportance.png"),shapplot1,height=10,width=6)
######## shapley dependence plot ########
# need dataframe of training data:
# note X above was a matrix but this has to be a df
Xdftrain <- prep(rand_forest_processing_recipe, training(split)) %>% 
  juice()  #%>%
  #select(-c(variable,fractionOfAllSegSites))

head(Xdftrain)

# note X above was a matrix but this has to be a df
Xdftest <- prep(rand_forest_processing_recipe, testing(split)) %>% 
  juice()  #%>%
#select(-c(variable,fractionOfAllSegSites))

head(Xdftest)



featuresToPlot=names(shap)
dir.create(paste0(outdir,"perFeatureDependencePlots/"),showWarnings = F)
for(feature in featuresToPlot){

shapdependenceplot1 <- autoplot(shap, 
         type = "dependence", 
         feature = feature, 
         X = Xdftest, # this hsould be test data if you ran shap with newdata=testing data! sholud be train if you just ran it without any new data
         smooth = TRUE, color_by = "fractionOfAllSegSites")+
  scale_color_viridis_c(trans="log10")+
  ggtitle(feature)+
  theme_bw()+
  labs(color="fractionOfAllSegSites (log10)")
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



# requires python and shap to be correctly pathed:
forceplot1 <- force_plot(object = shap[rowNumberForForce,], 
           feature_values = select(Xdftest[rowNumberForForce,],-c(variable,fractionOfAllSegSites)), 
           display = "html")
write.table(forceplot1,paste0(outdir,mutationType,".random_forest.shapForcePlot.",motifOfInterest,".html"),row.names=F,quote=F,col.names=F)

# multiple obs at once:
forceplot2 <- force_plot(object = shap, 
           feature_values = select(Xdftest,-c(variable,fractionOfAllSegSites)), 
           display = "html") # , 
           #link = "logit") 
write.table(forceplot2,paste0(outdir,mutationType,".random_forest.shapForcePlot.AllSamples.html"),row.names=F,quote=F,col.names=F)
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


