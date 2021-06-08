require(tidymodels)
require(tidyverse) # instead of caret going to use tidymodels
require(workflows)
require(tune)
require(vip) # for importance and shap values
require(ranger) # for random forest model
require(glmnet) # for lasso and linear modeling (maybe)
require(ggplot2)
require(ggrepel)
require(fastshap)
tidymodels_prefer() # use this 'tidymodels_prefer()' uses the 'conflicted' package to handle common conflicts with tidymodels and other packages. <-- should add this to my script
set.seed(42) #so results are reproducible 
######### try for just one mutation type
mutationType="A.T.no_ancestralCpG"
################## modeling of mouse spectrum ##############
correlationThreshold=0.75
CVcount=8

############ >>> read in train and test data (NOTE MAY BE MULTIPLE POPS PRESENT) ###################
traindata <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/training_test_datasets/mouse/mutyperResults_20210317_NOSTRICT_7mer/training_and_test_spectra_SUMMEDOVERCHROMOSOMES/mouse.TRAINING.odd_bpOnly.allPops.summedOverChr.spectra.NOSTRICT.txt",header=T)

testdata <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/training_test_datasets/mouse/mutyperResults_20210317_NOSTRICT_7mer/training_and_test_spectra_SUMMEDOVERCHROMOSOMES/mouse.TESTING.even_bpOnly.allPops.summedOverChr.spectra.NOSTRICT.txt",header=T)

######## >>> restrict to a single species/population! Important! ######################
population="Mmd"
traindata_onepop <- traindata[traindata$population==population & traindata$mutationClassLabel==mutationType,]
testdata_onepop <- testdata[testdata$population==population & testdata$mutationClassLabel==mutationType,]

######### >>> read in shape table for all possible 7mers ############
# using normalized at least for now:
shapes <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/firstOrder_featureTypes_allPossible7mers.FeatureValues.min-maxNormalized.WillWorkForAnySpecies.usethis.txt",header=T,sep="\t")
dim(shapes)
# need to assign row names:
rownames(shapes) <- shapes$motif
######### CONSIDER: reading in non-normalized version for making CV object 
######## >>> make training feature matrix by combining the shape features with the training data by ancestral and derived motifs ########
# get ancestral derived:
traindf_intermediateMerge_ancestral <- merge(traindata_onepop,shapes,by.x="ancestralkmer",by.y="motif")
# get the derived features (and set suffixes as ancestral and .derived):
traindf <- merge(traindf_intermediateMerge_ancestral,shapes,by.x="derivedkmer",by.y="motif",suffixes=c(".ancestral",".derived")) # sweet, I checked that this works and it does! 
rownames(traindf) <- traindf$variable # important so you can keep track of what motif is in what row without a column for it 

head(traindf)

# note here that fractionOfAllSegSites will be the dependent Y variable. and you need to exclude all the other columns that aren't shape features 
columnsToExclude = names(select(traindata_onepop,-fractionOfAllSegSites))
columnsToExclude
# exclude columns that aren't shape features:
traindf <- select(traindf,-all_of(columnsToExclude))
head(traindf)


######### process test data set #########
# get ancestral derived:
testdf_intermediateMerge_ancestral <- merge(testdata_onepop,shapes,by.x="ancestralkmer",by.y="motif")
# get the derived features (and set suffixes as ancestral and .derived):
testdf <- merge(testdf_intermediateMerge_ancestral,shapes,by.x="derivedkmer",by.y="motif",suffixes=c(".ancestral",".derived")) # sweet, I checked that this works and it does! 
rownames(testdf) <- testdf$variable # important so you can keep track of what motif is in what row without a column for it 

head(testdf)

# note here that fractionOfAllSegSites will be the dependent Y variable. and you need to exclude all the other columns that aren't shape features 
columnsToExclude = names(select(testdata_onepop,-fractionOfAllSegSites))
columnsToExclude
# exclude columns that aren't shape features:
testdf <- select(testdf,-all_of(columnsToExclude))
head(testdf)



####### >>> SANDBOX. try some modeling with tidy models ########

# THIS IS JUST AN EXPLORATION
# want to make sure i've filtered for nvz ; intercorrelation
# following and quoting from: http://www.rebeccabarter.com/blog/2020-03-25_machine_learning/#getting-set-up


######## >>> set up a recipe ##########
# list of preprocessing steps: https://www.tidymodels.org/find/recipes/
# "Recipes allow you to specify the role of each variable as an outcome or predictor variable (using a “formula”), and any pre-processing steps you want to conduct (such as normalization, imputation, PCA, etc).
# Creating a recipe has two parts (layered on top of one another using pipes %>%):"
#  1. Specify the formula (recipe()): specify the outcome variable and predictor variables
# 2. specify pre-processing steps (step_zzz()): define the pre-processing steps, such as imputation, creating dummy variables, scaling, and more
# which consists of the formula (outcome ~ predictors)
my_recipe <- 
  # which consists of the formula (outcome ~ predictors)
  recipe(fractionOfAllSegSites~.,data=traindf) %>%
  # and some pre-processing steps
  # check for near zero variance 
  # specify all_numeric() to specify that's what I want nzv and corr to be run on 
  step_nzv(all_numeric_predictors()) %>%
  step_corr(all_numeric_predictors(),threshold =correlationThreshold)
my_recipe
summary(my_recipe)
# this hasn't done any modeling yet! is just setting up the recipe!

######## >>> optional step: juice your preprocessed data: YOU DON'T NEED TO DO THIS; THIS IS JUST TO LOOK AT YOUR PREPROCESSED DATA. ######## 
# "If you want to extract the pre-processed dataset itself, you can first prep() the recipe for a specific dataset and juice() the prepped recipe to extract the pre-processed data. It turns out that extracting the pre-processed data isn’t actually necessary for the pipeline, since this will be done under the hood when the model is fit, but sometimes it’s useful anyway." ###
traindf_preprocessed <- my_recipe %>%
  # apply the recipe to the training data
  prep(traindf) %>%
  # extract the pre-processed training dataset
  juice()

traindf_preprocessed
dim(traindf)
dim(traindf_preprocessed) # some variables have been removed by the above 

cat(toString(length(setdiff(names(traindf),names(traindf_preprocessed)))), " variables removed during preprocessing:\n", toString(setdiff(names(traindf),names(traindf_preprocessed))))
# note in other tutorials you might see prep/bake steps but don't need those here

###### >>> specify the model ####
# https://www.tidymodels.org/find/parsnip/
# linear model (has lasso options) https://parsnip.tidymodels.org/reference/linear_reg.html
# random forest: https://parsnip.tidymodels.org/reference/rand_forest.html
# using the parsnip package
# only have to learn one way of specifying model and will work for lots of models


#The model type: what kind of model you want to fit, set using a different function depending on the model, such as rand_forest() for random forest, logistic_reg() for logistic regression, svm_poly() for a polynomial SVM model etc. The full list of models available via parsnip can be found here.

#The arguments: the model parameter values (now consistently named across different models), set using set_args().

#The engine: the underlying package the model should come from (e.g. “ranger” for the ranger implementation of Random Forest), set using set_engine().

#The mode: the type of prediction - since several packages can do both classification (binary/categorical prediction) and regression (continuous prediction), set using set_mode().
# Another thing to note is that nothing about this model specification is specific to the dataset.

# EXAMPLE:
rf_model <- 
  # specify that the model is a random forest
  rand_forest() %>%
  # specify that the `mtry` parameter needs to be tuned
  set_args(mtry = tune()) %>%
  # select the engine/package that underlies the model
  set_engine("ranger", importance = "impurity") %>%
  # choose either the continuous regression or binary classification mode
  set_mode("regression")

# Note that this code doesn’t actually fit the model. Like the recipe, it just outlines a description of the model. Moreover, setting a parameter to tune() means that it will be tuned later in the tune stage of the pipeline (i.e. the value of the parameter that yields the best performance will be chosen). You could also just specify a particular value of the parameter if you don’t want to tune it e.g. using set_args(mtry = 4).

# what is impurity doing? need to read up on random forests
# If you want to be able to examine the variable importance of your final model later, you will need to set importance argument when setting the engine. For ranger, the importance options are "impurity" or "permutation".



######### >>> put it together in a workflow ##############
# make a workflow
# set the workflow
rf_workflow <- workflow() %>%
  # add the recipe
  add_recipe(my_recipe) %>%
  # add the model
  add_model(rf_model)
rf_workflow
# Note that we still haven’t yet implemented the pre-processing steps in the recipe nor have we fit the model. We’ve just written the framework. It is only when we tune the parameters or fit the model that the recipe and model frameworks are actually implemented.


############ >>> tune parameters (make CV object and tune parameters (can tune multiple at once)) ############
# Since we had a parameter that we designated to be tuned (mtry), we need to tune it (i.e. choose the value that leads to the best performance) before fitting our model. If you don’t have any parameters to tune, you can skip this step.


# Note that we will do our tuning using the cross-validation object (diabetes_cv). To do this, we specify the range of mtry values we want to try, and then we add a tuning layer to our workflow using tune_grid() (from the tune package). Note that we focus on two metrics: accuracy and roc_auc (from the yardstick package).

# make a cross validated object for tuning (NOTE: this was normalized pre CV splitting which may not be ideal; in future maybe do this on the non-noramlized version)
traindf_CV <- vfold_cv(traindf,v = CVcount) # split into folds 
traindf_CV


#To do this, we specify the range of mtry values we want to try, and then we add a tuning layer to our workflow using tune_grid() (from the tune package). Note that we focus on two metrics: accuracy and roc_auc (from the yardstick package).
# specify which values eant to try
rf_grid <- expand.grid(mtry = c(3, 4, 5))
rf_grid
#   mtry
#1    3
#2    4
#3    5

####### TUNE: slow! 
rf_tune_results <- rf_workflow %>%
  tune_grid(resamples = traindf_CV, #CV object
            grid = rf_grid, # grid of values to try
            metrics = metric_set(rmse, rsq) # metrics we care about: RMSE and R^2
  )

# You can tune multiple parameters at once by providing multiple parameters to the expand.grid() function, e.g. expand.grid(mtry = c(3, 4, 5), trees = c(100, 500)).

# want to try this with the lasso parameters too


#It’s always a good idea to explore the results of the cross-validation. collect_metrics() is a really handy function that can be used in a variety of circumstances to extract any metrics that have been calculated within the object it’s being used on. In this case, the metrics come from the cross-validation performance across the different values of the parameters.
rf_tune_results %>%
  collect_metrics()

# see which metric gives best performance and set that value
# We can extract the best value for the accuracy metric by applying the select_best() function to the tune object.

param_final <- rf_tune_results %>%
  select_best(metric = "rsq")
param_final
# in this case best value is 5 for mtry 
######## >>> add tuned parameter to workflow ############
rf_workflow <- rf_workflow %>%
  finalize_workflow(param_final)

########## >>> make a split object for fitting the model ###########
# they expected to do the splitting for you, so when you've already split your data it's a bit annoying:
# need to do what's in here: https://github.com/tidymodels/rsample/issues/158
# put them together and then you can split based on row indicies into a split object (annoying)
testdf$label <- "TEST"
traindf$label <- "TRAIN"
ALLDATA_dontusetotrain <- bind_rows(traindf,testdf)

# make an indicator of what is what: have to be labeled analysis (training) and assessment (test)
indices <-
  list(analysis   = which(ALLDATA_dontusetotrain$label=="TRAIN"), 
       assessment = which(ALLDATA_dontusetotrain$label=="TEST"))

split <- make_splits(indices,select(ALLDATA_dontusetotrain,-label)) # get rid of label column ; this separates train and test by index
# and makes the object that tidy models wants for model training and testing
# optional at this stage: get rid of label column:
# traindf <- select(traindf,-label)


############ >>> fit the model using the split (train/test) object #########
#Now we’ve defined our recipe, our model, and tuned the model’s parameters, we’re ready to actually fit the final model. Since all of this information is contained within the workflow object, we will apply the last_fit() function to our workflow and our train/test split object. This will automatically train the model specified by the workflow using the training data, and produce evaluations based on the test set.

rf_fit <- rf_workflow %>%
  # fit on the training set and evaluate on test set
  last_fit(split) # this is the split object we created above. it has the training and test spectra (they were bound together and then split up by index into a custom split object. represents spectrum measured from even/odd bp)
### SPEND TIME HERE MAKING SURE YOU GET EXACTLY WHAT'S GOING ON
rf_fit # Note that the fit object that is created is a data-frame-like object; specifically, it is a tibble with list columns.
# This is a really nice feature of tidymodels (and is what makes it work so nicely with the tidyverse) since you can do all of your tidyverse operations to the model object. While truly taking advantage of this flexibility requires proficiency with purrr, if you don’t want to deal with purrr and list-columns, there are functions that can extract the relevant information from the fit object that remove the need for purrr as we will see below.
rf_fit$.workflow
# Since we supplied the train/test object when we fit the workflow, the metrics are evaluated on the test set. Now when we use the collect_metrics() function (recall we used this when tuning our parameters), it extracts the performance of the final model (since rf_fit now consists of a single final model) applied to the test set.
test_performance <- rf_fit %>% collect_metrics()
test_performance
# jesus christ!
# You can also extract the test set predictions themselves using the collect_predictions() function. Note that there are 192 rows in the predictions object below which matches the number of test set observations (just to give you some evidence that these are based on the test set rather than the training set).

test_predictions <- rf_fit %>% collect_predictions()
test_predictions
test_predictions$motif <- rownames(testdf)

# how the hell does this fit so well!?!? is there some leakage somewhere? 
# what's actually being predicted? trying to predict fraction of seg sites based on shape features.
# model gets trained on odd bp and tested on even bp to predict their fraction of mutation rate
# need to make sure it's behaving as expected
####### make some assessment plots ########
test_predictions$abs_difference <- abs(test_predictions$fractionOfAllSegSites - test_predictions$.pred)
ggplot(test_predictions,aes(x=.pred,y=fractionOfAllSegSites))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  geom_label_repel(data=subset(test_predictions,test_predictions$abs_difference>0.00005),aes(label=motif))

##### >>> WHAT DOES THIS DO WHEN DONE CORRECTLY doing for now: can use model to run on new data ########
# final_model <- fit(rf_workflow, diabetes_clean) <- could fit on entire dataset 
# then predict based on new data: predict(final_model, new_data = new_woman)
final_model <- fit(rf_workflow,select(ALLDATA_dontusetotrain,-label)) ### NOTE SURE THIS STEP IS RIGHT. ??  DOING IT FOR NOW  ??? maybe use whole dataset? What does this do 
# they wanted met o do this on whole dataset? why
final_model # OOB Rsquared is better -- why is that different than Rsq above?? makes more sense thoguh. saw something about that online
#### NEED TO READ A LOT HERE ... ok so if you use 'whole' dataset the OOB is good, but just test was not good. what's going on??

# don't trust this Rsq! need the out of bag Rsq (maybe?)
# confused at this stage 
##### get variable importance ########
ranger_obj <- pull_workflow_fit(final_model)$fit
ranger_obj
ranger_obj$variable.importance # what are units here?

######### SHAP values ################
shap()
# https://koalaverse.github.io/vip/reference/vi_shap.html
# run vip and vip_shap(?) to get shapley values and importance 
vip(ranger_obj,include_type = T) # what is importance here?? not sure what impurity meansneed to read up
vi_shap(ranger_obj) # need fast shap 


########## >>> try to predict a different mouse species using final model ##############
otherSpecies="Ms"
# other test data : 
testdata_pop2 <- testdata[testdata$population==otherSpecies & testdata$mutationClassLabel==mutationType,]
testdata_pop2_intermediateMerge_ancestral <- merge(testdata_pop2,shapes,by.x="ancestralkmer",by.y="motif")
# get the derived features (and set suffixes as ancestral and .derived):
testdf_pop2 <- merge(testdata_pop2_intermediateMerge_ancestral,shapes,by.x="derivedkmer",by.y="motif",suffixes=c(".ancestral",".derived")) # sweet, I checked that this works and it does! 
testdf_pop2 <- select(testdf_pop2,-all_of(columnsToExclude))
dim(testdf_pop2)
names(testdf_pop2)
pop2prediction <- predict(final_model,new_data = testdf_pop2)

pop2prediction$fractionOfAllSegSites <- testdf_pop2$fractionOfAllSegSites
pop2prediction$abs_difference <- abs(pop2prediction$fractionOfAllSegSites - pop2prediction$.pred)
pop2prediction$motif <- testdata_pop2$variable # is this correct? does it get resorted at any point? 

ggplot(pop2prediction,aes(x=.pred,y=fractionOfAllSegSites))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  geom_label_repel(data=subset(pop2prediction,pop2prediction$abs_difference>0.00005),aes(label=motif))+
  ggtitle(paste0("model trained on ",population," predicting values in ",otherSpecies))
# get rmse too ; maybe more useful than rsq?

