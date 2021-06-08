###### want to try xgboost ##########
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

parsnip_addin()
####################### set up recipe for data processing: xgboost processing ###################
# see here for recommended processing for diff models: https://www.tmwr.org/pre-proc-table.html#pre-proc-table
# based on the table at that link, doing the non-zero variance removal and the decorrelation *may* help but also aren't necessary. Can compare with another recipe -- not doing for now.

#xgboost recommendations : https://www.r-bloggers.com/2020/05/using-xgboost-with-tidymodels/
  
  
# to get the basic code you need 
xgboost_processing_recipe <- 
  # which consists of the formula (outcome ~ predictors) (don't want to include the 'variable' column)
  recipe(fractionOfAllSegSites ~ .,data=train_data_unprocessed) %>% 
  step_rm(label) %>%
  update_role(variable, new_role="7mer mutation type label")  %>%
  step_dummy(all_nominal_predictors())  # don't currently have any categorial vars but need to transform if I get some
xgboost_processing_recipe
### need to recalculate things actually based on the totalSitesAllIntervals rather than frac of seg sites - because should be recalc'd for each fold ###
  
  # add a role for the variable (motif) label that isn't part of the model but keeps track of the label (so it was included in the "." above but now we make it not an outcome or predictor variable so it's not going to be part of the model but will be a label. you see that this works )

boost_tree_xgboost_spec <-
  boost_tree(tree_depth = tune(), trees = 1000, learn_rate = tune(), min_n = tune(), loss_reduction = tune()) %>%
  set_engine('xgboost') %>%
  set_mode('regression')
boost_tree_xgboost_spec
# other things you could tune: sample_size = tune()

######## make a workflow with recipe and model specs ########
xgboost_workflow <- workflow() %>%
  # add the recipe
  add_recipe(xgboost_processing_recipe) %>%
  # add the model
  add_model(boost_tree_xgboost_spec)
xgboost_workflow

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





######### set up cross-validation folds from training dataset: NEED TO FIX THIS FOR GENOMIC REGIONS NOT 7mers#########
train_data_CVfolds_unprocessed <- vfold_cv(train_data_unprocessed,v = CVcount) # split into folds (unprocessed data because processing is going to occur during the workflow per-split)
train_data_CVfolds_unprocessed # has been split into 8 folds 

###### tune parameters (will be slow) #######
xgboost_tune_results <- xgboost_workflow %>%
  tune_grid(resamples = train_data_CVfolds_unprocessed, #CV object
            grid = xgboost_grid, # grid of values to try
            metrics = metric_set(rmse, rsq) # metrics we care about: RMSE and R^2
  )
# started at 1:41 on 6/3
xgboost_tune_results # finished 6 pm (or earlier)
tune_metrics <- xgboost_tune_results %>%
  collect_metrics()
tune_metrics
