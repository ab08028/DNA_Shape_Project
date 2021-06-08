######### Following along to tutorial in this book
# https://www.tmwr.org/ames.html
require(modeldata)
require(tidymodels)
data(ames) # using the ames data
head(ames)
dim(ames)
require(purrr)
require(usemodels)
# to use add in that lets you browse models:
require(shiny)
require(miniUI)
rand_forest_randomForest_spec <-
  rand_forest(mtry = tune(), min_n = tune()) %>%
  set_engine('randomForest') %>%
  set_mode('regression')

tidymodels_prefer() # use this 'tidymodels_prefer()' uses the 'conflicted' package to handle common conflicts with tidymodels and other packages. <-- should add this to my script
ggplot(ames, aes(x = Sale_Price)) + 
  geom_histogram(bins = 50)

# log transform 
ggplot(ames, aes(x = Sale_Price)) + 
  geom_histogram(bins = 50) +
  scale_x_log10()

ggplot(ames, aes(x = Sale_Price)) + 
  geom_histogram(bins = 50) +
  scale_x_log10()
# While not perfect, this will probably result in better models than using the untransformed data.

#From this point on, the outcome column is pre-logged in the ames data frame:
  
ames <- ames %>% mutate(Sale_Price = log10(Sale_Price))

######### splitting the data ##############
# This is held in reserve until one or two models are chosen as the methods that are most likely to succeed. The test set is then used as the final arbiter to determine the efficacy of the model. It is critical to only look at the test set once; otherwise, it becomes part of the modeling process.
# Set the random number stream using `set.seed()` so that the results can be 
# reproduced later. 
set.seed(123)

# Save the split information for an 80/20 split of the data
ames_split <- initial_split(ames, prop = 0.80)
ames_split
# The printed information denotes the amount of data in the training set (n=2,344), the amount in the test set (n=586), and the size of the original pool of samples (n=2,930).
# The object ames_split is an rsplit object and only contains the partitioning information; to get the resulting data sets, we apply two more functions:
ames_train <- training(ames_split)
ames_test  <-  testing(ames_split)
# okay so with my data I could split them back up from the split to make sure it's right.

dim(ames_train)
#> [1] 2344   74

# These objects are data frames with the same columns as the original data but only the appropriate rows for each set.
# may be better to split within quartiles so data isn't skewed
# I think my even/odd bp split is fine
# can use strata to stratify within data subsets by sale price
set.seed(123)
ames_split <- initial_split(ames, prop = 0.80, strata = Sale_Price)
ames_train <- training(ames_split)
ames_test  <-  testing(ames_split)

dim(ames_train)
#> [1] 2342   74
# whole different set of things to do for time series


# Throughout this book, notice which data are exposed to the model at any given time. Remember that it is critical to quarantine the test set from any model building activities.

# Keeping the training data in a separate data frame from the test set is one small check to make sure that information leakage does not occur by accident.
# but then why do we have to use a split object which seems fishy?

############ chapter 6 : recipes ###########################
# combine different feature engingeering and preproc steps to apply to different datasets (nice!)
# focusing on subset of ames data: neighborhood (qual), living area (cont), year built, type of building (qual)
# normal lm:
lm(Sale_Price ~ Neighborhood + log10(Gr_Liv_Area) + Year_Built + Bldg_Type, data = ames)
# predicting sales price 
#When this function is executed, the data are converted from a data frame to a numeric design matrix (also called a model matrix) and #then the least squares method is used to estimate parameters. 
library(tidymodels) # Includes the recipes package
tidymodels_prefer()
# put the train data into the recipe:
simple_ames <- 
  recipe(Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type,
         data = ames_train) %>%
  step_log(Gr_Liv_Area, base = 10) %>% 
  step_dummy(all_nominal_predictors())
simple_ames
#> Data Recipe
#> 
#> Inputs:
#> 
#>       role #variables
#>    outcome          1
#>  predictor          4
#> 
#> Operations:
#> 
#> Log transformation on Gr_Liv_Area
#> Dummy variables from all_nominal_predictors()

# Let’s break this down:
#   
#   The call to recipe() with a formula tells the recipe the roles of the variables (e.g., predictor, outcome). It only uses the data to determine the data types for the columns.
# 
# step_log() declares that Gr_Liv_Area should be log transformed.
# 
# step_dummy() is used to specify which variables should be converted from a qualitative format to a quantitative format, in this case, using dummy or indicator variables. An indicator or dummy variable is a binary numeric variable (a column of ones and zeroes) that encodes qualitative information; we will dig deeper into these kinds of variables in Section 6.3.
# 
# The function all_nominal_predictors() captures the names of any predictor columns that are currently factor or character (i.e., nominal) in nature. This is a dplyr selector function similar to starts_with() or matches() but can only be used inside of a recipe.
# 
# Other selectors specific to the recipes package are: all_numeric_predictors(), all_numeric(), all_predictors(), and all_outcomes(). As with dplyr, one or more unquoted expressions, separated by commas, can be used to select which columns are affected by each step.
# 
# What is the advantage to using a recipe? There are a few, including:
#   

# cool: could be used across any model that has these same outcome/predictors and same pre-proc
# compact 
# all data processing is in a single script 

# key thing is just setting up a recipe doesn't actually run or change anything about your training dataset
# have to use prep:
# i don't thikn you need to do this with a workflow but if you want to do stuff outside of workflow you need to use prep()
simple_ames <- prep(simple_ames, training = ames_train) # this carres out our recipe on the train dataset
simple_ames # prep calculates stats from the data that are needed (so then if you add new data does it not recalc those stats? so if normalizing test would it normalize based on training stats?) confusing!
# can use # now it's actually done the transformations
# this keeps data in memory, can use retain=F to not do that 
# this actually does the calculations but doesn't actually return the dataset to you
# to get the dataset back you can use bake with newdata=NULL if you retained the data
# this gives you back train data:
bake(simple_ames,new_data=NULL)
# but you don't need to do it
# but if you wanted to run it on test data you could:
test_ex <- bake(simple_ames, new_data = ames_test) # so 
# very unclear about prep vs bake! need more info. 
# high level functions do prep and bake automatically so let's not fuss too much for now. come back if still confused.

# qualitative data
# step_unknown() : chang mising values to a factor level ? ooh step_novel() can allot a new level if you anticipate coming across a new factor level in the future (may be useful for new species)
# can use step_other() to group infrequently observed factor levels into 'other' category (handy but not for my model)
# example for neighborhood:
ggplot(ames_train, aes(y = Neighborhood)) + 
  geom_bar() + 
  labs(y = NULL)
# Here there are two neighborhoods that have less than five properties in the training data; in this case, no houses at all in the Landmark neighborhood were included in the training set.  step_other(Neighborhood, threshold = 0.01) to our recipe, the bottom 1% of the neighborhoods will be lumped into a new level called “other”

simple_ames <- 
  recipe(Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type,
         data = ames_train) %>%
  step_log(Gr_Liv_Area, base = 10) %>% 
  step_other(Neighborhood, threshold = 0.01) %>% 
  step_dummy(all_nominal_predictors())
#when making dummy vars each factor level with the qual column will become its own column with a 1 and 0s 
# in r the convention is to leave out the first factor since if you know all others then you know which it is by process of elimination
# don't want linear dependence between columns (would add up to intercept ? not sure about this )
# full set of encodings can be used for some models = "one_hot" encoding and can be applied in step_dummy
# may want that for RF for species perhaps
# can use starts_with("feature_") to select particular shapes or features of the model cool
########### interactions: #############

ggplot(ames_train, aes(x = Gr_Liv_Area, y = 10^Sale_Price)) + 
  geom_point(alpha = .2) + 
  facet_wrap(~ Bldg_Type) + 
  geom_smooth(method = lm, formula = y ~ x, se = FALSE, col = "red") + 
  scale_x_log10() + 
  scale_y_log10() + 
  labs(x = "General Living Area", y = "Sale Price (USD)")

# see that regression line for price vs living area are different when faceted by bldg type

# in base R how would we specify this?
# Sale_Price ~ Neighborhood + log10(Gr_Liv_Area) + Bldg_Type + log10(Gr_Liv_Area):Bldg_Type # use the ":"
# Sale_Price ~ Neighborhood + log10(Gr_Liv_Area) * Bldg_Type # or * --> * expands the columns to the main effects *AND* the interaction term 
# does a lot of steps like first converting to dummy vars and then interacting for bldg type
# recipes are more explitic and sequential
# want to add an interaction term:
# step_interact(~ interaction terms) where terms on the right hand of tilde are interactions and can include selectors
# for example
step_interact( ~ Gr_Liv_Area:starts_with("Bldg_Type_") )
step_interact( ~ Gr_Liv_Area:starts_with("Bldg_Type_") )
# if you try to interact variables that haven't been made into dummys yet, a warning will be given!  
# note that quadratic terms are not generated  var1:var1 would be ignored by the recipe (important!)
simple_ames <- 
  recipe(Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type,
         data = ames_train) %>%
  step_log(Gr_Liv_Area, base = 10) %>% 
  step_other(Neighborhood, threshold = 0.01) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  # Gr_Liv_Area is on the log scale from a previous step
  step_interact( ~ Gr_Liv_Area:starts_with("Bldg_Type_")) # using starts with because of dummy var names all starting with that

prep(simple_ames,data=ames_train)
# if you skipped the step_dummy step then this would be thrown: 
#Warning message:
#  Categorical variables used in `step_interact` should probably be avoided;  This can lead to differences in dummy variable values #that are produced by `step_dummy`. Please convert all involved variables to dummy variables first. simple_ames
# if you didn't set up the "starts with" for the interaction then you'll get an error: 
#Warning message:
#  Interaction specification failed for: ~Gr_Liv_Area:Bldg_Type. No interactions will be created. 
# good to note that errors come at prep stage not when you make the recipe! 

###### TRANSOFRMATIONS OF OUTCOME COLUMN SHOULD OCCUR PRIOR TO (OUTSIDE) THE RECIPE! important!

#exception: downsampling 
# so if you want a step of the recipe only to be aplied ot the initial training data but not any subsequent data
# you can use skip=T for the step -- it'll be done when newdata=NULL but will be skipped all other times#
 ## THIS SEEMS DANGEROUS USE WITH CAUTION
# seems most important to downsampling which we aren't doing
######## feature extraction: pca ########
# create new features from the predictors that capture the information in the broader set as a whole.  
step_pca(matches("(SF$)|(Gr_Liv)")) # this is saying to do pca on ay of the square footage or general living area features (must all be on same scale!!!)
# this would remove the original features from the data and replace them with PCs - wow cools
# other examples:
#There are existing recipe steps for other extraction methods, such as: independent component analysis (ICA), non-negative matrix #factorization (NNMF), multidimensional scaling (MDS), uniform manifold approximation and projection (UMAP), and others.
# so does this replace those features with pca? 
# The argument num_comp controls the number of components that will be retained (the original variables that are used to derive the components are removed from the data). The new components will have names that begin with prefix and a sequence of numbers. The variable names are padded with zeros. For example, if num_comp < 10, their names will be PC1 - PC9. If num_comp = 101, the names would be PC001 - PC101.

########## using recipe with traditional baes R modeling #############
ames_rec <- 
  recipe(Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type + 
           Latitude + Longitude, data = ames_train) %>%
  step_log(Gr_Liv_Area, base = 10) %>%  # log transform
  step_other(Neighborhood, threshold = 0.01) %>% # convert all lt 0.1 to 'other'
  step_dummy(all_nominal_predictors()) %>%  # make dummies
  step_interact( ~ Gr_Liv_Area:starts_with("Bldg_Type_") ) %>% # interactions
  step_ns(Latitude, Longitude, deg_free = 20) # splines

# get it ready with prep: 
ames_rec_prepped <- prep(ames_rec)
# extract the prepped traiing set with bake and newdata=NULL
ames_train_prepped <- bake(ames_rec_prepped, new_data = NULL)
# get the test dataset with new ata:
ames_test_prepped <- bake(ames_rec_prepped, new_data =ames_test)
# Fit the model; Note that the column Sale_Price has already been
# log transformed.
lm_fit <- lm(Sale_Price ~ ., data = ames_train_prepped)
glance(lm_fit)
# get coefs with tidy:
tidy(lm_fit)
# to make predictions using the test set, use predict():
predict(lm_fit, ames_test_prepped)
# what do you do with these predictions then? get Rsq?
# tidy a recipe
tidy(ames_rec_prepped)
# can get tidy info about a particular step number:
tidy(ames_rec_prepped,number=3) # gives info about which columns were convete to dummy! useful!
tidy(ames_rec_prepped,number=2)
# tells you 

######## can change roles of columns ########
# what if you want to keep a column in as a label but don't want it to be a predictor or outcome
# ooh my motif label would be great for this!!! 
# you can add a 'new role' that isn't outcome or predictor  ames_rec %>% update_role(address, new_role = "street address")
ames_rec %>% update_role(address, new_role = "street address")
ames_rec
# okay so I want to ha
# update_role is when a variable doesn't currently have a role in the recipe of the model
# I want to do this with motif somehow though my ~ . would currently include motif if I left it in the df which isn't right. hm maybe I could do ~ (names- motif) or something like that.
#could do this: lm(y ~ . - excluded_1 - excluded_2, data = myFrame)
# sow hat I can do is fracSegSites ~ . -motif 
# and then add motif as a role to be the "motif label" or whatever and not be part of the model
# that is cool.

################ chapter 7: fitting models with parsnip ###########
# example ordinary vs regularized linear regression
# regularized linear regression adds a pentaly to the elast sq method to encourage smplicity by shirnking oceffs toward zero. 
# can use the glmnet model for regularized regression (non bayesian)
# in base R , different models have different input requirements. some want y ~ x formulae
# others want numeric matrices
# need to memorize each package's syntax -a hassle
# tidy models: instead:
# you specify the type of model (linear regression, random forest, etc.)
# then the 'enginge' (the software package to be used)
# then sometimes need the mode of the model :  regression or classification if those are options. for linear regression it's already set as regression, but RF could be either!
# these specificiations all have no refernce to the data set yet
linear_reg() %>% set_engine("lm")
linear_reg() %>% set_engine("glmnet") 
linear_reg() %>% set_engine("stan")
# once you've set details of the model then the model estimation is done using fit() to use a formula
# or fit_xy() when data already preprocesed (???)
# parsnip takes care of getting the data into the package
# translate() gives details about how parsnip converts the code
linear_reg() %>% set_engine("lm") %>% translate()
linear_reg(penalty = 1) %>% set_engine("glmnet") %>% translate()
linear_reg() %>% set_engine("stan") %>% translate()
# Note that missing_arg() is just a placeholder for the data that has yet to be provided.

##### predicting house price by long/lat:
lm_model <- 
  linear_reg() %>% 
  set_engine("lm")

lm_form_fit <- 
  lm_model %>% 
  # Recall that Sale_Price has been pre-logged
  fit(Sale_Price ~ Longitude + Latitude, data = ames_train)

lm_xy_fit <- 
  lm_model %>% 
  fit_xy(
    x = ames_train %>% select(Longitude, Latitude),
    y = ames_train %>% pull(Sale_Price)
  )

lm_form_fit
lm_xy_fit
# The differences between fit() and fit_xy() may not be obvious
#When fit() is used with a model specification, this almost always means that dummy variables will be created from qualitative predictors. If the underlying function requires a matrix (like glmnet), it will make them. However, if the underlying function uses a formula, fit() just passes the formula to that function. We estimate that 99% of modeling functions using formulas make dummy variables. The other 1% include tree-based methods that do not require purely numeric predictors.12

#The fit_xy() function always passes the data as-is to the underlying model function. It will not create dummy variables before doing so.
# so that is a bit dangerous! fix_xy won't change the data at all - so stick with fit for most of the time?
# random forest examples: Three commonly used arguments are the number of trees in the ensemble, the number of predictors to randomly sample with each split within a tree, and the number of data points required to make a split. these have diff arg nanes for diff engings, so parsnip uses a common set of names for ease of use
# parsnip calls lambda 'penalty' to be more explicit
# number of neighbors is called neighbors rather than k
# can use the help fie for the model or translate to see what things mean
?rand_forest
# mtry	
# An integer for the number of predictors that will be randomly sampled at each split when creating the tree models.
# trees	
# An integer for the number of trees contained in the ensemble.
# min_n	
# An integer for the minimum number of data points in a node that are required for the node to be split further.


rand_forest(trees = 1000, min_n = 5) %>% 
  set_engine("ranger") %>% 
  set_mode("regression") %>% 
  translate()
# see how it translates to a ranger command

##### use model results ###########
# once you have your fit model you can use the results. 
# a parsnip model object stores the fitted model and other outputs. the fitted model is found in element called fit which can be found with "purrr" --pluck
lm_form_fit %>% pluck("fit")
lm_form_fit$fit
# these do same thing
#***** Never pass fit element to predict function! you can do predict(m_form_fit) but not predict(lm_form_fit$fit) !!!!! if any preprocessing was done then incorrect preds can be generated #
# wow that seems dangerous. did I ever do that - no . but be careful!

# so you could use summary and coef to get info (old way in R)
# save results:
model_res <- 
  lm_form_fit %>% 
  pluck("fit") %>% 
  summary()
# but this isn't a very useful data structure
# use broom::tidy to tidy it up
tidy(lm_form_fit)

######## predict ###########
# results always a tibble, predictable column names, as many rows in tibble as in input datase
# The row order of the predictions are always the same as the original data.! 
# These three rules make it easier to merge predictions with the original data:
ames_test_small <- ames_test %>% slice(1:5)
predict(lm_form_fit, new_data = ames_test_small)
# 

# why are there leading dot in some of the column names? Some tidyverse and tidymodels arguments and return values contain periods. This is to protect against merging data with duplicate names. There are some data sets that contain predictors names pred!

## merge wtih original sale price test data:
ames_test_small %>% 
  select(Sale_Price) %>% # select the obs data 
  bind_cols(predict(lm_form_fit, ames_test_small)) %>%  # bind columns together 
  # Add 95% prediction intervals to the results: 
  bind_cols(predict(lm_form_fit, ames_test_small, type = "pred_int")) 

# usemodels package takes a dataframe and formula then write code to use it. creates a recipe. uses a workflow object to bind together. 
# example:
use_xgboost(Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type + 
              Latitude + Longitude, 
            data = ames_train,
            # Don't create the model tuning code:
            tune = FALSE,
            # Add comments explaining some of the code:
            verbose = TRUE)
# WOW this is cool. it gives you tips on the model! what happens when I do RF?
use_ranger(Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type + 
              Latitude + Longitude, 
            data = ames_train,
            # Don't create the model tuning code (but can still have some tune parameters?)
            tune = FALSE,
            # Add comments explaining some of the code:
            verbose = TRUE)
# writes out code for you:
ranger_recipe <- 
  recipe(formula = Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type + 
           Latitude + Longitude, data = ames_train) 

ranger_spec <- 
  rand_forest(trees = 1000) %>% 
  set_mode("regression") %>% 
  set_engine("ranger") 

# to explore other models can have an addin 
parsnip_addin() # wow

######### Chapter 8: workflows ##################

# the model up to now has meant the structural equation 
lm_model <- 
  linear_reg() %>% 
  set_engine("lm")

lm_wflow <- 
  workflow() %>% 
  add_model(lm_model)
# needs a parsnip model obj
lm_wflow
# no preproc specified yet
# If our model were very simple, a standard R formula can be used as a preprocessor:
lm_wflow <- 
  lm_wflow %>% 
  add_formula(Sale_Price ~ Longitude + Latitude)
lm_wflow

# Workflows have a fit() method that can be used to create the model. 
lm_fit <- fit(lm_wflow, ames_train)
lm_fit

# We can also predict() on the fitted workflow:
predict(lm_fit,ames_test) # predition on test dataset

# The predict() method follows all of the same rules and naming conventions that we described for the parsnip package in Section 7.3.

#Both the model and preprocessor can be removed or updated:
 
lm_fit %>% update_formula(Sale_Price~Longitude)  # gettint rid of latitude 

# but instead of using formulae directly you can add a recipe to do the preprocessing
lm_wflow %>% 
  add_recipe(ames_rec)
# That did not work! We can only have one preprocessing method at a time, so we need to remove the formula before adding the recipe.
lm_wflow <- 
  lm_wflow %>% 
  remove_formula() %>% 
  add_recipe(ames_rec)
lm_wflow
# We described the prep() and bake() functions in Section 6.8 for using the recipe with a modeling function. This can be onerous, so the fit() method for workflow objects automates this process:
# Does `prep()`, `bake()`, and `fit()` in one step:
lm_fit <- fit(lm_wflow, ames_train)

# # Does `bake()` and `predict()` automatically:
predict(lm_fit, ames_test %>% slice(1:3)) 

# If we need the bare model object or recipe, there are pull_* functions that can retrieve them:
lm_fit %>% 
  pull_workflow_prepped_recipe() %>% 
  tidy()
# gives the info on the workflow 

# to get model fit:
# To tidy the model fit: 
lm_fit %>% 
  # This returns the parsnip object:
  pull_workflow_fit() %>% 
  # Now tidy the linear model object:
  tidy() %>% 
  slice(1:5)

####### can make a set of workflows using
library(workflowsets)
# set of possible models:
location <- list(
  longitude = Sale_Price ~ Longitude,
  latitude = Sale_Price ~ Latitude,
  coords = Sale_Price ~ Longitude + Latitude,
  neighborhood = Sale_Price ~ Neighborhood
)
location_models <- workflow_set(preproc = location, models = list(lm = lm_model))
location_models
# then can pull out a particular workflow
pull_workflow(location_models, id = "coords_lm")

# this might be useful! I could have different models with species etc. in it
# Workflow sets are mostly designed to work with resampling, which is discussed in Chapter 10. In the object above, the columns option and result must be populated with specific types of objects that result from resampling. We will demonstrate this in more detail in Chapters 11 and 15.

# create model fits for every formula:
location_models <-
  location_models %>%
  mutate(fit = map(info, ~ fit(.x$workflow[[1]], ames_train)))
location_models
location_models$fit
location_models$fit[[1]]
# There’s a lot more to workflow sets. Their nuances and advantages won’t be illustrated until Chapter 15.

########## chapter 9 judging model effectiveness #########
# RMSE measures accuracy; Rsq measures correlation (coefficient of determination)
# a model optimizing rmse has more variability but uniform accuracy across range of outcomes (!)
# a model for rsq shows tighter correlation between obs and pred but may be bad at tails (huH!)

# this chapter uses yardstick package
#  This chapter focuses on functions that can be used to measure predictive strength.
# seemsˆmportnat
#One missing piece of information in this approach is how closely this model fits the actual data. Using resampling methods, discussed in Chapter 10, we can estimate the accuracy of this model to be about 73.3%. Accuracy is often a poor measure of model performance; we use it here because it is commonly understood. If the model has 73.3% fidelity to the data, should we trust the conclusions produced by the model? We might think so until we realize that the baseline rate of non-impaired patients in the data is 72.7%. This means that, despite our statistical analysis, the two-factor model appears to be only 0.6% better than a simple heuristic that always predicts patients to be unimpaired, irregardless of the observed data.
# hmmm 
# *** optimization of statistical characteristics of the model does not imply that the model fits the data well. **** (just by getting significant coeffs and interaction terms)

# #  Even for purely inferential models, some measure of fidelity to the data should accompany the inferential results. Using this, the consumers of the analyses can calibrate their expectations of the results of the statistical analysis. 

# come back to this chapter after ch 10

######## chapter 10 ###############
####### so want to know the efficacy of model BEFORE running test set 

# random forests: require very little pre-processing
rf_model <- 
  rand_forest(trees = 1000) %>% 
  set_engine("ranger") %>% 
  set_mode("regression")

rf_wflow <- 
  workflow() %>% 
  add_formula(
    Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type + 
      Latitude + Longitude) %>% 
  add_model(rf_model) 


rf_fit <- rf_wflow %>% fit(data = ames_train) # get model fit
rf_fit

# predict the training set to generate the 'apparent error rate' or resub error rate":

estimate_perf <- function(model, dat) {
  # Capture the names of the objects used
  cl <- match.call()
  obj_name <- as.character(cl$model)
  data_name <- as.character(cl$dat)
  data_name <- gsub("ames_", "", data_name)
  
  # Estimate these metrics:
  reg_metrics <- metric_set(rmse, rsq)
  
  model %>% 
    predict(dat) %>% 
    bind_cols(dat %>% select(Sale_Price)) %>% 
    reg_metrics(Sale_Price, .pred) %>% 
    select(-.estimator) %>% 
    mutate(object = obj_name, data = data_name)
}
  
estimate_perf(rf_fit, ames_train)
# compare to lm:
estimate_perf(lm_fit, ames_train)
# Based on these results, the random forest is much more capable of predicting the sale prices; the RMSE estimate is 2.06-fold better than linear regression. If these two models were under consideration for this prediction problem, the random forest would probably be chosen. The next step applies the random forest model to the test set for final verification:
estimate_perf(rf_fit, ames_test)
# The test set RMSE estimate, 0.0695, is much worse than the training set value of 0.0364! Why did this happen? !!! <---
# For a low-bias model, the high degree of predictive capacity can sometimes result in the model nearly memorizing the training set data. As an obvious example, consider a 1-nearest neighbor model. It will always provide perfect predictions for the training set no matter how well it truly works for other data sets. Random forest models are similar; re-predicting the training set will always result in an artificially optimistic estimate of performance.

# resample just the train data and then do analysis/assessment (named diff to separate from train/test)
# These are somewhat analogous to training and test sets. Our language of analysis and assessment avoids confusion with initial split of the data. # 

# Suppose twenty iterations of resampling are conducted. This means that twenty separate models are fit on the analysis sets and the corresponding assessment sets produce twenty sets of performance statistics. The final estimate of performance for a model is the average of the twenty replicates of the statistics. This average has very good generalization properties and is far better than the resubstituion estimates.
############ cross validation #############
# Cross-validation is a well established resampling method. While there are a number of variations, the most common cross-validation method is V-fold cross-validation. 
# 
set.seed(55)
ames_folds <- vfold_cv(ames_train, v = 10)
ames_folds
# The column named splits contains the information on how to split the data (similar to the object used to create the initial training/test partition). While each row of splits has an embedded copy of the entire training set, R is smart enough not to make copies of the data in memory14. The print method inside of the tibble shows the frequency of each: [2K/220] indicates that roughly two thousand samples are in the analysis set and 220 are in that particular assessment set.

# To manually retrieve the partitioned data, the analysis() and assessment() functions return the corresponding data frames:
# For the first fold:
# first split:
ames_folds$splits[[1]]

ames_folds$splits[[1]] %>% analysis() # get the analysis data of the fold
# usually don't need to extract manually though
# each have an id fold01 etc
# can repeat the CV multiple times to reduce noise
# you can add "repeats"
vfold_cv(ames_train, v = 10, repeats = 5)
# 
# one other kind: Leave-one-out cross-validation: an early verse of CV are deficient for anything -- don't use it. focus on V-fold CV
# mone carlo CV (MCCV) proportion of data is hanged every time (picked at random) results in assessment sets that are not mutally exclusive

# validation sets: a partition set aside to estimate performance, before test set

# keep out part of the training data 
# set.seed(12)
val_set <- validation_split(ames_train, prop = 3/4)
val_set
# but if you're doing CV I don't think you need to this?
# should I make a validation set by splitting the 

# can get boostraps (assessments sometimes called out of bag)
bootstraps(ames_train, times = 5)
# bootstraps are very pessimistic (why?); very low variance . random forests use bootstraps internally -- 1000 individual decision trees : each one is a different bootstrap sample of the training set (ah so that is where the out of bag error comes from). hm so there is replacement. do I need more trees?

############### estimate performace on folds ########
keep_pred <- control_resamples(save_pred = TRUE, save_workflow = TRUE)
# ^ save predictions
set.seed(130)
rf_res <- 
  rf_wflow %>% 
  fit_resamples(resamples = ames_folds, control = keep_pred)
rf_res
# fit_resamples fits your folds 

# to tidy up:
collect_metrics(rf_res) 
collect_metrics(rf_res,summarize = F) # to get predictions per fold 

# get the assessment predictions:
assess_res <- collect_predictions(rf_res) # from rf_res$.predictions
assess_res

# note the .row column is the row number from the original data


# can plot results:
assess_res %>% 
  ggplot(aes(x = Sale_Price, y = .pred)) + 
  geom_point(alpha = .15) +
  geom_abline(col = "red") + 
  coord_obs_pred() + 
  ylab("Predicted")

# find the one overpredicted house:
over_predicted <- 
  assess_res %>% 
  mutate(residual = Sale_Price - .pred) %>% 
  arrange(desc(abs(residual))) %>% 
  slice(1)
over_predicted
# with its row info you can pull it out 
ames_train %>% 
  slice(over_predicted$.row) %>% 
  select(Gr_Liv_Area, Neighborhood, Year_Built, Bedroom_AbvGr, Full_Bath)

###### options for parallel processing ######

####### saving resampled objects ##########
# models created during resamp are not retained. usually don't need them after computing perforamcnce stats. 
# if particularl approahch ends up being best option, then best to fit the whole model again tot he full training set ** oh okay this interesting
# # but if yo uwant to keep them you can use extract 
# example: linear model:
ames_rec <- 
  recipe(Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type + 
           Latitude + Longitude, data = ames_train) %>%
  step_other(Neighborhood, threshold = 0.01) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_interact( ~ Gr_Liv_Area:starts_with("Bldg_Type_") ) %>% 
  step_ns(Latitude, Longitude, deg_free = 20)

lm_wflow <-  
  workflow() %>% 
  add_recipe(ames_rec) %>% 
  add_model(linear_reg() %>% set_engine("lm")) 

lm_fit <- lm_wflow %>% fit(data = ames_train)

# can save the linear models coeffs:
get_model <- function(x) {
  pull_workflow_fit(x) %>% tidy()
}
get_model(lm_fit)

# now try with resampled fits (10 fold CV)
ctrl <- control_resamples(extract = get_model)
lm_res <- lm_wflow %>%  fit_resamples(resamples = ames_folds, control = ctrl)
lm_res # there's a new column: .extracts
# which contains what?
lm_res$.extracts
lm_res$.extracts[[1]]
lm_res$.extracts[[1]][[1]] # getting the model out for one of the folds
# flatten and get all
all_coef <- map_dfr(lm_res$.extracts, ~ .x[[1]][[1]])
all_coef # this just stacks all the models together (so there are repeated "Intercept" rows etc. one for each of the 10 mdoels)
# you can then show all replicates for a single predictor
filter(all_coef, term == "Year_Built")

############### multiple models ##############

# make multiple recipes to compare different lms
# first recipe :
basic_rec <- 
  recipe(Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type + 
           Latitude + Longitude, data = ames_train) %>%
  step_log(Gr_Liv_Area, base = 10) %>% 
  step_other(Neighborhood, threshold = 0.01) %>% 
  step_dummy(all_nominal_predictors())

# recipe with interaction terms:
interaction_rec <- 
  basic_rec %>% 
  step_interact( ~ Gr_Liv_Area:starts_with("Bldg_Type_") ) 

# recipe with splines
spline_rec <- 
  interaction_rec %>% 
  step_ns(Latitude, Longitude, deg_free = 50)

# list of all possible ways to preprocess (all recipes)
preproc <- 
  list(basic = basic_rec, 
       interact = interaction_rec, 
       splines = spline_rec
  )
preproc
# makes a workflow set with the preproc list 
lm_models <- workflow_set(preproc, list(lm = lm_model), cross = FALSE)
lm_models

# want to do resampling so going to use map to map the fit_resamples function to each of them: 
lm_models <- 
  lm_models %>% 
  workflow_map("fit_resamples", 
               # Options to `workflow_map()`: 
               seed = 1101, verbose = TRUE,
               # Options to `fit_resamples()`: 
               resamples = ames_folds, control = keep_pred)
lm_models


# can add old RF model : 
# convert it to a workflow set (n=1) and add it in
four_models <- 
  as_workflow_set(random_forest = rf_res) %>% 
  bind_rows(lm_models)
four_models

# plot it
autoplot(four_models, metric = "rsq")

# collect rsq metrics
rsq_indiv_estimates <- 
  collect_metrics(four_models, summarize = FALSE) %>% 
  filter(.metric == "rsq") 
rsq_indiv_estimates
# make dataset wide 
rsq_wider <- 
  rsq_indiv_estimates %>% 
  select(wflow_id, .estimate, id) %>% 
  pivot_wider(id_cols = "id", names_from = "wflow_id", values_from = ".estimate")
rsq_wider

# looka t correlation
corrr::correlate(rsq_wider %>% select(-id), quiet = TRUE)

# indicates that each fold is very correlated across models because same folds were used

# you can plot it:
rsq_indiv_estimates %>% 
  mutate(wflow_id = reorder(wflow_id, .estimate)) %>% 
  ggplot(aes(x = wflow_id, y = .estimate, group = id, col = id)) + 
  geom_line(alpha = .5, lwd = 1.25) + 
  theme(legend.position = "none") + 
  labs(x = NULL, y = expression(paste(R^2, "statistics")))
 # If the resample-to-resample effect was not real, there would not be any parallel lines. A statistical test for the correlations evaluates whether the magnitudes of these correlations are not simply noise
rsq_wider %>% 
  with( cor.test(basic_lm, splines_lm) ) %>% 
  tidy() %>% 
  select(estimate, starts_with("conf"))
# highly correlated

# In other words, ignoring the resample-to-resample effect would bias our model comparisons towards finding no differences between models.! 

# try to come up with a practical effective size. eg two models not practically different if diff in rsq < 2%
# Practical significance is subjective; two people can have very different ideas on the threshold for importance. However, as shown later, this consideration can be very helpful when deciding between models.

######### simple hypothesis testing methods ##########
############ ANOVA
# ANOVA: a version of regression. predictrs are binary dummy vars for different groups, betas estimate whetehr two or more groups are different 
# can be used for model comparisons
# so if rsqs are outcome data (ys) and models are the predictors 
# so then you can have Xs are dummy vars for model labels and the folds 
# b0 is the estimate of mean rsq for basic linear models (no splines/interactions)
# b1 is change in mean rsq when interactions added
# b2 is change in rsq when you change to RF
# b3 is chane betweeen basic inear and one with ints and splies
# generates p values fo differences . resample groups would be a 'block' effect that all have in common
# and added as a term to th model. or it could be a random effect. but we don't care which, we just want to adjust for them. 
# to compare two models at a time is to use diffs in Rsq as outcome data in the anova. then since outcomes are matches by resample, teh diffs don't contain the resample effect

####### compare just a pair of models ##########
# so pair them by the diff in rsqs:
compare_lm <- 
  rsq_wider %>% 
  mutate(difference = splines_lm - basic_lm)
compare_lm

# anova using lm:
lm(difference ~ 1, data = compare_lm) %>% 
  tidy(conf.int = TRUE) %>% 
  select(estimate, p.value, starts_with("conf"))
#> # A tibble: 1 x 4
#>   estimate p.value conf.low conf.high
#>      <dbl>   <dbl>    <dbl>     <dbl>
#> 1  0.00797 0.00660  0.00284    0.0131

# Alternatively, a paired t-test could also be used: 
rsq_wider %>% 
  with( t.test(splines_lm, basic_lm, paired = TRUE) ) %>%
  tidy() %>% 
  select(estimate, p.value, starts_with("conf"))

########### Bayesian comparison methods #########
# have prior distributions on the parameters eg normal, exponential etc. 
# see book 1..4

######## model tuning and avoiding overfitting! ################
# boosting: an ensemble method that combines a series of base models each of which is sequential and depnds on preivous. number of boostering iterations needs to be tuned
# nearal network 
# should never tune a prior
# ## interesting take on Random Forest:
# should not tune the number of trees in RF or baggin model. Should be chosen to be large enough to ensure numerical stability in the reuslts (?); tuning it can't improve perforamnce as long as value is large enoguh for reliable results. Typically need thousands of trees , and for bagging 50-100
# that's interesting!!  

## okay so random forest examples: main args: trees, min_n, mtry (need to understand all these)
# then engine specific things you can tune
# like for random forest in ranger you can tune: regularization.factor (modulate trade off between number of predictors) - regularization. is only in ranger. maybe I want to tune this?
# is it bad to tune too many things?? 
rand_forest(trees = 2000, min_n = 10) %>%                   # <- main arguments
  set_engine("ranger", regularization.factor = 0.5)         # <- engine-specific 

# can mark a param for tunin by setting it =tune():
neural_net_spec <- 
  mlp(hidden_units = tune()) %>%  # setting hiddnen units to be tuned
  set_engine("keras")
parameters(neural_net_spec)
# it specifies which are for tuning:

rf_experiment <-
  rand_forest(min_n=tune(),mtry=tune()) %>%
  set_engine("ranger")
rf_experiment

# if ou are tuning different things you could give the parameter a name inside tune:
# mtry = tune("mtry1")
ames_rec <- 
  recipe(Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type + 
           Latitude + Longitude, data = ames_train)  %>%
  step_log(Gr_Liv_Area, base = 10) %>% 
  step_other(Neighborhood, threshold = tune()) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_interact( ~ Gr_Liv_Area:starts_with("Bldg_Type_") ) %>% 
  step_ns(Longitude, deg_free = tune("longitude df")) %>% 
  step_ns(Latitude,  deg_free = tune("latitude df"))

recipes_param <- parameters(ames_rec)
recipes_param
# shows up in the identifier column

# can combine recipe and model:
wflow_param <- 
  workflow() %>% 
  add_recipe(ames_rec) %>% 
  add_model(neural_net_spec) %>% 
  parameters()
wflow_param

# see all the paramters for tuning (threshold, long and lat tuning in recipe pre-processing step and hidden units tuning in the model step)
# range of possible values for tuning:
# primary tuning param for RF is mtre : the number of predictor columns that re randomly sampled for each split in the tree (?). requires finalizeation because don't know number of predcitors:
########### random forest tuning regularization and mtry! useful! ############
rf_spec <- 
  rand_forest(mtry = tune()) %>% 
  set_engine("ranger", regularization.factor = tune("regularization"))
# going to tune regularization and mtry ^

rf_spec

rf_param <- parameters(rf_spec)
rf_param

## aha! if a parameter objective is complete it will say nparam[+] but if it is missing at least one end of the range 
# so you could use update to add the range based on what you know about the data dimensions:
rf_param %>% 
  update(mtry = mtry(c(1, 70)))
# but maybe columns get added or subtracted or you don't know how many cols there are
# can use finalize to execute the recipe to obtain dimensions:
pca_rec <- 
  recipe(Sale_Price ~ ., data = ames_train) %>% 
  # Select the square-footage predictors and extract their PCA components:
  step_normalize(contains("SF")) %>% 
  # Select the number of components needed to capture 95% of
  # the variance in the predictors. 
  step_pca(contains("SF"), threshold = .95)


updated_param <- 
  workflow() %>% 
  add_model(rf_spec) %>% 
  add_recipe(pca_rec) %>% 
  parameters() %>% 
  finalize(ames_train)
updated_param # now they are complete 

updated_param %>% pull_dials_object("mtry")
# When the recipe is prepared, the finalize() function learns to set the upper range of mtry to 74 predictors.
rf_param
# some parameters (like lasso penalty) are in log units. need to make sure range is in the same units
# correct:
penalty(c(-1, 0)) %>% value_sample(1000) %>% summary()
# incorrect (not log units)
penalty(c(0.1, 1.0)) %>% value_sample(1000) %>% summary()
# THAT IS TRICKY!

# to pull out param information:
rf_spec %>% pull_dials_object("mtry")
########## regular grids with grid_regular ############
# can create parameter grids with grid_regular and levels -- levels lets you set how many you want per parameter
updated_param %>% grid_regular(levels=c(mtry=5,regularization=3))
# can be computationally expensive


#One advantage to using a regular grid is that the relationships and patterns between the tuning parameters and the model metrics are easily understood. The factorial nature of these designs allows for examination of each parameter separately with little confounding between parameters.

########## irregular grids  ############
set.seed(10)
updated_param %>% 
  grid_random(size = 1000) %>% # 'size' is the number of combinations
  summary()
 
# randomly samples unfirom fandom numbers across parameter range 

# but can result in overlapping parameter combos; need more grid values for good coverage
# better to use 'space filling designs' 
# if I need them: Latin hypercube and maximum entropy designs 

######## evaluating the grid ##########

####### finalizing the model########
# note tuning doesn't pick final model
# manually pick values or use select : 
# select_best() will choose parameters with numerically best results

select_best()
select_best(mlp_reg_tune, metric = "roc_auc") # example 
# make a tibble manually with these
logistic_param <- 
  tibble(
    num_comp = 0,
    epochs = 125,
    hidden_units = 1,
    penalty = 1
  )
# and feed them to workflow
final_mlp_wflow <- 
  mlp_wflow %>% 
  finalize_workflow(logistic_param)
final_mlp_wflow

######## there's much more advanced tuning you can do - chapte r14 #########

######## screening many models: use if i want to do lots of models ##########