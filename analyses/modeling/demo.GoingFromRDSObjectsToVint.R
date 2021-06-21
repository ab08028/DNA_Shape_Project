########## this demonstrates a way to read in RDS objects (eg on hoffman) and run them through vint ############
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
require(devtools)

# for doing this on hoffman you could set 
indir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210616_RF.Log10MutationRate.Multispecies.DummyPopVar/" # set this to the date you want in hoffman and then run vint with a LOT of cores and maybe less mem per core 
oneFoldSetToTrainAndAssessOn <- readRDS(paste0(indir,"oneFoldSetToTrainAndAssessOn.rds"))
model <- readRDS(paste0(indir,"modelTrainedOnOneFold.rds"))

ranger_obj <- pull_workflow_fit(model)$fit
ranger_obj

# need the correct recipe:
recipe <- readRDS(paste0(indir,"recipe.rds"))

analysisdf_processed <- prep(recipe, analysis(oneFoldSetToTrainAndAssessOn)) %>% 
  juice() %>% 
  select(-c(mutationType)) 
head(analysisdf_processed)

interact <- vint(ranger_obj,feature_names=c("feature_1_Shear_2.derived","feature_1_Shear_2.ancestral",progress=T),train=analysisdf_processed,progress = "text",parallel=T) # CRASHED
