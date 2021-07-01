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
# try adding this:
library(foreach)
library(doParallel)
#registerDoParallel()
# for doing this on hoffman you could set 
args <- commandArgs(trailingOnly = TRUE)
indir <- paste0(args[1],"/") # doing this to make sure slash comes after
topXParamsToCompare <- args[2] # number of top params to compare (will affect size of array) == start with 5.
parameter_pair_index <- as.numeric(args[3]) # give the index when you submit the script so you can do parameter pairs in parallel!  

oneFoldSetToTrainAndAssessOn <- readRDS(paste0(indir,"oneFoldSetToTrainAndAssessOn.rds"))

windowOfAssessment=toString(unique(assessment(oneFoldSetToTrainAndAssessOn)$window))

model <- readRDS(paste0(indir,"modelTrainedOnOneFold.rds"))

ranger_obj <- pull_workflow_fit(model)$fit
ranger_obj

# need the correct recipe:
recipe <- readRDS(paste0(indir,"recipe.rds"))

analysisdf_processed <- prep(recipe, analysis(oneFoldSetToTrainAndAssessOn)) %>% 
  juice() %>% 
  select(-c(mutationType)) 

head(analysisdf_processed)

#interact <- vint(ranger_obj,feature_names=c("feature_1_Shear_2.derived","feature_1_Shear_2.ancestral",progress=T),train=analysisdf_processed,progress = "text",parallel=T) # CRASHED

#vi_scores <- vip::vi(ranger_obj) ### I'm being dumb here -- this may not order parameters in same way across runs. Need to be more thoughtful
#write.table(vi_scores,paste0(outdir,"modelTrainedOnOneFold.VIPScores.PerMutatationType.AssessedOnChr",toString(unique(truth_prediction_df$window)),".txt"),quote = F,row.names=F,sep="\t")
# just doing top 5 scores *for now* for interactions
vi_scores <- read.table(paste0(indir,"modelTrainedOnOneFold.VIPScores.PerMutatationType.AssessedOnChr",toString(windowOfAssessment),".txt"),sep="\t",header=T)

top_vi_scores <- 
  vi_scores %>%
  filter(rank(desc(Importance))<=topXParamsToCompare)

#vip_plot <- vip(vi_scores,include_type = T,num_features = 100)
#vip_plot

#ggsave(paste0(outdir,"modelTrainedOnOneFold.VIP.Plot.AssessedOnChr",toString(unique(truth_prediction_df$window)),".png"),vip_plot,height=12,width=5)
#vip(ranger_obj,include_type = T,num_features = 20)

########## try to find interactions? ############
##### juice the training data:
# maybe save this as an rds object? save the recipe as an object? what do I need?

# get all pairs of top scores:
pairsOfTopScores <- combn(top_vi_scores$Variable,m=2,simplify = T) # get all pairs of top scores  == I checked that combn will return things in the same order, so I can use SGE_TASK_ID to pull different rows of combn out in parallel 
#allInteractionResults <- NULL
#for (i in seq(1,dim(pairsOfTopScores)[2])){
  #print(pairsOfTopScores[,i])

print('starting parameter pair:', toString(pairsOfTopScores[,parameter_pair_index]))
interact <- vint(ranger_obj,feature_names=pairsOfTopScores[,parameter_pair_index],train=analysisdf_processed,progress = "text",parallel=T) 

write.table(data.frame(interact),paste0(indir,"vint.parallel.",parameter_pair_index,".txt"),rownames=F,quote=F,sep="\t")
#allInteractionResults <- rbind(allInteractionResults,interact)
#  }
#saveRDS(allInteractionResults,paste0(outdir,"vint.rds"))
#This function quantifies the strength of interaction between features $X_1$ and $X_2$ by measuring the change in variance along slices of the partial dependence of $X_1$ and $X_2$ on the target $Y$. See Greenwell et al. (2018) for details and examples.
