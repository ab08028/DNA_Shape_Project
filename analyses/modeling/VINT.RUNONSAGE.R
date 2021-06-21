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
args <- commandArgs(trailingOnly = TRUE)
indir <- paste0(args[1],"/") # doing this to make sure slash comes after
 
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

#interact <- vint(ranger_obj,feature_names=c("feature_1_Shear_2.derived","feature_1_Shear_2.ancestral",progress=T),train=analysisdf_processed,progress = "text",parallel=T) # CRASHED

vi_scores <- vip::vi(ranger_obj)
#write.table(vi_scores,paste0(outdir,"modelTrainedOnOneFold.VIPScores.PerMutatationType.AssessedOnChr",toString(unique(truth_prediction_df$window)),".txt"),quote = F,row.names=F,sep="\t")
# just doing top 5 scores *for now* for interactions
top_vi_scores <- 
  vi_scores %>%
  filter(rank(desc(Importance))<=5)

#vip_plot <- vip(vi_scores,include_type = T,num_features = 100)
#vip_plot

#ggsave(paste0(outdir,"modelTrainedOnOneFold.VIP.Plot.AssessedOnChr",toString(unique(truth_prediction_df$window)),".png"),vip_plot,height=12,width=5)
#vip(ranger_obj,include_type = T,num_features = 20)

########## try to find interactions? ############
##### juice the training data:
# maybe save this as an rds object? save the recipe as an object? what do I need?

# get all pairs of top scores:
pairsOfTopScores <- combn(top_vi_scores$Variable,m=2,simplify = T) # get all pairs of top scores 
allInteractionResults <- NULL
for (i in seq(1,dim(pairsOfTopScores)[2])){
  #print(pairsOfTopScores[,i])
  interact <- vint(ranger_obj,feature_names=pairsOfTopScores[,i],train=analysisdf_processed,progress = "text",parallel=T) 
  saveRDS(interact,paste0(outdir,"vint.intermediate.",i,".rds"))
  allInteractionResults <- rbind(allInteractionResults,interact)
  }
saveRDS(allInteractionResults,paste0(outdir,"vint.rds"))

