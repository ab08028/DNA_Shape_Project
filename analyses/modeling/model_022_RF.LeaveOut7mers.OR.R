############# leave out 7mer model with Odds Ratio ##########
# make a list of train and test 7mers
# want it to be random each time, 3/4 1/4
# but want all ancestral 7mers to be present in same group
# so need to make sure those group together
# so really what you want to split is ancestral 7mers into 3/4 and 1/4
# and then have all their >X mutations 
# aha. so I could come up with train/test list of ancestral 7mers
# and then label df by anc 7mer %in% train or %in% test and split on that
# have a seed so I can repeat, but eventually figure out bootstrapping (run as an array eventually?)
# eventually need to parallelize and save seeds et c
# will need to save train/test list

# split into train test
# rescale mutation rates by train / test 
# (rank?)
# issue will be that test will be 3x bigger but hopefully OR will deal with that?
# then divide sp A / sp B to get OR
# what about dividing by avg number of mutations? 
# then for pair of species get OR
# test cases: mmd and ms; mmd and EUR? human? (maybe?) 
# input: genome-wide mutation spectrum per population 
# # 

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

outdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210903_model_022_modelingUnseen7mers_ORs/"
dir.create(outdir,showWarnings = F)
############## read in shapes ##############
# for testing on home computer: 
shapes <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/METHYLATEDfirstOrder_featureTypes_allPossible7mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.METHYLATED.txt",header=T,sep="\t")
rownames(shapes) <- shapes$motif


##### read in spectra and sum up fully over all train and test chromosomes #######
#allData_multipop <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_spectrum_target_merged_PerChr/TEMPORARY.NEWGROUPS.MULTIPOPULATION_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.INCLUDES0Entries.txt",header=T)
# combining genome wide spectra
# note these aren't repeat masked
genomeWideSpectrum <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/summaryTables/Mouse7merSpectrum.SummedUpOverAllChrs.NOSTRICT.txt",header=T)

head(genomeWideSpectrum)
# get anc 7mer:
genomeWideSpectrum$ancestral7mer <- unlist(lapply(strsplit(genomeWideSpectrum$variable,"\\."),"[",1))
genomeWideSpectrum$derived7mer <- unlist(lapply(strsplit(genomeWideSpectrum$variable,"\\."),"[",2))
head(genomeWideSpectrum)
tail(genomeWideSpectrum)
genomeWideTargets <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/summaryTables/ALLINTERVALS.mouse.mutyper.targets.7mer.nostrict.txt",header=T)

# merge them:
head(genomeWideTargets)
colnames(genomeWideTargets) <- c("ancestral7mer","targetCount_allIntervals")

allData_multipop <- merge(genomeWideSpectrum,genomeWideTargets,by="ancestral7mer")
# note that these mice have same target regions because mapped to same ref genome
#dim(allData_multipop)
#dim(genomeWideSpectrum)
head(allData_multipop)
# okay potential issue here: what about 0 mutation rates that OR will be 0 or inf.
allData_multipop$mutationRate_unscaled <- allData_multipop$totalSitesAllIntervals/allData_multipop$targetCount_allIntervals
# missing in Ms:
#Ms CGACGCG.CGAGGCG    
#Ms CGACGGT.CGAGGGT  
#Ms GCGACGA.GCGCCGA  
# so at least FOR NOW: exclude the 
head(allData_multipop)

###### add in shape info: 
allData_multipop_intermediate <-merge(allData_multipop,shapes,by.x="ancestral7mer",by.y="motif")

allData_withShapes_unprocessed <- merge(allData_multipop_intermediate,shapes,by.x="derived7mer",by.y="motif",suffixes=c(".ancestral",".derived"))

head(allData_withShapes_unprocessed)

######## add seq info (not yet 1-hot encoded, will encode with tidy models )##########
# just want ancestral 7mer 
allData_withShapes_unprocessed <- allData_withShapes_unprocessed %>% separate(ancestral7mer,into=c("Pos_1","Pos_2","Pos_3","Pos_4.ancestral","Pos_5","Pos_6","Pos_7"),remove=F,sep=c(1,2,3,4,5,6))
# also get derived bp:
allData_withShapes_unprocessed$Pos_4.derived <- substr(allData_withShapes_unprocessed$derived7mer,4,4)


###### DO NOT RESCALE OUTCOME AT THIS STAGE : WOULD RESULT IN DATA LEAKAGE #######
################# subset species ###############
# just want two species for now.
speciesForModel = c("Mmd","Ms")

allData_withShapes_unprocessed <- allData_withShapes_unprocessed[allData_withShapes_unprocessed$population %in% speciesForModel,]

################# split the data ##############
# keep 3/4 of 7mers in model
# want to keep 3/4 of ANCESTRAL 7mers in train and 1/4 in test and all their associated mutations
# this is to  keep any info about a test ancestral 7mer out of the train set so taht the denominator isn't known
# not doing any grouping by central 5mer or 3mer though (at least for now)

#### *want* to randomly pick 1/4 of 7mers but have it be the *same* 7mers in both populations
length(unique(allData_withShapes_unprocessed$ancestral7mer)) # just the ancestral 7mers
# pick out random set of 7mers in one species and assign the mutation types to train/test for both species
sampleFracForTraining=0.75 # this will be mulitplied by total number of unique anc 7mers (8192) to get 75%

uniqueAnc7mers <- unique(allData_withShapes_unprocessed$ancestral7mer)

# sample 75% for training, rest will be testing.
trainingAncestral7mers <- sample(uniqueAnc7mers,sampleFracForTraining*length(uniqueAnc7mers))

# now set labels: first label all as test, then a subset as train
allData_withShapes_unprocessed$label <- "TEST" 
allData_withShapes_unprocessed[allData_withShapes_unprocessed$ancestral7mer %in% trainingAncestral7mers,]$label <- "TRAIN"
## remember: the reason you are doing this by ancestral 7mer is so taht all mutation types with the same ancestral 7mer are all eitehr in train or in test -- so you wouldn't split up AAAAAAAA>T and AAAAAAA>G, they would both be either in train and test
# and note here that now all species ahve the same 7mers in train or in test (good)
# make sure same for both pops:
length(allData_withShapes_unprocessed[allData_withShapes_unprocessed$population=="Mmd" & allData_withShapes_unprocessed$label=="TRAIN",]$variable) # 18432 - may change slightly based on seed? I don't think so 

setequal(allData_withShapes_unprocessed[allData_withShapes_unprocessed$population=="Mmd" & allData_withShapes_unprocessed$label=="TRAIN",]$variable , allData_withShapes_unprocessed[allData_withShapes_unprocessed$population=="Ms" & allData_withShapes_unprocessed$label=="TRAIN",]$variable) # checks if are equal 


setequal(allData_withShapes_unprocessed[allData_withShapes_unprocessed$population=="Mmd" & allData_withShapes_unprocessed$label=="TEST",]$variable , allData_withShapes_unprocessed[allData_withShapes_unprocessed$population=="Ms" & allData_withShapes_unprocessed$label=="TEST",]$variable) # checks if are equal 

############# rescale data within TRAIN and TEST and within species and get odds ratios ###########
allData_withShapes_rescaled <- allData_withShapes_unprocessed %>%
  group_by(population,label) %>%
  mutate(mutationRate_scaledWithinTrainTest = mutationRate_unscaled / sum(mutationRate_unscaled))

# make sure that worked:
head(allData_withShapes_rescaled)
sum(allData_withShapes_rescaled[allData_withShapes_rescaled$label=="TEST"& allData_withShapes_rescaled$population=="Ms",]$mutationRate_scaledWithinTrainTest) #1 sums to one within a population and a TEST or TRAIN dataset
sum(allData_withShapes_rescaled[allData_withShapes_rescaled$label=="TRAIN"& allData_withShapes_rescaled$population=="Ms",]$mutationRate_scaledWithinTrainTest) #1 sums to one within a population and a TEST or TRAIN dataset

head(allData_withShapes_rescaled)


# now get the odds ratios and make some plots # 

# need to pivot it 
# want to keep all Pos and feature columns 
# give a generic name to each species (A, B)
allData_withShapes_rescaled$speciesGenericLabel <- ""
allData_withShapes_rescaled[allData_withShapes_rescaled$population==speciesForModel[1],]$speciesGenericLabel <- "SpeciesA"
allData_withShapes_rescaled[allData_withShapes_rescaled$population==speciesForModel[2],]$speciesGenericLabel <- "SpeciesB"

allData_withShapes_rescaled_wide <- pivot_wider(allData_withShapes_rescaled,id_cols = c("variable","label",starts_with("Pos_"),starts_with("feature_"),"speciesGenericLabel"),values_from = "mutationRate_scaledWithinTrainTest",names_from="speciesGenericLabel") # dropping here the metadata cols that I don't need like centralThreemer 
# this pivots so now the columns Mmd and Ms are the column names (this will be non ideal when going to other species, may want to relabel as species A and B )
dim(allData_withShapes_rescaled_wide)
### CALLING OddsRatio_AoverB_mutationRateScaledWithinTrainOrTest = "outcome"
# NOTE: for 3 mutations in Ms there is a mutation rate of 0. this leads to Inf as the OR -- not good! Could add a little epsilon or could exclude
# going to exclude those for now (????) but need to make note of that.

allData_withShapes_rescaled_wide$outcome  <- allData_withShapes_rescaled_wide$SpeciesA / allData_withShapes_rescaled_wide$SpeciesB


######### NEED to get rid of odds ratios of 0 or inf (where one species is totally missing the mutation) ########### NOT SURE ABOUT THIS
allData_withShapes_rescaled_wide_no0NoInf <- allData_withShapes_rescaled_wide[allData_withShapes_rescaled_wide$outcome!=Inf & allData_withShapes_rescaled_wide$outcome!=0,]

dim(allData_withShapes_rescaled_wide)
dim(allData_withShapes_rescaled_wide_no0NoInf)

# missing 3. 
###### do some exploring of OR distrubtions ############
ggplot(allData_withShapes_rescaled_wide_no0NoInf,aes(x=outcome,color=label))+
  geom_density()+
  theme_bw()+
  ggtitle("comparing density of ORs based on rescaled mutation rates in train vs test")


ggplot(allData_withShapes_rescaled_wide_no0NoInf,aes(x=label,y=outcome,color=label))+
  geom_boxplot()+
  theme_bw()+
  ggtitle("comparing dist of ORs based on rescaled mutation rates in train vs test")


# want to compare rescaled mutation rates per species to make sure on are not systematically different:

ggplot(allData_withShapes_rescaled_wide_no0NoInf,aes(x=SpeciesA,y=SpeciesB,color=label))+
  geom_point()+
  theme_bw()+
  ggtitle("comparing rescaled mutation rates in train vs test\ntest is shifted upward due to larger denominator when rescaling, but that should be wiped out when converted to OR")+
  scale_y_log10()+
  scale_x_log10()+
  facet_wrap(~label)+
  geom_abline()

ggplot(allData_withShapes_unprocessed,aes(x=label, y=mutationRate_unscaled,color=population))+
  geom_violin()

ggplot(allData_withShapes_rescaled,aes(x=label, y=mutationRate_scaledWithinTrainTest,color=population))+
  geom_violin()
#### you are here, plotting various distributions####
#not sure about best way to rescale things now... having doubts
# no longer going to have population membership in the model (not a feature anymore maybe that will be weird with SHAP? tbd)
# now to make the split, I am going to split on TRAIN/TEST
indices <-
  list(analysis   = which(allData_withShapes_rescaled_wide_no0NoInf$label=="TRAIN"), 
       assessment = which(allData_withShapes_rescaled_wide_no0NoInf$label=="TEST"))

split <- make_splits(indices,allData_withShapes_rescaled_wide_no0NoInf)
saveRDS(split, file = paste0(outdir,"split.rds"))

# if you want to see what's what:
names(training(split))

head(training(split),4)
head(testing(split),4)
# note will need to get rid of SpeciesA and SPecies B columns and "variable" and "label"
# not going to use cross validation (yet?) -- maybe will use for model tuning. probably won't be amazing because too small? we'll see. 

# not using yet:
train_data_cv <- vfold_cv(training(split),v=5) 
train_data_cv


############ recipe #########

rand_forest_processing_recipe <- 
  # which consists of the formula (outcome ~ predictors) (don't want to include the 'variable' column)
  
  recipe(outcome ~ .,data=training(split)) %>% # 
  update_role(variable, new_role="7mer mutation type label") %>% # this isn't a feature used in the model 
  #step_rm(derived7mer,ancestral7mer, mutationCount_AllChrs,ancestral7merCount_AllChrs,label) %>% # already removed a lot of these so just needs a simpler one
  step_rm(SpeciesA,SpeciesB,label) %>% # removes these columns 
  step_dummy(all_nominal_predictors())

rand_forest_processing_recipe %>% summary()

rand_forest_processing_recipe %>% prep() %>% juice() %>% names()
# juice training data so I can look at it:

######### MODEL SPECIFICATION #########
rand_forest_ranger_model_specs <-
  rand_forest(trees = 1000, mtry = 32, min_n = 1) %>% # I added in tree number = 1000
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


####### LAST FIT ON SPLIT ############ 
last_rf_fit <- 
  rand_forest_workflow %>% 
  last_fit(split)
saveRDS(last_rf_fit, file = paste0(outdir,"FitOfFullModelToHeldOut7mers.rds"))


last_rf_fit$.workflow

predictions <- last_rf_fit %>%
  collect_predictions()

metrics <- last_rf_fit %>%
  collect_metrics() 
metrics
write.table(metrics,paste0(outdir,'metrics.txt'),row.names=F,quote=F,sep="\t")

predictions_withData <- bind_cols(predictions,testing(split)[,c("variable")])
predictions_withData$centralMutationType <- paste0(substr(predictions_withData$variable,4,4),".",substr(predictions_withData$variable,12,12))


saveRDS(predictions_withData, file = paste0(outdir,"ModelPredictions.WithData.rds"))
# read back in if you need it:
#predictions_withData <- readRDS(paste0(outdir,"ModelPredictions.WithData.rds"))
dim(predictions_withData)
head(predictions_withData)

#### add in Cpg Labels ###########

predictions_withData$ancestralDimer <- substr(predictions_withData$variable,4,5)

predictions_withData$CpGLabel <- ""

predictions_withData[predictions_withData$ancestralDimer=="CG",]$CpGLabel <- "CpG_"


predictions_withData$MutationTypePlusCpGLabel <- paste0(predictions_withData$CpGLabel,gsub("\\.",">",predictions_withData$centralMutationType))

head(predictions_withData[,c("variable","MutationTypePlusCpGLabel")])
tail(predictions_withData[,c("variable","MutationTypePlusCpGLabel")])

#### get rsq per type ######
rsqsPerSpeciesAndMutationType <- predictions_withData %>%
  group_by(MutationTypePlusCpGLabel) %>%
  rsq(truth=outcome,estimate=.pred)
rsqsPerSpeciesAndMutationType
write.table(rsqsPerSpeciesAndMutationType,paste0(outdir,"rsq.percentraltype.CpGLabels.txt"),quote = F,row.names=F,sep="\t")


# make a label that has rsq in the central mutation type
predictions_withData_WithPerMutationTypeMetrics <- merge(predictions_withData,rsqsPerSpeciesAndMutationType,by="MutationTypePlusCpGLabel")

head(predictions_withData_WithPerMutationTypeMetrics)
########## combine pop and mutaiton type and rsq into one label to facet over 
predictions_withData_WithPerMutationTypeMetrics$comboLabel <- paste0(predictions_withData_WithPerMutationTypeMetrics$MutationTypePlusCpGLabel,"\nrsq = ",round(predictions_withData_WithPerMutationTypeMetrics$.estimate,2))

head(predictions_withData_WithPerMutationTypeMetrics)

plotWithRsq <-  ggplot(predictions_withData_WithPerMutationTypeMetrics, aes(y=.pred,x=outcome,color=MutationTypePlusCpGLabel))+
  geom_point()+
  geom_abline()+
  #geom_text(data=rsqsPerSpeciesAndMutationType,aes(x=5e-3,y=4e-3,label=round(.estimate,4)),color="black")+
  facet_wrap(~comboLabel,scales="free",ncol=3,dir = "v")+
  #ggtitle(paste0(description,"\ntrained on Fold01\n(all odd Windows but one, tested on just Window",toString(unique(truth_prediction_df$newGroup)),")"))+
  theme_bw() +
  theme(text=element_text(size=14))+
  ggtitle(paste0("trained on 3/4 of 7mers across species, tested on 1/4 of 7mers\nOdds Ratio of ",speciesForModel[1],"/",speciesForModel[2]))+
  theme(legend.position="none")+
  xlab("Odds Ratio (truth)")+
  ylab("Odds Ratio (model prediction)")+
  theme(text=element_text(size=14))
#scale_x_log10()+
#scale_y_log10()
plotWithRsq

ggsave(paste0(outdir,"plotWithRsq.png"),plotWithRsq,height=6,width=9)

############### try some other ways to plot ##############
epsilon=0.1
# categorize ORs as within 1+- eps vs > and <
predictions_withData_WithPerMutationTypeMetrics$empirical_ORCategory <- ""
# label ones that are within epsilon of 1
predictions_withData_WithPerMutationTypeMetrics[predictions_withData_WithPerMutationTypeMetrics$outcome <= (1+epsilon) & predictions_withData_WithPerMutationTypeMetrics$outcome >= (1-epsilon),]$empirical_ORCategory <- paste0(1-epsilon," <= OR <= ",1+epsilon)

predictions_withData_WithPerMutationTypeMetrics[predictions_withData_WithPerMutationTypeMetrics$outcome > 1+epsilon,]$empirical_ORCategory <- paste0("OR > ",1+epsilon)

predictions_withData_WithPerMutationTypeMetrics[predictions_withData_WithPerMutationTypeMetrics$outcome < 1-epsilon,]$empirical_ORCategory <- paste0("OR < ",1-epsilon)

# how often does it categorize correctly?
predictions_withData_WithPerMutationTypeMetrics$model_ORCategory <- ""
# label ones that are within epsilon of 1
predictions_withData_WithPerMutationTypeMetrics[predictions_withData_WithPerMutationTypeMetrics$.pred <= (1+epsilon) & predictions_withData_WithPerMutationTypeMetrics$.pred >= (1-epsilon),]$model_ORCategory <- paste0(1-epsilon," <= OR <= ",1+epsilon)

predictions_withData_WithPerMutationTypeMetrics[predictions_withData_WithPerMutationTypeMetrics$.pred > 1+epsilon,]$model_ORCategory <- paste0("OR > ",1+epsilon)

predictions_withData_WithPerMutationTypeMetrics[predictions_withData_WithPerMutationTypeMetrics$.pred < 1-epsilon,]$model_ORCategory <- paste0("OR < ",1-epsilon)

# how often do categories match?
sum(predictions_withData_WithPerMutationTypeMetrics$model_ORCategory==predictions_withData_WithPerMutationTypeMetrics$empirical_ORCategory)
# 4591
sum(predictions_withData_WithPerMutationTypeMetrics$model_ORCategory!=predictions_withData_WithPerMutationTypeMetrics$empirical_ORCategory)
# 1552 
# so it's doing decently at categorizing them
# how to plot that
predictions_withData_WithPerMutationTypeMetrics$match <- ""
predictions_withData_WithPerMutationTypeMetrics[predictions_withData_WithPerMutationTypeMetrics$empirical_ORCategory==predictions_withData_WithPerMutationTypeMetrics$model_ORCategory,]$match <- "yes"
predictions_withData_WithPerMutationTypeMetrics[predictions_withData_WithPerMutationTypeMetrics$empirical_ORCategory!=predictions_withData_WithPerMutationTypeMetrics$model_ORCategory,]$match <- "no"

# maybe we want the match rate (?)
# per central mutation type? I dunno

ggplot(predictions_withData_WithPerMutationTypeMetrics,aes(x=empirical_ORCategory,y=.pred))+
  geom_violin()
### get a confusion matrix #####
# arrange levels:
predictions_withData_WithPerMutationTypeMetrics$empirical_ORCategory <-  factor(predictions_withData_WithPerMutationTypeMetrics$empirical_ORCategory,levels=c(paste0("OR < ",1-epsilon),paste0(1-epsilon," <= OR <= ",1+epsilon),paste0("OR > ",1+epsilon)))

predictions_withData_WithPerMutationTypeMetrics$model_ORCategory <-  factor(predictions_withData_WithPerMutationTypeMetrics$model_ORCategory,levels=c(paste0("OR < ",1-epsilon),paste0(1-epsilon," <= OR <= ",1+epsilon),paste0("OR > ",1+epsilon)))
                                                                                
confusion_matrix <- conf_mat(predictions_withData_WithPerMutationTypeMetrics,truth="empirical_ORCategory",estimate="model_ORCategory")

confusion_matrix
confusion_matrix_plot <- autoplot(confusion_matrix,type="heatmap")+theme(axis.ticks = element_blank()) 
confusion_matrix_plot
ggsave(paste0(outdir,"confusion_matrix_plot.epsilon.",epsilon,".png"),confusion_matrix_plot,height=7,width=9)

############## VIP #####################
# need it to just be fit rather than last fit? annoying
rf_fit_notlastfit <- fit(rand_forest_workflow,data = training(split))

saveRDS(rf_fit_notlastfit, file = paste0(outdir,"rf_fit_notlastfit.rds"))

ranger_obj <- pull_workflow_fit(rf_fit_notlastfit)$fit
ranger_obj
#OOB prediction error (MSE):       3.130142e-06 
# R squared (OOB):                  0.9729254
vi_scores <- vip::vi(ranger_obj)
write.table(vi_scores,paste0(outdir,"VIPScores.txt"),quote = F,row.names=F,sep="\t")

vip_plot <- vip(vi_scores,include_type = T,num_features = 100)
vip_plot

ggsave(paste0(outdir,"vipplot.png"),vip_plot,height=12,width=5)
#vip(ranger_obj,include_type = T,num_features = 20)

# eventually do shap?

################# Design a couple naive models ##############
### okay I talked to KH and we want to calculate 5mer mutation rates separately for train/test even if they have the same ancestral 5mer (just calculate targets separately as well) --- do some plotting to make sure this makes sense.
### then algorithmically if a 5mer is missing from training that's in test then you need to look for the closest 5mer with overlapping 4mers on either end and use those 
# idea for a naive model : want to collapse down to 5mers and get ORs
allData_withShapes_unprocessed$central5mer <- paste0(substr(allData_withShapes_unprocessed$variable,2,6),".",substr(allData_withShapes_unprocessed$variable,10,14))

tail(allData_withShapes_unprocessed[,c("variable","central5mer")])
head(tail(allData_withShapes_unprocessed[,c("variable","central5mer")])
)

# CHANGING this: to separate by train test. this leads to different target counts in train/test but that is okay because the mutation rate numerators are also different -- it's disjoint (talked this over with KH)
countsOf5mers <- allData_withShapes_unprocessed %>% 
  group_by(central5mer,population,label) %>%
  summarise(total5mers=sum(totalSitesAllIntervals),total5merTargets_inTrainOrTest=sum(targetCount_allIntervals)) # as long as you're grouping by different central 5mer mutations you won't be ouble-counting targets because :
countsOf5mers$mutationRate_5mer <- countsOf5mers$total5mers/countsOf5mers$total5merTargets_inTrainOrTest


# get avg of ORs of 5mers that have the same central 4mer 
head(countsOf5mers)

# CHANGE: rescaling within just train or test
countsOf5mers <- countsOf5mers %>%
  group_by(population,label) %>%
  mutate(mutationRate_5mer_rescaled = mutationRate_5mer / sum(mutationRate_5mer))

head(countsOf5mers)
# make it wide:
countsOf5mers$speciesGenericLabel <- ""
countsOf5mers[countsOf5mers$population==speciesForModel[1],]$speciesGenericLabel <- "SpeciesA"
countsOf5mers[countsOf5mers$population==speciesForModel[2],]$speciesGenericLabel <- "SpeciesB"

# need to separate train/test
countsOf5mers_wide <- pivot_wider(countsOf5mers,id_cols = c("central5mer","speciesGenericLabel"),values_from = "mutationRate_5mer_rescaled",names_from=c("speciesGenericLabel","label")) # does this work?

head(countsOf5mers_wide)
# odds ratio of a /b 
countsOf5mers_wide$OddsRatio_5mer_TRAIN <- countsOf5mers_wide$SpeciesA_TRAIN / countsOf5mers_wide$SpeciesB_TRAIN # use this for naive model predictions

# and for interest:
# odds ratio of a /b 
countsOf5mers_wide$OddsRatio_5mer_TEST <- countsOf5mers_wide$SpeciesA_TEST / countsOf5mers_wide$SpeciesB_TEST # don't use this for predictions tho, this is more like the 'truth' but you don't need it for anything

head(countsOf5mers_wide[,c("central5mer","OddsRatio_5mer_TRAIN","OddsRatio_5mer_TEST")])


########### get central 4mers on left and right ############
countsOf5mers_wide$left4mer <- paste0(substr(countsOf5mers_wide$central5mer,1,4),".",substr(countsOf5mers_wide$central5mer,7,10))

countsOf5mers_wide$right4mer <- paste0(substr(countsOf5mers_wide$central5mer,2,5),".",substr(countsOf5mers_wide$central5mer,8,11))

head(countsOf5mers_wide[,c("central5mer","left4mer","right4mer")])

countsOf5mers_wide$left4mer <- paste0(substr(countsOf5mers_wide$central5mer,2,5),".",substr(countsOf5mers_wide$central5mer,8,11))


########### deal with missing 5mers ###########

# need to find 5mers that are missing from train or test.
all5mersinTRAIN = unique(countsOf5mers[countsOf5mers$label=="TRAIN",]$central5mer)
all5mersinTEST = unique(countsOf5mers[countsOf5mers$label=="TEST",]$central5mer)
# if there's a 5mer that's entirely not present in TRAIN we need to find it's close equivalents and get their combined odds ratio (rescaled? no? get 4mer rates? get average of odds ratio?)

missing5mers <- all5mersinTEST[!(all5mersinTEST %in% all5mersinTRAIN)] # 5mers that are in train not test (or vice versa) -- don't technically need to do something with these since all of test are in train thankfully. but that won't always be the case. so!
missing5mers # okay so no 5mers are missing from train at this point. but if they are
#allMissing5mersdf = data.frame()
countsOf5mers_wide$OddsRatio_basedOn4mers_TRAIN <- 0
# set naive prediction but update it if NAs exist:
countsOf5mers_wide$naivePrediction <- countsOf5mers_wide$OddsRatio_5mer_TRAIN
if(length(missing5mers)>0){
  # okay can't get vector version working so just going to have to loop over missing ones for now I'm stuck
  for(missing5mer in missing5mers){
  #eg if missing5mers=c("AGACG.AGCCG","ACGGG.ACTGG")
  df = data.frame(missing5mer)
  df$left4mer <- paste0(substr(df$missing5mer,1,4),".",substr(df$missing5mer,7,10))
  
  df$right4mer <- paste0(substr(df$missing5mer,2,5),".",substr(df$missing5mer,8,11))
  # merge with data that has left and right 4mers 
  left4merMatch <- merge(countsOf5mers_wide,df,by="left4mer")
  # gather up list of ORs for averaging per missing 5mer's left 4mer
  left4merMatch_GatherORs <- left4merMatch %>%
    group_by(left4mer) %>%
    summarise(left4merORList=list(OddsRatio_5mer_TRAIN))
  left4merMatch_GatherORs
  right4merMatch <- merge(countsOf5mers_wide,df,by="right4mer")
  right4merMatch_GatherORs <- right4merMatch %>%
    group_by(right4mer) %>%
    summarise(right4merORList=list(OddsRatio_5mer_TRAIN))
 
  # put back into df
  # combine above:
  mean = mean(c(left4merMatch$OddsRatio_5mer_TRAIN,right4merMatch$OddsRatio_5mer_TRAIN),na.rm = T)
  # each produces 4 hits, one+ of which will be NA in the TRAIN dataset 
  # so want to average over the other hits, omitting NA
  # for now just averaging over Odds Ratio -- that's not mathematically perfect (?) but I think close to what the RF is doing since it doesn't see hte underlying raw data. Maybe more proper would be to totally recalculate mutation rates for 4mers but it's not clear how to properly rescale
  meanOf4merORs = mean(c(left4merMatch$OddsRatio_5mer_TRAIN,right4merMatch$OddsRatio_5mer_TRAIN),na.rm = T) # okay this gets the average of the odds ratios of all ; I think this is okay. if this becomes critical go back through this.
  #df$meanOf4merORs <- meanOf4merORs 
  #allMissing5mersdf <- rbind(allMissing5mersdf,df)
  
  countsOf5mers_wide[countsOf5mers_wide$central5mer==missing5mer,]$OddsRatio_basedOn4mers_TRAIN <- meanOf4merORs # update prediction with this avg
  }
  # update all missing 5mers:
  countsOf5mers_wide[is.na(countsOf5mers_wide$OddsRatio_5mer_TRAIN),]$naivePrediction <- countsOf5mers_wide[is.na(countsOf5mers_wide$OddsRatio_5mer_TRAIN),]$OddsRatio_basedOn4mers_TRAIN
# so if you need this you can fill stuff in   
}



ggplot(countsOf5mers_wide,aes(x=naivePrediction,y=OddsRatio_5mer_TEST))+
  geom_point()+
  geom_abline()+
  ggtitle("got odds ratios of same 5mers based on 7mers contained in eithe train/test -- this is how similar those odds ratios are. do worse at high ORs when data is missing")

# okay now want to predict the test set by the odds ratio of the full spectrum (which includes them) -- some data leakage but want this is the naive model so ok if it cheats to some degree I think -- since it's the model to beat.

predictions_withData$central5mer <- paste0(substr(predictions_withData$variable,2,6),".",substr(predictions_withData$variable,10,14))
head(predictions_withData)


# so if an ancestral 5mer from predictions_withData isn't present in the TRAIN dataset
# what do you do? (a, need to discover it)
# (b, need to find a next-closest)
# so for example AAACCCC>T ; want to look for AAACCX>T and XAACCC>T and get their odds ratios from train? or not their odds ratios but their actual counts? 



# want to merge these with 5mer spectrum
head(predictions_withData)

### so here need to deal with what happens for missing 5mers 
naiveModel_5merPredictions <- merge(predictions_withData[,c(".pred","outcome","variable","central5mer","centralMutationType")],countsOf5mers_wide[,c('central5mer','naivePrediction')],by="central5mer")

colnames(naiveModel_5merPredictions) <- c("central5mer","rf.pred","OR_7merempirical","MutationType_7mer","centralMutationType","naivePrediction")

 
# so it's just predicting the odds ratio based on the odds ratio of the central 5mer based on data in train (with the paricular 7mer as part)

# note: are any rescaling effects at play here? maybe? 


#### get rsq per type ######

naiveModel_rsq <- naiveModel_5merPredictions %>%
  group_by(centralMutationType) %>%
  rsq(truth=OR_7merempirical,estimate=naivePrediction)
naiveModel_rsq # where are NAs comign from?

# make a label that has rsq in the central mutation type
naiveModel_5merPredictions_WithPerMutationTypeRsq <- merge(naiveModel_5merPredictions,naiveModel_rsq,by=c("centralMutationType"))


########## combine pop and mutaiton type and rsq into one label to facet over 
naiveModel_5merPredictions_WithPerMutationTypeRsq$comboLabel <- paste0(naiveModel_5merPredictions_WithPerMutationTypeRsq$centralMutationType,"\nnaive model rsq = ",round(naiveModel_5merPredictions_WithPerMutationTypeRsq$.estimate,3))

naiveplot1 <- ggplot(naiveModel_5merPredictions_WithPerMutationTypeRsq,aes(x=OR_7merempirical,y=naivePrediction,color="naive 5mer model"))+
  geom_point()+
  geom_abline(slope = 1,intercept = 0)+
  theme_bw()+
  geom_point(aes(x=OR_7merempirical,y=rf.pred,color="rf model"))+
  facet_wrap(~comboLabel,scales="free",ncol=2,dir = "v")+
  #ggtitle(paste0(description,"\ntrained on Fold01\n(all odd Windows but one, tested on just Window",toString(unique(truth_prediction_df$newGroup)),")"))+
  theme_bw() +
  theme(text=element_text(size=14))+
  ggtitle(paste0("naive model is 5mer OR across all 7mers, predicting every sub-7mer to have the same odds ratio. \npoints are just the 'test' 7mers from the rf trained model, but note that in the naive model these were not held out - so naive model has unfair advantage.")) 

naiveplot1
ggsave(paste0(outdir,"naive5merModel.vs.rfmodel.png"),naiveplot1,height=12,width=9)


write.table(naiveModel_rsq,paste0(outdir,"NAIVEMODEL.rsq.percentraltype.txt"),quote = F,row.names=F,sep="\t")

# compare:
allrsqs <- merge(naiveModel_rsq,rsqsPerSpeciesAndMutationType,by="centralMutationType",suffixes = c(".naive",".rf"))
head(allrsqs)
allrsqs_melt <- melt(allrsqs[,c("centralMutationType",".estimate.naive",".estimate.rf")])

rsqsNaiveRFComparisonPlot <- ggplot(allrsqs_melt,aes(x=centralMutationType,y=value,fill=variable))+
  geom_col(position="dodge")+
  theme_bw()
rsqsNaiveRFComparisonPlot
ggsave(paste0(outdir,"naive.vs.rf.rsqs.png"),rsqsNaiveRFComparisonPlot,height=6,width=8)


# make a label that has rsq in the central mutation type
