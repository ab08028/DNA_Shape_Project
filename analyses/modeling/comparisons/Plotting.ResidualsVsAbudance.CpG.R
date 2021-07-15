############ trying to stratify by abundance to see what CpGs are doing ############
########## variability of CpG containing 7mers ###########
allData_multipop <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_spectrum_target_merged_PerChr/TEMPORARY.NEWGROUPS.MULTIPOPULATION_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.INCLUDES0Entries.txt",header=T) ####### NOTE make sure this is same as model that you ran! needs to be correct grouping; previously was using per-chr but this is per larger group

# need to remove NAs
allData_multipop <- na.omit(allData_multipop)
# rescale by population/wdinow:
allData_multipop <- allData_multipop %>%
  group_by(population,newGroup,label) %>%
  mutate(outcome=mutationCount_divByTargetCount/sum(mutationCount_divByTargetCount)) 

# label according to CpG counts
allData_multipop$centralMutationType <- paste0(substr(allData_multipop$mutationType,4,4),".",substr(allData_multipop$mutationType,12,12))

allData_multipop$ancestral_central3mer <- substr(allData_multipop$mutationType,start=3,stop=5)
allData_multipop$centralCpGLabel <- "NoCentralCpG"
allData_multipop[grep('CG',allData_multipop$ancestral_central3mer),]$centralCpGLabel <- "YesCentralCpG"

allData_multipop$centralMutationType <- paste0(substr(allData_multipop$mutationType,4,4),".",substr(allData_multipop$mutationType,12,12))

# make labels: 
allData_multipop$mutationLabel <- paste0(allData_multipop$centralMutationType,"_",allData_multipop$centralCpGLabel)

allData_multipop$ancestral_CpG_Count <- str_count(allData_multipop$ancestral7mer,"CG")

##### restrict to similar abundances #######
ggplot(allData_multipop[allData_multipop$population=="Mmd",],aes(x=ancestral7merCount,color=as.factor(ancestral_CpG_Count)))+
  geom_density()+
  geom_vline(xintercept = 2.5e4)+
  geom_vline(xintercept = 1e4)+
  scale_x_log10()

# restrict to window 7 which tests were run on for model 12 (truth_pred_df below)
restricted <- allData_multipop[allData_multipop$ancestral7merCount < 2.5e4 & allData_multipop$ancestral7merCount > 1e4 & allData_multipop$newGroup==7,] # trying to keep similar abundances 
dim(restricted)
head(restricted)

########## model 12 predictions ##########

truth_prediction_df <- readRDS("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210704_model_012_RF.plusMethylInfo/modelTrainedOnOneFold.PREDICTIONS.onWindow7.rds")
head(truth_prediction_df)

truth_prediction_df$centralMutationType <- paste0(substr(truth_prediction_df$mutationType,4,4),".",substr(truth_prediction_df$mutationType,12,12))


truth_prediction_df$ancestral_central3mer <- substr(truth_prediction_df$mutationType,start=3,stop=5)
truth_prediction_df$centralCpGLabel <- "NoCentralCpG"
truth_prediction_df[grep('CG',truth_prediction_df$ancestral_central3mer),]$centralCpGLabel <- "YesCentralCpG"

truth_prediction_df$centralMutationType <- paste0(substr(truth_prediction_df$mutationType,4,4),".",substr(truth_prediction_df$mutationType,12,12))

# make labels: 
truth_prediction_df$mutationLabel <- paste0(truth_prediction_df$centralMutationType,"_",truth_prediction_df$centralCpGLabel)

truth_prediction_df$ancestral_CpG_Count <- str_count(truth_prediction_df$ancestral7mer,"CG")


truth_prediction_df$residual <- truth_prediction_df$.pred - truth_prediction_df$outcome


### plot residuals vs 7mer abundance #########
ggplot(truth_prediction_df,aes(y=residual,x=ancestral7merCount,color=as.factor(ancestral_CpG_Count)))+
  geom_point()+
  facet_wrap(~mutationLabel,scales="free")+
  scale_x_log10()+
  geom_hline(yintercept = 0)

ggplot(truth_prediction_df,aes(y=residual,x=ancestral7merCount,color=as.factor(ancestral_CpG_Count)))+
  geom_point()+
  scale_x_log10()+
  geom_hline(yintercept = 0)

ggplot(truth_prediction_df,aes(y=residual,x=mutationCount,color=as.factor(ancestral_CpG_Count)))+
  geom_point()+
  facet_wrap(~mutationLabel,scales="free")
  #scale_x_log10()


# try plotting mutation rate vs abundance
ggplot(truth_prediction_df,aes(y=outcome,x=mutationCount,color=as.factor(ancestral_CpG_Count)))+
  geom_point()+
  geom_point(data=truth_prediction_df,aes(x=mutationCount,y=.pred),color="gray")+
  facet_wrap(~mutationLabel,scales="free")
#scale_x_log10()
############### restrict to just those similar abundance types ########
truth_prediction_df_restricted <- truth_prediction_df[truth_prediction_df$mutationType %in% restricted$mutationType,]

dim(truth_prediction_df_restricted)
head(truth_prediction_df_restricted)

ggplot(truth_prediction_df_restricted,aes(y=residual,x=ancestral7merCount,color=as.factor(ancestral_CpG_Count)))+
  geom_point()+
  facet_wrap(~mutationLabel,scales="free")+
  scale_x_log10()+
  geom_hline(yintercept = 0)+
  scale_color_manual(values=c("tomato","#7CAE00"))
########## try getting poisson LLs ###########
# LL from dpois
head(truth_prediction_df)

sum(truth_prediction_df$mutationCount/truth_prediction_df$ancestral7merCount) # this was denominator for scaling rates in observed data 
truth_prediction_df <- truth_prediction_df %>%
  group_by(population,newGroup) %>%
  mutate(empirical_scalingFactor=sum(mutationCount/ancestral7merCount))
truth_prediction_df # what rate was divided by

truth_prediction_df <- truth_prediction_df %>%
  group_by(population,newGroup) %>%
  mutate(.predicted_count_UsingObsScalingFactor_TEMPORARY = .pred*empirical_scalingFactor*ancestral7merCount)

head(truth_prediction_df[,c("mutationCount",".predicted_count_UsingObsScalingFactor_TEMPORARY",".pred","outcome")]) # off by order of mag because of how I rescaled rates (rate = rate/sum(rates)) so I need to undo that somehow hm.
# okay this works okay for now but eventually need to get the scaling factor from the modelfor the rates? 
truth_prediction_df <- truth_prediction_df %>% 
  mutate(poissonLLs_stillUsingEmpiricalScalingFactorForNow = dpois(mutationCount,lambda = .predicted_count_UsingObsScalingFactor_TEMPORARY,log=T)) # from kelley: mean should be predicted rate * ancestral 7mer count (but I need to add in scaling factor which this contains)

head(truth_prediction_df$poissonLLs_stillUsingEmpiricalScalingFactorForNow)

ggplot(truth_prediction_df,aes(x=mutationLabel,y=poissonLLs_stillUsingEmpiricalScalingFactorForNow,fill=mutationLabel))+
  geom_boxplot()
# so the problem here is that I actually modeled the rescaled rate (rate/sum(rates)) and that's what model puts out, so it's not quite as simple as multiplying predret * ancestral7mer count. also need some scaling factor (that's ~100) to scale predictions totally up to counts . for now I'm just using the empirical per-species per-window scaling factor
# did I calculate poisson LLs right? I think so 
# see how this related to residuals

View(truth_prediction_df[,c("poissonLLs_stillUsingEmpiricalScalingFactorForNow","mutationCount",".predicted_count_UsingObsScalingFactor_TEMPORARY","mutationType","mutationLabel")])
ggplot(truth_prediction_df,aes(x=abs(residual),y=poissonLLs_stillUsingEmpiricalScalingFactorForNow))+
  geom_point()+
  facet_wrap(~mutationLabel,scales="free")+
  geom_smooth(method="lm")

########### want to plot central 5mer average abundance (?) and see how 7mer predictions relate ###########
