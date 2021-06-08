############ trying out caret on mouse data ##########################
######## TO USE GLMNET YOU NEED TO BE IN THE ANACONDA R ENVIRONMENT! ###############
# don't try to install within R studio -- actually go to the anaconda navigator and install from there
# copied from other script, has nothing to do with caret #######################

#install.packages("BiocManager") --- in
# dont do it this way BiocManager::install("DNAshapeR") -- install on anaconda navigator instead (not in R studio)
require(DNAshapeR)
require(ggplot2) # done
require(ggfortify) # for pca autoplot - done
require(Biostrings) # needed to go with DNAShapeR - done
require(dplyr) # done
require(reshape2) # done
require(phylotools) # to get names of sequences from fasta  # done
require(caret) # for linear modeling # done
require(broom) # for linear modeling tidying (tidy, augment to combine with original data, glance) # done
require(GGally) # done 
require(vip) # for variable importance
# don't need library(devtools)
#### install on anaconda navigator instead of here install_version("glmnet", version = "3.0") # >3.0 is not available for R 3.5 (needs R 3.6) so specify version 3 here
require(glmnet) # done
### starting with mouse as a case study. Want to make a fasta in R so no bash is needed.  
### meaning of 'Rounding' sampler message:
#R 3.6.0 changes the method used to generate random integers in the sample
#function. In prior versions, the probability of generating each integer could vary from equal by up to 0.04% (or possibly more if generating more than a million different integers). This change will mainly be relevant to people who do large-scale simulations, and it also means that scripts using the sample
#function will generate different results in R 3.6.0 than they did in prior versions of R. If you need to keep the results the same (for reproducibility or for automated testing), you can revert to the old behavior by adding RNGkind(sample.kind="Rounding"))
#to the top of your script. https://www.r-bloggers.com/2019/05/whats-new-in-r-3-6-0/

############# >>> set some vars ##############
correlationThreshold=0.75
CVNumber=8


######### >>> set some dirs #############

fastadir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/fastasOf7mers/"
dir.create(fastadir,recursive = T,showWarnings = F)

################################### THIS SECTION GETS YOUR FEATURE MATRIX READY ###########################################
# this is identical to the lead in for the pca plotting script too
############ >>> Read in the 7mer spectrum data (not targets yet) and add some metadata ########################## 

# NOTE THAT THIS SPECTRUM HAS MULTIPLE MOUSE SPECIES CONTAINED IN IT (that's ok for when we are making fastas because we are doing "unique")
spectrum_multipop <- read.table('/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/summaryTables/Mouse7merSpectrum.SummedUpOverAllChrs.NOSTRICT.txt',header=T) # don't need to do this for every species.

###### IMPORTANT: RESTRICT TO ONE SPECIES/POPULATION HERE!!!! #### 
population="Mmd"
spectrum_onepop <- spectrum_multipop[spectrum_multipop$population==population,] # to prevent accidents

if(length(unique(spectrum_onepop$population))>1){
  print("More than one population!!!! Stop here!")
}
#just making a lookup table. unless a species has one that mouse doesn't.
head(spectrum_onepop)
speciesLabel="mouse"
plotdir=paste0("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/plots/",speciesLabel,"/modeling_plots/")
dir.create(plotdir,recursive = T,showWarnings = F)

spectrum_onepop$ancestralkmer <- unlist(lapply(strsplit(spectrum_onepop$variable,"\\."),"[",1))

spectrum_onepop$derivedkmer <- unlist(lapply(strsplit(spectrum_onepop$variable,"\\."),"[",2))
# add in central 5mer:

spectrum_onepop$central5mer.ancestral <- substr(spectrum_onepop$ancestralkmer,start=2,stop=6)


############ >>> label central mutation type by whether it is CpG or not and what nucleotides are involved (ala L&S) #########
head(spectrum_onepop)
spectrum_onepop$ancestral_central3mer <- unlist(lapply(strsplit(spectrum_onepop$ThreemerMutationType,"\\."),"[",1))
spectrum_onepop$derived_central3mer <- unlist(lapply(strsplit(spectrum_onepop$ThreemerMutationType,"\\."),"[",2))

# identify presence of CpG (CG) [not GC] in central 3mer
spectrum_onepop$ancestral_central3merCpG <- "no"
spectrum_onepop$derived_central3merCpG <- "no"
spectrum_onepop[grepl("CG",spectrum_onepop$ancestral_central3mer),]$ancestral_central3merCpG <- "yes"
spectrum_onepop[grepl("CG",spectrum_onepop$derived_central3mer),]$derived_central3merCpG <- "yes"

spectrum_onepop$mutationClassLabel <- paste0(spectrum_onepop$centralMutationType,".",spectrum_onepop$ancestral_central3mer,"_ancestralCpG") 
unique(spectrum_onepop$mutationClassLabel) # 9 mutation types (same as L&S)
# "A.C.no_ancestralCpG"  "A.G.no_ancestralCpG"  "A.T.no_ancestralCpG"  "C.A.no_ancestralCpG"  "C.G.no_ancestralCpG" 
# "C.T.no_ancestralCpG"  "C.A.yes_ancestralCpG" "C.G.yes_ancestralCpG" "C.T.yes_ancestralCpG"
############ >>> Make a fasta containing all ancestral and derived 7mers. These should be applicable to every species, though some species may be missing a certain 7mer (have a 0 count), but mouse probably has them all - we will see. ################


allUniqueAncestralPlusDerivedkmers_Unique <- unique(c(spectrum_onepop$derivedkmer,spectrum_onepop$ancestralkmer)) # aha so all the ancestral kmers  are contained within the derived kmer list so could actually just use the derived kmer list. doign this to be exra thorough. End up with 16384 kmers 
# make fasta -- just making one with all (unique) 7mers. Don't need anc/derived separated because just going to use this as a lookup table for my mutation spectrum df. But do need them both to be present because derived aren't collapsed to A.X and C.X the way ancestral targets are. # derived is 2x bigger because central mutation type in ancestral 7mers is collapsed to A>X and C>X whereas derived can be ATCG.

fastaContainingAllKmers=paste0(fastadir,"allAncestralPlusDerived7mers.unique.basedOnMouseSpectrum.ShouldWorkForAllSpp.fasta")

sink(fastaContainingAllKmers)
# using unique to not do repeats (there are repeats because of X>Y and X>Z mutations. just want to get shape features of each one once)
for(kmer in allUniqueAncestralPlusDerivedkmers_Unique){
  cat(paste0(">",kmer,"\n"))
  cat(paste0(kmer,"\n"))
}
sink()

# note if you do wc -l on fasta files that they are 2x the number of 7mers because each entry is:
# >name of 7mer
# 7mer 

################ >>> get sequence names from the fasta files ##########################
# get sequence names in order from the fasta: 
kmer_SeqNames_fromFasta=get.fasta.name(fastaContainingAllKmers, clean_name = TRUE)


############## >>> prepare dna shape R parameters ###########################
all14ShapeTypes=c("HelT", "Rise", "Roll", "Shift", "Slide", "Tilt", "Buckle", "Opening", "ProT", "Shear", "Stagger", "Stretch", "MGW", "EP") # default is just HelT proT EP MGW and Roll -- need to manually list all these others


############## >>> make DNAShapeR shape prediction for each fasta sequence. ############

DNAShapeR_Prediction <- getShape(fastaContainingAllKmers,parse=T,shapeType = all14ShapeTypes) # runs fast (a minute or two)

# this yields a massive list of dataframes that has all entries for every 7mer in each fasta. Note that it is done in sliding 5mers so you only get values for the center 3bp and their flanking steps which are at the center of each 5mer. but those data are informed by the whole sequence. Some measurements are per bp and some are per bp 'step' between bps depending on whether its an intra or inter-bp measurement. 
# you can access particular dataframes using :
# FOR EXAMPLE:
head(DNAShapeR_Prediction$HelT) # this would give you the HelT measurements between bp steps
head(DNAShapeR_Prediction$ProT) # this would give you Propeller twist measurements per bp (central 3 bp only)
# you can also plot measurements using:
# FOR EXAMPLE:
plotShape(DNAShapeR_Prediction$MGW)
plotShape(DNAShapeR_Prediction$ProT)
######## NOTE! At this stage the measurements are in units (eg degrees, angstroms), but when they go into your feature 
# matrix they must be normalized so that one with big units doesn't dominate all the others. DNAShapeR helpfully does this for you below

####### >>> Choose the features you want to include in your feature matrix  -- this step will probably be most modified!#############

# choose the variables you want. You could do first order, and/or second order shape features. You could also add sequence features (encoding the nucleotides) and lots of other things that you want in your model. 

firstOrder_featureTypes <- paste0("1-",all14ShapeTypes) # get names of all 14 first order feature types
secondOrder_featureTypes <- paste0("2-",all14ShapeTypes)
sequence_featureTypes_forNucleotideModels= c("1-mer","2-mer","3-mer") # NOT USING THESE FOR NOW BUT CAN INCLUDE IN FUTURE IF I WANT TO DO SEQUENCE+ORDER FEATURES

########### MODIFY DESIRED FEATURES HERE: ###################
desiredFeatures=c(firstOrder_featureTypes) # secondOrder_featureTypes) # for now going to use first order only. then can try  and second order
orderLabel="firstOrder_featureTypes" # ADJUST THIS LABEL IF YOU CHANGE THE ABOVE SETTING ^^^^ 
# the sequence features are very slow
# if you are using second order features toggle this label:

####### >>> Put together a feature dataframe with labelled columns containing all your normalized measurements #############
# note instead of just selecting all at once, I am looping over them so I can label them. Otherwise you lose the labels which makes things less informative downstream!
# this includes all ancestral and derived 
allFeatureVectors_labelled <- data.frame(motif=kmer_SeqNames_fromFasta) # set up with seqNames -- these *must* be the same as the order of sequences from the fasta.
for(feature in desiredFeatures){
  feature_noDash=gsub("-","_",feature) # need to get rid of dash in the name for the columns otherwise R hates it.
  featVec <- data.frame(encodeSeqShape(fastaContainingAllKmers,DNAShapeR_Prediction,feature, normalize=T))
  # label columns based on how many columsn there are for each feature: (adding _1 _2 _3 to the feature name)
  featVecNum=as.numeric(unlist(dim(featVec)[2]))
  # note that for sequence vectors they encode each bp as four columns (1000 0100 etc). So this numbering system doesn't work great for that. so each 7mer gets 28 columns for the 1mer column. Want to rename them somehow feature_1_mer_1 _2 etc is okay I guess. 
  colnames(featVec) <- unlist(paste0("feature_",feature_noDash,"_",as.character(seq(1,featVecNum))))
  
  # then want to bind columns:
  allFeatureVectors_labelled <- cbind(allFeatureVectors_labelled,featVec)
}
# so a first order feature will be: feature_1_ and a second order will be feature_2. 
head(allFeatureVectors_labelled)


###### WRITE THIS OUT TO USE AS FEATURE TABLE ###########
# note this may contain multiple entries because of multiple populations (eg many mice spp)
write.table(mutationSpectrum_featureMatrix,paste0("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/","FEATURE.TABLE.",orderLabel,".ancestral.derived.feautres.AllShapes.ShouldBeGoodForAllSpecies.txt"),row.names=T,quote=F)

########### >>> Make a heatmap of all correlations of the features #################
# https://briatte.github.io/ggcorr/
correlationOfFeaturesAcrossAll7mers_plot <- ggcorr(subset(allFeatureVectors_labelled,select=-c(motif)),layout.exp = 20,hjust=1,vjust=.5,size=1.5)
correlationOfFeaturesAcrossAll7mers_plot
ggsave(paste0(plotdir,"correlationOfFeaturesAcrossAll7mers.",orderLabel,".png"),correlationOfFeaturesAcrossAll7mers_plot,height=8,width=9)

#################### >>> Make data frame feature matrix for both ancestral and derived motifs based on species' mutation spectrum file (note may be multiple pops in spectrum) #################

# NOTE THIS MAY INCLUDE MULTIPLE POPULATIONS
# get the ancestral features into an intermediate file
# make  sure you restricted spectrum to just one pop! 

intermediateMerge_ancestral <- merge(spectrum_onepop,allFeatureVectors_labelled,by.x="ancestralkmer",by.y="motif")
# get the derived features (and set suffixes as ancestral and .derived)

mutationSpectrum_featureMatrix <- merge(intermediateMerge_ancestral,allFeatureVectors_labelled,by.x="derivedkmer",by.y="motif",suffixes=c("_ancestral","_derived")) # sweet, I checked that this works and it does! 


###################### >>> ACTUAL CARET STUFF STARTS HERE #######################


# okay working froma couple tutorials
# https://bioconductor.org/packages/release/bioc/vignettes/DNAshapeR/inst/doc/DNAshapeR.html
# https://bioconductor.org/packages/release/bioc/vignettes/DNAshapeR/inst/doc/DNAshapeR.html

############# >>> make a df just of the predictors and dependent variable: #############
columnsToExclude=names(spectrum_onepop)
columnsToExclude=columnsToExclude[columnsToExclude!="fractionOfAllSegSites"] # want to KEEP frac of seg sites in so this is how I'm doing that for now. 
# but want to keep in the frac of seg sites

# get rid of columns that aren't predictors or the dependent variable


featuredf <- select(mutationSpectrum_featureMatrix,-c(columnsToExclude))
# making rownames the mutation type:
rownames(featuredf) <- mutationSpectrum_featureMatrix$variable
# check if list of names is good:
names(featuredf) # okay now this has frac of all seg sites and all the features, but nothing else (can add back in from mutationSpectrum_featureMatrix as-needed)

################## >>> check for features that have near zero variance #######################
# remove the columns that are not my predictors or dependent variable (?)
# The function preProcess doesn’t actually pre-process the data.  predict.preProcess is used to pre-process. 
# this results in a data.frame with two metrics: freqRatio which indicates the ratio of 
# - the frequency of the most prevalent value over the second most frequent value (called the “frequency ratio’’), which would be near one for well-behaved predictors and very large for highly-unbalanced data 
# so if the most common ROLL-1 value was 0.8 with 500 mutation types having that value and then the next most frequent observation was 0.5 with only 2 mutation types, and the rest are 1, this would be very unbalanced and the freq ratio would be 500/2 
# percentUnique is the percent of the values that are totally unique (not shared by any other mutation types). so "variable" which is the mutation types is 100% unique. 
nzv_withInfo <- nearZeroVar(featuredf, saveMetrics= TRUE) # what is threshold?    freqCut = 95/5, uniqueCut = 10, <- cutoffs default. ok.
nzv_withInfo
# to just get a list to exclude, run it again without saveMetric , then can just subtract because it'll be column 
nzvToExclude <- rownames(nzv_withInfo[nzv_withInfo$nzv==T,])
cat(paste0("non-zero variance predictors to exclude: ",nzvToExclude))
# the ones with "TRUE" are the ones with nzv so the ones you want to exclude
# none of my variables are "TRUE" so none have near zero variance so we are good to go.
featuredf_nzvFilt <- select(featuredf,-all_of(nzvToExclude)) # this doesn't exclude anything because no nzv variables. 

########### >>> look for correlated predictors and remove ###############
# set corr threshold up top
# get correlation matix: 
predCor <-  cor(featuredf_nzvFilt)
# This function searches through a correlation matrix and returns a vector of integers corresponding to columns to remove to reduce pair-wise correlations. 
#The absolute values of pair-wise correlations are considered. If two variables have a high correlation, the function looks at the mean absolute correlation of each variable and removes the variable with the largest mean absolute correlation. 
# doesn't remove both of a pair? 
# exact tells it to recalculate avg correlation after each removal (slows it down but may remove fewer features)
highCorlist <- findCorrelation(predCor,cutoff=correlationThreshold,names = T,exact = T) # gives list to exclude, it chooses one of the two variables to remove based on how correlated it is with all the other variables. recalcs avg cor after each removal with exact=T
cat(paste0("correlation threshold: ", correlationThreshold,".\n\nhighly correlated predictors chosen to be excluded (note: one of pair gets chosen to be removed based on avg corr with other variables; avg corr recalculated after every removal):\n\n",toString(highCorlist)))

featuredf_nzvFilt_corFilt <- select(featuredf_nzvFilt,-all_of(highCorlist))
dim(featuredf_nzvFilt)
dim(featuredf_nzvFilt_corFilt)

########### >>> experiment with lm linear model ####################
# set up controls for the model; cross-validation with 8-fold validation
trainControl <- trainControl(method = "cv", number = CVNumber, 
                             savePredictions = TRUE)
lm_model <- train (fractionOfAllSegSites~ ., data = featuredf_nzvFilt_corFilt, 
                trControl=trainControl, method="lm", preProcess=NULL) # could do some preprocessing here, but I did it already above

lm_model

# useful things:
 # model $ has lots of different interesting things. most imporant is finalModel 
lm_model$resample # gives resampled information of Rsquared.

# final model (best fit out of the CVs?)
lm_finalModel <- lm_model$finalModel # the coefficients and intercept! 
summary(lm_finalModel)
# finalModel has lots of its own columns:
lm_finalModel$coefficients
# residuals (distance between estimate of what frac of seg sites should be and what it actually is)
lm_finalModel$residuals
lm_finalModel$fitted.values
# lots of other goodies in here, not sure what to do with them all!

# broom works on these ouputs! glance tells you fit summary and tidy gives you a nice dataframe that is easy to plot  
lm_modelFit <- glance(lm_finalModel)
lm_modelFit
lm_coefs <- tidy(lm_finalModel)
lm_coefs
# can use augment to combine it with the original data 
lm_coefPlot1 <- ggplot(lm_coefs,aes(y=term,x=estimate,color=-log10(p.value)))+
  geom_point()+
  theme_bw()+
  ggtitle(paste0(lm_model$modelInfo$label," coefficients\n",CVNumber,"-fold Cross Validation\nR^2 = ",round(modelFit$r.squared,4)))+
  scale_color_viridis_c(option="A")
lm_coefPlot1

ggsave(paste0(plotdir,"lm_allMutationClasses",orderLabel,".coefficients.png"),lm_coefPlot1,height=8,width=9)
# eventually write out model fit and coeffs somewhere/how


# variable importance using vip
vipplot <- vip(lm_finalModel,num_features = 100,include_type = T,geom = "point")+theme_bw()
vipplot
ggsave(paste0(plotdir,"correlationOfFeaturesAcrossAll7mers.",orderLabel,".png"),vipplot,height=8,width=9)
# note that t stat is coeff/std error. so if a coeff is much bigger than its std error then it probably doesn't equal zero --> big t stat -> low p value (highly significant)
# plot predictions and obs?

head(lm_finalModel$fitted.values)

# base R does some nice easy plots:
png(paste0(plotdir,"lm.plots.png"))
par(mfrow = c(2, 2))
plot(lm_finalModel)
dev.off()
par(mfrow=c(1,1))



################ >>> experiment with lasso #####################
# https://www.datacareer.ch/blog/ridge-and-lasso-in-r/
# "The “glmnet” method in caret has an alpha argument that determines what type of model is fit. If alpha = 0 then a ridge regression model is fit, and if alpha = 1 then a lasso model is fit. Here we’ll use caret as a wrapper for glment."
# "for now we won’t set the regularization parameter lambda. It is set to 1."
trainControl <- trainControl(method = "cv", number = CVNumber, 
                             savePredictions = TRUE)
# set alpha to  1 
# this requires glmnet
lasso_glm <- train(fractionOfAllSegSites~ ., data = featuredf_nzvFilt_corFilt, method = 'glmnet', 
                            tuneGrid = expand.grid(alpha = 1, lambda = 1), trainControl=trainControl) # should I choose tuning parameters ahead of time or not? lasso? hm.
summary(lasso_glm)                            
lasso_finalModel <- lasso_glm$finalModel
glance(lasso_finalModel)
tidy(lasso_finalModel)
# ehy does best model only have intercept -- gotta spend more time here
plot(lasso_finalModel,xvar="lambda") # this plots the coeffs https://glmnet.stanford.edu/articles/glmnet.html 
#  It shows the path of its coefficient against the ℓ1-norm of the whole coefficient vector as 𝜆 varies.  ??? what does this mean

# have to get coeffs at a value of 'lambda' - what is lambda??


############ >>> lm for just one mutaiton class at a time ##############
all_lmModelsPerClass_modelFit <- data.frame()
all_lmModelsPerClass_modelCoefficients <- data.frame()
for(mutationClass in unique(mutationSpectrum_featureMatrix$mutationClassLabel)){
   # get mutation class
    mutationSpectrum_featureMatrix_oneClassOnly <- mutationSpectrum_featureMatrix[mutationSpectrum_featureMatrix$mutationClassLabel==mutationClass,]
    # get numerical columns only plus frac seg sites
  featuredf <- select(mutationSpectrum_featureMatrix_oneClassOnly,-c(columnsToExclude))
  # set rownames
  rownames(featuredf) <- mutationSpectrum_featureMatrix_oneClassOnly$variable
  # check for near zero var
  nzv_withInfo <- nearZeroVar(featuredf, saveMetrics= TRUE) # what is threshold?    freqCut = 95/5, uniqueCut = 10, <- cutoffs default. ok.
  nzvToExclude <- rownames(nzv_withInfo[nzv_withInfo$nzv==T,])
  featuredf_nzvFilt <- select(featuredf,-all_of(nzvToExclude)) # this doesn't exclude anything because no nzv variables. 
  # check for inter correlation 
  predCor <-  cor(featuredf_nzvFilt)
  highCorlist <- findCorrelation(predCor,cutoff=correlationThreshold,names = T,exact = T) # gives list to exclude, it chooses one of the two variables to remove based on how correlated it is with all the other variables. recalcs avg cor after each removal with exact=T
  featuredf_nzvFilt_corFilt <- select(featuredf_nzvFilt,-all_of(highCorlist))
  trainControl <- trainControl(method = "cv", number = CVNumber, 
                               savePredictions = TRUE)
  lm_model <- train (fractionOfAllSegSites~ ., data = featuredf_nzvFilt_corFilt, 
                     trControl=trainControl, method="lm", preProcess=NULL) # could do some preprocessing here, but I did it already 
  # broom works on these ouputs! glance tells you fit summary and tidy gives you a nice dataframe that is easy to plot  
  lm_finalModel <- lm_model$finalModel 
  lm_modelFit <- glance(lm_finalModel)
  lm_modelFit$mutationClassLabel <- mutationClass

  lm_coefs <- tidy(lm_finalModel)
  lm_coefs$mutationClassLabel <- mutationClass 
  lm_coefs$minYvariable <- min(featuredf_nzvFilt$fractionOfAllSegSites)
  lm_coefs$maxYvariable <- max(featuredf_nzvFilt$fractionOfAllSegSites)
  
  # make some plots: 
  # variable importance using vip
  vipplot <- vip(lm_finalModel,num_features = 100,include_type = T,geom = "point")+theme_bw()+ggtitle(paste0(lm_model$modelInfo$label," coefficients\n",CVNumber,"-fold Cross Validation\nR^2 = ",round(lm_modelFit$r.squared,4)))
  vipplot
  ggsave(paste0(plotdir,mutationClass,".",orderLabel,".vlm.ipPlot.png"),vipplot,height=8,width=5)
  # note that t stat is coeff/std error. so if a coeff is much bigger than its std error then it probably doesn't equal zero --> big t stat -> low p value (highly significant)

  # base R does some nice easy plots:
  png(paste0(plotdir,mutationClass,".",orderLabel,".baseR.lm.plots.png"))
  par(mfrow = c(2, 2))
  plot(lm_finalModel)
  dev.off()
  par(mfrow=c(1,1))

  all_lmModelsPerClass_modelFit <- rbind(all_lmModelsPerClass_modelFit,lm_modelFit)
  all_lmModelsPerClass_modelCoefficients <- rbind(all_lmModelsPerClass_modelCoefficients,lm_coefs)
}

lm_rsq_plot <- ggplot(all_lmModelsPerClass_modelFit,aes(y=mutationClassLabel,x=r.squared))+
  geom_col()+
  ggtitle(paste0("linear model (no lasso) R^2\n",orderLabel))+
  theme_bw()

ggsave(paste0(plotdir,"allClasses.lm.",orderLabel,".R2Plot.png"),lm_rsq_plot,height=8,width=5)

# to compare coefficients between models, L & S divided by range of Y variable per class (so I'm saving ymin and ymax for each model with model coefs)

########## >>> sandbox get second order terms for adjacent positions #########
test=select(mutationSpectrum_featureMatrix,-names(spectrum_onepop))
test=test[1:10,]
test
pairs <- data.frame(t(combn(names(test),2)))
pairs
pairs$X1.position <- as.numeric(unlist(lapply(strsplit(pairs$X1,"_"),"[",4)))
pairs$X2.position <- as.numeric(unlist(lapply(strsplit(pairs$X2,"_"),"[",4)))
pairs$X1.ancOrDerived <- unlist(lapply(strsplit(pairs$X1,"_"),"[",5))
pairs$X2.ancOrDerived <- unlist(lapply(strsplit(pairs$X2,"_"),"[",5))

pairs$positionDistance <- pairs$X1.position - pairs$X2.position # so can select features that are adjacent (distance = 1 or -1? or plus 0?)
head(pairs)
# does it make sense to have interactions between ancestral and derived features? I guess? 
dim(pairs)
## aha so L&s had 4752 features: 4752= 4560 pairs, plus the original 96 (first order), plus an extra 96 #that I think are each feature sqaured. what is the extra 96? includes the feature *2 squared? Seems like #it. does that make sense though? why would it interact with itself? Ask Will/Kelley?
dim(pairs[pairs$X1.ancOrDerived ==pairs$X2.ancOrDerived & pairs$positionDistance %in% c(-1,0,1),]) # can be interaction between two shapes at same site, or one site/step over
# or should it just be within a 7mer -- so just ancestral 
# what terms actually make sense? interactions among adjacent positions only? what counts as adjacent? 
pairs$X1.feature <- unlist(lapply(strsplit(pairs$X1,"_"),"[",3))
pairs$X2.feature <- unlist(lapply(strsplit(pairs$X2,"_"),"[",3))

View(pairs[pairs$X1.ancOrDerived ==pairs$X2.ancOrDerived & pairs$positionDistance %in% c(-1,1) & pairs$X1.feature!=pairs$X2.feature,]) # can be interaction between two shapes at same site, or one site/step over


