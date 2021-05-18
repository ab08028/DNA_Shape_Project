################# DNA shapeR sandbox ##################
# to install
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/sandbox_plots/"
dir.create(plotdir,recursive = T,showWarnings = F)
#BiocManager::install("DNAshapeR")
library(DNAshapeR)
require(ggplot2)
require(ggfortify) # for pca autoplot
library(Biostrings) # needed to go with DNAShapeR
require(dplyr)
require(reshape2)
require(phylotools) # to get names of sequences from fasta
require(caret) # for linear modeling
require(broom) # for linear modeling tidying (tidy, augment to combine with original data, glance)
require(GGally)
# testFasta <- "/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/DNAShapeR/sandbox.all7mers.uniq.FromBear.someMayBeMissing.RevCompCollapsed.fasta" # this is based on mutation types
# this is based on targets:
testFasta='/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/DNAShapeR/sandbox.7mertargets.fromBeluga.someMayBeMissing.RevCompCollapsed.fasta'
# eventually get from mice /Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/summaryTables/ALLINTERVALS.mouse.mutyper.targets.7mer.nostrict.txt # doesn't matter wht species, I just need the 7mers.
# following this vignette https://bioconductor.org/packages/release/bioc/vignettes/DNAshapeR/inst/doc/DNAshapeR.html

seqNames=get.fasta.name(testFasta, clean_name = TRUE)

#x <- expand.grid(rep(list(c('A', 'G', 'T', 'C')), 7))
#allPossible7mers <- do.call(paste0, x)
# not quite right -- need to collapse by reverse complement though! 
# just getting from an output file for now;
# get from bear output: splitting 7mers and sort/uniq and then adding a >seq header to each to make it a fasta
# some could be missing 


### DO THIS IN THE SHELL: 
# testFasta="/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/DNAShapeR/sandbox.all7mers.uniq.FromBear.someMayBeMissing.RevCompCollapsed.fasta"
# > $testFasta

# awk '{print $2}' /Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210331_7mer_NoStrict/mapped_to_brown_bear/summaryTables/mapped_to_brown_bear.AllPops.bears.7merspectrum.SummedOverAllIntervals.ALLFREQS.NOSTRICT.txt | grep -v "variable" | sed 's/\./\n/g' | sort | uniq | while read sequence  ; do  echo -e ">"$sequence"\n"$sequence >> $testFasta; done
### this isn't quite rev-comp-collapsed right because of separating mutation types

######### doing it from targets instead
#testFasta="/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/DNAShapeR/sandbox.7mertargets.fromBeluga.someMayBeMissing.RevCompCollapsed.fasta"
#awk '{print $1}' /Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_targets_files/beluga_whale_GCF_002288925.2_ASM228892v3_dleu1.mutyper.targets.NOSTRICT.7mers.txt | sort | uniq | while read sequence  ; do  echo -e ">"$sequence"\n"$sequence >> $testFasta; done
# get from targets isntead
# okay to get additional shapes you need to add them in manually (no default -- annoying!)
all14ShapeTypes=c("HelT", "Rise", "Roll", "Shift", "Slide", "Tilt", "Buckle", "Opening", "ProT", "Shear", "Stagger", "Stretch", "MGW", "EP")
pred <- getShape(testFasta,parse=T,shapeType = all14ShapeTypes)
pred

# ok I checked that this keeps order (except when making heatmap I think)
# can get names 
# make some plots
# https://www.codecademy.com/articles/normalization ? 
for(shapeType in all14ShapeTypes){
  plotShape(pred$as.character(shapeType))
}
plotShape(pred$MGW)
plotShape(pred$ProT)
plotShape(pred$Roll)
plotShape(pred$HelT)
plotShape(pred$Shift)
# # etc.

heatShape(pred$HelT, nBins=1) # what is this showing???


####### okay this is cool want to have first (or second order) features and explore simialrity
# okay this makes a feature vector of all the things you're itnerested in in order, and normalizes so different units
# don't dominate (cool)
# some features only give you a few bp worht of info (eg Roll gives 4, MGW gives 3) so it's giving 3 values of 3pb for MGW 
# each row is a 7mer 
# and then each 
# can look at second order interactions too using 2-MGW etc (I'm not yet)
#featureTypes <- c("1-Roll","1-MGW","1-ProT","1-HelT","1-EP")
firstOrder_feautureTypes <- paste0("1-",all14ShapeTypes) # get names of all 14 first order feature types
#featureTypes="1-shape" # this is all four of the above 
# or could just do 1-shape for all 1st order interaciotns maybe?
#testFeatureVector <- encodeSeqShape(testFasta,pred,featureTypes, normalize=T) # aha! this will by default normalize the numbers. I guess that makes sense actually

# I want to keep the labels so maybe I can loop over them and label them? let's try that
# this more lines of code than just using "1-shape" for all but it lets me label them! 
# aha the - is a problem in the names. Rename them !
allFeatureVectors_labelled <- data.frame(motif=seqNames)
for(feature in firstOrder_feautureTypes){
  feature_noDash=gsub("-","_",feature)
  featVec <- data.frame(encodeSeqShape(testFasta,pred,feature, normalize=T))
  # label columns: 
  featVecNum=as.numeric(unlist(dim(featVec)[2]))
  colnames(featVec) <- unlist(paste0("feature_",feature_noDash,"_",as.character(seq(1,featVecNum))))
  # then want to bind columns:
  allFeatureVectors_labelled <- cbind(allFeatureVectors_labelled,featVec)
}

######### NOTE: the output of pred *MUST* be in the same order as in the input fasta. I checked that it is but if anything disturbs this things will be mislabelled!

allFeatureVectors_labelled$centralBP <- substr(allFeatureVectors_labelled$motif,start=4,stop=4)
allFeatureVectors_labelled$central3mer <- substr(allFeatureVectors_labelled$motif,start=3,stop=5)
allFeatureVectors_labelled$central5mer <- substr(allFeatureVectors_labelled$motif,start=2,stop=6)

# some optional labels:
allFeatureVectors_labelled$label <- "-"
allFeatureVectors_labelled[allFeatureVectors_labelled$motif=="TTTAAAA",]$label <- "TTTAAAA"
allFeatureVectors_labelled$InterestingXpYInCentral3merLabel <- "None"
allFeatureVectors_labelled[grep('TA',allFeatureVectors_labelled$central3mer),]$InterestingXpYInCentral3merLabel <- "TpA in central 3mer"
allFeatureVectors_labelled[grep('AT',allFeatureVectors_labelled$central3mer),]$InterestingXpYInCentral3merLabel <- "ApT in central 3mer"
allFeatureVectors_labelled[grep('CG',allFeatureVectors_labelled$central3mer),]$InterestingXpYInCentral3merLabel <- "CpG in central 3mer"
allFeatureVectors_labelled[grep('GC',allFeatureVectors_labelled$central3mer),]$InterestingXpYInCentral3merLabel <- "GpC in central 3mer"
#allFeatureVectors_labelled[allFeatureVectors_labelled$motif=="TTTTAAAA",]$label <- "TTTTAAA" # this one is collapsed into the above one in the rev complement collapse

head(allFeatureVectors_labelled) # how can I get sequence names in there as in the example?
# okay this is making a vector of central 3mer values 
# without norm (is the same as pred):  3.38 3.38 3.38 
#                                     4.65 4.97 5.21
# with norm:

# goal do a PCA based on normalized features to see what 7mers are similar?


# need the non-numeric columns to not be in here (using subset)
pca <- prcomp(subset(allFeatureVectors_labelled,select=-c(motif,centralBP,central3mer,label,InterestingXpYInCentral3merLabel,central5mer)), scale=T,center=T) # Michael says to scale/center; similar with and without


# maybe I also want to try to subset each feature and plot it?


############## a few different pca plots ####################
# am making our fave 7mer a big triangle
# color by central 3mer
threemerPCAPlot <- autoplot(pca,data=allFeatureVectors_labelled,colour="central3mer",shape="label",size="label",loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black")+ggtitle(toString(all14ShapeTypes))+theme_bw()+scale_size_manual(values=c(1,8))
threemerPCAPlot
ggsave(paste0(plotdir,"PCA.7mertargets.someMayBeMissing.ColoredByCentral3mer.All14Features.png"),threemerPCAPlot,height=8,width=15)


fivemerPCAPlot <- autoplot(pca,data=allFeatureVectors_labelled,colour="central5mer",shape="label",size="label",loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black")+ggtitle(toString(all14ShapeTypes))+theme_bw()+scale_size_manual(values=c(1,8))+theme(legend.position = "none")
# fivemerPCAPlot -- slow!
ggsave(paste0(plotdir,"PCA.7mertargets.someMayBeMissing.ColoredByCentral5mer.All14Features.png"),fivemerPCAPlot,height=8,width=15)


# color by containing TpA or ApT step or not:
TpA_CpG_pcaPlot <- autoplot(pca,data=allFeatureVectors_labelled,colour="InterestingXpYInCentral3merLabel",shape="label",size="label",loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black")+ggtitle(toString(all14ShapeTypes))+theme_bw()+scale_size_manual(values=c(1,8))
TpA_CpG_pcaPlot
ggsave(paste0(plotdir,"PCA.7mertargets.someMayBeMissing.ColoredByTpACpGContent.All14Features.png"),TpA_CpG_pcaPlot,height=8,width=15)

# color by central mutaiton type
centralBPPCAPlot <- autoplot(pca,data=allFeatureVectors_labelled,colour="centralBP",shape="label",size="label",loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black")+ggtitle(toString(all14ShapeTypes))+theme_bw()+scale_size_manual(values=c(1,8))
centralBPPCAPlot
ggsave(paste0(plotdir,"PCA.7mertargets.someMayBeMissing.ColoredByCentralBP.All14Features.png"),centralBPPCAPlot,height=8,width=15)



################## looking to see which features are different at each bp by doing a pca based on the _1 _2 _3 feautres ################
for(feature in all14ShapeTypes){
  feature_noDash=gsub("-","_",feature)
  featureColumns=names(allFeatureVectors_labelled)[grep(feature_noDash,names(allFeatureVectors_labelled))] # eg  "feature_1_EP_1" "feature_1_EP_2" "feature_1_EP_3"
  pca_oneFeature <- prcomp(subset(allFeatureVectors_labelled,select=featureColumns), scale=T,center=T) # Michael says to scale/center; similar with and without
  feature_PCAPlot <- autoplot(pca_oneFeature,data=allFeatureVectors_labelled,colour="central3mer",shape="label",size="label",loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black")+ggtitle(toString(feature))+theme_bw()+scale_size_manual(values=c(1,8))
  feature_PCAPlot
  ggsave(paste0(plotdir,"PCA.",feature_noDash,"MeasurementsForEachBP.7mertargets.someMayBeMissing.central3merColor.png"),feature_PCAPlot,height=8,width=15)
  
  
}


################# do some basic plotting of each parameter ###########
head(allFeatureVectors_labelled)
allFeatureVectors_labelled_melt <- melt(allFeatureVectors_labelled)
allFeatureVectors_labelled_melt$overallFeature <- unlist(lapply(strsplit(as.character(allFeatureVectors_labelled_melt$variable),"_"),"[",3))

featurePlot1 <- ggplot(allFeatureVectors_labelled_melt,aes(x=variable,y=value,group=motif,color=central3mer))+
  facet_grid(~central3mer~overallFeature,scales="free")+
  geom_line(size=0.01,alpha=0.6)+
  geom_point(size=0.1)+
  theme_bw()+
  ggtitle("How 7mer shape features change over central bps (note: NORMALIZED!)\nseparated by central 3mer")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90),text=element_text(size=14))
featurePlot1
ggsave(paste0(plotdir,"FeaturesVaryignOverBp.7mertargets.someMayBeMissing.central3merColor.all14Features.png"),featurePlot1,height=24,width=30)

################################# second order feature types ########################################
#SECONDORDER_featureTypes=c("2-Roll","2-MGW","2-ProT","2-HelT","2-EP")
SECONDORDER_featureTypes=paste0("2-",all14ShapeTypes)
SECONDORDER_allFeatureVectors_labelled <- data.frame(motif=seqNames)
for(feature in SECONDORDER_featureTypes){
  feature_noDash=gsub("-","_",feature)
  featVec <- data.frame(encodeSeqShape(testFasta,pred,feature, normalize=T))
  # label columns: 
  featVecNum=as.numeric(unlist(dim(featVec)[2]))
  colnames(featVec) <- unlist(paste0("feature_",feature_noDash,"_",as.character(seq(1,featVecNum))))
  # then want to bind columns:
  SECONDORDER_allFeatureVectors_labelled <- cbind(SECONDORDER_allFeatureVectors_labelled,featVec)
}



SECONDORDER_allFeatureVectors_labelled$centralBP <- substr(SECONDORDER_allFeatureVectors_labelled$motif,start=4,stop=4)
SECONDORDER_allFeatureVectors_labelled$central3mer <- substr(SECONDORDER_allFeatureVectors_labelled$motif,start=3,stop=5)
SECONDORDER_allFeatureVectors_labelled$central5mer <- substr(SECONDORDER_allFeatureVectors_labelled$motif,start=2,stop=6)

# some optional labels:
SECONDORDER_allFeatureVectors_labelled$label <- "-"
SECONDORDER_allFeatureVectors_labelled[SECONDORDER_allFeatureVectors_labelled$motif=="TTTAAAA",]$label <- "TTTAAAA"
SECONDORDER_allFeatureVectors_labelled$InterestingXpYInCentral3merLabel <- "None"
SECONDORDER_allFeatureVectors_labelled[grep('TA',SECONDORDER_allFeatureVectors_labelled$central3mer),]$InterestingXpYInCentral3merLabel <- "TpA in central 3mer"
SECONDORDER_allFeatureVectors_labelled[grep('AT',SECONDORDER_allFeatureVectors_labelled$central3mer),]$InterestingXpYInCentral3merLabel <- "ApT in central 3mer"
SECONDORDER_allFeatureVectors_labelled[grep('CG',SECONDORDER_allFeatureVectors_labelled$central3mer),]$InterestingXpYInCentral3merLabel <- "CpG in central 3mer"
SECONDORDER_allFeatureVectors_labelled[grep('GC',SECONDORDER_allFeatureVectors_labelled$central3mer),]$InterestingXpYInCentral3merLabel <- "GpC in central 3mer"
#SECONDORDER_allFeatureVectors_labelled[SECONDORDER_allFeatureVectors_labelled$motif=="TTTTAAAA",]$label <- "TTTTAAA" # this one is collapsed into the above one in the rev complement collapse

head(SECONDORDER_allFeatureVectors_labelled) # how can I get sequence names in there as in the example?

############## plot second order interactions ################
head(allFeatureVectors_labelled)
SECONDORDER_allFeatureVectors_labelled_melt <- melt(SECONDORDER_allFeatureVectors_labelled)
SECONDORDER_allFeatureVectors_labelled_melt$overallFeature <- unlist(lapply(strsplit(as.character(SECONDORDER_allFeatureVectors_labelled_melt$variable),"_"),"[",3))

SECONDORDER_featurePlot1 <- ggplot(SECONDORDER_allFeatureVectors_labelled_melt,aes(x=variable,y=value,group=motif,color=central3mer))+
  facet_grid(~central3mer~overallFeature,scales="free")+
  geom_line(size=0.01,alpha=0.6)+
  geom_point(size=0.1)+
  theme_bw()+
  ggtitle("How SECOND ORDER 7mer shape features change over central bps (note: NORMALIZED!)\nseparated by central 3mer")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90),text=element_text(size=14))
# SECONDORDER_featurePlot1
ggsave(paste0(plotdir,"FeaturesVaryignOverBp.SECONDORDER.7mertargets.someMayBeMissing.central3merColor.png"),SECONDORDER_featurePlot1,height=24,width=15)


##################### second order PCA ##########     

SECONDORDER_pca <- prcomp(subset(SECONDORDER_allFeatureVectors_labelled,select=-c(motif,centralBP,central3mer,label,InterestingXpYInCentral3merLabel,central5mer)), scale=T,center=T) # Michael says to scale/center; similar with and without


SECONDORDER_threemerPCAPlot <- autoplot(SECONDORDER_pca,data=SECONDORDER_allFeatureVectors_labelled,colour="central3mer",shape="label",size="label",loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black")+ggtitle(toString(SECONDORER_featureTypes))+theme_bw()+scale_size_manual(values=c(1,8))
SECONDORDER_threemerPCAPlot
ggsave(paste0(plotdir,"PCA.7mertargets.someMayBeMissing.ColoredByCentral3mer.SECONDORDERFEATURES.all14Features.png"),SECONDORDER_threemerPCAPlot,height=8,width=15)

######### don't need so many plots ##########
#SECONDORDER_fivemerPCAPlot <- autoplot(SECONDORDER_pca,data=SECONDORDER_allFeatureVectors_labelled,colour="central5mer",shape="label",size="label",loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black")+ggtitle(toString(SECONDORER_featureTypes))+theme_bw()+scale_size_manual(values=c(1,8))
#fivemerPCAPlot
#ggsave(paste0(plotdir,"PCA.7mertargets.someMayBeMissing.ColoredByCentral5mer.SECONDORDERFEATURES.all14Features.png"),fivemerPCAPlot,height=8,width=15)


# color by containing TpA or ApT step or not:
#SECONDORDER_TpA_CpG_pcaPlot <- autoplot(SECONDORDER_pca,data=SECONDORDER_allFeatureVectors_labelled,colour="InterestingXpYInCentral3merLabel",shape="label",size="label",loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black")+ggtitle(toString(SECONDORER_featureTypes))+theme_bw()+scale_size_manual(values=c(1,8))
#TpA_CpG_pcaPlot
#ggsave(paste0(plotdir,"PCA.7mertargets.someMayBeMissing.ColoredByTpACpGContent.SECONDORDERFEATURES.all14Features.png"),TpA_CpG_pcaPlot,height=8,width=15)

# color by central mutaiton type
#centralBPPCAPlot <- autoplot(pca,data=allFeatureVectors_labelled,colour="centralBP",shape="label",size="label",loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black")+ggtitle(toString(SECONDORER_featureTypes))+theme_bw()+scale_size_manual(values=c(1,8))
#centralBPPCAPlot
#ggsave(paste0(plotdir,"PCA.7mertargets.someMayBeMissing.ColoredByCentralBP.all14Features.png"),centralBPPCAPlot,height=8,width=15)



#### SANDBOX: from vignette: training a linear model ##############
# need to read up on:
#Then, a machine learning package (which can be any learning tools) is used to train a multiple linear regression (MLR) model based on 3-fold cross-validation. In this example, we used the caret package (see http://caret.r-forge.r-project.org/ for more information).
###################################### mouse and vaquita mutation spectr ####################
################# before doing a linear model you must check :  #######
# 1. linearity of relationships between every variable and the dependent variable
# 2. autocorrelation between independent varibles
# 3. whether dependent variable is ~normally distributed
# 4. after you run the model, check for Step 4: Check for homoscedasticity
mouseSpectrumData <- read.table('/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/summaryTables/Mouse7merSpectrum.SummedUpOverAllChrs.NOSTRICT.txt',header=T)
# select just the mmd spectrum vector to be the data
# need just the ancestral 7mer
mouseSpectrumData$ancestral7mer <- unlist(lapply(strsplit(as.character(mouseSpectrumData$variable),"\\."),"[",1))
# combine (merge) with feature Vector:
mouse_merged <- merge(allFeatureVectors_labelled,mouseSpectrumData[mouseSpectrumData$population=="Mmd",c("variable","fractionOfAllSegSites","ancestral7mer")],by.x="motif",by.y="ancestral7mer")
head(mouse_merged)

# 1. / 2. use a correlogram to check if there's a linear relationship between variables / check for ~ normal dist of dependent variable / 
## trying to plot correlation between variables
# this is slow!
correlogram <- ggpairs(subset(mouse_merged,select=-c(motif,centralBP,central3mer,label,InterestingXpYInCentral3merLabel,variable)))
ggsave(paste0(plotdir,"mouse.Correlogram.AllShapeParameters.PlusMouseFracSegSites.CheckForLinearity.Normality.png"),correlogram,height=20,width=20)

mouse_model_df <- subset(mouse_merged,select=-c(motif,centralBP,central3mer,label,InterestingXpYInCentral3merLabel,variable))
# note tht I'm taking out all the labels but will put them bck in??
# now want to make  model in caret based on the vignette: https://bioconductor.org/packages/release/bioc/vignettes/DNAshapeR/inst/doc/DNAshapeR.html




############ modeling with caret (skipping got same coeffs from lm )##############
# sets up parameters: cv=??? 
##### USE IF USING CARET: trainControl <- trainControl(method = "cv", number = 3, 
 #                            savePredictions = TRUE)
# then make the model: this is trying to predict the frac of seg sites for each mutation type based on the shape of the ancestral 7mer -- not sure what is going on here!! this is from the "caret" package as suggested by the DNAShapeR tutorial but I don't actually know what I"m doing here
##### USE IF USING CARET: mouse_model <- train (fractionOfAllSegSites~ ., data = mouse_model_df, 
               # trControl=trainControl, method="lm", preProcess=NULL) # okay this produces the same thing as just simple lm !!! # so this is currently adding unnecessary complexity. 
# summary(mouse_model)
# can use mouse_model$finalModel to explore .. $coefficients?  $residuals?
# hm. okay what next.
# then I want those coefficients
# am I doing this remotely right?? find a tutorial
# then want to figure out how to use it to predict rtes in aother species

############## model mouse with lm which I actually understand better #################
# gives same coeffs as tht complicated caret 'train' thing: frac~. means verses all variables. 
mouse_model_simple_lm = lm(fractionOfAllSegSites~., data=mouse_model_df)
summary(mouse_model_simple_lm)
############ tidy up with broom  ###########
glance(mouse_model_simple_lm)  # gives a summary of model fit with r.squared 
mouse_model_simple_lm_tidy <- tidy(mouse_model_simple_lm)  # makes a nice little dataframe 
# use augment to add data back in: this is super nice:
mouse_model_augmentedWithData <- augment(mouse_model_simple_lm,mouse_merged)

# check tht model is good fit to data
png(paste0(plotdir,"mouseMouse.lm.residuals.png"))
par(mfrow=c(2,2))
plot(mouse_model_simple_lm)
dev.off()
par(mfrow=c(1,1)) # resets to normal plotting
# how to interpret these residuals? looks like there's a bias.


############# vaquita model ##############
vaquitaSpectrumData=read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210225_7mer/vaquita.7merspectrum.SummedOverAllIntervals.txt",header=T,sep="\t")
head(vaquitaSpectrumData)
vaquitaSpectrumData$ancestral7mer <- unlist(lapply(strsplit(as.character(vaquitaSpectrumData$variable),"\\."),"[",1))
head(vaquitaSpectrumData)

vaquita_merge <- merge(allFeatureVectors_labelled,vaquitaSpectrumData[,c("variable","fractionOfAllSegSites","ancestral7mer")],by.x="motif",by.y="ancestral7mer")
head(vaquita_merge)

vaquita_model_df <- subset(vaquita_merge,select=-c(motif,centralBP,central3mer,label,InterestingXpYInCentral3merLabel,variable))

# skipping caret version in favor of regular lm
#vaquita_model <- train (fractionOfAllSegSites~ ., data = vaquita_model_df, 
#                      trControl=trainControl, method="lm", preProcess=NULL)
#summary(vaquita_model)

vaquita_model_simple_lm = lm(fractionOfAllSegSites~., data=vaquita_model_df)
summary(vaquita_model_simple_lm)
############ tidy up with broom  ###########
glance(vaquita_model_simple_lm)  # gives a summary of model fit with r.squared 
vaquita_model_simple_lm_tidy <- tidy(vaquita_model_simple_lm)  # makes a nice little dataframe 
# use augment to add data back in: this is super nice:
vaquita_model_augmentedWithData <- augment(vaquita_model_simple_lm,vaquita_merge)

####################### compare coefficients #####################
# can I plot correlation of coefficients?
comboCoeffs <- merge(mouse_model_simple_lm_tidy,vaquita_model_simple_lm_tidy,by="term",suffixes=c(".mouse",".vaquita"))

mouse_vaquita_coefficientsPlot <- ggplot(comboCoeffs,aes(x=estimate.mouse,y=estimate.vaquita,color=term))+
  geom_point()+
  geom_abline(intercept=0,slope=1)
mouse_vaquita_coefficientsPlot
ggsave(paste0(plotdir,"mouse.vaquita.compareLmCoefficients.png"),mouse_vaquita_coefficientsPlot,height=4,width=6)
# okay so now I have these mouse coefficients. Want to see how well they can predict some other species fract of seg sites?
########## OKAY! stop here and read about caret and make sure you UDNERSTAND WHAT IS GOING ON! can then try to use mouse coeffs in vaquita model or something but do NOT PROCEED UNTIL YOU"VE TAKEN TIME TO LEARN SOME STUFF ###########
################ SLOW. DOWN. #################

# not sure if it makes sense to make a model based on all these features??

######## SANDBOX: experiment with deltas idea- want to subtract one type from another to get deltas  ######


# motif pairs: this is just first two rows to get deltas. 
deltaTest <- subset(allFeatureVectors_labelled,select=-c(motif,centralBP,central3mer,label,InterestingXpYInCentral3merLabel))[1,] - subset(allFeatureVectors_labelled,select=-c(motif,centralBP,central3mer,label,InterestingXpYInCentral3merLabel))[2,]
colnames(deltaTest) <- paste0("delta_",names(deltaTest))
head(deltaTest)
############## sandbox plotting features ###########
