################## modeling of mouse spectrum ##############


############ >>> read in train and test data ###################
traindata <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/training_test_datasets/mouse/mutyperResults_20210317_NOSTRICT_7mer/training_and_test_spectra_SUMMEDOVERCHROMOSOMES/mouse.TRAINING.odd_bpOnly.allPops.summedOverChr.spectra.NOSTRICT.txt",header=T)
  
testdata <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/training_test_datasets/mouse/mutyperResults_20210317_NOSTRICT_7mer/training_and_test_spectra_SUMMEDOVERCHROMOSOMES/mouse.TESTING.even_bpOnly.allPops.summedOverChr.spectra.NOSTRICT.txt",header=T)

######## >>> restrict to a single species/population ######################
population="Mmd"
traindata_onepop <- traindata[traindata$population==population,]
testdata_onepop <- testdata[testdata$population==population,]

######### >>> read in shape table for all possible 7mers ############
# using normalized at least for now:
shapes <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/firstOrder_featureTypes_allPossible7mers.FeatureValues.min-maxNormalized.WillWorkForAnySpecies.usethis.txt",header=T,sep="\t")
dim(shapes)
# need to assign row names:
rownames(shapes) <- shapes$motif

######## >>> make training feature matrix by combining the shape features with the training data by ancestral and derived motifs ########
# get ancestral derived:
intermediateMerge_ancestral <- merge(traindata_onepop,shapes,by.x="ancestralkmer",by.y="motif")
# get the derived features (and set suffixes as ancestral and .derived):
training_fullFeature_df <- merge(intermediateMerge_ancestral,shapes,by.x="derivedkmer",by.y="motif",suffixes=c(".ancestral",".derived")) # sweet, I checked that this works and it does! 
rownames(training_fullFeature_df) <- training_fullFeature_df$variable # important so you can keep track of what motif is in what row without a column for it 

head(training_fullFeature_df)
# note here that fractionOfAllSegSites will be the dependent Y variable. and you need to exclude all the other columns that aren't shape features 
columnsToExclude = names(select(traindata_onepop,-fractionOfAllSegSites))
columnsToExclude

training_fullFeature_df <- select(training_fullFeature_df,-all_of(columnsToExclude))
head(training_fullFeature_df)
