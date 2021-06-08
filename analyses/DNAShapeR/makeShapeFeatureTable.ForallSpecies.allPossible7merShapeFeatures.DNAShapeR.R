###### Get features of all possible 7mers  (not rev-comped; will be once I merge with a spectrum ###########
require(tidyr)
require(DNAshapeR)
require(phylotools) # to get names of sequences from fasta

fastadir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/fastasOf7mers/"
plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/plots/"
dir.create(fastadir,recursive = T,showWarnings = F)
#### >>> get  all possible 7mers (not rev comped because don't need to at this stage) #####
bases=c("A","C","G","T")
kmerlength=7
allkmers <- expand.grid((rep(list(bases),kmerlength)))
dim(allkmers)

# paste all columns together:
allkmers_united <- allkmers %>%
                          unite(allkmers,1:ncol(allkmers),sep="")
# should all be unique

#### >>> write to fasta ####
fastaContainingAllKmers=paste0(fastadir,"allPossible7mers.NotRevCompedCollapsedBecauseNotNeeded.notSpeciesSpecific.WillWorkForAnySpecies.usethis.fasta")

sink(fastaContainingAllKmers)
# using unique to not do repeats (there are repeats because of X>Y and X>Z mutations. just want to get shape features of each one once)
for(kmer in allkmers_united$allkmers){
  cat(paste0(">",kmer,"\n"))
  cat(paste0(kmer,"\n"))
}
sink()

################ >>> get sequence names from the fasta files ##########################
# get sequence names in order from the fasta: 
kmer_SeqNames_fromFasta=get.fasta.name(fastaContainingAllKmers, clean_name = TRUE)


############## >>> prepare dna shape R parameters ###########################
all14ShapeTypes=c("HelT", "Rise", "Roll", "Shift", "Slide", "Tilt", "Buckle", "Opening", "ProT", "Shear", "Stagger", "Stretch", "MGW", "EP") # default is just HelT proT EP MGW and Roll -- need to manually list all these others

############## >>> make DNAShapeR shape prediction for each fasta sequence. ############

DNAShapeR_Prediction <- getShape(fastaContainingAllKmers,parse=T,shapeType = all14ShapeTypes) # runs fast (a minute or two)

# this yields a massive list of dataframes that has all entries for every 7mer in each fasta. Note that it is done in sliding 5mers so you only get values for the center 3bp and their flanking steps which are at the center of each 5mer. but those data are informed by the whole sequence. Some measurements are per bp and some are per bp 'step' between bps depending on whether its an intra or inter-bp measurement. 
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
dim(allFeatureVectors_labelled)

# write out the table:
write.table(allFeatureVectors_labelled,paste0("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/",orderLabel,"_allPossible7mers.FeatureValues.min-maxNormalized.WillWorkForAnySpecies.goodForPCAOrAnythingWhereDontWantUnits.txt"),row.names = F,quote=F,sep = "\t")

# note when you read it in you need to set rownames
######## >>> also write out a version that isn't normalized ##########
allFeatureVectors_labelled_NOTnormalized <- data.frame(motif=kmer_SeqNames_fromFasta) # set up with seqNames -- these *must* be the same as the order of sequences from the fasta.

for(feature in desiredFeatures){
  feature_noDash=gsub("-","_",feature) # need to get rid of dash in the name for the columns otherwise R hates it.
  featVec <- data.frame(encodeSeqShape(fastaContainingAllKmers,DNAShapeR_Prediction,feature, normalize=F))
  # label columns based on how many columsn there are for each feature: (adding _1 _2 _3 to the feature name)
  featVecNum=as.numeric(unlist(dim(featVec)[2]))
  # note that for sequence vectors they encode each bp as four columns (1000 0100 etc). So this numbering system doesn't work great for that. so each 7mer gets 28 columns for the 1mer column. Want to rename them somehow feature_1_mer_1 _2 etc is okay I guess. 
  colnames(featVec) <- unlist(paste0("feature_",feature_noDash,"_",as.character(seq(1,featVecNum))))
  
  # then want to bind columns:
  allFeatureVectors_labelled_NOTnormalized <- cbind(allFeatureVectors_labelled_NOTnormalized,featVec)
}
# so a first order feature will be: feature_1_ and a second order will be feature_2. 
head(allFeatureVectors_labelled_NOTnormalized)
dim(allFeatureVectors_labelled_NOTnormalized)

# write out the table:
write.table(allFeatureVectors_labelled_NOTnormalized,paste0("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/",orderLabel,"_allPossible7mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.txt"),row.names = F,quote=F,sep = "\t")

########## >>> make a few summary plots #############
########### >>> plot correlations of the features #################
# https://briatte.github.io/ggcorr/
correlationOfFeaturesAcrossAll7mers_plot <- ggcorr(subset(allFeatureVectors_labelled,select=-c(motif)),layout.exp = 20,hjust=1,vjust=.5,size=1.5)
correlationOfFeaturesAcrossAll7mers_plot
ggsave(paste0(plotdir,"correlationOfFeaturesAcrossAll7mers.",orderLabel,".png"),correlationOfFeaturesAcrossAll7mers_plot,height=8,width=9)

####### >>> plot heatmap of mutation types by their shape features ###########
rownames(allFeatureVectors_labelled) <- allFeatureVectors_labelled$motif
png(paste0(plotdir,"heatmap.shapefeatures.allPossibleShapeFeatures.png"),height=20,width = 20,units="in",res=300)
heatmap(as.matrix(select(allFeatureVectors_labelled,-c("motif"))))
dev.off()
