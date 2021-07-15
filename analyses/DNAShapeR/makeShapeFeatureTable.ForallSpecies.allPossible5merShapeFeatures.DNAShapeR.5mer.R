###### Get features of all possible 5mers  (not rev-comped; will be once I merge with a spectrum ###########
require(tidyr)
require(DNAshapeR)
require(phylotools) # to get names of sequences from fasta

fastadir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/fastasOf5mers/"
plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/plots/"
dir.create(fastadir,recursive = T,showWarnings = F)
#### >>> get  all possible 5mers (not rev comped because don't need to at this stage) #####
bases=c("A","C","G","T")
kmerlength=5
allkmers <- expand.grid((rep(list(bases),kmerlength)))
dim(allkmers) # 1024 -- way fewer than 7mers, nice

# paste all columns together:
allkmers_united <- allkmers %>%
                          unite(allkmers,1:ncol(allkmers),sep="")
# should all be unique

#### >>> write to fasta ####
fastaContainingAllKmers=paste0(fastadir,"allPossible5mers.NotRevCompedCollapsedBecauseNotNeeded.notSpeciesSpecific.WillWorkForAnySpecies.usethis.fasta")

sink(fastaContainingAllKmers)
# using unique to not do repeats (there are repeats because of X>Y and X>Z mutations. just want to get shape features of each one once)
for(kmer in allkmers_united$allkmers){
  cat(paste0(">",kmer,"\n"))
  cat(paste0(kmer,"\n"))
}
sink()


############# make a CpG position file ###############
# make an optional methylatedPosFile file with locations of CpGs so you can use methylate=T (note that we don't actually know which are methylated and which aren't so just going to try with all CpGs methylated or none. 

# need to get numeric positions within sequence of Cs of CpGs:
# this works!
#test="CGCGCGCCTT"
#unlist(gregexpr(pattern ='CG',test)) # gives CpG positions - yes! 
CpGPositionFile=paste0(fastadir,"allPossible5mers.CpGPositions.fasta")

sink(CpGPositionFile)
for(kmer in allkmers_united$allkmers){
  CpGPositions= unlist(gregexpr(pattern ='CG',kmer)) # this gives positions in the string of CpGs
  # this returns -1 if not present
  # so make sure it != -1
  if(!(-1 %in% CpGPositions)){
    cat(paste0(">",kmer,"\n"))
    cat(paste0(gsub("\\s+","",(toString(CpGPositions))),"\n")) # the gsub gets rid of spaces so that positions are like : 1,3,7 not 1, 3, 7 
  }

}
sink()




################ >>> get sequence names from the fasta files ##########################
# get sequence names in order from the fasta: 
kmer_SeqNames_fromFasta=get.fasta.name(fastaContainingAllKmers, clean_name = TRUE)


############## >>> prepare dna shape R parameters ###########################
all14ShapeTypes=c("HelT", "Rise", "Roll", "Shift", "Slide", "Tilt", "Buckle", "Opening", "ProT", "Shear", "Stagger", "Stretch", "MGW", "EP") # default is just HelT proT EP MGW and Roll -- need to manually list all these others

############## >>> make DNAShapeR shape prediction for each fasta sequence. ############

DNAShapeR_Prediction <- getShape(fastaContainingAllKmers,parse=T,shapeType = all14ShapeTypes) # runs fast (a minute or two)

# way faster with 5mers 
# this yields a massive list of dataframes that has all entries for every 5mer in each fasta. Note that it is done in sliding 5mers so you only get values for the center 3bp and their flanking steps which are at the center of each 5mer. but those data are informed by the whole sequence. Some measurements are per bp and some are per bp 'step' between bps depending on whether its an intra or inter-bp measurement. 

########## GET METHYLATED VERSION (Only works for basic shape features, not extended list) ###################
# #from manual:
# " To achieve a better mechanistic understanding of the effect of CpG methylation on local DNA structure, we developed a high-throughput method, named methyl-DNAshape, for predicting the impact of cytosine methylation on DNA shape features. In analogy to the DNAshape method (Zhou, et al., 2013), the method predicts DNA shape features (ProT, HelT, Roll, and MGW) in the context of CpG methylation based on methyl-DNAshape Pentamer Query Table (mPQT) derived from the results of all-atom Monte Carlo simulations on a total of 3,518 DNA fragments of lengths varying from 13 to 24 bp."  
# To predict DNA shape features in the context of CpG methylation, in addition to providing regular FASTA file (without symbolizing ‘Mg’) one can provide an additional input file identifying methylated positions. For example,

#>seq1

#4,16 # notifying the cytosine at position 4th and 16th on leading strand, and 5th and 17th on lagging strand are methylated.

# okay I like the second way better.
# so I want to get the CpG positions for every sequence and write out a secondary file with CpG positions. 

DNAShapeR_Prediction_METHYLATED <- getShape(fastaContainingAllKmers,parse=T,methylate = T,methylatedPosFile = CpGPositionFile) # runs fast (a minute or two)
# shapeType = all14ShapeTypes
# compare
#kmer_SeqNames_fromFasta[794] # "CGCATAA"
#DNAShapeR_Prediction_METHYLATED$Roll[794,] #NA -3.63  4.51 -3.08  7.92    NA
#DNAShapeR_Prediction$Roll[794,] #    NA -0.76  5.14 -3.12  7.92    NA
# changes things when it thinkst the CpG is methylated

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
  # note that for sequence vectors they encode each bp as four columns (1000 0100 etc). So this numbering system doesn't work great for that. 
  colnames(featVec) <- unlist(paste0("feature_",feature_noDash,"_",as.character(seq(1,featVecNum))))
  
  # then want to bind columns:
  allFeatureVectors_labelled <- cbind(allFeatureVectors_labelled,featVec)
}
# so a first order feature will be: feature_1_ and a second order will be feature_2. 
head(allFeatureVectors_labelled)
dim(allFeatureVectors_labelled)

# write out the table:
write.table(allFeatureVectors_labelled,paste0("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/",orderLabel,"_allPossible5mers.FeatureValues.min-maxNormalized.WillWorkForAnySpecies.goodForPCAOrAnythingWhereDontWantUnits.txt"),row.names = F,quote=F,sep = "\t")

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
write.table(allFeatureVectors_labelled_NOTnormalized,paste0("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/",orderLabel,"_allPossible5mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.txt"),row.names = F,quote=F,sep = "\t")


######### write a version with METHYLATED features that is NOT normalized ###########

allFeatureVectors_labelled_NOTnormalized_METHYLATED <- data.frame(motif=kmer_SeqNames_fromFasta) # set up with seqNames -- these *must* be the same as the order of sequences from the fasta.
### the only feautres that have methylated versions are: 
methylatedFeatureNames=names(DNAShapeR_Prediction_METHYLATED)
methylatedFeatureNames
# note feature is like 1-EP 1-Prot etc so doesn't quite match the methylated feature names
# don't want to label things differently becaues I want them to behave in my models the same way just note that  "MGW"  "HelT" "ProT" "Roll" are what differ due to methylation
for(feature in desiredFeatures){
  # check if feature is in the methylated list:
  featureGroup_stripOrder=unlist(lapply(strsplit(feature,"-"),"[",2)) # (grouping with out the 1-)
  if(featureGroup_stripOrder %in% methylatedFeatureNames){
    print(paste0('using methylated version of ',feature))
    featVec <- data.frame(encodeSeqShape(fastaContainingAllKmers,DNAShapeR_Prediction_METHYLATED,feature, normalize=F))
  } else {
    print(paste0("no methylation status for feature ", feature))
    # for non methyl status just use the existing DNAShapeR_Prediction prediction 
  featVec <- data.frame(encodeSeqShape(fastaContainingAllKmers,DNAShapeR_Prediction,feature, normalize=F))
  }
  # label columns based on how many columsn there are for each feature: (adding _1 _2 _3 to the feature name)
  feature_noDash=gsub("-","_",feature) # need to get rid of dash in the name for the columns otherwise R hates it.
  featVecNum=as.numeric(unlist(dim(featVec)[2]))
  # note that for sequence vectors they encode each bp as four columns (1000 0100 etc). So this numbering system doesn't work great for that. so each 7mer gets 28 columns for the 1mer column. Want to rename them somehow feature_1_mer_1 _2 etc is okay I guess. 
  colnames(featVec) <- unlist(paste0("feature_",feature_noDash,"_",as.character(seq(1,featVecNum))))
  
  # then want to bind columns:
  allFeatureVectors_labelled_NOTnormalized_METHYLATED <- cbind(allFeatureVectors_labelled_NOTnormalized_METHYLATED,featVec)
}
# so a first order feature will be: feature_1_ and a second order will be feature_2. ED)


# write out the table:
write.table(allFeatureVectors_labelled_NOTnormalized_METHYLATED,paste0("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/","METHYLATED",orderLabel,"_allPossible5mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.METHYLATED.txt"),row.names = F,quote=F,sep = "\t")
# writing methylated at beg and end of filename so I can't miss it 

########## >>> make a few summary plots #############
########### >>> plot correlations of the features #################
# https://briatte.github.io/ggcorr/
require(GGally)
correlationOfFeaturesAcrossAll5mers_plot <- ggcorr(subset(allFeatureVectors_labelled,select=-c(motif)),layout.exp = 20,hjust=1,vjust=.5,size=1.5)
correlationOfFeaturesAcrossAll5mers_plot
ggsave(paste0(plotdir,"correlationOfFeaturesAcrossAll5mers.",orderLabel,".png"),correlationOfFeaturesAcrossAll5mers_plot,height=8,width=9)

####### >>> plot heatmap of mutation types by their shape features ###########
rownames(allFeatureVectors_labelled) <- allFeatureVectors_labelled$motif
png(paste0(plotdir,"heatmap.shapefeatures.allPossibleShapeFeatures.5mers.png"),height=20,width = 20,units="in",res=300)
heatmap(as.matrix(select(allFeatureVectors_labelled,-c("motif"))))
dev.off()
