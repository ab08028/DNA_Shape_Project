###### Get features of all possible 7mers  (not rev-comped; will be once I merge with a spectrum ###########
require(tidyr)
require(DNAshapeR)
require(phylotools) # to get names of sequences from fasta
#install.packages("BiocManager")
#BiocManager::install("DNAshapeR")
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


############# >>> get sequence features #############
test <- data.frame(encodeSeqShape(fastaContainingAllKmers,feature="7-mer"))

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
sequenceDataFrame_binary <- data.frame(encodeSeqShape(fastaContainingAllKmers,feature="1-mer")) # this encodes the base of each position in the 7mer <- this is what we want ; end up with 28 features (7*4 per 7mer because it encodes A like 1 0 0 0 )
#colnames(test) 
#*** note: a typo in the DNA SHapeR Manual
# correct: (checked their github code)
# A = 1000, C = 0100, G = 0010, T = 0001
dim(sequenceDataFrame_binary) # 28 features (4 per position ACGT)
# checking with their github code:
#k=1
#lookupTable <- diag( 4**k )
#lookupTable
#row.names( lookupTable ) <- mkAllStrings( c( "A", "C", "G", "T" ), k )
# okay the 1-mer encodes the variable at every position (4 columns )
# whereas 2-mer and 3-mer encode the threemer of every position
# so I want to just encode the 1-mer sequence at each position but I also want to name it 

# want to label feature names 
head(sequenceDataFrame_binary)

sequenceDataFrame_binary$motif <- rownames(sequenceDataFrame_binary)
# careful here
assignColNames
listOfColNames <- list()
for(position in seq(1,7)){

  listOfColNames <- c(listOfColNames,unlist(paste0("Pos_",position,"_",c("A","C","G","T"))))
  
}
listOfColNames 
colnames(sequenceDataFrame_binary) <- listOfColNames
head(sequenceDataFrame_binary)
View(sequenceDataFrame_binary)
# works _ I checked carefully! 
# want them to go like Pos_1_A Pos_1_C Pos1_G Pos1_T
#         Pos_1_A Pos_1_C Pos_1_G Pos_1_T Pos_2_A Pos_2_C Pos_2_G Pos_2_T Pos_3_A Pos_3_C Pos_3_G Pos_3_T
#TCCCGAA       0       0       0       1       0       1       0       0       0       1       0       0
#Pos_4_A Pos_4_C Pos_4_G Pos_4_T Pos_5_A Pos_5_C Pos_5_G Pos_5_T Pos_6_A Pos_6_C Pos_6_G Pos_6_T
#TCCCGAA       0       1       0       0       0       0       1       0       1       0       0       0
#Pos_7_A Pos_7_C Pos_7_G Pos_7_T
#TCCCGAA       1       0       0       0

# write out the table:
write.table(sequenceDataFrame_binary,paste0("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/BinarySequenceEncoding_1bpPerPosition_allPossible7mers.WillWorkForAnySpecies.txt"),row.names = F,quote=F,sep = "\t")


