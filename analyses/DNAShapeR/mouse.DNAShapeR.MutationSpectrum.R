##################### DNA shape of full mutation spectrum #############
# following this vignette https://bioconductor.org/packages/release/bioc/vignettes/DNAshapeR/inst/doc/DNAshapeR.html
#  install.packages("BiocManager")
fastadir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/fastasOf7mers/"
dir.create(fastadir,recursive = T,showWarnings = F)

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
# for upmap:
require(umap)
# for plotting clusters:
require(ggdendro)
require(dendextend)
### starting with mouse as a case study. Want to make a fasta in R so no bash is needed.  

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
plotdir=paste0("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/plots/",speciesLabel,"/")
dir.create(plotdir,recursive = T,showWarnings = F)

spectrum_onepop$ancestralkmer <- unlist(lapply(strsplit(spectrum_onepop$variable,"\\."),"[",1))

spectrum_onepop$derivedkmer <- unlist(lapply(strsplit(spectrum_onepop$variable,"\\."),"[",2))
# add in central 5mer:

spectrum_onepop$central5mer.ancestral <- substr(spectrum_onepop$ancestralkmer,start=2,stop=6)


############ >>> label central mutation type by whether it is CpG or not and what nucleotides are involved (ala L&S) #########
head(spectrum_onepop)
spectrum_onepop$ancestral3mer <- unlist(lapply(strsplit(spectrum_onepop$ThreemerMutationType,"\\."),"[",1))
spectrum_onepop$derived3mer <- unlist(lapply(strsplit(spectrum_onepop$ThreemerMutationType,"\\."),"[",2))

# identify presence of CpG (CG) [not GC] in central 3mer
spectrum_onepop$ancestral3merCpG <- "no"
spectrum_onepop$derived3merCpG <- "no"
spectrum_onepop[grepl("CG",spectrum_onepop$ancestral3mer),]$ancestral3merCpG <- "yes"
spectrum_onepop[grepl("CG",spectrum_onepop$derivedkmer),]$derived3merCpG <- "yes"

spectrum_onepop$mutationClassLabel <- paste0(spectrum_onepop$centralMutationType,".",spectrum_onepop$ancestral3merCpG,"_ancestralCpG") 
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

mutationSpectrum_featureMatrix <- merge(intermediateMerge_ancestral,allFeatureVectors_labelled,by.x="derivedkmer",by.y="motif",suffixes=c(".ancestral",".derived")) # sweet, I checked that this works and it does! 

# note this may contain multiple entries because of multiple populations (eg many mice spp)

############ >>> PCA -- all mutation types ##############
# make  sure you restricted spectrum to just one pop/species
if(length(unique(mutationSpectrum_featureMatrix$population))>1){
  print("More than one population!!!! Stop here!")
}
rownames(mutationSpectrum_featureMatrix) <- mutationSpectrum_featureMatrix$variable

# set size limits
scaleSizeLimits=c(min(mutationSpectrum_featureMatrix$fractionOfAllSegSites),max(mutationSpectrum_featureMatrix$fractionOfAllSegSites))

pca <- prcomp(select(mutationSpectrum_featureMatrix,-c(names(spectrum_onepop))), scale=T,center=T) # Michael says to scale/center; similar with and without
PCAPlot_all7mers1 <- autoplot(pca,data=mutationSpectrum_featureMatrix,colour="fractionOfAllSegSites",size="fractionOfAllSegSites",alpha="fractionOfAllSegSites")+ggtitle(paste0(orderLabel,"\n",toString(desiredFeatures)))+theme_bw()+scale_size_continuous(limits=scaleSizeLimits,range=c(0.01,8)) +scale_color_viridis_c(option = "B",begin = 0.1) +scale_alpha_continuous(range=c(0.1,1)) #loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black" +scale_color_viridis_d(option = "B",begin = 0.2)
PCAPlot_all7mers1
ggsave(paste0(plotdir,"ALL.7mers.PCA.coloredByFracOfSegSites.",orderLabel,".png"),PCAPlot_all7mers1,height=10,width=14)

# color by 3mer
PCAPlot_all7mers2 <- autoplot(pca,data=mutationSpectrum_featureMatrix,colour="ThreemerMutationType",size="fractionOfAllSegSites",alpha="fractionOfAllSegSites")+ggtitle(paste0(orderLabel,"\n",toString(desiredFeatures)))+theme_bw()+scale_size_continuous(limits=scaleSizeLimits,range=c(0.01,8)) +scale_alpha_continuous(range=c(0.1,1)) #loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black" 
PCAPlot_all7mers2
ggsave(paste0(plotdir,"ALL.7mers.PCA.coloredByCentral3Mer.",orderLabel,".png"),PCAPlot_all7mers2,height=10,width=16)

# color by mutation class
PCAPlot_all7mers3 <- autoplot(pca,data=mutationSpectrum_featureMatrix,colour="mutationClassLabel",size="fractionOfAllSegSites",alpha="fractionOfAllSegSites")+ggtitle(paste0(orderLabel,"\n",toString(desiredFeatures)))+theme_bw()+scale_size_continuous(limits=scaleSizeLimits,range=c(0.01,8))+scale_alpha_continuous(range=c(0.1,1)) #loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black" 
PCAPlot_all7mers3
ggsave(paste0(plotdir,"ALL.7mers.PCA.coloredByMutationClass.",orderLabel,".png"),PCAPlot_all7mers3,height=8,width=16)


########### >>> screeplot ###########
# fyi this is how you calc var explained
# https://datavizpyr.com/how-to-make-scree-plot-in-r-with-ggplot2/ var_explained_df <- data.frame(PC= paste0("PC",1:5), var_explained=(pca_obj$sdev)^2/sum((pca_obj$sdev)^2)) # so tis' std-dev squared divided by sum of all std devs squared

png(paste0(plotdir,"ALL.7mers.PCA.coloredByMutationClass.",orderLabel,"screeplot.png"))
screeplot(pca)+title(paste0(orderLabel,"\n",toString(desiredFeatures)))
dev.off()
# okay so PC3 still explains a lot of variance.  need to plot that too

###### >>>> plot PC 2&3 and 1&3 ###
PCAPlot_all7mers4_PC2_PC3 <- autoplot(pca,x=2,y=3,data=mutationSpectrum_featureMatrix,colour="ThreemerMutationType",size="fractionOfAllSegSites",alpha="fractionOfAllSegSites")+ggtitle(paste0(orderLabel,"\n",toString(desiredFeatures)))+theme_bw()+scale_size_continuous(limits=scaleSizeLimits,range=c(0.01,8)) +scale_alpha_continuous(range=c(0.1,1)) #loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black" +scale_color_viridis_d(option = "B",begin = 0.2)
PCAPlot_all7mers4_PC2_PC3
ggsave(paste0(plotdir,"ALL.7mers.PCA.coloredByFracOfSegSites.",orderLabel,".PC2.PC3.png"),PCAPlot_all7mers4_PC2_PC3,height=10,width=16)

PCAPlot_all7mers4_PC1_PC3 <- autoplot(pca,x=1,y=3,z=2,data=mutationSpectrum_featureMatrix,colour="ThreemerMutationType",size="fractionOfAllSegSites",alpha="fractionOfAllSegSites")+ggtitle(paste0(orderLabel,"\n",toString(desiredFeatures)))+theme_bw()+scale_size_continuous(limits=scaleSizeLimits,range=c(0.01,8))+scale_alpha_continuous(range=c(0.1,1)) #loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black" +scale_color_viridis_d(option = "B",begin = 0.2)
PCAPlot_all7mers4_PC1_PC3
ggsave(paste0(plotdir,"ALL.7mers.PCA.coloredByFracOfSegSites.",orderLabel,".PC1.PC3.png"),PCAPlot_all7mers4_PC1_PC3,height=10,width=16)

########## >>> combine PC combos into one plot ############
comboplot <- grid.arrange(PCAPlot_all7mers2,PCAPlot_all7mers4_PC1_PC3,PCAPlot_all7mers4_PC2_PC3)
comboplot
ggsave(paste0(plotdir,"ALL.7mers.PCA.coloredByFracOfSegSites.",orderLabel,".PCCOMBO.png"),comboplot,height=30,width=20)
# okay so PC3 is important. need to add it in to other plots
######### 3d pca plot? ################




########## >>> SLOW!!!! plot PCA per 3mer (commented out bc slow) ##################
# for(threemer in unique(mutationSpectrum_featureMatrix$ThreemerMutationType)){
#   # this population selection is redundant but extra safe
#   mutationSpectrum_featureMatrix_onePopOnly_one3merOnly <- mutationSpectrum_featureMatrix[mutationSpectrum_featureMatrix$population==population & mutationSpectrum_featureMatrix$ThreemerMutationType==threemer,]
# # want to color by central 5mer
#   mutationSpectrum_featureMatrix_onePopOnly_one3merOnly$central5mer.ancestral <- substr(mutationSpectrum_featureMatrix_onePopOnly_one3merOnly$ancestralkmer,start=2,stop=6)
# # need the non-numeric columns to not be in here (using subset) and want dependent variables likt eh fraction of seg sites to not be in there -- only indepedendent shape/sequence
# # names(spectrum)
# 
# # figure out if symmetric or not
# #TEMP_forPCA$symmetric <- "non-symmetric" # need to make this work better
# #TEMP_forPCA[(TEMP_forPCA$derivedkmer %in% TEMP_forPCA$ancestralkmer) & (TEMP_forPCA$ancestralkmer %in% TEMP_forPCA$derivedkmer),]$symmetric <- "symmetric"
# # use select from dyplyr to get rid of the non-feature columns that came in from the spectrum 
#   rownames(mutationSpectrum_featureMatrix_onePopOnly_one3merOnly) <- mutationSpectrum_featureMatrix_onePopOnly_one3merOnly$variable
#   pca <- prcomp(select(mutationSpectrum_featureMatrix_onePopOnly_one3merOnly,-c(names(spectrum_onepop))), scale=T,center=T) # Michael says to scale/center; similar with and without
# 
#   ############# plot PCA colored by fraction of seg sites
#   PCAPlot <- autoplot(pca,data=mutationSpectrum_featureMatrix_onePopOnly_one3merOnly,colour="central5mer.ancestral",size="fractionOfAllSegSites",alpha=0.6,label.size=3,label.colour="black",label=T)+ggtitle(paste0(threemer,"\n",orderLabel,"\n",toString(desiredFeatures)))+theme_bw()+scale_size_continuous(limits=scaleSizeLimits,range=c(0.01,8)) #)+scale_color_viridis_d(option = "B",begin = 0.2) # loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black" scale_alpha_continuous(range = c(0.1,1))+
#   PCAPlot
#   
#   threemerplotdir=paste0(plotdir,"/pca_Per3mer/")
#   dir.create(threemerplotdir,recursive = T,showWarnings = F)
#   ggsave(paste0(threemerplotdir,threemer,".PCA.coloredByCentral5Mer.",orderLabel,".png"),PCAPlot,height=8,width=10)
# 
# }



############### >>> PCA per mutation type (nine types) #############
for(mutationClass in unique(mutationSpectrum_featureMatrix$mutationClassLabel)){
  # this population selection is redundant but extra safe
  mutationSpectrum_featureMatrix_oneClassOnly <- mutationSpectrum_featureMatrix[mutationSpectrum_featureMatrix$mutationClassLabel==mutationClass,]

    # use select from dyplyr to get rid of the non-feature columns that came in from the spectrum 
  rownames(mutationSpectrum_featureMatrix_oneClassOnly) <- mutationSpectrum_featureMatrix_oneClassOnly$variable
  pca <- prcomp(select(mutationSpectrum_featureMatrix_oneClassOnly,-c(names(spectrum_onepop))), scale=T,center=T) # Michael says to scale/center; similar with and without
  
  ############# plot PCA colored by fraction of seg sites
  PCAPlot <- autoplot(pca,data=mutationSpectrum_featureMatrix_oneClassOnly,colour="ThreemerMutationType",size="fractionOfAllSegSites",alpha="fractionOfAllSegSites")+ggtitle(paste0(mutationClass,"\n",orderLabel,"\n",toString(desiredFeatures)))+theme_bw()+scale_size_continuous(limits=scaleSizeLimits,range=c(0.01,8))+scale_alpha_continuous(range=c(0.1,1))
  PCAPlot
  

  PCAPlot13 <- autoplot(pca,x=1,y=3,data=mutationSpectrum_featureMatrix_oneClassOnly,colour="ThreemerMutationType",size="fractionOfAllSegSites",alpha="fractionOfAllSegSites")+ggtitle(paste0(mutationClass,"\n",orderLabel,"\n",toString(desiredFeatures)))+theme_bw()+scale_size_continuous(limits=scaleSizeLimits,range=c(0.01,8))+scale_alpha_continuous(range=c(0.1,1))
  PCAPlot13
  
  
  PCAPlot23 <-autoplot(pca,x=2,y=3,data=mutationSpectrum_featureMatrix_oneClassOnly,colour="ThreemerMutationType",size="fractionOfAllSegSites",alpha="fractionOfAllSegSites")+ggtitle(paste0(mutationClass,"\n",orderLabel,"\n",toString(desiredFeatures)))+theme_bw()+scale_size_continuous(limits=scaleSizeLimits,range=c(0.01,8))+scale_alpha_continuous(range=c(0.1,1))
  PCAPlot23
  
  
  comboplot_mutationType <- grid.arrange(PCAPlot,PCAPlot13,PCAPlot23)
  
  classplotdir=paste0(plotdir,"/pca_PerMutationClass/")
  dir.create(classplotdir,recursive = T,showWarnings = F)
  ggsave(paste0(classplotdir,mutationClass,".PCA.coloredByCentral3Mer.",orderLabel,".png"),PCAPlot,height=8,width=10)
  ggsave(paste0(classplotdir,mutationClass,".PCA.coloredByCentral3Mer.",orderLabel,"ALLPCs.png"),comboplot_mutationType,height=10,width=8)
  
  
}





########### >>> PCA of targets only (ancestral kmers shape features ) colored by mutation frequency #############
head(spectrum_onepop)
#### IMPORTANT this is assuming a single population -- won't always necessarily be the case so by careful! I filtered spectrum to be one pop up above, but am coding defensively here to not overcount sites.
# restric to one population: 
ancestralKmerTargets_spectrum <- spectrum_onepop %>%
  group_by(ancestralkmer,population) %>%
  summarise(totalMutationsInvolvingAncestralKmer=sum(totalSitesAllIntervals))

ancestralKmerTargets_spectrum$fractionOfAllSegSitesInvolvingAncestralKmer <- ancestralKmerTargets_spectrum$totalMutationsInvolvingAncestralKmer/sum(ancestralKmerTargets_spectrum$totalMutationsInvolvingAncestralKmer)

# combine this with features:
ancestralKmerTargets_spectrum_merged_withFeatures <- merge(ancestralKmerTargets_spectrum,allFeatureVectors_labelled,by.x="ancestralkmer",by.y="motif")
rownames(ancestralKmerTargets_spectrum_merged_withFeatures) <- ancestralKmerTargets_spectrum_merged_withFeatures$ancestralkmer


pca <- prcomp(select(ancestralKmerTargets_spectrum_merged_withFeatures,-c(names(ancestralKmerTargets_spectrum))), scale=T,center=T) # Michael says to scale/center; similar with and without

PCAPlot_targets <- autoplot(pca,data=ancestralKmerTargets_spectrum_merged_withFeatures,colour="fractionOfAllSegSitesInvolvingAncestralKmer",size="fractionOfAllSegSitesInvolvingAncestralKmer",alpha="fractionOfAllSegSitesInvolvingAncestralKmer")+ggtitle(paste0("TARGETS ONLY; ", orderLabel,"\n",toString(desiredFeatures)))+theme_bw()+scale_size_continuous(range=c(0.01,8)) +scale_color_viridis_c(option = "B",begin = 0.1) +scale_alpha_continuous(range=c(0.1,1)) #loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black" +scale_color_viridis_d(option = "B",begin = 0.2)
PCAPlot_targets
ggsave(paste0(plotdir,"ALL.TargetsOnly.PCA.coloredByFracOfSegSites.",orderLabel,".png"),PCAPlot_targets,height=10,width=16)

# color by 3mer:
ancestralKmerTargets_spectrum_merged_withFeatures$central3mer <- substr(ancestralKmerTargets_spectrum_merged_withFeatures$ancestralkmer,start=3,stop=5)
PCAPlot_targets2 <- autoplot(pca,data=ancestralKmerTargets_spectrum_merged_withFeatures,colour="central3mer",size="fractionOfAllSegSitesInvolvingAncestralKmer",alpha="fractionOfAllSegSitesInvolvingAncestralKmer",label=T,label.size=0.1)+ggtitle(paste0("TARGETS ONLY; ", orderLabel,"\n",toString(desiredFeatures)))+theme_bw()+scale_size_continuous(range=c(0.01,8))+scale_alpha_continuous(range=c(0.1,1)) #loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black" +scale_color_viridis_d(option = "B",begin = 0.2)
PCAPlot_targets2
ggsave(paste0(plotdir,"ALL.TargetsOnly.PCA.coloredBy3mer",orderLabel,".png"),PCAPlot_targets2,height=10,width=14)


# PC23

PCAPlot_targets2_PC23 <- autoplot(pca,x=2,y=3,data=ancestralKmerTargets_spectrum_merged_withFeatures,colour="central3mer",size="fractionOfAllSegSitesInvolvingAncestralKmer",alpha="fractionOfAllSegSitesInvolvingAncestralKmer",label=T,label.size=0.1)+ggtitle(paste0("TARGETS ONLY; ", orderLabel,"\n",toString(desiredFeatures)))+theme_bw()+scale_size_continuous(range=c(0.01,8))+scale_alpha_continuous(range=c(0.1,1)) #loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black" +scale_color_viridis_d(option = "B",begin = 0.2)
PCAPlot_targets2_PC23
ggsave(paste0(plotdir,"ALL.TargetsOnly.PCA.coloredBy3mer",orderLabel,"PC23.png"),PCAPlot_targets2_PC23,height=10,width=14)

PCAPlot_targets2_PC13 <- autoplot(pca,x=1,y=3,data=ancestralKmerTargets_spectrum_merged_withFeatures,colour="central3mer",size="fractionOfAllSegSitesInvolvingAncestralKmer",alpha="fractionOfAllSegSitesInvolvingAncestralKmer",label=T,label.size=0.1)+ggtitle(paste0("TARGETS ONLY; ", orderLabel,"\n",toString(desiredFeatures)))+theme_bw()+scale_size_continuous(range=c(0.01,8))+scale_alpha_continuous(range=c(0.1,1)) #loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black" +scale_color_viridis_d(option = "B",begin = 0.2)
PCAPlot_targets2_PC23
ggsave(paste0(plotdir,"ALL.TargetsOnly.PCA.coloredBy3mer",orderLabel,"PC13.png"),PCAPlot_targets2_PC13,height=10,width=14)

######## >>> PCA on derived mutations only ###################
head(spectrum_onepop)
derivedKmer_spectrum <- spectrum_onepop %>%
  group_by(derivedkmer) %>%
  summarise(totalMutationsToDerivedKmer_fromAnyAncestralTarget=sum(totalSitesAllIntervals))
derivedKmer_spectrum
# get the proportion of mutations that involve the derived kmer:
derivedKmer_spectrum$fractionOfAllSegSitesInvolvingDerivedKmer <- derivedKmer_spectrum$totalMutationsToDerivedKmer_fromAnyAncestralTarget/sum(derivedKmer_spectrum$totalMutationsToDerivedKmer_fromAnyAncestralTarget)

# combine this with features:
derivedKmer_spectrum_merged_withFeatures <- merge(derivedKmer_spectrum,allFeatureVectors_labelled,by.x="derivedkmer",by.y="motif")
rownames(derivedKmer_spectrum_merged_withFeatures) <- derivedKmer_spectrum_merged_withFeatures$derivedkmer

pca <- prcomp(select(derivedKmer_spectrum_merged_withFeatures,-c(names(derivedKmer_spectrum))), scale=T,center=T) # Michael says to scale/center; similar with and without

PCAPlot_derivedMutations <- autoplot(pca,data=derivedKmer_spectrum_merged_withFeatures,colour="fractionOfAllSegSitesInvolvingDerivedKmer",size="fractionOfAllSegSitesInvolvingDerivedKmer",alpha="fractionOfAllSegSitesInvolvingDerivedKmer")+ggtitle(paste0("Derived kmers ONLY; ", orderLabel,"\n",toString(desiredFeatures)))+theme_bw()+scale_size_continuous(range=c(0.01,8)) +scale_color_viridis_c(option = "B",begin = 0.1) +scale_alpha_continuous(range=c(0.1,1)) #loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black" +scale_color_viridis_d(option = "B",begin = 0.2)
PCAPlot_derivedMutations
ggsave(paste0(plotdir,"ALL.DerivedOnly.PCA.coloredByFracOfSegSites.",orderLabel,".png"),PCAPlot_derivedMutations,height=10,width=14)

# color by 3mer:
derivedKmer_spectrum_merged_withFeatures$central3mer <- substr(derivedKmer_spectrum_merged_withFeatures$derivedkmer,start=3,stop=5)

PCAPlot_derivedMutations2 <- autoplot(pca,data=derivedKmer_spectrum_merged_withFeatures,colour="central3mer",size="fractionOfAllSegSitesInvolvingDerivedKmer",alpha="fractionOfAllSegSitesInvolvingDerivedKmer")+ggtitle(paste0("Derived kmers ONLY; ", orderLabel,"\n",toString(desiredFeatures)))+theme_bw()+scale_size_continuous(range=c(0.01,8)) +scale_alpha_continuous(range=c(0.1,1)) #loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black" +scale_color_viridis_d(option = "B",begin = 0.2)
PCAPlot_derivedMutations2
ggsave(paste0(plotdir,"ALL.DerivedOnly.PCA.coloredBy3mer.",orderLabel,".png"),PCAPlot_derivedMutations2,height=10,width=14)

############## >>> to focus on TCC pulse --- a PCA of just CC>T mutations #############
# Kelley says to also add TCT>TTT because it's second most abundant. 
mutationSpectrum_featureMatrix_XCC <- mutationSpectrum_featureMatrix
mutationSpectrum_featureMatrix_XCC$ancestral3mer <- unlist(lapply(strsplit(mutationSpectrum_featureMatrix_XCC$ThreemerMutationType,"\\."),"[",1))
# only 3mers that contain CC in ancestral state and are C>T mutating?
mutationSpectrum_featureMatrix_XCC <- mutationSpectrum_featureMatrix_XCC[(grepl("CC",mutationSpectrum_featureMatrix_XCC$ancestral3mer) & mutationSpectrum_featureMatrix_XCC$centralMutationType=="C.T") | mutationSpectrum_featureMatrix_XCC$ancestral3mer=="TCT",]

pca <- prcomp(select(mutationSpectrum_featureMatrix_XCC,-c("ancestral3mer",names(spectrum_onepop))), scale=T,center=T) # Michael says to scale/center; similar with and without
PCAPlot_XCC <- autoplot(pca,data=mutationSpectrum_featureMatrix_XCC,colour="ThreemerMutationType")+theme_bw()+ggtitle("C>T Mutation 7mers that contain CC in ancestral 3mer, plus TCT>TTT (second most abundant)") # +scale_color_brewer(palette="Set2") #loadings=T,loadings.label=T,loadings.colour="black
PCAPlot_XCC
ggsave(paste0(plotdir,"TCCPulseComparisonPlusTCT.TTT.PCA.",orderLabel,".png"),PCAPlot_XCC,height=8,width=10)

############ try this just with ancestral / derived kmers related to TCC as well? ###########

########## want to try delta shape instead of one column for ancestral and one for derived ###########

##################### >>>  try UMAP  ###################################
# this is still v ery weird!!!
# set up umap config:
custom.config = umap.defaults
######## very slow!

# trying to subset UMAP to see why it's forming rows/columns
# what do rows/columsn have in common?
umapselection=c("TAA.TTA","CAA.CTA","GCC.GTC","GAA.GCA")
test_umap <- umap(select(mutationSpectrum_featureMatrix[mutationSpectrum_featureMatrix$ThreemerMutationType %in% umapselection,],-c("ancestral3mer",names(spectrum_onepop)))) # slow! 

test_umap_df <- data.frame(test_umap$layout)
colnames(test_umap_df) <- c("v1","v2")

head(test_umap_df)
test_umap_df$fractionOfAllSegSites <- mutationSpectrum_featureMatrix[mutationSpectrum_featureMatrix$ThreemerMutationType %in% umapselection,]$fractionOfAllSegSites
test_umap_df$ThreemerMutationType <- mutationSpectrum_featureMatrix[mutationSpectrum_featureMatrix$ThreemerMutationType  %in% umapselection,]$ThreemerMutationType
test_umap_df$variable <- mutationSpectrum_featureMatrix[mutationSpectrum_featureMatrix$ThreemerMutationType %in% umapselection,]$variable
test_umap_df$central5mer.ancestral <- mutationSpectrum_featureMatrix[mutationSpectrum_featureMatrix$ThreemerMutationType %in% umapselection,]$central5mer.ancestral

# using sample_n to randomly thin it out for now:
test_umap_plot <- ggplot(sample_n(test_umap_df,50),aes(x=v1,y=v2,colour=central5mer.ancestral,label=variable))+
  geom_point(shape=1,alpha=0.8)+
  geom_text(color="black")+
  ggtitle(toString(umapselection)) #+
  #theme(legend.position = "none")


test_umap_plot
ggsave(paste0(plotdir,"TEST.UMAP.",toString(umapselection),".",orderLabel,".png"),test_umap_plot,height=10,width=18)

############ >>> PCA colored by ratio between mouse species ##########
secondPop="Ms" # spretus
spectrum_secondpop = spectrum_multipop[spectrum_multipop$population==secondPop,c("variable","fractionOfAllSegSites","population")]
head(spectrum_secondpop)

mutationSpectrum_featureMatrix_twoPops <- merge(spectrum_secondpop,mutationSpectrum_featureMatrix,by="variable",suffixes=c(".pop2",".pop1")) # note order of populations here 
head(mutationSpectrum_featureMatrix_twoPops)

mutationSpectrum_featureMatrix_twoPops$Population_fracSegSitesRatio <- mutationSpectrum_featureMatrix_twoPops$fractionOfAllSegSites.pop1 / mutationSpectrum_featureMatrix_twoPops$fractionOfAllSegSites.pop2

rownames(mutationSpectrum_featureMatrix_twoPops) <- mutationSpectrum_featureMatrix_twoPops$variable

columnsToExclude=c("population.pop1","population.pop2","fractionOfAllSegSites.pop1","fractionOfAllSegSites.pop2","Population_fracSegSitesRatio",names(select(spectrum_onepop,-c("fractionOfAllSegSites","population"))))
#columnsToExclude


pca <- prcomp(select(mutationSpectrum_featureMatrix_twoPops,-all_of(columnsToExclude)), scale=T,center=T) 

# transform to log2 scale
mutationSpectrum_featureMatrix_twoPops$Population_fracSegSitesRatio_log2Scale <- log2(mutationSpectrum_featureMatrix_twoPops$Population_fracSegSitesRatio)

PCAPlot_ColoredByRatioOfTwoPops <- autoplot(pca,data=mutationSpectrum_featureMatrix_twoPops,colour="Population_fracSegSitesRatio_log2Scale",size=2,shape=1)+scale_color_gradient2() +ggtitle(paste0("Ratio of ",population,"/",secondPop," frac of SegSites"))
PCAPlot_ColoredByRatioOfTwoPops

ggsave(paste0(plotdir,population,".",secondPop,".PCA.coloredByRatioOfFracOfSegSites.",orderLabel,".PC12.png"),PCAPlot_ColoredByRatioOfTwoPops,height=10,width=14)


PCAPlot_ColoredByRatioOfTwoPops_23 <- autoplot(pca,x=2,y=3,data=mutationSpectrum_featureMatrix_twoPops,colour="Population_fracSegSitesRatio_log2Scale",size=2,shape=1)+scale_color_gradient2() +ggtitle(paste0("Ratio of ",population,"/",secondPop," frac of SegSites"))
PCAPlot_ColoredByRatioOfTwoPops_23
ggsave(paste0(plotdir,population,".",secondPop,".PCA.coloredByRatioOfFracOfSegSites.",orderLabel,".PC23.png"),PCAPlot_ColoredByRatioOfTwoPops_23,height=10,width=14)     

PCAPlot_ColoredByRatioOfTwoPops_13 <- autoplot(pca,x=1,y=3,data=mutationSpectrum_featureMatrix_twoPops,colour="Population_fracSegSitesRatio_log2Scale",size=2,shape=1)+scale_color_gradient2() +ggtitle(paste0("Ratio of ",population,"/",secondPop," frac of SegSites"))
PCAPlot_ColoredByRatioOfTwoPops_13
ggsave(paste0(plotdir,population,".",secondPop,".PCA.coloredByRatioOfFracOfSegSites.",orderLabel,".PC13.png"),PCAPlot_ColoredByRatioOfTwoPops_13,height=10,width=14) 

################ >>> trying to form clusters : hierarchical clustering #####################
# need a dissimilairty matrix from dist
?dist # This function computes and returns the distance matrix computed by using the specified distance measure to compute the distances between the rows of a data matrix. (default euclidean)
dist_mutationSpectrum_featureMatrix <- dist(select(mutationSpectrum_featureMatrix,-names(spectrum_onepop))) # slow
head(dist_mutationSpectrum_featureMatrix)
?hclust 
##### NEED to try other cluster methods
clusterCountk=12
hclust_results <- hclust(dist_mutationSpectrum_featureMatrix,method = "ward.D") # slow ; trying with Ward method?
dend <- as.dendrogram(hclust_results) # slow; from the dendextend package

dend <- color_branches(dend, k = clusterCountk)
labels(dend) <- ""
png(paste0(plotdir,"hclust.AllFeatures.ByShape.tree.ColoredByCluster.png"),units="cm",height=20,width=20,res = 300)
plot(dend)
dev.off()
#dend <- color_branches(dend,k=12) # this calls cutree with k = whatever to color the branches and assigns colors

# http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning#visualize-dendrogram-using-ggdendrogram-function <- nicer plotting of hclust dendrograms
# too slow # test_hclust_dendroPlot1 <- ggdendrogram(test_hclust,rotate=T,leaf_labels = F)

# can cut a tree into k groups using cutree 
?cutree
sample_clusters <- data.frame(groupID=as.factor(cutree(hclust_results,k = clusterCountk)))
head(sample_clusters)
sample_clusters$variable <- rownames(sample_clusters)
# merge with dataset: 
mutationSpectrum_featureMatrix_Clusters <- merge(sample_clusters,mutationSpectrum_featureMatrix,by.x="variable",by.y="variable")


# now do pca and color by group ID

pca <- prcomp(select(mutationSpectrum_featureMatrix_Clusters,-c("groupID",names(spectrum_onepop))), scale=T,center=T) # Michael says to scale/center; similar with and without
PCAPlot_ColoredByClusters <- autoplot(pca,data=mutationSpectrum_featureMatrix_Clusters,colour="groupID",size="fractionOfAllSegSites")+ggtitle(paste0(orderLabel,"\n",toString(desiredFeatures)))+theme_bw()+scale_size_continuous(limits=scaleSizeLimits,range=c(0.01,8))  +scale_alpha_continuous(range=c(0.1,1)) #+scale_color_brewer(palette="Set3") #loadings=T,loadings.label=T,loadings.colour="black",loadings.label.colour="black" +scale_color_viridis_d(option = "B",begin = 0.2)
PCAPlot_ColoredByClusters
#ggsave(paste0(plotdir,"ALL.7mers.PCA.coloredByFracOfSegSites.",orderLabel,".png"),PCAPlot_all7mers1,height=10,width=14)



ggplot(mutationSpectrum_featureMatrix_Clusters,aes(x=groupID,y=fractionOfAllSegSites,color=mutationClassLabel))+
  geom_point(size=2,alpha=0.8,shape=1)+
  scale_y_log10()


##########  >>> cluster shape features with heat map to form clusters? ########
rownames(mutationSpectrum_featureMatrix_Clusters) <- mutationSpectrum_featureMatrix_Clusters$variable
png(paste0(plotdir,"heatmap.shapefeatures.png"),height=20,width = 20,units="in",res=300)
heatmap(as.matrix(select(mutationSpectrum_featureMatrix_Clusters,-c("groupID",names(spectrum_onepop)))))
dev.off()

########## try two pops
mutationSpectrum_featureMatrix_twoPops_Clusters <- merge(sample_clusters, mutationSpectrum_featureMatrix_twoPops,by="variable")
head(mutationSpectrum_featureMatrix_twoPops_Clusters)

mutationSpectrum_featureMatrix_twoPops_Clusters$mutationClassLabel <- factor(mutationSpectrum_featureMatrix_twoPops_Clusters$mutationClassLabel,levels=c("A.C.no_ancestralCpG" , "A.G.no_ancestralCpG"  ,"A.T.no_ancestralCpG",  "C.A.no_ancestralCpG",  "C.G.no_ancestralCpG" , "C.T.no_ancestralCpG","C.A.yes_ancestralCpG", "C.G.yes_ancestralCpG", "C.T.yes_ancestralCpG" ))




boxplot1a <- ggplot(mutationSpectrum_featureMatrix_twoPops_Clusters,aes(x=groupID,color=groupID,y=Population_fracSegSitesRatio))+
  geom_boxplot(size=0.5)+
  geom_hline(yintercept=1,linetype="dotted")+
  
  #geom_point(size=0.8,alpha=0.8)+
  facet_wrap(~mutationClassLabel,dir = "v",ncol=1,scales="free_y")+
  theme_bw()+
  theme(text=element_text(size=14))
boxplot1a
ggsave(paste0(plotdir,population,".",secondPop,".",clusterCountk,"hclustgroups.boxplot.good.png"),boxplot1a,height=20,width=8)


# exclude cpgs: 
boxplot1b <- ggplot(mutationSpectrum_featureMatrix_twoPops_Clusters[!(mutationSpectrum_featureMatrix_twoPops_Clusters$mutationClassLabel %in% c("C.A.yes_ancestralCpG", "C.G.yes_ancestralCpG", "C.T.yes_ancestralCpG")),],aes(x=ThreemerMutationType,color=mutationClassLabel,y=Population_fracSegSitesRatio))+
  geom_boxplot()+
  geom_point(size=0.1,alpha=0.8)+
  facet_wrap(~groupID,scales="free")+
  geom_hline(yintercept=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("hclust groups -- no CpG categories")
boxplot1b
ggsave(paste0(plotdir,population,".",secondPop,".",clusterCountk,"hclustgroups.boxplot.PerGroup.Per3mer.good.png"),boxplot1b,height=14,width=24)




ggplot(mutationSpectrum_featureMatrix_twoPops_Clusters,aes(x=Population_fracSegSitesRatio,color=groupID))+
  geom_density()

#### look at 3mer composition of groups : 
total3merPerGroup <- mutationSpectrum_featureMatrix_twoPops_Clusters %>%
  group_by(groupID,ThreemerMutationType) %>%
  summarise(total_3mertype=n())
total3merPerGroup$frac_3mertype = total3merPerGroup$total_3mertype/sum(total3merPerGroup$total_3mertype)

ggplot(total3merPerGroup,aes(x=groupID,y=frac_3mertype,label=ThreemerMutationType,fill=ThreemerMutationType))+
  geom_col(position="stack")+
  geom_text(position = "stack")

######## >>> clustering: targets only ################
dist_mutationSpectrum_featureMatrix_ancOnly <- dist(select(ancestralKmerTargets_spectrum_merged_withFeatures,-c("ancestralkmer", "population","totalMutationsInvolvingAncestralKmer","fractionOfAllSegSitesInvolvingAncestralKmer"))) # slow

head(dist_mutationSpectrum_featureMatrix_ancOnly)
hclust_results_ancOnly <- hclust(dist_mutationSpectrum_featureMatrix_ancOnly,method = "ward.D") # slow ; trying with Ward method?
dend_ancOnly <- as.dendrogram(hclust_results_ancOnly) # slow; from the dendextend package
dend_ancOnly <- color_branches(dend_ancOnly,k=8) # this calls cutree with k = whatever to color the branches and assigns colors
labels(dend_ancOnly) <- ""
plot(dend_ancOnly)

