########## try looking at shapes of 5mers ################
require(tidyr)
require(dplyr)
require(ggplot2)
require(ggfortify)
require(reshape2)
wd="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/comparing_spectrum_and_phylo_divergence_WITHPROJECTIONS/20211216.plots.projectedCounts.epsilon.1.sepCpGs.no.addingApes_maskedMouseVaq.addingMultinomialDownsampling.MorePlotsForPAG.newTimeTree.MantelTest/"
plotdir=paste0(wd,"shapeExperiments/")
dir.create(plotdir)
spectra=read.table(paste0(wd,"AllProcessedSpectra.AllConditions.AllVariables.AllSpecies.txt"),header=T)
# note this has lots of conditions! don't double up!
head(spectra)
# subset to jsut 5mer spectra and not transformed and just the fully normalized mutation rate 
transform="none"
variable="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon"
spectra_5mer_notransform <- spectra[spectra$id=="spectra_5mer" & spectra$variable==variable & spectra$transform_id==transform,]
dim(spectra_5mer_notransform) # should be 1536 (correct)

###### add shape info to spectra table #######
# get ancestral/derived 
spectra_5mer_notransform$ancestral5mer <- unlist(lapply(strsplit(spectra_5mer_notransform$mutation_label,"\\."),"[",1))
spectra_5mer_notransform$derived5mer <- unlist(lapply(strsplit(spectra_5mer_notransform$mutation_label,"\\."),"[",2))


# get shapes:
shapes_5mers <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/DNAShapeR/firstOrder_featureTypes_allPossible5mers.FeatureValues.NOTNormalized.WillWorkForAnySpecies.notnormalized.UseForRandomForest.andTidyModels.5mer.txt",header=T)

######## merge with shape info #########
spectra_5mer_notransform_merge1 <- merge(spectra_5mer_notransform,shapes_5mers,by.x="ancestral5mer",by.y="motif")
spectra_5mer_notransform_merge2 <- merge(spectra_5mer_notransform_merge1,shapes_5mers,by.x="derived5mer",by.y="motif",suffixes = c(".ancestral",".derived"))

spectra_5mer_notransform_shapeInfo <- spectra_5mer_notransform_merge2
############### pca on shapes ###############
head(shapes_5mers)
shapePCA <- prcomp(select(shapes_5mers,-motif),center = T,scale=T)
rownames(shapes_5mers) <- shapes_5mers$motif
shapePCAplot <- autoplot(shapePCA,loadings=T,loadings.label=T,label=T,data=shapes_5mers,label.size=1,shape=F) + theme_bw()   +
  ggtitle("PCA based on shape characteristics\nreverse complements lead to symmetry\ncenter=T, scale=T")
shapePCAplot          
ggsave(paste0(plotdir,"shape.pca.5mer.symmetryDueToRevComps.png"),shapePCAplot)

####### pca just on ancestral shapes (rev comps collapsed) ########
shapes_5mers_collapsedRevComps <- shapes_5mers[shapes_5mers$motif %in% unique(spectra_5mer_notransform_shapeInfo$ancestral5mer),] # just pick motifs that are in the ancestral set so rev comps aren't included ; note that 1536 is total paired mutation types X>Y so there are fewer "X"s (rev comped); should be 1/2 of total size. (1024 motifs --> 512)
shapePCA_revCompCollapsed <- prcomp(select(shapes_5mers_collapsedRevComps,-motif),center=T,scale=T)
# get for mutation types #
shapePCAplot_RevCompsExcluded <- autoplot(shapePCA_revCompCollapsed,loadings=T,loadings.label=T,label=T,data=shapes_5mers_collapsedRevComps,label.size=1,shape=F) + theme_bw()   +
  ggtitle("PCA based on shape characteristics\nno rev comps\ncenter=T, scale=T")
shapePCAplot_RevCompsExcluded          
ggsave(paste0(plotdir,"shape.pca.5mer.CollapseRevComps.png"),shapePCAplot_RevCompsExcluded)

############## pca on all mutation types ################
rownames(spectra_5mer_notransform_shapeInfo) <- spectra_5mer_notransform_shapeInfo$mutation_label

p3 <- prcomp(select(spectra_5mer_notransform_shapeInfo,starts_with('feature_')),center=T,scale=T)

allMutTypesPCA_mouse <- autoplot(p3,data=spectra_5mer_notransform_shapeInfo,label=T,label.size=1,shape=T,loadings=T,loadings.label=T,colour="mice_Mmd",size="mice_Mmd")+
  theme_bw()+
  ggtitle("Colored by mutatbility in mice")+
  scale_color_viridis_c()

allMutTypesPCA_mouse
ggsave(paste0(plotdir,"allMutTypesPCA.mutabilityinmouse.png"),allMutTypesPCA)


allMutTypesPCA_bear <- autoplot(p3,data=spectra_5mer_notransform_shapeInfo,label=T,label.size=1,shape=T,loadings=T,loadings.label=T,colour="bears_ABC",size="bears_ABC")+
  theme_bw()+
  ggtitle("Colored by mutatbility in bears_ABC")+
  scale_color_viridis_c()

allMutTypesPCA_bear
ggsave(paste0(plotdir,"allMutTypesPCA.mutabilityinbear.png"),allMutTypesPCA_bear)


######### try to facet by central type ######### 
spectra_5mer_notransform_shapeInfo$centralBP <- paste0(substr(spectra_5mer_notransform_shapeInfo$mutation_label,3,3),".",substr(spectra_5mer_notransform_shapeInfo$mutation_label,9,9))
# add cpgs in 
for(mutationType in c("A.C","A.T","A.G","C.A","C.T","C.G")){
  subset= spectra_5mer_notransform_shapeInfo[spectra_5mer_notransform_shapeInfo$centralBP==mutationType,]
  pca_specificType <- prcomp(select(subset,starts_with('feature_')),center=T,scale=T)
  specificPCA <-  autoplot(pca_specificType,loadings=T,loadings.label=T,label=T,data=subset,label.size="mice_Mmd",shape=F,size="mice_Mmd",colour="mice_Mmd",loadings.label.size=2) + theme_bw()   +
    ggtitle(paste0(mutationType,"mice_Mmd"))+
    scale_color_viridis_c()
  specificPCA
  ggsave(paste0(plotdir,mutationType,"pca.mice_Mmd.png"),specificPCA)
}
########## keep exploring this #############
# would like to find that a) the most and least mutable 5mers across species have shape features in common? b) the most variable across species have shapes in common? c) pairs of species have differences? signatures? 
head(spectra_5mer_notransform_shapeInfo)
spectra_5mer_notransform_shapeInfo_melt <- melt(spectra_5mer_notransform_shapeInfo,id.vars = c("mutation_label","ancestral5mer","derived5mer","id","variable","transform_id",names(spectra_5mer_notransform_shapeInfo)[ grep("feature",names(spectra_5mer_notransform_shapeInfo))]),variable.name = "species",value.name = "mutrate")

head(spectra_5mer_notransform_shapeInfo_melt)                
buckleplot <- ggplot(spectra_5mer_notransform_shapeInfo_melt,aes(x=feature_1_Buckle_1.ancestral,y=mutrate,color=species))+
  geom_point()
buckleplot
ggsave(paste0(plotdir,"example.BuckleVsPerSpecies.png"),buckleplot)
require(matrixStats) # for rowVars
species=c( "bears_EUR", "bears_ABC","bears_PB","fin_whale_GOC","fin_whale_ENP" , "apes_Gorilla"  ,"humans_AMR","humans_EUR", "humans_EAS", "humans_AFR" , "humans_SAS" ,"mice_Ms" ,"mice_Mmd"      ,"mice_Mmm" , "mice_Mmc" , "apes_Pan_paniscus","apes_Pan_troglodytes","apes_Pongo_abelii"   , "apes_Pongo_pygmaeus" , "vaquita" )
spectra_5mer_notransform_shapeInfo$variance <- rowVars(as.matrix(select(spectra_5mer_notransform_shapeInfo,c(species) ))) # 
### something wrong here 
variance <- spectra_5mer_notransform_shapeInfo_melt %>%
  group_by(mutation_label) %>%
  summarise(varianceAcrossSpecies=var(mutrate))
varianceHistogram <- ggplot(variance,aes(x=varianceAcrossSpecies))+
  geom_histogram()+
  scale_x_log10()+
  geom_text(data=filter(variance,varianceAcrossSpecies>5e-6),aes(x=varianceAcrossSpecies,y=50,label=mutation_label),size=1)
varianceHistogram
ggsave(paste0(plotdir,"varianceHistogram.CGrichraremotifshavehighvar.png"),varianceHistogram)
# the variance depends on the scale , so most mutable may have most variance? figure out
# also get highest and lowest ? and look for similarity in shapes? how? 
# what about another distance metric? distance between shapes? could calculate distance between shapes and get the most different and most similar sequences based on shape? 
# somehow need to condense information. don't need modeling.
# want to find spikes . maybe like shap values? maybe for each feature, plot feature values and y axis is mut rate across all species? avg? all pts? variance? 
spectra_5mer_notransform_shapeInfo_variance <- merge(spectra_5mer_notransform_shapeInfo,variance,by="mutation_label")
allMutTypesPCA_colorByVariance <- autoplot(p3,data=spectra_5mer_notransform_shapeInfo_variance,label=T,shape=F,loadings=F,loadings.label=F,label.size="variance",label.colour="variance",colour="variance",size="variance")+
  theme_bw()+
  scale_color_viridis_c()

allMutTypesPCA_colorByVariance
ggsave(paste0(plotdir,"allMutTypesPCA_colorByVariance.png"),allMutTypesPCA_colorByVariance)
