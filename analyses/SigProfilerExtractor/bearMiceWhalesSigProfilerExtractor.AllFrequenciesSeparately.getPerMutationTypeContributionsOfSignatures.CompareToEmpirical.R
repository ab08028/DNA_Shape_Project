######## trying to plot actual mutation types in predicted models from sig profiler extractor ##########

require(dplyr)
require(tidyverse)
require(reshape2)
wd="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/SigProfilerExtractor/20210825_SigProfilerExtractorResults_bearsMiceVaquitaFinWhales_kSFSAllFreqsSep_mm10GenomeRef/SBS96/"
plotdir=paste0(wd,"myplots/")
dir.create(plotdir,showWarnings = F)
# not precisely sure how to do this, but have a few pieces:
# activities: these are the number of sites per individual per signature; look like: SBS96_S2_NMF_Activities.txt
#Samples	SBS96A	SBS96B
#ABC	407948	467188
#EUR	910756	937447
#PB	255031	120660

# then you also have a mutation probabilities file:
# SBS96_S2_Signatures.txt
# MutationType	SBS96A	SBS96B
# T[T>G]T	0.007084022636525333	0.009425224520266055
# T[T>C]T	0.012327333338558673	0.02368069513887167
# T[T>A]T	0.008128965561278165	0.004811476237606257
# G[T>G]T	0.003798518561758101	0.0052798039438202975
# G[T>C]T	0.013695589473471045	0.016363554568961262
# G[T>A]T	0.002367109747370705	0.002344377912348136
# C[T>G]T	0.004684840297792107	0.007032036006450653
# C[T>C]T	0.008242251921445131	0.019093077009543775
# C[T>A]T	0.003055387395899743	0.004241104795131833
# which gives proportion of each mutation type within each signature 
# So I think I could mulitply each column of each signature by the totals for each species to get counts

### and you ahve your actual data in the Samples.txt dir

# there's one more file I don't know if I need or not which shows probability of each mutation type coming from each signature:
# De_Novo_Mutation_Probabilities_refit.txt
# Sample Names	MutationTypes	SBS96A	SBS96B
# ABC	T[T>G]T	0.3962440692018421	0.6037559307981579
# ABC	T[T>C]T	0.3125051957708169	0.6874948042291832
# ABC	T[T>A]T	0.5960028679840764	0.40399713201592363
# ABC	G[T>G]T	0.3858311960272073	0.6141688039727926
# ABC	G[T>C]T	0.4222424197914	0.5777575802085999
# ABC	G[T>A]T	0.4685559464817049	0.5314440535182952
# ABC	C[T>G]T	0.36778372423712363	0.6322162757628763

# ^^^ what is the point of this one? understanding odds of an individual getting aha so it's the "and the probability for each signature to cause a specific mutation type in a cancer sample."  so it says that for a given sample it's 60% TTG>T is from SBS96B . okay so I don't need that for this.

############### read in files ##########3
########## this is the actual data: #########
empirical = read_tsv(paste0(wd,"Samples.txt"),col_names=T)
head(empirical)

#######3 have both cosmic and de novo activities (total sites contributed)  ########
cosmic_activities = read_tsv(paste0(wd,"Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities_refit.txt"),col_names=T)
head(cosmic_activities)

denovo_activities = read_tsv(paste0(wd,"Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt"),col_names=T)
head(denovo_activities)

####### and then get the make up of the signatures ######
denovo_signatures = read_tsv(paste0(wd,"Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt"),col_names=T)
head(denovo_signatures)

cosmic_signatures = read_tsv(paste0(wd,"Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Signatures/COSMIC_SBS96_Signatures.txt"),col_names=T)
head(cosmic_signatures)

### okay so for each sample, I want to mulitply the signature by total contributed sites for each signature
#### for now just focusing on the cosmic activities because any sites that can't be attributed to a cosmic are represented as the de novos, so it still represents them all ### 


head(cosmic_activities)
head(cosmic_signatures)
# trying something:
cosmic_signatures_melt <- melt(cosmic_signatures)
colnames(cosmic_signatures_melt) <- c("MutationType","signature","fractionOfSignature")
denovo_signatures_melt <- melt(denovo_signatures)
colnames(denovo_signatures_melt) <- c("MutationType","signature","fractionOfSignature")

cosmic_activities_melt <- melt(cosmic_activities)
colnames(cosmic_activities_melt) <- c("Samples","signature","totalSitesOfAllTypesContibutedBySignatureToSample")

denovo_activities_melt <- melt(denovo_activities)
colnames(denovo_activities_melt) <- c("Samples","signature","totalSitesOfAllTypesContibutedBySignatureToSample")

head(cosmic_signatures_melt)
head(cosmic_activities_melt)

cosmic_combined <- merge(cosmic_signatures_melt,cosmic_activities_melt)
head(cosmic_combined)

cosmic_combined$totalSites_specificType_ContributedBySignature <- cosmic_combined$fractionOfSignature*cosmic_combined$totalSitesOfAllTypesContibutedBySignatureToSample

head(cosmic_combined)
tail(cosmic_combined)

cosmic_combined$centralMutationType <- substr(cosmic_combined$MutationType,3,5)
head(cosmic_combined)

cosmic_combined$frequency <- as.numeric(unlist(lapply(strsplit(cosmic_combined$Samples,"_"),"[",2)))
head(cosmic_combined)

cosmic_combined$species <- unlist(lapply(strsplit(cosmic_combined$Samples,"_"),"[",1))
head(cosmic_combined)

# summing 3mers over central mutation type:
cosmic_combined_summed <- cosmic_combined %>%
  group_by(centralMutationType,frequency,species,signature) %>%
  summarise(totalMutatationsOfThisTypeContributedBySignature=sum(totalSites_specificType_ContributedBySignature)) 

# get proportions # 
cosmic_combined_summed <- cosmic_combined_summed %>%
  group_by(frequency,species) %>%
  mutate(proportionOfSites = totalMutatationsOfThisTypeContributedBySignature/sum(totalMutatationsOfThisTypeContributedBySignature))

# check:
sum(cosmic_combined_summed[cosmic_combined_summed$frequency==1 & cosmic_combined_summed$species=="ABC",]$totalMutatationsOfThisTypeContributedBySignature) 
sum(empirical$ABC_1) # nice these match (within 1 site)
sum(cosmic_combined_summed[cosmic_combined_summed$frequency==1 & cosmic_combined_summed$species=="ABC",]$proportionOfSites)  #1 


plot1 <- ggplot(cosmic_combined_summed,aes(x=frequency,y=proportionOfSites,color=species,group=species))+
  facet_wrap(~centralMutationType~signature,scales="free",ncol=4)+
  geom_line()+
  theme_bw()
plot1
ggsave(paste0(plotdir,"model.MutationsContributedPerSignature.freqsAsCounts.png"),plot1,height=30,width=30)
# want to get total (regardless of signature)

cosmic_combined_summed <- cosmic_combined_summed %>%
  group_by(species) %>%
  mutate(alleleFrequency= frequency/(max(frequency)+1)) # have to add one because fixed sites 16/16 isn't a category but AF is 1/16 not 1/15 <- 
# View(cosmic_combined_summed)
plot2 <- ggplot(cosmic_combined_summed,aes(x=alleleFrequency,y=proportionOfSites,color=species,group=species))+
  facet_wrap(~centralMutationType~signature,scales="free")+
  geom_line()+
  theme_bw()
plot2
ggsave(paste0(plotdir,"model.MutationsContributedPerSignature.alleleFreqs.png"),plot2,height=30,width=20)

######## now group overall and see how well it explains this #####

cosmic_combined_summed_OverSignatures <- cosmic_combined_summed %>%
  group_by(species,frequency,centralMutationType) %>%
  summarise(totalSites_acrossSigs = sum(totalMutatationsOfThisTypeContributedBySignature))

cosmic_combined_summed_OverSignatures <- cosmic_combined_summed_OverSignatures %>%
  group_by(species,frequency) %>%
  mutate(proportionOfSites=totalSites_acrossSigs/sum(totalSites_acrossSigs))


# get allele freqs:
cosmic_combined_summed_OverSignatures <- cosmic_combined_summed_OverSignatures %>%
  group_by(species) %>%
  mutate(alleleFrequency= frequency/(max(frequency)+1)) # have to add one because fixed sites 16/16 isn't a category but AF is 1/16 not 1/15 <- 
# View(cosmic_combined_summed)
head(cosmic_combined_summed_OverSignatures)


plot3 <- ggplot(cosmic_combined_summed_OverSignatures,aes(x=alleleFrequency,y=proportionOfSites,color=species,group=species))+
  facet_wrap(~centralMutationType,scales="free") +
  geom_line()+
  ggtitle("all signatures combined")+
  theme_bw()
plot3
ggsave(paste0(plotdir,"model.MutationsContributedPerCentralType.allSigsCombined.alleleFreqs.png"),plot3,width=12,height=7)
####### want to add empirical ######

empirical_melt <- melt(empirical)
head(empirical_melt)
empirical_melt$centralMutationType <- substr(empirical_melt$`Mutation Types`,3,5)
head(empirical_melt)
empirical_melt$frequency <- as.numeric(unlist(lapply(strsplit(as.character(empirical_melt$variable),"_"),"[",2)))
empirical_melt$species <- unlist(lapply(strsplit(as.character(empirical_melt$variable),"_"),"[",1))
head(empirical_melt)

# get proportion of sites:
empirical_melt_summed <- empirical_melt %>% 
  group_by(frequency,species,centralMutationType) %>% 
  summarise(totalSites=sum(value))

# get proportions
empirical_melt_summed <- empirical_melt_summed %>%
  group_by(species,frequency) %>%
  mutate(proportionOfSites =totalSites / sum(totalSites) )

# get alelle freqs
empirical_melt_summed <- empirical_melt_summed %>%
  group_by(species) %>%
  mutate(alleleFrequency= frequency/(max(frequency)+1)) 

head(empirical_melt_summed) 

###### plot with empirical ######
plot4 <- ggplot(cosmic_combined_summed_OverSignatures,aes(x=alleleFrequency,y=proportionOfSites,color=species,group=species,linetype="sigprofilerextractor"))+
  facet_wrap(~centralMutationType,scales="free") +
  geom_line()+
  geom_line(data=empirical_melt_summed,aes(x=alleleFrequency,y=proportionOfSites,color=species,group=species,linetype="empirical"))+
  ggtitle("all signatures combined, with empirical")+
  theme_bw()
plot4
ggsave(paste0(plotdir,"model.MutationsContributedPerCentralType.allSigsCombined.alleleFreqs.WITHEMPIRICAL.png"),plot4,width=12,height=7)


plot5 <- ggplot(cosmic_combined_summed_OverSignatures,aes(x=alleleFrequency,y=proportionOfSites,color=species,group=species,linetype="sigprofilerextractor"))+
  facet_wrap(~centralMutationType~species,scales="free") +
  geom_line()+
  geom_line(data=empirical_melt_summed,aes(x=alleleFrequency,y=proportionOfSites,color=species,group=species,linetype="empirical"))+
  ggtitle("all signatures combined, with empirical")+
  theme_bw()
plot5
ggsave(paste0(plotdir,"model.MutationsContributedPerCentralType.allSigsCombined.alleleFreqs.WITHEMPIRICAL.facetedPerSpecies.png"),plot5,width=16,height=16)

# I was worried that these results looked slightly different than what I remembered -- but actually it's just the use of nostrict or strict that makes a difference 
sum(cosmic_combined_summed_OverSignatures[cosmic_combined_summed_OverSignatures$species=="PB" & cosmic_combined_summed_OverSignatures$frequency==1,]$proportionOfSites)


######### make de novo signature plots #####
denovo_combined <- merge(denovo_signatures_melt,denovo_activities_melt)
head(denovo_combined)

denovo_combined$totalSites_specificType_ContributedBySignature <- denovo_combined$fractionOfSignature*denovo_combined$totalSitesOfAllTypesContibutedBySignatureToSample

head(denovo_combined)
tail(denovo_combined)

denovo_combined$centralMutationType <- substr(denovo_combined$MutationType,3,5)
head(denovo_combined)

denovo_combined$frequency <- as.numeric(unlist(lapply(strsplit(denovo_combined$Samples,"_"),"[",2)))
head(denovo_combined)

denovo_combined$species <- unlist(lapply(strsplit(denovo_combined$Samples,"_"),"[",1))
head(denovo_combined)

denovo_combined_summed <- denovo_combined %>%
  group_by(centralMutationType,frequency,species,signature) %>%
  summarise(totalMutatationsOfThisTypeContributedBySignature=sum(totalSites_specificType_ContributedBySignature)) 

# get proportions # 
denovo_combined_summed <- denovo_combined_summed %>%
  group_by(frequency,species) %>%
  mutate(proportionOfSites = totalMutatationsOfThisTypeContributedBySignature/sum(totalMutatationsOfThisTypeContributedBySignature))




denovo_combined_summed <- denovo_combined_summed %>%
  group_by(species) %>%
  mutate(alleleFrequency= frequency/(max(frequency)+1)) # have to add one because fixed sites 16/16 isn't a category but AF is 1/16 not 1/15 <- 

######## now group overall and see how well it explains this #####

denovo_combined_summed_OverSignatures <- denovo_combined_summed %>%
  group_by(species,frequency,centralMutationType) %>%
  summarise(totalSites_acrossSigs = sum(totalMutatationsOfThisTypeContributedBySignature))

denovo_combined_summed_OverSignatures <- denovo_combined_summed_OverSignatures %>%
  group_by(species,frequency) %>%
  mutate(proportionOfSites=totalSites_acrossSigs/sum(totalSites_acrossSigs))

head(denovo_combined_summed_OverSignatures)

# get allele freqs:
denovo_combined_summed_OverSignatures <- denovo_combined_summed_OverSignatures %>%
  group_by(species) %>%
  mutate(alleleFrequency= frequency/(max(frequency)+1)) # have to add one because fixed sites 16/16 isn't a category but AF is 1/16 not 1/15 <- 
# View(cosmic_combined_summed)
head(denovo_combined_summed_OverSignatures)


########## plot with de novo and empirical together ##########
plot4_denovo <- ggplot(denovo_combined_summed_OverSignatures,aes(x=alleleFrequency,y=proportionOfSites,color=species,group=species,linetype="sigprofilerextractor"))+
  facet_wrap(~centralMutationType,scales="free") +
  geom_line()+
  geom_line(data=empirical_melt_summed,aes(x=alleleFrequency,y=proportionOfSites,color=species,group=species,linetype="empirical"))+
  ggtitle("all signatures combined, with empirical")+
  theme_bw()
plot4_denovo
ggsave(paste0(plotdir,"model.MutationsContributedPerCentralType.allSigsCombined.alleleFreqs.WITHEMPIRICAL.DeNovoOnly.NoCosmicDecomposition.png"),plot4_denovo,width=30,height=20)

plot5_denovo <- ggplot(denovo_combined_summed_OverSignatures,aes(x=alleleFrequency,y=proportionOfSites,color=species,group=species,linetype="sigprofilerextractor"))+
  facet_wrap(~centralMutationType~species,scales="free") +
  geom_line()+
  geom_line(data=empirical_melt_summed,aes(x=alleleFrequency,y=proportionOfSites,color=species,group=species,linetype="empirical"))+
  ggtitle("all signatures combined, with empirical")+
  theme_bw()
plot5_denovo
ggsave(paste0(plotdir,"model.MutationsContributedPerCentralType.allSigsCombined.alleleFreqs.WITHEMPIRICAL.DeNovoOnly.NoCosmicDecomposition.facetedPerSpp.png"),plot5_denovo,width=30,height=16)

########## TROUBLESHOOTING: YOU ARE HERE look into what's going on with decomposition ##########
head(cosmic_activities)
head(denovo_activities)
head(cosmic_activities_melt)

totalPredSitesPerSample_cosmic <- cosmic_combined_summed_OverSignatures %>%
  group_by(species,frequency) %>%
  summarise(total=sum(totalSites_acrossSigs))

head(totalPredSitesPerSample_cosmic)

totalPredSitesPerSample_denovo <- denovo_combined_summed_OverSignatures %>%
  group_by(species,frequency) %>%
  summarise(total=sum(totalSites_acrossSigs))

head(totalPredSitesPerSample_denovo)
# okay so both methods are yielding the same total amount of sites predicted. so what differs is how they break down per category?
# try to plot another way?

# plot all together?
ggplot(denovo_combined_summed_OverSignatures[denovo_combined_summed_OverSignatures$species %in% "Mmd",],aes(x=frequency,y=totalSites_acrossSigs,fill="denovo prediction",color="denovo prediction"))+
  geom_point()+
  facet_wrap(~species~centralMutationType,scales="free")+
  geom_point(data=cosmic_combined_summed_OverSignatures[cosmic_combined_summed_OverSignatures$species %in% "Mmd",],aes(x=frequency,y=totalSites_acrossSigs,fill="cosmic prediction",color="cosmic prediction"))+
  geom_point(data=empirical_melt_summed[empirical_melt_summed$species %in% "Mmd",],aes(x=frequency,y=totalSites,fill="empirical",color="empirical"))+
  scale_y_log10()


######### YOU ARE HERE: TRYING TO UNDERSTAND WHY CONTRIBUTIONS DIFFER SO MUCH ########### 
ggplot(denovo_combined_summed_OverSignatures[denovo_combined_summed_OverSignatures$species %in% "Mmd",],aes(x=frequency,y=proportionOfSites,fill="denovo prediction",color="denovo prediction"))+
  geom_point()+
  facet_wrap(~species~centralMutationType,scales="free")+
  geom_point(data=cosmic_combined_summed_OverSignatures[cosmic_combined_summed_OverSignatures$species %in% "Mmd",],aes(x=frequency,y=proportionOfSites,fill="cosmic prediction",color="cosmic prediction"))+
  geom_point(data=empirical_melt_summed[empirical_melt_summed$species %in% "Mmd",],aes(x=frequency,y=proportionOfSites,fill="empirical",color="empirical"))

sum(cosmic_combined_summed[cosmic_combined_summed$frequency==1 & cosmic_combined_summed$species=="Mmd",]$totalMutatationsOfThisTypeContributedBySignature)
# 4141244
sum(empirical$Mmd_1)
# 4141244
test = cosmic_combined_summed[cosmic_combined_summed$frequency==1 & cosmic_combined_summed$species=="Mmd",]

test$newPropCalc <- test$totalMutatationsOfThisTypeContributedBySignature / 4141244

# okay so it's getting total sites right.
# why are proportions so odd? I think it's because cosmic signatures don't fit well.
