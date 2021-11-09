######## trying to plot actual mutation types in predicted models from sig profiler extractor ##########

require(dplyr)
require(tidyverse)
require(reshape2)
wd="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/SigProfilerExtractor/20210809_SigProfilerExtractorResults_bears_kSFSAllFreqsSep_GRCh37GenomeRef/SBS96/"
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

### okay so for each sample, I want to mulitply the signature by total contributed sites for each signature somehow...
# merge them ? per sample and per signature?

head(cosmic_activities)
head(cosmic_signatures)
# trying something:
cosmic_signatures_melt <- melt(cosmic_signatures)
colnames(cosmic_signatures_melt) <- c("MutationType","signature","fractionOfSignature")

cosmic_activities_melt <- melt(cosmic_activities)
colnames(cosmic_activities_melt) <- c("Samples","signature","totalSitesOfAllTypesContibutedBySignatureToSample")

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

cosmic_cominbed_summed <- cosmic_combined %>%
  group_by(centralMutationType,frequency,species,signature) %>%
  summarise(totalMutatationsOfThisTypeContributedBySignature=sum(totalSites_specificType_ContributedBySignature)) 

# get proportions # 
cosmic_cominbed_summed <- cosmic_cominbed_summed %>%
  group_by(frequency,species) %>%
  mutate(proportionOfSites = totalMutatationsOfThisTypeContributedBySignature/sum(totalMutatationsOfThisTypeContributedBySignature))

# check:
sum(cosmic_cominbed_summed[cosmic_cominbed_summed$frequency==1 & cosmic_cominbed_summed$species=="ABC",]$totalMutatationsOfThisTypeContributedBySignature) 
sum(empirical$ABC_1) # nice these match (within 1 site)
sum(cosmic_cominbed_summed[cosmic_cominbed_summed$frequency==1 & cosmic_cominbed_summed$species=="ABC",]$proportionOfSites)  #1 


plot1 <- ggplot(cosmic_cominbed_summed[cosmic_cominbed_summed$signature %in% c("SBS96A","SBS96B"),],aes(x=frequency,y=proportionOfSites,color=species,group=species))+
  facet_wrap(~centralMutationType~signature,scales="free",ncol=2)+
  geom_line()+
  theme_bw()
plot1
ggsave(paste0(plotdir,"model.MutationsContributedPerSignature.freqsAsCounts.png"),plot1,height=12,width=5)
# want to get total (regardless of signature)

cosmic_cominbed_summed <- cosmic_cominbed_summed %>%
  group_by(species) %>%
  mutate(alleleFrequency= frequency/(max(frequency)+1)) # have to add one because fixed sites 16/16 isn't a category but AF is 1/16 not 1/15 <- 
# View(cosmic_cominbed_summed)
plot2 <- ggplot(cosmic_cominbed_summed[cosmic_cominbed_summed$signature %in% c("SBS96A","SBS96B"),],aes(x=alleleFrequency,y=proportionOfSites,color=species,group=species))+
  facet_wrap(~centralMutationType~signature,scales="free", ncol=2)+
  geom_line()+
  theme_bw()
plot2
ggsave(paste0(plotdir,"model.MutationsContributedPerSignature.alleleFreqs.png"),plot2,height=12,width=5)

######## now group overall and see how well it explains this #####

cosmic_cominbed_summed_OverSignatures <- cosmic_cominbed_summed %>%
  group_by(species,frequency,centralMutationType) %>%
  summarise(totalSites_acrossSigs = sum(totalMutatationsOfThisTypeContributedBySignature))

cosmic_cominbed_summed_OverSignatures <- cosmic_cominbed_summed_OverSignatures %>%
  group_by(species,frequency) %>%
  mutate(proportionOfSites=totalSites_acrossSigs/sum(totalSites_acrossSigs))
  
head(cosmic_cominbed_summed_OverSignatures)

# get allele freqs:
cosmic_cominbed_summed_OverSignatures <- cosmic_cominbed_summed_OverSignatures %>%
  group_by(species) %>%
  mutate(alleleFrequency= frequency/(max(frequency)+1)) # have to add one because fixed sites 16/16 isn't a category but AF is 1/16 not 1/15 <- 
# View(cosmic_cominbed_summed)
head(cosmic_cominbed_summed_OverSignatures)


plot3 <- ggplot(cosmic_cominbed_summed_OverSignatures,aes(x=alleleFrequency,y=proportionOfSites,color=species,group=species))+
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
plot4 <- ggplot(cosmic_cominbed_summed_OverSignatures,aes(x=alleleFrequency,y=proportionOfSites,color=species,group=species,linetype="sigprofilerextractor"))+
  facet_wrap(~centralMutationType,scales="free") +
  geom_line()+
  geom_line(data=empirical_melt_summed,aes(x=alleleFrequency,y=proportionOfSites,color=species,group=species,linetype="empirical"))+
  ggtitle("all signatures combined, with empirical")+
  theme_bw()
plot4
ggsave(paste0(plotdir,"model.MutationsContributedPerCentralType.allSigsCombined.alleleFreqs.WITHEMPIRICAL.png"),plot4,width=12,height=7)

# I was worried that these results looked slightly different than what I remembered -- but actually it's just the use of nostrict or strict that makes a difference 
sum(cosmic_cominbed_summed_OverSignatures[cosmic_cominbed_summed_OverSignatures$species=="PB" & cosmic_cominbed_summed_OverSignatures$frequency==1,]$proportionOfSites)

#empirical_melt_summed %>%
#  group_by(species,frequency) %>%
#  summarise(tot=sum(proportionOfSites)) %>%
#  View()
# sums to one correctly, taht's good
############ SANDBOX: what the heck is going on with neutral sites that they look kind of different for me from when I did this as 1mers for BGC?? ###########
### first troubleshooting thing: see if this is something weird that sigprofilerextractor did #####
# origdata=read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210401_3mer_nostrict/mapped_to_brown_bear/summaryTables/All_Intervals_kSFS.ABC.NonABC.NotSeparatedYet.brown_bear.ALLFREQS.nostrict.txt",header=T)
# head(origdata)
# 
# origdata$centralMutationType <- paste0(substr(origdata$variable,2,2),">",substr(origdata$variable,6,6))
# head(origdata)
# tail(origdata)
# 
# origdata_summed <- origdata %>%
#   group_by(population,sample_frequency,centralMutationType) %>%
#   summarise(totalSites=sum(totalSites))
# head(origdata_summed)  
# origdata_summed <- origdata_summed %>% group_by(population,sample_frequency) %>%
#   mutate(proportionOfSites=totalSites/sum(totalSites)) 
# # get allele freq
# origdata_summed <- origdata_summed %>%
#   group_by(population) %>%
#   mutate(alleleFrequency= sample_frequency/(max(sample_frequency)+1)) # have to add one because fixed sites 16/16 isn't a category but AF is 1/16 not 1/15 <- 
# head(origdata_summed)
# # rename population to species:
# colnames(origdata_summed)[1] <- "species"
# head(origdata_summed)
# 
# ### add this original this to the plot:
# 
# sandbox_plot5  <- ggplot(cosmic_cominbed_summed_OverSignatures,aes(x=alleleFrequency,y=proportionOfSites,color=species,group=species,linetype="sigprofilerextractor"))+
#   facet_wrap(~centralMutationType,scales="free") +
#   geom_line()+
#   geom_line(data=empirical_melt_summed,aes(x=alleleFrequency,y=proportionOfSites,color=species,group=species,linetype="empirical"))+
#   geom_line(data=origdata_summed,aes(x=alleleFrequency,y=proportionOfSites,color=species,linetype="origdata_3mer"))+
#   ggtitle("all signatures combined, with empirical")
# sandbox_plot5
# 
# ggsave(paste0(plotdir,"model.MutationsContributedPerCentralType.allSigsCombined.alleleFreqs.WITHEMPIRICAL.withOrigDataEmpiricalToo.png"),sandbox_plot5,width=12,height=7)

####### okay now let's try with 1mer BGC ######
