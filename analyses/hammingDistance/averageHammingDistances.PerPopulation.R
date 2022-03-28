####### want total targets and want population information #########

require(viridis)
require(ggplot2)
require(dplyr)
require(reshape2)
require(scales)
# human pop files #
popfiledir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/populationFiles/"
# note not doing hamming for vaquita; adding apes on 20211206
hammingindir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/hamming_distance/"
# update to make sure this is correct !! targets must be masked in same way as hamming distances
# so if rep masking occurred for hamming dist, rep masking must occur for targets (can stay consisten by using targets and variants files from same mutyper runs)
### note th
tableOfTargetsToInclude=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/spectra/tableOfSpectraToIncludeInComparisons.SpeciesWithMultipleIndividualsOnly.PlusApes.MaskedMouse.Vaquita.20211123.txt",header=T,sep="\t") 
labels=c("merged_apes","mice","humans","bears","fin_whale")
# read in pop files 
allPopFiles=data.frame()
for(label in labels){
  popfiles=list.files(paste0(popfiledir,label),pattern=".txt")
  for(popfilename in popfiles){
    pop=unlist(lapply(strsplit(popfilename,"\\."),"[",1)) # get pop name from filename
    popfile=read.table(paste0(popfiledir,label,"/",popfilename),header=F)
    colnames(popfile) <- "sample_id"
    popfile$pop <- pop
    popfile$label <- label
    allPopFiles <- bind_rows(allPopFiles,popfile) 
  }
  
}
write.table(allPopFiles,paste0(popfiledir,"allPopFilesCombined.txt"),row.names=F,quote=F,sep="\t")
# want to assign ind1 and ind2 to pops ; note some won't be assigned -- keep those as "unassigned"

# process hamming distances
allAveragedHamming = data.frame()
for(label in labels){
  print(label)
  hamming = read.table(paste0(hammingindir,label,".totalAlleleCounts.HammingDistance.allIntervals.txt"),header=T)
  speciesPopFile = allPopFiles[allPopFiles$label==label,]
  hamming_merge_intermediate <- merge(hamming,speciesPopFile,by.x=c("label","ind1_alphabetical"),by.y=c("label","sample_id"))
  hamming_merge <- merge(hamming_merge_intermediate,speciesPopFile,by.x=c("label","ind2_alphabetical"),by.y=c("label","sample_id"),suffixes=c(".ind1",".ind2"))
  # create an alphabetical label of populations
  hamming_merge$populationComparisonLabel_alphabetical <- paste0(pmin(hamming_merge$pop.ind1,hamming_merge$pop.ind2),".",pmax(hamming_merge$pop.ind1,hamming_merge$pop.ind2))
  # and then want to average over each type
  write.table(hamming_merge,paste0(hammingindir,label,".totalAlleleCounts.HammingDistance.allIntervals.WITHPOPINFO.txt"))
  hamming_averaged <- hamming_merge %>% 
    group_by(label,populationComparisonLabel_alphabetical) %>%
    summarise(average_HammingDistance=mean(totalHammingDistance))
  write.table(hamming_merge,paste0(hammingindir,label,".totalAlleleCounts.HammingDistance.allIntervals.WITHPOPINFO.AVERAGEDPERCOMPARISONTYPE.txt"))
  allAveragedHamming <- bind_rows(allAveragedHamming,hamming_averaged)
}

# for the apes, change "merged_apes" to just "apes"
allAveragedHamming[allAveragedHamming$label=="merged_apes",]$label <- "apes"
allAveragedHamming$pop1 <- paste0(allAveragedHamming$label,"_",unlist(lapply(strsplit(allAveragedHamming$populationComparisonLabel_alphabetical,"\\."),"[",1)))
allAveragedHamming$pop2 <- paste0(allAveragedHamming$label,"_",unlist(lapply(strsplit(allAveragedHamming$populationComparisonLabel_alphabetical,"\\."),"[",2)))


write.table(allAveragedHamming,paste0(hammingindir,"ALLSPECIESCOMBINED.totalAlleleCounts.HammingDistance.allIntervals.WITHPOPINFO.AVERAGEDPERCOMPARISONTYPE.INCLUDESAPES.txt"),quote=F,row.names=F,sep="\t")

plot1 <- ggplot(allAveragedHamming,aes(x=pop1,y=pop2,fill=log10(average_HammingDistance)))+
  geom_tile()+
  theme(axis.text.x = element_text(angle=90))+
  #geom_text(aes(label=average_HammingDistance))+
  scale_fill_viridis()+
  ggtitle("raw hamming distance (not yet corrected by 2*callable genome size)")
plot1
ggsave(paste0(hammingindir,"hammingDistance.raw.notYetCorrectedByTargetSize.includesapes.png"),plot1,height=8,width = 11)


######## need to gather/read in targets ##########
# this also has info on spectra but don't use that. 
all7merTargets=data.frame()
# exclude vaquita because not needed for hamming 
for(sample in tableOfTargetsToInclude$sample){
  # some 'samples' may contain multipops like mice and bears
  sampleinfo = tableOfTargetsToInclude[tableOfTargetsToInclude$sample==sample,]
  targetsfilename = paste0(sampleinfo$targetsindir,"/",sampleinfo$targetsFile)
  
  targets=read.table(targetsfilename,header=sampleinfo$targetsHeader) # assigns whether has header or not
  colnames(targets) <- c("target_7mer","countOverAllIntervals")
  targets$sample <- sample
  all7merTargets = bind_rows(all7merTargets,targets)

}
# note: vaquita not needed for hamming
########## you are here -- something is odd with halbels
# change "human" sample label to "humans"
all7merTargets$label <- as.character(all7merTargets$sample)
#all7merTargets[all7merTargets$sample=="human",]$label <- "humans"
# change merged_apes to apes


# get total callable genome size
totalCallableGenome <- all7merTargets %>%
  group_by(sample,label) %>%
  summarise(totalCallableGenomeSize_sumOverAllTargets =sum(countOverAllIntervals))
totalCallableGenome

# apes have same callable sequence so can restrict just to one and call it merged_apes
# keep non-apes and add in one ape (Pan_paniscus )
totalCallableGenome_oneApe <- totalCallableGenome[!totalCallableGenome$label %in% c("Gorilla","Pan_troglodytes","Pongo_abelii","Pongo_pygmaeus"),]

# keeping in Pan_paniscus as rep (all have same amnt callabl genome tho)
totalCallableGenome_oneApe[totalCallableGenome_oneApe$sample=="Pan_paniscus",]$label <- "apes"

# keeping in Pan_paniscus as rep (all have same amnt callabl genome tho)
totalCallableGenome_oneApe[totalCallableGenome_oneApe$label=="merged_apes",]$sample <- "apes"

#### convert one ape to merged_apes because they are all the same

# okay so all the apes are the same and need them to be labelled as merged_apes (599074701)

write.table(totalCallableGenome_oneApe,paste0(hammingindir,"totalCallableGenomeSize.SummedOver7merTargets.weird_fin_whale.includesOneApeOnlyBCAllThesame.txt"),quote=F,sep="\t",row.names=F)
# fin whale remains weird -- very little of the genome 

# merge it:
hamming_with_callableGenome <- merge(allAveragedHamming,totalCallableGenome_oneApe,by="label")

hamming_with_callableGenome$totalCallableGenomeSize_sumOverAllTargets_x2 <- 2*hamming_with_callableGenome$totalCallableGenomeSize_sumOverAllTargets

hamming_with_callableGenome$average_HammingDistance_divBy2xGenomeSize <- hamming_with_callableGenome$average_HammingDistance/hamming_with_callableGenome$totalCallableGenomeSize_sumOverAllTargets_x2

head(hamming_with_callableGenome)
write.table(hamming_with_callableGenome,paste0(hammingindir,"ALLSPECIESCOMBINED.totalAlleleCounts.HammingDistance.allIntervals.WITHPOPINFO.AVERAGEDPERCOMPARISONTYPE.DIVIDEDBY2xGenomeSIZE.INCLUDESAPES.usethis.txt"),quote=F,sep="\t",row.names=F)

plot2 <- ggplot(hamming_with_callableGenome,aes(x=pop1,y=pop2,fill=log10(average_HammingDistance_divBy2xGenomeSize)))+
  geom_tile()+
  theme(axis.text.x = element_text(angle=90))+
  geom_text(aes(label=scientific(average_HammingDistance_divBy2xGenomeSize,digits = 1)))+
  scale_fill_viridis()+
  ggtitle("per-site hamming distance (divided by 2*callable genome size)")
plot2
ggsave(paste0(hammingindir,"hammingDistance.perSite.CorrectedBy2xGenomeSizee.log10.includesapes.png"),plot2,height=8,width = 15)


########### plot bear frac seg sites to compare to Hahn paper ########
# head(all1merSpectra)
# bear1merFrac <- all1merSpectra %>%
#   group_by(label) %>%
#   mutate(fracSegSites = total_mutations_projected/sum(total_mutations_projected)) %>%
#   filter(label %in% c("bears_ABC","bears_EUR","bears_PB","humans_AFR")) %>%
#   select(species,mutation_label,population,total_mutations_projected,fracSegSites )
# 
# bearplot1 <- ggplot(bear1merFrac,aes(x=mutation_label,y=fracSegSites, fill=label))+
#   geom_col(position="dodge")+
#   ggtitle("bear fraction of seg sites \nnot target corrected\nbased on projected number of mutations, to reduce sample size, but not yet downsampled to match vaquita snp count ")+
#   theme_bw()+
#   scale_fill_manual(values=c("sienna4","sienna2","dodgerblue","darkred"))+
#   scale_y_continuous(breaks=c(seq(0,0.4,0.05)))
# bearplot1
# ggsave("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/comparing_spectrum_and_phylo_divergence_WITHPROJECTIONS/20220222.plots.projectedCounts.epsilon.1.sepCpGs.no.addingApes_maskedMouseVaq.addingMultinomialDownsampling.MorePlotsForPAG.newTimeTree.MantelTest.hammingApes/BEARCOMPARISONPLOT.HAHN.png",bearplot1)
# 
# ####### make plot look more similar to hahn bear plot ########
# bear1merFrac$overallType <- bear1merFrac$mutation_label
# bear1merFrac[bear1merFrac$mutation_label %in% c("C.T","C.T_CpG"),]$overallType <- "C.T"
# bear1merFrac$mutation_label <- factor(bear1merFrac$mutation_label,levels=c("A.C","A.G","A.T","C.A","C.G","C.T_CpG","C.T"))
# bearplot2 <- ggplot(bear1merFrac[bear1merFrac$label=="bears_ABC",],aes(x=overallType,y=fracSegSites, fill=mutation_label))+
#   geom_col(position="stack")+
#   ggtitle("Bears_ABC only\nbear fraction of seg sites \nnot target corrected\nbased on projected number of mutations, to reduce sample size, but not yet downsampled to match vaquita snp count ")+
#   theme_bw()+
#   scale_fill_manual(values=c("grey","grey","grey","grey","grey","black","grey"))+
#   scale_y_continuous(breaks=c(seq(0,0.5,0.1)))
# bearplot2
# ggsave("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/comparing_spectrum_and_phylo_divergence_WITHPROJECTIONS/20220222.plots.projectedCounts.epsilon.1.sepCpGs.no.addingApes_maskedMouseVaq.addingMultinomialDownsampling.MorePlotsForPAG.newTimeTree.MantelTest.hammingApes/BEARCOMPARISONPLOT.HAHN.ONEBEARONLY.png",bearplot2)

############ TROUBLESHOOTING APES: ###########
# apes_missingness <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/hamming_distance/Troubleshooting_Apes/apes.plink.imiss",header=T)
# apes_missingness$species <- apes_missingness$IID
# ggplot(apes_missingness,aes(y=IID,x=F_MISS))+
#   geom_col()
# apes_missingness_popInfo <- merge(apes_missingness,allPopFiles[allPopFiles$label=="merged_apes",],by.x="IID",by.y="sample_id")
# # get average missingness per species:
# apes_missingness_popInfo_avg <- apes_missingness_popInfo %>%
#   group_by(pop) %>%
#   summarise(N_MISS_mean = mean(N_MISS),N_GENO_mean= mean(N_GENO),F_MISS_mean=mean(F_MISS))
# apes_missingness_popInfo_avg$pop <- paste0("apes_",apes_missingness_popInfo_avg$pop)
# ## maybe I can multiply by 1-FMiss to get genome size? 
# 
# # so what if I multiplied callable genome by 1-Fmiss_avg for the species (?) -- but which species? the one with the most missingness? minimum? 
# allAveragedHamming_apes <- allAveragedHamming[allAveragedHamming$label=="apes",]
# 
# allAveragedHamming_apes_m1 <- merge(hamming_with_callableGenome,apes_missingness_popInfo_avg,by.x="pop1",by.y="pop")
# allAveragedHamming_apes_m2 <- merge(allAveragedHamming_apes_m1,apes_missingness_popInfo_avg,by.x="pop2",by.y="pop",suffixes = c(".pop1",".pop2"))
# head(allAveragedHamming_apes_m2)
# 
# # okay so could try 1. add up missingness for pop 1 and 2; 2x min
# allAveragedHamming_apes_m2$pop1_genomeSizeTimes1minusFmiss <- allAveragedHamming_apes_m2$totalCallableGenomeSize_sumOverAllTargets*(1-allAveragedHamming_apes_m2$F_MISS_mean.pop1)
# 
# allAveragedHamming_apes_m2$pop2_genomeSizeTimes1minusFmiss <- allAveragedHamming_apes_m2$totalCallableGenomeSize_sumOverAllTargets*(1-allAveragedHamming_apes_m2$F_MISS_mean.pop2)
# 
# allAveragedHamming_apes_m2$pop1_plus_pop2_correctedGenomeSizes <- allAveragedHamming_apes_m2$pop1_genomeSizeTimes1minusFmiss+allAveragedHamming_apes_m2$pop2_genomeSizeTimes1minusFmiss
# 
# allAveragedHamming_apes_m2$average_HammingDistance_divBySumOfCorrectedGenomeSizes <- allAveragedHamming_apes_m2$average_HammingDistance / allAveragedHamming_apes_m2$pop1_plus_pop2_correctedGenomeSizes
# 
# ggplot(allAveragedHamming_apes_m2,aes(y=populationComparisonLabel_alphabetical,x=average_HammingDistance_divBySumOfCorrectedGenomeSizes))+
#   geom_col() # this looks a bit better 
# 
# # another option: take the minimum x2 of the two callable genome sizes:
# 
# allAveragedHamming_apes_m2$min_x2_of_pop1_plus_pop2_correctedGenomeSizes <- 2*pmin(allAveragedHamming_apes_m2$pop1_genomeSizeTimes1minusFmiss,allAveragedHamming_apes_m2$pop2_genomeSizeTimes1minusFmiss)
# 
# 
# allAveragedHamming_apes_m2$average_HammingDistance_divBy2xMinOfCorrectedGenomeSizes <- allAveragedHamming_apes_m2$average_HammingDistance / allAveragedHamming_apes_m2$min_x2_of_pop1_plus_pop2_correctedGenomeSizes
# 
# ggplot(allAveragedHamming_apes_m2,aes(y=populationComparisonLabel_alphabetical,x=average_HammingDistance_divBy2xMinOfCorrectedGenomeSizes))+
#   geom_col() # okay this looks worse 
# 
# ggplot(allAveragedHamming_apes_m2,aes(x=pop1,y=pop2,fill=average_HammingDistance_divBySumOfCorrectedGenomeSizes))+
#   geom_tile()
# 
# allAveragedHamming_apes_m2_melt <- melt(allAveragedHamming_apes_m2[,c("populationComparisonLabel_alphabetical","average_HammingDistance_divBy2xGenomeSize","average_HammingDistance_divBySumOfCorrectedGenomeSizes","average_HammingDistance_divBy2xMinOfCorrectedGenomeSizes")])
# 
# ggplot(allAveragedHamming_apes_m2_melt,aes(y=populationComparisonLabel_alphabetical,x=value,fill=variable))+
#   geom_col(position="dodge")
# 
# ############## add in phylo distances to compare ##########
# # process raxml tree here instead of in previous 
# raxmlTree=ape::read.tree("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/raxml_tree/RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick")
# 
# # want to get rid of the family/order info at end of tip. just keep species : 
# raxmlTree_renamedTips <- raxmlTree
# raxmlTree_renamedTips$tip.label <- paste0(unlist(lapply(strsplit(raxmlTree_renamedTips$tip.label,"_"),"[",1)),"_",unlist(lapply(strsplit(raxmlTree_renamedTips$tip.label,"_"),"[",2)))
# 
# # restrict to just the species to include
# raxmlTree_renamedTips_subset <- keep.tip(raxmlTree_renamedTips,speciesToInclude_in_raxmlTree)
# 
# raxml_cophenetic_dist <- cophenetic(raxmlTree_renamedTips_subset)
# # save this as an object
# saveRDS(raxml_cophenetic_dist,file=paste0(plotdir,"raxml_cophenetic_dist.dist"))
# 
# melt_cophenetic_func <- function(cophenetic_dist){
#   cophenetic_dist_melt <- melt(cophenetic_dist)
#   colnames(cophenetic_dist_melt) <- c("Sp1","Sp2","cophenetic_distance")
#   # get species names # ASSUMES ARE SEPARATED BY "_" eg mus_musculus
#   cophenetic_dist_melt$Sp1_species <- paste0(unlist(lapply(strsplit(as.character(cophenetic_dist_melt$Sp1),"_"),"[",1)),"_",unlist(lapply(strsplit(as.character(cophenetic_dist_melt$Sp1),"_"),"[",2)))
#   
#   cophenetic_dist_melt$Sp2_species <- paste0(unlist(lapply(strsplit(as.character(cophenetic_dist_melt$Sp2),"_"),"[",1)),"_",unlist(lapply(strsplit(as.character(cophenetic_dist_melt$Sp2),"_"),"[",2)))
#   
#   # need to deal with reciprocal duplicates
#   # put in alphabetical order so I can get rid of dups:
#   cophenetic_dist_melt$comparisonLabel <- paste0(pmin(cophenetic_dist_melt$Sp1_species,cophenetic_dist_melt$Sp2_species),".",pmax(cophenetic_dist_melt$Sp1_species,cophenetic_dist_melt$Sp2_species))  # # need to use pmin and pmax to get row by row min max. otherwise gives min max of whoel vector
#   # get rid of self to self comparisons:
#   cophenetic_dist_melt_distinct <- cophenetic_dist_melt[cophenetic_dist_melt$Sp1!=cophenetic_dist_melt$Sp2,]
#   # and get rid of reciprocal dups (sp A - B and sp B- A)
#   cophenetic_dist_melt_distinct <- cophenetic_dist_melt_distinct %>%
#     ungroup() %>% # adding ungroup as a safety measure on 20211201 not sure if ottally necessary though
#     distinct(comparisonLabel,.keep_all=T)  # restrict to just distinct comparison labels. 
#   dim(cophenetic_dist_melt_distinct)
#   return(cophenetic_dist_melt_distinct)
#   
#   
# }
# 
# phyloDistances = melt_cophenetic_func(raxml_cophenetic_dist)
# 
# ############## compare hamming and phylo #######
# apes=c("Pan_paniscus",    "Pan_troglodytes" ,"Pongo_abelii"  ,  "Gorilla_gorilla","Pongo_pygmaeus" )
# 
# phyloDistances_apes <- phyloDistances[phyloDistances$Sp1 %in%apes & phyloDistances$Sp2 %in% apes, ]
# phyloDistances_apes$Sp1 <- paste0("apes_",phyloDistances_apes$Sp1)
# phyloDistances_apes$Sp2 <- paste0("apes_",phyloDistances_apes$Sp2)
# phyloDistances_apes$comparisonLabel <-  paste0(pmin(phyloDistances_apes$Sp1,phyloDistances_apes$Sp2),".",pmax(phyloDistances_apes$Sp2,phyloDistances_apes$Sp2))
# 
# allAveragedHamming_apes_m2[allAveragedHamming_apes_m2$pop1=="apes_Gorilla",]$pop1 <- "apes_Gorilla_gorilla"
# allAveragedHamming_apes_m2[allAveragedHamming_apes_m2$pop2=="apes_Gorilla",]$pop2 <- "apes_Gorilla_gorilla"
# 
# allAveragedHamming_apes_m2$populationComparisonLabel_alphabetical <- paste0(pmin(allAveragedHamming_apes_m2$pop1,allAveragedHamming_apes_m2$pop2),".",pmax(allAveragedHamming_apes_m2$pop1,allAveragedHamming_apes_m2$pop2))
# 
# allAveragedHamming_apes_withPHYLO <- merge(allAveragedHamming_apes_m2,phyloDistances_apes,by.x="populationComparisonLabel_alphabetical",by.y="comparisonLabel")
# 
# 
# ggplot(allAveragedHamming_apes_withPHYLO,aes(x=cophenetic_distance,y=average_HammingDistance_divBy2xMinOfCorrectedGenomeSizes,color="2xMinCorrectedGenomes"))+
#   geom_point()+
#   geom_point(data=allAveragedHamming_apes_withPHYLO,aes(x=cophenetic_distance,y=average_HammingDistance_divBySumOfCorrectedGenomeSizes,color="sumOfFMissCorrectedGenomeSize"))+
#   geom_point(data=allAveragedHamming_apes_withPHYLO,aes(x=cophenetic_distance,y=average_HammingDistance_divBy2xGenomeSize,color="2xUncorrectedCallableGenomeSize"))+
#   geom_abline()+
#   geom_text(aes(label=populationComparisonLabel_alphabetical))
#   
# ###########OOH MISSING A LOT DUE TO MERGE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# #  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# #  !!!!!!!!!!!!!!!!!!!!! # when I merge, it's not necessarily that the species isn't called at that site. but that it doesn't have a snp. so may be missing or monormorphic. so lots and and lots of missing and will be correlated with 
# # so not necessarily missing data. aha. # ahhhha. hmmmm what to do about THAT ...
# ###### AHA AHA AHA 