####### want total targets and want population information #########

require(viridis)
require(ggplot2)
require(dplyr)
require(reshape2)
require(scales)
# human pop files #
popfiledir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/populationFiles/"
# note not doing hamming for apes or vaquita
hammingindir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/hamming_distance/"
# update to make sure this is correct !! targets must be masked in same way as hamming distances
# so if rep masking occurred for hamming dist, rep masking must occur for targets (can stay consisten by using targets and variants files from same mutyper runs)
### note th
tableOfTargetsToInclude=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/spectra/tableOfSpectraToIncludeInComparisons.SpeciesWithMultipleIndividualsOnly.PlusApes.MaskedMouse.Vaquita.20211123.txt",header=T,sep="\t") 
labels=c("mice","humans","bears","fin_whale")
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

write.table(allAveragedHamming,paste0(hammingindir,"ALLSPECIESCOMBINED.totalAlleleCounts.HammingDistance.allIntervals.WITHPOPINFO.AVERAGEDPERCOMPARISONTYPE.txt"))

allAveragedHamming$pop1 <- paste0(allAveragedHamming$label,"_",unlist(lapply(strsplit(allAveragedHamming$populationComparisonLabel_alphabetical,"\\."),"[",1)))
allAveragedHamming$pop2 <- paste0(allAveragedHamming$label,"_",unlist(lapply(strsplit(allAveragedHamming$populationComparisonLabel_alphabetical,"\\."),"[",2)))

plot1 <- ggplot(allAveragedHamming,aes(x=pop1,y=pop2,fill=log10(average_HammingDistance)))+
  geom_tile()+
  theme(axis.text.x = element_text(angle=90))+
  #geom_text(aes(label=average_HammingDistance))+
  scale_fill_viridis()+
  ggtitle("raw hamming distance (not yet corrected by 2*callable genome size)")
plot1
ggsave(paste0(hammingindir,"hammingDistance.raw.notYetCorrectedByTargetSize.png"),plot1,height=8,width = 11)


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
# get total callable genome size
totalCallableGenome <- all7merTargets %>%
  group_by(sample,label) %>%
  summarise(totalCallableGenomeSize_sumOverAllTargets =sum(countOverAllIntervals))
totalCallableGenome
write.table(totalCallableGenome,paste0(hammingindir,"totalCallableGenomeSize.SummedOver7merTargets.weird_fin_whale.txt"),quote=F,sep="\t",row.names=F)
# fin whale remains weird -- very little of the genome 

# merge it:
hamming_with_callableGenome <- merge(allAveragedHamming,totalCallableGenome,by="label")

hamming_with_callableGenome$totalCallableGenomeSize_sumOverAllTargets_x2 <- 2*hamming_with_callableGenome$totalCallableGenomeSize_sumOverAllTargets

hamming_with_callableGenome$average_HammingDistance_divBy2xGenomeSize <- hamming_with_callableGenome$average_HammingDistance/hamming_with_callableGenome$totalCallableGenomeSize_sumOverAllTargets_x2

head(hamming_with_callableGenome)
write.table(hamming_with_callableGenome,paste0(hammingindir,"ALLSPECIESCOMBINED.totalAlleleCounts.HammingDistance.allIntervals.WITHPOPINFO.AVERAGEDPERCOMPARISONTYPE.DIVIDEDBY2xGenomeSIZE.usethis.txt"),quote=F,sep="\t",row.names=F)

plot2 <- ggplot(hamming_with_callableGenome,aes(x=pop1,y=pop2,fill=log10(average_HammingDistance_divBy2xGenomeSize)))+
  geom_tile()+
  theme(axis.text.x = element_text(angle=90))+
  geom_text(aes(label=scientific(average_HammingDistance_divBy2xGenomeSize,digits = 1)))+
  scale_fill_viridis()+
  ggtitle("per-site hamming distance (divided by 2*callable genome size)")
plot2
ggsave(paste0(hammingindir,"hammingDistance.perSite.CorrectedBy2xGenomeSizee.log10.png"),plot2,height=8,width = 15)

