########## troubleshooting targets ###############
require(reshape2)
require(dplyr)
require(ggplot2)
require(tidyr)
require(compositions)  # for clr
require(broom) # for tidy
require(ggrepel)
require(ggfortify) # for pca autoplot

# plotdir
plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/troubleshooting_target_sizes/"
########### note: dplyr is a bit dangerous when dividing by sum(X) -- if the df has ever been grouped by things, it maintains that grouping. This manifests as sneaky bugs for fraction of segregating sites.
# Be very cautious ! !!! need to go back and check some resutls 
separateCpGs="no" # yes or no

########## table of target locations ########
# restricting to just species with >1 individual:
#tableOfTargetsToInclude=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/spectra/tableOfSpectraToIncludeInComparisons.SpeciesWithMultipleIndividualsOnly.20211026.txt",header=T,sep="\t") # this also has info on spectra but don't use that. 
tableOfTargetsToInclude=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/spectra/tableOfSpectraToIncludeInComparisons.SpeciesWithMultipleIndividualsOnly.PlusApes.MaskedMouse.Vaquita.20211123.txt",header=T,sep="\t")

# okay do need to process the targets #
########## read in targets ##############
all7merTargets=data.frame()
all5merTargets=data.frame()
all3merTargets = data.frame()
all1merTargets = data.frame()
########## YOU ARE HERE -- deal with cpg targets ##########
for(sample in tableOfTargetsToInclude$sample){
  # some 'samples' may contain multipops like mice and bears
  sampleinfo = tableOfTargetsToInclude[tableOfTargetsToInclude$sample==sample,]
  targetsfilename = paste0(sampleinfo$targetsindir,"/",sampleinfo$targetsFile)
  
  targets=read.table(targetsfilename,header=sampleinfo$targetsHeader) # assigns whether has header or not
  colnames(targets) <- c("target_7mer","countOverAllIntervals")
  
  # get central 5mer, 3mer and 1mer:
  targets <- targets %>%
    mutate(target_5mer=substr(target_7mer,2,6),target_3mer=substr(target_7mer,3,5),target_1mer= substr(target_7mer,4,4)) # get substrings 
  
  ######### relabel target_1mer with CpGs if you are separating out CpGs ###########
  if(separateCpGs=="yes"){
    targets$CpGLabel <- ""
    targets[grepl("CG",targets$target_3mer),]$CpGLabel <- "_CpG"
    targets$target_1mer_CpGsNotLabeled <- targets$target_1mer # update with this 
    targets$target_1mer <- paste0(targets$target_1mer_CpGsNotLabeled,targets$CpGLabel)
    
  } # don't need an else statement, just keep going with targets target$target_1mer as is.
  
  
  
  # TARGETS per 3mer type #
  targets_7mer <- targets %>% group_by(target_7mer) %>% summarise(total_target_count=sum(countOverAllIntervals)) # original targets file, just doing to be consistent but it's identical to original targets
  targets_3mer <- targets %>% group_by(target_3mer) %>% summarise(total_target_count=sum(countOverAllIntervals))
  targets_5mer <- targets %>% group_by(target_5mer) %>% summarise(total_target_count=sum(countOverAllIntervals))
  targets_1mer <- targets %>% group_by(target_1mer) %>% summarise(total_target_count=sum(countOverAllIntervals))
  targets_7mer$species <- sample
  targets_5mer$species <- sample
  targets_3mer$species <- sample
  targets_1mer$species <- sample
  # make all species targets:
  
  
  # combine  all targets:
  all7merTargets = bind_rows(all7merTargets,targets_7mer)
  all5merTargets = bind_rows(all5merTargets,targets_5mer)
  all3merTargets = bind_rows(all3merTargets,targets_3mer)
  all1merTargets = bind_rows(all1merTargets,targets_1mer)
  
  
}

# plot all targets
p0 <- ggplot(all1merTargets,aes(y=species,x=total_target_count,fill=target_1mer))+
  geom_col(position="dodge")
p0
ggsave(paste0(plotdir,"plot0.comparing1merTargets.counts.png"),p0)

head(all7merTargets)
p1 <- ggplot(all7merTargets[all7merTargets$target_7mer %in% c("AAAAAAA","AAAAAAC","TTTAAAA"),],aes(x=species,y=total_target_count,fill=target_7mer))+
  geom_col(position="dodge")
ggsave(paste0(plotdir,"plot1.comparingAFew7merTargets.counts.png"),p1)

all7merTargets_proportions <- all7merTargets %>%
  group_by(species) %>% 
  mutate(total_target_proportion=total_target_count/sum(total_target_count))

p2 <- ggplot(all7merTargets_proportions[all7merTargets_proportions$target_7mer %in% c("AAAAAAA","AAAAAAC","TTTAAAA"),],aes(x=species,y=total_target_proportion,fill=target_7mer))+
  geom_col(position="dodge")
p2
ggsave(paste0(plotdir,"plot2.comparingAFew7merTargets.proportions.png"),p2)



all7merTargets_proportions_wide <- pivot_wider(all7merTargets_proportions,id_cols = target_7mer,names_from = species,values_from = total_target_proportion)

head(all7merTargets_proportions_wide)

p3 <- ggplot(all7merTargets_proportions_wide,aes(x=vaquita,y=humans))+
  geom_point()+
  geom_abline() +
  geom_label_repel(data=all7merTargets_proportions_wide[abs(all7merTargets_proportions_wide$vaquita-all7merTargets_proportions_wide$humans)>0.0003,],aes(label=target_7mer))
p3
ggsave(paste0(plotdir,"plot3.comparingVaquitaHumanTargets.png"),p3)

## look at those outliers across species: AAAAAAAA, TATATAT, ACACACA ,CACACAC

p4<- ggplot(all7merTargets_proportions[all7merTargets_proportions$target_7mer %in% c("AAAAAAA","TATATAT","CACACAC","ACACACA"),],aes(x=species,y=total_target_proportion,fill=target_7mer))+
  geom_col(position="dodge")
p4
ggsave(paste0(plotdir,"plot4.comparingOutlierTargets.proportions.png"),p4)

######### does this translate to weird mutation rates? ###########
allProjectedSpectra <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/allKSFSes_and_projected_Spectra/All_SpectrumCounts_ProjectedDownTo.16.Haploids.txt",header=T)
allProjectedSpectra$ancestral7mer <- substr(allProjectedSpectra$variable,1,7)

merge2 <- merge(allProjectedSpectra,all7merTargets_proportions,by.x = c("species","ancestral7mer"),by.y=c("species","target_7mer"))

merge2$mutRate <- merge2$totalSites_projectedDown_allFreqsSummed / merge2$total_target_count

merge2 <- merge2 %>%
  group_by(species,population,label) %>%
  mutate(rescaled_mut_rate=mutRate/sum(mutRate))

head(merge2)
sum(merge2[merge2$label=="bears_ABC",]$rescaled_mut_rate)
merge2$species <- factor(merge2$species,levels=c('vaquita','fin_whale','bears','humans','mice'))
merge2$label <- factor(merge2$label,levels=c('vaquita','fin_whale_ENP','fin_whale_GOC','bears_ABC','bears_PB','humans_AFR','humans_AMR','humans_EUR','humans_SAS','humans_EAS','mice_Mmc','mice_Mmm','mice_Mmd','mice_Ms'))
p5 <- ggplot(merge2[merge2$ancestral7mer %in% c("AAAAAAA","TATATAT","CACACAC","ACACACA","TTTAAAA"),] ,aes(y=variable,fill=species,x=rescaled_mut_rate,group=label) )+
  geom_col(position="dodge")
p5
ggsave(paste0(plotdir,"plot5.selectionOf7mer.normalizedMutationRates.AcrossSpecies.png"),height=6,width=10)
