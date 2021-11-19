############# sum up mouse 7mer targets over chromosomes ##########
# now using data without strict 
require(reshape2)
require(dplyr)
wd="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/mutyper/apes/mutyperResults_20211118/mutyper_targets_files/" ## 
tabledir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/mutyper/apes/mutyperResults_20211118/summaryTables/"
dir.create(tabledir)
intervalCount=22
speciesList=c('Gorilla', 'Pan_paniscus', 'Pongo_abelii', 'Pongo_pygmaeus', 'Pan_troglodytes')

for(species in speciesList){
  allTargets=data.frame() # must clear between species 
for(i in seq(1,intervalCount)){
  input=read.table(paste(wd,species,".mutyper.targets.chr.",i,".NOSTRICT.BedMasked.7mers.txt",sep=""),header=F)
  colnames(input) <- c("target","count")
  input$chr <- paste0("chr",i) # this fixes bug from bear script
  input_melt <- melt(input,id.vars = c("target","chr"))
  input_melt$species <- species # safety net
  allTargets <- bind_rows(allTargets,data.frame(input_melt))
  
}
# sum up over chrs:
allTargetsSummedOverIntervals <- allTargets %>%
  group_by(target,species) %>% # added insurance; should just be one species though
  summarise(countOverAllIntervals=sum(value)) %>%
  select(target,countOverAllIntervals)

write.table(data.frame(allTargetsSummedOverIntervals),paste(tabledir,"/",species,".ALLINTERVALS.mutyper.targets.7mer.nostrict.bedmasked.txt",sep=""),row.names = F,quote=F,sep="\t")
}
