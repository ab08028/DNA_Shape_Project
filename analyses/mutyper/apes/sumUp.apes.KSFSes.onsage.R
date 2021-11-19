require(reshape2)
require(dplyr)
############ NEED TO MODIFY FOR APES ###########
######### human info ###############
indir="/net/harris/vol1/home/beichman/apes/analyses/mutyper/mutyperResults_20211118/mutyper_ksfs_files/"
outdir=paste0(indir,"summedOverAllIntervals/")
dir.create(outdir)
intervals=paste0("chr",seq(1,22)) 
populationList=c('Gorilla', 'Pan_paniscus', 'Pongo_abelii', 'Pongo_pygmaeus', 'Pan_troglodytes')  # treating as pops even though are diff species - just a naming thing for convenience of merging with other datasets

speciesLabel="apes"
# adding population; very big so have to sum as you go.
ksfsFunction_apeSpecific <- function(indir,vectorOfIntervals,speciesLabel,pop){
  allKSFSes_summedUp = data.frame()
  # need to add as I go 
  #for(pop in ListOfPopulations){
  #print(paste0("starting ",pop))
  for(interval in vectorOfIntervals){
    print(interval)
    ksfs=read.table(paste0(indir,pop,".",interval,".mutyper.ksfs.nostrict.PASSONLY.BEDMASKED.txt"),header=T)

    ksfs_melt <- melt(ksfs,id.vars = c("sample_frequency"))
    colnames(ksfs_melt) <- c("sample_frequency","variable","totalSites")
    #ksfs_melt$interval <- as.character(interval)
    #ksfs_melt$species <- speciesLabel
    #ksfs_melt$population <- pop 
    allKSFSes_summedUp=bind_rows(allKSFSes_summedUp,ksfs_melt)
    # re sum up: (adding them up as you go)
    allKSFSes_summedUp <- allKSFSes_summedUp %>% group_by(variable,sample_frequency) %>% # don't need to group by pop because this is all pop specific already
      summarise(totalSites=sum(totalSites))
  }
  # allKSFSes_summed <- allKSFSes %>% group_by(species,population,variable,sample_frequency) %>%
  #  summarise(totalSites=sum(value))
  allKSFSes_summedUp$species <- speciesLabel
  allKSFSes_summedUp$population <- pop
  return(allKSFSes_summedUp)
}

for(pop in populationList){
  print(pop)
  ksfs_onepop = ksfsFunction_apeSpecific(indir,intervals,speciesLabel,pop)
  write.table(ksfs_onepop,paste0(outdir,"ksfs.7mer.",speciesLabel,".SUMMEDUP.",pop,".txt"),row.names = F,quote=F,sep="\t")

}