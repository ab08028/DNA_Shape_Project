require(reshape2)
require(dplyr)

######### human info ###############
human_indir="/net/harris/vol1/home/beichman/humans/analyses/mutyper/mutyperResults_20211101_STRICT_PASSONLY_BEDMASK/mutyper_ksfs_files/"
human_outdir="/net/harris/vol1/home/beichman/humans/analyses/mutyper/mutyperResults_20211101_STRICT_PASSONLY_BEDMASK/mutyper_ksfs_files/summedOverAllIntervals/"
dir.create(human_outdir)
human_intervals=paste0("chr",seq(1,22)) 
human_pops=c("AFR","AMR", "EAS", "EUR", "SAS") # addd more?
#human_pops=c("AFR")
human_speciesLabel="human"
# adding population; very big so have to sum as you go.
ksfsFunction_humanSpecific <- function(indir,vectorOfIntervals,ListOfPopulations,speciesLabel,pop){
  allKSFSes_summedUp = data.frame()
  # need to add as I go 
  #for(pop in ListOfPopulations){
  #print(paste0("starting ",pop))
  for(interval in vectorOfIntervals){
    print(interval)
    ksfs=read.table(paste0(indir,pop,"/",pop,".",interval,".mutyper.ksfs.StrictButThatsFineForHumans.PASSONLY.BEDMASKED.txt"),header=T)
    ksfs_melt <- melt(ksfs,id.vars = c("sample_frequency"))
    colnames(ksfs_melt) <- c("sample_frequency","variable","totalSites")
    #ksfs_melt$interval <- as.character(interval)
    #ksfs_melt$species <- speciesLabel
    #ksfs_melt$population <- pop 
    allKSFSes_summedUp=bind_rows(allKSFSes_summedUp,ksfs_melt)
    # re sum up:
    allKSFSes_summedUp <- allKSFSes_summedUp %>% group_by(variable,sample_frequency) %>% # don't need to group by pop because this is all pop specific already
      summarise(totalSites=sum(totalSites))
  }
  # allKSFSes_summed <- allKSFSes %>% group_by(species,population,variable,sample_frequency) %>%
  #  summarise(totalSites=sum(value))
  allKSFSes_summed$species <- speciesLabel
  allKSFSes_summed$population <- pop
  return(allKSFSes_summed)
}

for(pop in human_pops){
  print(pop)
  human_ksfs_onePop = ksfsFunction_humanSpecific(human_indir,human_intervals,human_pops,human_speciesLabel,pop)
  write.table(human_ksfs,paste0(human_outdir,"ksfs.7mer.humans.SUMMEDUP.",pop,".txt"),row.names = F,quote=F,sep="\t")

}