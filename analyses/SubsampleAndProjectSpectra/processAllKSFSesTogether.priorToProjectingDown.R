############ Process all ksfses and project ###########
require(reshape2)
require(dplyr)
# for some species need to sum up over intervals
# and some just are fine by themselves
# 
outdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/allKSFSes_and_projected_Spectra/"



######### human info ###############
# NOTE : humans occur in their own script because it needs to be run on sage (too big/slow for laptop)


######## bear info ##########
bear_indir="/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210331_7mer_NoStrict/mapped_to_brown_bear/mutyper_ksfs_files/" # just brown bear ref
bear_intervals=seq(1,31) 
bear_pops=c("EUR","ABC","PB")
bear_speciesLabel="bears"

ksfsFunction_bearSpecific <- function(indir,vectorOfIntervals,ListOfPopulations,speciesLabel){
  allKSFSes = data.frame()
  for(pop in ListOfPopulations){
    print(paste0("starting ",pop))
    for(interval in vectorOfIntervals){
      print(interval)
      intpopdir=paste0(indir,"interval_",interval,"/",pop,"/") # this will differ for other species annoying
      ksfs=read.table(paste0(intpopdir,pop,"_samples_interval_",interval,".ref_brown_bear.mutyper.7mer.ksfs.nostrict.txt"),header=T) # just brown bear ref.
      ksfs_melt <- melt(ksfs,id.vars = c("sample_frequency"))
      ksfs_melt$interval <- as.character(interval)
      ksfs_melt$species <- speciesLabel
      ksfs_melt$population <- pop 
      allKSFSes=bind_rows(allKSFSes,ksfs_melt)
    }
    #   # a check for uniqueness -- make sure the values of each sfs aren't the same (a sign that you accidentally read in the same sfs many times ) so am making sure that the freq=1 bins aren't all the same.
    if(length(unique(allKSFSes[allKSFSes$sample_frequency==1,]$value))==1){ # if this reduces down to one then you've added the same SFS together a bunch of times. 
      print("are you reading the same SFS in??")
      break 
    }
  }
  
  
  allKSFSes_summed <- allKSFSes %>% group_by(species,population,variable,sample_frequency) %>%
    summarise(totalSites=sum(value))
  return(allKSFSes_summed)
}

bear_ksfs = ksfsFunction_bearSpecific(bear_indir,bear_intervals,bear_pops,bear_speciesLabel)

write.table(bear_ksfs,paste0(outdir,"ksfs.7mer.bears.txt"),row.names = F,quote=F,sep="\t")
# write another copy out in original dir
write.table(bear_ksfs,paste0(bear_indir,"ksfs.7mer.bears.txt"),row.names = F,quote=F,sep="\t")


########## mice info #############
# as of 20211123 am using new rep masked ksfs:

#mice_indir="/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_ksfs_files/"
mice_indir="/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20211123_NOSTRICT_7mer_REPEATMASKED/mutyper_ksfs_files/"
mice_intervals=paste0("chr",seq(1,19))
mice_pops=c("Mmc","Mmd","Mmm","Ms")
mice_speciesLabel="mice"

ksfsFunction_miceSpecific <- function(indir,vectorOfIntervals,ListOfPopulations,speciesLabel){
  allKSFSes = data.frame()
  for(pop in ListOfPopulations){
    print(paste0("starting ",pop))
    for(interval in vectorOfIntervals){
      print(interval)
      ksfs=read.table(paste0(indir,pop,"/",pop,"_samples_interval_",interval,".mutyper.ksfs.7mer.nostrict.REPEATMASKED.txt"),header=T)
      ksfs_melt <- melt(ksfs,id.vars = c("sample_frequency"))
      ksfs_melt$interval <- as.character(interval)
      ksfs_melt$species <- speciesLabel
      ksfs_melt$population <- pop 
      allKSFSes=bind_rows(allKSFSes,ksfs_melt)
    }
    ## a check for uniqueness -- make sure the values of each sfs aren't the same (a sign that you accidentally read in the same sfs many times ) so am making sure that the freq=1 bins aren't all the same.
    if(length(unique(allKSFSes[allKSFSes$sample_frequency==1,]$value))==1){ # if this reduces down to one then you've added the same SFS together a bunch of times. 
      print("are you reading the same SFS in??")
      break 
    }
  }
  allKSFSes_summed <- allKSFSes %>% group_by(species,population,variable,sample_frequency) %>%
    summarise(totalSites=sum(value))
  return(allKSFSes_summed)
}

mice_ksfs = ksfsFunction_miceSpecific(mice_indir,mice_intervals,mice_pops,mice_speciesLabel)

write.table(mice_ksfs,paste0(outdir,"ksfs.7mer.mice.RepeatMasked.AllIntervals.summed.txt"),row.names = F,quote=F,sep="\t")
# write another copy out in original dir
write.table(mice_ksfs,paste0(mice_indir,"ksfs.7mer.mice.RepeatMasked.AllIntervals.summed.txt"),row.names = F,quote=F,sep="\t")


############ fin whale info ############
finwhale_indir="/Users/annabelbeichman/Documents/UW/WhaleProject/fin_whales/results/mutyper/downloaddate_20211015_whole_genomes_plusOutgroups/mutyperResults_20211029_NOSTRICT_CALLABLESITESONLY_POLARIZED_7mer/mutyper_ksfs_files/" # this is improved version with more sites and 2 low cov inds removed 
finwhale_intervals=sprintf("%02d", seq(1,96))
finwhale_pops=c("ENP","GOC")
finwhale_speciesLabel="fin_whale"

################### fin whale sum up #####################
ksfsFunction_finWhaleSpecific <- function(indir,vectorOfIntervals,ListOfPopulations,speciesLabel){
  allKSFSes = data.frame()
  for(pop in ListOfPopulations){
    print(paste0("starting ",pop))
    for(interval in vectorOfIntervals){
      print(interval)
      intpopdir=paste0(indir,"interval_",interval,"/",pop,"/") # this will differ for other species annoying
      ksfs=read.table(paste0(intpopdir,pop,"_samples_interval_",interval,".mutyper.7mer.ksfs.txt"),header=T)
      ksfs_melt <- melt(ksfs,id.vars = c("sample_frequency"))
      ksfs_melt$interval <- as.character(interval)
      ksfs_melt$species <- speciesLabel
      ksfs_melt$population <- pop 
      allKSFSes=bind_rows(allKSFSes,ksfs_melt)
    }
    #   # a check for uniqueness -- make sure the values of each sfs aren't the same (a sign that you accidentally read in the same sfs many times ) so am making sure that the freq=1 bins aren't all the same.
    if(length(unique(allKSFSes[allKSFSes$sample_frequency==1,]$value))==1){ # if this reduces down to one then you've added the same SFS together a bunch of times. 
      print("are you reading the same SFS in??")
      break 
    }
  }
  
  
  allKSFSes_summed <- allKSFSes %>% group_by(species,population,variable,sample_frequency) %>%
    summarise(totalSites=sum(value))
  return(allKSFSes_summed)
}

finwhale_ksfs = ksfsFunction_finWhaleSpecific(finwhale_indir,finwhale_intervals,finwhale_pops,finwhale_speciesLabel)

write.table(finwhale_ksfs,paste0(outdir,"ksfs.7mer.finwhale.WholeGenome.PASSOnly.allIntervals.txt"),row.names = F,quote=F,sep="\t")
# write another copy out in original dir
write.table(finwhale_ksfs,paste0(finwhale_indir,"ksfs.7mer.finwhale.WholeGenome.PASSOnly.allIntervals.txt"),row.names = F,quote=F,sep="\t")

########### vaquita info #################
vaquita_indir="/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210225_7mer/mutyper_ksfs_files/"
vaquita_infile="vaquita.mutyper.ksfs.7mer.txt"
#vaquita_intervals=NA
#vaquita_pops=NA
vaquita_speciesLabel="vaquita"

vaquita_unmelted = read.table(paste0(vaquita_indir,vaquita_infile),header=T)
vaquita_ksfs <- melt(vaquita_unmelted,id.vars =c("sample_frequency"),value.name = "totalSites")
vaquita_ksfs$species <- vaquita_speciesLabel
vaquita_ksfs$population <- vaquita_speciesLabel
head(vaquita_ksfs)

write.table(vaquita_ksfs,paste0(outdir,"ksfs.7mer.vaquita.txt"),row.names = F,quote=F,sep="\t")
# write another copy out in original dir
write.table(vaquita_ksfs,paste0(vaquita_indir,"ksfs.7mer.vaquita.melted.txt"),row.names = F,quote=F,sep="\t")


############ write a note on where the ksfses came from ###########
sink(paste0(outdir,"mutyperRunsUsedToMakeTheseKSFes.txt"))
print(Sys.Date())
print(bear_indir)
print(mice_indir)
print(vaquita_indir)
print(finwhale_indir)
sink()
