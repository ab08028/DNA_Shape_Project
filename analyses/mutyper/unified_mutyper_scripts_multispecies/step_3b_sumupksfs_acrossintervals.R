# want to run for each species separately, but do pops within the script

args = commandArgs(trailingOnly=TRUE) # supply args
# first arg will be species
# second arg will be ? 

indir=as.character(args[1]) 
outdir=as.character(args[2]) 
species=as.character(args[3])
intervalCount=as.character(args[4])
prepend0=as.character(args[5])

######## going to get the pops from the dirs that are in the indir (if no dirs then it's no pops) #######

######## get the number of intervals ########
getIntervals <- function(intervalCount,prepend0) {
  
  if(intervalCount=="NA"){
    intervals="allautos"
  } else {
    
    if(prepend0=="FALSE") {
      intervals=as.character(c(seq(1,intervalCount))) # careful to make these characters when looping 
    } else if(prepend0=="TRUE") {
      intervals=as.character(sprintf("%02d", seq(1,intervalCount))) # if you need leading 0s 
    } else {
      stop("not a valid prepend0 option")
    }
  }
  return(intervals)
}

########## get the populations #######
getPopsFromIndir <- function(indir) {
  poplist=list.dirs(indir,full.names=F,recursive = F) 
  if(length(poplist)==0){
    print("no populations, using all inds")
    poplist=""
  } else {
    print("using pops:")
    print(poplist)
  }
  return(poplist)
}


######### sum up ksfs ########
sumupksfsacrossintervalsANDpopulations <- function(species,poplist,intervals,indir) {
  allKSFSes = data.frame()
  if(poplist!="") {
  for(pop in poplist){
    print(paste0("starting ",pop))
    for(interval in intervals){
      print(interval)
      popdir=paste0(indir,"/",pop,"/") # this will differ for other species annoying
      ksfs=read.table(paste0(popdir,species,".",pop,".int_or_chr_",interval,".mutyper.ksfs.SeeLogForFilters.maskALL.7mer.txt",header=T))
      ksfs_melt <- melt(ksfs,id.vars = c("sample_frequency"))
      ksfs_melt$interval <- as.character(interval)
      ksfs_melt$species <- species
      ksfs_melt$population <- pop 
      allKSFSes=bind_rows(allKSFSes,ksfs_melt)
    }
    #   # a check for uniqueness -- make sure the values of each sfs aren't the same (a sign that you accidentally read in the same sfs many times ) so am making sure that the freq=1 bins aren't all the same.
    if(length(unique(allKSFSes[allKSFSes$sample_frequency==1,]$value))==1){ # if this reduces down to one then you've added the same SFS together a bunch of times. 
      stop("are you reading the same SFS in??")
    }
  }
  } else {
    for(interval in intervals){
      print(interval)
      ksfs=read.table(paste0(indir,species,".int_or_chr_",interval,".mutyper.ksfs.SeeLogForFilters.maskALL.7mer.txt",header=T))
      ksfs_melt <- melt(ksfs,id.vars = c("sample_frequency"))
      ksfs_melt$interval <- as.character(interval)
      ksfs_melt$species <- species
      ksfs_melt$population <- species # just put in species for pop
      allKSFSes=bind_rows(allKSFSes,ksfs_melt)
    }
  }
  
  allKSFSes_summed <- allKSFSes %>% group_by(species,population,variable,sample_frequency) %>%
    summarise(totalSites=sum(value))
  return(allKSFSes_summed)
}

##### unified function ####
unifiedFunction_sumupksfs <- function(species,intervalCount,prepend0,indir){
  intervals=getIntervals(intervalCount,prepend0)
  poplist=getPopsFromIndir(indir)
  allKSFSes_summed = sumupksfsacrossintervalsANDpopulations(species,poplist,intervals,indir)
  return(allKSFSes_summed)
}

##### run it #######
allKSFSes_summed=unifiedFunction_sumupksfs(species,intervalCount,prepend0,indir)

##### write it out #######
write.table(allKSFSes_summed,paste0(outdir,species,".summedup.ksfs.allpops.txt"),row.names = F,quote=F,sep="\t")
