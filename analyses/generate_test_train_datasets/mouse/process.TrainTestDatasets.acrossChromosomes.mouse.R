############## Process training/test datasets across chrosomes/scaffolds #########################################
require(ggplot2)
require(reshape2)
require(ggrepel)
require(dplyr)

############### add in mice next ##############
species="mouse"
populations=c("Mmd", "Mmc" ,"Mmm", "Ms")
chromosomecount=19
testSpectradf=data.frame()
trainSpectradf=data.frame()

outdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/training_test_datasets/mouse/mutyperResults_20210317_NOSTRICT_7mer/training_and_test_spectra_SUMMEDOVERCHROMOSOMES/"
dir.create(outdir,showWarnings = F,recursive = T)
indir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/training_test_datasets/mouse/mutyperResults_20210317_NOSTRICT_7mer/training_and_test_spectra_perChr/"
for(state in c("TESTING.even_bpOnly","TRAINING.odd_bpOnly")){
  allIntervalSpectra <- data.frame()
  
  for(i in seq(1,chromosomecount)){
    spectra=data.frame()
    chromosome=paste("chr",i,sep="")
    for(pop in populations){
      input=read.table(paste(indir,chromosome,"_",pop,"_",state,".mutyper.spectra.PERPOPULATION.ALLFREQS.NOSTRICT.txt",sep=""),header=T,stringsAsFactors = F)
      input$population <- pop
      spectra <- bind_rows(spectra, input) # in the only singletons category some are missing, so having NA put in using dplyr:bind_rows
      # ksfs:
    }
    spectra_melt <- melt(spectra)
    spectra_melt$interval <- chromosome
    spectra_melt$state <- state
    allIntervalSpectra <- bind_rows(allIntervalSpectra,data.frame(spectra_melt))
    
  }
  allIntervalSpectra[is.na(allIntervalSpectra)] <- 0
  
  
  allIntervalSpectra_totals_SummedUp <- data.frame(allIntervalSpectra %>% 
                                                           group_by(population,variable,state) %>%
                                                           summarise(totalSitesAllIntervals=sum(value)))
  # just grouping by state to kep it in as a variable and to be extra safe -- should only be train or test separated
  ### get fraction of all segregating sites 
  allIntervalSpectra_totals_SummedUp <- allIntervalSpectra_totals_SummedUp %>%
    group_by(population,state) %>%
    mutate(fractionOfAllSegSites=totalSitesAllIntervals/sum(totalSitesAllIntervals))  
  
  write.table(allIntervalSpectra_totals_SummedUp,paste0(outdir,species,".",state,".allPops.summedOverChr.spectra.NOSTRICT.txt"),row.names = F,quote=F,sep="\t")
  
}


