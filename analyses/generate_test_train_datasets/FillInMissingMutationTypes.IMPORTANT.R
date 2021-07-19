########### goal: process mouse data with 0s for missing mutation types ############
# and bigger windows (sum up some chrs)
# chr 1 and 2 are on their own 
# but the rest get summed up 

#### want to use tidyr complete to try to fill in observations 
spectrumdir="/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_spectrum_target_merged_PerChr/" # eventually save in resutls dir of dnashape project 
allData_multipop <- read.table(paste0(spectrumdir,"MULTIPOPULATION_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.txt"),header=T) # 

head(allData_multipop)
require(tidyr)
require(dplyr)
# want to try to add in empty mutation types 


# 24573 mutation types -- some must be missing per chr? 
#sum(length(unique(allData_multipop[allData_multipop$window==9,]$mutationType)))
# 24372 # okay so you are missing some. so why isn't complete working?? 

#setdiff(unique(allData_multipop$mutationType),unique(allData_multipop[allData_multipop$window==9,]$mutationType)) # note set diff is order dependent. 

######### this works > except doesn't fill in ancestral count ######

## start by filling in the missing targets:
# this data frame gives you the count of each ancestral 7mer per window/population
# and inserts 0s for missing 7mer targets (rather than just having them  not represented)
ancestral7mercounts_with0s <- allData_multipop %>%
  complete(ancestral7mer,nesting(window,population,label),fill=list(ancestral7merCount=0)) %>%
  select(window,population,label,ancestral7mer,ancestral7merCount) %>%
  unique()
dim(ancestral7mercounts_with0s) # 311258 

# okay now i want to fill in all combinations of missing mutations/targets
allData_multipop_fillInMissingRows <- allData_multipop %>%
  complete(nesting(mutationType,ancestral7mer),nesting(window,population,label),fill=list(mutationCount=0,mutationCount_divByTargetCount=0)) # don't want combinations of ancestral7mer and window that don't exist in the data so don't actually want it to be full length. or want to fill in targetcount = 0 if that is appropriate.
dim(df)

# now want to update ancestral counts: and take out old counts and old mutation rates
allData_multipop_fillInMissingRows_merged <- merge(select(allData_multipop_fillInMissingRows,-c(ancestral7merCount,mutationCount_divByTargetCount)),ancestral7mercounts_with0s,by=c("ancestral7mer","window","population","label"))

dim(allData_multipop_fillInMissingRows_merged)
head(allData_multipop_fillInMissingRows_merged)
#dim(ancestral7mercounts_with0s)
View(allData_multipop_fillInMissingRows_merged)

# okay this works! there are some rows with 0 targets and therefore 0 mutations and they are marked as such
# now get rid of old designations

# recalculate rate: (NOTE GENERATES SOME N/As) 
allData_multipop_fillInMissingRows_merged$mutationCount_divByTargetCount <- allData_multipop_fillInMissingRows_merged$mutationCount / allData_multipop_fillInMissingRows_merged$ancestral7merCount

sum(is.na(allData_multipop_fillInMissingRows_merged)) # 1080 NA's
sum(allData_multipop_fillInMissingRows_merged$ancestral7merCount==0) # 1080
# okay cool so this generates NAs for rows that have 0 and 0 

### write.table(allData_multipop_fillInMissingRows_merged,paste0(spectrumoutdir,"MULTIPOPULATION_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.INCLUDES0Entries.USETHIS.txt"),row.names=F,quote=F,sep="\t")
