########### goal: process mouse data with 0s for missing mutation types ############
# and bigger windows (sum up some chrs)
# chr 1 and 2 are on their own 
# but the rest get summed up 
require(dplyr)
spectrumdir="/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/mutyper_spectrum_target_merged_PerChr/" # eventually save in resutls dir of dnashape project 

#datadf <- read.table(paste0(spectrumdir,"MULTIPOPULATION_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.INCLUDES0Entries.USETHIS.txt"),header=T,sep="\t")
# as of 20210719 changing to this one that has 0 entires dealt with *correctly*
datadf <- read.table(paste0(spectrumdir,"Mmd_Ms_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.FixedMissingMutationTypes.usethis.txt"),header=T,sep="\t")


dim(datadf) # 933774 good

# okay I want to sum things up across multiple chroms
# want 1+2, 2+4 to be on their own but then want 5+7 and 6+8 and 9+11+13 and 10+12 and 15+17+19 14+16+18
datadf$newGroup <- NA
datadf[datadf$window %in% c(1,3),]$newGroup <- 1                  
datadf[datadf$window %in% c(2,4),]$newGroup <- 2                  
datadf[datadf$window %in% c(5,7),]$newGroup <- 3                  
datadf[datadf$window %in% c(6,8),]$newGroup <- 4        
datadf[datadf$window %in% c(9,11,13),]$newGroup <- 5
datadf[datadf$window %in% c(10,12),]$newGroup <- 6
datadf[datadf$window %in% c(15,17,19),]$newGroup <- 7
datadf[datadf$window %in% c(14,16,18),]$newGroup <- 8

df2 <- datadf %>%
  group_by(newGroup,population,ancestral7mer, mutationType) %>%
  summarise(ancestral7merCount=sum(ancestral7merCount),mutationCount=sum(mutationCount)) ### actually think it may be ok maybe ? this isn't right! shouldn't sum up 7mer count per group -- that'll lead to a big over estimation I think ! uhoh! fix this!
# so issue I think is that it's somehow missing some 7mers?
# oh! maybe I actually did this ok???? seems to add up correctly actually.... let's chill out and find out later.
head(df2)
View(df2)


df2$mutationCount_divByTargetCount <- df2$mutationCount / df2$ancestral7merCount
head(df2)
# df2 %>%
#   group_by(newGroup) %>%
#   summarise(totalMutations=sum(mutationCount))

# need to label as train/test
df2$label <- ""
# odd is train even is test
df2[df2$newGroup %% 2 ==0,]$label <- "TEST"
df2[df2$newGroup %% 2 !=0,,]$label <- "TRAIN"

write.table(df2,paste0(spectrumdir,"TEMPORARY.NEWGROUPS.MULTIPOPULATION_spectrumCountsAndTargetCounts_perChromosome.allChrs.Labelled.INCLUDESCorrected0Entries.txt"),row.names=F,quote=F,sep="\t")
# okay this works. 
#df2[df2$mutationType=="AAAAAAC.AAAGAAC",]

#datadf[datadf$mutationType=="AAAAAAC.AAAGAAC" & datadf$population=="Mmd" & datadf$window %in% c(14,16,18),]

#sum(datadf[datadf$mutationType=="AAAAAAC.AAAGAAC" & datadf$population=="Mmd" & datadf$window %in% c(14,16,18),]$mutationCount)
