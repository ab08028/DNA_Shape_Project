########## Explore some prelim figures/analyses for the dna shape paper ##############
require(ggplot2)
require(reshape2)
require(dplyr)
require(GGally)
############# 7 mer targets : first want to compare genome targets ##################
# going to start with vaquita, mouse, bear, condor; these have targets in low complexity regions and in repeats (nothing masked from genome)
vaquita <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210225_7mer/mutyper_target_files/mutyper.targets.7mer.txt",header=F) # these are whole genome no mask
colnames(vaquita) <- c("target","countOverAllIntervals")
vaquita$species <- "vaquita"

mouse <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/summaryTables/ALLINTERVALS.mouse.mutyper.targets.7mer.nostrict.txt",header=T)
mouse$species <- "mouse"

brown_bear <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210331_7mer_NoStrict/mapped_to_brown_bear/summaryTables/ALLINTERVALS.brown_bear.mutyper.targets.7mer.nostrict.txt",header=T)
brown_bear$species <- "brown_bear"

condor <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/birds/mutyperResults_20210419/mutyper_targets_files/CA_condor.mutyper.targets.NOSTRICT.NOMASK.7mers.txt",header=F)
colnames(condor) <- c("target","countOverAllIntervals")
condor$species <- "condor"

########## adding other species targets #############


########################### combine ####################
all <- bind_rows(vaquita,mouse,brown_bear,condor)

# get proportions
all <- all %>%
  group_by(species) %>%
  mutate(proportion=countOverAllIntervals/sum(countOverAllIntervals))
# spread them:
head(all)
all_proportions_spread <- data.frame(spread(all[,c("target","species","proportion")],key=species,value=proportion))
head(all_proportions_spread)
########################## want to get regression between them ########################
species=c("vaquita","mouse","brown_bear","condor")
speciesAssignments=c(vaquita,mouse,bear,condor)
pairs=data.frame(combn(species,2))
dim(pairs)
numberOfPairs=dim(pairs)[2]

ggpairs(all_proportions_spread,columns=species,title="Genome Target Proportions",)
# can add diag="blank" if dont' want distribution in middle 
# Scatterplots of each pair of numeric variable are drawn on the left part of the figure. Pearson correlation is displayed on the right. Variable distribution is available on the diagonal.

ggcorr(all_proportions_spread[,species])
# to do ggplot stuff you can ggplot2::aes(colour=species)) add inside the call to ggpairs

# want to do this for 3mer targets as well. and see if correlations are similar. 
###############