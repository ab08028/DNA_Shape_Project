############# want to get all KSFSes and subsample down to some min sample size (haploids) #############
wd="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/allKSFSes_and_projected_Spectra/"
bears=read.table(paste0(wd,"ksfs.7mer.bears.txt"),header=T,sep="\t")
vaquita=read.table(paste0(wd,"ksfs.7mer.vaquita.WholeGenome.PASSOnly.allIntervals.txt"),header=T,sep="\t")
mice=read.table(paste0(wd,"ksfs.7mer.mice.NotRepeatMasked.AllIntervals.summed.txt"),header=T,sep="\t")
#humans
finwhale=read.table(paste0(wd,"ksfs.7mer.finwhale.WholeGenome.PASSOnly.allIntervals.txt"),header=T,sep="\t")
# humans=""

# combine:
allSFSes <- bind_rows(bears,vaquita,mice,finwhale) # add humans

allSFSes$label <- paste0(allSFSes$species,"_",allSFSes$population)
# reset vaquita so it's not vaquita_vaquita
allSFSes[allSFSes$species=="vaquita",]$label <- "vaquita"

####### restrict to: 
Species_Pops_ToProject =c("fin_whale_ENP","mice_Mmd","mice_Ms","human_AFR","vaquita","bears_PB","bears_ABC")

speciesIWant_SFSes <-
  allSFSes %>%
  filter(label %in% Species_Pops_ToProject )
# note sample sizes are in haploids 
##### all SFSes:

# add a column for each species/pop that has its max freq (will be the same for all rows within a species)
speciesIWant_SFSes <- speciesIWant_SFSes %>%
  group_by(species,population,label) %>%
  mutate(n_maxFreq_plus1_haploids=max(sample_frequency)+1) 

tail(speciesIWant_SFSes)

# print a summary as well
speciesIWant_SFSes %>%
  group_by(species,population,label) %>%
  summarise(n_maxFreq_plus1_haploids=max(sample_frequency)+1) 
# so Ms and ABC are limited at n = 16
# get projection size
k=min(speciesIWant_SFSes$n_maxFreq_plus1_haploids)
print(paste0("projection size is ",k," haploids"))
# choose pops to project (don't want to do all - too restrictive)
######################## get min sample size for projection #############
# GET THESE LABELS Species_Pops_ToProject =c("fin_whale_ENP","mouse_Mmc","mouse_Ms","human_AFR","vaquita","polar_bear","brown_bear_ABC")
# note sample sizes are in haploids 
head(speciesIWant_SFSes)

# subset to just the pops I want
######## then let's do the projection #######
 # 1. use dplyr to add full sample size (n ) in haploids (including fixed) to each population (so max(freq)+1)
  # 2. use equation to get prob of observing mutation  at least once (so don't need binomial sampling because not getting ksfs) : for each mutation type/species SUM_OVER_ALL_i[(count_mutationType_Freq_i) * (1 - (i/n)^k - ((n-i)/n)^k) ] <- this multiplies count # by probability of if you drew k haploids out of the n haploids, you don't want to get a fixed site (i/n ^k) and you want to draw your site at least once (n-i/n)^k -- but you don't care how many times it appears (like --population in mutyper spectrum. just gets counted once. so don't need binomial coefficients etc.). Then sum up over all frequencies. 

# 3 . sum up counts for each mutation type across all frequencies
# that's your new adjusted count! 
# you set k above 

# for checking
#intermediateFile <- speciesIWant_SFSes %>%
#  mutate(totalSites_projected_intermediate = totalSites*(1 - (sample_frequency/n_maxFreq_plus1_haploids)^k - #((n_maxFreq_plus1_haploids - sample_frequency)/n_maxFreq_plus1_haploids)^k))


projectedSFSes <- speciesIWant_SFSes %>%
  mutate(totalSites_projected_intermediate = totalSites*(1 - (sample_frequency/n_maxFreq_plus1_haploids)^k - ((n_maxFreq_plus1_haploids - sample_frequency)/n_maxFreq_plus1_haploids)^k)) %>% # did some spot checking by hand, this works 
  group_by(species,population,label,variable) %>%
  summarise(totalSites_projectedDown_allFreqsSummed=sum(totalSites_projected_intermediate))

  # this makes sense because at intermediate frequencies you see mostly same count, but at high or low freqs
  # you either miss it (low freqs) or it becomes fixed (high freqs) so you get slightly lower counts
  # cool.
  
# some checks
#intermediateFile[intermediateFile$species=="vaquita" & intermediateFile$variable=="AAAAAAA.AAACAAA",]
#sum(intermediateFile[intermediateFile$species=="vaquita" & intermediateFile$variable=="AAAAAAA.AAACAAA",]$totalSites)
#sum(intermediateFile[intermediateFile$species=="vaquita" & intermediateFile$variable=="AAAAAAA.AAACAAA",]$totalSites_projected_intermediate) # should observe fewer sites
#projectedSFSes[projectedSFSes$species=="vaquita" & projectedSFSes$variable=="AAAAAAA.AAACAAA",]
head(projectedSFSes)

write.table(projectedSFSes,paste0(wd,"All_SpectrumCounts_ProjectedDownTo.",k,".Haploids.txt"),col.names =T,row.names = F,quote=F)
