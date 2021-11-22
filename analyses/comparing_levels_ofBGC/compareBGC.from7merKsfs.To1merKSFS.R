#################### want to plot BGC from 7mer ksfs for all species ########
# note: not projecting so sample sizes will be very different
require(ggrepel)
require(dplyr)
require(ggplot2)
require(reshape2)
require(viridis)
wd="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/allKSFSes_and_projected_Spectra/"
outdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/comparing_levels_ofBGC/"
bears=read.table(paste0(wd,"ksfs.7mer.bears.txt"),header=T,sep="\t")
vaquita=read.table(paste0(wd,"ksfs.7mer.vaquita.WholeGenome.PASSOnly.allIntervals.txt"),header=T,sep="\t")
mice=read.table(paste0(wd,"ksfs.7mer.mice.NotRepeatMasked.AllIntervals.summed.txt"),header=T,sep="\t")
finwhale=read.table(paste0(wd,"ksfs.7mer.finwhale.WholeGenome.PASSOnly.allIntervals.txt"),header=T,sep="\t")
# humans="" # getting from sage

# these are very big: (gzip them)

humans_AFR=read.table(paste0(wd,"ksfs.7mer.humans.SUMMEDUP.AFR.txt.gz"),header=T,sep="\t")
#humans_EUR=read.table(paste0(wd,"ksfs.7mer.humans.SUMMEDUP.EUR.txt.gz"),header=T,sep="\t")
#humans_SAS=read.table(paste0(wd,"ksfs.7mer.humans.SUMMEDUP.SAS.txt.gz"),header=T,sep="\t")
#humans_EAS=read.table(paste0(wd,"ksfs.7mer.humans.SUMMEDUP.EAS.txt.gz"),header=T,sep="\t")
#humans_AMR=read.table(paste0(wd,"ksfs.7mer.humans.SUMMEDUP.AMR.txt.gz"),header=T,sep="\t")

# combine:
allSFSes <- bind_rows(bears,vaquita,mice,finwhale,humans_AFR) #,humans_AMR,humans_EUR,humans_SAS,humans_EAS) # add humans
# add s to end of humans
allSFSes[allSFSes$species=="human",]$species <- "humans"
# make a label that is species_pop so it'd be bears_EUR and humans_EUR
allSFSes$label <- paste0(allSFSes$species,"_",allSFSes$population)
# reset vaquita so it's not vaquita_vaquita
allSFSes[allSFSes$species=="vaquita",]$label <- "vaquita"

######### want to get 1mer spectra ########
allSFSes$mutation_type_1mer <- paste0(substr(allSFSes$variable,4,4),".",substr(allSFSes$variable,12,12))

# sum up over 1mer spectrum
allSFSes_1mer <- allSFSes %>%
  group_by(label,species,population,sample_frequency,mutation_type_1mer) %>%
  summarise(total_mutations_1mer = sum(totalSites))

# then want to get allele freqs
allSFSes_1mer <- allSFSes_1mer %>%
  group_by(species,label,population) %>%
  mutate(alleleFrequency=sample_frequency/(max(sample_frequency)+1)) %>% # get alelle freq
  ungroup() %>%
  group_by(species,label,population,sample_frequency) %>% # get fraction of mutations at a given frequency
  mutate(mutationFractionPerFrequency=total_mutations_1mer/sum(total_mutations_1mer)) # divide by max sfs +1 
# note that this is mutation Fraction per frequency not overall

bgcPlot1 <- ggplot(allSFSes_1mer,aes(x=alleleFrequency,y=mutationFractionPerFrequency,color=label))+
  geom_line()+
  facet_wrap(~mutation_type_1mer,scales="free_y")+
  theme_bw()
bgcPlot1
ggsave(paste0(outdir,"1merFractions.PerMutationType.ReducedFrom7merKSFS.SomethingWrongWithFinWhale.png"),bgcPlot1,height=7,width=12)
######### want to categorize as SW and WS as well ########
allSFSes$BGC_category <- ""
allSFSes[allSFSes$mutation_type_1mer %in% c("A.T","C.G"),]$BGC_category <- "BGC-conservative"
allSFSes[allSFSes$mutation_type_1mer %in% c("C.A","C.T"),]$BGC_category <- "Strong to Weak"
allSFSes[allSFSes$mutation_type_1mer %in% c("A.C","A.G"),]$BGC_category <- "Weak to Strong"

allSFSes_BGC <- allSFSes %>%
  group_by(label,species,population,sample_frequency,BGC_category) %>%
  summarise(total_mutations_BGC = sum(totalSites))


# then want to get allele freqs
allSFSes_BGC <- allSFSes_BGC %>%
  group_by(species,label,population) %>%
  mutate(alleleFrequency=sample_frequency/(max(sample_frequency)+1)) %>% # get alelle freq
  ungroup() %>%
  group_by(species,label,population,sample_frequency) %>% # get fraction of mutations at a given frequency
  mutate(mutationFractionPerFrequency=total_mutations_BGC/sum(total_mutations_BGC)) # divide by max sfs +1 
# note that this is mutation Fraction per frequency not overall

bgcPlot2 <- ggplot(allSFSes_BGC,aes(x=alleleFrequency,y=mutationFractionPerFrequency,color=label))+
  geom_line()+
  facet_wrap(~BGC_category,scales="free_y")+
  theme_bw()
bgcPlot2
ggsave(paste0(outdir,"1merFractions.PerBGCCategory.ReducedFrom7merKSFS.SomethingWrongWithFinWhale.png"),bgcPlot2,height=7,width=12)

######### make nicer plots #############
# exclude fin whale ENP and GOC too 
# add labels 
bgcPlot1b <- ggplot(allSFSes_1mer[allSFSes_1mer$species!="fin_whale",],aes(x=alleleFrequency,y=mutationFractionPerFrequency,color=label,label=label))+
  geom_line()+
  facet_wrap(~mutation_type_1mer,scales="free_y")+
  theme_bw()# +
  #geom_label_repel() 
bgcPlot1b

ggsave(paste0(outdir,"1merFractions.PerMutationType.ReducedFrom7merKSFS.ExcludeFinWhale.png"),bgcPlot1b,height=7,width=12)

bgcPlot2b <- ggplot(allSFSes_BGC[allSFSes_BGC$species!="fin_whale",],aes(x=alleleFrequency,y=mutationFractionPerFrequency,color=label))+
  geom_line()+
  facet_wrap(~BGC_category,scales="free_y")+
  theme_bw()
bgcPlot2b
ggsave(paste0(outdir,"1merFractions.PerBGCCategory.ReducedFrom7merKSFS.ExcludeFinWhaleENP.png"),bgcPlot2b,height=7,width=12)

############## make even nicer : for presentations ###########


bgcPlot2c <- ggplot(allSFSes_BGC[allSFSes_BGC$label %in% c("humans_AFR","bears_ABC","bears_PB","mice_Mmd","vaquita"),],aes(x=alleleFrequency,y=mutationFractionPerFrequency,color=label))+
  geom_line(size=1.5)+
  facet_wrap(~BGC_category,scales="free_y")+
  theme_bw()+
  scale_color_manual(values= c("#D55E00", "#0072B2","#999999",  "#009E73", "#CC79A7")) # this is a CB friendly pal from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
bgcPlot2c
ggsave(paste0(outdir,"1merFractions.PerBGCCategory.ReducedFrom7merKSFS.SubSetOfSpecies.png"),bgcPlot2c,height=4,width=12)

# add in fin whale GOC for bmha talk:
bgcPlot2d <- ggplot(allSFSes_BGC[allSFSes_BGC$label %in% c("humans_AFR","bears_ABC","bears_PB","mice_Mmd","vaquita","fin_whale_GOC"),],aes(x=alleleFrequency,y=mutationFractionPerFrequency,color=label))+
  geom_line(size=1.5)+
  facet_wrap(~BGC_category,scales="free_y")+
  theme_bw()+
  scale_color_manual(values= c("#D55E00", "#0072B2","#F0E442","#999999",  "#009E73", "#CC79A7")) # this is a CB friendly pal from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
bgcPlot2d
ggsave(paste0(outdir,"1merFractions.PerBGCCategory.ReducedFrom7merKSFS.SubSetOfSpecies.plusGOC.png"),bgcPlot2d,height=4,width=12)


########### troubleshoot fin whale ########
fw_plot1 <- ggplot(allSFSes_1mer[allSFSes_1mer$label=="fin_whale_ENP",],aes(x=sample_frequency,y=total_mutations_1mer,fill=mutation_type_1mer))+
  geom_col()
fw_plot1
ggsave(paste0(outdir,"troubleShoot.FinWhaleEnp.1.png"),fw_plot1,height=7,width=12)

fw_plot1b <- ggplot(allSFSes_1mer[allSFSes_1mer$label=="fin_whale_ENP",],aes(x=sample_frequency,y=mutationFractionPerFrequency,fill=mutation_type_1mer))+
  geom_col()
fw_plot1b
ggsave(paste0(outdir,"troubleShoot.FinWhaleEnp.1.png"),fw_plot1,height=7,width=12)

######## want to look for this weird periodicity in others species #########
labels=unique(allSFSes_1mer$label)
for(label in labels){
  spSpecPlot <- ggplot(allSFSes_1mer[allSFSes_1mer$label==label,],aes(x=alleleFrequency,y=mutationFractionPerFrequency,fill=mutation_type_1mer))+
    geom_col()
  ggsave(paste0(outdir,"speciesSpecificPlot.",label,".png"),spSpecPlot,height=7,width=12)
  
}
