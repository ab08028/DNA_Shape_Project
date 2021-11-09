####### plotting proportional mutaiton signatures:

require(dplyr)
require(ggplot2)
require(gtools) # for mixedsort
require(reshape2)

sigprofilerRunLabel="20210825_SigProfilerExtractorResults_bearsMiceVaquitaFinWhales_kSFSAllFreqsSep_mm10GenomeRef"
wd=paste0("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/SigProfilerExtractor/",sigprofilerRunLabel,"/SBS96/")
plotdir=paste0(wd,"myplots/")
suggestedSolutionDir=paste0(wd,"/Suggested_Solution/")
dir.create(plotdir,showWarnings = F)


cosmicCounts <- read.table(paste0(suggestedSolutionDir,"COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities_refit.txt"),header=T)
head(cosmicCounts)
# temoproary

denovoSBSCounts <- read.table(paste0(suggestedSolutionDir,"SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt"),header=T)
head(denovoSBSCounts)



# want to plot them proportionally instead of counts

denovoSBSCounts_melt <- melt(denovoSBSCounts)
head(denovoSBSCounts_melt)

# convert coutns to proportions of total sites
denovoSBSCounts_melt <- denovoSBSCounts_melt %>% 
  group_by(Samples) %>%
  mutate(proportionOfSites=value/sum(value))


cosmicCounts_melt <- melt(cosmicCounts)
cosmicCounts_melt <- cosmicCounts_melt %>%
  group_by(Samples) %>%
  mutate(proportionOfSites=value/sum(value))



head(denovoSBSCounts_melt)
head(cosmicCounts_melt)
# these orderings are specific 
cosmicCounts_melt$Samples <- factor(cosmicCounts_melt$Samples,levels=mixedsort(unique(cosmicCounts_melt$Samples)))

denovoSBSCounts_melt$Samples <- factor(denovoSBSCounts_melt$Samples,levels=mixedsort(unique(denovoSBSCounts_melt$Samples)))

denovoSBSCounts_melt$species <- unlist(lapply(strsplit(as.character(denovoSBSCounts_melt$Samples),"_"),"[",1))

denovoSBSCounts_melt$sample_frequency <- as.numeric(unlist(lapply(strsplit(as.character(denovoSBSCounts_melt$Samples),"_"),"[",2))) ## be careful here to make numeric otherwise it messes up ordering of plot and wuld give wrong 'max'

# convert to allele frequency -- did I do this right???? singleton is 1/15 if unfolded or 1/30 if folded? 
denovoSBSCounts_melt <- denovoSBSCounts_melt %>%
  group_by(species,variable) %>%
  mutate(sample_freq_rescaled=sample_frequency/(max(sample_frequency)+1))
# adding 1 to max sample freq because it doesn't include fixed sites so max of 15 actually represents 16 haploids (8 diploids) so you want allele freq to be sample_freq/total haploids aka 1/16 not 1/15 (took this code from the bear BGC plotting code)


cosmicCounts_melt$species <- unlist(lapply(strsplit(as.character(cosmicCounts_melt$Samples),"_"),"[",1))

cosmicCounts_melt$sample_frequency <- as.numeric(unlist(lapply(strsplit(as.character(cosmicCounts_melt$Samples),"_"),"[",2))) ## be careful here to make numeric otherwise it messes up ordering of plot and wuld give wrong 'max'

# convert to allele frequency -- did I do this right???? singleton is 1/15 if unfolded or 1/30 if folded? 
cosmicCounts_melt <- cosmicCounts_melt %>%
  group_by(species,variable) %>%
  mutate(sample_freq_rescaled=sample_frequency/(max(sample_frequency)+1))
# adding 1 to max sample freq because it doesn't include fixed sites so max of 15 actually represents 16 haploids (8 diploids) so you want allele freq to be sample_freq/total haploids aka 1/16 not 1/15 (took this code from the bear BGC plotting code)


head(denovoSBSCounts_melt)
head(cosmicCounts_melt)

plot1a <- ggplot(denovoSBSCounts_melt,aes(fill=variable,y=proportionOfSites,x=sample_frequency))+
  geom_col(position="stack")+
  theme_bw()+
  facet_wrap(~species,scales = "free",ncol=1)+
  ggtitle("proportion of sites made up of each de novo signature across the SFS")
plot1a
ggsave(paste0(plotdir,"ProportionalSignatures.deNovo.png"),plot1a,height=16,width=10)

plot2a <- ggplot(cosmicCounts_melt,aes(fill=variable,y=proportionOfSites,x=sample_frequency))+
  geom_col(position="stack")+
  theme_bw()+
  facet_wrap(~species,scales = "free",ncol=1)+
  ggtitle("proportion of sites made up of each COSMIC signature")
plot2a
ggsave(paste0(plotdir,"ProportionalSignatures.COSMIC.png"),plot2a,height=16,width=10)


######### plot how signatures change in their proportion of sites across the SFS -- all species grouped together ###########
plot3 <- ggplot(denovoSBSCounts_melt,aes(color=species,y=proportionOfSites,x=sample_freq_rescaled))+
  geom_line()+
  theme_bw()+
  facet_wrap(~variable,scales = "free")+
  ggtitle("Changes in proportion of sites from each de novo signature per species across allele frequencies")
plot3
ggsave(paste0(plotdir,"ProportionalSignatures.deNovo.Lines.AlleleFreqs.PerDeNovoSignature.png"),plot3,height=5,width=10)


plot4 <- ggplot(cosmicCounts_melt,aes(color=species,y=proportionOfSites,x=sample_freq_rescaled))+
  geom_line()+
  theme_bw()+
  facet_wrap(~variable,scales = "free")+
  ggtitle("Changes in proportion of sites from each de novo signature per species across allele frequencies")
plot4
ggsave(paste0(plotdir,"ProportionalSignatures.deNovo.Lines.AlleleFreqs.PerCOSMICSignature.png"),plot4,height=5,width=10)
######## nice! this is a cool plot ala mt dna ########
