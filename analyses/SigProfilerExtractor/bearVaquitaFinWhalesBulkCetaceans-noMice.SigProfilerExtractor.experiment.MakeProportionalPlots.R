####### plotting proportional mutaiton signatures:

require(dplyr)
require(ggplot2)
require(gtools) # for mixedsort
require(reshape2)

sigprofilerRunLabel="20210805_SigProfilerExtractorResults_bearsVaquitaFinWhaleBulkCetaceans_NoMice_GRCh37GenomeRef"
wd=paste0("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/SigProfilerExtractor/",sigprofilerRunLabel,"/SBS96/")
plotdir=paste0(wd,"myplots/")
suggestedSolutionDir=paste0(wd,"/Suggested_Solution/")
dir.create(plotdir,showWarnings = F)


cosmicCounts <- read.table(paste0(suggestedSolutionDir,"COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities_refit.txt"),header=T)
head(cosmicCounts)

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
cosmicCounts_melt$Samples <- factor(cosmicCounts_melt$Samples,levels=c("ABC","EUR","PB","ENP","GOC","Vaquita","bacu2","bacu3","bmus1","dleu1","dleu2","gmel1","lobl1","mmon1","nasi1","npho1","oorc1","oorc2","pmac1","pmac2","ppho1","ttru1","ttru2"))

denovoSBSCounts_melt$Samples <- factor(denovoSBSCounts_melt$Samples,levels=c("ABC","EUR","PB","ENP","GOC","Vaquita","bacu2","bacu3","bmus1","dleu1","dleu2","gmel1","lobl1","mmon1","nasi1","npho1","oorc1","oorc2","pmac1","pmac2","ppho1","ttru1","ttru2"))

plot1 <- ggplot(denovoSBSCounts_melt,aes(fill=variable,y=proportionOfSites,x=Samples))+
  geom_col(position="stack")+
  theme_bw()+
  ggtitle("proportion of sites made up of each de novo signature")
plot1
ggsave(paste0(plotdir,"ProportionalSignatures.deNovo.png"),plot1,height=5,width=10)

plot2 <- ggplot(cosmicCounts_melt,aes(fill=variable,y=proportionOfSites,x=Samples))+
  geom_col(position="stack")+
  theme_bw()+
  ggtitle("proportion of sites made up of each COSMIC signature")
plot2
ggsave(paste0(plotdir,"ProportionalSignatures.COSMIC.png"),plot2,height=5,width=10)
