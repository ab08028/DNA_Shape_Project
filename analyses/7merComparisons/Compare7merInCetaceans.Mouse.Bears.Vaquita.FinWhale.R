require(ggplot2)
require(reshape2)
require(dplyr)
require(ggrepel)
require(gtools)
############ these are now iwthout the --strict flag (good ) ###########
plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/"
########### combine cetacean spectra with vaquita and mouse ###########
cetaceans <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_spectrum_files/noSTRICTflag/AllCetaceans_melted.IncludesTwoConditions.HetsOnly.and.HetsPlusHoms.NOSTRICTFLAG.txt",header=T,sep="\t")
# adding cols with diff names for convenience to bind rows with others:
cetaceans$population <- cetaceans$Species_common # for convenience
cetaceans$fractionOfAllSegSites <- cetaceans$proportionOfSites

# note samples that were mapped to diff reference genome
cetaceans$MappedTo <- "own species" # most were mapped to their own ref genome with exceptions:
# some were mapped to vaquita instead of to themselves
cetaceans[cetaceans$Reference=="mPhoSin1.pri.cur.20190723",]$MappedTo <- "vaquita (no softmask)"
cetaceans[cetaceans$Species_common=="Indo-Pacific finless porpoise" & cetaceans$Reference=="GCF_003031525.1_Neophocaena_asiaeorientalis_V1",]$MappedTo <- "Yangtze finless porpoise"

### need to add vaquita and fin whale with no strict to this 
cetaceans$centralMutationType <- NA
# "AAAAAAA.AAACAAA" so the position of the central mutations are 4 and 12 (remember ther's a dot in between')
# and 3-5 is teh 3mer context and 11-13
cetaceans$centralMutationType <- paste(substr(cetaceans$variable,4,4),substr(cetaceans$variable,12,12),sep=".")
cetaceans$ThreemerMutationType <- NA
cetaceans$ThreemerMutationType <- paste(substr(cetaceans$variable,3,5),substr(cetaceans$variable,11,13),sep=".")

### beluga is weird and buggy! has a 7mer that is cutoff somehow.
#### NOTE contains hetsonly or hetsplushomalt -- just want hets for now... 
########## adding in vaquita 
vaquita <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210225_7mer/vaquita.7merspectrum.SummedOverAllIntervals.txt",header=T,sep="\t")
vaquita$MappedTo <- "vaquita (no softmask)"
vaquita$population <- "Vaquita (population)"
vaquita$sample <- "all"
##### adding in mouse:
mouse <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/summaryTables/Mouse7merSpectrum.SummedUpOverAllChrs.NOSTRICT.txt",header=T,sep="\t")
mouse$MappedTo <- "mm10 (nostrict)"
## this is the NOSTRICT version
mouse$sample <- "all"
####### humans (no strict) -- all populations #######
humans <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/humans/mutyperResults_20210329/mutyper_spectrum_files/ALLINTERVALS.human.PerPopulation.mutyper.spectra.7mer.txt",header=T,sep="\t")
humans$MappedTo <- "human (nostrict_nomask)"
humans$sample <- "all"


### bears (no strict) ####
bears <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210331_7mer_NoStrict/mapped_to_brown_bear/summaryTables/mapped_to_brown_bear.AllPops.bears.7merspectrum.SummedOverAllIntervals.ALLFREQS.NOSTRICT.txt",header=T,sep="\t")
bears$MappedTo <- "brown_bear"
bears$sample <- "all"
bears$population <- as.character(bears$population)
# need to rename because 
bears[bears$population=="ABC",]$population <- "ABC_Bear"
bears[bears$population=="EUR",]$population <- "EUR_Bear"
bears[bears$population=="PB",]$population <- "Polar_Bear"

######## combine ##########
combo <- bind_rows(mouse[,c("population","variable","centralMutationType","ThreemerMutationType","fractionOfAllSegSites","MappedTo","sample")],vaquita[,c("population","variable","centralMutationType","ThreemerMutationType","fractionOfAllSegSites","MappedTo","sample")],cetaceans[cetaceans$condition=="hetsonly",c("population","variable","centralMutationType","ThreemerMutationType","fractionOfAllSegSites","MappedTo","sample")],humans[,c("population","variable","centralMutationType","ThreemerMutationType","fractionOfAllSegSites","MappedTo","sample")],bears[,c("population","variable","centralMutationType","ThreemerMutationType","fractionOfAllSegSites","MappedTo","sample")])

# want to order them approx phylogenetically
combo$population <- factor(combo$population,levels=c("Polar_Bear","EUR_Bear","ABC_Bear","AFR","EUR","EAS","Ms","Mmm","Mmd","Mmc","Minke whale","Blue whale","Sperm whale","Bottlenose dolphin","Pacific white-sided dolphin","Killer whale/orca","Long-finned pilot whale","Vaquita (population)","Harbor porpoise","Indo-Pacific finless porpoise","Yangtze finless porpoise","Beluga whale","Narwhal"))


########### need to add: ######


#### fin whale -- need genome wide snps with no strict #######

####### plot ##########
### Note that some individuals (npho1) are mapped to two different references ! and therefore appears twice (only one I think )

# compare just that 7mer
p1 <- ggplot(combo[combo$variable=="TTTAAAA.TTTTAAA",],aes(y=population,x=fractionOfAllSegSites,color=MappedTo))+
  geom_point()+
  theme_bw()+
  ggtitle("Proportion of segregating sites that are TTTAAAA.TTTTAAA")
p1

ggsave(paste(plotdir,"ProportionOfHets.TTTAAAA.TTTTAAA.InCetaceans.PlusMice.PlusHumans.PlusBears.png",sep=""),height=6,width=8)

############## want to calculate the amount of variance for each variable (7mer) across species)
varianceOf7mers <- combo %>%
  group_by(variable) %>%
  summarise(varianceOf7merAcrossSpecies=var(fractionOfAllSegSites),MeanAbundanceOf7MerAcrossSpecies=mean(fractionOfAllSegSites))

varplot <- ggplot(varianceOf7mers,aes(x=varianceOf7merAcrossSpecies))+
  geom_histogram(bins=1000)+
  scale_x_log10()+
  #geom_vline(data=varianceOf7mers[varianceOf7mers$varianceOf7merAcrossSpecies>1e-7,],aes(xintercept=varianceOf7merAcrossSpecies))+
  geom_label_repel(data=varianceOf7mers[varianceOf7mers$varianceOf7merAcrossSpecies>2e-7 | varianceOf7mers$varianceOf7merAcrossSpecies<1e-13,], aes(y=1,label=variable), size = 3, box.padding = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  ggtitle("variance of each 7mer across cetacean species and mouse and humans")+
  theme_bw()
varplot

ggsave(paste(plotdir,"Variance7mers.CetaceansPlusMousePlusHumans.png",sep=""),varplot,height=5,width=7)
######### plotting mean abundance ######
meanabundanceplot <- ggplot(varianceOf7mers,aes(x=MeanAbundanceOf7MerAcrossSpecies))+
  geom_histogram(bins=1000)+
  scale_x_log10()+
  geom_label_repel(data=varianceOf7mers[varianceOf7mers$MeanAbundanceOf7MerAcrossSpecies>1e-3 | varianceOf7mers$MeanAbundanceOf7MerAcrossSpecies<5e-7,], aes(y=1,label=variable), size = 3, box.padding = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  ggtitle("mean abundance of each 7mer across cetacean species and mouse and human")+
  theme_bw()

meanabundanceplot
ggsave(paste(plotdir,"MeanAbundance7mers.CetaceansPlusMouse.PlusHuman.png",sep=""),meanabundanceplot,height=5,width=7)

####### abundance vs variance ########
meanvarplot <- ggplot(varianceOf7mers,aes(x=MeanAbundanceOf7MerAcrossSpecies,y=sqrt(varianceOf7merAcrossSpecies)))+
  geom_point()+
  geom_label_repel(data=varianceOf7mers[varianceOf7mers$MeanAbundanceOf7MerAcrossSpecies>1e-3 | varianceOf7mers$MeanAbundanceOf7MerAcrossSpecies<5e-7,], aes(y=sqrt(varianceOf7merAcrossSpecies),label=variable), size = 3, box.padding = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  theme_bw()+
  ggtitle("std deviation and mean abundance (across cetacean species and mouse and humans)")

meanvarplot
ggsave(paste(plotdir,"MeanVariance7mers.CetaceansPlusMouse.png",sep=""),meanvarplot,height=5,width=7)

########## want to just plot cetaceans or just mice #############
varianceOf7mers_mouseOnly <- combo %>%
  subset(population %in% c("Mmc","Mmd","Ms","Mmm")) %>%
  group_by(variable) %>%
  summarise(varianceOf7merAcrossSpecies=var(fractionOfAllSegSites),MeanAbundanceOf7MerAcrossSpecies=mean(fractionOfAllSegSites))
  
varianceOf7mers_cetaceansOnly <- combo %>%
  subset(!(population %in% c("Mmc","Mmd","Ms","Mmm"))) %>%
  group_by(variable) %>%
  summarise(varianceOf7merAcrossSpecies=var(fractionOfAllSegSites),MeanAbundanceOf7MerAcrossSpecies=mean(fractionOfAllSegSites))

# don't have SAS currently
varianceOf7mers_humansOnly <- combo %>%
  subset(!(population %in% c("EAS","AFR","EUR","SAS"))) %>%
  group_by(variable) %>%
  summarise(varianceOf7merAcrossSpecies=var(fractionOfAllSegSites),MeanAbundanceOf7MerAcrossSpecies=mean(fractionOfAllSegSites))

############## plot mouse only #############
meanvarplot_mouseOnly <- ggplot(varianceOf7mers_mouseOnly,aes(x=MeanAbundanceOf7MerAcrossSpecies,y=sqrt(varianceOf7merAcrossSpecies)))+
  geom_point()+
  geom_label_repel(data=varianceOf7mers_mouseOnly[sqrt(varianceOf7mers_mouseOnly$varianceOf7merAcrossSpecies)> .00015 | varianceOf7mers_mouseOnly$MeanAbundanceOf7MerAcrossSpecies>0.0005,], aes(y=sqrt(varianceOf7merAcrossSpecies),label=variable), size = 3, box.padding = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  theme_bw()+
  ggtitle("std deviation and mean abundance (mouse only)")

meanvarplot_mouseOnly
ggsave(paste(plotdir,"MeanVariance7mers.MOUSEONLY.png",sep=""),meanvarplot_mouseOnly,height=5,width=7)


########## plot cetaceans only ###########
meanvarplot_cetaceansOnly <- ggplot(varianceOf7mers_cetaceansOnly,aes(x=MeanAbundanceOf7MerAcrossSpecies,y=sqrt(varianceOf7merAcrossSpecies)))+
  geom_point()+
  geom_label_repel(data=varianceOf7mers_cetaceansOnly[sqrt(varianceOf7mers_cetaceansOnly$varianceOf7merAcrossSpecies)> .001 | varianceOf7mers_cetaceansOnly$MeanAbundanceOf7MerAcrossSpecies>0.001,], aes(y=sqrt(varianceOf7merAcrossSpecies),label=variable), size = 3, box.padding = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  theme_bw()+
  ggtitle("std deviation and mean abundance (cetaceans only)")

meanvarplot_cetaceansOnly
ggsave(paste(plotdir,"MeanVariance7mers.CETACEANSONLY.png",sep=""),meanvarplot_cetaceansOnly,height=5,width=7)


######### plot humans only ##############
meanvarplot_humanOnly <- ggplot(varianceOf7mers_humansOnly,aes(x=MeanAbundanceOf7MerAcrossSpecies,y=sqrt(varianceOf7merAcrossSpecies)))+
  geom_point()+
  geom_label_repel(data=varianceOf7mers_humansOnly[sqrt(varianceOf7mers_humansOnly$varianceOf7merAcrossSpecies)> .00075 | varianceOf7mers_humansOnly$MeanAbundanceOf7MerAcrossSpecies>0.00075,], aes(y=sqrt(varianceOf7merAcrossSpecies),label=variable), size = 3, box.padding = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  theme_bw()+
  ggtitle("std deviation and mean abundance (human populations only)")+
  xlab("mean frequency (out of all seg sites) in human populations")+
  ylab("std deviation of frequency among human populations")

meanvarplot_humanOnly
ggsave(paste(plotdir,"MeanVariance7mers.HUMANONLY.png",sep=""),meanvarplot_humanOnly,height=5,width=7)

############# make new diagonal plot with whale vs mouse ###############
## just going to use Mmd for now 
# choosing just one minke whale (bacu2) for now
Mmd_mouse_plus_minke_bacu2 <- merge(combo[combo$sample=="bacu2",c("variable", "ThreemerMutationType","centralMutationType","fractionOfAllSegSites")],combo[combo$population=="Mmd",c("variable", "ThreemerMutationType","centralMutationType" ,"fractionOfAllSegSites")],by=c("variable","centralMutationType","ThreemerMutationType"),suffixes = c(".minke",".Mmd"))

Mmd_minke_diagonalPlot <- ggplot(Mmd_mouse_plus_minke_bacu2,aes(x=fractionOfAllSegSites.Mmd,y=fractionOfAllSegSites.minke,color=centralMutationType,label=variable))+
  geom_point(size=1.2,alpha=0.7)+
  geom_abline(intercept=0,slope=1)+
  theme_bw() +
  geom_label_repel(data=subset(Mmd_mouse_plus_minke_bacu2,fractionOfAllSegSites.minke>=0.002),size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  xlab("Mouse 7mer proportion of all segregating sites (M. musculus domesticus)")+
  ylab("Minke whale 7mer proportion of all segregating sites")+
  ggtitle("Comparing mouse and minke whale 7mer proportions")+
  theme(text=element_text(size=14))+
  guides(colour = guide_legend(override.aes = list(size=3))) # this makes the dots in the legend larger
Mmd_minke_diagonalPlot
ggsave(paste(plotdir,"Diagonal.Mmdmouse.vs.Minke.bacu2.Comparison.7mers.NOSTRICT.png",sep=""),Mmd_minke_diagonalPlot,height=7,width=12)


############ comparing minke whale and dolphin ###############

bottlenose_dolphin_plus_minke_bacu2 <- merge(combo[combo$sample=="bacu2",c("variable", "ThreemerMutationType","centralMutationType","fractionOfAllSegSites")],combo[combo$sample=="ttru1",c("variable", "ThreemerMutationType","centralMutationType" ,"fractionOfAllSegSites")],by=c("variable","centralMutationType","ThreemerMutationType"),suffixes = c(".minke",".dolphin"))

dolphin_minke_diagonalPlot <- ggplot(bottlenose_dolphin_plus_minke_bacu2,aes(x=fractionOfAllSegSites.dolphin,y=fractionOfAllSegSites.minke,color=centralMutationType,label=variable))+
  geom_point(size=1.2,alpha=0.7)+
  geom_abline(intercept=0,slope=1)+
  theme_bw() +
  geom_label_repel(data=subset(bottlenose_dolphin_plus_minke_bacu2,fractionOfAllSegSites.minke>=0.002),size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  xlab("Bottlenose dolphin 7mer proportion of all segregating sites")+
  ylab("Minke whale 7mer proportion of all segregating sites")+
  ggtitle("Comparing bottlenose dolphin and minke whale 7mer proportions")+
  theme(text=element_text(size=14))+
  guides(colour = guide_legend(override.aes = list(size=3))) # this makes the dots in the legend larger
dolphin_minke_diagonalPlot

ggsave(paste(plotdir,"Diagonal.BottlenoseDolphin.vs.Minke.bacu2.Comparison.7mers.png",sep=""),dolphin_minke_diagonalPlot,height=7,width=12)
 
############## vaquita vs mouse (no strict) ############
Mmd_mouse_plus_vaquita_pop <- merge(combo[combo$population=="Vaquita (population)",c("variable", "ThreemerMutationType","centralMutationType","fractionOfAllSegSites")],combo[combo$population=="Mmd",c("variable", "ThreemerMutationType","centralMutationType" ,"fractionOfAllSegSites")],by=c("variable","centralMutationType","ThreemerMutationType"),suffixes = c(".vaquita",".Mmd"))

Mmd_vaquita_diagonalPlot <- ggplot(Mmd_mouse_plus_vaquita_pop,aes(x=fractionOfAllSegSites.Mmd,y=fractionOfAllSegSites.vaquita,color=centralMutationType,label=variable))+
  geom_point(size=1.2,alpha=0.7)+
  geom_abline(intercept=0,slope=1)+
  theme_bw() +
  geom_label_repel(data=subset(Mmd_mouse_plus_vaquita_pop,fractionOfAllSegSites.vaquita>=0.002),size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  xlab("Mouse 7mer proportion of all segregating sites (M. musculus domesticus)")+
  ylab("Vaquita 7mer proportion of all segregating sites")+
  ggtitle("Comparing mouse and vaquita 7mer proportions")+
  theme(text=element_text(size=14))+
  guides(colour = guide_legend(override.aes = list(size=3))) # this makes the dots in the legend larger
Mmd_vaquita_diagonalPlot
ggsave(paste(plotdir,"Diagonal.Mmdmouse.vs.Minke.bacu2.Comparison.7mers.NOSTRICT.png",sep=""),Mmd_minke_diagonalPlot,height=7,width=12)
