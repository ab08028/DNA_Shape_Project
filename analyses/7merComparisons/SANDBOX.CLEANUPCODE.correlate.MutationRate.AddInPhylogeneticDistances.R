########## want to show similarities between mutation spectra across species ############
require(dplyr)
require(reshape2)
require(tidyr)
require(GGally)
require(broom)
require(ggrepel)
plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/diagonal_comparison_plots/"
# want 3mer 5mer 7mer and want to show correlations 
# want it to be corrected for targets ideally
# proportion of segregating sites ?
# or rescaled mut rate <-- do this way

# whales, mice, bears, (condor leave out for now), humans 

#allCetaceans # /Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_spectrum_files/noSTRICTflag/AllCetaceans_melted.IncludesTwoConditions.HetsOnly.and.HetsPlusHoms.NOSTRICTFLAG.txt

# cetacean targets /Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_targets_files

# maybe just pick a couple cetaceans?
# let's do

# 4 wild mice, (fin whale -- has weirdness with being neutral so let's skip for now), minke whale, vaquita, polar bear, brown bear (eur) [condor? human?]

## maybe for now we do a case study -- just show that two mouse species and two whale species and mouse + whale -- I think that might be good.

##################### mice ################
mouseSpectrum <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/summaryTables/Mouse7merSpectrum.SummedUpOverAllChrs.NOSTRICT.txt",header=T)











######### set of processing steps that are in common (some may differ) #######

# get anc 7mer:
mouseSpectrum$ancestral7mer <- unlist(lapply(strsplit(mouseSpectrum$variable,"\\."),"[",1))
mouseSpectrum$derived7mer <- unlist(lapply(strsplit(mouseSpectrum$variable,"\\."),"[",2))
head(mouseSpectrum)
tail(mouseSpectrum)
mouseTargets <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/summaryTables/ALLINTERVALS.mouse.mutyper.targets.7mer.nostrict.txt",header=T)

# merge them:
head(mouseTargets)
colnames(mouseTargets) <- c("ancestral7mer","targetCount_allIntervals")

mouseAll <- merge(mouseSpectrum,mouseTargets,by="ancestral7mer")



mouseAll$mutationRate_unscaled <- mouseAll$totalSitesAllIntervals / mouseAll$targetCount_allIntervals

mouseAll <- mouseAll %>%
  group_by(population) %>%
  mutate(mutationRate_rescaled = mutationRate_unscaled/sum(mutationRate_unscaled))

mouseAll_wide <- pivot_wider(mouseAll[,c("population","variable","mutationRate_rescaled")],id_cols = c("variable"),values_from = "mutationRate_rescaled",names_from = "population")

head(mouseAll_wide)
############# pick a whale #############
# ok let's just do minke what the hell ; for now? maybe human would be better but don't have targets on hand. maybe want to keep vaquita private?
whale1Choice="bacu2" # minke whale 
whale1ChoiceNiceLabel="minke whale"
allCetaceansSpectrum <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_spectrum_files/noSTRICTflag/AllCetaceans_melted.IncludesTwoConditions.HetsOnly.and.HetsPlusHoms.NOSTRICTFLAG.txt",sep="\t",header=T)

# restrict to just one whale for now
whale1Spectrum <- allCetaceansSpectrum[allCetaceansSpectrum$sample==whale1Choice & allCetaceansSpectrum$condition=="hetsonly",]

whale1Spectrum$ancestral7mer <- unlist(lapply(strsplit(whale1Spectrum$variable,"\\."),"[",1))


# species specific targets: 
whale1targetfilename=list.files("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_targets_files", pattern=whale1Choice,full.names = T)
whale1targetfilename #make sure this is right
whale1targets <- read.table(whale1targetfilename ,header=F)
# careful here; make sure it's the right file 
head(whale1targets)
colnames(whale1targets) <- c("ancestral7mer","targetCount_allIntervals")

whale1All <- merge(whale1Spectrum,whale1targets,by="ancestral7mer")

whale1All$mutationRate_unscaled <- whale1All$value / whale1All$targetCount_allIntervals

whale1All <- whale1All %>%
  group_by(sample) %>%
  mutate(mutationRate_rescaled = mutationRate_unscaled/sum(mutationRate_unscaled))

whale1All_wide = pivot_wider(whale1All[,c("sample","variable","mutationRate_rescaled")],id_cols = c("variable"),values_from = "mutationRate_rescaled",names_from = "sample")
head(whale1All_wide)

######## pick a second whale #######
# bottlenose dolphin
whale2Choice="ttru1" # bottlenosedolphin 
whale2ChoiceNiceLabel="dolphin" # bottlenosedolphin 

# restrict to just one whale for now
whale2Spectrum <- allCetaceansSpectrum[allCetaceansSpectrum$sample==whale2Choice & allCetaceansSpectrum$condition=="hetsonly",]

whale2Spectrum$ancestral7mer <- unlist(lapply(strsplit(whale2Spectrum$variable,"\\."),"[",1))


# species specific targets: 
whale2targetfilename=list.files("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_targets_files", pattern=whale2Choice,full.names = T)
whale2targetfilename #make sure this is right
whale2targets <- read.table(whale2targetfilename ,header=F)
# careful here; make sure it's the right file 
head(whale2targets)
colnames(whale2targets) <- c("ancestral7mer","targetCount_allIntervals")

whale2All <- merge(whale2Spectrum,whale2targets,by="ancestral7mer")

whale2All$mutationRate_unscaled <- whale2All$value / whale2All$targetCount_allIntervals

whale2All <- whale2All %>%
  group_by(sample) %>%
  mutate(mutationRate_rescaled = mutationRate_unscaled/sum(mutationRate_unscaled))

whale2All_wide = pivot_wider(whale2All[,c("sample","variable","mutationRate_rescaled")],id_cols = c("variable"),values_from = "mutationRate_rescaled",names_from = "sample")
head(whale2All_wide)
############# add in bears ###########
bearSpectrum <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210331_7mer_NoStrict/mapped_to_brown_bear/summaryTables/mapped_to_brown_bear.AllPops.bears.7merspectrum.SummedOverAllIntervals.ALLFREQS.NOSTRICT.txt",header=T,sep="\t")
head(bearSpectrum)
bearSpectrum$ancestral7mer <- unlist(lapply(strsplit(bearSpectrum$variable,"\\."),"[",1))
head(bearSpectrum)
bearTargets <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210331_7mer_NoStrict/mapped_to_brown_bear/summaryTables/ALLINTERVALS.brown_bear.mutyper.targets.7mer.nostrict.txt",header=T)
head(bearTargets)
colnames(bearTargets) <- c("ancestral7mer","targetCount_allIntervals")

bearAll <- merge(bearSpectrum,bearTargets,by="ancestral7mer")

bearAll$mutationRate_unscaled <- bearAll$totalSitesAllIntervals / bearAll$targetCount_allIntervals


bearAll <- bearAll %>%
  group_by(population) %>%
  mutate(mutationRate_rescaled = mutationRate_unscaled/sum(mutationRate_unscaled))

head(bearAll)


bearAll_wide = pivot_wider(bearAll[,c("population","variable","mutationRate_rescaled")],id_cols = c("variable"),values_from = "mutationRate_rescaled",names_from = "population")
head(bearAll_wide)

########## okay let's do some plotting ###########
merge1 <- merge(mouseAll_wide,whale1All_wide,by="variable")
head(merge1)

merge2 <- merge(merge1,whale2All_wide,by="variable")

merge3 <- merge(merge2,bearAll_wide,by="variable")
# can add more merges as needed
head(merge3)

mergeFinal <- merge3

species=c("Mmd","Ms","EUR","PB","bacu2","ttru1")
# https://r-graphics.org/recipe-scatter-splom # #
# these are functions from the above cookbook ^ 
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex =  cex.cor * (1 + r) / 2)
}

# would fit a linear model
panel.lm <- function (x, y, col = par("col"), bg = NA, pch = par("pch"),
                      cex = 1, col.smooth = "red", ...) {
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  abline(stats::lm(y ~ x),  col = col.smooth, ...)
}

# fits the y=x line in red
panel.abline <- function (x, y, col = par("col"), bg = NA, pch = par("pch"),
                      cex = 1, col.smooth = "red", ...) {
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  abline(a=0,b=1,  col = col.smooth, ...)
}

png("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/diagonal_comparison_plots/pairsOfSpecies.20210908.forposter.png",width=9,height=8,res=300) 
pairs(select(mergeFinal,species),
      upper.panel =panel.cor,
      lower.panel= panel.abline,
      pch=".") # excludes points that are 0 
dev.off()

#### label CpG sites ########
mergeFinal$centralMutationType <- paste0(substr(mergeFinal$variable,4,4),".",substr(mergeFinal$variable,12,12))
mergeFinal$central3mer <- paste0(substr(mergeFinal$variable,3,5),".",substr(mergeFinal$variable,11,13))
head(mergeFinal)
tail(mergeFinal)

mergeFinal$ancestralDimer_forCpGLabeling <- substr(mergeFinal$central3mer,2,3) # get YZ from XYZ to check for the central bp being in a CpG (not just a CpG present anywhere in 3mer whic is how I used to do it -- I think this is more precise because specifically we expceted the center bp to mutate more if it's a C in a CpG not if it's a G in a CpG (I think))
head(mergeFinal)
mergeFinal$CpGPresentOrNot <- ""
mergeFinal[mergeFinal$ancestralDimer_forCpGLabeling=="CG",]$CpGPresentOrNot <- "CpG_"

mergeFinal$centralMutationType_includingCpG <- paste0(mergeFinal$CpGPresentOrNot,gsub("\\.",">",mergeFinal$centralMutationType))


####### want to plot mutation types separately? ########
MutTypes=unique(mergeFinal$centralMutationType_includingCpG)
allCorrelations_PerMutType_7mers = data.frame()
for(MutType in MutTypes){
  subset=filter(mergeFinal,centralMutationType_includingCpG==MutType)
  png(paste0(plotdir,"pairsOfSpecies.20210908.forposter.",MutType,".png"),width=7,height=6,units = "in",res=300)  # must specify res for png
  pairs(select(subset,species),
      upper.panel =panel.cor,
      lower.panel= panel.abline,
      pch=".",main=paste0(MutType," 7mers")) # excludes points that are 0 
  dev.off()

  # correlation matrix for this species type::
  cormat = cor(select(subset,species))
  ### Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  lower_tri_cormat <- get_lower_tri(cormat)
  cormat_melt <- melt(lower_tri_cormat,na.rm=T)
  colnames(cormat_melt) <- c("SpeciesA","SpeciesB","PearsonRsq")
  cormat_melt$label <- MutType
  allCorrelations_PerMutType_7mers <- rbind(allCorrelations_PerMutType_7mers,cormat_melt)
  
}



################ let's plot correlation coeffs by mutation type as a bar chart #######
# make nicer labels:

allCorrelations_PerMutType_7mers$nicerSpeciesLabelA <- ""

allCorrelations_PerMutType_7mers[allCorrelations_PerMutType_7mers$SpeciesA=="EUR",]$nicerSpeciesLabelA  <- "brown bear"
allCorrelations_PerMutType_7mers[allCorrelations_PerMutType_7mers$SpeciesA=="PB",]$nicerSpeciesLabelA  <- "polar bear"
allCorrelations_PerMutType_7mers[allCorrelations_PerMutType_7mers$SpeciesA==whale1Choice,]$nicerSpeciesLabelA  <- whale1ChoiceNiceLabel
allCorrelations_PerMutType_7mers[allCorrelations_PerMutType_7mers$SpeciesA==whale2Choice,]$nicerSpeciesLabelA  <- whale2ChoiceNiceLabel

allCorrelations_PerMutType_7mers[allCorrelations_PerMutType_7mers$SpeciesA=="Mmd",]$nicerSpeciesLabelA  <- "Mmd"

allCorrelations_PerMutType_7mers[allCorrelations_PerMutType_7mers$SpeciesA=="Ms",]$nicerSpeciesLabelA  <- "Ms"

allCorrelations_PerMutType_7mers$nicerSpeciesLabelB <- ""

allCorrelations_PerMutType_7mers[allCorrelations_PerMutType_7mers$SpeciesB=="EUR",]$nicerSpeciesLabelB  <- "brown bear"
allCorrelations_PerMutType_7mers[allCorrelations_PerMutType_7mers$SpeciesB=="PB",]$nicerSpeciesLabelB  <- "polar bear"
allCorrelations_PerMutType_7mers[allCorrelations_PerMutType_7mers$SpeciesB==whale1Choice,]$nicerSpeciesLabelB  <- whale1ChoiceNiceLabel
allCorrelations_PerMutType_7mers[allCorrelations_PerMutType_7mers$SpeciesB==whale2Choice,]$nicerSpeciesLabelB  <- whale2ChoiceNiceLabel
allCorrelations_PerMutType_7mers[allCorrelations_PerMutType_7mers$SpeciesB=="Mmd",]$nicerSpeciesLabelB  <- "Mmd"

allCorrelations_PerMutType_7mers[allCorrelations_PerMutType_7mers$SpeciesB=="Ms",]$nicerSpeciesLabelB  <- "Ms"

allCorrelations_PerMutType_7mers$nicerSpeciesLabelA <- factor(allCorrelations_PerMutType_7mers$nicerSpeciesLabelA,levels=c("Mmd","Ms","brown bear","polar bear","minke whale","dolphin"))


allCorrelations_PerMutType_7mers$nicerSpeciesLabelB <- factor(allCorrelations_PerMutType_7mers$nicerSpeciesLabelB,levels=c("Mmd","Ms","brown bear","polar bear","minke whale","dolphin"))

# nicer labels for mutation types
allCorrelations_PerMutType_7mers$NicerMutationlabel <- gsub("\\.",">",allCorrelations_PerMutType_7mers$label)
# specify if nonCpG 

########## all mut type correlation plot 7mers ########
corPlot7mer <- ggplot(data = allCorrelations_PerMutType_7mers[allCorrelations_PerMutType_7mers$SpeciesA!=allCorrelations_PerMutType_7mers$SpeciesB,], aes(x=nicerSpeciesLabelA, y=nicerSpeciesLabelB, fill = PearsonRsq)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(-.2,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  #scale_fill_gradient(low = "blue", high = "red", limit = c(-.2,1), space = "Lab", 
  #name="Pearson\nCorrelation") +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  facet_wrap(~NicerMutationlabel)+
  xlab("")+
  ylab("")+
  geom_text(aes(label=round(PearsonRsq,2)))+
  ggtitle("7mers")
corPlot7mer
ggsave(paste0(plotdir,"Correlations.AllMutTypes.Heatmap.7mers.png"),corPlot7mer,width=11,height=8)  # must specify res for png
# should I make a correlogram?

####### make a plot of a single species pair colored by mutation type ######
##### make a single plot colored by central mut type ######
mergeFinal$nicerLabel <- gsub("\\.",">",mergeFinal$centralMutationType_includingCpG)

# wow these plots are nice -- need to make a lot of them ## YOU ARE HERE
examplePlotColoredByMutationType <- ggplot(mergeFinal,aes(x=Mmd,y=Ms,color=centralMutationType_includingCpG))+
  geom_point()+
  facet_wrap(~nicerLabel,scales="free")+
  geom_abline(intercept=0,slope=1)+
  theme_bw()
# the log scaling isn't great here because you miss how different the high mut ones are 
examplePlotColoredByMutationType

# let's make a plot of all C>T 
##### need to plot 5mers tooo ######
# get 5mers:
##### bears
bearAll$FivemerMutationType <- paste0(substr(bearAll$variable,2,6),".",substr(bearAll$variable,10,14))

bearAll_5mer <- bearAll %>% group_by(population,FivemerMutationType) %>%
  summarise(FivemerMutationCount=sum(totalSitesAllIntervals),FivemerTargetCount=sum(targetCount_allIntervals),FivemerMutationRate_unscaled=FivemerMutationCount/FivemerTargetCount)

bearAll_5mer <- bearAll_5mer %>%
  group_by(population) %>%
  mutate(FivemerMutationRate_rescaled=FivemerMutationRate_unscaled/sum(FivemerMutationRate_unscaled))

bearAll_5mer_wide = pivot_wider(bearAll_5mer[,c("population","FivemerMutationType","FivemerMutationRate_rescaled")],id_cols = c("FivemerMutationType"),values_from = "FivemerMutationRate_rescaled",names_from = "population")
head(bearAll_5mer_wide)
dim(bearAll_5mer_wide)
####### mouse
mouseAll$FivemerMutationType <- paste0(substr(mouseAll$variable,2,6),".",substr(mouseAll$variable,10,14))

mouseAll_5mer <- mouseAll %>% group_by(population,FivemerMutationType) %>%
  summarise(FivemerMutationCount=sum(totalSitesAllIntervals),FivemerTargetCount=sum(targetCount_allIntervals),FivemerMutationRate_unscaled=FivemerMutationCount/FivemerTargetCount)

mouseAll_5mer <- mouseAll_5mer %>%
  group_by(population) %>%
  mutate(FivemerMutationRate_rescaled=FivemerMutationRate_unscaled/sum(FivemerMutationRate_unscaled))


mouseAll_5mer_wide = pivot_wider(mouseAll_5mer[,c("population","FivemerMutationType","FivemerMutationRate_rescaled")],id_cols = c("FivemerMutationType"),values_from = "FivemerMutationRate_rescaled",names_from = "population")
head(mouseAll_5mer_wide)

###### whale 1
whale1All$FivemerMutationType <- paste0(substr(whale1All$variable,2,6),".",substr(whale1All$variable,10,14))


whale1All_5mer <- whale1All %>% group_by(sample,FivemerMutationType) %>%
  summarise(FivemerMutationCount=sum(value),FivemerTargetCount=sum(targetCount_allIntervals),FivemerMutationRate_unscaled=FivemerMutationCount/FivemerTargetCount)

whale1All_5mer <- whale1All_5mer %>%
  group_by(sample) %>%
  mutate(FivemerMutationRate_rescaled=FivemerMutationRate_unscaled/sum(FivemerMutationRate_unscaled))

whale1All_5mer_wide = pivot_wider(whale1All_5mer[,c("sample","FivemerMutationType","FivemerMutationRate_rescaled")],id_cols = c("FivemerMutationType"),values_from = "FivemerMutationRate_rescaled",names_from = "sample")
head(mouseAll_5mer_wide)

###### whale2
whale2All$FivemerMutationType <- paste0(substr(whale2All$variable,2,6),".",substr(whale2All$variable,10,14))

whale2All_5mer <- whale2All %>% group_by(sample,FivemerMutationType) %>%
  summarise(FivemerMutationCount=sum(value),FivemerTargetCount=sum(targetCount_allIntervals),FivemerMutationRate_unscaled=FivemerMutationCount/FivemerTargetCount)

whale2All_5mer <- whale2All_5mer %>%
  group_by(sample) %>%
  mutate(FivemerMutationRate_rescaled=FivemerMutationRate_unscaled/sum(FivemerMutationRate_unscaled))

whale2All_5mer_wide = pivot_wider(whale2All_5mer[,c("sample","FivemerMutationType","FivemerMutationRate_rescaled")],id_cols = c("FivemerMutationType"),values_from = "FivemerMutationRate_rescaled",names_from = "sample")
head(whale2All_5mer_wide)
############# merge the 5mers ############
merge5mer1 <- merge(mouseAll_5mer_wide,bearAll_5mer_wide,by="FivemerMutationType")
merge5mer2 <- merge(merge5mer1,whale1All_5mer_wide,by="FivemerMutationType")
merge5mer3 <- merge(merge5mer2,whale2All_5mer_wide,by="FivemerMutationType")

merge5mer_final <- merge5mer3
head(merge5mer_final)


######### plot 5mers #########
png("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/diagonal_comparison_plots/pairsOfSpecies.20210908.forposter.5mers.png",width=9,height=8,units = "in",res=300) 
pairs(select(merge5mer_final,species),
      upper.panel =panel.cor,
      lower.panel= panel.abline,
      pch=".") # excludes points that are 0 
dev.off()

#### label CpG sites ########
merge5mer_final$centralMutationType <- paste0(substr(merge5mer_final$FivemerMutationType,3,3),".",substr(merge5mer_final$FivemerMutationType,9,9))

merge5mer_final$central3mer <- paste0(substr(merge5mer_final$FivemerMutationType,2,4),".",substr(merge5mer_final$FivemerMutationType,8,10))
head(merge5mer_final)
tail(merge5mer_final)

merge5mer_final$ancestralDimer_forCpGLabeling <- substr(merge5mer_final$central3mer,2,3) # get YZ from XYZ to check for the central bp being in a CpG (not just a CpG present anywhere in 3mer whic is how I used to do it -- I think this is more precise because specifically we expceted the center bp to mutate more if it's a C in a CpG not if it's a G in a CpG (I think))
head(merge5mer_final)
merge5mer_final$CpGPresentOrNot <- ""
merge5mer_final[merge5mer_final$ancestralDimer_forCpGLabeling=="CG",]$CpGPresentOrNot <- "CpG_"

merge5mer_final$centralMutationType_includingCpG <- paste0(merge5mer_final$CpGPresentOrNot,gsub("\\.",">",merge5mer_final$centralMutationType))

head(merge5mer_final)
####### want to plot mutation types separately? ########
MutTypes=unique(merge5mer_final$centralMutationType_includingCpG)
allCorrelations_PerMutType_5mers = data.frame()
for(MutType in MutTypes){
  subset=filter(merge5mer_final,centralMutationType_includingCpG==MutType)
  png(paste0(plotdir,"pairsOfSpecies.20210908.forposter.",MutType,"5mers.png"),width=7,height=6,units = "in",res=300)  # must specify res for png
  pairs(select(subset,species),
        upper.panel =panel.cor,
        lower.panel= panel.abline,
        pch=".",main=paste0(MutType," 5mers")) # excludes points that are 0 
  dev.off()
  
  # correlation matrix for this species type::
  cormat = cor(select(subset,species))
  ### Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  lower_tri_cormat <- get_lower_tri(cormat)
  cormat_melt <- melt(lower_tri_cormat,na.rm=T)
  colnames(cormat_melt) <- c("SpeciesA","SpeciesB","PearsonRsq")
  cormat_melt$label <- MutType
  allCorrelations_PerMutType_5mers <- rbind(allCorrelations_PerMutType_5mers,cormat_melt)
  
}





################ let's plot correlation coeffs by mutation type as a bar chart #######
# make nicer labels:

allCorrelations_PerMutType_5mers$nicerSpeciesLabelA <- ""

allCorrelations_PerMutType_5mers[allCorrelations_PerMutType_5mers$SpeciesA=="EUR",]$nicerSpeciesLabelA  <- "brown bear"
allCorrelations_PerMutType_5mers[allCorrelations_PerMutType_5mers$SpeciesA=="PB",]$nicerSpeciesLabelA  <- "polar bear"
allCorrelations_PerMutType_5mers[allCorrelations_PerMutType_5mers$SpeciesA==whale1Choice,]$nicerSpeciesLabelA  <- whale1ChoiceNiceLabel
allCorrelations_PerMutType_5mers[allCorrelations_PerMutType_5mers$SpeciesA==whale2Choice,]$nicerSpeciesLabelA  <- whale2ChoiceNiceLabel

allCorrelations_PerMutType_5mers[allCorrelations_PerMutType_5mers$SpeciesA=="Mmd",]$nicerSpeciesLabelA  <- "Mmd"

allCorrelations_PerMutType_5mers[allCorrelations_PerMutType_5mers$SpeciesA=="Ms",]$nicerSpeciesLabelA  <- "Ms"

allCorrelations_PerMutType_5mers$nicerSpeciesLabelB <- ""

allCorrelations_PerMutType_5mers[allCorrelations_PerMutType_5mers$SpeciesB=="EUR",]$nicerSpeciesLabelB  <- "brown bear"
allCorrelations_PerMutType_5mers[allCorrelations_PerMutType_5mers$SpeciesB=="PB",]$nicerSpeciesLabelB  <- "polar bear"
allCorrelations_PerMutType_5mers[allCorrelations_PerMutType_5mers$SpeciesB==whale1Choice,]$nicerSpeciesLabelB  <- whale1ChoiceNiceLabel
allCorrelations_PerMutType_5mers[allCorrelations_PerMutType_5mers$SpeciesB==whale2Choice,]$nicerSpeciesLabelB  <- whale2ChoiceNiceLabel
allCorrelations_PerMutType_5mers[allCorrelations_PerMutType_5mers$SpeciesB=="Mmd",]$nicerSpeciesLabelB  <- "Mmd"

allCorrelations_PerMutType_5mers[allCorrelations_PerMutType_5mers$SpeciesB=="Ms",]$nicerSpeciesLabelB  <- "Ms"

allCorrelations_PerMutType_5mers$nicerSpeciesLabelA <- factor(allCorrelations_PerMutType_5mers$nicerSpeciesLabelA,levels=c("Mmd","Ms","brown bear","polar bear","minke whale","dolphin"))


allCorrelations_PerMutType_5mers$nicerSpeciesLabelB <- factor(allCorrelations_PerMutType_5mers$nicerSpeciesLabelB,levels=c("Mmd","Ms","brown bear","polar bear","minke whale","dolphin"))


allCorrelations_PerMutType_5mers$NicerMutationlabel <- gsub("\\.",">",allCorrelations_PerMutType_5mers$label)

########## all mut types correlation plot 5mers ###############
corPlot5mer <- ggplot(data = allCorrelations_PerMutType_5mers[allCorrelations_PerMutType_5mers$SpeciesA!=allCorrelations_PerMutType_5mers$SpeciesB,], aes(x=nicerSpeciesLabelA, y=nicerSpeciesLabelB, fill = PearsonRsq)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(-.2,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  #scale_fill_gradient(low = "blue", high = "red", limit = c(-.2,1), space = "Lab", 
                       #name="Pearson\nCorrelation") +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  facet_wrap(~label)+
  xlab("")+
  ylab("")+
  geom_text(aes(label=round(PearsonRsq,2)))+
  ggtitle("5mers")
corPlot5mer
ggsave(paste0(plotdir,"Correlations.AllMutTypes.Heatmap.5mers.png"),corPlot5mer,width=9,height=7)  # must specify res for png



##### quickplot for KH #######
# are the C>A mouse outliers in SBS18?

# let's see:
ggplot(merge5mer_final[merge5mer_final$centralMutationType=="C.A",],aes(x=Mmd,y=EUR,label=FivemerMutationType))+
  geom_point()+
  theme_bw()+
  geom_text()

ggplot(merge5mer_final[merge5mer_final$centralMutationType=="C.A",],aes(x=Mmd,y=Ms,label=FivemerMutationType))+
  geom_point()+
  theme_bw()+
  geom_text()

############################ Adding in phylogenetic distance #################
# need to get proper distance with clr transform soon
# start with euclidean distance first 
head(merge5mer_final)
  
?dist # uses euclidean distance by default. doing this for now then move on to clr once I understand it better

# okay so dist works on ROWS of a matrix but my data are in columns
# so I could transpose?
matrix7mer <- t(select(mergeFinal,-c(variable,centralMutationType,central3mer,ancestralDimer_forCpGLabeling,CpGPresentOrNot,centralMutationType_includingCpG,nicerLabel)))

distances7mer <- tidy(dist(matrix7mer)) # tidy with broom packge 
head(distances7mer)
colnames(distances7mer) <- c("sp1","sp2","spectrum_euclidean_distance_7mer")
#distances7mer$label <- "7mer"
distances7mer$comparisonLabel1 <- paste0(pmin(as.character(distances7mer$sp1),as.character(distances7mer$sp2)),".",pmax(as.character(distances7mer$sp1),as.character(distances7mer$sp2)))

matrix5mer <- t(select(merge5mer_final,-c(FivemerMutationType,centralMutationType,central3mer,ancestralDimer_forCpGLabeling,CpGPresentOrNot,centralMutationType_includingCpG)))


distances5mer <- tidy(dist(matrix5mer)) # tidy with broom packge 
head(distances5mer)
colnames(distances5mer) <- c("sp1","sp2","spectrum_euclidean_distance_5mer")
#distances5mer$label <- "5mer"
distances5mer$comparisonLabel1 <- paste0(pmin(as.character(distances5mer$sp1),as.character(distances5mer$sp2)),".",pmax(as.character(distances5mer$sp1),as.character(distances5mer$sp2)))



distances_7mer_5mer <- merge(distances7mer,distances5mer,by=c("comparisonLabel1"))
head(distances_7mer_5mer)
# why did some drop out. 

ggplot(distances_7mer_5mer,aes(x=spectrum_euclidean_distance_7mer,y=spectrum_euclidean_distance_5mer,label=comparisonLabel1))+
  geom_point()+
  geom_abline()+
  geom_text_repel()
# need to give them proper species names to match with 
###### read in phylo distances: 
phyloDistances = read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/raxml_tree/raxml.cophenetic.pairwisedistances.fromUphametal.MYSPECIES.usethis.txt",header=T)
head(phyloDistances)

# way to get species names:
# THESE ARE TEMPORARY AND IMPERFECT 
speciesCodes=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/scripts/DNA_Shape_Project/analyses/PhylogeneticDistances/speciesCodes.txt",header=T)
head(speciesCodes)


distances_sp_code1 <- merge(distances_7mer_5mer,speciesCodes,by.x="sp1.x",by.y="code")
# will end up with NAs if comparisons aren't preseent in dataset
distances_sp_code2 <- merge(distances_sp_code1,speciesCodes,by.x="sp2.x",by.y="code",suffixes=c(".sp1",".sp2"))
head(distances_sp_code2)
# alphabetical comparison label
distances_sp_code2$comparisonLabel <- paste0(pmin(distances_sp_code2$species.sp1,distances_sp_code2$species.sp2),".",pmax(distances_sp_code2$species.sp1,distances_sp_code2$species.sp2))

# sweet all are in there: because of alphabetical sorting
distances_sp_code2$comparisonLabel %in% phyloDistances$comparisonLabel

distances_sp_code_phylo <- merge(distances_sp_code2,phyloDistances,by="comparisonLabel")

head(distances_sp_code_phylo)

ggplot(distances_sp_code_phylo,aes(x=cophenetic_distance,y=spectrum_euclidean_distance_7mer))+
  geom_point(aes(color="7mer"))+
  geom_point(aes(x=cophenetic_distance,y=spectrum_euclidean_distance_5mer,color="5mer"))+
  geom_text_repel(aes(label=comparisonLabel))+
  geom_text_repel(aes(label=comparisonLabel,x=cophenetic_distance,y=spectrum_euclidean_distance_5mer))
########## going to add all the comparisons 

############ make a function or something to make this faster -- this block code is too clunky!!!! ############## 
######## need to add 3mer and 1 mer in intelligent ways -- stop with the clunky code! 