################ want to plot de novo spectrum from S.P.E with 3mers between pairs of species and see if it captures species specific diffs in spectrum ########

# rescale in some fashion as well (figure out next)
require(tidyverse)
require(dplyr)
require(reshape2)
wd="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/SigProfilerExtractor/20210805_SigProfilerExtractorResults_bearsMiceVaquitaFinWhaleBulkCetaceans__GRCh37GenomeRef/SBS96/"
plotdir=paste0(wd,"myplots/")
dir.create(plotdir,showWarnings = F)
# empirical spectra:
empirical=read_tsv(paste0(wd,"Samples.txt"),col_names=T)
head(empirical)
empirical_melt <- melt(empirical)
head(empirical_melt)
colnames(empirical_melt) <- c("MutationType","Samples","totalEmpiricalMutations")

# how many sites are contributed by each de novo signature (sticking to de novo because cosmic compositions don't fit data well already)
denovo_activities = read_tsv(paste0(wd,"Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt"),col_names=T)
head(denovo_activities)

# makeup of each signature (each column sums to 1)
denovo_signatures = read_tsv(paste0(wd,"Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt"),col_names=T)
head(denovo_signatures)

# melt them:
denovo_signatures_melt <- melt(denovo_signatures)
colnames(denovo_signatures_melt) <- c("MutationType","signature","fractionOfSignature")
head(denovo_signatures_melt)

denovo_activities_melt <- melt(denovo_activities)
colnames(denovo_activities_melt) <- c("Samples","signature","totalSitesOfAllTypesContibutedBySignatureToSample")
head(denovo_activities_melt)

denovo_combined <- merge(denovo_signatures_melt,denovo_activities_melt)
head(denovo_combined)

denovo_combined$totalSites_specificType_ContributedBySignature <- denovo_combined$fractionOfSignature*denovo_combined$totalSitesOfAllTypesContibutedBySignatureToSample

head(denovo_combined)
tail(denovo_combined)
# okay then want to get 3 mer spectrum for each species that is predicted

# want to sum up total mutations of each Mutation Type that are contributed per species to get the predicted spectrum:
predicted_spectra <- denovo_combined %>%
  group_by(MutationType, Samples) %>%
  summarise(totalPredictedMutations=sum(totalSites_specificType_ContributedBySignature))

head(predicted_spectra)


predicted_plus_empirical <- merge(predicted_spectra,empirical_melt,by=c("MutationType","Samples"))
head(predicted_plus_empirical)
# okay now I have the 3mer spectrum for each species/sample within species


# rescale to fraction of segregating sites per species:
predicted_plus_empirical <- predicted_plus_empirical %>%
  group_by(Samples) %>%
  mutate(fractionSegSitesPredicted=totalPredictedMutations/sum(totalPredictedMutations),fractionSegSitesEmpirical=totalEmpiricalMutations/sum(totalEmpiricalMutations))

# check:
sum(predicted_plus_empirical[predicted_plus_empirical$Samples=="Mmd",]$fractionSegSitesEmpirical)
sum(predicted_plus_empirical[predicted_plus_empirical$Samples=="Mmd",]$fractionSegSitesPredicted)
# so let's plot a pair of species (loop over in some smarter way momentarily)


plot1 <- ggplot(predicted_plus_empirical,aes(x=totalEmpiricalMutations,y=totalPredictedMutations,color=Samples))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  theme_bw()+
  ggtitle("every dot is a 3mer of the spectrum")+
  scale_y_log10()+
  scale_x_log10()
plot1
ggsave(paste0(plotdir,"predicted.empirical.counts.log10.png"),plot1,height=6,width=8)

# plot frac of seg sites:

plot2 <- ggplot(predicted_plus_empirical,aes(x=fractionSegSitesEmpirical,y=fractionSegSitesPredicted,color=Samples))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  theme_bw()+
  ggtitle("every dot is a 3mer of the spectrum")+
  scale_y_log10()+
  scale_x_log10()
plot2
ggsave(paste0(plotdir,"predicted.empirical.counts.FracOfSegSites.png"),plot2,height=6,width=8)


# try pairs of species:
# remove some species that have duplicates
#keySpecies=c("EUR","bacu2","bmus1","dleu1","ENP","gmel1","lobl1","Mmc","Mmd","Mmm","mmon1","Ms","nasi1","npho1","oorc1","PB","pmac1","ppho1","ttru1","Vaquita") #
#leaving out duplicates of species
keySpecies=c("ENP","Ms","Mmd","Mmc","Mmm","PB","EUR","Vaquita","EUR")
allpairs=combn(unique(keySpecies),2)
allpairs
dim(allpairs)

subplotdir=paste0(plotdir,"allSppPairs_predictedEmpirical_3merSpectra/")
dir.create(subplotdir,showWarnings = F)

for(i in seq(1,dim(allpairs)[2])) {
  pair=c(allpairs[,i])
  #pair=c("Mmd","Ms")
  # add a generic label:
  paireddf = predicted_plus_empirical[predicted_plus_empirical$Samples %in% pair,]
  paireddf$label <- ""
  paireddf[paireddf$Samples == pair[1],]$label <- "1"
  paireddf[paireddf$Samples == pair[2],]$label <- "2"
  
  paireddf_wide = pivot_wider(paireddf,id_cols=c(MutationType,label),names_from = label,values_from=c(totalPredictedMutations,totalEmpiricalMutations,fractionSegSitesEmpirical,fractionSegSitesPredicted))
  head(paireddf_wide)
  
  pairplot1 <- ggplot(paireddf_wide,aes(x=fractionSegSitesEmpirical_1,y=fractionSegSitesEmpirical_2,color="empirical",shape="empirical"))+
    geom_point()+
    geom_abline(slope=1,intercept = 0)+
    geom_point(aes(x=fractionSegSitesPredicted_1,y=fractionSegSitesPredicted_2,color="predicted by model",shape="predicted by model"))+
    xlab(pair[1])+
    ylab(pair[2])+
    ggtitle("Comparing empirical (black) and predicted by SPE (red) fraction of segregating sites")+
    theme_bw()+
    #scale_y_log10()+
    #scale_x_log10()+
    scale_color_manual(values=c("black","red"))+
    geom_segment(aes(x=fractionSegSitesEmpirical_1,xend=fractionSegSitesPredicted_1,y=fractionSegSitesEmpirical_2,yend=fractionSegSitesPredicted_2),alpha=0.5,color="darkgrey")
  pairplot1
  

  ggsave(paste0(subplotdir,pair[1],".",pair[2],".png"),pairplot1,height=6,width=8)
  
  pairplot2 <- pairplot1 + scale_x_log10()+scale_y_log10()
  
  ggsave(paste0(subplotdir,pair[1],".",pair[2],".log10.png"),pairplot2,height=6,width=8)
  

  }
# trial of plotly:
# need to add  text=MutationType to the aes()
#require(plotly)
# require(htmlwidgets)
#test <- ggplotly(pairplot1,tooltip="text")
# htmlwidgets::saveWidget(test, file = "test.html")

