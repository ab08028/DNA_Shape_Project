########## Explore some prelim figures/analyses for the dna shape paper ##############
require(ggplot2)
require(reshape2)
require(dplyr)
require(GGally)
require(tidyr)
require(ggrepel)
require(corrr)

#require(corrplot) $ another plotting option in addtion to ggally, but doesn't use ggplot

plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/CorreloGramsAcrossSpeciesSpectra/"
dir.create(paste0(plotdir,"/perCentral3merCorrelations/"),showWarnings = F,recursive = T)
############# 7 mer mutatiton spectra ##################
# going to start with vaquita, mouse, bear, condor; these have targets in low complexity regions and in repeats (nothing masked from genome)
vaquita <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210225_7mer/vaquita.7merspectrum.SummedOverAllIntervals.txt",header=T,sep="\t") # these are whole genome no mask


mouse_allSpp <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/summaryTables/Mouse7merSpectrum.SummedUpOverAllChrs.NOSTRICT.txt",header=T,sep="\t") # NOTE multiple mouse species!!
# select just one mouse species:
mouse_Mmd <- mouse_allSpp[mouse_allSpp$population=="Mmd",]
mouse_Mmd$species <- "mouse_Mmd"

bears_all <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210331_7mer_NoStrict/mapped_to_brown_bear/summaryTables/mapped_to_brown_bear.AllPops.bears.7merspectrum.SummedOverAllIntervals.ALLFREQS.NOSTRICT.txt",header=T,sep="\t") # all bears 
# restrict to one population:
bears_EUR <- bears_all[bears_all$population=="EUR",]
bears_EUR$species <- "brown_bear_EUR"


######### STARTING with unfiltered CONDOR because of the low complexity issues
# only once
#condor_wide <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/birds/mutyperResults_20210419/UNFILTEREDDATA_mutyper_spectrum_files/CA_condor_CRW1112.vcf.gz.mutyper.spectrum.NORANDOMIZE.hetsonly.NOSTRICT.UNFILTEREDDATA.txt",header=T) # note just one individual
#head(condor_wide)
#condor_melt <- melt(condor_wide)
#condor_melt$fractionOfAllSegSites <- condor_melt$value/sum(condor_melt$value)
#head(condor_melt)
#condor_melt$species <- "condor_unfiltered"
#write.table(condor_melt,"/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/birds/mutyperResults_20210419/UNFILTEREDDATA_mutyper_spectrum_files/CA_condor_CRW1112.vcf.gz.mutyper.spectrum.NORANDOMIZE.hetsonly.NOSTRICT.UNFILTEREDDATA.LONG-FORMAT.usethis.txt",quote=F,sep="\t",row.names=F)
condor_melt <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/birds/mutyperResults_20210419/UNFILTEREDDATA_mutyper_spectrum_files/CA_condor_CRW1112.vcf.gz.mutyper.spectrum.NORANDOMIZE.hetsonly.NOSTRICT.UNFILTEREDDATA.LONG-FORMAT.usethis.txt",header=T,sep="\t")

########## adding other species targets #############


########################### combine ####################
commonColumns=c("species","variable","fractionOfAllSegSites")
all <- bind_rows(vaquita[,commonColumns],mouse_Mmd[,commonColumns],bears_EUR[,commonColumns],condor_melt[,commonColumns])

# spread them:
head(all)
all_proportions_spread <- data.frame(spread(all[,c("variable","species","fractionOfAllSegSites")],key=species,value=fractionOfAllSegSites))
head(all_proportions_spread) # some spp don't have certain 7mers! need to make it 0 not NA!

# turn NAs to 0s
all_proportions_spread[is.na(all_proportions_spread)] <- 0
########################## want to get regression between them ########################
species=c("vaquita","mouse_Mmd","brown_bear_EUR","condor_unfiltered")



################ get correlation matrix ################
# some nice ways of working with corr mats: https://www.datanovia.com/en/blog/easy-correlation-matrix-analysis-in-r-using-corrr-package/
pearson_correlationMatrix <- cor(subset(all_proportions_spread,select=-c(variable,central3mer)),method = "pearson")
pearson_correlationMatrix
write.table(pearson_correlationMatrix,"/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/CorreloGramsAcrossSpeciesSpectra/correlationMatrix.7merSpectra.txt",sep="",row.names = F,quote=F)

PlotCorrelationWithoutLabels <- ggpairs(all_proportions_spread,columns=c(2,3,4,5))+ggtitle("Correlogram of 7mer mutation spectra \n(lower = scatter; diag = dist of 7mers per spp; upper=pearson corr.)")+
  theme_bw()+
  theme(text=element_text(size=14))
PlotCorrelationWithoutLabels
ggsave(paste0(plotdir,"7mer.MutationSpectra.CorrelogramPlot.SubsetOfSpecies.WithoutLabels.png"),PlotCorrelationWithoutLabels,height=7,width=9)
# do columns= to get rid of the "motif" label from the main plot (hits cardinality threshold -- but want it as part of geom_label_repel data -- this works! otherwise is insanely slow.)
dataToLabel <- subset(all_proportions_spread,(all_proportions_spread$vaquita>= thresholdForLabel | all_proportions_spread$condor_unfiltered >= thresholdForLabel  | all_proportions_spread$brown_bear_EUR >= thresholdForLabel  | all_proportions_spread$mouse_Mmd>=thresholdForLabel))


PlotCorrelationWithLabels <- ggpairs(all_proportions_spread,columns=c(2,3,4,5),lower = list(continuous = function(data, mapping, ...) ggally_points(data, mapping, ...) + geom_text(data=dataToLabel,aes(label = variable))))+ggtitle("Correlogram of 7mer mutation spectra\n(lower = scatter; diag = dist of 7mers per spp; upper=pearson corr.)")+
  theme_bw()+
  theme(text=element_text(size=14))

PlotCorrelationWithLabels
#ok this labelling isn't ideal -- labels them all on all facets. can'nt fix it!!! 
ggsave(paste0(plotdir,"7mer.MutationSpectra.CorrelogramPlot.SubsetOfSpecies.WithLabels.labelsNotGreat.png"),PlotCorrelationWithLabels,height=7,width=9)


ggcorr(all_proportions_spread[,species])


############ try to get per-3mer correlation matrix #############

all_proportions_spread$central3mer <- paste0(substr(all_proportions_spread$variable,start=3,stop=5),".",substr(all_proportions_spread$variable,start=11,stop=13))
#  okay time to do it the 'dumb' way; didn't wrok with dplyr
central3mers <- unique(all_proportions_spread$central3mer)
# make plots for every 3mer (SLOW!)
for(central3mer in central3mers){
  # get per 3mer data
  data=all_proportions_spread[all_proportions_spread$central3mer==central3mer,]
  # this part is slow only do once if you can -->
  plot <- ggpairs(data,columns=c(2,3,4,5))+ggtitle(central3mer)+theme_bw()
  ggsave(paste0(plotdir,"/perCentral3merCorrelations/",central3mer,".Correlations.png"),plot,height=5,width=8)
}

# get correlation matrix for each one  using nice formatting in the corrr package
allCorMats <- data.frame()
for(central3mer in central3mers){
  # get per 3mer data
  data=all_proportions_spread[all_proportions_spread$central3mer==central3mer,]
  cormatrix <- correlate(subset(data,select=-c(central3mer,variable))) # use correlate in stead of cor to use corrr tools
  # from corr: to make correlation matrices nice (removing duplicate values etc)
  # NICE! very nice looking matrix.
  niceCor <- cormatrix %>%
    shave()
  niceCor$central3mer <- central3mer
  allCorMats <- rbind(allCorMats,niceCor)
  
}
# can use fashion() for nice printing, but I want the NAs
allCorMats
View(allCorMats) 
allCorMats_melt <- melt(allCorMats,id.vars = c("term","central3mer"))
head(allCorMats_melt)
# get rid of Nas (they are duplicates
allCorMats_melt_noNA <- na.omit(allCorMats_melt)
allCorMats_melt_noNA$comparison <- paste0(allCorMats_melt_noNA$term,".",allCorMats_melt_noNA$variable)

corrCoeffsPlot <- ggplot(allCorMats_melt_noNA,aes(y=central3mer,x=value,color=comparison))+
  geom_point()+
  ggtitle("Comparison of Correlation Coefficients between species in 7mer frequencies as fraction of segregating sites, analyzed per central 3mer, between species\nNOTE This is currently the corr in fraction of ALL segregating sites, \nnot just fraction of sites with central3mer -- not sure if that's good or not.")
corrCoeffsPlot
ggsave(paste0(plotdir,"CorrelationCoeffs.7mersWithinCentral3mers.FractionofALLSegSites.NotSureIfThatMakesSenseStatisticallyORNOT.png"),corrCoeffsPlot,height=16,width=7)
############ try to get per-3mer covariance matrices #############


# get variance-covariance matrix for each 3mer 
AllCovarianceMats <- data.frame()
for(central3mer in central3mers){
  # get per 3mer data
  data=all_proportions_spread[all_proportions_spread$central3mer==central3mer,]
  covariancematrix <- colpair_map(subset(data,select=-c(central3mer,variable)), stats::cov) # use this approach instead of cov so that you can tidy it with corrr package https://www.tidyverse.org/blog/2020/12/corrr-0-4-3/
  # from corr: to make correlation matrices nice (removing duplicate values etc)
  # NICE! very nice looking matrix.
  niceCovarianceMatrix <- covariancematrix %>%
    shave()
  niceCovarianceMatrix$central3mer <- central3mer
  AllCovarianceMats <- rbind(AllCovarianceMats,niceCovarianceMatrix)
  
}
# plot covariance 

head(AllCovarianceMats) 
AllCovarianceMats_melt <- melt(AllCovarianceMats,id.vars = c("term","central3mer"))
head(allCorMats_melt)
# get rid of Nas (they are duplicates
AllCovarianceMats_melt_noNA <- na.omit(AllCovarianceMats_melt)
AllCovarianceMats_melt_noNA$comparison <- paste0(AllCovarianceMats_melt_noNA$term,".",AllCovarianceMats_melt_noNA$variable)
covarPlot <- ggplot(AllCovarianceMats_melt_noNA,aes(y=central3mer,x=value,color=comparison))+
  geom_point()+
  ggtitle("Comparison of Covariance between species in 7mer frequencies as fraction of segregating sites, analyzed per central 3mer, between species\nNOTE This is currently the covar in fraction of all segregating sites, \nnot just fraction of sites with central3mer -- not sure if that's good or not.")
ggsave(paste0(plotdir,"CovarianceBetweenSpecies.7mersWithinCentral3mers.FractionofALLSegSites.NotSureIfThatMakesSenseStatisticallyORNOT.png"),covarPlot,height=16,width=7)
############### HANG ON!!!!! Do I want to do this for covariance of 7mers with
