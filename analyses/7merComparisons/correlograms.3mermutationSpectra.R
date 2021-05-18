########## Explore some prelim figures/analyses for the dna shape paper ##############
require(ggplot2)
require(reshape2)
require(dplyr)
require(GGally)
require(tidyr)
require(ggrepel)

plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/CorreloGramsAcrossSpeciesSpectra/"
############# 3mer mutation spectra ##################
# going to start with vaquita, mouse, bear, condor; these have targets in low complexity regions and in repeats (nothing masked from genome)
vaquita <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210216/summaryTables/All_Intervals_CountsProps.ALLFREQS.txt",header=T,sep="\t") # these are whole genome no mask
head(vaquita)
vaquita$species <- "vaquita"
######### don't have 3mer for mouse yet! (only 1mer and 7mer!!)
#mouse_allSpp <- read.table("",header=T,sep="\t") # NOTE multiple mouse species!!
# select just one mouse species:
#mouse_Mmd <- mouse_allSpp[mouse_allSpp$population=="Mmd",]
#mouse_Mmd$species <- "mouse_Mmd"

bears_all <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210401_3mer_nostrict/mapped_to_brown_bear/summaryTables/All_Intervals_CountsProps.ref.brown_bear.ALLFREQS.nostrict.txt",header=T,sep="\t") # all bears 
head(bears_all)
# restrict to one population:
bears_EUR <- bears_all[bears_all$population=="EUR",]
bears_EUR$species <- "brown_bear_EUR"


######### don't have 3mer for condor yet!!!
# only once
#condor_wide <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/birds/mutyperResults_20210419/UNFILTEREDDATA_mutyper_spectrum_files/CA_condor_CRW1112.vcf.gz.mutyper.spectrum.NORANDOMIZE.hetsonly.NOSTRICT.UNFILTEREDDATA.txt",header=T) # note just one individual
#head(condor_wide)
#condor_melt <- melt(condor_wide)
#condor_melt$fractionOfAllSegSites <- condor_melt$value/sum(condor_melt$value)
#head(condor_melt)
#condor_melt$species <- "condor_unfiltered"
#write.table(condor_melt,"/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/birds/mutyperResults_20210419/UNFILTEREDDATA_mutyper_spectrum_files/CA_condor_CRW1112.vcf.gz.mutyper.spectrum.NORANDOMIZE.hetsonly.NOSTRICT.UNFILTEREDDATA.LONG-FORMAT.usethis.txt",quote=F,sep="\t",row.names=F)
#condor_melt <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/birds/mutyperResults_20210419/UNFILTEREDDATA_mutyper_spectrum_files/CA_condor_CRW1112.vcf.gz.mutyper.spectrum.NORANDOMIZE.hetsonly.NOSTRICT.UNFILTEREDDATA.LONG-FORMAT.usethis.txt",header=T,sep="\t")

########## adding other species targets #############


########################### combine ####################
commonColumns=c("species","variable","fractionOfAllSegSites")
#all <- bind_rows(vaquita[,commonColumns],mouse_Mmd[,commonColumns],bears_EUR[,commonColumns],condor_melt[,commonColumns])
all <- bind_rows(bears_EUR,vaquita)
# spread them:
head(all)
all_proportions_spread <- data.frame(spread(all[,c("variable","species","proportionOfSites")],key=species,value=proportionOfSites))
head(all_proportions_spread) # some spp don't have certain 7mers! need to make it 0 not NA!

# turn NAs to 0s
all_proportions_spread[is.na(all_proportions_spread)] <- 0
########################## want to get regression between them ########################
#species=c("vaquita","mouse_Mmd","brown_bear_EUR","condor_unfiltered")
species=c("vaquita","brown_bear_EUR")
pairs=data.frame(combn(species,2))
dim(pairs)
numberOfPairs=dim(pairs)[2]
numberOfPairs


PlotCorrelationWithoutLabels <- ggpairs(all_proportions_spread,columns=c(2,3))+ggtitle("Correlogram of 3mer mutation spectra (lower = scatter; diag = dist of 7mers per spp; upper=pearson corr.")+
  theme_bw()+
  theme(text=element_text(size=14))
PlotCorrelationWithoutLabels

ggsave(paste0(plotdir,"3mer.MutationSpectra.CorrelogramPlot.SubsetOfSpecies.WithoutLabels.png"),PlotCorrelationWithoutLabels,height=4,width=6)




ggcorr(all_proportions_spread[,species])
