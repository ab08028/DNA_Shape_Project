############ comparing vaquita and fin whale 7mer context ##############

########## NOTE THIS SCRIPT IS VERY SLOW AND DOES A LOT OF PROCESSING.
#If you want to replot anything, use the output files instead which are here:
  
#fw = read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/fin_whales/results/mutyper/downloaddate_20210111/mutyperResults_20210225_7mer/finwhale.7merspectrum.SummedOverAllIntervals.txt",header=T,sep="\t")
#vaq = read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210225_7mer/vaquita.7merspectrum.SummedOverAllIntervals.txt",header = T,sep="\t")
#mouse = read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyer_results_20210226_7mer/mouse.7merspectrum.SummedOverAllChrs.txt",header = T,sep="\t")


require(ggplot2)
require(reshape2)
require(ggrepel)
plotdir="/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210225_7mer/"
require(dplyr)
vaquita7mer <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210225_7mer/mutyper_spectrum_files/mutyper.spectra.PERPOPULATION.ALLFREQS.7mer.txt",header=T)
vaquita7mer_melt <- melt(vaquita7mer)
vaquita7mer_melt$species <- "vaquita"
dim(vaquita7mer) # 22423 7mers
# need to condense 
# want to condense down to
vaquita7mer_melt$centralMutationType <- NA
# "AAAAAAA.AAACAAA" so the position of the central mutations are 4 and 12 (remember ther's a dot in between')
# and 3-5 is teh 3mer context and 11-13
vaquita7mer_melt$centralMutationType <- paste(substr(vaquita7mer_melt$variable,4,4),substr(vaquita7mer_melt$variable,12,12),sep=".")
vaquita7mer_melt$ThreemerMutationType <- NA
vaquita7mer_melt$ThreemerMutationType <- paste(substr(vaquita7mer_melt$variable,3,5),substr(vaquita7mer_melt$variable,11,13),sep=".")

vaquita7mer_melt$labelOfInterest <- "Not the motif"
### want to label any that contain the 6mer of interest: NTTAAAA>NTTTAAA
vaquita7mer_melt[grep("*TTAAAA.*TTTAAA",vaquita7mer_melt$variable),]$labelOfInterest <- "LINE1 receptor motif NTTAAAA.NTTTAAA"
# 6mer of interest from Jed's paper (LINE1 binding domain) NTTAAAA>NTTTAAA (how does sequence cntext matter?)
# I am also seeing this one show up a lot: AAAATTT.AAATTTT
vaquita7mer_melt[vaquita7mer_melt$variable=="AAAATTT.AAATTTT",]$labelOfInterest <- "AAAATTT.AAATTTT possible reverse of motif?"

### get fraction of A>T mutations 
vaquita7mer_melt <- vaquita7mer_melt %>%
  group_by(centralMutationType) %>%
  mutate(fractionOfCentralMutationType=value/sum(value))

### get fraction of all vaquita7mer_melt sites 
vaquita7mer_melt$fractionOfAllSegSites <- vaquita7mer_melt$value/sum(vaquita7mer_melt$value) 



# plot1 <- ggplot(vaquita7mer_melt[vaquita7mer_melt$centralMutationType=="A>T",],aes(x=variable,y=value,color=labelOfInterest))+
#   geom_point()+
#   coord_flip()+
#   ggtitle("A>T mutation types (7mer context)")+
#   xlab("SNP count (genome-wide")
# plot1
# ggsave(paste(plotdir,"vaquita.A.T.7MerContext.pdf",sep=""),plot1,height=12,width=18)
# also want to plot fraction of A>T mutation types?

plot2 <- ggplot(vaquita7mer_melt,aes(y=centralMutationType,x=fractionOfAllSegSites,color=labelOfInterest))+
  geom_point()+
  theme_bw()+
  ggtitle("Vaquita\nComparison of different 7mers. Each dot is a different 7mer")+
  scale_color_manual(values=c("tomato","orange","gray"))

plot2
ggsave(paste(plotdir,"vaquita.7mersAsFractionOfSegSites.pdf",sep=""),plot2,height=12,width=18)







############## fin whale ################### GOC only for now
# need to gather up spectra an populations 

genotypeDate="downloaddate_20210111"
mutyperDate="20210225_7mer" # date you ran mutyperl 20210118 is the first run of mutyper after the bug fix! 
intervalCount=96
intervals <- sprintf("%02d", 1:intervalCount) # need to have leading 0 for things <10.

populations=c("GOC","ENP")

wd=paste("/Users/annabelbeichman/Documents/UW/WhaleProject/fin_whales/results/mutyper/",genotypeDate,"/mutyperResults_",mutyperDate,"/",sep="")


allIntervalSpectra=data.frame()
for(interval in intervals){
  spectra=data.frame()
  for(pop in populations){
    input=read.table(paste(wd,"/mutyper_spectrum_files/interval_",interval,"/",pop,"/",pop,"_samples_interval_",interval,".mutyper.spectra.PERPOPULATION.ALLFREQS.txt",sep=""),header=T,stringsAsFactors = F)
    input$population <- pop
    spectra <- bind_rows(spectra, input) # in the only singletons category some are missing, so having NA put in using dplyr:bind_rows
    # ksfs:
  }
  spectra_melt <- melt(spectra)
  spectra_melt$interval <- interval
  allIntervalSpectra <- bind_rows(allIntervalSpectra,data.frame(spectra_melt))

}


# when mutyper finds no sites it doesn't output column, instead of just outputting 0.
# so when rbinding there are NAs filled in. Need to replace these with 0s before summing up otherwise introduces lots of Nas into sums

allIntervalSpectra[is.na(allIntervalSpectra)] <- 0


finwhale_allIntervalSpectra_totals_SummedUp <- data.frame(allIntervalSpectra %>% 
  group_by(population,variable) %>%
  summarise(totalSitesAllIntervals=sum(value)))


finwhale_allIntervalSpectra_totals_SummedUp$centralMutationType <- NA
# "AAAAAAA.AAACAAA" so the position of the central mutations are 4 and 12 (remember ther's a dot in between')
# and 3-5 is teh 3mer context and 11-13
finwhale_allIntervalSpectra_totals_SummedUp$centralMutationType <- paste(substr(finwhale_allIntervalSpectra_totals_SummedUp$variable,4,4),substr(finwhale_allIntervalSpectra_totals_SummedUp$variable,12,12),sep=".")
finwhale_allIntervalSpectra_totals_SummedUp$ThreemerMutationType <- NA
finwhale_allIntervalSpectra_totals_SummedUp$ThreemerMutationType <- paste(substr(finwhale_allIntervalSpectra_totals_SummedUp$variable,3,5),substr(finwhale_allIntervalSpectra_totals_SummedUp$variable,11,13),sep=".")

finwhale_allIntervalSpectra_totals_SummedUp$labelOfInterest <- "Not the motif"
### want to label any that contain the 6mer of interest: NTTAAAA>NTTTAAA
finwhale_allIntervalSpectra_totals_SummedUp[grep("*TTAAAA.*TTTAAA",finwhale_allIntervalSpectra_totals_SummedUp$variable),]$labelOfInterest <- "LINE1 receptor motif NTTAAAA.NTTTAAA"
# 6mer of interest from Jed's paper (LINE1 binding domain) NTTAAAA>NTTTAAA (how does sequence cntext matter?)
finwhale_allIntervalSpectra_totals_SummedUp[finwhale_allIntervalSpectra_totals_SummedUp$variable=="AAAATTT.AAATTTT",]$labelOfInterest <- "AAAATTT.AAATTTT possible reverse of motif?"


### get fraction of 1mer mutation type ### NEED TO GROUP BY POPULATION OTHERWISE IS WRONG!! 
finwhale_allIntervalSpectra_totals_SummedUp <- finwhale_allIntervalSpectra_totals_SummedUp %>%
  group_by(centralMutationType,population) %>%
  mutate(fractionOfCentralMutationType=totalSitesAllIntervals/sum(totalSitesAllIntervals))

### get fraction of all segregating sites 
finwhale_allIntervalSpectra_totals_SummedUp <- finwhale_allIntervalSpectra_totals_SummedUp %>%
  group_by(population) %>%
  mutate(fractionOfAllSegSites=totalSitesAllIntervals/sum(totalSitesAllIntervals))

######### plot comparing fw and vaquita #########

vaquita7mer_melt$population <- "vaquita"
plotfv1 <- ggplot(finwhale_allIntervalSpectra_totals_SummedUp,aes(y=centralMutationType,x=fractionOfAllSegSites,color=labelOfInterest))+
  geom_point()+
  geom_point(data=vaquita7mer_melt,aes(y=centralMutationType,x=fractionOfAllSegSites,color=labelOfInterest))+
  theme_bw()+
  facet_wrap(~population,scales="free_x")+
  ggtitle("Comparison of different 7mers. Each dot is a different 7mer")+
  scale_color_manual(values=c("tomato","orange","gray"))

plotfv1
ggsave(paste(plotdir,"FINwhale.Vaquita.Comparison.7mers.pdf",sep=""),plotfv1,height=7,width=12)

#################### plot GOC FW and vaquita on same plot along diagonal ###############

# need to merge based on variable and just one FW population
GOCfw_plus_vaq <- merge(vaquita7mer_melt[,c("variable", "ThreemerMutationType","centralMutationType" ,"labelOfInterest","fractionOfAllSegSites")],finwhale_allIntervalSpectra_totals_SummedUp[finwhale_allIntervalSpectra_totals_SummedUp$population=="GOC",c("variable","labelOfInterest", "ThreemerMutationType","centralMutationType" ,"fractionOfAllSegSites")],by=c("variable","labelOfInterest","centralMutationType","ThreemerMutationType"),suffixes = c(".vaq",".fw"))

ENPfw_plus_vaq <- merge(vaquita7mer_melt[,c("variable", "ThreemerMutationType","centralMutationType" ,"labelOfInterest","fractionOfAllSegSites")],finwhale_allIntervalSpectra_totals_SummedUp[finwhale_allIntervalSpectra_totals_SummedUp$population=="ENP",c("variable","labelOfInterest", "ThreemerMutationType","centralMutationType" ,"fractionOfAllSegSites")],by=c("variable","labelOfInterest","centralMutationType","ThreemerMutationType"),suffixes = c(".vaq",".fw"))

GOC_vaquita_diagonalPlot <- ggplot(GOCfw_plus_vaq,aes(x=fractionOfAllSegSites.fw,y=fractionOfAllSegSites.vaq,color=centralMutationType,label=variable))+
  geom_point(size=1.2,alpha=0.7)+
  geom_abline(intercept=0,slope=1)+
  theme_bw() +
  geom_label_repel(data=subset(GOCfw_plus_vaq,fractionOfAllSegSites.vaq>=0.00075), nudge_y=0.003,size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",direction="x",show.legend = F)+
  xlab("Fin Whale 7mer proportion of all segregating sites (Gulf of California)")+
  ylab("Vaquita 7mer proportion of all segregating sites")+
  ggtitle("Comparing fin whale and vaquita 7mer proportions")+
  theme(text=element_text(size=14))+
  guides(colour = guide_legend(override.aes = list(size=3))) # this makes the dots in the legend larger
GOC_vaquita_diagonalPlot
ggsave(paste(plotdir,"Diagonal.GOCfinwhale.vs.Vaquita.Comparison.7mers.pdf",sep=""),GOC_vaquita_diagonalPlot,height=7,width=12)

ENP_vaquita_diagonalPlot <- ggplot(ENPfw_plus_vaq,aes(x=fractionOfAllSegSites.fw,y=fractionOfAllSegSites.vaq,color=centralMutationType,label=variable))+
  geom_point(size=1.2,alpha=0.7)+
  geom_abline(intercept=0,slope=1)+
  theme_bw() +
  geom_label_repel(data=subset(ENPfw_plus_vaq,fractionOfAllSegSites.vaq>=0.00075), nudge_y=0.003,size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",direction="x",show.legend = F)+
  xlab("Fin Whale 7mer proportion of all segregating sites (ENP)")+
  ylab("Vaquita 7mer proportion of all segregating sites")+
  ggtitle("Comparing fin whale and vaquita 7mer proportions")+
  theme(text=element_text(size=14))+
  guides(colour = guide_legend(override.aes = list(size=3))) # this makes the dots in the legend larger
ENP_vaquita_diagonalPlot
ggsave(paste(plotdir,"Diagonal.ENPfinwhale.vs.Vaquita.Comparison.7mers.pdf",sep=""),ENP_vaquita_diagonalPlot,height=7,width=12)



############### add in mice next ##############
mpopulations=c("Mmd", "Mmc" ,"Mmm", "Ms")
mchromosomecount=19
mallIntervalSpectra=data.frame()
### can skip this and just read in table since I did this once

for(i in seq(1,mchromosomecount)){
  mspectra=data.frame()
  mchromosome=paste("chr",i,sep="")
  for(mpop in mpopulations){
    minput=read.table(paste("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyer_results_20210226_7mer/mutyper_spectrum_files/",mpop,"/",mchromosome,"_",mpop,"_samples.mutyper.spectra.PERPOPULATION.ALLFREQS.txt",sep=""),header=T,stringsAsFactors = F)
    minput$population <- mpop
    mspectra <- bind_rows(mspectra, minput) # in the only singletons category some are missing, so having NA put in using dplyr:bind_rows
    # ksfs:
  }
  mspectra_melt <- melt(mspectra)
  mspectra_melt$interval <- mchromosome
  mallIntervalSpectra <- bind_rows(mallIntervalSpectra,data.frame(mspectra_melt))
  
}

mallIntervalSpectra[is.na(mallIntervalSpectra)] <- 0


mouse_allIntervalSpectra_totals_SummedUp <- data.frame(mallIntervalSpectra %>% 
                                                            group_by(population,variable) %>%
                                                            summarise(totalSitesAllIntervals=sum(value)))


mouse_allIntervalSpectra_totals_SummedUp$centralMutationType <- NA
# "AAAAAAA.AAACAAA" so the position of the central mutations are 4 and 12 (remember ther's a dot in between')
# and 3-5 is teh 3mer context and 11-13
mouse_allIntervalSpectra_totals_SummedUp$centralMutationType <- paste(substr(mouse_allIntervalSpectra_totals_SummedUp$variable,4,4),substr(mouse_allIntervalSpectra_totals_SummedUp$variable,12,12),sep=".")
mouse_allIntervalSpectra_totals_SummedUp$ThreemerMutationType <- NA
mouse_allIntervalSpectra_totals_SummedUp$ThreemerMutationType <- paste(substr(mouse_allIntervalSpectra_totals_SummedUp$variable,3,5),substr(mouse_allIntervalSpectra_totals_SummedUp$variable,11,13),sep=".")

mouse_allIntervalSpectra_totals_SummedUp$labelOfInterest <- "Not the motif"
### want to label any that contain the 6mer of interest: NTTAAAA>NTTTAAA
mouse_allIntervalSpectra_totals_SummedUp[grep("*TTAAAA.*TTTAAA",mouse_allIntervalSpectra_totals_SummedUp$variable),]$labelOfInterest <- "LINE1 receptor motif NTTAAAA.NTTTAAA"
# 6mer of interest from Jed's paper (LINE1 binding domain) NTTAAAA>NTTTAAA (how does sequence cntext matter?)
mouse_allIntervalSpectra_totals_SummedUp[mouse_allIntervalSpectra_totals_SummedUp$variable=="AAAATTT.AAATTTT",]$labelOfInterest <- "AAAATTT.AAATTTT possible reverse of motif?"


### get fraction of 1mer mutation type ### NEED TO GROUP BY POPULATION OTHERWISE IS WRONG!! 
mouse_allIntervalSpectra_totals_SummedUp <- mouse_allIntervalSpectra_totals_SummedUp %>%
  group_by(centralMutationType,population) %>%
  mutate(fractionOfCentralMutationType=totalSitesAllIntervals/sum(totalSitesAllIntervals))

### get fraction of all segregating sites 
mouse_allIntervalSpectra_totals_SummedUp <- mouse_allIntervalSpectra_totals_SummedUp %>%
  group_by(population) %>%
  mutate(fractionOfAllSegSites=totalSitesAllIntervals/sum(totalSitesAllIntervals))


######### plot comparing fw and vaquita and mouse #########

plotfvm2 <- ggplot(finwhale_allIntervalSpectra_totals_SummedUp,aes(y=centralMutationType,x=fractionOfAllSegSites,color=labelOfInterest))+
  geom_point()+
  geom_point(data=vaquita7mer_melt,aes(y=centralMutationType,x=fractionOfAllSegSites,color=labelOfInterest))+
  geom_point(data=mouse_allIntervalSpectra_totals_SummedUp,aes(y=centralMutationType,x=fractionOfAllSegSites,color=labelOfInterest))+
  theme_bw()+
  facet_wrap(~population)+
  ggtitle("Comparison of different 7mers. Each dot is a different 7mer")+
  scale_color_manual(values=c("tomato","orange","gray"))

plotfvm2
ggsave(paste(plotdir,"Mouse.Finwhale.Vaquita.Comparison.7mers.pdf",sep=""),plotfvm2,height=7,width=12)

########### write out tables: ###########
write.table(mouse_allIntervalSpectra_totals_SummedUp,paste("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyer_results_20210226_7mer/mouse.7merspectrum.SummedOverAllChrs.txt",sep=""),row.names=F,quote=F,sep="\t")

fgenotypeDate="downloaddate_20210111"
fmutyperDate="20210225_7mer" # date you ran mutyperl 20210118 is the first run of mutyper after the bug fix! 
fwd=paste("/Users/annabelbeichman/Documents/UW/WhaleProject/fin_whales/results/mutyper/",genotypeDate,"/mutyperResults_",mutyperDate,"/",sep="")
write.table(finwhale_allIntervalSpectra_totals_SummedUp,paste(fwd,"/finwhale.7merspectrum.SummedOverAllIntervals.txt",sep=""),row.names=F,quote=F,sep="\t")

write.table(vaquita7mer_melt,paste(plotdir,"/vaquita.7merspectrum.SummedOverAllIntervals.txt",sep=""),row.names=F,quote=F,sep="\t")
