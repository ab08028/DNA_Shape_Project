require(ggplot2)
require(ggrepel)
plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/diagonal_comparison_plots/targetsInGenomes/"

## okay there remains an issue: for cetaceans I used a bed mask to get just where JAR's calls where (no repeat regions) whereas for mouse/bear I allowed targets to be called in whole genome , so simple repeat regions were included
# so this makes the mouse have an excess of simple repeat 7mers like ACACACAC that aren't in the cetaceans.
# how do I want to deal with this?
######## diagonal plots of targets 7mer
# get proprtions
#vaquita <- read.table(" haven't yet run vaquita targets iwthout strict ",header=F)
#colnames(vaquita) <- c("target","count.vaq")
#vaquita$proportion.vaq <- vaquita$count.vaq/sum(vaquita$count.vaq)

#fw_minke_wholegenome <- read.table(" haven't yet redone fin whale with out strict",header=F)
#colnames(fw_minke_wholegenome) <- c("target","count.minkeWG")
#fw_minke_wholegenome$proportion.minkeWG<- fw_minke_wholegenome$count.minkeWG/sum(fw_minke_wholegenome$count.minkeWG)

#fw_minke_neutralonly <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/fin_whales/results/mutyper/downloaddate_20210111/mutyperResults_20210225_7mer/mutyper_target_files/neutral_only/mutyper.targets.7mer.NEUTRALBEDCOORDSONLY.NoteAreVeryChoppyCoords.strict.txt",header=F)
#colnames(fw_minke_neutralonly) <- c("target","count.minkeNeutral")
#fw_minke_neutralonly$proportion.minkeNeutral <- fw_minke_neutralonly$count.minkeNeutral/sum(fw_minke_neutralonly$count.minkeNeutral)

# note for the following there *is* a header, but the above there is *not* 
# 20210407 these are new with no strict (good) for brown bear mouse and polar bear (good !)
# note they didnn't use any sort of bedmask of called sites yet though 
brownbear=read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210331_7mer_NoStrict/mapped_to_brown_bear/summaryTables/ALLINTERVALS.brown_bear.mutyper.targets.7mer.nostrict.txt",header=T)
colnames(brownbear) <- c("target","count")
brownbear$proportion <- brownbear$count/sum(brownbear$count)
brownbear$species <- "brown_bear"
brownbear$id <- "brown_bear"
brownbear$reference <- "brown_bear"

polarbear=read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210331_7mer_NoStrict/mapped_to_polar_bear/summaryTables/ALLINTERVALS.polar_bear.mutyper.targets.7mer.nostrict.txt",header=T)
colnames(polarbear) <- c("target","count")
polarbear$proportion <- polarbear$count/sum(polarbear$count)
polarbear$species <- "polar_bear"
polarbear$id <- "polar_bear"
polarbear$reference <- "polar_bear"

mouse=read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/summaryTables/ALLINTERVALS.mouse.mutyper.targets.7mer.nostrict.txt",header=T)
colnames(mouse) <- c("target","count")
mouse$proportion <- mouse$count/sum(mouse$count)
mouse$id <- "mm10_nostrict"
mouse$reference <- "mm10"
mouse$species <- "mouse"
# need to read in cetaceans #

cetaceansampleinfo <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/information/Cetacean.SampleInfo.txt",header = T,sep="\t")
#### note: one individual 
allCetaceanTargets=data.frame() # each spectra has an individual name attached 
# skip minke for now
ids=cetaceansampleinfo$ID
ids=c("GCF_002837175.2_ASM283717v2_pmac2", "GCF_009873245.2_mBalMus1.pri.v3_bmus1", "GCF_011762595.1_mTurTru1.mat.Y_ttru1", "GCF_000493695.1_BalAcu1.0_bacu3", "GCF_000331955.2_Oorc_1.1_oorc1", "GCF_006547405.1_ASM654740v1_gmel1", "GCF_002837175.2_ASM283717v2_pmac1")
for(id in ids){
  file=list.files("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_targets_files_incomplete/",pattern=paste(id,".mutyper.targets.NOSTRICT.7mers.txt",sep=""),full.names = T)
  targets <- read.table(file,header=F,sep="\t")
  colnames(targets) <- c("target","count")
  targets$id <- id 
  targets$species <- cetaceansampleinfo[cetaceansampleinfo$ID==id,]$Species_common
  targets$reference <- cetaceansampleinfo[cetaceansampleinfo$ID==id,]$Reference
  targets$proportion <- targets$count/sum(targets$count)
  allCetaceanTargets = bind_rows(allCetaceanTargets,targets) # use bind-rows to fill in some columns that may not exist in some spectra (if 0 observations)
  
}

  
mergedTargets <- bind_rows(allCetaceanTargets,brownbear,polarbear,mouse)


# want to make it wide for diagonal plots (other ways to plot?)
mergedTargets_wide <- spread(mergedTargets[,c("id","proportion","target")],key=id,value=c(proportion))
head(mergedTargets_wide)


mergedTargets_wide$centralBP <- substr(mergedTargets_wide$target,4,4)
mergedTargets_wide$centralthreemer <- substr(mergedTargets_wide$target,3,5)
head(mergedTargets_wide)
# some diagonal plots
# minke vs bottlenose?
minke_v_bottlenose <- ggplot(mergedTargets_wide,aes(x=GCF_011762595.1_mTurTru1.mat.Y_ttru1,y=GCF_000493695.1_BalAcu1.0_bacu3,label=target))+
  geom_point()+
  geom_abline(slope=1,intercept=0)+
  geom_label_repel(data=subset(mergedTargets_wide,GCF_011762595.1_mTurTru1.mat.Y_ttru1>0.001))+
  theme_bw()+
  ggtitle("note that bed mask was used to avoid repeat masker regions/uncalled sites")
minke_v_bottlenose
ggsave(paste(plotdir,"targets.minke_v_bottlenose.nostrict.png",sep = ""),minke_v_bottlenose,height=5,width=7)

brown_v_polar <-  ggplot(mergedTargets_wide,aes(x=brown_bear,y=polar_bear,label=target))+
  geom_point()+
  geom_abline(slope=1,intercept=0)+
  geom_label_repel(data=subset(mergedTargets_wide,abs(polar_bear-brown_bear)>0.00005 | brown_bear>0.0015))+
  theme_bw()
brown_v_polar
ggsave(paste(plotdir,"targets.brown_v_polar.nostrict.png",sep = ""),brown_v_polar,height=5,width=7)

# minke vs mouse 
minke_v_mouse <- ggplot(mergedTargets_wide,aes(x=mm10_nostrict,y=GCF_000493695.1_BalAcu1.0_bacu3,label=target))+
  geom_point()+
  geom_abline(slope=1,intercept=0)+
  geom_label_repel(data=subset(mergedTargets_wide,mm10_nostrict>0.001 | GCF_000493695.1_BalAcu1.0_bacu3>0.001))+
  theme_bw()+
  ggtitle("note that bed mask was used in MINKE (not mouse) to avoid repeat masker regions/uncalled sites")

minke_v_mouse
ggsave(paste(plotdir,"targets.minke_v_mouse.IssueThatMouseHasNoBedMask.soLotsOfSimpleRepeats.png",sep = ""),minke_v_mouse,height=5,width=7)


######### experimental target correction ########
# (count in minke / minke target proportion ) * mouse target proportion
# doing this manually just to think about it first 

cetaceanMutationProportions <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_spectrum_files/noSTRICTflag/AllCetaceans_melted.IncludesTwoConditions.HetsOnly.and.HetsPlusHoms.NOSTRICTFLAG.txt",header=T,sep="\t") # NOTE THIS HAD TWO CONDITIONS (hets only and hets plus homs)

# note this does have the simple repeats issue still!!! 
# want to restrict this just to one condition and one species
speciesID="GCF_000493695.1_BalAcu1.0_bacu3"
condition="hetsonly"
cetaceanMutationProportions_Restricted=cetaceanMutationProportions[cetaceanMutationProportions$ID==speciesID & condition==condition,c("variable","proportionOfSites")]
colnames(cetaceanMutationProportions_Restricted) <- c("variable","fractionOfAllSegSites")

# note this does have the simple repeats issue still!!! 
mouseMutationProportions <- read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/summaryTables/Mouse7merSpectrum.SummedUpOverAllChrs.NOSTRICT.txt",header=T,sep="\t")
# restrict just to Mmd:
mouseMutationProportions_Restricted <- mouseMutationProportions[mouseMutationProportions$population=="Mmd",c("variable","fractionOfAllSegSites")]


# combine them
merge1_restricted <- merge(mouseMutationProportions_Restricted,cetaceanMutationProportions_Restricted,by="variable",suffixes = c(".mouse",".minke"))


head(merge1_restricted)
# get targets:
merge1_restricted$target <- unlist(lapply(strsplit(as.character(merge1_restricted$variable),"\\."),"[",1))
### now it's time to rescale using targets: this adds mouse proportion
merge2_restricted_addTargets <- merge(merge1_restricted,mouse[,c("target","proportion")],by="target")
merge3_restricted_addTargets <- merge(merge2_restricted_addTargets,allCetaceanTargets[allCetaceanTargets$id==speciesID,c("target","proportion")],by="target",suffixes = c(".MouseTargets",".MinkeTargets"))


# okay time to do the rescaling
# am i doing this right????
merge3_restricted_addTargets$fractionOfAllSegSites.minke.RESCALED.BY.MOUSE <- (merge3_restricted_addTargets$fractionOfAllSegSites.minke / merge3_restricted_addTargets$proportion.MinkeTargets) * merge3_restricted_addTargets$proportion.MouseTargets


# okay let's see what this looks like
diagonalPlot_effectOfRescaling <- ggplot(merge3_restricted_addTargets,aes(x=fractionOfAllSegSites.minke,y=fractionOfAllSegSites.minke.RESCALED.BY.MOUSE,label=variable))+
  geom_point()+
  geom_abline(slope=1,intercept=0)+
  theme_bw()+
  geom_label_repel(data=subset(merge3_restricted_addTargets,fractionOfAllSegSites.minke.RESCALED.BY.MOUSE>0.001 | fractionOfAllSegSites.minke>0.001))+
  ggtitle("effect of rescaling minke proportions by relative abundance of target in mouse\n[(minke fraction/proportion of minke targets) * proportion of mouse targets]\nnote cetaceans had positive bed mask to avoid uncalled sites (repeats)")
diagonalPlot_effectOfRescaling
ggsave(paste(plotdir,"RESCALINGTEST.minke.EffectOfRescalingByMouseTargets.HasSimpleRepeatIssue.png",sep=""),diagonalPlot_effectOfRescaling,height=5,width=7)
# okay so the motif becomes less common

diagonalPlot_mouse_v_RescaledMinke <-  ggplot(merge3_restricted_addTargets,aes(x=fractionOfAllSegSites.mouse,y=fractionOfAllSegSites.minke.RESCALED.BY.MOUSE,label=variable))+
  geom_point()+
  geom_abline(slope=1,intercept=0)+
  theme_bw()+
  geom_label_repel(data=subset(merge3_restricted_addTargets,fractionOfAllSegSites.minke.RESCALED.BY.MOUSE>0.001 | fractionOfAllSegSites.mouse>0.001))+
  ggtitle("mouse vs rescaled minke\n[(minke fraction/proportion of minke targets) * proportion of mouse targets]\nnote cetaceans had positive bed mask to avoid uncalled sites (repeats)")
diagonalPlot_mouse_v_RescaledMinke
ggsave(paste(plotdir,"RESCALINGTEST.minke.vs.mouse.EffectOfRescalingByMouseTargets.HasSimpleRepeatIssue.png",sep=""),diagonalPlot_mouse_v_RescaledMinke,height=5,width=7)



############ sandbox: rescaling by hand ########
# minkeFracOfAllSegSites_TTTAAAA_TTTTAAA_manual=0.0154037511507136 # bacu2 manual from from /Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_spectrum_files/noSTRICTflag/AllCetaceans_melted.IncludesTwoConditions.HetsOnly.and.HetsPlusHoms.NOSTRICTFLAG.txt
# mouseFracOfAllSegSites_TTTAAAA_TTTTAAA_manual=0.00142284945490138  # mmd. from /Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210317_NOSTRICT_7mer/summaryTables/Mouse7merSpectrum.SummedUpOverAllChrs.NOSTRICT.txt
# minkeID="GCF_000493695.1_BalAcu1.0_bacu3"
# minkeProportion_TTTAAAA=mergedTargets_wide[mergedTargets_wide$target=="TTTAAAA",]$GCF_000493695.1_BalAcu1.0_bacu3
# # THIs proportion has an issue: mouse contains lots of repeat regions which are bedmasked away in the cetaceans; so that takes up more of the denominator and so the TTTAAAA motif then seems less but may have a similar count? 
# mouseProportion_TTTAAAA=mergedTargets_wide[mergedTargets_wide$target=="TTTAAAA",]$mm10_nostrict
# 
# rescaledMinkeFraction = (minkeFracOfAllSegSites_TTTAAAA_TTTTAAA_manual/minkeProportion_TTTAAAA)*mouseProportion_TTTAAAA # okay so now 
# rescaledMinkeFraction
# mouseFracOfAllSegSites_TTTAAAA_TTTTAAA_manual
# # okay so the minke rate is still like 6x higher. cool
# # so need to do this on a bigger scale. actually read in the correct files (listed above) 
# # and do this for every mutation type
