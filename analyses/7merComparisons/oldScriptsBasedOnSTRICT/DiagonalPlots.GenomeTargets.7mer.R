require(ggplot2)
require(ggrepel)
plotdir="/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/diagonal_comparison_plots/targetsInGenomes/"
######## diagonal plots of targets 7mer
# get proprtions
vaquita <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210225_7mer/mutyper_target_files/mutyper.targets.7mer.txt",header=F)
colnames(vaquita) <- c("target","count.vaq")
vaquita$proportion.vaq <- vaquita$count.vaq/sum(vaquita$count.vaq)

fw_minke_wholegenome <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/fin_whales/results/mutyper/downloaddate_20210111/mutyperResults_20210225_7mer/mutyper_target_files/whole_genome/mutyper.targets.7mer.WHOLEMINKEGENOME.strict.txt",header=F)
colnames(fw_minke_wholegenome) <- c("target","count.minkeWG")
fw_minke_wholegenome$proportion.minkeWG<- fw_minke_wholegenome$count.minkeWG/sum(fw_minke_wholegenome$count.minkeWG)

fw_minke_neutralonly <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/fin_whales/results/mutyper/downloaddate_20210111/mutyperResults_20210225_7mer/mutyper_target_files/neutral_only/mutyper.targets.7mer.NEUTRALBEDCOORDSONLY.NoteAreVeryChoppyCoords.strict.txt",header=F)
colnames(fw_minke_neutralonly) <- c("target","count.minkeNeutral")
fw_minke_neutralonly$proportion.minkeNeutral <- fw_minke_neutralonly$count.minkeNeutral/sum(fw_minke_neutralonly$count.minkeNeutral)

# note for the following there *is* a header, but the above there is *not* 
brownbear=read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210304_7mer/mapped_to_brown_bear/mutyper_target_files/ALLINTERVALS.brown_bear.mutyper.targets.7mer.txt",header=T)
colnames(brownbear) <- c("target","count.brown_bear")
brownbear$proportion.brown_bear <- brownbear$count.brown_bear/sum(brownbear$count.brown_bear)


polarbear=read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210304_7mer/mapped_to_polar_bear/mutyper_target_files/ALLINTERVALS.polar_bear.mutyper.targets.7mer.txt",header=T)
colnames(polarbear) <- c("target","count.polar_bear")
polarbear$proportion.polar_bear <- polarbear$count.polar_bear/sum(polarbear$count.polar_bear)

mouse=read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyer_results_20210226_7mer/mutyper_target_files/ALLINTERVALS.mouse.mutyper.targets.7mer.txt",header=T)
colnames(mouse) <- c("target","count.mouse")
mouse$proportion.mouse <- mouse$count.mouse/sum(mouse$count.mouse)


interestingmotifs=c("TTTAAAA","TTTTAAA","AAAATTT","AAATTTT")
interestingmotifs2=c("CCACCAC","CCAACAC") # this is in fin whale CCACGTG CCATGTG are most abundant 7mers so these translate to CACCAC and CAACAC when AT --> A and CG --> C
# merge them:
m1 <- merge(vaquita,fw_minke_wholegenome,by="target")
m2 <- merge(m1,fw_minke_neutralonly,by="target")
m3 <- merge(m2,brownbear,by="target")
m4 <- merge(m3,polarbear,by="target")
m5 <- merge(m4,mouse,by="target")
finalmerge = m5

finalmerge$centralBP <- substr(finalmerge$target,4,4)
finalmerge$centralthreemer <- substr(finalmerge$target,3,5)
head(finalmerge)
########## compare minke (fin whale) whole genome vs putatively neutral regions which we got rid of CpGs etc #########
p1 <- ggplot(finalmerge,aes(x=proportion.minkeWG,y=proportion.minkeNeutral,color=centralthreemer,label=target ))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  theme_bw()+  
  geom_label_repel(data=subset(finalmerge,target %in% interestingmotifs), nudge_y=0.003,size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",direction="x",show.legend = F)+
  ggtitle("Minke whale neutral vs whole genome (losing CpGs because we exclude Cpg Islands")

p1
ggsave(paste(plotdir,"minkeNeutral.vs.minkeWG.7mer.targets.png",sep=""),p1,height=5,width=8)
############### vaquita vs minke whole genome ###############
p2 <- ggplot(finalmerge,aes(x=proportion.vaq,y=proportion.minkeWG,color=centralthreemer,label=target ))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  theme_bw()+  
  geom_label_repel(data=subset(finalmerge,proportion.vaq>=0.001),size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  ggtitle("minke whole genome vs vaquita targets (note AT are collapsed to A, CG to C)")

p2
ggsave(paste(plotdir,"minkeWG.vs.vaquita.7mer.targets.png",sep=""),p2,height=5,width=8)

######## vaquita vs minke neutral regions ##########
p3 <- ggplot(finalmerge,aes(x=proportion.vaq,y=proportion.minkeNeutral,color=centralthreemer,label=target ))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  theme_bw()+  
  geom_label_repel(data=subset(finalmerge,proportion.vaq>=0.001),size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  ggtitle("minke neutral regions vs vaquita 7mer targets (note AT are collapsed to A, CG to C)")


p3
ggsave(paste(plotdir,"minkeNeutral.vs.vaquita.7mer.targets.png",sep=""),p3,height=5,width=8)


########## brown bear vs polar bear #############
p4a <- ggplot(finalmerge,aes(x=proportion.polar_bear,y=proportion.brown_bear,color=centralthreemer,label=target ))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  theme_bw()+  
  geom_label_repel(data=subset(finalmerge,abs(proportion.brown_bear-proportion.polar_bear)>=1e-5),size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  ggtitle("brown bear vs polar bear genome targets (note AT are collapsed to A, CG to C)")


p4a
ggsave(paste(plotdir,"brownbear.vs.polarbear.7mer.targets.WithLabels.png",sep=""),p4a,height=5,width=8)

p4b <- ggplot(finalmerge,aes(x=proportion.polar_bear,y=proportion.brown_bear,color=centralthreemer,label=target ))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  theme_bw()+  
  ggtitle("brown bear vs polar bear genome targets (note AT are collapsed to A, CG to C)")


p4b
ggsave(paste(plotdir,"brownbear.vs.polarbear.7mer.targets.NoLabels.png",sep=""),p4b,height=5,width=8)

################ mouse vs vaquita ##################

p5 <- ggplot(finalmerge,aes(x=proportion.mouse,y=proportion.vaq,color=centralthreemer,label=target ))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  theme_bw()+  
  geom_label_repel(data=subset(finalmerge,proportion.vaq>=0.001),size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  ggtitle("vaquita vs mouse genome targets (note AT are collapsed to A, CG to C)")


p5
ggsave(paste(plotdir,"vaquita.vs.mouse.7mer.targets.WithLabels.png",sep=""),p5,height=5,width=8)

################ mouse vs minke whale (fin) ##################
# also want to annotate the CACGTG and CATGTG domains which would be CACCAC and CAACAC this lingo in which AT -> A and CG --> C
# 
p6 <- ggplot(finalmerge,aes(x=proportion.mouse,y=proportion.minkeWG,color=centralthreemer,label=target ))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  theme_bw()+  
  geom_label_repel(data=subset(finalmerge,proportion.mouse>=0.00075),size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  geom_label_repel(data=subset(finalmerge,target %in% c("CCACCAC","CCAACAC")),size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = .8,segment.color = "yellow",show.legend = F)+
  ggtitle("minke whale (fin whale) vs mouse genome targets (note AT are collapsed to A, CG to C)")


p6
ggsave(paste(plotdir,"minkewhale.vs.mouse.7mer.targets.WithLabels.png",sep=""),p6,height=5,width=8)

########### brown bear vs mouse ##########
p7 <- ggplot(finalmerge,aes(x=proportion.mouse,y=proportion.brown_bear,color=centralthreemer,label=target ))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  theme_bw()+  
  geom_label_repel(data=subset(finalmerge,proportion.mouse>=0.00075),size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  ggtitle("brown bear vs mouse genome targets (note AT are collapsed to A, CG to C)")


p7
ggsave(paste(plotdir,"brownbear.vs.mouse.7mer.targets.WithLabels.png",sep=""),p6,height=5,width=8)

