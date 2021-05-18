############# make diagonal plots
# this uses the results of the CompareVaquitaFinWhale7mers.R script which is slow because it needs to process a lot. this has the processing done already
plotdir="/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210225_7mer/"
fw = read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/fin_whales/results/mutyper/downloaddate_20210111/mutyperResults_20210225_7mer/finwhale.7merspectrum.SummedOverAllIntervals.txt",header=T,sep="\t")
vaq = read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210225_7mer/vaquita.7merspectrum.SummedOverAllIntervals.txt",header = T,sep="\t")
mouse = read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyer_results_20210226_7mer/mouse.7merspectrum.SummedOverAllChrs.txt",header = T,sep="\t")
bears=""
######### be very careful with these to check carefully! some sites had bad write out due to ">" in label! Messed up some things. But now it's fixed.

############ merge vaquita and fw ###############
# need to merge based on variable and just one FW population
GOCfw_plus_vaq <- merge(vaq[,c("variable", "ThreemerMutationType","centralMutationType" ,"labelOfInterest","fractionOfAllSegSites")],fw[fw$population=="GOC",c("variable","labelOfInterest", "ThreemerMutationType","centralMutationType" ,"fractionOfAllSegSites")],by=c("variable","labelOfInterest","centralMutationType","ThreemerMutationType"),suffixes = c(".vaq",".fw"))

ENPfw_plus_vaq <- merge(vaq[,c("variable", "ThreemerMutationType","centralMutationType" ,"labelOfInterest","fractionOfAllSegSites")],fw[fw$population=="ENP",c("variable","labelOfInterest", "ThreemerMutationType","centralMutationType" ,"fractionOfAllSegSites")],by=c("variable","labelOfInterest","centralMutationType","ThreemerMutationType"),suffixes = c(".vaq",".fw"))

############ vaquita and fin whale plots ###########
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
ggsave(paste(plotdir,"Diagonal.GOCfinwhale.vs.Vaquita.Comparison.7mers.png",sep=""),GOC_vaquita_diagonalPlot,height=7,width=12)

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


############### merge mouse and vaquita ##################
## just going to use Mmd for now 
Mmd_mouse_plus_vaq <- merge(vaq[,c("variable", "ThreemerMutationType","centralMutationType" ,"labelOfInterest","fractionOfAllSegSites")],mouse[mouse$population=="Mmd",c("variable","labelOfInterest", "ThreemerMutationType","centralMutationType" ,"fractionOfAllSegSites")],by=c("variable","labelOfInterest","centralMutationType","ThreemerMutationType"),suffixes = c(".vaq",".Mmd"))

Mmd_vaquita_diagonalPlot <- ggplot(Mmd_mouse_plus_vaq,aes(x=fractionOfAllSegSites.Mmd,y=fractionOfAllSegSites.vaq,color=centralMutationType,label=variable))+
  geom_point(size=1.2,alpha=0.7)+
  geom_abline(intercept=0,slope=1)+
  theme_bw() +
  geom_label_repel(data=subset(Mmd_mouse_plus_vaq,fractionOfAllSegSites.vaq>=0.00075), nudge_y=0.003,size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",direction="x",show.legend = F)+
  xlab("Mouse 7mer proportion of all segregating sites (M. musculus domesticus)")+
  ylab("Vaquita 7mer proportion of all segregating sites")+
  ggtitle("Comparing mouse and vaquita 7mer proportions")+
  theme(text=element_text(size=14))+
  guides(colour = guide_legend(override.aes = list(size=3))) # this makes the dots in the legend larger
Mmd_vaquita_diagonalPlot
ggsave(paste(plotdir,"Diagonal.Mmdmouse.vs.Vaquita.Comparison.7mers.pdf",sep=""),Mmd_vaquita_diagonalPlot,height=7,width=12)


############### merge mouse and fin whale #####################
## just going to use Mmd for now 
Mmd_mouse_plus_GOCfw <- merge(fw[fw$population=="GOC",c("variable", "ThreemerMutationType","centralMutationType" ,"labelOfInterest","fractionOfAllSegSites")],mouse[mouse$population=="Mmd",c("variable","labelOfInterest", "ThreemerMutationType","centralMutationType" ,"fractionOfAllSegSites")],by=c("variable","labelOfInterest","centralMutationType","ThreemerMutationType"),suffixes = c(".fw",".Mmd"))

Mmd_GOCfw_diagonalPlot <- ggplot(Mmd_mouse_plus_GOCfw,aes(y=fractionOfAllSegSites.Mmd,x=fractionOfAllSegSites.fw,color=centralMutationType,label=variable))+
  geom_point(size=1.2,alpha=0.7)+
  geom_abline(intercept=0,slope=1)+
  theme_bw() +
  geom_label_repel(data=subset(Mmd_mouse_plus_GOCfw,fractionOfAllSegSites.Mmd>=0.0005), nudge_y=0.003,size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",direction="x",show.legend = F)+
  geom_label_repel(data=subset(Mmd_mouse_plus_GOCfw,fractionOfAllSegSites.fw>=5e-4), nudge_y=-0.003,size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",direction="x",show.legend = F)+
  xlab("Fin whale 7mer proportion of all segregating sites (Gulf of California)")+
  ylab("Mouse 7mer proportion of all segregating sites (M. musculus domesticus)")+
  ggtitle("Comparing mouse and fin whale 7mer proportions")+
  theme(text=element_text(size=14))+
  guides(colour = guide_legend(override.aes = list(size=3))) # this makes the dots in the legend larger
Mmd_GOCfw_diagonalPlot
ggsave(paste(plotdir,"Diagonal.Mmdmouse.vs.GOCfinwhale.Comparison.7mers.pdf",sep=""),Mmd_GOCfw_diagonalPlot,height=7,width=12)


Mmd_mouse_plus_ENPfw <- merge(fw[fw$population=="ENP",c("variable", "ThreemerMutationType","centralMutationType" ,"labelOfInterest","fractionOfAllSegSites")],mouse[mouse$population=="Mmd",c("variable","labelOfInterest", "ThreemerMutationType","centralMutationType" ,"fractionOfAllSegSites")],by=c("variable","labelOfInterest","centralMutationType","ThreemerMutationType"),suffixes = c(".fw",".Mmd"))

Mmd_ENPfw_diagonalPlot <- ggplot(Mmd_mouse_plus_ENPfw,aes(y=fractionOfAllSegSites.Mmd,x=fractionOfAllSegSites.fw,color=centralMutationType,label=variable))+
  geom_point(size=1.2,alpha=0.7)+
  geom_abline(intercept=0,slope=1)+
  theme_bw() +
  geom_label_repel(data=subset(Mmd_mouse_plus_ENPfw,fractionOfAllSegSites.Mmd>=0.0005), nudge_y=0.003,size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",direction="x",show.legend = F)+  
  geom_label_repel(data=subset(Mmd_mouse_plus_ENPfw,fractionOfAllSegSites.fw>=5e-4), nudge_y=-0.003,size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",direction="x",show.legend = F)+

  xlab("Fin whale 7mer proportion of all segregating sites (ENP)")+
  ylab("Mouse 7mer proportion of all segregating sites (M. musculus domesticus)")+
  ggtitle("Comparing mouse and fin whale 7mer proportions")+
  theme(text=element_text(size=14))+
  guides(colour = guide_legend(override.aes = list(size=3))) # this makes the dots in the legend larger
Mmd_ENPfw_diagonalPlot
ggsave(paste(plotdir,"Diagonal.Mmdmouse.vs.ENPfinwhale.Comparison.7mers.pdf",sep=""),Mmd_ENPfw_diagonalPlot,height=7,width=12)

############# adding in bears ##############
