require(ggplot2)
require(ggrepel)
############# these use old strict calls --d ont want ######
### want to modify so using without the --strict flag! Don't have this yet for fin whale or bear
############# make diagonal plots
# this uses the results of the CompareVaquitaFinWhale7mers.R script which is slow because it needs to process a lot. this has the processing done already
plotdir="/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/diagonal_comparison_plots/"
#USED STRICT: fw = read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/fin_whales/results/mutyper/downloaddate_20210111/mutyperResults_20210225_7mer/finwhale.7merspectrum.SummedOverAllIntervals.txt",header=T,sep="\t") ## STILL USES 
vaq = read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210225_7mer/vaquita.7merspectrum.SummedOverAllIntervals.txt",header = T,sep="\t")
## USED STRICT : mouse = read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyer_results_20210226_7mer/mouse.7merspectrum.SummedOverAllChrs.txt",header = T,sep="\t")
## USED STRICT :  bears=read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210304_7mer/mapped_to_brown_bear/summaryTables/mapped_to_brown_bear.AllPops.bears.7merspectrum.SummedOverAllIntervals.ALLFREQS.txt",header=T,sep="\t") # just doing mapped to bb for now but can later compare ref genomes

#bears_mappedToPolarBear = read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210304_7mer/mapped_to_polar_bear/summaryTables/mapped_to_polar_bear.AllPops.bears.7merspectrum.SummedOverAllIntervals.ALLFREQS.txt",header=T,sep="\t")

######### be very careful with these to check carefully! some sites had bad write out due to ">" in label! Messed up some things. But now it's fixed.

############ merge vaquita and fw ###############
# need to merge based on variable and just one FW population
GOCfw_plus_vaq <- merge(vaq[,c("variable", "ThreemerMutationType","centralMutationType" ,"labelOfInterest","fractionOfAllSegSites")],fw[fw$population=="GOC",c("variable","labelOfInterest", "ThreemerMutationType","centralMutationType" ,"fractionOfAllSegSites")],by=c("variable","labelOfInterest","centralMutationType","ThreemerMutationType"),suffixes = c(".vaq",".fw"))

########### fin whale vs vaq #########
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
ggsave(paste(plotdir,"Diagonal.GOCfinwhale.vs.Vaquita.Comparison.7mers.png",sep=""),GOC_vaquita_diagonalPlot,height=7,width=12)
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
ggsave(paste(plotdir,"Diagonal.ENPfinwhale.vs.Vaquita.Comparison.7mers.png",sep=""),ENP_vaquita_diagonalPlot,height=7,width=12)

############## ENP vs GOC fin whales ##############
ENPfw_plus_GOCfw <- merge(fw[fw$population=="GOC",c("variable", "ThreemerMutationType","centralMutationType" ,"labelOfInterest","fractionOfAllSegSites")],fw[fw$population=="ENP",c("variable","labelOfInterest", "ThreemerMutationType","centralMutationType" ,"fractionOfAllSegSites")],by=c("variable","labelOfInterest","centralMutationType","ThreemerMutationType"),suffixes = c(".GOCfw",".ENPfw"))

ENPfw_vsGOCfw_diagonalPlot <- ggplot(ENPfw_plus_GOCfw,aes(x=fractionOfAllSegSites.ENPfw,y=fractionOfAllSegSites.GOCfw,color=centralMutationType,label=variable))+
  geom_point(size=1.2,alpha=0.7)+
  geom_abline(intercept=0,slope=1)+
  theme_bw() +
  geom_label_repel(data=subset(ENPfw_plus_GOCfw,fractionOfAllSegSites.GOCfw>=5e-4), size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  xlab("ENP fin whale")+
  ylab("GOC fin whale")+
  ggtitle("Comparing fin whale populations")+
  theme(text=element_text(size=14))+
  guides(colour = guide_legend(override.aes = list(size=3))) # this makes the dots in the legend larger
ENPfw_vsGOCfw_diagonalPlot
ggsave(paste(plotdir,"Diagonal.ENPfinwhale.vs.GOCfinwhale.Comparison.7mers.png",sep=""),ENPfw_vsGOCfw_diagonalPlot,height=7,width=12)

############### merge mouse and vaquita ##################
## just going to use Mmd for now 
Mmd_mouse_plus_vaq <- merge(vaq[,c("variable", "ThreemerMutationType","centralMutationType" ,"labelOfInterest","fractionOfAllSegSites")],mouse[mouse$population=="Mmd",c("variable","labelOfInterest", "ThreemerMutationType","centralMutationType" ,"fractionOfAllSegSites")],by=c("variable","labelOfInterest","centralMutationType","ThreemerMutationType"),suffixes = c(".vaq",".Mmd"))

Mmd_vaquita_diagonalPlot <- ggplot(Mmd_mouse_plus_vaq,aes(x=fractionOfAllSegSites.Mmd,y=fractionOfAllSegSites.vaq,color=centralMutationType,label=variable))+
  geom_point(size=1.2,alpha=0.7)+
  geom_abline(intercept=0,slope=1)+
  theme_bw() +
  geom_label_repel(data=subset(Mmd_mouse_plus_vaq,fractionOfAllSegSites.vaq>=0.00075 | fractionOfAllSegSites.Mmd>=0.0005 |variable =="CCACGTG.CCATGTG"),size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  #geom_label_repel(data=subset(Mmd_mouse_plus_vaq,fractionOfAllSegSites.Mmd>=0.0005), nudge_y=-0.003,size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",direction="x",show.legend = F)+
  #geom_label_repel(data=subset(Mmd_mouse_plus_vaq,variable =="CCACGTG.CCATGTG"), nudge_y=-0.006,size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",direction="x",show.legend = F)+
  xlab("Mouse 7mer proportion of all segregating sites (M. musculus domesticus)")+
  ylab("Vaquita 7mer proportion of all segregating sites")+
  ggtitle("Comparing mouse and vaquita 7mer proportions")+
  theme(text=element_text(size=14))+
  guides(colour = guide_legend(override.aes = list(size=3))) # this makes the dots in the legend larger
Mmd_vaquita_diagonalPlot
ggsave(paste(plotdir,"Diagonal.Mmdmouse.vs.Vaquita.Comparison.7mers.png",sep=""),Mmd_vaquita_diagonalPlot,height=7,width=12)


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
ggsave(paste(plotdir,"Diagonal.Mmdmouse.vs.GOCfinwhale.Comparison.7mers.png",sep=""),Mmd_GOCfw_diagonalPlot,height=7,width=12)


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
ggsave(paste(plotdir,"Diagonal.Mmdmouse.vs.ENPfinwhale.Comparison.7mers.png",sep=""),Mmd_ENPfw_diagonalPlot,height=7,width=12)

# compare mouse species?
############# adding in bears ##############
## compare ABC vs EUR and ABC vs PB and EUR vs PB?
# compare to Mmd 

####### EUR BB vs ABC BB ###############

abc_plus_eurBB <- merge(bears[bears$population=="ABC",c("variable", "ThreemerMutationType","centralMutationType" ,"labelOfInterest","fractionOfAllSegSites")],bears[bears$population=="EUR",c("variable","labelOfInterest", "ThreemerMutationType","centralMutationType" ,"fractionOfAllSegSites")],by=c("variable","labelOfInterest","centralMutationType","ThreemerMutationType"),suffixes = c(".ABCBB",".EURBB"))

abc_plus_eurBB_diagonalPlot <- ggplot(abc_plus_eurBB,aes(x=fractionOfAllSegSites.EURBB,y=fractionOfAllSegSites.ABCBB,color=centralMutationType,label=variable))+
  geom_point(size=1.2,alpha=0.7)+
  geom_abline(intercept=0,slope=1)+
  theme_bw() +
  xlab("Brown Bear (EUR) 7mer proportion of all segregating sites")+
  ylab("Brown Bear (ABC) 7mer proportion of all segregating sites")+
  ggtitle("Comparing brown bear 7mer proportions")+
  theme(text=element_text(size=14))+
  guides(colour = guide_legend(override.aes = list(size=3))) # this makes the dots in the legend larger
abc_plus_eurBB_diagonalPlot
ggsave(paste(plotdir,"Diagonal.ABCBrownBear.vs.EURBrownBear.mappedtobrownbear.Comparison.7mers.png",sep=""),abc_plus_eurBB_diagonalPlot,height=7,width=12)

abc_plus_PB <- merge(bears[bears$population=="ABC",c("variable", "ThreemerMutationType","centralMutationType" ,"labelOfInterest","fractionOfAllSegSites")],bears[bears$population=="PB",c("variable","labelOfInterest", "ThreemerMutationType","centralMutationType" ,"fractionOfAllSegSites")],by=c("variable","labelOfInterest","centralMutationType","ThreemerMutationType"),suffixes = c(".ABCBB",".PB"))


abc_plus_PB_diagonalPlot <- ggplot(abc_plus_PB,aes(x=fractionOfAllSegSites.ABCBB,y=fractionOfAllSegSites.PB,color=centralMutationType,label=variable))+
  geom_point(size=1.2,alpha=0.7)+
  geom_abline(intercept=0,slope=1)+
  theme_bw() +
  xlab("Brown Bear (ABC) 7mer proportion of all segregating sites")+
  ylab("Polar Bear 7mer proportion of all segregating sites")+
  ggtitle("Comparing brown bear and polar bear 7mer proportions")+
  theme(text=element_text(size=14))+
  guides(colour = guide_legend(override.aes = list(size=3))) # this makes the dots in the legend larger
abc_plus_PB_diagonalPlot
ggsave(paste(plotdir,"Diagonal.ABCBrownBear.vs.PolarBear.mappedtobrownbear.Comparison.7mers.png",sep=""),abc_plus_PB_diagonalPlot,height=7,width=12)

########## facet the plot:
abc_plus_PB_diagonalPlot_facet <- ggplot(abc_plus_PB,aes(x=fractionOfAllSegSites.ABCBB,y=fractionOfAllSegSites.PB,color=ThreemerMutationType,label=variable))+
  geom_point(size=1.2,alpha=0.7)+
  geom_abline(intercept=0,slope=1)+
  theme_bw() +
  facet_wrap(~centralMutationType,scales="free")+
  xlab("Brown Bear (ABC) 7mer proportion of all segregating sites")+
  ylab("Polar Bear 7mer proportion of all segregating sites")+
  ggtitle("Comparing brown bear and polar bear 7mer proportions")+
  theme(text=element_text(size=14))+
  guides(colour = guide_legend(override.aes = list(size=3))) # this makes the dots in the legend larger
abc_plus_PB_diagonalPlot_facet

ggsave(paste(plotdir,"Diagonal.ABCBrownBear.vs.PolarBear.mappedtobrownbear.Comparison.7mers.FACETED.png",sep=""),abc_plus_PB_diagonalPlot_facet,height=7,width=12)

eur_plus_PB <- merge(bears[bears$population=="EUR",c("variable", "ThreemerMutationType","centralMutationType" ,"labelOfInterest","fractionOfAllSegSites")],bears[bears$population=="PB",c("variable","labelOfInterest", "ThreemerMutationType","centralMutationType" ,"fractionOfAllSegSites")],by=c("variable","labelOfInterest","centralMutationType","ThreemerMutationType"),suffixes = c(".EURBB",".PB"))

eur_plus_PB_diagonalPlot <- ggplot(eur_plus_PB,aes(x=fractionOfAllSegSites.EURBB,y=fractionOfAllSegSites.PB,color=centralMutationType,label=variable))+
  geom_point(size=1.2,alpha=0.7)+
  geom_abline(intercept=0,slope=1)+
  theme_bw() +
  xlab("Brown Bear (EUR) 7mer proportion of all segregating sites")+
  ylab("Polar Bear 7mer proportion of all segregating sites")+
  ggtitle("Comparing brown bear and polar bear 7mer proportions")+
  theme(text=element_text(size=14))+
  guides(colour = guide_legend(override.aes = list(size=3))) # this makes the dots in the legend larger
eur_plus_PB_diagonalPlot
ggsave(paste(plotdir,"Diagonal.EURBrownBear.vs.PolarBear.mappedtobrownbear.Comparison.7mers.png",sep=""),eur_plus_PB_diagonalPlot,height=7,width=12)

######## mmd plus EUR BB ########
Mmd_mouse_plus_EURBB <- merge(bears[bears$population=="EUR",c("variable", "ThreemerMutationType","centralMutationType" ,"labelOfInterest","fractionOfAllSegSites")],mouse[mouse$population=="Mmd",c("variable","labelOfInterest", "ThreemerMutationType","centralMutationType" ,"fractionOfAllSegSites")],by=c("variable","labelOfInterest","centralMutationType","ThreemerMutationType"),suffixes = c(".EURBB",".Mmd"))

Mmd_EURBB_diagonalPlot <- ggplot(Mmd_mouse_plus_EURBB,aes(x=fractionOfAllSegSites.EURBB,y=fractionOfAllSegSites.Mmd,color=centralMutationType,label=variable))+
  geom_point(size=1.2,alpha=0.7)+
  geom_abline(intercept=0,slope=1)+
  theme_bw() +
  geom_label_repel(data=subset(Mmd_mouse_plus_EURBB,fractionOfAllSegSites.Mmd>=0.0005), nudge_y=0.003,size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",direction="x",show.legend = F)+
  geom_label_repel(data=subset(Mmd_mouse_plus_EURBB,fractionOfAllSegSites.EURBB>=5e-4), nudge_y=-0.003,size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",direction="x",show.legend = F)+ # adding in label for sites that are at high freq for BB
  
  xlab("Brown Bear (EUR) 7mer proportion of all segregating sites")+
  ylab("Mouse 7mer proportion of all segregating sites (M. musculus domesticus)")+
  ggtitle("Comparing bear and mouse 7mer proportions")+
  theme(text=element_text(size=14))+
  guides(colour = guide_legend(override.aes = list(size=3))) # this makes the dots in the legend larger
Mmd_EURBB_diagonalPlot
ggsave(paste(plotdir,"Diagonal.Mmdmouse.vs.EURBrownBear.mappedtobrownbear.Comparison.7mers.png",sep=""),Mmd_EURBB_diagonalPlot,height=7,width=12)


############ focusing just on E box domain which is CANNTG (a hexamer so anythng on either side ok too) ########
#### E box mouse and bears ######
Mmd_mouse_plus_EURBB_EBOXONLY_CANNTG <- Mmd_mouse_plus_EURBB[grep("CA[ACGT][ACGT]TG",Mmd_mouse_plus_EURBB$variable),]

EBOX_Mmd_EURBB_diagonalPlot <- ggplot(Mmd_mouse_plus_EURBB_EBOXONLY_CANNTG,aes(x=fractionOfAllSegSites.EURBB,y=fractionOfAllSegSites.Mmd,color=centralMutationType,label=variable))+
  geom_point(size=1.2,alpha=0.7)+
  geom_abline(intercept=0,slope=1)+
  theme_bw() +
  #geom_label_repel(data=subset(Mmd_mouse_plus_EURBB,fractionOfAllSegSites.Mmd>=0.0005), nudge_y=0.003,size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",direction="x",show.legend = F)+
  #geom_label_repel(data=subset(Mmd_mouse_plus_EURBB,fractionOfAllSegSites.EURBB>=5e-4), nudge_y=-0.003,size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",direction="x",show.legend = F)+ # adding in label for sites that are at high freq for BB
  
  xlab("Brown Bear (EUR) 7mer proportion of all segregating sites")+
  ylab("Mouse 7mer proportion of all segregating sites (M. musculus domesticus)")+
  ggtitle("EBOX containing 7mers only (CANNCG) Comparing bear and mouse 7mer proportions")+
  theme(text=element_text(size=14))+
  guides(colour = guide_legend(override.aes = list(size=3))) # this makes the dots in the legend larger
EBOX_Mmd_EURBB_diagonalPlot
# huh! not sure what to make of this yet.