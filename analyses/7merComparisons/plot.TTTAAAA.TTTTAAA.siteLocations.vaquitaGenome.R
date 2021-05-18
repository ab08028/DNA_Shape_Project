################ plot TTTAAAA>TTTTAAA sites by position
require(ggplot2)
require(gtools)
plotdir="/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210225_7mer/"

vaqSites=read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210225_7mer/mutyper_variant_files/TTTAAAA.TTTTAAA.Variants.vcf",comment.char="#")
head(vaqSites)
colnames(vaqSites) <- c("chr","pos")
vaqSites$chr <- factor(vaqSites$chr,levels=mixedsort(unique(as.character(vaqSites$chr))))
sitesplot <- ggplot(vaqSites)+
  geom_segment(aes(x=pos,xend=pos+1,y=0,yend=1),size=0.1)+
  facet_wrap(~chr,scales="free_x",ncol=3)+
  ggtitle("Location of TTTAAAA>TTTTAAA variants in vaquita genome")+
  theme_classic()+
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank())+
  xlab("")+
  ylab("")

sitesplot
ggsave(paste(plotdir,"TTTAAAA.TTTTAAA.SitesAlongVaquitaGenome.pdf",sep=""),sitesplot,height=10,width=12)
