########### plot TTTAAAA>TTTTAAA ksfs ##########
plotdir="/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210225_7mer/"

ksfs <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210225_7mer/mutyper_ksfs_files/TTTAAAA.TTTTAAA.ksfs.txt",header=T)

require(ggplot2)
sfsplot <- ggplot(ksfs,aes(x=sample_frequency,y=TTTAAAA.TTTTAAA))+
  geom_col()+ylab("count of TTTAAAA>TTTTAAA mutations")+
  theme_bw()+
  ggtitle(paste("total sites = ",sum(ksfs$TTTAAAA.TTTTAAA)))
sfsplot
ggsave(paste(plotdir,"TTTAAAA.TTTTAAA.ksfs.pdf",sep = ""),sfsplot,height=4,width=8)
