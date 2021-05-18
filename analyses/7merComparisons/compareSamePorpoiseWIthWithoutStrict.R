require(ggplot2)
require(reshape2)
require(gtools)
require(ggrepel)
######### compare just two porpoises with and without scrict
# mapped to finless porpoise genome
finlessWithStrict=read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_spectrum_files/indopacific_finless_porpoise_GCF_003031525.1_Neophocaena_asiaeorientalis_V1_npho1.mutyper.spectrum.NORANDOMIZE.hetsonly.txt",header=T)
finlessWithStrict_melt <- melt(finlessWithStrict)
finlessWithStrict_melt$proportion <- finlessWithStrict_melt$value/sum(finlessWithStrict_melt$value)

finlessWithoutStrict=read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_spectrum_files/noSTRICTflag/indopacific_finless_porpoise_GCF_003031525.1_Neophocaena_asiaeorientalis_V1_npho1.mutyper.spectrum.NORANDOMIZE.hetsonly.NOSTRICT.txt",header=T)
finlessWithoutStrict_melt <- melt(finlessWithoutStrict)
finlessWithoutStrict_melt$proportion <- finlessWithoutStrict_melt$value/sum(finlessWithoutStrict_melt$value)

# mapped to vaquita
finlessMappedToVaq=read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_spectrum_files/indopacific_porpoise_vaquita_outgroup_npho1.vcf.gz.mutyper.spectrum.NORANDOMIZE.hetsonly.txt",header=T)
finlessMappedToVaq_melt <- melt(finlessMappedToVaq)
finlessMappedToVaq_melt$proportion <- finlessMappedToVaq_melt$value/sum(finlessMappedToVaq_melt$value)
####################### merge ############
finless_withWithoutStrict <- merge(finlessWithStrict_melt,finlessWithoutStrict_melt,by=c("variable","sample"),suffixes=c(".WithStrictNoRepeats",".WithoutStrict.HasRepeats"))

finlessWithoutStrict_vsMappedToVaquita <- merge(finlessWithoutStrict_melt,finlessMappedToVaq_melt,by=c("variable","sample"),suffixes=c(".WithoutStrict.HasRepeats",".MappedToVaquita.HasRepeats"))


sameIndividualWithWithoutStrict <- ggplot(finless_withWithoutStrict,aes(x=proportion.WithStrictNoRepeats,y=proportion.WithoutStrict.HasRepeats,label=variable))+
  geom_point()+
  geom_abline(slope=1,intercept=0)+
  theme_bw()+
  xlab("Proportion when soft-masked repeats excluded")+
  ylab("Proportion when soft-masked repeats are included")+
  geom_label_repel(data=subset(finless_withWithoutStrict,proportion.WithoutStrict.HasRepeats>=0.002), size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  ggtitle("Finless porpoise, mapped to finless porpoise with and without soft-masked repeats included/excluded")
  
finlessWithAndWithoutRepeatsPlot <- ggplot(finlessWithoutStrict_vsMappedToVaquita,aes(x=proportion.WithoutStrict.HasRepeats,y=proportion.MappedToVaquita.HasRepeats,label=variable))+
  geom_point()+
  geom_abline(slope=1,intercept=0)+
  theme_bw()+
  xlab("Proportion when mapped to the finless porpoise (soft-masked repeats INCLUDED)")+
  ylab("Proportion when mapped to the vaquita (soft-masked repeats INCLUDED)")+
  geom_label_repel(data=subset(finlessWithoutStrict_vsMappedToVaquita,proportion.WithoutStrict.HasRepeats>=0.002), size = 3, box.padding   = 1.5,point.padding = 0.5, force  = 100,segment.size  = 0.2,segment.color = "grey50",show.legend = F)+
  ggtitle("Finless porpoise, mapped to finless porpoise and vaquita with and without soft-masked repeats included/excluded")
finlessWithAndWithoutRepeatsPlot
