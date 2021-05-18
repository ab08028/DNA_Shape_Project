#### beluga find 
require(reshape2)
beluga1 <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_spectrum_files/noSTRICTflag/beluga_whale_GCF_002288925.2_ASM228892v3_dleu1.mutyper.spectrum.NORANDOMIZE.hetsonly.NOSTRICT.txt",header=T)
beluga1_melt <- melt(beluga1)

beluga2 <- read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyper_20210311/mutyper_spectrum_files/noSTRICTflag/beluga_whale_GCF_002288925.2_ASM228892v3_dleu2.mutyper.spectrum.NORANDOMIZE.hetsonly.NOSTRICT.txt",header=T)
beluga2_melt <- melt(beluga2)
### okay so beluga 2 has a single 5 mer but that's it -- its in the variants file too
