########## combine bear and mouse sfses for sig profiler extractor ###########
# have to use tibble version of stuff
require(tidyverse)
require(dplyr)

todaysdate=format(Sys.time(), "%Y%m%d")
outdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/SigProfilerExtractor/input_data/"

bears <- read_tsv("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210401_3mer_nostrict/mapped_to_brown_bear/summaryTables/All_Intervals_kmerCountsForSigProfilerExtractor.ALLFREQUENCIESSEPARATED.bears.txt",col_names =T)
head(bears)


mice <- read_tsv("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210706_NOSTRICT_3mer_REPEATMASKED/summaryTables/All_Intervals_kmerCountsForSigProfilerExtractor.ALLFREQUENCIESSEPARATED.mice.txt",col_names = T)
head(mice)

vaquita <- read_tsv("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210216/summaryTables/All_Intervals_kmerCountsForSigProfilerExtractor.ALLFREQUENCIESSEPARATED.vaquita.txt",col_names=T)

fin_whale <- read_tsv("/Users/annabelbeichman/Documents/UW/WhaleProject/fin_whales/results/mutyper/downloaddate_20210111/mutyperResults_20210804_NOSTRICT_3mer/summaryTables/All_Intervals_kmerCountsForSigProfilerExtractor.NOTEONLYNEUTRALREGIONS.HIGHLYFILTERED.NOTCOMPARABLETOOTHERCETACEANSASIS.ALLFREQUENCIESSEPARATED.fin_whale.txt",col_names=T)


# need to merge:
merge1 <- full_join(bears,mice,by="Mutation Types")
merge2 <- full_join(merge1,vaquita,by="Mutation Types")
merge3 <- full_join(merge2,fin_whale,by="Mutation Types")
names(merge3)
merged <- merge3
#merged <- full_join(bears,mice,by="Mutation Types")
head(merged)
dim(merged)

write.table(merged,paste0(outdir,todaysdate,".combined.Bear.Mice.Vaquita.FinWhales.SpectrumPerFrequency.ForSigProfilerExtractor.txt"),col.names = T,row.names = F,quote=F,sep="\t")
