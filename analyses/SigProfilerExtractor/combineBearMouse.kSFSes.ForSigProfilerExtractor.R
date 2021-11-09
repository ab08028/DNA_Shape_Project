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


# need to merge:
merged <- full_join(bears,mice,by="Mutation Types")
head(merged)
dim(merged)

write.table(merged,paste0(outdir,todaysdate,".combined.Bear.Mice.SpectrumPerFrequency.ForSigProfilerExtractor.txt"),col.names = T,row.names = F,quote=F,sep="\t")
