# combine into one matrix for sigprofilerextractor

todaysdate=format(Sys.time(), "%Y%m%d")
todaysdate
outdir=paste0("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/SigProfilerExtractor/input_data/")
dir.create(outdir)
# use read_delim in setad of read.table() because it makes a tibble and treats the `Mutation types` header way better 
mice=read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210706_NOSTRICT_3mer_REPEATMASKED/summaryTables/mouse.All_Intervals_Counts.perPopulation.ALLFREQS.nostrict.repeatmasked.CentralAsRevCompedToTs.SIGPROFILEREXTRACTORFORMAT.txt",header=T,sep="\t") # repeat masked
head(mice)
bears=read.table("/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210401_3mer_nostrict/mapped_to_brown_bear/summaryTables/All_Intervals_Counts.perPopulation.ref.brown_bear.ALLFREQS.nostrict.CentralAsRevCompedToTs.SIGPROFILEREXTRACTORFORMAT.txt",header=T,sep="\t") # not repeat masked I don't think

vaquita=read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/downloadDate_20210216/mutyperResults_20210216/summaryTables/vaquita.All_Intervals_Counts.perPopulation.ALLFREQS.CentralAsRevCompedToTs.SIGPROFILEREXTRACTORFORMAT.txt",header=T,sep="\t") # repeat masked (though targets aren't)

# WAITING ON:
bulk_cetaceans=read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyperResults_20210804_3mer/mutyper_spectrum_files/bulkCetaceans.All_Intervals_Counts.HetsOnly.PerIndividual.nostrict.CentralAsRevCompedToTs.SIGPROFILEREXTRACTORFORMAT.txt",header=T,sep="\t") # running for 3mers on sage
# need to combine all species ^^
  # repeat masked

finwhales=read.table("/Users/annabelbeichman/Documents/UW/WhaleProject/fin_whales/results/mutyper/downloaddate_20210111/mutyperResults_20210804_NOSTRICT_3mer/summaryTables/finWhale.All_Intervals_Counts.perPopulation.nostrict.NOTEONLYNEUTRALREGIONS.HIGHLYFILTERED.ALLFREQS.CentralAsRevCompedToTs.SIGPROFILEREXTRACTORFORMAT.txt",header=T,sep="\t")

# need to process 96 ints and convert 
# neutral regions

# humans? # need to run with 3mers if I want them or sum over 7mers

speciesToInclude=list(mice,bears,vaquita,finwhales,bulk_cetaceans)

# note it converts Mutation Types to Mutatation.Types so need to deal with that before I write out
# okay so next step is to combine:
combinedDF  <- Reduce(function(x,y) merge(x,y,by="Mutation.Types",all=TRUE) ,speciesToInclude)

# need to rename the column
colnames(combinedDF)[colnames(combinedDF) == 'Mutation.Types'] <- '`Mutation Types`'

head(combinedDF)

write.table(combinedDF,paste0(outdir,todaysdate,"MultispeciesMatrix.Mouse.Bear.Vaquita.FinWhale.BulkCetaceans.txt"),quote=F,sep="\t",row.names = F)
# need some way to detail what sorts of data went into this? I guess keeping this script 
