######## convert mutyper output to input format for sigprofiler matrix #######
require(dplyr)
require(reshape2)
require(tidyr)
require(spgs)  # for reverse complements

# example matrix is: /Users/annabelbeichman/Documents/UW/DNAShapeProject/scripts/DNA_Shape_Project/analyses/SigProfilerExtractor/sandbox/21BRCA

#Mutation Types	PD4199a	PD4005a	PD3851a	PD4116a	PD4086a	PD4194a	PD4248a	PD4120a	PD4198a	PD3904a	PD3945a	PD4107a	PD3905a	PD4192a	PD4109a	PD4103a	PD4115a	PD4085a	PD3890a	PD4006a	PD4088a
#A[C>A]A	58	74	31	128	31	28	64	165	58	122	243	210	94	48	122	112	228	61	110	198	34
#A[C>A]C	36	66	34	138	19	24	43	55	72	112	163	176	69	42	133	50	169	51	91	173	18
#A[C>A]G	13	12	9	15	8	6	7	20	12	13	24	33	11	16	12	7	21	7	33	3
#A[C>A]T	37	64	21	142	23	13	35	65	45	107	155	176	65	40	96	44	158	49	87	192	26
# this is one that I've processed (so not raw from mutyper)

# bulk cetaceans: don't need to add up (whole genome is in one vcf)
# but do need to transpose and combine 

# okay so the indir here contains lots of individuals so I need to get them all and tranpose them


indir="/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyperResults_20210804_3mer/mutyper_spectrum_files/" # mouse 3mers with repeat masking
outdir="/Users/annabelbeichman/Documents/UW/WhaleProject/vaquita/results/mutyper/bulk_cetaceans/mutyperResults_20210804_3mer/summaryTables"
dir.create(outdir)


filelist <- list.files(indir,full.names = T,pattern=".txt")
# need to make sure no species exist
allSpectra <- data.frame()
runningPopList <- list("") # list of species so far 
for(file in filelist){
  input <- read.table(file,header=T)
  # need to transpose it
  if(unique(input$sample) %in% runningPopList){
    print(paste0("duplicate sample? ",unique(input$sample)," look into it"))
    break
  }
  runningPopList <- c(runningPopList,unique(input$sample))
  input_longer <- pivot_longer(input,cols = -sample,names_to='variable')
  head(input_longer)
  colnames(input_longer) <- c("population",'variable',"totalSites")

  head(input_longer)
  allSpectra <- bind_rows(allSpectra,input_longer)
}
runningPopList
dim(allSpectra)
head(allSpectra)
# have to rev comp those with central -A- ancestral state to be ancestral T state
# have to fully rev comp both mutations 
# so A[A>G]C (AAC > AGC) is going to be come GTT > GCT so that the center position is T
# but don't want to change direction of mutation (ancestral is still ancestral, derived is still derived)
# example of rev comping:
#toupper(reverseComplement("AAC"))
#toupper(reverseComplement("AGC"))

# separate ancestral and derived 3mers:
allSpectra$ancestral3mer <- substr(allSpectra$variable,1,3)
head(allSpectra)
allSpectra$derived3mer <- substr(allSpectra$variable,5,7)
head(allSpectra)

# get ancestral central BP: 
allSpectra$ancestralCentralBP <- substr(allSpectra$variable,2,2)
head(allSpectra)
allSpectra$derivedCentralBP <- substr(allSpectra$variable,6,6)
head(allSpectra)
# get the reverse complement of all of them (then going to put together)

allSpectra$ancestral3mer.RevComplement <- toupper(reverseComplement(allSpectra$ancestral3mer))

allSpectra$derived3mer.RevComplement <- toupper(reverseComplement(allSpectra$derived3mer))

# so to be in sigprofiler format, if the ancestral central BP is an "A" then it needs to be rev complemented so that central BP is a "T"
# if central bp is a C then it's already good to go for sigprofiler
# note here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6717374/
# 
allSpectra$`Mutation Types` <- ""

# so if the center bp is a "C" it is good to go (no rev comping needed)
allSpectra[allSpectra$ancestralCentralBP=="C",]$`Mutation Types` <-  paste0(substr(allSpectra[allSpectra$ancestralCentralBP=="C",]$variable,1,1),"[",substr(allSpectra[allSpectra$ancestralCentralBP=="C",]$variable,2,2),">",substr(allSpectra[allSpectra$ancestralCentralBP=="C",]$variable,6,6),"]",substr(allSpectra[allSpectra$ancestralCentralBP=="C",]$variable,3,3))

# but if the center bp is an "A" then you need to form the rev complement:
allSpectra[allSpectra$ancestralCentralBP=="A",]$`Mutation Types` <-  paste0(substr(allSpectra[allSpectra$ancestralCentralBP=="A",]$ancestral3mer.RevComplement,1,1),"[",substr(allSpectra[allSpectra$ancestralCentralBP=="A",]$ancestral3mer.RevComplement,2,2),">",substr(allSpectra[allSpectra$ancestralCentralBP=="A",]$derived3mer.RevComplement,2,2),"]",substr(allSpectra[allSpectra$ancestralCentralBP=="A",]$ancestral3mer.RevComplement,3,3))


#View(input)
head(allSpectra)


head(input)
names(allSpectra)
allSpectra_wider <- pivot_wider(select(allSpectra,-c(variable)),id_cols = c(`Mutation Types`,population),names_from =  population,values_from = totalSites)
head(allSpectra_wider)                    


# 
write.table(allSpectra_wider,paste0(indir,"bulkCetaceans.All_Intervals_Counts.HetsOnly.PerIndividual.nostrict.CentralAsRevCompedToTs.SIGPROFILEREXTRACTORFORMAT.txt"),row.names=F,quote=F,sep="\t")
