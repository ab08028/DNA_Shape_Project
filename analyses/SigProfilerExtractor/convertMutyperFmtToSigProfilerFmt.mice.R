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
indir="/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/wild_mouse_data/mutyperResults_20210706_NOSTRICT_3mer_REPEATMASKED/summaryTables/" # mouse 3mers with repeat masking
input <- read.table(paste0(indir,"Mouse3merSpectrum.SummedUpOverAllChrs.NOSTRICT.RepeatMasked.txt"),header=T)
head(input)


# have to rev comp those with central -A- ancestral state to be ancestral T state
# have to fully rev comp both mutations 
# so A[A>G]C (AAC > AGC) is going to be come GTT > GCT so that the center position is T
# but don't want to change direction of mutation (ancestral is still ancestral, derived is still derived)
# example of rev comping:
toupper(reverseComplement("AAC"))
toupper(reverseComplement("AGC"))

# separate ancestral and derived 3mers:
input$ancestral3mer <- substr(input$variable,1,3)
head(input)
input$derived3mer <- substr(input$variable,5,7)
head(input)

# get ancestral central BP: 
input$ancestralCentralBP <- substr(input$variable,2,2)
head(input)
input$derivedCentralBP <- substr(input$variable,6,6)
head(input)
# get the reverse complement of all of them (then going to put together)

input$ancestral3mer.RevComplement <- toupper(reverseComplement(input$ancestral3mer))

input$derived3mer.RevComplement <- toupper(reverseComplement(input$derived3mer))

# so to be in sigprofiler format, if the ancestral central BP is an "A" then it needs to be rev complemented so that central BP is a "T"
# if central bp is a C then it's already good to go for sigprofiler
# note here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6717374/
# 
input$`Mutation Types` <- ""

# so if the center bp is a "C" it is good to go (no rev comping needed)
input[input$ancestralCentralBP=="C",]$`Mutation Types` <-  paste0(substr(input[input$ancestralCentralBP=="C",]$variable,1,1),"[",substr(input[input$ancestralCentralBP=="C",]$variable,2,2),">",substr(input[input$ancestralCentralBP=="C",]$variable,6,6),"]",substr(input[input$ancestralCentralBP=="C",]$variable,3,3))

# but if the center bp is an "A" then you need to form the rev complement:
input[input$ancestralCentralBP=="A",]$`Mutation Types` <-  paste0(substr(input[input$ancestralCentralBP=="A",]$ancestral3mer.RevComplement,1,1),"[",substr(input[input$ancestralCentralBP=="A",]$ancestral3mer.RevComplement,2,2),">",substr(input[input$ancestralCentralBP=="A",]$derived3mer.RevComplement,2,2),"]",substr(input[input$ancestralCentralBP=="A",]$ancestral3mer.RevComplement,3,3))


#View(input)
head(input)
#input %>%
#  group_by(`Mutation Types`,population) %>%
#  summarise(sum(totalSitesAllIntervals)) %>%
#  View()

#input %>%
#  group_by(variable,population) %>%
#  summarise(sum(totalSitesAllIntervals)) 
# this contains multiple populations and I need to fmt it like ^^
# this isn't per individual (should it be? maybe / maybe not. should )


head(input)
input_wider <- pivot_wider(select(input,-c(variable,fractionOfAllSegSites,fractionOfCentralMutationType)),id_cols = c(`Mutation Types`,population),names_from =  population,values_from = totalSitesAllIntervals)
head(input_wider)                    


# 
write.table(input_wider,paste0(indir,"mouse.All_Intervals_Counts.perPopulation.ALLFREQS.nostrict.repeatmasked.CentralAsRevCompedToTs.SIGPROFILEREXTRACTORFORMAT.txt"),row.names=F,quote=F,sep="\t")
