require(dplyr)
require(reshape2)
require(tidyr)
require(spgs)  # for reverse complements

############ convert ksfs to sigprofilerextractor format ##########
## restricting to signletons only for now #######
# contains multiple populations
## bears:
indir="/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210401_3mer_nostrict/mapped_to_brown_bear/summaryTables/"
ksfs = read.table(paste0(indir,"/All_Intervals_kSFS.ABC.NonABC.NotSeparatedYet.brown_bear.ALLFREQS.nostrict.txt"),header=T,sep="\t")
head(ksfs)
# so I want to subset to singletons only:

singletonsOnly <- ksfs %>%
  subset(sample_frequency==1) # singletons only

head(singletonsOnly)  
# convert to sigprofilerextractr format which looks like 
#Mutation Types	PD4199a	PD4005a	PD3851a	PD4116a	PD4086a	PD4194a	PD4248a	PD4120a	PD4198a	PD3904a	PD3945a	PD4107a	PD3905a	PD4192a	PD4109a	PD4103a	PD4115a	PD4085a	PD3890a	PD4006a	PD4088a
#A[C>A]A	58	74	31	128	31	28	64	165	58	122	243	210	94	48	122	112	228	61	110	198	34
#A[C>A]C	36	66	34	138	19	24	43	55	72	112	163	176	69	42	133	50	169	51	91	173	18
#A[C>A]G	13	12	9	15	8	6	7	20	12	13	24	33	11	16	12	7	21	7	33	3
#A[C>A]T	37	64	21	142	23	13	35	65	45	107	155	176	65	40	96	44	158	49	87	192	26
# and has rev comped mutations with central A

# separate ancestral and derived 3mers:
singletonsOnly$ancestral3mer <- substr(singletonsOnly$variable,1,3)
head(singletonsOnly)
singletonsOnly$derived3mer <- substr(singletonsOnly$variable,5,7)
head(singletonsOnly)

# get ancestral central BP: 
singletonsOnly$ancestralCentralBP <- substr(singletonsOnly$variable,2,2)
head(singletonsOnly)
singletonsOnly$derivedCentralBP <- substr(singletonsOnly$variable,6,6)
head(singletonsOnly)
# get the reverse complement of all of them (then going to put together)

singletonsOnly$ancestral3mer.RevComplement <- toupper(reverseComplement(singletonsOnly$ancestral3mer))

singletonsOnly$derived3mer.RevComplement <- toupper(reverseComplement(singletonsOnly$derived3mer))

# so to be in sigprofiler format, if the ancestral central BP is an "A" then it needs to be rev complemented so that central BP is a "T"
# if central bp is a C then it's already good to go for sigprofiler
# note here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6717374/
# 
singletonsOnly$`Mutation Types` <- ""

# so if the center bp is a "C" it is good to go (no rev comping needed)
singletonsOnly[singletonsOnly$ancestralCentralBP=="C",]$`Mutation Types` <-  paste0(substr(singletonsOnly[singletonsOnly$ancestralCentralBP=="C",]$variable,1,1),"[",substr(singletonsOnly[singletonsOnly$ancestralCentralBP=="C",]$variable,2,2),">",substr(singletonsOnly[singletonsOnly$ancestralCentralBP=="C",]$variable,6,6),"]",substr(singletonsOnly[singletonsOnly$ancestralCentralBP=="C",]$variable,3,3))

# but if the center bp is an "A" then you need to form the rev complement:
singletonsOnly[singletonsOnly$ancestralCentralBP=="A",]$`Mutation Types` <-  paste0(substr(singletonsOnly[singletonsOnly$ancestralCentralBP=="A",]$ancestral3mer.RevComplement,1,1),"[",substr(singletonsOnly[singletonsOnly$ancestralCentralBP=="A",]$ancestral3mer.RevComplement,2,2),">",substr(singletonsOnly[singletonsOnly$ancestralCentralBP=="A",]$derived3mer.RevComplement,2,2),"]",substr(singletonsOnly[singletonsOnly$ancestralCentralBP=="A",]$ancestral3mer.RevComplement,3,3))


#View(input)
head(singletonsOnly)



names(singletonsOnly)
singletonsOnly_wider <- pivot_wider(select(singletonsOnly,-c(variable)),id_cols = c(`Mutation Types`,population),names_from =  population,values_from = totalSites)
head(singletonsOnly_wider)                    


# 
write.table(singletonsOnly_wider,paste0(indir,"All_Intervals_kmerCountsForSigProfilerExtractor.SINGLETONSONLY.bears.txt"),row.names=F,quote=F,sep="\t")

####### also want to make a full SFS (?) version? might be innovative?? ##########
# separate ancestral and derived 3mers:
ksfs$ancestral3mer <- substr(ksfs$variable,1,3)
head(ksfs)
ksfs$derived3mer <- substr(ksfs$variable,5,7)
head(ksfs)

# get ancestral central BP: 
ksfs$ancestralCentralBP <- substr(ksfs$variable,2,2)
head(ksfs)
ksfs$derivedCentralBP <- substr(ksfs$variable,6,6)
head(ksfs)
# get the reverse complement of all of them (then going to put together)

ksfs$ancestral3mer.RevComplement <- toupper(reverseComplement(ksfs$ancestral3mer))

ksfs$derived3mer.RevComplement <- toupper(reverseComplement(ksfs$derived3mer))

# so to be in sigprofiler format, if the ancestral central BP is an "A" then it needs to be rev complemented so that central BP is a "T"
# if central bp is a C then it's already good to go for sigprofiler
# note here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6717374/
# 
ksfs$`Mutation Types` <- ""

# so if the center bp is a "C" it is good to go (no rev comping needed)
ksfs[ksfs$ancestralCentralBP=="C",]$`Mutation Types` <-  paste0(substr(ksfs[ksfs$ancestralCentralBP=="C",]$variable,1,1),"[",substr(ksfs[ksfs$ancestralCentralBP=="C",]$variable,2,2),">",substr(ksfs[ksfs$ancestralCentralBP=="C",]$variable,6,6),"]",substr(ksfs[ksfs$ancestralCentralBP=="C",]$variable,3,3))

# but if the center bp is an "A" then you need to form the rev complement:
ksfs[ksfs$ancestralCentralBP=="A",]$`Mutation Types` <-  paste0(substr(ksfs[ksfs$ancestralCentralBP=="A",]$ancestral3mer.RevComplement,1,1),"[",substr(ksfs[ksfs$ancestralCentralBP=="A",]$ancestral3mer.RevComplement,2,2),">",substr(ksfs[ksfs$ancestralCentralBP=="A",]$derived3mer.RevComplement,2,2),"]",substr(ksfs[ksfs$ancestralCentralBP=="A",]$ancestral3mer.RevComplement,3,3))


#View(input)
head(ksfs)



names(ksfs)
ksfs_wider <- pivot_wider(select(ksfs,-c(variable)),id_cols = c(`Mutation Types`,population),names_from =  c(population,sample_frequency),values_from = totalSites)
head(ksfs_wider)                    
write.table(ksfs_wider,paste0(indir,"All_Intervals_kmerCountsForSigProfilerExtractor.ALLFREQUENCIESSEPARATED.bears.txt"),row.names=F,quote=F,sep="\t")
