require(tidyverse)
require(reshape2)
require(spgs)  # for reverse complements
require(sigfit)
# load in your data
# SAMPLE: df=data.frame(sample="bear","AAA>ACA"=50,"AAA>AGA"=80,"TCC.TTC"=85)
#colnames(select(df,-sample))

# reshape data from wide to long: 
df_melt = melt(df)
df_melt
# set a new variable for sigfit that you can modify 
df_melt$variable_forSigfit <- as.character(df_melt$variable)
df_melt$ancestral1mer <- substr(df_melt$variable,2,2) # 
df_melt$derived1mer <- substr(df_melt$variable,6,6) # 
df_melt$ancestral3mer <- substr(df_melt$variable,1,3) # 
df_melt$derived3mer <- substr(df_melt$variable,5,7) # 

# for any mutation type with ancestral1mer == "A", reverse complement the mutation 
df_melt[df_melt$ancestral1mer=="A",]$variable_forSigfit <- paste0(toupper(reverseComplement(df_melt[df_melt$ancestral1mer=="A",]$ancestral3mer)),".",toupper(reverseComplement(df_melt[df_melt$ancestral1mer=="A",]$derived3mer)))
df_melt

head(df_melt)
# at this stage View() the df and make sure everything looks ok!!!

# make it wide with mutation types becoming col names
df_wide <- pivot_wider(df_melt,id_cols=c("sample"),names_from = "variable_forSigfit",values_from ="value" )

# convert . to > in col names:
colnames(df_wide) <- gsub("\\.",">",colnames(df_wide))

head(df_wide)
# change the 'sample' column to rownames (don't want it in final dataframe as its own column)

rownames(df_wide) <- df_wide$sample
# make it in same order as cosmic database in sigfit
data(cosmic_signatures_v3.2) # load cosmic data

df_readyforsigfit <- df_wide[colnames(cosmic_signatures_v3.2)] # this selects jsut columns that are the colnames of the cosmic signature, and *critically* reorders them to match the cosmic dataframe

# ** it is critical that your data and the cosmic signatures are in the exact same column order!!
# check:
colnames(df_readyforsigfit) == colnames(cosmic_signatures_v3.2) # should be true for all 
