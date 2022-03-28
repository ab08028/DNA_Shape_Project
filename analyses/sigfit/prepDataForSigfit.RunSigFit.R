############# gather all spectra and calculate rescaled mutation rates  ###########
require(tidyverse)
require(reshape2)
require(dplyr)
require(ggplot2)
require(tidyr)
require(compositions)  # for clr
require(broom) # for tidy
require(ggrepel)
require(ggfortify) # for pca autoplot
require(gridExtra)
require(sjPlot) # for grids of plots
require(ggpmisc) # for nicer plotting of linear models
require(ape)
require(sigfit)
########### NOTE FOR SIGFIT: 
# COUNTS, SIGNATURES, TARGETS MUST BE IN SAME ORDER. NOT SUFFICIENT TO HAVE SAME NAMES
# SO IF WORKING WITH TRIPLETS YOU MUST GET INTO SAME ORDER AS COSMIC IF YOU WANT TO COMPARE TO COSMIC 
niter=5000 #default is 2000 ; but want more 
model="multinomial" # (NMF)


epsilon=1 # not using epsilon for now; don't get confused # amount to add to 0 bins
separateCpGs="yes" # yes or no
todaysdate=format(Sys.Date(),"%Y%m%d")
flagOfAnalysis=paste0(todaysdate,".sigfitResults") # a label that will be in title of plots and outdir 
plotdir=paste0("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/sigfit/prepSpectraForSigfit/",flagOfAnalysis,"/")
dir.create(plotdir,showWarnings = F,recursive = T)

####### signatures ##########
agingSignatures <- read.table("/Users/annabelbeichman/jonsson.agingsignatures.tables9.convertedtoproportionsinexcel.REVCOMP.REORDEREDTOMATCHMYDATA.txt",header=T,sep="\t",row.names = "signature") # USE VERSION THAT HAS BEEN REVCOMPED 
# NEED TO GET THESE IN SAME ORDER AND NAMING CONVENTIONS AS MY DATA 

agingSignatures
#THESE MUST BE ORDERED IN SAME WAY AS COUNTS (this is a big deal)
if(sum(names(spectrum_1mer_forSigfit_downsampledCounts)!=names(agingSignatures))!=0){
  print("NAMES ARENT IN SAME ORDER BETWEEN SIGNATURES AND DATA -- THIS IS FATALLY BAD")
  # do this for all signatures and stuff.
  
}


########## test:adding a made up cpg signature ##########
CpGSignature = data.frame(row.names = "CpG.TpG signature",A.C=0,A.G=0,	A.T=0,	C.A=0	,C.G=0,	C.T=0,	C.T_CpG=1.0) # this signature is just CpG>T (100% of mutations from this signature)
# Can I use it to sop up CpGs? and then look more closely at other signatures?
agingPlusCpGSignatures <- bind_rows(CpGSignature,agingSignatures)

######### NOTE: If you use the cosmic signatures, need to normalize them away from genome targets (see github notes for sigfit about this!)
########## table of target locations ########

tableOfTargetsToInclude=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/spectra/tableOfSpectraToIncludeInComparisons.SpeciesWithMultipleIndividualsOnly.PlusApes.MaskedMouse.Vaquita.20211123.txt",header=T,sep="\t")  # UPDATED this file; it now has bedmasked vaquita targets and repmasked mouse targets

############## all projected spectra #############
# okay have all the spectra together already in:
# THESE have been projected down to k haploids:
#allProjectedSpectra <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/allKSFSes_and_projected_Spectra/All_SpectrumCounts_ProjectedDownTo.16.Haploids.txt",header=T)
allProjectedSpectra <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/allKSFSes_and_projected_Spectra/20211123_All_SpectrumCounts_ProjectedDownTo.10.Haploids.includesApes.RepMaskedMice.txt",header=T) # updated this file 20211123 to include repmasked mouse 
# this now includes bears_EUR and all apes (projectd down to 5 diploids, 10 haps)
speciesCodes=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/scripts/DNA_Shape_Project/analyses/PhylogeneticDistances/speciesCodes.txt",header=T) # for now just using AFR as human ; added apes
head(speciesCodes) # for getting phylo distances 



######### 1mers are dealt with in the if statement below this (dealing with CpGs) ########
# get ancestral kmers:
allProjectedSpectra$ancestral7mer <- substr(allProjectedSpectra$variable,1,7)
allProjectedSpectra$ancestral5mer <- substr(allProjectedSpectra$variable,2,6)
allProjectedSpectra$ancestral3mer <- substr(allProjectedSpectra$variable,3,5)

# get central mutation types:
allProjectedSpectra$mutation_7mer <- allProjectedSpectra$variable
allProjectedSpectra$mutation_5mer <- paste0(substr(allProjectedSpectra$variable,2,6),".",substr(allProjectedSpectra$variable,10,14))
allProjectedSpectra$mutation_3mer <- paste0(substr(allProjectedSpectra$variable,3,5),".",substr(allProjectedSpectra$variable,11,13))

# want to get ancestral 1mers for all data frames so I can plot things with central mutattion type separated later:

######## SEPARATE OUT  -- OPTIONAL ############
# just want to separate out C>Ts not other CpG types here (different in other script)
if(separateCpGs=="yes"){
  #### label CpGs and specify them as ancestral target as well
  print("separating out CpG sites!")
  allProjectedSpectra$mutation_1mer_CpGNotLabeled <- paste0(substr(allProjectedSpectra$variable,4,4),".",substr(allProjectedSpectra$variable,12,12))
  allProjectedSpectra$ancestral1mer_CpGNotLabeled <- paste0(substr(allProjectedSpectra$variable,4,4)) # adding this on 20220112 
  #### label CpGs and specify them as ancestral target as well
  allProjectedSpectra$CpGLabel <- ""
  # find all CpG sites:
  allProjectedSpectra[grepl("CG",allProjectedSpectra$ancestral3mer) & allProjectedSpectra$mutation_1mer_CpGNotLabeled=="C.T",]$CpGLabel <- "_CpG" # on 20220112 am changing this to just do C.T mtuations 
  allProjectedSpectra$mutation_1mer <- paste0(allProjectedSpectra$mutation_1mer_CpGNotLabeled,allProjectedSpectra$CpGLabel)
  # label ancestral state with CpG as well
  allProjectedSpectra$ancestral1mer <- allProjectedSpectra$ancestral1mer_CpGNotLabeled
  allProjectedSpectra[allProjectedSpectra$CpGLabel=="_CpG",]$ancestral1mer <- "C_CpG"
  
  allProjectedSpectra$mutation_1mer <- paste0(allProjectedSpectra$mutation_1mer_CpGNotLabeled,allProjectedSpectra$CpGLabel)
  
  allProjectedSpectra$ancestral1mer_CpGNotLabeled <- substr(allProjectedSpectra$variable,4,4) # not yet CpG identified 
  
  
} else if(separateCpGs=="no"){
  print("note that CpGs are not separated out!")
  allProjectedSpectra$mutation_1mer <- paste0(substr(allProjectedSpectra$variable,4,4),".",substr(allProjectedSpectra$variable,12,12)) # if you don't want to separate out CpGs then this is just the regular center bp
  # still get ancestral1mer (but don't specify if CpG or not)
  allProjectedSpectra$ancestral1mer <- substr(allProjectedSpectra$variable,4,4) # if you don't want to separate out CpGs then this is just the regular center bp
  
  
  
} else {
  print("invalid choice for separateCpGs")
  break
}



# get spectrum summed up per type
all7merSpectraOnly <- allProjectedSpectra %>%
  # adding groupbing by central 1mer but that doesn't actually change the sums since 7mers already group by those 
  # I just want to keep it in the dataframe 
  group_by(species,population,label,ancestral7mer,mutation_7mer,ancestral1mer,mutation_1mer) %>%
  summarise(total_mutations_projected = sum(totalSites_projectedDown_allFreqsSummed)) # this will be identical to original spectrum (7mer), just doing for consistency of formatting

all5merSpectraOnly <- allProjectedSpectra %>%
  group_by(species,population,label,ancestral5mer,mutation_5mer,ancestral1mer,mutation_1mer) %>%
  summarise(total_mutations_projected = sum(totalSites_projectedDown_allFreqsSummed)) 

all3merSpectraOnly <- allProjectedSpectra %>%
  group_by(species,population,label,ancestral3mer,mutation_3mer,ancestral1mer,mutation_1mer) %>%
  summarise(total_mutations_projected = sum(totalSites_projectedDown_allFreqsSummed)) 

all1merSpectraOnly <- allProjectedSpectra %>%
  group_by(species,population,label,ancestral1mer,mutation_1mer) %>%
  summarise(total_mutations_projected = sum(totalSites_projectedDown_allFreqsSummed)) 
# if you don't want to partition over CpGs you can use mutation_1mer_CpGNotLabeled instead
# okay NOTE HERE: "C" counts are nonCpG and then C_CpG are CpG ancestral targets


########## read in targets ##############
all7merTargets=data.frame()
all5merTargets=data.frame()
all3merTargets = data.frame()
all1merTargets = data.frame()

for(sample in tableOfTargetsToInclude$sample){
  # some 'samples' may contain multipops like mice and bears
  sampleinfo = tableOfTargetsToInclude[tableOfTargetsToInclude$sample==sample,]
  targetsfilename = paste0(sampleinfo$targetsindir,"/",sampleinfo$targetsFile)
  
  targets=read.table(targetsfilename,header=sampleinfo$targetsHeader) # assigns whether has header or not
  colnames(targets) <- c("target_7mer","countOverAllIntervals")
  
  # get central 5mer, 3mer and 1mer:
  targets <- targets %>%
    mutate(target_5mer=substr(target_7mer,2,6),target_3mer=substr(target_7mer,3,5),target_1mer= substr(target_7mer,4,4)) # get substrings 
  
  ######### relabel target_1mer with CpGs if you are separating out CpGs ###########
  if(separateCpGs=="yes"){
    targets$CpGLabel <- ""
    targets[grepl("CG",targets$target_3mer),]$CpGLabel <- "_CpG"
    targets$target_1mer_CpGsNotLabeled <- targets$target_1mer # update with this 
    targets$target_1mer <- paste0(targets$target_1mer_CpGsNotLabeled,targets$CpGLabel)
    
  } # don't need an else statement, just keep going with targets target$target_1mer as is.
  
  
  
  # TARGETS per 3mer type #
  targets_7mer <- targets %>% group_by(target_7mer) %>% summarise(total_target_count=sum(countOverAllIntervals)) # original targets file, just doing to be consistent but it's identical to original targets
  targets_3mer <- targets %>% group_by(target_3mer) %>% summarise(total_target_count=sum(countOverAllIntervals))
  targets_5mer <- targets %>% group_by(target_5mer) %>% summarise(total_target_count=sum(countOverAllIntervals))
  targets_1mer <- targets %>% group_by(target_1mer) %>% summarise(total_target_count=sum(countOverAllIntervals))
  targets_7mer$species <- sample
  targets_5mer$species <- sample
  targets_3mer$species <- sample
  targets_1mer$species <- sample
  # make all species targets:
  
  
  # combine  all targets:
  all7merTargets = bind_rows(all7merTargets,targets_7mer)
  all5merTargets = bind_rows(all5merTargets,targets_5mer)
  all3merTargets = bind_rows(all3merTargets,targets_3mer)
  all1merTargets = bind_rows(all1merTargets,targets_1mer)
  
  
}
# plot targets:
target1merPlot <- ggplot(all1merTargets,aes(y=species,x=total_target_count,fill=target_1mer))+
  geom_col(position="dodge")

target1merPlot2 <- all1merTargets %>%
  group_by(species) %>%
  mutate(target_proportion = total_target_count / sum(total_target_count)) %>%
  ungroup() %>%
  ggplot(.,aes(y=species,x=target_proportion,fill=target_1mer))+geom_col()+facet_wrap(~target_1mer,scales="free")+ggtitle("target proportions per genome")
target1merPlot2

ggsave(paste0(plotdir,"targetPlot.1mers.png"),target1merPlot,height=9,width=5)
ggsave(paste0(plotdir,"targetPlot.1mers.proportions.png"),target1merPlot2,height=9,width=12)

#ggplot(all1merSpectraOnly,aes(y=species,x=total_mutations_projected,fill=mutation_1mer))+
#  geom_col(position="dodge")
# complete so that missing mutation types are filled in 
# okay this is kind of sinister, depending on how it was grouped above this may or may not work (but doesn't throw an error. Frustrating!)
# so put in a check
# keep an eye on this in other code. Don't think it's made things go wrong before but be alert.
all7merSpectraOnly_filledin <- all7merSpectraOnly %>%
  ungroup() %>%
  complete(nesting(mutation_7mer,ancestral7mer,ancestral1mer,mutation_1mer),nesting(species,population,label),fill=list(total_mutations_projected=0))%>%
  mutate(mutation_label=mutation_7mer)

if(dim(all7merSpectraOnly_filledin)[1]==dim(all7merSpectraOnly)[1]){
  print("filling in didn't work!")
}
all5merSpectraOnly_filledin <- all5merSpectraOnly %>%
  ungroup() %>% # need to pre-ungroup ; got grouped up above
  complete(nesting(mutation_5mer,ancestral5mer,ancestral1mer,mutation_1mer),nesting(species,population,label),fill=list(total_mutations_projected=0))%>%
  mutate(mutation_label=mutation_5mer) # adding mutation label for use in functions

all3merSpectraOnly_filledin <- all3merSpectraOnly %>%
  ungroup() %>%
  complete(nesting(mutation_3mer,ancestral3mer,ancestral1mer,mutation_1mer),nesting(species,population,label),fill=list(total_mutations_projected=0))%>%
  mutate(mutation_label=mutation_3mer)

all1merSpectraOnly_filledin <- all1merSpectraOnly %>%
  ungroup() %>%
  complete(nesting(mutation_1mer,ancestral1mer),nesting(species,population,label),fill=list(total_mutations_projected=0)) %>%
  mutate(mutation_label=mutation_1mer)
# fill in should only be needed for 7mer (won't change dim of others)

# so when I was merging before all species were combined there was a bug where missing mtuation types wouldn't get targets filled in which led to them being NA downstream 
# okay now fill in missing mutation types : as 0s PRIOR TO MERGING WITH TARGEST
#  note that targets are the same within species (mouse1 and mouse2 have same targets; bear 1 and bear2 have targets)

#### to merge with ancestral targets, need apes species not to be 'apes' (but still want label to be apes_Species)
# so let's just update them here: # making the species equal to the popualtion (which is where I stored ape species. not ideal from a gneralizing point of view but works fine)
############### fill in ape species info
all7merSpectraOnly_filledin[all7merSpectraOnly_filledin$species=="apes",]$species <- all7merSpectraOnly_filledin[all7merSpectraOnly_filledin$species=="apes",]$population


all5merSpectraOnly_filledin[all5merSpectraOnly_filledin$species=="apes",]$species <- all5merSpectraOnly_filledin[all5merSpectraOnly_filledin$species=="apes",]$population


all3merSpectraOnly_filledin[all3merSpectraOnly_filledin$species=="apes",]$species <- all3merSpectraOnly_filledin[all3merSpectraOnly_filledin$species=="apes",]$population

all1merSpectraOnly_filledin[all1merSpectraOnly_filledin$species=="apes",]$species <- all1merSpectraOnly_filledin[all1merSpectraOnly_filledin$species=="apes",]$population

############## merge with targets ###############

all7merSpectra <- merge(all7merSpectraOnly_filledin, all7merTargets,by.x=c("species","ancestral7mer"),by.y=c("species","target_7mer"))

if(dim(all7merSpectra)[1] != dim(all7merSpectraOnly_filledin)[1]){
  print("something went wrong with merge!")
  
}


all5merSpectra <- merge(all5merSpectraOnly_filledin, all5merTargets,by.x=c("species","ancestral5mer"),by.y=c("species","target_5mer"))
all3merSpectra <- merge(all3merSpectraOnly_filledin, all3merTargets,by.x=c("species","ancestral3mer"),by.y=c("species","target_3mer"))
all1merSpectra <- merge(all1merSpectraOnly_filledin, all1merTargets,by.x=c("species","ancestral1mer"),by.y=c("species","target_1mer"))

################ multinomial downsample ######################


totalSegSites <- all1merSpectra %>%
  group_by(species,label) %>%
  summarise(totalSegSites=sum(total_mutations_projected))

multinomlogfile <- file(paste0(plotdir,"multinomialLogFile.txt"))

writeLines(paste0("downsampling to ",min(totalSegSites$totalSegSites)," sites (",totalSegSites[totalSegSites$totalSegSites==min(totalSegSites$totalSegSites),]$species,")"),paste0(plotdir,"multinomialDownsampling.info.txt"))
close(multinomlogfile)

# add multinomial variables to this function.
processSpectra <- function(spectradf){
  # get mutation rates and fraction of segregating sites: also adding a column in where I add +epsilon (1) to every count to deal with 0 entries ; taking out adding epsilon to target counts since I fixed the empty target bug
  spectradf <- spectradf %>%
    mutate(total_mutations_projected_plusEpsilon = (total_mutations_projected+epsilon),mutation_rate_projected = (total_mutations_projected/total_target_count), mutation_rate_projected_plusEpsilon=(total_mutations_projected_plusEpsilon/total_target_count)) %>%
    group_by(species,population,label) %>%
    mutate(mutation_rate_projected_normalized = (mutation_rate_projected/sum(mutation_rate_projected)),fractionOfSegregatingSites_projected=(total_mutations_projected/sum(total_mutations_projected)), mutation_rate_projected_normalized_plusEpsilon=(mutation_rate_projected_plusEpsilon/sum(mutation_rate_projected_plusEpsilon)),fractionOfSegregatingSites_projected_plusEpsilon=(total_mutations_projected_plusEpsilon/sum(total_mutations_projected_plusEpsilon))) # adding plus1 versions.
  
  # add multinomial sampling:
  # calculate min number of seg sites:
  totalSegSites <- spectradf %>%
    group_by(species,population,label) %>%
    summarise(totalSegSites=sum(total_mutations_projected))
  
  subsampleValue=min(totalSegSites$totalSegSites)
  
  subsampleSpecies=totalSegSites[totalSegSites$totalSegSites==min(totalSegSites$totalSegSites),]$species
  
  print(paste0("downsampling to ",subsampleValue," sites (",subsampleSpecies,")"))
  # carry out multinomial sampling:
  
  spectradf <- spectradf %>%
    group_by(species,label,population) %>%
    mutate(total_mutations_projected_multinom_downsampled=as.numeric(rmultinom(n=1,size=subsampleValue,prob=fractionOfSegregatingSites_projected))) %>%  # n=1 is how many times to sample (reps); size is how many sites (~230,000), probs is vector of probabilities aka the raw fraction (no +epsilon, no target correction)
    ungroup() %>% # ungroup just in case here 
    # now am adding epsilon and doing all transforms to the downsampled counts 
    mutate(total_mutations_projected_multinom_downsampled_plusEpsilon = (total_mutations_projected_multinom_downsampled+epsilon),mutation_rate_projected_multinom_downsampled = (total_mutations_projected_multinom_downsampled/total_target_count), mutation_rate_projected_multinom_downsampled_plusEpsilon=(total_mutations_projected_multinom_downsampled_plusEpsilon/total_target_count)) %>%
    group_by(species,population,label) %>%
    mutate(mutation_rate_projected_multinom_downsampled_normalized = (mutation_rate_projected_multinom_downsampled/sum(mutation_rate_projected_multinom_downsampled)), mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon=(mutation_rate_projected_multinom_downsampled_plusEpsilon/sum(mutation_rate_projected_multinom_downsampled_plusEpsilon)),fractionOfSegregatingSites_projected_multinom_downsampled_plusEpsilon=(total_mutations_projected_multinom_downsampled/sum(total_mutations_projected_multinom_downsampled))) # adding plus1 versions.
  # get mutation rates and fraction of segregating sites: also adding a column in where I add +epsilon (1) to every count to deal with 0 entries ; taking out adding epsilon to target counts since I fixed the empty target bug
  # adding write out.
  return(spectradf)
  
}

pivotSpectra_perVariable <- function(processedspectradf,variableToPivot) {
  processedspectradf_wider <-  pivot_wider(processedspectradf,id_cols=c(mutation_label),names_from=c("label"),values_from =all_of(variableToPivot ))
  processedspectradf_wider$variable <- variableToPivot
  # pivots on a variable 
  return(processedspectradf_wider)
}


################# process spectra #################
# for 1mers this doesn't matter, but for 3mers need to rev comp ala SPE for sigfit to compare to cosmic signatures!!! (pull code from ProjSpectra...script for writing out for SPE)
# also need to write out processed spectra and oppos here
# need to get opportunities (oppos) to work -- how do they get matched? is ">" essential? does it work for non triplets? 
spectrum_1mer_forSigfit = processSpectra(all1merSpectra)


# okay now want to format for sigfit:
# want sample column and then want mtuaiton labels as column names

######## get counts for sigfit: 
spectrum_1mer_forSigfit_downsampledCounts <- pivot_wider(spectrum_1mer_forSigfit[,c("label","total_mutations_projected_multinom_downsampled","mutation_label")],id_cols=c("label"),names_from = "mutation_label",values_from ="total_mutations_projected_multinom_downsampled" )

########## plot spectra #########
plot1merspectra <- ggplot(melt(spectrum_1mer_forSigfit_downsampledCounts),aes(y=label,x=value,fill=variable))+
  geom_col()+
  facet_wrap(~variable,scales="free")
ggsave(paste0(plotdir,"plot1merspectra.png"),plot1merspectra,height=7,width=12)
# change 'label' column name to 'sample' for sigfit format 
# rename C.T_CpG to CpG.TpG
spectrum_1mer_forSigfit_downsampledCounts <- spectrum_1mer_forSigfit_downsampledCounts %>%
  remove_rownames %>% # have to remove any preexisting rownames 
  column_to_rownames(var="label") # adds row names

write.table(spectrum_1mer_forSigfit_downsampledCounts,paste0(plotdir,"spectrum_1mer.forsigfit.downsampledCounts.txt"),row.names=T,sep="\t",quote=F)


# need unique because only care about 1 appearance of each target type: 
spectrum_1mer_forSigfit_opportunities <- pivot_wider(unique(spectrum_1mer_forSigfit[,c("label","total_target_count","mutation_label")]),id_cols=c("label"),names_from = "mutation_label",values_from ="total_target_count" )

# okay actually I want this to be for every mutation type 
# so C>A and C>T would have same oppos; C>T_CpG would have CpG oppos; A>T and A>C woudl have same oppos etc. 
# okay the https://github.com/kgori/sigfit/issues/58 they say to make opportunities normalize to 1 -- that's intriguing... let's try it! 
spectrum_1mer_forSigfit_opportunities <- spectrum_1mer_forSigfit_opportunities %>%
  remove_rownames %>% # have to remove any preexisting rownames 
  column_to_rownames(var="label") # adds row names




if(sum(names(spectrum_1mer_forSigfit_opportunities)!=names(spectrum_1mer_forSigfit_downsampledCounts))!=0){
  print("opportunities and counts are not in same order! this is fatally bad!")
  exit()
}

write.table(spectrum_1mer_forSigfit_opportunities,paste0(plotdir,"opportunities_1mer.forsigfit.noteThatTargetsAreRepeatedPerMutType.txt"),row.names=T,sep="\t",quote=F)

########### normalize opportunities -- this gets done internally with sigfit don't have to do it #############
#spectrum_1mer_forSigfit_opportunities_RowSums <- rowSums(spectrum_1mer_forSigfit_opportunities)
# normalize (note these are not GENOME FRACTIONS -- they are just normalized counts; will be the same for A>X and C>X because it's just the A, C counts. but have to add them up to get genome fractions -- weird. not sure about this.)
#spectrum_1mer_forSigfit_opportunities_normalized <- spectrum_1mer_forSigfit_opportunities %>%
#  mutate_all(~./spectrum_1mer_forSigfit_opportunities_RowSums)
# this is waht the sigfit dev wants but note that fractions DO NOT REPRESENT GENOME FRACTIONS. BECAUSE TARGETS ARE MULTI-COUNTED. I'M NOT SURE ABOUT THIS. maintains right relationships between them (e.g. A / C targets is still correct ratio ) ; but doesn't rep genome fraction. 
# does let you see that fin whales are weird. 



######## think about contextualizing these

# to try:
# figure out how to plot cosines per species and extract more info! 
# try sigfit with normalized oppos
# try triplets
# try other models 
############################# SIGFIT STUFF STARTS HERE ########################################


######### function to get cossine similarities in a dataframe instead of as plot titles ######
GetCosineSimilarity_allSamples <- function(mcmc_results,runlabel){
  counts=mcmc_results$data$counts_int
  samples=rownames(counts)
  reconstructions=retrieve_pars(mcmc_results,"reconstructions")$mean # get mean
  NSAMP=length(samples)
  allCos=data.frame()
  for(i in seq(1,NSAMP)){
    samplename=samples[i]
    cossimilarity=cosine_sim(counts[i,],reconstructions[i,]) # get cosine similarity (this is how sigfit does it in plotting code -- not very efficient, but oh well )
    cosdf = data.frame(species=samplename,cossine_similarity=cossimilarity)
    allCos=bind_rows(allCos,cosdf)
  }
  allCos$runlabel <- runlabel
  return(allCos)
}

########## for extracting de novo signatures ###############
run_plot_sigfit_EXTRACTING <- function(counts,opportunities=NULL,nsig_min,nsig_max,niter,model,wd,note){
  
  runlabel=paste("sigfit.extract_signatures",note,model,"iter",niter,sep=".")
  outdir=paste0(wd,runlabel,"/")
  
  dir.create(outdir,showWarnings = F,recursive=T)
  # first need to check that counts and oppos have same names in order (same species and same signatures)
  if(sum(rownames(counts)!=rownames(opportunities))!=0){
    print("different rownames for counts and opportunities")
    break
  }
  if(sum(names(counts)!=names(opportunities))!=0){
    print("different column names/order for counts and opportunities")
    break
  }
  sigfit_results <- extract_signatures(counts = counts,opportunities = opportunities,nsignatures=nsig_min:nsig_max,iter=niter,model=model)
  saveRDS(sigfit_results,paste0(outdir,runlabel,".rds"))
  # what is going on with R?
  #png(filename=paste0(outdir,runlabel,".gof.png"))
  #plot_gof(sigfit_results) # see what fits best -- still trying to save this plot 
  #dev.off()
  best=plot_gof(sigfit_results) # see what fits best -- still trying to save this plot 
  plot_all(sigfit_results[[best]],out_path = outdir,prefix=runlabel)
  # pull out and plot cossine similarities
  cossine_sims = GetCosineSimilarity_allSamples(sigfit_results[[best]],runlabel)
  write.table(cossine_sims,file=paste0(outdir,runlabel,".cossineSimilarities.txt"),row.names = F,quote=F,sep="\t")
  return(cossine_sims) # maybe return this so you can save it 
  # pull out and plot exposures (? maybe later)
  
}


####### for fitting a prior signatures ###############

run_plot_sigfit_FITTING <- function(counts,opportunities=NULL,signaturesToFit,signatureslabel,niter,model,wd,note){
  runlabel=paste("sigfit.fit",signatureslabel,"signatures",note,model,"iter",niter,sep=".")
  
  outdir=paste0(wd,runlabel,"/")
  dir.create(outdir,showWarnings = F,recursive=T)
  
  # if(sum(rownames(counts)!=rownames(opportunities)!=rownames(signaturesToFit))!=0){
  #   print("different rownames for counts and opportunities")
  #   break
  # }
  # if(sum(names(counts)!=names(opportunities)!=names())!=0){
  #   print("different column names/order for counts and opportunities")
  #   break
  # }
  # if(sum(names(counts)!=names(opportunities))!=0){
  #   print("different column names/order for counts and opportunities")
  #   break
  # }
  
  sigfit_results <- fit_signatures(counts,signatures = signaturesToFit,opportunities = opportunities,iter=niter,model=model)
  saveRDS(sigfit_results,paste0(outdir,runlabel,".rds"))
  plot_all(sigfit_results,out_path = outdir,prefix=runlabel)
  # pull out and plot cossine similarities
  cossine_sims = GetCosineSimilarity_allSamples(sigfit_results,runlabel)
  write.table(cossine_sims,file=paste0(outdir,runlabel,".cossineSimilarities.txt"),row.names = F,quote=F,sep="\t")
  return(cossine_sims) # maybe return this so you can save it 
  # pull out and plot exposures (? maybe later)
  
}


##### for fitting a priori signatures then extracting additional ones: 
run_plot_sigfit_FIT_EXT  <- function(counts,opportunities=NULL,signaturesToFit,signatureslabel,ndenovo, niter,model,wd,note) {
  runlabel=paste("sigfit.fit",signatureslabel,"signatures",note,model,"iter",niter,sep=".")
  
  outdir=paste0(wd,runlabel,"/")
  dir.create(outdir,showWarnings =F,recursive=T)
  
  sigfit_results <-  fit_extract_signatures(counts,signatures = signaturesToFit,num_extra_sigs = ndenovo,opportunities = opportunities,iter=niter,model=model)
  
  saveRDS(sigfit_results,paste0(outdir,runlabel,".rds"))
  plot_all(sigfit_results,out_path = outdir,prefix=runlabel)
  # pull out and plot cossine similarities
  cossine_sims = GetCosineSimilarity_allSamples(sigfit_results,runlabel)
  write.table(cossine_sims,file=paste0(outdir,runlabel,".cossineSimilarities.txt"),row.names = F,quote=F,sep="\t")
  return(cossine_sims) # maybe return this so you can save it 
  # pull out and plot exposures (? maybe later)
  
}

#################### OKAY TRY TO TEST THESE FUNCTIONS OUT ####################
############ YOU ARE HERE!  TO DO: add in checks to function about names/orders of rows!!! DO THIS!!
# eventually loop this better 
allCossineSims = data.frame()
####### extract sigs ###########
cossinesims1 = run_plot_sigfit_EXTRACTING(counts=spectrum_1mer_forSigfit_downsampledCounts,opportunities = spectrum_1mer_forSigfit_opportunities,nsig_min = 2,nsig_max = 8,model = "multinomial",niter = niter,wd=plotdir,note = "extract_sigs") 

allCossineSims = bind_rows(allCossineSims,cossinesims1)
###### fit aging sigs: ##########
cossinesims2 <- run_plot_sigfit_FITTING(counts=spectrum_1mer_forSigfit_downsampledCounts,signaturesToFit = agingSignatures,signatureslabel = "agingSignatures",opportunities = spectrum_1mer_forSigfit_opportunities,wd=plotdir,niter=niter,model=model,note="fittingAgingSignatures")
allCossineSims = bind_rows(allCossineSims,cossinesims2)


########## fit aging + CpG ########## 
cossinesims3 <- run_plot_sigfit_FITTING(counts=spectrum_1mer_forSigfit_downsampledCounts,signaturesToFit = agingPlusCpGSignatures,signatureslabel = "agingPlusCpGSignatures",opportunities = spectrum_1mer_forSigfit_opportunities,wd=plotdir,niter=niter,model=model,note="fittingAgingSignatures")
allCossineSims = bind_rows(allCossineSims,cossinesims3)

######## FIT+Extract: aging +1 denovo #######
cossinesims4 <- run_plot_sigfit_FIT_EXT(counts=spectrum_1mer_forSigfit_downsampledCounts,opportunities =spectrum_1mer_forSigfit_opportunities,signaturesToFit = agingSignatures,signatureslabel = "agingSignatures",ndenovo = 1,niter = niter,model = model,note = "agingPlus1DeNovo" ,wd=plotdir)
allCossineSims = bind_rows(allCossineSims,cossinesims4)

######## FIT+Extract: aging +2 denovo #######
cossinesims5 <- run_plot_sigfit_FIT_EXT(counts=spectrum_1mer_forSigfit_downsampledCounts,opportunities =spectrum_1mer_forSigfit_opportunities,signaturesToFit = agingSignatures,signatureslabel = "agingSignatures",ndenovo = 2,niter = niter,model = model,note = "agingPlus2DeNovo",wd=plotdir )
allCossineSims = bind_rows(allCossineSims,cossinesims5)

####### FIT+ extract: aging+Cpg+1 de novo #############
cossinesims6 <- run_plot_sigfit_FIT_EXT(counts=spectrum_1mer_forSigfit_downsampledCounts,opportunities =spectrum_1mer_forSigfit_opportunities,signaturesToFit = agingPlusCpGSignatures,signatureslabel = "agingPlusCpGSignatures",ndenovo = 1,niter = niter,model = model,note = "agingPlusCpGPlus1DeNovo",wd=plotdir )
allCossineSims = bind_rows(allCossineSims,cossinesims6)



############ try without correcting for opportunities ############
########## fit aging + CpG ########## 
cossinesims7 <- run_plot_sigfit_FITTING(counts=spectrum_1mer_forSigfit_downsampledCounts,signaturesToFit = agingPlusCpGSignatures,signatureslabel = "agingPlusCpGSignatures",opportunities = NULL,wd=plotdir,niter=niter,model=model,note="fittingAgingSignatures_NO_OPPORTUNITYCORRECTION")
allCossineSims = bind_rows(allCossineSims,cossinesims7)

########## plot all cossime similarities ########
write.table(allCossineSims,paste0(plotdir,"AllCossineSims.txt"),row.names = F,quote=F,sep="\t")



allCossineSims$runlabel <- factor(allCossineSims$runlabel,levels=c("sigfit.extract_signatures.extract_sigs.multinomial.iter.5000" , "sigfit.fit.agingSignatures.signatures.fittingAgingSignatures.multinomial.iter.5000" ,"sigfit.fit.agingPlusCpGSignatures.signatures.fittingAgingSignatures.multinomial.iter.5000","sigfit.fit.agingPlusCpGSignatures.signatures.fittingAgingSignatures_NO_OPPORTUNITYCORRECTION.multinomial.iter.5000", "sigfit.fit.agingSignatures.signatures.agingPlus1DeNovo.multinomial.iter.5000"  , "sigfit.fit.agingSignatures.signatures.agingPlus2DeNovo.multinomial.iter.5000" ,"sigfit.fit.agingPlusCpGSignatures.signatures.agingPlusCpGPlus1DeNovo.multinomial.iter.5000"   ))

cossineSimPlot <- ggplot(allCossineSims,aes(y=species,x=cossine_similarity,color=runlabel))+
  geom_point()
cossineSimPlot
ggsave(paste0(plotdir,"cossineSimilarityPlot.comparingRuns.png"),cossineSimPlot,height=7,width=18)

cossineSimPlot2 <- ggplot(allCossineSims,aes(x=cossine_similarity,y=reorder(runlabel,desc(runlabel))))+
  geom_boxplot()+
  ylab("")
cossineSimPlot2 
ggsave(paste0(plotdir,"cossineSimilarityPlot.comparingRuns2.png"),cossineSimPlot2,height=7,width=15)
  
cossineSimPlot3 <- ggplot(allCossineSims,aes(y=runlabel,x=cossine_similarity))+
geom_col()+
  facet_wrap(~species)
cossineSimPlot3


#



#

############ OLD STUFF
########### FIT JUST AGING SIGNATURES ###########
test_fit1_aging <- fit_signatures(spectrum_1mer_forSigfit_downsampledCounts,signatures = agingSignatures,opportunities = spectrum_1mer_forSigfit_opportunities,iter=niter,model=model)

plot_all(test_fit1_aging,out_path = plotdir,prefix="testing_fittingAgingSignatures")
# okay so this gets to about .9ish
cosine_similarities_test_fit1_aging <- GetCosineSimilarity_allSamples(test_fit1_aging,"fit aging signatures")

allCosineSimilarities
############## FIT AGING + 1 DENOVO #########
test_fit1_aging_plusdenovo <- fit_extract_signatures(spectrum_1mer_forSigfit_downsampledCounts,signatures = agingSignatures,num_extra_sigs = 1,opportunities = spectrum_1mer_forSigfit_opportunities,iter=niter,model=model)


plot_all(test_fit1_aging_plusdenovo,out_path = plotdir,prefix="testing_fittingAgingSignatures_plusOneDeNovo")
cosine_similarities_test_fit1_aging_plusdenovo <- GetCosineSimilarity_allSamples(test_fit1_aging_plusdenovo,"fit aging signatures + 1 denovo")

# TO DO: make a combo function that does all these things (fit, plot etc. ) -- need one function for fitting; one for extracting (with a range); and one for fit-extractign 
# need to figure out a way to save that gof plot nicely (later)

############## FIT AGING + 2DENOVO #########
test_fit1_aging_plusdenovo2 <- fit_extract_signatures(spectrum_1mer_forSigfit_downsampledCounts,signatures = agingSignatures,num_extra_sigs = 2,opportunities = spectrum_1mer_forSigfit_opportunities,iter=niter,model=model)

plot_all(test_fit1_aging_plusdenovo2,out_path = plotdir,prefix="testing_fittingAgingSignatures_plusTwoDeNovo")

cosine_similarities_test_fit1_aging_plusdenovo2 <- GetCosineSimilarity_allSamples(test_fit1_aging_plusdenovo2,"fit aging signatures + 2 denovo")
# idea: can I make a CpG signature (all 0 except CpGs?)
# make a bgc signature?
# make sure C>T counts with CpGs is good 
# try to control for sbs1?

########### Try adding a CpG Signature ############

# fit with no extra signatures and see if fit improves. 
names(agingPlusCpGSignatures)==names(spectrum_1mer_forSigfit_downsampledCounts) #Make sure this is all true 

test_fit_agingPlusCpGOnly <- fit_signatures(spectrum_1mer_forSigfit_downsampledCounts,signatures = agingPlusCpGSignatures,opportunities = spectrum_1mer_forSigfit_opportunities,iter=niter,model=model)

plot_all(test_fit_agingPlusCpGOnly,out_path = plotdir,prefix="testing_fittingAgingSignaturesPlusCpGSignature")
# okay so this gets to about .9ish
cosine_similarities_test_fit_agingPlusCpGOnly <- GetCosineSimilarity_allSamples(test_fit_agingPlusCpGOnly,"fit aging signatures + CpG signature")

############ try without correcting for opportunities (so I can understand why CpGs are such a big fraction ) ############
test_fit_agingPlusCpGOnly_NoCorrection <- fit_signatures(spectrum_1mer_forSigfit_downsampledCounts,signatures = agingPlusCpGSignatures,iter=niter,model=model)
plot_all(test_fit_agingPlusCpGOnly_NoCorrection,out_path = plotdir,prefix="testing_fittingAgingSignaturesPlusCpGSignature_NOOPPORTUNITYCORRECTION")
cosine_similarities_test_fit_agingPlusCpGOnly_NoCorrection <- GetCosineSimilarity_allSamples(test_fit_agingPlusCpGOnly_NoCorrection,"fit aging signatures + CpG signature (no opportunity correction)")
exposures_test_fit_agingPlusCpGOnly_NoCorrection <- retrieve_pars(test_fit_agingPlusCpGOnly_NoCorrection,"exposures")$mean # signatures and exposures get you to same place
exposures_test_fit_agingPlusCpGOnly_NoCorrection$sample <- rownames(exposures_test_fit_agingPlusCpGOnly_NoCorrection)
 
exposureplot_test_fit_agingPlusCpGOnly_NoCorrection <- ggplot(melt(exposures_test_fit_agingPlusCpGOnly_NoCorrection),aes(y=sample,x=value,color=variable,fill=variable))+
  geom_col()+
  ggtitle("exposures")
exposureplot_test_fit_agingPlusCpGOnly_NoCorrection
########### OKAY WANT TO MAKE A FUNCTION THAT MAKES ALL THESE PLOTS AND LABELS THEM #########
############ Try with non-downsampled data ##############
######## get counts for sigfit: 
spectrum_1mer_forSigfit_actualCounts <- pivot_wider(spectrum_1mer_forSigfit[,c("label","total_mutations_projected","mutation_label")],id_cols=c("label"),names_from = "mutation_label",values_from ="total_mutations_projected" )

spectrum_1mer_forSigfit_actualCounts <- spectrum_1mer_forSigfit_actualCounts %>%
  remove_rownames %>% # have to remove any preexisting rownames 
  column_to_rownames(var="label") # adds row names
# opportunities is the same
# I'm wondering if having different exposures based on opportunity correction is bad since I've downsampled to same number of snps. Hm. Not sure let's see.
test_extract_actualCounts <- extract_signatures(spectrum_1mer_forSigfit_actualCounts,nsignatures = 2:5 ,opportunities = spectrum_1mer_forSigfit_opportunities,iter=niter,model=model)

plot_gof(test_extract_actualCounts)
best_actualCounts <- plot_gof(test_extract_actualCounts)
plot_all(test_extract_actualCounts[[best_actualCounts]],out_path = plotdir,prefix="test_extract_actualCounts_notDownsampled")
# testing with other signatures 


########### try with non-downsampled Data and adding CpG signature and de novo (no aging) ######

test_fit_CpG_PlusDenovo <- fit_extract_signatures(spectrum_1mer_forSigfit_downsampledCounts,signatures = CpGSignature,num_extra_sigs = 3,opportunities = spectrum_1mer_forSigfit_opportunities,iter=niter,model=model)
plot_all(test_fit_CpG_PlusDenovo,out_path = plotdir,prefix="test_fit_CpG_PlusDenovo")

######################## Make some better summary plots #############
# want to summarise cosine per species
# want to summarise signatures in a grid of some sort 

# okay so the sigfit code loops this -- what a pain



############### try with triplets ################
# do this to better understand/compare to cosmic. 
