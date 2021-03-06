######### get branch lengths from raxml tree (not time calibrated) #########
require(ape)
require(reshape2)
require(dplyr)
require(ggplot2)
require(ggtree)
raxmlMLTree <- read.tree("/Users/annabelbeichman/Documents/UW/DNAShapeProject/information/vertlife/Data_S3_globalRAxML_files/RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick") # this is newick tree from https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000494#sec030 

# okay this tree has bootstraps (whole integers around 50-100)
# but also has branch lengths. 
# the branch lengths is what I'm itnerested in.
# they say they just did raxml to get topology
# but I want substitutions. is that what this is???
# I think the branch lengths are what I want (are related to subtitutions from Uphams 31 gene alignment (which I also have: .phy file in indir.))
# SI data s3
# I don't think time calibrate.d think represents diverfgence/substitutions, but make sure!
raxmlMLTree

# use keep.tip
speciesList=c("Balaenoptera_acutorostrata" ,"Balaenoptera_musculus" ,"Balaenoptera_physalus" ,"Tursiops_truncatus" ,"Physeter_macrocephalus" ,"Neophocaena_phocaenoides" ,"Phocoena_phocoena" ,"Orcinus_orca" ,"Neophocaena_asiaeorientalis" ,"Monodon_monoceros" ,"Lagenorhynchus_obliquidens" ,"Globicephala_melas" ,"Delphinapterus_leucas" ,"Phocoena_sinus" ,"Ursus_maritimus" ,"Ursus_arctos" ,"Mus_musculus" ,"Mus_spretus" ,"Homo_sapiens","Pongo_pygmaeus","Pongo_abelii","Pan_troglodytes","Pan_paniscus","Gorilla_gorilla","Microcebus_murinus")
# 20211206 adding mouse lemur Microcebus murinus

speciesList2=c("Balaenoptera_physalus" ,"Phocoena_sinus" ,"Ursus_maritimus" ,"Ursus_arctos" ,"Mus_musculus" ,"Mus_spretus" ,"Homo_sapiens")

speciesList3=c("Balaenoptera_physalus" ,"Phocoena_sinus" ,"Ursus_maritimus" ,"Ursus_arctos" ,"Mus_musculus" ,"Mus_spretus","Homo_sapiens")

speciesList4=c("Balaenoptera_physalus" ,"Phocoena_sinus" ,"Ursus_maritimus" ,"Ursus_arctos" ,"Mus_musculus" ,"Mus_spretus","Homo_sapiens","Pan_paniscus") # add in human bonobo for another pair

speciesList5=c("Balaenoptera_physalus" ,"Phocoena_sinus" ,"Ursus_maritimus" ,"Ursus_arctos" ,"Mus_musculus" ,"Mus_spretus","Homo_sapiens","Pongo_pygmaeus","Pongo_abelii","Pan_troglodytes","Pan_paniscus","Gorilla_gorilla") # add in human bonobo for another pair

speciesList6=c("Balaenoptera_physalus" ,"Phocoena_sinus" ,"Ursus_maritimus" ,"Ursus_arctos" ,"Mus_musculus" ,"Mus_spretus","Homo_sapiens","Pongo_pygmaeus","Pongo_abelii","Pan_troglodytes","Pan_paniscus","Gorilla_gorilla","Microcebus_murinus") # added lemur

# no fin whale:
speciesList7=c("Phocoena_sinus" ,"Ursus_maritimus" ,"Ursus_arctos" ,"Mus_musculus" ,"Mus_spretus","Homo_sapiens","Pongo_pygmaeus","Pongo_abelii","Pan_troglodytes","Pan_paniscus","Gorilla_gorilla")

# for mantel test etc 
speciesList_PAG =speciesList5
########## plot subset of tips ##########
listOftips <- raxmlMLTree$tip.label

tipsToKeep <- listOftips[grepl(paste(speciesList, collapse="|"),listOftips)] # tips I want
tipsToKeep2 <- listOftips[grepl(paste(speciesList2, collapse="|"),listOftips)] # tips I want

tipsToKeep3 <- listOftips[grepl(paste(speciesList3, collapse="|"),listOftips)] # tips I want

tipsToKeep4 <- listOftips[grepl(paste(speciesList4, collapse="|"),listOftips)]
tipsToKeep5 <- listOftips[grepl(paste(speciesList5, collapse="|"),listOftips)]
tipsToKeep6 <- listOftips[grepl(paste(speciesList6, collapse="|"),listOftips)]
tipsToKeep7 <- listOftips[grepl(paste(speciesList7, collapse="|"),listOftips)]

tipsToKeep_PAG <- listOftips[grepl(paste(speciesList_PAG, collapse="|"),listOftips)]
outdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/raxml_tree/"

# okay this is good. maybe slightly fewer species but will work well and is direclty based on DNA!


png(paste0(outdir,"plotOfRaxmlTree.subsetofTips.WithTipLabels.Allwhalesincluded.png"),width = 15,height=5,units="in",res=300)
plot(keep.tip(raxmlMLTree,tipsToKeep),show.tip.label = T)
axisPhylo(side = 1, root.time = NULL, backward = FALSE)

dev.off()

png(paste0(outdir,"plotOfRaxmlTree.subsetofTips.WithoutTipLabels.Allwhalesincluded.png"),width = 7,height=5,units="in",res=300)
plot(keep.tip(raxmlMLTree,tipsToKeep),show.tip.label = F)
axisPhylo(side = 1, root.time = NULL, backward = FALSE)

dev.off()


png(paste0(outdir,"plotOfRaxmlTree.subsetofTips.WithTipLabels.png"),width = 12,height=5,units="in",res=300)
plot(keep.tip(raxmlMLTree,tipsToKeep2),show.tip.label = T)
axisPhylo(side = 1, root.time = NULL, backward = FALSE)

dev.off()

png(paste0(outdir,"plotOfRaxmlTree.subsetofTips.WithoutTipLabels.png"),width = 7,height=5,units="in",res=300)
plot(keep.tip(raxmlMLTree,tipsToKeep2),show.tip.label = F)
axisPhylo(side = 1, root.time = NULL, backward = FALSE)

dev.off()

# plot without humans

png(paste0(outdir,"plotOfRaxmlTree.subsetofTips.NoHumans.WithTipLabels.png"),width = 12,height=5,units="in",res=300)
plot(keep.tip(raxmlMLTree,tipsToKeep3),show.tip.label = T)
axisPhylo(side = 1, root.time = NULL, backward = FALSE)

dev.off()

png(paste0(outdir,"plotOfRaxmlTree.subsetofTips.NoHumans.WithoutTipLabels.png"),width = 7,height=5,units="in",res=300)
plot(keep.tip(raxmlMLTree,tipsToKeep3),show.tip.label = F)
axisPhylo(side = 1, root.time = NULL, backward = FALSE)

dev.off()

### add in apes:
png(paste0(outdir,"plotOfRaxmlTree.subsetofTips.HumansAndAllApes.WithTipLabels.png"),width = 14,height=5,units="in",res=300)
plot(keep.tip(raxmlMLTree,tipsToKeep5),show.tip.label =T)
axisPhylo(side = 1, root.time = NULL, backward = FALSE)

dev.off()

png(paste0(outdir,"plotOfRaxmlTree.subsetofTips.HumansAndAllApes.WithoutTipLabels.png"),width = 7,height=5,units="in",res=300)
plot(keep.tip(raxmlMLTree,tipsToKeep5),show.tip.label =F)
axisPhylo(side = 1, root.time = NULL, backward = FALSE)

dev.off()

### add in lemur:
png(paste0(outdir,"plotOfRaxmlTree.subsetofTips.HumansAndAllApes.Lemur.WithTipLabels.png"),width = 14,height=5,units="in",res=300)
plot(keep.tip(raxmlMLTree,tipsToKeep6),show.tip.label =T)
axisPhylo(side = 1, root.time = NULL, backward = FALSE)

dev.off()

png(paste0(outdir,"plotOfRaxmlTree.subsetofTips.HumansAndAllApes.Lemur.WithoutTipLabels.png"),width = 7,height=5,units="in",res=300)
plot(keep.tip(raxmlMLTree,tipsToKeep6),show.tip.label =F)
axisPhylo(side = 1, root.time = NULL, backward = FALSE)

dev.off()

# remove fin whale:
png(paste0(outdir,"plotOfRaxmlTree.subsetofTips.HumansAndAllApes.v.WithTipLabels.png"),width = 14,height=5,units="in",res=300)
plot(keep.tip(raxmlMLTree,tipsToKeep7),show.tip.label =T)
axisPhylo(side = 1, root.time = NULL, backward = FALSE)

dev.off()

png(paste0(outdir,"plotOfRaxmlTree.subsetofTips.HumansAndAllApes.NoFinWhale.WithoutTipLabels.png"),width = 7,height=5,units="in",res=300)
plot(keep.tip(raxmlMLTree,tipsToKeep7),show.tip.label =F)
axisPhylo(side = 1, root.time = NULL, backward = FALSE)

dev.off()


########## get distances #########


pairwiseDistances <- cophenetic.phylo(raxmlMLTree) # okay this gets pairwise distances between all pairs for all trees

pairwiseDistances_melt <- melt(pairwiseDistances)
head(pairwiseDistances_melt)
colnames(pairwiseDistances_melt) <- c("Sp1","Sp2","cophenetic_distance")
dim(pairwiseDistances_melt)

pairwiseDistances_melt$Sp1_species <- paste0(unlist(lapply(strsplit(as.character(pairwiseDistances_melt$Sp1),"_"),"[",1)),"_",unlist(lapply(strsplit(as.character(pairwiseDistances_melt$Sp1),"_"),"[",2)))

pairwiseDistances_melt$Sp1_family <- unlist(lapply(strsplit(as.character(pairwiseDistances_melt$Sp1),"_"),"[",3))

pairwiseDistances_melt$Sp1_order <- unlist(lapply(strsplit(as.character(pairwiseDistances_melt$Sp1),"_"),"[",4))
            

pairwiseDistances_melt$Sp2_species <- paste0(unlist(lapply(strsplit(as.character(pairwiseDistances_melt$Sp2),"_"),"[",1)),"_",unlist(lapply(strsplit(as.character(pairwiseDistances_melt$Sp2),"_"),"[",2)))

pairwiseDistances_melt$Sp2_family <- unlist(lapply(strsplit(as.character(pairwiseDistances_melt$Sp2),"_"),"[",3))

pairwiseDistances_melt$Sp2_order <- unlist(lapply(strsplit(as.character(pairwiseDistances_melt$Sp2),"_"),"[",4))

head(pairwiseDistances_melt)


######## okay then I want to.... restric tto my species 


write.table(pairwiseDistances_melt,paste0(outdir,"raxml.cophenetic.pairwisedistances.fromUphametal.ALLMAMMALS.txt"),quote=F,row.names = F,sep="\t")
# okay want to pull my species out 

##### want to restrict just to my species list #######
# speciesList defined above speciesList=c("Balaenoptera_acutorostrata" ,"Balaenoptera_musculus" ,"Balaenoptera_physalus" ,"Tursiops_truncatus" ,"Physeter_macrocephalus" ,"Neophocaena_phocaenoides" ,"Phocoena_phocoena" ,"Orcinus_orca" ,"Neophocaena_asiaeorientalis" ,"Monodon_monoceros" ,"Lagenorhynchus_obliquidens" ,"Globicephala_melas" ,"Delphinapterus_leucas" ,"Phocoena_sinus" ,"Ursus_maritimus" ,"Ursus_arctos" ,"Mus_musculus" ,"Mus_spretus" ,"Homo_sapiens")
# want to grep these 

pairwiseDistances_melt_mySpecies <- pairwiseDistances_melt[pairwiseDistances_melt$Sp1_species %in% speciesList & pairwiseDistances_melt$Sp2_species %in% speciesList,]

# NICE!! works!!
 # need to deal with reciprocal duplicates
# put in alphabetical order so I can get rid of dups:
pairwiseDistances_melt_mySpecies$comparisonLabel <- paste0(pmin(pairwiseDistances_melt_mySpecies$Sp1_species,pairwiseDistances_melt_mySpecies$Sp2_species),".",pmax(pairwiseDistances_melt_mySpecies$Sp1_species,pairwiseDistances_melt_mySpecies$Sp2_species))

# need to use pmin and pmax to get row by row min max. otherwise gives min max of whoel vector
# restrict to unique comp labels:
# get rid of self to self comparisons:
pairwiseDistances_melt_mySpecies_distinct <- pairwiseDistances_melt_mySpecies[pairwiseDistances_melt_mySpecies$Sp1!=pairwiseDistances_melt_mySpecies$Sp2,]
# and get rid of reciprocal dups (sp A - B and sp B- A)
pairwiseDistances_melt_mySpecies_distinct <- pairwiseDistances_melt_mySpecies_distinct %>%
  ungroup() %>% # adding ungroup as a safety measure on 20211201 not sure if ottally necessary though
  distinct(comparisonLabel,.keep_all=T)  # restrict to just distinct comparison labels. 
dim(pairwiseDistances_melt_mySpecies_distinct)
write.table(pairwiseDistances_melt_mySpecies_distinct,paste0(outdir,"raxml.cophenetic.pairwisedistances.fromUphametal.MYSPECIES.usethis.txt"),quote=F,row.names = F,sep="\t")

