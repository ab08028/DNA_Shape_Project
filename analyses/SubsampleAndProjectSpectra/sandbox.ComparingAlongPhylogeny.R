require(ape)
require(castor)
require(phylosignal) 
require(ade4) # for mantel test
spectra <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/comparing_spectrum_and_phylo_divergence_WITHPROJECTIONS/20211214.plots.projectedCounts.epsilon.1.sepCpGs.no.addingApes_maskedMouseVaq.addingMultinomialDownsampling.MorePlotsForPAG/AllProcessedSpectra.AllConditions.AllVariables.AllSpecies.txt",header=T,sep='\t')

speciesCodes=read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/scripts/DNA_Shape_Project/analyses/PhylogeneticDistances/speciesCodes.txt",header=T) 

head(spectra)
# restrict to just multinom downsampled variable and clr transform and just 1mer spectrum:
spectra_1mer_clr <- spectra[spectra$variable=="mutation_rate_projected_multinom_downsampled_normalized_plusEpsilon" & spectra$transform_id=="clr" & spectra$id=="spectra_1mer",]
head(spectra_1mer_clr)

spectra_1mer_clr_melt <- melt(spectra_1mer_clr,variable.name = "species_label")
spectra_1mer_clr_melt <- merge(spectra_1mer_clr_melt,speciesCodes[,c("code","species")],by.x="species_label",by.y="code") # this adds it in to match the tree

timeTree=ape::read.tree("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/timeTreePhyloTree/listOfSpecies.ForTimeTree.20211210(1).nwk")


# make a vector of just one mutation type:
# and name it same way as tree

speciesToExclude = c("bears_EUR","fin_whale_ENP") # only want one per species 

spectra_1mer_clr_A_G <- spectra_1mer_clr_melt[spectra_1mer_clr_melt$mutation_label=="A.G" & !spectra_1mer_clr_melt$species_label %in% speciesToExclude,]$value

names(spectra_1mer_clr_A_G) <-spectra_1mer_clr_melt[spectra_1mer_clr_melt$mutation_label=="A.G" & !spectra_1mer_clr_melt$species_label %in% speciesToExclude,]$species

spectra_1mer_clr_A_G


spectra_1mer_clr_C_T <- spectra_1mer_clr_melt[spectra_1mer_clr_melt$mutation_label=="C.T" & !spectra_1mer_clr_melt$species_label %in% speciesToExclude,]$value

names(spectra_1mer_clr_C_T) <-spectra_1mer_clr_melt[spectra_1mer_clr_melt$mutation_label=="C.T" & !spectra_1mer_clr_melt$species_label %in% speciesToExclude,]$species

spectra_1mer_clr_C_T

all_ACF_df <- data.frame()
for(mutation_type in c("C.T","C.G","C.A","A.C","A.G","A.T")){
  # get just the one mut type 
  vectorOfMutationType <- spectra_1mer_clr_melt[spectra_1mer_clr_melt$mutation_label==mutation_type & !spectra_1mer_clr_melt$species_label %in% speciesToExclude,]$value
  # add names
  names(vectorOfMutationType) <-spectra_1mer_clr_melt[spectra_1mer_clr_melt$mutation_label==mutation_type & !spectra_1mer_clr_melt$species_label %in% speciesToExclude,]$species
  # could do pic or osmethign esle
  ACF = get_trait_acf(timeTree, vectorOfMutationType, Npairs=1e7, Nbins=4) # seems like Npairs is overkill 
  ACF_df <- data.frame(distance=ACF$phylodistances,autocor=ACF$autocorrelations,label=mutation_type)
  all_ACF_df <- bind_rows(all_ACF_df,ACF_df)

}

ggplot(all_ACF_df,aes(x=distance,y=autocor,color=label))+
  geom_line(size=2) 


# add a simulation: 

#+
 # geom_line()

# spectra_1mer_clr_A_G.pic <- pic(spectra_1mer_clr_A_G,timeTree)  # what does this even mean/do??!?!


# simulate continuous trait evolution on the tree on my own time tree
set.seed(123)
tip_states_simulated = simulate_bm_model(timeTree, diffusivity=1,Nsimulations = 1000)$tip_states
# brownian motion 

# calculate autocorrelation function
ACF_simulated = get_trait_acf(timeTree, tip_states_simulated, Npairs=1e7, Nbins=4)

ACF_simulated_df <- data.frame(distance=ACF_simulated$phylodistances,autocor=ACF_simulated$autocorrelations)

ggplot(ACF_AG_df,aes(x=distance,y=autocor))+
  geom_line(aes(color="A.G"))+
  geom_line(data=ACF_CT_df,aes(x=distance,y=autocor,color="C.T"))+
  geom_line(data=ACF_simulated_df,aes(x=distance,y=autocor,color="simulated trait diffusing across phylo"))


# plot ACF (autocorrelation vs phylogenetic distance)
#plot(ACF$phylodistances, ACF$autocorrelations, type="l", xlab="distance", ylab="ACF")
# still needs singular values #######


################ try with just one mutation type at a time #############
distances[distances$transform_label=="clr"]


############# try mantel test ##############

distances <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/comparing_spectrum_and_phylo_divergence_WITHPROJECTIONS/20211212.plots.projectedCounts.epsilon.1.sepCpGs.no.addingApes_maskedMouseVaq.addingMultinomialDownsampling.MorePlotsForPAG/all_distances_all_conditions_phylo.txt",header=T)

head(distances)
kmer="spectrum_1mer"
# start without downsampled version 
variable="mutation_rate_projected_normalized_plusEpsilon"
spectrum_dist_mat_1mer <- pivot_wider(distances[distances$id=="spectra_1mer" & distances$transform_label=="clr" & distances$variable==variable,c("item1","item2","spectrum_distance")],id_cols = "item1",names_from = "item2",values_from = spectrum_distance)

# fill in

tree_dist_mat_TIME <- cophenetic.phylo(timeTree)

# made this in other script ; need to code nicely

# uses lower triangles (so okay to have upper tri filled in by maybe not necessary?)
spectra_1mer_spectrumDistMat_ForMantel

tree_dist_mat_TIME_reordered <- as.matrix(tree_dist_mat_TIME)[rownames(as.matrix(spectra_1mer_spectrumDistMat_ForMantel)),colnames(as.matrix(spectra_1mer_spectrumDistMat_ForMantel))] # reorder to be in same order at 1mer spectrum 

mantel.test(as.matrix(spectra_1mer_spectrumDistMat_ForMantel),tree_dist_mat_TIME_reordered,nperm = 1000)


spectra_5mer_spectrumDistMat_ForMantel

mantel.test(as.matrix(spectra_5mer_spectrumDistMat_ForMantel),tree_dist_mat_TIME_reordered,nperm = 1000,graph=T)

spectra_7mer_spectrumDistMat_ForMantel
# this is from ape:
mt7 <- mantel.test(as.matrix(spectra_7mer_spectrumDistMat_ForMantel),tree_dist_mat_TIME_reordered,nperm = 999,graph=T)
mt7

# this is from ade4:
mt7 <- mantel.rtest(as.dist(spectra_7mer_spectrumDistMat_ForMantel),as.dist(tree_dist_mat_TIME_reordered),nrepet = 999)
mt7

plot(mt7)


mt5 <- mantel.rtest(as.dist(spectra_5mer_spectrumDistMat_ForMantel),as.dist(tree_dist_mat_TIME_reordered),nrepet = 999,)
mt5

plot(mt5)

# mantel coefficient 

####### vegan uses pearson #######
# and mantel to figure out significance 
require(vegan)
# https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/mantel
mt5_vegan <- mantel(as.dist(tree_dist_mat_TIME_reordered),as.dist(spectra_5mer_spectrumDistMat_ForMantel),permutations = 999,method="pearson") # careful with x and y? though does it  matter?
mt5_vegan
# okay I think vegan is way to go 

# maybe if time: 