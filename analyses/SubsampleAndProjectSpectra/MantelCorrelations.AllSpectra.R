############ New script. 
# what do I need
# 1, 3, 5, 7 mer distance matrices (for a few variables)
# time tree; raxml tree; upham time tree (?) that one is a bit whack so maybe let's skip that one? 
# upham tree is kind of weird. Let's go with time tree

require(ape)
require(vegan) # for mantel test 
############ read in trees #############
timeTree=ape::read.tree("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/timeTreePhyloTree/listOfSpecies.ForTimeTree.20211210(1).nwk")

raxmlTree=ape::read.tree("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/raxml_tree/RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick")

# want to get rid of the family/order info at end of tip. just keep species : 
raxmlTree_renamedTips <- raxmlTree
raxmlTree_renamedTips$tip.label <- paste0(unlist(lapply(strsplit(raxmlTree_renamedTips$tip.label,"_"),"[",1)),"_",unlist(lapply(strsplit(raxmlTree_renamedTips$tip.label,"_"),"[",2)))

speciesList=c("Balaenoptera_physalus" ,"Phocoena_sinus" ,"Ursus_maritimus" ,"Ursus_arctos" ,"Mus_musculus" ,"Mus_spretus","Homo_sapiens","Pongo_pygmaeus","Pongo_abelii","Pan_troglodytes","Pan_paniscus","Gorilla_gorilla") # add in human bonobo for another pair
# restric raxml tree to subset of species:
raxmlTree_renamedTips_subset <- keep.tip(raxmlTree_renamedTips,speciesList)
plot(raxmlTree_renamedTips_subset)
plot(timeTree)
# yay okay this is nice! 

################### get cophenetic distance matrices ###########
timeTree_dist <- cophenetic(timeTree)
raxml_dist <- cophenetic(raxmlTree_renamedTips_subset)

######## read in species codes #########

############ read in spectra distance matrices #########