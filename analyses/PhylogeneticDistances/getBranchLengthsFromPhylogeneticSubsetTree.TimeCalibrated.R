############### getting branch lengths of phylogenetic tree #################
require(ape)
require(plyr)
require(dplyr)
require(reshape2)
# this is

# these are from vertlife # go to vert life: http://vertlife.org/phylosubsets/
# and get pruned trees (get mlutiple so you can average distances)

# I think branch lengths are time-corrected to be MYA (?) -- so not ideal for mice and things with different generation times. But at least is based on phylogenetic distances -- branch lengths.... broad strokes. Talk to KH bout this. I would prefer just substitutions rather than time... can I get that out?
# I need to compare node and tip scaled trees and I need to average over 100+ trees (1000s probably better)
# for now I'm using node scaled and I'm using 100 trees
date="20211011"
datatypes=c("All_Species","DNA_only")
datingtypes=c("node_dated","tip_dated")
############## YOU ARE HERE --- only did it for sandbox set of trees, need to get for all.
# finish up
wd=paste0("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/vertlife_phylogeneticTrees/time_calibrated_trees/",date,"/")
#taskid="tree-pruner-02091bed-8b42-4a93-bf74-0e7fe22fe48c"
for(datatype in datatypes){
  for(datingtype in datingtypes){
    indir=paste0(wd,datatype,"/",datingtype,"/")
    jobid=list.files(indir)
    nexusFile = read.nexus(paste0(indir,jobid,"/output.nex"))
    
    allTreesPairwiseDistances <- lapply(nexusFile,cophenetic.phylo) # okay this gets pairwise distances between all pairs for all trees
    totalTrees=length(allTreesPairwiseDistances)
    # https://www.rdocumentation.org/packages/ape/versions/5.5/topics/cophenetic.phylo
    # "Values in the matrix are the sum of the branch lengths separating each pair of species."  sum of branch lengths . so not just the divergence time. it's the sum of branch lengths between tips (so includes amnt of evoution on both branches.) so sholud it be ~2x divergence time?
    allTreesPairwiseDistances_melt <- melt(allTreesPairwiseDistances) # oh wow melt does so well! puts them all together and labels them
    # so easy!
    
    colnames(allTreesPairwiseDistances_melt) <- c("Sp1","Sp2","branchLength","treeLabel")
    
    allTreesPairwiseDistances_melt$comparisonLabel <- paste0(allTreesPairwiseDistances_melt$Sp1,".",allTreesPairwiseDistances_melt$Sp2)
    
    head(allTreesPairwiseDistances_melt)
    ## get avg and std 
    
    allTreesPairwiseDistances_melt$datingType <- datingtype
    allTreesPairwiseDistances_melt$dataType <- datatype
    
    write.table(allTreesPairwiseDistances_melt,paste0(indir,"cophenetic.Pairwise.BranchLengthSums.PerTree.txt"),row.names = F,quote=F,sep="\t")
    
    averages <- allTreesPairwiseDistances_melt %>%
      group_by(Sp1,Sp2,comparisonLabel,datingType,dataType) %>%
      summarize(avgBranchLen=mean(branchLength),sdBranchLen=sd(branchLength)) 
    
    write.table(averages,paste0(indir,"cophenetic.Pairwise.BranchLengthSums.AveragedOverTrees.avg.sd.txt"),row.names = F,quote=F,sep="\t")
  }
}