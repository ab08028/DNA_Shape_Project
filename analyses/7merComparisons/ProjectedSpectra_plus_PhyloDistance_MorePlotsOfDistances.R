############ making more plots ###########
wd="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/comparing_spectrum_and_phylo_divergence_WITHPROJECTIONS_CpGSitesSeparatedOut/"
# these distances all are based on sample-size corrected (projected) counts with different epsilons added


fracSegSitesDistances_epsilon1=read.table(paste0(wd,"ALLDISTANCETYPES.includesNonTransformed.distances.FracSegSites.epsilon.1.PROJECTEDCOUNTS.txt"),header=T)

fracSegSitesDistances_epsilon0.1=read.table(paste0(wd,"ALLDISTANCETYPES.includesNonTransformed.distances.FracSegSites.epsilon.0.1.PROJECTEDCOUNTS.txt"),header=T)
                                       
fracSegSitesDistances_epsilon10=read.table(paste0(wd,"ALLDISTANCETYPES.includesNonTransformed.distances.FracSegSites.epsilon.10.PROJECTEDCOUNTS.txt"),header=T)

mutRateDistances_epsilon1=read.table(paste0(wd,"ALLDISTANCETYPES.includesNonTransformed.distances.NormMutRate.epsilon.1.PROJECTEDCOUNTS.txt"),header=T)

mutRateDistances_epsilon0.1=read.table(paste0(wd,"ALLDISTANCETYPES.includesNonTransformed.distances.NormMutRate.epsilon.0.1.PROJECTEDCOUNTS.txt"),header=T)

mutRateDistances_epsilon10=read.table(paste0(wd,"ALLDISTANCETYPES.includesNonTransformed.distances.NormMutRate.epsilon.10.PROJECTEDCOUNTS.txt"),header=T)

# want to plot spectrum distance between different 7mer types

# so want to pivot (?) distances (?)

# convert distances to ranks 

rank_and_pivot_function <- function(df) {
  # add rank:
  df <- df %>%
    group_by(label,variable) %>%
    mutate(rankLowestToHighestDistance = rank(spectrum_distance, ties.method = "first"))
  df_wide <-  pivot_wider(df,id_cols=c("comparisonLabel","comparisonLabel_broad_alphabetical","comparisonLabel_common_alphabetical"),names_from=c("label"),values_from =rankLowestToHighestDistance)
  
  return(df_wide)
}

fracSegSitesDistances_epsilon1_rank <- rank_and_pivot_function(fracSegSitesDistances_epsilon1)
fracSegSitesDistances_epsilon0.1_rank <- rank_and_pivot_function(fracSegSitesDistances_epsilon0.1)
fracSegSitesDistances_epsilon10_rank <- rank_and_pivot_function(fracSegSitesDistances_epsilon10)

mutRateDistances_epsilon1_rank <- rank_and_pivot_function(mutRateDistances_epsilon1)
mutRateDistances_epsilon0.1_rank <- rank_and_pivot_function(mutRateDistances_epsilon0.1)
mutRateDistances_epsilon10_rank <- rank_and_pivot_function(mutRateDistances_epsilon10)


############ plot mut rate #########
png(paste0(wd,"pairsPlot.clr.mutationRate.epsilon1.png"),height=5,width=9,units="in",res=300)
pairs(select(mutRateDistances_epsilon1_rank,c("1mer_clr","3mer_clr","5mer_clr","7mer_clr")),upper.panel=NULL,cex.labels = 1.4)
dev.off()

png(paste0(wd,"pairsPlot.clr.mutationRate.epsilon0.1.png"),height=5,width=9,units="in",res=300)
pairs(select(mutRateDistances_epsilon0.1_rank,c("1mer_clr","3mer_clr","5mer_clr","7mer_clr")),upper.panel=NULL,cex.labels = 1.4)
dev.off()

########## plot frac seg sites ########
png(paste0(wd,"pairsPlot.clr.fracSegSites.epsilon1.png"),height=5,width=9,units="in",res=300)
pairs(select(fracSegSitesDistances_epsilon1_rank,c("1mer_clr","3mer_clr","5mer_clr","7mer_clr")),upper.panel=NULL,cex.labels = 1.4)
dev.off()

png(paste0(wd,"pairsPlot.clr.fracSegSites.epsilon0.1.png"),height=5,width=9,units="in",res=300)
pairs(select(fracSegSitesDistances_epsilon0.1_rank,c("1mer_clr","3mer_clr","5mer_clr","7mer_clr")),upper.panel=NULL,cex.labels = 1.4)
dev.off()

########### plot epsilon 1 vs 0.1 ###########
mutRateDistances_epsilon0.1_rank$epsilon <- "epsilon_0.1"
mutRateDistances_epsilon1_rank$epsilon <- "epsilon_1"
mutRateDistances_epsilon10_rank$epsilon <- "epsilon_10"

mutRate_bothEpsilons <- bind_rows(mutRateDistances_epsilon0.1_rank,mutRateDistances_epsilon1_rank,mutRateDistances_epsilon10_rank)
mutRate_bothEpsilons_wide <- pivot_wider(mutRate_bothEpsilons,id_cols = c("comparisonLabel"),values_from=contains("clr"),names_from="epsilon")

png(paste0(wd,"pairsPlot.clr.MutationRate.ComparingEpsilonValues.png"),height=10,width=14,units="in",res=300)
pairs(select(mutRate_bothEpsilons_wide,contains("epsilon")),upper.panel=NULL)
dev.off()
######### impact of transform: 
png(paste0(wd,"pairsPlot.clr.plusNoTransform.mutationRate.epsilon1.png"),height=12,width=16,units="in",res=300)
pairs(select(mutRateDistances_epsilon1_rank,c(contains("ilr"),contains("clr"),contains("Transform"))),upper.panel=NULL,cex.labels = .8)
dev.off()



