########## sum up hamming distances across intervals/chromosomes ###########
# want this to work on sage eventually
# 
require(reshape2)
require(ggplot2)
require(dplyr)
#labels=c('humans','bears','mice','fin_whale')
labels=c("merged_apes")
for(label in labels){
wd=paste0("/net/harris/vol1/home/beichman/DNAShape/analyses/hamming_distance/",label,"/")
indir=paste0(wd,"/perInterval/")
list.files(indir)
######## count up dist files ########
distFiles=list.files(indir,pattern=".dist$") # list dist files
intervalCount=length(distFiles)
intervals=unlist(lapply(strsplit(distFiles,"\\."),"[",4))
#intervals # want to keep as strings because some will be '1' and others will be 'chr1'
# note this misses if some chromosomes have been left off
# maybe code up better in future (missed that humans 20-22 were missing but i've fixed that now -20211122)
allDistances=data.frame()
for(interval in intervals){
  distdf=read.table(paste0(indir,"plink.",label,".interval.",interval,".FromMutyperVariantsVCF.dist"),header=F)
  ids=read.table(paste0(indir,"plink.",label,".interval.",interval,".FromMutyperVariantsVCF.dist.id"),header=F)
  colnames(distdf) <- ids$V2 # colnames 
  head(distdf)
  distdf$ind2 <- as.character(ids$V2) # this works because its a sqaure with 0s on diagonal 
  distdf_melt <- melt(distdf,variable.name = "ind1") # note ind2 is a character, ind1 is factor - leads to some annoyances
  # get rid of the empty part of the triangle (0 entries):
  distdf_melt <- distdf_melt[distdf_melt$value!=0,]
  distdf_melt$interval <- as.character(interval)
  distdf_melt$label <- label
  # make an alphabetical comparison label so that even if are in different orders this label will always label the same comparisons the same way: 
  distdf_melt$ind1_alphabetical <- pmin(as.character(distdf_melt$ind1),distdf_melt$ind2) # note his is ind1 and ind2 not just ind1 -- is selecting the first alphabetically
  distdf_melt$ind2_alphabetical <- pmax(as.character(distdf_melt$ind1),distdf_melt$ind2) # note this should be ind1 and ind2
  distdf_melt$comparisonLabel <-  paste0(distdf_melt$ind1_alphabetical,".",distdf_melt$ind2_alphabetical)
  allDistances <- bind_rows(allDistances,distdf_melt)
  
}
head(allDistances)
if(length(unique(allDistances$interval))!=length(intervals)){
  print('something is wrong!')
  break
}
allDistances_summedUp <- allDistances %>%
  group_by(label,comparisonLabel,ind1_alphabetical,ind2_alphabetical) %>% # don't gorup by ind1 and ind2 any more
  summarise(totalHammingDistance=sum(value))
# dim of this should be sample size choose 2:
dim(allDistances_summedUp)[1]==choose(dim(ids)[1],2) # good
#allDistances_summedUp$ind2 <- factor(allDistances_summedUp$ind2,levels=sort(unique(as.character(allDistances_summedUp$ind2))))
#allDistances_summedUp$ind1 <- factor(allDistances_summedUp$ind1,levels=sort(unique(as.character(allDistances_summedUp$ind1),decreasing = T)))
distPlot <- ggplot(allDistances_summedUp,aes(x=ind1_alphabetical,y=ind2_alphabetical,fill=totalHammingDistance))+
  geom_tile()+
  theme(axis.text.x = element_text(angle=45))+
  scale_fill_gradient(low="blue",high="red")+
  ylab("")+
  xlab("")
distPlot
ggsave(paste0(wd,label,"totalAlleleCounts.HammingDistance.allIntervals.png"),distPlot,height=8,width=10)
write.table(allDistances_summedUp,paste0(wd,label,".totalAlleleCounts.HammingDistance.allIntervals.txt"),row.names=F,quote=F,sep="\t")


}
