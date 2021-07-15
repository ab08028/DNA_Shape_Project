######## model 11 additional tuning plots ########3
wd="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210710_model_011a_tune_RF/"
tuningMetrics <- read.table(paste0(wd,"tuningMetrics.txt"),header=T)
head(tuningMetrics)

plot1 <- ggplot(tuningMetrics[tuningMetrics$.metric=="rmse",],aes(x=mtry,y=min_n,color=mean))+
  geom_point()+
  facet_wrap(~.metric)+
  scale_color_viridis_c(direction = -1,option = "D")
plot1
ggsave(paste0(wd,"tuningMetricsplot.png"),plot1,height=5,width=7)

# plot with log scale:
plot2a <- ggplot(tuningMetrics[tuningMetrics$.metric=="rsq",],aes(x=mtry,y=mean))+
  geom_point()+
  ggtitle("rsq mtry")+
  theme(text=element_text(size=14))
plot2a
ggsave(paste0(wd,"tuningMetricsplot.rsq.mtry.png"),plot2a,height=5,width=7)


plot2b <- ggplot(tuningMetrics[tuningMetrics$.metric=="rsq",],aes(x=,min_n,y=mean))+
  geom_point()+
  ggtitle("rsq min_n")+
  theme(text=element_text(size=14))
plot2b
ggsave(paste0(wd,"tuningMetricsplot.rsq.min_n.png"),plot2b,height=5,width=7)
