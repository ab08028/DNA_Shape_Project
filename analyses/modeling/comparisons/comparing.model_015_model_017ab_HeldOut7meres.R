##### comparing fits of models 17a, 17b and 15
plotdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/comparisonPlots/"
model15 <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210713_model_015_modelingUnseen7mers.LeaveOutSameSetBothSpp_ranOnLaptop/rsq.percentraltype.txt",header=T)
model15$model <- "model_015 (seq+shape)"

model17a <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210714_model_017a_modelUnseen7mers.SequenceFeatsOnly_ranOnLaptop/rsq.percentraltype.txt",header=T)
model17a$model <- "model_017a (seq only)"

model17b <- read.table("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210714_model_017b_modelUnseen7mers.ShapeFeatsOnly_ranOnLaptop/rsq.percentraltype.txt",header=T)
model17b$model <- "model_017b (shape only)"

allmodels <- bind_rows(model15,model17a,model17b)
head(allmodels)
rsqplot <- ggplot(allmodels,aes(x=centralMutationType,y=.estimate,fill=model))+
  facet_wrap(~population)+
  geom_col(position="dodge")+
  theme_bw()+
  ggtitle("comparing rsq of held-out 7mer models")+
  theme(text=element_text(size=14))
rsqplot
ggsave(paste(plotdir,"comparingRsq.model_015.model_017a.model_017b.holdingout7mers.png"),rsqplot,height=8,width=14)
