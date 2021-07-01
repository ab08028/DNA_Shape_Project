######### sandbox: experimenting with vint output

####3 set up a tuned RF model
####### try a tuned xgboost model 
vint1 <-readRDS("/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/modeling/experiments/20210618_model_002_RF.MutationRate.Log10.Rescaled.Multispecies.DummyPopVar/vint/vint.intermediate.1.rds")
vint1
summary(vint1)
vint1
# okay this yields interaction strength 
# A tibble: 1 x 2
#Variables                                             Interaction

# feature_1_Shear_2.derived*feature_1_Shear_2.ancestral      0.0143
# could maybe calc interaction between spp and all others? or get shap interactions?