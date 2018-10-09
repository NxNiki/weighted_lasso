#!/usr/bin/env Rscript

behav.data = read.table("data_a01_behav/Output01_Pnc.csv", sep = ",", header = T)

match.id = read.table("data_a01_behav/match_id.txt", sep = " ", header = T)

behav.data = merge(match.id, behav.data, by = "bblid")
print(behav.data)

#write.table(behav.data, file = "data_a01_behav/resulta01_behav_data_matchid.csv", sep = ",", row.names = F)

age.pred = read.table("result_d03_glmnet_age_folds_5_5_alpha_0_none_logtransform_0_2018-05-01/age_predicted_age.csv", sep = ",", header = T)

behav.data = merge(age.pred, behav.data, by  = "SUBJID")

write.table(behav.data, file = "data_a01_behav/resulta01_behav_data_matchid.csv", sep = ",", row.names = F)

