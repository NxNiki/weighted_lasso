#!/usr/bin/env Rscript
slash = .Platform$file.sep

# alpha: 0 ridge regression, type.measure: class
file_dir1 = "result_d03_glmnet_folds_5_5_alpha_0_none_logtransform_0_2018-09-18"
file_name1 = "glmnetfolds_5_5_alpha_0_none_logtransform_02018_Sep_18_23_04_17.csv" 

# alpha: 1 lasso, type.measure: class
file_dir2 = "result_d03_glmnet_folds_5_5_alpha_1_none_logtransform_0_2018-09-18"
file_name2 = "glmnetfolds_5_5_alpha_1_none_logtransform_02018_Sep_18_23_32_00.csv" 

file1 = read.table(paste(file_dir1, slash, file_name1, sep = ""), sep = ",", header = T)
file2 = read.table(paste(file_dir2, slash, file_name2, sep = ""), sep = ",", header = T)

#file_dir3 = "result_d03_glmnet_folds_5_5_alpha_tuning_none_logtransform_0_2018-05-01"
#file_name1 = "glmnetfolds_5_5_alpha_tuning_none_logtransform_02018_May_01_17_56_28.csv"

figure_name = "class"
label_name = c("ridge", "lasso")

print(file1)
print(file2)

library(ggplot2)
library(reshape2)


library(ggpubr)
# select column of two files:
column.idx = c(2,6)
y.label = c("Gender_Acc", "Gender_reproducibility", "Gender_robustness")

for (i.plot in 1:length(column.idx)){
  
  data.plot = rbind(file1[,c(1, column.idx[i.plot])], file2[,c(1, column.idx[i.plot])]) 
  names(data.plot)[1:2] = c("feature", "value")
  label = c(rep(label_name[1], dim(file1)[1]), rep(label_name[2], dim(file2)[1]))
  data.plot$method = label
  
  print(data.plot)
  print(dim(data.plot))
  
  plot.name = paste(figure_name, y.label[i.plot], label_name[1], "vs", label_name[2],".png", sep = "_")
  
  #png(plot.name)
  ggline(data.plot, "feature", "value",
         linetype = "method", shape = "method",
         color = "method", palette = c("#00AFBB", "#E7B800")) + rotate_x_text(45)
  #dev.off()
}