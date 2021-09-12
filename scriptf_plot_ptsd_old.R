#!/usr/bin/env Rscript

######################
## this is OLD results use scriptf_plot_result.Rmd in the same folder.
######################

slash = .Platform$file.sep

#file_dir1 = "result_d03_glmnet_folds_5_5_alpha_1_none_logtransform_0_2018-04-27"
#file_name1 = "glmnet_folds_5_5_alpha_1_none_logtransform_02018_Apr_27_15_52_36.csv" #class
#file_name2 = "glmnet_folds_5_5_alpha_1_none_logtransform_02018_Apr_27_17_27_40.csv" #mse

file_dir1 = "result_d03_glmnet_folds_5_5_alpha_1_mean.diff.boot_logtransform_0_2018-05-11"
file_name1 = "glmnetfolds_5_5_alpha_1_mean.diff.boot_logtransform_02018_May_11_18_14_55.csv"
file_name2 = "glmnetfolds_5_5_alpha_1_mean.diff.boot_logtransform_02018_May_11_18_15_48.csv"
file_name3 = "glmnetfolds_5_5_alpha_1_mean.diff.boot_logtransform_02018_May_11_18_16_08.csv"

file1 = read.table(paste(file_dir1, slash, file_name2, sep = ""), sep = ",", header = T)

#file_dir2 = "result_d03_glmnet_folds_5_5_alpha_1_mean.diff.boot_logtransform_0_2018-05-01"
#file_name1 = "glmnetfolds_5_5_alpha_1_mean.diff.boot_logtransform_02018_May_01_16_29_14.csv" #mse
#file_name2 = "glmnetfolds_5_5_alpha_1_mean.diff.boot_logtransform_02018_May_01_17_10_40.csv" #class

file_dir2 = "result_d03_glmnet_folds_5_5_alpha_1_none_logtransform_0_2018-05-11"
file_name1 = "glmnetfolds_5_5_alpha_1_none_logtransform_02018_May_11_17_52_13.csv"
file_name2 = "glmnetfolds_5_5_alpha_1_none_logtransform_02018_May_11_17_53_15.csv"
file_name3 = "glmnetfolds_5_5_alpha_1_none_logtransform_02018_May_11_17_53_57.csv"

#file_dir3 = "result_d03_glmnet_folds_5_5_alpha_tuning_none_logtransform_0_2018-05-01"
#file_name1 = "glmnetfolds_5_5_alpha_tuning_none_logtransform_02018_May_01_17_56_28.csv"

figure_name = "auc"
label_name = c("weighted_lasso", "lasso")

file2 = read.table(paste(file_dir2, slash, file_name2, sep = ""), sep = ",", header = T)

library(ggplot2)
library(reshape2)

column.idx = c(2,6,10,5,9,13)
y.label = c("HC_Trauma_Acc", "HC_PTSD_Acc", "PTSD_Trauma_Acc",
			"HC_Trauma_Reproducibility", "HC_PTSD_Reproducibility", "PTSD_Trauma_Reproducibility")

for (i.plot in 1:length(column.idx)){

data.plot = rbind(file1[,c(1, column.idx[i.plot])], file2[,c(1, column.idx[i.plot])]) 
names(data.plot)[1:2] = c("feature", "value")
label = c(rep(label_name[1], dim(file1)[1]), rep(label_name[2], dim(file2)[1]))
data.plot$method = label 

print(data.plot)
print(dim(data.plot))

#png(paste(figure_name, y.label[i.plot], label_name[1], "vs_", label_name[2],".png", sep = "_"), width = 1100, height = 400)
plot.name = paste(figure_name, y.label[i.plot], label_name[1], "vs", label_name[2],".pdf", sep = "_")

p = ggplot(data = data.plot, aes(y = value, x = feature, group = method, color = method)) +
		geom_point(size = 1) + geom_line(size = 1) +
		geom_text(data = data.plot[data.plot$method==label_name[1],], aes(label=sprintf("%0.2f", value), color=method), hjust = 1, position = position_dodge(width = .9), size=5) +
		geom_text(data = data.plot[data.plot$method==label_name[2],], aes(label=sprintf("%0.2f", value), color=method), hjust = 0, position = position_dodge(width = .9), size=5) +
    	#scale_fill_brewer(palette="Paired") +
		theme_minimal() +
		labs(y = y.label[i.plot])+
		coord_cartesian(ylim = c(min(data.plot$value)-.03, max(data.plot$value)+.03)) +
		theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 15, angle = 90, hjust = 0), legend.title = element_blank(), legend.position = c(.1, .9), legend.key.size = unit(3, "cm"))	
		#coord_flip()

print(p)
ggsave(plot.name, device = "pdf", width = 11, height = 4, dpi = 300)
dev.off()

}

