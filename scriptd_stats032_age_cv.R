#!/usr/bin/env Rscript
# Xin Niu May.4.2017

rm(list=ls())
#graphics.off()

data_dir = "data_d01_features"
slash = .Platform$file.sep
date = Sys.Date()

# svm parameters:
##feature.selection: "coef.variation" "cv.wilcox" "kendall"

#report.name = "svm_"
#cost.seq = 10^seq(-4, 3, length = 30)
#
#svm.para = data.frame(nfolds=5)
#svm.para$feature.selection1 = "wilcox"
#svm.para$num.feature1 = 40
#
#svm.para$feature.selection2 = "cv.wilcox"
#svm.para$num.feature2.start = 5
#svm.para$num.feature2.end = 30
#svm.para$num.feature2.step = 3
#
#file.name.tail = paste("nfolds", toString(svm.para$nfolds), "num.feature1", toString(svm.para$num.feature1), svm.para$feature.selection1, "num.feature2", toString(svm.para$num.feature2.start), toString(svm.para$num.feature2.end), toString(svm.para$num.feature2.step), svm.para$feature.selection2, sep = "_")

# glmnet parameters:
report.name="glmnet" # do not change as it is called in script. change file.name.tail instead.

glmnet.para = list(nfolds = 5, nfolds.inner = 5, family = "gaussian")
glmnet.para$type.measure = "mse"
glmnet.para$predict.type = "response"
glmnet.para$lambda.seq = 10^seq(-6,3, length = 90)

#"mean.diff.boot" #"mean.diff" #"cv.wilcox" #"none"
#glmnet.para$penalty.weight = "mean.diff.boot"
#glmnet.para$penalty.weight = "pearson.boot"
#glmnet.para$penalty.weight = "glmnet.coef.boot"
glmnet.para$penalty.weight = "none"

glmnet.para$log.penalty.weight = 0 # log transform the penalty weight
glmnet.para$alpha = 0 # 1: lasso, 0: ridge regression, 0.5 elastic net
#glmnet.para$alpha = seq(0,1,0.1) # input a sequence of alpha to run cross validation. 
glmnet.para$quantile.thresh = 1
glmnet.para$return.mod = T # return regression model trained on all the data

if (length(glmnet.para$alpha)==1){
	alpha.name = toString(glmnet.para$alpha)
}else{
	alpha.name = "tuning"
}

train.group = 0
test.group = c(1,2) 

file.name.tail = paste("age_folds", toString(glmnet.para$nfolds), toString(glmnet.para$nfolds.inner), "alpha", alpha.name, glmnet.para$penalty.weight, 'logtransform', toString(glmnet.para$log.penalty.weight), sep = "_")

result_dir = paste("result_d03",report.name, file.name.tail, Sys.Date(), sep = "_")
dir.create(result_dir)

library(reshape2)
library(ggplot2)
#library(easyGgplot2)
library(ggpubr)

source("scriptd_stats01_1_read_all_feature2.R")
#source("scriptd_stats01_1_read_feature_fc.R")
source("scriptd_stats02_cv_functions.R")

# construct multimodal features:
# feature.name and feature.list will be constructed:
#source("scriptd_stats01_2_multimodal_feature_bn246.R")
source("scriptd_stats01_2_multimodal_feature.R")
#source("scriptd_stats01_2_multimodal_feature_fc.R")

## add feautre.list with empty brain.feature to examine result of using only gender and age:
#feature.list = c(list(data.frame()), feature.list)
#feature.name = c("demographic", feature.name)

# delete subjects with severe head motion or other problems:
subject.delete = read.table(paste(data_dir, slash, "subject_delete.txt", sep = ""), header=F)
print("deleted subjects:")
print(subject.delete)
#delete.idx = subject.info$SUBJID%in%subject.delete$V1
delete.idx = rep(FALSE, dim(subject.info)[1])

# save deleted subjects id as log:
glmnet.para$delete.subjects = delete.idx

print("ptsd:")
print(table(subject.info[delete.idx,]$ptsd))
print("Sex:")
print(table(subject.info[delete.idx,]$Sex))
report.rows = length(feature.list)
print(report.rows)
print(length(feature.name))

subject.info = subject.info[!delete.idx, ] # do not use -delete.idx here, it will just remove the 1st row!
#subject.info$age = scale(subject.info$age)
print(dim(subject.info))
print(subject.info)

report  =  data.frame(group = feature.name, 
	accuracy.hc = rep(NA,report.rows), std = rep(NA, report.rows), rmse.hc = rep(NA, report.rows), Reproducibility = rep(NA, report.rows), accuracy.trauma = rep(NA, report.rows), rmse.trauma = rep(NA, report.rows),  accuracy.ptsd = rep(NA, report.rows), rmse.ptsd = rep(NA, report.rows) 
	)

list.prediction = list()

# use healthy control to train the model (also run nested CV to evalutate performance) and predict age on ptsd group
ptsd.name = c("hc", "trauma", "ptsd")

	plot.name = paste(result_dir, slash, "hist_age.png", sep = "")
	png(plot.name)
	
	plot.age = subset(subject.info, select = c(ptsd, age))
	plot.age$ptsd[plot.age$ptsd==0] = ptsd.name[1] 
	plot.age$ptsd[plot.age$ptsd==1] = ptsd.name[2] 
	plot.age$ptsd[plot.age$ptsd==2] = ptsd.name[3] 

	p = gghistogram(plot.age, x = "age",
		add = "mean", rug = TRUE,
		color = "ptsd", fill = "ptsd",
		palette = "jco")
	print(p)
	dev.off()


for (i.feature in 1:length(feature.list)) {
#for (i.feature in 1:1) {
	
	print(" running on feature:")
	print(feature.name[i.feature])
	
	brain.feature = scale(feature.list[[i.feature]][!delete.idx,])
	
	if (dim(brain.feature)[1]!=0){
		print("adding brain feature to subject.info")
		df.all = cbind(subset(subject.info, select = -c(SUBJID)), brain.feature)
	} else {
		# do not add brain feature:
		df.all = subset(subject.info, select = -c(SUBJID))
	}
			
	# remove NAN columns:
	df.all = df.all[ , apply(df.all, 2, function(x) !any(is.na(x)))]
	
	subset.idx = df.all$ptsd==train.group
	df.train = subset(df.all[subset.idx,], select = -ptsd)
	factor.train = df.train$age	

	print("dimension for subset dataset")
	print(dim(df.train))
	
	if (report.name=="glmnet"){
		#x = as.matrix(subset(df.subset, select = -c(ptsd, Sex)))
		#x1 = as.matrix(subset(df.subset, select = -c(ptsd, age)))
		x = model.matrix(df.train$age~., df.train)[,-1] # remove 1st column (intercept)
		
		y = as.numeric(factor.train)
		#print(head(x))
		#print(head(y))
		cv.result = glmnet.cv.fun(x, y, glmnet.para)
		
		#coefs.name = paste("coefs", file.name.tail, feature.name[i.feature], '.csv', sep="_")
		#write.table(cv.result$coefs, coefs.name, sep=",", row.names=F)
		
		cv.mod = cv.result$mod
	
		print("test on trauma:")
		idx.trauma = df.all$ptsd==test.group[1]
		df.test.trauma = subset(df.all[idx.trauma,], select = -ptsd)
		# apply model on ptsd samples:	
		x.trauma = model.matrix(df.test.trauma$age~., df.test.trauma)[,-1] # remove 1st column (intercept)
		
		y.trauma = as.numeric(df.test.trauma$age)
		#print(head(x))
		#print(head(y))
		y.pred.trauma = predict(cv.mod, x.trauma, s = cv.mod$lambda.min, type=toString(glmnet.para$predict.type))
		result.cor = cor.test(y.trauma, y.pred.trauma, method = "pearson")
	
		print("test on ptsd:")
		## test on ptsd:
		idx.ptsd = df.all$ptsd==test.group[2]
		df.test.ptsd = subset(df.all[idx.ptsd,], select = -ptsd)
		# apply model on ptsd samples:	
		x.ptsd = model.matrix(df.test.ptsd$age~., df.test.ptsd)[,-1] # remove 1st column (intercept)
		
		y.ptsd = as.numeric(df.test.ptsd$age)
		#print(head(x))
		#print(head(y))
		y.pred.ptsd = predict(cv.mod, x.ptsd, s = cv.mod$lambda.min, type=toString(glmnet.para$predict.type))
		result.cor2 = cor.test(y.ptsd, y.pred.ptsd, method = "pearson")
		

	}
	else if (report.name=="svm_"){
		subject.info.subset = subset(subject.info, select = -c(SUBJID, ptsd))
		# change the column name of the group index to meet convention of svm.cv.fun:
		#colnames(subject.info.subset)[names(subject.info.subset)=="ptsd"]="factor"
		cv.result = svm.cv.fun(brain.feature[subset.idx,], subject.info.subset, cost.seq, svm.para)
		num.feature.tune.name = paste("num.feature.tune", file.name.tail, feature.name[i.feature], '.csv', sep="_")
		write.table(cv.result[[4]], num.feature.tune.name, sep=",", row.names=F)
	}
	# save result:
	test.result = cv.result$test.result
	train.result = cv.result$train.result
	coefs.result = cv.result$coefs
	print(test.result)
	print(train.result)
	#print(coefs.result)

	report[i.feature, 2] = mean(test.result[,1], na.rm=T)
	report[i.feature, 3] = sd(test.result[,1], na.rm=T)
	# rmse hc:
	report[i.feature, 4] = mean(test.result[,2], na.rm=T)
	# reproducibility index:
	report[i.feature, 5] = mean(as.numeric(coefs.result[, ncol(coefs.result)]), na.rm=T)
	# prediction correlation on test.group:
	report[i.feature, 6] = result.cor$estimate
	report[i.feature, 7] = sqrt(mean((y.trauma - y.pred.trauma)^2))
	report[i.feature, 8] = result.cor2$estimate
	report[i.feature, 9] = sqrt(mean((y.ptsd - y.pred.ptsd)^2))


	# scatter plot:
	plot.name = paste(result_dir, slash, "scatter_", file.name.tail, feature.name[i.feature], '.pdf', sep="")

	df.trauma = data.frame(Imaging = rep(feature.name[i.feature], length(y.trauma)), SUBJID = subject.info[idx.trauma,1], age = y.trauma, predicted_age = y.pred.trauma, group = rep("trauma", length(y.trauma)))
	df.ptsd = data.frame(Imaging = rep(feature.name[i.feature], length(y.ptsd)), SUBJID = subject.info[idx.ptsd,1], age = y.ptsd, predicted_age = y.pred.ptsd, group = rep("ptsd", length(y.ptsd)))

	plot.data = rbind(df.trauma, df.ptsd)
	list.prediction[[i.feature]] = plot.data

	# it is werid the name of the 2nd column is X1 rather than what we defined, so we change it:
#	colnames(plot.data)[2] = "age"
	colnames(plot.data)[4] = "predicted_age"
#	colnames(plot.data)[4] = "group"
	print(head(plot.data))
	#png(plot.name)
	
	legend1 = sprintf(result.cor$estimate, fmt = 'r = %#.2f (trauma)') 
	legend2 = sprintf(result.cor2$estimate, fmt = 'r = %#.2f (ptsd)')
	#p = ggplot(data = plot.data, aes(x=age, y=predicted_age, color = group, shape=group, palette =c("#00AFBB", "#E7B800"))) +
	p = ggplot(data = plot.data, aes(x=age, y=predicted_age, color = group, shape=group)) +
			geom_point(size = 3) +	
			geom_smooth(method = lm, fullrange = T) +
			theme(text = element_text(size=25), legend.position = c(0.2, .9), legend.title = element_blank()) + 
			scale_color_manual(values = c("#00AFBB", "#E7B800"), 
				breaks=c("trauma", "ptsd"),
			    labels=c(legend1, legend2)) +
			scale_shape_discrete(
				breaks=c("trauma", "ptsd"),
			    labels=c(legend1, legend2)) 
		#	annotate("text", x = min(plot.data$age) + (max(plot.data$age)-min(plot.data$age))/10, y=max(plot.data$predicted_age), 
		#		label = paste(sprintf(result.cor$estimate, fmt = 'r = %#.2f (trauma)\n'), 
		#		sprintf(result.cor2$estimate, fmt = 'r = %#.2f (ptsd)')))
	print(p)
	ggsave(plot.name, device = "pdf")
	dev.off()

	# permutation test:
	perm.t = perm.t.test(x1 = y.pred.trauma - y.trauma, x2 = y.pred.ptsd - y.ptsd, n = 5000, paired = F)

	# plot age difference for trauma and ptsd:
	plot.data$age.diff = plot.data$predicted_age - plot.data$age
	plot.name = paste(result_dir, slash, "boxplot_age_difference_", file.name.tail, feature.name[i.feature], '.eps', sep="")
	#png(plot.name)

	p <- ggboxplot(plot.data, y = "age.diff", x = "group",
	        color = "group", palette =c("#00AFBB", "#E7B800"),
			add = "jitter", shape = "group") +
			theme(text = element_text(size=25), legend.position = "none", axis.title.x = element_blank()) + 
			annotate(geom = "text", x = 1.5, y = max(plot.data$age.diff), size = 7,
				label = (sprintf(perm.t$p, fmt = 'permutation t-test: %#.2f'))) +
			labs(y = "The age difference")

	print(p)
	ggsave(plot.name, device = "eps")
	dev.off()

	# paired bar plot of age and predict age for trauma and ptsd separately:
	plot.data.trauma = data.frame(age = y.trauma, predicted_age = y.pred.trauma)
	colnames(plot.data.trauma)[2] = "predicted_age"
	plot.data.ptsd = data.frame(age = y.ptsd, predicted_age = y.pred.ptsd)
	colnames(plot.data.ptsd)[2] = "predicted_age"
	
	perm.t = perm.t.test(x1 = y.trauma, x2 = y.pred.trauma, n = 5000, paired = T)
	
	plot.name = paste(result_dir, slash, "boxplot_trauma_", file.name.tail, feature.name[i.feature], '.png', sep="")
	png(plot.name)
	p = ggpaired(plot.data.trauma, cond1 = "age", cond2 = "predicted_age",
			line.color = "gray", line.size = 0.4,
			color = "condition", palette = "jco") + 
			#stat_compare_means(method = "t.test", paired = TRUE, label.x = 1.5) 
			annotate("text", x = 1.5, y = max(plot.data.trauma),
				label = (sprintf(perm.t$p, fmt = 'permutation t-test: %#.2f')))

	
	print(p)
	dev.off()
	
	perm.t = perm.t.test(x1 = y.ptsd, x2 = y.pred.ptsd, n = 5000, paired = T)
	
	plot.name = paste(result_dir, slash, "boxplot_ptsd_", file.name.tail, feature.name[i.feature], '.eps', sep="")
	#png(plot.name)
	
	p = ggpaired(plot.data.ptsd, cond1 = "age", cond2 = "predicted_age",
			line.color = "gray", line.size = 0.4,
			color = "condition", palette = "jco") + 
			theme(text = element_text(size=20))+ 
			#stat_compare_means(method = "t.test", paired = TRUE, label.x = 1.5) 
			annotate("text", x = 1.5, y = max(plot.data.ptsd),
				label = (sprintf(perm.t$p, fmt = 'permutation t-test: %#.2f')))
	
	print(p)
	ggsave(plot.name, device = "eps")
	dev.off()

#	plot.name = paste(result_dir, slash, "hist_residual_", file.name.tail, feature.name[i.feature], '.png', sep="")
#	png(plot.name)
#
#	plot.residual = data.frame(age = y.trauma, predicted_age = y.pred.trauma, residual = y.pred.trauma - y.trauma)
#	colnames(plot.residual) = c("age", "predicted_age", "residual")
#	print(head(plot.residual))
#	histogram.df = melt(plot.residual)
#	p = gghistogram(histogram.df, x = "value",
#		add = "mean", rug = TRUE,
#		color = "variable", fill = "variable",
#		palette = "jco")
#	print(p)
#	dev.off()

}
# -------------------- save results : -------------------------------------
print(report)
#print(report.sd)


time = format(Sys.time(), '%Y_%b_%d_%H_%M_%S')
filename = paste(result_dir, slash, report.name, file.name.tail, time, ".csv", sep = "")
write.table(report, filename, sep = ",", row.names = F)

df.prediction = do.call("rbind", list.prediction)
filename = paste(result_dir, slash, "raw_prediction", report.name, file.name.tail, time, ".csv", sep = "")
write.table(df.prediction, filename, sep = ",", row.names = F)

# save parameter to text file:
sink(paste(result_dir, slash, slash, "para", time, ".txt", sep = ""))
print(glmnet.para)
sink()



