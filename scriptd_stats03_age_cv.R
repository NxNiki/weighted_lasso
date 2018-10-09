#!/usr/bin/env Rscript
# Xin Niu May.4.2017

# notes: adding penalty.weights (mean.diff.boot) causes error on no.brain features. remove the no.brain elements (NULL) in feature.list (in scriptd_stats01_2_...) when using penalty.weights.

rm(list=ls())
graphics.off()


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
report.name="glmnet_" # do not change as it is called in script. change file.name.tail instead.

glmnet.para = list(nfolds = 5, nfolds.inner = 5, family = "gaussian")
glmnet.para$type.measure = "mse"
glmnet.para$predict.type = "response"
glmnet.para$lambda.seq = 10^seq(-6,3, length = 90)

#"mean.diff.boot" #"mean.diff" #"cv.wilcox" #"none"
glmnet.para$penalty.weight = "mean.diff.boot"
#glmnet.para$penalty.weight = "none"
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

file.name.tail = paste("nfold_age", toString(glmnet.para$nfolds), toString(glmnet.para$nfolds.inner), "alpha", alpha.name, glmnet.para$penalty.weight, 'logtransform', toString(glmnet.para$log.penalty.weight), sep = "_")

#library(reshape2)
source("scriptd_stats01_1_read_feature2.R")
#source("scriptd_stats01_1_read_feature_fc.R")
source("scriptd_stats02_cv_functions.R")

# construct multimodal features:
# feature.name and feature.list will be constructed:
#source("scriptd_stats01_2_multimodal_feature_bn246.R")
source("scriptd_stats01_2_multimodal_feature.R")
#source("scriptd_stats01_2_multimodal_feature_fc.R")

## add feautre.list with empty brain.feature to examine result of using only gender and age:
#feature.list = c(list(data.frame()), feature.list)
#feature.name = c("gender", feature.name)


report.rows = length(feature.list)
print(report.rows)
print(length(feature.name))
report  =  data.frame(group = feature.name, 
	accuracy = rep(NA,report.rows), std = rep(NA, report.rows), Reproducibility = rep(NA, report.rows)
	)
report = cbind(report, report)

# select ptsd type separately to predict age: 0: hc, 1: trauma, 2: ptsd
ptsd.type = c(0,2) 
ptsd.name = c("hc", "trauma", "ptsd")

for (i.ptsd in 1:length(ptsd.type)){

	colnames(report)[(i.ptsd-1)*4+1] = paste(ptsd.name[ptsd.type[i.ptsd]+1])
	
	for (i.feature in 1:length(feature.list)) {
		
		print(" running on feature:")
		print(feature.name[i.feature])
		
		brain.feature = scale(feature.list[[i.feature]])
		
		if (dim(brain.feature)[1]!=0){
			df.all = cbind(subject.info[,-1], brain.feature)
		} else {
			# do not add brain feature:
			df.all = subject.info[, -1]
		}
				
		subset.idx = df.all$ptsd==ptsd.type[i.ptsd]
		df.subset = subset(df.all[subset.idx,], select = -ptsd)
		factor = df.subset$age	
		
		print("dimension for subset dataset")
		print(dim(df.subset))
		ptsd.comparision.name = paste("ptsd", toString(i.ptsd), sep = "_")
	
		#print(head(df.subset))
		if (report.name=="glmnet_"){
			#x = as.matrix(subset(df.subset, select = -c(ptsd, Sex)))
			#x1 = as.matrix(subset(df.subset, select = -c(ptsd, age)))
			x = model.matrix(df.subset$age~., df.subset)[,-1] # remove 1st column (intercept)
			
			y = as.numeric(factor)
			#print(head(x))
			#print(head(y))
			cv.result = glmnet.cv.fun(x, y, glmnet.para)
			
			coefs.name = paste("coefs", ptsd.comparision.name, file.name.tail, feature.name[i.feature], '.csv', sep="_")
			write.table(cv.result[[3]], coefs.name, sep=",", row.names=F)
		}
		else if (report.name=="svm_"){
			subject.info.subset = subset(subject.info, select = -c(SUBJID, ptsd))
			# change the column name of the group index to meet convention of svm.cv.fun:
			#colnames(subject.info.subset)[names(subject.info.subset)=="ptsd"]="factor"
			cv.result = svm.cv.fun(brain.feature[subset.idx,], subject.info.subset, cost.seq, svm.para)
			num.feature.tune.name = paste("num.feature.tune_", ptsd.comparision.name, file.name.tail, feature.name[i.feature], '.csv', sep="_")
			write.table(cv.result[[4]], num.feature.tune.name, sep=",", row.names=F)
		}
		# save result:
		test.result = cv.result[[1]]
		train.result = cv.result[[2]]
		coefs.result = cv.result[[3]]
		print(test.result)
		print(train.result)
		#print(coefs.result)
	
		report[i.feature, (i.ptsd-1)*4 + 2] = mean(test.result[,1], na.rm=T)
		report[i.feature, (i.ptsd-1)*4 + 3] = sd(test.result[,2], na.rm=T)
		report[i.feature, (i.ptsd-1)*4 + 4] = mean(as.numeric(coefs.result[, ncol(coefs.result)]), na.rm=T)
	}
}
# -------------------- save results : -------------------------------------
print(report)
#print(report.sd)

filename = paste("result_m_",report.name, file.name.tail, Sys.time(), ".csv", sep = "")
write.table(report, filename, sep = ",", row.names = F)
#filename = paste("result_sd_", file.name.tail, Sys.time(), ".csv", sep = "")
#write.table(report.sd, filename, sep = ",", row.names = F)

