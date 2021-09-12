#!/usr/bin/env Rscript
# Xin Niu May.4.2017


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
report.name="glmnet_"
glmnet.para = data.frame(nfolds = 5)

#"mean.diff" #"wilcox" #"none"
glmnet.para$penalty.weight = "pearson.boot"

glmnet.para$alpha = 1 # 1: lasso, 0: ridge regression, 0.5 elastic net
glmnet.para$quantile.thresh = 1

lambda.seq = 10^seq(-6,3, length = 90)

file.name.tail = paste("Sex_nfold", toString(glmnet.para$nfolds), "alpha", toString(glmnet.para$alpha), glmnet.para$penalty.weight, sep = "_")
nfolds = 5


#library(reshape2)
#source("scriptd_stats01_read_all_feature2.R")
#source("scriptd_stats01_read_feature2.R")
source("scriptd_stats01_read_feature_fc.R")
#source("scriptd_stats01_read_feature_fc.R")
source("scriptd_stats02_cv_functions.R")

subject.info$age = scale(subject.info$age)
subject.info$Sex = as.factor(subject.info$Sex)

#multimodal.feature = scale(cbind(spm.vbm, alff, falff, reho, label.fa, tract.fa, label.md, tract.md))
#multimodal.feature.aal = scale(cbind(spm.vbm.aal, alff.aal, falff.aal, reho.aal, label.fa, tract.fa, label.md, tract.md))
#multimodal.feature.bn246 = scale(cbind(spm.vbm.bn246, alff.bn246, falff.bn246, reho.bn246, label.fa, tract.fa, label.md, tract.md))
#
#mf.dti = scale(cbind(label.fa, label.md, tract.fa, tract.md))
#dti.fa = scale(cbind(label.fa, tract.fa))
#
#mf.vbm.dti = scale(cbind(spm.vbm, label.fa, label.md, tract.fa, tract.md))
#mf.vbm.aal.dti = scale(cbind(spm.vbm.aal, label.fa, label.md, tract.fa, tract.md))
#mf.vbm.bn246.dti = scale(cbind(spm.vbm.bn246, label.fa, label.md, tract.fa, tract.md))
#
#mf.resting = scale(cbind(alff, reho))
#mf.resting.aal = scale(cbind(alff.aal, reho.aal, falff.aal))
#mf.resting.bn246 = scale(cbind(alff.bn246, reho.bn246, falff.bn246))
#
#mf.vbm.fa = scale(cbind(spm.vbm, label.fa, tract.fa))
#mf.vbm.aal.fa = scale(cbind(spm.vbm.aal, label.fa, tract.fa))
#mf.vbm.bn246.fa = scale(cbind(spm.vbm.bn246, label.fa, tract.fa))
#
#mf.vbm.resting = scale(cbind(spm.vbm, alff, reho))
#mf.vbm.resting.aal = scale(cbind(spm.vbm.aal, alff.aal, reho.aal))
#mf.vbm.resting.bn246 = scale(cbind(spm.vbm.bn246, alff.aal, reho.aal))
#
#feature.name = c("spm.vbm", "spm.vbm.aal", "spm.vbm.bn246", "alff", "alff.aal", "alff.bn246", "falff", "falff.aal", "falff.bn246", "reho", "reho.aal", "reho.bn246", "label.fa", "tract.fa", "dti.fa", "multimodal.feature", "multimodal.feature.aal", "multimodal.feature.bn246", "mf.dti", "mf.vbm.dti", "mf.vbm.aal.dti", "fm.vbm.bn246.dti", "mf.resting", "mf.resting.aal", "mf.resting.bn246", "mf.vbm.fa", "mf.vbm.aal.fa", "mf.vbm.bn246.fa", "mf.vbm.resting", "mf.vbm.resting.aal", "mf.vbm.resting.bn246")
#feature.list = list(spm.vbm, spm.vbm.aal, spm.vbm.bn246, alff, alff.aal, alff.bn246, falff, falff.aal, falff.bn246, reho, reho.aal, reho.bn246, label.fa, tract.fa, dti.fa, multimodal.feature, multimodal.feature.aal, multimodal.feature.bn246, mf.dti, mf.vbm.dti, mf.vbm.aal.dti, mf.vbm.bn246.dti, mf.resting, mf.resting.aal, mf.resting.bn246, mf.vbm.fa, mf.vbm.aal.fa, mf.vbm.bn246.fa, mf.vbm.resting, mf.vbm.resting.aal, mf.vbm.resting.bn246)
#

#feature.name = feature.name[c(14,16,20)]
#feature.list = feature.list[c(14,16,20)]

feature.name = "fc"
feature.list = list(fc)


report.rows = length(feature.list)

report  =  data.frame(group = feature.name, 
	Sex.acc = rep(NA, report.rows), Sex.sensi = rep(NA, report.rows), Sex.speci = rep(NA, report.rows)
	)

report.sd  =  data.frame(group = feature.name, 
	Sex.acc = rep(NA, report.rows), Sex.sensi = rep(NA, report.rows), Sex.speci = rep(NA, report.rows)
	)

for (i.feature in 1:length(feature.list)) {
	
	print(" running on feature:")
	print(feature.name[i.feature])

	brain.feature = scale(feature.list[[i.feature]])
	df.all = cbind(subject.info[,-1], brain.feature)
	
	df.subset = subset(df.all, select=-ptsd)
	print("dimension for subset dataset")
	print(dim(df.subset))
	print(table(df.subset$Sex))
	
	if (report.name=="glmnet_"){
		x = model.matrix(df.subset$Sex~., df.subset)
		y = df.subset$Sex
		cv.result = glmnet.cv.fun(x[,-1], y, lambda.seq = lambda.seq, glmnet.para)
		num.feature.tune.name = paste("feature.weights_Sex0vs1", file.name.tail, feature.name[i.feature], '.csv', sep="_")
		write.table(cv.result[[3]], num.feature.tune.name, sep=",", row.names=F)
	}
	else if (report.name=="svm_"){
		subject.info.subset = subject.info[subset.idx,c(-1)]
		# change the column name of the group index to meet convention of svm.cv.fun:
		colnames(subject.info.subset)[names(subject.info.subset)=="Sex"]="factor"
		cv.result = svm.cv.fun(brain.feature[subset.idx,], subject.info.subset, cost.seq, svm.para)
		num.feature.tune.name = paste("num.feature.tune_Sex0vs1", file.name.tail, feature.name[i.feature], '.csv', sep="_")
		write.table(cv.result[[4]], num.feature.tune.name, sep=",", row.names=F)
	}
	# save result:
	result = cv.result[[1]]
	train.result = cv.result[[2]]
	
	print(result)
	print(train.result)
	
	report[i.feature,2:4] = colMeans(result,na.rm=T)
	report.sd[i.feature,2:4] = apply(result, 2, function(x)sd(na.omit(x)))
	
	
}

# -------------------- save results : -------------------------------------
print(report)
#print(report.sd)

filename = paste("result_m_",report.name, file.name.tail, Sys.time(), ".csv", sep = "")
write.table(report, filename, sep = ",", row.names = F)
#filename = paste("result_sd_", file.name.tail, Sys.time(), ".csv", sep = "")
#write.table(report.sd, filename, sep = ",", row.names = F)

