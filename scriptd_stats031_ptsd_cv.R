#!/usr/bin/env Rscript
# Xin Niu May.4.2017
options(error = traceback)

rm(list = ls())
graphics.off()

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
report.name = "glmnet" # this can only be glmnet or svm, DO NOT COSTOMIZE THE NAME.

# glmnet parameters: all the parameters will be saved in text file as log.
glmnet.para = list(nfolds = 5,
                   nfolds.inner = 5,
                   family = "binomial")

#glmnet.para$type.measure = "class"
glmnet.para$type.measure = "auc"
#glmnet.para$type.measure = "mse"
glmnet.para$predict.type = "class"
#glmnet.para$lambda.seq = 10^seq(-3,3, length = 90)
glmnet.para$lambda.seq = "default" # default: default lambda sequence, 0: no regularization.

#"mean.diff.boot" #"mean.diff" #"cv.wilcox" #"none"
glmnet.para$penalty.weight = "mean.diff"
#glmnet.para$penalty.weight = "pearson.boot"
#glmnet.para$penalty.weight = "glmnet.coef.boot"
#glmnet.para$penalty.weight = "none"

glmnet.para$log.penalty.weight = 0 # log transform the penalty weight
glmnet.para$alpha = 1 # 1: lasso, 0: ridge regression, 0.5 elastic net
#glmnet.para$alpha = seq(0,1,0.1) # input a sequence of alpha to run cross validation.
glmnet.para$quantile.thresh = 1
glmnet.para$return.mod = F # return regression model trained on all the data
# note is not used in the functions but useful for keeping information about changes of other parts of the codes.
#glmnet.para$note = "running standard lasso without setting penalty weights. rather than (all 1 for penalty weights."

if (length(glmnet.para$alpha) == 1) {
    alpha.name = toString(glmnet.para$alpha)
} else{
    alpha.name = "tuning"
}

file.name.tail = paste(
    "no_reg_folds",
    toString(glmnet.para$nfolds),
    toString(glmnet.para$nfolds.inner),
    "alpha",
    alpha.name,
    glmnet.para$penalty.weight,
    'logtransform',
    toString(glmnet.para$log.penalty.weight),
    sep = "_"
)


result_dir = paste("result_d03", report.name, file.name.tail, Sys.Date(), sep = "_")
dir.create(result_dir)

library(ggplot2)
library(reshape2)

#### read features:
source("scriptd_stats01_1_read_all_feature2.R")
#source("scriptd_stats01_1_read_feature2.R")
#source("scriptd_stats01_read_feature_fc.R")
#source("scriptd_stats01_read_feature_fc.R")
source("scriptd_stats02_cv_functions.R")

#### construct multimodal features:
# feature.name and feature.list will be constructed:
#source("scriptd_stats01_2_multimodal_feature_bn246.R")
source("scriptd_stats01_2_multimodal_feature_cat12.R")

# delete subjects with severe head motion or other problems:
subject.delete = read.table(paste(data_dir, slash, "subject_delete.txt", sep = ""), 
                            header = F)
print("deleted subjects:")
print(subject.delete)
delete.idx = subject.info$SUBJID %in% subject.delete$V1
#delete.idx = rep(FALSE, length(delete.idx))

# save deleted subjects id as log:
glmnet.para$delete.subjects = subject.delete

print("ptsd:")
print(table(subject.info[delete.idx, ]$ptsd))
print("Sex:")
print(table(subject.info[delete.idx, ]$Sex))

# add no brain option in feature list, this does not work for glmnet model with weighted penalty factors.
#feature.name = c("no.brain", feature.name)
#feature.list = c(list(NULL), feature.list)


# prepare reports to fill:
report.rows = length(feature.name)

report  =  data.frame(
    group = feature.name,
    Hc_Trauma.acc = rep(NA, report.rows),
    Hc_Trauma.sensi = rep(NA, report.rows),
    Hc_Trauma.speci = rep(NA, report.rows),
    Hc_Trauma.R.index = rep(NA, report.rows),
    Hc_Ptsd.acc = rep(NA, report.rows),
    Hc_Ptsd.sensi = rep(NA, report.rows),
    Hc_Ptsd.speci = rep(NA, report.rows),
    Hc_ptsd.R.index = rep(NA, report.rows),
    Trauma_Ptsd.acc = rep(NA, report.rows),
    Trauma_Ptsd.sensi = rep(NA, report.rows),
    Trauma_Ptsd.speci = rep(NA, report.rows),
    Trauma_ptsd.R.index = rep(NA, report.rows)
)

report.sd  =  data.frame(
    group = feature.name,
    Hc_Trauma.acc = rep(NA, report.rows),
    Hc_Trauma.sensi = rep(NA, report.rows),
    Hc_Trauma.speci = rep(NA, report.rows),
    Hc_Ptsd.acc = rep(NA, report.rows),
    Hc_Ptsd.sensi = rep(NA, report.rows),
    Hc_Ptsd.speci = rep(NA, report.rows),
    Trauma_Ptsd.acc = rep(NA, report.rows),
    Trauma_Ptsd.sensi = rep(NA, report.rows),
    Trauma_Ptsd.speci = rep(NA, report.rows)
)

subject.info = subject.info[!delete.idx,] # do not use -delete.idx here, it will just remove the 1st row!
subject.info$age = scale(subject.info$age)
print("subject info:")
print(dim(subject.info))
print(table(subject.info[, c("Sex", "ptsd")]))

for (i.feature in 1:length(feature.list)) {
    print(" running on feature:")
    print(feature.name[i.feature])
    
    if (is.null(feature.list[[i.feature]])) {
        print("null")
        df.all = subject.info[, -1]
    } else {
        brain.feature = scale(feature.list[[i.feature]][!delete.idx, ])
        #brain.feature = scale(feature.list[[i.feature]])
        df.all = cbind(subject.info[, -1], brain.feature)
    }
    
    # remove NAN columns:
    df.all = df.all[, apply(df.all, 2, function(x)
        ! any(is.na(x)))]
    
    report.col.idx = c(2:5)
    shift = 0
    
    for (i.ptsd1 in 0:1) {
        for (i.ptsd2 in (i.ptsd1 + 1):2) {
            subset.idx = df.all$ptsd == i.ptsd1 | df.all$ptsd == i.ptsd2
            df.subset = df.all[subset.idx, ]
            print("dimension for subset dataset")
            print(dim(df.subset))
            print(table(df.subset$ptsd))
            ptsd.comparision.name = paste("ptsd",
                                          toString(i.ptsd1),
                                          "vs",
                                          toString(i.ptsd2),
                                          sep = "_")
            
            if (report.name == "glmnet") {
                x = model.matrix(df.subset$ptsd ~ ., df.subset)
                y = df.subset$ptsd
                cv.result = glmnet.nested.cv(x[, -1], y, glmnet.para)
                
                coefs.name = paste(
                    result_dir,
                    slash,
                    "coefs_",
                    ptsd.comparision.name,
                    "_",
                    file.name.tail,
                    "_",
                    feature.name[i.feature],
                    '.csv',
                    sep = ""
                )
                
                write.table(cv.result$coefs,
                            coefs.name,
                            sep = ",",
                            row.names = F)
            }
            else if (report.name == "svm") {
                subject.info.subset = subject.info[subset.idx, c(-1)]
                
                # change the column name of the group index to meet convention of svm.cv.fun:
                colnames(subject.info.subset)[names(subject.info.subset) == "ptsd"] =
                    "factor"
                
                cv.result = svm.cv.fun(brain.feature[subset.idx, ],
                                       subject.info.subset,
                                       cost.seq,
                                       svm.para)
                
                num.feature.tune.name = paste(
                    "num.feature.tune_",
                    ptsd.comparision.name,
                    file.name.tail,
                    feature.name[i.feature],
                    '.csv',
                    sep = "_"
                )
                write.table(
                    cv.result[[4]],
                    num.feature.tune.name,
                    sep = ",",
                    row.names = F
                )
            }
            # save result:
            test.result = cv.result$test.result
            train.result = cv.result$train.result
            R.index = cv.result$reproducibility
            print(test.result)
            print(train.result)
            
            report[i.feature, report.col.idx + shift] = c(colMeans(test.result, na.rm =
                                                                       T), R.index)
            #report.sd[i.feature, report.col.idx + shift] = apply(test.result, 2, function(x)sd(na.omit(x)))
            shift = shift + 4
            
        }
    }
}

# -------------------- save results : -------------------------------------
print(report)
#print(report.sd)

time = format(Sys.time(), '%Y_%b_%d_%H_%M_%S')
filename = paste(result_dir,
                 slash,
                 report.name,
                 file.name.tail,
                 time,
                 ".csv",
                 sep = "")

write.table(report, filename, sep = ",", row.names = F)

# save parameter to text file:
sink(paste(result_dir, slash, "para", time, ".txt", sep = ""))
print(glmnet.para)
sink()
