#!/usr/bin/env Rscript
# Xin Niu May.4.2017


rm(list = ls())
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
report.name = "glmnet" # this can only be glmnet or svm, DO NOT COSTOMIZE THE NAME.

# glmnet parameters: all the parameters will be saved in text file as log.
glmnet.para = list(nfolds = 5,
                   nfolds.inner = 5,
                   family = "binomial")

glmnet.para$type.measure = "class"
#glmnet.para$type.measure = "auc"
#glmnet.para$type.measure = "mse"
glmnet.para$predict.type = "class"
glmnet.para$lambda.seq = 10 ^ seq(-3, 3, length = 90)
#glmnet.para$lambda.seq = "default" # default: default lambda sequence, 0: no regularization.

glmnet.para$alpha = 1 # 1: lasso, 0: ridge regression, 0.5 elastic net
#glmnet.para$alpha = seq(0,1,0.1) # input a sequence of alpha to run cross validation.
glmnet.para$return.mod = F # return regression model trained on all the data

#"mean.diff.boot" #"mean.diff" #"cv.wilcox" #"none"
glmnet.para$penalty.weight = "mean.diff.boot"
#glmnet.para$penalty.weight = "pearson.boot"
#glmnet.para$penalty.weight = "glmnet.coef.boot"
#glmnet.para$penalty.weight = "none"
glmnet.para$cutoff = c(0,1)
glmnet.para$log.penalty.weight = 1 # log transform the penalty weight
glmnet.para$quantile.thresh = 1 # feature select according to penalty weight.




# note is not used in the functions but useful for keeping information about changes of other parts of the codes.
#glmnet.para$note = "running standard lasso without setting penalty weights. rather than (all 1 for penalty weights."

if (length(glmnet.para$alpha) == 1) {
    alpha.name = toString(glmnet.para$alpha)
} else{
    alpha.name = "tuning"
}

file.name.tail = paste(
    "folds",
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

#library(reshape2)
source("scriptd_stats01_1_read_all_feature.R")
source("scriptd_stats01_2_multimodal_feature.R")
source("scriptd_stats02_cv_functions.R")

subject.info = subject.info[, c("SUBJID", "Gender", "age_at_cnb")]
subject.info$age_at_cnb = scale(subject.info$age_at_cnb)

report.rows = length(feature.list)

report  =  data.frame(
    group = rep(NA, report.rows),
    Gender.acc = rep(NA, report.rows),
    Gender.sensi = rep(NA, report.rows),
    Gender.speci = rep(NA, report.rows),
    Gender.reproducibility = rep(NA, report.rows),
    Gender.robustness = rep(NA, report.rows)
)

report.sd  =  data.frame(
    group = rep(NA, report.rows),
    Gender.acc = rep(NA, report.rows),
    Gender.sensi = rep(NA, report.rows),
    Gender.speci = rep(NA, report.rows)
)

for (i.feature in 1:length(feature.list)) {
    print(" running on feature:")
    feature.name = feature.list[[i.feature]][[2]]
    print(feature.name)
    report[i.feature, "group"] = feature.name
    report.sd[i.feature, "group"] = feature.name
    
    brain.feature = cbind(cat.vbm.subid, feature.list[[i.feature]][[1]])
    df.all = merge(subject.info, brain.feature, by = "SUBJID")
    
    df.subset = subset(df.all, select = -SUBJID)
    # model.matrix remove rows with Nans, which makes x has fewer rows than y.
    #df.subset = scale(df.subset[complete.cases(df.subset),])
    df.subset = df.subset[complete.cases(df.subset), ]
    print("dimension for subset dataset")
    print(dim(df.subset))
    print(table(df.subset$Gender))
    
    #print(head(df.subset))
    
    if (report.name == "glmnet") {
        x = model.matrix(df.subset$Gender ~ ., df.subset)
        #x = scale(subset(df.subset, select = -Gender))
        y = df.subset$Gender
        print(any(is.nan(y)))
        print(dim(x))
        print(length(y))
        
        cv.result = glmnet.nested.cv(scale(x[, -1]), y, glmnet.para)
        
        coefs.name = paste(
            result_dir,
            slash,
            "coefs_gender",
            "_",
            file.name.tail,
            "_",
            feature.name,
            '.csv',
            sep = ""
        )
        
        write.table(cv.result$coefs,
                    coefs.name,
                    sep = ",",
                    row.names = F)
        
        coefs.ind.name = paste(
            result_dir,
            slash,
            "coefs_ind_gender",
            "_",
            file.name.tail,
            "_",
            feature.name,
            '.csv',
            sep = ""
        )
        
        write.table(cv.result$coefs.ind,
                    coefs.ind.name,
                    sep = ",",
                    row.names = F)
    }
    
    else if (report.name == "svm") {
        
        subject.info.subset = subject.info[subset.idx, c(-1)]
        
        # change the column name of the group index to meet convention of svm.cv.fun:
        colnames(subject.info.subset)[names(subject.info.subset) == "Gender"] =
            "factor"
        
        cv.result = svm.cv.fun(brain.feature[subset.idx, ],
                               subject.info.subset,
                               cost.seq,
                               svm.para)
        
        num.feature.tune.name = paste("num.feature.tune_Gender0vs1",
                                      file.name.tail,
                                      feature.name,
                                      '.csv',
                                      sep = "_")
        
        write.table(cv.result[[4]],
                    num.feature.tune.name,
                    sep = ",",
                    row.names = F)
        
    }
    # save result:
    result = cv.result$test.result
    train.result = cv.result$train.result
    
    print(result)
    print(train.result)
    
    report[i.feature, 2:4] = colMeans(result, na.rm = T)
    report[i.feature, 5] = cv.result$reproducibility
    report[i.feature, 6] = cv.result$coefs.robustness
    report[i.feature, 6] = cv.result$coefs.robustness
    
    report.sd[i.feature, 2:4] = apply(result, 2, function(x)
        sd(na.omit(x)))
    
    
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
