#!/usr/bin/env Rscript
# Xin Niu May.4.2017

# udpates 1/23/2021 add age of trauma to features:

# udpates 2/26/2021 save orignal accuracy for each folds

rm(list=ls())
graphics.off()
library(plyr) #ldply

#setwd("home/xin/Dropbox/weighted_lasso_ptsd_brain")
setwd("C:/Users/Xin/Dropbox/weighted_lasso_ptsd_brain")


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

# glmnet parameters:
report.name = "glmnet" # this can only be glmnet or svm, DO NOT COSTOMIZE THE NAME.

# glmnet parameters: all the parameters will be saved in text file as log.
glmnet.para = list(nfolds = 10,
                   nfolds.inner = 5,
                   family = "binomial")

#glmnet.para$type.measure = "class"
glmnet.para$type.measure = "auc"
#glmnet.para$type.measure = "mse"

glmnet.para$predict.type = "class"
glmnet.para$lambda.seq = 10 ^ seq(-3, 3, length = 90)
#glmnet.para$lambda.seq = "default" # default: default lambda sequence, 0: no regularization.

glmnet.para$alpha = .5 # 1: lasso, 0: ridge regression, 0.5 elastic net
#glmnet.para$alpha = seq(0,1,0.1) # input a sequence of alpha to run cross validation.
glmnet.para$return.mod = F # return regression model trained on all the data

#"mean.diff.boot" #"mean.diff" #"cv.wilcox" #"none"
#glmnet.para$penalty.weight = "mean.diff.boot"
#glmnet.para$penalty.weight = "kendall"
#glmnet.para$penalty.weight = "pearson.boot"
#glmnet.para$penalty.weight = "glmnet.coef.boot"
#glmnet.para$penalty.weight = "none"

glmnet.para$cutoff = c(0.2,.8) # threshhold value below 1st element to 0 and larger than 2nd element to 1.
glmnet.para$log.penalty.weight = 1 # log transform the penalty weight
#glmnet.para$quantile.thresh = 1 # feature select according to penalty weight.


glmnet.para1 = glmnet.para
glmnet.para1$penalty.weight = "none"
glmnet.para1$alpha = 1

glmnet.para2 = glmnet.para
glmnet.para2$penalty.weight = "none"
glmnet.para2$alpha = .5

glmnet.para3 = glmnet.para
glmnet.para3$penalty.weight = "mean.diff.boot"
glmnet.para3$alpha = 1

glmnet.para4 = glmnet.para
glmnet.para4$penalty.weight = "mean.diff.boot"
glmnet.para4$alpha = .5


batch.list = list(glmnet.para1, glmnet.para2, glmnet.para3, glmnet.para4)
#batch.list = list(glmnet.para1)

for (i.batch in 1:length(batch.list)) {
    
    glmnet.para = batch.list[[i.batch]]

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
    
    
    result_dir = paste("result_d03_trauma_age_cv", report.name, file.name.tail, Sys.Date(), sep = "_")
    
    dir.create(result_dir)
    
    #library(reshape2)
    #source("scriptd_stats01_1_read_feature2.R")
    source("scriptd_stats01_1_read_all_feature2.R")
    #source("scriptd_stats01_1_read_feature_fc.R")
    source("scriptd_stats02_cv_functions2.R")
    
    # construct multimodal features:
    # feature.name and feature.list will be constructed:
    #source("scriptd_stats01_2_multimodal_feature_bn246.R")
    #source("scriptd_stats01_2_multimodal_feature.R")
    source("scriptd_stats01_2_multimodal_feature_cat12.R")
    #source("scriptd_stats01_2_multimodal_feature_fc.R")
    
    report.rows = length(feature.list)
    print(report.rows)
    report  =  data.frame(group = rep(NA, report.rows), 
                          accuracy = rep(NA,report.rows), sensitivity = rep(NA, report.rows), specificity = rep(NA, report.rows), 
                          std = rep(NA, report.rows), Reproducibility = rep(NA, report.rows), robustness = rep(NA, report.rows),
                          Reproducibility.nonzero = rep(NA, report.rows)
    )
    
    
    report.acc.nfolds = list()
    # choose ptsd type to classify: 0: hc, 1: trauma, 2: ptsd
    #i.ptsd1 = c(0,1) # healthy control and trauma
    i.ptsd1 = 1 # trauma
    i.ptsd2 = 2 # ptsd
    
    colnames(report)[1] = paste(toString(i.ptsd1), "vs", toString(i.ptsd2))
    
    for (i.feature in 1:length(feature.list)) {
        
        
        
        print(" running on feature:")
        feature.name = feature.list[[i.feature]][[2]]
        print(feature.name)
        report[i.feature, 1] = feature.name
        
        acc.nfolds = data.frame(group = rep(feature.name, glmnet.para$nfolds), 
                                cv = 1:glmnet.para$nfolds, 
                                acc = rep(NA, glmnet.para$nfolds))
        
        brain.feature = scale(feature.list[[i.feature]][[1]])
        
        if (dim(brain.feature)[1]!=0){
            df.all = cbind(subject.info[,-1], brain.feature)
        } else {
            df.all = subject.info[, -1]
        }
        
        subset.idx = df.all$ptsd %in% i.ptsd1|df.all$ptsd %in% i.ptsd2
        df.subset = df.all[subset.idx,]
        
        
        # remove columns with NaN values
        nan.col = apply(df.subset, 2, function(x) any(is.nan(x)))
        df.subset = df.subset[, !nan.col]
        
        df.subset = df.subset[complete.cases(df.subset),]
        
        #df.subset$factor = df.subset$ptsd %in% i.ptsd1	
        factor = df.subset$ptsd %in% i.ptsd1	
        
        print("dimension for subset dataset")
        print(dim(df.subset))
        print(table(df.subset$ptsd))
        ptsd.comparision.name = paste("ptsd", toString(i.ptsd1), "vs", toString(i.ptsd2), sep = "_")
        
        #print(head(df.subset))
        if (report.name=="glmnet"){
            #x = as.matrix(subset(df.subset, select = -c(ptsd, Sex)))
            #x = as.matrix(subset(df.subset, select = -c(ptsd)))
            #x = data.matrix(subset(df.subset, select = -c(ptsd)))
            #print(head(x))
            
            x = model.matrix(df.subset$ptsd~., df.subset)[,-1]
            y = as.numeric(factor)
            #print(head(x))
            #print(head(y))
            

            cv.result = glmnet.nested.cv(x, y, glmnet.para)
            
            coefs.name = paste0(result_dir, slash, "coefsidx_", ptsd.comparision.name, "_", file.name.tail, "_", feature.name, '.csv')
            write.table(cv.result$coefs.ind, coefs.name, sep = ",", row.names = F)
            coefs.name = paste0(result_dir, slash, "coefs_", ptsd.comparision.name, "_", file.name.tail, "_", feature.name, '.csv')
            write.table(cv.result$coefs, coefs.name, sep = ",", row.names = F)
        }
        else if (report.name=="svm_"){
            subject.info.subset = subset(subject.info, select = -c(SUBJID, ptsd))
            # change the column name of the group index to meet convention of svm.cv.fun:
            #colnames(subject.info.subset)[names(subject.info.subset)=="ptsd"]="factor"
            cv.result = svm.cv.fun(brain.feature[subset.idx,], subject.info.subset, cost.seq, svm.para)
            num.feature.tune.name = paste("num.feature.tune_", ptsd.comparision.name, file.name.tail, feature.name[i.feature], '.csv', sep="_")
            write.table(cv.result[[5]], num.feature.tune.name, sep=",", row.names=F)
        }
        # save result:
        test.result = cv.result$test.result
        train.result = cv.result$train.result
        print(test.result)
        print(train.result)
        #print(coefs.result)
        
        report[i.feature, 2:4] = apply(test.result, 2, mean)
        report[i.feature, 5] = sd(test.result[,1], na.rm=T)
        report[i.feature, 6] = cv.result$reproducibility
        report[i.feature, 7] = cv.result$coefs.robustness
        report[i.feature, 8] = cv.result$reproducibility.nonzero
        
        acc.nfolds[, 'acc'] = test.result[,1] # first column of test.result is accuracy.
        report.acc.nfolds = c(report.acc.nfolds, list(acc.nfolds))
        
        
    }
    
    # -------------------- save results : -------------------------------------
    print(report)
    #print(report.sd)
    
    time = format(Sys.time(), '%Y_%b_%d_%H_%M_%S')
    filename = paste(result_dir,
                     slash,
                     report.name,
                     file.name.tail,
                     #time,
                     ".csv",
                     sep = "")
    
    write.table(report, filename, sep = ",", row.names = F)
    
    filename = paste(result_dir,
                     slash,
                     report.name,
                     file.name.tail,
                     #time,
                     "acc.csv",
                     sep = "")
    
    report.acc = ldply(report.acc.nfolds, rbind)
    write.table(report.acc, filename, sep = ',', row.names = F)
}



