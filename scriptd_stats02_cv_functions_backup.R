#!/usr/bin/env Rscript

library(caret)
# createFolds()

scale.0.1 = function(dat) {
    # the output will be coerced to matrix.
    
    dat = as.matrix(dat)
    
    mins = apply(dat, 2, min)
    maxs = apply(dat, 2, max)
    
    scaled.dat = scale(dat, center = mins, scale = maxs - mins)
    return(scaled.dat)
}

compute.acc = function(y, yhat) {
    y = as.numeric(y)
    yhat = as.numeric(yhat)
    
    acc <- sum(y == yhat) / length(y)
    
    ylevel = sort(unique(y), decreasing = F)
    
    if (length(ylevel) == 2) {
        # asuming ylevel = c(0, 1)
        sensi <-
            sum(y == yhat & y == ylevel[2]) / sum(y == ylevel[2])
        speci <-
            sum(y == yhat & y == ylevel[1]) / sum(y == ylevel[1])
    }
    else if (max(yhat) == max(y)) {
        sensi <- sum(y == yhat & y == ylevel) / sum(y == ylevel)
        speci <- NaN
    } else{
        speci <- sum(y == yhat & y == ylevel) / sum(y == ylevel)
        sensi <- NaN
    }
    
    temp <- c(acc, sensi, speci)
    return(temp)
}

remap.factor = function(f, min = 0, max = 1) {
    f = as.numeric(f)
    
    f.remap = f
    f.remap[f == max(f)] = max
    f.remap[f == min(f)] = min
    return(f.remap)
}


perm.t.test = function(x1, x2, n = 5000, paired = F) {
    # permutation t test for paired and independent samples:
    
    mean.diff = mean(x1) - mean(x2)
    perm.stats = rep(NA, n)
    test.out = list()
    
    if (paired) {
        # paired permutation t-test: corresponding pairs in x1 and x2 were randomly
        # shuffled across the columns of cbind(x1, x2), and the mean differenced were
        # computed.
        x = cbind(x1, x2)
        n.obs = dim(x)[1]
        for (i in 1:n) {
            x1.perm.idx = sample(1:2, n.obs, replace = T)
            x2.perm.idx = 2 - x1.perm.idx
            perm.stats[i] = mean(x[cbind(1:n.obs, x1.perm.idx)] - x[cbind(1:n.obs, x2.perm.idx)])
        }
    } else{
        # independent permutation t-test:
        x = c(x1, x2)
        n.obs = length(x)
        for (i in 1:n) {
            x.perm.idx = sample(1:2, n.obs, replace = T)
            perm.stats[i] = mean(x[x.perm.idx == 1]) - mean(x[x.perm.idx == 2])
        }
    }
    
    test.out$density = perm.stats
    test.out$p.greater = sum(perm.stats > mean.diff) / n
    test.out$p.smaller = sum(perm.stats < mean.diff) / n
    test.out$p = sum(perm.stats < -abs(mean.diff) |
                         perm.stats > abs(mean.diff)) / n
    
    return(test.out)
}


library(VGAM)
# t.test()
# cor.test()

library(ltm)
# biserial.cor()

#library(matrixcalc)
## matrix.inverse()

#library(MASS)
# ginv()

feature.weight.mean.diff.boot = function(x, y, nboots, cutoff = 0){
    
    f.w.boot = matrix(NA, nboots, ncol(x))
    
    for (i in 1:nboots) {
        set.seed(i)
        boot.idx = sample(1:nrow(x), nrow(x), replace = T)
        x.boot = x[boot.idx,]
        y.boot = y[boot.idx]
        
        #f.w.boot[i, ] = abs(apply(x.boot[y.boot==0,], 2, mean) - apply(x.boot[y.boot==1,], 2, mean))
        f.w.boot[i, ] = apply(x.boot[y.boot==0,], 2, mean) - apply(x.boot[y.boot==1,], 2, mean)
    }
    
    feature.weight = scale.0.1(apply(f.w.boot, 2, sd) / apply(f.w.boot, 2, mean))
    
    if (cutoff>0){
        cutoff.value = quantile(feature.weight, cutoff)
        # when we set this to 0, we have convergence problem, wo we set it to the min value of feature.weight.
        feature.weight[feature.weight<cutoff.value] = min(feature.weight) 
    }
    return(feature.weight)
    
}

feature.cv.test.pca = function(feature.in,
                               factor,
                               k = 10,
                               ncomp = 10,
                               method = "wilcox") {
    
    pca.out = prcomp(feature.in, scale = F, center = F)
    feature.pca = pca.out$x
    
    feature.weight.pca = feature.cv.test(feature.pca, factor, k, method)
    #inv.rotation = ginv(pca.out$rotation[, 1:length(feature.weight.pca)])
    inv.rotation = ginv(pca.out$rotation[, 1:ncomp])
    feature.weight = feature.weight.pca[1:ncomp] %*% inv.rotation
    
    return(feature.weight)
}

feature.cv.test = function(feature.in,
                           factor,
                           k,
                           method = "wilcox",
                           glmnet.para,
                           seed = 111) {
    # statistic tests to compute feature weights: if k == 1, run test on the whole
    # sample without cross validation. the output is a vector with length same as
    # the number of features i.e ncol(feature.in) this is useful in case we do boot
    # strapping and compute coefficient of variation (CV) across boot strap samples.
    
    num.feature = ncol(feature.in)
    num.sample = nrow(feature.in)
    set.seed(seed)
    
    factor = as.numeric(factor)
    idx.factor.0 = factor == min(factor)
    idx.factor.1 = factor == max(factor)
    
    if (k > 1) {
        cv.k = createFolds(factor, k, list = F)
        print("feature.cv.test: number of samples in each CV, for factor 0 and 1:")
        print(table(cv.k[idx.factor.0]))
        print(table(cv.k[idx.factor.1]))
    } else{
        # k is set 1 in bootstrap, in which all bootstrapped sample are selected
        # setting cross-validation index as all 0 will enable all samples being selected
        # in the for loop over 1 to k. feature.in[cv.k != i,]
        cv.k = rep(0, num.sample)
    }
    
    test.out = matrix(NA, k, num.feature)
    
    for (i in 1:k) {
        idx = cv.k != i
        
        if (method == "mean.diff") {
            test.out[i, ] = abs(apply(feature.in[cv.k != i &
                                                     idx.factor.0,], 2, mean) -
                                    apply(feature.in[cv.k != i &
                                                         idx.factor.1,], 2, mean))
        } else if (method == "glmnet.coef") {
            # run cv.glmnet to get the coefficients and use the CV of them as feature weights.
            set.seed(444)
            cv.fit = cv.glmnet(
                feature.in[cv.k != i, ],
                factor[cv.k != i],
                nfold = glmnet.para$nfolds.inner,
                alpha = glmnet.para$alpha,
                family = glmnet.para$family,
                standardize = F
            )
            
            #test.out[i,] = abs(coef(cv.fit, s="lambda.min"))[-1]
            test.out[i, ] = coef(cv.fit, s = "lambda.min")[-1]
        } else{
            for (i.feature in 1:num.feature) {
                if (method == "kendall") {
                    test.result = cor.test(feature.in[idx, i.feature], factor[dix], method = "kendall")
                    value = abs(test.result$estimate)
                } else if (method == "wilcox") {
                    test.result = wilcox.test(feature.in[idx, i.feature], factor[idx])
                    value = abs(test.result$estimate)
                } else if (method == "spearman") {
                    test.result = cor.test(feature.in[idx, i.feature], factor[idx], method = "spearman")
                    value = abs(test.result$estimate)
                } else if (method == "pearson") {
                    test.result = cor.test(feature.in[idx, i.feature], factor[idx])
                    value = abs(test.result$estimate)
                } else if (method == "biserial") {
                    value = abs(biserial.cor(feature.in[idx, i.feature], factor[idx]))
                }
                test.out[i, i.feature] = value
            }
        }
    }
    return(test.out)
}

feature.cv.boot = function(feature.in,
                           factor,
                           n = 100,
                           method,
                           glmnet.para,
                           pca = F) {
    
    f.cv.boot = matrix(NA, n, ncol(feature.in))
    factor = remap.factor(factor)
    
    for (i in 1:n) {
        set.seed(i)
        boot.idx = sample(1:nrow(feature.in), nrow(feature.in), replace = T)
        feature.boot = feature.in[boot.idx,]
        factor.boot = factor[boot.idx]
        if (pca) {
            f.cv.boot[i, ] = feature.cv.test.pca(feature.boot, factor.boot, k = 1, method)
        } else{
            f.cv.boot[i, ] = feature.cv.test(feature.boot, factor.boot, k = 1, method, glmnet.para)
        }
    }
    return(f.cv.boot)
}

compute.feature.weight = function(x.train, y.train, method = "mean.diff", glmnet.para) {
    penalty.weight = method
    
    if (penalty.weight == "none") {
        f.weight.cv = rep(1, dim(x.train)[2])
    } else {
        if (penalty.weight == "mean.diff") {
            f.weight = feature.cv.test(x.train, y.train, 10, "mean.diff", glmnet.para)
        } else if (penalty.weight == "wilcox") {
            f.weight = feature.cv.test(x.train, y.train, 10, "wilcox", glmnet.para)
        } else if (penalty.weight == "pearson") {
            f.weight = feature.cv.test(x.train, y.train, 10, "pearson", glmnet.para)
        } else if (penalty.weight == "mean.diff.boot") {
            f.weight = feature.cv.boot(x.train,
                                       y.train,
                                       n = 500,
                                       method = "mean.diff",
                                       glmnet.para)
        } else if (penalty.weight == "pearson.boot") {
            f.weight = feature.cv.boot(x.train,
                                       y.train,
                                       n = 500,
                                       method = "pearson",
                                       glmnet.para)
        } else if (penalty.weight == "biserial.boot") {
            f.weight = feature.cv.boot(x.train,
                                       y.train,
                                       n = 500,
                                       method = "biserial",
                                       glmnet.para)
        } else if (penalty.weight == "glmnet.coef.boot") {
            f.weight = feature.cv.boot(x.train,
                                       y.train,
                                       n = 500,
                                       method = "glmnet.coef",
                                       glmnet.para)
        } else if (penalty.weight == "mean.diff.boot.pca") {
            f.weight = feature.cv.boot(
                x.train,
                y.train,
                n = 500,
                method = "mean.diff",
                glmnet.para,
                pca = T
            )
        }
        # compute coefficient of variation:
        f.weight.cv = scale.0.1(apply(f.weight, 2, sd) / apply(f.weight, 2, mean))
        # log transform the weights:
        
        if (glmnet.para$log.penalty.weight) {
            print('log transform weight:')
            f.weight.cv = log(f.weight.cv + 1)
        }
        
    }
    return(f.weight.cv)
}


library(glmnet)
#library(SIS)

glmnet.tune = function(x,
                       y,
                       k,
                       alpha = 1,
                       lambda.seq = 10 ^ seq(-4, 3, length = 70),
                       f.weights) {
    # use balanced accuracy (mean of sensitivity and specificity) rather than accuracy to run cross-validation:
    
    #set.seed(222)
    set.seed(123)
    
    if (missing(f.weights)) {
        f.weights = rep(1, dim(x)[2])
    }
    
    cv.k = createFolds(y, k, list = F)
    test.acc = matrix(NA, k, length(lambda.seq))
    
    glmnet.control(mxit = 1000000)
    
    k = max(c(5, k))
    
    for (i in 1:k) {
        x.train = x[cv.k != i, ]
        x.test = x[cv.k == i, ]
        y.train = y[cv.k != i]
        y.test = y[cv.k == i]
        
        for (j in 1:length(lambda.seq)) {
            mod = glmnet(
                x.train,
                y.train,
                family = "binomial",
                alpha = alpha,
                lambda = lambda.seq[j],
                standardize = F,
                penalty.factor = f.weights
            )
            y.pred = predict(mod, x.test, s = lambda.seq[j], type = "class")
            acc = compute.acc(y.pred, y.test)
            
            test.acc[i, j] = (acc[2] + acc[3]) / 2
            #test.acc[i,j] = acc[1]
        }
    }
    
    acc.cv = apply(test.acc, 2, sd) / apply(test.acc, 2, mean)
    min.idx = which.min(acc.cv)
    return(lambda.seq[min.idx])
}


glmnet.nested.cv = function(x, y, glmnet.para) {
    k = glmnet.para$nfolds
    set.seed(222)
    cv.k = createFolds(y, k, list = F)
    
    if (glmnet.para$family == "binomial") {
        y.num = as.numeric(y)
        idx.factor.0 = y.num == min(y.num)
        idx.factor.1 = y.num == max(y.num)
        print("glmnet.nested.cv: number of samples in each CV, for factor 0 and 1:")
        print(table(cv.k[idx.factor.0]))
        print(table(cv.k[idx.factor.1]))
        
    }
    
    penalty.weight = glmnet.para$penalty.weight
    quantile.thresh = glmnet.para$quantile.thresh
    alpha = glmnet.para$alpha
    
    #print(t(rbind(cv.k, y)))
    
    test.result = data.frame(
        acc = rep(NA, k),
        sensi = rep(NA, k),
        speci = rep(NA, k)
    )
    train.result = data.frame(
        acc = rep(NA, k),
        sensi = rep(NA, k),
        speci = rep(NA, k)
    )
    
    result.coefs = matrix(NA, ncol(x) + 1, k)
    result.feature.weights = matrix(NA, ncol(x), k)
    
    # outer CV:
    for (i in 1:k) {
        print(length(cv.k))
        print(dim(x))
        x.train = x[cv.k != i, ]
        x.test = x[cv.k == i, ]
        y.train = y[cv.k != i]
        y.test = y[cv.k == i]
        
        glmnet.control(mxit = 1000000)
        
        if (length(glmnet.para$lambda.seq) == 1 &
            glmnet.para$lambda.seq == 0) {
            # no regularization, in this case, alpha is not necessary.
            glmnet.fit = glmnet(
                x.train,
                y.train,
                family = toString(glmnet.para$family),
                lambda = glmnet.para$lambda.seq,
                standardize = F
            )
            
            result.feature.weights = NULL
            
        } else {
            # with regularization.
            
            #print("computing feature weights:")
            #print(penalty.weight)
            
            f.weight.cv = compute.feature.weight(x.train, y.train, penalty.weight, glmnet.para)
            #print("feature penalty weights:")
            #print(f.weight.cv)
            
            error.list = rep(NA, length(alpha))
            fit.list = vector("list", length(alpha))
            
            if (glmnet.para$lambda.seq == "default") {
                # default lambda sequence of cv.glmnet:
                lambda.seq = formals(cv.glmnet)$lambda
                
            } else {
                # tune on lambda and alpha:
                if (length(glmnet.para$lambda.seq) == 1) {
                    # fixed lambda, just tune on alpha:
                    lambda.min = glmnet.para$lambda.seq
                } else if (length(glmnet.para$lambda.seq) > 1) {
                    # costomized lambda sequence:
                    lambda.seq = glmnet.para$lambda.seq
                }
            }
            #print(lambda.seq)
            
            
            # inner CV to optimized alpha:
            for (i.alpha in 1:length(alpha)) {
                print("alpha:")
                print(alpha[i.alpha])
                print(toString(glmnet.para$family))
                
                set.seed(123)
                # inner CV to optimize lambda:
                fit.list[[i.alpha]] = cv.glmnet(
                    x.train,
                    y.train,
                    nfolds = glmnet.para$nfolds.inner,
                    family = toString(glmnet.para$family),
                    # added on July 25 2018. all analysis before used the default lambda.
                    lambda = lambda.seq,
                    alpha = alpha[i.alpha],
                    penalty.factor = f.weight.cv,
                    type.measure = toString(glmnet.para$type.measure),
                    standardize = F
                )
            }
            
            error.list[i.alpha] = min(fit.list[[i.alpha]]$cvm)
            #print("tuning on alpha")
            #print(cbind(error.list, alpha))
            
            min.i.alpha = which.min(error.list)
            cv.fit = fit.list[[min.i.alpha]]
            lambda.min = cv.fit$lambda.min
            alpha.min = alpha[min.i.alpha]
            print(dim(result.feature.weights))
            print(length(f.weight.cv))
            result.feature.weights[, i] <- f.weight.cv
            
            #lambda.min = glmnet.tune(x.train, y.train, k, alpha, lambda.seq, f.cv)
            print("glmnet.cv.fun: min lambda:")
            print(lambda.min)
            
            # fit with all the training data, this can be obtained directly from output of cv.glmnet (out$glmnet.fit)
            # glmnet.fit = glmnet(x.train, y.train,
            #              	family = glmnet.para$family,
            #              	alpha=alpha[min.i.alpha],
            #              	penalty.factor = f.weight.cv,
            #              	lambda = lambda.min,
            #              	standardize = F)
            
            glmnet.fit = cv.fit$glmnet.fit
            coefs <- coef(cv.fit, s = "lambda.min")
            
        }# end with regularization.
        
        print(length(as.vector(coefs)))
        print(dim(result.coefs))
        result.coefs[, i] <- as.vector(coefs)
        
        #f.idx = f.weight.cv<=quantile(f.weight.cv, quantile.thresh)
        #cv.fit = cv.glmnet(x.train[, f.idx], y.train, nfolds = k, family = "binomial", penalty.factor = f.weight.cv[f.idx], alpha=alpha, type.measure = "class", standardize = F)
        
        
        #print("number of predictors included in lasso:")
        #print(length(coefs!=0))
        
        # balanced accuracy:
        #lambda.min = glmnet.tune(x.train, y.train, k, alpha, lambda.seq, f.cv)
        
        #figure.name = paste("cvfit", feature.name[i.feature], toString(k),".png", sep = "_")
        #png(figure.name)
        #plot(cv.fit)
        #dev.off()
        
        #y.pred = predict(cv.fit$glmnet.fit, x.test, s = lambda.min, type = glmnet.para$predict.type)
        #y.pred.train = predict(cv.fit$glmnet.fit, x.train, s = lambda.min, type = glmnet.para$predict.type)
        
        y.pred = predict(glmnet.fit,
                         x.test,
                         s = lambda.min,
                         type = glmnet.para$predict.type)
        
        y.pred.train = predict(glmnet.fit,
                               x.train,
                               s = lambda.min,
                               type = glmnet.para$predict.type)
        
        if (glmnet.para$family == "gaussian") {
            result.test = cor.test(y.test, y.pred, method = "pearson")
            result.train = cor.test(y.train, y.pred.train, method = "pearson")
            
            rmse.test  = sqrt(mean((y.test - y.pred) ^ 2))
            
            print("prediction of testing data:")
            print(result.test)
            print("prediction of training data:")
            print(result.train)
            
            test.result[i, 1] = result.test$estimate
            test.result[i, 2] = rmse.test
            colnames(test.result)[2] = "rmse"
            
            train.result[i, 1] = result.train$estimate
        } else {
            print("prediction of testing data:")
            print(table(y.test, y.pred))
            print("prediction of training data:")
            print(table(y.train, y.pred.train))
            
            result.test = compute.acc(y.test, y.pred)
            result.train = compute.acc(y.train, y.pred.train)
            test.result[i, ] = result.test
            train.result[i, ] = result.train
        }
    } # end outer CV
    
    # compute robustness of coefs across CV:
    coefs.robustness = abs(apply(result.coefs, 1, mean) / apply(result.coefs, 1, sd))
    coefs.sd = apply(result.coefs, 1, sd)
    #print("robustness")
    #print(coefs.robustness)
    result.coefs.df = data.frame(c("intercept", colnames(x)),
                                 result.coefs,
                                 coefs.sd,
                                 coefs.robustness,
                                 stringsAsFactors = F)
    
    # convert coefs to 1 or 0 as indicator of feature selection. 1st row as intercept is removed:
    result.coefs.ind = matrix(
        as.numeric(result.coefs != 0),
        nrow = nrow(result.coefs),
        ncol = ncol(result.coefs)
    )[-1, ]
    
    # compute the number of times each feature is selected in the cross validation:
    coefs.ind.sum = apply(result.coefs.ind, 1, sum)
    #reproducibility.index = (abs(coefs.ind.sum - k/2)-(k/2)%%1)/floor(k/2)
    reproducibility.index = rep(NA, length(coefs.ind.sum))
    none.zero.index = coefs.ind.sum > 0
    #print(result.coefs)
    reproducibility.index[none.zero.index] = (abs(coefs.ind.sum[none.zero.index] - k /
                                                      2) - (k / 2) %% 1) / floor(k / 2)
    # it is so weird that when reproducibility.index is cbind to other matrix and vector and changed to data frame, it become a factor!!!
    # and it won't happen if we use data.frame(...) directly rather than as.data.frame(cbind(...))
    result.coefs.ind = data.frame(
        colnames(x),
        result.coefs.ind,
        coefs.ind.sum,
        reproducibility.index,
        stringsAsFactors = F
    )
    
    #print(result.coefs.df)
    
    glmnet.out = list()
    glmnet.out$test.result = test.result
    glmnet.out$train.result = train.result
    glmnet.out$coefs.ind = result.coefs.ind
    glmnet.out$coefs = result.coefs.df
    glmnet.out$coefs.robustness = mean(coefs.robustness, na.rm = T)
    glmnet.out$reproducibility = mean(reproducibility.index, na.rm = T)
    glmnet.out$penalty.weights = result.feature.weights
    
    if (glmnet.para$return.mod) {
        mod.all = glmnet(
            x,
            y,
            family = glmnet.para$family,
            lambda = lambda.min,
            alpha = alpha[min.i.alpha],
            penalty.factor = f.weight.cv,
            standardize = F
        )
        
        glmnet.out$mod = mod.all
    }
    return(glmnet.out)
    
    print("glmnet.nested.cv: finished")
}

library(e1071)
# svm tune.svm
library(plyr)
# rbind.fill()

svm.weights <- function(model) {
    # This function gives the weights of the hiperplane
    w = 0
    if (model$nclasses == 2) {
        w = t(model$coefs) %*% model$SV
    } else{
        #when we deal with OVO svm classification
        ## compute start-index
        start <- c(1, cumsum(model$nSV) + 1)
        start <- start[-length(start)]
        
        calcw <- function (i, j) {
            ## ranges for class i and j:
            ri <- start[i]:(start[i] + model$nSV[i] - 1)
            rj <- start[j]:(start[j] + model$nSV[j] - 1)
            
            ## coefs for (i,j):
            coef1 <- model$coefs[ri, j - 1]
            coef2 <- model$coefs[rj, i]
            ## return w values:
            w = t(coef1) %*% model$SV[ri, ] + t(coef2) %*% model$SV[rj, ]
            return(w)
        }
        
        W = NULL
        for (i in 1:(model$nclasses - 1)) {
            for (j in (i + 1):model$nclasses) {
                wi = calcw(i, j)
                W = rbind(W, wi)
            }
        }
        w = W
    }
    return(w)
}

select.feature = function(feature.in,
                          factor,
                          p,
                          method = "ttest",
                          k = 5) {
    factor = as.numeric(factor)
    f.idx = rep(NA, ncol(feature.in))
    p = min(p, ncol(feature.in))
    
    for (i.feature in 1:ncol(feature.in)) {
        if (method == "ttest") {
            test.result = t.test(feature.in[, i.feature] ~ factor, var.equal = TRUE)
            f.idx[i.feature] = test.result$p.value
        } else if (method == "wilcox") {
            test.result = wilcox.test(feature.in[, i.feature], factor)
            f.idx[i.feature] = test.result$p.value
        } else if (method == "kendall") {
            test.result = cor.test(feature.in[, i.feature], factor, method = "kendall")
            f.idx[i.feature] = test.result$p.value
        } else if (method == "spearman") {
            test.result = cor.test(feature.in[, i.feature], factor, method = "spearman")
            f.idx[i.feature] = test.result$p.value
        } else if (method == "coef.variation") {
            #cv.mean = aggregate(feature.in, list(cv.k), FUN=mean)
            f.idx = feature.cv(feature.in, factor, k)
            break
        } else if (method == "boot.var") {
            f.idx = feature.cv.boot(feature.in, factor, k)
            break
        } else if (method == "cv.wilcox") {
            f.idx = feature.cv.test(feature.in, factor, method = "wilcox", k)
            break
        }
    }
    #test.result = lapply(feature.in, function(x) t.test(x ~ factor, var.equal = TRUE))
    #p.value = test.result$p.value;
    if (p < 1) {
        # select feature with p less than p threshold:
        feature.idx = which(f.idx < p)
    }
    else{
        # select features with most significant p values
        p.sort = sort(f.idx,
                      index.return = T,
                      decreasing = F)
        feature.idx = p.sort$ix[1:p]
    }
    #print("selected features:")
    #print(feature.idx)
    return(feature.idx)
}

svm.cv.fun = function(brain.feature,
                      subject.info,
                      cost.seq,
                      svm.para) {
    # num.feature: a sequence of number of selected features or threshold of p values: see function select.feature
    # the best value is selected based on the training data:
    subject.info$factor = as.factor(subject.info$factor)
    set.seed(333)
    cv.k = createFolds(subject.info$factor, svm.para$nfolds, list = F)
    
    num.feature2 = seq(
        svm.para$num.feature2.start,
        svm.para$num.feature2.end,
        by = svm.para$num.feature2.step
    )
    
    test.result = data.frame(
        acc = rep(NA, svm.para$nfolds),
        sensi = rep(NA, svm.para$nfolds),
        speci = rep(NA, svm.para$nfolds)
    )
    train.result = data.frame(
        acc = rep(NA, svm.para$nfolds),
        sensi = rep(NA, svm.para$nfolds),
        speci = rep(NA, svm.para$nfolds)
    )
    
    tune.result = matrix(NA, svm.para$nfolds + 1, length(num.feature2))
    tune.result[1,] = num.feature2
    
    tune.control = tune.control(cross = svm.para$nfolds)
    
    for (i in 1:svm.para$nfolds) {
        brain.feature.train = brain.feature[cv.k != i, ]
        brain.feature.test = brain.feature[which(cv.k == i), , drop = F]
        subject.info.train = subject.info[cv.k != i, ]
        subject.info.test = subject.info[which(cv.k == i), , drop = F]
        
        length.num.feature = length(num.feature2)
        feature.idx = list()
        best.cost = vector()
        
        # select features by runing t test or kendall's tau:
        f.idx1 = select.feature(
            brain.feature.train,
            subject.info.train$factor,
            svm.para$num.feature1,
            svm.para$feature.selection1,
            svm.para$nfolds
        )
        
        for (j in 1:length.num.feature) {
            # select feature by coefficient of variation:
            f.idx2 = select.feature(
                brain.feature.train[, f.idx1],
                subject.info.train$factor,
                num.feature2[j],
                svm.para$feature.selection2,
                svm.para$nfolds
            )
            feature.idx[[j]] = f.idx2[!is.na(f.idx2)]
            
            #feature.idx = svmrfe(brain.feature.train, subject.info.train$factor, num.feature)
            
            print("selecting features for svm:")
            print(feature.idx[[j]])
            feature.train.j = cbind(subject.info.train, brain.feature.train[, feature.idx[[j]], drop =
                                                                                F])
            set.seed(123)
            tune.svm = tune(
                svm,
                factor ~ .,
                data = feature.train.j,
                kernel = "linear",
                ranges = list(cost = cost.seq),
                tunecontrol = tune.control
            )
            
            tune.result[i + 1, j] = tune.svm$best.performance
            best.cost[j] = tune.svm$best.parameters$cost
            
            print(summary(tune.svm))
            #print(tune.svm$best.parameters$cost)
            #print(tune.svm$best.modal$coefs)
            #print(svm.weights(tune.svm$best.modal))
        }
        
        best.j = which.min(tune.result[i + 1, ])
        feature.train = cbind(subject.info.train, brain.feature.train[, feature.idx[[best.j]], drop =
                                                                          F])
        feature.test = cbind(subject.info.test, brain.feature.test[, feature.idx[[best.j]], drop =
                                                                       F])
        svm.mod = svm(factor ~ .,
                      data = feature.train,
                      kernel = "linear",
                      cost = best.cost[best.j])
        
        weights.i = data.frame(svm.weights(svm.mod))
        #colnames(weights.i) = colnames(subset(feature.train, select = -factor))
        
        if (i == 1) {
            feature.weights = weights.i
        } else{
            feature.weights = rbind.fill(feature.weights, weights.i)
        }
        
        svm.pred = predict(svm.mod, subset(feature.test, select = -factor))
        svm.pred.train = predict(svm.mod, subset(feature.train, select = -factor))
        
        #print("prediction of testing data:")
        #print(cbind(feature.test$factor, svm.pred))
        #print(table(feature.test$factor, svm.pred))
        #print("prediction of training data:")
        #print(table(feature.train$factor, svm.pred.train))
        
        result.test = compute.acc(feature.test$factor, svm.pred)
        result.train = compute.acc(feature.train$factor, svm.pred.train)
        
        test.result[i, ] = result.test
        train.result[i, ] = result.train
    }
    
    print("tune number of features:")
    print(tune.result)
    
    return(list(test.result, train.result, feature.weights, tune.result))
}
