# called by simulation scripts.
# compute feature weights and generate simulation data:

# Xin Niu 2019.06.05


library(caret)

scale.0.1 = function(dat) {
    # the output will be coerced to matrix.
    
    dat = as.matrix(dat)
    
    mins = apply(dat, 2, min)
    maxs = apply(dat, 2, max)
    
    scaled.dat = scale(dat, center = mins, scale = maxs - mins)
    return(scaled.dat)
}


default.cutoff = c(.05, .95)
default.base = 10


weight.cutoff = function(feature.weight, cutoff = default.cutoff) {
    #print(feature.weight)
    
    cutoff.value = quantile(feature.weight, cutoff[1], na.rm = T)
    # when we set this to 0, we have convergence problem
    feature.weight[feature.weight<=cutoff.value] = 1e-2
    
    cutoff.value = quantile(feature.weight, cutoff[2], na.rm = T)
    feature.weight[feature.weight>cutoff.value] = 1
    
    return(feature.weight)
    
}

weight.transform = function(feature.weight, log.base = default.base) {
    #print(feature.weight)
    
    feature.weight = scale(log(feature.weight, base = log.base))
    feature.weight = 1/(1 + exp(-1*feature.weight))
    
    return(feature.weight)
}


compute.acc = function(y, yhat, type = 'class') {
  
  
  if (type == 'response'){
    print('compute.acc: type = response')
    
    cor = cor(y,yhat, method = 'pearson')
    mae = mean(abs(y-yhat))
    #rmse = sqrt(mean((y-yhat)**2))
    
    bias = cor(y, yhat-y, method = 'pearson')
    
    out = c(mae, cor, bias)
    return(out)
    
  }else if (type == 'class') {
    print('compute.acc: type = class')
    # caution: y and yhat should not be swapped.
    # 
    y = as.numeric(as.character(y))
    yhat = as.numeric(as.character(yhat))
    
    acc <- sum(y == yhat, na.rm = TRUE) / length(y)
    
    ylevel = sort(unique(y), decreasing = F)
    
    if (length(ylevel) == 2) {
      # asuming ylevel = c(0, 1)
      sensi <- sum(y == yhat & y == ylevel[2], na.rm = TRUE) / sum(y == ylevel[2], na.rm = TRUE)
      speci <- sum(y == yhat & y == ylevel[1], na.rm = TRUE) / sum(y == ylevel[1], na.rm = TRUE)
    }
    else if (max(yhat) == max(y)) {
      print('compute.acc: levels of y is not 2, may(y)==max(yhat)')
      sensi <- sum(y == yhat & y == ylevel, na.rm = TRUE) / sum(y == ylevel, na.rm = TRUE)
      speci <- NaN
    } else{
      print('compute.acc: levels of y is not 2')
      speci <- sum(y == yhat & y == ylevel, na.rm = TRUE) / sum(y == ylevel, na.rm = TRUE)
      sensi <- NaN
    }
    
    temp <- c(acc, sensi, speci)
    return(temp)
  }
}



reproducibility.index = function(x, k){
  r =  (abs(x - k / 2) - (k / 2) %% 1) / floor(k / 2)
  return(r)
}


feature.weight.test = function(x, y, method = "ttest", cutoff = default.cutoff, log.base = default.base){
    
    feature.weight = matrix(NA, ncol(x),1)
    for (i in 1:ncol(x)) {
        if (method == "ttest") {
            # ttest is not eligible for regression problems.
            test.result = t.test(x[, i]~y)
            #p = plot(y, x[,i], main = toString(test.result$p.value))
            #print(p)
            value = test.result$p.value
        } else if (method == "kendall") {
            test.result = cor.test(x[, i], y, method = "kendall")
            value = abs(test.result$estimate)
        } else if (method == "wilcox") {
            test.result = wilcox.test(x[, i], y)
            value = abs(test.result$estimate)
        } else if (method == "spearman") {
            test.result = cor.test(x[, i], y, method = "spearman")
            value = abs(test.result$estimate)
        } else if (method == "pearson") {
            test.result = cor.test(x[, i], y)
            value = abs(test.result$estimate)
        } else if (method == "biserial") {
            value = abs(biserial.cor(x[, i], y))
        }
        feature.weight[i] = value
    }
    
    feature.weight = weight.transform(feature.weight, log.base)
    feature.weight = weight.cutoff(feature.weight, cutoff)
    
    return(feature.weight)
    
}


# t-test with bootstrap: this method is extremely slow:
feature.weight.test.boot = function(x, y, method = "ttest", nboots = 50, 
                                    cutoff = default.cutoff, log.base = default.base){
    
    f.w.boot = matrix(NA, nboots, ncol(x))
    
    for (i.boot in 1:nboots){
        boot.idx = sample(1:nrow(x), nrow(x), replace = T)
        x.boot = x[boot.idx,]
        y.boot = y[boot.idx]
        
        for (i in 1:ncol(x)) {
            f.w.boot[i.boot,] = feature.weight.test(x.boot, y.boot, method)
        }
    }
    
    feature.weight = apply(f.w.boot, 2, sd)/apply(f.w.boot, 2, mean)
    feature.weight = weight.transform(feature.weight, log.base)
    feature.weight = weight.cutoff(feature.weight, cutoff)
    
    return(feature.weight)
    
}

# mean.diff.boot:
feature.weight.mean.diff.boot = function(x, y, nboots = 500, cutoff = default.cutoff, log.base = default.base){
    
    f.w.boot = matrix(NA, nboots, ncol(x))
    
    for (i in 1:nboots) {
        set.seed(i)
        boot.idx = sample(1:nrow(x), nrow(x), replace = T)
        
        x.boot = x[boot.idx,]
        y.boot = y[boot.idx]
        
        mean1 = apply(x.boot[y.boot==0,], 2, mean)
        mean2 = apply(x.boot[y.boot==1,], 2, mean)
        
        # move abs() to abs(apply(f.w.boot, 2, mean))
        #f.w.boot[i, ] = abs(mean1 - mean2)
        f.w.boot[i, ] = mean1 - mean2
    }
    
    feature.weight = apply(f.w.boot, 2, sd)/abs(apply(f.w.boot, 2, mean))
    feature.weight = weight.transform(feature.weight, log.base)
    feature.weight = weight.cutoff(feature.weight, cutoff)
    
    return(feature.weight)
    
}

fisher_z = function(r){
  
  z = .5*(log(1+r)-log(1-r))
  return(z)
  
}


feature.weight.corr.boot = function(x, y, nboots = 500, cutoff = default.cutoff, log.base = default.base){
  
  f.w.boot = matrix(NA, nboots, ncol(x))
  
  for (i in 1:nboots) {
    set.seed(i)
    boot.idx = sample(1:nrow(x), nrow(x), replace = T)
    
    x.boot = x[boot.idx,]
    y.boot = y[boot.idx]
    f.w.boot[i, ] = fisher_z(cor(x.boot, y.boot, method = 'pearson'))
    
  }
  
  feature.weight = apply(f.w.boot, 2, sd)/abs(apply(f.w.boot, 2, mean))
  feature.weight = weight.transform(feature.weight, log.base)
  feature.weight = weight.cutoff(feature.weight, cutoff)
  
  return(feature.weight)
  
}

# sd and mean diff boot:
feature.weight.sd.mean.diff.boot = function(x, y, nboots, cutoff = default.cutoff, log.base = default.base){
    
    f.w.boot = matrix(NA, nboots, ncol(x))
    
    for (i in 1:nboots) {
        set.seed(i)
        boot.idx = sample(1:nrow(x), nrow(x), replace = T)
        x.boot = x[boot.idx,]
        y.boot = y[boot.idx]
        
        mean1 = apply(x.boot[y.boot==0,], 2, mean)
        mean2 = apply(x.boot[y.boot==1,], 2, mean)
        std1 = apply(x.boot[y.boot==0,], 2, sd)
        std2 = apply(x.boot[y.boot==1,], 2, sd)
        
        f.w.boot[i, ] = abs(mean1 - mean2) + abs(std1 - std2)
    }
    
    feature.weight = apply(f.w.boot, 2, sd)/apply(f.w.boot, 2, mean)
    feature.weight = weight.transform(feature.weight, log.base)
    feature.weight = weight.cutoff(feature.weight, cutoff)
    
    return(feature.weight)
    
}

# mean.diff kfold:
feature.weight.mean.diff.kfold = function(x, y, k = 10, cutoff = default.cutoff, log.base = default.base){
    
    f.w.boot = matrix(NA, k, ncol(x))
    set.seed(111)
    cv.k = createFolds(y, k, list = F)
    
    for (i in 1:k) {
        x.boot = x[cv.k ==i,]
        y.boot = y[cv.k ==i]
        mean1 = apply(x.boot[y.boot==0,], 2, mean)
        mean2 = apply(x.boot[y.boot==1,], 2, mean)
        #std1 = apply(x.boot[y.boot==0,], 2, sd)
        #std2 = apply(x.boot[y.boot==1,], 2, sd)
        #f.w.boot[i, ] = ((mean1 + mean2)*(std1 + std2))/((mean1 - mean2)*(std1 - std2))
        f.w.boot[i, ] = mean1 - mean2
    }
    
    feature.weight = apply(f.w.boot, 2, sd)/abs(apply(f.w.boot, 2, mean))
    feature.weight = weight.transform(feature.weight, log.base)
    feature.weight = weight.cutoff(feature.weight, cutoff)
    
    return(feature.weight)
    
}

#meshgrid funtion copied form:
#https://github.com/cran/pracma/blob/master/R/meshgrid.R

meshgrid <- function(x, y = x) {
  if (!is.numeric(x) || !is.numeric(y))
    stop("Arguments 'x' and 'y' must be numeric vectors.")
  
  x <- c(x); y <- c(y)
  n <- length(x)
  m <- length(y)
  
  X <- matrix(rep(x, each = m),  nrow = m, ncol = n)
  Y <- matrix(rep(y, times = n), nrow = m, ncol = n)
  
  return(list(X = X, Y = Y))
}


library(MASS) # mvrnorm
simulate_data = function(num.samples, beta, x.sd, type = 'class', seed = 111) {
   
    num.predictors = length(beta)
    
    set.seed(seed)
    x.err = runif(num.predictors ,0.5, 5)
    
    set.seed(seed)
    # create x compound gaussian distribution:
    #x = matrix(rnorm(num.predictors* num.samples, 0, x.sd), nrow = num.samples, ncol = num.info.predictors)
    
    # multinormal distribution:
    sigma = diag(x = 1, nrow = num.predictors, ncol = num.predictors)
    sigmamesh = meshgrid(1:num.predictors)
    
    for (i.pred in 1: num.predictors-1){
        sigma[abs(sigmamesh[[1]]-sigmamesh[[2]])==i.pred] = .5^i.pred
    }
    x = mvrnorm(n = num.samples, mu = rep(0, num.predictors), Sigma = x.sd*sigma)
    
    for (ix in 1:num.predictors){
        set.seed(ix+seed)
        x[,ix] = x[,ix] + rnorm(num.samples, 0, x.err[ix])
    }
    
    #print(dim(x))
    set.seed(seed)
    
    if (type == 'class'){
      linpred = x %*% beta + rnorm(num.samples, 0, 5) 
      prob <- 1 / (1 + exp(-linpred)) 
      set.seed(seed+222)
      y <- as.factor(rbinom(num.samples, 1, prob))
    } else {
      linpred = x %*% beta + rnorm(num.samples, 0, 1) 
      y = linpred
    }

    
    #x = scale.0.1(x)
    
    return(list(X = x, y = y))
}
