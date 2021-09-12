#!/usr/bin/env Rscript

library(caret) 
# createFolds()

scale.0.1 = function(dat){
	# the output will be coerced to matrix.
	
	dat = as.matrix(dat)
	
	mins = apply(dat, 2, min)
	maxs = apply(dat, 2, max)
	
	scaled.dat = scale(dat, center = mins, scale = maxs-mins)
	return(scaled.dat)
}

compute.acc = function(y,yhat){
	
	y = as.numeric(y)
	yhat = as.numeric(yhat)

	acc<-sum(y==yhat)/length(y)
	
	ylevel = sort(unique(y), decreasing = F)
	
	if (length(ylevel)==2){
		# asuming ylevel = c(0, 1)
		sensi<-sum(y==yhat & y==ylevel[2])/sum(y==ylevel[2])
		speci<-sum(y==yhat & y==ylevel[1])/sum(y==ylevel[1])
	}
	else if (max(yhat)==max(y)){
		sensi<-sum(y==yhat & y==ylevel)/sum(y==ylevel)
		speci<-NaN
	} else{
		speci<-sum(y==yhat & y==ylevel)/sum(y==ylevel)
		sensi<-NaN
	}

	temp<-c(acc, sensi, speci)
	return(temp)
}

remap.factor = function(f,min=0,max=1){
	f = as.numeric(f)
	
	f.remap = f
	f.remap[f==max(f)]=max
	f.remap[f==min(f)]=min
	return(f.remap)
}

library(VGAM)
# t.test()
# cor.test()
library(ltm)
# biserial.cor()

feature.cv.test = function(feature.in, factor, k, method="wilcox", seed = 111){
	
	num.feature = ncol(feature.in)
	num.sample = nrow(feature.in)
	set.seed(seed)
	
	factor = as.numeric(factor)
	idx.factor.0 = factor==min(factor)
	idx.factor.1 = factor==max(factor)
	
	if (k>1){
		cv.k = createFolds(factor, k, list=F)
		print("number of samples in each CV, for factor 0 and 1:")
		print(table(cv.k[idx.factor.0]))
		print(table(cv.k[idx.factor.1]))
	}else{
		cv.k = rep(0, num.sample)
	}


	
	test.out = matrix(NA, k, num.feature)
	
	for ( i in 1:k) {
		
		idx = cv.k!=i

		if (method=="mean.diff"){
			test.out[i,] = abs(apply(feature.in[cv.k!=i&idx.factor.0, ], 2, mean) - apply(feature.in[cv.k!=i&idx.factor.1, ], 2, mean))
		}else{

			for (i.feature in 1:num.feature){
				if (method=="kendall"){
					test.result = cor.test(feature.in[idx, i.feature], factor[dix], method = "kendall")
				}else if (method == "wilcox"){
					test.result = wilcox.test(feature.in[idx, i.feature], factor[idx])
				}else if (method == "spearman"){
					test.result = cor.test(feature.in[idx,i.feature], factor[idx], method = "spearman")
				}else if (method == "pearson"){
					test.result = cor.test(feature.in[idx,i.feature], factor[idx])
				}else if (method == "biserial"){
					test.out[i, i.feature] = abs(biserial.cor(feature.in[idx,i.feature], factor[idx]))
					next
				}
				#test.out[i, i.feature] = test.result$p.value
				test.out[i, i.feature] = abs(test.result$estimate)
			}
		}
	}
	return(test.out)
}

#feature.cv.mean.diff = function(feature.in, factor, k=10, seed = 111){
#	#feature.in = feature.in[index,]
#	
#	num.sample = nrow(feature.in)
#	set.seed(seed)
#	cv.k = createFolds(factor, k, list=F)
#
#	factor = as.numeric(factor)
#
#	idx.factor.0 = factor==min(factor)
#	idx.factor.1 = factor==max(factor)
#
#	#print("number of samples in each CV, for factor 0 and 1:")
#	#print(table(cv.k[idx.factor.0]))
#	#print(table(cv.k[idx.factor.1]))
#	
#	# create NA matrix to save mean value for each folds, in case some folds may missing thus cause the dimension of cv.mean.0 and cv.mean.1 are not match.
#	cv.mean.0 = matrix(NA, k, ncol(feature.in))
#	cv.mean.1 = matrix(NA, k, ncol(feature.in))
#	
#	#agg.mean.0 = aggregate(feature.in[idx.factor.0,], list(cv.k[idx.factor.0]), FUN=mean)
#	#agg.mean.1 = aggregate(feature.in[idx.factor.1,], list(cv.k[idx.factor.1]), FUN=mean)
#	
#	#cv.idx.0 = sort(unique(cv.k[idx.factor.0]))
#	#cv.idx.1 = sort(unique(cv.k[idx.factor.1]))
#
#	#cv.mean.0[cv.idx.0,] = as.matrix(agg.mean.0[,-1])
#	#cv.mean.1[cv.idx.1,] = as.matrix(agg.mean.1[,-1])
#	
#	for ( i in 1:k){
#		cv.mean.0[i,] = apply(feature.in[cv.k!=i&idx.factor.0,], 2, mean)
#		cv.mean.1[i,] = apply(feature.in[cv.k!=i&idx.factor.1,], 2, mean)
#	}
#
#	cv.mean.diff = abs(cv.mean.0 - cv.mean.1)
#	
#	return(cv.mean.diff)
#}

glmnet.coef.cv = function(x, y, nfold, seed = 111) {
  num.sample = nrow(x)
  set.seed(seed)
  cv.k = createFolds(y, nfold, list=F)
  coefs = matrix(NA, nfold, ncol(x))
  
  for (i in 1:nfold) {
    idx.train = cv.k!=i
    x.train = x[idx.train,]
    y.train = y[idx.train]
    cv.fit = cv.glmnet(x.train, y.train, nfolds = nfold-1, family = "binomial", standardize = F)
    coefs.i = coef(cv.fit, s="lambda.min")
	coefs[i,] = coefs.i[-1]
  }
  
  coef.cv = apply(coefs, 2, sd, na.rm=T)/apply(coefs, 2, mean, na.rm=T)
  # replace NaN with inf:
  coef.cv[is.nan(coef.cv)]=10
  return(coef.cv)
}

    
feature.cv.boot = function(feature.in, factor, n=100, method){
	f.cv.boot = matrix(NA, n, ncol(feature.in))
	factor = remap.factor(factor)
		
	for (i in 1:n) {
		set.seed(i)
    	boot.idx = sample(1:nrow(feature.in), nrow(feature.in), replace = T)
    	feature.boot = feature.in[boot.idx, ]
    	factor.boot = factor[boot.idx]
		f.cv.boot[i,] = feature.cv.test(feature.boot, factor.boot, k=1, method)
	}
  	
	return(f.cv.boot)
}


library(glmnet)
#library(SIS)

glmnet.tune = function(x, y, k, alpha = 1, lambda.seq = 10^seq(-4, 3, length=70), f.weights){
	#set.seed(222)
	set.seed(123)
	
	if (missing(f.weights)){
		f.weights = rep(1, dim(x)[2])
	}

	cv.k = createFolds(y, k, list=F)
	test.acc = matrix(NA, k, length(lambda.seq))

	glmnet.control(mxit = 1000000)
	
	k = max(c(5,k))
	
	for (i in 1:k){	
		x.train = x[cv.k!=i,]
		x.test = x[cv.k==i,]
		y.train = y[cv.k!=i]
		y.test = y[cv.k==i]
		
		for (j in 1:length(lambda.seq)){
	
			mod = glmnet(x.train, y.train, family = "binomial", alpha=alpha, lambda = lambda.seq[j], standardize = F, penalty.factor = f.weights)
			y.pred = predict(mod, x.test, s = lambda.seq[j], type="class")
			acc = compute.acc(y.pred, y.test)

			test.acc[i,j] = (acc[2]+acc[3])/2
			#test.acc[i,j] = acc[1]
		}
	}

	acc.cv = apply(test.acc, 2, sd)/apply(test.acc, 2, mean)
	min.idx = which.min(acc.cv)
	return(lambda.seq[min.idx])
}


glmnet.cv.fun = function(x, y, lambda.seq, glmnet.para){
	set.seed(222)
	
	k = glmnet.para$nfolds
	penalty.weight = glmnet.para$penalty.weight
	quantile.thresh = glmnet.para$quantile.thresh
	alpha = glmnet.para$alpha

	cv.k = createFolds(y, k, list=F)
	y.num = as.numeric(y)
	idx.factor.0 = y.num==min(y.num)
	idx.factor.1 = y.num==max(y.num)
	print("number of samples in each CV, for factor 0 and 1:")
	print(table(cv.k[idx.factor.0]))
	print(table(cv.k[idx.factor.1]))

#	print(t(rbind(cv.k, y)))

	test.result=data.frame(acc=rep(NA,k), sensi=rep(NA,k), speci=rep(NA,k))
	train.result=data.frame(acc=rep(NA,k), sensi=rep(NA,k), speci=rep(NA,k))
	result.coefs = matrix(NA, ncol(x)+1, k)
	result.feature.weights= matrix(NA, ncol(x), k)

	for (i in 1:k){	
		x.train = x[cv.k!=i,]
		x.test = x[cv.k==i,]
		y.train = y[cv.k!=i]
		y.test = y[cv.k==i]
	
		#cv.fit = cv.ncvreg(x.train, y.train, seed=123, family="binomial")
		#cv.fit<-SIS(x.train, y.train, family = "binomial", penalty = p, tune = "cv", nfolds = 10, type.measure = "deviance", seed = 123, iter.max = 5000)
		print("computing feature weights:")
		print(penalty.weight)
		
		set.seed(123)
		glmnet.control(mxit = 1000000)		
		
		if (penalty.weight == "none"){
			f.weight.cv = rep(1, dim(x)[2])
		}else { 
			if (penalty.weight == "mean.diff"){
				f.weight = feature.cv.test(x.train, y.train, 10, "mean.diff")
			}else if (penalty.weight == "wilcox"){
				f.weight = feature.cv.test(x.train, y.train, 10, "wilcox")
			}else if (penalty.weight == "pearson"){
				f.weight = feature.cv.test(x.train, y.train, 10, "pearson")
			}else if (penalty.weight == "mean.diff.boot"){
				f.weight = feature.cv.boot(x.train, y.train, n=500, method = "mean.diff")
			}else if (penalty.weight == "pearson.boot"){
				f.weight = feature.cv.boot(x.train, y.train, n=500, method = "pearson")
			}	
			f.weight.cv = scale.0.1(apply(f.weight, 2, sd)/apply(f.weight, 2, mean))
		}
	

		f.idx = f.weight.cv<=quantile(f.weight.cv, quantile.thresh)
	
		#print(f.weight.cv)
		#print(dim(x.train))
		#print(y.train)
		cv.fit = cv.glmnet(x.train, y.train, nfolds = k, family = "binomial", alpha=alpha, penalty.factor = f.weight.cv, type.measure = "class", standardize = F)
		
		#cv.fit = cv.glmnet(x.train[, f.idx], y.train, nfolds = k, family = "binomial", penalty.factor = f.weight.cv[f.idx], alpha=alpha, type.measure = "class", standardize = F)
		lambda.min = cv.fit$lambda.min
		
		coefs <- coef(cv.fit, s="lambda.min")
		result.coefs[,i] <- as.vector(coefs)
		result.feature.weights[,i] <- f.weight.cv

		#print("number of predictors included in lasso:")
		#print(length(coefs!=0))
		
		#lambda.min = glmnet.tune(x.train, y.train, k, alpha, lambda.seq, f.cv)

		#figure.name = paste("cvfit", feature.name[i.feature], toString(k),".png", sep = "_")
		#png(figure.name)
		#plot(cv.fit)
		#dev.off()

		#lambda.min = glmnet.tune(x.train, y.train, k, alpha, lambda.seq, f.cv)
		print(lambda.min)

		mod = glmnet(x.train, y.train, family = "binomial", alpha=alpha, penalty.factor = f.weight.cv, lambda = lambda.min, standardize = F)
		y.pred = predict(mod, x.test, s = lambda.min, type="class")
		y.pred.train = predict(mod, x.train, s = lambda.min, type="class")
		
		print("prediction of testing data:")
		print(table(y.test, y.pred))
		print("prediction of training data:")
		print(table(y.train, y.pred.train))
		
		result.test = compute.acc(y.test, y.pred)
		result.train = compute.acc(y.train, y.pred.train)

		test.result[i,] = result.test
		train.result[i,] = result.train
	}
	result.coefs = matrix(as.numeric(result.coefs!=0), nrow = nrow(result.coefs), ncol = ncol(result.coefs))

	return(list(test.result, train.result, result.coefs, result.feature.weights))
}

library(e1071)
# svm tune.svm
library(plyr)
# rbind.fill()

svm.weights<-function(model){
	# This function gives the weights of the hiperplane
	w=0
	if(model$nclasses==2){
		w=t(model$coefs)%*%model$SV
	}else{ #when we deal with OVO svm classification
  		## compute start-index
  		start <- c(1, cumsum(model$nSV)+1)
  		start <- start[-length(start)]

  		calcw <- function (i,j) {
    	## ranges for class i and j:
    	ri <- start[i] : (start[i] + model$nSV[i] - 1)
    	rj <- start[j] : (start[j] + model$nSV[j] - 1)

  		## coefs for (i,j):
		coef1 <- model$coefs[ri, j-1]
		coef2 <- model$coefs[rj, i]
		## return w values:
		w=t(coef1)%*%model$SV[ri,]+t(coef2)%*%model$SV[rj,]
		return(w)
  }

	W=NULL
	for (i in 1 : (model$nclasses - 1)){
		for (j in (i + 1) : model$nclasses){
			wi=calcw(i,j)
			W=rbind(W,wi)
		}
	}
	w=W
	}
	return(w)
}

select.feature = function(feature.in, factor, p, method = "ttest", k=5){
	
	factor = as.numeric(factor)
	f.idx = rep(NA, ncol(feature.in))
	p = min(p, ncol(feature.in))

	for (i.feature in 1:ncol(feature.in)) {
		if (method == "ttest"){
			test.result = t.test(feature.in[,i.feature] ~ factor, var.equal=TRUE)
			f.idx[i.feature] = test.result$p.value
		}else if (method == "wilcox") {
			test.result = wilcox.test(feature.in[,i.feature], factor)
			f.idx[i.feature] = test.result$p.value
		}else if (method == "kendall") {
			test.result = cor.test(feature.in[,i.feature], factor, method = "kendall")
			f.idx[i.feature] = test.result$p.value
		}else if (method == "spearman") {
			test.result = cor.test(feature.in[, i.feature], factor, method = "spearman")
			f.idx[i.feature] = test.result$p.value
		}else if (method == "coef.variation") {
			#cv.mean = aggregate(feature.in, list(cv.k), FUN=mean)
			f.idx = feature.cv(feature.in, factor, k)
			break
		}else if (method =="boot.var") {
			f.idx = feature.cv.boot(feature.in, factor, k)
			break
		}else if (method =="cv.wilcox") {
			f.idx = feature.cv.test(feature.in, factor, method = "wilcox", k)
			break
		}
	}
	#test.result = lapply(feature.in, function(x) t.test(x ~ factor, var.equal = TRUE))
	#p.value = test.result$p.value;
	if (p<1){
		# select feature with p less than p threshold:
		feature.idx = which(f.idx<p)
		}
	else{
		# select features with most significant p values
		p.sort = sort(f.idx, index.return=T, decreasing=F)
		feature.idx = p.sort$ix[1:p]
	}
	#print("selected features:")
	#print(feature.idx)
	return(feature.idx)
}

svm.cv.fun = function(brain.feature, subject.info, cost.seq, svm.para){
	# num.feature: a sequence of number of selected features or threshold of p values: see function select.feature
	# the best value is selected based on the training data:
	subject.info$factor = as.factor(subject.info$factor)
	set.seed(333)
	cv.k = createFolds(subject.info$factor, svm.para$nfolds, list=F)
		
	num.feature2 = seq(svm.para$num.feature2.start, svm.para$num.feature2.end, by = svm.para$num.feature2.step)
	
	test.result=data.frame(acc=rep(NA,svm.para$nfolds), sensi=rep(NA,svm.para$nfolds), speci=rep(NA,svm.para$nfolds))
	train.result=data.frame(acc=rep(NA,svm.para$nfolds), sensi=rep(NA,svm.para$nfolds), speci=rep(NA,svm.para$nfolds))
	
	tune.result = matrix(NA, svm.para$nfolds+1, length(num.feature2))
	tune.result[1, ] = num.feature2

	tune.control = tune.control(cross = svm.para$nfolds)

	for (i in 1:svm.para$nfolds){	
		
		brain.feature.train = brain.feature[cv.k!=i,]
		brain.feature.test = brain.feature[which(cv.k==i), , drop = F]
		subject.info.train = subject.info[cv.k!=i,]
		subject.info.test = subject.info[which(cv.k==i), , drop = F]

		length.num.feature = length(num.feature2)
		feature.idx = list()
		best.cost = vector()

		# select features by runing t test or kendall's tau:
		f.idx1 = select.feature(brain.feature.train, subject.info.train$factor, svm.para$num.feature1, svm.para$feature.selection1, svm.para$nfolds)

		for (j in 1:length.num.feature) {
			
			# select feature by coefficient of variation:
			f.idx2 = select.feature(brain.feature.train[,f.idx1], subject.info.train$factor, num.feature2[j], svm.para$feature.selection2, svm.para$nfolds)
			feature.idx[[j]] = f.idx2[!is.na(f.idx2)]

			#feature.idx = svmrfe(brain.feature.train, subject.info.train$factor, num.feature)
			
			print("selecting features for svm:")
			print(feature.idx[[j]])	
			feature.train.j = cbind(subject.info.train, brain.feature.train[, feature.idx[[j]], drop=F])
			set.seed(123)
			tune.svm = tune(svm, factor ~., data = feature.train.j, kernel = "linear", ranges = list(cost=cost.seq), tunecontrol = tune.control)
			
			tune.result[i+1, j] = tune.svm$best.performance
			best.cost[j] = tune.svm$best.parameters$cost
			
			print(summary(tune.svm))
			#print(tune.svm$best.parameters$cost)
			#print(tune.svm$best.modal$coefs)		
			#print(svm.weights(tune.svm$best.modal))	
		}
		
		best.j = which.min(tune.result[i+1,])
		feature.train = cbind(subject.info.train, brain.feature.train[, feature.idx[[best.j]], drop=F])
		feature.test = cbind(subject.info.test, brain.feature.test[, feature.idx[[best.j]], drop=F])
		svm.mod = svm(factor~., data = feature.train, kernel = "linear", cost = best.cost[best.j])
		
		weights.i = data.frame(svm.weights(svm.mod))
		#colnames(weights.i) = colnames(subset(feature.train, select = -factor))
		
		if (i==1){
			feature.weights = weights.i
		}else{
			feature.weights = rbind.fill(feature.weights, weights.i)
		}

		svm.pred = predict(svm.mod, subset(feature.test, select=-factor))
		svm.pred.train = predict(svm.mod, subset(feature.train, select=-factor))
		
		#print("prediction of testing data:")
		#print(cbind(feature.test$factor, svm.pred))
		#print(table(feature.test$factor, svm.pred))
		#print("prediction of training data:")
		#print(table(feature.train$factor, svm.pred.train))

		result.test = compute.acc(feature.test$factor, svm.pred)
		result.train = compute.acc(feature.train$factor, svm.pred.train)

		test.result[i,] = result.test
		train.result[i,] = result.train
	}
	
	print("tune number of features:")
	print(tune.result)
	
	return(list(test.result, train.result, feature.weights, tune.result))
}


