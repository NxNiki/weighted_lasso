

######################
##generate simulation data
############################
#r=1
#niu=0.5
beta0=0.5
n=sn=1023
#ftao=0.33
#fc=5

##sample size need to be larege enough
m=600
tempx<-read.table(file="tempx1000.txt")
tempx<-as.matrix(tempx)


dis=dist(tempx,method = "euclidean",diag=TRUE,upper=TRUE)
#subdis=dist(subtempx,method = "euclidean",diag=TRUE,upper=TRUE)
dismatrix=as.matrix(dis)
#subdismatrix=as.matrix(subdis)

#112,223,334,445,556,667,778,889,990
set.seed(889)
err<-rnorm(n=1023*600,mean=0,sd=1)
err<-matrix(err,1023,600)
corr.err<-matrix(0,1023,600)

for (i in 1:1023){
indx.temp=which(dismatrix[i,]<=1)
corr.err[i,]=colSums(err[indx.temp,])/sqrt(length(indx.temp))
}
corr.x<-corr.err


xtest=t(corr.x)
xtest=scale(xtest,center=TRUE,scale=TRUE)

betas<-read.table(file="trbetas1000.txt")
betas<-as.matrix(betas)

linpred <- beta0 + (xtest%*%betas)
#prob1 <- exp(linpred)/(1 + exp(linpred))
prob<-pnorm(linpred,0,1)
set.seed(38)
runis <- runif(m,0,1)
ytest <- ifelse(runis < prob,1,0)
sum(ytest==1)
sum(ytest==0)
write.table(ytest,file="ytest1000a600new9.txt",sep="\t",row.names=FALSE,col.names=FALSE)
write.table(xtest,file="xscaled1000a600new9.txt",sep="\t",row.names=FALSE,col.names=FALSE)



