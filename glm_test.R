n = 50
p = 100


x = matrix(rnorm(n*p, 0, 1), n, p)
beta = rnorm(p)

y = x %*% beta + rnorm(n)

library(glmnet)

glmnet.fit = glmnet(x, y, lambda = 0)

coef(glmnet.fit)

glmfit = glm.fit(x, y, family = gaussian())
glmfit$coef
