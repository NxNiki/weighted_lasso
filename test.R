x = matrix(rnorm(500*50, 0, 1), nrow = 500, ncol = 50)
y = rnorm(num.samples, 0, 1)
y[y>0] = 1
y[y<=0] = 0
p = 
for (i in 1:50){
    y = sample(y)
    test.result = t.test(x[, i]~y)
    p = plot(y, x[,i], main = toString(test.result$p.value))
    print(p)
    value = test.result$p.value
}