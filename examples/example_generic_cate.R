library(drl.cate)

n <- 2000
d <- 20
x <- matrix(rnorm(n * d), ncol = d, nrow = n)
a <- rbinom(n, size = 1, prob = 0.5)
v <- matrix(rnorm(n), ncol = 2, nrow = n)
v0 <- matrix(rep(seq(-0.5, 0.5, 0.1), 2), ncol = 2,
             dimnames = list(NULL, colnames(v)))
y <- rnorm(n)
cate(v0 = v0, learner = "dr", y = y, a = a, x = x, v = v,
     nsplits = 2, pi.x.method = "lasso", mu1.x.method = "lasso",
     mu0.x.method = "lasso", drl.method = "lasso")

plot()