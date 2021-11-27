# Examples for the DR-Learner and Lp-R-Learner.
set.seed(1234)
require(orthopolynom)
require(SuperLearner)
# some helper functions
expit <- function(x) exp(x) / (1 + exp(x))
logit <- function(x) log( x / (1 - x))

# functions to estimate nuisance regressions
mu11.x <- function(y.tr, a.tr, x.tr, new.x, sl.lib = c("SL.mean", "SL.lm", "SL.gam",
                                                       "SL.polymars", "SL.rpart")){
  y1.tr <- y.tr[a.tr == 1]
  x1.tr <- as.data.frame(x.tr[a.tr == 1, , drop = FALSE])
  fit <- SuperLearner::SuperLearner(Y = y1.tr, X = x1.tr, SL.library = sl.lib,
                                    newX = as.data.frame(new.x))
  return(fit$SL.predict)
}

mu0.x <- function(y.tr, a.tr, x.tr, new.x, sl.lib = c("SL.mean", "SL.lm", "SL.gam",
                                                      "SL.polymars", "SL.rpart")){
  y0.tr <- y.tr[a.tr == 0]
  x0.tr <- as.data.frame(x.tr[a.tr == 0, , drop = FALSE])
  fit <- SuperLearner::SuperLearner(Y = y0.tr, X = x0.tr, SL.library = sl.lib,
                                    newX = as.data.frame(new.x))
  return(fit$SL.predict)
}

mu.x <- function(y.tr, a.tr, x.tr, new.x, sl.lib = c("SL.mean", "SL.lm", "SL.gam",
                                                      "SL.polymars", "SL.rpart")){
  fit <- SuperLearner::SuperLearner(Y = y.tr, X = as.data.frame(x.tr), SL.library = sl.lib,
                                    newX = as.data.frame(new.x))
  return(fit$SL.predict)
}

pi.x <- function(a.tr, x.tr, new.x, sl.lib = c("SL.mean", "SL.glm", "SL.ranger",
                                         "SL.rpart")){
  fit <- SuperLearner::SuperLearner(Y = a.tr, X = as.data.frame(x.tr), SL.library = sl.lib,
                                    newX = as.data.frame(new.x))
  return(fit$SL.predict)
}

drl.x <- function(y.tr, x.tr, new.x, sl.lib = c("SL.mean", "SL.lm", "SL.gam",
                                          "SL.loess", "SL.glm", "SL.polymars",
                                          "SL.rpart")){
  fit <- SuperLearner::SuperLearner(Y = y.tr, X = as.data.frame(x.tr),
                                    SL.library = sl.lib,
                                    newX = as.data.frame(new.x))
  return(fit$SL.predict)
}

order_basis <- 10
lp <- legendre.polynomials(n = order_basis, normalized = TRUE)
lp.fns <- lapply(1:(order_basis + 1), function(u) as.function(lp[[u]]))

basis.legendre <- function(x, j) {
  return(lp.fns[[j]](x))
}

h <- 0.1
kernel.gaussian <- function(x, x0) {
  tmp <- function(u) prod(dnorm( as.matrix((u - x0)) / h) / h)
  out <- apply(x, 1, tmp)
  return(out)
}

# some generic parameters
n <- 2000
#################################
# Example 1 from Kennedy (2020) #
#################################
x0 <- seq(-0.8, 0.8, length.out = 101)
mu0_true <- function(y.tr, a.tr, x.tr, new.x) {
  u <- new.x[, 1]
  (u <= -0.5) * 0.5 * (u + 2)^2 +
    (u * 0.5 + 0.875) * (u > -0.5 & u < 0) +
    (u > 0 & u < 0.5) * (-5 * (u - 0.2)^2 + 1.075) +
    (u > .5) * (u + 0.125)
}
mu1_true <- function(y.tr, a.tr, x.tr, new.x) {
  mu0_true(y.tr, a.tr, x.tr, new.x)
}
pix_true <- function(a.tr, x.tr, new.x){
  u <- new.x[, 1]
  0.1 + 0.8 * (u > 0)
}

gen_data <- function(n){
  x <- data.frame(x = runif(n, -1, 1))
  a <- rbinom(n, 1, pix_true(a = NULL, x = NULL, new.x = x))
  y <- a * mu1_true(y = NULL, x = NULL, new.x = x) +
    (1 - a) * mu0_true(y = NULL, x = NULL, new.x = x) +
    rnorm(n, sd = (0.2 - 0.1 * cos(2 * pi * x$x)))

  dat <- data.frame(y = y, a = a, x1 = x$x)
  return(dat)
}

dat <- gen_data(n)
x0 <- matrix(x0, ncol = 1, dimnames = list(NULL, "x1"))
res_drl <- cate(v0 = x0, v = dat[, -c(1:2), drop = FALSE], learner = "dr", y = dat$y,
                a = dat$a, x = dat[, -c(1:2), drop = FALSE], drl = drl.x,
                mu1.x = mu1.x, mu0.x = mu0.x, pi.x = pi.x, nsplits = 5)

res_lprl <- lp_r_learner(x0 = x0, y = dat$y, a = dat$a, x = dat[, -c(1:2), drop = FALSE],
                         mu.x = mu.x, pi.x = pi.x, basis = basis.legendre,
                         order_basis = order_basis, kernel = kernel.gaussian)

oracle <- cate(v0 = x0, v = dat[, -c(1:2), drop = FALSE], learner = "dr", y = dat$y,
               a = dat$a, x = dat[, -c(1:2), drop = FALSE], drl = drl.x,
               mu1.x = mu1_true, mu0.x = mu0_true, pi.x = pix_true, nsplits = 1)

plot(x0, y = rep(0, length(x0)), type = "l", col = "red", ylim = c(-1.5, 1.5),
     ylab = "E(Y^1 - Y^0)", xlab = "X0")
lines(x0, y = oracle$est, col = "blue")
lines(x0, y = res_drl$est, col = "black")
lines(x0, y = res_lprl$est, col = "pink")
legend(0.15, 1.5, legend = c("Truth", "Oracle", "DR-Learner", "Lp-R-Learner"),
       col = c("red", "blue", "black", "pink"), cex = 0.8, lty = 1)



