# Examples for the DR-Learner and Lp-R-Learner.
set.seed(1234)
# some helper functions
expit <- function(x) exp(x) / (1 + exp(x))
logit <- function(x) log( x / (1 - x))

# functions to estimate nuisance regressions
mu1.x <- function(y, x, new.x, sl.lib = c("SL.mean", "SL.lm", "SL.gam",
                                          "SL.polymars", "SL.rpart")){
  fit <- SuperLearner::SuperLearner(Y = y, X = x, SL.library = sl.lib,
                                    newX = new.x)
  return(fit$SL.predict)
}

mu0.x <- mu1.x

pi.x <- function(a, x, new.x, sl.lib = c("SL.mean", "SL.glm", "SL.ranger",
                                         "SL.rpart")){
  fit <- SuperLearner::SuperLearner(Y = a, X = x, SL.library = sl.lib,
                                    newX = new.x, family = binomial())
  return(fit$SL.predict)
}

drl.x <- function(y, x, new.x, sl.lib = c("SL.mean", "SL.lm", "SL.gam",
                                          "SL.loess", "SL.glm", "SL.polymars",
                                          "SL.rpart")){
  fit <- SuperLearner::SuperLearner(Y = y, X = x, SL.library = sl.lib,
                                    newX = new.x)
  return(fit$SL.predict)
}

# some generic parameters
n <- 2000
#################################
# Example 1 from Kennedy (2020) #
#################################
x0 <- seq(-0.8, 0.8, length.out = 101)
mu0_true <- function(y, x, new.x) {
  u <- new.x[, 1]
  (u <= -0.5) * 0.5 * (u + 2)^2 +
    (u * 0.5 + 0.875) * (u > -0.5 & u < 0) +
    (u > 0 & u < 0.5) * (-5 * (u - 0.2)^2 + 1.075) +
    (u > .5) * (u + 0.125)
}
mu1_true <- function(y, x, new.x) {
  mu0_true(y, x, new.x)
}
pix_true <- function(a, x, new.x){
  u <- new.x[, 1]
  0.1 + 0.8 * (u > 0)
}

gen_data <- function(n){
  x <- data.frame(x = runif(n, -1, 1))
  a <- rbinom(n, 1, pix_true(a = NULL, x = NULL, new.x = x))
  y <- a * mu1_true(y = NULL, x = NULL, new.x = x) +
    (1 - a) * mu0_true(y = NULL, x = NULL, new.x = x) +
    rnorm(n, sd = (0.2 - 0.1 * cos(2 * pi * x$x)))

  dat <- data.frame(y = y, a = a, x = x$x)
  return(dat)
}

dat <- gen_data(n)
res <- dr_learner(x0 = x0, y = dat$y, a = dat$a, x = as.data.frame(dat$x),
                  drl.x = drl.x, mu1.x = mu1.x, mu0.x = mu0.x, pi.x = pi.x,
                  nsplits = 5)
oracle <- dr_learner(x0 = x0, y = dat$y, a = dat$a, x = as.data.frame(dat$x),
                     drl.x = drl.x, mu1.x = mu1_true, mu0.x = mu0_true,
                     pi.x = pix_true, nsplits = 5)

plot(x0, y = rep(0, length(x0)), type = "l", col = "red", ylim = c(-1.5, 1.5),
     ylab = "E(Y^1 - Y^0)", xlab = "X0")
lines(x0, y = oracle$est, col = "blue")
lines(x0, y = res$est, col = "black")
legend(0.25, 1.5, legend = c("Truth", "Oracle", "DR-Learner"),
       col = c("red", "blue", "black"), cex = 0.8, lty = 1)
