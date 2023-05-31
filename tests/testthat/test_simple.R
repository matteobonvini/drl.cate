# context("tests with simple simulation setups")
# library(drl.cate)
#
# simple_sim <- function(n) {
#
#   x1 <- rnorm(n)
#   x2 <- rnorm(n)
#   y <- rnorm(n)
#   a <- rbinom(n, 1, 0.5)
#
#   dat <- cbind(y= y, a = a, X1 = x1, X2 = x2)
#
#   return(dat)
#
# }
#
# test_that("lasso based tests with no confounding and no effect", {
#
#   dat <- simple_sim(1000)
#   y <- dat[, "y"]
#   a <- dat[, "a"]
#   x <- dat[, c("X1", "X2"), drop = FALSE]
#   v <- x
#   v0 <- matrix(rep(seq(-0.5, 0.5, length.out = 100), 2), ncol = 2,
#                dimnames = list(NULL, c("X1", "X2")))
#
#   order_basis <- 0
#   lp <- orthopolynom::legendre.polynomials(n = order_basis + 1, normalized = TRUE)
#   lp.fns <- lapply(1:(order_basis + 1), function(u) as.function(lp[[u]]))
#
#   basis.legendre <- function(x, j) {
#     return(lp.fns[[j]](x))
#   }
#
#   h <- 0.2
#   kernel.gaussian <- function(x, x0) {
#     tmp <- function(u) prod(dnorm( as.matrix((u - x0)) / h) / h)
#     out <- apply(x, 1, tmp)
#     return(out)
#   }
#
#   dr.fit <- cate(v0 = v0, v = v, learner = c("dr", "u", "t", "lp-r"), y = y,
#                  a = a, x = x, nsplits = 1, pi.x.method = "lasso",
#                  mu1.x.method = "lasso", mu0.x.method = "lasso",
#                  mu.x.method = "lasso", drl.method = "lasso", ul.method = "lasso",
#                  order_basis = order_basis, kernel = kernel.gaussian,
#                  basis = basis.legendre)
#   # naive estimate is "true" because no confounding
#   oracle_truth <- mean(y[a == 1]) - mean(y[a == 0])
#   ests <- dr.fit$est
#   index.v0 <- 1:nrow(v0)
#   # Because we are using a linear model, estimated CATE needs to be linear.
#   expect_true(max(abs(resid(lm(ests[, "dr"] ~ index.v0)))) < 1e-15)
#   expect_true(max(abs(resid(lm(ests[, "u"] ~ index.v0)))) < 1e-15)
#   expect_true(max(abs(resid(lm(ests[, "t"] ~ index.v0)))) < 1e-15)
#
#   # predictions need to be ballpark correct
#   expect_true(max(abs(ests[, "dr"] - oracle_truth)) < 0.2)
#   expect_true(max(abs(ests[, "u"] - oracle_truth)) < 0.2)
#   expect_true(max(abs(ests[, "t"] - oracle_truth)) < 0.2)
#   expect_true(max(abs(ests[, "lp-r"] - oracle_truth)) < 10) # why lp-r performs bad?! need to find the bug.
# })
