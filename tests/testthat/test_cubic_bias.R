test_that("cubic bias fit and influence function describe the same estimator", {
  set.seed(607)
  x <- runif(80, -1, 2)
  y <- x^3 + rnorm(80, sd=0.1)
  n <- length(x)
  h <- 1.4
  b <- 1.7
  # Independent weighted-empirical-law implementation; bandwidths are fixed.
  estimate <- function(p, a, debias=TRUE) {
    u <- (x-a)/h
    z <- (x-a)/b
    X <- cbind(1, u)
    B <- cbind(1, z, z^2, z^3)
    wh <- p * dnorm(u)/h
    wb <- p * dnorm(z)/b
    Dh <- crossprod(X, X*wh)
    beta <- solve(Dh, crossprod(X, y*wh))
    c2 <- solve(Dh, crossprod(X, u^2*wh))[1]
    bias <- solve(crossprod(B, B*wb), crossprod(B, y*wb))[3]
    unname(beta[1] - as.numeric(debias)*(h/b)^2*c2*bias)
  }
  p <- rep(1/n, n)
  for (a in c(0.2, 1.8)) {
    for (debias in c(FALSE, TRUE)) {
      target <- estimate(p, a, debias)
      actual <- .lprobust(x, y, h, b, debias, eval.pt=a,
                         kernel.type="gau")[1, "theta.hat"]
      expect_equal(as.numeric(actual), as.numeric(target), tolerance=1e-10)
      eps <- 1e-6
      numeric.if <- vapply(seq_len(n), function(i) {
        delta <- as.numeric(seq_len(n)==i)-p
        (estimate(p+eps*delta, a, debias) -
           estimate(p-eps*delta, a, debias))/(2*eps)
      }, numeric(1))
      actual.if <- .compute.rinfl.func(y, x, a, h, b, dnorm, debias)$est
      expect_equal(actual.if, numeric.if, tolerance=1e-6)
    }
  }
  # Test hat diagonals away from the separately audited bandwidth-floor bug.
  h <- b <- 5
  diagonal <- vapply(seq_len(n), function(i) {
    .lprobust(x, as.numeric(seq_len(n)==i), h, b, TRUE,
              eval.pt=x[i], kernel.type="gau")[1, "theta.hat"]
  }, numeric(1))
  actual.hat <- vapply(x, function(a)
    .local.fit(x, y, a, h, b, kern=dnorm, debias=TRUE)$hat, numeric(1))
  expect_equal(actual.hat, diagonal, tolerance=1e-10)
})
