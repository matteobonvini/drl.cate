test_that("PD integration uses aligned fold rows and disjoint columns", {
  n <- 45
  v <- seq(-1, 2, length.out=n)
  w <- data.frame(w=sin(seq_len(n)))
  s <- rep(1:3, length.out=n)
  fits <- lapply(1:3, function(k) {
    force(k)
    list(fit=function(z) z$v1j*z$w + k*z$v1j^2)
  })
  blocks <- .get.muhat(s, fits, v, w, max.n.integral=4, stream=FALSE)
  full <- .get.muhat(s, fits, v, w, max.n.integral=n, stream=FALSE)
  expect_equal(Reduce(`+`, lapply(blocks, function(z) as.integer(z$sub.idx))),
               rep(1L, n))
  for (block in blocks) {
    rows <- block$row.idx
    cols <- which(block$sub.idx)
    expect_lte(length(cols), 4)
    expect_equal(rows, which(s==s[cols[1]]))
    expected <- outer(v[rows], w$w[cols], `*`) + s[cols[1]]*v[rows]^2
    expect_equal(block$muhat.mat, expected)
  }
  m <- numeric(n)
  for (k in 1:3) {
    idx <- which(s==k)
    m[idx] <- v[idx]*mean(w$w[idx]) + k*v[idx]^2
  }
  y <- cos(v) + w$w
  h <- 2
  b <- 2.5
  a <- 0.3
  for (debias in c(FALSE, TRUE)) {
    calc <- function(grid) .compute.rinfl.func(y, v, a, h, b, dnorm,
                                              debias, grid, m)$est
    expect_equal(calc(blocks), calc(full), tolerance=1e-12)
    # Direct per-observation version of the manuscript's integral.
    u <- (v-a)/h
    z <- (v-a)/b
    X <- cbind(1, u)
    B <- cbind(1, z, z^2, z^3)
    kh <- dnorm(u)/h
    kb <- dnorm(z)/b
    Dh <- crossprod(X, X*kh)/n
    Db <- crossprod(B, B*kb)/n
    c2 <- solve(Dh, crossprod(X, u^2*kh)/n)[1]
    gamma <- drop(X %*% solve(Dh)[,1])*kh
    if(debias) gamma <- gamma - (h/b)^2*c2*drop(B %*% solve(Db)[,3])*kb
    direct <- vapply(seq_len(n), function(j) {
      idx <- which(s==s[j])
      mean(gamma[idx]*(v[idx]*w$w[j]+s[j]*v[idx]^2-m[idx]))
    }, numeric(1))
    baseline <- .compute.rinfl.func(y, v, a, h, b, dnorm, debias)$est
    expect_equal(calc(blocks)-baseline, direct, tolerance=1e-10)
  }
})

test_that("PD correction vanishes for a regression independent of W", {
  v <- seq(-1, 2, length.out=60)
  s <- rep(1:3, length.out=length(v))
  fits <- rep(list(list(fit=function(z) z$v1j^2)), 3)
  grid <- .get.muhat(s, fits, v, data.frame(w=cos(v)), 7)
  for(debias in c(FALSE, TRUE)) {
    baseline <- .compute.rinfl.func(sin(v), v, 0, 2, 2, dnorm, debias)
    corrected <- .compute.rinfl.func(sin(v), v, 0, 2, 2, dnorm,
                                     debias, grid, v^2)
    expect_equal(corrected, baseline, tolerance=1e-12)
  }
})
