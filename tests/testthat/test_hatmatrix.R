context("non-debiased inference helpers")

test_that(".hatmatrix is correct when the bandwidth is huge (linear regression)", {
  n <- 100
  nsim <- 500
  for(i in 1:nsim) {
    x <- rnorm(n)
    y <- 2*x + rnorm(n)
    fit <- lm(y ~ x)
    hatvals1 <- .hatmatrix(x, y, h=1000, b=1000, eval.pt=x, kernel.type="uni",
                           debias=FALSE)
    hatvals2 <- hatvalues(fit)
    expect_true(max(abs(hatvals1 - hatvals2)) < 1e-12)

    loocv.risk1 <- .robust.loocv(x, y, h=1000, b=1000, FALSE, x, "uni")
    loocv.risk2 <- mean((resid(fit) / (1-hatvals2))^2)
    expect_true(max(abs(loocv.risk1 - loocv.risk2)) < 1e-12)
    print(i)
  }

})


test_that(".hatmatrix is correct when the bandwidth is small", {

  n <- 1000
  nsim <- 2
  h <- b <- 0.5
  for(i in 1:nsim) {
    x <- rnorm(n)
    y <- 2*x + rnorm(n)
    fit <- lm(y ~ x)

    bw.min <- Vectorize(function(u){
      # count the unique points nearby u
      # length(unique(A[.kern((A-u)/min(h.seq), "gau") > 1e-6]))
      sort(unique(abs(x - u)))[10]
    })
    h <- max(h, max(bw.min(x)))
    b <- max(b, max(bw.min(x)))

    # choose gaussian kernel for simplicity so that weighted lm does not discard hatvalues
    # for weights = 0.
    hatvals1 <- .hatmatrix(x, y, h=h, b=h, eval.pt=x, kernel.type="gau",
                           debias=FALSE)
    hatvals2 <- hatvalues(fit)
    # because the bandwidth is small we should have different hatvals
    expect_true(max(abs(hatvals1 - hatvals2)) > 1e-2)
    hatvals3 <- err <- rep(NA, n)
    for(j in 1:n){

      ww <- .kern((x-x[j])/h, "gau")/h
      fit.j <- lm(y ~ x, weights = ww)
      fit.not.j <- lm(y[-j] ~ x[-j], weights = ww[-j])
      hatvals3[j] <- hatvalues(fit.j)[j]
      err[j] <- y[j] - sum(c(1, x[j]) * coef(fit.not.j))

    }

    loocv.risk1 <- .robust.loocv(x, y, h=h, b=h, FALSE, x, "gau")
    loocv.risk3 <- mean(err^2)
    expect_true(max(abs(hatvals1 - hatvals3)) < 1e-10)
    expect_true(abs(loocv.risk1 - loocv.risk3) < 1e-10)

    print(i)
  }

})

test_that("we get the right *optimal* bw between 0.1 and 1000", {

  n <- 500
  nsim <- 5
  for(i in 1:nsim) {
    x <- runif(n, -1, 1)

    # the optimal bandwidth should be small in this model
    y <- cos(2*pi*x) + rnorm(n, sd = 0.1)

    # check that this is the case
    fit <- lm(y ~ x)
    loocv.risk2 <- mean((resid(fit) / (1-hatvalues(fit)))^2)

    hatvals3 <- err <- preds.wlm <- rep(NA, n)
    h <- b <- 0.1
    bw.min <- Vectorize(function(u){
      # count the unique points nearby u
      # length(unique(A[.kern((A-u)/min(h.seq), "gau") > 1e-6]))
      sort(unique(abs(x - u)))[10]
    })
    h <- max(h, max(bw.min(x)))
    b <- max(b, max(bw.min(x)))
    for(j in 1:n){

      ww <- .kern((x-x[j])/h, "gau")/h
      fit.j <- lm(y ~ x, weights = ww)
      preds.wlm[j] <- sum(c(1, x[j]) * coef(fit.j))
      fit.not.j <- lm(y[-j] ~ x[-j], weights = ww[-j])
      hatvals3[j] <- hatvalues(fit.j)[j]
      err[j] <- y[j] - sum(c(1, x[j]) * coef(fit.not.j))

    }
    loocv.risk3 <- mean(err^2)
    expect_true(loocv.risk3 < loocv.risk2)
    loocv.risk1 <- .robust.loocv(x, y, h=h, b=h, FALSE, x, "gau")

    expect_true(max(abs(loocv.risk1$loocv.risk - loocv.risk3)) < 1e-10)

    # en passant, check the predictions are correct
    preds <- .lprobust(x, y, h, h, FALSE, x, "gau")
    expect_true(max(abs(preds[, 2] - preds.wlm))<1e-10)

    loocv.risk45 <- debiased_inference(A=x, pseudo.out=y, debias=FALSE,
                                       eval.pts=x, bw.seq=c(h, 1000),
                                       bandwidth.method="LOOCV",
                                       kernel.type="gau")$risk
    loocv.risk4 <- loocv.risk45$loocv.risk[which(loocv.risk45$h==h)]
    loocv.risk5 <- loocv.risk45$loocv.risk[which(loocv.risk45$h==1000)]

    expect_true(max(abs(loocv.risk3 - loocv.risk4)) < 1e-10)
    expect_true(max(abs(loocv.risk2 - loocv.risk5)) < 1e-7)
    print(i)
  }
})


test_that("we get the right *optimal* bw more refined", {

  n <- 500
  nsim <- 5
  for(i in 1:nsim) {
    x <- runif(n, -1, 1)

    # the optimal bandwidth should be small in this model
    y <- cos(pi*x) + rnorm(n, sd = 0.2)

    # check that this is the case
    fit <- lm(y ~ x)
    loocv.risk2 <- mean((resid(fit) / (1-hatvalues(fit)))^2)
    h.seq <- c(0.01, 0.05, 0.1, 0.15, 0.2)
    hatvals3 <- err <- preds.wlm <- matrix(NA, nrow=n, ncol=length(h.seq))
    bw.min <- Vectorize(function(u){
      # count the unique points nearby u
      # length(unique(A[.kern((A-u)/min(h.seq), "gau") > 1e-6]))
      sort(unique(abs(x - u)))[10]
    })
    h.seq <- pmax(h.seq, max(bw.min(x)))
    for(k in 1:length(h.seq)){
      h <- h.seq[k]
      for(j in 1:n){
        ww <- .kern((x-x[j])/h, "gau")/h
        fit.j <- lm(y ~ x, weights = ww)
        preds.wlm[j, k] <- sum(c(1, x[j]) * coef(fit.j))
        fit.not.j <- lm(y[-j] ~ x[-j], weights = ww[-j])
        hatvals3[j, k] <- hatvalues(fit.j)[j]
        err[j, k] <- y[j] - sum(c(1, x[j]) * coef(fit.not.j))
      }
    }
    loocv.risk3 <- apply(err, 2, function(x) mean(x^2))

    loocv.risk4 <- debiased_inference(A=x, pseudo.out=y, debias=FALSE, eval.pts=x,
                                      bw.seq=h.seq, bandwidth.method="LOOCV",
                                      kernel.type="gau")$risk

    expect_true(max(abs(loocv.risk3 - loocv.risk4$loocv.risk)) < 1e-10)
    print(i)
  }
})
