context("second-stage variance")

test_that("expected behavior of second stage variance estimator", {
  for(i in 1:30) {
    n <- 200
    A <- rnorm(n, sd = 0.1)
    pseudo.out <- cos(2*pi*A) + rnorm(n)
    h <- b <- max(A) - min(A) + 10 # bandwidth is so large that we are doing linear regression
    res <- debiased_inference(A = A, pseudo.out = pseudo.out,
                              bw.seq = h, eval.pts=A,
                              debias=FALSE,
                              kernel.type="uni",
                              bandwidth.method = "LOOCV")
    sigmahat <- (res$res$ci.ul.pts - res$res$theta) / qnorm(0.975) *
      sqrt(n - 1) / sqrt(n)
    fit <- lm(y ~ x, data = data.frame(y = pseudo.out, x = A))
    expect_true(max(abs(fitted(fit) - res$res$theta)) < 1e-10)
    beta.vcov <- sandwich::vcovHC(fit, type = "HC")
    design.mat <- cbind(1, A)
    sigmahat2 <- sqrt(diag(design.mat %*% beta.vcov %*% t(design.mat)))
    expect_true(max(abs(sigmahat - sigmahat2)) < 1e-10)
  }
})
