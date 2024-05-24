context("second-stage loc linear regression")

test_that("expected behavior of second stage local linear smoother", {

  for(bw in c(0.05, 0.1, 1, 100)) {
    for(i in 1:30) {
      n <- 500
      A <- rnorm(n)
      pseudo.out <- cos(2*pi*A) + rnorm(n)
      eval.pts <- seq(min(A), max(A), length.out = 100)
      require(locpol)
      fit2 <- suppressWarnings({
        locpol(y ~ x, data=data.frame(x=A, y=pseudo.out), bw=bw, xeval=eval.pts,
               kernel=gaussK)
        })
      est <- rep(NA, length(eval.pts))
      for(j in 1:length(eval.pts)) {
        A.std <- (A - eval.pts[j])/bw
        W <- dnorm(A.std)
        fit4 <- lm(y ~ x, data.frame(x = A.std, y = pseudo.out), weights = W)
        est[j] <- coef(fit4)[1]

      }

      fit <- debiased_inference(A, pseudo.out, eval.pts=eval.pts,
                                bandwidth.method="LOOCV",
                                bw.seq=bw,
                                debias=FALSE,
                                kernel.type="gau")
      expect_true(max(abs(fit2$lpFit$y - fit$res$theta), na.rm = TRUE) < 1e-10)
      expect_true(max(abs(est - fit$res$theta), na.rm = TRUE) < 1e-10)
      print(i)
    }
  }
})
