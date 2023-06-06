context("second-stage loc linear regression")

test_that("expected behavior of second stage local linear smoother", {
  
  for(bw in c(0.1, 0.5, 100)) {
    for(i in 1:30) {
      n <- 100
      A <- rnorm(n)
      pseudo.out <- cos(2*pi*A) + rnorm(n)
      
      # fit2 <- KernSmooth::locpoly(x = A, y = pseudo.out,
      #                             bandwidth = 1,
      #                             gridsize = 100,
      #                             range.x = c(-2, 2),
      #                             truncate = FALSE, 
      #                             bwdisc = 100)
      eval.pts <- seq(min(A), max(A), length.out = 100)
      require(locpol)
      fit2 <- suppressWarnings({locpol(y ~ x, data = data.frame(x = A, y = pseudo.out),
                                               bw = bw, xeval = eval.pts, kernel = gaussK)})
      # fit3 <- tryCatch(list(res = locpol(y ~ x, data = data.frame(x = A, y = pseudo.out), 
      #                         bw = bw, xeval = eval.pts, kernel = gaussK)),
      #                  error=function(e) e, 
      #                  warning=function(w) {
      #                    out = list(w = w, 
      #                               res = locpol(y ~ x, 
      #                                            data = data.frame(x = A, 
      #                                                              y = pseudo.out), 
      #                                            bw = bw, xeval = eval.pts, 
      #                                            kernel = gaussK))
      #                    return(out)
      #                  }
      #                    )
      # print(is(fit3$w, "warning"))
      # fit2 <- fit3$res
      est <- rep(NA, length(eval.pts))
      for(j in 1:length(eval.pts)) {
        A.std <- (A - eval.pts[j])/bw
        W <- dnorm(A.std)
        fit4 <- lm(y ~ x, data.frame(x = A.std, y = pseudo.out), weights = W)
        est[j] <- coef(fit4)[1]
        
      }
      
      fit <- debiased_inference(A, pseudo.out, eval.pts = eval.pts,
                                bandwidth.method = "LOOCV",
                                bw.seq = bw,
                                kernel.type = "gau")
      # print(i)
      # lines(fit$res$eval.pts, fit$res$theta, col = "pink")
      expect_true(max(abs(fit2$lpFit$y - fit$res$theta), na.rm = TRUE) < 1e-10)
      # print(sum(is.na(fit2$lpFit$y)))
      # print(sum(max(abs(fit2$lpFit$y - est), na.rm = TRUE) > 1e-6))
      # expect_true(max(abs(fit2$lpFit$y - est), na.rm = TRUE) < 1e-6 || 
      #               is(fit3, "warning") || 
      #               (which(abs(fit2$lpFit$y - est) > 1e-6) 
      #                %in% c(1, 2, 99, 100)))
    }
  }
})
