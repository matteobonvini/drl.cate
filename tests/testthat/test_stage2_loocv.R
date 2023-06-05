context("second-stage LOOCV")

test_that("expected behavior of second stage bw selector", {
  for(i in 1:30) {
    n <- 200
    A <- rnorm(n, sd = 0.1)
    pseudo.out <- cos(2*pi*A) + rnorm(n)
    h <- b <- 0.95
    res <- debiased_inference(A = A, pseudo.out = pseudo.out, 
                              bw.seq = h, eval.pts=A, 
                              kernel.type="epa", 
                              bandwidth.method = "LOOCV")
    idx <- which(!is.na(res$res$theta))
    if(sum(is.na(res$res$theta)) > 0) {
      print("Skipping")
      next()
    }
    preds <- rep(NA, n)
    for(j in 1:n) {
      if(is.na(res$res$theta[j])) next()
      preds[j] <- .lprobust(x = A[-j], 
                            y = pseudo.out[-j], 
                            h = h, b = b, 
                            eval.pt = A[j], kernel.type="epa")[, "mu.hat"]
      
    }
    loocv.risk <- ifelse(mean(!is.na(preds)) <= 0.8,
                         Inf, mean((pseudo.out[idx] - preds[idx])^2))
   expect_true(abs(loocv.risk - res$risk$loocv.risk) < 1e-10 ||
               is.infinite(loocv.risk))
   # print(loocv.risk)
   # print(abs(loocv.risk - res$risk$loocv.risk))
   # print(i)
  }
  
  
})
