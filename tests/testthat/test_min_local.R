test_that("componentwise deletion matches independent local QR refits", {
  ref <- function(x,y,a,h,b,debias,kernel,k) {
    bw <- function(v) {
      if(is.null(k)) return(v)
      r <- sort(abs(x-a))[k]
      if(r>=v) r*(1+sqrt(.Machine$double.eps)) else v
    }
    h <- bw(h); b <- bw(b)
    u <- (x-a)/h; z <- (x-a)/b
    w <- .kern(u,kernel)/h; wb <- .kern(z,kernel)/b
    X <- cbind(1,u); B <- cbind(1,z,z^2,z^3)
    lin <- stats::lm.wfit(X,y,w)$coefficients
    if(anyNA(lin)) return(NA_real_)
    if(!debias) return(unname(lin[1]))
    cub <- stats::lm.wfit(B,y,wb)$coefficients
    c2 <- stats::lm.wfit(X,u^2,w)$coefficients[1]
    if(anyNA(cub)) return(NA_real_)
    unname(lin[1]-(h/b)^2*c2*cub[3])
  }
  set.seed(142)
  for(ties in c(FALSE,TRUE)) {
    x <- sort(runif(45,-1,1)); if(ties) x <- round(x,1)
    y <- sin(x)+rnorm(length(x),sd=.2)
    for(kernel in c("epa","tri","uni","gau"))
      for(k in list(NULL,15L)) for(debias in c(FALSE,TRUE))
        for(hb in list(c(.9,1.2),c(1.2,.9),c(.05,.05))) {
          actual <- .local.loo(x,y,hb[1],hb[2],debias,kernel,k)
          expected <- vapply(seq_along(x),function(i)
            ref(x[-i],y[-i],x[i],hb[1],hb[2],debias,kernel,k),numeric(1))
          # Numerically unresolved Gaussian systems may be rejected rather
          # than extrapolated from effectively zero weights.
          valid <- is.finite(actual$pred) & is.finite(expected)
          expect_equal(actual$pred[valid],expected[valid],tolerance=1e-6)
          if(hb[1]>.1 || !is.null(k)) expect_equal(is.na(actual$pred),is.na(expected))
          expect_equal(actual$risk,if(anyNA(actual$pred)) Inf else mean((y-actual$pred)^2))
        }
  }
})

test_that("local failure is pointwise and minimum counts observations, not ranks", {
  x <- rep(c(0,1,2),each=5); y <- x^2
  expect_equal(.local.bandwidth(x,0,.01,10),1+sqrt(.Machine$double.eps))
  expect_warning(fit <- .lprobust(x,y,.01,.01,FALSE,c(0,1)),"2 local")
  expect_true(all(is.na(fit[,2])))
  expect_true(all(is.finite(.lprobust(x,y,.01,.01,FALSE,c(0,1),min.local=10)[,2])))
  expect_warning(fit <- .lprobust(x,y,3,3,TRUE,1,min.local=10),"1 local")
  expect_true(is.na(fit[1,2]))
  for(k in list(0,-1,1.5,NA,Inf,c(1,2),"10"))
    expect_error(.validate.min.local(k),"positive integer")
})

test_that("inference preserves valid points and rejects invalid CV candidates", {
  set.seed(11); x <- seq(-1,1,length.out=40); y <- x^2+rnorm(40)
  expect_warning(out <- debiased_inference(x,y,TRUE,eval.pts=c(0,100),
    bandwidth.method="LOOCV(h=b)",bw.seq=c(.001,1),bootstrap=20),"1 evaluation")
  expect_true(is.finite(out$res$theta[1]))
  expect_true(is.finite(out$res$ci.ll.pts[1]))
  expect_true(is.na(out$res$theta[2]))
  expect_true(all(is.na(out$res$ci.ll.unif)))
  expect_equal(out$risk$loocv.risk[1],Inf)
  expect_error(debiased_inference(x,y,TRUE,bw.seq=.001,
    bandwidth.method="LOOCV"),"All bandwidth")
  a <- debiased_inference(x,y,TRUE,eval.pts=0,min.local=15,
    bandwidth.method="LOOCV(h=b)",bw.seq=.01,unif=FALSE)
  b <- debiased_inference(x,y,TRUE,eval.pts=c(-.5,0,.5),min.local=15,
    bandwidth.method="LOOCV(h=b)",bw.seq=.01,unif=FALSE)
  expect_equal(a$risk,b$risk)
  expect_gt(a$res$h.effective,.01)
  expect_equal(a$res$theta,b$res$theta[2])
})
