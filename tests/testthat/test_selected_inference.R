test_that("selected-only inference preserves estimates and avoids repeated PD prediction", {
  set.seed(814)
  n <- 40
  v <- sort(runif(n,-1,1))
  w <- data.frame(w=rnorm(n))
  s <- rep(1:2,length.out=n)
  calls <- 0L
  fit <- list(fit=function(z) {
    calls <<- calls+1L
    z$v1j^2+z$w
  })
  src <- .get.muhat(s,rep(list(fit),2),v,w,max.n.integral=7)
  m <- v^2+ave(w$w,s)
  y <- sin(v)+w$w
  run <- function(debias, all, unif=FALSE) debiased_inference(
    v,y,debias,eval.pts=c(-.7,0,.7),bw.seq=c(1,1.5),
    bandwidth.method="LOOCV",kernel.type="gau",unif=unif,
    bootstrap=20,muhat.vals=src,mhat.obs=m,inference.all=all,cv.eval.size=NULL)
  for(debias in c(FALSE,TRUE)) {
    calls <- 0L
    slow <- run(debias,TRUE)
    slow.calls <- calls
    calls <- 0L
    fast <- run(debias,FALSE)
    candidates <- if(debias) 4L else 2L
    expect_equal(calls,2L*ceiling(20/7)^2)
    expect_equal(slow.calls,calls*(candidates+1L))
    expect_equal(fast$risk,slow$risk,tolerance=1e-12)
    expect_equal(fast$res,slow$res,tolerance=1e-12)
    for(i in seq_along(fast$res.list)) {
      cols <- c("eval.pts","theta","h","b","loocv.risk")
      expect_equal(fast$res.list[[i]][cols],slow$res.list[[i]][cols])
      selected <- fast$res.list[[i]]$h[1]==fast$res$h[1] &&
                  fast$res.list[[i]]$b[1]==fast$res$b[1]
      expect_equal(all(is.na(fast$res.list[[i]]$ci.ul.pts)),!selected)
    }
    # Match the Gaussian draws used for the final selected fit in old mode.
    set.seed(123)
    slow <- run(debias,TRUE,TRUE)
    set.seed(123)
    invisible(rnorm(n*20*candidates))
    fast <- run(debias,FALSE,TRUE)
    expect_equal(fast$res,slow$res,tolerance=1e-12)
  }
})
