test_that("streamed PD integrals agree with dense inference", {
  set.seed(815)
  n <- 48
  v <- sort(runif(n,-1,2))
  w <- data.frame(w=rnorm(n))
  s <- rep(1:3, length.out=n)
  calls <- 0L
  largest <- 0L
  fits <- lapply(1:3, function(k) {
    force(k)
    list(fit=function(z) {
      calls <<- calls+1L
      largest <<- max(largest,nrow(z))
      sin(z$v1j)*z$w + k*z$v1j^2
    })
  })
  dense <- .get.muhat(s,fits,v,w,stream=FALSE)
  m <- numeric(n)
  for(k in 1:3) {
    idx <- which(s==k)
    m[idx] <- sin(v[idx])*mean(w$w[idx])+k*v[idx]^2
  }
  y <- cos(v)+w$w
  for(size in c(1,5,100)) {
    calls <- largest <- 0L
    source <- .get.muhat(s,fits,v,w,max.n.integral=size)
    expect_equal(calls,0L)
    expect_s3_class(source,"pd_integration_source")
    for(debias in c(FALSE,TRUE)) {
      pts <- c(-.7,.3,1.8)
      batched <- .pd.moments(source,v,pts,1.2,1.5,dnorm,debias,m)
      for(i in seq_along(pts)) {
        ref <- .compute.rinfl.func(y,v,pts[i],1.2,1.5,dnorm,debias,dense,m)
        got <- .compute.rinfl.func(y,v,pts[i],1.2,1.5,dnorm,debias,batched[[i]],m)
        expect_equal(got,ref,tolerance=1e-10)
      }
    }
    expect_lte(largest,size^2)
  }
  # Full inference path, including risk tables and randomized uniform bands.
  source <- .get.muhat(s,fits,v,w,max.n.integral=5)
  for(debias in c(FALSE,TRUE)) {
    run <- function(input) {
      set.seed(512)
      debiased_inference(v,y,debias,eval.pts=c(-.7,.3,1.8),
                        bw.seq=c(1.2,1.5),bandwidth.method="LOOCV",
                        kernel.type="gau",bootstrap=30,
                        muhat.vals=input,mhat.obs=m)
    }
    expect_equal(run(source),run(dense),tolerance=1e-9)
  }
  calls <- 0L
  .pd.moments(source,v,c(-.7,.3,1.8),1.2,1.5,dnorm,TRUE,m)
  expect_equal(calls,3L*ceiling(16/5)^2) # Not multiplied by evaluation points.
  # The same enlarged bandwidth must enter the estimator, IF and PD weights.
  for(debias in c(FALSE,TRUE)) {
    pts <- c(-.7,.3,1.8)
    batch <- .pd.moments(source,v,pts,.01,.02,dnorm,debias,m,min.local=15)
    for(i in seq_along(pts)) {
      hh <- .local.bandwidth(v,pts[i],.01,15)
      bb <- if(debias) .local.bandwidth(v,pts[i],.02,15) else .02
      ref <- .compute.rinfl.func(y,v,pts[i],hh,bb,dnorm,debias,dense,m)
      got <- .compute.rinfl.func(y,v,pts[i],.01,.02,dnorm,debias,batch[[i]],m,min.local=15)
      expect_equal(got,ref,tolerance=1e-9)
    }
  }
})

test_that("PD source creation is lazy even for large folds", {
  n <- 120000L
  s <- rep(1:4,length.out=n)
  fits <- rep(list(list(fit=function(z) stop("must not predict"))),4)
  source <- .get.muhat(s,fits,seq_len(n),data.frame(w=seq_len(n)))
  expect_s3_class(source,"pd_integration_source")
  expect_lt(as.numeric(object.size(source)),10*1024^2)
})
