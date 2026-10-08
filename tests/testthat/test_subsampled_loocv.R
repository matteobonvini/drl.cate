test_that("subset losses equal full LOOCV and explicit deleted fits", {
  set.seed(98)
  x <- round(runif(50,-1,1),1); y <- sin(x)+rnorm(50,.0,.1)
  ids <- c(30L,1L,45L,12L,17L)
  for(db in c(FALSE,TRUE)) for(k in list(NULL,15L)) {
    full <- .local.loo(x,y,.8,1,db,"gau",k)
    sub <- .local.loo(x,y,.8,1,db,"gau",k,ids)
    brute <- vapply(ids,function(i)
      .local.fit(x[-i],y[-i],x[i],.8,1,function(z).kern(z,"gau"),db,k)$theta,
      numeric(1))
    expect_equal(sub$pred,full$pred[ids],tolerance=1e-10)
    expect_equal(sub$pred,unname(brute),tolerance=1e-8)
    expect_equal(sub$risk,mean((y[ids]-sub$pred)^2))
  }
  expect_error(.local.loo(x,y,1,1,FALSE,eval.indices=c(1,1)),"indices")
  expect_error(.local.loo(x,y,1,1,FALSE,eval.indices=51),"indices")
})

test_that("repeated subsets are reproducible and use median winners", {
  set.seed(21); x <- runif(70,-1,1); y <- sin(4*x)+rnorm(70,sd=.4)
  g <- expand.grid(h=c(.3,.6,1),b=c(.3,.6,1))
  run <- function() .local.cv(x,y,g$h,g$b,TRUE,"gau",NULL,20,5)
  set.seed(24); a <- run()
  set.seed(24); expect_identical(run(),a)
  expect_equal(lengths(a$diagnostics$indices),rep(20L,5))
  expect_true(all(vapply(a$diagnostics$indices,function(i)!anyDuplicated(i),logical(1))))
  ref <- vapply(a$diagnostics$indices,function(ids)
    vapply(seq_len(nrow(g)),function(k)
      .local.loo(x,y,g$h[k],g$b[k],TRUE,"gau",eval.indices=ids)$risk,
      numeric(1)),numeric(nrow(g)))
  expect_equal(a$diagnostics$risk,ref)
  expect_equal(a$risk,rowMeans(ref))
  winners <- apply(ref,2,function(v) {
    tied <- which(abs(v-min(v))<1e-6)
    tied[order(-g$h[tied],-g$b[tied])[1]]
  })
  expect_equal(a$diagnostics$winners$candidate,unname(winners))
  med <- c(h=median(g$h[winners]),b=median(g$b[winners]))
  expect_equal(a$diagnostics$median,med)
  distance <- log(g$h/med[1])^2+log(g$b/med[2])^2
  expect_equal(distance[a$selected],min(distance))
})

test_that("explicit exhaustive CV is unchanged and does not consume sampling RNG", {
  x <- seq(-1,1,length.out=40); y <- sin(x)+cos(5*x)
  for(db in c(FALSE,TRUE)) {
    run <- function(size) debiased_inference(x,y,db,eval.pts=c(-.5,0,.5),
      bw.seq=c(.5,1),kernel.type="gau",bandwidth.method="LOOCV",unif=FALSE,cv.eval.size=size)
    set.seed(5); state <- .Random.seed
    full <- run(NULL)
    expect_identical(.Random.seed,state)
    expect_equal(run(40),full)
    expect_equal(run(100),full)
    expect_true(full$cv$exhaustive)
    set.seed(4); sub <- run(12)
    selected <- debiased_inference(x,y,db,eval.pts=c(-.5,0,.5),
      bw.seq=unique(c(sub$res$h[1],sub$res$b[1])),kernel.type="gau",unif=FALSE,
      bandwidth.method="LOOCV",inference.all=TRUE)
    matches <- Filter(function(z) z$h[1]==sub$res$h[1] && z$b[1]==sub$res$b[1],selected$res.list)
    cols <- setdiff(names(sub$res),"loocv.risk")
    expect_equal(sub$res[cols],matches[[1]][cols],tolerance=1e-10)
  }
})

test_that("default CV uses five subsets of 5000 except on small samples", {
  for(fn in list(cate,debiased_inference)) {
    expect_identical(formals(fn)$cv.eval.size,5000L)
    expect_identical(formals(fn)$cv.repeats,5L)
  }
  # Test the default sampling policy without thousands of local fits.
  # Numerical deletion and full-data inference are exercised separately above.
  local_mocked_bindings(.local.loo=function(x,y,h,b,debias,kernel.type,
                                           min.local,eval.indices) {
    list(pred=y[eval.indices],risk=0)
  })
  run <- function(n,...) {
    x <- seq(-1,1,length.out=n)
    debiased_inference(x,sin(x),TRUE,eval.pts=0,bw.seq=1,
      bandwidth.method="LOOCV",kernel.type="gau",unif=FALSE,...)
  }
  set.seed(321); sampled <- run(5200)
  expect_false(sampled$cv$exhaustive)
  expect_equal(lengths(sampled$cv$indices),rep(5000L,5))
  set.seed(321); expect_identical(run(5200),sampled)
  expect_true(run(5200,cv.eval.size=NULL)$cv$exhaustive)
  set.seed(321); state <- .Random.seed
  small <- run(5000)
  expect_true(small$cv$exhaustive)
  expect_equal(lengths(small$cv$indices),5000L)
  expect_identical(.Random.seed,state)
  expect_equal(small,run(5000,cv.eval.size=NULL))
})

test_that("only requested local systems are built", {
  original <- .local.fit; calls <- 0L
  local_mocked_bindings(.local.fit=function(...) {
    calls <<- calls+1L
    original(...)
  })
  x <- seq(-1,1,length.out=2000)
  z <- .local.loo(x,sin(x),1,1,TRUE,"gau",eval.indices=c(1,200,1000,1400,2000))
  expect_equal(calls,5L)
  expect_true(all(is.finite(z$pred)))
})

test_that("cate forwards subset controls without extra nuisance fits", {
  d <- expand.grid(v=seq(-1,1,length.out=20),w=c(-1,1),a=0:1)
  d$y <- d$a*(d$v+d$w)
  calls <- 0L
  signal <- function(x) rowSums(x)
  stage <- function(pseudo,x,new.x) list(res=cbind(signal(new.x),NA,NA))
  set.seed(171)
  fit <- cate(d,"dr",c("v","w"),"y","a",c("v","w"),
    expand.grid(v=c(-.5,.5),w=c(-1,1)),
    mu1.x=function(y,a,x,new.x) {calls <<- calls+1L;list(res=signal(new.x))},
    mu0.x=function(y,a,x,new.x) list(res=rep(0,nrow(new.x))),
    pi.x=function(a,x,new.x) list(res=rep(.5,nrow(new.x))),
    drl.v=function(pseudo,v,new.v) stage(pseudo,v,new.v), drl.x=stage,
    foldid=rep(1:2,length.out=nrow(d)),univariate_reg=TRUE,partial_dependence=TRUE,
    bw.stage2=list(c(.5,1),NULL),cv.eval.size=12,cv.repeats=5,
    inference.method="influence-function",
    density.ratio=rep(list(function(v1,v2)
      list(predict=function(new.v1,new.v2) rep(1,length(new.v1)))),2),
    cate.w=rep(list(function(tau,w,new.w) list(fit=function(new.w) signal(new.w))),2))
  expect_equal(calls,2L)
  for(z in list(fit$univariate.res$dr[[1]],fit$pd.res$dr[[1]])) {
    expect_equal(lengths(z$cv$indices),rep(12L,5))
    expect_equal(lengths(z$cv.debias$indices),rep(12L,5))
    expect_true(all(is.finite(z$res$theta)))
  }
})

test_that("invalid controls and failed candidates are explicit", {
  for(v in list(0,-1,1.5,NA,Inf,"10",c(1,2))) {
    expect_error(.validate.local.cv(v,5),"cv.eval.size")
    expect_error(.validate.local.cv(10,v),"cv.repeats")
  }
  expect_error(.local.cv(rep(1,20),1:20,1,1,TRUE,"gau",NULL,5,5),
               "All bandwidth candidates")
  x <- rep(c(0,1,2),each=5)
  ids <- c(1L,8L,15L)
  expect_true(all(is.na(.local.loo(x,x,.01,.01,TRUE,"epa",eval.indices=ids)$pred)))
})

test_that("CV repeats do not multiply PD inference predictions", {
  x <- seq(-1,1,length.out=40); y <- sin(x)
  calls <- 0L
  fit <- list(fit=function(z) {calls <<- calls+1L; z$v1j+z$w})
  src <- .get.muhat(rep(1:2,each=20),rep(list(fit),2),x,
                    data.frame(w=cos(x)),max.n.integral=10)
  set.seed(15)
  z <- debiased_inference(x,y,TRUE,eval.pts=c(-.5,0,.5),bw.seq=c(.8,1),
    bandwidth.method="LOOCV",kernel.type="gau",unif=FALSE,
    muhat.vals=src,mhat.obs=x,cv.eval.size=8,cv.repeats=5)
  expect_equal(calls,8L)
  expect_equal(length(z$cv$indices),5L)
})
