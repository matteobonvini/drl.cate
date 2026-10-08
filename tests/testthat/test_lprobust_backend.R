native_fixture <- function(...) {
  set.seed(319)
  n <- 160L
  d <- data.frame(v=runif(n,-2,2),w=rep(0:1,n/2),a=rep(c(0,1,1,0),n/4))
  d$y <- sin(d$v)+d$a*(d$v+d$w)+rnorm(n)
  signal <- function(x) rowSums(x)
  stage <- function(pseudo,x,new.x) list(res=cbind(signal(new.x),NA,NA))
  cate(d,"dr",c("v","w"),"y","a",c("v","w"),
    expand.grid(v=c(-1,0,1),w=0:1),
    mu1.x=function(y,a,x,new.x) list(res=signal(new.x)),
    mu0.x=function(y,a,x,new.x) list(res=rep(0,nrow(new.x))),
    pi.x=function(a,x,new.x) list(res=rep(.5,nrow(new.x))),
    drl.v=function(pseudo,v,new.v) stage(pseudo,v,new.v),drl.x=stage,
    foldid=rep(c(1,1,2,2),n/4),univariate_reg=TRUE,partial_dependence=TRUE,
    density.ratio=rep(list(function(v1,v2)
      list(predict=function(new.v1,new.v2) rep(1,length(new.v1)))),2),
    cate.w=rep(list(function(tau,w,new.w)
      list(fit=function(new.w) signal(new.w))),2), ...)
}

test_that("native helper matches nprobust and calls it only once", {
  set.seed(51)
  x <- runif(160,-2,2); y <- sin(x)+rnorm(160); grid <- c(-1,0,1)
  reference <- nprobust::lprobust(y,x,eval=grid,p=1,kernel="gau",
    bwselect="imse-dpi",rho=1,bwcheck=21,vce="nn")$Estimate
  calls <- 0L; native <- nprobust::lprobust
  local_mocked_bindings(lprobust=function(...) {
    calls <<- calls+1L; native(...)
  },.package="nprobust")
  pair <- .continuous.inference(x,y,grid,.resolve.inference("lprobust",NULL))
  expect_equal(calls,1L)
  expect_equal(pair$regular$res$theta,reference[,"tau.us"])
  expect_equal(pair$debiased$res$theta,reference[,"tau.bc"])
  for(debias in c(FALSE,TRUE)) {
    z <- pair[[if(debias) "debiased" else "regular"]]
    se <- reference[,if(debias) "se.rb" else "se.us"]
    expect_equal(z$res$if.val.sd,se)
    expect_equal(z$res$ci.ll.pts,z$res$theta-qnorm(.975)*se)
    expect_equal(z$res$ci.ul.pts,z$res$theta+qnorm(.975)*se)
    expect_equal(z$res$h.effective,reference[,"h"])
    expect_equal(z$res$b.effective,reference[,"b"])
    expect_true(all(is.na(z$res$ci.ll.unif)))
    expect_true(all(is.na(z$res$ci.ul.unif)))
    expect_true(all(is.na(z$res$unif.quantile)))
    expect_null(z$cv); expect_null(z$risk)
  }
  expect_equal(pair$metadata$q,2L)
  expect_false(pair$metadata$uniform.bands)
})

test_that("cate defaults to native inference without CV or extra PD integration", {
  local_mocked_bindings(.local.cv=function(...) stop("unexpected CV"),
    .get.muhat=function(...) stop("unexpected PD integration source"),
    .pd.integrate=function(...) stop("unexpected PD integration"))
  fit <- native_fixture()
  expect_equal(fit$inference,list(method="lprobust",bandwidth.method="imse-dpi"))
  for(z in list(fit$univariate.res$dr[[1]],fit$pd.res$dr[[1]])) {
    reference <- nprobust::lprobust(z$data$pseudo,z$data$exposure,
      eval=c(-1,0,1),p=1,kernel="gau",bwselect="imse-dpi",rho=1,
      bwcheck=21,vce="nn")$Estimate
    expect_equal(z$res$theta,reference[,"tau.us"])
    expect_equal(z$res$theta.debias,reference[,"tau.bc"])
    expect_equal(z$res$bias,z$res$theta-z$res$theta.debias)
    expect_equal(z$inference$q,2L)
    expect_null(z$cv); expect_null(z$cv.debias)
    expect_null(z$risk$risk); expect_null(z$risk$risk.debias)
    expect_true(all(is.finite(z$res$ci.ll.pts.debias)))
    expect_true(all(is.na(z$res$ci.ll.unif.debias)))
  }
  expect_false(fit$pd.res$dr[[1]]$inference$pd.integration.if)
  # Discrete inference still returns the existing one-step/subgroup outputs.
  for(z in list(fit$univariate.res$dr[[2]],fit$pd.res$dr[[2]])) {
    expect_equal(z$res,z$res.empVar)
    expect_true(all(is.finite(z$res$theta)))
    expect_null(z$inference)
  }
})

test_that("previous backend preserves direct inference and CV diagnostics", {
  set.seed(52)
  x <- seq(-2,2,length.out=60); y <- sin(x)+rnorm(60)
  args <- list(A=x,pseudo.out=y,eval.pts=c(-.5,.5),bw.seq=c(.8,1.2),
               cv.eval.size=10L,cv.repeats=2L)
  set.seed(15)
  pair <- do.call(.continuous.inference,c(args,
    list(inference=.resolve.inference("influence-function",NULL))))
  set.seed(15)
  old <- lapply(c(FALSE,TRUE),function(debias)
    do.call(debiased_inference,c(args,list(debias=debias,kernel.type="gau",
                                         bandwidth.method="LOOCV"))))
  expect_equal(pair$regular,old[[1]])
  expect_equal(pair$debiased,old[[2]])
  expect_equal(pair$metadata$q,3L)
  expect_true(pair$metadata$uniform.bands)
  expect_true(all(is.finite(pair$debiased$res$ci.ll.unif)))
})

test_that("backend settings are explicit and incompatible choices fail early", {
  expect_equal(formals(cate)$inference.method,"lprobust")
  expect_error(cate(inference.method="other"),"arg")
  expect_error(cate(bandwidth.method="LOOCV"),"bandwidth.method")
  expect_error(cate(inference.method="influence-function",bandwidth.method="imse-dpi"),
               "bandwidth.method")
  for(method in c("imse-dpi","mse-dpi","imse-rot","mse-rot","ce-dpi","ce-rot"))
    expect_equal(.resolve.inference("lprobust",method)$bandwidth.method,method)
  expect_equal(.resolve.inference("influence-function","LOOCV(h=b)")$bandwidth.method,
               "LOOCV(h=b)")
  fit <- native_fixture(min.local=30,bandwidth.method="imse-rot")
  expect_equal(fit$pd.res$dr[[1]]$inference$bwcheck,30)
  expect_equal(fit$pd.res$dr[[1]]$inference$bandwidth.method,"imse-rot")
  local_mocked_bindings(lprobust=function(...) stop("native fit failed"),.package="nprobust")
  expect_error(.continuous.inference(1:30,1:30,15,
    .resolve.inference("lprobust",NULL)),"native fit failed")
})
