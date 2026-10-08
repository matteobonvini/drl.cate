ratio_fixture <- function(ratios=NULL, split=TRUE, min.local=NULL) {
  set.seed(501) # Uniform bands use a randomized Gaussian-process calculation.
  n <- 80
  d <- data.frame(v=seq(-2,2,length.out=n), w=rep(0:1,n/2),
                  a=rep(c(0,1,1,0),n/4))
  d$y <- sin(d$v) + d$a
  zero_stage <- function(pseudo,x,new.x) list(res=cbind(rep(0,nrow(new.x)),NA,NA))
  cond <- function(v1,v2) list(predict.cond.dens=function(v1,v2,new.v1,new.v2)
    rep(.5,length(new.v1)))
  cate(d,"dr",c("v","w"),"y","a",c("v","w"),
       expand.grid(v=c(-1,0,1),w=0:1),
       mu1.x=function(y,a,x,new.x) list(res=rep(0,nrow(new.x))),
       mu0.x=function(y,a,x,new.x) list(res=rep(0,nrow(new.x))),
       pi.x=function(a,x,new.x) list(res=rep(.5,nrow(new.x))),
       drl.v=function(pseudo,v,new.v) zero_stage(pseudo,v,new.v),
       drl.x=zero_stage, foldid=rep(c(1,1,2,2),n/4),
       partial_dependence=TRUE, sample.split.cond.dens=split,
       cond.dens=rep(list(cond),2),
       cate.w=rep(list(function(tau,w,new.w) list(fit=function(new.w) rep(0,nrow(new.w)))),2),
       bw.stage2=list(.7,NULL), density.ratio=ratios, min.local=min.local,
       inference.method="influence-function")
}

test_that("cate forwards min.local to continuous univariate and PD fits", {
  fit <- ratio_fixture(min.local=50)
  for(res in list(fit$univariate.res$dr[[1]]$res,fit$pd.res$dr[[1]]$res)) {
    expect_true(all(res$h.effective>.7))
    expect_true(all(res$h.debias.effective>.7))
    expect_true(all(is.finite(res$theta.debias)))
  }
})

test_that("custom ratios use training folds and bypass density floors", {
  seen <- list()
  factory <- function(v1,v2) {
    seen[[length(seen)+1L]] <<- v1
    list(predict=function(new.v1,new.v2) {
      expect_length(intersect(v1,new.v1),0)
      rep(1e-4,length(new.v1))
    })
  }
  fit <- ratio_fixture(list(factory,NULL))
  expect_length(seen,2)
  expect_equal(lengths(seen),c(40L,40L))
  pd <- fit$pd.res$dr[[1]]$data
  expect_equal(pd$density.ratio,rep(1e-4,80))
  expect_equal(pd$pseudo,fit$cate.x.res$pseudo$dr*1e-4)
  expect_true(all(is.na(pd$cond.dens.vals)))
  expect_true(all(is.finite(fit$pd.res$dr[[1]]$res$theta.debias)))
  expect_false("density.ratio" %in% names(fit$pd.res$dr[[2]]$data))
})

test_that("custom discrete ratios agree with legacy weights in both split modes", {
  for (split in c(TRUE,FALSE)) {
    calls <- 0L
    factory <- function(v1,v2) {
      calls <<- calls+1L
      expect_length(v1,if (split) 40L else 80L)
      list(predict=function(new.v1,new.v2) rep(1,length(new.v1)))
    }
    old <- ratio_fixture(split=split)
    new <- ratio_fixture(list(NULL,factory),split=split)
    expect_equal(calls,if (split) 2L else 1L)
    expect_equal(new$pd.res$dr[[2]]$res,old$pd.res$dr[[2]]$res)
    expect_equal(new$pd.res$dr[[2]]$res.empVar,old$pd.res$dr[[2]]$res.empVar)
    expect_equal(new$pd.res$dr[[1]]$res,old$pd.res$dr[[1]]$res)
  }
})

test_that("invalid custom ratios fail clearly instead of returning placeholders", {
  expect_error(ratio_fixture(list(NULL)),"one function or NULL")
  for (bad in list(-1,NA_real_,Inf,"bad",numeric(0),matrix(1,80,1))) {
    factory <- function(v1,v2) list(predict=function(new.v1,new.v2) bad)
    expect_error(ratio_fixture(list(factory,NULL)),"finite, nonnegative numeric vector")
    expect_error(ratio_fixture(list(factory,NULL),split=FALSE),"finite, nonnegative numeric vector")
  }
})
