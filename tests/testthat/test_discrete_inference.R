discrete_fixture <- function(constant=FALSE, factor.v=FALSE, ratio=FALSE,
                             noise=FALSE, split=TRUE) {
  d <- expand.grid(v=0:2,w=c(-1,1),a=0:1,rep=1:8)
  signal <- function(x) if(constant) rep(2,nrow(x)) else {
    as.numeric(as.character(x$v))^2+x$w
  }
  if(factor.v) d$v <- factor(d$v)
  eps <- if(noise) sin(seq_len(nrow(d))) else rep(0,nrow(d))
  d$y <- d$a*signal(d)+eps
  stage <- function(pseudo,x,new.x) list(res=cbind(signal(new.x),NA,NA))
  fit <- cate(d,"dr",c("v","w"),"y","a",c("v","w"),
    expand.grid(v=unique(d$v),w=c(-1,1)),
    mu1.x=function(y,a,x,new.x) list(res=signal(new.x)),
    mu0.x=function(y,a,x,new.x) list(res=rep(0,nrow(new.x))),
    pi.x=function(a,x,new.x) list(res=rep(.5,nrow(new.x))),
    drl.v=function(pseudo,v,new.v) stage(pseudo,v,new.v), drl.x=stage,
    foldid=rep(c(1,1,2,3),length.out=nrow(d)),
    univariate_reg=TRUE, partial_dependence=TRUE,
    sample.split.cond.dens=split,
    cond.dens=lapply(c(1/3,1/2),function(p) {
      force(p)
      function(v1,v2) list(predict.cond.dens=function(v1,v2,new.v1,new.v2)
        rep(p,length(new.v1)))
    }),
    density.ratio=if(ratio) rep(list(function(v1,v2)
      list(predict=function(new.v1,new.v2) rep(1,length(new.v1)))),2) else NULL,
    cate.w=list(
      function(tau,w,new.w) list(fit=function(new.w) {
        names(new.w)[1] <- "v"
        signal(new.w)
      }),
      function(tau,w,new.w) list(fit=function(new.w) {
        names(new.w)[1] <- "w"
        signal(new.w)
      })))
  list(fit=fit,d=d,signal=signal(d))
}

test_that("discrete PD retains uncertainty from averaging other modifiers", {
  for(factor.v in c(FALSE,TRUE)) {
    z <- discrete_fixture(factor.v=factor.v)
    r <- z$fit$pd.res$dr[[1]]$res
    expect_equal(r,z$fit$pd.res$dr[[1]]$res.empVar)
    expect_equal(r$theta,c(0,1,4),tolerance=1e-12)
    se <- (r$ci.ul.pts-r$theta)/qnorm(.975)
    expect_equal(se,rep(sd(z$d$w)/sqrt(nrow(z$d)),3),tolerance=1e-12)
    # Constant oracle scores have zero IF for both estimands.
    cst <- discrete_fixture(constant=TRUE,factor.v=factor.v)$fit
    for(r in list(cst$univariate.res$dr[[1]]$res.empVar,
                  cst$pd.res$dr[[1]]$res)) {
      expect_equal(r$theta,rep(2,3))
      expect_equal(r$ci.ll.pts,r$theta)
      expect_equal(r$ci.ul.pts,r$theta)
    }
  }
})

test_that("primary discrete univariate results use subgroup means and IF intervals", {
  for(factor.v in c(FALSE,TRUE)) {
    z <- discrete_fixture(factor.v=factor.v)
    r <- z$fit$univariate.res$dr[[1]]$res
    expect_equal(r,z$fit$univariate.res$dr[[1]]$res.empVar)
    expect_equal(r$theta,c(0,1,4),tolerance=1e-12)
    phi <- z$fit$cate.x.res$pseudo$dr
    for(i in seq_len(nrow(r))) {
      idx <- z$d$v==r$eval.pts[i]
      influence <- idx/mean(idx)*(phi-r$theta[i])
      se <- sd(influence)/sqrt(length(phi))
      expect_equal(r$ci.ll.pts[i],r$theta[i]-qnorm(.975)*se)
      expect_equal(r$ci.ul.pts[i],r$theta[i]+qnorm(.975)*se)
    }
  }
})

test_that("discrete intervals agree with directly calculated influence functions", {
  for(split in c(FALSE,TRUE)) {
    old <- discrete_fixture(noise=TRUE,split=split)
    new <- discrete_fixture(noise=TRUE,ratio=TRUE,split=split)
    for(j in 1:2) {
      expect_equal(old$fit$pd.res$dr[[j]]$res,new$fit$pd.res$dr[[j]]$res)
      expect_equal(old$fit$pd.res$dr[[j]]$res.empVar,new$fit$pd.res$dr[[j]]$res.empVar)
    }
    phi <- old$fit$cate.x.res$pseudo$dr
    d <- old$d
    n <- nrow(d)
    for(level in 0:2) {
      idx <- d$v==level
      theta <- mean(phi[idx])
      influence <- idx/mean(idx)*(phi-theta)
      r <- old$fit$univariate.res$dr[[1]]$res.empVar[level+1,]
      expect_equal(r$theta,theta)
      expect_equal((r$ci.ul.pts-r$theta)/qnorm(.975),sd(influence)/sqrt(n))
      score <- idx*3*(phi-old$signal)+level^2+d$w
      r <- old$fit$pd.res$dr[[1]]$res[level+1,]
      expect_equal(r$theta,mean(score))
      expect_equal((r$ci.ul.pts-r$theta)/qnorm(.975),sd(score)/sqrt(n))
    }
  }
})
