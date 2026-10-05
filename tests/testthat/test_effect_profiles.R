test_that("additive levels are marginal predictions, including factor modifiers", {
  set.seed(401)
  d <- data.frame(age=rnorm(160), age_extra=rnorm(160),
                  group=factor(rep(c("a", "b"), 80)))
  y <- 2 + d$age - .4*d$age_extra + (d$group == "b") + rnorm(160)
  fit <- lm(y ~ splines::bs(age, df=4) + age_extra + group, data=d)
  for (nm in names(d)) {
    pts <- if (is.factor(d[[nm]])) levels(d[[nm]]) else c(-1, 0, 1)
    out <- additive_effect_profile(fit, y, d, nm, pts)
    brute <- vapply(pts, function(pt) {
      nd <- d
      nd[[nm]] <- if (is.factor(d[[nm]])) factor(pt, levels=levels(d[[nm]])) else pt
      mean(predict(fit, nd))
    }, numeric(1))
    expect_equal(out$res$theta, unname(brute), tolerance=1e-10)
    at.sample <- additive_effect_profile(fit, y, d, nm, d[[nm]])
    expect_equal(mean(at.sample$res$theta), mean(y), tolerance=1e-10)
    expect_equal(mean(at.sample$component.res$theta), 0, tolerance=1e-10)
    expect_true(all(out$res$se > 0))
  }
})

test_that("level IF intervals match numerical perturbations of the empirical law", {
  set.seed(402)
  n <- 70
  d <- data.frame(v=rnorm(n), w=rnorm(n))
  d$v <- d$v + d$w
  y <- 2 + d$v - .2*d$v^2 + d$w + rnorm(n)
  pts <- c(-1, 0, 1)
  fit <- lm(y ~ v + I(v^2) + w, data=d)
  gam <- additive_effect_profile(fit, y, d, "v", pts)
  # Independent finite-difference calculation: perturb the empirical law,
  # refit, and marginalize over that same perturbed reference distribution.
  eps <- 1e-6
  perturbed <- vapply(seq_len(n), function(i) {
    weights <- rep((1-eps)/n, n)
    weights[i] <- weights[i] + eps
    f <- lm(y ~ v + I(v^2) + w, data=d, weights=weights)
    vapply(pts, function(pt) {
      nd <- d; nd$v <- pt
      sum(weights * predict(f, nd))
    }, numeric(1))
  }, numeric(length(pts)))
  numerical.if <- (perturbed - gam$res$theta) / eps
  expected.se <- sqrt(apply(numerical.if, 1, var)/n)
  expect_equal(gam$res$se, expected.se, tolerance=2e-5)

  regress <- function(y, x, new.x) predict(lm(y ~ ., data=cbind(y=y, x)), new.x)
  rob <- robinson(y, d["w"], d$v, pts, rep(1, n), regress, regress, dfs=2)
  out <- robinson_effect_profile(rob$model, y, d$v, pts)
  expect_equal(out$res$theta, gam$res$theta, tolerance=1e-10)
  expect_equal(out$res$se, expected.se, tolerance=2e-5)
  expect_gt(out$res$se[2], 0)
  expect_equal(out$res$ci.ul.pts-out$res$theta, qnorm(.975)*out$res$se)
  # Translation of the polynomial origin cannot change the level or its SE.
  shifted <- robinson(y, d["w"], d$v+5, pts+5, rep(1, n), regress, regress, dfs=2)
  translated <- robinson_effect_profile(shifted$model, y, d$v+5, pts+5)
  expect_equal(translated$res$theta, out$res$theta, tolerance=1e-9)
  expect_equal(translated$res$se, out$res$se, tolerance=1e-9)
})

test_that("cate returns levels and retains components for both approximations", {
  set.seed(403)
  n <- 180
  d <- data.frame(v=rnorm(n), w=rep(0:1, n/2), a=rbinom(n, 1, .5))
  d$y <- d$a*(2 + d$v + d$w) + rnorm(n)
  regress <- function(y, x, new.x) predict(lm(y ~ ., data=cbind(y=y, x)), new.x)
  stage <- function(pseudo, x, new.x) {
    fit <- lm(y ~ ., data=cbind(y=pseudo, x))
    list(res=cbind(predict(fit, new.x), NA, NA), model=fit)
  }
  fit <- cate(d, "dr", c("v", "w"), "y", "a", c("v", "w"),
              expand.grid(v=c(-1, 0, 1), w=0:1),
              mu1.x=function(y,a,x,new.x) list(res=2+new.x$v+new.x$w),
              mu0.x=function(y,a,x,new.x) list(res=rep(0,nrow(new.x))),
              pi.x=function(a,x,new.x) list(res=rep(.5,nrow(new.x))),
              drl.v=function(pseudo,v,new.v) stage(pseudo,v,new.v), drl.x=stage,
              nsplits=2, additive_approx=TRUE, partially_linear=TRUE,
              cate.not.j=rep(list(regress),2), reg.basis.not.j=rep(list(regress),2),
              pl.dfs=list(1,1),
              fit.basis.additive=function(y,x,new.x) list(model=lm(y ~ .,data=cbind(y=y,x))))
  psi <- mean(fit$cate.x.res$pseudo$dr)
  for (method in c("additive.res", "robinson.res")) {
    for (j in 1:2) {
      out <- fit[[method]]$dr[[j]]
      expect_equal(out$res$theta, psi + out$component.res$theta)
      expect_true(all(is.finite(out$res$se) & out$res$se > 0))
      expect_equal(out$ate, psi)
      if(method == "additive.res") expect_true(is.data.frame(out$legacy.res))
      else expect_false("legacy.res" %in% names(out))
    }
    expect_gt(fit[[method]]$dr[[2]]$res$theta[1], 1)
  }
})

test_that("unidentified and nonadditive profiles fail explicitly", {
  set.seed(404)
  d <- data.frame(v=rnorm(30), w=rnorm(30))
  y <- rnorm(30)
  bad <- lm(y ~ v*w, data=d)
  expect_error(additive_effect_profile(bad,y,d,"v",0), "additive lm")
  d$w <- d$v
  bad <- lm(y ~ v+w, data=d)
  expect_error(additive_effect_profile(bad,y,d,"v",0), "full-rank")
})
