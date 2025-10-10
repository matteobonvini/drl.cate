context("non-debiased inference helpers")

test_that("Partial dependence with independent covariates and effect modifiers", {
  n <- 5000
  X <- matrix(rnorm(n*3), n, 3); colnames(X) <- paste0("X", 1:3)
  V <- cbind(X[,1], X[, 2])                        # effect modifier
  A <- rbinom(n, 1, plogis(0.4*X[,2]-0.3*X[,3]))
  tau_true <- 1 + V%*%c(0.5, -0.8)               # linear CATE in V
  Y <- tau_true*A - 0.5*X[,3] + rnorm(n, 0, 1)
  
  df <- data.frame(Y=Y, A=A, X)
  
  # Use simple nuisance fits that are reasonably well-specified
  # (or your package's default learners). Keep this minimal & fast.
  v0 <- matrix(c(seq(-2, 2, length.out = 9),seq(-2, 2, length.out = 9)), nrow=9, ncol = 2)
  colnames(v0) <- c("X1", "X2")
  
  pi.x <- function(a, x, new.x) {
    SL.lib <- c("SL.mean", "SL.lm", "SL.glm")
    model <- SuperLearner::SuperLearner(Y = a,
                                        X = x,
                                        newX = new.x,
                                        family = binomial(),
                                        SL.library = SL.lib)
    fit <- function(new.x) {
      out <- predict.SuperLearner(model, newdata = new.x)$pred[, 1]
      out[out < 0.01] <- 0.01
      out[out > 0.99] <- 0.99
      return(out)
    }
    preds.tr <- model$SL.predict[, 1]
    preds.tr[preds.tr < 0.01] <- 0.01
    preds.tr[preds.tr > 0.99] <- 0.99
    return(list(res = preds.tr, model = model, fit = fit))
  }
  
  mu1.x <- function(y, a, x, new.x) {
    SL.lib <- c("SL.mean", "SL.lm", "SL.glm")
    model <- SuperLearner::SuperLearner(Y = y[a == 1],
                                        X = x[a == 1, , drop = FALSE],
                                        newX = new.x,
                                        SL.library = SL.lib)
    fit <- function(new.x) predict.SuperLearner(model, newdata = new.x)$pred[, 1]
    return(list(res = model$SL.predict[, 1], model = model, fit = fit))
  }
  
  mu0.x <- function(y, a, x, new.x) {
    SL.lib <- c("SL.mean", "SL.lm", "SL.glm")
    model <- SuperLearner::SuperLearner(Y = y[a == 0],
                                        X = x[a == 0, , drop = FALSE],
                                        newX = new.x,
                                        SL.library = SL.lib)
    fit <- function(new.x) predict.SuperLearner(model, newdata = new.x)$pred[, 1]
    return(list(res = model$SL.predict[, 1], model = model, fit = fit))
  }
  
  drl.x <- function(pseudo, x, new.x) {
    SL.lib <- c("SL.mean", "SL.lm", "SL.glm")
    model <- SuperLearner::SuperLearner(Y = pseudo,
                                        X = x,
                                        newX = new.x,
                                        family = gaussian(),
                                        SL.library = SL.lib)
    fit <- function(new.x) predict.SuperLearner(model, newdata = new.x)$pred[, 1]
    out <- cbind(model$SL.predict[, 1], NA, NA)
    return(list(res = out, model = model, fit = fit))
  }
  
  drl.v <- function(pseudo, v, new.v) {
    SL.lib <- c("SL.mean", "SL.lm", "SL.glm")
    model <- SuperLearner::SuperLearner(Y = pseudo,
                                        X = v,
                                        newX = new.v,
                                        family = gaussian(),
                                        SL.library = SL.lib)
    fit <- function(new.x) predict.SuperLearner(model, newdata = new.v)$pred[, 1]
    out <- cbind(model$SL.predict[, 1], NA, NA)
    return(list(res = out, model = model, fit = fit))
  }
  
  cond.dens.V1 <- function(v1, v2) {
    # v1 = age; v2 = all V (effect modifiers) except for age
    # semiparametric model: V1 = E(V1 | V2) + \sigma(V2) * eps
    # E(V1 | V2) and \sigma^2(V2) are estimated via GAM
    # eps is univariate and its density is estimated by KDE.
    dat <- data.frame(X1=v1, X2=v2$X2)
    mean.fit <- mgcv::gam(X1 ~ s(X2), data=dat)
    preds.means <- mgcv::predict.gam(mean.fit, newdata=dat)
    
    dat.var <- dat
    dat.var$eps2 <- resid(mean.fit)^2
    
    var.fit <- mgcv::gam(eps2 ~ s(X2),
                         data=dat.var, family=gaussian(link=log))
    preds.vars <- mgcv::predict.gam(var.fit, type="response", newdata=dat.var)
    
    v1.std.tr <- resid(mean.fit)/sqrt(preds.vars)
    dens.est <- ks::kde(x=as.matrix(v1.std.tr), eval.points=v1.std.tr, density=TRUE)
    res.tr <- dens.est$estimate/sqrt(preds.vars)
    
    predict.cond.dens <- function(v1, v2, new.v1, new.v2) {
      # evaluate the estimated conditional density f(v1 | v2) at f(new.v1 | new.v2)
      new.dat <- data.frame(X2=new.v2$X2)
      new.preds.means <- mgcv::predict.gam(mean.fit,
                                           newdata=as.data.frame(new.dat),
                                           type="response")
      new.preds.vars <- mgcv::predict.gam(var.fit,
                                          newdata=as.data.frame(new.dat),
                                          type="response")
      v1.std.te <- (new.v1-new.preds.means)/sqrt(new.preds.vars)
      dens.est <- ks::kde(x=as.matrix(v1.std.te), eval.points=v1.std.te, density=TRUE)
      res <- dens.est$estimate/sqrt(new.preds.vars)
      return(res)
    }
    out <- list(predict.cond.dens=predict.cond.dens)
    return(out)
  }
  
  cond.dens.V2 <- function(v1, v2) {
    # v1 = age; v2 = all V (effect modifiers) except for age
    # semiparametric model: V1 = E(V1 | V2) + \sigma(V2) * eps
    # E(V1 | V2) and \sigma^2(V2) are estimated via GAM
    # eps is univariate and its density is estimated by KDE.
    dat <- data.frame(X2=v1, X1=v2$X1)
    mean.fit <- mgcv::gam(X2 ~ s(X1), data=dat)
    preds.means <- mgcv::predict.gam(mean.fit, newdata=dat)
    
    dat.var <- dat
    dat.var$eps2 <- resid(mean.fit)^2
    
    var.fit <- mgcv::gam(eps2 ~ s(X1),
                         data=dat.var, family=gaussian(link=log))
    preds.vars <- mgcv::predict.gam(var.fit, type="response", newdata=dat.var)
    
    v1.std.tr <- resid(mean.fit)/sqrt(preds.vars)
    dens.est <- ks::kde(x=as.matrix(v1.std.tr), eval.points=v1.std.tr, density=TRUE)
    res.tr <- dens.est$estimate/sqrt(preds.vars)
    
    predict.cond.dens <- function(v1, v2, new.v1, new.v2) {
      # evaluate the estimated conditional density f(v1 | v2) at f(new.v1 | new.v2)
      new.dat <- data.frame(X1=new.v2$X1)
      new.preds.means <- mgcv::predict.gam(mean.fit,
                                           newdata=as.data.frame(new.dat),
                                           type="response")
      new.preds.vars <- mgcv::predict.gam(var.fit,
                                          newdata=as.data.frame(new.dat),
                                          type="response")
      v1.std.te <- (new.v1-new.preds.means)/sqrt(new.preds.vars)
      dens.est <- ks::kde(x=as.matrix(v1.std.te), eval.points=v1.std.te, density=TRUE)
      res <- dens.est$estimate/sqrt(new.preds.vars)
      return(res)
    }
    out <- list(predict.cond.dens=predict.cond.dens)
    return(out)
  }
  
  cate.w1 <- function(tau, w, new.w) {
    SL.lib <- c("SL.mean", "SL.lm", "SL.glm")
    
    model <- SuperLearner::SuperLearner(Y = tau,
                                        X = w,
                                        newX = new.w,
                                        family = gaussian(),
                                        SL.library = SL.lib)

    fit <- function(new.w) predict.SuperLearner(model, newdata = new.w,
                                                verbose = TRUE)$pred[, 1]
    return(list(res = model$SL.predict, model = model, cate.w.form = NULL,
                fit = fit))
  }
  
  cate.not.j <- function(y, x, new.x) {
    SL.lib <- c("SL.mean", "SL.lm", "SL.glm")
    
    fit <- SuperLearner::SuperLearner(Y = y,
                                      X = x,
                                      newX = new.x,
                                      family = gaussian(),
                                      SL.library = SL.lib)
    
    return(fit$SL.predict)
    
  }
  
  reg.basis.not.j <- cate.not.j
  
  fit <- suppressWarnings(cate(
    data=df, x_names=c("X1","X2","X3"), y_name="Y", a_name="A",
    v_names=c("X1", "X2"), v0=v0, learner="dr", pi.x=pi.x, mu1.x=mu1.x, mu0.x=mu0.x, 
    drl.x=drl.x, drl.v=drl.v, nsplits=1, partial_dependence = TRUE,
    cond.dens = list(cond.dens.V1, cond.dens.V2), cate.w = list(cate.w1, cate.w1),
    bw.stage2 = list(c(0.15, 0.5, 1), c(0.15, 0.5, 1))
  ))
  
  pd.theta1 <- fit$pd.res$dr[[1]]$res$theta
  pd.theta2 <- fit$pd.res$dr[[2]]$res$theta
  m1 <- lm(pd.theta1 ~ v0[, 1])
  expect_equal(unname(coef(m1)[2]), 0.5, tolerance = 0.1)
  m2 <- lm(pd.theta2 ~ v0[, 2])
  expect_equal(unname(coef(m2)[2]), -0.8, tolerance = 0.1)
})
