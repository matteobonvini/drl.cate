context("gaps in the support of eff modif")

test_that("expected NAs when there are gaps", {
  
  n <- 500
  y <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  x <- as.data.frame(matrix(rnorm(n * 3), ncol = 3))
  x[1:2, 1] <- rep(min(x[, 1]) - 10, 2)
  x$x4 <- factor(sample(2:3, n, replace = TRUE), levels = c(2, 3))
  x$x5 <- factor(sample(4:6, n, replace = TRUE), levels = c(4,5,6))
  colnames(x) <- paste0("x", 1:5)
  
  pi.x <- function(a, x, new.x) return(rep(0.5, nrow(new.x)))
  mu1.x <- function(y, a, x, new.x) {
    mm <- model.matrix(as.formula("~."), new.x)
    beta <- c(0, rep(1, ncol(mm) -1))
    return(mm %*% beta)
  }
  mu0.x <- function(y, a, x, new.x) {
    mm <- model.matrix(as.formula("~."), new.x)
    beta <- c(0, rep(0.5, ncol(mm) - 1))
    return(mm %*% beta)
  }
  drl.v <- function(y, x, new.x) {
    fit <- lm(y ~ ., data = cbind(data.frame(y = y), x))
    res <- cbind(predict.lm(fit, newdata = as.data.frame(new.x)), NA, NA)
    return(list(drl.form = "~ .", res = res, model = fit))
  }
  
  drl.x <- drl.v
  
  cond.dens <- function(v1, v2) {
    
    if(is.factor(v1)) {
      probs.fit <- nnet::multinom("y ~ .", data = cbind(data.frame(y = v1), v2))
      n.lvl <- length(levels(v1))
      predict.cond.dens <- function(v1, v2, new.v1, new.v2) {
        preds <- predict(probs.fit, newdata = cbind(data.frame(y = new.v1), new.v2), 
                         type = "probs")
        preds.t <- t(preds)
        
        if(n.lvl == 2) preds.t <- rbind(1-preds, preds.t)
        
        tmp <- matrix(rep(new.v1, each = n.lvl), byrow = FALSE, ncol = length(new.v1), nrow = n.lvl)
        tmp2 <- matrix(levels(v1), byrow = FALSE, ncol = length(new.v1), nrow = n.lvl)
        return(preds.t[which(tmp == tmp2, arr.ind = TRUE)])
      }
      out <- list(predict.cond.dens = predict.cond.dens)
    } else {
      
      mean.fit <- drl.basis(y = v1, x = v2, new.x = v2, kmax = 10)
      preds.means <- mean.fit$res[, 1]
      # new.preds.means <- mean.fit$res[(nrow(v2)+1:length(preds.means)), 1]
      
      var.fit <- drl.basis(y = (v1 - preds.means)^2, x = v2, new.x = v2, kmax = 10)
      preds.vars <- var.fit$res[, 1]
      preds.vars[preds.vars < 0.05] <- 0.05
      # new.preds.vars <- preds.vars.all[(nrow(v2)+1:length(preds.vars)), 1]
      
      v1.std.tr <- (v1 - preds.means) / sqrt(preds.vars)
      # v1.std.te <- (new.v1 - new.preds.means) / sqrt(new.preds.vars)
      
      res.tr <- approx(density(v1.std.tr)$x, density(v1.std.tr)$y,
                       xout = v1.std.tr)$y / sqrt(preds.vars)
      
      predict.cond.dens <- function(v1, v2, new.v1, new.v2) {
        
        new.preds.means <- predict.lm(mean.fit$model, 
                                      newdata = as.data.frame(new.v2))
        new.preds.vars <- predict.lm(var.fit$model, 
                                     newdata = as.data.frame(new.v2))
        new.preds.vars[new.preds.vars < 0.05] <- 0.05
        v1.std.te <- (new.v1 - new.preds.means) / sqrt(new.preds.vars)
        
        res <- approx(density(v1.std.tr)$x, density(v1.std.tr)$y, xout = v1.std.te,
                      rule = 2)$y / 
          sqrt(new.preds.vars)
        return(res)
        
      }
    }
    out <- list(predict.cond.dens = predict.cond.dens)
    return(out)
  }
  
  cate.w <- function(tau, w, new.w) {
    
    model <- lm(y ~ ., data = cbind(data.frame(y = tau), w))
    
    res <- predict(model, newdata = as.data.frame(new.w))
    fit <- function(new.w) predict.lm(model, newdata = as.data.frame(new.w))
    return(list(res = res, model = model, cate.w.form = "~ .", fit = fit))
  }
  
  data <- cbind(data.frame(y = y, a = a), x)
  
  v0.long <- expand.grid(c(min(x$x1) + 5, 
                           quantile(x$x1, c(0.30, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65))),
                         seq(quantile(x$x3, 0.05),
                             quantile(x$x3, 0.95), length.out = 10),
                         levels(x$x5))
  colnames(v0.long) <- c("x1", "x3", "x5")
  cate.fit <- cate(data_frame = data, learner = "dr",
                   x_names = paste0("x", 1:5),
                   y_name = "y",
                   a_name = "a",
                   v_names = c("x1", "x3", "x5"),
                   univariate_reg = TRUE,
                   partial_dependence = TRUE,
                   additive_approx = TRUE,
                   nsplits = 2,
                   v0.long = v0.long,
                   mu1.x = mu1.x,
                   mu0.x = mu0.x,
                   pi.x = pi.x,
                   drl.v = drl.v,
                   drl.x = drl.x,
                   cond.dens = cond.dens,
                   cate.w = cate.w,
                   bw.stage2 = NULL)
  
  cate.fit$univariate_res$dr[[1]]$res
  expect_true(is.na(cate.fit$univariate_res$dr[[1]]$res[1, "theta"]))
  expect_true(is.na(cate.fit$univariate_res$dr[[1]]$res[1, "theta.debias"]))
  expect_true(is.na(cate.fit$univariate_res$dr[[1]]$res[1, "ci.ul.pts"]))
  expect_true(!is.na(cate.fit$univariate_res$dr[[1]]$res[2, "theta"]))
})

