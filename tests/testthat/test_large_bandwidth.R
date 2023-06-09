context("large bandwidth")

test_that("expected results in simple linear second-stage model", {

  n <- 500
  y <- rnorm(n)
  a <- rbinom(n, 1, 0.5)
  x <- as.data.frame(matrix(rnorm(n * 3), ncol = 3))
  x$x4 <- factor(sample(2:3, n, replace = TRUE), levels = c(2, 3))
  x$x5 <- factor(sample(4:6, n, replace = TRUE), levels = c(4,5,6))
  colnames(x) <- paste0("x", 1:5)

  pi.x <- function(a, x, new.x) {
    return(list(res = rep(0.5, nrow(new.x))))
  }
  mu1.x <- function(y, a, x, new.x) {
    mm <- model.matrix(as.formula("~."), new.x)
    beta <- c(0, rep(1, ncol(mm) -1))
    return(list(res = mm %*% beta))
  }
  mu0.x <- function(y, a, x, new.x) {
    mm <- model.matrix(as.formula("~."), new.x)
    beta <- c(0, rep(0.5, ncol(mm) - 1))
    return(list(res = mm %*% beta))
  }
  drl.v <- function(y, x, new.x) {
    fit <- lm(y ~ ., data = cbind(data.frame(y = y), x))
    res <- cbind(predict.lm(fit, newdata = as.data.frame(new.x)), NA, NA)
    return(list(drl.form = "~ .", res = res, model = fit))
  }

  drl.x <- drl.v

  # cond.dens <- function(v1, v2) {
  #   probs.fit <- nnet::multinom("y ~ .", data = cbind(data.frame(y = v1), v2))
  #   n.lvl <- length(levels(v1))
  #   fit <- function(v1, v2, new.v1, new.v2) {
  #     preds <- predict(probs.fit, newdata = cbind(data.frame(y = new.v1), new.v2),
  #                      type = "probs")
  #     preds.t <- t(preds)
  #
  #     if(n.lvl == 2) preds.t <- rbind(1-preds, preds.t)
  #
  #     tmp <- matrix(rep(new.v1, each = n.lvl), byrow = FALSE, ncol = length(new.v1), nrow = n.lvl)
  #     tmp2 <- matrix(levels(v1), byrow = FALSE, ncol = length(new.v1), nrow = n.lvl)
  #     return(preds.t[which(tmp == tmp2, arr.ind = TRUE)])
  #   }
  #   out <- list(predict.cond.dens = fit)
  #   return(out)
  # }

  cond.dens <- function(v1, v2) {
    fit <- function(v1, v2, new.v1, new.v2) {
      preds <- rep(NA, length(new.v1))
      for(u in unique(new.v1)) preds[new.v1 == u] <- mean(new.v1 == u)
      return(preds)
    }
    out <- list(predict.cond.dens = fit)
    return(out)
  }

  cate.w <- function(tau, w, new.w) {

    model <- lm(y ~ ., data = cbind(data.frame(y = tau), w))

    res <- predict(model, newdata = as.data.frame(new.w))
    fit <- function(new.w) predict.lm(model, newdata = as.data.frame(new.w))
    return(list(res = res, model = model, cate.w.form = "~ .", fit = fit))
  }

  data <- cbind(data.frame(y = y, a = a), x)

  cate.fit <- cate(data_frame = data, learner = "dr",
                   x_names = paste0("x", 1:5),
                   y_name = "y",
                   a_name = "a",
                   v_names = "x4",
                   univariate_reg = TRUE,
                   partial_dependence = TRUE,
                   additive_approx = TRUE,
                   nsplits = 2,
                   mu1.x = mu1.x,
                   mu0.x = mu0.x,
                   pi.x = pi.x,
                   drl.v = drl.v,
                   drl.x = drl.x,
                   cond.dens = cond.dens,
                   cate.w = cate.w,
                   bw.stage2 = NULL)

  foldid <- cate.fit$foldid
  mu1hat <- mu1.x(y, a, x, x)$res
  mu0hat <- mu0.x(y, a, x, x)$res
  true.pseudo <- a / 0.5 * (y - mu1hat) + mu1hat -
    (1-a) / 0.5 * (y - mu0hat) - mu0hat

  expect_true(max(abs(true.pseudo - cate.fit$pd_res$dr[[1]]$data$pseudo)) < 1e-12)

  expect_true(max(abs(cate.fit$univariate_res$dr[[1]]$data$pseudo - true.pseudo)) < 1e-12)

  res <- cate.fit$est[[1]][, 1]

  true.additive_approx1 <- mean(true.pseudo[x$x4 == 2])
  true.additive_approx2 <- mean(true.pseudo[x$x4 == 3])
  expect_true(max(abs(cate.fit$additive_res$dr[[1]]$theta[1] - true.additive_approx1)) < 1e-12)
  expect_true(max(abs(cate.fit$additive_res$dr[[1]]$theta[2] - true.additive_approx2)) < 1e-12)

  true.res1 <- (mean(true.pseudo[foldid == 1 & x$x4 == 2]) +
                  mean(true.pseudo[foldid == 2 & x$x4 == 2]))/2
  true.res2 <- (mean(true.pseudo[foldid == 1 & x$x4 == 3]) +
                  mean(true.pseudo[foldid == 2 & x$x4 == 3]))/2

  expect_true(abs(res[1] - true.res1) < 1e-12)
  expect_true(abs(res[2] - true.res2) < 1e-12)

  univariate_res <- cate.fit$univariate_res$dr[[1]]$res
  true.univariate_res1 <- mean(true.pseudo[x$x4 == 2])
  true.univariate_res2 <- mean(true.pseudo[x$x4 == 3])
  expect_true(abs(univariate_res$theta[1] - true.univariate_res1) < 1e-12)
  expect_true(abs(univariate_res$theta[2] - true.univariate_res2) < 1e-12)

  bw.stage2 <- list(100, 100, NULL)

  cate.fit2 <- cate(data_frame = data, learner = "dr",
                    x_names = paste0("x", 1:5),
                    y_name = "y",
                    a_name = "a",
                    v_names = c("x1", "x3", "x5"),
                    univariate_reg = TRUE,
                    partial_dependence = FALSE,
                    additive_approx = TRUE,
                    nsplits = 2,
                    mu1.x = mu1.x,
                    mu0.x = mu0.x,
                    pi.x = pi.x,
                    drl.v = drl.v,
                    drl.x = drl.x,
                    cond.dens = NULL,
                    cate.w = NULL,
                    bw.stage2 = bw.stage2)

  expect_true(max(abs(cate.fit2$univariate_res$dr[[1]]$data$pseudo - true.pseudo)) < 1e-12)
  expect_true(max(abs(cate.fit2$univariate_res$dr[[2]]$data$pseudo - true.pseudo)) < 1e-12)
  expect_true(max(abs(cate.fit2$univariate_res$dr[[3]]$data$pseudo - true.pseudo)) < 1e-12)

  true.univariate_res1 <- predict.lm(lm(y ~ x, data = data.frame(y = true.pseudo,
                                                                 x = x$x1)),
                                     newdata = data.frame(x = cate.fit2$univariate_res$dr[[1]]$res$eval.pts))

  # expect_true(max(abs(cate.fit2$univariate_res$dr[[1]]$res$theta - true.univariate_res1)) < 1e-12)

  true.univariate_res2 <- predict.lm(lm(y ~ x, data = data.frame(y = true.pseudo,
                                                                 x = x$x3)),
                                     newdata = data.frame(x = cate.fit2$univariate_res$dr[[2]]$res$eval.pts))

  # expect_true(max(abs(cate.fit2$univariate_res$dr[[2]]$res$theta - true.univariate_res2)) < 1e-12)

  true.v0.long <- expand.grid(seq(quantile(x$x1, 0.05),
                                  quantile(x$x1, 0.95), length.out = 100),
                              seq(quantile(x$x3, 0.05),
                                  quantile(x$x3, 0.95), length.out = 100),
                              levels(x$x5))
  colnames(true.v0.long) <- c("x1", "x3", "x5")
  expect_true(max(abs(true.v0.long$x1 - cate.fit2$v0.long$v1)) < 1e-12)
  expect_true(all(true.v0.long$x5 == cate.fit2$v0.long$v3))
  fit1 <- lm(y ~ . , data = cbind(data.frame(y = true.pseudo[cate.fit2$foldid == 1]),
                                       x[cate.fit2$foldid == 1, c(1, 3, 5), drop = FALSE]))

  fit2 <- lm(y ~ x1 + x3 + x5, data = cbind(data.frame(y = true.pseudo[cate.fit2$foldid == 2]),
                                       x[cate.fit2$foldid == 2, c(1, 3, 5), drop = FALSE]))

  true.res <- (predict(fit1, newdata = true.v0.long) +
                 predict(fit2, newdata = true.v0.long)) / 2

  expect_true(max(abs(cate.fit2$est[[1]][, 1] - true.res)) < 1e-12)

})

