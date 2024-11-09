context("second-stage LOOCV")

test_that("expected behavior of second stage bw selector", {
  for(i in 1:2) {

    n <- 200
    x <- rnorm(n, sd = 0.1)
    y <- cos(2*pi*x) + rnorm(n)
    h <- b <- 0.95

    loo.preds <- loo.preds.db <- rep(NA, n)
    for(j in 1:n) {
      loo.preds[j] <- .lprobust(x=x[-j], y=y[-j], h=h, b=b, eval.pt=x[j],
                                debias=FALSE, kernel.type="epa")[, "theta.hat"]
      loo.preds.db[j] <- .lprobust(x=x[-j], y=y[-j], h=h, b=b, eval.pt=x[j],
                                   debias=TRUE, kernel.type="epa")[, "theta.hat"]
    }

    loocv.risk <- mean((y-loo.preds)^2)
    loocv.risk.db <- mean((y-loo.preds.db)^2)

    res <- debiased_inference(A=x, pseudo.out=y,
                              bw.seq=h, eval.pts=x,
                              kernel.type="epa",
                              debias=FALSE,
                              bandwidth.method="LOOCV")
    res.db <- debiased_inference(A=x, pseudo.out=y,
                                 bw.seq=h, eval.pts=x,
                                 kernel.type="epa",
                                 debias=TRUE,
                                 bandwidth.method="LOOCV")

   expect_true(abs(loocv.risk - res$risk$loocv.risk) < 1e-10 ||
               is.infinite(loocv.risk))
   expect_true(abs(loocv.risk.db - res.db$risk$loocv.risk) < 1e-10 ||
                 is.infinite(loocv.risk.db))
  }
})

test_that("expected behavior of LOOCV additive basis", {

  drl.basis.additive2 <- function(y, x, new.x, kmin = 1, kmax = 10) {
    # require(splines)
    x <- as.data.frame(x)
    n.vals <- apply(x, 2, function(u) length(unique(u)))
    var.type <- unlist(lapply(x, function(u) paste0(class(u), collapse = " ")))
    factor.boolean <- (n.vals <= 10) | (var.type %in% c("factor", "ordered factor"))
    x.cont <- x[, which(!factor.boolean), drop = FALSE]
    x.disc <- x[, which(factor.boolean), drop = FALSE]

    n.basis <- expand.grid(rep(list(kmin:kmax), ncol(x.cont)))
    risk <- models <- rep(NA, nrow(n.basis))
    for(i in 1:nrow(n.basis)){
      if(ncol(x.cont) > 0) {
        # lm.form <- paste0("~ ", paste0("ns(", colnames(x.cont)[1], ", degree = ", n.basis[i, 1], ")"))
        lm.form <- paste0("~ ", paste0("poly(", colnames(x.cont)[1], ", degree = ", n.basis[i, 1], ")"))
        if(ncol(x.cont) > 1) {
          for(k in 2:ncol(x.cont)) {
            # lm.form <- c(lm.form, paste0("ns(", colnames(x.cont)[k], ", degree = ", n.basis[i, k], ")"))
            lm.form <- c(lm.form, paste0("poly(", colnames(x.cont)[k], ", degree = ", n.basis[i, k], ")"))
          }
        }
      }
      if(ncol(x.disc) > 0) {
        for(k in 1:ncol(x.disc)) {
          if(ncol(x.cont) == 0 & k == 1) {
            lm.form <- paste0("~ ", colnames(x.disc)[k])
          } else {
            lm.form <- c(lm.form, colnames(x.disc)[k])
          }
        }
      }

      lm.form <- paste0(lm.form, collapse = " + ")
      fit <- lm(as.formula(paste0("y", lm.form)), data = cbind(data.frame(y = y), x))
      # x.mat <- model.matrix(as.formula(lm.form), data = x)
      # hat.mat <- x.mat %*% solve(crossprod(x.mat, x.mat)) %*% t(x.mat)
      diag.hat.mat <- lm.influence(fit, do.coef = FALSE)$hat
      # diag.hat.mat <- diag(hat.mat)
      risk[i] <- mean((resid(fit) / (1 - diag.hat.mat))^2)
      models[i] <- lm.form
    }

    best.model <- lm(as.formula(paste0("y", models[which.min(risk)])),
                     data = cbind(data.frame(y = y), x))

    out <- predict(best.model, newdata = as.data.frame(new.x))
    res <- cbind(out, NA, NA)
    return((list(drl.form = models[which.min(risk)], res = res, model =  best.model,
                 risk = risk)))

  }

  for(i in 1:2) {

    n <- 250
    x1 <- rnorm(n)
    x2 <- runif(n)
    x3 <- rbinom(n, 1, 0.5)
    x <- data.frame(x1 = x1, x2 = x2, x3 = x3)
    y <- rnorm(n, x1^2 + x2  + x3, 0.5)
    dat <- cbind(y = y, x)
    new.x <- x
    fit <- drl.basis.additive(y, x, new.x, kmin = 1, kmax = 5)
    fit3 <- drl.basis.additive2(y, x, new.x, kmin = 1, kmax = 2)
    fit2 <- lm(as.formula(paste0("y ", fit$drl.form)),
               data = dat)
    resids <- resids3 <- rep(NA, n)
    x.mat <- model.matrix(as.formula(fit$drl.form), data = x)
    x.mat3 <- model.matrix(as.formula(fit3$drl.form), data = x)
    for(j in 1:n) {

      fit.not.j3 <- lm(as.formula(paste0("y ", fit3$drl.form)),
                      data = dat[-j, ])
      resids3[j] <- y[j] - predict(fit.not.j3,
                                  newdata = x[j, , drop = FALSE])
      x.mat.j <- x.mat[-j, ]
      y.j <- y[-j]
      betahat <- solve(crossprod(x.mat.j, x.mat.j)) %*% t(x.mat.j) %*% y.j
      resids[j] <- y[j] - x.mat[j, ] %*% betahat

    }
    expect_true(max(abs(fit$res[, 1] - fitted(fit2))) < 1e-10)
    expect_true(abs(min(fit$risk) - mean(resids^2)) < 1e-10)
    expect_true(abs(min(fit3$risk) - mean(resids3^2)) < 1e-10)
  }
})
