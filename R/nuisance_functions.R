get_input <- function(data, x_names, y_name, a_name, v_names, v0.long, num_grid = 100){

  # a function that return sanitized input according to covariate names

  if ((!is.data.frame(data)&!is.matrix(data))| any(is.na(data))) {
    stop("input data need to be a dataframe/matrix with no missing data")
  }
  # check whether names in x,y,a,v in colnames(data)
  if (!all(x_names %in% colnames(data))|!y_name %in% colnames(data)|!a_name %in% colnames(data)|!all(v_names %in% x_names)){
    stop("variable names do not match")
  }

  a <- data[,a_name]
  y <- data[,y_name]
  x <- data[,x_names, drop = FALSE]
  v <- x[,v_names, drop = FALSE]

  colnames(v) <- stringr::str_c("v", 1:ncol(v))

  v0 <- list()

  for (i in 1:ncol(v)){

    if(!is.null(v0.long)) {

      v0[[paste0("v",i)]] <- unique(v0.long[,i])

    } else {
      if (length(unique(v[,i])) <= 20 | class(v[,i]) == 'factor'){
        tmp.vals <- sort(unique(v[, i]))
      } else {
        tmp.vals <- seq(quantile(v[, i], 0.05),
                        quantile(v[, i], 0.95), length.out = num_grid)
      }
      v0[[paste0("v", i)]] <- tmp.vals
    }
  }

  v0.long <- expand.grid(v0)

  res <- list(a = a, y = y, x = x, v = v, v0 = v0, v0.long = v0.long)
  return(res)
}

pi.x.lasso <- function(a, x, new.x) {
  fit <-  glmnet::cv.glmnet(as.matrix(x), a, family = "binomial",
                            type.measure = "deviance")
  out <- predict(fit, newx = as.matrix(new.x), s = fit$lambda.min, type = "response")
  return(as.vector(out))
}

mu1.x.lasso <- function(y, a, x, new.x) {
  y1 <- y[a == 1]
  x1 <- x[a == 1, , drop = FALSE]
  fit <-  glmnet::cv.glmnet(y = y1, x = as.matrix(x1))
  out <- predict(fit, newx = as.matrix(new.x), s = fit$lambda.min)
  return(as.vector(out))
}

mu0.x.lasso <- function(y, a, x, new.x) {
  y0 <- y[a == 0]
  x0 <- x[a == 0, , drop = FALSE]
  fit <-  glmnet::cv.glmnet(y = y0, x = as.matrix(x0))
  out <- predict(fit, newx = as.matrix(new.x), s = fit$lambda.min)
  return(as.vector(out))
}

# update to also return the model
drl.lasso <- function(y, x, new.x) {
  fit <-  glmnet::cv.glmnet(y = y, x = as.matrix(x))
  out <- predict(fit, newx = as.matrix(new.x), s = fit$lambda.min)
  # to align with drl function
  res <- cbind(as.vector(out), NA, NA)
  return(list(res = res, model = fit))
}

ul.lasso <- function(y.tr, x.tr, new.x) {
  fit <-  glmnet::cv.glmnet(y = y.tr, x = x.tr)
  out <- predict(fit, newx = new.x, s = fit$lambda.min)
  return(as.vector(out))
}

mu.x.lasso <- function(y.tr, x.tr, new.x) {
  fit <-  glmnet::cv.glmnet(y = y.tr, x = x.tr)
  return(predict(fit, newx = new.x, s = fit$lambda.min))
}

mu1.x.SL <- function(y.tr, a.tr, x.tr, new.x, sl.lib = NULL){
  if(!is.null(sl.lib)) {
    sl.lib <- c("SL.mean", "SL.lm", "SL.gam", "SL.polymars", "SL.rpart")
  }
  y1.tr <- y.tr[a.tr == 1]
  x1.tr <- x.tr[a.tr == 1, , drop = FALSE]
  fit <- SuperLearner::SuperLearner(Y = y1.tr, X = x1.tr, SL.library = sl.lib,
                                    newX = new.x)
  return(fit$SL.predict)
}

mu0.x.SL <- function(y.tr, a.tr, x.tr, new.x, sl.lib = NULL){
  if(!is.null(sl.lib)) {
    sl.lib <- c("SL.mean", "SL.lm", "SL.gam", "SL.polymars", "SL.rpart")
  }
  y0.tr <- y.tr[a.tr == 0]
  x0.tr <- x.tr[a.tr == 0, , drop = FALSE]
  fit <- SuperLearner::SuperLearner(Y = y0.tr, X = x0.tr, SL.library = sl.lib,
                                    newX = new.x)
  return(fit$SL.predict)
}

mu.x.SL <- function(y.tr, a.tr, x.tr, new.x, sl.lib = NULL){
  if(!is.null(sl.lib)) {
    sl.lib <- c("SL.mean", "SL.lm", "SL.gam", "SL.polymars", "SL.rpart")
  }
  fit <- SuperLearner::SuperLearner(Y = y.tr, X = x.tr, SL.library = sl.lib,
                                    newX = new.x)
  return(fit$SL.predict)
}

pi.x.SL <- function(a.tr, x.tr, new.x, sl.lib = NULL){
  if(!is.null(sl.lib)) {
    sl.lib <- c("SL.mean", "SL.lm", "SL.gam", "SL.polymars", "SL.rpart")
  }
  fit <- SuperLearner::SuperLearner(Y = a.tr, X = a.tr, SL.library = sl.lib,
                                    newX = new.x)
  return(fit$SL.predict)
}

mu1.x.lm <- function(y, a, x, new.x) {
  fit <- lm(y ~ ., data = cbind(y = y[a == 1], x[a == 1, , drop = FALSE]))
  return(predict(fit, newdata = new.x))
}

mu0.x.lm <- function(y, a, x, new.x) {
  fit <- lm(y ~ ., data = cbind(y = y[a == 0], x[a == 0, , drop = FALSE]))
  return(predict(fit, newdata = new.x))
}

pi.x.glm <- function(a, x, new.x) {
  fit <- glm(a ~ ., data = cbind(a = a, x), family = binomial(link = "logit"))
  return(predict(fit, newdata = new.x, type = "response"))
}

drl.lm <- function(y, x, new.x) {
  dataset <- as.data.frame(cbind(y = y, x))

  # caution ncol(x)
  if (ncol(x) == 1){
    fit <- lm(y ~ poly(v1, 3, raw = TRUE),
              # data = as.data.frame(cbind(y = y, x)))
              data = dataset)
  } else {
    fit <- lm(y ~ poly(v1, 3, raw = TRUE) + poly(v2, 3, raw = TRUE),
              # data = as.data.frame(cbind(y = y, x)))
              data = dataset)
  }
  out <- predict(fit, newdata = as.data.frame(new.x))
  res <- cbind(out, NA, NA)
  return((list(res = res, model = fit)))
}

drl.ite.lm <- function(y, x, new.x) {

  fit <- lm(y ~ poly(x[, 1], 3, raw = TRUE) + poly(x[, 2], 3, raw = TRUE),
            data = cbind(y, x))
  return(predict(fit, newdata = new.x))

}

drl.gam <- function(y, x, new.x) {
  x1 <- colnames(x)[1]
  x2 <- colnames(x)[2]
  data <- cbind(y, x)

  fit <- mgcv::gam(as.formula(paste('y ~ s(', x1, ')+', x2)), data = as.data.frame(data))
  out <- predict(fit, newdata = as.data.frame(new.x))
  res <- cbind(out, NA, NA)
  return((list(res = res, model = fit)))
}

drl.glm <- function(y, x, new.x) {
  dataset <- as.data.frame(cbind(y = y, x))
  fit <- glm(y ~ ., data = dataset, family = binomial(link = "logit"))
  out <- predict(fit, newdata = as.data.frame(new.x), type = "response")
  res <- cbind(out, NA, NA)
  return((list(res = res, model = fit)))
}

drl.rf <- function(y, x, new.x) {
  fit <- randomForest::randomForest(y~., data = cbind(y, x), ntree = 50)
  out <- predict(fit, as.matrix(new.x))
  res <- cbind(out, NA, NA)
  return((list(res = res, model = fit)))
}

# drl.basis <- function(y, x, new.x, kmin = 3, kmax = 10) {
#   require(splines)
#   x <- as.data.frame(x)
#   n.vals <- apply(x, 2, function(u) length(unique(u)))
#   var.type <- unlist(lapply(x, function(u) paste0(class(u), collapse = " ")))
#   factor.boolean <- n.vals <= 10 | var.type == "factor" | var.type == "ordered factor"
#   x.cont <- x[, which(!factor.boolean), drop = FALSE]
#   x.disc <- x[, which(factor.boolean), drop = FALSE]
#
#   n.basis <- expand.grid(rep(list(kmin:kmax), ncol(x.cont)))
#   risk <- models <- rep(NA, nrow(n.basis))
#   for(i in 1:nrow(n.basis)){
#     if(ncol(x.cont) > 0) {
#       lm.form <- paste0("~ ", paste0("bs(", colnames(x.cont)[1], ", df = ", n.basis[i, 1], ")"))
#       if(ncol(x.cont) > 1) {
#         for(k in 2:ncol(x.cont)) {
#           lm.form <- c(lm.form, paste0("bs(", colnames(x.cont)[k], ", df = ", n.basis[i, k], ")"))
#         }
#       }
#     }
#     if(ncol(x.disc) > 0) {
#       for(k in 1:ncol(x.disc)) {
#         if(ncol(x.cont) == 0 & k == 1) {
#           lm.form <- paste0("~ as.factor(", colnames(x.disc)[k], ")")
#         } else {
#           lm.form <- c(lm.form, paste0("as.factor(", colnames(x.disc)[k], ")"))
#         }
#       }
#     }
#
#     lm.form <- paste0(lm.form, collapse = "*")
#     fit <- lm(as.formula(paste0("y", lm.form)), data = cbind(data.frame(y = y), x))
#     # x.mat <- model.matrix(as.formula(lm.form), data = x)
#     # hat.mat <- x.mat %*% solve(crossprod(x.mat, x.mat)) %*% t(x.mat)
#     diag.hat.mat <- lm.influence(fit, do.coef = FALSE)$hat
#     risk[i] <- mean((resid(fit) / (1 - diag.hat.mat))^2)
#     models[i] <- lm.form
#   }
#
#   best.model <- lm(as.formula(paste0("y", models[which.min(risk)])),
#                    data = cbind(data.frame(y = y), x))
#
#   out <- predict(best.model, newdata = as.data.frame(new.x))
#   res <- cbind(out, NA, NA)
#   return((list(drl.form = models[which.min(risk)], res = res, model =  best.model)))
#
# }


lm.discrete.v <- function(y, x, new.x) {
  # univariate
  fit <- lm(y ~ x, data = data.frame(y = y, x = x))
  new.dat <- data.frame(x = new.x)
  preds <- predict.lm(fit, newdata = new.dat)
  tt <- delete.response(terms(fit))
  m <- model.frame(tt, new.dat)
  design.mat <- model.matrix(tt, m)
  beta.vcov <- sandwich::vcovHC(fit)
  sigma2hat <- diag(design.mat %*% beta.vcov %*% t(design.mat))

  ci.l <- preds - 1.96 * sqrt(sigma2hat)
  ci.u <- preds + 1.96 * sqrt(sigma2hat)

  out <- data.frame(
    eval.pts=new.x,
    theta=preds,
    ci.ll.pts=ci.l,
    ci.ul.pts=ci.u,
    ci.ul.unif= NA,
    ci.ll.unif= NA
  )

  return(out)

}

.parse.cate <- function(learner, ...) {

  params <- list(...)
  arg <- list()

  pi.x <- params[["pi.x"]]
  pi.x.method <- params[["pi.x.method"]]

  if(is.null(pi.x) & !is.null(pi.x.method)) {
    if(pi.x.method == "lasso") pi.x <- pi.x.lasso
    else if(pi.x.method == "glm") pi.x <- pi.x.glm
    else if(pi.x.method == "SL") pi.x <- pi.x.SL(params[["sl.lib.pi"]])
    else stop("Provide valid method for estimating the propensity score.")
  }

  arg$pi.x <- pi.x

  if(any(learner %in% c("dr", "t"))) {

    mu1.x <- params[["mu1.x"]]
    mu0.x <- params[["mu0.x"]]
    drl.v <- params[["drl.v"]]
    drl.x <- params[["drl.x"]]
    cate.w <- params[["cate.w"]]
    cond.dens <- params[["cond.dens"]]

    if(is.null(mu1.x)) {
      mu1.x.method <- params[["mu1.x.method"]]
      if(mu1.x.method == "lasso") {
        mu1.x <- mu1.x.lasso
      } else if(mu1.x.method == "lm"){
        mu1.x <- mu1.x.lm
      } else if(mu1.x.method == "SL") {
        sl.lib <- params[["SL.lib"]]
        mu1.x <- mu1.x.SL(sl.lib)
      } else stop("Provide valid method for estimating the E(Y|A=1,X).")

    }

    if(is.null(mu0.x)) {
      mu0.x.method <- params[["mu1.x.method"]]
      if(mu0.x.method == "lasso") {
        mu0.x <- mu0.x.lasso
      } else if(mu0.x.method == "lm"){
        mu0.x <- mu0.x.lm
      } else if(mu0.x.method == "SL") {
        sl.lib <- params[["SL.lib"]]
        mu0.x <- mu0.x.SL
      } else stop("Provide valid method for estimating the E(Y|A=0,X).")

    }

    if(is.null(drl.v) & any(learner == "dr")) {
      drl.v.method <- params[["drl.v.method"]]
      if(drl.v.method == "lasso") {
        drl.v <- drl.lasso
      } else if(drl.v.method == "lm"){
        drl.v <- drl.lm
      } else if(drl.v.method == "glm"){
        drl.v <- drl.glm
      } else if(drl.v.method == "gam") {
        drl.v <- drl.gam
      } else if(drl.v.method == "rf") {
        drl.v <- drl.rf
      } else stop("Provide valid method for second-stage regression.")
    }

    if(is.null(drl.x) & any(learner == "dr")) {
      drl.x.method <- params[["drl.x.method"]]
      if(is.null(drl.x.method)) {
        # by default use lasso
        drl.x <- drl.lasso
      } else if (drl.x.method == "lm"){
        drl.x <- drl.ite.lm
      } else if (drl.x.method == "lasso"){
        drl.x <- drl.lasso
      } else stop("Provide valid method for second-stage regression.")
    }

    arg$mu1.x <- mu1.x
    arg$mu0.x <- mu0.x
    arg$drl.v <- drl.v
    arg$drl.x <- drl.x
    arg$cate.w <- cate.w
    arg$cond.dens <- cond.dens

  }

  if(any(learner %in% c("u", "r", "lp-r"))) {

    mu.x.method <- params[["mu.x.method"]]
    ul.method <- params[["ul.method"]]
    ul <- params[["ul"]]

    if(is.null(mu.x)) {
      if(mu.x.method == "lasso") {
        mu.x <- mu.x.lasso
      } else if(mu.x.method == "SL") {
        mu.x <- mu.x.SL
      } else stop("Provide valid method for estimating the E(Y|X).")
    }

    if(is.null(ul) & any(learner == "u")) {
      if(ul.method == "lasso") {
        ul <- ul.lasso
      } else if(ul.method == "SL") {
        ul <- ul.SL
      } else stop("Provide valid method for second-stage regression.")
    }

    arg$mu.x <- mu.x
    arg$ul <- ul

  }

  return(arg)
}


drl.basis.additive <- function(y, x, new.x, kmin = 3, kmax = 20) {
  require(splines)
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
      lm.form <- paste0("~ ", paste0("bs(", colnames(x.cont)[1], ", df = ", n.basis[i, 1], ")"))
      if(ncol(x.cont) > 1) {
        for(k in 2:ncol(x.cont)) {
          lm.form <- c(lm.form, paste0("bs(", colnames(x.cont)[k], ", df = ", n.basis[i, k], ")"))
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
    risk[i] <- mean((resid(fit) / (1 - diag.hat.mat))^2)
    models[i] <- lm.form
  }

  best.model <- lm(as.formula(paste0("y", models[which.min(risk)])),
                   data = cbind(data.frame(y = y), x))

  out <- predict(best.model, newdata = as.data.frame(new.x))
  res <- cbind(out, NA, NA)
  return((list(drl.form = models[which.min(risk)], res = res, model =  best.model)))

}

