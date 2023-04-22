get_input <- function(data, x_names, y_name, a_name, v_names, num_grid = 100){

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

  # for now v0.long only works for 2 dimensional v
  colnames(v) <- stringr::str_c("v", 1:ncol(v))
  v0 <- list()

  for (i in 1:ncol(v)){
    if (length(unique(v[,i])) <= 20|class(v[,i]) == 'factor'){
      tmp.vals <- sort(unique(v[,i]))
    } else {
      tmp.vals <- seq(min(v[, i]), max(v[, i]), length.out = num_grid)
    }

    # v0[[paste0("v", i, ".vals")]] <- tmp.vals
    v0[[paste0("v", i)]] <- tmp.vals
  }
  v0.long <- expand.grid(v0)

  res <- list(a = a, y = y, x = x, v = v, v0 = v0, v0.long= v0.long)
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

drl.basis <- function(y, x, new.x, kmax = 10) {
  x <- as.data.frame(x)
  n.vals <- apply(x, 2, function(u) length(unique(u)))
  x.cont <- x[, which(n.vals > kmax), drop = FALSE]
  x.disc <- x[, which(n.vals <= kmax), drop = FALSE]

   n.basis <- expand.grid(rep(list(1:kmax), ncol(x.cont)))
   risk <- models <- rep(NA, nrow(n.basis))
   for(i in 1:nrow(n.basis)){
     if(ncol(x.cont) > 0) {
     lm.form <- paste0("y ~ ", paste0("poly(", colnames(x.cont)[1], ",", n.basis[i, 1], ")"))
     if(ncol(x.cont) > 1) {
       for(k in 2:ncol(x.cont)) {
         lm.form <- c(lm.form, paste0("poly(", colnames(x.cont)[k], ",", n.basis[i, k], ")"))
       }
     }
     }
     if(ncol(x.disc) > 0) {
       for(k in 1:ncol(x.disc)) {
         lm.form <- c(lm.form, paste0(colnames(x.disc)[k]))
       }
     }

     lm.form <- paste0(lm.form, collapse = "+")
     fit <- lm(as.formula(lm.form), data = cbind(data.frame(y = y), x))
     x.mat <- model.matrix(as.formula(lm.form), data = x)
     hat.mat <- x.mat %*% solve(t(x.mat) %*% x.mat) %*% t(x.mat)
     risk[i] <- mean((resid(fit) / (1 - diag(hat.mat)))^2)
     models[i] <- lm.form
   }

  best.model <- lm(as.formula(models[which.min(risk)]), data = cbind(data.frame(y = y), x))

  out <- predict(best.model, newdata = as.data.frame(new.x))
  res <- cbind(out, NA, NA)
  return((list(drl.form = models[which.min(risk)], res = res, model =  best.model)))

}


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

