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

.get.pseudo.y <- function(y, a, x, v1, v2, mu0.x, mu1.x, pi.x, cond.dens,
                          cate.w, nsplits = 1) {

  n <- length(y)
  s <- sample(rep(1:nsplits, ceiling(n / nsplits))[1:n])
  mu0hat <- mu1hat <- pihat <- tauhat <- cate.out <- rep(NA, n)
  cate.out.folds <- pd.out.folds <- vector("list", nsplits)

  if(!is.null(v2)) {
    theta.bar <- ghat <- tauhat.w <- rep(NA, n)
    theta.mat <- matrix(NA, ncol = n, nrow = n)
    pd.out <- matrix(NA, ncol = ncol(v1), nrow = n,
                     dimnames = list(NULL, colnames(v1)))
  }
  else pd.out <- theta.bar <- theta.mat <- NULL

  for(k in 1:nsplits) {

    test.idx <- k == s
    train.idx <- k != s

    if(all(!train.idx)) train.idx <- test.idx

    n.te <- sum(test.idx)
    n.tr <- sum(train.idx)

    a.tr <- a[train.idx]
    y.tr <- y[train.idx]
    x.tr <- x[train.idx, , drop = FALSE]

    a.te <- a[test.idx]
    y.te <- y[test.idx]
    x.te <- x[test.idx, , drop = FALSE]

    mu0hat.vals <- mu0.x(y = y.tr, a = a.tr, x = x.tr,
                         new.x = rbind(x.te, x.tr))
    mu0hat[test.idx] <- mu0hat.vals[1:n.te]

    mu1hat.vals <- mu1.x(y = y.tr, a = a.tr, x = x.tr,
                         new.x = rbind(x.te, x.tr))
    mu1hat[test.idx] <- mu1hat.vals[1:n.te]

    pihat[test.idx] <- pi.x(a = a.tr, x = x.tr, new.x = x.te)

    tauhat[test.idx] <- mu1hat[test.idx] - mu0hat[test.idx]

    cate.out[test.idx] <- a.te / pihat[test.idx] * (y.te - mu1hat[test.idx]) -
      (1-a.te) / (1 - pihat[test.idx]) * (y.te - mu0hat[test.idx]) +
      tauhat[test.idx]

    cate.out.folds[[k]] <- cate.out[test.idx]

    if(!is.null(v2)) {

      for(j in 1:ncol(v1)) {
        v1.j <- v1[, j]
        v2.not.v1.j <- v2[, colnames(v2) != colnames(v1)[j], drop = FALSE]
        v1.tr <- v1.j[train.idx]
        v2.tr <- v2[train.idx, , drop = FALSE]
        v2.not.v1.j.tr <- v2.not.v1.j[train.idx, , drop = FALSE]

        v1.te <- v1.j[test.idx]
        v2.te <- v2[test.idx, , drop = FALSE]
        v2.not.v1.j.te <- v2.not.v1.j[test.idx, , drop = FALSE]

        w.long.test <- cbind(v1[test.idx, j, drop = FALSE],
                             v2.not.v1.j.te[rep(1:n.te, n.te), , drop = FALSE])

        cond.dens.vals <- cond.dens(v1 = v1.tr,
                                    v2 = v2.not.v1.j.tr,
                                    new.v1 = c(v1.te, w.long.test[, 1]),
                                    new.v2 = rbind(v2.not.v1.j.te,
                                                   w.long.test[, -1, drop = FALSE]))
        marg.dens <- colMeans(matrix(cond.dens.vals[-c(1:n.te)], ncol = n.te,
                                     nrow = n.te))
        ghat[test.idx] <- marg.dens / cond.dens.vals[1:n.te]

        w.long <- cbind(v1[test.idx, j, drop = FALSE],
                        v2.not.v1.j.te[rep(1:n.te, n), , drop = FALSE])

        cate.preds <- cate.w(tau = mu1hat.vals[-c(1:n.te)] - mu0hat.vals[-c(1:n.te)],
                             w = v2.tr,
                             new.w = rbind(v2.te, w.long.test, w.long))

        tauhat.w[test.idx] <- cate.preds$res[1:n.te]
        theta.mat.test <- matrix(cate.preds$res[(n.te+1):(n.te + n.te^2)],
                                 nrow = n.te, ncol = n.te)
        theta.mat[, test.idx] <- matrix(cate.preds$res[-c(1:(n.te + n.te^2))],
                                        nrow = n, ncol = n.te)
        theta.bar[test.idx] <- apply(theta.mat.test, 2, mean)

        pd.out[test.idx, j] <-
          (cate.out[test.idx] - tauhat.w[test.idx]) * ghat[test.idx] +
          theta.bar[test.idx]

      }
      pd.out.folds[[k]] <- pd.out[test.idx, ]
    }
  }

  return(list(cate.out = cate.out, pd.out = pd.out, theta.mat = theta.mat,
              theta.bar = theta.bar,
              cate.out.folds = cate.out.folds,
              pd.out.folds = pd.out.folds))
}

