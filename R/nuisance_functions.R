get_input <- function(data, x_names, y_name, a_name, v_names){
  
  # a function that return sanitized input according to covariate names
  
  if (!is.data.frame(data)| any(is.na(data))) {
    stop("input data need to be a dataframe with no missing data")
  }
  # check whether names in x,y,a,v in colnames(data)
  if (!all(x_names %in% colnames(data))|!y_name %in% colnames(data)|!a_name %in% colnames(data)|!all(v_names %in% x_names)){
    stop("variable names do not match")
  }
  
  data.est <- data %>% rename(a = a_name, y = y_name) %>%
    dplyr::select(all_of(c("y", "a", x_names)))
  
  a <- pull(data.est, "a")
  y <- pull(data.est, "y")
  x <- data.est %>% dplyr::select(all_of(x_names))
  v <- data.est %>% dplyr::select(all_of(v_names))
  
  # for now v0.long only works for 2 dimensional v
  colnames(v) <- stringr::str_c("v", 1:ncol(v))
  
  v1.vals <- seq(min(v[, 1]), max(v[, 1]), length.out = 100)
  v2.vals <- seq(min(v[, 2]), max(v[, 2]), length.out = 100)
  v0 <- cbind(v1 = v1.vals, v2 = v2.vals)
  v0.long <- expand.grid(v1 = v1.vals, v2 = v2.vals)
  
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

drl.lasso <- function(y, x, new.x) {
  fit <-  glmnet::cv.glmnet(y = y, x = as.matrix(x))
  out <- predict(fit, newx = as.matrix(new.x), s = fit$lambda.min)
  # to align with drl function
  res <- cbind(as.vector(out), NA, NA)
  return(res)
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
  fit <- lm(y ~ poly(v1, 3, raw = TRUE) + poly(v2, 3, raw = TRUE),
            data = as.data.frame(cbind(y = y, x)))
  ret <- cbind(predict(fit, newdata = as.data.frame(new.x)), NA, NA)
  return(ret)
  
}

drl.ite.lm <- function(y, x, new.x) {
  
  fit <- lm(y ~ poly(x[, 1], 3, raw = TRUE) + poly(x[, 2], 3, raw = TRUE),
            data = cbind(y, x))
  return(predict(fit, newdata = new.x))
  
}
