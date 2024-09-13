get_input <- function(data, x_names, y_name, a_name, v_names, v0){

  # a function that return sanitized input according to covariate names
  if ((!is.data.frame(data)&!is.matrix(data))|any(is.na(data))) {
    stop("input data need to be a dataframe/matrix with no missing data")
  }
  # check whether names in x,y,a,v in colnames(data)
  if (!all(x_names %in% colnames(data))|!y_name %in% colnames(data)|!a_name %in% colnames(data)|!all(v_names %in% x_names)){
    stop("variable names do not match")
  }

  a <- data[,a_name]
  y <- data[,y_name]
  x <- data[,x_names, drop=FALSE]
  v <- x[,v_names, drop=FALSE]

  unique.v0 <- list()
  for(i in 1:ncol(v0)) unique.v0[[colnames(v0)[i]]] <- unique(v0[,i])

  res <- list(a=a, y=y, x=x, v=v, unique.v0=unique.v0, v0=v0)
  return(res)
}

# pi.x.lasso <- function(a, x, new.x) {
#   fit <-  glmnet::cv.glmnet(as.matrix(x), a, family = "binomial",
#                             type.measure = "deviance")
#   out <- predict(fit, newx = as.matrix(new.x), s = fit$lambda.min, type = "response")
#   return(as.vector(out))
# }
#
# mu1.x.lasso <- function(y, a, x, new.x) {
#   y1 <- y[a == 1]
#   x1 <- x[a == 1, , drop = FALSE]
#   fit <-  glmnet::cv.glmnet(y = y1, x = as.matrix(x1))
#   out <- predict(fit, newx = as.matrix(new.x), s = fit$lambda.min)
#   return(as.vector(out))
# }
#
# mu0.x.lasso <- function(y, a, x, new.x) {
#   y0 <- y[a == 0]
#   x0 <- x[a == 0, , drop = FALSE]
#   fit <-  glmnet::cv.glmnet(y = y0, x = as.matrix(x0))
#   out <- predict(fit, newx = as.matrix(new.x), s = fit$lambda.min)
#   return(as.vector(out))
# }
#
# # update to also return the model
# drl.lasso <- function(y, x, new.x) {
#   fit <-  glmnet::cv.glmnet(y = y, x = as.matrix(x))
#   out <- predict(fit, newx = as.matrix(new.x), s = fit$lambda.min)
#   # to align with drl function
#   res <- cbind(as.vector(out), NA, NA)
#   return(list(res = res, model = fit))
# }
#
# ul.lasso <- function(y.tr, x.tr, new.x) {
#   fit <-  glmnet::cv.glmnet(y = y.tr, x = x.tr)
#   out <- predict(fit, newx = new.x, s = fit$lambda.min)
#   return(as.vector(out))
# }
#
# mu.x.lasso <- function(y.tr, x.tr, new.x) {
#   fit <-  glmnet::cv.glmnet(y = y.tr, x = x.tr)
#   return(predict(fit, newx = new.x, s = fit$lambda.min))
# }
#
# mu1.x.SL <- function(y.tr, a.tr, x.tr, new.x, sl.lib = NULL){
#   if(!is.null(sl.lib)) {
#     sl.lib <- c("SL.mean", "SL.lm", "SL.gam", "SL.polymars", "SL.rpart")
#   }
#   y1.tr <- y.tr[a.tr == 1]
#   x1.tr <- x.tr[a.tr == 1, , drop = FALSE]
#   fit <- SuperLearner::SuperLearner(Y = y1.tr, X = x1.tr, SL.library = sl.lib,
#                                     newX = new.x)
#   return(fit$SL.predict)
# }
#
# mu0.x.SL <- function(y.tr, a.tr, x.tr, new.x, sl.lib = NULL){
#   if(!is.null(sl.lib)) {
#     sl.lib <- c("SL.mean", "SL.lm", "SL.gam", "SL.polymars", "SL.rpart")
#   }
#   y0.tr <- y.tr[a.tr == 0]
#   x0.tr <- x.tr[a.tr == 0, , drop = FALSE]
#   fit <- SuperLearner::SuperLearner(Y = y0.tr, X = x0.tr, SL.library = sl.lib,
#                                     newX = new.x)
#   return(fit$SL.predict)
# }
#
# mu.x.SL <- function(y.tr, a.tr, x.tr, new.x, sl.lib = NULL){
#   if(!is.null(sl.lib)) {
#     sl.lib <- c("SL.mean", "SL.lm", "SL.gam", "SL.polymars", "SL.rpart")
#   }
#   fit <- SuperLearner::SuperLearner(Y = y.tr, X = x.tr, SL.library = sl.lib,
#                                     newX = new.x)
#   return(fit$SL.predict)
# }
#
# pi.x.SL <- function(a.tr, x.tr, new.x, sl.lib = NULL){
#   if(!is.null(sl.lib)) {
#     sl.lib <- c("SL.mean", "SL.lm", "SL.gam", "SL.polymars", "SL.rpart")
#   }
#   fit <- SuperLearner::SuperLearner(Y = a.tr, X = a.tr, SL.library = sl.lib,
#                                     newX = new.x)
#   return(fit$SL.predict)
# }
#
# mu1.x.lm <- function(y, a, x, new.x) {
#   fit <- lm(y ~ ., data = cbind(y = y[a == 1], x[a == 1, , drop = FALSE]))
#   return(predict(fit, newdata = new.x))
# }
#
# mu0.x.lm <- function(y, a, x, new.x) {
#   fit <- lm(y ~ ., data = cbind(y = y[a == 0], x[a == 0, , drop = FALSE]))
#   return(predict(fit, newdata = new.x))
# }
#
# pi.x.glm <- function(a, x, new.x) {
#   fit <- glm(a ~ ., data = cbind(a = a, x), family = binomial(link = "logit"))
#   return(predict(fit, newdata = new.x, type = "response"))
# }

lm.discrete.v <- function(y, x, new.x) {
  fit <- lm(y ~ x, data = data.frame(y = y, x = x))
  new.dat <- data.frame(x = new.x)
  preds <- predict.lm(fit, newdata = new.dat)
  tt <- delete.response(terms(fit))
  m <- model.frame(tt, new.dat)
  design.mat <- model.matrix(tt, m)
  beta.vcov <- sandwich::vcovHC(fit, type = "HC")
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
  # pi.x.method <- params[["pi.x.method"]]

  # if(is.null(pi.x) & !is.null(pi.x.method)) {
  #   if(pi.x.method == "lasso") pi.x <- pi.x.lasso
  #   else if(pi.x.method == "glm") pi.x <- pi.x.glm
  #   else if(pi.x.method == "SL") pi.x <- pi.x.SL(params[["sl.lib.pi"]])
  #   else stop("Provide valid method for estimating the propensity score.")
  # }

  arg$pi.x <- pi.x

  if(any(learner %in% c("dr", "t"))) {

    mu1.x <- params[["mu1.x"]]
    mu0.x <- params[["mu0.x"]]
    drl.v <- params[["drl.v"]]
    drl.x <- params[["drl.x"]]
    cate.w <- params[["cate.w"]]
    cond.dens <- params[["cond.dens"]]
    cate.not.j <- params[["cate.not.j"]]
    reg.basis.not.j <- params[["reg.basis.not.j"]]
    dfs.rob <- params[["pl.dfs"]]

    # if(is.null(mu1.x)) {
    #   mu1.x.method <- params[["mu1.x.method"]]
    #   if(mu1.x.method == "lasso") {
    #     mu1.x <- mu1.x.lasso
    #   } else if(mu1.x.method == "lm"){
    #     mu1.x <- mu1.x.lm
    #   } else if(mu1.x.method == "SL") {
    #     sl.lib <- params[["SL.lib"]]
    #     mu1.x <- mu1.x.SL(sl.lib)
    #   } else stop("Provide valid method for estimating the E(Y|A=1,X).")
    #
    # }

    # if(is.null(mu0.x)) {
    #   mu0.x.method <- params[["mu1.x.method"]]
    #   if(mu0.x.method == "lasso") {
    #     mu0.x <- mu0.x.lasso
    #   } else if(mu0.x.method == "lm"){
    #     mu0.x <- mu0.x.lm
    #   } else if(mu0.x.method == "SL") {
    #     sl.lib <- params[["SL.lib"]]
    #     mu0.x <- mu0.x.SL
    #   } else stop("Provide valid method for estimating the E(Y|A=0,X).")
    #
    # }

    # if(is.null(drl.v) & any(learner == "dr")) {
    #   drl.v.method <- params[["drl.v.method"]]
    #   if(drl.v.method == "lasso") {
    #     drl.v <- drl.lasso
    #   } else if(drl.v.method == "lm"){
    #     drl.v <- drl.lm
    #   } else if(drl.v.method == "glm"){
    #     drl.v <- drl.glm
    #   } else if(drl.v.method == "gam") {
    #     drl.v <- drl.gam
    #   } else if(drl.v.method == "rf") {
    #     drl.v <- drl.rf
    #   } else stop("Provide valid method for second-stage regression.")
    # }

    # if(is.null(drl.x) & any(learner == "dr")) {
    #   drl.x.method <- params[["drl.x.method"]]
    #   if(is.null(drl.x.method)) {
    #     # by default use lasso
    #     drl.x <- drl.lasso
    #   } else if (drl.x.method == "lm"){
    #     drl.x <- drl.ite.lm
    #   } else if (drl.x.method == "lasso"){
    #     drl.x <- drl.lasso
    #   } else stop("Provide valid method for second-stage regression.")
    # }

    arg$mu1.x <- mu1.x
    arg$mu0.x <- mu0.x
    arg$drl.v <- drl.v
    arg$drl.x <- drl.x
    arg$cate.w <- cate.w
    arg$cond.dens <- cond.dens
    arg$cate.not.j <- cate.not.j
    arg$reg.basis.not.j <- reg.basis.not.j
    arg$dfs.rob <- dfs.rob

  }

  # if(any(learner %in% c("u", "r", "lp-r"))) {
  #
  #   mu.x.method <- params[["mu.x.method"]]
  #   ul.method <- params[["ul.method"]]
  #   ul <- params[["ul"]]
  #
  #   if(is.null(mu.x)) {
  #     if(mu.x.method == "lasso") {
  #       mu.x <- mu.x.lasso
  #     } else if(mu.x.method == "SL") {
  #       mu.x <- mu.x.SL
  #     } else stop("Provide valid method for estimating the E(Y|X).")
  #   }
  #
  #   if(is.null(ul) & any(learner == "u")) {
  #     if(ul.method == "lasso") {
  #       ul <- ul.lasso
  #     } else if(ul.method == "SL") {
  #       ul <- ul.SL
  #     } else stop("Provide valid method for second-stage regression.")
  #   }
  #
  #   arg$mu.x <- mu.x
  #   arg$ul <- ul
  #
  # }

  return(arg)
}
robinson <- function(pseudo, w, v, new.v, s, cate.not.j, reg.basis.not.j, dfs) {
  nsplits <- length(unique(s))
  risk <- rep(NA, length(dfs))
  fits <- vector("list", length=length(dfs))
  # preds <- matrix(NA, ncol = length(dfs), nrow = length(new.v))
  # todo: make the code below faster by only estimating E(p_j(V) | V_{-j}) once.
  for(k in 1:length(dfs)) {

    res.v <- matrix(NA, nrow=length(pseudo), ncol=dfs[k])
    res.y <- rep(NA, length(pseudo))

    for(i in 1:nsplits){
      test.idx <- i==s
      train.idx <- i!=s
      if(all(!train.idx)) train.idx <- test.idx
      w.tr <- w[train.idx, , drop = FALSE]
      w.te <- w[test.idx, , drop = FALSE]
      pseudo.tr <- pseudo[train.idx]
      pseudo.te <- pseudo[test.idx]
      v.tr <- v[train.idx]
      v.te <- v[test.idx]

      p.v.tr <- poly(v.tr, degree=dfs[k], raw=TRUE)
      p.v.te <- poly(v.te, degree=dfs[k], raw=TRUE)

      for(j in 1:dfs[k]){
        res.v[test.idx, j] <- p.v.te[, j] - reg.basis.not.j(y=p.v.tr[, j],
                                                            x=w.tr, new.x=w.te)
      }
      res.y[test.idx] <- pseudo.te - cate.not.j(y=pseudo.tr, x=w.tr, new.x=w.te)
    }

    fit.k <-  lm(res.y ~ -1 + res.v)
    fits[[k]] <- fit.k
    diag.hat.mat <- lm.influence(fit.k, do.coef=FALSE)$hat
    risk[k] <- mean((resid(fit.k)/(1-diag.hat.mat))^2)

  }
  fit.star <- fits[[which.min(risk)]]
  risk.dat <- data.frame(dfs = dfs, risk = risk)
  design.mat <- poly(new.v, degree=dfs[which.min(risk)], raw=TRUE)
  preds <- design.mat%*%coef(fit.star)
  beta.vcov <- sandwich::vcovHC(fit.star)
  sigma2hat <- diag(design.mat%*%beta.vcov%*%t(design.mat))
  ci.ll <- preds-1.96*sqrt(sigma2hat)
  ci.uu <- preds+1.96*sqrt(sigma2hat)
  res <- data.frame(preds=preds, ci.ll=ci.ll, ci.uu=ci.uu)
  out <- list(res=res, model=fit.star, risk=risk.dat, fits=fits)
  return(out)

}


drl.basis.additive <- function(y, x, new.x, kmin = 1, kmax = 10) {
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
      lm.form <- paste0("~ ", paste0("poly(", colnames(x.cont)[1], ", raw = TRUE, degree = ", n.basis[i, 1], ")"))
      # lm.form <- paste0("~ ", paste0("ns(", colnames(x.cont)[1], ", df = ", n.basis[i, 1], ")"))
      if(ncol(x.cont) > 1) {
        for(k in 2:ncol(x.cont)) {
          lm.form <- c(lm.form, paste0("poly(", colnames(x.cont)[k], ", raw = TRUE, degree = ", n.basis[i, k], ")"))
          # lm.form <- c(lm.form, paste0("ns(", colnames(x.cont)[k], ", df = ", n.basis[i, k], ")"))
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
    fit <- lm(as.formula(paste0("y", lm.form)), data = cbind(data.frame(y = y), x.cont, x.disc))
    # x.mat <- model.matrix(as.formula(lm.form), data = x)
    # hat.mat <- x.mat %*% solve(crossprod(x.mat, x.mat)) %*% t(x.mat)
    diag.hat.mat <- lm.influence(fit, do.coef=FALSE)$hat
    # diag.hat.mat <- diag(hat.mat)
    risk[i] <- mean((resid(fit)/(1-diag.hat.mat))^2)
    models[i] <- lm.form
  }
  risk.dat <- cbind(n.basis, risk)
  colnames(risk.dat) <- c(colnames(x.cont), "loocv.risk")
  best.model <- lm(as.formula(paste0("y", models[which.min(risk)])),
                   data = cbind(data.frame(y=y), x))

  out <- predict(best.model, newdata = as.data.frame(new.x))
  res <- cbind(out, NA, NA)
  return((list(drl.form = models[which.min(risk)], res = res, model =  best.model,
               risk = risk.dat)))

}

draw_VIMP <- function(vimp_df){
  # a helper function to illustrate vimp

  vimp_df <- cbind(as.matrix(row.names(vimp_df), nrow(vimp_df), 1), vimp_df)
  colnames(vimp_df) <- c('variable', 'DR', 'Lower_Bound', 'Upper_Bound')
  vimp_df <- vimp_df[order(vimp_df$DR, decreasing = FALSE),]
  order_2b <- vimp_df$variable

  fig_2b <- ggplot(vimp_df, aes(x = factor(variable, level = order_2b), y = DR))+
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Lower_Bound, ymax = Upper_Bound), width=.1, position=position_dodge(0.2))+
    theme(axis.title.y=element_blank(),
          axis.text.y = element_text(color = "grey20", size = 10),
          axis.text.x = element_text(color = "grey20", size = 10)) +
    ylim(-0.1, 1) +
    coord_flip()
  fig_2b
}
