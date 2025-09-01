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

lm.discrete.v <- function(y, x, new.x) {
  fit <- lm(y ~ x, data = data.frame(y = y, x = x))
  new.dat <- data.frame(x = new.x)
  preds <- predict.lm(fit, newdata = new.dat)
  tt <- delete.response(terms(fit))
  m <- model.frame(tt, new.dat)
  design.mat <- model.matrix(tt, m)
  beta.vcov <- sandwich::vcovHC(fit, type = "HC")
  sigma2hat <- diag(design.mat %*% beta.vcov %*% t(design.mat))

  ci.l <- preds-1.96*sqrt(sigma2hat)
  ci.u <- preds+1.96*sqrt(sigma2hat)

  out <- data.frame(
    eval.pts=new.x,
    theta=preds,
    ci.ll.pts=ci.l,
    ci.ul.pts=ci.u,
    ci.ul.unif=NA,
    ci.ll.unif=NA
    )
  return(out)
}

.parse.cate <- function(learner, ...) {

  params <- list(...)
  arg <- list()

  pi.x <- params[["pi.x"]]
  arg$pi.x <- pi.x

  if(any(learner %in% c("dr"))) {

    mu1.x <- params[["mu1.x"]]
    mu0.x <- params[["mu0.x"]]
    drl.v <- params[["drl.v"]]
    drl.x <- params[["drl.x"]]
    cate.w <- params[["cate.w"]]
    cond.dens <- params[["cond.dens"]]
    cate.not.j <- params[["cate.not.j"]]
    reg.basis.not.j <- params[["reg.basis.not.j"]]
    dfs.rob <- params[["pl.dfs"]]

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
  return(arg)
}
robinson <- function(pseudo, w, v, new.v, s, cate.not.j, reg.basis.not.j, dfs) {
  # Estimate \tau(V) = \rho(V_j)^T \beta + m(V_{-j})
  # using Robinson's tranformation:
  # lm(\tau(V) - \E(\tau(V)| V_{-j}) ~ -1 + \rho(V_j) - \E(\rho(V_j) | V_{-j}))
  nsplits <- length(unique(s))
  risk <- rep(NA, length(dfs))
  fits <- vector("list", length=length(dfs))
  # todo: make the code below faster by only estimating E(\rho(V) | V_{-j}) once.
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
  risk.dat <- data.frame(dfs=dfs, risk=risk)
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


#' drl.basis.additive
#' This function fits a low dimensional additive model
#' @param y a numeric vector of outcomes
#' @param x a matrix or data frame of covariates' values
#' @param new.x a matrix or data frame of evaluation points
#' @param kmin minimum number of basis terms (for each covariate) to try in LOOCV step
#' @param kmax maximum number of basis terms (for each covariate) to try in LOOCV step.
#' A total of kmin*kmax model will be evaluated by LOOCV.
#' @return All the fits as well as the best model
#' @export
drl.basis.additive <- function(y, x, new.x, kmin=3, kmax=10) {
  bsc <- function(x, ..., center = TRUE) {
    B <- splines::bs(x, ...)
    if (center) B <- sweep(B, 2, colMeans(B), "-")
    B
  }
  x <- as.data.frame(x)
  n.vals <- apply(x, 2, function(u) length(unique(u)))
  var.type <- unlist(lapply(x, function(u) paste0(class(u), collapse=" ")))
  factor.boolean <- (n.vals <= 10) | (var.type %in% c("factor", "ordered factor"))
  x.cont <- x[, which(!factor.boolean), drop=FALSE]
  x.disc <- x[, which(factor.boolean), drop=FALSE]

  n.basis <- expand.grid(rep(list(kmin:kmax), ncol(x.cont)))
  if(ncol(x.cont)==0) n.basis <- expand.grid(rep(list(1), ncol(x.disc)))
  risk <- models <- rep(NA, nrow(n.basis))
  fits <- vector("list", length=max(nrow(n.basis), 1))
  for(i in 1:nrow(n.basis)){
    if(ncol(x.cont) > 0) {
      # lm.form <- paste0("~ ", paste0("poly(", colnames(x.cont)[1], ", raw = TRUE, degree = ", n.basis[i, 1], ")"))
      lm.form <- paste0("~ ", paste0("bsc(", colnames(x.cont)[1], ", df = ", n.basis[i, 1], ")"))
      if(ncol(x.cont) > 1) {
        for(k in 2:ncol(x.cont)) {
          # lm.form <- c(lm.form, paste0("poly(", colnames(x.cont)[k], ", raw = TRUE, degree = ", n.basis[i, k], ")"))
          lm.form <- c(lm.form, paste0("bsc(", colnames(x.cont)[k], ", df = ", n.basis[i, k], ")"))
          }
      }
    }
    if(ncol(x.disc) > 0) {
      for(k in 1:ncol(x.disc)) {
        if(ncol(x.cont)==0 & k==1) {
          lm.form <- paste0("~ ", colnames(x.disc)[k])
        } else {
          lm.form <- c(lm.form, colnames(x.disc)[k])
        }
      }
    }
    lm.form <- paste0(lm.form, collapse = " + ")
    fits[[i]] <- lm(as.formula(paste0("y", lm.form)),
                    data=cbind(data.frame(y=y), x.cont, x.disc))
    risk[i] <- mean((resid(fits[[i]])/(1-hatvalues(fits[[i]])))^2)
    models[i] <- lm.form
  }
    risk.dat <- cbind(n.basis, risk)
    if(ncol(x.cont)>0) colnames(risk.dat) <- c(colnames(x.cont), "loocv.risk")
    if(ncol(x.cont)==0) colnames(risk.dat) <- c(colnames(x.disc), "loocv.risk")
    # other choices are possible, always plot the estimates risks!
    best.model <- lm(as.formula(paste0("y", models[which.min(risk)])),
                     data=cbind(data.frame(y=y), x))

  out <- predict(best.model, newdata=as.data.frame(new.x))
  res <- cbind(out, NA, NA)
  return((list(drl.form=formula(best.model), res=res, model=best.model,
               risk=risk.dat, fits=fits)))
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

#' get.smooth.fit.gam
#' This function isolates the invdividual smooth components in low dimensional gam fit
#' @param fit output from a low dimensional lm fit where the GAM is stored
#' @param eval.pts the evaluation points for the additive component of interst
#' @param eff.modif.name the name of the effect modifier
#' @param v the matrix of effect modifiers values
#' (the original data subsetted to the effect modifiers considered enetering the GAM)
#' @return matrix of results: estimates and pointwise CIs.
#' @export
get.smooth.fit.gam <- function(fit, eval.pts, eff.modif.name, v) {
  new.dat.additive <- as.data.frame(matrix(0, nrow=length(eval.pts),
                                           ncol=ncol(v),
                                           dimnames=list(NULL, colnames(v))))
  j <- which(colnames(v)==eff.modif.name)
  for(l in 1:ncol(v)) {
    if(l==j) {
      new.dat.additive[, l] <- eval.pts
    } else {
      if(is.factor(v[, l])) {
        new.dat.additive[, l] <- factor(levels(v[, l])[1], levels=levels(v[, l]))
      } else {
        new.dat.additive[, l] <- min(v[, l])
      }
    }
  }
  form <- formula(fit)
  mm <- model.matrix(fit)
  coefs <- coef(fit)

  bs_term_for_vj <- grep(eff.modif.name, colnames(mm), value=TRUE)
  coefs.names.vj <- grep(eff.modif.name, names(coefs), value=TRUE)

  coefs.vj <- coefs[coefs.names.vj]
  tt <- delete.response(terms(fit))
  new.design.mat <- as.matrix(model.matrix(tt, new.dat.additive)[, bs_term_for_vj])

  if(!is.factor(v[, j])) {
    design.mat <- as.matrix(model.matrix(tt, v)[, bs_term_for_vj])
    mean.point <- apply(design.mat, 2, mean)
    new.design.mat <- sweep(new.design.mat, 2, mean.point, FUN = "-")
  }

  preds.j.additive <- new.design.mat %*% coefs.vj
  beta.vcov <- sandwich::vcovHC(fit, type="HC")[coefs.names.vj, coefs.names.vj]
  sigma2hat <- diag(new.design.mat %*% beta.vcov %*% t(new.design.mat))
  ci.l <- preds.j.additive-1.96*sqrt(sigma2hat)
  ci.u <- preds.j.additive+1.96*sqrt(sigma2hat)
  return(data.frame(eval.pts=eval.pts, theta=preds.j.additive, se=sigma2hat,
                    ci.ll.pts=ci.l, ci.ul.pts=ci.u, ci.ll.unif=NA, ci.ul.unif=NA))
}
