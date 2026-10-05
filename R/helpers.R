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

  if ((!is.data.frame(v0) && !is.matrix(v0)) || nrow(v0)==0L ||
      is.null(colnames(v0)) || anyDuplicated(colnames(v0)) ||
      anyDuplicated(v_names) || !setequal(colnames(v0), v_names) || anyNA(v0)) {
    stop("v0 must be a nonempty matrix/data frame with one uniquely named column per v_names and no missing values")
  }
  # Downstream modifier loops use positions: align names before extracting grids.
  v0 <- v0[, v_names, drop=FALSE]
  unique.v0 <- list()
  for(i in 1:ncol(v0)) unique.v0[[colnames(v0)[i]]] <- unique(v0[,i])

  res <- list(a=a, y=y, x=x, v=v, unique.v0=unique.v0, v0=v0)
  return(res)
}

# Effect-level inference for an unweighted, full-rank second-stage lm.
# B and B0 use the full coefficient layout, with zero columns for coefficients
# outside the modifier's component. QR avoids forming ill-conditioned X'X.
effect_profile <- function(fit, pseudo, B, B0, eval.pts) {
  n <- length(pseudo)
  q <- fit$qr
  p <- length(stats::coef(fit))
  if (n < 2L || stats::nobs(fit) != n || nrow(B) != n ||
      !is.null(fit$weights) || is.null(q) || q$rank != p ||
      any(!is.finite(stats::coef(fit)))) {
    stop("Effect profiles require an unweighted, full-rank lm on all pseudo-outcomes.")
  }
  C <- sweep(B0, 2, colMeans(B), "-")
  centered <- drop(sweep(B, 2, colMeans(B), "-") %*% stats::coef(fit))
  component <- drop(C %*% stats::coef(fit))
  psi <- mean(pseudo)
  piv <- q$pivot[seq_len(p)]
  U <- qr.R(q)[seq_len(p), seq_len(p), drop=FALSE]
  Q <- qr.Q(q)[, seq_len(p), drop=FALSE]
  se <- component.se <- numeric(nrow(C))
  # Bound memory when evaluation points include every observation.
  for (start in seq.int(1L, nrow(C), by=64L)) {
    ii <- start:min(start + 63L, nrow(C))
    projection <- n * Q %*% backsolve(U, t(C[ii, piv, drop=FALSE]),
                                     transpose=TRUE)
    component.if <- projection * stats::resid(fit) - centered
    level.if <- component.if + (pseudo - psi)
    component.se[ii] <- sqrt(apply(component.if, 2, stats::var) / n)
    se[ii] <- sqrt(apply(level.if, 2, stats::var) / n)
  }
  result <- function(theta, se) {
    data.frame(eval.pts=eval.pts, theta=theta, se=se,
               ci.ll.pts=theta-stats::qnorm(.975)*se,
               ci.ul.pts=theta+stats::qnorm(.975)*se,
               ci.ll.unif=NA_real_, ci.ul.unif=NA_real_)
  }
  list(res=result(psi + component, se),
       component.res=result(component, component.se), ate=psi)
}

robinson_effect_profile <- function(fit, pseudo, v, eval.pts) {
  powers <- seq_along(stats::coef(fit))
  effect_profile(fit, pseudo, outer(as.numeric(v), powers, `^`),
                 outer(as.numeric(eval.pts), powers, `^`), eval.pts)
}

additive_effect_profile <- function(fit, pseudo, v, modifier, eval.pts) {
  tt <- stats::delete.response(stats::terms(fit))
  X <- stats::model.matrix(fit)
  term.vars <- lapply(attr(tt, "term.labels"), function(term)
    all.vars(stats::as.formula(paste("~", term))))
  if (any(lengths(term.vars) > 1L) || attr(tt, "intercept") != 1L ||
      !is.null(stats::model.offset(stats::model.frame(fit))) ||
      !isTRUE(all.equal(as.numeric(stats::model.response(stats::model.frame(fit))),
                        as.numeric(pseudo)))) {
    stop("Additive effect profiles require an additive lm with an intercept and the pseudo-outcome as response.")
  }
  terms.j <- which(vapply(term.vars, function(vars) modifier %in% vars, logical(1)))
  cols <- which(attr(X, "assign") %in% terms.j)
  nd <- as.data.frame(v)[rep(1L, length(eval.pts)), , drop=FALSE]
  nd[[modifier]] <- if (is.factor(v[[modifier]])) {
    factor(eval.pts, levels=levels(v[[modifier]]), ordered=is.ordered(v[[modifier]]))
  } else eval.pts
  X0 <- stats::model.matrix(tt, nd, contrasts.arg=fit$contrasts, xlev=fit$xlevels)
  B <- X * 0
  B0 <- X0 * 0
  B[, cols] <- X[, cols, drop=FALSE]
  B0[, cols] <- X0[, cols, drop=FALSE]
  effect_profile(fit, pseudo, B, B0, eval.pts)
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

      p.v.tr <- stats::poly(v.tr, degree=dfs[k], raw=TRUE)
      p.v.te <- stats::poly(v.te, degree=dfs[k], raw=TRUE)

      for(j in 1:dfs[k]){
        res.v[test.idx, j] <- p.v.te[, j] - reg.basis.not.j(y=p.v.tr[, j],
                                                            x=w.tr, new.x=w.te)
      }
      res.y[test.idx] <- pseudo.te - cate.not.j(y=pseudo.tr, x=w.tr, new.x=w.te)
    }

    fit.k <- stats::lm(res.y ~ -1 + res.v)
    fits[[k]] <- fit.k
    diag.hat.mat <- stats::lm.influence(fit.k, do.coef=FALSE)$hat
    risk[k] <- mean((stats::resid(fit.k)/(1-diag.hat.mat))^2)

  }
  fit.star <- fits[[which.min(risk)]]
  risk.dat <- data.frame(dfs=dfs, risk=risk)
  # Effect levels and centered-component inference are computed by
  # robinson_effect_profile(); no zero-anchored prediction/band is needed.
  out <- list(model=fit.star, risk=risk.dat, fits=fits)
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
  # bsc <- function(x, ..., center = TRUE) {
  #   B <- splines::bs(x, ...)
  #   if (center) B <- sweep(B, 2, colMeans(B), "-")
  #   B
  # }
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
      lm.form <- paste0("~ ", paste0("splines::bs(", colnames(x.cont)[1], ", df = ", n.basis[i, 1], ")"))
      if(ncol(x.cont) > 1) {
        for(k in 2:ncol(x.cont)) {
          # lm.form <- c(lm.form, paste0("poly(", colnames(x.cont)[k], ", raw = TRUE, degree = ", n.basis[i, k], ")"))
          lm.form <- c(lm.form, paste0("splines::bs(", colnames(x.cont)[k], ", df = ", n.basis[i, k], ")"))
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
    fits[[i]] <- stats::lm(stats::as.formula(paste0("y", lm.form)),
                           data=cbind(data.frame(y=y), x.cont, x.disc))
    risk[i] <- mean((stats::resid(fits[[i]])/(1-stats::hatvalues(fits[[i]])))^2)
    models[i] <- lm.form
  }
    risk.dat <- cbind(n.basis, risk)
    if(ncol(x.cont)>0) colnames(risk.dat) <- c(colnames(x.cont), "loocv.risk")
    if(ncol(x.cont)==0) colnames(risk.dat) <- c(colnames(x.disc), "loocv.risk")
    # other choices are possible, always plot the estimates risks!
    best.model <- stats::lm(stats::as.formula(paste0("y", models[which.min(risk)])),
                            data=cbind(data.frame(y=y), x))

  out <- stats::predict(best.model, newdata=as.data.frame(new.x))
  res <- cbind(out, NA, NA)
  return((list(drl.form=stats::formula(best.model), res=res, model=best.model,
               risk=risk.dat, fits=fits)))
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
  form <- stats::formula(fit)
  mm <- stats::model.matrix(fit)
  coefs <- stats::coef(fit)

  bs_term_for_vj <- grep(eff.modif.name, colnames(mm), value=TRUE)
  coefs.names.vj <- grep(eff.modif.name, names(coefs), value=TRUE)

  coefs.vj <- coefs[coefs.names.vj]
  tt <- stats::delete.response(stats::terms(fit))
  new.design.mat <- as.matrix(stats::model.matrix(tt, new.dat.additive)[, bs_term_for_vj])

  if(!is.factor(v[, j])) {
    design.mat <- as.matrix(stats::model.matrix(tt, v)[, bs_term_for_vj])
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
