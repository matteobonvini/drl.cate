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


drl.basis.additive <- function(y, x, new.x, kmin=1, kmax=10) {
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

  best.model <- lm(as.formula(paste0("y", models[which.min(risk)])),
                   data=cbind(data.frame(y=y), x))

  out <- predict(best.model, newdata=as.data.frame(new.x))
  res <- cbind(out, NA, NA)
  return((list(drl.form=models[which.min(risk)], res=res, model=best.model,
               risk=risk.dat, fits)))
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
