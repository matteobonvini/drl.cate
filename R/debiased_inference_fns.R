# These functions are slight modifications of those contained in the R package
# for dose-response debiased inference available at
# https://github.com/Kenta426/DebiasedDoseResponse
# Main reference is https://arxiv.org/abs/2210.06448
# Article and repository authored by Kenta Takatsu and Ted Westling.

.get.muhat <- function(splits.id, cate.w.fit, v1, v2, max.n.integral=1000) {

  nsplits <- length(cate.w.fit)
  res <- list()
  counter <- 1

  for(w in 1:nsplits) {

    n.k <- sum(splits.id==w)
    n.subsplits <- max(round(n.k / max.n.integral), 1)
    sub.s <- sample(rep(1:n.subsplits, ceiling(n.k/nsplits))[1:n.k])

    for(u in 1:n.subsplits){
      sub.idx <- w==splits.id
      sub.idx[!u==sub.s] <- FALSE
      n.small <- sum(sub.idx)
      v2.idx <- v2[sub.idx, , drop=FALSE]
      v2.idx.rep <- v2.idx[rep(1:n.small, n.small), , drop=FALSE]
      new.w <- cbind(v1j=sample(v1, nrow(v2.idx.rep), replace=TRUE),
                     v2.idx.rep)
      muhat.mat <- matrix(cate.w.fit[[w]]$fit(new.w), nrow=n.small,
                          ncol=n.small)
      res[[counter]] <- list(muhat.mat=muhat.mat, sub.idx=sub.idx)
      counter <- counter + 1
    }
  }
  return(res)
}

.compute.rinfl.func <- function(Y, A, a, h, b, kern, debias, muhat.vals=NULL,
                                mhat.vals=NULL){
  #####################################################
  ## ATTN: check with Kenta the definition of if.c2! ##
  #####################################################
  n <- length(A)
  bw.min <- sort(abs(unique(A - a)))[10]
  h <- max(h, bw.min); b <- max(b, bw.min)
  a.std.h <- (A-a)/h
  kern.std.h <- kern(a.std.h)/h
  a.std.b <- (A-a)/b
  kern.std.b <- kern(a.std.b)/b
  if(!is.null(muhat.vals)) {
    int1.h <- int2.h <- int1.b <- int2.b <- int3.b <- rep(NA, n)
    nsubsplits <- length(muhat.vals)
    for(i in 1:nsubsplits) {

      sub.idx <- muhat.vals[[i]]$sub.idx
      muhat.mat <- muhat.vals[[i]]$muhat.mat
      mhat <- mhat.vals[sub.idx]

      int1.h[sub.idx] <- colMeans(kern.std.h[sub.idx]*(muhat.mat-mhat))
      int2.h[sub.idx] <- colMeans(a.std.h[sub.idx]*kern.std.h[sub.idx]*(muhat.mat-mhat))
      if(debias) {
        int1.b[sub.idx] <- colMeans(kern.std.b[sub.idx]*(muhat.mat-mhat))
        int2.b[sub.idx] <- colMeans(a.std.b[sub.idx]*kern.std.b[sub.idx]*(muhat.mat-mhat))
        int3.b[sub.idx] <- colMeans((a.std.b[sub.idx])^2*kern.std.b[sub.idx]*(muhat.mat-mhat))
      }
    }
  } else {
    int1.h <- int2.h <- int1.b <- int2.b <- int3.b <- rep(0, n)
  }

  c0.h <- mean(kern.std.h)
  c1.h <- mean(kern.std.h*a.std.h)
  c2.h <- mean(kern.std.h*a.std.h^2)
  Dh <- matrix(c(c0.h, c1.h,
                 c1.h, c2.h), nrow=2)
  Dh.inv <- solve(Dh)
  gamma.h <- coef(lm(Y ~ a.std.h, weights=kern.std.h))
  res.h <- Y - (gamma.h[1] + gamma.h[2]*a.std.h)
  inf.fn <- t(Dh.inv %*% rbind(res.h*kern.std.h + int1.h,
                               a.std.h*res.h*kern.std.h + int2.h))

  if(debias){
    c0.b <- mean(kern.std.b)
    c1.b <- mean(kern.std.b*a.std.b)
    c2.b <- mean(kern.std.b*a.std.b^2)
    c3.b <- mean(kern.std.b*a.std.b^3)
    c4.b <- mean(kern.std.b*a.std.b^4)

    Db <- matrix(c(c0.b, c1.b, c2.b,
                   c1.b, c2.b, c3.b,
                   c2.b, c3.b, c4.b), nrow=3)
    Db.inv <- solve(Db)

    # c2 <- integrate(function(u){u^2 * kern(u)}, -Inf, Inf)$value # old c2
    w.vec <- cbind(kern.std.h, kern.std.h*a.std.h)
    c2.vec <- ((Dh.inv %*% crossprod(w.vec, a.std.h^2))/n)
    c2 <- c2.vec[1]
    # the EIF for the fixed-h c2 -----------------------------------------------
    w.1 <- cbind(1, a.std.h)
    w.1.tilde <- cbind(kern.std.h*a.std.h^2, kern.std.h*a.std.h^3)
    term1 <- (w.1%*%Dh.inv)[,1]*c(w.1 %*% c2.vec)*kern.std.h
    term2 <- (w.1.tilde %*% Dh.inv)[,1]
    deriv2 <- (Db.inv%*%crossprod(cbind(1, a.std.b, a.std.b^2)*kern.std.b, Y)/n)[3,]
    # if.c2 <- h^2/2*deriv2*(term2-term1)
    if.c2 <- (h/b)^2*deriv2*(term2-term1)

    model.b <- lm(Y ~ poly(a.std.b, 3), weights=kern.std.b)
    res.b <- Y - predict(model.b)
    inf.fn.robust <- t(Db.inv %*% rbind(res.b*kern.std.b + int1.b,
                                        a.std.b*res.b*kern.std.b + int2.b,
                                        a.std.b^2*res.b*kern.std.b + int3.b))
  }
  if(debias){
    out <- data.frame(est=inf.fn[,1] - (h/b)^2*c2*inf.fn.robust[,3] - if.c2)
  } else {
    out <- data.frame(est=inf.fn[,1])
  }
  return(out)
}


.parse.debiased_inference <- function(...){
  option <- list(...); arg <- list()
  if (is.null(option$kernel.type)){
    kernel.type <- "epa"
  }
  else{
    kernel.type <- option$kernel.type
  }
  if (is.null(option$bandwidth.method)){
    bandwidth.method <- "DPI"
  }
  else{
    bandwidth.method <- option$bandwidth.method
  }
  if (is.null(option$alpha.pts)){
    alpha.pts <- 0.05
  }
  else{
    alpha.pts <- option$alpha.pts
  }
  if (is.null(option$unif)){
    unif <- TRUE
  }
  else{
    unif <- option$unif
  }
  if (is.null(option$alpha.unif)){
    alpha.unif <- 0.05
  }
  else{
    alpha.unif <- option$alpha.unif
  }
  if (is.null(option$bootstrap)){
    bootstrap <- 1e4
  }
  else{
    bootstrap <- option$bootstrap
  }
  arg$kernel.type <- kernel.type
  arg$bandwidth.method <- bandwidth.method
  arg$alpha.pts <- alpha.pts
  arg$unif <- unif
  arg$alpha.unif <- alpha.unif
  arg$bw.seq <- option$bw.seq
  arg$bootstrap <- bootstrap
  arg$mu <- option$mu
  arg$g <- option$g
  return(arg)
}

#' plot_debiased_curve
#'
#' @param res.df ADD
#' @param ci ADD
#' @param unif ADD
#' @export
plot_debiased_curve <- function(pseudo, exposure, res.df, ci=TRUE, unif=TRUE,
                                add.pseudo=TRUE){
  p <-  ggplot2::ggplot() + ggplot2::xlab("Exposure") +
    ggplot2::ylab("Covariate-adjusted outcome") +
    ggplot2::theme_minimal()
  if(class(res.df$eval.pts) == "factor") {
    if(add.pseudo) {
      p <- p + ggplot2::geom_point(data = NULL, aes(x=as.factor(exposure), y=pseudo), col = "gray")
    }
    p <- p + ggplot2::geom_point(data = res.df, aes(x=eval.pts, y=theta))
  } else {
    if(add.pseudo) {
      p <- p + ggplot2::geom_point(data = NULL, aes(x=exposure, y=pseudo), col = "gray")
    }
    p <- p + ggplot2::geom_line(data = res.df, aes(x=eval.pts, y=theta))
  }
  if(ci){
    p <- p + ggplot2::geom_pointrange(data = res.df, aes(x=eval.pts, y=theta,
                                                         ymin=ci.ll.pts,
                                                         ymax=ci.ul.pts,
                                                         size="Pointwise CIs"), col = "black") +
      ggplot2::scale_size_manual("",values=c("Pointwise CIs"=0.2))
  }
  if(unif){
    p <- p + ggplot2::geom_line(data = res.df, aes(x=eval.pts,
                                                   y=ci.ll.unif,
                                                   linetype="Uniform band"), col = "red") +
      ggplot2::geom_line(data = res.df, aes(x=eval.pts,
                                            y=ci.ul.unif,
                                            linetype="Uniform band"),
                         col = "red")+
      ggplot2::scale_linetype_manual("",values=c("Uniform band"=2))
  }
  p <- p + ggplot2::theme(legend.position = "bottom") +
    ggplot2::geom_hline(yintercept = 0, col = "blue") +
    ggplot2::geom_hline(yintercept = mean(pseudo), col = "orange")
  return(p)
}

#' debiased_inference
#' @param A ADD
#' @param pseudo.out ADD
#' @param muhat.mat ADD
#' @param mhat.obs ADD
#' @param debias boolean ADD
#' @param tau ratio bandwidths ADD
#' @param eval.pts ADD
#' @param ... ADD
#' @export
debiased_inference <- function(A, pseudo.out, debias, tau=1, eval.pts=NULL,
                               muhat.vals=NULL, mhat.obs=NULL, ...){
  # Parse control inputs ------------------------------------------------------
  # control <- .parse.debiased_inference(alpha=0.05, unif=FALSE, kernel.type = "gau",
  #                                      eval.pts = x, A = x, pseudo.out=y, debias=FALSE,
  #                                      muhat.vals = NULL, mahat.obs=NULL, bw.seq = c(0.1, 0.5))
  control <- .parse.debiased_inference(...)
  kernel.type <- control$kernel.type
  # Compute an estimated pseudo-outcome sequence ------------------------------
  # ord <- order(A)
  # A <- A[ord]
  # pseudo.out <- pseudo.out[ord]
  n <- length(A)
  kern <- function(t){.kern(t, kernel=kernel.type)}
  if (is.null(eval.pts)){
    eval.pts <- seq(quantile(A, 0.05), quantile(A, 0.95),
                    length.out=min(30, length(unique(A))))
  }
  # Compute bandwidth ---------------------------------------------------------
  bw.seq <- control$bw.seq
  if(debias) {
    if(control$bandwidth.method=="LOOCV"){
      bw.seq.h <- rep(bw.seq, length(bw.seq))
      bw.seq.b <- rep(bw.seq, each=length(bw.seq))
    } else if(control$bandwidth.method=="LOOCV(h=b)"){
      bw.seq.h <- bw.seq.b <- bw.seq
    } else stop("Specify a valid bandwidth method.")
  } else {
    bw.seq.h <-  bw.seq.b <- bw.seq
  }


  est.proc <- function(h, b){
    est.res <- matrix(NA, ncol=2, nrow=length(eval.pts),
                      dimnames=list(NULL, c("eval", "theta.hat")))

    # good.pts <- sapply(eval.pts,
    #                    function(u) {
    #                      length(unique(A[.kern((A-u)/min(h, b), kernel.type)>1e-6]))>4
    #                          })
    # if(mean(good.pts) < 0.8) {
    #   res <- data.frame(
    #     eval.pts=eval.pts,
    #     theta=NA,
    #     ci.ul.pts=NA,
    #     ci.ll.pts=NA,
    #     ci.ul.unif=NA,
    #     ci.ll.unif=NA,
    #     if.val.sd=NA,
    #     unif.quantile=NA,
    #     h=h,
    #     b=b,
    #     loocv.risk=Inf)
    #   return(res)
    # }
    # good.eval.pts <- eval.pts[good.pts]


    est.res <- .lprobust(x=A, y=pseudo.out, h=h, b=b,
                         debias=debias, eval=eval.pts,
                         kernel.type=kernel.type)

    hat.val <- .hatmatrix(x=A, y=pseudo.out, h=h, b=b, debias=debias,
                          eval.pt=eval.pts, kernel.type=kernel.type)

    est.fn <- approx(eval.pts, est.res[,"theta.hat"], xout=A, rule=2)$y
    sq.resid <- ((pseudo.out-est.fn)/(1-hat.val))^2
    loocv.risk <- mean(sq.resid)

    # Estimate influence function sequence --------------------------------------
    rinf.fns <- mapply(function(a, h.val, b.val){
      .compute.rinfl.func(Y=pseudo.out, A=A, a=a, h=h.val, b=b.val, kern=kern,
                          debias=debias, muhat.vals=muhat.vals, mhat.vals=mhat.obs)
    }, eval.pts, h, b, SIMPLIFY=FALSE)
    rif.se <- matrix(NA, ncol=1, nrow=length(eval.pts), dimnames=list(NULL, "est"))
    rif.se <- do.call(rbind, lapply(rinf.fns, function(u) {
                                                 apply(u, 2, sd)/sqrt(n)
                                               }))
    alpha <- control$alpha.pts
    z.val <- qnorm(1-alpha/2)
    ci.ll.p <- ci.ul.p <- rep(NA, length(eval.pts))
    ci.ll.p <- est.res[,"theta.hat"] - z.val*rif.se[,"est"]
    ci.ul.p <- est.res[,"theta.hat"] + z.val*rif.se[,"est"]

    # Compute uniform bands by simulating GP ------------------------------------
    if(control$unif){
      get.unif.ep <- function(alpha){
        std.inf.vals <- do.call(cbind, lapply(rinf.fns, function(u) scale(u)))
        boot.samples <- control$bootstrap
        ep.maxes<- replicate(boot.samples,
                             max(abs(rbind(rnorm(n)/sqrt(n)) %*% std.inf.vals)))
        quantile(ep.maxes, alpha)
      }
      alpha <- control$alpha.unif
      unif.quantile <- get.unif.ep(1-alpha)
      names(unif.quantile) <- NULL
      ep.unif.quant <- rep(NA, length(eval.pts))
      ep.unif.quant <- unif.quantile*rif.se[, "est"]
      ci.ll.u <- ci.ul.u <- rep(NA, length(eval.pts))
      ci.ll.u <- est.res[,"theta.hat"] - ep.unif.quant
      ci.ul.u <- est.res[,"theta.hat"] + ep.unif.quant
    }
    else{
      unif.quantile <- ci.ll.u <- ci.ul.u <- NA
    }

    # Output data.frame ---------------------------------------------------------
    res <- data.frame(
      eval.pts=eval.pts,
      theta=est.res[,"theta.hat"],
      ci.ul.pts=ci.ul.p,
      ci.ll.pts=ci.ll.p,
      ci.ul.unif=ci.ul.u,
      ci.ll.unif=ci.ll.u,
      if.val.sd=rif.se[,"est"],
      unif.quantile=unif.quantile,
      h=h,
      b=b,
      loocv.risk=loocv.risk)
  }

  res.list <- mapply(est.proc, bw.seq.h, bw.seq.b, SIMPLIFY=FALSE)
  loocv.vals <- matrix(NA, ncol = 3, nrow = length(res.list),
                       dimnames=list(NULL, c("loocv.risk", "h", "b")))
  for(k in 1:length(res.list)) {
    loocv.vals[k, ] <- c(res.list[[k]]$loocv.risk[1],
                         res.list[[k]]$h[1],
                         res.list[[k]]$b[1])
  }
  best.loocv.idx <- abs(loocv.vals[, 1] - min(loocv.vals[, 1])) < 1e-6
  h.opt <- max(loocv.vals[best.loocv.idx, 2])
  b.opt <- max(loocv.vals[best.loocv.idx, 3])
  res <- est.proc(h=h.opt, b=b.opt)
  return(list(res=res, risk=as.data.frame(loocv.vals), res.list=res.list))
}

# Compute hat matrix for leave-one-out cross validation for a linear smoother
.hatmatrix <- function(x, y, h, b, eval.pt=NULL, kernel.type="epa", debias=TRUE){
  # ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval.pt)){
    eval.pt <- seq(quantile(x, 0.05), quantile(x, 0.95), length.out=30)
  }

  neval <- length(eval.pt)

  e1 <- matrix(c(1, 0), ncol=2); e3 <- matrix(c(1, 0, 0), ncol=3)
  w.h.zero   <- .kern(0, kernel.type)/h
  w.b.zero   <- .kern(0, kernel.type)/b
  hat.mat <- rep(0, neval)
  for (i in 1:neval){
    bw.min   <- sort(abs(unique(x-eval.pt[i])))[10]
    h0     <- max(h, bw.min); b0 <- max(b, bw.min)
    k.h   <- .kern((x-eval.pt[i])/h0, kernel.type)/h0
    k.b   <- .kern((x-eval.pt[i])/b0, kernel.type)/b0
    ind.h <- k.h>0;  ind.b <- k.b>0
    N.h   <- sum(ind.h);  N.b <- sum(ind.b); ind <- ind.b
    if (h>b) ind <- ind.h
    eY  <- y[ind]; eX  <- x[ind]
    K.h <- k.h[ind]; K.b <- k.b[ind]
    W.b <- matrix(NA, sum(ind), 3) # (1, u, u^2)
    for (j in 1:3)  W.b[,j] <- ((eX-eval.pt[i])/b0)^(j-1)
    W.h  <- matrix(NA, sum(ind), 2) # (1, u)
    for (j in 1:2)  W.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)
    invD.b  <- .qrXXinv((sqrt(K.b)*W.b))
    invD.h  <- .qrXXinv((sqrt(K.h)*W.h))
    if(debias) {
      # c.2 <- integrate(function(t){t^2*.kern(t, kernel.type)}, -Inf, Inf)$value
      u <- (eX-eval.pt[i])/h0
      c.2 <- (invD.h %*% crossprod(W.h*K.h,u^(2)))[1]
    } else {
      c.2 <- 0
    }
    hat.mat[i] <- (invD.h%*%t(e1*w.h.zero))[1,] -
      ((h/b)^2*c.2*invD.b%*%t(e3*w.b.zero))[3,]
  }
  return(approx(eval.pt, hat.mat, xout=x, rule=2)$y)
}


# Perform leave-one-out cross-validation based on the "hat matrix trick"
.robust.loocv <- function(x, y, h, b, debias, eval.pt=NULL, kernel.type="epa"){
  # ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval.pt)){
    eval.pt <- seq(quantile(x, 0.05), quantile(x, 0.95), length.out=30)
  }
  hat.val <- .hatmatrix(x=x, y=y, h=h, b=b, debias=debias, eval.pt=eval.pt,
                        kernel.type=kernel.type)

  est <- .lprobust(x=x, y=y, h=h, b=b, debias=debias, eval.pt=eval.pt,
                   kernel.type=kernel.type)
  est.fn <- approx(eval.pt, est[,"theta.hat"], xout=x, rule=2)$y
  sq.resid <- ((y-est.fn)/(1-hat.val))^2
  loocv.risk <- mean(sq.resid)
  return(data.frame(loocv.risk=loocv.risk))
}


# Debiased local linear regression
.lprobust <- function(x, y, h, b, debias, eval.pt=NULL, kernel.type="epa"){
  # ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval.pt)){
    eval.pt <- seq(quantile(x, 0.05), quantile(x, 0.95),  length.out=30)
  }
  neval <- length(eval.pt)
  estimate.mat <- matrix(NA, neval, 2)
  colnames(estimate.mat) <- c("eval", "theta.hat")
  for (i in 1:neval){
    bw.min <- sort(abs(unique(x-eval.pt[i])))[10]
    h0 <- max(h, bw.min); b0 <- max(b, bw.min)
    k.h <- .kern((x-eval.pt[i])/h0, kernel.type)/h0
    k.b <- .kern((x-eval.pt[i])/b0, kernel.type)/b0
    ind.h <- k.h>0;  ind.b <- k.b>0; ind <- ind.b
    if (h>b) ind <- ind.h
    eN <- sum(ind); eY <- y[ind]; eX <- x[ind]
    K.h <- k.h[ind]; K.b <- k.b[ind]
    W.b <- matrix(NA, eN, 3)
    for (j in 1:3)  W.b[, j] <- ((eX-eval.pt[i])/b0)^(j-1)
    W.h <- matrix(NA, eN, 2)
    for (j in 1:2)  W.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)
    invD.b  <- .qrXXinv((sqrt(K.b)*W.b))
    invD.h  <- .qrXXinv((sqrt(K.h)*W.h))
    beta.ll  <- invD.h%*%crossprod(W.h*K.h, eY) #R.p^T W.h Y

    if(debias) {
      # c.2 <- integrate(function(t){t^2*.kern(t, kernel.type)}, -Inf, Inf)$value
      u <- (eX-eval.pt[i])/h0
      c.2 <- (invD.h %*% crossprod(W.h*K.h,u^2))[1]
    } else {
      c.2 <- 0
    }

    beta.bias <- (h0/b0)^2*invD.b%*%crossprod(W.b*K.b, eY)*c.2

    estimate.mat[i,] <- c(eval.pt[i], beta.ll[1,1]-beta.bias[3,1])
  }
  estimate.mat
}

# Generate kernel function
.kern = function(u, kernel="epa"){
  if (kernel=="epa") w <- 0.75*(1-u^2)*(abs(u)<=1)
  if (kernel=="uni") w <- 0.5*(abs(u)<=1)
  if (kernel=="tri") w <- (1-abs(u))*(abs(u)<=1)
  if (kernel=="gau") w <- dnorm(u)
  return(w)
}

# Matrix inverse via Cholesky decomposition
.qrXXinv = function(x, ...) {
  chol2inv(chol(crossprod(x)))
}

# # Local polynomial regression
# .locpoly <- function(x, y, h, eval.pt=NULL, kernel.type="epa", degree=1){
#   ind <- order(x); x <- x[ind]; y <- y[ind]
#   if (is.null(eval.pt)){
#     eval.pt <- seq(quantile(x, 0.05), quantile(x, 0.95),  length.out=30)
#   }
#   neval <- length(eval.pt)
#   estimate.mat <- matrix(NA, neval, 2)
#   colnames(estimate.mat) <- c("eval", "theta.hat")
#   # Compute for each evaluation point
#   for (i in 1:neval){
#     bw.min   <- sort(abs(x-eval.pt[i]))[10]
#     h0     <- max(h, bw.min)
#     k.h   <- .kern((x-eval.pt[i])/h0, kernel=kernel.type)/h0
#     ind <- k.h>0
#     eY  <- y[ind]; eX  <- x[ind]; K.h <- k.h[ind]
#     W.h <- matrix(NA, sum(ind), degree+1)
#     for (j in 1:(degree+1))  W.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)
#     invD.h <- .qrXXinv((sqrt(K.h)*W.h))
#     beta <- invD.h%*%crossprod(W.h*K.h, eY)
#     estimate.mat[i,] <- c(eval.pt[i], beta[1,1])
#   }
#   estimate.mat
# }

# # Local linear regression
# .loclinear <- function(x, y, h, eval.pt=NULL, kernel.type="epa"){
#   ind <- order(x); x <- x[ind]; y <- y[ind]
#   if (is.null(eval.pt)){
#     eval.pt <- seq(quantile(x, 0.05), quantile(x, 0.95),  length.out=30)
#   }
#   neval <- length(eval.pt)
#   estimate.mat <- matrix(NA, neval, 2)
#   colnames(estimate.mat) <- c("eval","theta.hat")
#   for (i in 1:neval){
#     bw.min   <- sort(abs(x-eval.pt[i]))[10]
#     h0 <- max(h, bw.min)
#     k.h <- .kern((x-eval.pt[i])/h0, kernel=kernel.type)/h0
#     ind <- k.h>0
#     eY <- y[ind]; eX <- x[ind]; K.h <- k.h[ind]
#     w.h <- matrix(NA, sum(ind), 2)
#     for (j in 1:2)  w.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)
#     invD.h  <- .qrXXinv((sqrt(K.h)*w.h))
#     beta.p  <- invD.h%*%crossprod(w.h*K.h, eY)
#     estimate.mat[i,] <- c(eval.pt[i], beta.p[1,1])
#   }
#   estimate.mat
# }
#
#
# .hatmatrix2 <- function(x, y, h, b, eval.pt=NULL, kernel.type="epa", debias=TRUE){
#   # ind <- order(x); x <- x[ind]; y <- y[ind]
#   if (is.null(eval.pt)){
#     eval.pt <- seq(quantile(x, 0.05), quantile(x, 0.95), length.out=30)
#   }
#   neval <- length(eval.pt)
#   if(debias) {
#     c.2 <- integrate(function(t){t^2*.kern(t, kernel.type)}, -Inf, Inf)$value
#   } else {
#     c.2 <- 0
#   }
#
#   e1 <- matrix(c(1, 0), ncol=2); e3 <- matrix(c(1, 0, 0, 0), ncol=4)
#   w.h.zero   <- .kern(0, kernel.type)/h
#   w.b.zero   <- .kern(0, kernel.type)/b
#   hat.mat <- rep(0, neval)
#   for (i in 1:neval){
#     bw.min   <- sort(abs(x-eval.pt[i]))[10]
#     bw.min <- 0
#     h0     <- max(h, bw.min); b0 <- max(b, bw.min)
#
#     k.h   <- .kern((x-eval.pt[i])/h0, kernel.type)/h0
#     k.b   <- .kern((x-eval.pt[i])/b0, kernel.type)/b0
#     ind.h <- k.h>0;  ind.b <- k.b>0
#     N.h   <- sum(ind.h);  N.b <- sum(ind.b); ind <- ind.b
#     if (h>b) ind <- ind.h
#     eY  <- y[ind]; eX  <- x[ind]
#     K.h <- k.h[ind]; K.b <- k.b[ind]
#     W.b <- matrix(NA, sum(ind), 4) # (1, u, u^2, u^3)
#     for (j in 1:4)  W.b[,j] <- ((eX-eval.pt[i])/b0)^(j-1)
#     W.h  <- matrix(NA, sum(ind), 2) # (1, u)
#     for (j in 1:2)  W.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)
#     # print(W.h)
#     # print(unique(sqrt(K.b)*W.b))
#     # solve(crossprod(unique(sqrt(K.b)*W.b)))
#     invD.b  <- .qrXXinv((sqrt(K.b)*W.b))
#     invD.h  <- .qrXXinv((sqrt(K.h)*W.h))
#     hat.mat[i] <- (invD.h%*%t(e1*w.h.zero))[1,] -
#       ((h/b)^2*c.2*invD.b%*%t(e3*w.b.zero))[3,]
#   }
#   return(hat.mat)
# }



# debiased_inference <- function(A, pseudo.out, debias, tau=1, eval.pts=NULL,
#                                muhat.vals=NULL, mhat.obs=NULL, ...){
#   # Parse control inputs ------------------------------------------------------
#   # control <- .parse.debiased_inference(alpha=0.05, unif=FALSE, kernel.type = "gau",
#   #                                      eval.pts = x, A = x, pseudo.out=y, debias=FALSE,
#   #                                      muhat.vals = NULL, mahat.obs=NULL)
#   control <- .parse.debiased_inference(...)
#   kernel.type <- control$kernel.type
#   # Compute an estimated pseudo-outcome sequence ------------------------------
#   # ord <- order(A)
#   # A <- A[ord]
#   # pseudo.out <- pseudo.out[ord]
#   n <- length(A)
#   kern <- function(t){.kern(t, kernel=kernel.type)}
#   if (is.null(eval.pts)){
#     eval.pts <- seq(quantile(A, 0.05), quantile(A, 0.95),
#                     length.out=min(30, length(unique(A))))
#   }
#   # Compute bandwidth ---------------------------------------------------------
#   bw.seq <- control$bw.seq
#   if(debias) {
#     if(control$bandwidth.method=="LOOCV"){
#       bw.seq.h <- rep(bw.seq, length(bw.seq))
#       bw.seq.b <- rep(bw.seq, each=length(bw.seq))
#     } else if(control$bandwidth.method=="LOOCV(h=b)"){
#       bw.seq.h <- bw.seq.b <- bw.seq
#     } else stop("Specify a valid bandwidth method.")
#   } else {
#     bw.seq.h <-  bw.seq.b <- bw.seq
#   }
#
#
#   risk <- mapply(function(h, b){
#     count.pts.nearby <- Vectorize(function(u){
#       # count the unique points nearby u
#       length(unique(A[.kern((A-u)/min(h, b), kernel.type) > 1e-6]))
#     })
#
#     # keep only those evaluation points for which we have at least 4 observations
#     # in a min(h, b)-neighborhood.
#     good.pts <- count.pts.nearby(eval.pts) > 4
#     good.eval.pts <- eval.pts[good.pts]
#
#     if(mean(good.pts) > 0.8) {
#       # If we can retain at least 80% of the evaluation points
#       # we compute the loocv estimated risk
#       loocv.res <- .robust.loocv(x=A, y=pseudo.out, h=h, b=b, debias=debias,
#                                  eval.pt=good.eval.pts, kernel.type=kernel.type)
#     } else {
#       # otherwise, we set the risk to infinity
#       loocv.res <- data.frame(loocv.risk=Inf)
#     }
#
#     return(cbind(prop.good.pts=mean(good.pts), loocv.res))
#
#   }, bw.seq.h, bw.seq.b, SIMPLIFY=FALSE)
#   risk <- do.call(rbind, risk)
#   risk$h <- bw.seq.h
#   risk$b <- bw.seq.b
#   h.opt <- bw.seq.h[which.min(risk[, "loocv.risk"])]
#   b.opt <- bw.seq.b[which.min(risk[, "loocv.risk"])]
#   est.res <- matrix(NA, ncol=2, nrow=length(eval.pts),
#                     dimnames=list(NULL, c("eval", "theta.hat")))
#   good.pts.h.opt <- sapply(eval.pts,
#                            function(u) {
#                              length(unique(A[.kern((A-u)/min(h.opt, b.opt), kernel.type)>1e-6]))>4
#                              })
#   good.eval.pts.h.opt <- eval.pts[good.pts.h.opt]
#
#   est.res[good.pts.h.opt, ] <- .lprobust(x=A, y=pseudo.out, h=h.opt, b=b.opt,
#                                          debias=debias, eval=good.eval.pts.h.opt,
#                                          kernel.type=kernel.type)
#
#   # Estimate influence function sequence --------------------------------------
#   rinf.fns <- mapply(function(a, h.val, b.val){
#     .compute.rinfl.func(Y=pseudo.out, A=A, a=a, h=h.val, b=b.val, kern=kern,
#                         debias=debias, muhat.vals=muhat.vals, mhat.vals=mhat.obs)
#   }, good.eval.pts.h.opt, h.opt, b.opt, SIMPLIFY=FALSE)
#   rif.se <- matrix(NA, ncol=1, nrow=length(eval.pts), dimnames=list(NULL, "est"))
#   rif.se[good.pts.h.opt, ] <- do.call(rbind,
#                                       lapply(rinf.fns,
#                                              function(u) {
#                                                apply(u, 2, sd)/sqrt(n)
#                                                }))
#   alpha <- control$alpha.pts
#   z.val <- qnorm(1-alpha/2)
#   ci.ll.p <- ci.ul.p <- rep(NA, length(eval.pts))
#   ci.ll.p[good.pts.h.opt] <- est.res[good.pts.h.opt,"theta.hat"] -
#     z.val*rif.se[good.pts.h.opt, "est"]
#   ci.ul.p[good.pts.h.opt] <- est.res[good.pts.h.opt,"theta.hat"] +
#     z.val*rif.se[good.pts.h.opt, "est"]
#
#   # Compute uniform bands by simulating GP ------------------------------------
#   if(control$unif){
#     get.unif.ep <- function(alpha){
#       std.inf.vals <- do.call(cbind, lapply(rinf.fns, function(u) scale(u)))
#       boot.samples <- control$bootstrap
#       ep.maxes<- replicate(boot.samples,
#                            max(abs(rbind(rnorm(n)/sqrt(n)) %*% std.inf.vals)))
#       quantile(ep.maxes, alpha)
#     }
#     alpha <- control$alpha.unif
#     unif.quantile <- get.unif.ep(1-alpha)
#     names(unif.quantile) <- NULL
#     ep.unif.quant <- rep(NA, length(eval.pts))
#     ep.unif.quant[good.pts.h.opt] <- unif.quantile*rif.se[good.pts.h.opt, "est"]
#     ci.ll.u <- ci.ul.u <- rep(NA, length(eval.pts))
#     ci.ll.u[good.pts.h.opt] <- est.res[good.pts.h.opt,"theta.hat"] - ep.unif.quant[good.pts.h.opt]
#     ci.ul.u[good.pts.h.opt] <- est.res[good.pts.h.opt,"theta.hat"] + ep.unif.quant[good.pts.h.opt]
#   }
#   else{
#     unif.quantile <- ci.ll.u <- ci.ul.u <- NA
#   }
#
#   # Output data.frame ---------------------------------------------------------
#   res <- data.frame(
#     eval.pts=eval.pts,
#     theta=est.res[, "theta.hat"],
#     ci.ul.pts=ci.ul.p,
#     ci.ll.pts=ci.ll.p,
#     ci.ul.unif=ci.ul.u,
#     ci.ll.unif=ci.ll.u,
#     if.val.sd=rif.se[, "est"],
#     unif.quantile=unif.quantile,
#     h=h.opt,
#     b=b.opt)
#   return(list(res=res, risk=risk))
# }
