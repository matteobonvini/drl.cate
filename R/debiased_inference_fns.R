# These functions are slight modifications of those contained in the R package
# for dose-response debiased inference available at
# https://github.com/Kenta426/DebiasedDoseResponse
# Main reference is https://arxiv.org/abs/2210.06448
# Article and repository authored by Kenta Takatsu and Ted Westling.

.get.muhat <- function(splits.id, cate.w.fit, v1, v2, max.n.integral = 1000) {

  nsplits <- length(cate.w.fit)
  res <- list()
  counter <- 1
  for(w in 1:nsplits) {

    n.k <- sum(splits.id == w)
    n.subsplits <- max(round(n.k / max.n.integral), 1)
    sub.s <- sample(rep(1:n.subsplits, ceiling(n.k / nsplits))[1:n.k])

    for(u in 1:n.subsplits){
      sub.idx <- w == splits.id
      sub.idx[!u == sub.s] <- FALSE
      n.small <- sum(sub.idx)
      v2.idx <- v2[sub.idx, , drop = FALSE]
      v2.idx.rep <- v2.idx[rep(1:n.small, n.small), , drop = FALSE]
      new.w <- cbind(v1j = sample(v1, nrow(v2.idx.rep), replace = TRUE),
                     v2.idx.rep)
      muhat.mat <- matrix(cate.w.fit[[w]]$fit(new.w), nrow = n.small,
                          ncol = n.small)
      res[[counter]] <- list(muhat.mat = muhat.mat,
                             sub.idx = sub.idx)
      counter <- counter + 1
    }
  }
  return(res)
}

.compute.rinfl.func <- function(Y, A, a, h, b, kern, muhat.vals = NULL,
                                mhat.vals = NULL){

  n <- length(A); bw.min <- sort(abs(A - a))[21]
  int1.h <- int2.h <- int1.b <- int2.b <- int3.b <- int4.b <-
    rep(0, n)
  h <- max(h, bw.min); b <- max(b, bw.min)
  a.std.h <- (A - a)/h; kern.std.h <- kern(a.std.h)/h
  a.std.b <- (A - a)/b; kern.std.b <- kern(a.std.b)/b

  # Compute 2x2 inverse matrix ------------------------------------------------
  c0.h <- mean(kern.std.h)
  c1.h <- mean(kern.std.h * a.std.h)
  c2.h <- mean(kern.std.h * a.std.h^2)
  Dh <- matrix(c(c0.h, c1.h,
                 c1.h, c2.h), nrow = 2)
  Dh.inv <- solve(Dh)

  # Compute 4x4 inverse matrix ------------------------------------------------
  c0.b <- mean(kern.std.b)
  c1.b <- mean(kern.std.b * a.std.b)
  c2.b <- mean(kern.std.b * a.std.b^2)
  c3.b <- mean(kern.std.b * a.std.b^3)
  c4.b <- mean(kern.std.b * a.std.b^4)
  c5.b <- mean(kern.std.b * a.std.b^5)
  c6.b <- mean(kern.std.b * a.std.b^6)
  Db <- matrix(c(c0.b, c1.b, c2.b, c3.b,
                 c1.b, c2.b, c3.b, c4.b,
                 c2.b, c3.b, c4.b, c5.b,
                 c3.b, c4.b, c5.b, c6.b), nrow = 4)
  Db.inv <- solve(Db)

  # Estimate local linear components ------------------------------------------
  g2.h <- (A - a)/h
  g2.b <- (A - a)/b
  g3.b <- ((A - a)/b)^2
  g4.b <- ((A - a)/b)^3

  if(is.null(muhat.vals)) {
    int1.h <- int2.h <- 0
  } else {

    nsubsplits <- length(muhat.vals)
    for(i in 1:nsubsplits) {

      sub.idx <- muhat.vals[[i]]$sub.idx
      muhat.mat <- muhat.vals[[i]]$muhat.mat
      mhat <- mhat.vals[sub.idx]

      int1.h[sub.idx] <-
        colMeans(kern.std.h[sub.idx] * (muhat.mat - mhat))
      int2.h[sub.idx] <-
        colMeans(g2.h[sub.idx] * kern.std.h[sub.idx] *
                   (muhat.mat - mhat))
        int1.b[sub.idx] <- colMeans(kern.std.b[sub.idx] *
                                      (muhat.mat - mhat))

        int2.b[sub.idx] <-
          colMeans(g2.b[sub.idx] * kern.std.b[sub.idx] *
                     (muhat.mat - mhat))
        int3.b[sub.idx] <-
          colMeans(g3.b[sub.idx] * kern.std.b[sub.idx] *
                     (muhat.mat - mhat))
        int4.b[sub.idx] <-
          colMeans(g4.b[sub.idx] * kern.std.b[sub.idx] *
                     (muhat.mat - mhat))
    }
  }
  gamma.h <- coef(lm(Y ~ a.std.h, weights = kern.std.h))
  res.h <- Y - (gamma.h[1] + gamma.h[2] * a.std.h)
  inf.fn <- t(Dh.inv %*% rbind(res.h * kern.std.h + int1.h,
                               g2.h * res.h * kern.std.h + int2.h))

  # Estimate local polynomial components --------------------------------------
  # g2.b <- (A - a)/b
  # g3.b <- ((A - a)/b)^2
  # g4.b <- ((A - a)/b)^3
  # if(all(muhat.mat == 0) & all(muhat.mat == 0)) {
  #   int1.b <- int2.b <- int3.b <- int4.b <- 0
  # } else{
  #   int1.b <- colMeans(kern.std.b * (muhat.mat - mhat.obs))
  #   int2.b <- colMeans(g2.b * kern.std.b * (muhat.mat - mhat.obs))
  #   int3.b <- colMeans(g3.b * kern.std.b * (muhat.mat - mhat.obs))
  #   int4.b <- colMeans(g4.b * kern.std.b * (muhat.mat - mhat.obs))
  # }
    model.b <- lm(Y ~ poly(a.std.b, 3), weights = kern.std.b)
    res.b <- Y - predict(model.b)
    inf.fn.robust <- t(Db.inv %*% rbind(res.b * kern.std.b + int1.b,
                                        g2.b * res.b * kern.std.b + int2.b,
                                        g3.b * res.b * kern.std.b + int3.b,
                                        g4.b * res.b * kern.std.b + int4.b))
    c2 <- integrate(function(u){u^2 * kern(u)}, -Inf, Inf)$value
    debiased_est <- inf.fn[,1] - (h/b)^2 * c2 * inf.fn.robust[,3]
    est <- inf.fn[,1]
    out <- data.frame(debiased_est = debiased_est, est = est)
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
debiased_inference <- function(A, pseudo.out, tau = 1, eval.pts = NULL,
                               muhat.vals = NULL, mhat.obs = NULL, ...){
  # Parse control inputs ------------------------------------------------------
  control <- .parse.debiased_inference(...)
  kernel.type <- control$kernel.type
  # Compute an estimated pseudo-outcome sequence ------------------------------
  ord <- order(A)
  A <- A[ord]
  pseudo.out <- pseudo.out[ord]
  n <- length(A)
  kern <- function(t){.kern(t, kernel=kernel.type)}
  if (is.null(eval.pts)){
    # eval.pts <- seq(quantile(A, 0.05), quantile(A, 0.95),  length.out=30)
    eval.pts <- seq(min(A), max(A),
                    length.out = min(30, length(unique(A))))
  }

  # Compute bandwidth ---------------------------------------------------------
  if(control$bandwidth.method == "LOOCV"){
    if (is.null(control$bw.seq)){
      bw.seq <- seq(0.1*n^{-1/5}, 1*n^{-1/5}, length.out = 10)
    }
    else{
      bw.seq <- control$bw.seq
    }
    bw.seq.h <- rep(bw.seq, length(bw.seq))
    bw.seq.b <- rep(bw.seq, each=length(bw.seq))
    risk <- mapply(function(h, b){
      .robust.loocv(A, pseudo.out, h, b, eval.pt=eval.pts,
                    kernel.type=kernel.type)}, bw.seq.h, bw.seq.b, SIMPLIFY = FALSE)
    risk <- do.call(rbind, risk)
    # risk <- mapply(function(h, b){
    #   .robust.loocv(A, pseudo.out, h, b, debias = FALSE, eval.pt=eval.pts,
    #                 kernel.type=kernel.type)}, bw.seq.h, bw.seq.b)
    h.opt.debias <- bw.seq.h[which.min(risk[, "loocv.risk.debias"])]
    b.opt.debias <- bw.seq.b[which.min(risk[, "loocv.risk.debias"])]
    h.opt <- bw.seq.h[which.min(risk[, "loocv.risk"])]
    b.opt <- bw.seq.b[which.min(risk[, "loocv.risk"])]
  }
  else if(control$bandwidth.method == "LOOCV(h=b)"){
    if (is.null(control$bw.seq)){
      bw.seq <- seq(0.1*n^{-1/5}, 1*n^{-1/5}, length.out = 50)
    }
    else{
      bw.seq <- control$bw.seq
    }
    risk <- mapply(function(h, b){
      .robust.loocv(A, pseudo.out, h, b, debias, eval.pt=eval.pts,
                    kernel.type=kernel.type)}, bw.seq, bw.seq)
    h.opt <- bw.seq[which.min(risk)]
    print(risk)
    b.opt <- bw.seq[which.min(risk)]
  }
  else{
    h.opt <- lpbwselect(pseudo.out, A, eval=eval.pts,
                        bwselect="imse-dpi")$bws[,2]
    b.opt <- h.opt/tau
  }

  # Fit debiased local linear regression --------------------------------------
  est.res <- .lprobust(A, pseudo.out, h.opt, b.opt, eval=eval.pts,
                       kernel.type=kernel.type)
  # if(!is.null(cate.w.fit)) {
  #   muhat.vals <- .get.muhat(splits.id, cate.w.fit, A, v2, max.n.integral)
  # } else muhat.vals <- NULL

  # Estimate influence function sequence --------------------------------------
  rinf.fns <- mapply(function(a, h.val, b.val){
    .compute.rinfl.func(Y = pseudo.out,
                        A = A,
                        a = a,
                        h = h.val,
                        b = b.val,
                        kern = kern,
                        muhat.vals = muhat.vals,
                        mhat.vals = mhat.obs)
  }, eval.pts, h.opt, b.opt, SIMPLIFY = FALSE)
  rif.se <- sapply(rinf.fns, function(u) apply(u, 2, sd)/sqrt(n))
  alpha <- control$alpha.pts
  ci.ll.p.debias <- est.res[,"theta.hat"]-qnorm(1-alpha/2)*rif.se["debiased_est", ]
  ci.ul.p.debias <- est.res[,"theta.hat"]+qnorm(1-alpha/2)*rif.se["debiased_est", ]

  ci.ll.p <- est.res[,"mu.hat"]-qnorm(1-alpha/2)*rif.se["est", ]
  ci.ul.p <- est.res[,"mu.hat"]+qnorm(1-alpha/2)*rif.se["est", ]

  # Compute uniform bands by simulating GP ------------------------------------
  if(control$unif){
    get.unif.ep <- function(alpha, debias){
      var.name <- ifelse(debias, "debiased_est", "est")
      std.inf.vals <- sapply(rinf.fns, function(u) scale(u[, var.name]))
      boot.samples <- control$bootstrap
      ep.maxes<- replicate(boot.samples,
                           max(abs(rbind(rnorm(n)/sqrt(n)) %*% std.inf.vals)))
      quantile(ep.maxes, alpha)
    }
    alpha <- control$alpha.unif
    ep.unif.quant <- get.unif.ep(1-alpha, FALSE) * rif.se["est", ]
    ep.unif.quant.debias <- get.unif.ep(1-alpha, TRUE) * rif.se["debiased_est", ]
    ci.ll.u <- est.res[,"theta.hat"] - ep.unif.quant
    ci.ul.u <- est.res[,"theta.hat"] + ep.unif.quant
    ci.ll.u.debias <- est.res[,"theta.hat"] - ep.unif.quant.debias
    ci.ul.u.debias <- est.res[,"theta.hat"] + ep.unif.quant.debias
  }
  else{
    ci.ll.u <- ci.ul.u <- 0
  }

  # Output data.frame ---------------------------------------------------------
  res <- data.frame(
    eval.pts=eval.pts,
    theta=est.res[, "mu.hat"],
    theta.debias=est.res[,"theta.hat"],
    ci.ul.pts=ci.ul.p,
    ci.ll.pts=ci.ll.p,
    ci.ul.pts.debias=ci.ul.p.debias,
    ci.ll.pts.debias=ci.ll.p.debias,
    ci.ul.unif=ci.ul.u,
    ci.ll.unif=ci.ll.u,
    ci.ul.unif.debias=ci.ul.u.debias,
    ci.ll.unif.debias=ci.ll.u.debias,
    bias=est.res[,"b.hat"],
    if.val=rif.se["est", ],
    if.val.debias=rif.se["debiased_est", ],
    unif.if.val.debias=ep.unif.quant,
    unif.if.val.debias=ep.unif.quant.debias,
    h=h.opt,
    b=b.opt,
    h.debias=h.opt.debias,
    b.bias=b.opt.debias)
  return(res)
}


# Compute hat matrix for leave-one-out cross validation for a linear smoother
.hatmatrix <- function(x, y, h, b, eval.pt=NULL, kernel.type="epa",
                       debias = TRUE){
  ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval.pt)){
    eval.pt <- seq(quantile(x, 0.05), quantile(x, 0.95), length.out=30)
  }
  neval <- length(eval.pt)
  if(debias) {
    c.2 <- integrate(function(t){t^2*.kern(t, kernel.type)}, -Inf, Inf)$value
  } else {
    c.2 <- 0
  }

  e1 <- matrix(c(1, 0), ncol=2); e3 <- matrix(c(1, 0, 0, 0), ncol=4)
  w.h.zero   <- .kern(0, kernel.type)/h
  w.b.zero   <- .kern(0, kernel.type)/b
  hat.mat <- rep(0, neval)
  for (i in 1:neval){
    bw.min   <- sort(abs(x-eval.pt[i]))[10]
    h0     <- max(h, bw.min); b0     <- max(b, bw.min)
    k.h   <- .kern((x-eval.pt[i])/h0, kernel.type)/h0
    k.b   <- .kern((x-eval.pt[i])/b0, kernel.type)/b0
    ind.h <- k.h>0;  ind.b <- k.b>0
    N.h   <- sum(ind.h);  N.b <- sum(ind.b); ind   <- ind.b
    if (h>b) ind <- ind.h
    eY  <- y[ind]; eX  <- x[ind]
    K.h <- k.h[ind]; K.b <- k.b[ind]
    W.b <- matrix(NA, sum(ind), 4) # (1, u, u^2, u^3)
    for (j in 1:4)  W.b[,j] <- ((eX-eval.pt[i])/b0)^(j-1)
    W.h  <- matrix(NA, sum(ind), 2) # (1, u)
    for (j in 1:2)  W.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)
    invD.b  <- .qrXXinv((sqrt(K.b)*W.b))
    invD.h  <- .qrXXinv((sqrt(K.h)*W.h))
    hat.mat[i] <- (invD.h%*%t(e1*w.h.zero))[1,] -
      ((h/b)^2*c.2*invD.b%*%t(e3*w.b.zero))[3,]
  }
  approx(eval.pt, hat.mat, xout = x)$y
}


# Perform leave-one-out cross-validation based on the "hat matrix trick"
.robust.loocv <- function(x, y, h, b, eval.pt=NULL, kernel.type="epa"){
  ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval.pt)){
    eval.pt <- seq(quantile(x, 0.05), quantile(x, 0.95), length.out=30)
  }

  # Compute hat matrix --------------------------------------------------------
  hat.val.debias <- .hatmatrix(x, y, h, b, debias = TRUE, eval.pt=eval.pt,
                               kernel.type=kernel.type)
  hat.val <- .hatmatrix(x, y, h, b, debias = FALSE, eval.pt=eval.pt,
                        kernel.type=kernel.type)
  # Compute debiased local linear estimator -----------------------------------
  est <- .lprobust(x, y, h, b, eval.pt=eval.pt, kernel.type=kernel.type)
  est.fn <- approx(eval.pt, est[,"mu.hat"], xout = x)$y
  est.fn.debias <- approx(eval.pt, est[,"theta.hat"], xout = x)$y

  # Compute the estimate of LOOCV risk ----------------------------------------
  loocv.risk.debias <- mean(((y - est.fn.debias) / (1 - hat.val.debias))^2, na.rm=TRUE)
  loocv.risk <- mean(((y - est.fn) / (1 - hat.val))^2, na.rm=TRUE)

  return(data.frame(loocv.risk.debias = loocv.risk.debias,
                    loocv.risk = loocv.risk))
}


# Debiased local linear regression
.lprobust <- function(x, y, h, b, eval.pt=NULL, kernel.type="epa"){
  ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval.pt)){
    eval.pt <- seq(quantile(x, 0.05), quantile(x, 0.95),  length.out=30)
  }
  neval <- length(eval.pt)
  estimate.mat <- matrix(NA, neval, 4)
  colnames(estimate.mat) <- c("eval", "theta.hat", "mu.hat", "b.hat")
  c.2 <- integrate(function(t){t^2*.kern(t, kernel.type)}, -Inf, Inf)$value
  # if(debias){
  #   c.2 <- integrate(function(t){t^2*.kern(t, kernel.type)}, -Inf, Inf)$value
  # } else {
  #   c.2 <- 0; b <- h
  # }
  for (i in 1:neval){
    bw.min <- sort(abs(x-eval.pt[i]))[10]
    h0 <- max(h, bw.min); b0 <- max(b, bw.min)
    k.h <- .kern((x-eval.pt[i])/h0, kernel.type)/h0
    k.b <- .kern((x-eval.pt[i])/b0, kernel.type)/b0
    ind.h <- k.h>0;  ind.b <- k.b>0; ind <- ind.b
    if (h>b) ind <- ind.h
    eN <- sum(ind); eY <- y[ind]; eX <- x[ind]
    K.h <- k.h[ind]; K.b <- k.b[ind]
    W.b <- matrix(NA, eN, 4)
    for (j in 1:4)  W.b[, j] <- ((eX-eval.pt[i])/b0)^(j-1)
    W.h <- matrix(NA, eN, 2)
    for (j in 1:2)  W.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)
    invD.b  <- .qrXXinv((sqrt(K.b)*W.b))
    invD.h  <- .qrXXinv((sqrt(K.h)*W.h))
    beta.ll  <- invD.h%*%crossprod(W.h*K.h, eY) #R.p^T W.h Y
    beta.bias <- (h/b)^2*invD.b%*%crossprod(W.b*K.b, eY)*c.2
    estimate.mat[i,] <- c(eval.pt[i], beta.ll[1,1]-beta.bias[3,1],
                          beta.ll[1,1], beta.bias[3,1])
  }
  estimate.mat
}

# Local polynomial regression
.locpoly <- function(x, y, h, eval.pt=NULL, kernel.type="epa", degree=1){
  ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval.pt)){
    eval.pt <- seq(quantile(x, 0.05), quantile(x, 0.95),  length.out=30)
  }
  neval <- length(eval.pt)
  estimate.mat <- matrix(NA, neval, 2)
  colnames(estimate.mat) <- c("eval", "theta.hat")
  # Compute for each evaluation point
  for (i in 1:neval){
    bw.min   <- sort(abs(x-eval.pt[i]))[10]
    h0     <- max(h, bw.min)
    k.h   <- .kern((x-eval.pt[i])/h0, kernel=kernel.type)/h0
    ind <- k.h>0
    eY  <- y[ind]; eX  <- x[ind]; K.h <- k.h[ind]
    W.h <- matrix(NA, sum(ind), degree+1)
    for (j in 1:(degree+1))  W.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)
    invD.h <- .qrXXinv((sqrt(K.h)*W.h))
    beta <- invD.h%*%crossprod(W.h*K.h, eY)
    estimate.mat[i,] <- c(eval.pt[i], beta[1,1])
  }
  estimate.mat
}

# Local linear regression
.loclinear <- function(x, y, h, eval.pt=NULL, kernel.type="epa"){
  ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval.pt)){
    eval.pt <- seq(quantile(x, 0.05), quantile(x, 0.95),  length.out=30)
  }
  neval <- length(eval.pt)
  estimate.mat <- matrix(NA, neval, 2)
  colnames(estimate.mat) <- c("eval","theta.hat")
  for (i in 1:neval){
    bw.min   <- sort(abs(x-eval.pt[i]))[10]
    h0 <- max(h, bw.min)
    k.h <- .kern((x-eval.pt[i])/h0, kernel=kernel.type)/h0
    ind <- k.h>0
    eY <- y[ind]; eX <- x[ind]; K.h <- k.h[ind]
    w.h <- matrix(NA, sum(ind), 2)
    for (j in 1:2)  w.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)
    invD.h  <- .qrXXinv((sqrt(K.h)*w.h))
    beta.p  <- invD.h%*%crossprod(w.h*K.h, eY)
    estimate.mat[i,] <- c(eval.pt[i], beta.p[1,1])
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
