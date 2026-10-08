# These functions are slight modifications of those contained in the R package
# for dose-response debiased inference available at
# https://github.com/Kenta426/DebiasedDoseResponse
# Main reference is https://arxiv.org/abs/2210.06448
# Article and repository authored by Kenta Takatsu and Ted Westling.

# Continuous cate() backend. Keep debiased_inference() itself unchanged so
# existing direct callers retain the original influence-function procedure.
.resolve.inference <- function(inference.method, bandwidth.method) {
  inference.method <- match.arg(inference.method, c("lprobust", "influence-function"))
  allowed <- if(inference.method=="lprobust")
    c("imse-dpi", "mse-dpi", "imse-rot", "mse-rot", "ce-dpi", "ce-rot") else
    c("LOOCV", "LOOCV(h=b)")
  if(is.null(bandwidth.method)) bandwidth.method <- allowed[1L]
  if(!is.character(bandwidth.method) || length(bandwidth.method)!=1L ||
     is.na(bandwidth.method) || !bandwidth.method %in% allowed)
    stop("bandwidth.method for inference.method = '", inference.method,
         "' must be one of: ", paste(allowed, collapse=", "))
  list(method=inference.method, bandwidth.method=bandwidth.method)
}

.continuous.inference <- function(A, pseudo.out, eval.pts, inference,
                                  bw.seq=NULL, min.local=NULL,
                                  cv.eval.size=5000L, cv.repeats=5L,
                                  muhat.vals=NULL, mhat.obs=NULL, pd=FALSE) {
  native <- inference$method=="lprobust"
  metadata <- c(inference, list(p=1L, q=if(native) 2L else 3L,
    kernel="gau", uniform.bands=!native,
    pd.integration.if=if(pd) !native else NA,
    variance=if(native) "nn" else "influence-function", min.local=min.local))
  if(!native) {
    run <- function(debias) debiased_inference(A=A, pseudo.out=pseudo.out,
      eval.pts=eval.pts, debias=debias, bandwidth.method=inference$bandwidth.method,
      kernel.type="gau", bw.seq=bw.seq, min.local=min.local,
      cv.eval.size=cv.eval.size, cv.repeats=cv.repeats,
      muhat.vals=muhat.vals, mhat.obs=mhat.obs)
    return(list(regular=run(FALSE), debiased=run(TRUE), metadata=metadata))
  }
  # Native default floor, capped at n for small samples. Explicit min.local
  # overrides the floor; it is not a guarantee of a full-rank local fit.
  bwcheck <- min(length(A), if(is.null(min.local)) 21L else min.local)
  fit <- nprobust::lprobust(y=pseudo.out, x=A, eval=eval.pts, p=1L,
    kernel="gau", bwselect=inference$bandwidth.method, rho=1,
    bwcheck=bwcheck, vce="nn", covgrid=FALSE, level=95)
  est <- fit$Estimate
  convert <- function(debias) {
    theta <- est[,if(debias) "tau.bc" else "tau.us"]
    se <- est[,if(debias) "se.rb" else "se.us"]
    res <- data.frame(eval.pts=est[,"eval"], theta=theta,
      ci.ul.pts=theta+stats::qnorm(.975)*se,
      ci.ll.pts=theta-stats::qnorm(.975)*se,
      ci.ul.unif=NA_real_, ci.ll.unif=NA_real_,
      if.val.sd=se, unif.quantile=NA_real_, h=est[,"h"], b=est[,"b"],
      h.effective=est[,"h"], b.effective=est[,"b"], loocv.risk=NA_real_)
    list(res=res, risk=NULL, cv=NULL, res.list=list(res))
  }
  metadata$bwcheck <- bwcheck
  metadata$rho <- 1
  list(regular=convert(FALSE), debiased=convert(TRUE), metadata=metadata)
}

# Cross-fitted empirical version of the supplement's PD marginalization IF:
# integrate over all modifier observations in the held-out fold, with the
# nuisance regression trained outside that fold. By default keep a lazy source,
# not pairwise predictions. Both axes of each prediction block are bounded.
# stream=FALSE retains the dense reference implementation for small checks.
.get.muhat <- function(splits.id, cate.w.fit, v1, v2, max.n.integral=1000,
                       stream=TRUE) {

  nsplits <- length(cate.w.fit)
  if(length(max.n.integral)!=1L || !is.finite(max.n.integral) ||
     max.n.integral < 1 || max.n.integral != floor(max.n.integral)) {
    stop("max.n.integral must be a positive integer")
  }
  if(length(splits.id)!=length(v1) || nrow(v2)!=length(v1) ||
     anyNA(splits.id) || !all(splits.id %in% seq_len(nsplits))) {
    stop("PD integration requires matching observations and valid fold indices")
  }
  if(any(tabulate(splits.id, nbins=nsplits)==0L)) {
    stop("PD integration requires nonempty folds")
  }
  if(stream) {
    return(structure(list(splits.id=splits.id, fits=cate.w.fit, v1=v1,
                          v2=v2, block.size=max.n.integral),
                     class="pd_integration_source"))
  }
  res <- list()
  counter <- 1

  for(w in seq_len(nsplits)) {
    row.idx <- which(splits.id==w)
    n.k <- length(row.idx)
    if(n.k==0L) stop("PD integration requires nonempty folds")
    blocks <- split(row.idx, ceiling(seq_len(n.k)/max.n.integral))
    for(col.idx in blocks){
      # Rows integrate over all V_i in the held-out fold; columns hold W_j
      # fixed. Column-major vectorization gives M[i,j] = tau(V_i, W_j).
      new.w <- cbind(v1j=rep(v1[row.idx], times=length(col.idx)),
                     v2[rep(col.idx, each=n.k), , drop=FALSE])
      predictions <- cate.w.fit[[w]]$fit(new.w)
      if(!is.numeric(predictions) || length(predictions)!=nrow(new.w) ||
         any(!is.finite(predictions))) {
        stop("PD integration predictions must be finite and match the prediction grid")
      }
      muhat.mat <- matrix(predictions, nrow=n.k, ncol=length(col.idx))
      sub.idx <- seq_along(v1) %in% col.idx
      res[[counter]] <- list(muhat.mat=muhat.mat, sub.idx=sub.idx,
                             row.idx=row.idx)
      counter <- counter + 1
    }
  }
  return(res)
}

# Accumulate all requested weighted integrals in one pass over prediction
# blocks. The output is n by number-of-weights; no n_k by n_k matrix is saved.
.pd.integrate <- function(source, weights, mhat) {
  n <- length(source$v1)
  if(length(mhat)!=n || any(!is.finite(mhat)) || nrow(weights)!=n) {
    stop("PD integration needs finite marginal predictions and aligned weights")
  }
  out <- matrix(0, n, ncol(weights))
  for(k in seq_along(source$fits)) {
    idx <- which(source$splits.id==k)
    blocks <- split(idx, ceiling(seq_along(idx)/source$block.size))
    for(cols in blocks) {
      accum <- matrix(0, length(cols), ncol(weights))
      for(rows in blocks) {
        new.w <- cbind(v1j=rep(source$v1[rows], times=length(cols)),
                       source$v2[rep(cols, each=length(rows)), , drop=FALSE])
        pred <- source$fits[[k]]$fit(new.w)
        if(!is.numeric(pred) || length(pred)!=nrow(new.w) ||
           any(!is.finite(pred))) {
          stop("PD integration predictions must be finite and match the prediction grid")
        }
        centered <- matrix(pred, nrow=length(rows), ncol=length(cols))-mhat[rows]
        accum <- accum + crossprod(centered, weights[rows, , drop=FALSE])
      }
      out[cols, ] <- accum/length(idx)
    }
  }
  out
}

# Batch all evaluation points so nuisance predictions are reused across them.
.pd.moments <- function(source, A, eval.pts, h, b, kern, debias, mhat, min.local=NULL) {
  width <- if(debias) 6L else 2L
  weights <- do.call(cbind, lapply(eval.pts, function(a) {
    hh <- .local.bandwidth(A, a, h, min.local)
    bb <- if(debias) .local.bandwidth(A, a, b, min.local) else b
    if(!is.finite(hh) || !is.finite(bb)) return(matrix(0, length(A), width))
    u <- (A-a)/hh; z <- (A-a)/bb
    ans <- cbind(kern(u)/hh, u*kern(u)/hh)
    if(debias) ans <- cbind(ans, cbind(1,z,z^2,z^3)*kern(z)/bb)
    ans
  }))
  values <- .pd.integrate(source, weights, mhat)
  lapply(seq_along(eval.pts), function(i) {
    structure(list(values=values[, (i-1L)*width+seq_len(width), drop=FALSE]),
              class="pd_integration_moments")
  })
}

.compute.rinfl.func <- function(Y, A, a, h, b, kern, debias, muhat.vals=NULL,
                                mhat.vals=NULL, min.local=NULL){
  #####################################################
  ## ATTN: check with Kenta the definition of if.c2! ##
  #####################################################
  n <- length(A)
  fit <- .local.fit(A, Y, a, h, b, kern, debias, min.local)
  if(is.null(fit)) return(data.frame(est=rep(NA_real_, n)))
  h <- fit$h; b <- fit$b
  a.std.h <- (A-a)/h
  kern.std.h <- kern(a.std.h)/h
  a.std.b <- (A-a)/b
  kern.std.b <- kern(a.std.b)/b
  if(inherits(muhat.vals, "pd_integration_source")) {
    muhat.vals <- .pd.moments(muhat.vals, A, a, h, b, kern, debias, mhat.vals)[[1]]
  }
  if(inherits(muhat.vals, "pd_integration_moments")) {
    vals <- muhat.vals$values
    int1.h <- vals[,1]; int2.h <- vals[,2]
    if(debias) {
      int1.b <- vals[,3]; int2.b <- vals[,4]
      int3.b <- vals[,5]; int4.b <- vals[,6]
    }
  } else if(!is.null(muhat.vals)) {
    int1.h <- int2.h <- int1.b <- int2.b <- int3.b <- int4.b <- rep(NA, n)
    nsubsplits <- length(muhat.vals)
    for(i in 1:nsubsplits) {

      sub.idx <- muhat.vals[[i]]$sub.idx
      row.idx <- muhat.vals[[i]]$row.idx
      # Support callers supplying the older, square-matrix representation.
      if(is.null(row.idx)) row.idx <- which(sub.idx)
      muhat.mat <- muhat.vals[[i]]$muhat.mat
      mhat <- mhat.vals[row.idx]

      int1.h[sub.idx] <- colMeans(kern.std.h[row.idx]*(muhat.mat-mhat))
      int2.h[sub.idx] <- colMeans(a.std.h[row.idx]*kern.std.h[row.idx]*(muhat.mat-mhat))
      if(debias) {
        int1.b[sub.idx] <- colMeans(kern.std.b[row.idx]*(muhat.mat-mhat))
        int2.b[sub.idx] <- colMeans(a.std.b[row.idx]*kern.std.b[row.idx]*(muhat.mat-mhat))
        int3.b[sub.idx] <- colMeans((a.std.b[row.idx])^2*kern.std.b[row.idx]*(muhat.mat-mhat))
        int4.b[sub.idx] <- colMeans((a.std.b[row.idx])^3*kern.std.b[row.idx]*(muhat.mat-mhat))
      }
    }
  } else {
    int1.h <- int2.h <- int1.b <- int2.b <- int3.b <- int4.b <- rep(0, n)
  }

  c0.h <- mean(kern.std.h)
  c1.h <- mean(kern.std.h*a.std.h)
  c2.h <- mean(kern.std.h*a.std.h^2)
  Dh <- matrix(c(c0.h, c1.h,
                 c1.h, c2.h), nrow=2)
  Dh.inv <- n*fit$ih
  gamma.h <- fit$alpha
  res.h <- Y - (gamma.h[1] + gamma.h[2]*a.std.h)
  inf.fn <- t(Dh.inv %*% rbind(res.h*kern.std.h + int1.h,
                               a.std.h*res.h*kern.std.h + int2.h))

  if(debias){
    # Use the same local cubic bias fit as the point estimator.
    X.b <- cbind(1, a.std.b, a.std.b^2, a.std.b^3)
    Db <- crossprod(X.b, X.b*kern.std.b)/n
    Db.inv <- n*fit$ib

    # c2 <- integrate(function(u){u^2 * kern(u)}, -Inf, Inf)$value # old c2
    w.vec <- cbind(kern.std.h, kern.std.h*a.std.h)
    c2.vec <- ((Dh.inv %*% crossprod(w.vec, a.std.h^2))/n)
    c2 <- c2.vec[1]
    # the EIF for the fixed-h c2 -----------------------------------------------
    w.1 <- cbind(1, a.std.h)
    w.1.tilde <- cbind(kern.std.h*a.std.h^2, kern.std.h*a.std.h^3)
    term1 <- (w.1%*%Dh.inv)[,1]*c(w.1 %*% c2.vec)*kern.std.h
    term2 <- (w.1.tilde %*% Dh.inv)[,1]
    beta.b <- matrix(fit$beta,ncol=1)
    deriv2 <- beta.b[3,]
    # if.c2 <- h^2/2*deriv2*(term2-term1)
    if.c2 <- (h/b)^2*deriv2*(term2-term1)

    res.b <- Y - drop(X.b %*% beta.b)
    inf.fn.robust <- t(Db.inv %*% rbind(res.b*kern.std.b + int1.b,
                                        a.std.b*res.b*kern.std.b + int2.b,
                                        a.std.b^2*res.b*kern.std.b + int3.b,
                                        a.std.b^3*res.b*kern.std.b + int4.b))
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
  arg$inference.all <- if(is.null(option$inference.all)) FALSE else option$inference.all
  if(!is.logical(arg$inference.all) || length(arg$inference.all)!=1L ||
     is.na(arg$inference.all)) stop("inference.all must be TRUE or FALSE")
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

# #' plot_debiased_curve
# #' @param res.df ADD
# #' @param ci ADD
# #' @param unif ADD
# #' @export
# plot_debiased_curve <- function(pseudo, exposure, res.df, ci=TRUE, unif=TRUE,
#                                 add.pseudo=TRUE){
#   p <-  ggplot2::ggplot() + ggplot2::xlab("Exposure") +
#     ggplot2::ylab("Covariate-adjusted outcome") +
#     ggplot2::theme_minimal()
#   if(class(res.df$eval.pts) == "factor") {
#     if(add.pseudo) {
#       p <- p + ggplot2::geom_point(data = NULL, ggplot2::aes(x=as.factor(exposure), y=pseudo), col = "gray")
#     }
#     p <- p + ggplot2::geom_point(data = res.df, ggplot2::aes(x=eval.pts, y=theta))
#   } else {
#     if(add.pseudo) {
#       p <- p + ggplot2::geom_point(data = NULL, ggplot2::aes(x=exposure, y=pseudo), col = "gray")
#     }
#     p <- p + ggplot2::geom_line(data = res.df, ggplot2::aes(x=eval.pts, y=theta))
#   }
#   if(ci){
#     p <- p + ggplot2::geom_pointrange(data = res.df,
#                                       ggplot2::aes(x=.data$eval.pts,
#                                                    y=.data$theta,
#                                                    ymin=.data$ci.ll.pts,
#                                                    ymax=.data$ci.ul.pts,
#                                                    size="Pointwise CIs"), col = "black") +
#       ggplot2::scale_size_manual("",values=c("Pointwise CIs"=0.2))
#   }
#   if(unif){
#     p <- p + ggplot2::geom_line(data = res.df, aes(x=eval.pts,
#                                                    y=ci.ll.unif,
#                                                    linetype="Uniform band"), col = "red") +
#       ggplot2::geom_line(data = res.df, aes(x=eval.pts,
#                                             y=ci.ul.unif,
#                                             linetype="Uniform band"),
#                          col = "red")+
#       ggplot2::scale_linetype_manual("",values=c("Uniform band"=2))
#   }
#   p <- p + ggplot2::theme(legend.position = "bottom") +
#     ggplot2::geom_hline(yintercept = 0, col = "blue") +
#     ggplot2::geom_hline(yintercept = mean(pseudo), col = "orange")
#   return(p)
# }

#' debiased_inference
#' @param A ADD
#' @param pseudo.out ADD
#' @param muhat.vals ADD
#' @param mhat.obs ADD
#' @param debias boolean ADD
#' @param tau ratio bandwidths ADD
#' @param eval.pts ADD
#' @param ... Smoothing and inference controls. \code{inference.all=FALSE}
#'   (default) computes intervals only after selecting bandwidths. Set
#'   \code{inference.all=TRUE} for the historical, expensive behavior of
#'   computing intervals for every candidate as well.
#' @param min.local Minimum local observation count, including ties. NULL
#'   (default) disables bandwidth enlargement. Failed local fits return NA with
#'   a warning; failed held-out predictions give infinite candidate risk.
#'   Uniform bands are unavailable if any requested point fails. Effective
#'   bandwidths are returned as h.effective and b.effective. Inference treats
#'   the selected effective bandwidths as fixed.
#' @param cv.eval.size Number of randomly sampled observations at which to
#'   evaluate leave-one-out loss (default 5000); NULL or at least n uses exhaustive
#'   LOOCV. All other observations remain in each deleted regression. Use
#'   set.seed() for reproducibility. This approximates bandwidth selection only.
#' @param cv.repeats Number of independent subsets (default 5). Each subset is
#'   reused across all candidates. Select the coordinatewise median of repeat
#'   winners, projected to the nearest candidate with finite risk in every
#'   repeat using squared log-bandwidth distance; ties use mean risk then larger
#'   bandwidths. Ignored for exhaustive CV. Final estimation uses all rows.
#' @return A list with \code{res} (selected estimates and intervals),
#'   \code{risk} (all candidate LOOCV risks), and \code{res.list} (candidate
#'   estimates; unselected interval fields are \code{NA} by default).
#' @details Bandwidth selection does not evaluate PD nuisance predictions or
#'   simulate confidence bands in the default mode. Selected PD inference
#'   reuses the supplied fitted nuisances; it does not retrain them. Pointwise
#'   results agree with exhaustive inference up to numerical rounding. Uniform
#'   bands may differ because skipping candidate simulations changes RNG use.
#'   Subsampled CV retains full-data estimation. Candidate risks are means over
#'   repeats, whereas bandwidth selection uses median repeat winners and need
#'   not minimize that mean. The returned cv list records indices, candidate
#'   bandwidths, the candidate-by-repeat risk matrix, repeat winners, median
#'   bandwidths and the selected candidate index. For fixed subset size and
#'   repeat count, CV loss evaluation scales linearly in n per candidate;
#'   separate full-fold PD integration may still require quadratic work.
#'   Intervals do not add uncertainty for the randomized selection step.
#' @export
debiased_inference <- function(A, pseudo.out, debias, tau=1, eval.pts=NULL,
                               muhat.vals=NULL, mhat.obs=NULL, ..., min.local=NULL,
                               cv.eval.size=5000L, cv.repeats=5L){
  .validate.min.local(min.local)
  .validate.local.cv(cv.eval.size,cv.repeats)
  # Parse control inputs ------------------------------------------------------
  # control <- .parse.debiased_inference(alpha=0.05, unif=FALSE, kernel.type = "gau",
  #                                      eval.pts = x, A = x, pseudo.out=y, debias=FALSE,
  #                                      muhat.vals = NULL, mahat.obs=NULL, bw.seq = c(0.1, 0.5))
  control <- .parse.debiased_inference(...)
  kernel.type <- control$kernel.type
  if(!kernel.type %in% c("epa","uni","tri","gau")) stop("Unknown kernel.type")
  if(!is.numeric(A) || !length(A) || any(!is.finite(A)) ||
     !is.numeric(pseudo.out) || length(A)!=length(pseudo.out) ||
     any(!is.finite(pseudo.out))) stop("A and pseudo.out must be finite numeric vectors of equal length")
  # Compute an estimated pseudo-outcome sequence ------------------------------
  # ord <- order(A)
  # A <- A[ord]
  # pseudo.out <- pseudo.out[ord]
  n <- length(A)
  kern <- function(t){.kern(t, kernel=kernel.type)}
  if (is.null(eval.pts)){
    eval.pts <- seq(stats::quantile(A, 0.05), stats::quantile(A, 0.95),
                    length.out=min(30, length(unique(A))))
  }
  # Compute bandwidth ---------------------------------------------------------
  bw.seq <- control$bw.seq
  if(!is.numeric(bw.seq) || !length(bw.seq) || any(!is.finite(bw.seq)) || any(bw.seq<=0))
    stop("bw.seq must contain positive finite bandwidths")
  if(!length(eval.pts) || any(!is.finite(eval.pts))) stop("eval.pts must be finite and nonempty")
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


  cv <- .local.cv(A,pseudo.out,bw.seq.h,bw.seq.b,debias,kernel.type,
                  min.local,cv.eval.size,cv.repeats)
  est.proc <- function(h, b, inference=TRUE, estimate=NULL, risk=NULL){
    h.effective <- vapply(eval.pts, function(a) .local.bandwidth(A,a,h,min.local), numeric(1))
    b.effective <- if(debias) vapply(eval.pts, function(a) .local.bandwidth(A,a,b,min.local), numeric(1)) else rep(b,length(eval.pts))
    if(is.null(estimate)) {
      est.res <- .lprobust(A, pseudo.out, h, b, debias, eval.pts,
                          kernel.type, min.local, warn=FALSE)
      loocv.risk <- risk
    } else {
      est.res <- cbind(eval=estimate$eval.pts, theta.hat=estimate$theta)
      loocv.risk <- estimate$loocv.risk[1]
    }
    if(!inference) {
      return(data.frame(eval.pts=eval.pts, theta=est.res[,"theta.hat"],
                        ci.ul.pts=NA_real_, ci.ll.pts=NA_real_,
                        ci.ul.unif=NA_real_, ci.ll.unif=NA_real_,
                        if.val.sd=NA_real_, unif.quantile=NA_real_,
                        h=h, b=b, h.effective=h.effective, b.effective=b.effective,
                        loocv.risk=loocv.risk))
    }

    # Estimate influence function sequence --------------------------------------
    pd.moments <- if(inherits(muhat.vals, "pd_integration_source")) {
      .pd.moments(muhat.vals, A, eval.pts, h, b, kern, debias, mhat.obs, min.local)
    } else NULL
    rinf.fns <- lapply(seq_along(eval.pts), function(i){
      .compute.rinfl.func(Y=pseudo.out, A=A, a=eval.pts[i], h=h, b=b, kern=kern,
                          debias=debias,
                          muhat.vals=if(is.null(pd.moments)) muhat.vals else pd.moments[[i]],
                          mhat.vals=mhat.obs, min.local=min.local)
    })
    rif.se <- matrix(NA, ncol=1, nrow=length(eval.pts), dimnames=list(NULL, "est"))
    rif.se <- do.call(rbind, lapply(rinf.fns, function(u) {
                                                 apply(u, 2, stats::sd)/sqrt(n)
                                               }))
    alpha <- control$alpha.pts
    z.val <- stats::qnorm(1-alpha/2)
    ci.ll.p <- ci.ul.p <- rep(NA, length(eval.pts))
    ci.ll.p <- est.res[,"theta.hat"] - z.val*rif.se[,"est"]
    ci.ul.p <- est.res[,"theta.hat"] + z.val*rif.se[,"est"]

    # Compute uniform bands by simulating GP ------------------------------------
    if(control$unif && all(is.finite(rif.se)) && all(rif.se>0)){
      get.unif.ep <- function(alpha){
        std.inf.vals <- do.call(cbind, lapply(rinf.fns, function(u) scale(u)))
        boot.samples <- control$bootstrap
        ep.maxes<- replicate(boot.samples,
                             max(abs(rbind(stats::rnorm(n)/sqrt(n)) %*% std.inf.vals)))
        stats::quantile(ep.maxes, alpha)
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
      unif.quantile=unname(unif.quantile),
      h=h,
      b=b,
      h.effective=h.effective,
      b.effective=b.effective,
      loocv.risk=loocv.risk)
  }

  res.list <- mapply(est.proc, bw.seq.h, bw.seq.b, risk=cv$risk,
                     MoreArgs=list(inference=control$inference.all), SIMPLIFY=FALSE)
  loocv.vals <- matrix(NA, ncol = 3, nrow = length(res.list),
                       dimnames=list(NULL, c("loocv.risk", "h", "b")))
  for(k in 1:length(res.list)) {
    loocv.vals[k, ] <- c(res.list[[k]]$loocv.risk[1],
                         res.list[[k]]$h[1],
                         res.list[[k]]$b[1])
  }
  if(!any(is.finite(loocv.vals[,1])))
    stop("All bandwidth candidates have failed leave-one-out predictions; increase bw.seq or min.local.")
  selected <- cv$selected
  h.opt <- unname(loocv.vals[selected,2]); b.opt <- unname(loocv.vals[selected,3])
  estimate <- res.list[[selected]]
  res <- est.proc(h=h.opt, b=b.opt, estimate=estimate)
  if(!control$inference.all && length(selected)) {
    for(i in selected) res.list[[i]] <- res
  }
  failed <- sum(!is.finite(res$theta) | !is.finite(res$if.val.sd))
  if(failed) warning(failed, " evaluation point(s) could not be estimated; affected results are NA.",
                     if(control$unif) " Uniform bands are unavailable over the requested grid." else "",
                     call.=FALSE)
  return(list(res=res, risk=as.data.frame(loocv.vals), res.list=res.list,
              cv=cv$diagnostics))
}


.lprobust <- function(x,y,h,b,debias,eval.pt=NULL,kernel.type="epa",
                      min.local=NULL,warn=TRUE) {
  .validate.min.local(min.local)
  if(is.null(eval.pt)) eval.pt <- seq(stats::quantile(x,.05),
                                     stats::quantile(x,.95),length.out=30)
  kern <- function(u) .kern(u,kernel.type)
  theta <- vapply(eval.pt, function(a) {
    fit <- .local.fit(x,y,a,h,b,kern,debias,min.local)
    if(is.null(fit)) NA_real_ else fit$theta
  }, numeric(1))
  if(warn && anyNA(theta))
    warning(sum(is.na(theta)), " local prediction(s) failed; returning NA.",call.=FALSE)
  cbind(eval=eval.pt,theta.hat=theta)
}

# Generate kernel function
.kern = function(u, kernel="epa"){
  if (kernel=="epa") w <- 0.75*(1-u^2)*(abs(u)<=1)
  if (kernel=="uni") w <- 0.5*(abs(u)<=1)
  if (kernel=="tri") w <- (1-abs(u))*(abs(u)<=1)
  if (kernel=="gau") w <- stats::dnorm(u)
  return(w)
}
