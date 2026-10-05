#' CATE
#'
#' Estimate heterogeneous treatment effects (HTEs) defined as \(E(Y^1 - Y^0 | V = v_0)\).
#'
#' @param data A data frame containing the dataset.
#' @param learner Character vector of learners to use (currently only `"dr"` is implemented).
#' @param x_names Character vector with the names of confounding variables \(X\).
#' @param y_name Character string: outcome variable name \(Y\).
#' @param a_name Character string: treatment variable name \(A\).
#' @param v_names Character vector with the names of the effect modifiers \(V\).
#' @param v0 Matrix of evaluation points; rows are values of \(V\) at which the CATE
#'   \eqn{E(Y^1 - Y^0 | V=v_0)} is estimated.
#'   A data frame is also accepted. Columns must be uniquely named with exactly
#'   the names in \code{v_names}; their order is aligned automatically.
#' @param mu1.x Function \code{function(y, a, x, new.x)} that trains a model for
#'   \eqn{E[Y | A=1, X]} and returns a list with elements \emph{res} (predictions at \code{new.x}),
#'   \emph{model} (fitted model object), and \emph{fit} (predictor function of \code{new.x}).
#' @param mu0.x Function \code{function(y, a, x, new.x)} that trains a model for
#'   \eqn{E[Y | A=0, X]} and returns a list with elements \emph{res}, \emph{model}, and \emph{fit}.
#' @param pi.x Function \code{function(a, x, new.x)} that trains a model for the propensity
#'   \eqn{P(A=1 | X)} and returns a list with elements \emph{res}, \emph{model}, and \emph{fit}.
#' @param drl.v Function \code{function(pseudo, v, new.v)} that regresses the pseudo-outcome on \(V\)
#'   to estimate \eqn{E(Y^1 - Y^0 | V)}, returning \emph{res}, \emph{model}, and \emph{fit}.
#' @param drl.x Function \code{function(pseudo, x, new.x)} that regresses the pseudo-outcome on \(X\)
#'   to estimate \eqn{E(Y^1 - Y^0 | X)}, returning \emph{res}, \emph{model}, and \emph{fit}.
#' @param nsplits Integer; number of splits used for cross-validation (ignored if \code{foldid} is given).
#' @param foldid Optional integer vector of fold assignments.
#' @param univariate_reg Logical; if \code{TRUE}, perform univariate regression of the CATE on each
#'   effect modifier separately. Default \code{FALSE}.
#' @param partial_dependence Logical; if \code{TRUE}, compute partial-dependence estimates. Default \code{FALSE}.
#'   Continuous-modifier inference integrates over every held-out-fold modifier
#'   value using streamed prediction blocks. Pairwise predictions are not retained;
#'   storage grows linearly with sample size for a fixed evaluation grid. Exact
#'   integration still requires quadratic work within folds, but is performed
#'   only at the selected ordinary and debiased bandwidths, not every candidate.
#' @param partially_linear Logical; if \code{TRUE}, compute partially linear approximations via Robinson's transformation.
#'   Default \code{FALSE}.
#' @param additive_approx Logical; if \code{TRUE}, compute an additive approximation to the CATE. Default \code{FALSE}.
#' @param bw.stage2 List of length equal to \code{length(v_names)}; each element is a vector of candidate
#'   bandwidths for second-stage regression used in univariate CATE or partial-dependence estimation. Required when
#'   \code{univariate_reg} or \code{partial_dependence} is \code{TRUE}. Default \code{NULL}.
#' @param sample.split.cond.dens Logical; if \code{TRUE}, use sample-splitting for conditional-density estimation.
#'   Default \code{FALSE}.
#' @param cond.dens List of functions for conditional-density estimation, one per effect modifier \eqn{V_j}.
#'   Each should return an object with a method \code{predict.cond.dens(new.v1, new.v2)}. Default \code{NULL}.
#' @param cate.w List of functions \code{function(tau, w, new.w)} (one per effect modifier) that fit \eqn{E[tau | W]}
#'   and return an object with a \code{fit(new.w)} method. Default \code{NULL}.
#' @param cate.not.j List of functions (one per \(j\)) used in the partially linear Robinson step for \eqn{V_{-j}}. Default \code{NULL}.
#' @param reg.basis.not.j List of basis/estimation helpers for the \eqn{V_{-j}} regression in the Robinson step. Default \code{NULL}.
#' @param pl.dfs List (length \code{length(v_names)}) where each element is a vector of candidate degrees of freedom
#'   for the partially linear approximation. Default \code{NULL}.
#' @param fit.basis.additive Optional function used to fit the additive model (GAM); if \code{NULL} and
#'   \code{additive_approx=TRUE}, \code{drl.additive.basis()} is used. Default \code{NULL}.
#'   Its \code{model} must be an unweighted, full-rank additive \code{lm} with an
#'   intercept, no offset, and all pseudo-outcomes as its response in input row order.
#' @param min.local Minimum number of observations (including ties) within each continuous local bandwidth. NULL disables enlargement (default); a positive integer enlarges bandwidths as needed, also in each deleted sample for LOOCV. Rank-deficient local fits return NA with a warning; candidates with any failed held-out prediction receive infinite risk. This does not guarantee full rank.
#'   For Gaussian kernels this counts observations within one bandwidth, not
#'   all nonzero weights. Effective bandwidths are returned in the continuous
#'   result tables; inference treats these selected bandwidths as fixed.
#'   Uniform bands are unavailable if any requested evaluation point fails.
#' @param density.ratio Optional list, in \code{v_names} order, of functions
#'   \code{function(v1, v2)} returning a list with
#'   \code{predict=function(new.v1, new.v2)}. Predictions must be a finite,
#'   nonnegative numeric vector estimating \eqn{f_j(v)/f_j(v\mid w)}.
#'   A \code{NULL} list entry uses the existing \code{cond.dens} pathway for
#'   that modifier; a supplied function overrides it. Uses training folds when
#'   \code{sample.split.cond.dens=TRUE}, otherwise fits on the whole sample.
#'   Custom ratios are used directly in PD estimation and inference, without
#'   a marginal-density calculation, density floor, or automatic weight cap.
#'   Apply any desired stabilization in the supplied predictor. The final weights
#'   are returned in \code{pd.res$dr[[j]]$data$density.ratio}; conditional-density
#'   values, when present in that data frame, are \code{NA} for this pathway.
#'   Ignored unless \code{partial_dependence=TRUE}.
#'
#' @details
#' For discrete modifiers, \code{pd.res$dr[[j]]$res} and \code{res.empVar}
#' both report the cross-fitted one-step PD estimate. Its score is the
#' inverse-conditional-probability weighted residual plus the fitted effect
#' at the requested level and each observation's other modifiers. The mean
#' score estimates PD; its sample standard deviation divided by \eqn{\sqrt n}
#' gives the standard error, including empirical-reference-distribution uncertainty.
#' These estimates may differ from the previous subgroup regression of PD
#' pseudo-outcomes; \code{data$pseudo} remains a diagnostic pseudo-outcome.
#' Discrete univariate \code{res} and \code{res.empVar} both report subgroup
#' means with subgroup-centered influence values, including for numeric
#' multi-category modifiers (no linear trend is imposed). Intervals are pointwise 95 percent;
#' validity requires the usual positivity and nuisance-rate conditions.
#'
#' For each modifier, \code{additive.res$dr[[j]]$res} and
#' \code{robinson.res$dr[[j]]$res} report treatment-effect levels:
#' \eqn{\widehat\psi + \{b_j(v)-\overline b_j\}^T\widehat\beta_j}, where
#' \eqn{\widehat\psi} is the mean cross-fitted DR pseudo-outcome. These are
#' model-based partial-dependence profiles, not generally the univariate CATE
#' when modifiers are dependent. Under model misspecification they are working
#' approximation profiles, not necessarily the true partial dependence.
#'
#' The \code{se} column is a standard error (not a variance). Pointwise 95 percent
#' intervals use the joint influence function of the ATE, empirical basis mean,
#' and regression coefficients, including their covariances. For Robinson,
#' write \eqn{R_i=b_j(V_{ij})-\widehat E[b_j(V_j)\mid V_{-j}]_i},
#' \eqn{e_i=\widehat\phi_i-\widehat E[\phi\mid V_{-j}]_i-R_i^T\widehat\beta_j},
#' and \eqn{\widehat Q=n^{-1}\sum_i R_iR_i^T}. The estimated influence value is
#' \deqn{\widehat\phi_i-\widehat\psi
#' -\{b_j(V_{ij})-\overline b_j\}^T\widehat\beta_j
#' +\{b_j(v)-\overline b_j\}^T\widehat Q^{-1}R_i e_i.}
#' For the additive fit, coefficient influence values come from the full joint
#' regression, including the intercept and all other components. Standard errors
#' are the sample standard deviation of influence values divided by \eqn{\sqrt n}.
#' Inference assumes a fixed basis and the usual nuisance-rate and regularity
#' conditions; it does not account for selecting the basis dimension. Uniform
#' bands are not supplied for these profiles.
#'
#' This changes the primary approximation curves from previous versions.
#' Each modifier also has \code{component.res} (a sample-centered component,
#' with centering uncertainty) and \code{ate} (the baseline added to the
#' centered component). Only additive results retain \code{legacy.res}, the
#' previous component output and intervals. Robinson's zero-anchored output
#' is no longer returned. Do not add another
#' baseline to \code{res}. Existing model and tuning-result fields are retained.
#'
#' @return A list with elements for V-based and X-based CATE results, per-fold estimates, pseudo-outcomes,
#'   univariate/partial-dependence/additive/Robinson outputs, and inputs used.
#'
#' @export
#' @references Kennedy, E. H. (2020). Optimal Doubly Robust Estimation of Heterogeneous Causal Effects.
#'   \emph{arXiv preprint} arXiv:2004.14497.

cate <- function(data, learner, x_names, y_name, a_name, v_names, v0,
                 mu1.x, mu0.x, pi.x, drl.v, drl.x,
                 nsplits=5,
                 foldid=NULL,
                 univariate_reg=FALSE,
                 partial_dependence=FALSE,
                 partially_linear=FALSE,
                 additive_approx=FALSE,
                 bw.stage2=NULL,
                 sample.split.cond.dens=FALSE,
                 cond.dens=NULL,
                 cate.w=NULL,
                 cate.not.j=NULL,
                 reg.basis.not.j=NULL,
                 pl.dfs=NULL,
                 fit.basis.additive=NULL,
                 density.ratio=NULL, min.local=NULL) {
  .validate.min.local(min.local)

  if(any(learner != "dr")) stop("Only learner = dr is currently implemented.")

  dta <- get_input(data=data, x_names=x_names, y_name=y_name,
                   a_name=a_name, v_names=v_names, v0=v0)

  a <- dta$a
  v <- dta$v
  v0.long <- dta$v0
  v0.short <- dta$unique.v0
  y <- dta$y
  x <- dta$x

  n <- length(y)
  n.eval.pts <- nrow(v0.long)
  n.eff.modif <- ncol(v)
  if (is.null(density.ratio)) density.ratio <- rep(list(NULL), n.eff.modif)
  if (!is.list(density.ratio) || length(density.ratio) != n.eff.modif ||
      !all(vapply(density.ratio, function(f) is.null(f) || is.function(f), logical(1)))) {
    stop("density.ratio must be NULL or a list with one function or NULL per effect modifier.")
  }
  # Fit on the same rows as the conditional-density nuisance. Custom ratios
  # are already final weights: no marginal multiplication or density floor.
  predict.ratio <- function(j, v1, v2, new.v1, new.v2) {
    tryCatch({
      fit <- density.ratio[[j]](v1=v1, v2=v2)
      out <- fit$predict(new.v1=new.v1, new.v2=new.v2)
      if (!is.numeric(out) || !is.null(dim(out)) || length(out) != length(new.v1) ||
          any(!is.finite(out)) || any(out < 0)) {
        stop("density.ratio predictions must be a finite, nonnegative numeric vector with one value per row.")
      }
      out
    }, error=function(e) stop(errorCondition(conditionMessage(e),
                                             class="cate_density_ratio_error")))
  }
  ratio.vals <- matrix(NA_real_, nrow=n, ncol=n.eff.modif)

  if(is.null(foldid)) {
    s <- sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  } else {
    s <- foldid
    nsplits <- length(unique(foldid))
  }

  est <- est.pi <- replicate(length(learner),
                             array(NA, dim=c(n.eval.pts, 3, nsplits)),
                             simplify=FALSE)
  univariate_res <- pd_res <- additive_res <- robinson_res <-
    replicate(length(learner), vector("list", n.eff.modif), simplify=FALSE)

  cate.w.fit <- stage2.reg.data.pd <-
    replicate(n.eff.modif, vector("list", nsplits), simplify=FALSE)

  pseudo.y <- replicate(length(learner), rep(NA, n), simplify=FALSE)
  pseudo.y.tr <- ites.x.tr <-
    replicate(length(learner), vector("list", nsplits), simplify=FALSE)

  pseudo.y.pd <- theta.bar <- cond.dens.vals <- cate.w.vals <-
    replicate(length(learner), matrix(NA, ncol=n.eff.modif, nrow=n),
              simplify=FALSE)

  ites_v <- ites_x <- replicate(length(learner),
                                matrix(NA, ncol=3, nrow=n), simplify=FALSE)

  names(est) <- names(est.pi) <- names(pseudo.y) <- names(ites_v) <-
    names(ites_x) <- names(pseudo.y.pd) <- names(theta.bar) <-
    names(cond.dens.vals) <- names(cate.w.vals) <- names(pseudo.y.tr) <- learner

  stage2.reg.data.v <- stage2.reg.data.x <- reg.model <-
    vector("list", nsplits)

  tmp <- tryCatch(
    {
      for(k in 1:nsplits) {

        print(paste0("Considering split # ", k, " out of ", nsplits))

        test.idx <- k == s
        train.idx <- k != s
        if(all(!train.idx)) train.idx <- test.idx
        n.te <- sum(test.idx)
        n.tr <- sum(train.idx)

        x.tr <- x[train.idx, , drop=FALSE]
        v.tr <- v[train.idx, , drop=FALSE]
        a.tr <- a[train.idx]
        y.tr <- y[train.idx]

        x.te <- x[test.idx, , drop=FALSE]
        v.te <- v[test.idx, , drop=FALSE]
        a.te <- a[test.idx]
        y.te <- y[test.idx]

        ## Estimate nuisance functions using all folds but k and predict on fold k ##
        pihat.vals <- pi.x(a=a.tr, x=x.tr, new.x=rbind(x.te, x.tr))$res
        pihat.te <- pihat.vals[1:n.te]
        pihat.tr <- pihat.vals[-c(1:n.te)]

        mu0hat.vals <- mu0.x(y=y.tr, a=a.tr, x=x.tr, new.x=rbind(x.te, x.tr))$res
        mu0hat.te <- mu0hat.vals[1:n.te]
        mu0hat.tr <-  mu0hat.vals[-c(1:n.te)]

        mu1hat.vals <- mu1.x(y=y.tr, a=a.tr, x=x.tr, new.x=rbind(x.te, x.tr))$res
        mu1hat.te <- mu1hat.vals[1:n.te]
        mu1hat.tr <-  mu1hat.vals[-c(1:n.te)]

        cate.tr <- mu1hat.tr-mu0hat.tr
        cate.te <- mu1hat.te-mu0hat.te

        for(alg in learner) {
          ## compute IF values, i.e., the pseudo-outcomes for the DR-Learner ##
          pseudo.te <- (a.te-pihat.te)/(pihat.te*(1-pihat.te)) *
            (y.te-a.te*mu1hat.te - (1-a.te)*mu0hat.te) + cate.te

          pseudo.tr <- (a.tr-pihat.tr)/(pihat.tr*(1-pihat.tr)) *
            (y.tr-a.tr*mu1hat.tr - (1-a.tr)*mu0hat.tr) + cate.tr

          drl.v.out <-  drl.v(pseudo=pseudo.te, v=v.te, new.v=rbind(v0.long, v.te))
          drl.v.out.pi <-  drl.v(pseudo=cate.te, v=v.te, new.v=rbind(v0.long, v.te)) # plug-in
          stage2.reg.data.v[[k]] <- cbind(data.frame(pseudo=pseudo.te,
                                                   mu1hat=mu1hat.te,
                                                   mu0hat=mu0hat.te,
                                                   pihat=pihat.te,
                                                   y=y.te,
                                                   a=a.te,
                                                   fold.id=k), v.te)
          stage2.reg.data.x[[k]] <- cbind(data.frame(pseudo=pseudo.te,
                                                     mu1hat=mu1hat.te,
                                                     mu0hat=mu0hat.te,
                                                     pihat=pihat.te,
                                                     y=y.te,
                                                     a=a.te,
                                                     fold.id=k), x.te)

          # drl.form[[k]] <- drl.res$drl.form
          reg.model[[k]] <- drl.v.out$model

          drl.v.res <- drl.v.out$res
          drl.v.res.pi <- drl.v.out.pi$res
          drl.x.res <- drl.x(pseudo=pseudo.tr, x=x.tr, new.x=rbind(x.te, x.tr))$res

          est[[alg]][, , k] <- drl.v.res[1:n.eval.pts, ]
          est.pi[[alg]][, , k] <- drl.v.res.pi[1:n.eval.pts, ]
          pseudo.y[[alg]][test.idx] <- pseudo.te
          pseudo.y.tr[[alg]][[k]] <- pseudo.tr
          ites.x.tr[[alg]][[k]] <- drl.x.res[-c(1:n.te), 1]
          ites_v[[alg]][test.idx, ] <- drl.v.res[-c(1:n.eval.pts), ]
          ites_x[[alg]][test.idx, ] <- drl.x.res[1:n.te, ]


          if(partial_dependence) {

            for(j in 1:ncol(v)) {

              v1.j.tr <- v.tr[, j]
              v1.j.te <- v.te[, j]
              not.v1.j.tr <- v.tr[, -j, drop=FALSE]
              not.v1.j.te <- v.te[, -j, drop=FALSE]

              w.tr <- cbind(v1j=v1.j.tr, not.v1.j.tr)
              w.te <-  cbind(v1j=v1.j.te, not.v1.j.te)

              if(sample.split.cond.dens && !is.null(density.ratio[[j]])) {
                ratio.vals[test.idx, j] <- predict.ratio(j, v1.j.tr, not.v1.j.tr,
                                                        v1.j.te, not.v1.j.te)
              } else if(sample.split.cond.dens){
                cond.dens.fit <- cond.dens[[j]](v1=v1.j.tr, v2=not.v1.j.tr)
                cond.dens.vals.te <-
                  cond.dens.fit$predict.cond.dens(v1=v1.j.tr, v2=not.v1.j.tr,
                                                  new.v1=v1.j.te, new.v2=not.v1.j.te)
                if(sum(cond.dens.vals.te < 0.001) > 0) {
                  warning(paste0("Effect modifier # ", j, ". There are ",
                                 sum(cond.dens.vals.te < 0.001),
                                 " conditional density values < 0.001. They will ",
                                 "truncated at 0.001."))
                  cond.dens.vals.te[cond.dens.vals.te < 0.001] <- 0.001
                }
                cond.dens.vals[[alg]][test.idx, j] <- cond.dens.vals.te
              }

              cate.w.fit[[j]][[k]] <- cate.w[[j]](tau=cate.tr, w=w.tr, new.w=w.tr)

              cate.w.te <- cate.w.fit[[j]][[k]]$fit(new.w=w.te)
              cate.w.vals[[alg]][test.idx, j] <- cate.w.te

              if(n.te > 1000) {

                if(is.factor(v1.j.te)) {
                  v1.j.seq <- factor(levels(v1.j.te), levels=levels(v1.j.te))
                }
                else {
                  v1.j.seq <- seq(min(v1.j.te), max(v1.j.te), length.out=100)
                }

                tmp.cate.w.fit.fn <- Vectorize(function(u) {
                  mean(cate.w.fit[[j]][[k]]$fit(new.w=cbind(v1j=u, not.v1.j.te)))
                }, vectorize.args = "u")

                cate.w.avg.vals <- tmp.cate.w.fit.fn(v1.j.seq)

                if(is.factor(v1.j.te)) {
                  theta.bar.vals <- rep(NA, length(v1.j.te))
                  for(u in levels(v1.j.te)) {
                    theta.bar.vals[v1.j.te==u] <- cate.w.avg.vals[v1.j.seq==u]
                  }
                } else {
                  theta.bar.vals <- stats::approx(x=v1.j.seq, y=cate.w.avg.vals,
                                                  xout=v1.j.te, rule=2)$y
                }
              }
              else {
                w.long.test <- cbind(v1j=rep(v1.j.te, each=n.te),
                                     not.v1.j.te[rep(1:n.te, n.te), , drop=FALSE])
                if(sample.split.cond.dens && is.null(density.ratio[[j]])) {
                  cond.dens.preds <-
                    cond.dens.fit$predict.cond.dens(v1=v1.j.tr, v2=not.v1.j.tr,
                                                    new.v1=w.long.test[, 1],
                                                    new.v2=w.long.test[, -1, drop=FALSE])

                  marg.dens <- colMeans(matrix(cond.dens.preds, ncol=n.te, nrow=n.te))
                }
                cate.preds <- cate.w.fit[[j]][[k]]$fit(new.w=w.long.test)
                theta.bar.vals <- colMeans(matrix(cate.preds, nrow=n.te, ncol=n.te))
              }

              theta.bar[[alg]][test.idx, j] <- theta.bar.vals
              data.pd <- data.frame(pseudo.cate=pseudo.te,
                                    mu1hat=mu1hat.te,
                                    mu0hat=mu0hat.te,
                                    pihat=pihat.te,
                                    tauhat.w=cate.w.te,
                                    theta.bar=theta.bar.vals,
                                    y=y.te,
                                    a=a.te,
                                    fold.id=k)
              stage2.reg.data.pd[[j]][[k]] <- cbind(data.pd, v.te)
            }
          }
        }
      }
      print("Done with fitting nuisance functions.")
    },
    error = function(cond) {
      if (inherits(cond, "cate_density_ratio_error")) stop(cond)
      message(conditionMessage(cond))
    }
  )
  if(is.null(tmp)){
    warning("Encountered error while fitting nuisance functions")
    univ.res <- pd.res <- add.res <- rob.res <-
      data.frame(theta=rep(NA, 10), theta.debias=rep(NA, 10))
    return(list(univariate_res=list(dr=list(list(res=univ.res))),
                pd_res=list(dr=list(list(res=pd.res))),
                additive_res=list(dr=list(list(res=add.res))),
                robinson_res=list(dr=list(list(res=rob.res)))))
  }
  for(alg in learner) {

    if(alg != "dr") stop("Only learner = dr is currently implemented.")

    # if(additive_approx) {
    #   additive_model <- drl.basis.additive(y=pseudo.y[[alg]], x=v, new.x=v)
    #   tt <- delete.response(terms(additive_model$model))
    # }
    if(univariate_reg | partial_dependence | additive_approx | partially_linear) {

      for(j in 1:ncol(v)){
        vj <- v[, j]
        is.var.factor <- paste0(class(vj), collapse=" ") %in% c("factor", "ordered factor")

        if(partially_linear) {

          j.robinson <- robinson(pseudo=pseudo.y[[alg]],
                                 w=v[, -j, drop=FALSE],
                                 v=v[, j],
                                 new.v=v0.short[[j]],
                                 s=s,
                                 cate.not.j=cate.not.j[[j]],
                                 reg.basis.not.j=reg.basis.not.j[[j]],
                                 dfs=pl.dfs[[j]])

          rob.level <- robinson_effect_profile(j.robinson$model, pseudo.y[[alg]],
                                                vj, v0.short[[j]])
          robinson_res[[alg]][[j]] <- list(res=rob.level$res,
                                           component.res=rob.level$component.res,
                                           ate=rob.level$ate,
                                           model=j.robinson$model,
                                           risk=j.robinson$risk,
                                           fits=j.robinson$fits)

        }

        if(partial_dependence) {

          if (!is.null(density.ratio[[j]])) {
            if (!sample.split.cond.dens) {
              ratio.vals[, j] <- predict.ratio(j, vj, v[, -j, drop=FALSE],
                                              vj, v[, -j, drop=FALSE])
            }
            ghat <- ratio.vals[, j]
          } else {
            if(!sample.split.cond.dens) {
              cond.dens.fit <- cond.dens[[j]](v1=v[, j], v2=v[, -j, drop=FALSE])
              cond.dens.vals[[alg]][, j] <-
                cond.dens.fit$predict.cond.dens(v1=v[, j], v2=v[, -j, drop=FALSE],
                                                new.v1=v[, j],
                                                new.v2=v[, -j, drop=FALSE])

              if(sum(cond.dens.vals[[alg]][, j] < 0.001) > 0) {
                warning(paste0("Effect modifier # ", j, ". There are ",
                               sum(cond.dens.vals[[alg]][, j] < 0.001),
                               " conditional density values < 0.001. They will ",
                               "truncated at 0.001."))
                cond.dens.vals[[alg]][cond.dens.vals[[alg]][, j] < 0.001, j] <- 0.001
              }
            }
            if(length(unique(vj)) < 15 | is.var.factor) {
              marg.dens <- rep(NA, length(vj))
              for(u in unique(vj)) marg.dens[vj==u] <- mean(vj==u)
            } else {
              marg.dens <- local({
                h <- tryCatch(bw.SJ(vj), error = function(e) bw.nrd0(vj))
                vapply(vj, function(z) mean(dnorm((z - vj) / h) / h), numeric(1))
              })
            }
            ghat <- marg.dens/cond.dens.vals[[alg]][, j]
          }
          pseudo.y.pd[[alg]][, j] <-
            (pseudo.y[[alg]]-cate.w.vals[[alg]][, j])*ghat + theta.bar[[alg]][, j]
        }
        if(length(unique(vj)) < 15 | is.var.factor) {
          if(univariate_reg) {
            res.empVar <- NULL
            for(ll in 1:length(unique(vj))) {
              pts.vj <- unique(vj)[ll]
              theta <- mean(pseudo.y[[alg]][vj==pts.vj])
              if.vals <- (pseudo.y[[alg]]-theta)*(vj==pts.vj) / mean(vj==pts.vj)
              tmp <- data.frame(eval.pts=pts.vj,
                                    theta=theta,
                                    ci.ll.pts=theta - stats::qnorm(.975)*sqrt(stats::var(if.vals)/n),
                                    ci.ul.pts=theta + stats::qnorm(.975)*sqrt(stats::var(if.vals)/n),
                                    ci.ul.unif=NA,
                                    ci.ll.unif=NA)
              res.empVar <- rbind(res.empVar, tmp)
            }
            univariate_res[[alg]][[j]] <-
              list(data=data.frame(pseudo=pseudo.y[[alg]], exposure=vj),
                   res=res.empVar,
                   res.empVar=res.empVar)
          }
          if(partial_dependence) {
            res.empVar <- NULL
            for(ll in 1:length(unique(vj))) {
              pts.vj <- unique(vj)[ll]
              correction.weight <- if (is.null(density.ratio[[j]])) {
                (vj==pts.vj)/cond.dens.vals[[alg]][, j]
              } else {
                (vj==pts.vj)*ghat/mean(vj==pts.vj)
              }
              # Each patient's W contributes to the empirical PD reference
              # distribution, regardless of that patient's observed V_j.
              tau.at.level <- rep(NA_real_, n)
              for(k in seq_len(nsplits)) {
                idx <- which(s==k)
                new.w <- cbind(v1j=rep(pts.vj, length(idx)),
                               v[idx, -j, drop=FALSE])
                predictions <- cate.w.fit[[j]][[k]]$fit(new.w)
                if(!is.numeric(predictions) || length(predictions)!=length(idx) ||
                   any(!is.finite(predictions))) {
                  stop("Discrete PD predictions must be finite and match the held-out rows")
                }
                tau.at.level[idx] <- predictions
              }
              score <- (pseudo.y[[alg]]-cate.w.vals[[alg]][, j]) * correction.weight +
                tau.at.level
              theta <- mean(score)
              if.vals <- score-theta
              tmp <- data.frame(eval.pts=pts.vj,
                                     theta=theta,
                                     ci.ll.pts=theta - stats::qnorm(.975)*sqrt(stats::var(if.vals)/n),
                                     ci.ul.pts=theta + stats::qnorm(.975)*sqrt(stats::var(if.vals)/n),
                                     ci.ul.unif=NA,
                                     ci.ll.unif=NA)
              res.empVar <- rbind(res.empVar, tmp)
            }
            pd_res[[alg]][[j]] <-
              list(data=data.frame(pseudo=pseudo.y.pd[[alg]][, j], exposure=vj),
                   stage2.reg.data.pd=stage2.reg.data.pd,
                   res=res.empVar,
                   res.empVar=res.empVar)
          }
        }
        else {
          if(univariate_reg) {
            univ.inf <- debiased_inference(A=vj,
                                           pseudo.out=pseudo.y[[alg]],
                                           eval.pts=v0.short[[j]],
                                           debias=FALSE,
                                           bandwidth.method="LOOCV",
                                           kernel.type="gau",
                                           bw.seq=bw.stage2[[j]], min.local=min.local)
            univ.debias.inf <- debiased_inference(A=vj,
                                                  pseudo.out=pseudo.y[[alg]],
                                                  eval.pts=v0.short[[j]],
                                                  debias=TRUE,
                                                  bandwidth.method="LOOCV",
                                                  kernel.type="gau",
                                                  bw.seq=bw.stage2[[j]], min.local=min.local)
            univ.res <- data.frame(eval.pts=univ.inf$res$eval.pts,
                                   theta=univ.inf$res$theta,
                                   theta.debias=univ.debias.inf$res$theta,
                                   ci.ul.pts=univ.inf$res$ci.ul.pts,
                                   ci.ll.pts=univ.inf$res$ci.ll.pts,
                                   ci.ul.pts.debias=univ.debias.inf$res$ci.ul.pts,
                                   ci.ll.pts.debias=univ.debias.inf$res$ci.ll.pts,
                                   ci.ul.unif=univ.inf$res$ci.ul.unif,
                                   ci.ll.unif=univ.inf$res$ci.ll.unif,
                                   ci.ul.unif.debias=univ.debias.inf$res$ci.ul.unif,
                                   ci.ll.unif.debias=univ.debias.inf$res$ci.ll.unif,
                                   bias=univ.inf$res$theta-univ.debias.inf$res$theta,
                                   if.val.sd=univ.inf$res$if.val.sd,
                                   if.val.sd.debias=univ.debias.inf$res$if.val.sd,
                                   unif.quantile=univ.inf$res$unif.quantile,
                                   unif.quantile.debias=univ.debias.inf$res$unif.quantile,
                                   h=univ.inf$res$h,
                                   b=univ.inf$res$b,
                                   h.debias=univ.debias.inf$res$h,
                                   b.debias=univ.debias.inf$res$b,
                                   h.effective=univ.inf$res$h.effective,
                                   b.effective=univ.inf$res$b.effective,
                                   h.debias.effective=univ.debias.inf$res$h.effective,
                                   b.debias.effective=univ.debias.inf$res$b.effective)

            univariate_res[[alg]][[j]] <-
              list(data=data.frame(pseudo=pseudo.y[[alg]], exposure=vj),
                   res=univ.res,
                   risk=list(risk=univ.inf$risk, risk.debias=univ.debias.inf$risk),
                   res.list=univ.inf$res.list,
                   res.list.debias=univ.debias.inf$res.list)
          }

          if(partial_dependence) {

            muhat.vals <- .get.muhat(splits.id=s, cate.w.fit=cate.w.fit[[j]],
                                     v1=vj, v2=v[, -j, drop=FALSE],
                                     max.n.integral=1000)

            pd.inf <- debiased_inference(A=vj, debias=FALSE,
                                         pseudo.out=pseudo.y.pd[[alg]][, j],
                                         eval.pts=v0.short[[j]],
                                         mhat.obs=theta.bar[[alg]][, j],
                                         muhat.vals=muhat.vals,
                                         bandwidth.method="LOOCV",
                                         kernel.type="gau",
                                         bw.seq=bw.stage2[[j]], min.local=min.local)

            pd.debias.inf <- debiased_inference(A=vj, debias=TRUE,
                                                pseudo.out=pseudo.y.pd[[alg]][, j],
                                                eval.pts=v0.short[[j]],
                                                mhat.obs=theta.bar[[alg]][, j],
                                                muhat.vals=muhat.vals,
                                                bandwidth.method="LOOCV",
                                                kernel.type="gau",
                                                bw.seq=bw.stage2[[j]], min.local=min.local)

            pd.inf.res <- data.frame(eval.pts=pd.inf$res$eval.pts,
                                     theta=pd.inf$res$theta,
                                     theta.debias=pd.debias.inf$res$theta,
                                     ci.ul.pts=pd.inf$res$ci.ul.pts,
                                     ci.ll.pts=pd.inf$res$ci.ll.pts,
                                     ci.ul.pts.debias=pd.debias.inf$res$ci.ul.pts,
                                     ci.ll.pts.debias=pd.debias.inf$res$ci.ll.pts,
                                     ci.ul.unif=pd.inf$res$ci.ul.unif,
                                     ci.ll.unif=pd.inf$res$ci.ll.unif,
                                     ci.ul.unif.debias=pd.debias.inf$res$ci.ul.unif,
                                     ci.ll.unif.debias=pd.debias.inf$res$ci.ll.unif,
                                     bias=pd.inf$res$theta-pd.debias.inf$res$theta,
                                     if.val.sd=pd.inf$res$if.val.sd,
                                     if.val.sd.debias=pd.debias.inf$res$if.val.sd,
                                     unif.quantile=pd.inf$res$unif.quantile,
                                     unif.quantile.debias=pd.debias.inf$res$unif.quantile,
                                     h=pd.inf$res$h,
                                     b=pd.inf$res$b,
                                     h.debias=pd.debias.inf$res$h,
                                     b.debias=pd.debias.inf$res$b,
                                   h.effective=pd.inf$res$h.effective,
                                   b.effective=pd.inf$res$b.effective,
                                   h.debias.effective=pd.debias.inf$res$h.effective,
                                   b.debias.effective=pd.debias.inf$res$b.effective)

            pd_res[[alg]][[j]] <-
              list(data=data.frame(pseudo=pseudo.y.pd[[alg]][, j],
                                   cond.dens.vals=cond.dens.vals[[alg]][, j],
                                   exposure=vj),
                   res=pd.inf.res,
                   risk=list(risk=pd.inf$risk, risk.debias=pd.debias.inf$risk),
                   res.list=pd.inf$res.list,
                   res.list.debias=pd.debias.inf$res.list)
          }
        }

        if (partial_dependence && !is.null(density.ratio[[j]])) {
          pd_res[[alg]][[j]]$data$density.ratio <- ghat
        }

        if(additive_approx){
          if(is.null(fit.basis.additive)) fit.basis.additive <- drl.basis.additive
          additive_model <- fit.basis.additive(y=pseudo.y[[alg]], x=v, new.x=v)

          res.gam <- get.smooth.fit.gam(fit=additive_model$model,
                                        eval.pts=v0.short[[j]],
                                        eff.modif.name=colnames(v)[j], v=v)

          # new.dat.additive <- as.data.frame(matrix(0, nrow=length(v0.short[[j]]),
          #                                          ncol=ncol(v),
          #                                          dimnames=list(NULL, colnames(v))))
          # for(l in 1:ncol(v)) {
          #   if(l == j) {
          #     new.dat.additive[, l] <- v0.short[[j]]
          #   } else {
          #     if(is.factor(v[, l])) {
          #       new.dat.additive[, l] <- factor(levels(v[, l])[1], levels=levels(v[, l]))
          #     } else {
          #       new.dat.additive[, l] <- min(v[, l])
          #     }
          #   }
          # }
          # form <- formula(additive_model$model)
          # # term_labels <- attr(terms(form), "term.labels")
          # mm <- model.matrix(additive_model$model)
          # coefs <- coef(additive_model$model)
          #
          # bs_term_for_vj <- grep(colnames(v)[j], colnames(mm), value=TRUE)
          # coefs.names.vj <- grep(colnames(v)[j], names(coefs), value=TRUE)
          #
          # coefs.vj <- coefs[coefs.names.vj]
          # new.design.mat <- as.matrix(model.matrix(tt, new.dat.additive)[, bs_term_for_vj])
          #
          # if(!is.factor(v[, j])) {
          #   design.mat <- as.matrix(model.matrix(tt, v)[, bs_term_for_vj])
          #   mean.point <- apply(design.mat, 2, mean)
          #   new.design.mat <- sweep(new.design.mat, 2, mean.point, FUN = "-")
          # }
          #
          # preds.j.additive <- new.design.mat %*% coefs.vj
          # beta.vcov <- sandwich::vcovHC(additive_model$model, type="HC")[coefs.names.vj, coefs.names.vj]
          # sigma2hat <- diag(new.design.mat %*% beta.vcov %*% t(new.design.mat))
          # ci.l <- preds.j.additive-1.96*sqrt(sigma2hat)
          # ci.u <- preds.j.additive+1.96*sqrt(sigma2hat)


          # preds.j.additive <- predict.lm(additive_model$model, newdata=new.dat.additive)
          # m <- model.frame(tt, new.dat.additive)
          # design.mat <- model.matrix(tt, m)
          # design.mat[, apply(design.mat, 2, function(u) length(unique(u))==1)] <- 0
          # preds.j.additive <- design.mat %*% coef(additive_model$model)
          # beta.vcov <- sandwich::vcovHC(additive_model$model)
          # sigma2hat <- diag(design.mat %*% beta.vcov %*% t(design.mat))
          # ci.l <- preds.j.additive-1.96*sqrt(sigma2hat)
          # ci.u <- preds.j.additive+1.96*sqrt(sigma2hat)
          gam.level <- additive_effect_profile(additive_model$model, pseudo.y[[alg]],
                                                v, colnames(v)[j], v0.short[[j]])
          additive_res[[alg]][[j]] <- list(res=gam.level$res,
                                           component.res=gam.level$component.res,
                                           legacy.res=res.gam,
                                           ate=gam.level$ate,
                                           drl.form=additive_model$drl.form,
                                           model=additive_model$model,
                                           risk=additive_model$risk,
                                           res.list=additive_model$fits)
        }
      }
    }

  }

  out <- lapply(learner, function(w) apply(est[[w]], c(1, 2), mean))

  cate.v.res <- list(est=out,
                     fold.est=est,
                     fold.est.pi=est.pi,
                     pseudo=pseudo.y,
                     pseudo.tr=pseudo.y.tr,
                     cate.v.sample=ites_v,
                     v0.long=v0.long,
                     v0.short=v0.short,
                     stage2.reg.data.v=stage2.reg.data.v,
                     stage2.reg.model=reg.model)

  cate.x.res <- list(pseudo=pseudo.y,
                     pseudo.tr=pseudo.y.tr,
                     cate.x.sample=ites_x,
                     cate.x.sample.tr=ites.x.tr,
                     stage2.reg.data.x=stage2.reg.data.x)

  ret <- list(cate.v.res=cate.v.res,
              cate.x.res=cate.x.res,
              univariate.res=univariate_res,
              pd.res=pd_res,
              additive.res=additive_res,
              robinson.res=robinson_res,
              v0.long=v0.long, v0.short=v0.short,
              foldid=s,
              x=x, y=y, a=a, v=v,
              drl.x=drl.x, drl.v=drl.v)
  return(ret)
}
