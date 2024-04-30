#' CATE
#'
#' This function estimates heterogeneous treatment effects (HTEs) defined as E(Y^1 - Y^0 | V = v0).
#' @param data A data frame containing the dataset.
#' @param learner A character string specifying which learner to use (e.g., "dr").
#' @param x_names A character vector specifying the names of the confouding variables.
#' @param y_name A character string specifying the outcome variable.
#' @param a_name A character string specifying the treatment variable.
#' @param v_names A character vector specifying the names of the effect modifiers.
#' @param v0 A matrix of evaluation points, i.e., values of V for which the CATE is estimated (E(Y^1 - Y^0 | V = v0)).
#' @param mu1.x A function taking arguments (y, a, x, new.x). It trains a model estimating
#' E(Y | A = 1, X) and returns a list of 3 elements: res, model and fit. \emph{res} is a vector of predictions of the model evaluated at
#' new.x, \emph{model} is the model object used to estimate E(Y | A = 1, X) and \emph{fit} is a function with argument new.x that
#' returns the predictions of the model. See examples.
#' @param mu0.x A function taking arguments (y, a, x, new.x). It trains a model estimating
#' E(Y | A = 0, X) and returns a list of 3 elements: res, model and fit. \emph{res} is a vector of predictions of the model evaluated at
#' new.x, \emph{model} is the model object used to estimate E(Y | A = 0, X) and \emph{fit} is a function with argument new.x that
#' returns the predictions of the model. See examples.
#' @param pi.x A function taking arguments (a, x, new.x). It trains a model estimating
#' P(A = 1 | X) and returns a list of 3 elements: res, model and fit. \emph{res} is a vector of predictions of the model evaluated at
#' new.x, \emph{model} is the model object used to estimate P(A = 1 | X) and \emph{fit} is a function with argument new.x that
#' returns the predictions of the model. See examples.
#' @param drl.v A function taking arguments (pseudo, v, new.v). It trains a model estimating E(Y^1 - Y^0 | V) by
#' regressing a pseudo-outcome \emph{pseudo} onto v and returns a list of 3 elements: res, model and fit.
#' \emph{res} is a vector of predictions of the model evaluated at new.v,
#' \emph{model} is the model object used to estimate E(Y^1 - Y^0 | V) (after possibly model selection)
#' and \emph{fit} is a function with argument new.v that returns the predictions of the model. See examples.
#' #' @param drl.x A function taking arguments (pseudo, x, new.x). It trains a model estimating E(Y^1 - Y^0 | X) by
#' regressing a pseudo-outcome \emph{pseudo} onto x and returns a list of 3 elements: res, model and fit.
#' \emph{res} is a vector of predictions of the model evaluated at new.x,
#' \emph{model} is the model object used to estimate E(Y^1 - Y^0 | X) (after possibly model selection)
#' and \emph{fit} is a function with argument new.v that returns the predictions of the model. See examples.
#' @param nsplits An integer indicating the number of splits used for cross-validation. Ignored if foldid is specified.
#' @param foldid An optional vector specifying fold assignments for cross-validation.
#' @param univariate_reg A logical indicating whether to perform univariate regression for estimating the CATE
#' as a function of each effect modifier separately (default: FALSE).
#' @param partial_dependence A logical indicating whether to compute partial dependence plots (default: FALSE).
#' @param partially_linear A logical indicating whether to compute partially linear approximations
#' via Robinson's transformation (default: FALSE).
#' @param additive_approx A logical indicating whether to compute the CATE assuming an additive structure (default: FALSE).
#' @param variable_importance A logical indicating whether to compute variable importance measures (default: FALSE).
#' @param vimp_num_splits An integer specifying the number of splits for variable importance computation (default: 1).
#' @param bw.stage2 A list of length equal to the number of effect modifiers considered, where each element if a vector of
#' candidate bandwidths for second-stage regression of the pseudo-outcome onto the effect modifier
#' that calculates either the univariate CATE or the Partial Dependence measure (default: NULL).
#' It needs to be provided if \emph{univariate_reg} or \emph{partial_dependence} is set to TRUE.
#' @param sample.split.cond.dens A logical indicating whether to do sample-splitting for conditional density estimation
#' (default: FALSE).
#' @param cond.dens A function
#' @param cate.w A function
#' @param cate.not.j A function
#' @param reg.basis.not.j A function
#' @param pl.dfs A list of length equal to the number of effect modifiers considered, where each element is a vector of
#' candidate number of basis elements for the partially linear approximation computed via Robinson trick.
#' @return A list containing the estimated CATE at v0 and per-fold estimates of the CATE at v0 for each learner.
#' @export
#' @references Kennedy, EH. (2020). Optimal Doubly Robust Estimation of
#' Heterogeneous Causal Effects. \emph{arXiv preprint arXiv:2004.14497}.

cate <- function(data, learner, x_names, y_name, a_name, v_names, v0,
                 mu1.x, mu0.x, pi.x, drl.v, drl.x,
                 nsplits=5,
                 foldid=NULL,
                 univariate_reg=FALSE,
                 partial_dependence=FALSE,
                 partially_linear=FALSE,
                 additive_approx=FALSE,
                 variable_importance=FALSE,
                 vimp_num_splits=1,
                 bw.stage2=NULL,
                 sample.split.cond.dens=FALSE,
                 cond.dens=NULL,
                 cate.w=NULL,
                 cate.not.j=NULL,
                 reg.basis.not.j=NULL,
                 pl.dfs=NULL) {

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

  if(is.null(foldid)) {
    s <- sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  } else {
    s <- foldid
    nsplits <- length(unique(foldid))
  }

  est <- est.pi <- replicate(length(learner), array(NA, dim=c(n.eval.pts, 3, nsplits)),
                             simplify=FALSE)
  univariate_res <- pd_res <- additive_res <- robinson_res <-
    replicate(length(learner), vector("list", ncol(v)), simplify=FALSE)

  cate.w.fit <- replicate(ncol(v), vector("list", nsplits), simplify=FALSE)

  pseudo.y <- replicate(length(learner), rep(NA, n), simplify=FALSE)
  pseudo.y.tr <- ites.x.tr <-
    replicate(length(learner), vector("list", nsplits), simplify=FALSE)
  pseudo.y.pd <- theta.bar <- cond.dens.vals <- cate.w.vals <-
    replicate(length(learner), matrix(NA, ncol=ncol(v), nrow=n), simplify=FALSE)
  ites_v <- ites_x <- replicate(length(learner), matrix(NA, ncol=3, nrow=n), simplify=FALSE)
  names(est) <- names(est.pi) <- names(pseudo.y) <- names(ites_v) <-
    names(ites_x) <- names(pseudo.y.pd) <- names(theta.bar) <-
    names(cond.dens.vals) <- names(cate.w.vals) <- names(pseudo.y.tr) <- learner

  stage2.reg.data.v <- stage2.reg.data.x <- reg.model <- vector("list", nsplits)
  stage2.reg.data.pd <- replicate(ncol(v), vector("list", nsplits), simplify=FALSE)

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

              if(sample.split.cond.dens){
                cond.dens.fit <- cond.dens[[j]](v1=v1.j.tr, v2=v2.not.v1.j.tr)
                cond.dens.vals[[alg]][test.idx, j] <-
                  cond.dens.fit$predict.cond.dens(v1=v1.j.tr, v2=not.v1.j.tr,
                                                  new.v1=v1.j.te, new.v2=not.v1.j.te)
                if(sum(cond.dens.vals.te < 0.001) > 0) {
                  warning(paste0("Effect modifier # ", j, ". There are ",
                                 sum(cond.dens.vals.te < 0.001),
                                 " conditional density values < 0.01. They will ",
                                 "truncated at 0.01."))
                  cond.dens.vals.te[cond.dens.vals.te < 0.001] <- 0.01
                }
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
                  theta.bar.vals <- approx(x=v1.j.seq, y=cate.w.avg.vals,
                                           xout=v1.j.te, rule=2)$y
                }
              }
              else {
                w.long.test <- cbind(v1j=rep(v1.j.te, each=n.te),
                                     not.v1.j.te[rep(1:n.te, n.te), , drop=FALSE])
                if(sample.split.cond.dens) {
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
      message(conditionMessage(cond))
      # each_filename <- paste0('sim_data_', as.character(n), '-',
      #                         as.character(iter), '.rda')
      # save(data, file = "/Users/matteobonvini/Dropbox/Hetero Effects/Simulations/debug/dataset_error.rda")
      # save(s, file = "/Users/matteobonvini/Dropbox/Hetero Effects/Simulations/debug/foldid_error.rda")
      # NULL
    }
  )
  if(is.null(tmp)){
    warning("Encountered error while fitting nuisance functions")
    # print("Encountered error while fitting nuisance functions, exiting procedure.")
    univ.res <- pd.res <- add.res <- rob.res <- data.frame(theta = rep(NA, 10),
                                                           theta.debias = rep(NA, 10))
    return(list(univariate_res = list(dr = list(list(res = univ.res))),
                pd_res = list(dr = list(list(res = pd.res))),
                additive_res = list(dr = list(list(res = add.res))),
                robinson_res = list(dr = list(list(res = rob.res)))))
  }
  for(alg in learner) {

    if(alg != "dr") stop("Only learner = df is currently implemented.")

    if(additive_approx) {
      additive_model <- drl.basis.additive(y=pseudo.y[[alg]], x=v, new.x=v)
      tt <- delete.response(terms(additive_model$model))
    }
    if(univariate_reg | partial_dependence | additive_approx | partially_linear) {

      for(j in 1:ncol(v)){
        vj <- v[, j]
        is.var.factor <- paste0(class(vj), collapse = " ") %in% c("factor", "ordered factor")

        if(partially_linear) {

          j.robinson <- robinson(pseudo=pseudo.y[[alg]],
                                 w=v[, -j, drop = FALSE],
                                 v=v[, j],
                                 new.v=v0.short[[j]],
                                 s=s,
                                 cate.not.j=cate.not.j[[j]],
                                 reg.basis.not.j=reg.basis.not.j[[j]],
                                 dfs=pl.dfs[[j]])

          rob.data <- data.frame(eval.pts=v0.short[[j]],
                                 theta=j.robinson$res$preds,
                                 ci.ll.pts=j.robinson$res$ci.ll,
                                 ci.ul.pts=j.robinson$res$ci.uu,
                                 ci.ll.unif=NA,
                                 ci.ul.unif=NA)

          robinson_res[[alg]][[j]] <- list(res=rob.data,
                                           model=j.robinson$model,
                                           risk=j.robinson$risk,
                                           fits=j.robinson$fits)

        }

        if(partial_dependence) {

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
            marg.dens <- ks::kde(x=vj, eval.points=vj, density=TRUE)$estimate
          }
          ghat <- marg.dens/cond.dens.vals[[alg]][, j]
          pseudo.y.pd[[alg]][, j] <-
            (pseudo.y[[alg]]-cate.w.vals[[alg]][, j])*ghat + theta.bar[[alg]][, j]
        }
        if(length(unique(vj)) < 15 | is.var.factor) {
          if(univariate_reg) {
            res.empVar <- NULL
            for(ll in 1:length(unique(vj))) {
              pts.vj <- unique(vj)[ll]
              if.vals <- pseudo.y[[alg]]*I(vj==pts.vj) / mean(I(vj==pts.vj))
              tmp <- data.frame(eval.pts=pts.vj,
                                    theta=mean(if.vals),
                                    ci.ll.pts=mean(if.vals) - 1.96*sqrt(var(if.vals)/n),
                                    ci.ul.pts=mean(if.vals) + 1.96*sqrt(var(if.vals)/n),
                                    ci.ul.unif=NA,
                                    ci.ll.unif=NA)
              res.empVar <- rbind(res.empVar, tmp)
            }
            univariate_res[[alg]][[j]] <-
              list(data=data.frame(pseudo=pseudo.y[[alg]], exposure=vj),
                   res=lm.discrete.v(y=pseudo.y[[alg]], x=vj, new.x=unique(vj)),
                   res.empVar=res.empVar)
          }
          if(partial_dependence) {
            res.empVar <- NULL
            for(ll in 1:length(unique(vj))) {
              pts.vj <- unique(vj)[ll]
              if.vals <- (pseudo.y[[alg]]-cate.w.vals[[alg]][, j]) *
                I(vj==pts.vj)/cond.dens.vals[[alg]][, j]  +
                mean(theta.bar[[alg]][I(vj==pts.vj), j])
              tmp <- data.frame(eval.pts=pts.vj,
                                     theta=mean(if.vals),
                                     ci.ll.pts=mean(if.vals) - 1.96*sqrt(var(if.vals)/n),
                                     ci.ul.pts=mean(if.vals) + 1.96*sqrt(var(if.vals)/n),
                                     ci.ul.unif=NA,
                                     ci.ll.unif=NA)
              res.empVar <- rbind(res.empVar, tmp)
            }
            pd_res[[alg]][[j]] <-
              list(data=data.frame(pseudo=pseudo.y.pd[[alg]][, j], exposure=vj),
                   stage2.reg.data.pd=stage2.reg.data.pd,
                   res=lm.discrete.v(y=pseudo.y.pd[[alg]][, j], x=vj, new.x=unique(vj)),
                   res.empVar=res.empVar)
          }
        }
        else {
          if(univariate_reg) {
            fit_debias_inf <- debiased_inference(A=vj,
                                                 pseudo.out=pseudo.y[[alg]],
                                                 eval.pts=v0.short[[j]],
                                                 bandwidth.method="LOOCV(h=b)",
                                                 kernel.type="gau",
                                                 bw.seq=bw.stage2[[j]])
            univariate_res[[alg]][[j]] <-
              list(data=data.frame(pseudo=pseudo.y[[alg]], exposure=vj),
                   res=fit_debias_inf$res,
                   risk=fit_debias_inf$risk)
          }

          if(partial_dependence) {

            muhat.vals <- .get.muhat(splits.id=s, cate.w.fit=cate.w.fit[[j]],
                                     v1=vj, v2=v[, -j, drop=FALSE],
                                     max.n.integral=1000)

            fit_debias_inf <- debiased_inference(A=vj,
                                                 pseudo.out=pseudo.y.pd[[alg]][, j],
                                                 eval.pts=v0.short[[j]],
                                                 mhat.obs=theta.bar[[alg]][, j],
                                                 muhat.vals=muhat.vals,
                                                 bandwidth.method="LOOCV(h=b)",
                                                 kernel.type="gau",
                                                 bw.seq=bw.stage2[[j]])

            pd_res[[alg]][[j]] <-
              list(data=data.frame(pseudo=pseudo.y.pd[[alg]][, j],
                                   cond.dens.vals=cond.dens.vals[[alg]][, j],
                                   exposure=vj),
                   res=fit_debias_inf$res,
                   risk=fit_debias_inf$risk)
          }
        }

        if(additive_approx){

          new.dat.additive <- as.data.frame(matrix(0, nrow=length(v0.short[[j]]),
                                                   ncol=ncol(v),
                                                   dimnames=list(NULL, colnames(v))))
          for(l in 1:ncol(v)) {
            if(l == j) {
              new.dat.additive[, l] <- v0.short[[j]]
            } else {
              if(is.factor(v[, l])) {
                new.dat.additive[, l] <- factor(levels(v[, l])[1], levels=levels(v[, l]))
              } else {
                new.dat.additive[, l] <- min(v[, l])
              }
            }
          }

          preds.j.additive <- predict.lm(additive_model$model, newdata=new.dat.additive)
          m <- model.frame(tt, new.dat.additive)
          design.mat <- model.matrix(tt, m)
          design.mat[, apply(design.mat, 2, function(u) length(unique(u))==1)] <- 0
          preds.j.additive <-  design.mat %*% coef(additive_model$model)
          beta.vcov <- sandwich::vcovHC(additive_model$model)
          sigma2hat <- diag(design.mat %*% beta.vcov %*% t(design.mat))

          ci.l <- preds.j.additive - 1.96*sqrt(sigma2hat)
          ci.u <- preds.j.additive + 1.96*sqrt(sigma2hat)
          additive_res[[alg]][[j]] <- list(res=data.frame(eval.pts=v0.short[[j]],
                                                          theta=preds.j.additive,
                                                          ci.ll.pts=ci.l,
                                                          ci.ul.pts=ci.u,
                                                          ci.ll.unif=NA,
                                                          ci.ul.unif=NA),
                                           drl.form=additive_model$drl.form,
                                           model=additive_model$model,
                                           risk=additive_model$risk)
        }
      }
    }

    if(variable_importance){
      print('start calculating vimp')
      tau_hat <- ites_x[[alg]][,1]
      pseudo_hat <- pseudo.y[[alg]]
      vimp_df <- get_VIMP(tau_hat, pseudo_hat, x, y, a, v_names,
                          vimp_num_splits=vimp_num_splits, option=option)
      draw_VIMP(vimp_df)
    } else{
      vimp_df <- data.frame(matrix(ncol = 4, nrow = length(v_names)))
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
              vimp.df=vimp_df,
              v0.long=v0.long, v0.short=v0.short,
              foldid=s,
              x=x, y=y, a=a, v=v,
              drl.x=drl.x, drl.v=drl.v)
  return(ret)
}

