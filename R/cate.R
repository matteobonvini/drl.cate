#' CATE
#'
#' This function is a wrapper to estimate HTEs defined as E(Y^1 - Y^0 | V = v0)
#' @param v0 evaluation points, i.e. E(Y^1 - Y^0 | V = v0)
#' @param learner string specifying which learners to use, e.g. DR-Learner.
#' @param y vector of outcomes
#' @param a binary vector of treatments
#' @param x matrix of confounders
#' @param v matrix of effect modifiers
#' @param nsplits a number indicating the number of splits used to do
#' cross-fitting. Ignored if foldid is specified.
#' @param foldid id vector specifying which observation belongs to which fold.
#' @return a list containing the following components:
#' \item{est}{ estimates of the CATE at v0}
#' \item{fold_est}{ list of length(v0)xnsplits matrices with estimates
#' of the CATE at v0 in each fold, one matrix for each learner.}
#' @export
#' @references Kennedy, EH. (2020). Optimal Doubly Robust Estimation of
#' Heterogeneous Causal Effects. \emph{arXiv preprint arXiv:2004.14497}.

cate <- function(data_frame, learner, x_names, y_name, a_name, v_names,
                 num_grid = 100, nsplits = 5, foldid = NULL, v0.long = NULL,
                 univariate_reg = FALSE, partial_dependence = FALSE,
                 additive_approx = FALSE, variable_importance = FALSE, bw.stage2 = NULL,
                 sample.split.cond.dens = FALSE, ...) {

  option <- .parse.cate(learner, ...)

  dta <- get_input(data = data_frame, x_names = x_names, y_name = y_name,
                   a_name = a_name, v_names = v_names, v0.long = v0.long,
                   num_grid = num_grid)

  a <- dta$a
  v <- dta$v
  v0.long <- dta$v0.long
  v0 <- dta$v0
  y <- dta$y
  x <- dta$x

  n <- length(y)
  n.eval.pts <- nrow(v0.long)

  if(is.null(foldid)) {
    s <- sample(rep(1:nsplits, ceiling(n / nsplits))[1:n])
  } else {
    s <- foldid
    nsplits <- length(unique(foldid))
  }

  est <- est.pi <- replicate(length(learner),
                             array(NA, dim = c(n.eval.pts, 3, nsplits)),
                             simplify = FALSE)
  univariate_res <- pd_res <- additive_res <-
    replicate(length(learner), vector("list", ncol(v)), simplify = FALSE)

  cate.w.fit <-
    replicate(ncol(v), vector("list", nsplits), simplify = FALSE)

  pseudo.y <- replicate(length(learner), matrix(NA, ncol = 1, nrow = n),
                        simplify = FALSE)
  pseudo.y.pd <- theta.bar <- cond.dens.vals <- cate.w.vals <-
    replicate(length(learner), matrix(NA, ncol = ncol(v), nrow = n),
              simplify = FALSE)
  ites_v <- replicate(length(learner), matrix(NA, ncol = 3, nrow = n),
                      simplify = FALSE)
  ites_x <- replicate(length(learner), matrix(NA, ncol = 1, nrow = n),
                      simplify = FALSE)
  names(est) <- names(est.pi) <- names( pseudo.y) <- names(ites_v) <-
    names(ites_x) <- names(pseudo.y.pd) <- names(theta.bar) <-
    names(cond.dens.vals) <- names(cate.w.vals) <- learner

  stage2.reg.data <- vector("list", nsplits)
  drl.form <- reg.model <- vector("list", nsplits)
  stage2.reg.data.pd <- replicate(ncol(v), vector("list", nsplits),
                                  simplify = FALSE)

  for(k in 1:nsplits) {
    print(paste0("Considering split # ", k, " out of ", nsplits))
    test.idx <- k == s
    train.idx <- k != s
    if(all(!train.idx)) train.idx <- test.idx
    n.te <- sum(test.idx)
    n.tr <- sum(train.idx)

    x.tr <- x[train.idx, , drop = FALSE]
    a.tr <- a[train.idx]
    y.tr <- y[train.idx]
    x.te <- x[test.idx, , drop = FALSE]
    a.te <- a[test.idx]
    y.te <- y[test.idx]

    if ((!is.matrix(v)) & (!is.data.frame(v))){
      v.te <- v[test.idx]
      v.tr <- v[train.idx]
    } else{
      v.te <- v[test.idx, , drop = FALSE]
      v.tr <- v[train.idx, , drop = FALSE]
    }
    # option <- list(pi.x = pi.x, mu0.x = mu0.x, mu1.x = mu1.x,
    #                drl.v = drl.v, drl.x = drl.x, cond.dens = cond.dens,
    #                cate.w = cate.w)

    pihat <- option$pi.x(a = a.tr, x = x.tr, new.x = x.te)$res

    if(any(learner %in% c("dr", "t"))) {

      mu0hat.vals <- option$mu0.x(y = y.tr, a = a.tr, x = x.tr,
                                  new.x = rbind(x.te, x.tr))$res
      mu0hat <- mu0hat.vals[1:n.te]

      mu1hat.vals <- option$mu1.x(y = y.tr, a = a.tr, x = x.tr,
                                  new.x = rbind(x.te, x.tr))$res
      mu1hat <- mu1hat.vals[1:n.te]

    }

    if(any(learner %in% c("u", "r"))) {

      muhat <- option$mu.x(y = y.tr, x = x.tr, new.x = x.te)

    }

    for(alg in learner) {

      if(alg == "dr") {

        pseudo <- (a.te - pihat) / (pihat * (1 - pihat)) *
          (y.te - a.te * mu1hat - (1 - a.te) * mu0hat) + mu1hat - mu0hat

        drl.res <-  option$drl.v(y = pseudo, x = v.te,
                                 new.x = rbind(v0.long, v.te))
        drl.res.pi <-  option$drl.v(y = mu1hat - mu0hat, x = v.te,
                                    new.x = rbind(v0.long, v.te)) # plug-in
        stage2.reg.data[[k]] <- cbind(data.frame(pseudo = pseudo,
                                                 mu1hat = mu1hat,
                                                 mu0hat = mu0hat,
                                                 pihat = pihat,
                                                 y = y.te,
                                                 a = a.te,
                                                 fold.id = k), v.te)

        drl.form[[k]] <- drl.res$drl.form
        reg.model[[k]] <- drl.res$model

        drl.vals <-  as.matrix(drl.res$res)
        drl.vals.pi <-  as.matrix(drl.res.pi$res)
        drl.vals.x <- as.matrix(option$drl.x(y = pseudo, x = x.te,
                                             new.x = x.te)$res)

        est[[alg]][, , k] <- drl.vals[1:n.eval.pts, ]
        est.pi[[alg]][, , k] <- drl.vals.pi[1:n.eval.pts, ]
        pseudo.y[[alg]][test.idx, 1] <- pseudo
        ites_v[[alg]][test.idx, ] <- drl.vals[-c(1:n.eval.pts), ]
        ites_x[[alg]][test.idx, 1] <- drl.vals.x[,1]


        if(partial_dependence) {

          for(j in 1:ncol(v)) {

            v1.j.tr <- v.tr[, j]
            v2.not.v1.j.tr <- v.tr[, -j, drop = FALSE]
            v1.j.te <- v.te[, j]
            v2.not.v1.j.te <- v.te[, -j, drop = FALSE]

            w.tr <- cbind(v1j = v1.j.tr, v2.not.v1.j.tr)
            w.te <-  cbind(v1j = v1.j.te, v2.not.v1.j.te)

            if(sample.split.cond.dens){
              cond.dens.fit <- option$cond.dens(v1 = v1.j.tr,
                                                v2 = v2.not.v1.j.tr)
              cond.dens.vals[[alg]][test.idx, j] <-
                cond.dens.fit$predict.cond.dens(v1 = v1.j.tr,
                                                v2 = v2.not.v1.j.tr,
                                                new.v1 = v1.j.te,
                                                new.v2 = v2.not.v1.j.te)
              if(sum(cond.dens.vals.te < 0.001) > 0) {
                warning(paste0("Effect modifier # ", j, ". There are ", sum(cond.dens.vals.te < 0.01),
                               " conditional density values < 0.01. They will ",
                               "truncated at 0.01."))
                cond.dens.vals.te[cond.dens.vals.te < 0.001] <- 0.01
              }
            }

            # print(range(1/cond.dens.vals.te))
            cate.tr <- mu1hat.vals[-c(1:n.te)] - mu0hat.vals[-c(1:n.te)]

            cate.w.fit[[j]][[k]] <- option$cate.w(tau = cate.tr, w = w.tr,
                                                  new.w = w.tr)

            cate.w.te <- cate.w.fit[[j]][[k]]$fit(new.w = w.te)
            cate.w.vals[[alg]][test.idx, j] <- cate.w.te

            if(n.te > 1000) {

              if(is.factor(v1.j.te)) {
                v1.j.seq <- factor(levels(v1.j.te), levels = levels(v1.j.te))
              } else {
                v1.j.seq <- seq(min(v1.j.te), max(v1.j.te), length.out = 100)
              }

              # tmp.cond.dens.fn <- Vectorize(function(u) {
              #   mean(cond.dens.fit$predict.cond.dens(v1 = v1.j.tr,
              #                                        v2 = v2.not.v1.j.tr,
              #                                        new.v1 = u,
              #                                        new.v2 = v2.not.v1.j.te))
              # }, vectorize.args = "u")

              # tmp.cond.dens.fn <- Vectorize(function(u) {
              #   mean(cond.dens.fit$predict.cond.dens(v1 = v[, j],
              #                                        v2 = v[, -j, drop = FALSE],
              #                                        new.v1 = u,
              #                                        new.v2 = v[, -j, drop = FALSE]))
              # }, vectorize.args = "u")

              tmp.cate.w.fit.fn <- Vectorize(function(u) {
                mean(cate.w.fit[[j]][[k]]$fit(new.w = cbind(v1j = u,
                                                            v2.not.v1.j.te)))
              }, vectorize.args = "u")

              # marg.dens.vals <- tmp.cond.dens.fn(v1.j.seq)
              cate.w.avg.vals <- tmp.cate.w.fit.fn(v1.j.seq)
              if(is.factor(v1.j.te)) {
                # names(marg.dens.vals) <- v1.j.seq
                # marg.dens <- theta.bar.vals <- rep(NA, length(v1.j.te))
                theta.bar.vals <- rep(NA, length(v1.j.te))
                for(u in levels(v1.j.te)) {
                  # marg.dens[v1.j.te == u] <-
                  #   marg.dens.vals[names(marg.dens.vals) == u]
                  theta.bar.vals[v1.j.te == u] <-
                    cate.w.avg.vals[v1.j.seq == u]
                }
              } else {
                # marg.dens <- approx(x = as.numeric(v1.j.seq), y = marg.dens.vals,
                #                     xout = v1.j.te)$y
                # marg.dens <- approx(x = as.numeric(v1.j.seq), y = marg.dens.vals,
                #                     xout = v[, j], rule = 2)$y
                # marg.dens2 <- ks::kde(x = v[, j],
                #                       eval.points = v[, j],
                #                       density = TRUE)
                theta.bar.vals <- approx(x = v1.j.seq, y = cate.w.avg.vals,
                                         xout = v1.j.te, rule = 2)$y
              }
            }
            else {
              w.long.test <- cbind(v1j = rep(v1.j.te, each = n.te),
                                   v2.not.v1.j.te[rep(1:n.te, n.te), ,
                                                  drop = FALSE])
              if(sample.split.cond.dens) {
                cond.dens.preds <-
                  cond.dens.fit$predict.cond.dens(v1 = v1.j.tr,
                                                  v2 = v2.not.v1.j.tr,
                                                  new.v1 = w.long.test[, 1],
                                                  new.v2 = w.long.test[, -1, drop = FALSE])

                marg.dens <- colMeans(matrix(cond.dens.preds, ncol = n.te,
                                             nrow = n.te))
              }


              cate.preds <- cate.w.fit[[j]][[k]]$fit(new.w = w.long.test)

              theta.bar.vals <-  colMeans(matrix(cate.preds,
                                                 nrow = n.te, ncol = n.te))

            }

            theta.bar[[alg]][test.idx, j] <- theta.bar.vals
            # ghat <- marg.dens / cond.dens.vals.te
            # pseudo.y.pd.vals <- (pseudo - cate.w.te) * ghat + theta.bar.vals

            # fit <- locpol::locpol(y ~ x, data = data.frame(y = pseudo.y.pd.vals,
            #                             x = v[, j]))
            # pseudo.y.pd[[alg]][test.idx, j] <- pseudo.y.pd.vals
            # pseudo.y.pd[[alg]][, j] <- pseudo.y.pd.vals

            # data.pd <- data.frame(pseudo.pd = pseudo.y.pd.vals,
            data.pd <- data.frame(pseudo.cate = pseudo,
                                  mu1hat = mu1hat,
                                  mu0hat = mu0hat,
                                  pihat = pihat,
                                  tauhat.w = cate.w.te,
                                  # marg.dens = marg.dens,
                                  # cond.dens = cond.dens.vals.te,
                                  theta.bar = theta.bar.vals,
                                  y = y.te,
                                  a = a.te,
                                  fold.id = k)

            stage2.reg.data.pd[[j]][[k]] <- cbind(data.pd, v.te)

          }
        }
      } else if(alg == "t"){
        if(all(colnames(x) %in% colnames(v))) {
          next # no need to smooth the t-learner if V = X
        } else {
          pseudo <- mu1hat - mu0hat
          est[[alg]][, , k] <- tl(y = pseudo, x = v.te, new.x = v0.long)
        }
      } else if(alg == "u") {

        pseudo <- (y.te - muhat) / (a.te - pihat)
        est[[alg]][, , k] <- ul(y = pseudo, x = v.te, new.x = v0.long)

      } else if(alg == "lp-r") next
    }
  }
  print("Done with fitting nuisance functions.")
  for(alg in learner) {

    if(alg != "dr") next

    if(additive_approx) {
      additive_model <- drl.basis.additive(y = pseudo.y[[alg]][, 1], x = v,
                                           new.x = v)
      tt <- delete.response(terms(additive_model$model))
    }
    if(univariate_reg | partial_dependence | additive_approx) {

      for(j in 1:ncol(v)){
        vj <- v[, j]
        is.var.factor <-
          paste0(class(vj), collapse = " ") %in% c("factor", "ordered factor")

        if(partial_dependence) {

          if(!sample.split.cond.dens) {
            cond.dens.fit <- option$cond.dens(v1 = v[, j],
                                              v2 = v[, -j, drop = FALSE])
            cond.dens.vals[[alg]][, j] <-
              cond.dens.fit$predict.cond.dens(v1 = v[, j],
                                              v2 = v[, -j, drop = FALSE],
                                              new.v1 = v[, j],
                                              new.v2 = v[, -j, drop = FALSE])

            if(sum(cond.dens.vals[[alg]][, j] < 0.001) > 0) {
              warning(paste0("Effect modifier # ", j, ". There are ",
                             sum(cond.dens.vals[[alg]][, j] < 0.001),
                             " conditional density values < 0.001. They will ",
                             "truncated at 0.001."))
              cond.dens.vals[[alg]][cond.dens.vals[[alg]][, j] < 0.001, j] <- 0.001
            }

          }
          if(is.var.factor) {
            marg.dens <- rep(NA, length(vj))
            for(u in levels(vj)) marg.dens[vj == u] <- mean(vj == u)
          } else {
            marg.dens <- ks::kde(x = v[, j], eval.points = v[, j],
                                 density = TRUE)$estimate
          }
          ghat <- marg.dens / cond.dens.vals[[alg]][, j]
          pseudo.y.pd[[alg]][, j] <-
            (pseudo.y[[alg]][, 1] - cate.w.vals[[alg]][, j]) * ghat +
            theta.bar[[alg]][, j]
        }
        if(length(unique(vj)) < 15 | is.var.factor) {
          if(univariate_reg) {
            univariate_res[[alg]][[j]] <-
              list(data = data.frame(pseudo = pseudo.y[[alg]][, 1],
                                     exposure = vj),
                   res = lm.discrete.v(y = pseudo.y[[alg]][, 1], x = vj,
                                       new.x = levels(vj))
              )
          }
          if(partial_dependence) {
            pd_res[[alg]][[j]] <-
              list(data = data.frame(pseudo = pseudo.y.pd[[alg]][, j],
                                     exposure = vj),
                   res = lm.discrete.v(y = pseudo.y.pd[[alg]][, j], x = vj,
                                       new.x = levels(vj))
              )

          }
        }
        else {
          if(univariate_reg) {
            fit_debias_inf <- debiased_inference(A = vj,
                                                 pseudo.out = pseudo.y[[alg]][, 1],
                                                 eval.pts = v0[[j]],
                                                 bandwidth.method = "LOOCV(h=b)",
                                                 kernel.type = "gau",
                                                 bw.seq = bw.stage2[[j]])
            univariate_res[[alg]][[j]] <-
              list(
                data = data.frame(pseudo = pseudo.y[[alg]][, 1],
                                  exposure = vj),
                res = fit_debias_inf$res,
                risk = fit_debias_inf$risk
              )
          }

          if(partial_dependence) {

            muhat.vals <- .get.muhat(splits.id = s, cate.w.fit = cate.w.fit[[j]],
                                     v1 = vj, v2 = v[, -j, drop = FALSE],
                                     max.n.integral = 1000)

            fit_debias_inf <- debiased_inference(A = vj,
                                                 pseudo.out = pseudo.y.pd[[alg]][, j],
                                                 eval.pts = v0[[j]],
                                                 mhat.obs = theta.bar[[alg]][, j],
                                                 muhat.vals = muhat.vals,
                                                 bandwidth.method = "LOOCV(h=b)",
                                                 kernel.type = "gau",
                                                 bw.seq = bw.stage2[[j]])

            pd_res[[alg]][[j]] <-
              list(
                data = data.frame(pseudo = pseudo.y.pd[[alg]][, j],
                                  cond.dens.vals = cond.dens.vals[[alg]][, j],
                                  exposure = vj),
                res = fit_debias_inf$res,
                risk = fit_debias_inf$risk
              )
          }
        }

        if(additive_approx){

          new.dat.additive <- as.data.frame(matrix(0, nrow = length(v0[[j]]),
                                                   ncol = ncol(v),
                                                   dimnames = list(NULL,
                                                                   colnames(v))))
          for(l in 1:ncol(v)) {
            if(l == j) {
              new.dat.additive[, l] <- v0[[j]]
            } else {
              if(is.factor(v[, l])) {
                new.dat.additive[, l] <- factor(levels(v[, l])[1],
                                                levels = levels(v[, l]))
              } else {
                new.dat.additive[, l] <- min(v[, l])
              }
            }
          }

          preds.j.additive <- predict.lm(additive_model$model,
                                         newdata = new.dat.additive)
          m <- model.frame(tt, new.dat.additive)
          design.mat <- model.matrix(tt, m)
          beta.vcov <- sandwich::vcovHC(additive_model$model)
          sigma2hat <- diag(design.mat %*% beta.vcov %*% t(design.mat))

          ci.l <- preds.j.additive - 1.96 * sqrt(sigma2hat)
          ci.u <- preds.j.additive + 1.96 * sqrt(sigma2hat)
          additive_res[[alg]][[j]] <- data.frame(eval.pts = v0[[j]],
                                                 theta = preds.j.additive,
                                                 ci.ll.pts = ci.l,
                                                 ci.ul.pts = ci.u,
                                                 ci.ll.unif = NA,
                                                 ci.ul.unif = NA)
        }
      }
    }

    if(variable_importance){
      print('start calculating vimp')
      tau_hat <- ites_x[[alg]][,1]
      pseudo_hat <- pseudo.y[[alg]]
      vimp_df <- get_VIMP(tau_hat, pseudo_hat, x, v_names)
      draw_VIMP(vimp_df)
    } else{
      vimp_df <- data.frame(matrix(ncol = 4, nrow = ncol(v_names)))
    }
  }


  if(any(learner == "lp-r")) {

    est[["lp-r"]] <- lp_r_learner(x0 = v0.long, y = y, a = a, x = x,
                                  mu.x = mu.x,
                                  pi.x = pi.x, basis = params[["basis"]],
                                  order_basis = params[["order_basis"]],
                                  kernel = params[["kernel"]])$fold_est

  }

  if(any(learner == "t") & all(colnames(x) %in% colnames(v))) {
    est[["t"]][, 1] <- mu1.x(y = y, a = a, x = x, new.x = v0.long) -
      mu0.x(y.tr = y, a.tr = a, x = x, new.x = v0.long)
  }

  out <- lapply(learner, function(w) apply(est[[w]], c(1, 2), mean))

  ret <- list(est = out, fold_est = est, fold_est_pi = est.pi,
              pseudo.y = pseudo.y, ites_v = ites_v,
              ites_x = ites_x, v0.long = v0.long, v0 = v0, v = v,
              drl.form = drl.form, reg_model = reg.model,
              stage2.reg.data = stage2.reg.data,
              stage2.reg.data.pd = stage2.reg.data.pd,
              univariate_res = univariate_res,
              pd_res = pd_res,
              additive_res = additive_res,
              foldid = s, vimp_df = vimp_df)
  return(ret)
}

