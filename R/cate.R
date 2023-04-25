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

cate <- function(data_frame, learner, x_names, y_name, a_name, v_names, num_grid = 100,
                 nsplits = 5, foldid = NULL, partial_dependence = TRUE,
                 additive_approx = TRUE, ...) {

  params <- list(...)

  dta <- get_input(data_frame, x_names, y_name, a_name, v_names, num_grid)

  a <- dta$a
  v <- dta$v
  v0.long <- dta$v0.long
  v0 <- dta$v0
  y <- dta$y
  x <- dta$x

  # input data can override assigned covariates?
  v0_input <- params[["v0"]]
  if(!is.null(v0_input)) {v0 <- v0_input}
  v0_long_input <- params[["v0.long"]]
  if(!is.null(v0_long_input)) {v0.long <- v0_long_input}

  n <- length(y)

  if(is.null(foldid)) {
    s <- sample(rep(1:nsplits, ceiling(n / nsplits))[1:n])
  } else {
    s <- foldid
    nsplits <- length(unique(foldid))
  }

  est <- est.pi <- replicate(length(learner),
                             array(NA, dim = c(nrow(v0.long), 3, nsplits)),
                             simplify = FALSE)
  univariate_res <- univariate_debias_res <- pd_res <- pd_debias_res <-
    additive_res <-
    replicate(length(learner), vector("list", ncol(v)), simplify = FALSE)

  pseudo.y <- replicate(length(learner), matrix(NA, ncol = 1, nrow = n),
                        simplify = FALSE)
  pseudo.y.pd <- theta.bar <- replicate(length(learner), matrix(NA, ncol = ncol(v), nrow = n),
                                        simplify = FALSE)
  theta.mat <- replicate(length(learner), array(NA, dim = c(n, n, ncol(v))),
                         simplify = FALSE)
  ites_v <- replicate(length(learner), matrix(NA, ncol = 3, nrow = n),
                      simplify = FALSE)
  ites_x <- replicate(length(learner), matrix(NA, ncol = 1, nrow = n),
                      simplify = FALSE)
  names(est) <- names(est.pi) <- names( pseudo.y) <- names(ites_v) <-
    names(ites_x) <- names(theta.mat) <- names(pseudo.y.pd) <-
    names(theta.bar) <- learner

  stage2.reg.data <- vector("list", nsplits)
  stage2.reg.data.pd <- replicate(ncol(v), vector("list", nsplits), simplify = FALSE)

  if(any(learner == "lp-r")) {
    est[["lp-r"]] <- matrix(NA, ncol = 3, nrow = nrow(v0.long))
  }
  if(any(learner == "t") & all(colnames(x) %in% colnames(v))) {
    est[["t"]] <- matrix(NA, ncol = 1, nrow = nrow(v0.long))
  }

  pi.x <- params[["pi.x"]]
  pi.x.method <- params[["pi.x.method"]]

  if(is.null(pi.x) & !is.null(pi.x.method)) {

    if(pi.x.method == "lasso") pi.x <- pi.x.lasso
    else if(pi.x.method == "glm") pi.x <- pi.x.glm
    else if(pi.x.method == "SL") pi.x <- pi.x.SL(params[["sl.lib.pi"]])
    else stop("Provide valid method for estimating the propensity score.")

  }

  sl.lib <- params[["sl.lib"]]

  if(any(learner %in% c("dr", "t"))) {

    mu1.x.method <- params[["mu1.x.method"]]
    mu0.x.method <- params[["mu1.x.method"]]
    drl.v.method <- params[["drl.v.method"]]
    drl.x.method <- params[["drl.x.method"]]
    mu1.x <- params[["mu1.x"]]
    mu0.x <- params[["mu0.x"]]
    drl.v <- params[["drl.v"]]
    drl.x <- params[["drl.x"]]
    cate.w <- params[["cate.w"]]
    cond.dens <- params[["cond.dens"]]

    if(is.null(mu1.x)) {

      if(mu1.x.method == "lasso") {
        mu1.x <- mu1.x.lasso
      } else if(mu1.x.method == "lm"){
        mu1.x <- mu1.x.lm
      } else if(mu1.x.method == "SL") {
        mu1.x <- mu1.x.SL(sl.lib)
      } else stop("Provide valid method for estimating the E(Y|A=1,X).")

    }

    if(is.null(mu0.x)) {

      if(mu0.x.method == "lasso") {
        mu0.x <- mu0.x.lasso
      } else if(mu0.x.method == "lm"){
        mu0.x <- mu0.x.lm
      } else if(mu0.x.method == "SL") {
        mu0.x <- mu0.x.SL
      } else stop("Provide valid method for estimating the E(Y|A=0,X).")

    }

    if(is.null(drl.v) & any(learner == "dr")) {

      if(drl.v.method == "lasso") {
        drl.v <- drl.lasso
      } else if(drl.v.method == "lm"){
        drl.v <- drl.lm
      } else if(drl.v.method == "glm"){
        drl.v <- drl.glm
      } else if(drl.v.method == "gam") {
        drl.v <- drl.gam
      } else if(drl.v.method == "rf") {
        drl.v <- drl.rf
      } else stop("Provide valid method for second-stage regression.")

    }

    if(is.null(drl.x) & any(learner == "dr")) {
      if(is.null(drl.x.method)) {
        # by default use lasso
        drl.x <- drl.lasso
      } else if (drl.x.method == "lm"){
        drl.x <- drl.ite.lm
      } else if (drl.x.method == "lasso"){
        drl.x <- drl.lasso
      } else stop("Provide valid method for second-stage regression.")
    }
  }

  if(any(learner %in% c("u", "r", "lp-r"))) {

    mu.x.method <- params[["mu.x.method"]]
    ul.method <- params[["ul.method"]]
    ul <- params[["ul"]]

    if(is.null(mu.x)) {

      if(mu.x.method == "lasso") {
        mu.x <- mu.x.lasso
      } else if(mu.x.method == "SL") {
        mu.x <- mu.x.SL
      } else stop("Provide valid method for estimating the E(Y|X).")

    }

    if(is.null(ul) & any(learner == "u")) {

      if(ul.method == "lasso") {
        ul <- ul.lasso
      } else if(ul.method == "SL") {
        ul <- ul.SL
      } else stop("Provide valid method for second-stage regression.")

    }


  }

  for(k in 1:nsplits) {

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

    pihat <- pi.x(a = a.tr, x = x.tr, new.x = x.te)

    if(any(learner %in% c("dr", "t"))) {

      mu0hat.vals <- mu0.x(y = y.tr, a = a.tr, x = x.tr,
                           new.x = rbind(x.te, x.tr))
      mu0hat <- mu0hat.vals[1:n.te]

      mu1hat.vals <- mu1.x(y = y.tr, a = a.tr, x = x.tr,
                           new.x = rbind(x.te, x.tr))
      mu1hat <- mu1hat.vals[1:n.te]

      # mu1hat <- mu1.x(y = y.tr, a = a.tr, x = x.tr, new.x = x.te)
      # mu0hat <- mu0.x(y = y.tr, a = a.tr, x = x.tr, new.x = x.te)

    }

    if(any(learner %in% c("u", "r"))) {

      muhat <- mu.x(y = y.tr, x = x.tr, new.x = x.te)

    }

    for(alg in learner) {

      if(alg == "dr") {

        pseudo <- (a.te - pihat) / (pihat * (1 - pihat)) *
          (y.te - a.te * mu1hat - (1 - a.te) * mu0hat) + mu1hat - mu0hat

        drl.res <-  drl.v(y = pseudo, x = v.te, new.x = rbind(v0.long, v.te))
        drl.res.pi <-  drl.v(y = mu1hat - mu0hat, x = v.te,
                             new.x = rbind(v0.long, v.te))
        stage2.reg.data[[k]] <- cbind(data.frame(pseudo = pseudo,
                                                 mu1hat = mu1hat,
                                                 mu0hat = mu0hat,
                                                 pihat = pihat,
                                                 y = y.te,
                                                 a = a.te,
                                                 fold.id = k), v.te)
        if(k == 1) {
          drl.form <- drl.res$drl.form
          reg.model <- drl.res$model
        }
        drl.vals <-  drl.res$res
        drl.vals.pi <-  drl.res.pi$res

        est[[alg]][, , k] <- drl.vals[1:nrow(v0.long), ]
        est.pi[[alg]][, , k] <- drl.vals.pi[1:nrow(v0.long), ]
        pseudo.y[[alg]][test.idx, 1] <- pseudo
        ites_v[[alg]][test.idx, ] <- drl.vals[(nrow(v0.long)+1):nrow(drl.vals), ]
        ites_x[[alg]][test.idx, 1] <- drl.x(y = pseudo, x = x.te, new.x = x.te)$res[,1]


        if(partial_dependence) {

          for(j in 1:ncol(v)) {

            v1.j.tr <- v.tr[, j]
            v2.not.v1.j.tr <- v.tr[, colnames(v) != colnames(v)[j], drop = FALSE]
            v1.j.te <- v.te[, j]
            v2.not.v1.j.te <- v.te[, colnames(v) != colnames(v)[j], drop = FALSE]


            w.long.test <- cbind(v1j = v1.j.te,
                                 v2.not.v1.j.te[rep(1:n.te, n.te), , drop = FALSE])
            w.long <- cbind(v1j = v1.j.te,
                            v2.not.v1.j.te[rep(1:n.te, n), , drop = FALSE])
            w.tr <- cbind(v1j = v1.j.tr, v2.not.v1.j.tr)
            w.te <-  cbind(v1j = v1.j.te, v2.not.v1.j.te)

            cond.dens.vals <- cond.dens(v1 = v1.j.tr,
                                        v2 = v2.not.v1.j.tr,
                                        new.v1 = c(v1.j.te, w.long.test[, 1]),
                                        new.v2 = rbind(v2.not.v1.j.te,
                                                       w.long.test[, -1, drop = FALSE]))
            marg.dens <- colMeans(matrix(cond.dens.vals[-c(1:n.te)], ncol = n.te,
                                         nrow = n.te))
            ghat <- marg.dens / cond.dens.vals[1:n.te]

            cate.preds <- cate.w(tau = mu1hat.vals[-c(1:n.te)] - mu0hat.vals[-c(1:n.te)],
                                 w = w.tr,
                                 new.w = rbind(w.te, w.long.test, w.long))

            tauhat.w <- cate.preds$res[1:n.te]
            theta.mat.test <- matrix(cate.preds$res[(n.te+1):(n.te + n.te^2)],
                                     nrow = n.te, ncol = n.te)

            theta.mat[[alg]][, test.idx, j] <-
              matrix(cate.preds$res[-c(1:(n.te + n.te^2))], nrow = n, ncol = n.te)
            theta.bar[[alg]][test.idx, j] <- apply(theta.mat.test, 2, mean)
            pseudo.y.pd[[alg]][test.idx, j] <- (pseudo - tauhat.w) * ghat +
              theta.bar[[alg]][test.idx, j]

            stage2.reg.data.pd[[j]][[k]] <- cbind(data.frame(pseudo.pd = pseudo.y.pd[[alg]][test.idx, j],
                                                             pseudo.cate = pseudo,
                                                             mu1hat = mu1hat,
                                                             mu0hat = mu0hat,
                                                             pihat = pihat,
                                                             tauhat.w = tauhat.w,
                                                             marg.dens = marg.dens,
                                                             cond.dens = cond.dens.vals[1:n.te],
                                                             theta.bar = theta.bar[[alg]][test.idx, j],
                                                             y = y.te,
                                                             a = a.te,
                                                             fold.id = k), v.te)

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
  for(alg in learner) {

    if(alg != "dr") next

    if(additive_approx) {
      additive_model <- drl.basis(y = pseudo.y[[alg]][, 1], x = v, new.x = v)
      tt <- delete.response(terms(additive_model$model))
    }

    for(j in 1:ncol(v)){
      vj <- v[, j]

      if(length(unique(vj)) < 10 | class(vj) == "factor") {
        univariate_debias_res[[alg]][[j]] <- univariate_res[[alg]][[j]] <-
          list(data = data.frame(pseudo = pseudo.y[[alg]][, 1],
                                 exposure = vj),
               res = lm.discrete.v(y = pseudo.y[[alg]][, 1], x = as.factor(vj),
                             new.x = unique(as.factor(vj)))
          )
        if(partial_dependence) {
          pd_debias_res[[alg]][[j]] <- pd_res[[alg]][[j]] <-
            list(data = data.frame(pseudo = pseudo.y.pd[[alg]][, j],
                                   exposure = vj),
                 res = lm.discrete.v(y = pseudo.y.pd[[alg]][, j], x = as.factor(vj),
                               new.x = unique(as.factor(vj)))
            )

        }
      }
      else {
        univariate_debias_res[[alg]][[j]] <-
          list(
            data = data.frame(pseudo = pseudo.y[[alg]][, 1],
                              exposure = vj),
            res = debiased_inference(A = vj, pseudo.out = pseudo.y[[alg]][, 1], tau = 1,
                               muhat.mat = matrix(0, nrow(data)^2),
                               mhat.obs = rep(0, nrow(data)),
                               debias = TRUE,
                               bandwidth.method = "LOOCV",
                               kernel.type = "epa",
                               bw.seq = params[["bw.stage2"]][[j]])
          )
        univariate_res[[alg]][[j]] <-
          list(data = data.frame(pseudo = pseudo.y[[alg]][, 1],
                          exposure = vj),
               res = debiased_inference(A = vj, pseudo.out = pseudo.y[[alg]][, 1], tau = 1,
                                  muhat.mat = matrix(0, nrow(data)^2),
                                  mhat.obs = rep(0, nrow(data)),
                                  debias = FALSE,
                                  bandwidth.method = "LOOCV",
                                  kernel.type = "epa",
                                  bw.seq = params[["bw.stage2"]][[j]])
          )

        if(partial_dependence) {
          pd_debias_res[[alg]][[j]] <-
            list(
              data = data.frame(pseudo = pseudo.y.pd[[alg]][, j],
                                exposure = vj),
              res = debiased_inference(A = vj, pseudo.out = pseudo.y.pd[[alg]][, j], tau = 1,
                                 muhat.mat = theta.mat[[alg]][, , j],
                                 mhat.obs = theta.bar[[alg]][, j],
                                 debias = TRUE,
                                 bandwidth.method = "LOOCV",
                                 kernel.type = "epa",
                                 bw.seq = params[["bw.stage2"]][[j]])
            )
          pd_res[[alg]][[j]] <-
            list(
              data = data.frame(pseudo = pseudo.y.pd[[alg]][, j],
                                exposure = vj),
              res = debiased_inference(A = vj, pseudo.out = pseudo.y.pd[[alg]][, j], tau = 1,
                                 muhat.mat = theta.mat[[alg]][, , j],
                                 mhat.obs = theta.bar[[alg]][, j],
                                 debias = FALSE,
                                 bandwidth.method = "LOOCV",
                                 kernel.type = "epa",
                                 bw.seq = params[["bw.stage2"]][[j]])
            )
        }
      }
      if(additive_approx){
        new.dat.additive <- as.data.frame(matrix(0, nrow = length(v0[[j]]), ncol = ncol(v),
                                                 dimnames = list(NULL, colnames(v))))
        new.dat.additive[, j] <- v0[[j]]
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
              univariate_debias_res = univariate_debias_res,
              pd_res = pd_res,
              pd_debias_res = pd_debias_res,
              additive_res = additive_res)
  return(ret)
}

