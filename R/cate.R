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
                 nsplits = 5, foldid = NULL, ...) {

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
  pseudo.y <- replicate(length(learner), matrix(NA, ncol = 1, nrow = n),
                        simplify = FALSE)
  ites_v <- replicate(length(learner), matrix(NA, ncol = 3, nrow = n),
                      simplify = FALSE)
  ites_x <- replicate(length(learner), matrix(NA, ncol = 1, nrow = n),
                      simplify = FALSE)
  names(est) <- names(est.pi) <- names( pseudo.y) <- names(ites_v) <-
    names(ites_x) <- learner

  stage2.reg.data <- vector("list", nsplits)
  stage2.reg.preds <- vector("list", nsplits)

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
    drl.method <- params[["drl.method"]]
    mu1.x <- params[["mu1.x"]]
    mu0.x <- params[["mu0.x"]]
    drl <- params[["drl"]]
    drl.ite <- params[["drl.ite"]]
    drl.ite.method <- params[["drl.ite.method"]]

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
    
    if(is.null(drl) & any(learner == "dr")) {

      if(drl.method == "lasso") {
        drl <- drl.lasso
      } else if(drl.method == "lm"){
        drl <- drl.lm
      } else if(drl.method == "glm"){
        drl <- drl.glm
      } else if(drl.method == "gam") {
        drl <- drl.gam
      } else if(drl.method == "SL") {
        drl <- drl.SL
      } else stop("Provide valid method for second-stage regression.")

    }
    
    if(is.null(drl.ite) & any(learner == "dr")) {
      if(is.null(drl.ite.method)) {
        # by default use lasso
        drl.ite <- drl.lasso 
      } else if (drl.ite.method == "lm"){
        drl.ite <- drl.ite.lm
      } else if (drl.ite.method == "lasso"){
        drl.ite <- drl.lasso
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

    x.tr <- x[train.idx, , drop = FALSE]
    a.tr <- a[train.idx]
    y.tr <- y[train.idx]
    x.te <- x[test.idx, , drop = FALSE]
    a.te <- a[test.idx]
    y.te <- y[test.idx]

    if ((!is.matrix(v)) & (!is.data.frame(v))){
      v.te <- v[test.idx]
    } else{
      v.te <- v[test.idx, , drop = FALSE]}

    pihat <- pi.x(a = a.tr, x = x.tr, new.x = x.te)

    if(any(learner %in% c("dr", "t"))) {

      mu1hat <- mu1.x(y = y.tr, a = a.tr, x = x.tr, new.x = x.te)
      mu0hat <- mu0.x(y = y.tr, a = a.tr, x = x.tr, new.x = x.te)

    }

    if(any(learner %in% c("u", "r"))) {

      muhat <- mu.x(y = y.tr, x = x.tr, new.x = x.te)

    }

    for(alg in learner) {

      if(alg == "dr") {

        pseudo <- (a.te - pihat) / (pihat * (1 - pihat)) *
          (y.te - a.te * mu1hat - (1 - a.te) * mu0hat) + mu1hat - mu0hat
        drl.res <-  drl(y = pseudo, x = v.te, new.x = rbind(v0.long, v.te))
        drl.res.pi <-  drl(y = mu1hat - mu0hat, x = v.te,
                           new.x = rbind(v0.long, v.te))
        stage2.reg.data[[k]] <- cbind(data.frame(pseudo = pseudo,
                                                 mu1hat = mu1hat,
                                                 mu0hat = mu0hat,
                                                 pihat = pihat,
                                                 y = y.te,
                                                 a = a.te), v.te)
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
        ites_x[[alg]][test.idx, 1] <- drl.ite(y = pseudo, x = x.te, new.x = x.te)$res[,1]

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
              stage2.reg.data = stage2.reg.data)
  return(ret)
}

