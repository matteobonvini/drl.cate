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

cate <- function(v0, learner, y, a, x, v, nsplits = 5, foldid = NULL, ...) {
  
  params <- list(...)
  
  n <- length(y)
  
  if(is.null(foldid)) {
    s <- sample(rep(1:nsplits, ceiling(n / nsplits))[1:n]) 
  } else {
    s <- foldid
    nsplits <- length(unique(foldid))
  }
  
  est <- replicate(length(learner), matrix(NA, ncol = nsplits, nrow = nrow(v0)),
                   simplify = FALSE)
  names(est) <- learner
  if(any(learner == "lp-r")) {
    est[["lp-r"]] <- matrix(NA, ncol = 3, nrow = nrow(v0))
  }
  if(any(learner == "t") & all(colnames(x) %in% colnames(v))) {
    est[["t"]] <- matrix(NA, ncol = 1, nrow = nrow(v0))
  }
  
  pi.x <- params[["pi.x"]]
  pi.x.method <- params[["pi.x.method"]]

  if(is.null(pi.x) & !is.null(pi.x.method)) {
    
    if(pi.x.method == "lasso") pi.x <- pi.x.lasso
    else if(pi.x.method == "SL") pi.x <- pi.x.SL(args[["sl.lib.pi"]])
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
    
    if(is.null(mu1.x)) {
      
      if(mu1.x.method == "lasso") {
        mu1.x <- mu1.x.lasso
      } else if(mu1.x.method == "SL") {
        mu1.x <- mu1.x.SL(sl.lib)
      } else stop("Provide valid method for estimating the E(Y|A=1,X).")
      
    }
    
    if(is.null(mu0.x)) {
      
      if(mu0.x.method == "lasso") {
        mu0.x <- mu0.x.lasso
      } else if(mu0.x.method == "SL") {
        mu0.x <- mu0.x.SL
      } else stop("Provide valid method for estimating the E(Y|A=0,X).")
      
    } 
    if(is.null(drl) & any(learner == "dr")) {
      
      if(drl.method == "lasso") {
        drl <- drl.lasso
      } else if(drl.method == "SL") {
        drl <- drl.SL
      } else stop("Provide valid method for second-stage regression.")
      
    } 
    
  }
  
  if(any(learner %in% c("u", "r", "lp-r"))) {
    
    mu.x.method <- args[["mu.x.method"]]
    
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
    v.te <- v[test.idx, , drop = FALSE]
    
    pihat <- pi.x(a = a.tr, x = x.tr, new.x = x.te)
    
    if(any(learner %in% c("dr", "t"))) {
      
      mu1hat <- mu1.x(y = y.tr, a.tr = a.tr, x = x.tr, new.x = x.te)
      mu0hat <- mu0.x(y = y.tr, a.tr = a.tr, x = x.tr, new.x = x.te)
      
    }
    
    if(any(learner %in% c("u", "r"))) {
      
      muhat <- mu.x(y = y.tr, x = x.tr, new.x = x.te)
      
    }
    
    for(alg in learner) {
      
      if(alg == "dr") {
        
        pseudo <- (a.te - pihat) / (pihat * (1 - pihat)) *
          (y.te - a.te * mu1hat - (1 - a.te) * mu0hat) + mu1hat - mu0hat
        est[[alg]][, k] <- drl(y.tr = pseudo, x.tr = v.te, new.x = v0)
        
      } else if(alg == "t"){
        if(all(colnames(x) %in% colnames(v))) {
          next # no need to smooth the t-learner if V = X
        } else {
          pseudo <- mu1hat - mu0hat
          est[[alg]][, k] <- tl(y.tr = pseudo, x.tr = v.te, new.x = v0)
        }
      } else if(alg == "u") {
        
        pseudo <- (y.te - muhat) / (a.te - pihat)
        est[[alg]][, k] <- ul(y.tr = pseudo, x.tr = v.te, new.x = v0)
        
      } else if(alg == "lp-r") next
    }
  }
  
  if(any(learner == "lp-r")) {
    
    est[["lp-r"]] <- lp_r_learner(x0 = v0, y = y, a = a, x = x, mu.x = mu.x,
                                  pi.x = pi.x, basis = args[["basis"]], 
                                  order_basis = args[["order_basis"]], 
                                  kernel = args[["kernel"]])$fold_est
    
  }
  
  if(any(learner == "t") & all(colnames(x) %in% colnames(v))) {
    est[["t"]][, 1] <- mu1.x(y.tr = y, a.tr = a, x = x, new.x = v0) -
      mu0.x(y.tr = y, a.tr = a, x = x, new.x = v0)
  }
  
  out <- sapply(learner, function(w) apply(est[[w]], 1, mean))
  ret <- list(est = out, fold_est = est)
  return(ret)
}

