#' unif_ctff_series
#'
#' This function calculates the pseudo-outcomes depending on the estimand
#' @param y outcome
#' @param a treatment
#' @param x confounders
#' @param v1 effect modifiers 1
#' @param v2 effect modifiers 2
#' @param mu0x function for fitting mu0
#' @param mu1x function for fitting mu1
#' @param pix function for fitting propensity score
#' @param cond.dens function for the conditional density of any of random variable
#' in v1 given all variables in v2
#' @param cate.w function for fitting E(tau(X) | v2)
#' @param nsplits number of splints
#' @return a list containing the pseudo outcomes
#' @export
get.pseudo.y <- function(y, a, x, v1, v2, mu0.x, mu1.x, pi.x, cond.dens,
                         cate.w, nsplits = 1) {
  
  n <- length(y)
  s <- sample(rep(1:nsplits, ceiling(n / nsplits))[1:n])
  mu0hat <- mu1hat <- pihat <- tauhat <- cate.out <- rep(NA, n)
  cate.out.folds <- pd.out.folds <- vector("list", nsplits)
  
  if(!is.null(v2)) {
    theta.bar <- ghat <- tauhat.w <- rep(NA, n)
    theta.mat <- matrix(NA, ncol = n, nrow = n)
    pd.out <- matrix(NA, ncol = ncol(v1), nrow = n,
                     dimnames = list(NULL, colnames(v1)))
  }
  else pd.out <- theta.bar <- theta.mat <- NULL
  
  for(k in 1:nsplits) {
    
    test.idx <- k == s
    train.idx <- k != s
    
    if(all(!train.idx)) train.idx <- test.idx
    
    n.te <- sum(test.idx)
    n.tr <- sum(train.idx)
    
    a.tr <- a[train.idx]
    y.tr <- y[train.idx]
    x.tr <- x[train.idx, , drop = FALSE]
    
    a.te <- a[test.idx]
    y.te <- y[test.idx]
    x.te <- x[test.idx, , drop = FALSE]
    
    mu0hat.vals <- mu0.x(y = y.tr, a = a.tr, x = x.tr,
                         new.x = rbind(x.te, x.tr))
    mu0hat[test.idx] <- mu0hat.vals[1:n.te]
    
    mu1hat.vals <- mu1.x(y = y.tr, a = a.tr, x = x.tr,
                         new.x = rbind(x.te, x.tr))
    mu1hat[test.idx] <- mu1hat.vals[1:n.te]
    
    pihat[test.idx] <- pi.x(a = a.tr, x = x.tr, new.x = x.te)
    
    tauhat[test.idx] <- mu1hat[test.idx] - mu0hat[test.idx]
    
    cate.out[test.idx] <- a.te / pihat[test.idx] * (y.te - mu1hat[test.idx]) -
      (1-a.te) / (1 - pihat[test.idx]) * (y.te - mu0hat[test.idx]) +
      tauhat[test.idx]
    
    cate.out.folds[[k]] <- cate.out[test.idx]
    
    if(!is.null(v2)) {
      
      for(j in 1:ncol(v1)) {
        v1.j <- v1[, j]
        v2.not.v1.j <- v2[, colnames(v2) != colnames(v1)[j], drop = FALSE]
        v1.tr <- v1.j[train.idx]
        v2.tr <- v2[train.idx, , drop = FALSE]
        v2.not.v1.j.tr <- v2.not.v1.j[train.idx, , drop = FALSE]
        
        v1.te <- v1.j[test.idx]
        v2.te <- v2[test.idx, , drop = FALSE]
        v2.not.v1.j.te <- v2.not.v1.j[test.idx, , drop = FALSE]
        
        w.long.test <- cbind(v1[test.idx, j, drop = FALSE],
                             v2.not.v1.j.te[rep(1:n.te, n.te), , drop = FALSE])
        
        cond.dens.vals <- cond.dens(v1 = v1.tr,
                                    v2 = v2.not.v1.j.tr,
                                    new.v1 = c(v1.te, w.long.test[, 1]),
                                    new.v2 = rbind(v2.not.v1.j.te,
                                                   w.long.test[, -1, drop = FALSE]))
        marg.dens <- colMeans(matrix(cond.dens.vals[-c(1:n.te)], ncol = n.te,
                                     nrow = n.te))
        ghat[test.idx] <- marg.dens / cond.dens.vals[1:n.te]
        
        w.long <- cbind(v1[test.idx, j, drop = FALSE],
                        v2.not.v1.j.te[rep(1:n.te, n), , drop = FALSE])
        
        cate.preds <- cate.w(tau = mu1hat.vals[-c(1:n.te)] - mu0hat.vals[-c(1:n.te)],
                             w = v2.tr,
                             new.w = rbind(v2.te, w.long.test, w.long))
        
        tauhat.w[test.idx] <- cate.preds$res[1:n.te]
        theta.mat.test <- matrix(cate.preds$res[(n.te+1):(n.te + n.te^2)],
                                 nrow = n.te, ncol = n.te)
        theta.mat[, test.idx] <- matrix(cate.preds$res[-c(1:(n.te + n.te^2))],
                                        nrow = n, ncol = n.te)
        theta.bar[test.idx] <- apply(theta.mat.test, 2, mean)
        
        pd.out[test.idx, j] <-
          (cate.out[test.idx] - tauhat.w[test.idx]) * ghat[test.idx] +
          theta.bar[test.idx]
        
      }
      pd.out.folds[[k]] <- pd.out[test.idx, ]
    }
  }
  
  return(list(cate.out = cate.out, pd.out = pd.out, theta.mat = theta.mat,
              theta.bar = theta.bar,
              cate.out.folds = cate.out.folds,
              pd.out.folds = pd.out.folds))
}