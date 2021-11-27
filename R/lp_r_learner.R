#' Lp-R-Learner
#'
#' This function computes the Lp-R-Learner approach to CATE estimation.
#' @param x0 evaluation points, i.e. E(Y^1 - Y^0 | x0)
#' @param y vector of outcomes
#' @param a vector of treatments
#' @param x matrix of covariates
#' @param mu.x a function with arguments y, x, new.x computing the regression
#' of y on x and evaluating it at new.x.
#' @param pi.x a function with arguments a, x, new.x computing
#' the propensity score and evaluating it at new.x.
#' @param basis a function with arguments x and j returning the j^th basis element
#' applied to x, e.g. x^j. It will be the building block to compute a tensor product basis.
#' @param order_basis the order of the basis
#' @param kernel a function with arguments x and x_0 returning K((x - x0) / h) / h
#'
#' @return a list containing the following components:
#' \item{est}{ estimate of the CATE at x0}
#' \item{fold_k_est}{ length(x0)xnsplits matrix with estimates of the CATE at x0 in each fold.}
#' @export

lp_r_learner <- function(x0, y, a, x, mu.x, pi.x, basis, order_basis, kernel) {

  n <- length(y)
  x <- as.data.frame(x)
  d <- ncol(x)
  colnames(x) <- paste0("X", 1:d)
  x0 <- as.data.frame(x0)
  colnames(x0) <- colnames(x)

  splits.mat <- matrix(c(1, 2, 3, 2, 3, 1, 3, 1, 2), ncol = 3, nrow = 3,
                       byrow = TRUE)
  nsplits <- 3

  s <- sample(rep(1:nsplits, ceiling(n / nsplits))[1:n])

  est <- matrix(NA, ncol = nsplits, nrow = nrow(x0))

  for(k in 1:nsplits) {

    tr.idx0 <- s == splits.mat[k, 1]
    tr.idx1 <- s == splits.mat[k, 2]
    te.idx <- s == splits.mat[k, 3]

    x.tr0 <- x[tr.idx0, , drop = FALSE]
    x.tr1 <- x[tr.idx1, , drop = FALSE]
    x.te <- x[te.idx, , drop = FALSE]

    y.tr0 <- y[tr.idx0]
    y.tr1 <- y[tr.idx1]
    y.te <- y[te.idx]

    a.tr0 <- a[tr.idx0]
    a.tr1 <- a[tr.idx1]
    a.te <- a[te.idx]

    n.te <- nrow(x.te)


    muhat <- mu.x(y = y.tr0, x = x.tr0, new.x = x.te)
    pi0hat <- pi.x(a = a.tr0, x = x.tr0, new.x = x.te)
    pi1hat <- pi.x(a = a.tr1, x = x.tr1, new.x = x.te)

    for(j in 1:nrow(x0)) {

      x0.mat <- matrix(rep(as.numeric(x0[j, ]), n.te), ncol = ncol(x0), nrow = n.te,
                       byrow = TRUE)

      ww <- (a.te - pi0hat) / (a.te - pi1hat) * kernel(x = x.te, x0 = x0[j, ])

      pseudo.y <- y.te - muhat

      # create tensor product basis
      new.expand.grid <- function(input, reps) {
        expand.grid(replicate(reps, input, simplify = FALSE))
      }

      tensor.idx <- new.expand.grid(1:(order_basis + 1), d)
      basis.x <- matrix(NA, nrow = n.te, ncol = (order_basis + 1)^d)
      basis.additive <- array(dim = c(n.te, d, order_basis + 1))
      basis.0 <- rep(NA, (order_basis + 1)^d)
      basis.additive.0 <- array(dim = c(1, d, order_basis + 1))
      for(b in 1:(order_basis + 1)) {
        basis.additive[, , b] <- apply(as.matrix(x.te) - x0.mat, 2, basis, j = b)
        basis.additive.0[1, , b] <- basis(0, j = b)
      }
      for(tb in 1:(order_basis + 1)^d) {
        basis.x[, tb] <- apply(sapply(1:d,
                                      function(l) basis.additive[, l, tensor.idx[tb, l]],
                                      simplify = TRUE), 1, prod)
        tmp <- sapply(1:d,
                      function(l) basis.additive.0[, l, tensor.idx[tb, l]],
                      simplify = TRUE)
        basis.0[tb] <- prod(as.matrix(tmp, nrow = 1, ncol = d))
      }

      pseudo.x <- matrix(rep((a.te - pi1hat), (order_basis + 1)^d), ncol = (order_basis + 1)^d) * basis.x
      est[j, k] <- sum(basis.0 * coef(lm(pseudo.y ~ -1 + pseudo.x, weights = ww)))

    }
  }

  out <- list(est = apply(est, 1, mean), fold_est = est)

  return(out)
}
