#' DR-Learner
#'
#' This functions computes the DR-Learner approach to CATE estimation.
#' @param x0 evaluation points, i.e. E(Y^1 - Y^0 | x0)
#' @param y vector of outcomes
#' @param a vector of treatments
#' @param x matrix of covariates
#' @param drl.x a function with arguments y, x, new.x computing the regression
#' of y on x and evaluating it at new.x. This is the model for the second-stage
#' regression.
#' @param mu1.x a function with arguments y, x, new.x computing the regression
#' of y on x among the treated and evaluating it at new.x.
#' @param mu0.x a function with arguments y, x, new.x computing the regression
#' of y on x among the un-treated and evaluating it at new.x.
#' @param pi.x a function with arguments a, x, new.x computing
#' the propensity score and evaluating it at new.x.
#' @param nsplits a number indicating the number of splits used to do cross-fitting
#'
#' @return a list containing the following components:
#' \item{est}{ estimate of the CATE at x0}
#' \item{fold_k_est}{ length(x0)xnsplits matrix with estimates of the CATE at x0 in each fold.}
#' @export
#' @examples
#' drl.x <- function(y, x, new.x) {
#' dat <- cbind(data.frame(y = y), x)
#' predict(lm(y ~ ., data = dat), newdata = new.x)
#' }
#' mu1.x <- mu0.x <- drl.x
#' pi.x <- function(a, x, new.x) {
#' dat = cbind(data.frame(a = a), x)
#' predict(glm(a ~ ., data = dat, family = binomial()),
#' newdata = new.x, type = "response")
#' }
#' dr_learner(x0 = cbind(seq(0, 1, 0.1), seq(0, 1, 0.1)),
#'            y = rnorm(100), a = rbinom(100, 1, 0.5),
#'            x = matrix(runif(200), ncol = 2, nrow = 100),
#'            drl.x = drl.x, mu1.x = mu1.x, mu0.x = mu0.x, pi.x = pi.x,
#'            nsplits = 2)
#'
#' @references Kennedy, EH. (2020). Optimal Doubly Robust Estimation of
#' Heterogeneous Causal Effects. \emph{arXiv preprint arXiv:2004.14497}.

dr_learner <- function(x0, y, a, x, drl.x, mu1.x, mu0.x, pi.x, nsplits = 5) {

  n <- length(y)
  x <- as.data.frame(x)
  colnames(x) <- paste0("X", 1:ncol(x))
  x0 <- as.data.frame(x0)
  colnames(x0) <- colnames(x)
  if(nsplits > 1) s <- sample(rep(1:nsplits, ceiling(n / nsplits))[1:n])

  est <- matrix(NA, ncol = nsplits, nrow = nrow(x0))

  for(k in 1:nsplits) {

    test.idx <- k == s
    train.idx <- k != s

    x.tr <- x[train.idx, , drop = FALSE]
    a.tr <- a[train.idx]
    y.tr <- y[train.idx]
    x.te <- x[test.idx, , drop = FALSE]
    a.te <- a[test.idx]
    y.te <- y[test.idx]

    mu1hat <- mu1.x(y = y.tr[a.tr == 1], x = x.tr[a.tr == 1, , drop = FALSE],
                    new.x = x.te)
    mu0hat <- mu0.x(y = y.tr[a.tr == 0], x = x.tr[a.tr == 0, , drop = FALSE],
                    new.x = x.te)
    pihat <- pi.x(a = a.tr, x = x.tr, new.x = x.te)

    pseudo <- (a.te - pihat) / (pihat * (1 - pihat)) *
      (y.te - a.te * mu1hat - (1 - a.te) * mu0hat) + mu1hat - mu0hat

    est[, k] <- drl.x(y = pseudo, x = x.te, new.x = x0)

  }

  out <- list(est = apply(est, 1, mean), fold_k_est = est)

  return(out)
}
