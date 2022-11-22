#' unif_ctff_series
#'
#' This function calculates the quantile of sup_v|tauhat(v) - tau(v)| / sigma(v)
#' when tauhat(x) is based on a orthogonal series second stage regression.
#' @param design.mat.v0 design matrix of evaluation points.
#' @param design.mat.v design matrix of observed v points.
#' @param residuals observed residuals, i.e. pseudo.y - fitted.
#' @param B number of bootstrap replications.
#' @param alpha 1-confidence level, e.g. alpha = 0.05.
#' @return a list containing the following components:
#' \item{cutoff}{ the estimated quantile}
#' \item{sigma.hat}{ vector of estimated standard devition of sigma(v)}
#' \item{Omega.hat.sqrt}{ Omega matrix}
#' @export
#' @references Chernozhukov et al (2014).
unif_ctff_series <- function(design.mat.v0, residuals, design.mat, B, alpha) {

  # estimate c_n(1-alpha) such that
  # P(sup_v|tauhat(v) - tau(v)| / sigma(v) <= c_n(1-alpha)) = 1 - alpha - o(1)
  # using the method by Chernozhukov et al (2014) where
  # tauhat is based on orthogonal series regression.
  norm.sq <- function(x) sqrt(sum(x^2))
  k <- ncol(design.mat)
  n <- nrow(design.mat)
  Q.hat <- t(design.mat) %*% design.mat / n
  Q.hat.inv <- solve(Q.hat)
  res.mat <- matrix(residuals, ncol = k, nrow = n, byrow = FALSE)
  Sigma.hat <- t(design.mat * res.mat) %*% (design.mat * res.mat) / n
  Omega.hat <- Q.hat.inv %*% Sigma.hat %*% Q.hat.inv
  Omega.hat.sqrt <- expm::sqrtm(Omega.hat)
  sigma.hat <- apply(Omega.hat.sqrt %*% t(design.mat.v0), 2, norm.sq) / sqrt(n)
  sigma.hat.mat <- matrix(sigma.hat, ncol = k, nrow = nrow(design.mat.v0),
                          byrow = FALSE)
  mult <- design.mat.v0 %*% Omega.hat.sqrt / (sigma.hat.mat * sqrt(n))
  gauss <- matrix(rnorm(k * B), ncol = B, nrow = k)
  res <- apply(abs(mult %*% gauss), 2, max)

  return(list(cutoff = quantile(res, 1-alpha), sigma.hat = sigma.hat,
              Omega.hat.sqrt = Omega.hat.sqrt))
}


#' cate_lvl_set
#'
#' This function estimates the upper level set of the cate
#'
#' @param theta the user-specified level
#' @param cate.obj a list of predictions of the cate at v0, fitted cate values,
#' residuals, design matrix at v0 and design matrix at oberved V.
#' @param set_type string indicating if upper level set, lower level set or level set
#' should be computed. Currently ignored.
#' @param se boolean for whether confidence sets should be returned.
#' @param B number of bootstrap replications, ignored if se = FALSE.
#' @param alpha 1-confidence level, ignored if se = FALSE.
#' @return a list containing the following components:
#' \item{level_set}{ the estimated level sets}
#' \item{chat.l}{ lower confidence set}
#' \item{chat.u}{ upper confidence set}
#' \item{cutoff}{ cutoff from multiplier bootstrap}
#' @export
cate_lvl_set <- function(theta, cate.obj, set_type = "upper",
                         se = TRUE, B = 1000, alpha = 0.05) {
  # it currently works only when effect modifiers are 2D
  cate.vals <- cate.obj$predictions
  v0 <- cate.obj$v0
  v0.size <- nrow(v0)
  cate.vals <- matrix(cate.vals, nrow = v0.size, ncol = v0.size, byrow = FALSE)

  res <- contourLines(x = v0$v1, y = v0$v2, z = cate.vals, levels = theta)

  is.there.lvl.set <- min(cate.vals) <= theta & theta <= max(cate.vals)

  if(is.there.lvl.set) res.binded <- as.matrix(bind_rows(res)[, -1])
  else res.binded <- cbind(x = -1, y = -1)
  if(se) {
  residuals <- cate.obj$residuals
  design.mat <- cate.obj$design.mat
  design.mat.v0 <- cate.obj$design.mat.v0

  res.boot <- unif_ctff_series(design.mat.v0, residuals = residuals,
                               design.mat = design.mat, B = B, alpha = alpha)

  ctf <- res.boot$cutoff
  sigma.hat <- res.boot$sigma.hat

  ctf.l <- matrix(theta + sigma.hat * ctf, byrow = FALSE, ncol = v0.size,
                  nrow = v0.size)
  ctf.u <- matrix(theta - sigma.hat * ctf, byrow = FALSE, ncol = v0.size,
                  nrow = v0.size)
  chat.l <- cate.vals > ctf.l
  chat.u <- cate.vals >= ctf.u

  } else chat.l <- chat.u <- ctf <- NA
  ret <- list(level.set = res.binded, chat.l = chat.l, chat.u = chat.u,
              cutoff = ctf)
  return(ret)

}

