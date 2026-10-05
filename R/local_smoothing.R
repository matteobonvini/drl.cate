# Shared bandwidth and numerical safeguards for local-polynomial estimation.
.validate.min.local <- function(min.local) {
  if(!is.null(min.local) && (!is.numeric(min.local) || length(min.local)!=1L ||
     !is.finite(min.local) || min.local<1 || min.local!=floor(min.local))) {
    stop("min.local must be NULL or a positive integer")
  }
}

.local.bandwidth <- function(x, a, bandwidth, min.local=NULL) {
  if(is.null(min.local)) return(bandwidth)
  if(length(x)<min.local) return(NA_real_)
  radius <- sort(abs(x-a), partial=min.local)[min.local]
  # Strictly inside the window: triangular/Epanechnikov weights vanish at 1.
  if(radius>=bandwidth) radius*(1+sqrt(.Machine$double.eps)) else bandwidth
}

.local.system <- function(X, w) {
  if(any(!is.finite(X)) || any(!is.finite(w)) || any(w<0)) return(NULL)
  keep <- w>0
  if(sum(keep)<ncol(X)) return(NULL)
  q <- qr(X[keep,,drop=FALSE]*sqrt(w[keep]), tol=1e-10)
  if(q$rank<ncol(X)) return(NULL)
  R <- qr.R(q)
  if(rcond(R)<1e-10) return(NULL)
  inv <- chol2inv(R)
  inv[q$pivot,q$pivot] <- inv
  inv
}

.local.fit <- function(x,y,a,h,b,kern,debias,min.local=NULL) {
  hh <- .local.bandwidth(x,a,h,min.local)
  bb <- if(debias) .local.bandwidth(x,a,b,min.local) else b
  if(!is.finite(hh) || !is.finite(bb)) return(NULL)
  u <- (x-a)/hh; wh <- kern(u)/hh; X <- cbind(1,u)
  ih <- .local.system(X,wh)
  if(is.null(ih)) return(NULL)
  alpha <- stats::lm.wfit(X,y,wh,tol=1e-10)$coefficients
  out <- list(theta=alpha[1],h=hh,b=bb,alpha=alpha,ih=ih,
              hat=kern(0)/hh*ih[1,1])
  if(debias) {
    z <- (x-a)/bb; wb <- kern(z)/bb; B <- cbind(1,z,z^2,z^3)
    ib <- .local.system(B,wb)
    if(is.null(ib)) return(NULL)
    beta <- stats::lm.wfit(B,y,wb,tol=1e-10)$coefficients
    c2 <- stats::lm.wfit(X,u^2,wh,tol=1e-10)$coefficients[1]
    out$theta <- alpha[1]-(hh/bb)^2*c2*beta[3]
    out$hat <- out$hat-(hh/bb)^2*c2*kern(0)/bb*ib[3,1]
    out$beta <- beta; out$c2 <- c2; out$ib <- ib
  }
  if(!is.finite(out$theta)) return(NULL)
  out
}

# Exact componentwise deletion; repeated modifier values share a system.
.local.loo <- function(x,y,h,b,debias,kernel.type="epa",min.local=NULL) {
  .validate.min.local(min.local)
  kern <- function(u) .kern(u,kernel.type)
  pred <- rep(NA_real_,length(x))
  for(a in unique(x)) {
    ids <- which(x==a)
    hh <- .local.bandwidth(x[-ids[1]],a,h,min.local)
    bb <- if(debias) .local.bandwidth(x[-ids[1]],a,b,min.local) else b
    if(!is.finite(hh) || !is.finite(bb)) next
    fit <- .local.fit(x,y,a,hh,bb,kern,debias)
    safe <- !is.null(fit)
    if(safe) {
      lh <- kern(0)/hh*fit$ih[1,1]
      safe <- is.finite(lh) && 1-lh>1e-9 && rcond(fit$ih)>1e-9
      if(debias) {
        lb <- kern(0)/bb*fit$ib[1,1]
        qb <- kern(0)/bb*fit$ib[3,1]
        safe <- safe && is.finite(lb) && 1-lb>1e-9 && rcond(fit$ib)>1e-9
      }
    }
    if(safe) {
      pred[ids] <- (fit$alpha[1]-lh*y[ids])/(1-lh)
      if(debias) pred[ids] <- pred[ids]-(hh/bb)^2*fit$c2/(1-lh)*
        (fit$beta[3]-qb*(y[ids]-fit$beta[1])/(1-lb))
    } else {
      for(i in ids) {
        # Rebuild the small weighted QR systems, not any nuisance model.
        deleted <- .local.fit(x[-i],y[-i],a,hh,bb,kern,debias)
        if(!is.null(deleted)) pred[i] <- deleted$theta
      }
    }
  }
  list(pred=pred, risk=if(any(!is.finite(pred))) Inf else mean((y-pred)^2))
}
