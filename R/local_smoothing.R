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
.local.loo <- function(x,y,h,b,debias,kernel.type="epa",min.local=NULL,
                       eval.indices=seq_along(x)) {
  .validate.min.local(min.local)
  if(!is.numeric(eval.indices) || !length(eval.indices) ||
     anyNA(eval.indices) || any(!is.finite(eval.indices)) ||
     any(eval.indices != floor(eval.indices)) ||
     any(eval.indices<1 | eval.indices>length(x)) || anyDuplicated(eval.indices))
    stop("eval.indices must be distinct valid observation indices")
  kern <- function(u) .kern(u,kernel.type)
  pred <- rep(NA_real_,length(eval.indices))
  for(a in unique(x[eval.indices])) {
    positions <- which(x[eval.indices]==a)
    ids <- eval.indices[positions]
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
      pred[positions] <- (fit$alpha[1]-lh*y[ids])/(1-lh)
      if(debias) pred[positions] <- pred[positions]-(hh/bb)^2*fit$c2/(1-lh)*
        (fit$beta[3]-qb*(y[ids]-fit$beta[1])/(1-lb))
    } else {
      for(pos in positions) {
        i <- eval.indices[pos]
        # Rebuild the small weighted QR systems, not any nuisance model.
        deleted <- .local.fit(x[-i],y[-i],a,hh,bb,kern,debias)
        if(!is.null(deleted)) pred[pos] <- deleted$theta
      }
    }
  }
  list(pred=pred, risk=if(any(!is.finite(pred))) Inf else mean((y[eval.indices]-pred)^2))
}

.validate.local.cv <- function(cv.eval.size, cv.repeats) {
  integer.scalar <- function(x) is.numeric(x) && length(x)==1L &&
    is.finite(x) && x>=1 && x==floor(x)
  if(!is.null(cv.eval.size) && !integer.scalar(cv.eval.size))
    stop("cv.eval.size must be NULL or a positive integer")
  if(!integer.scalar(cv.repeats)) stop("cv.repeats must be a positive integer")
}

# Sample losses, never the training observations. Cache predictions for the
# union of repeat indices; repeats may overlap, but each is sampled without replacement.
.local.cv <- function(x,y,h,b,debias,kernel,min.local,cv.eval.size,cv.repeats) {
  .validate.local.cv(cv.eval.size,cv.repeats)
  exhaustive <- is.null(cv.eval.size) || cv.eval.size>=length(x)
  indices <- if(exhaustive) list(seq_along(x)) else
    replicate(cv.repeats,sample.int(length(x),cv.eval.size),simplify=FALSE)
  union.indices <- unique(unlist(indices,use.names=FALSE))
  positions <- lapply(indices,match,table=union.indices)
  risks <- matrix(Inf,length(h),length(indices))
  for(k in seq_along(h)) {
    prediction <- .local.loo(x,y,h[k],b[k],debias,kernel,min.local,
                             eval.indices=union.indices)$pred
    loss <- (y[union.indices]-prediction)^2
    for(r in seq_along(indices)) {
      v <- loss[positions[[r]]]
      if(all(is.finite(v))) risks[k,r] <- mean(v)
    }
  }
  choose <- function(risk) {
    if(!any(is.finite(risk)))
      stop("All bandwidth candidates have failed leave-one-out predictions in a CV repeat; increase bw.seq or min.local.")
    tied <- which(is.finite(risk) & abs(risk-min(risk))<1e-6)
    tied[order(-h[tied],-b[tied])[1]]
  }
  winners <- vapply(seq_along(indices),function(r) choose(risks[,r]),integer(1))
  pooled <- rowMeans(risks)
  median.h <- stats::median(h[winners]); median.b <- stats::median(b[winners])
  if(exhaustive || length(indices)==1L) selected <- winners[1] else {
    valid <- which(is.finite(pooled))
    if(!length(valid)) stop("No bandwidth candidate has finite risk in every CV repeat; increase bw.seq or min.local.")
    # Log distances make the projection invariant to the modifier's units.
    distance <- log(h[valid]/median.h)^2 + log(b[valid]/median.b)^2
    nearest <- valid[abs(distance-min(distance))<1e-12]
    selected <- nearest[order(pooled[nearest],-h[nearest],-b[nearest])[1]]
  }
  list(risk=pooled,selected=selected,diagnostics=list(
    exhaustive=exhaustive,indices=indices,risk=risks,
    candidates=data.frame(h=h,b=b),
    winners=data.frame(candidate=winners,h=h[winners],b=b[winners]),
    median=c(h=median.h,b=median.b),selected=selected,
    selection=if(exhaustive) "exhaustive" else "median-bandwidth"))
}
