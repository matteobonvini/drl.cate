# Historical proposal validation using the former ten-distinct-value floor.
# Current NULL/numeric min.local behavior is tested in test_min_local.R.
pkgload::load_all('.',quiet=TRUE) # Run from the repository root.

proposed_loo <- function(x,y,h,b,debias=TRUE,kernel='epa',guard=1e-9) {
  stopifnot(length(x)==length(y),all(is.finite(x)),all(is.finite(y)),h>0,b>0)
  pred <- rep(NA_real_,length(x)); fallback <- floor.changed <- 0L
  for(a in unique(x)) {
    ids <- which(x==a)
    # Use the bandwidth that an ACTUAL deleted-sample refit would use.
    floor.full <- sort(abs(unique(x-a)))[10]
    floor.loo <- sort(abs(unique(x[-ids[1]]-a)))[10]
    if(!is.finite(floor.loo)) next
    hh <- max(h,floor.loo); bb <- max(b,floor.loo)
    floor.changed <- floor.changed+length(ids)*as.integer(
      hh!=max(h,floor.full) || (debias && bb!=max(b,floor.full)))
    u <- (x-a)/hh; z <- (x-a)/bb
    wh <- .kern(u,kernel)/hh; wb <- .kern(z,kernel)/bb
    X <- cbind(1,u); B <- cbind(1,z,z^2,z^3)
    Mh <- crossprod(X,X*wh); Mb <- crossprod(B,B*wb)
    stable <- rcond(Mh)>guard && (!debias || rcond(Mb)>guard)
    if(stable) {
      ih <- chol2inv(chol(Mh))
      alpha <- drop(ih %*% crossprod(X,y*wh))
      lambda.h <- wh[ids[1]]*ih[1,1]
      stable <- is.finite(lambda.h) && 1-lambda.h>guard
      if(debias) {
        ib <- chol2inv(chol(Mb))
        beta <- drop(ib %*% crossprod(B,y*wb))
        c2 <- (ih %*% crossprod(X,u^2*wh))[1]
        lambda.b <- wb[ids[1]]*ib[1,1]
        q.b <- wb[ids[1]]*ib[3,1]
        stable <- stable && is.finite(lambda.b) && 1-lambda.b>guard
      }
    }
    for(i in ids) {
      if(stable) {
        pred[i] <- (alpha[1]-lambda.h*y[i])/(1-lambda.h)
        if(debias) pred[i] <- pred[i]-(hh/bb)^2*c2/(1-lambda.h)*
          (beta[3]-q.b*(y[i]-beta[1])/(1-lambda.b))
      } else {
        fallback <- fallback+1L
        # Numerical fallback: direct weighted QR on the deleted sample.
        fit.h <- lm.wfit(X[-i,,drop=FALSE],cbind(y,u^2)[-i,,drop=FALSE],wh[-i])
        if(fit.h$rank<2L) next
        pred[i] <- fit.h$coefficients[1,1]
        if(debias) {
          fit.b <- lm.wfit(B[-i,,drop=FALSE],y[-i],wb[-i])
          if(fit.b$rank<4L) {pred[i]<-NA_real_; next}
          pred[i] <- pred[i]-(hh/bb)^2*fit.h$coefficients[1,2]*fit.b$coefficients[3]
        }
      }
    }
  }
  list(pred=pred,risk=if(any(!is.finite(pred))) Inf else mean((y-pred)^2),
       fallback=fallback,floor.changed=floor.changed)
}

explicit_loo <- function(x,y,h,b,debias,kernel) vapply(seq_along(x),function(i) {
  floor <- sort(abs(unique(x[-i]-x[i])))[10]
  unname(.lprobust(x[-i],y[-i],max(h,floor),max(b,floor),debias,
                  eval.pt=x[i],kernel.type=kernel)[1,'theta.hat'])
},numeric(1))

results <- list(); counter <- 0L
for(seed in 1:5) {
  set.seed(seed)
  designs <- list(normal=sort(rnorm(60)), skewed=rexp(60),
                  ties=rep(seq(-2,2,length.out=15),each=4))
  for(design in names(designs)) for(kernel in c('epa','uni','tri','gau'))
    for(debias in c(FALSE,TRUE)) for(pair in list(c(2,2),c(.8,2),c(2,.8),c(.001,.001))) {
      x <- designs[[design]]; y <- sin(x)+.1*x^3+rnorm(length(x),sd=.3)
      p <- proposed_loo(x,y,pair[1],pair[2],debias,kernel)
      ref <- explicit_loo(x,y,pair[1],pair[2],debias,kernel)
      err <- max(abs(p$pred-ref)/pmax(1,abs(ref)))
      risk.err <- abs(p$risk-mean((y-ref)^2))/max(1,mean((y-ref)^2))
      stopifnot(err<1e-7,risk.err<1e-7)
      counter <- counter+1L
      results[[counter]] <- data.frame(seed,design,kernel,debias,h=pair[1],b=pair[2],err,risk.err,
                                      fallback=p$fallback,floor.changed=p$floor.changed)
    }
}
results <- do.call(rbind,results)
cat('CASES',nrow(results),'MAX PRED ERROR',max(results$err),'MAX RISK ERROR',max(results$risk.err),
    'FLOOR CHANGES',sum(results$floor.changed),'FALLBACKS',sum(results$fallback),'\n')
saveRDS(results,file.path(tempdir(),'cate-loocv-proposal-results.rds'))

# Check every candidate and the selected pair, not just isolated predictions.
set.seed(913)
x <- rnorm(100); y <- sin(2*x)+rnorm(100,sd=.3)
grid <- expand.grid(h=c(.2,.7,1.5),b=c(.2,.7,1.5))
new <- ref <- numeric(nrow(grid))
for(i in seq_len(nrow(grid))) {
  new[i] <- proposed_loo(x,y,grid$h[i],grid$b[i])$risk
  ref[i] <- mean((y-explicit_loo(x,y,grid$h[i],grid$b[i],TRUE,'epa'))^2)
}
stopifnot(max(abs(new-ref))<1e-7,which.min(new)==which.min(ref))
cat('SELECTED PAIR',unlist(grid[which.min(new),]),'\n')

# Extreme leverage exercises the direct-QR fallback.
x <- c(seq(-.1,.1,length.out=59),3)
y <- sin(x)
p <- proposed_loo(x,y,6,6,TRUE,'gau')
stopifnot(p$fallback>0L,all(is.finite(p$pred)))
qr.ref <- vapply(seq_along(x),function(i) {
  u <- (x[-i]-x[i])/6; w <- dnorm(u)/6
  linear <- lm.wfit(cbind(1,u),cbind(y[-i],u^2),w)$coefficients
  cubic <- lm.wfit(cbind(1,u,u^2,u^3),y[-i],w)$coefficients
  linear[1,1]-linear[1,2]*cubic[3]
},numeric(1))
stopifnot(max(abs(p$pred-qr.ref))<1e-6)
cat('HIGH LEVERAGE FALLBACKS',p$fallback,'\n')

# Numerically rank-deficient deleted cubic systems must not get a finite score.
x <- c(seq(-.1,.1,length.out=59),100); y <- sin(x)
p <- proposed_loo(x,y,200,200,TRUE,'gau')
stopifnot(is.infinite(p$risk),p$fallback>0L)
cat('NUMERICALLY RANK-DEFICIENT CANDIDATE REJECTED\n')

# Repeated x values, as in age data: reuse systems for each distinct value.
set.seed(71)
x <- sample(18:100,1000,replace=TRUE); y <- sin(x/15)+rnorm(1000)
t.new <- system.time(p <- proposed_loo(x,y,15,10,TRUE,'gau'))[['elapsed']]
t.ref <- system.time(ref <- explicit_loo(x,y,15,10,TRUE,'gau'))[['elapsed']]
stopifnot(max(abs(p$pred-ref))<1e-7)
cat('N1000 TIMING PROPOSED',t.new,'EXPLICIT',t.ref,'SECONDS\n')

# The two repetitions in the original regression test, fixed seeds.
for(seed in c(1024,1025)) {
  set.seed(seed); x <- rnorm(200,sd=.1); y <- cos(2*pi*x)+rnorm(200)
  p <- proposed_loo(x,y,.95,.95,TRUE,'epa')
  ref <- explicit_loo(x,y,.95,.95,TRUE,'epa')
  stopifnot(abs(p$risk-mean((y-ref)^2))<1e-10)
}
cat('ORIGINAL LOOCV ASSERTIONS PASS WITH PROPOSED CRITERION\n')

# Permutation, response linearity, and modifier translation/rescaling.
set.seed(48); x <- rnorm(80); y <- rnorm(80); y2 <- rnorm(80)
ref <- proposed_loo(x,y,.7,1.2,TRUE,'gau')$pred
ix <- sample(length(x))
permuted <- proposed_loo(x[ix],y[ix],.7,1.2,TRUE,'gau')$pred
linear <- proposed_loo(x,2*y+3*y2,.7,1.2,TRUE,'gau')$pred
scaled <- proposed_loo(5*x+9,y,3.5,6,TRUE,'gau')$pred
stopifnot(max(abs(permuted-ref[ix]))<1e-9,
          max(abs(linear-2*ref-3*proposed_loo(x,y2,.7,1.2,TRUE,'gau')$pred))<1e-9,
          max(abs(scaled-ref))<1e-9)
forced <- proposed_loo(x,y,.7,1.2,TRUE,'gau',guard=1)
stopifnot(forced$fallback==length(x),max(abs(forced$pred-ref))<1e-9)
cat('PERMUTATION, LINEARITY, RESCALING, AND FORCED QR FALLBACK PASS\n')
