# nuisance functions are now all using cv.glmnet
# drlasso function: dr-learner with lasso

mu.x.new <- function(y.tr, a.tr, x.tr, new.x, nfolds, alpha = 1, treatment = NULL, ...){
  # estimate E(Y|X)
  
  # y.tr: vector of outcomes
  # a.tr: vector of binary treatments
  # x.tr: matrix of confounders
  # new.x: estimation point
  # nfolds: a number indicating the number of splits used to do cross-fitting.
  # alpha: hyper-parameter in glmnet. alpha = 1 indicates lasso regression
  # treatment: categorical variable. 
  #   E(Y|X, W=1) is estimated when treatment = 1, E(Y|X, W=0) is estimated when 0,  E(Y|X) when NULL
  
  params <- list(...)
  
  if(!is.null(treatment)){
    if(treatment == 0){
      y.tr <- y.tr[a.tr == 0]
      x.tr <- x.tr[a.tr == 0, , drop = FALSE]
    }else if(treatment == 1){
      y.tr <- y.tr[a.tr == 1]
      x.tr <- x.tr[a.tr == 1, , drop = FALSE]
    }
    else stop("invalid treatment")
  }
  
  n <- length(y.tr)
  
  family <- params[["family"]]
  foldid <- params[["foldid"]]
  lambda <- params[["lambda"]]
  if(is.null(family)) family <- 'gaussian'
  if(is.null(lambda)) lambda <- 'lambda.min'
  
  tmp <- as.integer(ceiling(n / nfolds))

  if(is.null(foldid)) foldid <- sample(rep(1:nfolds, tmp)[1:n])
  
  fit <- cv.glmnet(x.tr, y.tr, family = family, nfolds = nfolds, foldid = foldid, alpha = alpha)
  res <- predict(fit, new.x, s= lambda)
  
  return(res)
}


mu.x.new2 <- mu.x.new(dta$y.tr, dta$w.tr, dta$x.tr, new.x = dta$x.te, 
                      nfolds = 3, alpha = 1, treatment = NULL)

mu.x.new1 <- mu.x.new(dta$y.tr, dta$w.tr, dta$x.tr, new.x = dta$x.te, 
                      nfolds = 3, alpha = 1, treatment = 1)

mu.x.new0 <- mu.x.new(dta$y.tr, dta$w.tr, dta$x.tr, new.x = dta$x.te, 
                      nfolds = 3, alpha = 1, treatment = 0)



pi.x.new <- function(a.tr, x.tr, new.x, nfolds, ...){
  # estimate E(W|X)
  
  # a.tr: vector of binary treatments
  # x.tr: matrix of confounders
  # new.x: covariates to be estimated
  # nfolds: a number indicating the number of splits used to do cross-fitting.
  
  params <- list(...)
  n <- length(a.tr)
  # print(n)
  family <- params[["family"]]
  foldid <- params[["foldid"]]
  lambda <- params[["lambda"]]
  
  if(is.null(family)) family <- 'binomial'
  if(is.null(lambda)) lambda <- 'lambda.min'
  if(is.null(foldid)) foldid <- sample(rep(1:nfolds, ceiling(n / nfolds))[1:n])
  
  w_fit <- cv.glmnet(x.tr, a.tr, nfolds = 3,
                     foldid = foldid,
                     family= family)
  # w_lambda_min = w_fit$lambda[which.min(w_fit$cvm[!is.na(colSums(w_fit$fit.preval))])]
  res <- predict(w_fit, new.x, s = lambda, type = "response")
  # p_hat = w_fit$fit.preval[,!is.na(colSums(w_fit$fit.preval))][, w_fit$lambda[!is.na(colSums(w_fit$fit.preval))] == w_lambda_min]
  return(res)
}
pi.x.test <- pi.x.new(a.tr = dta$w.tr, x.tr = dta$x.tr, new.x = dta$x.te, nfolds = 3)
hist(pi.x.test)



drl.ite.new <- function(y.tr, x.tr, new.x, nfolds, ...){
  # used for second stage regression:
  
  # y.tr
  # x.tr
  # new.x
  # nfolds
  
  params <- list(...)
  
  n <- length(y.tr)
  family <- params[["family"]]
  foldid <- params[["foldid"]]
  lambda <- params[["lambda"]]
  
  if(is.null(family)) family <- 'gaussian'
  if(is.null(foldid)) foldid <- sample(rep(1:nfolds, ceiling(n / nfolds))[1:n])
  if(is.null(lambda)) lambda <- 'lambda.min'
  
  fit <- cv.glmnet(x.tr, y.tr, family = family, nfolds = nfolds, foldid = foldid, alpha = 1)
  res <- predict(fit, new.x, s = lambda)
  
  return(res)
}


drl.ite.new(y.tr, x.tr, new.x = x.te, nfolds = 3)



# create a new dr-learner function
drlasso <- function(v0, y, a, x, v, nsplits = 5, nfolds = 3, alpha = 1, ...){
  
  params <- list(...)
  foldid <- params[["foldid"]]
  oracle <- params[["oracle"]]

  n <- length(y)
  
  if(is.null(foldid)) {
    s <- sample(rep(1:nsplits, ceiling(n / nsplits))[1:n])
  } else {
    s <- foldid
    nsplits <- length(unique(foldid))
  }
  # print(s)
  if(is.null(oracle)){oracle <- FALSE}
  
  foldid_res <- matrix(nrow = length(y)/nsplits, ncol = nsplits)
  pred_res <- matrix(nrow = length(y), ncol = nsplits)
  mu1_res <- matrix(nrow = length(y)/nsplits, ncol = nsplits)
  mu0_res <- matrix(nrow = length(y)/nsplits, ncol = nsplits)
  pi_res <- matrix(nrow = length(y)/nsplits, ncol = nsplits)
  pseudo_res <- matrix(nrow = length(y)/nsplits, ncol = nsplits)
  
  for(k in 1:nsplits){
    
    # assign test and train set to calculate mu and pi
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
    
    if(oracle == 'FALSE'){
      mu1hat <- mu.x.new(y.tr = y.tr, a.tr = a.tr, x.tr = x.tr, new.x = x.te, 
                         alpha = 1, nfolds = 3, treatment = 1)
      mu0hat <- mu.x.new(y.tr = y.tr, a.tr = a.tr, x.tr = x.tr, new.x = x.te, 
                         alpha = 1, nfolds = 3, treatment = 0)
      pihat <- pi.x.new(a = a.tr, x = x.tr, new.x = x.te, nfolds = 3, alpha = 1)
    }else
    {
      mu1hat <- true_mu1(x.te)
      mu0hat <- true_mu0(x.te)
      pihat <- e(x.te)
      print('oracle mu and pi!')
      # print(mu1hat)
    }

    pseudo <- (a.te - pihat) / (pihat * (1 - pihat)) *
      (y.te - a.te * mu1hat - (1 - a.te) * mu0hat) + mu1hat - mu0hat

     
    foldid_res[,k] <- as.numeric(which(test.idx == TRUE))
    pseudo_res[,k] <- pseudo
    mu1_res[,k] <- mu1hat
    mu0_res[,k] <- mu0hat
    pi_res[,k] <- pihat
    pred_res[,k] <- drl.ite(y.tr = pseudo, x.tr = v.te, new.x = v0)
  }
  res = list(pred = pred_res, pseudo = pseudo_res, foldid = foldid_res,
             mu1 = mu1_res, mu0 = mu0_res, pi = pi_res)
  
  return(res)
}

# drlasso.test2 <- drlasso(v0 = dta$x.te, y = dta$y.tr, a = dta$w.tr, x = dta$x.tr, 
#                         v = dta$x.tr, nsplits = 5, nfolds = 3, oracle = TRUE)



# calculate MSE of oracle, average MSE ~ 0.01, looks good
tmp_res_oracle <- vector("numeric", 40L)

for (i in 1:40){
  drlasso.test2 <- drlasso(v0 = dta$x.te, y = dta$y.tr, a = dta$w.tr, x = dta$x.tr, 
                           v = dta$x.tr, nsplits = 4, nfolds = 3, oracle = TRUE)
  pred.drl <-rowMeans(drlasso.test2[["pred"]])
  true_est <- true_mu1(dta$x.te) - true_mu0(dta$x.te)
  
  tmp_res_oracle[i] <- mean((pred.drl - true_est)^2)
  print(paste('Iteration', i))
}

hist(tmp_res_oracle)


# calculate MSE of dr-learner
tmp_res <- vector("numeric", 40L)

for (i in 1:40){
  drlasso.test2 <- drlasso(v0 = dta$x.te, y = dta$y.tr, a = dta$w.tr, x = dta$x.tr, 
                           v = dta$x.tr, nsplits = 4, nfolds = 3, oracle = FALSE)
  pred.drl <-rowMeans(drlasso.test2[["pred"]])
  true_est <- true_mu1(dta$x.te) - true_mu0(dta$x.te)
  
  tmp_res[i] <- mean((pred.drl - true_est)^2)
  print(paste('Iteration', i))
  # print(tmp_res)
}

hist(tmp_res)


