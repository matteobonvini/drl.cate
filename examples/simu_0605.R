
# some simulations to compare the performance of r-learner, dr-learner, and other learners.
rm(list = ls())
#import packages
library(rlearner)
library(drl.cate)
require(orthopolynom)
require(SuperLearner)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(reshape2)

# nuisance functions:
sl.lib =  c("SL.glm", "SL.randomForest", "SL.polymars", "SL.mean")

mu1.x <- function(y.tr, a.tr, x.tr, new.x,
                  sl.lib = c("SL.glmnet")){
  y1.tr <- y.tr[a.tr == 1]
  x1.tr <- as.data.frame(x.tr[a.tr == 1, , drop = FALSE])
  fit <- SuperLearner::SuperLearner(Y = y1.tr, X = x1.tr, SL.library = sl.lib,
                                    newX = as.data.frame(new.x))
  return(fit$SL.predict)
}

# mu1.x <- function(y.tr, a.tr, x.tr, new.x, nfolds){
#   y1.tr <- y.tr[a.tr == 1]
#   x1.tr <- x.tr[a.tr == 1, , drop = FALSE]
#   foldid <- sample(rep(seq(nfolds), length.out = nrow(x.tr)))
#   model <- cv.glmnet(Y = y1.tr, X = x1.tr, nfolds = nfolds, foldid = foldid)
#   
#   res <- predict(model, newx = new.x, s = "lambda.min")
#   return(res)
# }
# 
# mu1.x(y.tr = y.tr, a.tr = w.tr, x.tr = x.tr, new.x = x.te, nfolds = 3)



mu0.x <- function(y.tr, a.tr, x.tr, new.x, 
                  sl.lib = c("SL.glmnet") ){
  y0.tr <- y.tr[a.tr == 0]
  x0.tr <- as.data.frame(x.tr[a.tr == 0, , drop = FALSE])
  fit <- SuperLearner::SuperLearner(Y = y0.tr, X = x0.tr, SL.library = sl.lib,
                                    newX = as.data.frame(new.x))
  return(fit$SL.predict)
}

mu.x <- function(y.tr, a.tr, x.tr, new.x, 
                 sl.lib = c("SL.glmnet")){
  fit <- SuperLearner::SuperLearner(Y = y.tr, X = as.data.frame(x.tr), SL.library = sl.lib,
                                    newX = as.data.frame(new.x))
  return(fit$SL.predict)
}

pi.x <- function(a.tr, x.tr, new.x, 
                 sl.lib = c("SL.glmnet")){
  fit <- SuperLearner::SuperLearner(Y = a.tr, X = as.data.frame(x.tr), SL.library = sl.lib,
                                    newX = as.data.frame(new.x))
  return(fit$SL.predict)
}



# drl functions:
drl.x <- function(y.tr, x.tr, new.x, 
                  sl.lib = c("SL.glmnet") ){
  fit <- SuperLearner::SuperLearner(Y = y.tr, X = as.data.frame(x.tr),
                                    SL.library = sl.lib,
                                    newX = as.data.frame(new.x))
  return(fit$SL.predict)
}

# individual treatment effect, same function as drl.x, for second stage regression
# note that y.tr could be pseudo \phi
drl.ite <- function(y.tr, x.tr, new.x, 
                    sl.lib = c("SL.glmnet")){
  fit <- SuperLearner::SuperLearner(Y = y.tr, X = as.data.frame(x.tr), 
                                    SL.library = sl.lib,
                                    newX = as.data.frame(new.x))
  return(fit$SL.predict)
}

# generate data
# n <- number of total data, including train and test set; p <- number of features
n <- 2000; p <- 10; eta <- 0.1;

b <- function(X){
  sin(pi * X[,1] * X[,2]) + 2 * (X[,3] - 0.5)^2 + X[,4] + 0.5 * X[,5]
}

tau <- function(X){
  (X[,1] + X[,2]) / 2
}

e <- function(X){
  pmax(eta, pmin(sin(pi * X[,1] * X[,2]), 1 - eta))
}

gen_data <- function(n, p){
  idx.test <- n/2; idx.train <- idx.test + 1
  
  X <- matrix(runif(n * p, min = 0, max = 1), n, p)
  X_train <- X[idx.train:n,]
  X_test <- X[1:idx.test,]
  
  W <- rbinom(n, 1, e(X))
  W_train <- W[idx.train: n]
  W_test <- W[1: idx.test]
  
  Y <- b(X) + (W - 0.5) * tau(X) + 0.5 * rnorm(n)
  Y_train <- Y[idx.train: n]
  Y_test <- Y[1: idx.test]
  
  dta <- list("y.tr" = Y_train, "y.te" = Y_test, 
              "w.tr" = W_train, "w.te" = W_test,
              "x.tr" = X_train, "x.te" = X_test)
  return(dta)
}
dta <- gen_data(n = n, p = p)

# oracle 
true_mu1 <- function(y, a, x, new.x) {
  b(new.x) + 0.5 * tau(new.x)
}

true_mu0 <- function(y, a, x, new.x){
  b(new.x) - 0.5 * tau(new.x)
}

true_mu <- function(y, a, x, new.x){
  b(new.x) + (e(new.x)-0.5) * tau(new.x)
}



get_MSEs <- function(dta, learner){
  
  # calculate MSE for a given learner and dataset
  # return a num
  
  # dta: list, combination of all X, Y, W data
  # learner: list of learner
  
  true_est <- true_mu1(new.x = dta$x.te) - true_mu0(new.x = dta$x.te)
  
  if (learner == 'R'){
    fit <- rlasso(dta$x.tr, dta$w.tr, dta$y.tr)
    tau_cv <- predict(fit, dta$x.te)
  }
  else if (learner == 'T'){
    fit = tlasso(dta$x.tr, dta$w.tr, dta$y.tr)
    tau_cv <- predict(fit, dta$x.te)
  }
  else if (learner == 'X'){
    fit = xlasso(dta$x.tr, dta$w.tr, dta$y.tr)
    tau_cv <- predict(fit, dta$x.te)
  }
  else if (learner == 'DR'){
    tau_cv <- cate(v0 = dta$x.te, v = dta$x.tr, learner = "dr", y = dta$y.tr,
                   a = dta$w.tr, x = dta$x.tr, drl = drl.x,
                   mu1.x = mu1.x, mu0.x = mu0.x, pi.x = pi.x, nsplits = 2)
  }
  else if (learner == 'R-Oracle'){
    fit <- rlasso(dta$x.tr, dta$w.tr, dta$y.tr, p_hat = e(dta$x.tr), 
                  m_hat = true_mu(new.x = dta$x.tr))
    tau_cv <- predict(fit, dta$x.te)
  }
  else if (learner == 'DR-Oracle'){
    tau_cv <- cate(v0 = dta$x.te, v = dta$x.tr, learner = "dr", y = dta$y.tr,
                   a = dta$w.tr, x = dta$x.tr, drl = drl.x,
                   mu1.x = true_mu1(new.x = dta$x.tr), 
                   mu0.x = true_mu0(new.x = dta$x.tr), 
                   pi.x = e(dta$x.tr), nsplits = 2)
  }
  else {
    stop("wrong learner.")
  }
  # calculate MSE iteration
  if(learner %in% c("R", "X", "T", "R-Oracle")){
    mse <- mean((tau_cv - true_est)^2)
  }
  else if(learner %in% c("DR", "DR-Oracle")){
    mse <- mean((tau_cv$est[[1]][, 1] - true_est)^2)
  }
  else{
    stop("wrong learner.")
  }
  return(mse)
}



MSE_matrix <- function(n, p, t, learner_list){
  
  # return a matrix of t*length(learner_list)
  # each row is a simulation using the same dataset in different learners
  
  # n: sample size
  # p: number of covariates
  # t: time of simulation
  # learner_list: a list of learner that to be simulated
  
  MSEs <- matrix(nrow = t, ncol = length(learner_list))
  rownames(MSEs) <- c(1:t)
  colnames(MSEs) <- learner_list
  
  for (i in 1:t){
    dta <- gen_data(n = n, p = p)
    
    for (j in 1:length(learner_list)){
      MSEs[i, j] <- get_MSEs(dta = dta, learner = learner_list[j]) 
    }
    print(paste("Iteration", i, "finished!"))
    print(paste("MSEs are", MSEs[i, ]))
  }
  return(MSEs)
}




# simulation 1, compare MSEs across different learners


learner_list <- list("T", "R", "X", "DR", "DR-Oracle", "R-Oracle")
# learner_list <- list("R","DR", "DR-Oracle", "R-Oracle")
MSEs <- MSE_matrix(n = 400, p = 10, t = 20, learner_list = learner_list)

# write.csv(MSEs,"C:\\Users\\ymiao\\Desktop\\test.csv",row.names=FALSE)

# draw a boxplot
MSEs_melt <- melt(MSEs, id = 'num')

ggplot(MSEs_melt, aes(Var2, value)) + geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim=c(0, 0.20)) +
  labs(title="Performance of different learners", x="Meta-Learners", y = "MSEs")


ggplot(MSEs_melt, aes(Var2, value)) + geom_boxplot() +
  labs(title="Performance of different learners", x="Meta-Learners", y = "MSEs")
  







# simulation 2, different n, different learner, train 100 times:
n_lst <- c(200*2, 400*2, 800*2, 1600*2, 3200*2)
learner_list <- list("T", "R", "X", "DR", "DR-Oracle", "R-Oracle")
MSEs_simu2 <- matrix(nrow = length(n_lst), ncol = length(learner_list))
colnames(MSEs_simu2) <- learner_list
rownames(MSEs_simu2) <- n_lst

t = 100

for (n in 1:length(n_lst)){
  start_time <- Sys.time()
  
  print(paste("Sample Size is", n_lst[n], ":"))
  
  temp_MSEs <- MSE_matrix(n = n_lst[n], t = t, learner_list = learner_list)
  mean_MSES <- colMeans(temp_MSEs)
  MSEs_simu2[n, ] <- mean_MSES
  
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  print(time_taken)
}


# grouped bar chart
MSE_sub2 <- melt(output, id = 'num')
MSE_sub2$Var1 <- as.factor(MSE_sub2$Var1)
pic <- ggplot(MSE_sub2, aes(x = Var2, y = value, fill = Var1)) +   
  geom_bar(stat="identity", position = "dodge2", linetype = 0 ) +
  coord_cartesian(ylim=c(0, 0.5)) +
  labs(title="Performance of different learners under different sample size n", x="Meta-Learners", y = "MSEs") + 
  scale_y_continuous(expand = c(0, 0)) 

pic



# simulation3: different p, different learner, train 100 times:
# p_lst <- c(10, 100, 500)
p_lst <- c(10, 100)
learner_list <- list("T", "R", "X", "DR", "DR-Oracle", "R-Oracle")
MSEs_simu3 <- matrix(nrow = length(p_lst), ncol = length(learner_list))
colnames(MSEs_simu3) <- learner_list
rownames(MSEs_simu3) <- p_lst

t = 20

for (n in 1:length(p_lst)){
  start_time <- Sys.time()
  
  print(paste("Covariate Space Dimension is", p_lst[n], ":"))
  
  temp_MSEs <- MSE_matrix(n = 400, p = p_lst[n], t = t, learner_list = learner_list)
  mean_MSES <- colMeans(temp_MSEs)
  MSEs_simu3[n, ] <- mean_MSES
  print(paste("mean MSES when p = ", p_lst[n], "is", mean_MSES))
  
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  print(time_taken)
}


# when p = 500, it takes a long time to run (~15min each loop)
tst3 <- MSE_matrix(n = 400, p = 500, t = 20, learner_list = c('R', 'DR', 'R-Oracle', 'DR-Oracle'))
tst3


remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}



# have a close look at n = 400, p = 10 situation
test_dr <- function(n, p, t){
  x.tr <- array(data = NA, dim = c(n/2, p, t))
  x.te <- array(data = NA, dim = c(n/2, p, t))
  y.tr <- matrix(nrow = n/2, ncol = t)
  y.te <- matrix(nrow = n/2, ncol = t)
  w.tr <- matrix(nrow = n/2, ncol = t)
  w.te <- matrix(nrow = n/2, ncol = t)
  true_dr <- matrix(nrow = n/2, ncol = t)
  est <- matrix(nrow = n/2, ncol = t)
  mse <- matrix(nrow = 1, ncol = t)
  
  for (i in 1:t){
    init_dta <- gen_data(n, p)
    true_est_dr <- true_mu1(new.x = init_dta$x.te) - true_mu0(new.x = init_dta$x.te)
    tst_drl <- cate(v0 = init_dta$x.te, v = init_dta$x.tr, learner = "dr", y = init_dta$y.tr,
                    a = init_dta$w.tr, x = init_dta$x.tr, drl = drl.x,
                    mu1.x = mu1.x, mu0.x = mu0.x, pi.x = pi.x, nsplits = 2)
    est_drl <- tst_drl$est[[1]][, 1]
    mse_drl <- mean((tst_drl$est[[1]][, 1] - true_est_dr)^2)
    print(mse_drl)
    
    x.tr[,,i] <- init_dta$x.tr
    x.te[,,i] <- init_dta$x.te
    y.tr[,i] <- init_dta$y.tr
    y.te[,i] <- init_dta$y.te
    w.tr[,i] <- init_dta$w.tr
    w.te[,i] <- init_dta$w.te
    true_dr[,i] <- true_est_dr
    est[,i] <- est_drl
    mse[,i] <- mse_drl
  }
  res = list(x.tr = x.tr, x.te = x.te, y.tr = y.tr, y.te = y.te, w.tr = w.tr, w.te = w.te,
             true_dr = true_dr, est = est, mse = mse)
  
  return(res)
}


# To capture the original x,y,w that causes large mse
test.1 <- test_dr(n = 400, p = 10, t = 20)

x.tr <- test.1[["x.tr"]][,,5]
x.te <- test.1[["x.te"]][,,5]
y.tr <- test.1[["y.tr"]][,5]
y.te <- test.1[["y.te"]][,5]
w.tr <- test.1[["w.tr"]][,5]
w.te <- test.1[["w.te"]][,5]


# run this iteration, then notice that outliers come from randomness of cate function
for (i in 1:50){
  true_est_dr <- test.1[["true_dr"]][,5]
  tst.res.7cate <- cate(v0 = x.te, v = x.tr, learner = "dr", y = y.tr,
                  a = w.tr, x = x.tr, drl = drl.x,
                  mu1.x = mu1.x, mu0.x = mu0.x, pi.x = pi.x, nsplits = 2)
  mse.7cate <- mean((tst.res.7cate$est[[1]][, 1] - true_est_dr)^2)
  print(mse.7cate)
  i = i + 1
}


# if manully run first and second stage regression: also get strange result occasionally
for (i in 1:20){
  mu1hat <- mu1.x(y = y.tr, a.tr = w.tr, x = x.tr, new.x = x.te)
  mu0hat <- mu0.x(y = y.tr, a.tr = w.tr, x = x.tr, new.x = x.te)
  pihat <- pi.x(a = w.tr, x = x.tr, new.x = x.te)
  pseudo <- (w.te - pihat) / (pihat * (1 - pihat)) *
    (y.te - w.te * mu1hat - (1 - w.te) * mu0hat) + mu1hat - mu0hat
  
  
  tst.res.7 <- drl.x(y.tr = pseudo, x.tr = x.te, new.x = x.te)
  res <- mean((tst.res.7 - true_est_dr)^2)
  
  print(res)
  i = i+1
}



# to do: revise nuisance functions with glmnet, record results in both first and second regression







# grouped bar chart
MSE_sub3 <- melt(MSEs_simu3, id = 'num')
MSE_sub3$Var1 <- as.factor(MSE_sub3$Var1)
pic <- ggplot(MSE_sub3, aes(x = Var2, y = value, fill = Var1)) +   
  geom_bar(stat="identity", position = "dodge2", linetype = 0 ) +
  coord_cartesian(ylim=c(0, 0.5)) +
  labs(title="Performance of different learners under different p", x="Meta-Learners", y = "MSEs") + 
  scale_y_continuous(expand = c(0, 0)) 

pic





