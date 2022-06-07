
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
mu1.x <- function(y.tr, a.tr, x.tr, new.x, sl.lib = c("SL.mean", "SL.lm", "SL.gam",
                                                      "SL.polymars", "SL.rpart")){
  y1.tr <- y.tr[a.tr == 1]
  x1.tr <- as.data.frame(x.tr[a.tr == 1, , drop = FALSE])
  fit <- SuperLearner::SuperLearner(Y = y1.tr, X = x1.tr, SL.library = sl.lib,
                                    newX = as.data.frame(new.x))
  return(fit$SL.predict)
}

mu0.x <- function(y.tr, a.tr, x.tr, new.x, sl.lib = c("SL.mean", "SL.lm", "SL.gam",
                                                      "SL.polymars", "SL.rpart")){
  y0.tr <- y.tr[a.tr == 0]
  x0.tr <- as.data.frame(x.tr[a.tr == 0, , drop = FALSE])
  fit <- SuperLearner::SuperLearner(Y = y0.tr, X = x0.tr, SL.library = sl.lib,
                                    newX = as.data.frame(new.x))
  return(fit$SL.predict)
}

mu.x <- function(y.tr, a.tr, x.tr, new.x, sl.lib = c("SL.mean", "SL.lm", "SL.gam",
                                                     "SL.polymars", "SL.rpart")){
  fit <- SuperLearner::SuperLearner(Y = y.tr, X = as.data.frame(x.tr), SL.library = sl.lib,
                                    newX = as.data.frame(new.x))
  return(fit$SL.predict)
}

pi.x <- function(a.tr, x.tr, new.x, sl.lib = c("SL.mean", 
                                               "SL.ranger", "SL.glmnet")){
  fit <- SuperLearner::SuperLearner(Y = a.tr, X = as.data.frame(x.tr), SL.library = sl.lib,
                                    newX = as.data.frame(new.x))
  return(fit$SL.predict)
}



# drl functions:
drl.x <- function(y.tr, x.tr, new.x, sl.lib = c("SL.mean", "SL.lm", "SL.gam",
                                                "SL.loess", "SL.glm", "SL.polymars",
                                                "SL.rpart")){
  fit <- SuperLearner::SuperLearner(Y = y.tr, X = as.data.frame(x.tr),
                                    SL.library = sl.lib,
                                    newX = as.data.frame(new.x))
  return(fit$SL.predict)
}

# individual treatment effect, same function as drl.x, for second stage regression
# note that y.tr could be pseudo \phi
drl.ite <- function(y.tr, x.tr, new.x, 
                    sl.lib = c("SL.mean", "SL.lm", "SL.gam",
                               "SL.loess", "SL.glm", "SL.polymars",
                               "SL.rpart")){
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



# run a learner return MSE
get_MSEs <- function(dta, learner){
  # data: list, combination of all X, Y, W data
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

# simulate several times and get a list of MSEs from each run
mse_list <- function(n, p, learner, t){
  # t: number of iterations 
  start_time <- Sys.time()
  mse <- numeric(length = length(t))
  for(i in 1:t){
    dta <- gen_data(n = n, p = p)
    mse[i] <- get_MSEs(dta = dta, learner = learner)
    print(paste(learner, "Learner,", "Iteration", i, ":MSE is", mse[i]))
    # print("MSE of iteration", i, "is", mse[i])
  }
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  print(time_taken)
  return(mse)
}

# first simulation, compare MSEs across different learners
# simulation is much faster with lasso than boost. The whole process takes ~50min
mse_t <- mse_list(n = 2000, p = 10, learner = 'T', t = 100)
mse_r <- mse_list(n = 2000, p = 10, learner = 'R', t = 100)
mse_x <- mse_list(n = 2000, p = 10, learner = 'X', t = 100)
mse_dr <- mse_list(n = 2000, p = 10, learner = 'DR', t = 100)
mse_dr_oracle <- mse_list(n = 2000, p = 10, learner = 'DR-Oracle', t = 100)
mse_r_oracle <- mse_list(n = 2000, p = 10, learner = 'R-Oracle', t = 100)

# box plot with and without outliers
x.ax <- c(1:100)
MSEs <- data.frame(
  "num" = x.ax,
  "DR-Learner" = mse_dr,
  "DR-Oracle" = mse_dr_oracle,
  "R-learner" = mse_r,
  "R-Oracle" = mse_r_oracle,
  "T-learner" = mse_t,
  "X-learner" = mse_x
)
write.csv(MSEs,"C:\\Users\\ymiao\\Desktop\\test.csv",row.names=FALSE)

MSEs_melt <- melt(MSEs, id = 'num')

ggplot(MSEs_melt, aes(variable, value)) + geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim=c(0, 0.10)) +
  labs(title="Performance of different learners", x="Meta-Learners", y = "MSEs")


ggplot(MSEs_melt, aes(variable, value)) + geom_boxplot() +
  labs(title="Performance of different learners", x="Meta-Learners", y = "MSEs")
  


# second simulation, with different sample size n:
n_lst <- c(200*2, 400*2, 800*2, 1600*2, 3200*2)
learner_list <- list("T", "R", "X", "DR", "DR-Oracle", "R-Oracle")
output <- matrix(nrow = length(n_lst), ncol = length(learner_list))
rownames(output) <- n_lst
colnames(output) <- learner_list

t = 100

# when t = 100, takes about 4~6h to run
for (i in 1:length(n_lst)){
  for (j in 1:length(learner_list)){
    print(n_lst[i])
    output[i, j] <- mean(mse_list(n = n_lst[i], p = p, learner = learner_list[j], t = t))
  }
}

# grouped bar chart
MSE_sub2 <- melt(output, id = 'num')
MSE_sub2$Var1 <- as.factor(MSE_sub2$Var1)
pic <- ggplot(MSE_sub2, aes(x = Var2, y = value, fill = Var1)) +   
  geom_bar(stat="identity", position = "dodge2", linetype = 0 ) +
  coord_cartesian(ylim=c(0, 0.15)) +
  labs(title="Performance of different learners under different sample size n", x="Meta-Learners", y = "MSEs") + 
  scale_y_continuous(expand = c(0, 0)) 

pic




# finally, need to plot true tau and x, to make sure that they are truly linear

