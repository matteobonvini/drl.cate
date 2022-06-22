

# some simulations to compare the performance of r-learner, dr-learner, and other learners.
rm(list = ls())
#import packages
library(rlearner)
library(drl.cate)
library(glmnet)

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
# settings are based on Okasa(2022)
n <- 2000; p <- 10;

tau <- function(X){
  1
}

true_mu0 <- function(X){
  sin(pi * X[,1] * X[,2]) + 2 * (X[,3] - 0.5)^2 + X[,4] + 0.5 * X[,5]
}

true_mu1 <- function(X){
  true_mu0(X) + tau(X)
}

true_mu <- function(X, W){
  true_mu0(X) + W 
}

e <- function(X){
  1/12 * (1 + dbeta(sin(pi * X[,1] * X[,2] * X[,3] * X[,4] * X[,5]), 2, 4))
}

gen_data <- function(n, p){
  idx.test <- n/2; idx.train <- idx.test + 1
  
  X <- matrix(runif(n*p, min=0, max=1), n, p)
  X_train <- X[idx.train:n,]
  X_test <- X[1:idx.test,]
  
  W <- rbinom(n, 1, e(X))
  W_train <- W[idx.train: n]
  W_test <- W[1: idx.test]
  
  Y <- W * (true_mu1(X) + rnorm(n)) + (1 - W) * (true_mu0(X) + rnorm(n))
  Y_train <- Y[idx.train: n]
  Y_test <- Y[1: idx.test]
  
  dta <- list("y.tr" = Y_train, "y.te" = Y_test, 
              "w.tr" = W_train, "w.te" = W_test,
              "x.tr" = X_train, "x.te" = X_test)
  return(dta)
}
dta <- gen_data(n = n, p = p)


get_MSEs <- function(dta, learner){
  
  # calculate MSE for a given learner and dataset
  # return a num
  
  # dta: list, combination of all X, Y, W data
  # learner: R, T, X, DR, R-Oracle or DR-Oracle
  
  true_est <- tau(dta$x.te)
  
  if (learner == 'R'){
    fit <- rlasso(dta$x.tr, dta$w.tr, dta$y.tr, lambda_choice = "lambda.min")
    tau_cv <- predict(fit, dta$x.te)
  }
  else if (learner == 'T'){
    fit = tlasso(dta$x.tr, dta$w.tr, dta$y.tr, lambda_choice = "lambda.min")
    tau_cv <- predict(fit, dta$x.te)
  }
  else if (learner == 'X'){
    fit = xlasso(dta$x.tr, dta$w.tr, dta$y.tr, lambda_choice = "lambda.min")
    tau_cv <- predict(fit, dta$x.te)
  }
  else if (learner == 'DR'){
    # tau_cv <- cate(v0 = dta$x.te, v = dta$x.tr, learner = "dr", y = dta$y.tr,
    #                a = dta$w.tr, x = dta$x.tr, drl = drl.x,
    #                mu1.x = mu1.x, mu0.x = mu0.x, pi.x = pi.x, nsplits = 2)
    tau_cv <- drlasso(v0 = dta$x.te, y = dta$y.tr, a = dta$w.tr, x = dta$x.tr,
                      v = dta$x.tr, nsplits = 5, nfolds = 3)
  }
  else if (learner == 'R-Oracle'){
    fit <- rlasso(dta$x.tr, dta$w.tr, dta$y.tr, lambda_choice = "lambda.min",
                  p_hat = e(dta$x.tr), 
                  m_hat = true_mu(X = dta$x.tr, W = dta$w.tr))
    tau_cv <- predict(fit, dta$x.te)
  }
  else if (learner == 'DR-Oracle'){
    # tau_cv <- cate(v0 = dta$x.te, v = dta$x.tr, learner = "dr", y = dta$y.tr,
    #                a = dta$w.tr, x = dta$x.tr, drl = drl.x,
    #                mu1.x = true_mu1(dta$x.tr),
    #                mu0.x = true_mu0(dta$x.tr),
    #                pi.x = e(dta$x.tr), nsplits = 2)
    # to do
    tau_cv <- drlasso(v0 = dta$x.te, y = dta$y.tr, a = dta$w.tr, x = dta$x.tr,
                      v = dta$x.tr, nsplits = 5, nfolds = 3, oracle =TRUE)
  }
  else {
    stop("wrong learner.")
  }
  # calculate MSE iteration
  if(learner %in% c("R", "X", "T", "R-Oracle")){
    mse <- mean((tau_cv - true_est)^2)
  }
  else if(learner %in% c("DR", "DR-Oracle")){
    dr.est <- rowMeans(tau_cv[["pred"]])
    # dr.est <- tau_cv[["est"]][[1]][,1]
    mse <- mean((dr.est - true_est)^2)
    print('dr')
    # print(tau_cv)
  }
  else{
    stop("wrong learner.")
  }
  return(mse)
}

get_MSEs(dta, "DR")

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
    # print(paste("MSEs are", MSEs[i, ]))
  }
  return(MSEs)
}




# simulation 1, compare MSEs across different learners


# learner_list <- list("T", "R", "X", "DR", "R-Oracle", "DR-Oracle")
learner_list <- list("T", "R", "X", "DR", "DR-Oracle")
MSEs_simu1 <- MSE_matrix(n = 2000, p = 10, t = 40, learner_list = learner_list)

# write.csv(MSEs,"C:\\Users\\ymiao\\Desktop\\test.csv",row.names=FALSE)

# draw a boxplot
MSEs100 <- MSEs_simu1 * 100
MSEs_melt <- melt(MSEs100, id = 'num')

ggplot(MSEs_melt, aes(Var2, value)) + geom_boxplot(outlier.shape = NA) +
  # coord_cartesian(ylim=c(0, 0.70)) +
  labs(title="Performance of different learners", x="Meta-Learners", y = "MSEs")


ggplot(MSEs_melt, aes(Var2, value)) + geom_boxplot() +
  labs(title="Performance of different learners", x="Meta-Learners", y = "MSEs")
  





# simulation 2, different n, different learner, train 100 times:
n_lst <- c(200*2, 400*2, 800*2, 1600*2, 3200*2)
# learner_list <- list("T", "R", "X", "DR", "DR-Oracle", "R-Oracle")
MSEs_simu2 <- matrix(nrow = length(n_lst), ncol = length(learner_list))
colnames(MSEs_simu2) <- learner_list
rownames(MSEs_simu2) <- n_lst

t = 40
p = 10

for (n in 1:length(n_lst)){
  start_time <- Sys.time()
  
  print(paste("Sample Size is", n_lst[n], ":"))
  
  temp_MSEs <- MSE_matrix(n = n_lst[n], p = p, t = t, learner_list = learner_list)
  mean_MSES <- colMeans(temp_MSEs)
  MSEs_simu2[n, ] <- mean_MSES
  
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  print(time_taken)
}


# grouped bar chart
MSE_sub2 <- melt(MSEs_simu2 * 100, id = 'num')

MSE_sub2$Var1 <- as.factor(MSE_sub2$Var1)
pic <- ggplot(MSE_sub2, aes(x = Var2, y = value, fill = Var1)) +   
  geom_bar(stat="identity", position = "dodge2", linetype = 0 ) +
  coord_cartesian(ylim=c(0, 30)) +
  labs(title="Performance of different learners under different sample size n", x="Meta-Learners", y = "MSEs*100") + 
  scale_y_continuous(expand = c(0, 0)) 

pic2 <- ggplot(MSE_sub2, aes(x = Var1, y = value, color = Var2)) +   
  geom_line(aes(group = Var2), size = 1) +
  coord_cartesian(ylim=c(0, 30)) +
  theme(legend.position = c(0.8, 0.8)) + 
  labs(title="Performance of different learners under different sample size n", x="Meta-Learners", y = "MSEs*100") + 
  scale_y_continuous(expand = c(0, 0)) 

pic2



# simulation3: different p, different learner, train 100 times:
# p_lst <- c(10, 100, 500)
p_lst <- c(10, 40, 100)
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
# tst3 <- MSE_matrix(n = 400, p = 500, t = 20, learner_list = c('R', 'DR', 'R-Oracle', 'DR-Oracle'))
# tst3
# 
# 
# remove_outliers <- function(x, na.rm = TRUE, ...) {
#   qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
#   H <- 1.5 * IQR(x, na.rm = na.rm)
#   y <- x
#   y[x < (qnt[1] - H)] <- NA
#   y[x > (qnt[2] + H)] <- NA
#   y
# }
# 


# grouped bar chart
MSE_sub3 <- melt(MSEs_simu3, id = 'num')
MSE_sub3$Var1 <- as.factor(MSE_sub3$Var1)
pic <- ggplot(MSE_sub3, aes(x = Var2, y = value, fill = Var1)) +   
  geom_bar(stat="identity", position = "dodge2", linetype = 0 ) +
  coord_cartesian(ylim=c(0, 0.5)) +
  labs(title="Performance of different learners under different p", x="Meta-Learners", y = "MSEs") + 
  scale_y_continuous(expand = c(0, 0)) 

pic





