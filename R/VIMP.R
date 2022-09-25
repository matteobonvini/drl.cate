
rm(list = ls())
library(speff2trial) # dataset for part 4
library(dplyr)
library(SuperLearner)
library(ggplot2)
library(mgcv)

################### Section 3 : Simulation Study #####################
################### no sample splitting ##############################

nseq <- c(500, 1000, 2000, 3000, 4000)
expit <- function(x){ exp(x)/(1+exp(x))}

tau <- function(x){
  x1 <- x[,1]
  x2 <- x[,2]
  res <- x1*x1*(x1 + 7/5) + 25*x2*x2/9
  return (res)
}

t <- 100
res <- data.frame(matrix(nrow = length(nseq), ncol = 8))
rownames(res) <- nseq
colnames(res) <- c('psi_1_hat', 'psi_1_hat_bias', 'variance_1', 'coverage_1', 'psi_2_hat', 'psi_2_hat_bias', 'variance_2', 'coverage_2')



for (j in 1:length(nseq)){
  
  n <- nseq[j]
  tmp_res <- data.frame(matrix(nrow = t, ncol = 8))
  colnames(tmp_res) <- c('psi_1_hat', 'psi_1_hat_bias', 'variance_1', 'coverage_1', 'psi_2_hat', 'psi_2_hat_bias', 'variance_2', 'coverage_2')
  
  
  for (i in 1:t){
    x <- as.data.frame(replicate(2, runif(n, -1, 1)))
    a <- rbinom(n, 1, expit(-0.4*x[,1] + 0.1*x[,1]*x[,2]))
    y <- rnorm(n, x[,1]*x[,2] + 2*x[,2]^2 - x[,1] + a*tau(x), 1)
    
    # true estimand values
    tau_p <- mean(tau(x))
    theta_p <- var(tau(x))
    
    tau_1 <- 25/9*x[,2]^2 + 7/15
    tau_2 <- x[,1]^2*(x[,1] + 7/5) + 25/27
    
    psi_1 <- 1 - var(tau_1)/var(tau(x))
    psi_2 <- 1 - var(tau_2)/var(tau(x))
    
    # estimated values
    data <- cbind(y, a, x)
    colnames(data) <- c('y', 'a', 'x1', 'x2')
    
    data_0 <- data[a == 0, , drop = FALSE]
    data_1 <- data[a == 1, , drop = FALSE]
    
    mu0hat <- predict(gam(y ~ s(x1) + s(x2) + te(x1, x2), data = data_0), newdata = data)
    mu1hat <- predict(gam(y ~ s(x1) + s(x2) + te(x1, x2), data = data_1), newdata = data)
    # pihat <- predict(gam(a ~ s(x1) + s(x2), family = 'binomial', data = data), newdata = data)
    pihat <- expit(predict(gam(a ~ x1 + x2 + x1*x2, family = 'binomial', data = data), newdata = data))
    
    pseudo_hat <- ((a-pihat)/(pihat*(1-pihat)))*(y-a*mu1hat-(1-a)*mu0hat) + mu1hat-mu0hat
    
    # T-learner, algorithm A
    # tau_hat <- mu1hat - mu0hat
    
    # DR-learner, algorithm B
    tau_hat <- predict(gam(pseudo_hat ~ s(x1) + s(x2), data = cbind(data, pseudo_hat)), newdata = cbind(data, pseudo_hat))
    
    
    # try mgcv
    tau_1_hat <- predict(gam(tau_hat ~ s(x2), data = cbind(data, tau_hat)), newdata = cbind(data, tau_hat))
    tau_2_hat <- predict(gam(tau_hat ~ s(x1), data = cbind(data, tau_hat)), newdata = cbind(data, tau_hat))
    
    
    tau_1_hat <- tau_1
    tau_2_hat <- tau_2
    tau_hat <- tau(x)
    pihat <- expit(-0.4*x[,1] + 0.1*x[,1]*x[,2])
    mu1hat <- x[,1]*x[,2] + 2*x[,2]^2 - x[,1] + tau(x)
    mu0hat <- x[,1]*x[,2] + 2*x[,2]^2 - x[,1]
    pseudo_hat <- ((a-pihat)/(pihat*(1-pihat)))*(y-a*mu1hat-(1-a)*mu0hat) + mu1hat-mu0hat
    
    tau_p_hat <- mean(tau_hat)
    
    theta_1_hat <- mean((pseudo_hat - tau_1_hat)^2 - (pseudo_hat - tau_hat)^2)
    theta_2_hat <- mean((pseudo_hat - tau_2_hat)^2 - (pseudo_hat - tau_hat)^2)
    theta_p_hat <- mean((pseudo_hat - tau_p_hat)^2 - (pseudo_hat - tau_hat)^2)
    
    psi_1_hat <- theta_1_hat / theta_p_hat
    psi_2_hat <- theta_2_hat / theta_p_hat
    
    # bias
    tmp_res[i, 1] <- psi_1_hat
    tmp_res[i, 2] <- (psi_1_hat - psi_1) * sqrt(n)
    tmp_res[i, 5] <- psi_2_hat
    tmp_res[i, 6] <- (psi_2_hat - psi_2) * sqrt(n)
    
    
    # variance
    phi_1_hat <- ((pseudo_hat - tau_1_hat)^2 - psi_1_hat*(pseudo_hat - tau_p_hat)^2 + (psi_1_hat - 1)*(pseudo_hat - tau_hat)^2) / theta_p_hat
    variance_1 <- 1/n^2 * sum((phi_1_hat)^2)
    
    phi_2_hat <- ((pseudo_hat - tau_2_hat)^2 - psi_2_hat*(pseudo_hat - tau_p_hat)^2 + (psi_2_hat - 1)*(pseudo_hat - tau_hat)^2) / theta_p_hat
    variance_2 <- 1/n^2 * sum((phi_2_hat)^2)
    
    tmp_res[i, 3] <- variance_1 * n
    tmp_res[i, 7] <- variance_2 * n
    
    # coverage
    if (psi_1 <= psi_1_hat + 1.96 * sqrt(variance_1) & psi_1 >= psi_1_hat - 1.96 * sqrt(variance_1)){
      coverage_1 <- 1
    } else {
      coverage_1 <- 0
    }
    
    if (psi_2 <= psi_2_hat + 1.96 * sqrt(variance_2) & psi_2 >= psi_2_hat - 1.96 * sqrt(variance_2)){
      coverage_2 <- 1
    } else {
      coverage_2 <- 0
    }
    
    tmp_res[i, 4] <- coverage_1
    tmp_res[i, 8] <- coverage_2
    
    print(paste('Iteration', i, 'finished!'))
  }
  print(tmp_res)
  res[j, ] <- colMeans(tmp_res)
}

res







################### Section 4: Real World Example #####################
dta_aids <- get(data(ACTG175))
dta2 <- filter(dta_aids, arms %in% c(1,3))
# add treatment
dta2 <- mutate(dta2, treatment = ifelse(arms == 3, 0, 1))
# dataset used in section 4:
dta2 <- dta2[, c('cd420', 'treatment', 
                 'age', 'wtkg', 'karnof', 'cd40', 'cd80', 'gender', 'homo', 'race', 'symptom', 'drugs', 'hemo', 'str2')]
y <- dta2[, 'cd420']
a <- dta2[, 'treatment']
x <- dta2[, 3:14]


sl.lib <- c("SL.glm", "SL.gam", "SL.glmnet", "SL.xgboost", "SL.ranger")
# sl.lib <- c("SL.glm", "SL.gam", "SL.glmnet")


# 1. without sample splitting. 

cate_estimate <- function(x, y, a, x.estimate, model){
  # function
  # x, y, a: standard input
  # model: A for t-learner, B for df-learner
  
  mu0hat <- SuperLearner(Y = y[a == 0], X = data.frame(x[a == 0, , drop = FALSE]), SL.library = sl.lib,
                         newX = data.frame(x))$SL.predict
  mu1hat <- SuperLearner(Y = y[a == 1], X = data.frame(x[a == 1, , drop = FALSE]), SL.library = sl.lib,
                         newX = data.frame(x))$SL.predict
  pihat <- mean(a)
  # pihat <- SuperLearner(Y = a, X = data.frame(x), SL.library = sl.lib)$SL.predict
  
  # pihat <- mean(a)
  pseudo <- ((a-pihat)/(pihat*(1-pihat)))*(y-a*mu1hat-(1-a)*mu0hat) + mu1hat-mu0hat
  
  if (model == 'A'){
    tau <- mu1hat - mu0hat
  }
  else if (model == 'B'){
    tau <- SuperLearner(Y = pseudo, X = as.data.frame(x), SL.library = sl.lib,
                        newX = as.data.frame(x))$SL.predict
  }
  else {
    stop('invalid type!')
  }
  out <- list(tau = tau, pseudo = pseudo)
  return (out)
  
}
## T learner
VIM_2a <- data.frame(matrix(nrow = length(x), ncol = 3))
rownames(VIM_2a) <- colnames(x)
colnames(VIM_2a) <- c('T', 'Lower_Bound', 'Upper_Bound')

cate <- cate_estimate(x, y, a, x, 'A')
tau_hat <- cate$tau
pseudo_hat <- cate$pseudo

tau_p_hat <- mean(tau_hat)

for (i in 1:length(x)){
  tau_s_hat <- SuperLearner(Y = tau_hat, X = x[,-i], 
                            SL.library = sl.lib,
                            newX = x[,-i])$SL.predict
  theta_s_hat <- mean((pseudo_hat - tau_s_hat)^2 - (pseudo_hat - tau_hat)^2)
  theta_p_hat <- mean((pseudo_hat - tau_p_hat)^2 - (pseudo_hat - tau_hat)^2)
  # theta_p_hat <- as.numeric(var(tau_hat))
  
  psi_hat <- theta_s_hat / theta_p_hat
  
  phi_hat <- ((pseudo_hat - tau_s_hat)^2 - psi_hat * (pseudo_hat - tau_p_hat)^2 + (psi_hat - 1)*(pseudo_hat - tau_hat)^2) / theta_p_hat
  var <- 1/length(phi_hat)^2 * sum((phi_hat)^2)
  
  VIM_2a[i, 1] <- psi_hat
  VIM_2a[i, 2] <- VIM_2a[i, 1] - 1.96 * sqrt(var)
  VIM_2a[i, 3] <- VIM_2a[i, 1] + 1.96 * sqrt(var)
  print(c(VIM_2a[i,1], VIM_2a[i,2], VIM_2a[i,3], theta_s_hat, theta_p_hat))
}

VIM_2a <- VIM_2a[order(VIM_2a$T, decreasing = FALSE),]
### draw picture
VIM_2a <- cbind(as.matrix(row.names(VIM_2a), nrow(VIM_2a), 1), VIM_2a)
colnames(VIM_2a) <- c('variable', 'T', 'Lower_Bound', 'Upper_Bound')
VIM_2a <- cbind(as.matrix(row.names(VIM_2a), nrow(VIM_2a), 1), VIM_2a)
VIM_2a <- VIM_2a[order(VIM_2a$T, decreasing = FALSE),]
order_2a <- VIM_2a$variable


fig_2a <- ggplot(VIM_2a, aes(x = factor(variable, level = order_2a), y = T))+
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower_Bound, ymax = Upper_Bound), width=.1, position=position_dodge(0.2))+
  theme(axis.title.y=element_blank(), 
        axis.text.y = element_text(color = "grey20", size = 10),
        axis.text.x = element_text(color = "grey20", size = 10)) +
  coord_flip() +
  ylim(0, 1)
fig_2a



## DR learner
VIM_2b <- data.frame(matrix(nrow = length(x), ncol = 3))
rownames(VIM_2b) <- colnames(x)
colnames(VIM_2b) <- c('DR', 'Lower_Bound', 'Upper_Bound')

cate <- cate_estimate(x, y, a, x, 'B')
tau_hat <- cate$tau
pseudo_hat <- cate$pseudo

tau_p_hat <- mean(tau_hat)

for (i in 1:length(x)){
  tau_s_hat <- SuperLearner(Y = tau_hat, X = x[,-i], 
                            SL.library = sl.lib,
                            newX = x[,-i])$SL.predict
  theta_s_hat <- mean((pseudo_hat - tau_s_hat)^2 - (pseudo_hat - tau_hat)^2)
  theta_p_hat <- mean((pseudo_hat - tau_p_hat)^2 - (pseudo_hat - tau_hat)^2)
  # theta_p_hat <- as.numeric(var(tau_hat))
  
  psi_hat <- theta_s_hat / theta_p_hat
  
  phi_hat <- ((pseudo_hat - tau_s_hat)^2 - psi_hat * (pseudo_hat - tau_p_hat)^2 + (psi_hat - 1)*(pseudo_hat - tau_hat)^2) / theta_p_hat
  var <- 1/length(phi_hat)^2 * sum((phi_hat)^2)
  
  VIM_2b[i, 1] <- psi_hat
  VIM_2b[i, 2] <- VIM_2b[i, 1] - 1.96 * sqrt(var)
  VIM_2b[i, 3] <- VIM_2b[i, 1] + 1.96 * sqrt(var)
  print(c(VIM_2b[i,1], VIM_2b[i,2], VIM_2b[i,3], theta_s_hat, theta_p_hat))
}
VIM_2b <- VIM_2b[order(VIM_2b$DR, decreasing = FALSE),]
# VIM_2b <- VIM_2b[order(VIM_2b$DR),]
### draw picture
VIM_2b <- cbind(as.matrix(row.names(VIM_2b), nrow(VIM_2b), 1), VIM_2b)
colnames(VIM_2b) <- c('variable', 'DR', 'Lower_Bound', 'Upper_Bound')
VIM_2b <- VIM_2b[order(VIM_2b$DR, decreasing = FALSE),]
order_2b <- VIM_2b$variable

fig_2b <- ggplot(VIM_2b, aes(x = factor(variable, level = order_2b), y = DR))+
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower_Bound, ymax = Upper_Bound), width=.1, position=position_dodge(0.2))+
  theme(axis.title.y=element_blank(), 
        axis.text.y = element_text(color = "grey20", size = 10),
        axis.text.x = element_text(color = "grey20", size = 10)) +
  coord_flip() +
  ylim(0, 1)
fig_2b

