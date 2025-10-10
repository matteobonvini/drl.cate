context("test VIMP code")
set.seed(100)
drl.v <- function(pseudo, v, new.v) {
  fit <- lm(y ~ ., data = cbind(data.frame(y = pseudo), v))
  out <- cbind(predict(fit, newdata = new.v), NA, NA)
  return(list(res = out))
}
drl.x <- function(pseudo, x, new.x) {
  fit <- lm(y ~ ., data = cbind(data.frame(y = pseudo), x))
  out <- cbind(predict(fit, newdata = new.x), NA, NA)
  return(list(res = out))
}
mu1.x <- function(y, a, x, new.x) {
  fit <- lm(y ~ ., data = cbind(data.frame(y =y[a == 1]), x[a == 1, , drop=FALSE]))
  out <- cbind(predict(fit, newdata = new.x), NA, NA)
  return(list(res = out[, 1]))
}
mu0.x <- function(y, a, x, new.x) {
  fit <- lm(y ~ ., data = cbind(data.frame(y =y[a == 0]), x[a == 0, , drop=FALSE]))
  out <- cbind(predict(fit, newdata = new.x), NA, NA)
  return(list(res = out[, 1]))
}
pi.x <- function(a, x, new.x){
  fit <- glm(a ~ ., data = cbind(data.frame(a = a), x), family = binomial())
  out <- cbind(predict(fit, newdata = new.x, type = "response"), NA, NA)
  return(list(res = out[, 1]))
}
test_that("expected behavior with pure noise and signal", {
  for(i in 1:10) {
    n <- 1000
    nsplits <- sample(1:3, 1)
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    x3 <- rnorm(n)
    a <- rbinom(n, 1, 0.5)
    y <- x3 + a * (x1 + x2) + rnorm(n, sd = 1)
    data <- data.frame(y = y, a = a, x1 = x1, x2 = x2, x3 = x3)
    v0 <-expand.grid(seq(min(x1), max(x1), length.out = 25),
                     seq(min(x2), max(x2), length.out = 25),
                     seq(min(x3), max(x3), length.out = 25))
    colnames(v0) <- c("x1", "x2", "x3")
    cate.fit <- cate(data = data, learner="dr",
                     x_names = c("x1", "x2", "x3"),
                     y_name = "y",
                     a_name = "a",
                     v_names = c("x1", "x2", "x3"),
                     v0 = v0,
                     nsplits = nsplits,
                     univariate_reg = FALSE,
                     partial_dependence = FALSE,
                     additive_approx = FALSE,
                     mu1.x = mu1.x,
                     mu0.x = mu0.x,
                     pi.x = pi.x,
                     drl.v = drl.v,
                     drl.x = drl.x)

    var.names <- list("x1", "x2", c("x1", "x2"), "x3")
    lab.var.names <- c("x1", "x2", "x1 and x2", "x3")
    VIM_3b <- get_vimp(cate.fit, var.names, lab.var.names)
    expect_true(which.max(VIM_3b$psi) == 3)
    expect_true(which.min(VIM_3b$psi) == 4)
    if(nsplits == 1) expect_true(all(0 <= VIM_3b$psi & VIM_3b$psi <= 1))
  }
})


test_that("expected behavior with only pure noise", {
  for(i in 1:10) {
    n <- 1000
    nsplits <- sample(1:3, 1)
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    x3 <- rnorm(n)
    a <- rbinom(n, 1, 0.5)
    y <- a * x1 + rnorm(n, sd = 0.1)
    data <- data.frame(y = y, a = a, x1 = x1, x2 = x2, x3 = x3)
    v0 <- expand.grid(seq(min(x1), max(x1), length.out=25),
                      seq(min(x2), max(x2), length.out = 25),
                      seq(min(x3), max(x3), length.out = 25))
    colnames(v0) <- c("x1", "x2", "x3")
    cate.fit <- cate(data = data, learner = "dr",
                     x_names = c("x1", "x2", "x3"),
                     y_name = "y",
                     a_name = "a",
                     v_names = c("x1", "x2", "x3"),
                     v0=v0,
                     nsplits = nsplits,
                     univariate_reg = FALSE,
                     partial_dependence = FALSE,
                     additive_approx = FALSE,
                     mu1.x = mu1.x,
                     mu0.x = mu0.x,
                     pi.x = pi.x,
                     drl.v = drl.v,
                     drl.x = drl.x)

    var.names <- list("x2", "x3")
    VIM_3b <- get_vimp(cate.fit, var.names, toupper(var.names))
    expect_true(all(abs(VIM_3b$psi) < 0.001))
    if(nsplits == 1) expect_true(all(0 <= VIM_3b$psi & VIM_3b$psi <= 1))
  }
})


