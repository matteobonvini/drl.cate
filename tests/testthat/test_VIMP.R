context("test VIMP code")
lm.mod <- function(y, x, new.x) {
  fit <- lm(y ~ ., data = cbind(data.frame(y = y), x))
  out <- cbind(predict(fit, newdata = new.x), NA, NA)
  return(list(res = out))
}
mu1.x <- function(y, a, x, new.x) {
  preds <- lm.mod(y[a == 1], x[a == 1, , drop = FALSE], new.x)
  return(list(res = preds$res[, 1]))
}
mu0.x <- function(y, a, x, new.x) {
  preds <- lm.mod(y[a == 0], x[a == 0, , drop = FALSE], new.x)
  return(list(res = preds$res[, 1]))
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
    y <- x3 + a * (x1 + x2) + rnorm(n, sd = 0.1)
    data <- data.frame(y = y, a = a, x1 = x1, x2 = x2, x3 = x3)
    cate.fit <- cate(data_frame = data, learner = "dr",
                     x_names = c("x1", "x2", "x3"),
                     y_name = "y",
                     a_name = "a",
                     v_names = c("x1", "x2", "x3"),
                     num_grid = 100,
                     v0.long = NULL,
                     nsplits = nsplits,
                     univariate_reg = FALSE,
                     partial_dependence = FALSE,
                     additive_approx = FALSE,
                     mu1.x = mu1.x,
                     mu0.x = mu0.x,
                     pi.x = pi.x,
                     drl.v = lm.mod,
                     drl.x = lm.mod)

    var.names <- list("x1", "x2", c("x1", "x2"), "x3")
    VIM_3b <- vimp(cate.fit, var.names, toupper(var.names))
    expect_true(which.max(VIM_3b$DR) == 3)
    expect_true(which.min(VIM_3b$DR) == 4)
    if(nsplits == 1) expect_true(all(0 <= VIM_3b$DR & VIM_3b$DR <= 1))
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
    cate.fit <- cate(data_frame = data, learner = "dr",
                     x_names = c("x1", "x2", "x3"),
                     y_name = "y",
                     a_name = "a",
                     v_names = c("x1", "x2", "x3"),
                     num_grid = 100,
                     v0.long = NULL,
                     nsplits = nsplits,
                     univariate_reg = FALSE,
                     partial_dependence = FALSE,
                     additive_approx = FALSE,
                     mu1.x = mu1.x,
                     mu0.x = mu0.x,
                     pi.x = pi.x,
                     drl.v = lm.mod,
                     drl.x = lm.mod)

    var.names <- list("x2", "x3")
    VIM_3b <- vimp(cate.fit, var.names, toupper(var.names))
    expect_true(all(abs(VIM_3b$DR) < 0.001))
    if(nsplits == 1) expect_true(all(0 <= VIM_3b$DR & VIM_3b$DR <= 1))
  }
})


