context("Robinson partially linear")

test_that("check we recover true beta if partially linear model is true", {
  cate.not.j <- function(y, x, new.x) {
    dat <- cbind(data.frame(y=y), x)
    fit <- lm(y~v2+v2^2, data=dat)
    new.x <- predict(fit, newdata=new.x)
    return(new.x)

  }
  reg.basis.not.j <- function(y, x, new.x) {
    dat <- cbind(data.frame(y=y), x)
    fit <- lm(y~-1+v2, data=dat)
    new.x <- predict(fit, newdata=new.x)
    return(new.x)

  }
  n <- 2000
  nsim <- 5000
  beta <- rep(NA, nsim)
  for(i in 1:nsim){
    v2 <- rnorm(n)
    v1 <- -2*v2 + rnorm(n, sd=0.1)
    y <- 2 + 5*v1 + v2^2 + rnorm(n, sd=0.1)

    fit <- robinson(pseudo=y, w=data.frame(v2=v2), v=v1, new.v=v1, s=rep(1, n),
                    cate.not.j, reg.basis.not.j, dfs=1)
    beta[i] <- coef(fit$model)
    print(i)
  }
  expect_true(abs(mean(beta) - (5)) < 1e-2)


  #### Sim setting 2
  cate.not.j <- function(y, x, new.x) {
    dat <- cbind(data.frame(y=y), x)
    fit <- lm(y~I(v2>0) + poly(v2, 4), data=dat)
    new.x <- predict(fit, newdata=new.x)
    return(new.x)

  }
  reg.basis.not.j <- function(y, x, new.x) {
    dat <- cbind(data.frame(y=y), x)
    fit <- lm(y~poly(v2, 4), data=dat)
    new.x <- predict(fit, newdata=new.x)
    return(new.x)

  }

  nsim <- 100
  beta1 <- beta2 <- rep(NA, nsim)
  n <- 5000
  for(i in 1:nsim){
    v2 <- rnorm(n)
    v1 <- v2+v2^2 + rnorm(n, sd=1)
    y <- 2 + 5*v1 - 7*v1^2 + v2 + v2^2 + 5*I(v2 > 0) + 2*v2^3 + rnorm(n, sd=0.1)
    res.y <- y - cate.not.j(y, data.frame(v2=v2), data.frame(v2=v2))
    res.v1 <- v1 - reg.basis.not.j(v1, data.frame(v2=v2), data.frame(v2=v2))
    res.v2 <- v1^2 - reg.basis.not.j(v1^2, data.frame(v2=v2), data.frame(v2=v2))

    fit <- robinson(pseudo=y, w=data.frame(v2=v2), v=v1, new.v=v1, s=rep(1, n),
                    cate.not.j, reg.basis.not.j, dfs=2)
    beta1[i] <- coef(fit$model)[1]
    beta2[i] <- coef(fit$model)[2]
    expect_true(max(abs(coef(lm(res.y~-1+res.v1+res.v2))-coef(fit$model))) < 1e-10)
    print(i)
  }

  expect_true(abs(mean(beta1) - (5)) < 1e-2)
  expect_true(abs(mean(beta2) - (-7)) < 1e-2)


})
