# Run from the package root. Synthetic data only; no nuisance fits or PD integrals.
pkgload::load_all(".",quiet=TRUE)
set.seed(20261008)
x <- runif(8000,-1,1)
y <- sin(3*x)+rnorm(length(x),sd=.5)
cat("Single debiased candidate: exhaustive versus 5 subsets of 250\n")
for(n in c(1000L,2000L,4000L,8000L)) {
  i <- seq_len(n)
  full <- system.time(.local.cv(x[i],y[i],.5,.5,TRUE,"gau",NULL,NULL,5))[['elapsed']]
  set.seed(123)
  sub <- system.time(.local.cv(x[i],y[i],.5,.5,TRUE,"gau",NULL,250,5))[['elapsed']]
  print(data.frame(n=n,exhaustive=full,subset=sub,speedup=full/sub),row.names=FALSE)
}
cat("\nBandwidth and curve stability on n=1000, relative to exhaustive CV\n")
i <- seq_len(1000)
g <- expand.grid(h=c(.2,.4,.8),b=c(.2,.4,.8))
full <- .local.cv(x[i],y[i],g$h,g$b,TRUE,"gau",NULL,NULL,5)
grid <- seq(-.9,.9,length.out=31)
curve <- function(k) .lprobust(x[i],y[i],g$h[k],g$b[k],TRUE,grid,"gau")[,2]
reference <- curve(full$selected)
print(g[full$selected,,drop=FALSE])
for(m in c(100L,250L,500L)) for(seed in 1:3) {
  set.seed(seed)
  elapsed <- system.time(z <- .local.cv(x[i],y[i],g$h,g$b,TRUE,"gau",NULL,m,5))[['elapsed']]
  delta <- curve(z$selected)-reference
  print(data.frame(size=m,seed=seed,h=g$h[z$selected],b=g$b[z$selected],
                   seconds=elapsed,curve.rmse=sqrt(mean(delta^2))),row.names=FALSE)
}
