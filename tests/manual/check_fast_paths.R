# Checks for the parallel loops (fork addition, rebased on 90ad039).
#
#   Rscript tests/manual/check_fast_paths.R [cores] [timing-n]
#
# 1. Parallel and serial leave-one-out agree exactly, on every row and on a
#    subsample of evaluation rows (the cv.eval.size path).
# 2. Parallel and serial PD integration agree exactly.
# 3. Neither touches the caller's RNG state.
# 4. .local.cv(): parallel == serial, exhaustive and subsampled (same seed).
# 5. debiased_inference(): parallel == serial, exhaustive and subsampled.
# 6. Timing: exact LOO serial and parallel, subsampled LOO parallel, at timing-n.
library(drl.cate)
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 8L
n.time <- if (length(args) >= 2) as.integer(args[2]) else 20000L
ns <- asNamespace("drl.cate")
loo <- get(".local.loo", ns); lcv <- get(".local.cv", ns); integ <- get(".pd.integrate", ns)
ok <- function(cond, what) { cat(if (cond) "PASS " else "FAIL ", what, "\n"); if (!cond) quit(status = 1) }

# Residual-risk-like modifier: mass near 0, long right tail, nearly all distinct,
# plus a block of ties so the shared-system path for repeated values is used.
make <- function(n, seed) {
  set.seed(seed)
  x <- c(rnorm(round(.9 * n), 0, .05), rexp(n - round(.9 * n), 8))
  x <- round(x, 6); x[sample.int(n, round(n / 20))] <- 0.01
  y <- -0.1 * x + 0.3 * x^2 + rnorm(n, 0, 1)
  list(x = x, y = y)
}
serial <- function(expr) { options(drl.cate.cores = 1L); force(expr) }
par <- function(expr) { options(drl.cate.cores = cores); force(expr) }

d <- make(3000, 1)
set.seed(11); sub <- sample.int(length(d$x), 700)
for (deb in c(FALSE, TRUE)) {
  s <- serial(loo(d$x, d$y, .05, .1, deb, "gau"))
  set.seed(99); before <- .Random.seed
  p <- par(loo(d$x, d$y, .05, .1, deb, "gau"))
  ok(identical(s, p), sprintf("1. LOO parallel == serial, all rows (debias=%s)", deb))
  ok(identical(before, .Random.seed), sprintf("3. LOO leaves RNG state (debias=%s)", deb))
  s <- serial(loo(d$x, d$y, .05, .1, deb, "gau", eval.indices = sub))
  p <- par(loo(d$x, d$y, .05, .1, deb, "gau", eval.indices = sub))
  ok(identical(s, p), sprintf("1. LOO parallel == serial, 700 sampled rows (debias=%s)", deb))
  full <- serial(loo(d$x, d$y, .05, .1, deb, "gau"))
  ok(identical(s$pred, full$pred[sub]), sprintf("1. sampled LOO == full LOO at those rows (debias=%s)", deb))
}

# PD integration with a stand-in nuisance fit.
set.seed(2); n <- 2500
src <- structure(list(splits.id = sample(rep(1:4, length.out = n)),
                      fits = lapply(1:4, function(k) list(fit = function(w) k * 0.01 + 0.1 * w[, 1] - 0.2 * w[, 2]^2)),
                      v1 = runif(n), v2 = data.frame(age = rnorm(n)), block.size = 300),
                 class = "pd_integration_source")
W <- matrix(rnorm(n * 6), n, 6); mh <- rnorm(n, 0, .01)
s <- serial(integ(src, W, mh))
set.seed(5); before <- .Random.seed
p <- par(integ(src, W, mh))
ok(identical(s, p), "2. PD integral parallel == serial")
ok(identical(before, .Random.seed), "3. PD integral leaves RNG state")

# .local.cv(): candidates as in debiased_inference(debias=TRUE).
h <- rep(c(.025, .05, .1), 3); b <- rep(c(.05, .1, .25), each = 3)
s <- serial(lcv(d$x, d$y, h, b, TRUE, "gau", NULL, NULL, 5))
p <- par(lcv(d$x, d$y, h, b, TRUE, "gau", NULL, NULL, 5))
ok(identical(s, p), "4. .local.cv exhaustive: parallel == serial")
set.seed(3); s <- serial(lcv(d$x, d$y, h, b, TRUE, "gau", NULL, 500, 5)); s.after <- .Random.seed
set.seed(3); p <- par(lcv(d$x, d$y, h, b, TRUE, "gau", NULL, 500, 5)); p.after <- .Random.seed
ok(identical(s, p), "4. .local.cv subsampled (seeded): parallel == serial")
ok(identical(s.after, p.after), "4. .local.cv subsampled: same RNG state afterwards")

# debiased_inference(), as cate() calls it for a univariate curve.
di <- function(deb, size) debiased_inference(A = d$x, pseudo.out = d$y, debias = deb,
  bandwidth.method = "LOOCV", kernel.type = "gau", bw.seq = c(.025, .05, .1, .25),
  cv.eval.size = size, cv.repeats = 5)
for (deb in c(FALSE, TRUE)) {
  set.seed(4); s <- serial(di(deb, NULL)); set.seed(4); p <- par(di(deb, NULL))
  ok(identical(s, p), sprintf("5. debiased_inference exhaustive: parallel == serial (debias=%s)", deb))
  set.seed(4); s <- serial(di(deb, 500)); set.seed(4); p <- par(di(deb, 500))
  ok(identical(s, p), sprintf("5. debiased_inference subsampled: parallel == serial (debias=%s)", deb))
}

# Timing.
d <- make(n.time, 3)
cat(sprintf("\n6. timing at n = %d (%d distinct), one debiased candidate (h = 0.05, b = 0.1)\n",
            n.time, length(unique(d$x))))
tm <- function(expr) round(unname(system.time(expr)["elapsed"]), 1)
t.loo1 <- tm(serial(loo(d$x, d$y, .05, .1, TRUE, "gau")))
t.looP <- tm(par(loo(d$x, d$y, .05, .1, TRUE, "gau")))
set.seed(1)
t.subP <- tm(par(lcv(d$x, d$y, .05, .1, TRUE, "gau", NULL, 5000, 5)))
cat(sprintf("   exact LOO, serial            : %7.1f s\n", t.loo1))
cat(sprintf("   exact LOO, %2d cores          : %7.1f s  (x%.1f)\n", cores, t.looP, t.loo1 / t.looP))
cat(sprintf("   5 x 5000 sampled, %2d cores   : %7.1f s  (x%.1f)\n", cores, t.subP, t.loo1 / t.subP))
cat("\nall checks passed\n")
