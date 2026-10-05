test_that("evaluation grids are aligned by name and invalid grids fail early", {
  d <- data.frame(v=seq(-1,1,length.out=40),w=seq(10,12,length.out=40),
                  a=rep(0:1,20),y=seq_len(40)/40)
  grid <- data.frame(v=c(-.5,.5),w=c(10.5,11.5))
  zero <- function(y,a,x,new.x) list(res=rep(0,nrow(new.x)))
  stage <- function(pseudo,x,new.x) list(res=cbind(rep(mean(pseudo),nrow(new.x)),NA,NA))
  run <- function(g) {
    set.seed(12)
    cate(d,"dr",c("v","w"),"y","a",c("v","w"),g,
      mu0.x=zero,mu1.x=zero,
      pi.x=function(a,x,new.x) list(res=rep(.5,nrow(new.x))),
      drl.v=function(pseudo,v,new.v) stage(pseudo,v,new.v),drl.x=stage,
      foldid=rep(1:2,20),univariate_reg=TRUE,partial_dependence=FALSE,
      additive_approx=FALSE,partially_linear=FALSE,bw.stage2=list(2,2))
  }
  original <- run(grid)
  expect_false("vimp.df" %in% names(original))
  reversed <- run(grid[c("w","v")])
  # Returned learner closures retain their invocation environments (including g).
  fields <- setdiff(names(original),c("drl.v","drl.x"))
  expect_equal(reversed[fields],original[fields])
  expect_equal(original$univariate.res$dr[[1]]$res$eval.pts,grid$v)
  expect_equal(original$univariate.res$dr[[2]]$res$eval.pts,grid$w)
  input <- function(g) get_input(d,c("v","w"),"y","a",c("v","w"),g)
  expect_equal(input(as.matrix(grid[c("w","v")]))$v0,as.matrix(grid))
  duplicate <- grid; names(duplicate) <- c("v","v")
  for(g in list(grid["v"],cbind(grid,extra=1),duplicate,
                unname(as.matrix(grid)),grid[FALSE,],transform(grid,v=NA_real_)))
    expect_error(input(g),"v0 must be")
})

test_that("the package API excludes VIMP and retains level sets", {
  ns <- asNamespace("drl.cate")
  expect_false(exists("get_vimp", envir=ns, inherits=FALSE))
  expect_false("get_vimp" %in% getNamespaceExports("drl.cate"))
  expect_false("vimp_num_splits" %in% names(formals(cate)))
  expect_true("cate_lvl_set" %in% getNamespaceExports("drl.cate"))
  expect_true(is.function(get("cate_lvl_set", envir=ns)))
})
