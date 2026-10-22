## Why L-HDCQR here is solved at level 0.5 of the augmented problem.
## The authors' routine quantproc() (reference/zheng2018_author_quantproc.R)
## calls rq(augy ~ 0 + augx, method = "lasso", tau = JJ[s]) at every grid
## point, but the augmented form of the Peng-Huang estimating equation is an
## absolute-loss (tau = 0.5) problem.  Without censoring the correct path must
## reproduce ordinary quantile regression; the author routine does not.
## Rscript reference/check_author_quantproc.R   (run from the bundle folder)
source("load.R")
suppressPackageStartupMessages(library(quantreg))
if (!file.exists(file.path(CQR_COMP_ROOT, "reference", "zheng2018_author_quantproc.R")))
  stop("the L-HDCQR author routine is not in the public repository; it is in the Supplementary Material of the paper")
author <- new.env(); author$rq <- quantreg::rq; author$H <- H
sys.source(file.path(CQR_COMP_ROOT, "reference", "zheng2018_author_quantproc.R"), envir = author)
set.seed(11)
n <- 3000; X <- matrix(rnorm(n * 3), n)
y <- drop(1 + X %*% c(1, -1, 0.5) + rnorm(n))
Z <- cbind(1, X); JJ <- seq(0.1, 0.5, length.out = 21)
beta0 <- coef(rq(y ~ X, tau = JJ[1]))
a_path <- author$quantproc(y, Z, rep(1, n), JJ, beta0, lambda = 0, tol = 0)
o_path <- hdcqr_path(y, rep(1, n), Z, JJ, 0)$B
for (t in c(6, 11, 16, 21))
  cat(sprintf("tau = %.2f | rq intercept %.3f | this solver %.3f | author quantproc %.3f\n",
              JJ[t], coef(rq(y ~ X, tau = JJ[t]))[1], o_path[1, t], a_path[t, 1]))
