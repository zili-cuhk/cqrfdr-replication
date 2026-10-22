## Usage: Rscript tests/check_author_lhdcqr.R
## Demonstrates the problem in the L-HDCQR author routine `quantproc` used by
## the draft's comparator: every grid step calls
##     rq(augy ~ 0 + augx, method = "lasso", tau = JJ[s], ...)
## but the augmented representation of the Peng-Huang estimating equation is
## an absolute-loss (tau = 0.5) problem.  With tau = JJ[s] the first-order
## condition becomes sum_i Z_i {d_i 1(y_i < Z_i'b) - 2 tau_s * Lambda_i} = 0,
## i.e. the compensator is multiplied by 2*tau_s.  Without censoring the
## correct path must reproduce ordinary quantile regression.
SIM_ROOT <- normalizePath(Sys.getenv("CQR_SIM_ROOT", "."))
source(file.path(SIM_ROOT, "R", "utils.R")); load_sources(SIM_ROOT)
suppressPackageStartupMessages(library(quantreg))
if (!file.exists(file.path(SIM_ROOT, "tests", "ref", "hdcqr_author_fns.R")))
  stop("the L-HDCQR author routine is not in the public repository; it is in the Supplementary Material of the paper")
author <- new.env(); author$rq <- quantreg::rq
sys.source(file.path(SIM_ROOT, "tests", "ref", "hdcqr_author_fns.R"), envir = author)

set.seed(11)
n <- 3000; X <- matrix(rnorm(n * 3), n)
y <- drop(1 + X %*% c(1, -1, 0.5) + rnorm(n))      # no censoring
Z <- cbind(1, X); JJ <- seq(0.1, 0.5, length.out = 21)
beta0 <- coef(rq(y ~ X, tau = JJ[1]))                # author code needs beta(tau_0)
a_path <- author$quantproc(y, Z, rep(1, n), JJ, beta0, lambda = 0, tol = 0)
o_path <- hdcqr_path(y, rep(1, n), Z, JJ, 0)$B
for (t in c(6, 11, 16, 21)) {
  truth <- coef(rq(y ~ X, tau = JJ[t]))
  cat(sprintf("tau = %.2f | rq intercept %.3f | corrected path %.3f | author quantproc %.3f\n",
              JJ[t], truth[1], o_path[1, t], a_path[t, 1]))
}
