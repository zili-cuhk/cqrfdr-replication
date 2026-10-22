## ---------------------------------------------------------------------------
## Censored quantile regression: variable selection methods (stand-alone)
## ---------------------------------------------------------------------------
## Usage:  source("path/to/cqr_comparators/load.R")
## Loads the solver and all methods into the global environment and defines
## one wrapper per method with the exact settings used in the paper
## (Yu and Liu, FDR-controlled variable selection for high-dimensional
## censored quantile regression).  All wrappers take
##   y     : log of the observed time, log(min(T, C))
##   delta : event indicator (1 = event, 0 = censored)
##   X     : n x p numeric covariate matrix (standardized internally)
##   taus  : increasing quantile grid tau_0 < ... < tau (target = last point),
##           e.g. seq(0.1, 0.5, length.out = 40)
## and return a list with status, selected (column indices of X) and failed.
## FDR-controlling methods take q (one value or a vector; with a vector the
## result has by_q, one selection per level).

CQR_COMP_ROOT <- local({
  f <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (is.null(f)) normalizePath(Sys.getenv("CQR_COMP_ROOT", ".")) else dirname(f)
})
for (f in c("dgp.R", "solver.R", "methods.R", "utils.R"))
  sys.source(file.path(CQR_COMP_ROOT, "R", f), envir = globalenv())

## --- proposed methods ---------------------------------------------------------
## SS-MS: one half-half split; lasso-penalized Peng-Huang fit on split 1 with the
## pivotal penalty scaled by a multiplier chosen by 5-fold CV; path-ordered
## screen capped at floor(d2/log d2); unpenalized refit on split 2; mirror
## statistics with the cutoff U+ (offset 1).
cqr_ssms <- function(y, delta, X, taus, q = c(0.1, 0.2), K = 5L, cores = 1L, foldid = NULL)
  ssms(y, delta, X, taus, q = q, offset = 1, lambda_rule = "cv", K = K, loss = "ph",
       screen_cap = "theory", screen_rule = "path", cores = cores, foldid = foldid)

## MS-MS: SS-MS on B random splits, multiplier tuned once, aggregated by the
## inclusion rates of Dai et al. (2023), ties broken at random.
cqr_msms <- function(y, delta, X, taus, q = c(0.1, 0.2), B = 50L, K = 5L, cores = 1L, foldid = NULL)
  mds(y, delta, X, taus, q = q, B = B, offset = 1, tune = "cv_once", ties = "random",
      K = K, screen_cap = "theory", screen_rule = "path", cores = cores, foldid = foldid)

## --- comparators, each following its original publication ----------------------
## DS (Dai, Lin, Xing and Liu, 2023, JASA): half-half split; lasso on split 1
## tuned by 10-fold CV within that half; candidate set = CV support (no cap);
## unpenalized Peng-Huang refit on split 2 (in place of their OLS refit);
## cutoff U: #{M < -t} / max(#{M > t}, 1) <= q.
cqr_ds <- function(y, delta, X, taus, q = c(0.1, 0.2), cores = 1L)
  ds(y, delta, X, taus, q = q, K = 10L, loss = "ph", cores = cores)

## MDS (Dai et al., 2023): B complete DS fits, each with its own 10-fold CV,
## aggregated by inclusion rates with their selection rule {j: I_j > I_(l)}.
cqr_mds <- function(y, delta, X, taus, q = c(0.1, 0.2), B = 50L, cores = 1L)
  mds(y, delta, X, taus, q = q, B = B, offset = 0, tune = "cv_each", ties = "dai",
      K = 10L, screen_cap = "none", screen_rule = "cv", foldid = NULL, cores = cores)

## L-HDCQR (Zheng, Peng and He, 2018, AoS, Section 5.1): sequential
## L1-penalized estimation on the grid; lambda_0 pivotal at tau_0, then at every
## grid point an increment chosen by 5-fold CV; each step solved as an
## augmented median (LAD) problem.  Selection = nonzero coefficients at tau.
cqr_lhdcqr <- function(y, delta, X, taus, K = 5L, cores = 1L, foldid = NULL)
  lhdcqr_z18(y, delta, X, taus, K = K, foldid = foldid, cores = cores)

## SCQR-Lasso / SCQR-SCAD (He, Pan, Tan and Zhou, 2022, AoS): the authors' C++
## code (vendor/hdscqr.cpp, from github.com/XiaoouPan/scqr), 50 log-spaced
## penalty levels in [0.01, 0.2], bandwidth 0.5 (log p / n)^(1/4), their CV.
## Needs Rcpp and RcppArmadillo; the first call compiles the C++ (about 1 min).
cqr_scqr_lasso <- function(y, delta, X, taus, K = 5L, foldid = NULL)
  scqr(y, delta, X, taus, penalty = "lasso", K = K, foldid = foldid)
cqr_scqr_scad <- function(y, delta, X, taus, K = 5L, foldid = NULL)
  scqr(y, delta, X, taus, penalty = "scad", K = K, foldid = foldid)

## Cox-lasso (glmnet), a non-quantile reference: cqr_cox_lasso(exp(y), delta, X)
cqr_cox_lasso <- function(time, delta, X, K = 5L, foldid = NULL) cox_lasso(time, delta, X, K, foldid)

## selected indices for one level q from any result
cqr_selected <- function(res, q = NULL) {
  if (!is.null(q) && !is.null(res$by_q)) {
    k <- which(abs(vapply(res$by_q, `[[`, 0, "q") - q) < 1e-12)
    return(res$by_q[[k]]$selected)
  }
  res$selected
}
invisible(NULL)
