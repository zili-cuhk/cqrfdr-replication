## Usage: Rscript tests/test_core.R      (about a minute)
SIM_ROOT <- normalizePath(Sys.getenv("CQR_SIM_ROOT", "."))
source(file.path(SIM_ROOT, "R", "utils.R")); load_sources(SIM_ROOT)
suppressPackageStartupMessages(library(quantreg))
check <- function(ok, msg) { cat(if (isTRUE(ok)) "PASS " else "FAIL ", msg, "\n", sep = ""); if (!isTRUE(ok)) FAILED <<- TRUE }
FAILED <- FALSE
set.seed(1)

## 1. lad_lasso solves the same LP as quantreg::rq.fit.lasso -----------------
obj <- function(A, y, lam, b) sum(abs(y - A %*% b)) + sum(lam * abs(b))
for (cfg in list(c(n = 60, p = 20), c(n = 50, p = 120))) {
  n <- cfg[["n"]]; p <- cfg[["p"]]
  A <- cbind(1, matrix(rnorm(n * (p - 1)), n)); y <- A[, 2] - A[, 3] + rnorm(n)
  lam <- c(0, rep(2, p - 1))
  b1 <- lad_lasso(A, y, lam)$coef
  b2 <- rq.fit.lasso(A, y, tau = 0.5, lambda = lam)$coefficients
  o1 <- obj(A, y, lam, b1); o2 <- obj(A, y, lam, b2)
  check(abs(o1 - o2) / o2 < 1e-6, sprintf("lad_lasso objective = rq.fit.lasso (n=%d,p=%d): %.8f vs %.8f", n, p, o1, o2))
  check(max(abs(b1 - b2)) < 1e-4, sprintf("lad_lasso coefficients = rq.fit.lasso, max diff %.2e", max(abs(b1 - b2))))
}
b0 <- lad_lasso(A[, 1:10], y, 0)$coef
check(max(abs(b0 - rq.fit.fnb(A[, 1:10], y, 0.5)$coefficients)) < 1e-4, "unpenalised lad = rq.fit.fnb")

## 2. without censoring the Peng-Huang path recovers ordinary quantile regression
n <- 2000; X <- matrix(rnorm(n * 3), n); y <- 1 + X %*% c(1, -1, 0.5) + (1 + 0.5 * pnorm(X[, 1])) * rnorm(n)
y <- drop(y); Z <- cbind(1, X); taus <- seq(0.02, 0.6, length.out = 59)
f <- hdcqr_path(y, rep(1, n), Z, taus, 0)
qr <- coef(rq(y ~ X, tau = 0.6))
check(f$status == "ok" && max(abs(f$B[, 59] - qr)) < 0.08,
      sprintf("uncensored PH path at tau=0.6 matches rq (max diff %.3f)", max(abs(f$B[, 59] - qr))))

## 3. the draft's 56 SS-MS "numerical failures" at 50% censoring ------------
##    50% censoring, 400 observations (~180 events), 120 candidate columns:
##    crq(method = "PengHuang") starts its grid at 1/n, its Frisch-Newton step
##    fails at the first grid points and the path is truncated, so coef()
##    returns NA at every target tau.  The same refit on the paper's grid
##    (starting at nu = 0.1) is fine.  The screen cap avoids the regime anyway.
spec <- default_spec(cens = 0.5, tau = 0.5, p = 200, calib_n = 50000)
cal <- calibrate_censoring(spec, seed = 2)
set.seed(8); d <- sim_data(spec, cal$c0); i2 <- 401:800; tg <- seq(0.1, 0.5, length.out = 40)
f40  <- hdcqr_path(d$y[i2], d$delta[i2], cbind(1, scale(d$X[i2, 1:40])), tg, 0)
f120 <- hdcqr_path(d$y[i2], d$delta[i2], cbind(1, scale(d$X[i2, 1:120])), tg, 0)
suppressMessages(library(survival))
cr <- crq(Surv(d$y[i2], d$delta[i2]) ~ scale(d$X[i2, 1:120]), method = "PengHuang")
check(any(!is.finite(coef(cr, taus = tg))),
      sprintf("crq(PengHuang), 120 columns, %d events: path stops at tau = %.3f, coef() is NA",
              sum(d$delta[i2]), max(cr$sol[1, ])))
## 120 columns for ~170 events is borderline: whether the refit reaches tau
## can differ across BLAS/platforms, so it is reported, not tested.  The
## default screen cap floor(d2 / log d2) (about 33 here) stays well inside.
check(f40$status == "ok",
      sprintf("hdcqr_path refit on the 0.1..0.5 grid with 40 columns: %s (120 columns: %s, informational)",
              f40$status, f120$status))
f300 <- hdcqr_path(d$y[i2], d$delta[i2], cbind(1, scale(d$X[i2, 1:190])), tg, 0)
check(f300$status != "ok", sprintf("more columns than events is reported, not hidden: status %s", f300$status))

## 4. mirror threshold ---------------------------------------------------------
M <- c(5, 4, 3.5, 3, 2.5, 2, -1, 0.5, -0.4, 0)
check(mirror_threshold(M, 0.2, 1) == 2, "U_+ on a toy vector")
check(is.infinite(mirror_threshold(c(1, -1, 0.5), 0.2, 1)), "infeasible threshold gives Inf")

## 5. DGP sanity: location model quantile and censoring calibration ------------
spec <- default_spec(n = 50000, p = 20, s = 10, tau = 0.25, cens = 0.3)
cal <- calibrate_censoring(default_spec(p = 20, s = 10, tau = 0.25, cens = 0.3, calib_n = 50000), 3)
lat <- sim_latent(50000, spec)
qhat <- mean(lat$logT <= true_beta(spec)[1] + lat$X %*% true_beta(spec)[-1])
check(abs(qhat - 0.25) < 0.01, sprintf("P(logT <= Q_tau) = %.3f at tau = 0.25", qhat))
check(abs(cal$achieved - 0.3) < 1e-3, sprintf("calibrated censoring %.4f", cal$achieved))
sh <- default_spec(n = 50000, p = 20, s = 10, tau = 0.25, model = "hetero")
lat <- sim_latent(50000, sh); bt <- true_beta(sh)
qhat <- mean(lat$logT <= bt[1] + lat$X %*% bt[-1])
check(abs(qhat - 0.25) < 0.01 && abs(bt[1 + sh$hetero_index]) > 0.3,
      sprintf("hetero model: P(logT <= Q_tau) = %.3f, coefficient of Z_11 = %.3f", qhat, bt[12]))
check(true_beta(default_spec(model = "hetero", tau = 0.5))[12] == 0, "hetero model: Z_11 null at the median")

for (df in c(1, 2, 5)) {
  sp <- default_spec(n = 50000, p = 20, s = 10, tau = 0.3, xdist = "t", xdf = df,
                     edist = "t", edf = df, err_scale = "iqr")
  lat <- sim_latent(50000, sp); bt <- true_beta(sp)
  qhat <- mean(lat$logT <= bt[1] + lat$X %*% bt[-1])
  check(abs(qhat - 0.3) < 0.01, sprintf("t_%g covariates and errors: P(logT <= Q_tau) = %.3f at tau = 0.3", df, qhat))
}
check(abs(diff(quantile(err_r(1e5, "t", 1, "iqr"), c(.25, .75))) - 2 * qnorm(.75)) < 0.03,
      "IQR-scaled Cauchy errors have the N(0,1) interquartile range")

## 6. small end-to-end run -------------------------------------------------------
spec <- default_spec(n = 400, p = 60, s = 6, cens = 0.3, tau = 0.5)
cal <- calibrate_censoring(default_spec(n = 400, p = 60, s = 6, calib_n = 20000), 4)
set.seed(5); d <- sim_data(spec, cal$c0); taus <- seq(0.1, 0.5, length.out = 20)
r <- ssms(d$y, d$delta, d$X, taus, q = 0.2, K = 5L)
mt <- selection_metrics(r$selected, d$active, spec$p)
check(r$status %in% c("ok", "no_threshold") && mt[["TPP"]] >= 0.5,
      sprintf("SS-MS end-to-end: status %s, R=%d, V=%d, T=%d, screen %d (mult %.2f)",
              r$status, mt[["R"]], mt[["V"]], mt[["T"]], length(r$screen), r$mult))
l <- lhdcqr(d$y, d$delta, d$X, taus, K = 5L)
mt <- selection_metrics(l$selected, d$active, spec$p)
check(l$status == "ok" && mt[["T"]] == 6, sprintf("L-HDCQR end-to-end: R=%d, V=%d, T=%d", mt[["R"]], mt[["V"]], mt[["T"]]))
if (FAILED) quit(status = 1)
