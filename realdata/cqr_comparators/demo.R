## Demo: every method on one small simulated data set (a few minutes).
## Rscript demo.R          (run from this folder)
source("load.R")
set.seed(2026)
spec <- default_spec(n = 300L, p = 100L, s = 10L, rho = 0.5, tau = 0.5, cens = 0.3)
c0 <- calibrate_censoring(utils::modifyList(spec, list(calib_n = 50000L)), seed = 1L)$c0
dat <- sim_data(spec, c0)
taus <- seq(0.1, 0.5, length.out = 15)
cat(sprintf("n = %d, p = %d, censoring = %.2f, active = %s\n\n", spec$n, spec$p,
            dat$cens_rate, paste(dat$active, collapse = " ")))
show <- function(name, res, q = NULL) {
  sel <- cqr_selected(res, q)
  m <- selection_metrics(sel, dat$active, spec$p)
  cat(sprintf("%-18s %-5s status %-12s selected %3d  FDP %.2f  power %.2f\n", name,
              if (is.null(q)) "" else format(q), res$status, length(sel), m[["FDP"]], m[["TPP"]]))
}
run <- function(name, expr, qs = NULL) {
  t0 <- proc.time()[["elapsed"]]; r <- tryCatch(expr, error = function(e) list(status = conditionMessage(e), selected = integer(0)))
  if (is.null(qs)) show(name, r) else for (q in qs) show(name, r, q)
  cat(sprintf("%-18s %.0f s\n", "", proc.time()[["elapsed"]] - t0))
}
run("SS-MS",      cqr_ssms(dat$y, dat$delta, dat$X, taus), c(0.1, 0.2))
run("MS-MS (B=5)", cqr_msms(dat$y, dat$delta, dat$X, taus, B = 5L), c(0.1, 0.2))
run("DS",         cqr_ds(dat$y, dat$delta, dat$X, taus), c(0.1, 0.2))
run("MDS (B=5)", cqr_mds(dat$y, dat$delta, dat$X, taus, B = 5L), c(0.1, 0.2))
run("L-HDCQR",    cqr_lhdcqr(dat$y, dat$delta, dat$X, taus))
run("SCQR-Lasso", cqr_scqr_lasso(dat$y, dat$delta, dat$X, taus))
run("SCQR-SCAD",  cqr_scqr_scad(dat$y, dat$delta, dat$X, taus))
if (requireNamespace("glmnet", quietly = TRUE))
  run("Cox-lasso", cqr_cox_lasso(exp(dat$y), dat$delta, dat$X))
