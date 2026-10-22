## Usage: Rscript sim/calibrate.R
## Computes the censoring location for every scenario (200,000 draws each)
## and two diagnostics at the target quantile: the at-risk probability
## P(C >= Q_T(tau|Z)) over Z, and the share of 20 oracle Peng-Huang refits
## (half sample, true active columns) whose path reaches tau.  Run once
## before the job array; writes sim/calibration.csv.

SIM_ROOT <- normalizePath(Sys.getenv("CQR_SIM_ROOT", "."))
source(file.path(SIM_ROOT, "R", "utils.R")); load_sources(SIM_ROOT)
source(file.path(SIM_ROOT, "sim", "config.R"))
cfg <- SIM_CONFIG
sc <- make_scenarios()
out <- do.call(rbind, lapply(seq_len(nrow(sc)), function(i) {
  spec <- do.call(default_spec, as.list(sc[i, -1L]))
  cal <- calibrate_censoring(spec, seed = make_seed(cfg$master_seed, sc$scenario[i], 0L, "calibration"))
  taus <- grid_taus(spec$tau, cfg); act <- true_active(spec)
  reach <- with_seed(make_seed(cfg$master_seed, sc$scenario[i], 0L, "oracle"), mean(replicate(20, {
    sp <- spec; sp$n <- spec$n %/% 2L; d <- sim_data(sp, cal$c0)
    hdcqr_path(d$y, d$delta, cbind(1, d$X[, act, drop = FALSE]), taus, 0)$status == "ok"
  })))
  data.frame(scenario = sc$scenario[i], c0 = cal$c0, achieved = cal$achieved,
             at_risk_mean = cal$at_risk_mean, at_risk_q05 = cal$at_risk_q05,
             oracle_refit_reaches_tau = reach, s_true = length(act))
}))
write.csv(out, calibration_file(SIM_ROOT), row.names = FALSE)
print(out, digits = 4)
if (any(out$oracle_refit_reaches_tau < 1))
  message("WARNING: even the oracle refit fails to reach tau in: ",
          paste(out$scenario[out$oracle_refit_reaches_tau < 1], collapse = ", "))
