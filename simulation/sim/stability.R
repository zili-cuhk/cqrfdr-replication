## Usage: Rscript sim/stability.R <task_id>
## Split-to-split stability of SS-MS (A2).  For each scenario in STAB below
## and each of STAB_DATA data sets, SS-MS is run with STAB_SPLITS different
## random splits of the SAME data set (each split tuned by the shared 5-fold
## CV, as in the main study), and MDS (corrected version, B = 50) is run twice
## with different split seeds.  Task t handles data set ((t-1) %% STAB_DATA) + 1
## of scenario ceiling(t / STAB_DATA).  One file per data set:
##   results/stability/<scenario>/data###.rds
## Re-submitting is safe: finished data sets are skipped.

SIM_ROOT <- normalizePath(Sys.getenv("CQR_SIM_ROOT", "."))
source(file.path(SIM_ROOT, "R", "utils.R")); load_sources(SIM_ROOT)
Sys.setenv(CQR_SIM_SET = "main")
source(file.path(SIM_ROOT, "sim", "config.R"))
cfg <- SIM_CONFIG

STAB <- list(                                   # scenario name -> spec arguments
  T_p1000_strong = list(p = 1000L),             # baseline of the main study
  T_p1000_weak   = list(p = 1000L, beta_level = 0.3, beta_slope = 0.3),
  T_p100_strong  = list(p = 100L))
STAB_DATA   <- 20L                              # data sets per scenario
STAB_SPLITS <- 50L                              # SS-MS splits per data set
## quick local check: CQR_STAB_TEST=1 shrinks n, p and the number of splits
STAB_TEST <- identical(Sys.getenv("CQR_STAB_TEST"), "1")
## first-stage screen of SS-MS: "cv" (default) or "path"; results of the two
## rules go to different folders
STAB_RULE <- Sys.getenv("CQR_STAB_RULE", "cv")
if (STAB_TEST) { STAB_SPLITS <- 4L; cfg$mds_cv$B <- 3L }

task <- as.integer(commandArgs(TRUE)[1])
s_id <- ceiling(task / STAB_DATA); d <- ((task - 1L) %% STAB_DATA) + 1L
if (is.na(s_id) || s_id > length(STAB)) stop("task id out of range: ", task)
name <- names(STAB)[s_id]
path <- file.path(SIM_ROOT, "results", paste0(if (STAB_TEST) "stability_test" else "stability",
                                               if (STAB_RULE == "cv") "" else paste0("_", STAB_RULE)), name,
                  sprintf("data%03d.rds", d))
if (file.exists(path)) { message("skip ", path); quit(save = "no") }

base <- list(n = 800L, rho = 0.5, tau = 0.5, xdist = "gaussian", edist = "normal",
             cens = 0.30, model = "location", xdf = XDF, edf = EDF)
spec <- do.call(default_spec, utils::modifyList(base, STAB[[name]]))
if (STAB_TEST) { spec$n <- 300L; spec$p <- min(spec$p, 80L) }
c0 <- calibrate_censoring(spec, seed = make_seed(cfg$master_seed, name, 0L, "calibration"))$c0
taus <- grid_taus(spec$tau, cfg)
cores <- cfg$cores
dat <- with_seed(make_seed(cfg$master_seed, name, d, "data"), sim_data(spec, c0))
foldid <- with_seed(make_seed(cfg$master_seed, name, d, "folds"), make_folds(dat$delta, cfg$cv_folds))

splits <- vector("list", STAB_SPLITS)
for (k in seq_len(STAB_SPLITS)) {
  t0 <- proc.time()[["elapsed"]]
  r <- with_seed(make_seed(cfg$master_seed, name, d, sprintf("split%03d", k)),
                 tryCatch(ssms(dat$y, dat$delta, dat$X, taus, q = cfg$q, offset = 1,
                               lambda_rule = "cv", loss = "ph", screen_cap = "theory",
                               K = cfg$cv_folds, foldid = foldid, cores = cores,
                               screen_rule = STAB_RULE),
                          error = function(e) fail("R_error", message = conditionMessage(e))))
  splits[[k]] <- list(status = r$status, failed = isTRUE(r$failed), screen = r$screen,
                      by_q = lapply(if (is.null(r$by_q)) lapply(cfg$q, function(qq) list(q = qq, selected = integer(0)))
                                    else r$by_q, function(l) list(q = l$q, selected = l$selected)),
                      seconds = proc.time()[["elapsed"]] - t0)
  message(sprintf("%s data %d split %d  %s  %.0fs", name, d, k, r$status, splits[[k]]$seconds))
}
mds_runs <- if (STAB_RULE != "cv") NULL else lapply(c("mds_a", "mds_b"), function(role) {
  r <- with_seed(make_seed(cfg$master_seed, name, d, role),
                 tryCatch(do.call(mds, c(list(y = dat$y, delta = dat$delta, X = dat$X, taus = taus,
                                              q = cfg$q, K = cfg$cv_folds, foldid = foldid, cores = cores),
                                         cfg$mds_cv)),
                          error = function(e) fail("R_error", message = conditionMessage(e))))
  list(failed = isTRUE(r$failed), by_q = lapply(r$by_q, function(l) list(q = l$q, selected = l$selected)))
})
save_atomic(list(scenario = name, data = d, spec = spec, active = dat$active, p = spec$p,
                 splits = splits, mds = mds_runs), path)
message("wrote ", path)
