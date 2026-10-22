## Usage: Rscript sim/run.R <task_id> [reps_per_task]
## Task t (1-based) runs replications ((t-1) %% chunks)*k + 1 ... of scenario
## ceiling(t / chunks), where k = reps_per_task and chunks = ceiling(reps/k).
## One result file per replication (all methods).  On re-submission a
## finished replication is skipped, except that methods that ended with a
## software error (status R_error / scqr_error) or are missing from the file
## are re-run and merged in; every method has its own seed stream, so the
## re-run gives exactly what an error-free first run would have given.

SIM_ROOT <- normalizePath(Sys.getenv("CQR_SIM_ROOT", "."))
source(file.path(SIM_ROOT, "R", "utils.R")); load_sources(SIM_ROOT)
source(file.path(SIM_ROOT, "sim", "config.R"))
cfg <- SIM_CONFIG
args <- commandArgs(TRUE)
task <- as.integer(args[1]); k <- if (length(args) > 1L) as.integer(args[2]) else 1L
sc <- make_scenarios()
cal <- read.csv(calibration_file(SIM_ROOT), stringsAsFactors = FALSE)
chunks <- ceiling(cfg$reps / k)
s_id <- ceiling(task / chunks)
if (is.na(s_id) || s_id > nrow(sc)) stop("task id out of range: ", task)
reps <- ((task - 1L) %% chunks) * k + seq_len(k)
reps <- reps[reps <= cfg$reps]
name <- sc$scenario[s_id]
spec <- do.call(default_spec, as.list(sc[s_id, -1L]))
c0 <- cal$c0[cal$scenario == name]
if (length(c0) != 1L) stop("no calibration for ", name, "; run sim/calibrate.R first")
taus <- grid_taus(spec$tau, cfg)
for (m in names(cfg)[names(cfg) %in% c("ssms", "ssms_size", "ssms_path", "ds", "msms", "mds", "mds_cv", "mds_dai", "lhdcqr", "lhdcqr_z18")]) cfg[[m]]$cores <- cfg$cores

for (r in reps) {
  path <- file.path(results_dir(SIM_ROOT), name, sprintf("rep%03d.rds", r))
  old <- if (file.exists(path)) readRDS(path) else NULL
  todo <- cfg$methods
  if (!is.null(old)) {
    bad <- unique(old$metrics$base_method[old$metrics$status %in% c("R_error", "scqr_error")])
    todo <- union(setdiff(cfg$methods, unique(old$metrics$base_method)), intersect(bad, cfg$methods))
    ## CQR_SIM_ONLY="ds" (comma-separated) restricts a re-run to these methods
    only <- strsplit(Sys.getenv("CQR_SIM_ONLY", ""), ",")[[1]]
    if (length(only) && nzchar(only[1])) todo <- intersect(todo, only)
    if (!length(todo)) { message("skip ", path); next }
    message("re-running ", paste(todo, collapse = ", "), " for ", path)
  }
  dat <- with_seed(make_seed(cfg$master_seed, name, r, "data"), sim_data(spec, c0))
  foldid <- with_seed(make_seed(cfg$master_seed, name, r, "folds"), make_folds(dat$delta, cfg$cv_folds))
  rows <- list(); details <- list()
  for (meth in todo) {
    t0 <- proc.time()[["elapsed"]]
    res <- with_seed(make_seed(cfg$master_seed, name, r, meth),
                     tryCatch(run_method(meth, dat, taus, cfg, foldid),
                              error = function(e) fail("R_error", message = conditionMessage(e))))
    secs <- proc.time()[["elapsed"]] - t0
    rec <- if (!is.null(res$screen)) length(intersect(res$screen, dat$active)) / length(dat$active) else NA
    ## q-dependent methods give one row per nominal level; the others one row (q = NA)
    levels <- if (!is.null(res$by_q)) res$by_q else if (meth %in% c("ssms", "ssms_size", "ssms_path", "ds", "msms", "mds", "mds_cv", "mds_dai"))
      lapply(cfg$q, function(qq) list(q = qq, selected = res$selected, status = res$status))
      else list(list(q = NA_real_, selected = res$selected, status = res$status))
    for (lv in levels) {
      met <- selection_metrics(lv$selected, dat$active, spec$p)
      label <- if (is.na(lv$q)) meth else sprintf("%s_q%s", meth, format(lv$q))
      rows[[label]] <- data.frame(scenario = name, rep = r, method = label, base_method = meth,
                                  q = lv$q, status = lv$status,
                                  failed = isTRUE(res$failed), R = met[["R"]], V = met[["V"]],
                                  T = met[["T"]], FDP = met[["FDP"]], TPP = met[["TPP"]],
                                  screen_size = if (is.null(res$screen)) NA else length(res$screen),
                                  screen_full = if (is.null(res$screen_full_size)) NA else res$screen_full_size,
                                  screen_recall = rec,
                                  mult = if (is.null(res$mult)) NA else res$mult,
                                  n_failed_splits = if (is.null(res$n_failed_splits)) NA else res$n_failed_splits,
                                  seconds = secs, cens_rate = dat$cens_rate,
                                  message = if (is.null(res$message)) "" else res$message,
                                  stringsAsFactors = FALSE)
      message(sprintf("%s rep %d %-14s %-14s R=%d V=%d T=%d  %.0fs", name, r, label, lv$status,
                      met[["R"]], met[["V"]], met[["T"]], secs))
    }
    details[[meth]] <- res
  }
  metrics <- do.call(rbind, rows)
  if (!is.null(old)) {
    metrics <- rbind(old$metrics[!old$metrics$base_method %in% todo, , drop = FALSE], metrics)
    details <- c(old$details[setdiff(names(old$details), todo)], details)
  }
  save_atomic(list(metrics = metrics, details = details, spec = spec,
                   taus = taus, active = dat$active, session = sessionInfo()$R.version$version.string),
              path)
}
